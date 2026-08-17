#!/usr/bin/env python3
"""Cluster and dissect the Orai × benchmark-ligand docked poses (no crystal).

Companion to ``pose_cluster_crystal_pocket_report.py``, but for the **Orai cross-
docking** run, where there is *no* experimental crystal pose to compare against.
The reference axes are therefore different:

  1. CLUSTERING — for every (Orai receptor frame, ligand) pair, cluster the docked
     poses from all tools (AutoDock Vina, DiffDock, EquiBind) by 3D centroid to find
     the candidate binding sites the tools collectively propose.
  2. CROSS-TOOL HOMOGENEITY — do the tools agree on *where* the ligand binds? For
     each pair we take each tool's own site (its largest cluster's centroid) and
     measure the pairwise inter-tool distances; small = the tools converge. The
     pair's **consensus site** is the geometric median of those per-tool sites
     (one point per tool), NOT the largest cluster's centre — so a high-throughput
     tool can no longer pull the reference toward wherever it sampled most. With
     only one tool present the consensus falls back to that tool's largest cluster
     (flagged ``consensus_source='single_tool'``).
  3. POSEBUSTERS SURVIVORS vs FAILURES — are the poses that pass PoseBusters
     spatially different from those that fail (tighter, nearer the consensus)?
     Evaluated on multi-tool pairs only (where the consensus is throughput-robust).
  4. DISTRIBUTION ACROSS THE ORAI MD SIMULATIONS — the four receptors are MD
     snapshots of the same channel (START-Fr0 + Fr300/400/499). We compare per-frame
     validity, cluster tightness and tool-agreement, and (with a coordinate-frame
     caveat) how much a ligand's consensus site moves between frames.
  5. OUTLIERS — poses far from their pair's consensus, and pairs where the tools
     most disagree.
  6. EXTRA INSIGHTS — binding position along the receptor's principal (≈ pore) axis
     (does a ligand sit at one end of the channel, the pore centre, …?) and a
     per-tool dispersion ("decisiveness") profile (how tight is each tool's own
     cloud of poses).

Pose source: the Orai PoseBusters per-pose results CSV (``docking_method``,
``protein`` = frame, ``ligand``, ``pose_file`` = SDF, plus the PoseBusters test
columns from which ``pb_valid`` is derived). The reader is **method-agnostic** —
it analyses whichever tools are present, so re-running it after the Orai PoseBusters
job finishes (which adds EquiBind) upgrades every panel to the full 3-tool view with
no code change.

Run in the analysis env (conda env ``vina`` — rdkit + sklearn + scipy + matplotlib):

    python Scripts/Analysis/orai_pose_cluster_report.py
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
# Reuse the clustering machinery from the crystal-pose report and the pocket parsers.
from pose_cluster_crystal_pocket_report import (   # noqa: E402
    load_heavy_atom_mol, centroid_from_mol, _extract_centroid,
    SimpleKMedoids, _select_k, cluster_sites, _centroid_dm,
    pockets_from_labels, _dist,
    _internal_indices, _bootstrap_stability,          # reference-free cluster-quality
)
from pose_topn import top_n_allowlist                 # noqa: E402  top-N pose cap

# Optional top-N pose allowlist (pose_file strings). When set (by --top-n-poses),
# load_poses keeps only these poses. Computed once in main() from the FULL pose set,
# so the cap is "top-N first, then the TM / PB-valid filters" for every load path.
_TOPN_ALLOW: Optional[set] = None

try:                                                # shared, unit-tested stats helpers
    import stats_utils as su                        # noqa: E402
except Exception:                                   # pragma: no cover
    su = None
try:
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                   # pragma: no cover
    def _label_panels(axes, fontsize=12):
        flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
        for i, ax in enumerate(flat):
            if ax is not None and getattr(ax, "get_visible", lambda: True)():
                ax.text(-0.08, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                        fontsize=fontsize, fontweight="bold", va="bottom", ha="right")

# tool key -> (pretty name, colour, 3D marker)
TOOL_STYLE = {
    "autodock": ("AutoDock Vina", "#4C72B0", "o"),
    "diffdock": ("DiffDock", "#55A868", "^"),
    "equibind": ("EquiBind", "#C44E52", "s"),
}
_FALLBACK = ("#8172B3", "#CCB974", "#64B5CD", "#937860")

# PoseBusters columns that are metadata or vacuously-failing (no cofactors/waters
# present), excluded when deriving the overall pass verdict — matches the
# validity-report logic in Master_Docking.ipynb.
_PB_META = {
    "file_path", "filepath", "file", "path", "sdf_file", "sdf_path", "method",
    "docking_method", "protein", "ligand", "pose_rank", "rank", "molecule",
    "mol_name", "name", "complex", "protein_path", "ligand_path", "mol_pred",
    "mol_true", "mol_cond", "pose_file", "pose_name", "file_format",
    "protein_file_used", "posebusters_mode", "diffdock_confidence",
}
_PB_EXCLUDE = {
    "mol_true_loaded", "mol_cond_loaded", "not_too_far_away_organic_cofactors",
    "not_too_far_away_inorganic_cofactors", "not_too_far_away_waters",
}
_BOOL_LIKE = frozenset({True, False, 1, 0, 1.0, 0.0, "True", "False", "true", "false"})
_BOOL_MAP = {"True": True, "true": True, "False": False, "false": False}


# ════════════════════════════════════════════════════════════════════════
# IO
# ════════════════════════════════════════════════════════════════════════

def _derive_pb_valid(df: pd.DataFrame) -> pd.Series:
    cols = []
    for c in df.columns:
        cl = c.lower().strip()
        if cl in _PB_META or c in _PB_EXCLUDE or cl in _PB_EXCLUDE:
            continue
        # ``number_*``/``num_*`` are counts; ``most_extreme_*`` are per-atom-pair
        # clash DIAGNOSTICS, not PoseBusters verdict tests. Critically,
        # ``most_extreme_clash_*`` is INVERTED-polarity (True == a clash is
        # present), so AND-ing it into an all-True pass verdict demands a clash
        # to exist — it forced every clash-free pose (all AutoDock Vina poses)
        # to "invalid" and collapsed the per-frame validity panel to a flat 0 %.
        if cl.startswith(("number_", "num_", "most_extreme_")):
            continue
        u = set(df[c].dropna().unique())
        if u and u.issubset(_BOOL_LIKE):
            cols.append(c)
    if not cols:
        return pd.Series(np.nan, index=df.index)
    T = df[cols].apply(lambda s: s.map(_BOOL_MAP) if s.dtype == object else s).astype(bool)
    return T.all(axis=1)


def load_poses(csv: Path, ids: Optional[set],
               dd_variant: Optional[str] = None) -> pd.DataFrame:
    df = pd.read_csv(csv, low_memory=False)
    tool_col = "docking_method" if "docking_method" in df.columns else "method"
    df["tool"] = df[tool_col].astype(str).str.lower()
    # Restrict DiffDock to a single optimizer variant (raw poses carry optimizer
    # "original"; refined poses "smina"/"gnina"). Keeps DiffDock represented by one
    # variant instead of blending raw+smina+gnina into a ~3x pose cloud.
    if dd_variant and dd_variant != "all" and "optimizer" in df.columns:
        _opt = df["optimizer"].astype(str).str.lower()
        _dd = df["tool"].str.startswith("diffdock")
        df = df[(~_dd) | (_opt == dd_variant)].copy()
    # collapse any equibind variant label (equibind_guided / equibind_unguided_*) to "equibind"
    df.loc[df["tool"].str.startswith("equibind"), "tool"] = "equibind"
    df["frame"] = df["protein"].astype(str)
    df["ligand"] = df["ligand"].astype(str)
    df["pb_valid"] = _derive_pb_valid(df)
    keep = ["tool", "frame", "ligand", "pose_file", "pb_valid"]
    df = df[[c for c in keep if c in df.columns]].copy()
    if _TOPN_ALLOW is not None:                       # top-N-poses cap (see main())
        df = df[df["pose_file"].astype(str).isin(_TOPN_ALLOW)].copy()
    if ids is not None:
        df = df[df["ligand"].isin(ids)].copy()
    df = df[df["pose_file"].astype(str).map(lambda p: Path(str(p)).exists())].copy()
    return df


def add_centroids(df: pd.DataFrame, workers: int, max_coord: float = 1000.0) -> pd.DataFrame:
    files = df["pose_file"].astype(str).unique().tolist()
    cents: Dict[str, Optional[np.ndarray]] = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(_extract_centroid, f): f for f in files}
        for fut in as_completed(futs):
            path, cx, cy, cz, _n = fut.result()
            cents[path] = None if cx is None else np.array([cx, cy, cz])
    df = df.copy()
    df["_cent"] = df["pose_file"].astype(str).map(cents)
    df = df[df["_cent"].notna()].copy()
    # Drop physically-impossible centroids (corrupt SDF coordinates, e.g. failed
    # poses written with astronomic values) — they otherwise create phantom
    # singleton clusters and wreck the axis-projection / outlier stats.
    ok = df["_cent"].map(lambda c: bool(np.all(np.isfinite(c)) and np.max(np.abs(c)) <= max_coord))
    n_bad = int((~ok).sum())
    if n_bad:
        print(f"  dropped {n_bad} poses with corrupt centroids (|coord| > {max_coord:g} Å)")
    return df[ok].copy()


# ════════════════════════════════════════════════════════════════════════
# Receptor principal (≈ pore) axis — for the binding-depth analysis
# ════════════════════════════════════════════════════════════════════════

def frame_axis(frame: str, fpocket_dir: Path) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """(origin, unit principal axis) from the frame's CA atoms; None if unavailable."""
    pdb = fpocket_dir / f"{frame}_cleaned.pdb"
    if not pdb.exists():
        return None
    xs = []
    for ln in pdb.read_text(errors="ignore").splitlines():
        if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
            try:
                xs.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            except ValueError:
                continue
    if len(xs) < 10:
        return None
    a = np.asarray(xs)
    origin = a.mean(axis=0)
    _, vecs = np.linalg.eigh(np.cov((a - origin).T))
    return origin, vecs[:, -1]            # largest-variance eigenvector


# ════════════════════════════════════════════════════════════════════════
# Per-(frame, ligand) analysis
# ════════════════════════════════════════════════════════════════════════

def _geometric_median(pts: np.ndarray, n_iter: int = 64, eps: float = 1e-6) -> np.ndarray:
    """Weiszfeld geometric median — a translation-robust centre that, unlike a
    pose-weighted mean, gives every input point equal pull (here: one point per
    tool). For <=2 points it is the mean/midpoint."""
    pts = np.asarray(pts, dtype=float)
    if len(pts) <= 2:
        return pts.mean(axis=0)
    g = pts.mean(axis=0)
    for _ in range(n_iter):
        d = np.linalg.norm(pts - g, axis=1)
        nz = d > eps
        if not nz.any():
            break
        w = 1.0 / d[nz]
        g_new = (pts[nz] * w[:, None]).sum(axis=0) / w.sum()
        if np.linalg.norm(g_new - g) < eps:
            g = g_new
            break
        g = g_new
    return g


def analyze_pair(frame: str, ligand: str, sub: pd.DataFrame, thr: float,
                 axis: Optional[Tuple[np.ndarray, np.ndarray]],
                 site_method: str = "threshold",
                 pocket_radius: float = 8.0) -> Tuple[dict, List[dict]]:
    C = np.vstack(sub["_cent"].to_numpy())
    tools = sub["tool"].to_numpy()
    pbv = sub["pb_valid"].to_numpy()
    files = sub["pose_file"].to_numpy()
    n = len(C)
    if n < 2:
        return {"frame": frame, "ligand": ligand, "n_poses": n, "skipped": True}, []

    # k=1-capable site clustering; largest cluster kept only as a COUNT-concentration
    # diagnostic (dominant_cluster_frac), NOT as the consensus reference.
    dm = _centroid_dm(C)
    labels, k, sil = cluster_sites(dm, C, site_method, pocket_radius)
    pockets = pockets_from_labels(C, labels, list(tools))   # size-ranked
    big = pockets[0]                                   # largest cluster (count only)
    dom_frac = big["size"] / n

    # per-tool site (largest cluster of that tool's poses)
    tool_site: Dict[str, Optional[np.ndarray]] = {}
    tool_spread: Dict[str, float] = {}
    for t in set(tools):
        ti = np.where(tools == t)[0]
        ct = C[ti]
        tool_spread[t] = float(np.linalg.norm(ct - ct.mean(0), axis=1).mean()) if len(ct) > 1 else 0.0
        if len(ct) < 2:
            tool_site[t] = ct.mean(0)
        else:
            tl, _, _ = cluster_sites(_centroid_dm(ct), ct, site_method, pocket_radius)
            lab = Counter(tl.tolist()).most_common(1)[0][0]
            tool_site[t] = ct[np.where(tl == lab)[0]].mean(0)

    inter = {}
    present = sorted(tool_site)
    for i in range(len(present)):
        for j in range(i + 1, len(present)):
            a, b = present[i], present[j]
            inter[f"{a[:2]}_{b[:2]}_dist"] = round(_dist(tool_site[a], tool_site[b]), 3)
    inter_vals = list(inter.values())

    # CONSENSUS SITE = geometric median of the per-tool sites (one point per tool,
    # throughput-invariant). Single-tool pairs fall back to the largest cluster.
    if len(present) >= 2:
        consensus = _geometric_median(np.vstack([tool_site[t] for t in present]))
        consensus_source = "tool_median"
    else:
        consensus = big["center"]
        consensus_source = "single_tool"
    consensus_multitool = consensus_source == "tool_median"

    # dist of each pose to the consensus (for outliers + pb comparison + depth)
    pose_rows = []
    dists = np.linalg.norm(C - consensus, axis=1)
    if axis is not None:
        origin, ax = axis
        depth = (C - origin) @ ax
    else:
        depth = np.full(n, np.nan)
    for i in range(n):
        pose_rows.append({
            "frame": frame, "ligand": ligand, "tool": tools[i],
            "pb_valid": (bool(pbv[i]) if pbv[i] == pbv[i] else np.nan),
            "cluster": int(labels[i]), "dist_to_consensus": round(float(dists[i]), 3),
            "axis_depth": (round(float(depth[i]), 3) if depth[i] == depth[i] else np.nan),
            "is_outlier": bool(dists[i] > thr * 2),
            "consensus_multitool": bool(consensus_multitool),
            "pose_file": files[i],
        })

    # pb-valid spatial: mean dist-to-consensus for valid vs invalid
    valid_mask = np.array([r["pb_valid"] is True for r in pose_rows])
    invalid_mask = np.array([r["pb_valid"] is False for r in pose_rows])
    pb_known = valid_mask | invalid_mask
    rec = {
        "frame": frame, "ligand": ligand, "skipped": False,
        "n_poses": n, "n_tools": len(present), "tools": ",".join(present),
        "n_clusters": len(pockets), "silhouette": round(sil, 3) if sil == sil else np.nan,
        "dominant_cluster_frac": round(dom_frac, 3),   # largest-cluster count share (diagnostic only)
        "consensus_source": consensus_source,
        "consensus_multitool": bool(consensus_multitool),
        "consensus_x": round(float(consensus[0]), 2),
        "consensus_y": round(float(consensus[1]), 2),
        "consensus_z": round(float(consensus[2]), 2),
        "mean_inter_tool_dist": round(float(np.mean(inter_vals)), 3) if inter_vals else np.nan,
        "max_inter_tool_dist": round(float(np.max(inter_vals)), 3) if inter_vals else np.nan,
        "tools_agree": bool(inter_vals and max(inter_vals) <= thr),
        "pb_valid_rate": round(float(valid_mask.sum()) / pb_known.sum(), 3) if pb_known.any() else np.nan,
        "mean_dist_valid": round(float(dists[valid_mask].mean()), 3) if valid_mask.any() else np.nan,
        "mean_dist_invalid": round(float(dists[invalid_mask].mean()), 3) if invalid_mask.any() else np.nan,
        "n_outliers": int(sum(r["is_outlier"] for r in pose_rows)),
        **inter,
        **{f"{t}_spread": round(tool_spread[t], 3) for t in present},
    }
    return rec, pose_rows


# ════════════════════════════════════════════════════════════════════════
# Statistical tests  (see STATISTICAL_VALIDATION_PLAN.md)
# ------------------------------------------------------------------------
# The independent n is the number of (frame, ligand) pairs, which varies by dataset
# (Orai × JKU is TINY — ~12 pairs, a handful of JKU ligands × 4 MD snapshots; the
# Orai × benchmark run has ~1200). Every test below is therefore GUARDED (skips
# gracefully when there are too few units) and clearly labelled EXPLORATORY, and
# we NEVER treat the hundreds of correlated poses as if they were independent:
# proportions/continuous comparisons aggregate to the pair, and the one place we
# look at pose-level distances uses a (frame,ligand) cluster bootstrap so the CI
# respects the pose clustering. Results are annotated on the figures and written
# in full to ``orai_cluster_stats.json``.
# ════════════════════════════════════════════════════════════════════════

_MIN_FRAMES_TREND = 3       # ordered MD frames needed for a Cochran–Armitage trend
_MIN_PAIRS_WILCOXON = 6     # per-pair distances needed for a signed-rank test
_MIN_CLUSTERS_MW = 4        # (frame,ligand) clusters for the pose-level bootstrap CI
_MW_SAMPLE_CAP = 1500       # cap per group for Cliff's-δ CI (its bootstrap is O(n²)/iter;
                            # keeps the always-on test cheap on big pose sets)


def _f(v):
    """nan/None-safe float for JSON sidecars."""
    try:
        if v is None or (isinstance(v, float) and v != v):
            return None
        return float(v)
    except Exception:
        return None


def _stats_pb_frame_trend(poses, frames, tools):
    """Cochran–Armitage trend of PoseBusters validity across the ORDERED MD
    frames, one test per tool. Unit = pose (correlated within a (frame,ligand)
    pair) so this is EXPLORATORY — it matches the pose-level rates the panel
    plots, and with only ~1–3 pairs per frame it cannot be more than suggestive."""
    out = {"test": "cochran_armitage_trend", "ordered_frames": list(frames),
           "unit": "pose (correlated within frame×ligand pair) — exploratory",
           "per_tool": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    if len(frames) < _MIN_FRAMES_TREND:
        out["skipped"] = "n too small — exploratory"
        return out
    d = poses[poses["pb_valid"].notna()]
    for t in tools:
        succ, tot = [], []
        for fr in frames:
            sel = d[(d["frame"] == fr) & (d["tool"] == t)]
            succ.append(int((sel["pb_valid"] == True).sum()))   # noqa: E712
            tot.append(int(len(sel)))
        if sum(tot) == 0:
            continue
        try:
            z, p, sign = su.cochran_armitage(succ, tot)
            out["per_tool"][t] = {"z": _f(z), "p": _f(p), "sign": int(sign),
                                  "successes": succ, "totals": tot,
                                  "star": su.p_stars(p)}
        except Exception as e:                          # pragma: no cover
            out["per_tool"][t] = {"error": str(e)}
    return out


def _stats_tool_pair_agreement(pairs, tools, thr):
    """Per tool-pair: is the per-(frame,ligand) inter-tool DOMINANT-site distance
    significantly BELOW the agreement threshold? One-sample Wilcoxon signed-rank
    of (distance − thr); rank-biserial < 0 ⇒ the tools tend to agree (sit within
    ``thr`` Å). Holm across the tool-pairs. Unit = (frame,ligand) pair (paired
    design). Guarded + exploratory given the tiny n."""
    out = {"test": "wilcoxon_signed_rank_vs_threshold", "thr_A": float(thr),
           "unit": "(frame,ligand) pair", "pairs": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    entries = []
    for i in range(len(tools)):
        for j in range(i + 1, len(tools)):
            a, b = tools[i], tools[j]
            col = f"{a[:2]}_{b[:2]}_dist"
            if col not in pairs.columns:
                col = f"{b[:2]}_{a[:2]}_dist"
            if col not in pairs.columns:
                continue
            v = pd.to_numeric(pairs[col], errors="coerce").dropna().to_numpy()
            key = f"{a}|{b}"
            rec = {"col": col, "n": int(v.size),
                   "median_dist": _f(np.median(v)) if v.size else None,
                   "agree_frac": _f(np.mean(v <= thr)) if v.size else None}
            if v.size >= _MIN_PAIRS_WILCOXON:
                try:
                    rb, p, npair = su.wilcoxon_rankbiserial(v, np.full(v.shape, float(thr)))
                    rec.update({"rank_biserial": _f(rb), "p_raw": _f(p),
                                "n_nonzero": int(npair)})
                    if rec["p_raw"] is not None:
                        entries.append((key, rec))
                except Exception as e:                  # pragma: no cover
                    rec["error"] = str(e)
            else:
                rec["note"] = "n too small — exploratory"
            out["pairs"][key] = rec
    if entries:
        adj = su.holm([rec["p_raw"] for _, rec in entries])
        for (key, rec), pa in zip(entries, adj):
            rec["p_holm"] = float(pa)
            rec["star"] = su.p_stars(pa)
    return out


def _tpa_by_col(stats):
    """Column-name -> per-tool-pair stat record, for figure annotation."""
    if not stats:
        return {}
    return {rec["col"]: rec for rec in stats.get("pairs", {}).values()
            if isinstance(rec, dict) and "col" in rec}


def _stats_valid_vs_consensus(poses):
    """Pose-level, UNPAIRED test: are PoseBusters survivors nearer their pair's
    consensus than failures? Mann–Whitney U + Cliff's δ (with a percentile-
    bootstrap δ CI); the group medians additionally get a CLUSTER bootstrap CI
    with clusters = (frame,ligand) pair, so the CI respects that poses are
    correlated within a pair rather than pretending each pose is independent.
    Multi-tool pairs only (single-tool consensus would be circular)."""
    out = {"test": "mann_whitney_u + cliffs_delta",
           "unit": "pose (clustered by frame×ligand pair) — exploratory"}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    dd = poses[poses["pb_valid"].notna()].copy()
    if "consensus_multitool" in dd:
        dd = dd[dd["consensus_multitool"] == True]      # noqa: E712
    if dd.empty:
        out["skipped"] = "no multi-tool poses"
        return out
    dd["_pair"] = dd["frame"].astype(str) + "|" + dd["ligand"].astype(str)
    vv = dd[dd["pb_valid"] == True]                      # noqa: E712
    iv = dd[dd["pb_valid"] == False]                     # noqa: E712
    vdist = pd.to_numeric(vv["dist_to_consensus"], errors="coerce").to_numpy()
    idist = pd.to_numeric(iv["dist_to_consensus"], errors="coerce").to_numpy()
    n_clusters = int(dd["_pair"].nunique())
    out.update({"n_valid": int(np.isfinite(vdist).sum()),
                "n_invalid": int(np.isfinite(idist).sum()),
                "n_clusters": n_clusters})
    if n_clusters < _MIN_CLUSTERS_MW:
        out["note"] = "n too small — exploratory"
    # Cliff's-δ CI inside mannwhitney_cliffs bootstraps at O(n_a·n_b) per iteration,
    # so cap each group to keep this always-on test cheap on large pose sets. The
    # U-test p-value and δ point estimate are near-invariant to this subsample.
    va = vdist[np.isfinite(vdist)]; ia = idist[np.isfinite(idist)]
    capped = va.size > _MW_SAMPLE_CAP or ia.size > _MW_SAMPLE_CAP
    if capped:
        rng = np.random.default_rng(0)
        if va.size > _MW_SAMPLE_CAP:
            va = rng.choice(va, _MW_SAMPLE_CAP, replace=False)
        if ia.size > _MW_SAMPLE_CAP:
            ia = rng.choice(ia, _MW_SAMPLE_CAP, replace=False)
        out["mw_sample_cap"] = _MW_SAMPLE_CAP
    try:
        mw = su.mannwhitney_cliffs(va, ia)
    except Exception as e:                              # pragma: no cover
        mw = None
        out["error"] = str(e)
    if mw:
        out.update({"U": _f(mw["U"]), "p": _f(mw["p"]),
                    "cliffs_delta": _f(mw["cliffs_delta"]),
                    "delta_ci": [_f(mw["delta_ci"][0]), _f(mw["delta_ci"][1])],
                    "median_valid": _f(mw["medians"][0]),
                    "median_invalid": _f(mw["medians"][1]),
                    "star": su.p_stars(mw["p"])})
    for name, sub in (("valid", vv), ("invalid", iv)):
        try:
            vals = pd.to_numeric(sub["dist_to_consensus"], errors="coerce").to_numpy()
            cl = sub["_pair"].to_numpy()
            est, lo, hi = su.cluster_bootstrap_ci(vals, cl, statistic=np.median)
            out[f"median_{name}_clusterCI"] = [_f(est), _f(lo), _f(hi)]
        except Exception:                               # pragma: no cover
            pass
    return out


def _paired_tool_test(complete):
    """Robust paired-across-tools comparison: Friedman + Kendall's W + Wilcoxon/Holm
    when ≥3 tools are present, a single signed-rank when exactly 2 (Friedman needs
    ≥3 conditions). ``complete`` is a DataFrame whose rows are the same units
    (one per (frame,ligand) pair) and columns are tools, already NaN-dropped."""
    labels = list(complete.columns)
    if len(labels) >= 3:
        return su.paired_continuous(complete, labels=labels)
    a, b = labels
    rb, p, npair = su.wilcoxon_rankbiserial(complete[a].to_numpy(float),
                                            complete[b].to_numpy(float))
    return {"omnibus": None,
            "pairwise": [{"a": a, "b": b, "rank_biserial": _f(rb), "p_raw": _f(p),
                          "p_holm": _f(p), "star": su.p_stars(p), "n": int(npair)}],
            "medians": {l: _f(complete[l].median()) for l in labels}}


def _paired_summary(res):
    """One-line caption for a _paired_tool_test / paired_continuous result."""
    if not res:
        return None
    om = res.get("omnibus")
    if om and om.get("p") is not None:
        return (f"Friedman {su.p_stars(om['p'])} {su.fmt_p(om['p'])}"
                f" (W={om['kendall_w']:.2f}, n={om['n']})")
    pw = res.get("pairwise") or []
    if pw and pw[0].get("p_raw") is not None:
        pr = pw[0]
        return (f"Wilcoxon {su.p_stars(pr.get('p_holm') or pr['p_raw'])}"
                f" {su.fmt_p(pr.get('p_holm') or pr['p_raw'])} (n={pr['n']})")
    return None


def _stats_cluster_structure_trend(pairs, frames):
    """Jonckheere–Terpstra trend of per-(frame,ligand) cluster structure across the
    ORDERED MD frames: does the number of distinct clusters / the dominant-cluster
    fraction rise or fall along the trajectory? Unit = (frame,ligand) pair (one
    value per pair — no pose pseudoreplication). Guarded + exploratory at this n."""
    out = {"test": "jonckheere_terpstra_trend", "ordered_frames": list(frames),
           "unit": "(frame,ligand) pair", "metrics": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    if len(frames) < _MIN_FRAMES_TREND:
        out["skipped"] = "n too small — exploratory"
        return out
    for metric in ("n_clusters", "dominant_cluster_frac"):
        if metric not in pairs.columns:
            continue
        groups = [pd.to_numeric(pairs[pairs.frame == fr][metric], errors="coerce")
                  .dropna().to_numpy() for fr in frames]
        if any(len(g) == 0 for g in groups) or sum(len(g) for g in groups) < _MIN_PAIRS_WILCOXON:
            out["metrics"][metric] = {"note": "n too small — exploratory"}
            continue
        try:
            JT, z, p = su.jonckheere(groups)
            out["metrics"][metric] = {"JT": _f(JT), "z": _f(z), "p": _f(p),
                                      "per_frame_median": [_f(np.median(g)) for g in groups],
                                      "star": su.p_stars(p)}
        except Exception as e:                          # pragma: no cover
            out["metrics"][metric] = {"error": str(e)}
    return out


def _stats_tool_cluster_contribution(poses, tools):
    """Paired-across-tools Friedman + Wilcoxon on the per-(frame,ligand) number of
    distinct clusters each tool occupies. Unit = (frame,ligand) pair, listwise-
    complete across the tools (paired design — no pose pseudoreplication).
    Exploratory at this n."""
    out = {"test": "friedman + wilcoxon (paired across tools)",
           "quantity": "distinct clusters occupied per pair",
           "unit": "(frame,ligand) pair"}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    if "cluster" not in poses.columns:
        out["skipped"] = "no cluster labels"
        return out
    rows = {}
    for (frame, ligand), g in poses.groupby(["frame", "ligand"]):
        rows[f"{frame}|{ligand}"] = {t: float(gt["cluster"].nunique())
                                     for t, gt in g.groupby("tool")}
    mat = pd.DataFrame.from_dict(rows, orient="index").reindex(columns=list(tools))
    complete = mat.dropna()
    out["n_pairs_complete"] = int(len(complete))
    out["medians"] = {t: _f(mat[t].median()) for t in mat.columns if mat[t].notna().any()}
    if len(complete) < _MIN_PAIRS_WILCOXON or complete.shape[1] < 2:
        out["note"] = "n too small — exploratory"
        return out
    try:
        out.update(_paired_tool_test(complete.loc[:, complete.notna().all()]))
    except Exception as e:                              # pragma: no cover
        out["error"] = str(e)
    return out


def _stats_tool_dispersion(pairs, tools):
    """Paired-across-tools Friedman + Wilcoxon on each tool's per-(frame,ligand)
    pose spread (decisiveness). Unit = (frame,ligand) pair, listwise-complete
    across the tools. Exploratory at this n."""
    out = {"test": "friedman + wilcoxon (paired across tools)",
           "quantity": "per-pair pose spread (Å)", "unit": "(frame,ligand) pair"}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    cols = {t: f"{t}_spread" for t in tools if f"{t}_spread" in pairs.columns}
    if len(cols) < 2:
        out["skipped"] = "fewer than 2 tools with a spread column"
        return out
    mat = pd.DataFrame({t: pd.to_numeric(pairs[c], errors="coerce") for t, c in cols.items()})
    complete = mat.dropna()
    out["n_pairs_complete"] = int(len(complete))
    out["medians"] = {t: _f(mat[t].median()) for t in mat.columns if mat[t].notna().any()}
    if len(complete) < _MIN_PAIRS_WILCOXON:
        out["note"] = "n too small — exploratory"
        return out
    try:
        out.update(_paired_tool_test(complete))
    except Exception as e:                              # pragma: no cover
        out["error"] = str(e)
    return out


def _stats_cross_pose_touch(pairs, poses_by_pair, thr):
    """Companion to ``_stats_tool_pair_agreement`` for the divergence-matrix RIGHT
    panel: per tool-pair, is the per-(frame,ligand) CLOSEST cross-tool pose
    distance (do the pose CLOUDS touch?) significantly below the ``thr`` agree
    line? One-sample Wilcoxon signed-rank of (distance − thr); rank-biserial < 0 ⇒
    the clouds tend to graze within ``thr`` Å. Holm across tool-pairs. Unit =
    (frame,ligand) pair. Guarded + exploratory at this n."""
    out = {"test": "wilcoxon_signed_rank_vs_threshold",
           "quantity": "closest cross-tool pose distance",
           "thr_A": float(thr), "unit": "(frame,ligand) pair", "pairs": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    dfp = pairs.sort_values(["frame", "ligand"]).reset_index(drop=True)
    entries = []
    for ta, tb, col, lab in _TOOL_PAIRS:
        vals = []
        for _, r in dfp.iterrows():
            d = _min_cross_pose_dist(poses_by_pair.get((r["frame"], r["ligand"]), []), ta, tb)
            if np.isfinite(d):
                vals.append(d)
        v = np.asarray(vals, float)
        key = f"{ta}|{tb}"
        rec = {"col": col, "label": lab.replace("\n", " "), "n": int(v.size),
               "median_dist": _f(np.median(v)) if v.size else None,
               "touch_frac": _f(np.mean(v <= thr)) if v.size else None}
        if v.size >= _MIN_PAIRS_WILCOXON:
            try:
                rb, p, npair = su.wilcoxon_rankbiserial(v, np.full(v.shape, float(thr)))
                rec.update({"rank_biserial": _f(rb), "p_raw": _f(p), "n_nonzero": int(npair)})
                if rec["p_raw"] is not None:
                    entries.append((key, rec))
            except Exception as e:                      # pragma: no cover
                rec["error"] = str(e)
        else:
            rec["note"] = "n too small — exploratory"
        out["pairs"][key] = rec
    if entries:
        adj = su.holm([rec["p_raw"] for _, rec in entries])
        for (key, rec), pa in zip(entries, adj):
            rec["p_holm"] = float(pa)
            rec["star"] = su.p_stars(pa)
    return out


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _style(tool):
    if tool in TOOL_STYLE:
        return TOOL_STYLE[tool]
    i = abs(hash(tool)) % len(_FALLBACK)
    return (tool, _FALLBACK[i], "o")


# ── individual panel draws ───────────────────────────────────────────────
# Each _draw_* renders one panel into a supplied Axes so the SAME code backs
# both the combined overview grid and the standalone per-panel PNGs.

def _draw_cross_tool_agreement(ax, pairs, thr, tools, stats=None):
    """(A) inter-tool consensus-site distance distribution. Legend labels carry a
    signed-rank star (distance vs the ``thr`` agree-line; see orai_cluster_stats.json)."""
    by_col = _tpa_by_col(stats)
    pair_cols = [c for c in pairs.columns if c.endswith("_dist")
                 and c not in ("mean_inter_tool_dist", "max_inter_tool_dist")]
    # Shared 2 Å bins spanning the full observed range so no pair's distances are
    # clipped off the axis — inter-tool sites can sit far past 40 Å apart (up to
    # ~93 Å on Orai), which a fixed [0,40] window silently dropped, hiding the
    # au_eq/di_eq bars and pushing their medians off the plot. A 40 Å floor keeps
    # tight runs from looking cramped.
    cols_v = [pd.to_numeric(pairs[c], errors="coerce").dropna().to_numpy() for c in pair_cols]
    allv = np.concatenate(cols_v) if cols_v else np.array([])
    top = 40.0 if allv.size == 0 else max(40.0, float(np.ceil(allv.max() / 2.0) * 2.0))
    bins = np.arange(0.0, top + 2.0, 2.0)
    for c, v in zip(pair_cols, cols_v):
        if v.size:
            med = float(np.median(v))
            lab = f"{c.replace('_dist','')} (med {med:.1f} Å)"
            rec = by_col.get(c)
            if rec and rec.get("star"):
                lab += f" {rec['star']}"
            _, _, patches = ax.hist(v, bins=bins, histtype="step", lw=2, label=lab)
            # Put each pair's median literally on the axis (matched colour).
            col = patches[0].get_edgecolor() if patches else None
            ax.axvline(med, color=col, ls=":", lw=1.2, alpha=0.75)
    ax.axvline(thr, color="k", ls="--", lw=0.8, label=f"agree ≤ {thr:g} Å")
    ax.set_title("Cross-tool agreement — distance between tools' consensus sites")
    ax.set_xlabel("Inter-tool consensus-site distance (Å)")
    ax.set_ylabel("Number of receptor-ligand pairs")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)
    if stats is not None and not stats.get("skipped"):
        ax.text(0.98, 0.02, "star: Wilcoxon dist vs threshold (exploratory)",
                transform=ax.transAxes, fontsize=6.5, color="0.4", ha="right", va="bottom")


def _draw_pb_validity_per_frame(ax, poses, frames, tools, stats=None):
    """(B) per-frame PoseBusters validity rate, one line per tool — frames are
    ordered MD snapshots, so the line traces validity along the trajectory. Each
    tool's legend label carries its Cochran–Armitage trend p across the ordered
    frames (exploratory; see orai_cluster_stats.json)."""
    per_tool = (stats or {}).get("per_tool", {})
    d = poses[poses["pb_valid"].notna()]
    x = np.arange(len(frames))
    for t in tools:
        rates = []
        for fr in frames:
            sel = d[(d["frame"] == fr) & (d["tool"] == t)]
            rates.append(float(sel["pb_valid"].mean()) * 100 if len(sel) else np.nan)
        lab = _style(t)[0]
        st = per_tool.get(t)
        if st and su is not None and st.get("p") is not None:
            lab += f" (trend {su.p_stars(st['p'])} {su.fmt_p(st['p'])})"
        ax.plot(x, rates, marker="o", lw=2, markersize=6,
                color=_style(t)[1], label=lab,
                markeredgecolor="black", markeredgewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("Orai1WT-", "") for f in frames], rotation=20, ha="right")
    ax.set_title("PoseBusters validity rate per Orai MD frame")
    ax.set_ylabel("Valid poses (percent of generated)")
    ax.set_ylim(0, 105); ax.legend(fontsize=8)
    ax.grid(alpha=0.3); ax.set_axisbelow(True)
    if stats is not None:
        note = ("trend across ordered MD frames: Cochran–Armitage (exploratory)"
                if not stats.get("skipped") else "trend: n too small — exploratory")
        ax.text(0.02, 0.04, note, transform=ax.transAxes, fontsize=6.5,
                color="0.4", ha="left", va="bottom")


def _draw_cluster_structure(ax, pairs, frames, stats=None):
    """(C) per-frame cluster tightness: mean #clusters + mean dominant fraction,
    pooled over all tools (the per-tool breakdown is the next panel). Each bar
    series carries a Jonckheere–Terpstra trend p across the ordered MD frames
    (per-(frame,ligand) pair; exploratory — see orai_cluster_stats.json)."""
    metrics = (stats or {}).get("metrics", {})

    def _trend_tag(metric, base):
        m = metrics.get(metric)
        if m and su is not None and m.get("p") is not None:
            return f"{base} (trend {su.p_stars(m['p'])} {su.fmt_p(m['p'])})"
        return base

    x = np.arange(len(frames))
    ncl = [pd.to_numeric(pairs[pairs.frame == fr]["n_clusters"], errors="coerce").mean() for fr in frames]
    dom = [pd.to_numeric(pairs[pairs.frame == fr]["dominant_cluster_frac"], errors="coerce").mean() for fr in frames]
    ax.bar(x - 0.2, ncl, 0.4, color="#8172B3", label=_trend_tag("n_clusters", "mean number of clusters"))
    ax2 = ax.twinx()
    ax2.bar(x + 0.2, dom, 0.4, color="#CCB974",
            label=_trend_tag("dominant_cluster_frac", "mean dominant-cluster fraction"))
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("Orai1WT-", "") for f in frames], rotation=20, ha="right")
    ax.set_title("Cluster structure per Orai MD frame")
    ax.set_ylabel("Mean number of distinct clusters")
    ax2.set_ylabel("Mean dominant-cluster fraction"); ax2.set_ylim(0, 1.05)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    ax.grid(alpha=0.25); ax.set_axisbelow(True)
    if stats is not None:
        note = ("trend across ordered MD frames: Jonckheere–Terpstra, "
                "unit (frame,ligand) pair (exploratory)"
                if not stats.get("skipped") else "trend: n too small — exploratory")
        ax.text(0.02, 0.02, note, transform=ax.transAxes, fontsize=6.5,
                color="0.4", ha="left", va="bottom")


def _per_tool_cluster_stats(poses, tools):
    """Per (frame, tool): mean number of distinct global clusters that tool's
    poses occupy, and mean fraction of that tool's poses sitting in the pair's
    DOMINANT (largest) cluster. A per-tool decomposition of the pooled
    'Cluster structure per frame' panel — it shows which tool multiplies the
    site count and which concentrates on the main site.
    Returns {(frame, tool): (mean_n_clusters, mean_dominant_share, n_pairs)}."""
    out: Dict[tuple, tuple] = {}
    if "cluster" not in poses.columns:
        return out
    acc: Dict[tuple, dict] = defaultdict(lambda: {"ncl": [], "dom": []})
    for (frame, _ligand), g in poses.groupby(["frame", "ligand"]):
        if g.empty:
            continue
        dom_label = g["cluster"].value_counts().idxmax()      # largest = dominant
        for t, gt in g.groupby("tool"):
            acc[(frame, t)]["ncl"].append(int(gt["cluster"].nunique()))
            acc[(frame, t)]["dom"].append(float((gt["cluster"] == dom_label).mean()))
    for key, d in acc.items():
        out[key] = (float(np.mean(d["ncl"])), float(np.mean(d["dom"])), len(d["ncl"]))
    return out


def _draw_tool_cluster_contribution(ax, poses, frames, tools, stats=None):
    """(new) per-tool contribution to the cluster structure of each MD frame:
    grouped bars = mean number of distinct clusters a tool occupies (left axis);
    markers = mean share of that tool's poses in the pair's dominant cluster
    (right axis). Together they show how much each docking tool drives the
    multi-cluster spread vs concentrating on the dominant site. A paired-across-
    tools Friedman/Wilcoxon on the per-(frame,ligand) cluster count is annotated
    (exploratory — see orai_cluster_stats.json)."""
    cstats = _per_tool_cluster_stats(poses, tools)
    x = np.arange(len(frames))
    nt = max(len(tools), 1)
    width = 0.8 / nt
    ax2 = ax.twinx()
    for k, t in enumerate(tools):
        off = (k - (nt - 1) / 2) * width
        ncl = [cstats.get((fr, t), (np.nan, np.nan, 0))[0] for fr in frames]
        dom = [cstats.get((fr, t), (np.nan, np.nan, 0))[1] for fr in frames]
        ax.bar(x + off, ncl, width, color=_style(t)[1], alpha=0.85,
               edgecolor="black", linewidth=0.3, label=_style(t)[0])
        ax2.plot(x + off, dom, marker="o", ls="none", markersize=6,
                 color=_style(t)[1], markeredgecolor="black", markeredgewidth=0.6)
    ax2.plot([], [], marker="o", ls="none", color="0.4", markeredgecolor="black",
             label="dominant-cluster share (right axis)")
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("Orai1WT-", "") for f in frames], rotation=20, ha="right")
    ax.set_title("Per-tool contribution to cluster structure per Orai MD frame")
    ax.set_ylabel("Mean number of distinct clusters a tool occupies")
    ax2.set_ylabel("Mean share of a tool's poses in the dominant cluster")
    ax2.set_ylim(0, 1.05)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
    ax.grid(alpha=0.25); ax.set_axisbelow(True)
    if stats is not None and su is not None:
        cap = _paired_summary(stats) if not stats.get("skipped") else None
        note = (f"clusters/pair across tools: {cap}\nunit (frame,ligand) pair; exploratory"
                if cap else "clusters/pair across tools: n too small — exploratory")
        ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=6.5, color="0.4",
                ha="left", va="top",
                bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.8))


def _tool_agreement_matrix(pairs, tools, thr):
    """Symmetric matrix of how often two tools land on the SAME site: fraction of
    the (frame, ligand) pairs where both tools are present and their per-tool
    consensus sites lie within ``thr`` Å. Built from the pre-computed inter-tool
    ``<ab>_<cd>_dist`` columns. Returns (agreement_fraction, n_shared_pairs)."""
    n = len(tools)
    frac = np.full((n, n), np.nan)
    cnt = np.zeros((n, n), dtype=int)
    for i in range(n):
        frac[i, i] = 1.0
        for j in range(i + 1, n):
            a, b = tools[i], tools[j]
            col = f"{a[:2]}_{b[:2]}_dist"
            if col not in pairs.columns:
                col = f"{b[:2]}_{a[:2]}_dist"
            if col not in pairs.columns:
                continue
            v = pd.to_numeric(pairs[col], errors="coerce").dropna()
            if len(v):
                f = float((v <= thr).mean())
                frac[i, j] = frac[j, i] = f
                cnt[i, j] = cnt[j, i] = len(v)
    return frac, cnt


def _draw_tool_agreement_matrix(ax, pairs, tools, thr, stats=None):
    """(new) heat-map of pairwise tool agreement on the binding site. Cells carry
    a Wilcoxon star (per-pair distance vs the ``thr`` agree-line; exploratory)."""
    by_col = _tpa_by_col(stats)
    frac, cnt = _tool_agreement_matrix(pairs, tools, thr)
    n = len(tools)
    im = ax.imshow(np.where(np.isnan(frac), 0.0, frac), cmap="YlGn", vmin=0, vmax=1)
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    labels = [_style(t)[0] for t in tools]
    ax.set_xticklabels(labels, rotation=20, ha="right"); ax.set_yticklabels(labels)
    for i in range(n):
        for j in range(n):
            if i == j:
                ax.text(j, i, "—", ha="center", va="center", fontsize=11, color="0.4")
            elif not np.isnan(frac[i, j]):
                col = f"{tools[i][:2]}_{tools[j][:2]}_dist"
                rec = by_col.get(col) or by_col.get(f"{tools[j][:2]}_{tools[i][:2]}_dist")
                star = f"\n{rec['star']}" if (rec and rec.get("star")) else ""
                ax.text(j, i, f"{frac[i, j]*100:.0f}%\n(n={cnt[i, j]}){star}",
                        ha="center", va="center", fontsize=9,
                        color="black" if frac[i, j] < 0.6 else "white")
            else:
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9, color="0.4")
    ax.set_title(f"How often do two tools agree on a site?\n(per-tool sites within {thr:g} Å; n = shared pairs)")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fraction of shared pairs in agreement")
    ax.grid(False)


def _draw_valid_vs_consensus(ax, poses, thr, stats=None):
    """(D) PoseBusters survivors vs failures — distance to the pair consensus.
    Multi-tool pairs only: for single-tool pairs the consensus is just that
    tool's own cluster centre, so 'nearer the consensus' would be circular.
    Mann–Whitney U + Cliff's δ (unit = pose, clustered by (frame,ligand) pair for
    the median CIs — see orai_cluster_stats.json)."""
    dd = poses[poses["pb_valid"].notna()]
    if "consensus_multitool" in dd:
        dd = dd[dd["consensus_multitool"] == True]  # noqa: E712
    vv = pd.to_numeric(dd[dd["pb_valid"] == True]["dist_to_consensus"], errors="coerce").dropna()  # noqa: E712
    iv = pd.to_numeric(dd[dd["pb_valid"] == False]["dist_to_consensus"], errors="coerce").dropna()  # noqa: E712
    bins = np.linspace(0, 30, 31)
    if len(vv):
        ax.hist(vv, bins=bins, density=True, alpha=0.55, color="#2C7FB8",
                label=f"survives PoseBusters (med {vv.median():.1f} Å)")
    if len(iv):
        ax.hist(iv, bins=bins, density=True, alpha=0.55, color="#D95F0E",
                label=f"fails PoseBusters (med {iv.median():.1f} Å)")
    ax.set_title("Are valid poses nearer the consensus site?\n(multi-tool pairs)")
    ax.set_xlabel("Pose distance to its pair's consensus site (Å)")
    ax.set_ylabel("Probability density"); ax.legend(fontsize=8); ax.grid(alpha=0.25)
    if stats and su is not None and stats.get("p") is not None:
        d, ci = stats.get("cliffs_delta"), stats.get("delta_ci") or [None, None]
        ci_txt = (f" [{ci[0]:.2f}, {ci[1]:.2f}]"
                  if ci and ci[0] is not None and ci[1] is not None else "")
        txt = (f"Mann–Whitney {su.p_stars(stats['p'])} {su.fmt_p(stats['p'])}"
               f"\nCliff's δ={d:.2f}{ci_txt}"
               f"\nunit: pose (clustered by frame×ligand pair); exploratory")
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=7,
                ha="right", va="top",
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))


def _draw_tool_dispersion(ax, pairs, tools, stats=None):
    """(E) per-tool homogeneity: how dispersed each tool's own poses are. A paired-
    across-tools Friedman/Wilcoxon on the per-(frame,ligand) spread is annotated
    (exploratory — see orai_cluster_stats.json)."""
    present = []
    for t in tools:
        col = pairs.get(f"{t}_spread")
        v = (pd.to_numeric(col, errors="coerce").dropna().to_numpy()
             if col is not None else np.array([]))
        if v.size:
            present.append((t, v))
    if present:
        bp = ax.boxplot([v for _, v in present], positions=np.arange(len(present)),
                        widths=0.6, patch_artist=True, showfliers=False,
                        medianprops=dict(color="black"))
        for patch, (t, _) in zip(bp["boxes"], present):
            patch.set_facecolor(_style(t)[1]); patch.set_alpha(0.75)
        ax.set_xticks(np.arange(len(present)))
        ax.set_xticklabels([_style(t)[0] for t, _ in present], rotation=15, ha="right")
    ax.set_title("Per-tool pose dispersion (decisiveness)")
    ax.set_ylabel("Mean spread of a tool's poses per pair (Å)"); ax.grid(alpha=0.25)
    if stats is not None and su is not None:
        cap = _paired_summary(stats) if not stats.get("skipped") else None
        note = (f"spread across tools: {cap}\nunit (frame,ligand) pair; exploratory"
                if cap else "spread across tools: n too small — exploratory")
        ax.text(0.02, 0.98, note, transform=ax.transAxes, fontsize=6.5, color="0.4",
                ha="left", va="top",
                bbox=dict(boxstyle="round", fc="white", ec="0.8", alpha=0.8))


def _draw_axis_depth(ax, poses, tools):
    """(F) binding position along the receptor principal (≈ pore) axis."""
    pd_ = poses.dropna(subset=["axis_depth"])
    for t in tools:
        v = pd.to_numeric(pd_[pd_.tool == t]["axis_depth"], errors="coerce").dropna().to_numpy()
        if v.size:
            ax.hist(v, bins=40, histtype="step", lw=2, color=_style(t)[1], label=_style(t)[0])
    ax.set_title("Binding position along the receptor principal (≈ pore) axis")
    ax.set_xlabel("Projected position along principal axis (Å, 0 = receptor centre)")
    ax.set_ylabel("Number of poses"); ax.legend(fontsize=8); ax.grid(alpha=0.25)
    allv = pd.to_numeric(pd_["axis_depth"], errors="coerce").dropna().to_numpy()
    if allv.size:                                     # robust window (ignore stray tails)
        lo, hi = np.percentile(allv, [0.5, 99.5])
        ax.set_xlim(lo - 5, hi + 5)


def fig_overview(pairs: pd.DataFrame, poses: pd.DataFrame, thr: float, out_dir: Path,
                 stats_out: Optional[dict] = None):
    """Combined 8-panel overview PLUS a standalone PNG per panel in ``panels/``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tools = sorted(poses["tool"].unique())
    frames = sorted(pairs["frame"].unique())

    # Compute the per-panel stats ONCE (guarded); they annotate both the grid
    # cell and the standalone panel, and are written to orai_cluster_stats.json.
    s_pb = s_tpa = s_vc = s_cs = s_tcc = s_td = None
    for name, fn in (("pb_validity_per_frame trend", lambda: _stats_pb_frame_trend(poses, frames, tools)),
                     ("tool-pair agreement", lambda: _stats_tool_pair_agreement(pairs, tools, thr)),
                     ("valid-vs-consensus", lambda: _stats_valid_vs_consensus(poses)),
                     ("cluster-structure trend", lambda: _stats_cluster_structure_trend(pairs, frames)),
                     ("tool cluster contribution", lambda: _stats_tool_cluster_contribution(poses, tools)),
                     ("tool dispersion", lambda: _stats_tool_dispersion(pairs, tools))):
        try:
            r = fn()
        except Exception as e:                          # pragma: no cover
            print(f"  [stats] {name} failed: {e}"); r = None
        if name.startswith("pb_validity"):   s_pb = r
        elif name.startswith("tool-pair"):    s_tpa = r
        elif name.startswith("valid-vs"):     s_vc = r
        elif name.startswith("cluster-str"):  s_cs = r
        elif name.startswith("tool cluster"): s_tcc = r
        elif name.startswith("tool disp"):    s_td = r
    if stats_out is not None:
        for key, val in (("pb_validity_per_frame", s_pb), ("tool_pair_agreement", s_tpa),
                         ("valid_vs_consensus", s_vc), ("cluster_structure_trend", s_cs),
                         ("tool_cluster_contribution", s_tcc), ("tool_dispersion", s_td)):
            if val is not None:
                stats_out[key] = val

    # (name, draw-fn) — one entry per panel; the draw-fn takes only an Axes so the
    # same call renders the grid cell and the standalone figure.
    panels = [
        ("cross_tool_agreement",     lambda a: _draw_cross_tool_agreement(a, pairs, thr, tools, s_tpa)),
        ("pb_validity_per_frame",    lambda a: _draw_pb_validity_per_frame(a, poses, frames, tools, s_pb)),
        ("cluster_structure",        lambda a: _draw_cluster_structure(a, pairs, frames, s_cs)),
        ("tool_cluster_contribution", lambda a: _draw_tool_cluster_contribution(a, poses, frames, tools, s_tcc)),
        ("tool_agreement_matrix",    lambda a: _draw_tool_agreement_matrix(a, pairs, tools, thr, s_tpa)),
        ("valid_vs_consensus",       lambda a: _draw_valid_vs_consensus(a, poses, thr, s_vc)),
        ("tool_dispersion",          lambda a: _draw_tool_dispersion(a, pairs, tools, s_td)),
        ("axis_depth",               lambda a: _draw_axis_depth(a, poses, tools)),
    ]

    # combined overview (2 x 4)
    fig, ax = plt.subplots(2, 4, figsize=(25, 12))
    for (name, fn), a in zip(panels, ax.ravel()):
        fn(a)
    _label_panels(ax)
    fig.suptitle("Orai × benchmark-ligand docked poses — clustering, cross-tool agreement, "
                 f"PoseBusters survival & MD-frame distribution (n={len(pairs)} pairs, "
                 f"tools: {', '.join(_style(t)[0] for t in tools)})", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "orai_pose_cluster_overview.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    # standalone per-panel PNGs (the split-out individual graphs)
    panel_dir = out_dir / "panels"; panel_dir.mkdir(parents=True, exist_ok=True)
    for name, fn in panels:
        f, a = plt.subplots(figsize=(8, 6.5))
        fn(a)
        f.tight_layout()
        f.savefig(panel_dir / f"orai_{name}.png", dpi=130, bbox_inches="tight")
        plt.close(f)
    print(f"  Individual panels: {panel_dir}/ ({len(panels)} PNGs)")
    return p


def _short_frame(f):
    return str(f).replace("Orai1WT-", "").replace("MDSnap-", "")


def _short_lig(l):
    return (str(l).replace("-OPT-Singlet", "").replace("-prot-OPT", "")
            .replace("-OPT", ""))


# tool-pairs, in the fixed display order, with the pre-computed site-distance
# column each maps to (see analyze_pair's `<ab>_<cd>_dist` naming).
_TOOL_PAIRS = [("autodock", "diffdock", "au_di_dist", "AutoDock\n↔ DiffDock"),
               ("autodock", "equibind", "au_eq_dist", "AutoDock\n↔ EquiBind"),
               ("diffdock", "equibind", "di_eq_dist", "DiffDock\n↔ EquiBind")]


def _min_cross_pose_dist(pose_recs, ta, tb):
    """Closest distance (Å) between ANY pose of tool ``ta`` and ANY pose of tool
    ``tb`` in one (frame, ligand) pair — do the two pose CLOUDS ever touch,
    regardless of where each tool's dominant cluster sits. NaN if either absent."""
    A = np.array([p["_cent"] for p in pose_recs if p["tool"] == ta and p.get("_cent") is not None])
    B = np.array([p["_cent"] for p in pose_recs if p["tool"] == tb and p.get("_cent") is not None])
    if len(A) == 0 or len(B) == 0:
        return np.nan
    return float(np.min(np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)))


def _draw_dist_matrix(ax, M, row_labels, col_labels, n_pairs, frames_seq,
                      norm, cmap, vmax, title, show_yticks):
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(M.shape[1])); ax.set_xticklabels(col_labels, fontsize=9)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels if show_yticks else [""] * len(row_labels), fontsize=8)
    ax.axvline(M.shape[1] - 1.5, color="black", lw=1.2)          # pairwise | mean divider
    for i in range(1, n_pairs):                                  # group rows by MD frame
        if frames_seq[i] != frames_seq[i - 1]:
            ax.axhline(i - 0.5, color="black", lw=0.8)
    ax.axhline(n_pairs - 0.5, color="black", lw=2.0)             # summary-row divider
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8,
                        fontweight="bold" if i == n_pairs else "normal",
                        color="white" if v >= 0.78 * vmax else "black")
            else:
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=7, color="0.4")
    ax.set_title(title, fontsize=10)
    return im


def fig_divergence_matrix(pairs: pd.DataFrame, poses_by_pair: dict, thr: float, out_dir: Path,
                          stats_out: Optional[dict] = None):
    """Two matched matrices over every (MD snapshot × ligand) pair (rows) and
    every tool-pair (cols), sharing one Å colour scale (green ≤ ``thr`` = agree,
    red = divergent, grey = a tool produced no pose):

      LEFT  — distance between the two tools' DOMINANT sites (largest-cluster
              centroid): do the tools *prefer* the same spot?
      RIGHT — closest distance between ANY two of their poses: do the pose
              *clouds* ever touch, even away from each tool's main site?

    Left red + right green means the clouds graze but the preferences differ;
    both red means the tools are genuinely working in different regions. A
    bottom row gives the all-pairs median; a trailing column the per-pair mean
    over present tool-pairs."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm

    if pairs.empty:
        return None
    dfp = pairs.sort_values(["frame", "ligand"]).reset_index(drop=True)
    n_pairs = len(dfp)

    # Tool-pair signed-rank stars annotate the LEFT (dominant-site) matrix, whose
    # columns are exactly the distances the Wilcoxon test uses. Reuse the stats
    # already computed by fig_overview when available, else compute here.
    s_tpa = (stats_out or {}).get("tool_pair_agreement")
    if s_tpa is None:
        try:
            tools_seen = sorted({p["tool"] for recs in poses_by_pair.values() for p in recs})
            s_tpa = _stats_tool_pair_agreement(pairs, tools_seen, thr)
            if stats_out is not None:
                stats_out.setdefault("tool_pair_agreement", s_tpa)
        except Exception:                              # pragma: no cover
            s_tpa = None
    _by_col = _tpa_by_col(s_tpa)
    left_labels = []
    for _ta, _tb, col, lab in _TOOL_PAIRS:
        rec = _by_col.get(col)
        left_labels.append(lab + (f"\n{rec['star']}" if (rec and rec.get("star")) else ""))
    left_labels.append("mean of\npresent pairs")

    # Tool-pair signed-rank stars for the RIGHT (cloud-touch) matrix, whose columns
    # are exactly the per-pair CLOSEST cross-tool pose distances the Wilcoxon test
    # uses (do the pose clouds graze within ``thr`` Å?).
    s_touch = (stats_out or {}).get("cross_pose_touch")
    if s_touch is None:
        try:
            s_touch = _stats_cross_pose_touch(pairs, poses_by_pair, thr)
            if stats_out is not None:
                stats_out.setdefault("cross_pose_touch", s_touch)
        except Exception:                              # pragma: no cover
            s_touch = None
    _touch_by_col = {rec.get("col"): rec for rec in (s_touch or {}).get("pairs", {}).values()
                     if isinstance(rec, dict) and rec.get("col")}
    right_labels = []
    for _ta, _tb, col, lab in _TOOL_PAIRS:
        rec = _touch_by_col.get(col)
        right_labels.append(lab + (f"\n{rec['star']}" if (rec and rec.get("star")) else ""))
    right_labels.append("mean of\npresent pairs")

    # LEFT: dominant-site distance (pre-computed columns) + per-row mean
    S = np.full((n_pairs, len(_TOOL_PAIRS)), np.nan)
    for j, (_ta, _tb, col, _lab) in enumerate(_TOOL_PAIRS):
        if col in dfp.columns:
            S[:, j] = pd.to_numeric(dfp[col], errors="coerce").to_numpy()
    S = np.column_stack([S, np.nanmean(S, axis=1)])

    # RIGHT: closest cross-tool pose distance (from raw centroids) + per-row mean
    D = np.full((n_pairs, len(_TOOL_PAIRS)), np.nan)
    for i, (_, r) in enumerate(dfp.iterrows()):
        recs = poses_by_pair.get((r["frame"], r["ligand"]), [])
        for j, (ta, tb, _col, _lab) in enumerate(_TOOL_PAIRS):
            D[i, j] = _min_cross_pose_dist(recs, ta, tb)
    D = np.column_stack([D, np.nanmean(D, axis=1)])

    row_labels = [f"{_short_frame(r.frame)} · {_short_lig(r.ligand)}" for _, r in dfp.iterrows()]
    frames_seq = list(dfp["frame"])
    with warnings.catch_warnings():                        # all-NaN column → NaN median
        warnings.simplefilter("ignore")
        S = np.vstack([S, np.nanmedian(S, axis=0)])
        D = np.vstack([D, np.nanmedian(D, axis=0)])
    row_labels.append("median (all pairs)")

    both = np.concatenate([S[np.isfinite(S)].ravel(), D[np.isfinite(D)].ravel()])
    vmax = max(float(np.nanmax(both)) if both.size else thr * 2, thr + 1.0)
    norm = TwoSlopeNorm(vmin=0.0, vcenter=float(thr), vmax=vmax)
    cmap = plt.get_cmap("RdYlGn_r").copy(); cmap.set_bad("#dddddd")

    fig, axes = plt.subplots(1, 2, figsize=(15.5, max(5.0, 0.5 * len(row_labels) + 2.2)))
    _draw_dist_matrix(axes[0], S, row_labels, left_labels, n_pairs, frames_seq, norm, cmap, vmax,
                      "Do the tools prefer the same site?\n"
                      "distance between each tool's DOMINANT (largest-cluster) site (Å)", True)
    im = _draw_dist_matrix(axes[1], D, row_labels, right_labels, n_pairs, frames_seq, norm, cmap, vmax,
                           "Do the tools' pose clouds ever touch?\n"
                           "closest distance between ANY two poses of the tools (Å)", False)
    for a, lab in zip(axes, "AB"):                         # panel labels clear of the 2-line titles
        a.text(-0.02, 1.17, f"({lab})", transform=a.transAxes,
               fontsize=13, fontweight="bold", va="bottom", ha="right")
    if _by_col:                                            # star meaning for the LEFT column headers
        axes[0].text(0.0, -0.02, "star: Wilcoxon per-pair distance vs threshold (exploratory)",
                     transform=axes[0].transAxes, fontsize=6.5, color="0.4", ha="left", va="top")
    if _touch_by_col:                                      # star meaning for the RIGHT column headers
        axes[1].text(0.0, -0.02, "star: Wilcoxon per-pair closest-pose distance vs threshold (exploratory)",
                     transform=axes[1].transAxes, fontsize=6.5, color="0.4", ha="left", va="top")
    fig.suptitle("Cross-tool divergence across every Orai MD-snapshot × ligand pair  "
                 f"(green ≤ {thr:g} Å = agree, red = divergent, grey = tool absent)",
                 fontsize=12, y=1.0)
    cbar = fig.colorbar(im, ax=axes, fraction=0.046, pad=0.04, extend="max")
    cbar.set_label("Distance (Å)")
    cbar.set_ticks([t for t in sorted({0.0, float(thr), 10.0, 20.0, 30.0, round(vmax)}) if t <= vmax])
    p = out_dir / "orai_tool_divergence_matrix.png"
    fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)
    return p


def fig_examples(records, poses_by_pair, thr, out_dir, n_ex=6):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    # pick pairs with the most tools + most poses for a clear picture
    ex = sorted([r for r in records if not r.get("skipped")],
                key=lambda r: (r["n_tools"], r["n_poses"]), reverse=True)[:n_ex]
    if not ex:
        return None
    ncol = 3
    nrow = (len(ex) + ncol - 1) // ncol
    fig = plt.figure(figsize=(6.2 * ncol, 5.4 * nrow))
    axes = []
    for i, r in enumerate(ex):
        ax = fig.add_subplot(nrow, ncol, i + 1, projection="3d")
        axes.append(ax)
        pr = poses_by_pair[(r["frame"], r["ligand"])]
        C = np.vstack([p["_cent"] for p in pr])
        labels = np.array([p["cluster"] for p in pr])
        clusters = sorted(set(labels.tolist()))
        cmap = cm.get_cmap("tab10", max(len(clusters), 1))
        for j, p in enumerate(pr):
            valid = p["pb_valid"]
            ax.scatter(*C[j], color=cmap(clusters.index(labels[j])),
                       marker=_style(p["tool"])[2], s=40,
                       edgecolors=("black" if valid is True else ("red" if valid is False else "grey")),
                       linewidths=1.0 if valid is not None else 0.3, alpha=0.85)
        ax.set_title(f"{r['ligand']} on {r['frame'].replace('Orai1WT-','')}\n"
                     f"{r['n_clusters']} clusters, dom {r['dominant_cluster_frac']:.0%}, "
                     f"agree {r['mean_inter_tool_dist']:.1f} Å", fontsize=8)
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    _label_panels(axes, fontsize=12)
    # Edge legend names only the validity categories actually plotted. Under
    # --pb-valid-only every marker is valid (black), so "red = invalid / grey =
    # unknown" would be a phantom legend describing categories that never appear.
    present = [p["pb_valid"] for r in ex
               for p in poses_by_pair[(r["frame"], r["ligand"])]]
    edge_bits = []
    if any(v is True for v in present):
        edge_bits.append("black = PoseBusters-valid")
    if any(v is False for v in present):
        edge_bits.append("red = invalid")
    if any(v is not True and v is not False for v in present):
        edge_bits.append("grey = unknown")
    edge_txt = ", ".join(edge_bits) if edge_bits else "black = PoseBusters-valid"
    fig.suptitle("Example pose clusterings — colour = cluster, marker = tool "
                 "(○ AutoDock / △ DiffDock / □ EquiBind), edge: " + edge_txt,
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "orai_example_clusters.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def fig_descriptor_quality(pairs: pd.DataFrame, features_csv: str, out_dir: Path,
                           stats_out: Optional[dict] = None):
    """Which ligand types dock CONSISTENTLY on Orai (no crystal)? Relate ligand
    descriptors / Lipinski Ro5 / PCA to the quality proxies we DO have: PoseBusters
    validity, cross-tool agreement, and cluster tightness. Reuses the same shared
    descriptor module + Lipinski/PCA as the crystal report, so findings transfer."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    try:
        from ligand_descriptors import load_descriptors, add_pca, PCA_FEATURES, nice
    except Exception:
        return None
    desc = load_descriptors(features_csv)
    if desc is None:
        print("  (ligand feature CSV not found; skipping descriptor analysis)")
        return None
    agg = pairs.groupby("ligand").agg(
        pb_valid_rate=("pb_valid_rate", "mean"),
        inter_tool_dist=("mean_inter_tool_dist", "mean"),
        dominant_frac=("dominant_cluster_frac", "mean"),
        n_frames=("frame", "nunique")).reset_index()
    d = agg.set_index("ligand").join(desc, how="inner")
    if len(d) < 12:
        return None
    d, loadings, evr, _ = add_pca(d)
    bins, labs = [-0.1, 2, 5, 8, 11, 100], ["0–2", "3–5", "6–8", "9–11", "12+"]
    d["_rb"] = pd.cut(pd.to_numeric(d["rot_bonds"], errors="coerce"), bins=bins, labels=labs)
    d.reset_index().to_csv(out_dir / "descriptor_vs_quality.csv", index=False)
    has_pb = d["pb_valid_rate"].notna()

    fig, ax = plt.subplots(2, 2, figsize=(14, 11))

    # (A) PB validity by Lipinski Ro5
    if has_pb.any():
        grp = [("Ro5 pass", d[d.ro5_pass & has_pb]["pb_valid_rate"].dropna() * 100),
               ("Ro5 fail", d[~d.ro5_pass & has_pb]["pb_valid_rate"].dropna() * 100)]
        ax[0, 0].boxplot([g for _, g in grp], labels=[f"{n}\n(n={len(g)})" for n, g in grp],
                         showfliers=False)
    ax[0, 0].set_title("PoseBusters validity vs Lipinski Rule-of-Five")
    ax[0, 0].set_ylabel("Per-ligand valid-pose rate (percent)")

    # (B) cross-tool agreement vs flexibility (lower distance = tools converge) —
    # a trend across ordinal rotatable-bond bins, so a line shows it directly.
    g = d.groupby("_rb")["inter_tool_dist"].agg(["median", "size"])
    ax[0, 1].plot(range(len(g)), g["median"], marker="o", color="#55A868",
                  lw=2, markersize=7, markeredgecolor="black", markeredgewidth=0.6)
    for i, (m, nn) in enumerate(zip(g["median"], g["size"])):
        if nn and m == m:
            ax[0, 1].annotate(f"n={int(nn)}", (i, m), textcoords="offset points",
                              xytext=(0, 8), ha="center", fontsize=8)
    ax[0, 1].set_xticks(range(len(g))); ax[0, 1].set_xticklabels(g.index)
    ax[0, 1].set_title("Cross-tool agreement vs ligand flexibility\n(lower = tools converge)")
    ax[0, 1].set_xlabel("Rotatable bonds")
    ax[0, 1].set_ylabel("Median inter-tool consensus distance (Å)")
    ax[0, 1].grid(alpha=0.3); ax[0, 1].set_axisbelow(True)

    # (C) Spearman of properties vs validity & vs agreement — now with p-values
    # (previously discarded), BH-FDR-corrected across the properties (per target),
    # and significance stars on the bars. n≈{len(d)} ligands ⇒ exploratory.
    feats = [f for f in PCA_FEATURES if f in d.columns] + [c for c in ("PC1", "PC2") if c in d.columns]
    rows = []                # [feature, rho_valid, rho_agree, p_valid, p_agree]
    for f in feats:
        v = pd.to_numeric(d[f], errors="coerce").to_numpy()
        rho, praw = {}, {}
        for tgt, kk in (("pb_valid_rate", "valid"), ("inter_tool_dist", "agree")):
            tv = pd.to_numeric(d[tgt], errors="coerce").to_numpy()
            m = ~(np.isnan(v) | np.isnan(tv))
            if m.sum() > 10 and su is not None:
                try:
                    r, p, _n = su.spearman(v, tv)
                except Exception:
                    r, p = np.nan, np.nan
            elif m.sum() > 10:
                r, p = spearmanr(v[m], tv[m])[0], np.nan
            else:
                r, p = np.nan, np.nan
            rho[kk], praw[kk] = r, p
        rows.append([f, rho["valid"], rho["agree"], praw["valid"], praw["agree"]])
    # BH-FDR across properties, separately per target
    if su is not None:
        try:
            qv = su.bh_fdr([r[3] for r in rows])
            qa = su.bh_fdr([r[4] for r in rows])
        except Exception:
            qv = qa = [np.nan] * len(rows)
    else:
        qv = qa = [np.nan] * len(rows)
    for r, q1, q2 in zip(rows, qv, qa):
        r.append(_f(q1)); r.append(_f(q2))     # -> [.., q_valid, q_agree]
    rows.sort(key=lambda x: (x[1] if x[1] == x[1] else 0))
    y = np.arange(len(rows)); w = 0.4
    ax[1, 0].barh(y - w / 2, [r[1] for r in rows], w, color="#2C7FB8", label="vs validity (↑ = more valid)")
    ax[1, 0].barh(y + w / 2, [r[2] for r in rows], w, color="#D95F0E", label="vs inter-tool dist (↑ = more disagreement)")
    # star bars whose BH-FDR q < 0.05
    if su is not None:
        for k, r in enumerate(rows):
            for val, q, dy in ((r[1], r[5], -w / 2), (r[2], r[6], +w / 2)):
                s = su.p_stars(q)
                if s and s not in ("", "ns") and val == val:
                    ax[1, 0].text(val + (0.02 if val >= 0 else -0.02), k + dy, s,
                                  va="center", ha="left" if val >= 0 else "right",
                                  fontsize=8, fontweight="bold")
    ax[1, 0].set_yticks(y); ax[1, 0].set_yticklabels([nice(r[0]) for r in rows], fontsize=7)
    ax[1, 0].axvline(0, color="k", lw=0.8); ax[1, 0].legend(fontsize=7)
    ax[1, 0].set_title("Which properties track Orai docking quality?\n"
                       "(Spearman ρ; star = BH-FDR q<0.05, exploratory)")
    ax[1, 0].set_xlabel("Spearman correlation")
    if stats_out is not None:
        stats_out["descriptor_quality"] = {
            "test": "spearman + BH-FDR across properties",
            "unit": f"per-ligand (n={len(d)}) — exploratory",
            "properties": [{"feature": r[0], "rho_valid": _f(r[1]), "p_valid": _f(r[3]),
                            "q_valid": r[5], "rho_agree": _f(r[2]), "p_agree": _f(r[4]),
                            "q_agree": r[6]} for r in rows]}

    # (D) PCA chemical space coloured by validity
    if "PC1" in d and "PC2" in d and has_pb.any():
        dd = d[has_pb]
        sc = ax[1, 1].scatter(dd["PC1"], dd["PC2"], c=dd["pb_valid_rate"] * 100,
                              cmap="RdYlGn", vmin=0, vmax=100, s=24, edgecolors="k", linewidths=0.2)
        plt.colorbar(sc, ax=ax[1, 1], label="Valid-pose rate (percent)")
        e1 = evr[0] * 100 if len(evr) > 0 else 0
        e2 = evr[1] * 100 if len(evr) > 1 else 0
        ax[1, 1].set_xlabel(f"Principal Component 1 ({e1:.0f}% of variance)")
        ax[1, 1].set_ylabel(f"Principal Component 2 ({e2:.0f}% of variance)")
        ax[1, 1].set_title("Ligand chemical space coloured by PoseBusters validity")

    for a in ax.ravel():
        a.grid(alpha=0.2)
    _label_panels(ax)
    fig.suptitle("Which ligand types dock consistently on Orai? Descriptors vs validity, "
                 f"tool-agreement & tightness (n={len(d)} ligands; no crystal reference)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "orai_descriptor_quality.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Reference-free per-tool cluster QUALITY  (+ raw-vs-filtered tightening)
# ------------------------------------------------------------------------
# The panels above answer WHERE the tools bind and whether they AGREE. This
# block answers a different question — for each tool ON ITS OWN, how well-formed
# is the cloud of poses it produces for a (frame, ligand) pair? With no crystal
# there is no RMSD-to-native, so "quality" is read from internal-validity indices
# on the 3D centroid cloud (silhouette / Calinski–Harabasz / Davies–Bouldin,
# compactness vs separation) plus a bootstrap-Jaccard STABILITY of the dominant
# cluster (does it survive resampling?). Each tool is clustered on its own poses
# so the tools compare head to head; ``n_poses`` is recorded because CH/DB/
# silhouette are pose-count sensitive (AutoDock ~10, DiffDock up to ~30, EquiBind
# was sampled ~30/pair but is thinned by the pore + PoseBusters filters). Two pose
# populations are scored so the effect of filtering is a first-class output:
#   raw       — every generated pose (incl. transmembrane + PB-invalid)
#   filtered  — the pore set (transmembrane-excluded) ∩ PoseBusters-valid
# ════════════════════════════════════════════════════════════════════════

# (key, short label, y-axis label, higher_is_better)
QUALITY_METRICS = [
    ("silhouette",  "Silhouette",          "Silhouette (−1…1)",                       True),
    ("ch_score",    "Calinski–Harabasz",   "Calinski–Harabasz index",                 True),
    ("db_score",    "Davies–Bouldin",      "Davies–Bouldin index",                    False),
    ("compactness", "Compactness",         "Within-cluster spread (Å)",               False),
    ("separation",  "Separation",          "Nearest-cluster distance (Å)",            True),
    ("stability",   "Bootstrap stability", "Dominant-cluster Jaccard stability (0…1)", True),
]
_QM_LABEL = {k: (short, ylab, hb) for k, short, ylab, hb in QUALITY_METRICS}


def _order_tools(values) -> List[str]:
    """Canonical tool order (AutoDock, DiffDock, EquiBind) then any extras."""
    seen = set(str(v) for v in values)
    ordered = [t for t in ("autodock", "diffdock", "equibind") if t in seen]
    ordered += [t for t in sorted(seen) if t not in ordered]
    return ordered


def compute_tool_cluster_quality(df: pd.DataFrame, pose_set: str,
                                 site_method: str, pocket_radius: float,
                                 min_poses: int = 3, n_boot: int = 25,
                                 cap_n: int = 0, seed: int = 0) -> List[dict]:
    """One row per (frame, ligand, tool): internal cluster-validity indices +
    bootstrap stability of that tool's own centroid cloud. ``pose_set`` labels the
    input population ('raw' or 'filtered') so raw-vs-filtered pairs on identity."""
    rows: List[dict] = []
    rng = np.random.RandomState(seed)
    keys = [k for k, *_ in QUALITY_METRICS]
    for (frame, ligand, tool), sub in df.groupby(["frame", "ligand", "tool"]):
        C = np.vstack(sub["_cent"].to_numpy())
        n = len(C)
        if cap_n and n > cap_n:                       # optional pose-count equalisation
            C = C[rng.choice(n, cap_n, replace=False)]; n = cap_n
        row = {"frame": frame, "ligand": ligand, "tool": tool,
               "pose_set": pose_set, "n_poses": int(n),
               "n_clusters": int(n >= 1), "dominant_frac": np.nan}
        if n < min_poses:                             # too few poses to form/score clusters
            row.update({k: np.nan for k in keys})
            rows.append(row)
            continue
        dm = _centroid_dm(C)
        labels, k, sil = cluster_sites(dm, C, site_method, pocket_radius)
        idx = _internal_indices(C, labels)
        stab, _ = _bootstrap_stability(C, labels, site_method, pocket_radius, n_boot=n_boot, seed=seed)
        sizes = np.bincount(labels)
        row.update({
            "n_clusters": int(k),
            "dominant_frac": round(float(sizes.max()) / n, 3),
            "silhouette": round(float(sil), 3) if sil == sil else np.nan,
            "ch_score": idx["ch_score"], "db_score": idx["db_score"],
            "compactness": idx["compactness"], "separation": idx["separation"],
            "stability": stab,
        })
        rows.append(row)
    return rows


def _raw_and_pore_csv(per_pose_csv: str, raw_override: Optional[str]) -> Tuple[Path, Path]:
    """Resolve the (raw, pore) pose CSVs for the raw-vs-filtered comparison.
    raw  = every generated pose (…filtered_results.csv, incl. transmembrane + PB-invalid).
    pore = the transmembrane-excluded set (…filtered_results.no_tm.csv)."""
    p = Path(per_pose_csv); name = p.name
    if name.endswith(".no_tm.csv"):
        pore = p
        raw = Path(raw_override) if raw_override else p.with_name(name[:-len(".no_tm.csv")] + ".csv")
    else:
        raw = Path(raw_override) if raw_override else p
        cand = (p.with_name(name[:-len(".csv")] + ".no_tm.csv")
                if name.endswith(".csv") else p)
        pore = cand if cand.exists() else p
    return raw, pore


# ── quality figures (descriptive only; inferential tests -> *_stats.txt) ────

def _quality_boxplot(ax, qs: pd.DataFrame, metric: str, tools: List[str]):
    short, ylab, hb = _QM_LABEL[metric]
    present = [(t, pd.to_numeric(qs[qs.tool == t][metric], errors="coerce").dropna().to_numpy())
               for t in tools]
    pos = np.arange(len(present))
    bp = ax.boxplot([v if v.size else np.array([np.nan]) for _, v in present],
                    positions=pos, widths=0.6, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black"))
    for patch, (t, _) in zip(bp["boxes"], present):
        patch.set_facecolor(_style(t)[1]); patch.set_alpha(0.75)
    for p_, (t, v) in zip(pos, present):              # jitter small-n so the box isn't a lie
        if 0 < v.size <= 80:
            jit = (np.linspace(-1, 1, v.size) if v.size > 1 else np.zeros(1)) * 0.18
            ax.scatter(np.full(v.size, p_) + jit, v, s=9, color=_style(t)[1],
                       edgecolor="white", lw=0.3, zorder=3)
    ax.set_xticks(pos)
    ax.set_xticklabels([f"{_style(t)[0]}\n(n={v.size})" for t, v in present], fontsize=8)
    ax.set_ylabel(ylab, fontsize=9)
    ax.set_title(f"{short}  ({'↑ better' if hb else '↓ better'})", fontsize=10)
    ax.grid(alpha=0.25, axis="y"); ax.set_axisbelow(True)


def fig_cluster_quality(q: pd.DataFrame, out_dir: Path, pose_set: str = "raw"):
    """Six per-tool intrinsic quality indices (silhouette, Calinski–Harabasz,
    Davies–Bouldin, compactness, separation, bootstrap stability)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    qs = q[q["pose_set"] == pose_set]
    if qs.empty:
        return None
    tools = _order_tools(qs["tool"])
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    for a, (key, *_rest) in zip(ax.ravel(), QUALITY_METRICS):
        _quality_boxplot(a, qs, key, tools)
    _label_panels(ax)
    fig.suptitle("Per-tool intrinsic cluster QUALITY on Orai (reference-free; "
                 f"'{pose_set}' pose cloud) — each tool clustered on its own poses; "
                 "sparse panels (CH/DB/separation) are the multi-cluster subset (k≥2)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "orai_cluster_quality_metrics.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def fig_cluster_compactness_separation(q: pd.DataFrame, out_dir: Path, pose_set: str = "raw"):
    """Compactness (tightness) vs separation for every multi-cluster (frame,ligand,
    tool): upper-left (tight + well-separated) is best."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    qs = q[q["pose_set"] == pose_set].copy()
    qs = qs[pd.to_numeric(qs["separation"], errors="coerce").notna()
            & pd.to_numeric(qs["compactness"], errors="coerce").notna()]
    if qs.empty:
        return None
    tools = _order_tools(qs["tool"])
    fig, ax = plt.subplots(figsize=(9, 7.6))
    for t in tools:
        g = qs[qs.tool == t]
        x = pd.to_numeric(g["compactness"], errors="coerce").to_numpy()
        y = pd.to_numeric(g["separation"], errors="coerce").to_numpy()
        ax.scatter(x, y, s=22, color=_style(t)[1], alpha=0.5, edgecolor="white",
                   lw=0.3, label=f"{_style(t)[0]} (n={len(g)})")
        ax.plot(np.median(x), np.median(y), marker="X", ms=15, color=_style(t)[1],
                mec="black", mew=1.3, zorder=5)
    ax.set_xlabel("Cluster compactness — within-cluster spread (Å, lower = tighter)")
    ax.set_ylabel("Cluster separation — nearest-cluster distance (Å, higher = better)")
    ax.set_title("Are each tool's clusters tight AND well-separated?\n"
                 f"(multi-cluster (frame,ligand,tool) cases, k≥2; '{pose_set}' cloud; "
                 "upper-left = best; X = per-tool median)", fontsize=11)
    ax.legend(fontsize=9, title="Tool (per-tool median = X)"); ax.grid(alpha=0.25)
    ax.set_axisbelow(True)
    fig.tight_layout()
    p = out_dir / "orai_cluster_compactness_separation.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


_FILTER_SET_STYLE = {
    "raw":      dict(color="#B0B0B0", hatch=None,  label="raw cloud (all generated poses)"),
    "filtered": dict(color="#2E7D32", hatch="///", label="filtered (outside-pore ∩ PB-valid)"),
}


def fig_cluster_quality_filtering(q: pd.DataFrame, out_dir: Path):
    """Does filtering to PB-valid poses outside the pore tighten each tool's clusters?
    Raw cloud vs filtered cloud, per tool, on the metrics with an unambiguous
    "better" direction (compactness ↓, silhouette ↑, stability ↑)."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not {"raw", "filtered"}.issubset(set(q["pose_set"])):
        return None
    tools = _order_tools(q["tool"])
    metrics = [("compactness", "Within-cluster spread (Å) — lower = tighter"),
               ("silhouette",  "Silhouette — higher = better"),
               ("stability",   "Bootstrap stability — higher = better")]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    width = 0.36
    for ax, (key, ylab) in zip(axes, metrics):
        for si, ps in enumerate(("raw", "filtered")):
            st = _FILTER_SET_STYLE[ps]
            data, pos = [], []
            for ti, t in enumerate(tools):
                v = pd.to_numeric(q[(q.tool == t) & (q.pose_set == ps)][key],
                                  errors="coerce").dropna().to_numpy()
                data.append(v if v.size else np.array([np.nan]))
                pos.append(ti + (si - 0.5) * (width + 0.04))
            bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True,
                            showfliers=False, medianprops=dict(color="black"))
            for box in bp["boxes"]:
                box.set(facecolor=st["color"], alpha=0.55, hatch=st["hatch"], edgecolor="0.3")
        ax.set_xticks(np.arange(len(tools)))
        ax.set_xticklabels([_style(t)[0] for t in tools], fontsize=9)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_title(key.capitalize(), fontsize=10)
        ax.grid(alpha=0.25, axis="y"); ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_FILTER_SET_STYLE[s]["color"],
                             alpha=0.55, hatch=_FILTER_SET_STYLE[s]["hatch"],
                             edgecolor="0.3", label=_FILTER_SET_STYLE[s]["label"])
               for s in ("raw", "filtered")]
    axes[0].legend(handles=handles, fontsize=8, loc="best")
    _label_panels(np.asarray(axes))
    fig.suptitle("Does filtering to PB-valid poses outside the pore tighten each tool's clusters? "
                 "Raw cloud vs filtered cloud (per tool; paired tests in "
                 "orai_cluster_quality_stats.txt)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "orai_cluster_quality_filtering.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ── quality stats (paired across tools; paired raw-vs-filtered) ─────────────

def _stats_cluster_quality_paired(q: pd.DataFrame, pose_set: str = "raw") -> dict:
    """Per metric: paired-across-tools Friedman + Wilcoxon on the per-(frame,ligand)
    value, listwise-complete across tools (paired design; unit = (frame,ligand) pair;
    no pose pseudoreplication). Exploratory at the JKU n."""
    out = {"test": "friedman + wilcoxon (paired across tools)",
           "pose_set": pose_set, "unit": "(frame,ligand) pair", "metrics": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    qs = q[q["pose_set"] == pose_set]
    tools = _order_tools(qs["tool"])
    for key, *_ in QUALITY_METRICS:
        piv = qs.pivot_table(index=["frame", "ligand"], columns="tool", values=key,
                             aggfunc="first")
        piv = piv.reindex(columns=[t for t in tools if t in piv.columns])
        complete = piv.dropna()
        rec = {"n_pairs_complete": int(len(complete)),
               "medians": {t: _f(pd.to_numeric(qs[qs.tool == t][key], errors="coerce").median())
                           for t in piv.columns}}
        if len(complete) >= _MIN_PAIRS_WILCOXON and complete.shape[1] >= 2:
            try:
                rec.update(_paired_tool_test(complete))
            except Exception as e:                      # pragma: no cover
                rec["error"] = str(e)
        else:
            rec["note"] = "n too small — exploratory"
        out["metrics"][key] = rec
    return out


def _stats_cluster_quality_filtering(q: pd.DataFrame) -> dict:
    """Per tool: paired raw-vs-filtered Wilcoxon signed-rank on compactness /
    silhouette / stability over the (frame,ligand) present in BOTH pose sets.
    rank-biserial sign is for (filtered − raw): compactness < 0 ⇒ filtering tightens;
    silhouette/stability > 0 ⇒ filtering improves. Holm across all (tool × metric)."""
    out = {"test": "wilcoxon signed-rank (paired filtered vs raw)",
           "unit": "(frame,ligand) present in BOTH pose sets", "per_tool": {}}
    if su is None:
        out["skipped"] = "stats_utils unavailable"
        return out
    metrics = ["compactness", "silhouette", "stability"]
    entries = []
    for t in _order_tools(q["tool"]):
        out["per_tool"][t] = {}
        raw = q[(q.tool == t) & (q.pose_set == "raw")].set_index(["frame", "ligand"])
        fil = q[(q.tool == t) & (q.pose_set == "filtered")].set_index(["frame", "ligand"])
        common = raw.index.intersection(fil.index)
        for m in metrics:
            if not len(common):
                out["per_tool"][t][m] = {"n": 0, "note": "no shared (frame,ligand)"}
                continue
            a = pd.to_numeric(raw.loc[common, m], errors="coerce")
            b = pd.to_numeric(fil.loc[common, m], errors="coerce")
            mask = a.notna().to_numpy() & b.notna().to_numpy()
            av, bv = a.to_numpy()[mask], b.to_numpy()[mask]
            rec = {"n": int(mask.sum()),
                   "median_raw": _f(np.median(av)) if av.size else None,
                   "median_filtered": _f(np.median(bv)) if bv.size else None}
            if av.size >= _MIN_PAIRS_WILCOXON:
                try:
                    rb, p, npair = su.wilcoxon_rankbiserial(bv, av)   # filtered vs raw
                    rec.update({"rank_biserial": _f(rb), "p_raw": _f(p), "n_nonzero": int(npair)})
                    if rec["p_raw"] is not None:
                        entries.append((t, m, rec))
                except Exception as e:                  # pragma: no cover
                    rec["error"] = str(e)
            else:
                rec["note"] = "n too small — exploratory"
            out["per_tool"][t][m] = rec
    if entries:
        adj = su.holm([rec["p_raw"] for *_ , rec in entries])
        for (t, m, rec), pa in zip(entries, adj):
            rec["p_holm"] = float(pa); rec["star"] = su.p_stars(pa)
    return out


def write_cluster_quality_stats_txt(out_dir: Path, q: pd.DataFrame,
                                    paired: dict, filtering: dict) -> Path:
    """Human-readable companion to the quality figures (tests kept OFF the panels)."""
    L: List[str] = []
    L.append("Orai per-tool cluster QUALITY — inferential tests (kept off the figure panels)")
    L.append("=" * 78)
    for ps in ("raw", "filtered"):
        sub = q[q["pose_set"] == ps]
        if len(sub):
            L.append(f"  {ps:<9} pose set: {len(sub)} (frame,ligand,tool) rows, "
                     f"{sub['ligand'].nunique()} ligands, tools "
                     f"{sorted(sub['tool'].unique())}")
    L.append("Unit for every test = the (frame,ligand) pair (poses within a pair are")
    L.append("NOT treated as independent). At the Orai×JKU n these are EXPLORATORY.")
    L.append("")

    L.append("-- Per-tool differences per metric (Friedman/Wilcoxon, paired across tools; "
             f"'{paired.get('pose_set','raw')}' cloud) --")
    for key, short, _yl, hb in QUALITY_METRICS:
        rec = (paired.get("metrics") or {}).get(key, {})
        meds = rec.get("medians", {})
        med_s = ", ".join(f"{_style(t)[0]}={meds[t]:.2f}" for t in meds
                          if meds.get(t) is not None) or "n/a"
        cap = _paired_summary(rec) if not rec.get("note") else rec.get("note")
        L.append(f"  {short:<20} ({'higher' if hb else 'lower'} = better)  "
                 f"n_complete={rec.get('n_pairs_complete', 0)}")
        L.append(f"      medians: {med_s}")
        L.append(f"      test: {cap or 'n/a'}")
    L.append("")

    L.append("-- Does filtering (pore ∩ PB-valid) tighten each tool's clusters? "
             "(Wilcoxon, paired filtered vs raw; Holm across tool×metric) --")
    L.append("   compactness rank-biserial<0 ⇒ tighter after filtering; "
             "silhouette/stability>0 ⇒ better.")
    for t, per_m in (filtering.get("per_tool") or {}).items():
        L.append(f"  {_style(t)[0]}:")
        for m in ("compactness", "silhouette", "stability"):
            rec = per_m.get(m, {})
            mr, mf = rec.get("median_raw"), rec.get("median_filtered")
            base = (f"raw={mr:.2f}→filt={mf:.2f}" if mr is not None and mf is not None
                    else "n/a")
            if rec.get("p_raw") is not None:
                pv = rec.get("p_holm", rec["p_raw"])
                tail = (f"  rb={rec['rank_biserial']:+.2f}  p_holm={su.fmt_p(pv)} "
                        f"{su.p_stars(pv)}  n={rec['n']}")
            else:
                tail = f"  ({rec.get('note', 'n/a')}; n={rec.get('n', 0)})"
            L.append(f"      {m:<12} {base}{tail}")
    L.append("")
    L.append("Panels carry descriptive content only (boxes, medians); the inferential")
    L.append("tests above are intentionally kept off the figures.")

    p = out_dir / "orai_cluster_quality_stats.txt"
    p.write_text("\n".join(L) + "\n")
    return p


def run_cluster_quality(args, ids, out_dir: Path, cluster_stats: dict) -> None:
    """Load the raw + pore pose CSVs, score reference-free per-tool cluster quality
    on both, write the table/summary/figures, and stash the tests in the sidecar.
    Deliberately independent of --pb-valid-only/--valid-ligands-only: the raw-vs-
    filtered comparison needs the UNFILTERED pose population."""
    raw_csv, pore_csv = _raw_and_pore_csv(args.per_pose_csv, args.raw_per_pose_csv)
    print("\n" + "=" * 86)
    print("CLUSTER-QUALITY BLOCK (reference-free intrinsic quality + raw-vs-filtered)")
    print("=" * 86)
    print(f"  raw  pose set : {raw_csv}{'' if raw_csv.exists() else '   (MISSING)'}")
    print(f"  pore pose set : {pore_csv}{'' if pore_csv.exists() else '   (MISSING)'}")

    def _load(csv: Path) -> pd.DataFrame:
        if not csv.exists():
            return pd.DataFrame()
        d = load_poses(csv, ids, args.diffdock_variant)
        return add_centroids(d, args.workers) if not d.empty else d

    rows: List[dict] = []
    raw_df = _load(raw_csv)
    if not raw_df.empty:
        print(f"  raw : {len(raw_df):,} poses | tools {sorted(raw_df['tool'].unique())}")
        rows += compute_tool_cluster_quality(
            raw_df, "raw", args.site_cluster, args.pocket_radius,
            args.quality_min_poses, args.quality_boot, args.quality_cap_n)
    pore_df = _load(pore_csv)
    filt_df = (pore_df[pore_df["pb_valid"] == True].copy()      # noqa: E712
               if not pore_df.empty else pd.DataFrame())
    if not filt_df.empty:
        print(f"  filt: {len(filt_df):,} poses (pore ∩ PB-valid) | "
              f"tools {sorted(filt_df['tool'].unique())}")
        rows += compute_tool_cluster_quality(
            filt_df, "filtered", args.site_cluster, args.pocket_radius,
            args.quality_min_poses, args.quality_boot, args.quality_cap_n)

    if not rows:
        print("  No poses available for the cluster-quality block — skipped.")
        return
    q = pd.DataFrame(rows)
    q.to_csv(out_dir / "cluster_quality_per_tool.csv", index=False)

    # per-(tool, pose_set) aggregate — quick reference + the cross-dataset overlay.
    agg = []
    for (tool, ps), g in q.groupby(["tool", "pose_set"]):
        rec = {"tool": tool, "pose_set": ps, "n_rows": int(len(g)),
               "mean_n_poses": round(float(g["n_poses"].mean()), 2)}
        for key, *_ in QUALITY_METRICS:
            vals = pd.to_numeric(g[key], errors="coerce").dropna()
            rec[f"median_{key}"] = round(float(vals.median()), 3) if len(vals) else np.nan
            rec[f"n_{key}"] = int(len(vals))
        agg.append(rec)
    pd.DataFrame(agg).to_csv(out_dir / "cluster_quality_summary.csv", index=False)

    # stats (guarded) + human-readable txt + sidecar keys
    paired = _stats_cluster_quality_paired(q, pose_set="raw")
    filtering = _stats_cluster_quality_filtering(q)
    cluster_stats["cluster_quality_per_tool"] = paired
    cluster_stats["cluster_quality_filtering"] = filtering
    txt = write_cluster_quality_stats_txt(out_dir, q, paired, filtering)

    if not args.no_plot:
        for f in (fig_cluster_quality(q, out_dir, "raw"),
                  fig_cluster_compactness_separation(q, out_dir, "raw"),
                  fig_cluster_quality_filtering(q, out_dir)):
            if f:
                print(f"  Figure: {f}")
    print(f"  Quality table : {out_dir / 'cluster_quality_per_tool.csv'}")
    print(f"  Quality stats : {txt}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/orai_benchmark/dock/posebusters_filtered_results.csv")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--ids-file", default="Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt")
    ap.add_argument("--match-thr", type=float, default=5.0,
                    help="Å threshold for tool agreement / pocket match.")
    ap.add_argument("--site-cluster", choices=("threshold", "gap", "silhouette"),
                    default="threshold",
                    help="Site clustering: 'threshold' (default, complete-linkage cut "
                         "at --pocket-radius, k>=1) | 'gap' | 'silhouette' (legacy).")
    ap.add_argument("--pocket-radius", type=float, default=8.0,
                    help="Distance-threshold (Å) for the 'threshold' site cut.")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="Cap #pairs (smoke test).")
    ap.add_argument("--out-dir", default="posebusters_results/orai_benchmark/pose_clusters")
    ap.add_argument("--features-csv",
                    default="PoseBusters_Benchmark_Analysis/ligand_protein_features.csv",
                    help="Per-ligand RDKit descriptors (from PoseBusters_DataSet_Analysis.ipynb); "
                         "enables the 'which ligand types dock consistently' analysis.")
    ap.add_argument("--diffdock-variant", default="all",
                    choices=("all", "original", "smina", "gnina"),
                    help="Restrict DiffDock to one optimizer variant (default 'all' "
                         "pools raw+smina+gnina). 'gnina' = keep only gnina-refined poses.")
    ap.add_argument("--valid-ligands-only", action="store_true",
                    help="Restrict the analysis to ligands that produced at least ONE "
                         "PoseBusters-valid pose anywhere (any Orai frame / any tool). "
                         "All poses of those ligands are kept (so the valid-vs-invalid "
                         "spatial split is preserved); ligands whose every docked pose "
                         "fails PoseBusters are dropped entirely before clustering.")
    ap.add_argument("--pb-valid-only", action="store_true",
                    help="Cluster ONLY the individual poses that pass PoseBusters "
                         "(pb_valid == True); every failing pose is dropped before "
                         "clustering. The site/agreement/descriptor axes then describe "
                         "the PB-valid pose cloud only. The valid-vs-invalid spatial "
                         "comparison and per-frame validity panels become degenerate "
                         "(everything is valid by construction) and are labelled as such.")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--cluster-quality", action="store_true",
                    help="Additionally compute reference-free per-tool cluster-QUALITY "
                         "metrics (silhouette / Calinski–Harabasz / Davies–Bouldin, "
                         "compactness vs separation, bootstrap stability) and a raw-vs-"
                         "filtered (all poses vs pore∩PB-valid) tightening analysis. "
                         "Writes cluster_quality_per_tool.csv, cluster_quality_summary.csv, "
                         "orai_cluster_quality_*.png and orai_cluster_quality_stats.txt.")
    ap.add_argument("--raw-per-pose-csv", default=None,
                    help="Unfiltered pose CSV (every generated pose incl. transmembrane + "
                         "PB-invalid) for the raw-vs-filtered comparison; default derives it "
                         "from --per-pose-csv by stripping the '.no_tm' suffix.")
    ap.add_argument("--quality-min-poses", type=int, default=3,
                    help="Min poses per (frame,ligand,tool) to score internal indices "
                         "(CH/DB need k>=2, so >=3 poses).")
    ap.add_argument("--quality-boot", type=int, default=25,
                    help="Bootstrap iterations for the dominant-cluster Jaccard stability.")
    ap.add_argument("--quality-cap-n", type=int, default=0,
                    help="If >0, subsample each tool's poses to this cap before scoring "
                         "quality, to remove the pose-count confound (CH/DB/silhouette are "
                         "n-sensitive). Default 0 = use every pose.")
    ap.add_argument("--top-n-poses", type=int, default=0,
                    help="If >0, restrict EVERY analysis to each tool's top-N ranked poses "
                         "per (frame, ligand) — AutoDock=Vina mode, DiffDock=confidence, "
                         "EquiBind=gnina-energy rank. The cap is computed on the FULL "
                         "(TM-inclusive) pose set, so it is top-N first THEN the TM/PB-valid "
                         "filters. Use 10 to equalise tools whose raw pose budget differs "
                         "(e.g. Orai×Experimental AutoDock/EquiBind at 30/unit).")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    global _TOPN_ALLOW
    if args.top_n_poses and args.top_n_poses > 0:
        raw_csv, _pore_csv = _raw_and_pore_csv(args.per_pose_csv, args.raw_per_pose_csv)
        allow_src = raw_csv if Path(raw_csv).exists() else Path(args.per_pose_csv)
        _TOPN_ALLOW = top_n_allowlist(allow_src, args.top_n_poses)
        print(f"  --top-n-poses {args.top_n_poses}: capped to {len(_TOPN_ALLOW)} pose files "
              f"(top-{args.top_n_poses}/tool per frame×ligand, ranked on {allow_src.name})")
    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    print(f"Loading poses from {args.per_pose_csv} ...")
    df = load_poses(Path(args.per_pose_csv), ids, args.diffdock_variant)
    if df.empty:
        print("No poses found (check CSV path / ids).")
        return 1
    print(f"  {len(df):,} poses | tools: {sorted(df['tool'].unique())} | "
          f"frames: {df['frame'].nunique()} | ligands: {df['ligand'].nunique()}")

    # Optional: keep only ligands with >=1 PoseBusters-valid pose (any frame/tool).
    # Every pose of a qualifying ligand is retained, so the per-pair clustering and
    # the valid-vs-invalid spatial comparison still run; ligands that never produced
    # a single valid pose are removed wholesale.
    if args.valid_ligands_only:
        valid_ligs = set(df.loc[df["pb_valid"] == True, "ligand"].unique())   # noqa: E712
        n_before = df["ligand"].nunique()
        df = df[df["ligand"].isin(valid_ligs)].copy()
        print(f"  --valid-ligands-only: kept {len(valid_ligs)}/{n_before} ligands "
              f"with >=1 PB-valid pose -> {len(df):,} poses remain.")
        if df.empty:
            print("No poses remain after the valid-ligand filter.")
            return 1

    # Pose-level PB-valid restriction: cluster ONLY passing poses. Applied after the
    # optional ligand filter. Poses with an unknown verdict (pb_valid is NaN, e.g. a
    # tool with no PoseBusters columns) are dropped too, since "valid-only" cannot
    # vouch for them.
    if args.pb_valid_only:
        n_poses_before, n_lig_before = len(df), df["ligand"].nunique()
        df = df[df["pb_valid"] == True].copy()   # noqa: E712
        print(f"  --pb-valid-only: kept {len(df):,}/{n_poses_before:,} PB-valid poses "
              f"across {df['ligand'].nunique()}/{n_lig_before} ligands "
              f"(valid-vs-invalid panels are degenerate under this filter).")
        if df.empty:
            print("No PB-valid poses remain after the filter.")
            return 1
    print("Computing pose centroids ...")
    df = add_centroids(df, args.workers)

    fp_dir = Path(args.fpocket_dir)
    axis_cache: Dict[str, Optional[tuple]] = {}

    records: List[dict] = []
    pose_rows: List[dict] = []
    poses_by_pair: Dict[tuple, list] = {}
    groups = list(df.groupby(["frame", "ligand"]))
    if args.limit:
        groups = groups[:args.limit]
    print(f"Analysing {len(groups)} (frame, ligand) pairs ...")
    for (frame, ligand), sub in groups:
        if frame not in axis_cache:
            axis_cache[frame] = frame_axis(frame, fp_dir)
        rec, prows = analyze_pair(frame, ligand, sub, args.match_thr, axis_cache[frame],
                                  args.site_cluster, args.pocket_radius)
        records.append(rec)
        if prows:
            # attach centroids for the 3D example figure
            cents = sub.set_index("pose_file")["_cent"].to_dict()
            for pr in prows:
                pr["_cent"] = cents.get(pr["pose_file"])
            poses_by_pair[(frame, ligand)] = prows
            pose_rows.extend(prows)

    pairs = pd.DataFrame([r for r in records if not r.get("skipped")])
    poses = pd.DataFrame([{k: v for k, v in r.items() if k != "_cent"} for r in pose_rows])
    pairs.to_csv(out_dir / "per_pair.csv", index=False)
    poses.to_csv(out_dir / "per_pose.csv", index=False)

    # ── per-frame summary ────────────────────────────────────────────────
    frows = []
    for fr, g in pairs.groupby("frame"):
        pg = poses[(poses.frame == fr) & poses["pb_valid"].notna()]
        frows.append({
            "frame": fr, "n_pairs": int(len(g)),
            "mean_n_clusters": round(float(g["n_clusters"].mean()), 3),
            "mean_dominant_frac": round(float(g["dominant_cluster_frac"].mean()), 3),
            "mean_inter_tool_dist": round(float(pd.to_numeric(g["mean_inter_tool_dist"], errors="coerce").mean()), 3),
            "tools_agree_rate": round(float(g["tools_agree"].mean()), 3),
            "pb_valid_rate": round(float(pg["pb_valid"].mean()), 3) if len(pg) else np.nan,
        })
    df_frame = pd.DataFrame(frows)
    df_frame.to_csv(out_dir / "per_frame_summary.csv", index=False)

    # ── outliers: pairs where tools most disagree + most-outlying poses ──
    out_pairs = pairs.sort_values("max_inter_tool_dist", ascending=False).head(25)
    out_pairs.to_csv(out_dir / "top_tool_disagreement.csv", index=False)
    out_poses = poses[poses["is_outlier"]].sort_values("dist_to_consensus", ascending=False).head(50)
    out_poses.to_csv(out_dir / "outlier_poses.csv", index=False)

    # ── console insights ─────────────────────────────────────────────────
    print("\n" + "=" * 86)
    print(f"ORAI POSE-CLUSTER REPORT  (n={len(pairs)} pairs, thr={args.match_thr:g} Å)")
    print("=" * 86)
    print(f"  Tools analysed: {', '.join(sorted(poses['tool'].unique()))}")
    multi = pairs[pairs["n_tools"] >= 2]
    if len(multi):
        agree = float(multi["tools_agree"].mean())
        print(f"  Cross-tool agreement (all tools' consensus within {args.match_thr:g} Å): "
              f"{agree:.0%} of {len(multi)} multi-tool pairs")
        print(f"  Median inter-tool consensus distance: "
              f"{pd.to_numeric(multi['mean_inter_tool_dist'], errors='coerce').median():.1f} Å")
    print(f"  Mean clusters/pair: {pairs['n_clusters'].mean():.2f} | "
          f"mean dominant-cluster fraction: {pairs['dominant_cluster_frac'].mean():.0%}")

    pk = poses[poses["pb_valid"].notna()]
    if "consensus_multitool" in pk:                    # throughput-robust pairs only
        pk = pk[pk["consensus_multitool"] == True]     # noqa: E712
    if len(pk):
        mv = pd.to_numeric(pk[pk.pb_valid == True]["dist_to_consensus"], errors="coerce").median()  # noqa: E712
        mi = pd.to_numeric(pk[pk.pb_valid == False]["dist_to_consensus"], errors="coerce").median()  # noqa: E712
        print(f"  PoseBusters (multi-tool pairs): valid poses sit {mv:.1f} Å from consensus vs "
              f"{mi:.1f} Å for invalid (median) → "
              f"{'valid are tighter' if mv < mi else 'no tightness advantage'}")

    print("\n  Per Orai MD frame:")
    print(f"  {'frame':<24}{'pairs':>6}{'clusters':>10}{'dom.frac':>10}"
          f"{'agree%':>8}{'pb-valid%':>11}")
    for _, r in df_frame.iterrows():
        pv = f"{r['pb_valid_rate']*100:.0f}%" if r['pb_valid_rate'] == r['pb_valid_rate'] else "n/a"
        print(f"  {r['frame']:<24}{r['n_pairs']:>6}{r['mean_n_clusters']:>10.2f}"
              f"{r['mean_dominant_frac']:>10.0%}{r['tools_agree_rate']:>7.0%}{pv:>11}")

    print(f"\n  Outliers: {int(poses['is_outlier'].sum())} poses > {2*args.match_thr:g} Å from "
          f"consensus; top tool-disagreement pairs in top_tool_disagreement.csv")
    print(f"\n  CSVs written to: {out_dir}/")

    summary = {
        "n_pairs": int(len(pairs)), "match_thr_A": args.match_thr,
        "tools": sorted(poses["tool"].unique()),
        "cross_tool_agree_rate": (float(multi["tools_agree"].mean()) if len(multi) else None),
        "mean_clusters_per_pair": float(pairs["n_clusters"].mean()),
        "mean_dominant_fraction": float(pairs["dominant_cluster_frac"].mean()),
        "per_frame": df_frame.to_dict(orient="records"),
        "n_outlier_poses": int(poses["is_outlier"].sum()),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # Statistical-test sidecar (see STATISTICAL_VALIDATION_PLAN.md). Every test is
    # guarded + labelled EXPLORATORY; a failure here never blocks the figures. The
    # independent n varies by dataset (Orai×JKU is tiny, ~12 (frame,ligand) pairs;
    # Orai×benchmark is ~1200), so the note reports the ACTUAL n rather than assuming.
    _n_pairs = int(len(pairs))
    _n_ligs = int(pairs["ligand"].nunique()) if len(pairs) else 0
    _tiny = " (a tiny independent n)" if _n_pairs < 60 else ""
    cluster_stats: Dict[str, dict] = {
        "n_pairs": _n_pairs, "match_thr_A": float(args.match_thr),
        "tools": sorted(poses["tool"].unique()),
        "note": (f"n={_n_pairs} independent (frame×ligand) pairs over {_n_ligs} "
                 f"ligands{_tiny}; all tests are exploratory and never treat "
                 f"correlated poses within a pair as independent."),
    }

    if not args.no_plot:
        f1 = fig_overview(pairs, poses, args.match_thr, out_dir, stats_out=cluster_stats)
        print(f"  Figure: {f1}")
        fdiv = fig_divergence_matrix(pairs, poses_by_pair, args.match_thr, out_dir,
                                     stats_out=cluster_stats)
        if fdiv:
            print(f"  Figure: {fdiv}")
        f2 = fig_examples(records, poses_by_pair, args.match_thr, out_dir)
        if f2:
            print(f"  Figure: {f2}")
        f3 = fig_descriptor_quality(pairs, args.features_csv, out_dir, stats_out=cluster_stats)
        if f3:
            print(f"  Figure: {f3}")

    # Reference-free per-tool cluster-QUALITY block (+ raw-vs-filtered tightening).
    # Loads its own pose populations, so it runs even under --pb-valid-only; a
    # failure here never blocks the main report.
    if args.cluster_quality:
        try:
            run_cluster_quality(args, ids, out_dir, cluster_stats)
        except Exception as e:                          # pragma: no cover
            print(f"  [cluster-quality] failed: {e}")

    try:
        (out_dir / "orai_cluster_stats.json").write_text(
            json.dumps(cluster_stats, indent=2, default=str))
        print(f"  Stats sidecar: {out_dir / 'orai_cluster_stats.json'}")
    except Exception as e:                              # pragma: no cover
        print(f"  [stats] failed to write orai_cluster_stats.json: {e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
