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
)
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


def load_poses(csv: Path, ids: Optional[set]) -> pd.DataFrame:
    df = pd.read_csv(csv, low_memory=False)
    tool_col = "docking_method" if "docking_method" in df.columns else "method"
    df["tool"] = df[tool_col].astype(str).str.lower()
    # collapse any equibind variant label (equibind_guided / equibind_unguided_*) to "equibind"
    df.loc[df["tool"].str.startswith("equibind"), "tool"] = "equibind"
    df["frame"] = df["protein"].astype(str)
    df["ligand"] = df["ligand"].astype(str)
    df["pb_valid"] = _derive_pb_valid(df)
    keep = ["tool", "frame", "ligand", "pose_file", "pb_valid"]
    df = df[[c for c in keep if c in df.columns]].copy()
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

def _draw_cross_tool_agreement(ax, pairs, thr, tools):
    """(A) inter-tool consensus-site distance distribution."""
    pair_cols = [c for c in pairs.columns if c.endswith("_dist")
                 and c not in ("mean_inter_tool_dist", "max_inter_tool_dist")]
    for c in pair_cols:
        v = pd.to_numeric(pairs[c], errors="coerce").dropna().to_numpy()
        if v.size:
            ax.hist(v, bins=np.linspace(0, 40, 21), histtype="step", lw=2,
                    label=f"{c.replace('_dist','')} (med {np.median(v):.1f} Å)")
    ax.axvline(thr, color="k", ls="--", lw=0.8, label=f"agree ≤ {thr:g} Å")
    ax.set_title("Cross-tool agreement — distance between tools' consensus sites")
    ax.set_xlabel("Inter-tool consensus-site distance (Å)")
    ax.set_ylabel("Number of receptor-ligand pairs")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)


def _draw_pb_validity_per_frame(ax, poses, frames, tools):
    """(B) per-frame PoseBusters validity rate, one line per tool — frames are
    ordered MD snapshots, so the line traces validity along the trajectory."""
    d = poses[poses["pb_valid"].notna()]
    x = np.arange(len(frames))
    for t in tools:
        rates = []
        for fr in frames:
            sel = d[(d["frame"] == fr) & (d["tool"] == t)]
            rates.append(float(sel["pb_valid"].mean()) * 100 if len(sel) else np.nan)
        ax.plot(x, rates, marker="o", lw=2, markersize=6,
                color=_style(t)[1], label=_style(t)[0],
                markeredgecolor="black", markeredgewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("Orai1WT-", "") for f in frames], rotation=20, ha="right")
    ax.set_title("PoseBusters validity rate per Orai MD frame")
    ax.set_ylabel("Valid poses (percent of generated)")
    ax.set_ylim(0, 105); ax.legend(fontsize=8)
    ax.grid(alpha=0.3); ax.set_axisbelow(True)


def _draw_cluster_structure(ax, pairs, frames):
    """(C) per-frame cluster tightness: mean #clusters + mean dominant fraction,
    pooled over all tools (the per-tool breakdown is the next panel)."""
    x = np.arange(len(frames))
    ncl = [pd.to_numeric(pairs[pairs.frame == fr]["n_clusters"], errors="coerce").mean() for fr in frames]
    dom = [pd.to_numeric(pairs[pairs.frame == fr]["dominant_cluster_frac"], errors="coerce").mean() for fr in frames]
    ax.bar(x - 0.2, ncl, 0.4, color="#8172B3", label="mean number of clusters")
    ax2 = ax.twinx()
    ax2.bar(x + 0.2, dom, 0.4, color="#CCB974", label="mean dominant-cluster fraction")
    ax.set_xticks(x)
    ax.set_xticklabels([f.replace("Orai1WT-", "") for f in frames], rotation=20, ha="right")
    ax.set_title("Cluster structure per Orai MD frame")
    ax.set_ylabel("Mean number of distinct clusters")
    ax2.set_ylabel("Mean dominant-cluster fraction"); ax2.set_ylim(0, 1.05)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    ax.grid(alpha=0.25); ax.set_axisbelow(True)


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


def _draw_tool_cluster_contribution(ax, poses, frames, tools):
    """(new) per-tool contribution to the cluster structure of each MD frame:
    grouped bars = mean number of distinct clusters a tool occupies (left axis);
    markers = mean share of that tool's poses in the pair's dominant cluster
    (right axis). Together they show how much each docking tool drives the
    multi-cluster spread vs concentrating on the dominant site."""
    stats = _per_tool_cluster_stats(poses, tools)
    x = np.arange(len(frames))
    nt = max(len(tools), 1)
    width = 0.8 / nt
    ax2 = ax.twinx()
    for k, t in enumerate(tools):
        off = (k - (nt - 1) / 2) * width
        ncl = [stats.get((fr, t), (np.nan, np.nan, 0))[0] for fr in frames]
        dom = [stats.get((fr, t), (np.nan, np.nan, 0))[1] for fr in frames]
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


def _draw_tool_agreement_matrix(ax, pairs, tools, thr):
    """(new) heat-map of pairwise tool agreement on the binding site."""
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
                ax.text(j, i, f"{frac[i, j]*100:.0f}%\n(n={cnt[i, j]})",
                        ha="center", va="center", fontsize=9,
                        color="black" if frac[i, j] < 0.6 else "white")
            else:
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9, color="0.4")
    ax.set_title(f"How often do two tools agree on a site?\n(per-tool sites within {thr:g} Å; n = shared pairs)")
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Fraction of shared pairs in agreement")
    ax.grid(False)


def _draw_valid_vs_consensus(ax, poses, thr):
    """(D) PoseBusters survivors vs failures — distance to the pair consensus.
    Multi-tool pairs only: for single-tool pairs the consensus is just that
    tool's own cluster centre, so 'nearer the consensus' would be circular."""
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


def _draw_tool_dispersion(ax, pairs, tools):
    """(E) per-tool homogeneity: how dispersed each tool's own poses are."""
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


def fig_overview(pairs: pd.DataFrame, poses: pd.DataFrame, thr: float, out_dir: Path):
    """Combined 8-panel overview PLUS a standalone PNG per panel in ``panels/``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tools = sorted(poses["tool"].unique())
    frames = sorted(pairs["frame"].unique())

    # (name, draw-fn) — one entry per panel; the draw-fn takes only an Axes so the
    # same call renders the grid cell and the standalone figure.
    panels = [
        ("cross_tool_agreement",     lambda a: _draw_cross_tool_agreement(a, pairs, thr, tools)),
        ("pb_validity_per_frame",    lambda a: _draw_pb_validity_per_frame(a, poses, frames, tools)),
        ("cluster_structure",        lambda a: _draw_cluster_structure(a, pairs, frames)),
        ("tool_cluster_contribution", lambda a: _draw_tool_cluster_contribution(a, poses, frames, tools)),
        ("tool_agreement_matrix",    lambda a: _draw_tool_agreement_matrix(a, pairs, tools, thr)),
        ("valid_vs_consensus",       lambda a: _draw_valid_vs_consensus(a, poses, thr)),
        ("tool_dispersion",          lambda a: _draw_tool_dispersion(a, pairs, tools)),
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


def fig_divergence_matrix(pairs: pd.DataFrame, poses_by_pair: dict, thr: float, out_dir: Path):
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
    col_labels = [lab for *_, lab in _TOOL_PAIRS] + ["mean of\npresent pairs"]

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
    _draw_dist_matrix(axes[0], S, row_labels, col_labels, n_pairs, frames_seq, norm, cmap, vmax,
                      "Do the tools prefer the same site?\n"
                      "distance between each tool's DOMINANT (largest-cluster) site (Å)", True)
    im = _draw_dist_matrix(axes[1], D, row_labels, col_labels, n_pairs, frames_seq, norm, cmap, vmax,
                           "Do the tools' pose clouds ever touch?\n"
                           "closest distance between ANY two poses of the tools (Å)", False)
    for a, lab in zip(axes, "AB"):                         # panel labels clear of the 2-line titles
        a.text(-0.02, 1.17, f"({lab})", transform=a.transAxes,
               fontsize=13, fontweight="bold", va="bottom", ha="right")
    fig.suptitle("Cross-tool divergence across all Orai MD snapshots × JKU ligands  "
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
    fig.suptitle("Example pose clusterings — colour = cluster, marker = tool "
                 "(○ AutoDock / △ DiffDock / □ EquiBind), edge: black = PoseBusters-valid, "
                 "red = invalid, grey = unknown", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "orai_example_clusters.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def fig_descriptor_quality(pairs: pd.DataFrame, features_csv: str, out_dir: Path):
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

    # (C) Spearman of properties vs validity & vs agreement
    feats = [f for f in PCA_FEATURES if f in d.columns] + [c for c in ("PC1", "PC2") if c in d.columns]
    rows = []
    for f in feats:
        v = pd.to_numeric(d[f], errors="coerce")
        out = {}
        for tgt, kk in (("pb_valid_rate", "valid"), ("inter_tool_dist", "agree")):
            m = v.notna() & d[tgt].notna()
            out[kk] = spearmanr(v[m], d[tgt][m])[0] if m.sum() > 10 else np.nan
        rows.append((f, out["valid"], out["agree"]))
    rows.sort(key=lambda x: (x[1] if x[1] == x[1] else 0))
    y = np.arange(len(rows)); w = 0.4
    ax[1, 0].barh(y - w / 2, [r[1] for r in rows], w, color="#2C7FB8", label="vs validity (↑ = more valid)")
    ax[1, 0].barh(y + w / 2, [r[2] for r in rows], w, color="#D95F0E", label="vs inter-tool dist (↑ = more disagreement)")
    ax[1, 0].set_yticks(y); ax[1, 0].set_yticklabels([nice(f) for f, _, _ in rows], fontsize=7)
    ax[1, 0].axvline(0, color="k", lw=0.8); ax[1, 0].legend(fontsize=7)
    ax[1, 0].set_title("Which properties track Orai docking quality?")
    ax[1, 0].set_xlabel("Spearman correlation")

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
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    print(f"Loading poses from {args.per_pose_csv} ...")
    df = load_poses(Path(args.per_pose_csv), ids)
    if df.empty:
        print("No poses found (check CSV path / ids).")
        return 1
    print(f"  {len(df):,} poses | tools: {sorted(df['tool'].unique())} | "
          f"frames: {df['frame'].nunique()} | ligands: {df['ligand'].nunique()}")
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

    if not args.no_plot:
        f1 = fig_overview(pairs, poses, args.match_thr, out_dir)
        print(f"  Figure: {f1}")
        fdiv = fig_divergence_matrix(pairs, poses_by_pair, args.match_thr, out_dir)
        if fdiv:
            print(f"  Figure: {fdiv}")
        f2 = fig_examples(records, poses_by_pair, args.match_thr, out_dir)
        if f2:
            print(f"  Figure: {f2}")
        f3 = fig_descriptor_quality(pairs, args.features_csv, out_dir)
        if f3:
            print(f"  Figure: {f3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
