#!/usr/bin/env python3
"""Tie docked-pose clusters, fpocket/p2rank pockets, and the crystal pose together.

For each benchmark complex this clusters the docked poses from the three docking
tools (AutoDock Vina, DiffDock, EquiBind) — in the spirit of
``pose_clustering_v3.ipynb`` — to get candidate binding pockets, then compares
those clusters, the fpocket/p2rank predicted pockets, and the experimental
crystal ligand in one 3D frame.

Key modeling fact: fpocket/p2rank pockets have **no ligand atoms**, so the only
quantity comparable across docked-pose clusters, predicted pockets, and the
crystal is the **3D center**. Centroid distance is therefore the primary, common
"how far off" axis; ligand RMSD-to-crystal (real poses only) is reported as a
complementary pose-quality axis.

Inputs (all already on disk — nothing re-docked):
  * per-pose CSV  : posebusters_results/.../pose_comparison_report/per_pose_metrics.csv
                    (method, protein, ligand, pose_file SDF, rank, rmsd, centroid_dist)
  * crystal pose  : Data/PoseBuster Benchmark Set/<id>/<id>_ligand.sdf
  * pockets       : pocket_results/{fpocket_results,p2rank_results}  (reuses the
                    parsers in pocket_comparison_report.py)

EquiBind is collapsed to one representative variant (``--equibind-variant``).

Run in the analysis env (conda env ``vina`` — rdkit + sklearn + scipy + matplotlib;
do NOT rely on sklearn_extra, it is broken under NumPy 2.x):

    python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
        --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
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
# Silence RDKit's harmless "tagged as 2D but has non-zero Z" warning before any
# molecule is read. Installed at import so it also covers the ProcessPool workers
# that compute pose centroids (fork inherits it; spawn re-imports this module).
# Also covers orai_pose_cluster_report.py, which imports this module.
from rdkit_quiet import silence_rdkit_2d3d_warning  # noqa: E402
silence_rdkit_2d3d_warning()
# Reuse the fpocket/p2rank parsers the user asked to compare against.
from pocket_comparison_report import (          # noqa: E402
    parse_fpocket, parse_p2rank, _load_ids, _label_panels,
)

from rdkit import Chem                            # noqa: E402
from rdkit.Chem import rdMolAlign                 # noqa: E402

TOOLS = ("autodock", "diffdock", "equibind")
TOOL_COLORS = {"autodock": "#4C72B0", "diffdock": "#55A868", "equibind": "#C44E52"}
TOOL_MARKERS = {"autodock": "o", "diffdock": "^", "equibind": "s"}


# ════════════════════════════════════════════════════════════════════════
# Ported pose_clustering_v3 helpers (notebook isn't importable)
# ════════════════════════════════════════════════════════════════════════

def load_heavy_atom_mol(sdf_path: str) -> Optional[Chem.Mol]:
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
            if mol is not None:
                mol = Chem.RemoveHs(mol)
        if mol is None or mol.GetNumConformers() == 0:
            return None
        return mol
    except Exception:
        return None


def heavy_atom_coords(mol: Chem.Mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                      conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])


def centroid_from_mol(mol: Chem.Mol) -> np.ndarray:
    return heavy_atom_coords(mol).mean(axis=0)


def compute_heavy_atom_rmsd(mol_a: Chem.Mol, mol_b: Chem.Mol,
                            max_matches: int = 1000) -> float:
    """Best-fit heavy-atom RMSD (symmetry-aware); centroid-distance fallback."""
    try:
        if mol_a.GetNumAtoms() == mol_b.GetNumAtoms():
            return float(rdMolAlign.GetBestRMS(mol_a, mol_b, maxMatches=max_matches))
        return float(np.linalg.norm(centroid_from_mol(mol_a) - centroid_from_mol(mol_b)))
    except Exception:
        try:
            return float(np.linalg.norm(centroid_from_mol(mol_a) - centroid_from_mol(mol_b)))
        except Exception:
            return float("nan")


class SimpleKMedoids:
    """PAM-style K-Medoids on a precomputed distance matrix (no sklearn_extra)."""

    def __init__(self, n_clusters=3, max_iter=300, random_state=42):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.random_state = random_state
        self.labels_ = None

    def fit(self, dm: np.ndarray):
        rng = np.random.RandomState(self.random_state)
        n = dm.shape[0]
        meds = rng.choice(n, self.n_clusters, replace=False)
        for _ in range(self.max_iter):
            labels = np.argmin(dm[:, meds], axis=1)
            new = meds.copy()
            for k in range(self.n_clusters):
                members = np.where(labels == k)[0]
                if len(members) == 0:
                    continue
                sub = dm[np.ix_(members, members)]
                new[k] = members[np.argmin(sub.sum(axis=1))]
            if np.array_equal(meds, new):
                break
            meds = new
        self.labels_ = np.argmin(dm[:, meds], axis=1)
        return self


def _normalise(dm: np.ndarray) -> np.ndarray:
    finite = dm[np.isfinite(dm)]
    if len(finite) == 0:
        return np.zeros_like(dm)
    lo, hi = finite.min(), finite.max()
    if hi - lo < 1e-12:
        return np.zeros_like(dm)
    out = (dm - lo) / (hi - lo)
    out[~np.isfinite(out)] = 1.0
    return out


# ════════════════════════════════════════════════════════════════════════
# Clustering (silhouette-selected k; KMedoids primary, Agglomerative cross-check)
# ════════════════════════════════════════════════════════════════════════

def _select_k(dm: np.ndarray, engine: str, max_k: int = 10):
    """Return (labels, k, silhouette) choosing k by max silhouette."""
    from sklearn.metrics import silhouette_score
    from sklearn.cluster import AgglomerativeClustering
    n = dm.shape[0]
    if n < 3:
        return np.zeros(n, dtype=int), 1, float("nan")
    best = (None, 1, -1.0)
    for k in range(2, min(max_k, n - 1) + 1):
        if engine == "kmedoids":
            labels = SimpleKMedoids(n_clusters=k).fit(dm).labels_
        else:
            labels = AgglomerativeClustering(
                n_clusters=k, metric="precomputed", linkage="average").fit_predict(dm)
        if len(set(labels)) < 2:
            continue
        try:
            sil = silhouette_score(dm, labels, metric="precomputed")
        except Exception:
            continue
        if sil > best[2]:
            best = (labels, k, float(sil))
    if best[0] is None:
        return np.zeros(n, dtype=int), 1, float("nan")
    return best


def _centroid_dm(cents: np.ndarray) -> np.ndarray:
    from scipy.spatial.distance import pdist, squareform
    return squareform(pdist(cents, metric="euclidean"))


def pockets_from_labels(cents: np.ndarray, labels: np.ndarray,
                        tools: List[str]) -> List[dict]:
    """Cluster -> candidate-pocket records, ranked by size (n poses)."""
    out = []
    for lab in sorted(set(labels.tolist())):
        idx = np.where(labels == lab)[0]
        c = cents[idx]
        center = c.mean(axis=0)
        radius = float(np.linalg.norm(c - center, axis=1).max()) if len(c) > 1 else 0.0
        tl = [tools[i] for i in idx]
        out.append({
            "center": center, "size": int(len(idx)), "radius": round(radius, 2),
            "tools": sorted(set(tl)), "n_tools": len(set(tl)),
            "tool_counts": dict(Counter(tl)),
            "label": int(lab), "members": idx,
        })
    out.sort(key=lambda p: p["size"], reverse=True)
    for i, p in enumerate(out, 1):
        p["rank"] = i
    return out


# Tool-combination ensembles to explore (singletons = baselines; combos = ensembles).
ENSEMBLES = {
    "AD": ("autodock",), "DD": ("diffdock",), "EB": ("equibind",),
    "AD+DD": ("autodock", "diffdock"), "AD+EB": ("autodock", "equibind"),
    "DD+EB": ("diffdock", "equibind"),
    "AD+DD+EB": ("autodock", "diffdock", "equibind"),
}


def _ensemble_stats(C: np.ndarray, tools: List[str], crystal: Optional[np.ndarray],
                    subset: Tuple[str, ...], thr: float) -> dict:
    """For one tool-combination: pose-level oracle + cluster top1 vs crystal."""
    idx = [i for i, t in enumerate(tools) if t in subset]
    if not idx or crystal is None:
        return dict(n_poses=len(idx), oracle_dist=np.nan, oracle_hit=np.nan,
                    top1_dist=np.nan, top1_hit=np.nan, n_clusters=np.nan)
    sc = C[idx]
    pose_d = np.linalg.norm(sc - crystal, axis=1)
    oracle = float(pose_d.min())
    # cluster the subset -> "top1" = largest cluster (what you'd actually pick)
    if len(sc) < 3:
        top1_center = sc.mean(axis=0)
        n_clusters = 1
    else:
        lab, _, _ = _select_k(_centroid_dm(sc), "kmedoids")
        big = Counter(lab.tolist()).most_common(1)[0][0]
        top1_center = sc[np.where(lab == big)[0]].mean(axis=0)
        n_clusters = len(set(lab.tolist()))
    top1 = _dist(top1_center, crystal)
    return dict(n_poses=len(idx), oracle_dist=round(oracle, 3),
                oracle_hit=bool(oracle <= thr), top1_dist=round(top1, 3),
                top1_hit=bool(top1 <= thr), n_clusters=int(n_clusters))


# ════════════════════════════════════════════════════════════════════════
# IO
# ════════════════════════════════════════════════════════════════════════

def _map_tool(method: str, eq_variant: str, dd_variant: str = "diffdock") -> Optional[str]:
    if method == "autodock":
        return "autodock"
    if method == dd_variant:                 # diffdock | diffdock_smina | diffdock_gnina
        return "diffdock"
    if method == eq_variant:
        return "equibind"
    return None


def load_poses(csv: Path, eq_variant: str, ids: Optional[set],
               pb_valid_only: bool, dd_variant: str = "diffdock") -> pd.DataFrame:
    df = pd.read_csv(csv, low_memory=False)
    df["tool"] = df["method"].map(lambda m: _map_tool(m, eq_variant, dd_variant))
    df = df[df["tool"].notna()].copy()
    # DiffDock writes its top pose twice — a bare ``rank1.sdf`` that duplicates
    # ``rank1_confidence-*.sdf`` (and their optimised copies rank1_<tool>.sdf).
    # Drop the bare copy so it doesn't double-count the rank-1 pose (which
    # otherwise makes rank1 == rank2 in per-rank stats).
    bare = (df["tool"].eq("diffdock")
            & df["pose_file"].astype(str).str.contains(
                r"rank\d+(?:_(?:smina|gnina))?\.sdf$", regex=True))
    df = df[~bare].copy()
    if pb_valid_only and "pb_valid" in df.columns:
        df = df[df["pb_valid"].astype(str).str.lower().isin(("true", "1", "1.0"))].copy()
    if ids is not None:
        df = df[df["protein"].isin(ids)].copy()
    for col in ("rank", "rmsd", "centroid_dist"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def _effective_ranks(sub: pd.DataFrame) -> Tuple[Dict[str, int], Dict[str, str]]:
    """Per-pose effective rank within each tool of a complex, + the source used.

    A tool's native ``rank`` is used when it varies (AutoDock = Vina mode,
    DiffDock = confidence rank). EquiBind has a constant sentinel rank (999) — it
    does not score poses — so we fall back to ascending ``smina_affinity`` (the
    only ordering the smina-refined variant has), else pose order. The source is
    recorded so the report can flag proxy rankings.
    """
    eff: Dict[str, int] = {}
    src: Dict[str, str] = {}
    for tool, g in sub.groupby("tool"):
        r = pd.to_numeric(g.get("rank"), errors="coerce")
        if r.notna().any() and r.nunique(dropna=True) > 1:
            order, how = r, "native"
        elif "smina_affinity" in g and pd.to_numeric(
                g["smina_affinity"], errors="coerce").nunique(dropna=True) > 1:
            order, how = pd.to_numeric(g["smina_affinity"], errors="coerce"), "smina_affinity"
        else:
            order, how = pd.Series(range(len(g)), index=g.index), "pose_order"
        er = order.rank(method="min", ascending=True)
        for f, e in zip(g["pose_file"], er):
            eff[f] = int(e) if e == e else None
        src[tool] = how
    return eff, src


def _extract_centroid(path: str):
    """Worker: SDF path -> (path, cx, cy, cz, n_atoms) or (path, None, ...)."""
    mol = load_heavy_atom_mol(path)
    if mol is None:
        return path, None, None, None, 0
    c = centroid_from_mol(mol)
    return path, float(c[0]), float(c[1]), float(c[2]), mol.GetNumAtoms()


def crystal_centroid(benchmark_dir: Path, cid: str) -> Optional[np.ndarray]:
    sdf = benchmark_dir / cid / f"{cid}_ligand.sdf"
    if not sdf.exists():
        return None
    mol = load_heavy_atom_mol(str(sdf))
    if mol is None:
        return None
    return centroid_from_mol(mol)


# ════════════════════════════════════════════════════════════════════════
# Per-complex analysis
# ════════════════════════════════════════════════════════════════════════

def _dist(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def _pocket_crystal_stats(pockets: List[dict], crystal: np.ndarray,
                          thr: float, key="center") -> dict:
    """top1 (rank-1) and oracle (best) center-to-crystal stats for an ordered list."""
    if not pockets or crystal is None:
        return dict(top1_dist=np.nan, top1_hit=np.nan,
                    oracle_dist=np.nan, oracle_rank=np.nan, oracle_hit=np.nan)
    dists = [_dist(p[key], crystal) for p in pockets]
    top1 = dists[0]
    o_i = int(np.argmin(dists))
    return dict(
        top1_dist=round(top1, 3), top1_hit=bool(top1 <= thr),
        oracle_dist=round(dists[o_i], 3),
        oracle_rank=int(pockets[o_i].get("rank", o_i + 1)),
        oracle_hit=bool(dists[o_i] <= thr),
    )


def _tool_pose_stats(sub: pd.DataFrame, crystal: Optional[np.ndarray],
                     cents: Dict[str, Tuple[float, float, float]],
                     thr: float) -> dict:
    """Per-tool pose-level stats vs crystal: top1/oracle centroid dist + RMSD."""
    if sub.empty:
        return {}
    sub = sub.sort_values("rank", na_position="last")
    # centroid distance to crystal from our own centroids (fallback to CSV col)
    def cdist(row):
        c = cents.get(row["pose_file"])
        if c is not None and crystal is not None:
            return _dist(c, crystal)
        return float(row.get("centroid_dist", np.nan))
    cd = sub.apply(cdist, axis=1).to_numpy(dtype=float)
    rm = pd.to_numeric(sub.get("rmsd"), errors="coerce").to_numpy(dtype=float) \
        if "rmsd" in sub else np.full(len(sub), np.nan)
    valid_cd = cd[np.isfinite(cd)]
    valid_rm = rm[np.isfinite(rm)]
    return dict(
        n_poses=int(len(sub)),
        top1_centroid_dist=round(float(cd[0]), 3) if np.isfinite(cd[0]) else np.nan,
        oracle_centroid_dist=round(float(valid_cd.min()), 3) if valid_cd.size else np.nan,
        oracle_centroid_hit=bool(valid_cd.min() <= thr) if valid_cd.size else np.nan,
        top1_rmsd=round(float(rm[0]), 3) if np.isfinite(rm[0]) else np.nan,
        oracle_rmsd=round(float(valid_rm.min()), 3) if valid_rm.size else np.nan,
        rmsd_le2=bool(valid_rm.min() <= 2.0) if valid_rm.size else np.nan,
    )


def _top5_analysis(sub: pd.DataFrame, eff_rank: Dict[str, int],
                   file_idx: Dict[str, int], labels: np.ndarray, C: np.ndarray,
                   crystal: Optional[np.ndarray], max_rank: int = 5):
    """Per-tool profile of the top-`max_rank` ranked poses.

    Returns (per_rank_rows, tool_summary):
      * per_rank_rows: (tool, rank_pos, centroid_dist, rmsd) for the k-th best
        pose of each tool — the "how far is the rank-k pose from crystal" profile.
      * tool_summary[tool]: how concentrated the top-5 are (n distinct clusters,
        modal-cluster fraction) and which of the top-5 is closest to the crystal
        (top5_best_rank, top5_best_dist).
    """
    file_rmsd = (dict(zip(sub["pose_file"],
                          pd.to_numeric(sub.get("rmsd"), errors="coerce")))
                 if "rmsd" in sub else {})
    per_rank: List[tuple] = []
    tool_summary: Dict[str, dict] = {}
    for tool, g in sub.groupby("tool"):
        ranked = sorted(((eff_rank[f], f) for f in g["pose_file"]
                         if eff_rank.get(f) is not None), key=lambda x: x[0])
        top = ranked[:max_rank]
        dists, clusters = [], []
        for pos, (_rnk, f) in enumerate(top, 1):
            idx = file_idx.get(f)
            cd = _dist(C[idx], crystal) if (idx is not None and crystal is not None) else np.nan
            rm = float(file_rmsd.get(f, np.nan))
            per_rank.append((tool, pos, cd, rm))
            dists.append(cd)
            if idx is not None:
                clusters.append(int(labels[idx]))
        if clusters:
            cc = Counter(clusters)
            modal_frac = max(cc.values()) / len(clusters)
            n_clusters = len(cc)
        else:
            modal_frac = n_clusters = np.nan
        valid = [(p, dd) for p, dd in enumerate(dists, 1) if dd == dd]
        if valid:
            best_pos, best_dist = min(valid, key=lambda x: x[1])
        else:
            best_pos = best_dist = np.nan
        tool_summary[tool] = dict(
            top5_n=len(top),
            top5_n_clusters=n_clusters,
            top5_modal_frac=round(modal_frac, 3) if modal_frac == modal_frac else np.nan,
            top5_best_rank=best_pos,
            top5_best_dist=round(best_dist, 3) if best_dist == best_dist else np.nan,
        )
    return per_rank, tool_summary


def analyze_complex(cid: str, sub: pd.DataFrame,
                    cents: Dict[str, Tuple[float, float, float]],
                    crystal: Optional[np.ndarray],
                    fp_pockets: List[dict], pr_pockets: List[dict],
                    thr: float, do_hybrid: bool, hybrid_w: float,
                    rmsd_cap: int) -> dict:
    # ── pose centroid array (drop poses with no centroid) ────────────────
    rows = [(r["pose_file"], r["tool"]) for _, r in sub.iterrows()
            if cents.get(r["pose_file"]) is not None]
    if len(rows) < 2:
        return {"protein": cid, "n_poses": len(rows), "skipped": True}
    files = [f for f, _ in rows]
    tools = [t for _, t in rows]
    C = np.asarray([cents[f] for f in files])

    # ── primary: centroid clustering -> candidate pockets ────────────────
    dm = _centroid_dm(C)
    labels, k, sil = _select_k(dm, "kmedoids")
    hi_labels, hi_k, hi_sil = _select_k(dm, "agglomerative")
    pockets = pockets_from_labels(C, labels, tools)

    # ── rank of each tool's poses in the crystal-closest cluster ─────────
    # The cluster whose center is nearest the crystal is the "correct" site.
    # For each tool we report the best (lowest) rank it assigns to a pose that
    # landed in that cluster — i.e. does the tool prioritize its near-native pose?
    eff_rank, rank_src = _effective_ranks(sub)
    correct_dist = correct_is_hit = np.nan
    rank_in_correct: Dict[str, Optional[int]] = {}
    n_in_correct: Dict[str, int] = {}
    if crystal is not None and pockets:
        cp = min(pockets, key=lambda p: _dist(p["center"], crystal))
        correct_dist = round(_dist(cp["center"], crystal), 3)
        correct_is_hit = bool(correct_dist <= thr)
        for t in TOOLS:
            rk = [eff_rank.get(files[i]) for i in cp["members"]
                  if tools[i] == t and eff_rank.get(files[i]) is not None]
            rank_in_correct[t] = int(min(rk)) if rk else None
            n_in_correct[t] = len(rk)

    # ── top-5 ranked poses: per-rank distance + cluster consistency ──────
    file_idx = {f: i for i, f in enumerate(files)}
    per_rank_rows, top5 = _top5_analysis(sub, eff_rank, file_idx, labels, C, crystal)

    # ── tool-combination ensemble exploration (which tools to combine) ───
    ensembles = {name: _ensemble_stats(C, tools, crystal, subset, thr)
                 for name, subset in ENSEMBLES.items()}

    # ── secondary: hybrid clustering (bounded) ───────────────────────────
    hy_k = hy_sil = ari = nmi = np.nan
    if do_hybrid and len(files) <= rmsd_cap:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
        mols = [load_heavy_atom_mol(f) for f in files]
        n = len(files)
        rdm = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                if mols[i] is None or mols[j] is None:
                    rdm[i, j] = rdm[j, i] = np.nan
                else:
                    r = compute_heavy_atom_rmsd(mols[i], mols[j])
                    rdm[i, j] = rdm[j, i] = r
        hyb = hybrid_w * _normalise(rdm) + (1 - hybrid_w) * _normalise(dm)
        hyb[~np.isfinite(hyb)] = np.nanmax(hyb[np.isfinite(hyb)]) if np.isfinite(hyb).any() else 0.0
        hyl, hy_k, hy_sil = _select_k(hyb, "kmedoids")
        if len(set(labels)) > 1 and len(set(hyl)) > 1:
            ari = float(adjusted_rand_score(labels, hyl))
            nmi = float(normalized_mutual_info_score(labels, hyl))

    # ── per-tool primary site (largest cluster of that tool's centroids) ──
    tool_site: Dict[str, Optional[np.ndarray]] = {}
    for t in TOOLS:
        ti = [i for i, tt in enumerate(tools) if tt == t]
        if not ti:
            tool_site[t] = None
            continue
        ct = C[ti]
        if len(ct) < 3:
            tool_site[t] = ct.mean(axis=0)
        else:
            tl, _, _ = _select_k(_centroid_dm(ct), "kmedoids")
            big = Counter(tl.tolist()).most_common(1)[0][0]
            tool_site[t] = ct[np.where(tl == big)[0]].mean(axis=0)

    # ── crystal comparison ───────────────────────────────────────────────
    ens = _pocket_crystal_stats(pockets, crystal, thr)            # ensemble clusters
    fp = _pocket_crystal_stats(fp_pockets, crystal, thr)
    pr = _pocket_crystal_stats(pr_pockets, crystal, thr)
    per_tool = {t: _tool_pose_stats(sub[sub["tool"] == t], crystal, cents, thr)
                for t in TOOLS}

    # ── cross-tool agreement (top sites within thr of each other) ────────
    inter = {}
    sites = {t: s for t, s in tool_site.items() if s is not None}
    for a in range(len(TOOLS)):
        for b in range(a + 1, len(TOOLS)):
            ta, tb = TOOLS[a], TOOLS[b]
            if ta in sites and tb in sites:
                inter[f"{ta[:2]}_{tb[:2]}_dist"] = round(_dist(sites[ta], sites[tb]), 3)
    n_agree_pairs = sum(1 for v in inter.values() if v <= thr)
    all3 = (len(sites) == 3 and all(v <= thr for v in inter.values()))

    # ── docked clusters vs predicted pockets (within thr) ────────────────
    def match_frac(pock):
        if not pockets or not pock:
            return np.nan
        from scipy.optimize import linear_sum_assignment
        A = np.asarray([p["center"] for p in pockets])
        B = np.asarray([p["center"] for p in pock])
        D = np.linalg.norm(A[:, None] - B[None], axis=2)
        ri, ci = linear_sum_assignment(D)
        m = sum(1 for i, j in zip(ri, ci) if D[i, j] <= thr)
        return round(m / min(len(pockets), len(pock)), 3)
    cluster_fp = match_frac(fp_pockets)
    cluster_pr = match_frac(pr_pockets)

    # ── triangulation: who finds the true site (oracle within thr) ───────
    tri = {
        "docking_hit": bool(ens["oracle_hit"]) if ens["oracle_hit"] == ens["oracle_hit"] else False,
        "fpocket_hit": bool(fp["oracle_hit"]) if fp["oracle_hit"] == fp["oracle_hit"] else False,
        "p2rank_hit": bool(pr["oracle_hit"]) if pr["oracle_hit"] == pr["oracle_hit"] else False,
    }

    return {
        "protein": cid, "skipped": False,
        "n_poses": len(files), "n_tools": len(set(tools)),
        "n_clusters": len(pockets), "silhouette": round(sil, 3) if sil == sil else np.nan,
        "hier_n_clusters": int(len(set(hi_labels))), "hier_silhouette": round(hi_sil, 3) if hi_sil == hi_sil else np.nan,
        "hybrid_n_clusters": hy_k, "hybrid_silhouette": round(hy_sil, 3) if hy_sil == hy_sil else np.nan,
        "ari_centroid_vs_hybrid": round(ari, 3) if ari == ari else np.nan,
        "nmi_centroid_vs_hybrid": round(nmi, 3) if nmi == nmi else np.nan,
        "has_crystal": crystal is not None,
        # ensemble-cluster vs crystal
        "ens_top1_dist": ens["top1_dist"], "ens_top1_hit": ens["top1_hit"],
        "ens_oracle_dist": ens["oracle_dist"], "ens_oracle_rank": ens["oracle_rank"],
        "ens_oracle_hit": ens["oracle_hit"],
        # predicted pockets vs crystal
        "fpocket_top1_dist": fp["top1_dist"], "fpocket_oracle_dist": fp["oracle_dist"],
        "fpocket_oracle_rank": fp["oracle_rank"], "fpocket_oracle_hit": fp["oracle_hit"],
        "p2rank_top1_dist": pr["top1_dist"], "p2rank_oracle_dist": pr["oracle_dist"],
        "p2rank_oracle_rank": pr["oracle_rank"], "p2rank_oracle_hit": pr["oracle_hit"],
        # cross-tool agreement
        **inter, "n_agree_tool_pairs": n_agree_pairs, "all_three_agree": all3,
        # clusters vs pockets
        "cluster_fpocket_match_frac": cluster_fp, "cluster_p2rank_match_frac": cluster_pr,
        # triangulation
        **{f"tri_{k2}": v for k2, v in tri.items()},
        # per tool (flattened)
        **{f"{t}_{k2}": v for t in TOOLS for k2, v in (per_tool[t] or {}).items()},
        # rank of each tool's poses in the crystal-closest cluster
        "correct_cluster_dist": correct_dist, "correct_cluster_is_hit": correct_is_hit,
        **{f"{t}_rank_in_correct_cluster": rank_in_correct.get(t) for t in TOOLS},
        **{f"{t}_n_in_correct_cluster": n_in_correct.get(t, 0) for t in TOOLS},
        **{f"{t}_rank_source": rank_src.get(t) for t in TOOLS},
        # tool-combination ensembles (which tools combined; oracle/top1 vs crystal)
        **{f"ens[{name}]_oracle_dist": s["oracle_dist"] for name, s in ensembles.items()},
        **{f"ens[{name}]_oracle_hit": s["oracle_hit"] for name, s in ensembles.items()},
        **{f"ens[{name}]_top1_dist": s["top1_dist"] for name, s in ensembles.items()},
        **{f"ens[{name}]_top1_hit": s["top1_hit"] for name, s in ensembles.items()},
        # top-5 ranked-pose consistency per tool
        **{f"{t}_{k2}": v for t in TOOLS for k2, v in (top5.get(t) or {}).items()},
        "_pockets": pockets, "_centroids": C, "_labels": labels, "_tools": tools,
        "_fp": fp_pockets, "_pr": pr_pockets, "_crystal": crystal,
        "_rank_in_correct": rank_in_correct, "_correct_is_hit": correct_is_hit,
        "_rank_src": rank_src, "_ensembles": ensembles,
        "_per_rank": [(cid, t, rp, cd, rm) for (t, rp, cd, rm) in per_rank_rows],
    }


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _fig_examples(results, thr, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm
    ex = [r for r in results if not r.get("skipped") and r.get("_crystal") is not None][:6]
    if not ex:
        return None
    n = len(ex)
    ncol = 3
    nrow = (n + ncol - 1) // ncol
    fig = plt.figure(figsize=(6.2 * ncol, 5.4 * nrow))
    panel_axes = []
    for i, r in enumerate(ex):
        ax = fig.add_subplot(nrow, ncol, i + 1, projection="3d")
        panel_axes.append(ax)
        C, labels, tools = r["_centroids"], r["_labels"], r["_tools"]
        clusters = sorted(set(labels.tolist()))
        cmap = cm.get_cmap("tab10", max(len(clusters), 1))
        for j in range(len(C)):
            ax.scatter(*C[j], color=cmap(clusters.index(labels[j])),
                       marker=TOOL_MARKERS.get(tools[j], "o"), s=34,
                       edgecolors="k", linewidths=0.2, alpha=0.85)
        cr = r["_crystal"]
        ax.scatter(*cr, marker="X", s=240, color="black", label="crystal", zorder=10)
        for p in r["_fp"][:3]:
            ax.scatter(*p["center"], marker="D", s=70, facecolors="none",
                       edgecolors="magenta", linewidths=1.6)
        for p in r["_pr"][:3]:
            ax.scatter(*p["center"], marker="D", s=70, facecolors="none",
                       edgecolors="cyan", linewidths=1.6)
        ax.set_title(f"{r['protein']}\nens oracle→crystal "
                     f"{r['ens_oracle_dist']:.1f} Å", fontsize=9)
        ax.set_xticklabels([]); ax.set_yticklabels([]); ax.set_zticklabels([])
    _label_panels(panel_axes, fontsize=12)
    fig.suptitle("Docked-pose clusters (color), tool (marker o/△/□), "
                 "crystal ✕, fpocket ◇ magenta, p2rank ◇ cyan", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "example_3d_clusters.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_summary(df, thr, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[df["has_crystal"]]
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))

    # (0,0) ECDF of oracle center-to-crystal distance per source
    src_oracle = {
        "autodock": "autodock_oracle_centroid_dist",
        "diffdock": "diffdock_oracle_centroid_dist",
        "equibind": "equibind_oracle_centroid_dist",
        "ensemble": "ens_oracle_dist",
        "fpocket": "fpocket_oracle_dist",
        "p2rank": "p2rank_oracle_dist",
    }
    ecdf = {}
    for name, col in src_oracle.items():
        if col not in d:
            continue
        v = np.sort(pd.to_numeric(d[col], errors="coerce").dropna().to_numpy())
        if v.size:
            ecdf[name] = (v, np.linspace(0, 1, v.size))
            ax[0, 0].plot(v, ecdf[name][1], label=name, lw=2)
    ax[0, 0].axvline(thr, color="k", ls="--", lw=0.8)
    ax[0, 0].set_xlim(0, 40)
    ax[0, 0].set_title("Oracle center-to-crystal distance (ECDF)")
    ax[0, 0].set_xlabel("distance (Å)"); ax[0, 0].set_ylabel("fraction of complexes")
    ax[0, 0].legend(fontsize=8)
    # Detail inset: the curves bunch up near 0 Å — zoom the 0–8 Å region so the
    # docking sources (which nearly all sit there) are distinguishable.
    axin = ax[0, 0].inset_axes([0.46, 0.12, 0.5, 0.55])
    for name, (vx, vy) in ecdf.items():
        axin.plot(vx, vy, lw=1.6)
    axin.axvline(thr, color="k", ls="--", lw=0.8)
    axin.set_xlim(0, 8); axin.set_ylim(0, 1)
    axin.set_title("detail: 0–8 Å", fontsize=8)
    axin.tick_params(labelsize=7)
    ax[0, 0].indicate_inset_zoom(axin, edgecolor="grey")

    # (0,1) hit-rate@thr per source (oracle; + top1 where defined)
    def rate(col):
        v = pd.to_numeric(d[col], errors="coerce").dropna()
        return float((v <= thr).mean()) if len(v) else np.nan
    labels = list(src_oracle)
    oracle_rates = [rate(src_oracle[s]) for s in labels]
    x = np.arange(len(labels))
    ax[0, 1].bar(x, oracle_rates, color=[TOOL_COLORS.get(s, "#8172B3") for s in labels])
    for xi, v in zip(x, oracle_rates):
        if v == v:
            ax[0, 1].text(xi, v + 0.01, f"{v:.0%}", ha="center", fontsize=9)
    ax[0, 1].set_xticks(x); ax[0, 1].set_xticklabels(labels, rotation=30, ha="right")
    ax[0, 1].set_ylim(0, 1.18); ax[0, 1].set_ylabel("fraction within thr")
    ax[0, 1].set_title(f"Finds the true site (oracle ≤ {thr:g} Å)")

    # (1,0) triangulation: who finds the true site
    combos = Counter()
    for _, r in d.iterrows():
        key = (bool(r.get("tri_docking_hit")), bool(r.get("tri_fpocket_hit")),
               bool(r.get("tri_p2rank_hit")))
        combos[key] += 1
    names = {(1, 1, 1): "all three", (1, 1, 0): "dock+fp", (1, 0, 1): "dock+p2r",
             (0, 1, 1): "fp+p2r", (1, 0, 0): "dock only", (0, 1, 0): "fp only",
             (0, 0, 1): "p2r only", (0, 0, 0): "none"}
    keys = sorted(names, key=lambda k: -combos[tuple(map(bool, k))])
    vals = [combos[tuple(map(bool, k))] for k in keys]
    tot = sum(vals) or 1
    ax[1, 0].bar(range(len(keys)), [v / tot for v in vals], color="#55A868")
    ax[1, 0].set_xticks(range(len(keys)))
    ax[1, 0].set_xticklabels([names[k] for k in keys], rotation=30, ha="right", fontsize=8)
    ax[1, 0].set_title(f"Triangulation — who finds the true site (≤ {thr:g} Å)")
    ax[1, 0].set_ylabel("fraction of complexes")
    for xi, v in enumerate(vals):
        ax[1, 0].text(xi, v / tot + 0.005, str(v), ha="center", fontsize=8)

    # (1,1) cross-tool agreement: inter-tool top-site distances
    for col, name, c in [("au_di_dist", "AD–DD", "#4C72B0"),
                         ("au_eq_dist", "AD–EB", "#55A868"),
                         ("di_eq_dist", "DD–EB", "#C44E52")]:
        if col in d:
            v = pd.to_numeric(d[col], errors="coerce").dropna().to_numpy()
            if v.size:
                ax[1, 1].hist(v, bins=np.linspace(0, 40, 21), histtype="step",
                              lw=2, label=f"{name} (med {np.median(v):.1f})", color=c)
    ax[1, 1].axvline(thr, color="k", ls="--", lw=0.8)
    ax[1, 1].set_title("Cross-tool agreement — top-site center distance")
    ax[1, 1].set_xlabel("distance (Å)"); ax[1, 1].set_ylabel("complexes")
    ax[1, 1].legend(fontsize=8)

    _label_panels(ax)
    fig.suptitle("Docked clusters vs fpocket/p2rank vs crystal "
                 f"(n={len(d)} complexes, thr={thr:g} Å)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "crystal_pocket_summary.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_rank_in_correct(df, thr, out_dir):
    """How each tool ranks the poses it puts in the crystal-closest cluster."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[(df.get("correct_cluster_is_hit") == True)]  # noqa: E712
    if d.empty:
        return None
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    # (left) ECDF of best rank in the correct cluster, per tool
    for t in TOOLS:
        col = f"{t}_rank_in_correct_cluster"
        if col not in d:
            continue
        v = np.sort(pd.to_numeric(d[col], errors="coerce").dropna().to_numpy())
        if v.size:
            ax[0].plot(v, np.linspace(0, 1, v.size), marker=".", lw=2,
                       color=TOOL_COLORS[t],
                       label=f"{t} (med {np.median(v):.0f}, n={v.size})")
    ax[0].axvline(1, color="k", ls=":", lw=0.8)
    ax[0].set_xlim(0, 31)
    ax[0].set_title("Best rank of a tool's pose in the crystal-closest cluster (ECDF)")
    ax[0].set_xlabel("rank (1 = tool's top pose)")
    ax[0].set_ylabel("fraction of true-site complexes")
    ax[0].grid(True, ls="--", lw=0.5, alpha=0.5)
    ax[0].set_axisbelow(True)
    ax[0].legend(fontsize=8)
    # (right) coverage + %top1 + %top<=5 bars
    metrics = ["coverage", "%top-1", "%top-5"]
    x = np.arange(len(TOOLS)); w = 0.25
    for mi, m in enumerate(metrics):
        vals = []
        for t in TOOLS:
            col = f"{t}_rank_in_correct_cluster"
            r = pd.to_numeric(d[col], errors="coerce") if col in d else pd.Series(dtype=float)
            present = r.dropna()   # complexes where the tool has a pose in the cluster
            if m == "coverage":
                vals.append(float(r.notna().mean()) if len(r) else np.nan)
            elif m == "%top-1":
                vals.append(float((present == 1).mean()) if len(present) else np.nan)
            else:
                vals.append(float((present <= 5).mean()) if len(present) else np.nan)
        ax[1].bar(x + (mi - 1) * w, vals, w, label=m)
    ax[1].set_xticks(x); ax[1].set_xticklabels(TOOLS)
    ax[1].set_ylim(0, 1.1)
    ax[1].set_ylabel("fraction (coverage: all; %top: given present)")
    ax[1].set_title("Coverage of, and ranking within, the crystal-closest cluster")
    ax[1].grid(True, axis="y", ls="--", lw=0.5, alpha=0.5)
    ax[1].set_axisbelow(True)
    ax[1].legend(fontsize=8)
    _label_panels(ax)
    fig.suptitle("How each tool ranks its poses in the crystal-closest cluster "
                 f"(true-site complexes, n={len(d)}; EquiBind rank = smina-affinity proxy)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "rank_in_correct_cluster.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_ensembles(df, thr, out_dir):
    """Which tool combination best covers the true site."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[df["has_crystal"]]
    names = list(ENSEMBLES)
    colors = ["#999999"] * 3 + ["#4C72B0", "#55A868", "#C44E52", "#8172B3"]
    oracle_hit, top1_hit, med_oracle = [], [], []
    for nm in names:
        oh = pd.to_numeric(d.get(f"ens[{nm}]_oracle_hit").map(
            {True: 1, False: 0}) if f"ens[{nm}]_oracle_hit" in d else pd.Series(dtype=float),
            errors="coerce")
        th = pd.to_numeric(d.get(f"ens[{nm}]_top1_hit").map(
            {True: 1, False: 0}) if f"ens[{nm}]_top1_hit" in d else pd.Series(dtype=float),
            errors="coerce")
        od = pd.to_numeric(d.get(f"ens[{nm}]_oracle_dist"), errors="coerce") \
            if f"ens[{nm}]_oracle_dist" in d else pd.Series(dtype=float)
        oracle_hit.append(float(oh.mean()) if len(oh.dropna()) else np.nan)
        top1_hit.append(float(th.mean()) if len(th.dropna()) else np.nan)
        med_oracle.append(float(od.median()) if len(od.dropna()) else np.nan)
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(names))
    # Vertical dumbbell: filled = oracle, open = top-1; the connector length is
    # the ranking loss (oracle ceiling vs the largest-cluster pick).
    for xi, oh, th in zip(x, oracle_hit, top1_hit):
        if oh == oh and th == th:
            ax[0].plot([xi, xi], [min(oh, th), max(oh, th)], color="#b8b8b8",
                       lw=2.6, zorder=1, solid_capstyle="round")
    ax[0].scatter(x, oracle_hit, s=90, c=colors, edgecolor="black",
                  linewidth=0.7, zorder=3)
    ax[0].scatter(x, top1_hit, s=90, facecolors="white", edgecolors=colors,
                  linewidths=1.8, zorder=3)
    for xi, oh, th in zip(x, oracle_hit, top1_hit):
        if oh == oh:
            ax[0].text(xi, oh + 0.025, f"{oh:.0%}", ha="center", va="bottom",
                       fontsize=7, fontweight="bold")
        if th == th:
            ax[0].text(xi + 0.12, th, f"{th:.0%}", ha="left", va="center",
                       fontsize=6.5, color="#444")
    ax[0].set_xticks(x); ax[0].set_xticklabels(names, rotation=30, ha="right")
    ax[0].set_xlim(-0.6, len(names) - 0.4)
    ax[0].set_ylim(0, 1.18); ax[0].set_ylabel("fraction of complexes")
    ax[0].set_title(f"Ensemble finds the true site (≤ {thr:g} Å)")
    ax[0].grid(axis="y", alpha=0.3); ax[0].set_axisbelow(True)
    ax[0].legend(handles=[
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555",
               markeredgecolor="black", markersize=9, label="oracle (best pose)"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="#555555", markeredgewidth=1.8, markersize=9,
               label="top-1 (largest cluster)"),
    ], fontsize=8, loc="upper left")
    ax[1].bar(x, med_oracle, color=colors)
    for xi, v in zip(x, med_oracle):
        if v == v:
            ax[1].text(xi, v + 0.05, f"{v:.1f}", ha="center", fontsize=7)
    ax[1].set_xticks(x); ax[1].set_xticklabels(names, rotation=30, ha="right")
    ax[1].set_ylabel("median oracle distance (Å)")
    ax[1].set_title("Ensemble median oracle center-to-crystal distance")
    _label_panels(ax)
    fig.suptitle("Tool-combination ensembles  "
                 "(AD=AutoDock, DD=DiffDock, EB=EquiBind)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "ensemble_exploration.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


_SRC_ORACLE = {
    "autodock": "autodock_oracle_centroid_dist", "diffdock": "diffdock_oracle_centroid_dist",
    "equibind": "equibind_oracle_centroid_dist", "ensemble": "ens_oracle_dist",
    "fpocket": "fpocket_oracle_dist", "p2rank": "p2rank_oracle_dist",
}
_SRC_COLOR = {"autodock": "#4C72B0", "diffdock": "#55A868", "equibind": "#C44E52",
              "ensemble": "#8172B3", "fpocket": "#CCB974", "p2rank": "#64B5CD"}


def _ecdf_xy(v):
    v = np.sort(np.asarray(v, dtype=float))
    return v, np.linspace(0, 1, len(v)) if len(v) else (v, v)


def _fig_ecdf_split(df, thr, out_dir):
    """Split the combined oracle ECDF into one panel per source + an AD/DD crossover."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[df["has_crystal"]]
    series = {name: pd.to_numeric(d[col], errors="coerce").dropna().to_numpy()
              for name, col in _SRC_ORACLE.items() if col in d}
    items = list(series.items())
    npan = len(items) + 1                      # one per source + AD/DD crossover
    ncol = 3
    nrow = (npan + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 3.8 * nrow))
    axes = axes.flatten()
    for ax, (name, v) in zip(axes, items):
        if v.size:
            vx, vy = _ecdf_xy(v)
            ax.plot(vx, vy, color=_SRC_COLOR.get(name, "k"), lw=2.2)
            ax.fill_between(vx, vy, alpha=0.12, color=_SRC_COLOR.get(name, "k"))
        ax.axvline(thr, color="k", ls="--", lw=0.8)
        ax.set_xlim(0, 15); ax.set_ylim(0, 1.02)
        med = np.median(v) if v.size else np.nan
        hit = np.mean(v <= thr) if v.size else np.nan
        ax.set_title(f"{name}  (median {med:.2f} Å, ≤{thr:g}Å {hit:.0%}, n={v.size})",
                     fontsize=10)
        ax.set_xlabel("oracle dist (Å)"); ax.set_ylabel("frac ≤ x")
        ax.grid(alpha=0.25)
    # AD vs DD crossover panel
    axc = axes[len(items)]
    ad = series.get("autodock", np.array([])); dd = series.get("diffdock", np.array([]))
    if ad.size and dd.size:
        ax_x = np.linspace(0, 8, 401)
        ea = np.mean(ad[:, None] <= ax_x[None, :], axis=0)
        ed = np.mean(dd[:, None] <= ax_x[None, :], axis=0)
        axc.plot(ax_x, ea, color=_SRC_COLOR["autodock"], lw=2.2, label="autodock")
        axc.plot(ax_x, ed, color=_SRC_COLOR["diffdock"], lw=2.2, label="diffdock")
        diff = ea - ed
        cross = [ax_x[i] for i in range(1, len(ax_x)) if diff[i - 1] * diff[i] < 0]
        for c in cross:
            axc.axvline(c, color="grey", ls=":", lw=1)
            axc.annotate(f"cross {c:.1f} Å", (c, 0.05), fontsize=8, rotation=90,
                         color="grey", va="bottom")
        axc.axvline(thr, color="k", ls="--", lw=0.8)
        axc.set_xlim(0, 8); axc.set_ylim(0, 1.02)
        axc.set_title("AutoDock vs DiffDock (crossover)", fontsize=10)
        axc.set_xlabel("oracle dist (Å)"); axc.set_ylabel("frac ≤ x")
        axc.legend(fontsize=8); axc.grid(alpha=0.25)
    for ax in axes[npan:]:
        ax.set_visible(False)
    _label_panels(axes)
    fig.suptitle("Oracle center-to-crystal distance — ECDF per source "
                 f"(n={len(d)} complexes, thr={thr:g} Å). Steeper/left = more accurate.",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "oracle_ecdf_split.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_top5(df_rank, df, out_dir):
    """Per-rank distance profile (top-5) + top-5 cluster consistency, per tool."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if df_rank.empty:
        return None
    d = df[df["has_crystal"]]
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    ranks = list(range(1, 6))
    # (0,0) centroid distance vs rank ; (0,1) RMSD vs rank.
    # Median is the primary (solid) line — robust to EquiBind's far/garbage poses
    # that otherwise blow the scale; mean (dashed) is overlaid and clipped.
    caps = {"centroid_dist": 15.0, "rmsd": 15.0}
    for col_name, axi, lab in [("centroid_dist", ax[0, 0], "centroid distance"),
                               ("rmsd", ax[0, 1], "RMSD")]:
        cap = caps[col_name]
        for t in TOOLS:
            g = df_rank[df_rank.tool == t]
            means, meds = [], []
            for k in ranks:
                vk = pd.to_numeric(g[g["rank"] == k][col_name], errors="coerce").dropna()
                means.append(vk.mean() if len(vk) else np.nan)
                meds.append(vk.median() if len(vk) else np.nan)
            clip = lambda xs: [x if (x == x and x <= cap) else np.nan for x in xs]
            axi.plot(ranks, clip(meds), "-o", color=TOOL_COLORS[t], lw=2.2,
                     label=f"{t} (median)")
            axi.plot(ranks, clip(means), "--", color=TOOL_COLORS[t], lw=1, alpha=0.55)
        axi.set_xticks(ranks); axi.set_ylim(0, cap)
        axi.set_xlabel("pose rank (1 = tool's top pose)")
        axi.set_ylabel(f"{lab} to crystal (Å)")
        axi.set_title(f"Top-5 {lab} from crystal (solid=median, dashed=mean; ≤{cap:g} Å)")
        axi.legend(fontsize=8); axi.grid(alpha=0.25)
    # (1,0) top-5 cluster consistency: mean modal fraction + mean #clusters
    x = np.arange(len(TOOLS)); w = 0.38
    mf = [pd.to_numeric(d.get(f"{t}_top5_modal_frac"), errors="coerce").mean()
          if f"{t}_top5_modal_frac" in d else np.nan for t in TOOLS]
    nc = [pd.to_numeric(d.get(f"{t}_top5_n_clusters"), errors="coerce").mean()
          if f"{t}_top5_n_clusters" in d else np.nan for t in TOOLS]
    ax[1, 0].bar(x - w / 2, mf, w, label="modal-cluster fraction",
                 color=[TOOL_COLORS[t] for t in TOOLS])
    ax[1, 0].bar(x + w / 2, [v / 5 for v in nc], w, alpha=0.4,
                 color=[TOOL_COLORS[t] for t in TOOLS],
                 label="#distinct clusters ÷ 5")
    for xi, v in zip(x, mf):
        if v == v:
            ax[1, 0].text(xi - w / 2, v + 0.01, f"{v:.0%}", ha="center", fontsize=8)
    ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels(TOOLS); ax[1, 0].set_ylim(0, 1.1)
    ax[1, 0].set_title("Are a tool's top-5 in one cluster? (modal fraction)")
    ax[1, 0].set_ylabel("fraction"); ax[1, 0].legend(fontsize=8)
    # (1,1) which of the top-5 is closest to crystal (distribution of best rank)
    width = 0.25
    for ti, t in enumerate(TOOLS):
        br = pd.to_numeric(d.get(f"{t}_top5_best_rank"), errors="coerce").dropna() \
            if f"{t}_top5_best_rank" in d else pd.Series(dtype=float)
        counts = [float((br == k).mean()) if len(br) else 0 for k in ranks]
        ax[1, 1].bar(np.array(ranks) + (ti - 1) * width, counts, width,
                     color=TOOL_COLORS[t], label=t)
    ax[1, 1].set_xticks(ranks)
    ax[1, 1].set_xlabel("rank of the top-5 pose closest to crystal")
    ax[1, 1].set_ylabel("fraction of complexes")
    ax[1, 1].set_title("Which of the top-5 is closest to the crystal?")
    ax[1, 1].legend(fontsize=8)
    _label_panels(ax)
    fig.suptitle("Top-5 ranked poses: per-rank accuracy and cluster consistency",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "top5_rank_distance.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Driver
# ════════════════════════════════════════════════════════════════════════

def _fig_descriptor_quality(df_complex, features_csv, out_dir):
    """Which ligand types dock well? Pose accuracy (oracle RMSD to crystal) vs the
    ligand's Lipinski Ro5 compliance, physicochemical descriptors, and PCA space."""
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
        print("  (ligand feature CSV not found; skipping descriptor-quality figure)")
        return None
    d = df_complex[df_complex["has_crystal"]].copy() if "has_crystal" in df_complex else df_complex.copy()
    rcols = [f"{t}_oracle_rmsd" for t in TOOLS if f"{t}_oracle_rmsd" in d.columns]
    if not rcols or d.empty:
        return None
    for c in rcols:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["best_oracle_rmsd"] = d[rcols].min(axis=1)        # closest any tool got to crystal
    d = d.set_index("protein").join(desc, how="inner")
    d = d[np.isfinite(d["best_oracle_rmsd"])]
    if len(d) < 12:
        return None
    d["success"] = d["best_oracle_rmsd"] <= 2.0
    d, loadings, evr, _ = add_pca(d)
    bins, labs = [-0.1, 2, 5, 8, 11, 100], ["0–2", "3–5", "6–8", "9–11", "12+"]
    d["_rb"] = pd.cut(pd.to_numeric(d["rot_bonds"], errors="coerce"), bins=bins, labels=labs)
    d.reset_index().to_csv(out_dir / "descriptor_vs_quality.csv", index=False)

    fig, ax = plt.subplots(2, 3, figsize=(18, 11))

    # (A) accuracy by Lipinski Rule-of-Five
    grp = [("Ro5 pass\n(≤1 viol.)", d[d.ro5_pass]["best_oracle_rmsd"].dropna()),
           ("Ro5 fail\n(≥2 viol.)", d[~d.ro5_pass]["best_oracle_rmsd"].dropna())]
    ax[0, 0].boxplot([g.clip(upper=20) for _, g in grp], labels=[n for n, _ in grp], showfliers=False)
    for i, (_n, g) in enumerate(grp, 1):
        sr = (g <= 2).mean() * 100 if len(g) else np.nan
        ax[0, 0].text(i, 0.4, f"{sr:.0f}% ≤2Å\nn={len(g)}", ha="center", va="bottom", fontsize=8)
    ax[0, 0].axhline(2, color="green", ls="--", lw=1)
    ax[0, 0].set_title("Pose accuracy vs Lipinski \n Rule-of-Five compliance")
    ax[0, 0].set_ylabel("Best oracle RMSD to crystal (Å)")

    # (B) success vs flexibility (rotatable bonds) — a trend across ordinal bins,
    # so a line reads the monotonic decline more directly than bars.
    sr = d.groupby("_rb")["success"].agg(["mean", "size"])
    ax[0, 1].plot(range(len(sr)), sr["mean"] * 100, marker="o", color="#4C72B0",
                  lw=2, markersize=7, markeredgecolor="black", markeredgewidth=0.6)
    for i, (m, nn) in enumerate(zip(sr["mean"], sr["size"])):
        if nn:
            ax[0, 1].text(i, m * 100 + 4, f"{m*100:.0f}%\nn={int(nn)}",
                          ha="center", va="bottom", fontsize=8)
    ax[0, 1].set_xticks(range(len(sr))); ax[0, 1].set_xticklabels(sr.index)
    ax[0, 1].set_title("Docking success vs ligand flexibility")
    ax[0, 1].set_xlabel("Rotatable bonds")
    ax[0, 1].set_ylabel("Complexes with a pose ≤2 Å (percent)"); ax[0, 1].set_ylim(0, 105)
    ax[0, 1].grid(alpha=0.3); ax[0, 1].set_axisbelow(True)

    # (C) Spearman correlation of each property (+ PCs) with pose error
    feats = [f for f in PCA_FEATURES if f in d.columns] + [c for c in ("PC1", "PC2", "PC3") if c in d.columns]
    rhos = []
    for f in feats:
        v = pd.to_numeric(d[f], errors="coerce")
        m = v.notna() & d["best_oracle_rmsd"].notna()
        if m.sum() > 10:
            rho, _ = spearmanr(v[m], d["best_oracle_rmsd"][m])
            rhos.append((f, rho))
    rhos.sort(key=lambda x: x[1])
    vals = [r for _, r in rhos]
    ax[0, 2].barh(range(len(vals)), vals, color=["#C44E52" if r > 0 else "#4C72B0" for r in vals])
    ax[0, 2].set_yticks(range(len(vals))); ax[0, 2].set_yticklabels([nice(f) for f, _ in rhos], fontsize=7)
    ax[0, 2].axvline(0, color="k", lw=0.8)
    ax[0, 2].set_title("Which properties track pose error?\n(Spearman ρ vs RMSD\n red = worse, blue = better)")
    ax[0, 2].set_xlabel("Spearman correlation with best oracle RMSD")

    # (D) PCA chemical space coloured by accuracy
    if "PC1" in d and "PC2" in d:
        sc = ax[1, 0].scatter(d["PC1"], d["PC2"], c=d["best_oracle_rmsd"].clip(upper=10),
                              cmap="viridis_r", s=22, edgecolors="k", linewidths=0.2)
        plt.colorbar(sc, ax=ax[1, 0], label="Best oracle RMSD (Å, clipped 10)")
        e1 = evr[0] * 100 if len(evr) > 0 else 0
        e2 = evr[1] * 100 if len(evr) > 1 else 0
        ax[1, 0].set_xlabel(f"Principal Component 1 ({e1:.0f}% of variance)")
        ax[1, 0].set_ylabel(f"Principal Component 2 ({e2:.0f}% of variance)")
        ax[1, 0].set_title("Ligand chemical space (PCA) \ncoloured by pose accuracy")

    # (E) success across PCA space (tertile heat-map)
    if "PC1" in d and "PC2" in d and d["PC1"].notna().sum() > 20:
        d["_p1"] = pd.qcut(d["PC1"], 3, labels=["low", "mid", "high"], duplicates="drop")
        d["_p2"] = pd.qcut(d["PC2"], 3, labels=["low", "mid", "high"], duplicates="drop")
        piv = d.pivot_table(index="_p2", columns="_p1", values="success", aggfunc="mean") * 100
        im = ax[1, 1].imshow(piv.values, cmap="RdYlGn", vmin=0, vmax=100, origin="lower", aspect="auto")
        ax[1, 1].set_xticks(range(piv.shape[1])); ax[1, 1].set_xticklabels(piv.columns)
        ax[1, 1].set_yticks(range(piv.shape[0])); ax[1, 1].set_yticklabels(piv.index)
        for i in range(piv.shape[0]):
            for j in range(piv.shape[1]):
                vv = piv.values[i, j]
                if vv == vv:
                    ax[1, 1].text(j, i, f"{vv:.0f}%", ha="center", va="center", fontsize=9)
        plt.colorbar(im, ax=ax[1, 1], label="Pose ≤2 Å (percent)")
        ax[1, 1].set_xlabel("Principal Component 1 tertile")
        ax[1, 1].set_ylabel("Principal Component 2 tertile")
        ax[1, 1].set_title("Docking success across\nligand chemical space")

    # (F) per-tool success vs flexibility
    for t in TOOLS:
        col = f"{t}_oracle_rmsd"
        if col not in d:
            continue
        rate = (pd.to_numeric(d[col], errors="coerce") <= 2.0).groupby(d["_rb"]).mean() * 100
        ax[1, 2].plot(range(len(rate)), rate.values, "-o", color=TOOL_COLORS[t], label=t)
    ax[1, 2].set_xticks(range(len(labs))); ax[1, 2].set_xticklabels(labs)
    ax[1, 2].set_title("Per-tool success vs ligand flexibility")
    ax[1, 2].set_xlabel("Rotatable bonds"); ax[1, 2].set_ylabel("Pose ≤2 Å (percent)")
    ax[1, 2].set_ylim(0, 105); ax[1, 2].legend(fontsize=8)

    for a in ax.ravel():
        a.grid(alpha=0.2)
    _label_panels(ax)
    fig.suptitle("Which ligands dock well? Pose accuracy vs ligand physicochemistry "
                 f"(n={len(d)} benchmark complexes with crystal). "
                 "Oracle RMSD = closest pose any tool produced.", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "descriptor_vs_quality.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")
    ap.add_argument("--benchmark-dir", default="Data/PoseBuster Benchmark Set")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--p2rank-dir", default="pocket_results/p2rank_results")
    ap.add_argument("--ids-file", default=None)
    ap.add_argument("--equibind-variant", default="equibind_unguided_smina_clampOFF")
    ap.add_argument("--diffdock-variant", default="diffdock",
                    help="Which DiffDock variant to cluster: 'diffdock' (raw) | "
                         "'diffdock_smina' | 'diffdock_gnina'.")
    ap.add_argument("--distance-mode", choices=("centroid", "hybrid", "both"),
                    default="both")
    ap.add_argument("--hybrid-weight", type=float, default=0.5)
    ap.add_argument("--match-thr", type=float, default=4.0)
    ap.add_argument("--top-n-pockets", type=int, default=3)
    ap.add_argument("--pb-valid-only", action="store_true")
    ap.add_argument("--outlier-dist", type=float, default=100.0)
    ap.add_argument("--rmsd-pose-cap", type=int, default=80,
                    help="Skip hybrid RMSD matrix for complexes with more poses.")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="Cap #complexes (smoke test).")
    ap.add_argument("--out-dir", default="posebusters_results/cluster_crystal_pocket")
    ap.add_argument("--features-csv",
                    default="PoseBusters_Benchmark_Analysis/ligand_protein_features.csv",
                    help="Per-ligand RDKit descriptors (from PoseBusters_DataSet_Analysis.ipynb).")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    csv = Path(args.per_pose_csv)
    if not csv.exists():
        print(f"per-pose CSV not found: {csv}")
        return 1
    ids = _load_ids(Path(args.ids_file)) if args.ids_file else None
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    df = load_poses(csv, args.equibind_variant, ids, args.pb_valid_only, args.diffdock_variant)
    if df.empty:
        print("No poses after filtering (check --equibind-variant / --diffdock-variant / --ids-file).")
        return 1
    complexes = sorted(df["protein"].unique())
    if args.limit:
        complexes = complexes[:args.limit]
        df = df[df["protein"].isin(complexes)].copy()
    print(f"Complexes: {len(complexes)} | poses: {len(df)} | "
          f"tools: {sorted(df['tool'].unique())} | equibind={args.equibind_variant} | "
          f"diffdock={args.diffdock_variant}")

    # ── parallel centroid extraction ─────────────────────────────────────
    files = sorted(df["pose_file"].unique())
    cents: Dict[str, Tuple[float, float, float]] = {}
    print(f"Extracting centroids for {len(files)} pose SDFs ({args.workers} workers)...")
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
        for path, cx, cy, cz, na in ex.map(_extract_centroid, files, chunksize=64):
            if cx is not None:
                cents[path] = (cx, cy, cz)
    print(f"  loaded {len(cents)}/{len(files)} centroids")

    # ── per-complex outlier removal (>outlier-dist from median) ──────────
    keep_idx = []
    for cid, sub in df.groupby("protein"):
        cc = np.asarray([cents[f] for f in sub["pose_file"] if f in cents])
        if len(cc) == 0:
            continue
        med = np.median(cc, axis=0)
        for i, f in zip(sub.index, sub["pose_file"]):
            c = cents.get(f)
            if c is not None and np.linalg.norm(np.asarray(c) - med) <= args.outlier_dist:
                keep_idx.append(i)
    df = df.loc[keep_idx].copy()

    do_hybrid = args.distance_mode in ("hybrid", "both")
    results = []
    for ci, cid in enumerate(complexes, 1):
        sub = df[df["protein"] == cid]
        if sub.empty:
            continue
        crystal = crystal_centroid(Path(args.benchmark_dir), cid)
        stem = f"{cid}_protein"
        fp_out = Path(args.fpocket_dir) / f"{stem}_out"
        pr_csv = Path(args.p2rank_dir) / f"{stem}.pdb_predictions.csv"
        fp_pockets = parse_fpocket(fp_out)[:args.top_n_pockets] if fp_out.exists() else []
        pr_pockets = parse_p2rank(pr_csv)[:args.top_n_pockets] if pr_csv.exists() else []
        res = analyze_complex(cid, sub, cents, crystal, fp_pockets, pr_pockets,
                              args.match_thr, do_hybrid, args.hybrid_weight,
                              args.rmsd_pose_cap)
        results.append(res)
        if ci % 50 == 0:
            print(f"  analyzed {ci}/{len(complexes)}")

    ok = [r for r in results if not r.get("skipped")]
    if not ok:
        print("No complex had >=2 usable poses.")
        return 1

    # ── per-complex CSV (drop figure payloads) ───────────────────────────
    df_complex = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")}
                               for r in ok])
    df_complex.to_csv(out_dir / "per_complex_summary.csv", index=False)

    # ── per-cluster CSV ──────────────────────────────────────────────────
    crows = []
    for r in ok:
        cr = r["_crystal"]
        for p in r["_pockets"]:
            crows.append({
                "protein": r["protein"], "cluster_rank": p["rank"], "size": p["size"],
                "n_tools": p["n_tools"], "tools": ",".join(p["tools"]),
                "radius_A": p["radius"],
                "center_x": round(float(p["center"][0]), 3),
                "center_y": round(float(p["center"][1]), 3),
                "center_z": round(float(p["center"][2]), 3),
                "dist_to_crystal": round(_dist(p["center"], cr), 3) if cr is not None else np.nan,
            })
    pd.DataFrame(crows).to_csv(out_dir / "per_cluster.csv", index=False)

    # ── long-form crystal distance by source ─────────────────────────────
    lrows = []
    for r in ok:
        if not r["has_crystal"]:
            continue
        for t in TOOLS:
            lrows.append({"protein": r["protein"], "source": t,
                          "top1_dist": r.get(f"{t}_top1_centroid_dist"),
                          "oracle_dist": r.get(f"{t}_oracle_centroid_dist"),
                          "oracle_hit": r.get(f"{t}_oracle_centroid_hit"),
                          "top1_rmsd": r.get(f"{t}_top1_rmsd"),
                          "oracle_rmsd": r.get(f"{t}_oracle_rmsd")})
        lrows.append({"protein": r["protein"], "source": "ensemble_clusters",
                      "top1_dist": r["ens_top1_dist"], "oracle_dist": r["ens_oracle_dist"],
                      "oracle_hit": r["ens_oracle_hit"], "oracle_rank": r["ens_oracle_rank"]})
        lrows.append({"protein": r["protein"], "source": "fpocket",
                      "top1_dist": r["fpocket_top1_dist"], "oracle_dist": r["fpocket_oracle_dist"],
                      "oracle_hit": r["fpocket_oracle_hit"], "oracle_rank": r["fpocket_oracle_rank"]})
        lrows.append({"protein": r["protein"], "source": "p2rank",
                      "top1_dist": r["p2rank_top1_dist"], "oracle_dist": r["p2rank_oracle_dist"],
                      "oracle_hit": r["p2rank_oracle_hit"], "oracle_rank": r["p2rank_oracle_rank"]})
    df_long = pd.DataFrame(lrows)
    df_long.to_csv(out_dir / "crystal_distance_by_source.csv", index=False)

    # ── rank-in-correct-cluster CSV (per complex × tool) ─────────────────
    rrows = []
    for r in ok:
        if not r["has_crystal"]:
            continue
        for t in TOOLS:
            rrows.append({
                "protein": r["protein"],
                "correct_cluster_dist": r.get("correct_cluster_dist"),
                "correct_cluster_is_hit": r.get("correct_cluster_is_hit"),
                "tool": t,
                "best_rank_in_correct_cluster": r.get(f"{t}_rank_in_correct_cluster"),
                "n_poses_in_correct_cluster": r.get(f"{t}_n_in_correct_cluster"),
                "rank_source": r.get(f"{t}_rank_source"),
            })
    pd.DataFrame(rrows).to_csv(out_dir / "rank_in_correct_cluster.csv", index=False)

    # ── ensemble summary CSV (per tool-combination) ──────────────────────
    dcx = df_complex[df_complex["has_crystal"]] if "has_crystal" in df_complex else df_complex
    def _hit_rate(series):
        if series is None:
            return np.nan
        v = series.map({True: 1.0, False: 0.0}).dropna()
        return round(float(v.mean()), 4) if len(v) else np.nan
    erows = []
    for nm, subset in ENSEMBLES.items():
        od = pd.to_numeric(dcx.get(f"ens[{nm}]_oracle_dist"), errors="coerce")
        erows.append({
            "ensemble": nm, "tools": "+".join(subset), "n_tools": len(subset),
            "oracle_hit_rate": _hit_rate(dcx.get(f"ens[{nm}]_oracle_hit")),
            "top1_hit_rate": _hit_rate(dcx.get(f"ens[{nm}]_top1_hit")),
            "median_oracle_dist": round(float(od.median()), 3) if od is not None and od.notna().any() else np.nan,
        })
    df_ens = pd.DataFrame(erows)
    df_ens.to_csv(out_dir / "ensembles_summary.csv", index=False)

    # ── per-rank (top-5) distance CSV (per complex × tool × rank) ─────────
    prrows = []
    for r in ok:
        if not r["has_crystal"]:
            continue
        for (cid_, t, rp, cd, rm) in r.get("_per_rank", []):
            prrows.append({"protein": cid_, "tool": t, "rank": rp,
                           "centroid_dist": None if cd != cd else round(cd, 3),
                           "rmsd": None if rm != rm else round(rm, 3)})
    df_rank = pd.DataFrame(prrows)
    df_rank.to_csv(out_dir / "per_rank_distance.csv", index=False)

    # ── aggregate report ─────────────────────────────────────────────────
    dc = df_complex[df_complex["has_crystal"]] if "has_crystal" in df_complex else df_complex
    n = len(dc)
    def med_hit(col):
        v = pd.to_numeric(dc[col], errors="coerce").dropna()
        return (float(v.median()), float((v <= args.match_thr).mean())) if len(v) else (np.nan, np.nan)

    print("\n" + "=" * 78)
    print(f"DOCKED CLUSTERS vs fpocket/p2rank vs CRYSTAL  "
          f"(n={n} complexes w/ crystal, thr={args.match_thr:g} Å)")
    print("=" * 78)
    print(f"  {'source':<20}{'median oracle Å':>16}{'oracle hit':>12}{'top1 hit':>12}")
    src_cols = {
        "autodock": ("autodock_oracle_centroid_dist", "autodock_top1_centroid_dist"),
        "diffdock": ("diffdock_oracle_centroid_dist", "diffdock_top1_centroid_dist"),
        "equibind": ("equibind_oracle_centroid_dist", "equibind_top1_centroid_dist"),
        "ensemble_clusters": ("ens_oracle_dist", "ens_top1_dist"),
        "fpocket": ("fpocket_oracle_dist", "fpocket_top1_dist"),
        "p2rank": ("p2rank_oracle_dist", "p2rank_top1_dist"),
    }
    for s, (oc, tc) in src_cols.items():
        if oc not in dc:
            continue
        med, ohit = med_hit(oc)
        _, thit = med_hit(tc)
        print(f"  {s:<20}{med:>16.2f}{ohit:>11.0%}{thit:>12.0%}")

    # pose-quality (RMSD) for the 3 tools
    print("\n  Pose quality vs crystal (RMSD, oracle over poses):")
    for t in TOOLS:
        oc = f"{t}_oracle_rmsd"
        if oc in dc:
            v = pd.to_numeric(dc[oc], errors="coerce").dropna()
            if len(v):
                print(f"    {t:<10} median {v.median():.2f} Å | "
                      f"≤2 Å: {(v <= 2).mean():.0%} | ≤5 Å: {(v <= 5).mean():.0%}")

    # triangulation
    tri = {k: float(dc[f"tri_{k}"].mean()) for k in ("docking_hit", "fpocket_hit", "p2rank_hit")
           if f"tri_{k}" in dc}
    print("\n  Triangulation — finds the true site (oracle ≤ thr):")
    for k, v in tri.items():
        print(f"    {k:<14} {v:.0%}")
    if {"tri_docking_hit", "tri_fpocket_hit", "tri_p2rank_hit"} <= set(dc.columns):
        all3 = float(((dc["tri_docking_hit"]) & (dc["tri_fpocket_hit"]) & (dc["tri_p2rank_hit"])).mean())
        none = float((~(dc["tri_docking_hit"]) & ~(dc["tri_fpocket_hit"]) & ~(dc["tri_p2rank_hit"])).mean())
        print(f"    all three: {all3:.0%} | none: {none:.0%}")

    # cross-tool agreement + cluster-vs-pocket
    if "n_agree_tool_pairs" in dc:
        print(f"\n  Cross-tool agreement: mean agreeing tool-pairs/complex "
              f"{dc['n_agree_tool_pairs'].mean():.2f}/3 | all-3-agree "
              f"{float(dc['all_three_agree'].mean()):.0%}")
    for col, lab in [("cluster_fpocket_match_frac", "fpocket"),
                     ("cluster_p2rank_match_frac", "p2rank")]:
        if col in dc:
            v = pd.to_numeric(dc[col], errors="coerce").dropna()
            if len(v):
                print(f"  Docked clusters matching {lab} pockets (≤thr): mean {v.mean():.2f}")

    # ── how each tool ranks poses in the crystal-closest cluster ─────────
    hit = dc[dc.get("correct_cluster_is_hit") == True] if "correct_cluster_is_hit" in dc else dc.iloc[0:0]  # noqa: E712
    print(f"\n  TOOL RANKING within the crystal-closest cluster "
          f"(true-site complexes only, n={len(hit)}):")
    print(f"  {'tool':<10}{'coverage':>10}{'median rank':>13}{'%top-1':>9}{'%top-5':>9}  rank source")
    for t in TOOLS:
        col = f"{t}_rank_in_correct_cluster"
        if col not in hit:
            continue
        r = pd.to_numeric(hit[col], errors="coerce")
        cov = float(r.notna().mean()) if len(hit) else np.nan
        rr = r.dropna()
        src_col = f"{t}_rank_source"
        srcs = hit[src_col].dropna().unique().tolist() if src_col in hit else []
        srctxt = "/".join(map(str, srcs)) if srcs else "?"
        if len(rr):
            print(f"  {t:<10}{cov:>9.0%}{rr.median():>13.0f}{(rr == 1).mean():>9.0%}"
                  f"{(rr <= 5).mean():>9.0%}  {srctxt}")
        else:
            print(f"  {t:<10}{cov:>9.0%}{'n/a':>13}{'n/a':>9}{'n/a':>9}  {srctxt}")

    # ── tool-combination ensembles ──────────────────────────────────────
    print(f"\n  ENSEMBLE EXPLORATION (which tools combined; thr={args.match_thr:g} Å):")
    print(f"  {'ensemble':<10}{'tools':<18}{'oracle hit':>12}{'top1 hit':>10}{'median Å':>11}")
    for _, row in df_ens.iterrows():
        oh = f"{row['oracle_hit_rate']:.0%}" if row['oracle_hit_rate'] == row['oracle_hit_rate'] else "n/a"
        th = f"{row['top1_hit_rate']:.0%}" if row['top1_hit_rate'] == row['top1_hit_rate'] else "n/a"
        md = f"{row['median_oracle_dist']:.2f}" if row['median_oracle_dist'] == row['median_oracle_dist'] else "n/a"
        print(f"  {row['ensemble']:<10}{row['tools']:<18}{oh:>12}{th:>10}{md:>11}")

    # ── per-rank (top-5) distance-from-crystal profile ───────────────────
    if not df_rank.empty:
        print("\n  TOP-5 RANKED POSES — MEDIAN distance from crystal (centroid Å | RMSD Å):")
        print("  (median is robust; means are tail-inflated by far poses — see figure/CSV)")
        print(f"  {'tool':<10}" + "".join(f"{'rank'+str(k):>16}" for k in range(1, 6)))
        for t in TOOLS:
            g = df_rank[df_rank.tool == t]
            cells = []
            for k in range(1, 6):
                gk = g[g["rank"] == k]
                cd = pd.to_numeric(gk["centroid_dist"], errors="coerce").median()
                rm = pd.to_numeric(gk["rmsd"], errors="coerce").median()
                cells.append(f"{cd:>6.1f}|{rm:>6.1f}" if cd == cd else f"{'n/a':>13}")
            print(f"  {t:<10}" + "".join(f"{c:>16}" for c in cells))

    # ── top-5 cluster consistency ────────────────────────────────────────
    print("\n  TOP-5 CLUSTER CONSISTENCY (are a tool's 5 best poses in one site?):")
    print(f"  {'tool':<10}{'mean modal frac':>16}{'mean #clusters':>16}"
          f"{'closest-is-rank1':>18}")
    for t in TOOLS:
        mf = pd.to_numeric(dc.get(f"{t}_top5_modal_frac"), errors="coerce").dropna() \
            if f"{t}_top5_modal_frac" in dc else pd.Series(dtype=float)
        nc = pd.to_numeric(dc.get(f"{t}_top5_n_clusters"), errors="coerce").dropna() \
            if f"{t}_top5_n_clusters" in dc else pd.Series(dtype=float)
        br = pd.to_numeric(dc.get(f"{t}_top5_best_rank"), errors="coerce").dropna() \
            if f"{t}_top5_best_rank" in dc else pd.Series(dtype=float)
        if len(mf):
            print(f"  {t:<10}{mf.mean():>16.2f}{nc.mean():>16.2f}"
                  f"{(br == 1).mean():>17.0%}")

    summary = {
        "n_complexes": int(n), "match_thr_A": args.match_thr,
        "equibind_variant": args.equibind_variant,
        "diffdock_variant": args.diffdock_variant,
        "median_oracle_dist": {s: (float(pd.to_numeric(dc[oc], errors="coerce").median())
                                    if oc in dc else None)
                               for s, (oc, _) in src_cols.items()},
        "oracle_hit_rate": {s: (float((pd.to_numeric(dc[oc], errors="coerce") <= args.match_thr).mean())
                                 if oc in dc else None)
                            for s, (oc, _) in src_cols.items()},
        "triangulation_hit_rate": tri,
        "ensembles": df_ens.to_dict(orient="records"),
        "rank_in_correct_cluster": {
            t: {
                "coverage": (float(pd.to_numeric(hit[f"{t}_rank_in_correct_cluster"],
                                                  errors="coerce").notna().mean())
                             if f"{t}_rank_in_correct_cluster" in hit and len(hit) else None),
                "median_best_rank": (float(pd.to_numeric(hit[f"{t}_rank_in_correct_cluster"],
                                                          errors="coerce").median())
                                     if f"{t}_rank_in_correct_cluster" in hit else None),
                "pct_top1": (float((pd.to_numeric(hit[f"{t}_rank_in_correct_cluster"],
                                                  errors="coerce") == 1).mean())
                             if f"{t}_rank_in_correct_cluster" in hit and len(hit) else None),
                "rank_source": (hit[f"{t}_rank_source"].dropna().unique().tolist()
                                if f"{t}_rank_source" in hit else []),
            } for t in TOOLS
        },
        "n_true_site_complexes": int(len(hit)),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(f"\n  CSVs + summary written to: {out_dir}/")

    if not args.no_plot:
        figs = [_fig_examples(ok, args.match_thr, out_dir),
                _fig_summary(df_complex, args.match_thr, out_dir),
                _fig_ecdf_split(df_complex, args.match_thr, out_dir),
                _fig_top5(df_rank, df_complex, out_dir),
                _fig_rank_in_correct(df_complex, args.match_thr, out_dir),
                _fig_ensembles(df_complex, args.match_thr, out_dir),
                _fig_descriptor_quality(df_complex, args.features_csv, out_dir)]
        for f in figs:
            if f:
                print(f"  Figure: {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
