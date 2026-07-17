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

Each of the three tool slots is filled by ONE docking-mode variant, chosen from
the CSV:
  * EquiBind  : pick the pocket source, refiner, and clamp independently
                (``--equibind-pocket``/``--equibind-refine``/``--equibind-clamp``),
                or pass a full method string with ``--equibind-variant``.
  * DiffDock  : pick the refiner (``--diffdock-refine {raw,smina,gnina}``) or a
                full method string with ``--diffdock-variant``.
The refiner axis is what lets you compare e.g. gnina- vs smina-refined poses.
Note the gnina EquiBind method names OMIT the ``gnina`` token (e.g.
``equibind_fpocket_clampON`` IS the gnina refine of fpocket/clampON); the
component selectors resolve that quirk from the CSV so you don't have to.

Run in the analysis env (conda env ``vina`` — rdkit + sklearn + scipy + matplotlib;
do NOT rely on sklearn_extra, it is broken under NumPy 2.x):

    # default: EquiBind unguided/smina/clampOFF + raw DiffDock
    python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
        --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"

    # gnina-refined EquiBind (fpocket-guided) vs gnina-refined DiffDock
    python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
        --equibind-pocket fpocket --equibind-refine gnina \
        --diffdock-refine gnina
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
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


def _equibind_rank_note(eq_variant: Optional[str]) -> str:
    """Figure-title note stating which docking score actually ranks EquiBind in
    THIS run (EquiBind has no native rank). The refiner token in the variant name
    determines the populated affinity column: gnina→gnina_affinity, smina→
    smina_affinity, raw→pose-generation order (no score). See _effective_ranks."""
    v = str(eq_variant).lower()
    if "gnina" in v:
        return "EquiBind rank = gnina-affinity proxy"
    if "smina" in v:
        return "EquiBind rank = smina-affinity proxy"
    return "EquiBind rank = pose-generation order"


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


def _safe_silhouette(dm: np.ndarray, labels: np.ndarray) -> float:
    """Silhouette only where it is defined (2 <= k <= n-1); NaN otherwise (e.g. k=1)."""
    from sklearn.metrics import silhouette_score
    k = len(set(labels.tolist()))
    if k < 2 or k >= len(labels):
        return float("nan")
    try:
        return float(silhouette_score(dm, labels, metric="precomputed"))
    except Exception:
        return float("nan")


def _gap_statistic(dm: np.ndarray, cents: np.ndarray, max_k: int = 10,
                   n_ref: int = 10, seed: int = 0) -> int:
    """Tibshirani gap statistic; returns the chosen k, WHICH MAY BE 1.

    Within-cluster dispersion W_k uses the sum of intra-cluster pairwise
    distances / (2 n_c); the null is a uniform cloud in the centroid bounding
    box. Picks the smallest k with Gap(k) >= Gap(k+1) - s_{k+1}.
    """
    n = cents.shape[0]

    def _Wk(d: np.ndarray, lab: np.ndarray) -> float:
        w = 0.0
        for l in set(lab.tolist()):
            idx = np.where(lab == l)[0]
            if len(idx) > 1:
                sub = d[np.ix_(idx, idx)]
                w += float(sub.sum()) / (2.0 * len(idx))
        return w

    ks = list(range(1, min(max_k, n - 1) + 1))
    if len(ks) <= 1:
        return 1
    logW = []
    for k in ks:
        lab = np.zeros(n, dtype=int) if k == 1 else SimpleKMedoids(n_clusters=k).fit(dm).labels_
        logW.append(np.log(_Wk(dm, lab) + 1e-12))
    lo, hi = cents.min(axis=0), cents.max(axis=0)
    rng = np.random.RandomState(seed)
    ref = np.zeros((n_ref, len(ks)))
    for b in range(n_ref):
        rc = rng.uniform(lo, hi, size=cents.shape)
        rdm = _centroid_dm(rc)
        for ik, k in enumerate(ks):
            lab = np.zeros(n, dtype=int) if k == 1 else SimpleKMedoids(n_clusters=k).fit(rdm).labels_
            ref[b, ik] = np.log(_Wk(rdm, lab) + 1e-12)
    gap = ref.mean(axis=0) - np.asarray(logW)
    sk = ref.std(axis=0) * np.sqrt(1.0 + 1.0 / n_ref)
    for ik in range(len(ks) - 1):
        if gap[ik] >= gap[ik + 1] - sk[ik + 1]:
            return ks[ik]
    return ks[-1]


def cluster_sites(dm: np.ndarray, cents: np.ndarray, method: str = "threshold",
                  pocket_radius: float = 8.0, max_k: int = 10):
    """Partition poses into spatial sites, ALLOWING k=1 (unlike _select_k).

    Returns (labels, k, silhouette-or-NaN).
      * 'threshold'  : complete-linkage cut at ``pocket_radius`` — k is discovered
                       (k>=1), every cluster's max internal spread <= pocket_radius,
                       deterministic, encodes the physical pocket scale.
      * 'gap'        : gap statistic picks k (may be 1), then KMedoids at that k.
      * 'silhouette' : legacy silhouette-max KMedoids (never returns k=1).
    """
    n = dm.shape[0]
    if n < 2:
        return np.zeros(n, dtype=int), 1, float("nan")
    if method == "silhouette":
        return _select_k(dm, "kmedoids", max_k)
    if method == "gap":
        k = _gap_statistic(dm, cents, max_k)
        if k <= 1:
            return np.zeros(n, dtype=int), 1, float("nan")
        labels = SimpleKMedoids(n_clusters=k).fit(dm).labels_
        _, labels = np.unique(labels, return_inverse=True)
        return labels, len(set(labels.tolist())), _safe_silhouette(dm, labels)
    # default: distance-threshold complete-linkage cut (k discovered, k>=1)
    from sklearn.cluster import AgglomerativeClustering
    labels = AgglomerativeClustering(
        n_clusters=None, distance_threshold=pocket_radius,
        metric="precomputed", linkage="complete").fit_predict(dm)
    _, labels = np.unique(labels, return_inverse=True)      # contiguous 0..k-1
    return labels, len(set(labels.tolist())), _safe_silhouette(dm, labels)


# ════════════════════════════════════════════════════════════════════════
# Placement-aware (in-place RMSD) binding-MODE clustering
# ════════════════════════════════════════════════════════════════════════
# The centroid axis answers WHERE a pose sits; this answers HOW it sits. In-place
# symmetry-aware RMSD (rdMolAlign.CalcRMS — NO re-superposition) folds translation +
# orientation + conformation into one Å distance, so two poses in different pockets,
# or flipped 180° in the same pocket, are correctly far apart (GetBestRMS would call
# them identical). All poses are the same ligand in the same receptor frame, so the
# raw coordinates are directly comparable and CalcRMS always applies.

def _inplace_rmsd(a: Optional[Chem.Mol], b: Optional[Chem.Mol]) -> float:
    if a is None or b is None:
        return float("nan")
    try:
        return float(rdMolAlign.CalcRMS(a, b))          # in-place, symmetry-aware
    except Exception:
        try:                                            # atom-count mismatch fallback
            return float(np.linalg.norm(centroid_from_mol(a) - centroid_from_mol(b)))
        except Exception:
            return float("nan")


def _inplace_rmsd_dm(mols: List[Optional[Chem.Mol]]) -> np.ndarray:
    n = len(mols)
    dm = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dm[i, j] = dm[j, i] = _inplace_rmsd(mols[i], mols[j])
    return dm


def cluster_modes(rmsd_dm: np.ndarray, thr: float = 2.0):
    """Placement-aware binding-mode partition: complete-linkage cut of the in-place
    RMSD matrix at ``thr`` Å (the near-native cutoff), k>=1 (a set of poses that are
    all one mode stays one cluster). Returns (labels, k, silhouette-or-NaN)."""
    n = rmsd_dm.shape[0]
    if n < 2:
        return np.zeros(n, dtype=int), 1, float("nan")
    dm = rmsd_dm.copy()
    if not np.isfinite(dm).all():                       # substitute failed pairs
        finite = dm[np.isfinite(dm)]
        dm[~np.isfinite(dm)] = float(finite.max()) if finite.size else 0.0
    from sklearn.cluster import AgglomerativeClustering
    labels = AgglomerativeClustering(
        n_clusters=None, distance_threshold=thr,
        metric="precomputed", linkage="complete").fit_predict(dm)
    _, labels = np.unique(labels, return_inverse=True)
    return labels, len(set(labels.tolist())), _safe_silhouette(dm, labels)


def _mode_quality(rmsd_dm: np.ndarray, labels: np.ndarray) -> dict:
    """Compactness (mean within-mode median RMSD) and separation (nearest inter-mode
    medoid RMSD) on the in-place RMSD matrix."""
    labs = sorted(set(labels.tolist()))
    spreads, medoids = [], []
    for l in labs:
        idx = np.where(labels == l)[0]
        if len(idx) > 1:
            sub = rmsd_dm[np.ix_(idx, idx)]
            spreads.append(float(np.nanmedian(sub[np.triu_indices(len(idx), 1)])))
            medoids.append(int(idx[np.argmin(np.nansum(sub, axis=1))]))
        else:
            spreads.append(0.0)
            medoids.append(int(idx[0]))
    out = dict(compactness=round(float(np.mean(spreads)), 3), separation=np.nan)
    if len(medoids) > 1:
        md = rmsd_dm[np.ix_(medoids, medoids)]
        vals = md[np.triu_indices(len(medoids), 1)]
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out["separation"] = round(float(vals.min()), 3)
    return out


# Ranking rules for candidate pockets. 'size' is the legacy throughput-biased key
# (kept for the ablation and for backward-compatible callers, e.g. the Orai script);
# 'consensus' ranks by distinct-tool agreement first, so the highest-throughput
# tool can no longer win a site by ballot-stuffing.
_POCKET_SORT_KEYS = {
    "size":       lambda p: (-p["size"],),
    "consensus":  lambda p: (-p["n_tools"], -p["conf_weight"], p["spread"], -p["size"]),
    "tight":      lambda p: (0 if p["size"] >= 2 else 1, p["spread"], -p["size"]),
    "confidence": lambda p: (-p["conf_weight"], -p["size"]),
}


def _sort_pockets(pockets: List[dict], rank_by: str) -> None:
    pockets.sort(key=_POCKET_SORT_KEYS.get(rank_by, _POCKET_SORT_KEYS["size"]))


def pockets_from_labels(cents: np.ndarray, labels: np.ndarray, tools: List[str],
                        files: Optional[List[str]] = None,
                        pose_weights: Optional[Dict[str, float]] = None,
                        rank_by: str = "size") -> List[dict]:
    """Cluster -> candidate-pocket records.

    Each pocket carries both a mean ``center`` and a robust ``medoid_center``
    (outlier-resistant), a ``radius`` (max) and a robust ``spread`` (median), the
    distinct-tool count, and a ``conf_weight`` (sum of member pose confidences,
    from ``pose_weights``; falls back to pose count). ``rank_by`` selects the
    ordering (default 'size' keeps legacy/Orai behaviour; the main pipeline
    passes 'consensus').
    """
    out = []
    for lab in sorted(set(labels.tolist())):
        idx = np.where(labels == lab)[0]
        c = cents[idx]
        center = c.mean(axis=0)
        d_to_center = np.linalg.norm(c - center, axis=1)
        radius = float(d_to_center.max()) if len(c) > 1 else 0.0
        spread = float(np.median(d_to_center)) if len(c) > 1 else 0.0
        if len(c) > 1:
            sub = _centroid_dm(c)
            medoid_center = c[int(np.argmin(sub.sum(axis=1)))]
        else:
            medoid_center = center
        tl = [tools[i] for i in idx]
        if files is not None and pose_weights is not None:
            conf = float(sum(pose_weights.get(files[i], 0.0) for i in idx))
        else:
            conf = float(len(idx))
        out.append({
            "center": center, "medoid_center": medoid_center,
            "size": int(len(idx)), "radius": round(radius, 2),
            "spread": round(spread, 2), "conf_weight": round(conf, 3),
            "tools": sorted(set(tl)), "n_tools": len(set(tl)),
            "tool_counts": dict(Counter(tl)),
            "label": int(lab), "members": idx,
        })
    _sort_pockets(out, rank_by)
    for i, p in enumerate(out, 1):
        p["rank"] = i
    return out


def _pose_weights(files: List[str], eff_rank: Dict[str, int]) -> Dict[str, float]:
    """Per-pose confidence weight = 1/effective-rank (rank-1 pose weighted most)."""
    w: Dict[str, float] = {}
    for f in files:
        r = eff_rank.get(f)
        w[f] = 1.0 / float(r) if (r is not None and r > 0) else 0.25
    return w


def _internal_indices(C: np.ndarray, labels: np.ndarray) -> dict:
    """Internal cluster-validity indices on the 3D centroid cloud.

    Calinski-Harabasz / Davies-Bouldin need k>=2 (NaN at k=1); compactness (mean
    within-cluster median spread, lower=tighter) and separation (nearest inter-
    cluster-center distance, higher=better) are defined for any k.
    """
    out = dict(ch_score=np.nan, db_score=np.nan, compactness=np.nan, separation=np.nan)
    labs = sorted(set(labels.tolist()))
    k = len(labs)
    centers, spreads = [], []
    for lab in labs:
        cc = C[np.where(labels == lab)[0]]
        ctr = cc.mean(axis=0)
        centers.append(ctr)
        spreads.append(float(np.median(np.linalg.norm(cc - ctr, axis=1))) if len(cc) > 1 else 0.0)
    out["compactness"] = round(float(np.mean(spreads)), 3)
    if k > 1:
        cen = np.asarray(centers)
        dc = _centroid_dm(cen)
        iu = np.triu_indices(k, 1)
        out["separation"] = round(float(dc[iu].min()), 3)
    if 2 <= k < len(labels):
        from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
        try:
            out["ch_score"] = round(float(calinski_harabasz_score(C, labels)), 3)
        except Exception:
            pass
        try:
            out["db_score"] = round(float(davies_bouldin_score(C, labels)), 3)
        except Exception:
            pass
    return out


def _bootstrap_stability(C: np.ndarray, base_labels: np.ndarray, method: str,
                         pocket_radius: float, n_boot: int = 25,
                         seed: int = 0) -> Tuple[float, int]:
    """Clusterboot-style stability: mean over base clusters of their best Jaccard
    to a bootstrap-resampled clustering (>0.75 stable, <0.5 dissolved)."""
    n = len(C)
    if n < 4 or n_boot <= 0:
        return float("nan"), int(len(set(base_labels.tolist())))
    base = [set(np.where(base_labels == l)[0].tolist()) for l in sorted(set(base_labels.tolist()))]
    rng = np.random.RandomState(seed)
    jacc = np.zeros(len(base))
    used = 0
    for _ in range(n_boot):
        uniq = np.unique(rng.choice(n, n, replace=True))
        if len(uniq) < 2:
            continue
        used += 1
        lab = cluster_sites(_centroid_dm(C[uniq]), C[uniq], method, pocket_radius)[0]
        boot = [set(uniq[np.where(lab == l)[0]].tolist()) for l in set(lab.tolist())]
        present = set(uniq.tolist())
        for ci, oc in enumerate(base):
            oc_p = oc & present
            best = 0.0
            for bc in boot:
                u = len(oc_p | bc)
                if u:
                    best = max(best, len(oc_p & bc) / u)
            jacc[ci] += best
    if used == 0:
        return float("nan"), len(base)
    return round(float((jacc / used).mean()), 3), len(base)


def _cluster_purity(cp: dict, C: np.ndarray, crystal: np.ndarray,
                    files: List[str], rmsd_map: Dict[str, float],
                    thr: float) -> Tuple[float, float, float]:
    """For the crystal-closest cluster: fraction of members within thr (centroid)
    and within 2 A RMSD of the crystal, plus the best (min) member RMSD. A tight
    correct cluster scores high; a diffuse one that only averages near the site
    scores low."""
    idx = cp["members"]
    d = np.linalg.norm(C[idx] - crystal, axis=1)
    purity_centroid = round(float(np.mean(d <= thr)), 3)
    rmsds = [rmsd_map.get(files[i]) for i in idx]
    rmsds = [float(x) for x in rmsds if x is not None and np.isfinite(x)]
    best_rmsd = round(float(min(rmsds)), 3) if rmsds else np.nan
    purity_rmsd = round(float(np.mean([x <= 2.0 for x in rmsds])), 3) if rmsds else np.nan
    return purity_centroid, purity_rmsd, best_rmsd


def _precision_at_1(pockets: List[dict], crystal: np.ndarray, thr: float) -> dict:
    """Rank-1 pocket under each ranking rule vs crystal — the ablation that shows
    whether a better ranking recovers the true site more often than raw size."""
    out: dict = {}

    def pick(key, eligible=None):
        pool = [p for p in pockets if eligible(p)] if eligible else list(pockets)
        pool = pool or list(pockets)
        return sorted(pool, key=key)[0]

    picks = {
        "size":       pick(_POCKET_SORT_KEYS["size"]),
        "ntools":     pick(_POCKET_SORT_KEYS["consensus"]),
        "tight":      pick(_POCKET_SORT_KEYS["tight"], eligible=lambda p: p["size"] >= 2),
        "confidence": pick(_POCKET_SORT_KEYS["confidence"]),
    }
    for name, p in picks.items():
        dd = _dist(p["center"], crystal)
        out[f"p1_{name}_dist"] = round(dd, 3)
        out[f"p1_{name}_hit"] = bool(dd <= thr)
    # medoid-centered: the consensus pick but measured from its robust medoid center
    pm = picks["ntools"]
    dmd = _dist(pm.get("medoid_center", pm["center"]), crystal)
    out["p1_medoid_dist"] = round(dmd, 3)
    out["p1_medoid_hit"] = bool(dmd <= thr)
    return out


# Tool-combination ensembles to explore (singletons = baselines; combos = ensembles).
ENSEMBLES = {
    "AD": ("autodock",), "DD": ("diffdock",), "EB": ("equibind",),
    "AD+DD": ("autodock", "diffdock"), "AD+EB": ("autodock", "equibind"),
    "DD+EB": ("diffdock", "equibind"),
    "AD+DD+EB": ("autodock", "diffdock", "equibind"),
}


def _ensemble_stats(C: np.ndarray, tools: List[str], crystal: Optional[np.ndarray],
                    subset: Tuple[str, ...], thr: float,
                    site_method: str = "threshold", pocket_radius: float = 8.0) -> dict:
    """For one tool-combination: pose-level oracle + cluster top1 vs crystal.

    The subset is clustered with the k=1-capable site method and 'top1' is the
    consensus-ranked pocket (distinct-tool agreement first), not the largest
    cluster — so a high-throughput tool no longer wins the pick by pose count.
    """
    idx = [i for i, t in enumerate(tools) if t in subset]
    if not idx or crystal is None:
        return dict(n_poses=len(idx), oracle_dist=np.nan, oracle_hit=np.nan,
                    top1_dist=np.nan, top1_hit=np.nan, n_clusters=np.nan)
    sc = C[idx]
    st = [tools[i] for i in idx]
    pose_d = np.linalg.norm(sc - crystal, axis=1)
    oracle = float(pose_d.min())
    if len(sc) < 2:
        top1_center = sc.mean(axis=0)
        n_clusters = 1
    else:
        lab, _, _ = cluster_sites(_centroid_dm(sc), sc, site_method, pocket_radius)
        pk = pockets_from_labels(sc, lab, st, rank_by="consensus")
        top1_center = pk[0]["center"]
        n_clusters = len(pk)
    top1 = _dist(top1_center, crystal)
    return dict(n_poses=len(idx), oracle_dist=round(oracle, 3),
                oracle_hit=bool(oracle <= thr), top1_dist=round(top1, 3),
                top1_hit=bool(top1 <= thr), n_clusters=int(n_clusters))


# ════════════════════════════════════════════════════════════════════════
# IO
# ════════════════════════════════════════════════════════════════════════

# Friendly refiner keywords accepted for both EquiBind and DiffDock. 'raw' = the
# unrefined docked pose; 'smina'/'gnina' = the smina-/gnina-minimised pose.
_REFINERS = ("raw", "smina", "gnina")
_DECOMP_COLS = ("pocket_source", "refine_variant", "clamp_variant")


def _variant_catalog(csv: Path) -> pd.DataFrame:
    """Light read of just the variant-defining columns (``method`` + the
    pocket/refine/clamp decomposition), de-duplicated. Lets the CLI resolve and
    validate a variant selection without parsing the whole ~190k-row per-pose CSV
    at full width first."""
    want = ("method",) + _DECOMP_COLS
    cat = pd.read_csv(csv, usecols=lambda c: c in want, low_memory=False)
    return cat.drop_duplicates().reset_index(drop=True)


def resolve_equibind_variant(cat: pd.DataFrame, variant: Optional[str],
                             pocket: str, refine: str, clamp: str) -> str:
    """Resolve the exact ``equibind_*`` method string for the EquiBind tool slot.

    An explicit full ``--equibind-variant`` method string wins (validated against
    the CSV); otherwise the (pocket, refine, clamp) components are matched. This
    is what makes ``--equibind-refine gnina`` selectable even though the gnina
    method names omit the ``gnina`` token (e.g. ``equibind_fpocket_clampON`` IS
    the gnina refine of fpocket/clampON) — the decomposition columns carry the
    true refiner, so we match on those rather than on the folded name.
    """
    eq = cat[cat["method"].astype(str).str.startswith("equibind")]
    methods = sorted(eq["method"].astype(str).unique())
    if not methods:
        raise SystemExit("No 'equibind' poses in the per-pose CSV.")
    if variant:
        if variant in methods:
            return variant
        raise SystemExit(
            f"--equibind-variant '{variant}' is not in the CSV. Available:\n    "
            + "\n    ".join(methods))
    if set(_DECOMP_COLS) <= set(cat.columns):
        mask = ((eq["pocket_source"].astype(str) == pocket)
                & (eq["refine_variant"].astype(str) == refine)
                & (eq["clamp_variant"].astype(str) == clamp))
        cand = sorted(eq[mask]["method"].astype(str).unique())
    else:                                    # older CSV: compose by naming rule
        mid = "" if refine == "gnina" else f"_{refine}"   # gnina omits the token
        cand = [m for m in methods if m == f"equibind_{pocket}{mid}_{clamp}"]
    sel = f"pocket={pocket}, refine={refine}, clamp={clamp}"
    if len(cand) == 1:
        return cand[0]
    if not cand:
        avail = ""
        if set(_DECOMP_COLS) <= set(cat.columns):
            avail = ("\n  pockets={} refiners={} clamps={}".format(
                sorted(eq["pocket_source"].dropna().unique()),
                sorted(eq["refine_variant"].dropna().unique()),
                sorted(eq["clamp_variant"].dropna().unique())))
        raise SystemExit(f"No equibind variant matches ({sel})." + avail)
    raise SystemExit(f"Ambiguous equibind selection ({sel}) -> {cand}. Narrow it down.")


def resolve_diffdock_variant(cat: pd.DataFrame, variant: Optional[str],
                             refine: str) -> str:
    """Resolve the exact DiffDock method string for the DiffDock tool slot.

    DiffDock encodes its refiner in the method NAME (unlike EquiBind's
    decomposition columns): raw -> ``diffdock``, smina -> ``diffdock_smina``,
    gnina -> ``diffdock_gnina``. ``--diffdock-variant`` accepts either a full
    method string or a bare refiner keyword; ``--diffdock-refine`` is the tidy
    way to say the same thing.
    """
    methods = sorted(m for m in cat["method"].astype(str).unique()
                     if m == "diffdock" or m.startswith("diffdock_"))
    if not methods:
        raise SystemExit("No 'diffdock' poses in the per-pose CSV.")
    if variant:
        if variant in methods:
            return variant
        if variant in _REFINERS:             # shorthand: --diffdock-variant gnina
            refine = variant
        else:
            raise SystemExit(
                f"--diffdock-variant '{variant}' not recognised (want a full method "
                f"name or one of {_REFINERS}). Available: {methods}")
    cand = "diffdock" if refine == "raw" else f"diffdock_{refine}"
    if cand not in methods:
        raise SystemExit(
            f"No DiffDock '{refine}' variant ('{cand}') in the CSV. Available: {methods}")
    return cand


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
    does not score poses — so we fall back to a docking-score ordering (ascending,
    most-negative = rank 1), matching the EquiBind ranking convention in
    ``posebusters_pose_comparison.py``.

    The gnina- and smina-refined EquiBind variants carry MUTUALLY EXCLUSIVE score
    columns (gnina poses populate only ``gnina_affinity``; smina poses only
    ``smina_affinity``), so both are tried — ``gnina_affinity`` FIRST, so the gnina
    variant is ranked by its gnina affinity instead of silently dropping through to
    pose-generation order (which is what happened when only smina was checked).
    Pose order is the last resort. The source is recorded so the report can flag
    proxy rankings.
    """
    eff: Dict[str, int] = {}
    src: Dict[str, str] = {}
    for tool, g in sub.groupby("tool"):
        r = pd.to_numeric(g.get("rank"), errors="coerce")
        if r.notna().any() and r.nunique(dropna=True) > 1:
            order, how = r, "native"
        else:
            order, how = None, None
            for col in ("gnina_affinity", "smina_affinity"):  # gnina first
                if col in g:
                    a = pd.to_numeric(g[col], errors="coerce")
                    if a.nunique(dropna=True) > 1:
                        order, how = a, col
                        break
            if order is None:
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


# Depth of the per-rank profile written to per_rank_distance.csv and used by the
# rank-by-rank figures (e.g. rank_in_crystal_cluster_by_rank.png). Independent of
# the top-5 concentration summary below, which stays fixed at 5.
_RANK_PROFILE_MAX = 15


def _top5_analysis(sub: pd.DataFrame, eff_rank: Dict[str, int],
                   file_idx: Dict[str, int], labels: np.ndarray, C: np.ndarray,
                   crystal: Optional[np.ndarray], correct_label: Optional[int] = None,
                   max_rank: int = 5, per_rank_max: int = _RANK_PROFILE_MAX):
    """Per-tool profile of the top ranked poses.

    Returns (per_rank_rows, tool_summary):
      * per_rank_rows: (tool, rank_pos, centroid_dist, rmsd, in_correct_cluster) for
        the k-th best pose of each tool, up to ``per_rank_max`` ranks — the "how far
        is the rank-k pose from crystal" profile, plus whether that pose landed in
        the crystal-closest cluster (``correct_label``; NaN when no crystal / cluster
        is defined).
      * tool_summary[tool]: how concentrated the top-``max_rank`` (5) poses are (n
        distinct clusters, modal-cluster fraction) and which of them is closest to
        the crystal (top5_best_rank, top5_best_dist). This summary is INDEPENDENT of
        ``per_rank_max`` — it always covers only the first ``max_rank`` positions.
    """
    file_rmsd = (dict(zip(sub["pose_file"],
                          pd.to_numeric(sub.get("rmsd"), errors="coerce")))
                 if "rmsd" in sub else {})
    per_rank: List[tuple] = []
    tool_summary: Dict[str, dict] = {}
    for tool, g in sub.groupby("tool"):
        ranked = sorted(((eff_rank[f], f) for f in g["pose_file"]
                         if eff_rank.get(f) is not None), key=lambda x: x[0])
        top = ranked[:per_rank_max]
        rows = []                                       # (centroid_dist, cluster|None)
        for pos, (_rnk, f) in enumerate(top, 1):
            idx = file_idx.get(f)
            cd = _dist(C[idx], crystal) if (idx is not None and crystal is not None) else np.nan
            rm = float(file_rmsd.get(f, np.nan))
            clab = int(labels[idx]) if idx is not None else None
            in_correct = (bool(clab == correct_label)
                          if (idx is not None and correct_label is not None) else np.nan)
            per_rank.append((tool, pos, cd, rm, in_correct, clab))
            rows.append((cd, clab))
        # ── top-5 concentration summary: first ``max_rank`` positions only ──
        dists = [cd for cd, _ in rows[:max_rank]]
        clusters = [cl for _, cl in rows[:max_rank] if cl is not None]
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
            top5_n=min(len(top), max_rank),
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
                    thr: float, do_placement: bool, mode_rmsd_thr: float,
                    rmsd_cap: int, site_method: str = "threshold",
                    pocket_radius: float = 8.0, rank_by: str = "consensus",
                    n_boot: int = 25) -> dict:
    # ── pose centroid array (drop poses with no centroid) ────────────────
    rows = [(r["pose_file"], r["tool"]) for _, r in sub.iterrows()
            if cents.get(r["pose_file"]) is not None]
    if len(rows) < 2:
        return {"protein": cid, "n_poses": len(rows), "skipped": True}
    files = [f for f, _ in rows]
    tools = [t for _, t in rows]
    C = np.asarray([cents[f] for f in files])
    rmsd_map = (dict(zip(sub["pose_file"], pd.to_numeric(sub.get("rmsd"), errors="coerce")))
                if "rmsd" in sub else {})

    # ── primary: centroid clustering (k=1-capable) -> candidate pockets ──
    # cluster_sites allows a single cluster when the tools converge on one site
    # (silhouette-max KMedoids never could), and pockets are ranked by cross-tool
    # consensus rather than raw pose count.
    eff_rank, rank_src = _effective_ranks(sub)
    pose_w = _pose_weights(files, eff_rank)
    dm = _centroid_dm(C)
    labels, k, sil = cluster_sites(dm, C, site_method, pocket_radius)
    km_labels, km_k, km_sil = _select_k(dm, "kmedoids")        # legacy, for agreement
    hi_labels, hi_k, hi_sil = _select_k(dm, "agglomerative")
    pockets = pockets_from_labels(C, labels, tools, files=files,
                                  pose_weights=pose_w, rank_by=rank_by)

    # ── cluster-quality metrics on the primary partition ─────────────────
    internal = _internal_indices(C, labels)
    boot_stab, boot_nc = _bootstrap_stability(C, labels, site_method, pocket_radius, n_boot)
    from sklearn.metrics import adjusted_rand_score as _ari
    ari_pk = (float(_ari(labels, km_labels))
              if len(set(labels)) > 1 and len(set(km_labels)) > 1 else np.nan)

    # ── rank of each tool's poses in the crystal-closest cluster ─────────
    # The cluster whose center is nearest the crystal is the "correct" site.
    # For each tool we report the best (lowest) rank it assigns to a pose that
    # landed in that cluster — i.e. does the tool prioritize its near-native pose?
    correct_dist = correct_is_hit = np.nan
    pur_centroid = pur_rmsd = best_rmsd_in_correct = np.nan
    p1 = {}
    correct_label: Optional[int] = None
    rank_in_correct: Dict[str, Optional[int]] = {}
    n_in_correct: Dict[str, int] = {}
    correct_members: List[Tuple[str, Optional[int]]] = []
    if crystal is not None and pockets:
        cp = min(pockets, key=lambda p: _dist(p["center"], crystal))
        correct_label = int(cp["label"])
        correct_dist = round(_dist(cp["center"], crystal), 3)
        correct_is_hit = bool(correct_dist <= thr)
        pur_centroid, pur_rmsd, best_rmsd_in_correct = _cluster_purity(
            cp, C, crystal, files, rmsd_map, thr)
        p1 = _precision_at_1(pockets, crystal, thr)
        for t in TOOLS:
            rk = [eff_rank.get(files[i]) for i in cp["members"]
                  if tools[i] == t and eff_rank.get(files[i]) is not None]
            rank_in_correct[t] = int(min(rk)) if rk else None
            n_in_correct[t] = len(rk)
        # Every member pose of the crystal-closest cluster as (tool, effective rank),
        # for the near-native-cluster composition figure (rank make-up + tool consensus).
        correct_members = [(tools[i],
                            (int(eff_rank[files[i]]) if eff_rank.get(files[i]) is not None else None))
                           for i in cp["members"]]

    # ── top-5 ranked poses: per-rank distance + cluster consistency ──────
    file_idx = {f: i for i, f in enumerate(files)}
    per_rank_rows, top5 = _top5_analysis(sub, eff_rank, file_idx, labels, C, crystal,
                                         correct_label=correct_label)

    # ── tool-combination ensemble exploration (which tools to combine) ───
    ensembles = {name: _ensemble_stats(C, tools, crystal, subset, thr,
                                       site_method, pocket_radius)
                 for name, subset in ENSEMBLES.items()}

    # ── secondary: placement-aware binding-MODE analysis (bounded) ───────
    # A SEPARATE clustering on the in-place RMSD matrix (translation+orientation+
    # conformation), reported alongside — not blended into — the centroid sites.
    # ari/nmi_site_vs_mode = do poses that share a spatial site also share a mode?
    # modes_in_top_site = how many distinct modes live inside the biggest site
    # (>1 means the centroid axis is hiding e.g. flipped poses within one pocket).
    n_modes = mode_sil = mode_compact = mode_sep = np.nan
    ari_sm = nmi_sm = modes_in_top_site = np.nan
    mode_oracle_rmsd = mode_top_site_best_rmsd = np.nan
    if do_placement and len(files) <= rmsd_cap:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
        mols = [load_heavy_atom_mol(f) for f in files]
        rdm = _inplace_rmsd_dm(mols)
        mode_labels, n_modes, mode_sil = cluster_modes(rdm, mode_rmsd_thr)
        mq = _mode_quality(rdm, mode_labels)
        mode_compact, mode_sep = mq["compactness"], mq["separation"]
        if len(set(labels)) > 1 and len(set(mode_labels)) > 1:
            ari_sm = float(adjusted_rand_score(labels, mode_labels))
            nmi_sm = float(normalized_mutual_info_score(labels, mode_labels))
        big_site = Counter(labels.tolist()).most_common(1)[0][0]
        smembers = np.where(labels == big_site)[0]
        if len(smembers):
            modes_in_top_site = int(len(set(mode_labels[smembers].tolist())))
        # crystal (RMSD axis): best pose RMSD overall, and within the biggest site
        rm_all = np.array([rmsd_map.get(files[i], np.nan) for i in range(len(files))], dtype=float)
        if np.isfinite(rm_all).any():
            mode_oracle_rmsd = round(float(np.nanmin(rm_all)), 3)
        if len(smembers) and np.isfinite(rm_all[smembers]).any():
            mode_top_site_best_rmsd = round(float(np.nanmin(rm_all[smembers])), 3)

    # ── per-tool primary site (largest cluster of that tool's centroids) ──
    tool_site: Dict[str, Optional[np.ndarray]] = {}
    for t in TOOLS:
        ti = [i for i, tt in enumerate(tools) if tt == t]
        if not ti:
            tool_site[t] = None
            continue
        ct = C[ti]
        if len(ct) < 2:
            tool_site[t] = ct.mean(axis=0)
        else:
            tl, _, _ = cluster_sites(_centroid_dm(ct), ct, site_method, pocket_radius)
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
        "site_method": site_method, "rank_by": rank_by, "k_is_one": bool(k == 1),
        "n_clusters": len(pockets), "silhouette": round(sil, 3) if sil == sil else np.nan,
        # legacy silhouette-KMedoids (never k=1) kept for method-agreement comparison
        "km_n_clusters": int(km_k), "km_silhouette": round(km_sil, 3) if km_sil == km_sil else np.nan,
        "ari_primary_vs_kmedoids": round(ari_pk, 3) if ari_pk == ari_pk else np.nan,
        "hier_n_clusters": int(len(set(hi_labels))), "hier_silhouette": round(hi_sil, 3) if hi_sil == hi_sil else np.nan,
        # placement-aware (in-place RMSD) binding-MODE analysis, alongside centroid sites
        "n_modes": n_modes, "mode_silhouette": round(mode_sil, 3) if mode_sil == mode_sil else np.nan,
        "mode_compactness": mode_compact, "mode_separation": mode_sep,
        "ari_site_vs_mode": round(ari_sm, 3) if ari_sm == ari_sm else np.nan,
        "nmi_site_vs_mode": round(nmi_sm, 3) if nmi_sm == nmi_sm else np.nan,
        "modes_in_top_site": modes_in_top_site,
        "mode_oracle_rmsd": mode_oracle_rmsd,
        "mode_top_site_best_rmsd": mode_top_site_best_rmsd,
        # internal validity indices on the primary partition
        "ch_score": internal["ch_score"], "db_score": internal["db_score"],
        "compactness": internal["compactness"], "separation": internal["separation"],
        # bootstrap cluster stability (mean best-Jaccard of base clusters)
        "boot_stability": boot_stab, "boot_n_clusters": boot_nc,
        "has_crystal": crystal is not None,
        # crystal-closest-cluster purity + best member RMSD (tight-vs-diffuse)
        "correct_cluster_purity_centroid": pur_centroid,
        "correct_cluster_purity_rmsd": pur_rmsd,
        "correct_cluster_best_rmsd": best_rmsd_in_correct,
        # precision@1: rank-1 pocket vs crystal under each ranking rule
        **p1,
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
        "_correct_members": correct_members, "_rank_src": rank_src, "_ensembles": ensembles,
        "_per_rank": [(cid, t, rp, cd, rm, ic, cl)
                      for (t, rp, cd, rm, ic, cl) in per_rank_rows],
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


# Oracle-distance columns feeding the ECDF / hit-rate panels, in draw order.
_SRC_ORACLE = {
    "autodock": "autodock_oracle_centroid_dist",
    "diffdock": "diffdock_oracle_centroid_dist",
    "equibind": "equibind_oracle_centroid_dist",
    "ensemble": "ens_oracle_dist",
    "fpocket": "fpocket_oracle_dist",
    "p2rank": "p2rank_oracle_dist",
}


def _panel_oracle_ecdf(ax, d, thr, detail_xlim=None):
    """ECDF of oracle center-to-crystal distance per source.

    ``detail_xlim=None`` draws the full 0–40 Å view with the 0–8 Å zoom inset
    (the combined-figure / standalone Panel-A behaviour). ``detail_xlim=(lo, hi)``
    draws ONLY the zoomed lo–hi Å view (the standalone detail graph). Returns the
    per-source ECDF dict so callers can reuse it.
    """
    ecdf = {}
    for name, col in _SRC_ORACLE.items():
        if col not in d:
            continue
        v = np.sort(pd.to_numeric(d[col], errors="coerce").dropna().to_numpy())
        if v.size:
            ecdf[name] = (v, np.linspace(0, 1, v.size))
    detail = detail_xlim is not None
    for name, (vx, vy) in ecdf.items():
        ax.plot(vx, vy, label=name, lw=2)
    ax.axvline(thr, color="k", ls="--", lw=0.8)
    ax.set_ylabel("fraction of complexes")
    ax.set_xlabel("center-to-crystal distance (Å)")
    ax.legend(fontsize=8)
    if detail:
        lo, hi = detail_xlim
        ax.set_xlim(lo, hi); ax.set_ylim(0, 1)
        ax.set_title(f"Oracle center-to-crystal distance (ECDF, detail {lo:g}–{hi:g} Å)")
        ax.grid(alpha=0.25); ax.set_axisbelow(True)
        return ecdf
    ax.set_xlim(0, 40)
    ax.set_title("Oracle center-to-crystal distance (ECDF)")
    # Detail inset: the curves bunch up near 0 Å — zoom the 0–8 Å region so the
    # docking sources (which nearly all sit there) are distinguishable.
    axin = ax.inset_axes([0.46, 0.12, 0.5, 0.55])
    for name, (vx, vy) in ecdf.items():
        axin.plot(vx, vy, lw=1.6)
    axin.axvline(thr, color="k", ls="--", lw=0.8)
    axin.set_xlim(0, 8); axin.set_ylim(0, 1)
    axin.set_title("detail: 0–8 Å", fontsize=8)
    axin.tick_params(labelsize=7)
    ax.indicate_inset_zoom(axin, edgecolor="grey")
    return ecdf


def _panel_hit_rate(ax, d, thr):
    """Hit-rate@thr per source (oracle center within thr Å of the crystal)."""
    def rate(col):
        v = pd.to_numeric(d[col], errors="coerce").dropna()
        return float((v <= thr).mean()) if len(v) else np.nan
    labels = list(_SRC_ORACLE)
    oracle_rates = [rate(_SRC_ORACLE[s]) for s in labels]
    x = np.arange(len(labels))
    ax.bar(x, oracle_rates, color=[TOOL_COLORS.get(s, "#8172B3") for s in labels])
    for xi, v in zip(x, oracle_rates):
        if v == v:
            ax.text(xi, v + 0.01, f"{v:.0%}", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylim(0, 1.18); ax.set_ylabel("fraction within thr")
    ax.set_title(f"Finds the true site (oracle ≤ {thr:g} Å)")


def _panel_triangulation(ax, d, thr):
    """Triangulation: which of docking / fpocket / p2rank finds the true site."""
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
    ax.bar(range(len(keys)), [v / tot for v in vals], color="#55A868")
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([names[k] for k in keys], rotation=30, ha="right", fontsize=8)
    ax.set_title(f"Triangulation — who finds the true site (≤ {thr:g} Å)")
    ax.set_ylabel("fraction of complexes")
    for xi, v in enumerate(vals):
        ax.text(xi, v / tot + 0.005, str(v), ha="center", fontsize=8)


def _panel_cross_tool(ax, d, thr):
    """Cross-tool agreement: histogram of inter-tool top-site center distances."""
    for col, name, c in [("au_di_dist", "AD–DD", "#4C72B0"),
                         ("au_eq_dist", "AD–EB", "#55A868"),
                         ("di_eq_dist", "DD–EB", "#C44E52")]:
        if col in d:
            v = pd.to_numeric(d[col], errors="coerce").dropna().to_numpy()
            if v.size:
                ax.hist(v, bins=np.linspace(0, 40, 21), histtype="step",
                        lw=2, label=f"{name} (med {np.median(v):.1f})", color=c)
    ax.axvline(thr, color="k", ls="--", lw=0.8)
    ax.set_title("Cross-tool agreement — top-site center distance")
    ax.set_xlabel("distance (Å)"); ax.set_ylabel("complexes")
    ax.legend(fontsize=8)


def _fig_summary(df, thr, out_dir):
    """Docked-clusters vs fpocket/p2rank vs crystal overview.

    Writes the combined 2×2 ``crystal_pocket_summary.png`` AND each of its four
    panels as its own standalone PNG (single-panel figures carry no (a)/(b) label),
    plus a 0–5 Å detail of the oracle-distance ECDF. Returns the list of paths.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = df[df["has_crystal"]]
    saved = []

    # ── combined 2×2 (kept for an at-a-glance overview) ──────────────────
    fig, ax = plt.subplots(2, 2, figsize=(14, 10))
    _panel_oracle_ecdf(ax[0, 0], d, thr)
    _panel_hit_rate(ax[0, 1], d, thr)
    _panel_triangulation(ax[1, 0], d, thr)
    _panel_cross_tool(ax[1, 1], d, thr)
    _label_panels(ax)
    fig.suptitle("Docked clusters vs fpocket/p2rank vs crystal "
                 f"(n={len(d)} complexes, thr={thr:g} Å)", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "crystal_pocket_summary.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    saved.append(p)

    # ── the four panels, each as its own standalone graph ────────────────
    panels = [
        ("summary_oracle_distance_ecdf.png", _panel_oracle_ecdf, (7.6, 5.6)),
        ("summary_hit_rate.png", _panel_hit_rate, (7.0, 5.2)),
        ("summary_triangulation.png", _panel_triangulation, (7.2, 5.2)),
        ("summary_cross_tool_agreement.png", _panel_cross_tool, (7.0, 5.2)),
    ]
    for fname, drawer, figsize in panels:
        fig, a = plt.subplots(figsize=figsize)
        drawer(a, d, thr)
        fig.tight_layout()
        p = out_dir / fname
        fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
        saved.append(p)

    # ── new graph: the Panel-A detail, zoomed to 0–5 Å (no inset) ────────
    fig, a = plt.subplots(figsize=(7.0, 5.2))
    _panel_oracle_ecdf(a, d, thr, detail_xlim=(0, 5))
    fig.tight_layout()
    p = out_dir / "oracle_distance_ecdf_detail_0-5A.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    saved.append(p)
    return saved


def _fig_rank_in_correct(df, thr, out_dir, eq_variant=None):
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
                 f"(true-site complexes, n={len(d)}; {_equibind_rank_note(eq_variant)})",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "rank_in_correct_cluster.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


_TOOL_DISPLAY = {"autodock": "AutoDock Vina", "diffdock": "DiffDock", "equibind": "EquiBind"}
# Rank bands for the near-native-cluster composition figure (ordinal -> sequential ramp).
_RANK_BANDS = [("rank 1", 1, 1), ("rank 2–5", 2, 5), ("rank 6–15", 6, 15), ("rank ≥ 16", 16, 10 ** 9)]
_RANK_BAND_COLORS = ["#08519c", "#3182bd", "#6baed6", "#c6dbef"]   # dark→light blue
_CONSENSUS_COLORS = {1: "#bdbdbd", 2: "#737373", 3: "#252525"}     # 1→3 tools, light→dark


def _fig_near_native_composition(results, thr, out_dir, eq_variant=None):
    """Composition of the near-native cluster.

    The near-native cluster is the crystal-closest cluster that actually holds a pose
    ≤ 2 Å from the crystal (i.e. contains near-native poses). Across those complexes:
      (left)  which RANKED poses of each tool land in that cluster — mean number of a
              tool's poses in it, stacked by rank band (does a tool put its TOP poses
              on the true site, or only deep ones?);
      (right) CONSENSUS — how many distinct tools have ≥ 1 pose in that cluster.
    EquiBind is ranked by smina affinity (see _effective_ranks) when that variant is used.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nn = [r for r in results
          if not r.get("skipped") and r.get("has_crystal") and r.get("_correct_members")
          and r.get("correct_cluster_best_rmsd") == r.get("correct_cluster_best_rmsd")
          and float(r["correct_cluster_best_rmsd"]) <= 2.0]
    if not nn:
        return None
    N = len(nn)

    band_sum = {t: [0.0] * len(_RANK_BANDS) for t in TOOLS}   # summed poses per band
    coverage = {t: 0 for t in TOOLS}                          # complexes where present
    consensus = {1: 0, 2: 0, 3: 0}
    for r in nn:
        present = set()
        per_tool = {t: [0] * len(_RANK_BANDS) for t in TOOLS}
        for tool, rk in r["_correct_members"]:
            if tool not in TOOLS:
                continue
            present.add(tool)
            if rk is None:
                continue
            for bi, (_lbl, lo, hi) in enumerate(_RANK_BANDS):
                if lo <= rk <= hi:
                    per_tool[tool][bi] += 1
                    break
        for t in TOOLS:
            for bi in range(len(_RANK_BANDS)):
                band_sum[t][bi] += per_tool[t][bi]
            if t in present:
                coverage[t] += 1
        nt = len(present & set(TOOLS))
        if nt in consensus:
            consensus[nt] += 1
    band_mean = {t: [s / N for s in band_sum[t]] for t in TOOLS}
    subtitle = (f"near-native cluster = crystal-closest cluster holding a ≤ 2 Å pose  "
                f"(n={N} complexes; {_equibind_rank_note(eq_variant)})")
    saved = []

    # ── (1) rank composition per tool (stacked by rank band) — own figure ──
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    x = np.arange(len(TOOLS))
    bottom = np.zeros(len(TOOLS))
    for bi, (lbl, _lo, _hi) in enumerate(_RANK_BANDS):
        vals = np.array([band_mean[t][bi] for t in TOOLS])
        ax.bar(x, vals, 0.6, bottom=bottom, color=_RANK_BAND_COLORS[bi],
               edgecolor="white", linewidth=0.6, label=lbl, zorder=3)
        bottom += vals
    for xi, t in zip(x, TOOLS):
        tot = float(sum(band_mean[t]))
        ax.text(xi, tot + 0.04,
                f"{tot:.2f} poses\n{100 * coverage[t] / N:.0f}% of complexes",
                ha="center", va="bottom", fontsize=8, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in TOOLS])
    ax.set_ylabel("mean # of a tool's poses in the near-native cluster")
    ax.set_ylim(0, max(0.5, float(bottom.max()) * 1.3))
    ax.set_title("Which ranked poses of each tool reach the near-native cluster")
    ax.legend(title="pose rank", fontsize=8, loc="upper right")
    ax.grid(axis="y", ls="--", lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    fig.suptitle(subtitle, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = out_dir / "near_native_cluster_rank_composition.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    saved.append(p)

    # ── (2) cross-tool consensus in the near-native cluster — own figure ──
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    cx = [1, 2, 3]
    cy = [100 * consensus[k] / N for k in cx]
    bars = ax.bar(cx, cy, 0.6, color=[_CONSENSUS_COLORS[k] for k in cx],
                  edgecolor="black", zorder=3)
    for b, val in zip(bars, cy):
        ax.text(b.get_x() + b.get_width() / 2, val + 1.2, f"{val:.0f}%",
                ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_xticks(cx); ax.set_xticklabels(["1 tool", "2 tools", "3 tools"])
    ax.set_ylim(0, max(cy) * 1.25 if cy else 100)
    ax.set_ylabel("% of near-native-cluster complexes")
    ax.set_title("How many tools have poses in the near-native cluster")
    ax.grid(axis="y", ls="--", lw=0.5, alpha=0.5); ax.set_axisbelow(True)
    fig.suptitle(subtitle, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = out_dir / "near_native_cluster_consensus.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    saved.append(p)
    return saved


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
    """Top-5 ranked poses — one standalone PNG per graph. Returns the list of
    written paths (skipping any figure that had no data):
      * top5_centroid_distance_by_rank.png — per-rank centroid distance to crystal
      * top5_rmsd_by_rank.png              — per-rank RMSD to crystal
      * top5_pose_concentration.png        — how concentrated the top-5 are (was a bar)
      * top5_distinct_sites.png            — # distinct sites in the top-5 (was a bar)
      * top5_best_rank.png                 — which of the top-5 is closest to crystal
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    saved = []
    ranks = list(range(1, 6))

    # ── per-rank distance / RMSD to the crystal ligand ───────────────────
    # Median is the primary (solid) line — robust to EquiBind's far/garbage
    # poses that otherwise blow the scale; mean (dashed) is overlaid and clipped.
    caps = {"centroid_dist": 15.0, "rmsd": 15.0}
    if not df_rank.empty:
        for col_name, lab, fname in (
                ("centroid_dist", "centroid distance",
                 "top5_centroid_distance_by_rank.png"),
                ("rmsd", "RMSD", "top5_rmsd_by_rank.png")):
            cap = caps[col_name]
            fig, ax = plt.subplots(figsize=(6.6, 4.7))
            for t in TOOLS:
                g = df_rank[df_rank.tool == t]
                means, meds = [], []
                for k in ranks:
                    vk = pd.to_numeric(g[g["rank"] == k][col_name],
                                       errors="coerce").dropna()
                    means.append(vk.mean() if len(vk) else np.nan)
                    meds.append(vk.median() if len(vk) else np.nan)
                clip = lambda xs: [x if (x == x and x <= cap) else np.nan for x in xs]
                ax.plot(ranks, clip(meds), "-o", color=TOOL_COLORS[t], lw=2.2,
                        label=f"{t} (median)")
                ax.plot(ranks, clip(means), "--", color=TOOL_COLORS[t], lw=1,
                        alpha=0.55)
            ax.set_xticks(ranks); ax.set_ylim(0, cap)
            ax.set_xlabel("Pose rank (1 = tool's top pose)")
            ax.set_ylabel(f"{lab} to crystal (Å)")
            ax.set_title(f"Top-5 {lab} from the crystal ligand\n"
                         f"(solid = median, dashed = mean; capped at {cap:g} Å)")
            ax.grid(alpha=0.25); ax.set_axisbelow(True); ax.legend(fontsize=8)
            fig.tight_layout()
            p = out_dir / fname
            fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
            saved.append(p)

    # ── pose concentration & distinct sites — box + strip per tool.
    # These replace the old grouped bar-of-means: showing the full per-complex
    # distribution (not just a mean) is the point of "how concentrated".
    d = df[df["has_crystal"]] if "has_crystal" in df else df
    rng = np.random.default_rng(0)

    def _dist_fig(col_suffix, fname, title, ylabel, ylim, pct):
        series = [(pd.to_numeric(d.get(f"{t}_{col_suffix}"), errors="coerce")
                   .dropna().to_numpy() if f"{t}_{col_suffix}" in d else np.array([]))
                  for t in TOOLS]
        if not any(len(s) for s in series):
            return
        fig, ax = plt.subplots(figsize=(6.4, 4.7))
        positions = np.arange(1, len(TOOLS) + 1)
        bp = ax.boxplot([s if len(s) else [np.nan] for s in series],
                        positions=positions, widths=0.55, showfliers=False,
                        patch_artist=True, medianprops=dict(color="k", lw=1.5),
                        whiskerprops=dict(color="#888"), capprops=dict(color="#888"))
        for patch, t in zip(bp["boxes"], TOOLS):
            patch.set_facecolor(TOOL_COLORS[t]); patch.set_alpha(0.35)
            patch.set_edgecolor(TOOL_COLORS[t])
        mean_labeled = False
        for pos, vals, t in zip(positions, series, TOOLS):
            if not len(vals):
                continue
            jit = rng.uniform(-0.16, 0.16, size=len(vals))
            ax.scatter(pos + jit, vals, s=10, color=TOOL_COLORS[t], alpha=0.45,
                       edgecolor="none")
            ax.scatter([pos], [np.mean(vals)], marker="D", s=42, color="k",
                       zorder=5, label=None if mean_labeled else "mean")
            mean_labeled = True
            ax.text(pos, ylim[1] * 0.97,
                    (f"{np.mean(vals):.0%}" if pct else f"{np.mean(vals):.1f}"),
                    ha="center", va="top", fontsize=9, fontweight="bold")
        ax.set_xticks(positions); ax.set_xticklabels(TOOLS)
        ax.set_ylim(*ylim); ax.set_ylabel(ylabel); ax.set_title(title)
        ax.grid(alpha=0.25, axis="y"); ax.set_axisbelow(True)
        ax.legend(fontsize=8, loc="lower right")
        fig.tight_layout()
        p = out_dir / fname
        fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
        saved.append(p)

    _dist_fig("top5_modal_frac", "top5_pose_concentration.png",
              "Pose concentration — fraction of a tool's top-5 in one cluster\n"
              "(higher = poses agree on one site; diamond = mean)",
              "Fraction of top-5 poses in the modal cluster", (0, 1.08), pct=True)
    _dist_fig("top5_n_clusters", "top5_distinct_sites.png",
              "Number of distinct sites among a tool's top-5 poses\n"
              "(higher = poses spread across more sites; diamond = mean)",
              "Distinct clusters among the top-5", (0, 5.4), pct=False)

    # ── which of the top-5 is closest to the crystal (a rank distribution,
    #    ordinal in rank -> a line per tool reads the trend more directly) ──
    if not d.empty:
        fig, ax = plt.subplots(figsize=(6.6, 4.7))
        any_data = False
        for t in TOOLS:
            br = (pd.to_numeric(d.get(f"{t}_top5_best_rank"), errors="coerce").dropna()
                  if f"{t}_top5_best_rank" in d else pd.Series(dtype=float))
            any_data = any_data or bool(len(br))
            fracs = [float((br == k).mean()) if len(br) else np.nan for k in ranks]
            ax.plot(ranks, fracs, "-o", color=TOOL_COLORS[t], lw=2.2, markersize=7,
                    markeredgecolor="black", markeredgewidth=0.5,
                    label=f"{t} (n={len(br)})")
        if any_data:
            ax.set_xticks(ranks)
            ax.set_xlabel("Rank of the top-5 pose closest to the crystal")
            ax.set_ylabel("Fraction of complexes")
            ax.set_ylim(0, None)
            ax.set_title("Which of the top-5 poses is closest to the crystal?")
            ax.grid(alpha=0.25); ax.set_axisbelow(True)
            ax.legend(fontsize=8)
            fig.tight_layout()
            p = out_dir / "top5_best_rank.png"
            fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
            saved.append(p)
        else:
            plt.close(fig)

    return saved


def _fig_rank_in_crystal_cluster(df_rank, out_dir, eq_variant=None):
    """Do a tool's n-th ranked pose tend to land in the crystal-closest cluster?

    For each tool and each rank position n (1..up to _RANK_PROFILE_MAX), the fraction
    of complexes whose n-th ranked pose falls in the crystal-closest cluster
    (``in_correct_cluster`` from the per-rank table). A high, rank-1-peaked line means
    the tool prioritises poses at the true site; a flat line means rank carries no
    site information. Note the per-rank denominator shrinks at deep ranks (a tool that
    emits fewer poses — e.g. AutoDock's Vina modes — has no pose there), so its line
    simply stops; the legend n is the rank-1 count.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if df_rank.empty or "in_correct_cluster" not in df_rank:
        return None
    max_present = pd.to_numeric(df_rank.get("rank"), errors="coerce").max()
    hi = int(min(_RANK_PROFILE_MAX, max_present)) if max_present == max_present else 5
    ranks = list(range(1, max(hi, 1) + 1))
    _truth = {True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    any_data = False
    for t in TOOLS:
        g = df_rank[df_rank["tool"] == t]
        fracs, n1 = [], 0
        for k in ranks:
            s = g[g["rank"] == k]["in_correct_cluster"].map(_truth).dropna()
            fracs.append(float(s.mean()) if len(s) else np.nan)
            if k == 1:
                n1 = len(s)
        if not np.isfinite(fracs).any():
            continue
        any_data = True
        ax.plot(ranks, fracs, "-o", color=TOOL_COLORS[t], lw=2.2, markersize=6,
                markeredgecolor="black", markeredgewidth=0.5,
                label=f"{t} (n={n1})")
    if not any_data:
        plt.close(fig)
        return None
    ax.set_xticks(ranks); ax.set_xlim(0.5, ranks[-1] + 0.5); ax.set_ylim(0, 1.05)
    ax.set_xlabel("Pose rank (1 = tool's top pose)")
    ax.set_ylabel("Fraction of complexes in the crystal-closest cluster")
    ax.set_title("Do a tool's top-ranked poses land in the crystal-closest cluster?\n"
                 f"({_equibind_rank_note(eq_variant)})")
    ax.grid(alpha=0.25); ax.set_axisbelow(True)
    ax.legend(fontsize=8)
    fig.tight_layout()
    p = out_dir / "rank_in_crystal_cluster_by_rank.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_rank_cluster_membership(df_rank, out_dir, eq_variant=None):
    """Per rank, where each tool's rank-k pose lands: IN the crystal-closest cluster
    vs a DIFFERENT cluster — the in/out split of rank_in_crystal_cluster_by_rank.

    One subplot per tool; each rank is a 100%-stacked bar over ALL complexes
    (n = crystal complexes), split into three shares that sum to 1:
      * in crystal cluster   (rank-k pose is in the crystal-closest cluster)
      * different cluster     (rank-k pose exists but is in another cluster)
      * no pose at this rank   (the tool emitted < k poses for that complex)
    The 'no pose' share makes the shrinking deep-rank denominator explicit (it grows
    for tools with few poses, e.g. AutoDock's Vina modes) instead of hiding it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    if df_rank.empty or "in_correct_cluster" not in df_rank:
        return None
    max_present = pd.to_numeric(df_rank.get("rank"), errors="coerce").max()
    hi = int(min(_RANK_PROFILE_MAX, max_present)) if max_present == max_present else 5
    ranks = list(range(1, max(hi, 1) + 1))
    total = int(df_rank["protein"].nunique())
    if total == 0:
        return None
    _truth = {True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}
    C_IN, C_OUT, C_ABS = "#55A868", "#C44E52", "#D9D9D9"
    fig, axes = plt.subplots(1, len(TOOLS), figsize=(5.0 * len(TOOLS), 4.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, t in zip(axes, TOOLS):
        g = df_rank[df_rank["tool"] == t]
        f_in, f_out, f_abs = [], [], []
        for k in ranks:
            m = g[g["rank"] == k]["in_correct_cluster"].map(_truth).dropna()
            nin = float((m == 1.0).sum()); nout = float((m == 0.0).sum())
            f_in.append(nin / total); f_out.append(nout / total)
            f_abs.append(max(0.0, 1.0 - (nin + nout) / total))
        f_in = np.array(f_in); f_out = np.array(f_out); f_abs = np.array(f_abs)
        ax.bar(ranks, f_in, width=0.82, color=C_IN)
        ax.bar(ranks, f_out, width=0.82, bottom=f_in, color=C_OUT)
        ax.bar(ranks, f_abs, width=0.82, bottom=f_in + f_out, color=C_ABS)
        ax.text(ranks[0], min(f_in[0] + 0.02, 0.98), f"{f_in[0]:.0%}", ha="center",
                va="bottom", fontsize=8, color="#2f6b45", fontweight="bold")
        ax.set_title(_TOOL_DISPLAY.get(t, t), color=TOOL_COLORS[t], fontweight="bold")
        ax.set_xlabel("Pose rank (1 = top pose)")
        ax.set_xticks(ranks); ax.set_xlim(0.4, ranks[-1] + 0.6); ax.set_ylim(0, 1.0)
        ax.grid(axis="y", alpha=0.25); ax.set_axisbelow(True)
    axes[0].set_ylabel(f"Fraction of complexes (n={total})")
    handles = [Patch(color=C_IN, label="in crystal-closest cluster"),
               Patch(color=C_OUT, label="in a different cluster"),
               Patch(color=C_ABS, label="no pose at this rank")]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, frameon=False)
    fig.suptitle("Where each tool's rank-k pose lands: crystal-closest cluster vs elsewhere\n"
                 f"({_equibind_rank_note(eq_variant)})", fontsize=12)
    fig.tight_layout(rect=(0, 0.055, 1, 0.93))
    p = out_dir / "rank_cluster_membership_stacked.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_rank_consensus(df_rank, out_dir, eq_variant=None):
    """Tool consensus across the top ranks: at each rank k, do AutoDock, DiffDock and
    EquiBind put their rank-k pose in the SAME cluster, and is it the crystal one?

    One 100%-stacked bar per rank k (1..up to _RANK_PROFILE_MAX), over the complexes
    where all three tools have a rank-k pose (n printed above each bar), split into four
    mutually exclusive outcomes:
      * all 3 → crystal cluster        (consensus on the true site)
      * all 3 → same other cluster     (consensus, wrong site)
      * split, >=1 in crystal cluster  (partial — some agree with the crystal site)
      * split, none in crystal cluster (scattered — all three miss)
    Needs the per-pose ``cluster`` label (added to per_rank_distance.csv) so that
    'same cluster' can be told apart from 'both merely out of the crystal cluster'.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = {"tool", "rank", "in_correct_cluster", "cluster"}
    if df_rank.empty or not need <= set(df_rank.columns):
        return None
    max_present = pd.to_numeric(df_rank.get("rank"), errors="coerce").max()
    hi = int(min(_RANK_PROFILE_MAX, max_present)) if max_present == max_present else 5
    ranks = list(range(1, max(hi, 1) + 1))
    _truth = {True: True, False: False, "True": True, "False": False}
    cats = ["all_crystal", "all_other", "split_some", "split_none"]
    colors = {"all_crystal": "#55A868", "all_other": "#8C8C8C",
              "split_some": "#DD8452", "split_none": "#C44E52"}
    lab = {"all_crystal": "all 3 → crystal cluster",
           "all_other": "all 3 → same other cluster",
           "split_some": "split — ≥1 in crystal cluster",
           "split_none": "split — none in crystal cluster"}
    d = df_rank[pd.to_numeric(df_rank["rank"], errors="coerce").isin(ranks)]
    frac = {c: [] for c in cats}
    ns = []
    for k in ranks:
        dk = d[d["rank"] == k]
        by_prot = {}
        for row in dk.itertuples(index=False):
            cl = getattr(row, "cluster", None)
            if cl is None or (isinstance(cl, float) and cl != cl):
                continue
            ic = _truth.get(getattr(row, "in_correct_cluster", None), False)
            by_prot.setdefault(getattr(row, "protein"), {})[getattr(row, "tool")] = (int(cl), bool(ic))
        counts = {c: 0 for c in cats}; n = 0
        for rec in by_prot.values():
            if not all(t in rec for t in TOOLS):
                continue
            n += 1
            clset = {rec[t][0] for t in TOOLS}
            n_in = sum(rec[t][1] for t in TOOLS)
            if len(clset) == 1 and n_in == len(TOOLS):
                counts["all_crystal"] += 1
            elif len(clset) == 1:
                counts["all_other"] += 1
            elif n_in >= 1:
                counts["split_some"] += 1
            else:
                counts["split_none"] += 1
        ns.append(n)
        for c in cats:
            frac[c].append(counts[c] / n if n else np.nan)
    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    bottom = np.zeros(len(ranks))
    for c in cats:
        vals = np.array([f if f == f else 0.0 for f in frac[c]])
        ax.bar(ranks, vals, width=0.82, bottom=bottom, color=colors[c], label=lab[c])
        bottom += vals
    for xi, k in enumerate(ranks):
        ax.text(k, 1.008, f"{ns[xi]}", ha="center", va="bottom", fontsize=7, color="#555")
    ax.set_xticks(ranks); ax.set_xlim(0.4, ranks[-1] + 0.6); ax.set_ylim(0, 1.0)
    ax.set_xlabel("Pose rank (each tool's k-th ranked pose)")
    ax.set_ylabel("Fraction of complexes (all 3 tools present at rank k)")
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12),
              ncol=2, frameon=False)
    ax.grid(axis="y", alpha=0.25); ax.set_axisbelow(True)
    fig.suptitle("Do the three tools' rank-k poses land in the SAME cluster — the crystal one or another?\n"
                 f"(number above each bar = n complexes; {_equibind_rank_note(eq_variant)})",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    p = out_dir / "rank_consensus_by_rank.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


_TOOL_ABBR = {"autodock": "AD", "diffdock": "DD", "equibind": "EB"}
_TRUTH_MAP = {True: True, False: False, "True": True, "False": False}


def _fig_rank1_cluster_matrix(df_rank, out_dir, eq_variant=None):
    """Pairwise matrix: how often two tools' RANK-1 poses land in the SAME cluster.

    cell[i][j] = fraction of complexes where tool i's rank-1 pose and tool j's rank-1
    pose share a cluster label (any cluster). Each cell also prints the count
    (same / both-present) and, off-diagonal, how many of those shared clusters are
    the crystal-closest one. The diagonal is trivially 100%.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = {"tool", "rank", "cluster", "in_correct_cluster"}
    if df_rank.empty or not need <= set(df_rank.columns):
        return None
    d1 = df_rank[pd.to_numeric(df_rank["rank"], errors="coerce") == 1]
    r1 = {}
    for row in d1.itertuples(index=False):
        cl = getattr(row, "cluster", None)
        if cl is None or (isinstance(cl, float) and cl != cl):
            continue
        r1.setdefault(getattr(row, "protein"), {})[getattr(row, "tool")] = (
            int(cl), bool(_TRUTH_MAP.get(getattr(row, "in_correct_cluster", None), False)))
    T = list(TOOLS); nT = len(T)
    same = np.zeros((nT, nT)); both = np.zeros((nT, nT)); samex = np.zeros((nT, nT))
    for rec in r1.values():
        for a in range(nT):
            for b in range(nT):
                ta, tb = T[a], T[b]
                if ta in rec and tb in rec:
                    both[a, b] += 1
                    if rec[ta][0] == rec[tb][0]:
                        same[a, b] += 1
                        if rec[ta][1] and rec[tb][1]:
                            samex[a, b] += 1
    frac = np.divide(same, both, out=np.full((nT, nT), np.nan), where=both > 0)
    # three-way agreement: rank-1 poses of ALL tools sharing one cluster (by cluster
    # transitivity this is exactly "EB shares a cluster with BOTH AD and DD", etc.)
    n3 = all3_same = all3_crystal = two_agree = all_diff = 0
    for rec in r1.values():
        if not all(t in rec for t in T):
            continue
        n3 += 1
        clset = {rec[t][0] for t in T}
        if len(clset) == 1:
            all3_same += 1
            if all(rec[t][1] for t in T):
                all3_crystal += 1
        elif len(clset) == 2:
            two_agree += 1
        else:
            all_diff += 1
    tri = np.tril(np.ones((nT, nT), dtype=bool))          # lower triangle + diagonal
    cmap = plt.cm.YlGnBu.copy(); cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(7.0, 6.2))
    im = ax.imshow(np.ma.masked_where(~tri, frac), cmap=cmap, vmin=0, vmax=1)
    ax.set_xticks(range(nT)); ax.set_yticks(range(nT))
    ax.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in T], rotation=45, ha="right",
                       rotation_mode="anchor")
    ax.set_yticklabels([_TOOL_DISPLAY.get(t, t) for t in T])
    for a in range(nT):
        for b in range(nT):
            if b > a or both[a, b] == 0:
                continue
            pct = frac[a, b]
            txt = f"{pct:.0%}\n{int(same[a, b])}/{int(both[a, b])}"
            if a != b:
                txt += f"\n(crystal {int(samex[a, b])})"
            ax.text(b, a, txt, ha="center", va="center",
                    color="white" if pct > 0.55 else "black", fontsize=9)
    summ = ("All 3 rank-1 in ONE cluster\n"
            f"  {all3_same}/{n3} ({all3_same / n3:.0%})   → crystal {all3_crystal}\n"
            f"Only 2 tools share:  {two_agree}\n"
            f"All 3 differ:  {all_diff}")
    ax.text(0.97, 0.97, summ, transform=ax.transAxes, ha="right", va="top", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.5", fc="#f5f5f5", ec="#bbbbbb"))
    n_total = int(df_rank["protein"].nunique())
    ax.set_title("Rank-1 poses in the SAME cluster — pairwise (any cluster)\n"
                 f"(n={n_total} complexes; {_equibind_rank_note(eq_variant)})", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="fraction of complexes")
    fig.tight_layout()
    p = out_dir / "rank1_cluster_agreement_matrix.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_topN_crystal_matrix(df_rank, out_dir, eq_variant=None):
    """Buckets: pairwise matrices of both tools reaching the CRYSTAL cluster within
    the top-N ranked poses, for N = 1, 5, 10, 15.

    cell[i][j] (i != j) = fraction of complexes where tool i AND tool j each place at
    least one of their top-N poses in the crystal-closest cluster (so both are in the
    same crystal cluster). Diagonal cell[i][i] = fraction where tool i alone reaches
    the crystal cluster within top-N.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = {"tool", "rank", "in_correct_cluster"}
    if df_rank.empty or not need <= set(df_rank.columns):
        return None
    maxr = int(pd.to_numeric(df_rank["rank"], errors="coerce").max())
    buckets = [b for b in (1, 5, 10, 15) if b <= maxr] or [maxr]
    T = list(TOOLS); nT = len(T)
    # protein -> tool -> best (min) rank landing in the crystal cluster (else inf)
    reach = {}; present = {}
    for row in df_rank.itertuples(index=False):
        prot = getattr(row, "protein"); t = getattr(row, "tool")
        r = pd.to_numeric(getattr(row, "rank"), errors="coerce")
        present.setdefault(prot, set()).add(t)
        if _TRUTH_MAP.get(getattr(row, "in_correct_cluster", None), False) and r == r:
            m = reach.setdefault(prot, {})
            m[t] = min(m.get(t, np.inf), float(r))
    proteins = list(present.keys())
    n_total = len(proteins)
    tri = np.tril(np.ones((nT, nT), dtype=bool))          # lower triangle + diagonal
    cmap = plt.cm.YlGnBu.copy(); cmap.set_bad("white")
    nb = len(buckets)
    ncols = 2 if nb > 2 else nb
    nrows = int(np.ceil(nb / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 5.0 * nrows),
                             constrained_layout=True, squeeze=False)
    flat = axes.ravel()
    im = None
    for i, N in enumerate(buckets):
        ax = flat[i]; col = i % ncols
        cnt = np.zeros((nT, nT)); both = np.zeros((nT, nT))
        all3 = n_all = 0                     # all three tools reach crystal within top-N
        for prot in proteins:
            pres = present.get(prot, set()); rc = reach.get(prot, {})
            hit = {t: (rc.get(t, np.inf) <= N) for t in T}
            for a in range(nT):
                for b in range(nT):
                    ta, tb = T[a], T[b]
                    if ta in pres and tb in pres:
                        both[a, b] += 1
                        if (hit[ta] if a == b else (hit[ta] and hit[tb])):
                            cnt[a, b] += 1
            if all(t in pres for t in T):
                n_all += 1
                if all(hit[t] for t in T):
                    all3 += 1
        M = np.divide(cnt, both, out=np.full((nT, nT), np.nan), where=both > 0)
        im = ax.imshow(np.ma.masked_where(~tri, M), cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(range(nT)); ax.set_yticks(range(nT))
        ax.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in T], rotation=45,
                           ha="right", rotation_mode="anchor", fontsize=8)
        ax.set_yticklabels([_TOOL_DISPLAY.get(t, t) for t in T] if col == 0
                           else [""] * nT)
        for a in range(nT):
            for b in range(nT):
                if b > a or both[a, b] == 0:
                    continue
                ax.text(b, a, f"{M[a, b]:.0%}\n{int(cnt[a, b])}", ha="center",
                        va="center", color="white" if M[a, b] > 0.55 else "black",
                        fontsize=9)
        if n_all:
            ax.text(0.96, 0.96, "All 3 tools reach crystal:\n"
                    f"{all3}/{n_all} ({all3 / n_all:.0%})", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#bbbbbb"))
        ax.set_title(f"Top-{N}")
    for j in range(nb, nrows * ncols):       # hide any unused grid cell
        flat[j].axis("off")
    if im is not None:
        fig.colorbar(im, ax=list(flat), fraction=0.046, pad=0.02,
                     label="fraction of complexes")
    fig.suptitle("Both tools reach the CRYSTAL cluster within top-N poses "
                 "(diagonal = one tool alone)\n"
                 f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}", fontsize=12)
    p = out_dir / "topN_crystal_cluster_matrix.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_topN_crystal_reach_curves(df_rank, out_dir, eq_variant=None):
    """How the topN_crystal_cluster_matrix values EVOLVE with ranking depth N.

    Continuous companion to the four discrete matrix buckets: fraction of complexes
    reaching the crystal-closest cluster within top-N, for N = 1..up to
    _RANK_PROFILE_MAX. Left panel = each tool alone (the matrix diagonals); right
    panel = each tool pair jointly (the off-diagonals) plus all three (the per-panel
    three-way box). Monotonically non-decreasing, so the slope is the marginal gain
    from allowing one more pose.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = {"tool", "rank", "in_correct_cluster"}
    if df_rank.empty or not need <= set(df_rank.columns):
        return None
    maxr = pd.to_numeric(df_rank["rank"], errors="coerce").max()
    hi = int(min(_RANK_PROFILE_MAX, maxr)) if maxr == maxr else 5
    Ns = list(range(1, max(hi, 1) + 1))
    T = list(TOOLS)
    reach = {}; present = {}
    for row in df_rank.itertuples(index=False):
        prot = getattr(row, "protein"); t = getattr(row, "tool")
        r = pd.to_numeric(getattr(row, "rank"), errors="coerce")
        present.setdefault(prot, set()).add(t)
        if _TRUTH_MAP.get(getattr(row, "in_correct_cluster", None), False) and r == r:
            m = reach.setdefault(prot, {}); m[t] = min(m.get(t, np.inf), float(r))
    proteins = [p for p in present if all(t in present[p] for t in T)]
    n_total = len(proteins)
    if n_total == 0:
        return None

    def frac(subset, N):
        c = sum(all(reach.get(p, {}).get(t, np.inf) <= N for t in subset)
                for p in proteins)
        return c / n_total

    pairs = [("autodock", "diffdock"), ("autodock", "equibind"),
             ("diffdock", "equibind")]
    pair_colors = {pairs[0]: "#8172B3", pairs[1]: "#DD8452", pairs[2]: "#937860"}
    mark = [r for r in (1, 5, 10, 15) if r in Ns]      # ranks to label with the %

    def _ann(ax, y, color, dy):
        for N in mark:
            v = y[Ns.index(N)]
            ax.annotate(f"{v:.0%}", (N, v), textcoords="offset points",
                        xytext=(0, dy), ha="center",
                        va="bottom" if dy >= 0 else "top", fontsize=7,
                        color=color, fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                                  alpha=0.55))

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13.0, 5.4), sharey=True)
    left_dy = {"autodock": -12, "diffdock": 9, "equibind": 9}   # AD near ceiling → below
    for t in T:
        y = [frac((t,), N) for N in Ns]
        axL.plot(Ns, y, "-o", color=TOOL_COLORS[t], lw=2.4, markersize=5,
                 label=_TOOL_DISPLAY.get(t, t))
        _ann(axL, y, TOOL_COLORS[t], left_dy[t])
    axL.set_title("Each tool alone (matrix diagonal)")
    axL.set_ylabel(f"Fraction of complexes reaching the crystal cluster (n={n_total})")
    pair_dy = {pairs[0]: 9, pairs[1]: 9, pairs[2]: -12}         # DD+EB below (≈ all-three)
    for pr in pairs:
        y = [frac(pr, N) for N in Ns]
        axR.plot(Ns, y, "--o", color=pair_colors[pr], lw=2, markersize=4,
                 label=f"{_TOOL_DISPLAY[pr[0]]} + {_TOOL_DISPLAY[pr[1]]}")
        _ann(axR, y, pair_colors[pr], pair_dy[pr])
    y3 = [frac(tuple(T), N) for N in Ns]
    axR.plot(Ns, y3, "-s", color="black", lw=2.6, markersize=5, label="All three tools")
    _ann(axR, y3, "black", -24)                                # further below DD+EB
    axR.set_title("Tool pairs and all three jointly (matrix off-diagonal)")
    for ax in (axL, axR):
        ax.set_xlabel("Ranking depth (number of top poses considered, N)")
        ax.set_xticks(Ns); ax.set_xlim(Ns[0] - 0.3, Ns[-1] + 0.3); ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.25); ax.set_axisbelow(True); ax.legend(fontsize=8)
    fig.suptitle("Reaching the crystal cluster vs ranking depth "
                 "(evolution of the top-N matrix values)\n"
                 f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    p = out_dir / "topN_crystal_reach_curves.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ── self-contained statistics for the homogeneity figure ─────────────────────
# No scikit-posthocs / statsmodels in the `vina` env, so the post-hoc tests
# (Dunn, Jonckheere-Terpstra, Cochran's Q, McNemar, Holm) are hand-rolled from
# scipy primitives. Each panel of crystal_cluster_homogeneity.png has a distinct
# data structure and therefore a distinct test; see _homogeneity_stats().

def _fmt_p(p) -> str:
    if p is None or p != p:
        return "n/a"
    if p < 1e-4:
        return "<1e-4"
    if p < 1e-3:
        return f"{p:.1e}"
    return f"{p:.3f}"


def _p_stars(p) -> str:
    if p is None or p != p:
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def _holm(pvals):
    """Holm-Bonferroni step-down adjusted p-values (order preserved)."""
    p = np.asarray(pvals, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (m - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    return adj


def _cochran_q(X):
    """Cochran's Q for k paired binary conditions (X is n x k of 0/1)."""
    from scipy.stats import chi2
    X = np.asarray(X, float)
    n, k = X.shape
    Cj = X.sum(axis=0)
    Ri = X.sum(axis=1)
    T = X.sum()
    denom = k * T - float(np.sum(Ri ** 2))
    if denom == 0:
        return np.nan, np.nan, k - 1
    Q = (k - 1) * (k * float(np.sum(Cj ** 2)) - T ** 2) / denom
    return float(Q), float(chi2.sf(Q, k - 1)), k - 1


def _mcnemar_exact(a, b):
    """Exact (binomial) McNemar for two paired binary vectors a, b."""
    from scipy.stats import binomtest
    a = np.asarray(a); b = np.asarray(b)
    n10 = int(np.sum((a == 1) & (b == 0)))
    n01 = int(np.sum((a == 0) & (b == 1)))
    disc = n10 + n01
    if disc == 0:
        return n10, n01, 1.0
    p = binomtest(min(n10, n01), disc, 0.5).pvalue
    return n10, n01, float(min(p, 1.0))


def _ablation_paired_stats(df):
    """Paired-complex significance for the precision@1 ranking ablation (panel F
    of cluster_quality_metrics.png).

    Every rule ranks the same clustering of the same complexes, so rank-1
    hit/miss is a fully paired binary design. Cochran's Q gives the omnibus (do
    the k rules differ at all?); an exact McNemar on every one of the C(k,2)
    pairs, Holm-adjusted across the whole family, gives the post-hoc — with the
    consensus-vs-legacy-size contrast pulled out for headline reporting. Each
    test uses listwise-complete complexes (a rule that could not produce a
    rank-1 pocket casts no vote; the omnibus needs all k present, each pair only
    its two). Returns ``None`` when the hit columns are absent (a no-crystal
    dataset) or too few complexes are scored to test.
    """
    if "has_crystal" in df:
        df = df[df["has_crystal"]]
    rules = ["size", "ntools", "tight", "confidence", "medoid"]
    cols = {r: f"p1_{r}_hit" for r in rules}
    if not all(c in df for c in cols.values()):
        return None

    def _bin(s):
        return pd.to_numeric(
            s.map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}),
            errors="coerce")

    M = pd.DataFrame({r: _bin(df[c]) for r, c in cols.items()})
    m_omni = M.dropna()
    if len(m_omni) < 3:
        return None
    Q, pQ, dfQ = _cochran_q(m_omni.to_numpy())
    prec = {r: float(m_omni[r].mean()) for r in rules}

    # every pairwise exact McNemar, oriented so 'a' is the higher-precision rule
    # of the pair (so a_wins >= b_wins reads as "a beats b"); Holm across all C(k,2).
    pairs = []
    for i in range(len(rules)):
        for j in range(i + 1, len(rules)):
            a, b = rules[i], rules[j]
            if prec[b] > prec[a]:
                a, b = b, a
            sub = M[[a, b]].dropna()
            n10, n01, p = _mcnemar_exact(sub[a].to_numpy(), sub[b].to_numpy())
            pairs.append({"a": a, "b": b, "a_wins": int(n10), "b_wins": int(n01),
                          "n_pairs": int(len(sub)), "p_raw": float(p)})
    for pr, padj in zip(pairs, _holm([pp["p_raw"] for pp in pairs])):
        pr["p_holm"] = float(padj)

    cvs = next((pp for pp in pairs if {pp["a"], pp["b"]} == {"ntools", "size"}), None)
    cons_vs_size = None
    if cvs is not None:
        c_wins = cvs["a_wins"] if cvs["a"] == "ntools" else cvs["b_wins"]
        s_wins = cvs["b_wins"] if cvs["a"] == "ntools" else cvs["a_wins"]
        cons_vs_size = {"n_pairs": cvs["n_pairs"], "consensus_wins": int(c_wins),
                        "size_wins": int(s_wins),
                        "p_raw": cvs["p_raw"], "p_holm": cvs["p_holm"]}
    return {
        "n_omnibus": int(len(m_omni)), "k": len(rules),
        "cochran_q": round(float(Q), 3), "df": int(dfQ), "p_omnibus": float(pQ),
        "precision_at_1": {r: round(prec[r], 4) for r in rules},
        "pairwise_mcnemar_holm": pairs,
        "consensus_vs_size_mcnemar": cons_vs_size,
    }


_ABLATION_RULE_ABBR = {"size": "size", "ntools": "cons", "tight": "tght",
                       "confidence": "conf", "medoid": "medo"}


def _ablation_pairwise_matrix_text(st):
    """Render the full Holm-adjusted pairwise McNemar family (``st`` from
    :func:`_ablation_paired_stats`) as a compact lower-triangular star matrix,
    rules ordered by precision@1 (best first). Same significance-star vocabulary
    as the sibling homogeneity figure's ``_pair_stars``."""
    pairs = st.get("pairwise_mcnemar_holm") or []
    prec = st.get("precision_at_1") or {}
    if not pairs or not prec:
        return None
    order = sorted(prec, key=lambda r: prec[r], reverse=True)   # best precision first
    star = {frozenset((p["a"], p["b"])): _p_stars(p["p_holm"]) for p in pairs}
    ab = _ABLATION_RULE_ABBR
    cols = order[:-1]
    lines = ["Holm-adj. pairwise McNemar",
             "     " + "".join(f"{ab.get(c, c):>5}" for c in cols)]
    for ri, r in enumerate(order[1:], start=1):
        cells = "".join(f"{star.get(frozenset((r, c)), ''):>5}" for c in order[:ri])
        lines.append(f"{ab.get(r, r):>4} {cells}")
    lines.append("*** p<.001  ** p<.01  * p<.05")
    return "\n".join(lines)


def _wilcoxon_rankbiserial(x, y):
    """Wilcoxon signed-rank p + matched-pairs rank-biserial effect size."""
    from scipy.stats import rankdata, wilcoxon
    d = np.asarray(x, float) - np.asarray(y, float)
    d = d[d != 0]
    n = len(d)
    if n < 1:
        return np.nan, np.nan, 0
    r = rankdata(np.abs(d))
    r_plus = float(np.sum(r[d > 0])); r_minus = float(np.sum(r[d < 0]))
    total = r_plus + r_minus
    rb = (r_plus - r_minus) / total if total else np.nan
    try:
        p = float(wilcoxon(x, y, zero_method="wilcox").pvalue)
    except Exception:
        p = np.nan
    return float(rb), p, n


def _dunn(groups):
    """Dunn's post-hoc (tie-corrected z, uncorrected two-sided p) for k groups."""
    from scipy.stats import rankdata, norm
    data = np.concatenate([np.asarray(g, float) for g in groups])
    N = len(data)
    ranks = rankdata(data)
    _, counts = np.unique(data, return_counts=True)
    ties = float(np.sum(counts ** 3 - counts))
    sigma2 = (N * (N + 1) / 12.0) - ties / (12.0 * (N - 1))
    sizes, meanranks, idx = [], [], 0
    for g in groups:
        m = len(g)
        meanranks.append(float(np.mean(ranks[idx:idx + m])))
        sizes.append(m); idx += m
    out = []
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            se = np.sqrt(sigma2 * (1.0 / sizes[i] + 1.0 / sizes[j]))
            z = (meanranks[i] - meanranks[j]) / se if se else np.nan
            out.append((i, j, float(z), float(2 * norm.sf(abs(z)))))
    return out


def _jonckheere(groups):
    """Jonckheere-Terpstra trend test across ordered groups (normal approx)."""
    from scipy.stats import norm
    sizes = [len(g) for g in groups]
    N = sum(sizes)
    JT = 0.0
    for i in range(len(groups)):
        gi = np.asarray(groups[i], float)
        for j in range(i + 1, len(groups)):
            gj = np.asarray(groups[j], float)
            for a in gi:
                JT += float(np.sum(gj > a) + 0.5 * np.sum(gj == a))
    mean = (N ** 2 - sum(s ** 2 for s in sizes)) / 4.0
    var = (N ** 2 * (2 * N + 3)
           - sum(s ** 2 * (2 * s + 3) for s in sizes)) / 72.0
    z = (JT - mean) / np.sqrt(var) if var > 0 else np.nan
    return float(JT), float(z), float(2 * norm.sf(abs(z)))


def _consensus_perm_test(R, n_perm=10000, seed=42):
    """Permutation test: do tools reach the crystal cluster on the SAME
    complexes more than if their per-complex successes were independent?
    Each tool's reach indicator is shuffled across complexes independently
    (marginals preserved), breaking only the co-occurrence structure.
    Returns observed vs expected counts of {0,1,..,k tools} and p-values."""
    R = np.asarray(R, int)
    n, k = R.shape
    rng = np.random.default_rng(seed)

    def dist(mat):
        s = mat.sum(axis=1)
        return np.array([np.sum(s == c) for c in range(k + 1)], float)

    obs = dist(R)
    cols = [R[:, j].copy() for j in range(k)]
    null = np.zeros((n_perm, k + 1))
    for b in range(n_perm):
        perm = np.column_stack([rng.permutation(c) for c in cols])
        null[b] = dist(perm)
    exp = null.mean(axis=0)
    mask = exp > 0

    def stat(o):
        return float(np.sum((o[mask] - exp[mask]) ** 2 / exp[mask]))

    obs_stat = stat(obs)
    null_stat = np.array([stat(null[b]) for b in range(n_perm)])
    p_omni = (np.sum(null_stat >= obs_stat) + 1) / (n_perm + 1)
    allk_null = null[:, k]
    dev = abs(obs[k] - allk_null.mean())
    p_allk = (np.sum(np.abs(allk_null - allk_null.mean()) >= dev) + 1) / (n_perm + 1)
    return {"observed": obs.tolist(), "expected": exp.round(2).tolist(),
            "chi2_like_stat": round(obs_stat, 2), "p_omnibus": float(p_omni),
            "all_tools_obs": int(obs[k]), "all_tools_exp": round(float(allk_null.mean()), 2),
            "all_tools_p": float(p_allk), "n_perm": n_perm}


def _homogeneity_stats(comp, reach_R, internal, internal_prot, rad_by_n,
                       nt_radius_pairs, tools, tool_disp):
    """Compute every recommended test for crystal_cluster_homogeneity.png.

    Panel A (paired counts, k=3)     -> Friedman + Wilcoxon post-hoc (Holm);
                                        reaches (paired binary) -> Cochran's Q + McNemar.
    Panel B (consensus co-occurrence)-> independence-null permutation test.
    Panel C (tightness, unpaired)    -> Kruskal-Wallis + Dunn (Holm);
                                        + common-subset Friedman sensitivity check.
    Panel D (radius vs #tools, ord.) -> Jonckheere-Terpstra trend + Spearman.
    """
    from scipy.stats import friedmanchisquare, kruskal, spearmanr
    T = list(tools)
    disp = [tool_disp.get(t, t) for t in T]
    rep = {}

    # ---- Panel A: counts of poses in the crystal cluster (fully paired) ----
    cols = [comp[t].to_numpy(float) for t in T]
    chi2A, pA = friedmanchisquare(*cols)
    nA = len(cols[0]); kA = len(T)
    W = float(chi2A) / (nA * (kA - 1)) if nA else np.nan       # Kendall's W
    posthocA = []
    for i in range(kA):
        for j in range(i + 1, kA):
            rb, p, npair = _wilcoxon_rankbiserial(cols[i], cols[j])
            posthocA.append({"pair": f"{disp[i]} vs {disp[j]}", "rank_biserial": rb,
                             "p_raw": p, "n_pairs": npair})
    padj = _holm([h["p_raw"] for h in posthocA])
    for h, pa in zip(posthocA, padj):
        h["p_holm"] = float(pa)
    rep["A_pose_counts"] = {
        "test": "Friedman (paired, k=3)", "chi2": round(float(chi2A), 3),
        "df": kA - 1, "p": float(pA), "kendall_w": round(W, 3), "n": nA,
        "medians": {disp[i]: float(np.median(cols[i])) for i in range(kA)},
        "posthoc_wilcoxon_holm": posthocA}

    # ---- Panel A companion: "reaches" >=1 pose (paired binary) ----
    Rmat = np.asarray(reach_R, int)
    Q, pQ, dfQ = _cochran_q(Rmat)
    posthocQ = []
    for i in range(kA):
        for j in range(i + 1, kA):
            n10, n01, p = _mcnemar_exact(Rmat[:, i], Rmat[:, j])
            posthocQ.append({"pair": f"{disp[i]} vs {disp[j]}",
                             "disc_i_only": n10, "disc_j_only": n01, "p_raw": p})
    padjQ = _holm([h["p_raw"] for h in posthocQ])
    for h, pa in zip(posthocQ, padjQ):
        h["p_holm"] = float(pa)
    rep["A_reaches"] = {
        "test": "Cochran's Q (paired binary, k=3)", "Q": round(float(Q), 3),
        "df": dfQ, "p": float(pQ),
        "reach_rate": {disp[i]: round(float(Rmat[:, i].mean()), 3) for i in range(kA)},
        "posthoc_mcnemar_holm": posthocQ}

    # ---- Panel B: consensus co-occurrence vs independence null ----
    rep["B_consensus"] = {"test": "independence-null permutation", **
                          _consensus_perm_test(Rmat)}

    # ---- Panel C: internal tightness std (unpaired; informative missingness) ----
    gC = [np.asarray(internal[t], float) for t in T]
    if all(len(g) >= 1 for g in gC) and sum(len(g) for g in gC) > len(T):
        Hc, pC = kruskal(*gC)
        NC = sum(len(g) for g in gC)
        eps2 = (float(Hc) - kA + 1) / (NC - kA) if NC > kA else np.nan
        posthocC = []
        for (i, j, z, p) in _dunn(gC):
            posthocC.append({"pair": f"{disp[i]} vs {disp[j]}", "z": round(z, 3),
                             "p_raw": p})
        padjC = _holm([h["p_raw"] for h in posthocC])
        for h, pa in zip(posthocC, padjC):
            h["p_holm"] = float(pa)
        # sensitivity: Friedman on complexes present for all three tools
        common = set.intersection(*[set(internal_prot[t]) for t in T]) if all(
            internal_prot[t] for t in T) else set()
        sens = None
        if len(common) >= 3:
            common = sorted(common)
            fcols = [np.array([internal_prot[t][pr] for pr in common]) for t in T]
            chi2s, ps = friedmanchisquare(*fcols)
            sens = {"n_common": len(common), "chi2": round(float(chi2s), 3),
                    "p": float(ps), "note": "paired subset present for all 3 tools"}
        rep["C_tightness"] = {
            "test": "Kruskal-Wallis (independent, unbalanced)",
            "H": round(float(Hc), 3), "df": kA - 1, "p": float(pC),
            "epsilon_sq": round(float(eps2), 3),
            "n": {disp[i]: int(len(gC[i])) for i in range(kA)},
            "medians": {disp[i]: (round(float(np.median(gC[i])), 3)
                        if len(gC[i]) else None) for i in range(kA)},
            "posthoc_dunn_holm": posthocC,
            "note": "n differs per tool because a >=2-pose in-cluster count is "
                    "required for a std; that missingness is informative "
                    "(a tool is absent when it places <2 poses). A complex-random-"
                    "effect mixed model would be more rigorous (needs statsmodels).",
            "sensitivity_friedman_common_subset": sens}

    # ---- Panel D: crystal-cluster radius vs number of contributing tools ----
    gD = [np.asarray(rad_by_n[k], float) for k in (1, 2, 3)]
    gD = [g[~np.isnan(g)] for g in gD]
    if all(len(g) >= 1 for g in gD):
        JT, zD, pD = _jonckheere(gD)
        nt = np.array([p[0] for p in nt_radius_pairs], float)
        rr = np.array([p[1] for p in nt_radius_pairs], float)
        m = ~np.isnan(rr)
        rho, prho = spearmanr(nt[m], rr[m]) if m.sum() > 2 else (np.nan, np.nan)
        rep["D_radius_trend"] = {
            "test": "Jonckheere-Terpstra (ordered 1<2<3 tools)",
            "JT": round(JT, 1), "z": round(zD, 3), "p": float(pD),
            "spearman_rho": round(float(rho), 3), "spearman_p": float(prho),
            "n": {f"{k}_tools": int(len(gD[k - 1])) for k in (1, 2, 3)},
            "medians": {f"{k}_tools": round(float(np.median(gD[k - 1])), 3)
                        for k in (1, 2, 3)},
            "caveat": "radius = max(pose->center) grows mechanically with the "
                      "number of poses in the cluster, and more-consensus "
                      "clusters tend to hold more poses; treat the trend as "
                      "confounded unless controlled for pose count."}
    return rep


def _print_homogeneity_stats(rep):
    def line(s=""):
        print("    " + s)
    print("\n  ── crystal_cluster_homogeneity.png statistics ──")
    a = rep.get("A_pose_counts")
    if a:
        line(f"(A) pose counts | Friedman chi2={a['chi2']} df={a['df']} "
             f"p={_fmt_p(a['p'])} {_p_stars(a['p'])}  Kendall W={a['kendall_w']}")
        for h in a["posthoc_wilcoxon_holm"]:
            line(f"      {h['pair']:<28} rank-biserial r={h['rank_biserial']:+.3f} "
                 f"p_Holm={_fmt_p(h['p_holm'])} {_p_stars(h['p_holm'])}")
    q = rep.get("A_reaches")
    if q:
        line(f"(A) reaches>=1  | Cochran's Q={q['Q']} df={q['df']} "
             f"p={_fmt_p(q['p'])} {_p_stars(q['p'])}  rates={q['reach_rate']}")
        for h in q["posthoc_mcnemar_holm"]:
            line(f"      {h['pair']:<28} McNemar p_Holm={_fmt_p(h['p_holm'])} "
                 f"{_p_stars(h['p_holm'])} (disc {h['disc_i_only']}/{h['disc_j_only']})")
    b = rep.get("B_consensus")
    if b:
        line(f"(B) consensus   | perm vs independence: chi2-like={b['chi2_like_stat']} "
             f"p={_fmt_p(b['p_omnibus'])} {_p_stars(b['p_omnibus'])}")
        line(f"      all-3 tools observed={b['all_tools_obs']} vs "
             f"expected={b['all_tools_exp']} (p={_fmt_p(b['all_tools_p'])})  "
             f"obs={b['observed']} exp={b['expected']} [0,1,2,3 tools]")
    c = rep.get("C_tightness")
    if c:
        line(f"(C) tightness   | Kruskal-Wallis H={c['H']} df={c['df']} "
             f"p={_fmt_p(c['p'])} {_p_stars(c['p'])}  eps^2={c['epsilon_sq']} n={c['n']}")
        for h in c["posthoc_dunn_holm"]:
            line(f"      {h['pair']:<28} Dunn z={h['z']:+.2f} "
                 f"p_Holm={_fmt_p(h['p_holm'])} {_p_stars(h['p_holm'])}")
        s = c.get("sensitivity_friedman_common_subset")
        if s:
            line(f"      sensitivity Friedman (paired, n={s['n_common']}) "
                 f"chi2={s['chi2']} p={_fmt_p(s['p'])}")
    d = rep.get("D_radius_trend")
    if d:
        line(f"(D) radius trend| Jonckheere-Terpstra z={d['z']} p={_fmt_p(d['p'])} "
             f"{_p_stars(d['p'])}  Spearman rho={d['spearman_rho']} "
             f"(p={_fmt_p(d['spearman_p'])})  n={d['n']}")
        line(f"      CAVEAT: {d['caveat']}")


def _fig_crystal_cluster_homogeneity(ok, df_rank, out_dir, eq_variant=None,
                                     stats=False):
    """How homogeneous is the crystal-closest cluster the tools land in?

    Four panels, all restricted to each tool's top-15 ranked poses (the scope used
    by the sibling rank figures):
      (A) Pose contribution — how many of a tool's top-15 poses fall in the crystal
          cluster (distribution over all complexes; 'reaches' = fraction with >=1).
      (B) Consensus richness — fraction of crystal clusters populated by 1 / 2 / 3
          tools.
      (C) Internal tightness — per tool, the std of that tool's own in-cluster
          poses' distances to the crystal centroid (needs >=2 in-cluster poses):
          how far apart a single tool's poses sit.
      (D) Spatial spread vs consensus — the crystal cluster's geometric radius
          (max pose->center, from the full clustering) stratified by how many tools
          contribute, i.e. does agreement make the cluster tighter or wider.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    need = {"protein", "tool", "in_correct_cluster", "centroid_dist"}
    if df_rank.empty or not need <= set(df_rank.columns):
        return None
    d = df_rank.copy()
    d["_ic"] = d["in_correct_cluster"].map(_TRUTH_MAP)
    proteins = list(d["protein"].unique())
    n_total = len(proteins)
    if n_total == 0:
        return None
    T = list(TOOLS)
    cc = d[d["_ic"] == True]                                             # noqa: E712
    comp = (cc.groupby(["protein", "tool"]).size().unstack(fill_value=0)
            .reindex(index=proteins, columns=T, fill_value=0))
    reach = {t: float((comp[t] > 0).mean()) for t in T}
    n_tools_per = (comp > 0).sum(axis=1)
    ntool_counts = {k: int((n_tools_per == k).sum()) for k in (1, 2, 3)}
    reach_R = (comp[T] > 0).astype(int).to_numpy()          # n x k reach indicator
    internal = {t: [] for t in T}
    internal_prot = {t: {} for t in T}                       # keep protein identity
    for (_prot, t), g in cc.groupby(["protein", "tool"]):
        if len(g) >= 2:
            v = pd.to_numeric(g["centroid_dist"], errors="coerce").to_numpy()
            s = float(np.std(v))
            internal[t].append(s)
            internal_prot[t][_prot] = s
    radius = {}
    for r in ok:
        if not r.get("has_crystal"):
            continue
        cr = r.get("_crystal"); pk = r.get("_pockets")
        if cr is None or not pk:
            continue
        cp = min(pk, key=lambda p: _dist(p["center"], cr))
        radius[r["protein"]] = float(cp.get("radius", np.nan))
    rad_by_n = {k: [] for k in (1, 2, 3)}
    nt_radius_pairs = []                     # (n_tools, radius) per complex for Spearman
    for prot in proteins:
        k = int(n_tools_per.get(prot, 0))
        rv = radius.get(prot, np.nan)
        if k in rad_by_n and rv == rv:
            rad_by_n[k].append(rv)
        nt_radius_pairs.append((k, rv))

    stats_rep = None
    if stats:
        try:
            stats_rep = _homogeneity_stats(comp, reach_R, internal, internal_prot,
                                           rad_by_n, nt_radius_pairs, T, _TOOL_DISPLAY)
        except Exception as e:                                  # never break the figure
            warnings.warn(f"homogeneity stats failed: {e}")
            stats_rep = None

    def _short(t):                              # AutoDock Vina -> AV, DiffDock -> DD
        return "".join(w[0] for w in _TOOL_DISPLAY.get(t, t).split())

    def _pair_stars(posthoc, key="p_holm"):
        # "AV-DD *** | AV-EB *** | DD-EB **" from a post-hoc list
        segs = []
        for h, (i, j) in zip(posthoc, [(0, 1), (0, 2), (1, 2)]):
            segs.append(f"{_short(T[i])}-{_short(T[j])} {_p_stars(h[key])}")
        return "Holm post-hoc:  " + "  |  ".join(segs)

    rng = np.random.default_rng(42)
    tri_col = ["#C44E52", "#DD8452", "#55A868"]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 10.4))
    axA, axB, axC, axD = axes.ravel()

    bp = axA.boxplot([comp[t].to_numpy() for t in T], patch_artist=True, widths=0.6,
                     showmeans=True, medianprops=dict(color="black"),
                     meanprops=dict(marker="D", markerfacecolor="white",
                                    markeredgecolor="black", markersize=6))
    for patch, t in zip(bp["boxes"], T):
        patch.set_facecolor(TOOL_COLORS[t]); patch.set_alpha(0.75)
    for i, t in enumerate(T, start=1):
        axA.scatter(rng.normal(i, 0.05, len(comp)), comp[t], s=6,
                    color=TOOL_COLORS[t], alpha=0.20, zorder=1)
        axA.text(i, 15.6, f"reaches\n{reach[t]:.0%}", ha="center", va="bottom",
                 fontsize=8, color=TOOL_COLORS[t], fontweight="bold")
    axA.set_xticks([1, 2, 3]); axA.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in T])
    axA.set_ylabel("Poses a tool places in the crystal cluster\n(top-15 scope)")
    axA.set_ylim(-0.5, 18)
    titleA = "Pose contribution per tool"
    if stats_rep and stats_rep.get("A_pose_counts"):
        a = stats_rep["A_pose_counts"]
        titleA += (f"\nFriedman χ²={a['chi2']:.1f}, p={_fmt_p(a['p'])}, "
                   f"Kendall W={a['kendall_w']:.2f}")
        axA.text(2, 17.0, _pair_stars(a["posthoc_wilcoxon_holm"]), ha="center",
                 va="center", fontsize=7.6, family="monospace", color="0.25")
    axA.set_title(titleA)
    axA.grid(axis="y", alpha=0.25); axA.set_axisbelow(True)

    xs = [1, 2, 3]
    fr = [ntool_counts[k] / n_total for k in xs]
    axB.bar(xs, fr, color=tri_col, width=0.7, alpha=0.85)
    for x in xs:
        axB.text(x, fr[x - 1] + 0.01, f"{fr[x - 1]:.0%}\n(n={ntool_counts[x]})",
                 ha="center", va="bottom", fontsize=9)
    axB.set_xticks(xs); axB.set_xticklabels(["1 tool", "2 tools", "all 3 tools"])
    axB.set_ylim(0, (max(fr) if fr else 1) + 0.12)
    axB.set_ylabel(f"Fraction of complexes (n={n_total})")
    axB.set_xlabel("Number of tools contributing ≥1 pose")
    titleB = "Consensus richness of the crystal cluster"
    if stats_rep and stats_rep.get("B_consensus"):
        b = stats_rep["B_consensus"]
        titleB += f"\nvs independence null: permutation p={_fmt_p(b['p_omnibus'])}"
    axB.set_title(titleB)
    axB.grid(axis="y", alpha=0.25); axB.set_axisbelow(True)

    bp = axC.boxplot([internal[t] if internal[t] else [np.nan] for t in T],
                     patch_artist=True, widths=0.6, showfliers=False,
                     medianprops=dict(color="black"))
    for patch, t in zip(bp["boxes"], T):
        patch.set_facecolor(TOOL_COLORS[t]); patch.set_alpha(0.75)
    for i, t in enumerate(T, start=1):
        yy = internal[t]
        if yy:
            axC.scatter(rng.normal(i, 0.05, len(yy)), yy, s=7,
                        color=TOOL_COLORS[t], alpha=0.25, zorder=1)
        axC.text(i, -0.08, f"n={len(yy)}", ha="center", va="top", fontsize=8,
                 color=TOOL_COLORS[t])
    axC.set_xticks([1, 2, 3]); axC.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in T])
    axC.set_ylabel("Std. of a tool's pose distances to the\ncrystal centroid, within the cluster (Å)")
    titleC = "Internal tightness of each tool's own poses"
    _c_top = None
    if stats_rep and stats_rep.get("C_tightness"):
        c = stats_rep["C_tightness"]
        titleC += (f"\nKruskal–Wallis H={c['H']:.1f}, p={_fmt_p(c['p'])}, "
                   f"ε²={c['epsilon_sq']:.2f}")
        _all = [v for t in T for v in internal[t]]
        _c_top = (max(_all) * 1.16 if _all else None)
        if _c_top:
            axC.text(2, max(_all) * 1.08, _pair_stars(c["posthoc_dunn_holm"]),
                     ha="center", va="center", fontsize=7.6, family="monospace",
                     color="0.25")
    axC.set_title(titleC)
    axC.set_ylim(-0.15, _c_top)
    axC.grid(axis="y", alpha=0.25); axC.set_axisbelow(True)

    dataD = [rad_by_n[k] if rad_by_n[k] else [np.nan] for k in xs]
    bp = axD.boxplot(dataD, patch_artist=True, widths=0.6, showmeans=True,
                     medianprops=dict(color="black"),
                     meanprops=dict(marker="D", markerfacecolor="white",
                                    markeredgecolor="black", markersize=6))
    for patch, c in zip(bp["boxes"], tri_col):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    for i, _k in enumerate(xs):
        vals = [v for v in dataD[i] if v == v]
        if vals:
            med = float(np.median(vals))
            axD.text(i + 1, med + 0.12, f"{med:.2f} Å", ha="center",
                     va="bottom", fontsize=9)
    axD.set_xticks(xs); axD.set_xticklabels(["1 tool", "2 tools", "all 3 tools"])
    axD.set_ylabel("Crystal-cluster radius (max pose->center, Å)")
    axD.set_xlabel("Number of tools contributing to the crystal cluster")
    titleD = "Spatial spread vs consensus richness"
    if stats_rep and stats_rep.get("D_radius_trend"):
        dd = stats_rep["D_radius_trend"]
        titleD += (f"\nJonckheere–Terpstra trend p={_fmt_p(dd['p'])}; "
                   f"ρ={dd['spearman_rho']:+.2f}")
    axD.set_title(titleD)
    axD.grid(axis="y", alpha=0.25); axD.set_axisbelow(True)

    _label_panels(axes)
    suptitle = ("How homogeneous is the crystal-closest cluster?  "
                "Composition, internal tightness and spatial spread\n"
                f"n={n_total} complexes; top-15 poses per tool; "
                f"{_equibind_rank_note(eq_variant)}")
    if stats_rep:
        suptitle += "   ·   stars: * p<.05  ** p<.01  *** p<.001 (Holm)"
    fig.suptitle(suptitle, fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out_dir / "crystal_cluster_homogeneity.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    if stats_rep:
        (out_dir / "crystal_cluster_homogeneity_stats.json").write_text(
            json.dumps(stats_rep, indent=2, default=str))
        _print_homogeneity_stats(stats_rep)
        print(f"  Stats JSON: {out_dir / 'crystal_cluster_homogeneity_stats.json'}")
    return p


def _fig_placement(df, mode_thr, out_dir):
    """Centroid SITES vs placement-aware MODES: how the two geometry axes compare."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if "n_modes" not in df:
        return None
    d = df[pd.to_numeric(df["n_modes"], errors="coerce").notna()].copy()
    if d.empty:
        return None
    ns = pd.to_numeric(d["n_clusters"], errors="coerce")
    nm = pd.to_numeric(d["n_modes"], errors="coerce")
    fig, ax = plt.subplots(2, 2, figsize=(13, 10))

    # (A) #sites vs #modes per complex — modes >= sites means a site holds >1 mode
    hi = int(np.nanmax([ns.max(), nm.max(), 1])) + 1
    ax[0, 0].scatter(ns + np.random.RandomState(0).uniform(-0.15, 0.15, len(ns)),
                     nm + np.random.RandomState(1).uniform(-0.15, 0.15, len(nm)),
                     s=18, alpha=0.5, color="#4C72B0", edgecolors="none")
    ax[0, 0].plot([0, hi], [0, hi], "k--", lw=0.8, label="modes = sites")
    ax[0, 0].set_xlim(0, hi); ax[0, 0].set_ylim(0, hi)
    ax[0, 0].set_xlabel("Number of centroid sites (where)")
    ax[0, 0].set_ylabel("Number of placement-aware modes (how)")
    ax[0, 0].set_title("Sites vs binding modes per complex\n(above line ⇒ a site splits into modes)")
    ax[0, 0].legend(fontsize=8); ax[0, 0].grid(alpha=0.25); ax[0, 0].set_axisbelow(True)

    # (B) modes inside the single biggest site
    mits = pd.to_numeric(d.get("modes_in_top_site"), errors="coerce").dropna()
    if len(mits):
        vc = mits.clip(upper=5).value_counts().sort_index()
        ax[0, 1].bar(vc.index, vc.values / vc.sum(), color="#55A868", width=0.8)
        ax[0, 1].set_xticks(sorted(vc.index))
        ax[0, 1].set_xlabel("Distinct binding modes inside the biggest site (≥5 binned)")
        ax[0, 1].set_ylabel("Fraction of complexes")
        ax[0, 1].set_title(f"Does one spatial site hold several modes?\n"
                           f">1 mode in {(mits > 1).mean():.0%} of complexes")
        ax[0, 1].grid(axis="y", alpha=0.25); ax[0, 1].set_axisbelow(True)

    # (C) site-vs-mode agreement (ARI) — low ⇒ the axes disagree
    ari = pd.to_numeric(d.get("ari_site_vs_mode"), errors="coerce").dropna()
    if len(ari):
        ax[1, 0].hist(ari, bins=np.linspace(0, 1, 21), color="#8172B3", alpha=0.85)
        ax[1, 0].axvline(float(ari.median()), color="k", ls="--", lw=1,
                         label=f"median {ari.median():.2f}")
        ax[1, 0].set_xlabel("ARI(centroid sites, placement modes)")
        ax[1, 0].set_ylabel("Number of complexes")
        ax[1, 0].set_title("Do co-located poses share a binding mode?\n(1 = identical partitions)")
        ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=0.25); ax[1, 0].set_axisbelow(True)

    # (D) mode compactness vs separation (are modes tight and well-separated?)
    comp = pd.to_numeric(d.get("mode_compactness"), errors="coerce")
    sep = pd.to_numeric(d.get("mode_separation"), errors="coerce")
    m = comp.notna() & sep.notna()
    if m.any():
        ax[1, 1].scatter(comp[m], sep[m], s=18, alpha=0.5, color="#C44E52", edgecolors="none")
        ax[1, 1].axhline(mode_thr, color="k", ls="--", lw=0.8,
                         label=f"mode threshold {mode_thr:g} Å")
        ax[1, 1].set_xlabel("Within-mode RMSD spread (Å, lower = tighter)")
        ax[1, 1].set_ylabel("Nearest inter-mode RMSD (Å, higher = better separated)")
        ax[1, 1].set_title("Binding-mode compactness vs separation")
        ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=0.25); ax[1, 1].set_axisbelow(True)

    _label_panels(ax)
    fig.suptitle("Placement-aware binding modes (in-place RMSD) vs centroid sites  "
                 f"(n={len(d)} complexes, mode-thr={mode_thr:g} Å)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "placement_vs_centroid.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def _fig_cluster_quality(df, ablation, ranking_rho, out_dir, match_thr,
                         stability_boot):
    """Visualise the cluster-VALIDITY numbers that otherwise live only in
    summary.json / per_complex_summary.csv / ranking_ablation.csv / console.

    Six panels on the primary (centroid) partition, one figure:
      (A-C) internal indices — silhouette, Calinski-Harabasz, Davies-Bouldin
            (all k>=2 only; k=1 complexes contribute NaN and are dropped);
      (D)   the k=1-safe pair compactness vs separation (tight & far-apart=good);
      (E)   bootstrap stability distribution (>0.75 stable, <0.5 dissolved);
      (F)   the precision@1 ranking ablation — does consensus ranking recover the
            true crystal site more often than legacy size, and how close to the
            oracle ceiling.
    Internal indices say "well-separated geometry", stability says "reproducible",
    precision@1 says "actually the right pocket" — three independent axes.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def _col(name):
        return (pd.to_numeric(df.get(name), errors="coerce").dropna()
                if name in df else pd.Series(dtype=float))

    def _hist(a, series, title, xlabel, color, ref=None):
        if len(series):
            a.hist(series, bins=20, color=color, alpha=0.85)
            a.axvline(float(series.median()), color="k", ls="--", lw=1,
                      label=f"median {series.median():.2f}")
            for x, lab, c in (ref or []):
                a.axvline(x, color=c, ls=":", lw=1.4, label=lab)
            a.legend(fontsize=8)
        else:
            a.text(0.5, 0.5, "no data (all k=1)", ha="center", va="center",
                   transform=a.transAxes, color="0.5")
        a.set_title(title)
        a.set_xlabel(xlabel)
        a.set_ylabel("Number of complexes")
        a.grid(alpha=0.25); a.set_axisbelow(True)

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.5))

    # (A) Silhouette
    _hist(ax[0, 0], _col("silhouette"),
          "Silhouette coefficient (k≥2 only)\nhigher = tighter, better-separated sites",
          "Silhouette coefficient (−1 to 1)", "#4C72B0")

    # (B) Calinski-Harabasz — long right tail; clip display at the 98th pctile
    ch = _col("ch_score")
    ch_disp = ch.clip(upper=float(ch.quantile(0.98))) if len(ch) else ch
    _hist(ax[0, 1], ch_disp,
          "Calinski-Harabasz index (k≥2 only)\nhigher = better (between/within variance)",
          "Calinski-Harabasz index (98th-pctile clipped)", "#55A868")

    # (C) Davies-Bouldin
    _hist(ax[0, 2], _col("db_score"),
          "Davies-Bouldin index (k≥2 only)\nlower = better",
          "Davies-Bouldin index", "#C44E52")

    # (D) compactness vs separation (both defined at k=1)
    comp = pd.to_numeric(df.get("compactness"), errors="coerce")
    sep = pd.to_numeric(df.get("separation"), errors="coerce")
    m = comp.notna() & sep.notna()
    a = ax[1, 0]
    if m.any():
        a.scatter(comp[m], sep[m], s=18, alpha=0.5, color="#8172B3", edgecolors="none")
        a.axhline(float(match_thr), color="k", ls="--", lw=0.8,
                  label=f"match threshold {match_thr:g} Å")
        a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "no data", ha="center", va="center",
               transform=a.transAxes, color="0.5")
    a.set_title("Site compactness vs separation\ntight (low x) and far-apart (high y) = good")
    a.set_xlabel("Within-site spread (Å, median distance-to-center; lower = tighter)")
    a.set_ylabel("Nearest inter-site distance (Å; higher = better)")
    a.grid(alpha=0.25); a.set_axisbelow(True)

    # (E) bootstrap stability
    _hist(ax[1, 1], _col("boot_stability"),
          f"Bootstrap stability ({stability_boot} resamples)\nreproducibility of the partition",
          "Mean best-cluster Jaccard overlap", "#CCB974",
          ref=[(0.75, "stable ≥0.75", "#2ca02c"),
               (0.5, "dissolved <0.5", "#d62728")])

    # (F) precision@1 ranking ablation
    a = ax[1, 2]
    rule_lab = {"size": "size (legacy count)", "ntools": "consensus (n tools)",
                "tight": "tightest spread", "confidence": "confidence-weighted",
                "medoid": "consensus + medoid"}
    bars = [(rule_lab[r["rule"]], float(r["precision_at_1"]), r["rule"])
            for r in (ablation or [])
            if r.get("rule") in rule_lab and r.get("precision_at_1") is not None]
    oracle = next((float(r["precision_at_1"]) for r in (ablation or [])
                   if r.get("rule") == "oracle" and r.get("precision_at_1") is not None), None)
    if bars:
        bars.sort(key=lambda x: x[1])
        labels = [b[0] for b in bars]
        vals = [b[1] for b in bars]
        cols = ["#DD8452" if b[2] == "ntools" else "#B0A08F" for b in bars]
        y = np.arange(len(bars))
        a.barh(y, vals, color=cols)
        a.set_yticks(y); a.set_yticklabels(labels, fontsize=9)
        for yi, v in zip(y, vals):
            a.text(min(v + 0.015, 1.0), yi, f"{v:.0%}", va="center", fontsize=8)
        if oracle is not None:
            a.axvline(oracle, color="k", ls="--", lw=1.3,
                      label=f"oracle ceiling {oracle:.0%}")
            a.legend(fontsize=8, loc="upper left", framealpha=0.92)
        a.set_xlim(0, 1.08)
    else:
        a.text(0.5, 0.5, "no data", ha="center", va="center",
               transform=a.transAxes, color="0.5")
    ttl = f"Precision@1 — rank-1 pocket within {match_thr:g} Å of crystal (default = consensus)"
    if ranking_rho is not None:
        ttl += f"\nranking enrichment Spearman ρ={ranking_rho:.2f}"
    st = _ablation_paired_stats(df)
    if st:
        # omnibus in the title; the full pairwise family (incl. consensus vs size)
        # is the star matrix below, and the console/summary.json carry the counts.
        ttl += (f"\nCochran Q={st['cochran_q']:.1f}, p={_fmt_p(st['p_omnibus'])} "
                f"(paired, k={st['k']}, n={st['n_omnibus']})")
    a.set_title(ttl)
    # complete Holm-adjusted pairwise McNemar family as a lower-triangular star matrix
    if bars and st:
        mtxt = _ablation_pairwise_matrix_text(st)
        if mtxt:
            a.text(0.985, 0.03, mtxt, transform=a.transAxes, ha="right", va="bottom",
                   family="monospace", fontsize=6.0, color="0.15", zorder=5,
                   linespacing=1.25,
                   bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                             edgecolor="0.75", alpha=0.92))
    a.set_xlabel("Precision@1 (fraction of complexes)")
    a.grid(axis="x", alpha=0.25); a.set_axisbelow(True)

    _label_panels(ax)
    fig.suptitle("Cluster-quality assessment — internal validity, stability, and "
                 f"ranking accuracy  (n={len(df)} complexes)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "cluster_quality_metrics.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Driver
# ════════════════════════════════════════════════════════════════════════

def _fig_descriptor_quality(df_complex, features_csv, out_dir,
                            per_pose_csv=None, dd_variant="diffdock",
                            eq_variant="equibind_unguided_smina"):
    """Which ligand types dock well? Pose accuracy (oracle RMSD to crystal) vs the
    ligand's Lipinski Ro5 compliance, physicochemical descriptors, and PCA space.

    Every success measure is reported for TWO categories: ``near-native`` (a pose
    RMSD <= 2 A, any validity) and ``PB-valid & <= 2 A`` (the near-native pose is
    also PoseBusters-valid) — the gap is the fraction of near-native poses lost to
    physical invalidity. Both oracles are recomputed from the raw per-pose CSV (over
    autodock + the chosen DiffDock/EquiBind variants) so the comparison holds
    regardless of any --pb-valid-only clustering filter. One PNG per measure; returns
    the list of written paths.
    """
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
    d = (df_complex[df_complex["has_crystal"]].copy()
         if "has_crystal" in df_complex else df_complex.copy())
    if d.empty:
        return None
    NEAR = 2.0
    tool_method = {"autodock": "autodock", "diffdock": dd_variant, "equibind": eq_variant}

    # ── per-complex oracle RMSD, near-native (any validity) vs PB-valid, recomputed
    #    from the raw per-pose CSV over the analysed methods (so it is independent of
    #    any --pb-valid-only clustering filter and can show BOTH categories) ──
    oracle_all = oracle_valid = None
    oracle_all_tool: dict = {}
    if per_pose_csv and Path(per_pose_csv).exists():
        pp = pd.read_csv(per_pose_csv, low_memory=False)
        pp = pp[pp["method"].astype(str).isin(set(tool_method.values()))].copy()
        pp["rmsd"] = pd.to_numeric(pp["rmsd"], errors="coerce")
        pp = pp.dropna(subset=["rmsd"])
        pbv_col = pp.get("pb_valid")
        pbv = (pbv_col.astype(str).str.lower().isin(("true", "1", "1.0"))
               if pbv_col is not None else pd.Series(False, index=pp.index))
        oracle_all = pp.groupby("protein")["rmsd"].min()
        oracle_valid = pp[pbv].groupby("protein")["rmsd"].min()
        for t, m in tool_method.items():
            oracle_all_tool[t] = pp[pp["method"].astype(str) == m].groupby("protein")["rmsd"].min()

    if oracle_all is None:
        rcols = [f"{t}_oracle_rmsd" for t in TOOLS if f"{t}_oracle_rmsd" in d.columns]
        if not rcols:
            return None
        for c in rcols:
            d[c] = pd.to_numeric(d[c], errors="coerce")
        d = d.set_index("protein")
        d["oracle_all"] = d[rcols].min(axis=1)
        d["oracle_valid"] = d["oracle_all"]          # no unfiltered baseline available
    else:
        d = d.set_index("protein")
        d["oracle_all"] = d.index.map(oracle_all)
        d["oracle_valid"] = d.index.map(oracle_valid)
        for t in TOOLS:
            d[f"{t}_oracle_all"] = d.index.map(oracle_all_tool.get(t, pd.Series(dtype=float)))

    d = d.join(desc, how="inner")
    d["oracle_all"] = pd.to_numeric(d["oracle_all"], errors="coerce")
    d["oracle_valid"] = pd.to_numeric(d["oracle_valid"], errors="coerce")
    d = d[np.isfinite(d["oracle_all"])]
    if len(d) < 12:
        return None
    d["best_oracle_rmsd"] = d["oracle_all"]          # correlations / PCA use geometry
    d["success_near"] = d["oracle_all"] <= NEAR
    d["success_valid"] = d["oracle_valid"] <= NEAR   # NaN oracle_valid -> False
    d, loadings, evr, _ = add_pca(d)
    bins, labs = [-0.1, 2, 5, 8, 11, 100], ["0–2", "3–5", "6–8", "9–11", "12+"]
    d["_rb"] = pd.cut(pd.to_numeric(d["rot_bonds"], errors="coerce"), bins=bins, labels=labs)
    d.reset_index().to_csv(out_dir / "descriptor_vs_quality.csv", index=False)

    NEAR_C, VALID_C = "#4C72B0", "#2ca02c"           # near-native (blue), PB-valid&2A (green)
    n = len(d)
    sub = (f"n={n} complexes with crystal · oracle = closest pose any tool produced · "
           "near-native = RMSD ≤ 2 Å; PB-valid & ≤ 2 Å = also PoseBusters-valid")
    saved = []

    def _save(fig, fname):
        fig.suptitle(sub, fontsize=8, y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        p = out_dir / fname
        fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
        saved.append(p)

    # (A) accuracy vs Lipinski Ro5 — boxplot of near-native RMSD + dual success rates
    fig, ax = plt.subplots(figsize=(6.8, 5.2))
    grp = [("Ro5 pass\n(≤1 viol.)", d[d.ro5_pass]),
           ("Ro5 fail\n(≥2 viol.)", d[~d.ro5_pass])]
    ax.boxplot([g["oracle_all"].dropna().clip(upper=20) for _, g in grp],
               labels=[nm for nm, _ in grp], showfliers=False)
    for i, (_nm, g) in enumerate(grp, 1):
        near = 100 * g["success_near"].mean() if len(g) else np.nan
        val = 100 * g["success_valid"].mean() if len(g) else np.nan
        ax.text(i, 0.4, f"≤2Å: {near:.0f}%\nPB-valid&≤2Å: {val:.0f}%\nn={len(g)}",
                ha="center", va="bottom", fontsize=8)
    ax.axhline(2, color="green", ls="--", lw=1)
    ax.set_title("Pose accuracy vs Lipinski Rule-of-Five compliance")
    ax.set_ylabel("Best oracle RMSD to crystal (Å)")
    ax.grid(alpha=0.2); ax.set_axisbelow(True)
    _save(fig, "descriptor_accuracy_vs_ro5.png")

    # (B) docking success vs flexibility — near-native + PB-valid&2A lines
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    gn = d.groupby("_rb")["success_near"].agg(["mean", "size"])
    gv = d.groupby("_rb")["success_valid"].agg(["mean", "size"])
    xs = range(len(gn))
    ax.plot(xs, gn["mean"] * 100, "-o", color=NEAR_C, lw=2, markersize=7,
            markeredgecolor="black", markeredgewidth=0.6, label="near-native (≤ 2 Å)")
    ax.plot(xs, gv["mean"] * 100, "--s", color=VALID_C, lw=2, markersize=6,
            label="PB-valid & ≤ 2 Å")
    for i, (m, nn) in enumerate(zip(gn["mean"], gn["size"])):
        if nn:
            ax.text(i, m * 100 + 3.5, f"{m*100:.0f}%\nn={int(nn)}", ha="center",
                    va="bottom", fontsize=8, color=NEAR_C)
    for i, m in enumerate(gv["mean"]):
        if m == m:
            ax.text(i, m * 100 - 5, f"{m*100:.0f}%", ha="center", va="top",
                    fontsize=8, color=VALID_C, fontweight="bold")
    ax.set_xticks(list(xs)); ax.set_xticklabels(gn.index)
    ax.set_title("Docking success vs ligand flexibility")
    ax.set_xlabel("Rotatable bonds")
    ax.set_ylabel("Complexes with a pose (percent)")
    ax.set_ylim(0, 105); ax.grid(alpha=0.3); ax.set_axisbelow(True); ax.legend(fontsize=9)
    _save(fig, "descriptor_success_vs_flexibility.png")

    # (C) Spearman correlation of each property (+ PCs) with pose error (geometry)
    feats = [f for f in PCA_FEATURES if f in d.columns] + [c for c in ("PC1", "PC2", "PC3") if c in d.columns]
    rhos = []
    for f in feats:
        v = pd.to_numeric(d[f], errors="coerce")
        m = v.notna() & d["best_oracle_rmsd"].notna()
        if m.sum() > 10:
            rho, _ = spearmanr(v[m], d["best_oracle_rmsd"][m])
            rhos.append((f, rho))
    if rhos:
        rhos.sort(key=lambda x: x[1])
        vals = [r for _, r in rhos]
        fig, ax = plt.subplots(figsize=(6.8, 5.4))
        ax.barh(range(len(vals)), vals, color=["#C44E52" if r > 0 else "#4C72B0" for r in vals])
        ax.set_yticks(range(len(vals))); ax.set_yticklabels([nice(f) for f, _ in rhos], fontsize=8)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_title("Which properties track pose error?\n(Spearman ρ vs RMSD; red = worse, blue = better)")
        ax.set_xlabel("Spearman correlation with best oracle RMSD")
        ax.grid(alpha=0.2, axis="x"); ax.set_axisbelow(True)
        _save(fig, "descriptor_property_correlations.png")

    # (D) PCA chemical space coloured by accuracy (geometry)
    if "PC1" in d and "PC2" in d:
        fig, ax = plt.subplots(figsize=(6.8, 5.4))
        sc = ax.scatter(d["PC1"], d["PC2"], c=d["best_oracle_rmsd"].clip(upper=10),
                        cmap="viridis_r", s=24, edgecolors="k", linewidths=0.2)
        plt.colorbar(sc, ax=ax, label="Best oracle RMSD (Å, clipped 10)")
        e1 = evr[0] * 100 if len(evr) > 0 else 0
        e2 = evr[1] * 100 if len(evr) > 1 else 0
        ax.set_xlabel(f"Principal Component 1 ({e1:.0f}% of variance)")
        ax.set_ylabel(f"Principal Component 2 ({e2:.0f}% of variance)")
        ax.set_title("Ligand chemical space (PCA) coloured by pose accuracy")
        ax.grid(alpha=0.2); ax.set_axisbelow(True)
        _save(fig, "descriptor_pca_accuracy.png")

    # (E) success across PCA space (tertile heat-map) — near-native shaded, each cell
    #     annotated 'near-native% / PB-valid&≤2Å%'
    if "PC1" in d and "PC2" in d and d["PC1"].notna().sum() > 20:
        d["_p1"] = pd.qcut(d["PC1"], 3, labels=["low", "mid", "high"], duplicates="drop")
        d["_p2"] = pd.qcut(d["PC2"], 3, labels=["low", "mid", "high"], duplicates="drop")
        piv_n = d.pivot_table(index="_p2", columns="_p1", values="success_near", aggfunc="mean") * 100
        piv_v = d.pivot_table(index="_p2", columns="_p1", values="success_valid", aggfunc="mean") * 100
        fig, ax = plt.subplots(figsize=(6.6, 5.4))
        im = ax.imshow(piv_n.values, cmap="RdYlGn", vmin=0, vmax=100, origin="lower", aspect="auto")
        ax.set_xticks(range(piv_n.shape[1])); ax.set_xticklabels(piv_n.columns)
        ax.set_yticks(range(piv_n.shape[0])); ax.set_yticklabels(piv_n.index)
        for i in range(piv_n.shape[0]):
            for j in range(piv_n.shape[1]):
                a = piv_n.values[i, j]
                b = piv_v.values[i, j] if (i < piv_v.shape[0] and j < piv_v.shape[1]) else np.nan
                if a == a:
                    ax.text(j, i, f"{a:.0f}%\n{b:.0f}%" if b == b else f"{a:.0f}%",
                            ha="center", va="center", fontsize=9)
        plt.colorbar(im, ax=ax, label="near-native ≤ 2 Å (percent)")
        ax.set_xlabel("Principal Component 1 tertile")
        ax.set_ylabel("Principal Component 2 tertile")
        ax.set_title("Docking success across ligand chemical space\n"
                     "(cell = near-native% / PB-valid & ≤2Å%)")
        _save(fig, "descriptor_success_pca.png")

    # (F) per-tool success vs flexibility (near-native) + PB-valid&≤2Å (any tool)
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for t in TOOLS:
        col = f"{t}_oracle_all" if f"{t}_oracle_all" in d else f"{t}_oracle_rmsd"
        if col not in d:
            continue
        rate = (pd.to_numeric(d[col], errors="coerce") <= NEAR).groupby(d["_rb"]).mean() * 100
        ax.plot(range(len(rate)), rate.values, "-o", color=TOOL_COLORS[t],
                label=f"{_TOOL_DISPLAY.get(t, t)} (near-native)")
    valrate = d.groupby("_rb")["success_valid"].mean() * 100
    ax.plot(range(len(valrate)), valrate.values, "--", color="black", lw=2,
            marker="s", markersize=5, label="any tool: PB-valid & ≤ 2 Å")
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs)
    ax.set_title("Per-tool success vs ligand flexibility")
    ax.set_xlabel("Rotatable bonds"); ax.set_ylabel("Pose ≤ 2 Å (percent)")
    ax.set_ylim(0, 105); ax.legend(fontsize=8); ax.grid(alpha=0.25); ax.set_axisbelow(True)
    _save(fig, "descriptor_per_tool_success_vs_flexibility.png")

    return saved

# ════════════════════════════════════════════════════════════════════════
# Analysis cache — the per-complex analysis (clustering / placement / bootstrap)
# is the expensive part; the CSVs + figures are cheap to regenerate from it. We
# pickle the analysed complexes keyed by a fingerprint of the inputs, so a re-run
# with unchanged inputs skips straight to (re)writing CSVs + figures.
# ════════════════════════════════════════════════════════════════════════
_CACHE_SCHEMA = 5


def _file_fp(path) -> Optional[dict]:
    """Content fingerprint (size + sha256) of a file, or None if it is absent."""
    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return {"size": p.stat().st_size, "sha256": h.hexdigest()}


def _analysis_signature(args, eq_variant: str, dd_variant: str) -> dict:
    """Fingerprint of everything the per-complex analysis (``ok``) depends on.

    Deliberately excludes figure-only inputs (``--features-csv``) and plotting code,
    so changing a figure and re-running still reuses the cached analysis. Pocket /
    crystal directories are keyed by path only (their contents are cheap to change
    but expensive to fingerprint) — pass ``--force`` if you regenerate those.
    """
    return {
        "schema": _CACHE_SCHEMA,
        "per_pose_csv": _file_fp(args.per_pose_csv),
        "ids_file": _file_fp(args.ids_file) if args.ids_file else None,
        "eq_variant": eq_variant, "dd_variant": dd_variant,
        "pb_valid_only": bool(args.pb_valid_only),
        "limit": int(args.limit),
        "outlier_dist": float(args.outlier_dist),
        "distance_mode": args.distance_mode,
        "match_thr": float(args.match_thr),
        "mode_rmsd_thr": float(args.mode_rmsd_thr),
        "rmsd_pose_cap": int(args.rmsd_pose_cap),
        "site_cluster": args.site_cluster,
        "pocket_radius": float(args.pocket_radius),
        "rank_by": args.rank_by,
        "stability_boot": int(args.stability_boot),
        "top_n_pockets": int(args.top_n_pockets),
        "benchmark_dir": str(args.benchmark_dir),
        "fpocket_dir": str(args.fpocket_dir),
        "p2rank_dir": str(args.p2rank_dir),
    }


def _load_analysis_cache(cache_path: Path, sig: dict):
    """Return the cached ``ok`` list if the cache matches ``sig``, else None."""
    if not cache_path.exists():
        return None
    try:
        with open(cache_path, "rb") as f:
            blob = pickle.load(f)
    except Exception as e:
        print(f"  (analysis cache unreadable — recomputing: {e})")
        return None
    if not isinstance(blob, dict) or blob.get("sig") != sig:
        return None
    ok = blob.get("ok")
    return ok if ok else None


def _save_analysis_cache(cache_path: Path, sig: dict, ok: list) -> None:
    tmp = cache_path.with_suffix(cache_path.suffix + ".tmp")
    try:
        with open(tmp, "wb") as f:
            pickle.dump({"schema": _CACHE_SCHEMA, "sig": sig, "ok": ok},
                        f, protocol=pickle.HIGHEST_PROTOCOL)
        tmp.replace(cache_path)
    except Exception as e:
        print(f"  (warning: could not write analysis cache: {e})")
        try:
            tmp.unlink()
        except OSError:
            pass


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")
    ap.add_argument("--benchmark-dir", default="Data/PoseBuster Benchmark Set")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--p2rank-dir", default="pocket_results/p2rank_results")
    ap.add_argument("--ids-file", default=None)
    # ── EquiBind tool slot: choose the variant by component, or by full name ──
    ap.add_argument("--equibind-variant", default=None,
                    help="Full equibind method string (e.g. "
                         "'equibind_fpocket_smina_clampON'). Overrides the "
                         "--equibind-pocket/-refine/-clamp components below. "
                         "Default: compose from the components.")
    ap.add_argument("--equibind-pocket", choices=("unguided", "fpocket", "p2rank"),
                    default="unguided",
                    help="EquiBind pocket source (default: unguided).")
    ap.add_argument("--equibind-refine", choices=_REFINERS, default="smina",
                    help="EquiBind pose refiner: 'raw' (unrefined) | 'smina' | "
                         "'gnina' (default: smina). The gnina method names omit "
                         "the token in the CSV; this resolves it for you.")
    ap.add_argument("--equibind-clamp", choices=("clampON", "clampOFF"),
                    default="clampOFF",
                    help="EquiBind clamp variant (default: clampOFF).")
    # ── DiffDock tool slot: choose the refiner, or a full name ────────────────
    ap.add_argument("--diffdock-variant", default=None,
                    help="Full DiffDock method string ('diffdock' | 'diffdock_smina' "
                         "| 'diffdock_gnina') or a bare refiner keyword. Overrides "
                         "--diffdock-refine. Default: compose from --diffdock-refine.")
    ap.add_argument("--diffdock-refine", choices=_REFINERS, default="raw",
                    help="DiffDock pose refiner: 'raw' -> diffdock | 'smina' -> "
                         "diffdock_smina | 'gnina' -> diffdock_gnina (default: raw).")
    ap.add_argument("--distance-mode", choices=("centroid", "hybrid", "both"),
                    default="both",
                    help="'centroid' = site clustering only; 'hybrid'/'both' also run "
                         "the placement-aware (in-place RMSD) binding-MODE analysis.")
    ap.add_argument("--mode-rmsd-thr", type=float, default=2.0,
                    help="RMSD (Å) threshold for the placement-aware mode cut "
                         "(2 Å = the near-native / same-pose cutoff).")
    ap.add_argument("--hybrid-weight", type=float, default=0.5,
                    help="(Deprecated; the old centroid+shape blend was replaced by "
                         "the placement-aware mode analysis. Unused.)")
    ap.add_argument("--match-thr", type=float, default=4.0)
    ap.add_argument("--site-cluster", choices=("threshold", "gap", "silhouette"),
                    default="threshold",
                    help="Primary site clustering. 'threshold' (default): complete-"
                         "linkage cut at --pocket-radius, k discovered incl. k=1. "
                         "'gap': gap statistic (admits k=1) then KMedoids. "
                         "'silhouette': legacy silhouette-max KMedoids (never k=1).")
    ap.add_argument("--pocket-radius", type=float, default=8.0,
                    help="Distance-threshold (Å) for the 'threshold' site cut "
                         "(≈ one pocket diameter; separate from --match-thr).")
    ap.add_argument("--rank-by", choices=("consensus", "size", "tight", "confidence"),
                    default="consensus",
                    help="How the rank-1 ('top1') pocket is chosen. 'consensus' "
                         "(default): distinct-tool agreement first; 'size': legacy "
                         "raw pose count.")
    ap.add_argument("--stability-boot", type=int, default=25,
                    help="Bootstrap resamples for cluster-stability Jaccard (0=off).")
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
    ap.add_argument("--stats", action="store_true",
                    help="Compute + annotate + export statistical tests for the "
                         "crystal_cluster_homogeneity figure (Friedman/Cochran's Q "
                         "for A, permutation for B, Kruskal-Wallis/Dunn for C, "
                         "Jonckheere-Terpstra for D); writes "
                         "crystal_cluster_homogeneity_stats.json.")
    ap.add_argument("--force", action="store_true",
                    help="Recompute the per-complex analysis even if a matching "
                         "cache (analysis_cache.pkl) exists for the current inputs.")
    args = ap.parse_args(argv)

    csv = Path(args.per_pose_csv)
    if not csv.exists():
        print(f"per-pose CSV not found: {csv}")
        return 1
    ids = _load_ids(Path(args.ids_file)) if args.ids_file else None
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve which docking-mode variant fills each tool slot (exits with a helpful
    # list of what's available if the selection can't be matched in the CSV).
    cat = _variant_catalog(csv)
    eq_variant = resolve_equibind_variant(
        cat, args.equibind_variant, args.equibind_pocket,
        args.equibind_refine, args.equibind_clamp)
    dd_variant = resolve_diffdock_variant(cat, args.diffdock_variant, args.diffdock_refine)

    # ── reuse the cached per-complex analysis when the inputs are unchanged ──
    sig = _analysis_signature(args, eq_variant, dd_variant)
    cache_path = out_dir / "analysis_cache.pkl"
    ok = None if args.force else _load_analysis_cache(cache_path, sig)
    if ok is not None:
        print(f"Reusing cached analysis: {len(ok)} complexes (inputs unchanged; "
              f"skipping re-analysis — use --force to recompute). "
              f"equibind={eq_variant} | diffdock={dd_variant}")
    else:
        df = load_poses(csv, eq_variant, ids, args.pb_valid_only, dd_variant)
        if df.empty:
            print("No poses after filtering (check --equibind-* / --diffdock-* / --ids-file).")
            return 1
        complexes = sorted(df["protein"].unique())
        if args.limit:
            complexes = complexes[:args.limit]
            df = df[df["protein"].isin(complexes)].copy()
        print(f"Complexes: {len(complexes)} | poses: {len(df)} | "
              f"tools: {sorted(df['tool'].unique())} | equibind={eq_variant} | "
              f"diffdock={dd_variant}")

        # ── parallel centroid extraction ─────────────────────────────────
        files = sorted(df["pose_file"].unique())
        cents: Dict[str, Tuple[float, float, float]] = {}
        print(f"Extracting centroids for {len(files)} pose SDFs ({args.workers} workers)...")
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
            for path, cx, cy, cz, na in ex.map(_extract_centroid, files, chunksize=64):
                if cx is not None:
                    cents[path] = (cx, cy, cz)
        print(f"  loaded {len(cents)}/{len(files)} centroids")

        # ── per-complex outlier removal (>outlier-dist from median) ──────
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

        do_placement = args.distance_mode in ("hybrid", "both")
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
                                  args.match_thr, do_placement, args.mode_rmsd_thr,
                                  args.rmsd_pose_cap, args.site_cluster,
                                  args.pocket_radius, args.rank_by, args.stability_boot)
            results.append(res)
            if ci % 50 == 0:
                print(f"  analyzed {ci}/{len(complexes)}")

        ok = [r for r in results if not r.get("skipped")]
        if not ok:
            print("No complex had >=2 usable poses.")
            return 1
        _save_analysis_cache(cache_path, sig, ok)
        print(f"  cached analysis → {cache_path} (re-run reuses it unless inputs change)")

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
                "radius_A": p["radius"], "spread_A": p.get("spread"),
                "conf_weight": p.get("conf_weight"),
                "center_x": round(float(p["center"][0]), 3),
                "center_y": round(float(p["center"][1]), 3),
                "center_z": round(float(p["center"][2]), 3),
                "dist_to_crystal": round(_dist(p["center"], cr), 3) if cr is not None else np.nan,
                "medoid_dist_to_crystal": (round(_dist(p.get("medoid_center", p["center"]), cr), 3)
                                           if cr is not None else np.nan),
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
        for (cid_, t, rp, cd, rm, ic, cl) in r.get("_per_rank", []):
            prrows.append({"protein": cid_, "tool": t, "rank": rp,
                           "centroid_dist": None if cd != cd else round(cd, 3),
                           "rmsd": None if rm != rm else round(rm, 3),
                           "in_correct_cluster": None if (isinstance(ic, float) and ic != ic) else bool(ic),
                           "cluster": None if cl is None else int(cl)})
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

    # ── clustering-quality metrics (primary partition) ───────────────────
    def _mean(df, col):
        v = pd.to_numeric(df.get(col), errors="coerce").dropna() if col in df else pd.Series(dtype=float)
        return float(v.mean()) if len(v) else np.nan
    print(f"\n  CLUSTERING QUALITY  (site-cluster={args.site_cluster}, "
          f"pocket-radius={args.pocket_radius:g} Å, rank-by={args.rank_by}):")
    k1_frac = float(pd.to_numeric(df_complex.get("k_is_one"), errors="coerce").mean()) \
        if "k_is_one" in df_complex else np.nan
    print(f"    mean #sites {_mean(df_complex,'n_clusters'):.2f}  "
          f"(legacy silhouette-KMedoids {_mean(df_complex,'km_n_clusters'):.2f}); "
          f"single-site (k=1) complexes {k1_frac:.0%}")
    print(f"    silhouette {_mean(df_complex,'silhouette'):.3f} (k≥2 only) | "
          f"Calinski-Harabasz {_mean(df_complex,'ch_score'):.1f} | "
          f"Davies-Bouldin {_mean(df_complex,'db_score'):.2f}")
    print(f"    compactness (median intra-spread) {_mean(df_complex,'compactness'):.2f} Å | "
          f"separation (nearest-site) {_mean(df_complex,'separation'):.2f} Å")
    if args.stability_boot > 0:
        print(f"    bootstrap stability (mean best-Jaccard, {args.stability_boot} resamples) "
              f"{_mean(df_complex,'boot_stability'):.2f}  (>0.75 stable, <0.5 dissolved)")
    print(f"    primary-vs-KMedoids agreement (ARI) {_mean(df_complex,'ari_primary_vs_kmedoids'):.2f}")

    # ── precision@1 ablation: does a better ranking recover the true site? ─
    print(f"\n  PRECISION@1 — rank-1 pocket within {args.match_thr:g} Å of crystal, per ranking rule:")
    print(f"    (oracle ceiling = ens_oracle_hit {_mean(dc,'ens_oracle_hit'):.0%})")
    for rule, lab in [("size", "size (legacy)"), ("ntools", "consensus (n_tools)"),
                      ("tight", "tightest"), ("confidence", "confidence-wtd"),
                      ("medoid", "consensus+medoid ctr")]:
        col = f"p1_{rule}_hit"
        if col in dc:
            print(f"    {lab:<24}{_mean(dc, col):>7.0%}")

    # ── crystal-closest-cluster purity (tight correct vs diffuse) ─────────
    print(f"\n  CRYSTAL-CLOSEST CLUSTER purity (are its members really near-native?):")
    print(f"    within {args.match_thr:g} Å (centroid) {_mean(dc,'correct_cluster_purity_centroid'):.0%} of members | "
          f"≤2 Å RMSD {_mean(dc,'correct_cluster_purity_rmsd'):.0%} of members | "
          f"best member RMSD median "
          f"{pd.to_numeric(dc.get('correct_cluster_best_rmsd'), errors='coerce').median():.2f} Å")

    # ── placement-aware (in-place RMSD) binding-MODE analysis, vs centroid sites ─
    if "n_modes" in df_complex and pd.to_numeric(df_complex["n_modes"], errors="coerce").notna().any():
        print(f"\n  PLACEMENT-AWARE binding modes  (in-place RMSD, mode-thr={args.mode_rmsd_thr:g} Å; "
              "translation+orientation+conformation):")
        print(f"    mean #modes/complex {_mean(df_complex,'n_modes'):.2f}  vs  "
              f"mean #sites {_mean(df_complex,'n_clusters'):.2f} (centroid)")
        mits = pd.to_numeric(df_complex.get('modes_in_top_site'), errors='coerce').dropna()
        if len(mits):
            print(f"    modes inside the biggest site: mean {mits.mean():.2f} | "
                  f">1 mode (site hides distinct orientations) in {(mits > 1).mean():.0%} of complexes")
        print(f"    mode compactness (median intra-mode RMSD) {_mean(df_complex,'mode_compactness'):.2f} Å | "
              f"separation {_mean(df_complex,'mode_separation'):.2f} Å")
        print(f"    site-vs-mode agreement: ARI {_mean(df_complex,'ari_site_vs_mode'):.2f} | "
              f"NMI {_mean(df_complex,'nmi_site_vs_mode'):.2f}  "
              "(low ⇒ one site holds several modes)")

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

    # ── ranking ablation (precision@1 per rule) + ranking enrichment ─────
    ablation = []
    for rule, desc in [("size", "size (legacy raw count)"),
                       ("ntools", "consensus (distinct tools)"),
                       ("tight", "tightest (robust spread)"),
                       ("confidence", "confidence-weighted"),
                       ("medoid", "consensus + medoid center"),
                       ("oracle", "oracle ceiling (best cluster)")]:
        if rule == "oracle":
            hcol, dcol = "ens_oracle_hit", "ens_oracle_dist"
        else:
            hcol, dcol = f"p1_{rule}_hit", f"p1_{rule}_dist"
        hv = (pd.to_numeric(dc[hcol].map({True: 1.0, False: 0.0}), errors="coerce")
              if hcol in dc else pd.Series(dtype=float))
        dv = pd.to_numeric(dc.get(dcol), errors="coerce") if dcol in dc else pd.Series(dtype=float)
        ablation.append({
            "rule": rule, "description": desc,
            "precision_at_1": round(float(hv.mean()), 4) if len(hv.dropna()) else None,
            "median_top1_dist": round(float(dv.median()), 3) if dv.notna().any() else None,
        })
    pd.DataFrame(ablation).to_csv(out_dir / "ranking_ablation.csv", index=False)

    # pooled ranking enrichment: Spearman(cluster rank, distance-to-crystal);
    # per complex there are too few clusters, so pool all clusters together.
    from scipy.stats import spearmanr
    rr, ddc = [], []
    for r in ok:
        cr = r.get("_crystal")
        if cr is None:
            continue
        for p in r["_pockets"]:
            rr.append(p["rank"]); ddc.append(_dist(p["center"], cr))
    ranking_rho = None
    if len(rr) > 10 and len(set(rr)) > 1:
        ranking_rho = float(spearmanr(rr, ddc)[0])
    print("\n  RANKING ENRICHMENT (pooled Spearman of cluster rank vs distance-to-crystal): "
          + (f"ρ={ranking_rho:.3f}  (positive = better-ranked clusters are nearer the crystal)"
             if ranking_rho is not None else "n/a (too few clusters)"))

    # paired-complex significance of the ablation (same complexes, k rules):
    #   Cochran's Q omnibus + Holm-adjusted pairwise exact McNemar family.
    ablation_sig = _ablation_paired_stats(dc)
    if ablation_sig is not None:
        mc = ablation_sig["consensus_vs_size_mcnemar"]
        print(f"  RANKING ABLATION SIGNIFICANCE (paired, n={ablation_sig['n_omnibus']}, "
              f"k={ablation_sig['k']} rules): Cochran Q={ablation_sig['cochran_q']:.2f} "
              f"p={_fmt_p(ablation_sig['p_omnibus'])}")
        if mc:
            print(f"    consensus vs size: {mc['consensus_wins']} wins / {mc['size_wins']} "
                  f"losses, McNemar p={_fmt_p(mc['p_raw'])} (Holm p={_fmt_p(mc['p_holm'])})")
        print("    pairwise McNemar (Holm-adjusted, a=better rule):")
        for pr in ablation_sig.get("pairwise_mcnemar_holm", []):
            print(f"      {pr['a']:>10} vs {pr['b']:<10} {pr['a_wins']:>3}↑/{pr['b_wins']:<3}↓"
                  f"  p={_fmt_p(pr['p_raw'])}  Holm={_fmt_p(pr['p_holm'])} {_p_stars(pr['p_holm'])}")

    summary = {
        "n_complexes": int(n), "match_thr_A": args.match_thr,
        "equibind_variant": eq_variant,
        "diffdock_variant": dd_variant,
        "clustering": {
            "site_method": args.site_cluster, "pocket_radius": args.pocket_radius,
            "rank_by": args.rank_by,
            "mean_n_clusters": _mean(df_complex, "n_clusters"),
            "mean_n_clusters_legacy_kmedoids": _mean(df_complex, "km_n_clusters"),
            "frac_single_site": (float(pd.to_numeric(df_complex["k_is_one"], errors="coerce").mean())
                                 if "k_is_one" in df_complex else None),
            "mean_silhouette": _mean(df_complex, "silhouette"),
            "mean_calinski_harabasz": _mean(df_complex, "ch_score"),
            "mean_davies_bouldin": _mean(df_complex, "db_score"),
            "mean_compactness_A": _mean(df_complex, "compactness"),
            "mean_separation_A": _mean(df_complex, "separation"),
            "mean_bootstrap_stability": _mean(df_complex, "boot_stability"),
            "mean_ari_primary_vs_kmedoids": _mean(df_complex, "ari_primary_vs_kmedoids"),
        },
        "ranking_ablation": ablation,
        "ranking_ablation_significance": ablation_sig,
        "ranking_enrichment_spearman": ranking_rho,
        "crystal_closest_cluster_purity": {
            "mean_purity_centroid": _mean(dc, "correct_cluster_purity_centroid"),
            "mean_purity_rmsd": _mean(dc, "correct_cluster_purity_rmsd"),
            "median_best_rmsd": (float(pd.to_numeric(dc["correct_cluster_best_rmsd"],
                                                     errors="coerce").median())
                                 if "correct_cluster_best_rmsd" in dc else None),
        },
        "placement_modes": {
            "mode_rmsd_thr_A": args.mode_rmsd_thr,
            "mean_n_modes": _mean(df_complex, "n_modes"),
            "mean_n_sites": _mean(df_complex, "n_clusters"),
            "mean_modes_in_top_site": _mean(df_complex, "modes_in_top_site"),
            "frac_top_site_multimode": (float((pd.to_numeric(df_complex["modes_in_top_site"],
                                                             errors="coerce") > 1).mean())
                                        if "modes_in_top_site" in df_complex else None),
            "mean_mode_compactness_A": _mean(df_complex, "mode_compactness"),
            "mean_mode_separation_A": _mean(df_complex, "mode_separation"),
            "mean_ari_site_vs_mode": _mean(df_complex, "ari_site_vs_mode"),
            "mean_nmi_site_vs_mode": _mean(df_complex, "nmi_site_vs_mode"),
        },
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
                _fig_top5(df_rank, df_complex, out_dir),   # returns a list of paths
                _fig_rank_in_crystal_cluster(df_rank, out_dir, eq_variant),
                _fig_rank_cluster_membership(df_rank, out_dir, eq_variant),
                _fig_rank_consensus(df_rank, out_dir, eq_variant),
                _fig_rank1_cluster_matrix(df_rank, out_dir, eq_variant),
                _fig_topN_crystal_matrix(df_rank, out_dir, eq_variant),
                _fig_topN_crystal_reach_curves(df_rank, out_dir, eq_variant),
                _fig_crystal_cluster_homogeneity(ok, df_rank, out_dir, eq_variant,
                                                 stats=args.stats),
                _fig_rank_in_correct(df_complex, args.match_thr, out_dir, eq_variant),
                _fig_near_native_composition(ok, args.match_thr, out_dir, eq_variant),
                _fig_ensembles(df_complex, args.match_thr, out_dir),
                _fig_cluster_quality(df_complex, ablation, ranking_rho, out_dir,
                                     args.match_thr, args.stability_boot),
                _fig_placement(df_complex, args.mode_rmsd_thr, out_dir),
                _fig_descriptor_quality(df_complex, args.features_csv, out_dir,
                                        per_pose_csv=args.per_pose_csv,
                                        dd_variant=dd_variant, eq_variant=eq_variant)]
        for f in figs:
            if not f:
                continue
            for pth in (f if isinstance(f, (list, tuple)) else [f]):
                print(f"  Figure: {pth}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
