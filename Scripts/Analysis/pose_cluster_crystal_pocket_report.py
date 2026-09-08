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
  * AutoDock  : pick the scoring function and the post-docking optimizer
                (``--autodock-scoring {vina,vinardo}`` /
                ``--autodock-optimize {none,gnina,gnina_refinement,smina}``), or
                pass a full method string with ``--autodock-variant``.
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

    # gnina-reranked AutoDock Vina (best AutoDock variant on the full-protein run)
    python Scripts/Analysis/pose_cluster_crystal_pocket_report.py \
        --autodock-variant autodock_gnina --diffdock-refine smina

Crystal-site convention (``--crystal-copies``, plan v2 decision D5):
  * ``reference`` (default, today's behaviour byte for byte): the crystal site is
    the heavy-atom centroid of the single ``<ID>_ligand.sdf`` reference instance.
  * ``any``: the crystal site is EVERY deposited copy of the ligand of interest
    (each record of ``<ID>_ligands.sdf``, record 0 = the reference instance, loaded
    on the same heavy-atom basis as the poses). Every crystal distance in the
    report becomes the minimum over copies; a cluster is correct through ANY
    copy — for every copy the cluster nearest to it counts when its centre lies
    within ``--match-thr`` of that copy (today's single-copy rule applied copy by
    copy); the crystal-closest cluster is the one nearest to the nearest copy. A
    cluster that is the nearest to two copies (six ids carry an adjacent copy
    4.9-7.8 A from the reference, inside the 8 A clustering radius) counts
    once: correctness is a SET of cluster labels, never a per-copy tally. The
    convention is part of the analysis-cache signature and, under ``any``, is
    printed in every stats ``*.txt`` header and in ``summary.json``; the
    reference outputs carry no such line and stay byte-identical.
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
import stats_utils as su                         # noqa: E402  (shared, unit-tested)
import method_filter as mf                       # noqa: E402  (shared exclusion filter)

from rdkit import Chem                            # noqa: E402
from rdkit.Chem import rdMolAlign                 # noqa: E402

# Cross-tool binding-site CONSENSUS here is defined over the three finalized
# tools (AutoDock, DiffDock, EquiBind): the ENSEMBLES enumeration, the
# distinct-tool agreement counters (values 1/2/3) and the pairwise site tables
# are all built for exactly these three. The new full-protein engines
# (autodock_vinardo, unidock, unidock2) are validated + compared per-tool by the
# PoseBusters validity, pose-comparison and PandaMap analyses; folding them into
# this spatial CONSENSUS would redefine "agreement" (a deliberate N-tool redesign
# of the counters/ensembles below), so they are intentionally not consensus tools
# here. Their palette entries are provided for forward-compatibility.
TOOLS = ("autodock", "diffdock", "equibind")
# Per-tool palette shared across every figure — matches 09f_pbvalid_yield_boxplot
# (posebusters_pose_comparison.py): AutoDock=blue, DiffDock=orange, EquiBind=green.
TOOL_COLORS = {
    "autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind": "#2ca02c",
    "autodock_vinardo": "#17becf", "unidock": "#9467bd", "unidock2": "#8c564b",
}
TOOL_MARKERS = {
    "autodock": "o", "diffdock": "^", "equibind": "s",
    "autodock_vinardo": "D", "unidock": "P", "unidock2": "X",
}


# ── crystal-site convention (D5) ─────────────────────────────────────────
# ``reference``: one crystal centroid (<ID>_ligand.sdf), a 1-D (3,) array.
# ``any``      : one centroid per deposited copy (<ID>_ligands.sdf), a 2-D (n, 3)
#                array whose row 0 is the reference instance. Every helper below
#                takes the minimum over rows when it meets a 2-D crystal, and
#                falls through to the exact reference arithmetic for a 1-D one.
_CRYSTAL_COPIES_CHOICES = ("reference", "any")
_CRYSTAL_COPIES = "reference"          # set once in main() from --crystal-copies
_TABLE_CONVENTION = None               # reference_convention of the per-pose table, if it carries one
_CRYSTAL_ANY_RULE = (
    "crystal site: any deposited copy (every record of <ID>_ligands.sdf, record 0 = "
    "reference instance; every crystal distance is the minimum over copies; a "
    "cluster is correct when it is the cluster nearest to ANY copy and lies within "
    "the threshold of it; the crystal-closest cluster is the one nearest to the "
    "nearest copy; a cluster that is nearest to two copies counts once)")


def _crystal_site_lines() -> List[str]:
    """Header line(s) naming the crystal-site convention for the stats sidecars.

    Empty under ``reference`` so that every existing header stays byte-identical
    (the absence of the line IS the reference-instance wording); one explicit
    rule line under ``any`` so the two conventions can never be confused."""
    if _CRYSTAL_COPIES != "any":
        return []
    return [_CRYSTAL_ANY_RULE,
            f"per-pose table reference convention: {_TABLE_CONVENTION or 'not recorded (pre-schema-6 table)'}"]


def _crystal_site_header(sep: str = "\n") -> str:
    """``_crystal_site_lines`` joined for f-string headers ('' under reference)."""
    lines = _crystal_site_lines()
    return (sep.join(lines) + sep) if lines else ""


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

def _strip_all_hs(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Return ``mol`` on the hub's heavy-atom basis, or None if RDKit cannot.

    Sanitising first is preferred (it also fixes up implicit-H counts), but the
    refined pose files carry geometries RDKit occasionally refuses to sanitise
    (boron ligands, odd valences from the refiner), so a ``sanitize=False``
    strip is the fallback.
    """
    # Mirror the hub (posebusters_pose_comparison.py: load_first_mol sanitises,
    # then Chem.RemoveHs) so this report and per_pose_metrics.csv share ONE
    # heavy-atom basis: RemoveHs keeps a stereo-defining hydrogen on two
    # ligands (5SAK_ZRY, 8D5D_5DK) and RemoveAllHs would not, which is why the
    # two generators disagreed by up to 0.22 A on those complexes. RemoveAllHs
    # stays as the last resort for a pose RDKit cannot sanitise.
    try:
        return Chem.RemoveHs(mol)
    except Exception:
        try:
            return Chem.RemoveHs(mol, sanitize=False)
        except Exception:
            try:
                return Chem.RemoveAllHs(mol, sanitize=False)
            except Exception:
                return None


def load_heavy_atom_mol(sdf_path: str) -> Optional[Chem.Mol]:
    # INVARIANT: every Mol returned here is heavy-atom only. The crystal ligand
    # (<ID>_ligand.sdf) carries no hydrogens while the refined docking poses do,
    # and every centroid / RMSD / atom-count comparison downstream assumes both
    # sides are on the same heavy-atom basis. NOTE ``removeHs=True`` on
    # MolFromMolFile is a NO-OP when ``sanitize=False`` (RDKit only strips Hs as
    # part of sanitisation), so the explicit strip below is what enforces it —
    # relying on the reader flag alone left 28-30 of 30 sampled refined poses
    # with their hydrogens (Phase 0.3 defect, plan v2 Section 1.7).
    try:
        mol = Chem.MolFromMolFile(sdf_path, removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)
        if mol is None or mol.GetNumConformers() == 0:
            return None
        if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
            mol = _strip_all_hs(mol)
            if mol is None or mol.GetNumConformers() == 0:
                return None
        return mol
    except Exception:
        return None


def load_heavy_atom_mols_all(sdf_path: str) -> List[Chem.Mol]:
    """Every record of a multi-record SDF on the SAME heavy-atom basis as
    ``load_heavy_atom_mol`` (which reads only the first record). Used for the
    deposited copies in ``<ID>_ligands.sdf`` under ``--crystal-copies any``; record
    order is preserved (record 0 = the reference instance). Unreadable or
    conformer-less records are skipped."""
    out: List[Chem.Mol] = []
    try:
        sup = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
    except Exception:
        return out
    for mol in sup:
        if mol is None or mol.GetNumConformers() == 0:
            continue
        if any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
            mol = _strip_all_hs(mol)
            if mol is None or mol.GetNumConformers() == 0:
                continue
        out.append(mol)
    return out


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
    d = _pose_cdists(C[idx], crystal)
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
        dd = _cdist(p["center"], crystal)
        out[f"p1_{name}_dist"] = round(dd, 3)
        out[f"p1_{name}_hit"] = bool(dd <= thr)
    # medoid-centered: the consensus pick but measured from its robust medoid center
    pm = picks["ntools"]
    dmd = _cdist(pm.get("medoid_center", pm["center"]), crystal)
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
    pose_d = _pose_cdists(sc, crystal)
    oracle = float(pose_d.min())
    if len(sc) < 2:
        top1_center = sc.mean(axis=0)
        n_clusters = 1
    else:
        lab, _, _ = cluster_sites(_centroid_dm(sc), sc, site_method, pocket_radius)
        pk = pockets_from_labels(sc, lab, st, rank_by="consensus")
        top1_center = pk[0]["center"]
        n_clusters = len(pk)
    top1 = _cdist(top1_center, crystal)
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

# AutoDock encodes its variant in the method NAME, as
# ``autodock[_<scoring>][_<optimizer>]``: the scoring token is present only for
# the non-default function (``vinardo``), and the optimizer token only when the
# Vina poses were post-processed (gnina/smina minimisation + re-ranking, and the
# CNN-gradient 'gnina_refinement' flavour). So the six full-protein methods are
# autodock, autodock_gnina, autodock_gnina_refinement, autodock_vinardo,
# autodock_vinardo_gnina, autodock_vinardo_gnina_refinement.
_AD_SCORINGS = ("vina", "vinardo")
_AD_OPTIMIZERS = ("none", "gnina", "gnina_refinement", "smina")


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


def resolve_autodock_variant(cat: pd.DataFrame, variant: Optional[str],
                             scoring: str, optimize: str) -> str:
    """Resolve the exact AutoDock method string for the AutoDock tool slot.

    ``--autodock-variant`` accepts either a full method string or a bare optimizer
    keyword (``gnina`` -> ``autodock_gnina``, kept alongside the current
    ``--autodock-scoring``); ``--autodock-scoring``/``--autodock-optimize`` are the
    tidy way to say the same thing. Composition follows the naming rule above, so
    a variant that is not in the CSV fails loudly with the available list instead
    of silently emptying the AutoDock slot.
    """
    methods = sorted(m for m in cat["method"].astype(str).unique()
                     if m == "autodock" or m.startswith("autodock_"))
    if not methods:
        raise SystemExit("No 'autodock' poses in the per-pose CSV.")
    if variant:
        if variant in methods:
            return variant
        if variant in _AD_OPTIMIZERS:        # shorthand: --autodock-variant gnina
            optimize = variant
        else:
            raise SystemExit(
                f"--autodock-variant '{variant}' not recognised (want a full method "
                f"name or one of {_AD_OPTIMIZERS}). Available:\n    "
                + "\n    ".join(methods))
    cand = "autodock" + ("_vinardo" if scoring == "vinardo" else "") \
                      + ("" if optimize == "none" else f"_{optimize}")
    if cand not in methods:
        raise SystemExit(
            f"No AutoDock variant (scoring={scoring}, optimize={optimize} -> '{cand}') "
            f"in the CSV. Available:\n    " + "\n    ".join(methods))
    return cand


# The thesis draws its three headline slots as AutoDock*, DiffDock* and EquiBind*,
# a convention defined once in the Results as AutoDock Vina + gnina, DiffDock +
# smina and unguided EquiBind + gnina. A starred label therefore means exactly
# that arm and nothing else; every other variant is spelled out so it can never be
# read as the starred default. Before 2026-08-21 only the AutoDock slot was
# relabelled here, which left it as the single spelled-out entry on axes where
# DiffDock and EquiBind were equally optimised but drawn bare.
_STARRED_VARIANT = {
    # The dominant AutoDock arm as of 2026-08-29: ADFRsuite-prepared ligands,
    # exhaustiveness 128, gnina-rescored. The Meeko arms it replaced are excluded
    # from the Benchmark analyses entirely (see Scripts/Analysis/method_filter.py).
    "autodock": "autodock_mgltools_exh128_gnina",
    "diffdock": "diffdock_smina",
    "equibind": "equibind_unguided_gnina",
}


def _autodock_display(ad_variant: str) -> str:
    """Figure label for the AutoDock slot, e.g. ``autodock_gnina`` -> 'AutoDock*'
    — so a plot never claims plain Vina while showing gnina-reranked poses."""
    name = str(ad_variant)
    if name == _STARRED_VARIANT["autodock"]:
        return "AutoDock*"
    base = "AutoDock Vinardo" if "vinardo" in name else "AutoDock Vina"
    # A non-dominant MGLTools/ADFRsuite arm must not be drawn as plain "AutoDock
    # Vina + gnina": that label belongs to the Meeko ligand prep, and the two are a
    # controlled A/B. Carry the prep and the search effort into the label instead.
    if name.startswith("autodock_mgltools"):
        base += " MGLTools-lig"
        for eff in ("18", "64", "92", "128"):
            if f"_exh{eff}" in name:
                base += f" exh{eff}"
                break
        else:
            base += " exh32"
    if name.endswith("_gnina_refinement"):
        return f"{base} + gnina (refine)"
    for opt in ("gnina", "smina"):
        if name.endswith(f"_{opt}"):
            return f"{base} + {opt}"
    return base


def _diffdock_display(dd_variant: str) -> str:
    """Figure label for the DiffDock slot, mirroring ``_autodock_display``."""
    name = str(dd_variant)
    if name == _STARRED_VARIANT["diffdock"]:
        return "DiffDock*"
    for opt in ("gnina", "smina"):
        if name.endswith(f"_{opt}"):
            return f"DiffDock + {opt}"
    return "DiffDock"


def _equibind_display(eq_variant: str) -> str:
    """Figure label for the EquiBind slot, mirroring ``_autodock_display``.

    EquiBind folds pocket source, refiner and clamp into the method name, so a
    non-starred arm is shown with that suffix rather than guessed at."""
    name = str(eq_variant)
    if name == _STARRED_VARIANT["equibind"]:
        return "EquiBind*"
    if name.startswith("equibind_"):
        return f"EquiBind ({name[len('equibind_'):].replace('_', ' ')})"
    return "EquiBind"


def _map_tool(method: str, eq_variant: str, dd_variant: str = "diffdock",
              ad_variant: str = "autodock") -> Optional[str]:
    if method == ad_variant:                 # autodock | autodock_gnina | autodock_vinardo | ...
        return "autodock"
    if method == dd_variant:                 # diffdock | diffdock_smina | diffdock_gnina
        return "diffdock"
    if method == eq_variant:
        return "equibind"
    return None


def load_poses(csv: Path, eq_variant: str, ids: Optional[set],
               pb_valid_only: bool, dd_variant: str = "diffdock",
               ad_variant: str = "autodock") -> pd.DataFrame:
    df = pd.read_csv(csv, low_memory=False)
    df["tool"] = df["method"].map(lambda m: _map_tool(m, eq_variant, dd_variant, ad_variant))
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

    A tool's native ``rank`` is used when it varies (AutoDock = Vina mode, or the
    optimizer's re-ranked order for the gnina/smina AutoDock variants, whose
    ``rank`` column already carries ``optimized_rank``;
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


def crystal_centroid(benchmark_dir: Path, cid: str,
                     copies: str = "reference") -> Optional[np.ndarray]:
    """Crystal site(s) of ``cid`` on the heavy-atom basis of the poses.

    ``copies="reference"``: the centroid of the single ``<ID>_ligand.sdf``
    reference instance as a 1-D (3,) array (today's behaviour, unchanged).
    ``copies="any"`` (D5): one centroid per record of ``<ID>_ligands.sdf`` as a
    2-D (n, 3) array, row 0 = the reference instance. Falls back to the single
    reference file (as a (1, 3) array) when the multi-record file is absent or
    unreadable, so a complex never silently changes basis inside one run.
    """
    sdf = benchmark_dir / cid / f"{cid}_ligand.sdf"
    if copies == "any":
        sdfs = benchmark_dir / cid / f"{cid}_ligands.sdf"
        if sdfs.exists():
            mols = load_heavy_atom_mols_all(str(sdfs))
            if mols:
                return np.asarray([centroid_from_mol(m) for m in mols], dtype=float)
        if not sdf.exists():
            return None
        mol = load_heavy_atom_mol(str(sdf))
        if mol is None:
            return None
        return np.asarray([centroid_from_mol(mol)], dtype=float)
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


def _cdist(a, crystal) -> float:
    """Distance from a point to the crystal site: the plain distance for a 1-D
    (reference) crystal, the MINIMUM over copies for a 2-D (any-copy) one."""
    crystal = np.asarray(crystal)
    if crystal.ndim == 1:
        return _dist(a, crystal)
    return float(np.linalg.norm(crystal - np.asarray(a, dtype=float), axis=1).min())


def _pose_cdists(P: np.ndarray, crystal) -> np.ndarray:
    """Vectorised ``_cdist`` for an (m, 3) block of pose centroids -> (m,)."""
    crystal = np.asarray(crystal)
    if crystal.ndim == 1:
        return np.linalg.norm(P - crystal, axis=1)
    return np.linalg.norm(P[:, None, :] - crystal[None, :, :], axis=2).min(axis=1)


def _nearest_copy_index(a, crystal) -> int:
    """Row of the copy nearest to ``a`` (0 for a 1-D reference crystal)."""
    crystal = np.asarray(crystal)
    if crystal.ndim == 1:
        return 0
    return int(np.argmin(np.linalg.norm(crystal - np.asarray(a, dtype=float), axis=1)))


def _pocket_crystal_stats(pockets: List[dict], crystal: np.ndarray,
                          thr: float, key="center") -> dict:
    """top1 (rank-1) and oracle (best) center-to-crystal stats for an ordered list."""
    if not pockets or crystal is None:
        return dict(top1_dist=np.nan, top1_hit=np.nan,
                    oracle_dist=np.nan, oracle_rank=np.nan, oracle_hit=np.nan)
    dists = [_cdist(p[key], crystal) for p in pockets]
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
    # centroid distance to crystal from our own centroids (fallback to CSV col;
    # the CSV column follows the hub's own reference convention)
    def cdist(row):
        c = cents.get(row["pose_file"])
        if c is not None and crystal is not None:
            return _cdist(c, crystal)
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
                   crystal: Optional[np.ndarray], correct_labels=None,
                   max_rank: int = 5, per_rank_max: int = _RANK_PROFILE_MAX):
    """Per-tool profile of the top ranked poses.

    Returns (per_rank_rows, tool_summary):
      * per_rank_rows: (tool, rank_pos, centroid_dist, rmsd, in_correct_cluster) for
        the k-th best pose of each tool, up to ``per_rank_max`` ranks — the "how far
        is the rank-k pose from crystal" profile, plus whether that pose landed in
        a correct cluster (``correct_labels``: the SET of correct cluster labels —
        exactly one under the reference convention, possibly several under D5
        ``--crystal-copies any``; NaN when no crystal / cluster is defined).
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
            cd = _cdist(C[idx], crystal) if (idx is not None and crystal is not None) else np.nan
            rm = float(file_rmsd.get(f, np.nan))
            clab = int(labels[idx]) if idx is not None else None
            in_correct = (bool(clab in correct_labels)
                          if (idx is not None and correct_labels) else np.nan)
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
    # The cluster whose center is nearest the crystal is the "correct" site, but
    # only when that center actually lands within `thr` of the crystal ligand. A
    # nearest cluster further away than that is not a recovery of the site, so the
    # complex carries NO correct cluster at all: `correct_label` stays None, every
    # tool's `in_correct_cluster` is NaN rather than False, and the complex
    # contributes neither reach nor cluster composition. `correct_cluster_is_hit`
    # is what gates this, and it is the same threshold the oracle ceiling uses.
    # For each tool we report the best (lowest) rank it assigns to a pose that
    # landed in that cluster — i.e. does the tool prioritize its near-native pose?
    #
    # D5 (``--crystal-copies any``, 2-D ``crystal`` of copy centroids): every
    # distance here is the minimum over copies (``_cdist``). The crystal-closest
    # cluster ``cp`` is the one nearest to the NEAREST copy, and the correct SET
    # holds, for every copy, the cluster nearest to that copy when it lies within
    # `thr` of it (a cluster is correct through ANY copy). Correctness is a set of
    # cluster LABELS, so a cluster that is the nearest to two copies (the six ids
    # whose adjacent copy sits 4.9-7.8 A from the reference, inside the 8 A
    # clustering radius) enters the set once and every pose in it is counted once
    # — never a per-copy tally. Under the reference convention the set is exactly
    # {cp}, so every quantity below reduces to today's arithmetic.
    correct_dist = correct_is_hit = np.nan
    pur_centroid = pur_rmsd = best_rmsd_in_correct = np.nan
    p1 = {}
    correct_labels: Optional[set] = None
    correct_pockets: List[dict] = []
    rank_in_correct: Dict[str, Optional[int]] = {}
    n_in_correct: Dict[str, int] = {}
    correct_members: List[Tuple[str, Optional[int]]] = []
    if crystal is not None and pockets:
        cp = min(pockets, key=lambda p: _cdist(p["center"], crystal))
        correct_dist = round(_cdist(cp["center"], crystal), 3)
        correct_is_hit = bool(correct_dist <= thr)
        p1 = _precision_at_1(pockets, crystal, thr)
    if correct_is_hit is True:
        # Today's rule per copy: for EVERY copy the cluster nearest to it is correct
        # when its centre lies within thr of that copy (any other cluster that also
        # happens to sit within thr is NOT promoted — exactly today's single-copy
        # rule, applied copy by copy). The labels form a SET, so a cluster that is
        # the nearest to two copies enters once. With one copy (the reference
        # convention, or a single-copy complex under ``any``) the per-copy nearest
        # cluster is cp itself and the set is {cp}: identical to today's arithmetic.
        copies = np.atleast_2d(np.asarray(crystal, dtype=float))
        correct_labels = set()
        for j in range(copies.shape[0]):
            pj = min(pockets, key=lambda p: _dist(p["center"], copies[j]))
            # rounded to 3 dp before the threshold test, exactly like correct_dist
            if (round(_dist(pj["center"], copies[j]), 3) <= thr
                    and int(pj["label"]) not in correct_labels):
                correct_labels.add(int(pj["label"]))
                correct_pockets.append(pj)
        correct_pockets.sort(key=lambda p: p["rank"])       # deterministic member order
        correct_idx = [i for p in correct_pockets for i in p["members"]]
        pur_centroid, pur_rmsd, best_rmsd_in_correct = _cluster_purity(
            cp, C, crystal, files, rmsd_map, thr)
        for t in TOOLS:
            rk = [eff_rank.get(files[i]) for i in correct_idx
                  if tools[i] == t and eff_rank.get(files[i]) is not None]
            rank_in_correct[t] = int(min(rk)) if rk else None
            n_in_correct[t] = len(rk)
        # Every member pose of the correct cluster(s) as (tool, effective rank),
        # for the near-native-cluster composition figure (rank make-up + tool consensus).
        correct_members = [(tools[i],
                            (int(eff_rank[files[i]]) if eff_rank.get(files[i]) is not None else None))
                           for i in correct_idx]

    # ── top-5 ranked poses: per-rank distance + cluster consistency ──────
    file_idx = {f: i for i, f in enumerate(files)}
    per_rank_rows, top5 = _top5_analysis(sub, eff_rank, file_idx, labels, C, crystal,
                                         correct_labels=correct_labels)

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
        # D5 bookkeeping — present ONLY under --crystal-copies any (2-D crystal),
        # so the reference per_complex_summary.csv keeps its exact column set
        **({"crystal_n_copies": int(np.asarray(crystal).shape[0]),
            "correct_cluster_copy_index": (_nearest_copy_index(cp["center"], crystal)
                                           if correct_is_hit is True else np.nan),
            "n_correct_clusters": len(correct_labels) if correct_labels else 0}
           if (crystal is not None and np.asarray(crystal).ndim == 2) else {}),
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
        cr = np.asarray(r["_crystal"])
        if cr.ndim == 2:
            # D5: draw every deposited copy — the reference instance (row 0) as
            # today's black X, every alternate copy as a grey X.
            ax.scatter(*cr[0], marker="X", s=240, color="black", label="crystal", zorder=10)
            for alt in cr[1:]:
                ax.scatter(*alt, marker="X", s=200, color="0.45", label="alt. copy",
                           zorder=10)
        else:
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


def _reach_stats_sink(out_dir, key, payload):
    """Merge one figure's numeric stats into ``descriptor_reach_stats.json``.

    Kept separate from the homogeneity figure's dedicated
    ``crystal_cluster_homogeneity_stats.json``; collects the lightweight per-figure
    tests added for the summary / top-N / descriptor figures so the annotated
    numbers stay recoverable. Best-effort: a write failure never breaks the run.
    """
    try:
        p = Path(out_dir) / "descriptor_reach_stats.json"
        blob = {}
        if p.exists():
            try:
                blob = json.loads(p.read_text())
            except Exception:
                blob = {}
        blob[key] = payload
        p.write_text(json.dumps(blob, indent=2, default=str))
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] could not write descriptor_reach_stats.json: {e}")


def _hit_rate_paired(d, thr):
    """Paired-binary comparison of 'oracle center within thr Å' across sources.

    Per-complex unit of analysis; a NaN oracle (source produced no pose) counts as
    a miss so every complex stays a complete row across all sources. Returns the
    ``su.paired_proportions`` dict (Cochran Q omnibus + pairwise McNemar/Holm +
    Wilson CIs) or None when there are too few sources / complexes.
    """
    labels = [s for s in _SRC_ORACLE if _SRC_ORACLE[s] in d]
    cols = {s: (pd.to_numeric(d[_SRC_ORACLE[s]], errors="coerce") <= thr)
                 .fillna(False).astype(int).to_numpy()
            for s in labels}
    if len(cols) < 2 or len(d) < 8:
        return None
    return su.paired_proportions(cols, labels=list(cols))


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
    labels = list(_SRC_ORACLE)
    hits, tots, oracle_rates = [], [], []
    for s in labels:
        v = pd.to_numeric(d[_SRC_ORACLE[s]], errors="coerce").dropna()
        n_s = len(v); k_s = int((v <= thr).sum())
        hits.append(k_s); tots.append(n_s)
        oracle_rates.append(k_s / n_s if n_s else np.nan)
    x = np.arange(len(labels))
    ax.bar(x, oracle_rates, color=[TOOL_COLORS.get(s, "#8172B3") for s in labels])
    # Wilson 95% CI error bars (same k/n as each displayed rate).
    tops = list(oracle_rates)
    try:
        yerr = np.zeros((2, len(labels)))
        for i, (r, k, n) in enumerate(zip(oracle_rates, hits, tots)):
            if r == r and n:
                lo, hi = su.wilson_ci(k, n)
                yerr[0, i] = max(0.0, r - lo); yerr[1, i] = max(0.0, hi - r)
                tops[i] = hi
        ax.errorbar(x, [r if r == r else 0 for r in oracle_rates], yerr=yerr,
                    fmt="none", ecolor="0.25", elinewidth=1.1, capsize=3, zorder=5)
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] hit-rate Wilson CI skipped: {e}")
    for xi, v, t in zip(x, oracle_rates, tops):
        if v == v:
            ax.text(xi, min(1.12, t + 0.02), f"{v:.0%}", ha="center", fontsize=9)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylim(0, 1.18); ax.set_ylabel("fraction within thr")
    # Paired-binary omnibus across sources (same complexes; per-complex unit).
    title = f"Finds the true site (oracle ≤ {thr:g} Å)"
    try:
        res = _hit_rate_paired(d, thr)
        if res is not None:
            o = res["omnibus"]
            title += (f"\nCochran Q={o['Q']:.1f}, {su.fmt_p(o['p'])} "
                      f"{su.p_stars(o['p'])} (paired across sources)")
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] hit-rate paired omnibus skipped: {e}")
    ax.set_title(title)


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

    # ── sidecar: full paired-binary stats behind summary_hit_rate.png ────
    try:
        res = _hit_rate_paired(d, thr)
        if res is not None:
            _reach_stats_sink(out_dir, "summary_hit_rate", {
                "thr_A": float(thr), "n_complexes": int(len(d)),
                "unit": "per-complex boolean (oracle center ≤ thr); NaN oracle = miss",
                "rates": res["rates"], "omnibus": res["omnibus"],
                "pairwise": res["pairwise"]})
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] summary_hit_rate sidecar skipped: {e}")

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
_SRC_COLOR = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind": "#2ca02c",
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


# Figure captions (footnotes) for the top-N matrix + its companion stats table.
# Kept off the figures themselves and written to topN_crystal_cluster_matrix_caption.txt.
_TOPN_MATRIX_CAPTION = (
    "Each cell: co-reach % (both tools land >=1 top-N pose in the crystal cluster) "
    "over the count of such complexes; diagonal = one tool alone. Inferential "
    "statistics (Wilson CIs, signed co-reach phi, permutation tests) are in the "
    "companion table topN_crystal_cluster_matrix_stats.png.")
_TOPN_STATS_TABLE_CAPTION = (
    "Reach rate = fraction of complexes with >=1 top-N pose in the crystal-closest "
    "cluster (Wilson 95% CI). Co-reach phi = signed association of the two tools' "
    "per-complex reach vs independence (two-sided Fisher; + redundant, - complementary; "
    "n/e when a tool sits at a reach ceiling). All three = observed vs expected joint "
    "reach, 10 000-shuffle permutation p. Depths are nested cumulative thresholds - not "
    "an independent test family.")


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
    # per-tool presence over the fixed complex order (reused every bucket)
    present_vec = {t: np.array([t in present.get(p, set()) for p in proteins], bool)
                   for t in T}
    matrix_stats = {}
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
                           ha="right", rotation_mode="anchor", fontsize=12)
        ax.set_yticklabels([_TOOL_DISPLAY.get(t, t) for t in T] if col == 0
                           else [""] * nT, fontsize=12)
        # ── per-cell inferential stats (same helpers as the reach-curves companion):
        #    diagonal  → Wilson 95% CI on the marginal reach rate;
        #    off-diag  → SIGNED phi / Fisher co-reach association vs independence
        #                (AutoDock sits at a ceiling → 'fails' margin too thin → n/e).
        #    Depths are NESTED cumulative thresholds, so they are not corrected as a
        #    family; the diagonal marginal-RATE comparison (Cochran Q/McNemar) lives
        #    in the topN_crystal_reach_curves companion, not duplicated here.
        #    These feed the companion table topN_crystal_cluster_matrix_stats.png; the
        #    heatmap cells themselves stay purely descriptive (co-reach % + count).
        bstat = {"diagonal": {}, "offdiagonal": {}}
        try:
            hitv = {t: np.array([(reach.get(p, {}).get(t, np.inf) <= N)
                                 for p in proteins], int) for t in T}
            for a in range(nT):
                ta = T[a]; mask = present_vec[ta]
                k = int(hitv[ta][mask].sum()); nn = int(mask.sum())
                lo, hi = su.wilson_ci(k, nn) if nn else (np.nan, np.nan)
                bstat["diagonal"][ta] = {"k": k, "n": nn,
                                         "rate": (k / nn if nn else None),
                                         "wilson_ci": [float(lo), float(hi)]}
            for a in range(nT):
                for b in range(a):                       # strictly lower triangle
                    ta, tb = T[a], T[b]; m2 = present_vec[ta] & present_vec[tb]
                    assoc = su.paired_2x2_association(hitv[ta][m2], hitv[tb][m2])
                    bstat["offdiagonal"][f"{ta}+{tb}"] = assoc
            if n_all:
                mask_all = present_vec[T[0]] & present_vec[T[1]] & present_vec[T[2]]
                allc = su.consensus_perm_test(
                    np.column_stack([hitv[t][mask_all] for t in T]), n_perm=10000)
                bstat["all_three"] = allc
        except Exception as e:                           # pragma: no cover
            print(f"  [stats] topN matrix cell stats skipped (top-{N}): {e}")
            allc = None
        matrix_stats[f"top{N}"] = bstat
        for a in range(nT):
            for b in range(nT):
                if b > a or both[a, b] == 0:
                    continue
                ax.text(b, a, f"{M[a, b]:.0%}\n{int(cnt[a, b])}", ha="center",
                        va="center",
                        color="white" if M[a, b] > 0.55 else "black", fontsize=12)
        if n_all:                            # descriptive three-way count (perm p → table)
            ax.text(0.96, 0.96, "All 3 tools reach crystal:\n"
                    f"{all3}/{n_all} ({all3 / n_all:.0%})",
                    transform=ax.transAxes, ha="right", va="top", fontsize=11,
                    bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbbbbb"))
        ax.set_title(f"Top-{N}", fontsize=14)
    for j in range(nb, nrows * ncols):       # hide any unused grid cell
        flat[j].axis("off")
    if im is not None:
        fig.colorbar(im, ax=list(flat), fraction=0.046, pad=0.02,
                     label="fraction of complexes")
    _reach_stats_sink(out_dir, "topN_crystal_cluster_matrix", {
        "n_complexes": int(n_total),
        "unit": "per-complex boolean (best rank in crystal cluster ≤ N)",
        "note": "diagonal = Wilson CI on marginal reach; off-diagonal = signed-phi / "
                "Fisher co-reach association vs independence (+ redundant, − complementary); "
                "depths are nested cumulative thresholds (not an independent family); the "
                "marginal-rate omnibus (Cochran Q/McNemar) is in topN_crystal_reach_curves.",
        "buckets": matrix_stats})
    fig.suptitle("Both tools reach the CRYSTAL cluster within top-N poses "
                 "(diagonal = one tool alone)\n"
                 f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}", fontsize=12)
    p = out_dir / "topN_crystal_cluster_matrix.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    tbl_path = _fig_topN_crystal_matrix_stats_table(
        matrix_stats, buckets, n_total, out_dir, eq_variant)
    # figure captions (footnotes) live in a txt sidecar, not on the figures
    cap = ("topN_crystal_cluster_matrix.png / topN_crystal_cluster_matrix_stats.png"
           " — figure captions\n"
           f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}\n"
           + _crystal_site_header()
           + "=" * 72 + "\n\n"
           "topN_crystal_cluster_matrix.png:\n" + _TOPN_MATRIX_CAPTION + "\n")
    if tbl_path:
        cap += ("\ntopN_crystal_cluster_matrix_stats.png:\n"
                + _TOPN_STATS_TABLE_CAPTION + "\n")
    cap_path = out_dir / "topN_crystal_cluster_matrix_caption.txt"
    cap_path.write_text(cap)
    print(f"  Caption TXT: {cap_path}")
    return [p, tbl_path] if tbl_path else p


def _fig_topN_crystal_matrix_stats_table(matrix_stats, buckets, n_total, out_dir,
                                         eq_variant=None):
    """Companion statistics table for topN_crystal_cluster_matrix.png.

    The matrix panels now carry only the descriptive co-reach %/counts; every
    inferential quantity that used to crowd the cells lives here instead, one
    column per depth bucket:
      · diagonal  → marginal crystal-reach rate + Wilson 95 % CI (one row per tool);
      · off-diag  → signed co-reach association φ + two-sided Fisher p (one row per
                    tool pair; 'n/e' when a tool sits at a reach ceiling);
      · all three → observed vs expected joint reach + permutation p.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if not matrix_stats:
        return None
    T = list(TOOLS)
    depth_keys = [f"top{N}" for N in buckets]
    col_labels = ["Statistic"] + [f"Top-{N}" for N in buckets]

    def _pair_assoc(bstat, t1, t2):
        od = bstat.get("offdiagonal", {})
        return od.get(f"{t1}+{t2}") or od.get(f"{t2}+{t1}")

    rows = []                                    # (label, [cell per depth], group_key)
    # diagonal — marginal crystal-reach rate + Wilson CI
    for t in T:
        cells = []
        for dk in depth_keys:
            d = matrix_stats.get(dk, {}).get("diagonal", {}).get(t)
            if d and d.get("n") and d.get("rate") is not None:
                lo, hi = d["wilson_ci"]
                cells.append(f"{d['rate']:.0%}  [{lo:.0%}–{hi:.0%}]")
            else:
                cells.append("—")
        rows.append((f"{_TOOL_DISPLAY.get(t, t)} — reach rate", cells, "reach"))
    # off-diagonal — signed co-reach association φ + Fisher p
    pairs = [("autodock", "diffdock"), ("autodock", "equibind"),
             ("diffdock", "equibind")]
    for t1, t2 in pairs:
        cells = []
        for dk in depth_keys:
            a = _pair_assoc(matrix_stats.get(dk, {}), t1, t2)
            if a and a.get("estimable"):
                cells.append(f"φ={a['phi']:+.2f}  {su.fmt_p(a['p'])} {su.p_stars(a['p'])}")
            elif a:
                cells.append("n/e (ceiling)")
            else:
                cells.append("—")
        rows.append((f"{_TOOL_DISPLAY[t1]} + {_TOOL_DISPLAY[t2]} — co-reach φ",
                     cells, "coreach"))
    # all three — observed vs expected joint reach + permutation p
    cells = []
    for dk in depth_keys:
        allc = matrix_stats.get(dk, {}).get("all_three")
        if allc:
            k = allc["k"]
            cells.append(f"{int(allc['observed'][k])} obs / {allc['expected'][k]:.0f} exp  "
                         f"perm {su.fmt_p(allc['p_all_agree_two_sided'])}")
        else:
            cells.append("—")
    rows.append(("All three reach — obs vs exp (perm)", cells, "all3"))

    group_fc = {"reach": "#eef3f8", "coreach": "#f4f0f7", "all3": "#f0f5ee"}
    cell_text = [[lab] + cs for lab, cs, _ in rows]
    fig_w = 3.4 + 2.7 * len(buckets)
    fig_h = 0.5 * (len(rows) + 1) + 1.7
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.axis("off")
    tbl = ax.table(cellText=cell_text, colLabels=col_labels, loc="upper center",
                   cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8.5); tbl.scale(1, 1.6)
    for (r_i, c_i), cell in tbl.get_celld().items():
        cell.set_edgecolor("0.85")
        if r_i == 0:                             # header row
            cell.set_text_props(fontweight="bold"); cell.set_facecolor("#e6e9ec")
        else:
            cell.set_facecolor(group_fc.get(rows[r_i - 1][2], "white"))
            if c_i == 0:                         # statistic label column
                cell.get_text().set_ha("left"); cell.PAD = 0.03
                cell.set_text_props(fontweight="bold")
    fig.suptitle("Top-N crystal-cluster reach — inferential statistics  "
                 "(companion to topN_crystal_cluster_matrix.png)\n"
                 f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}",
                 fontsize=11)   # caption/footnote → topN_crystal_cluster_matrix_caption.txt
    fig.subplots_adjust(left=0.02, right=0.98, top=0.82, bottom=0.05)
    p = out_dir / "topN_crystal_cluster_matrix_stats.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    # plain-text sidecar of the same table (monospace, aligned) so the numbers
    # are copy-pasteable without transcribing the PNG.
    header = (col_labels, *cell_text)
    widths = [max(len(str(r[c])) for r in header) for c in range(len(col_labels))]
    def _line(vals):
        return "  ".join(str(v).ljust(widths[c]) for c, v in enumerate(vals)).rstrip()
    lines = ["Top-N crystal-cluster reach — inferential statistics",
             "(companion to topN_crystal_cluster_matrix.png)",
             f"n={n_total} complexes; {_equibind_rank_note(eq_variant)}",
             *_crystal_site_lines(),
             "=" * (sum(widths) + 2 * (len(widths) - 1)),
             _line(col_labels),
             "-" * (sum(widths) + 2 * (len(widths) - 1))]
    lines += [_line(r) for r in cell_text]
    lines += ["", _TOPN_STATS_TABLE_CAPTION]
    txt_path = out_dir / "topN_crystal_cluster_matrix_stats.txt"
    txt_path.write_text("\n".join(lines) + "\n")
    print(f"  Stats TXT: {txt_path}")
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
    axR.legend(fontsize=8, loc="upper left")   # keep lower-right clear for the stat box
    # ── paired-binary comparison across the three tools at each labelled depth ──
    #    ('reaches the crystal cluster within top-N' per complex per tool; the same
    #    `proteins` are the paired units). Annotate the top-1 omnibus; sidecar full.
    try:
        if n_total >= 8:
            depth_stats = {}
            reach_bool = {}
            for N in mark:
                cols = {t: np.array([reach.get(p, {}).get(t, np.inf) <= N
                                     for p in proteins], int) for t in T}
                reach_bool[N] = cols
                depth_stats[f"top{N}"] = su.paired_proportions(cols, labels=list(T))
            o1 = depth_stats.get(f"top{mark[0]}", {}).get("omnibus") if mark else None
            if o1:
                axL.text(0.02, 0.02,
                         f"top-{mark[0]} across tools: Cochran Q={o1['Q']:.1f}, "
                         f"{su.fmt_p(o1['p'])} {su.p_stars(o1['p'])}",
                         transform=axL.transAxes, fontsize=8, va="bottom",
                         bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.7",
                                   alpha=0.85))
            # ── co-reach ASSOCIATION at top-1: the off-diagonal's unique content ──
            #    Do two tools reach the crystal cluster on the SAME complexes more
            #    than their marginal rates predict? Signed phi / Haldane odds ratio +
            #    two-sided Fisher (vs McNemar, which only tests equal marginal RATES
            #    and ignores the joint-reach cell). AutoDock sits at a top-1 ceiling,
            #    so its pairs have a near-empty 'fails' margin → not estimable.
            cooc = None
            if mark:
                N1 = mark[0]; rb = reach_bool[N1]
                pair_assoc = {f"{pr[0]}+{pr[1]}": su.paired_2x2_association(
                    rb[pr[0]], rb[pr[1]]) for pr in pairs}
                est = [(pr, pair_assoc[f"{pr[0]}+{pr[1]}"]) for pr in pairs
                       if pair_assoc[f"{pr[0]}+{pr[1]}"]["estimable"]]
                if est:                            # Holm across the estimable pairs only
                    for (pr, aa), ph in zip(est, su.holm([a2["p"] for _, a2 in est])):
                        aa["p_holm"] = float(ph); aa["star"] = su.p_stars(ph)
                allc = su.consensus_perm_test(
                    np.column_stack([rb[t] for t in T]), n_perm=10000)
                lines = [f"Top-{N1} co-reach vs independence (n={n_total}):"]
                for pr in pairs:
                    a2 = pair_assoc[f"{pr[0]}+{pr[1]}"]
                    nm = f"{_TOOL_DISPLAY[pr[0]]}+{_TOOL_DISPLAY[pr[1]]}"
                    if a2["estimable"]:
                        lines.append(
                            f"{nm}: {a2['observed']} obs / {a2['expected']:.0f} exp  "
                            f"φ={a2['phi']:+.2f}  OR={a2['odds_ratio']:.1f}"
                            f"[{a2['or_ci'][0]:.1f}–{a2['or_ci'][1]:.1f}]  "
                            f"{su.fmt_p(a2.get('p_holm', a2['p']))} {a2.get('star', '')}")
                    else:
                        lines.append(f"{nm}: AutoDock at ceiling → not estimable")
                ka = allc["k"]
                lines.append(
                    f"All 3: {int(allc['observed'][ka])} obs / "
                    f"{allc['expected'][ka]:.0f} exp  "
                    f"(2-sided perm {su.fmt_p(allc['p_all_agree_two_sided'])})")
                axR.text(0.98, 0.03, "\n".join(lines), transform=axR.transAxes,
                         ha="right", va="bottom", fontsize=6.5,
                         bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.7",
                                   alpha=0.92))
                cooc = {"depth": int(N1),
                        "null": "independence of the two per-complex reach indicators",
                        "effect": "signed phi / Haldane OR (+ redundant, − complementary)",
                        "pairs": pair_assoc, "all_three": allc}
            _reach_stats_sink(out_dir, "topN_crystal_reach_curves", {
                "n_complexes": int(n_total),
                "unit": "per-complex boolean (best rank in crystal cluster ≤ N)",
                "depths": depth_stats,
                "cooccurrence_top1": cooc})
        else:
            axL.text(0.02, 0.02, f"n={n_total} too small — exploratory",
                     transform=axL.transAxes, fontsize=8, va="bottom", color="0.4")
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] topN reach paired test skipped: {e}")
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
        # unrounded twin for the text formatter: printing ``:.1%`` from the 4-dp
        # value double-rounds (168/303 = 0.55446 -> 0.5545 -> "55.5%"; correct is
        # 55.4%). Underscore-prefixed so the summary.json writer can drop it and
        # the JSON keeps exactly the 4-dp field it always had.
        "_precision_at_1_exact": {r: float(prec[r]) for r in rules},
        "pairwise_mcnemar_holm": pairs,
        "consensus_vs_size_mcnemar": cons_vs_size,
    }


_ABLATION_RULE_ABBR = {"size": "size", "ntools": "cons", "tight": "tght",
                       "confidence": "conf", "medoid": "medo"}


def _ablation_pairwise_matrix_text(st):
    """Render the full Holm-adjusted pairwise McNemar family (``st`` from
    :func:`_ablation_paired_stats`) as a compact lower-triangular star matrix,
    rules ordered by precision@1 (best first). Same significance-star vocabulary
    (:func:`_p_stars`) as the other stats sidecars."""
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


_ABLATION_RULE_NAME = {"size": "size (legacy count)", "ntools": "consensus (n tools)",
                       "tight": "tightest spread", "confidence": "confidence-weighted",
                       "medoid": "consensus + medoid"}


def _format_ablation_stats(st, ranking_rho, match_thr, n_total):
    """Plain-text rendering of the precision@1 ranking-ablation statistics (panel F
    of cluster_quality_metrics.png) for the *_stats.txt sidecar and the console.

    Carries every test that used to be annotated on the panel: the Spearman ranking-
    enrichment, the Cochran-Q omnibus, the per-rule precision@1, and the full Holm-
    adjusted pairwise McNemar family (as a readable list plus the compact star matrix).
    """
    nm = _ABLATION_RULE_NAME
    out = ["cluster_quality_metrics.png — precision@1 ranking-ablation statistics",
           f"n={n_total} complexes; rank-1 pocket within {match_thr:g} Å of crystal",
           *_crystal_site_lines(),
           "stars: * p<.05  ** p<.01  *** p<.001 (Holm-adjusted McNemar)",
           "=" * 72]
    if ranking_rho is not None:
        out.append(f"ranking enrichment: Spearman rho={ranking_rho:+.3f}")
    if not st:
        out.append("paired ablation tests unavailable (no-crystal dataset or too few "
                   "scored complexes).")
        return "\n".join(out)
    out.append(f"omnibus (paired): Cochran Q={st['cochran_q']} df={st['df']} "
               f"p={_fmt_p(st['p_omnibus'])} {_p_stars(st['p_omnibus'])} "
               f"(k={st['k']} rules, n={st['n_omnibus']} complexes)")
    # print from the unrounded fraction (the 4-dp ``precision_at_1`` field is
    # what the JSON carries; formatting it with ``:.1%`` double-rounds).
    prec = st.get("_precision_at_1_exact") or st.get("precision_at_1", {})
    out.append("")
    out.append("precision@1 by ranking rule (best first):")
    for r in sorted(prec, key=lambda k: prec[k], reverse=True):
        out.append(f"    {nm.get(r, r):<24} {prec[r]:.1%}")
    out.append("")
    out.append("pairwise McNemar (a beats b; Holm-adjusted across the full family):")
    for pr in st.get("pairwise_mcnemar_holm", []):
        out.append(f"    {nm.get(pr['a'], pr['a']):<22} beats "
                   f"{nm.get(pr['b'], pr['b']):<22} "
                   f"wins {pr['a_wins']}/{pr['b_wins']} (n={pr['n_pairs']})  "
                   f"p_raw={_fmt_p(pr['p_raw'])}  p_Holm={_fmt_p(pr['p_holm'])} "
                   f"{_p_stars(pr['p_holm'])}")
    cvs = st.get("consensus_vs_size_mcnemar")
    if cvs:
        out.append("")
        out.append("headline — consensus vs legacy size: consensus wins "
                   f"{cvs['consensus_wins']}, size wins {cvs['size_wins']} "
                   f"(n={cvs['n_pairs']})  p_raw={_fmt_p(cvs['p_raw'])}  "
                   f"p_Holm={_fmt_p(cvs['p_holm'])} {_p_stars(cvs['p_holm'])}")
    mtxt = _ablation_pairwise_matrix_text(st)
    if mtxt:
        out.append("")
        out.append(mtxt)
    return "\n".join(out)


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


def _partial_spearman(x, y, z):
    """First-order partial Spearman corr(x, y | z): the rank association between
    x and y after removing what each shares with the control variable z. Uses the
    standard partial-correlation formula on Spearman coefficients (no extra deps).
    Returns (partial_rho, p, components{rho_xy, rho_xz, rho_yz, n})."""
    from scipy.stats import spearmanr, t as tdist
    x = np.asarray(x, float); y = np.asarray(y, float); z = np.asarray(z, float)
    m = ~(np.isnan(x) | np.isnan(y) | np.isnan(z))
    x, y, z = x[m], y[m], z[m]
    n = len(x)
    if n < 5:
        return np.nan, np.nan, {}
    rxy = float(spearmanr(x, y)[0]); rxz = float(spearmanr(x, z)[0])
    ryz = float(spearmanr(y, z)[0])
    comps = {"rho_xy": round(rxy, 3), "rho_xz": round(rxz, 3),
             "rho_yz": round(ryz, 3), "n": int(n)}
    denom = np.sqrt(max((1 - rxz ** 2) * (1 - ryz ** 2), 0.0))
    if denom == 0:
        return np.nan, np.nan, comps
    rp = float(np.clip((rxy - rxz * ryz) / denom, -1.0, 1.0))
    df = n - 3                                          # n - 2 - (1 control var)
    if df <= 0 or abs(rp) >= 1:
        p = np.nan
    else:
        tstat = rp * np.sqrt(df / (1 - rp ** 2))
        p = float(2 * tdist.sf(abs(tstat), df))
    return rp, p, comps


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
        cs = np.array([p[2] for p in nt_radius_pairs], float)
        m = ~np.isnan(rr)
        rho, prho = spearmanr(nt[m], rr[m]) if m.sum() > 2 else (np.nan, np.nan)
        # Control for the mechanical confound: radius (max pose->center) grows with
        # the cluster's pose count, and more-consensus clusters hold more poses.
        # Partial Spearman(radius, n_tools | pose_count) removes that.
        rho_p, prho_p, comps = _partial_spearman(rr, nt, cs)
        rep["D_radius_trend"] = {
            "test": "Jonckheere-Terpstra (ordered 1<2<3 tools)",
            "JT": round(JT, 1), "z": round(zD, 3), "p": float(pD),
            "spearman_rho": round(float(rho), 3), "spearman_p": float(prho),
            "n": {f"{k}_tools": int(len(gD[k - 1])) for k in (1, 2, 3)},
            "medians": {f"{k}_tools": round(float(np.median(gD[k - 1])), 3)
                        for k in (1, 2, 3)},
            "posecount_control": {
                "method": "partial Spearman(radius, n_tools | crystal-cluster "
                          "pose count)",
                "partial_rho": (round(rho_p, 3) if rho_p == rho_p else None),
                "partial_p": prho_p,
                "rho_radius_vs_posecount": comps.get("rho_xz"),
                "rho_ntools_vs_posecount": comps.get("rho_yz"),
                "n": comps.get("n"),
                "reading": "if partial_rho collapses toward 0 / loses "
                           "significance, the raw trend is a pose-count artifact; "
                           "if it survives, agreement widens the site beyond mere "
                           "pose count."},
            "caveat": "radius = max(pose->center) grows mechanically with the "
                      "number of poses in the cluster, and more-consensus "
                      "clusters tend to hold more poses; see posecount_control "
                      "for the confound-adjusted trend."}
    return rep


def _format_homogeneity_stats(rep):
    """Render the homogeneity statistical tests as a plain-text block (shared by
    the console printout and the crystal_cluster_homogeneity_stats.txt sidecar)."""
    out = ["  ── crystal_cluster_homogeneity.png statistics ──"]

    def line(s=""):
        out.append("    " + s)
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
        pc = d.get("posecount_control", {})
        if pc:
            line(f"      pose-count control: partial rho={pc.get('partial_rho')} "
                 f"(p={_fmt_p(pc.get('partial_p'))} {_p_stars(pc.get('partial_p'))}) "
                 f"| radius~pose# rho={pc.get('rho_radius_vs_posecount')}, "
                 f"n_tools~pose# rho={pc.get('rho_ntools_vs_posecount')}")
            line(f"      → {pc.get('reading')}")
        line(f"      CAVEAT: {d['caveat']}")
    return "\n".join(out)


def _print_homogeneity_stats(rep):
    print("\n" + _format_homogeneity_stats(rep))


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
    cluster_size = {}                        # #poses in the crystal-closest pocket
    for r in ok:
        if not r.get("has_crystal"):
            continue
        cr = r.get("_crystal"); pk = r.get("_pockets")
        if cr is None or not pk:
            continue
        cp = min(pk, key=lambda p: _cdist(p["center"], cr))   # nearest to the nearest copy
        radius[r["protein"]] = float(cp.get("radius", np.nan))
        cluster_size[r["protein"]] = int(cp.get("size", 0))
    rad_by_n = {k: [] for k in (1, 2, 3)}
    # (n_tools, radius, crystal-cluster pose count) per complex for the trend +
    # its pose-count-controlled partial correlation.
    nt_radius_pairs = []
    for prot in proteins:
        k = int(n_tools_per.get(prot, 0))
        rv = radius.get(prot, np.nan)
        if k in rad_by_n and rv == rv:
            rad_by_n[k].append(rv)
        nt_radius_pairs.append((k, rv, cluster_size.get(prot, np.nan)))

    stats_rep = None
    if stats:
        try:
            stats_rep = _homogeneity_stats(comp, reach_R, internal, internal_prot,
                                           rad_by_n, nt_radius_pairs, T, _TOOL_DISPLAY)
        except Exception as e:                                  # never break the figure
            warnings.warn(f"homogeneity stats failed: {e}")
            stats_rep = None

    rng = np.random.default_rng(42)
    # consensus-count palette for panels B & D (1 / 2 / all-3 tools): a purple
    # sequential ramp (light→dark = more consensus), deliberately distinct from the
    # per-tool colours (AutoDock blue / DiffDock green / EquiBind red) so a bar/box is
    # never confused with a tool.
    tri_col = ["#CBB9E0", "#9B72C0", "#6A3D9A"]
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
    axA.set_title("Pose contribution per tool")  # tests → *_stats.txt sidecar
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
    axB.set_title("Consensus richness of the crystal cluster")  # tests → *_stats.txt
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
        # n= over each box (blended transform: x in data, y in axes fraction) so the
        # label never falls below the axis and gets clipped
        axC.text(i, 0.98, f"n={len(yy)}", transform=axC.get_xaxis_transform(),
                 ha="center", va="top", fontsize=8, color=TOOL_COLORS[t])
    axC.set_xticks([1, 2, 3]); axC.set_xticklabels([_TOOL_DISPLAY.get(t, t) for t in T])
    axC.set_ylabel("Std. of a tool's pose distances to the\ncrystal centroid, within the cluster (Å)")
    axC.set_title("Internal tightness of each tool's own poses")  # tests → *_stats.txt
    _all_c = [v for t in T for v in internal[t]]
    axC.set_ylim(-0.02, (max(_all_c) * 1.14 if _all_c else None))   # headroom for n= labels
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
            # offset to the right of the box: the white mean diamond sits on the
            # box centre and would otherwise print through the median label
            axD.text(i + 1 + 0.34, med, f"{med:.2f} Å", ha="left",
                     va="center", fontsize=9)
    axD.set_xticks(xs); axD.set_xticklabels(["1 tool", "2 tools", "all 3 tools"])
    axD.set_ylabel("Crystal-cluster radius (max pose->center, Å)")
    axD.set_xlabel("Number of tools contributing to the crystal cluster")
    axD.set_title("Spatial spread vs consensus richness")  # tests → *_stats.txt
    axD.grid(axis="y", alpha=0.25); axD.set_axisbelow(True)

    _label_panels(axes)
    fig.suptitle("How homogeneous is the crystal-closest cluster?  "
                 "Composition, internal tightness and spatial spread", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    p = out_dir / "crystal_cluster_homogeneity.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    if stats_rep:
        (out_dir / "crystal_cluster_homogeneity_stats.json").write_text(
            json.dumps(stats_rep, indent=2, default=str))
        header = ("crystal_cluster_homogeneity.png — statistical tests\n"
                  f"n={n_total} complexes; top-15 poses per tool; "
                  f"{_equibind_rank_note(eq_variant)}\n"
                  + _crystal_site_header()
                  + "stars: * p<.05  ** p<.01  *** p<.001 (Holm-adjusted)\n"
                  + "=" * 72 + "\n")
        txt_path = out_dir / "crystal_cluster_homogeneity_stats.txt"
        txt_path.write_text(header + _format_homogeneity_stats(stats_rep) + "\n")
        _print_homogeneity_stats(stats_rep)
        print(f"  Stats TXT:  {txt_path}")
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

    fig, ax = plt.subplots(3, 2, figsize=(11, 14.25))

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
    _hist(ax[1, 0], _col("db_score"),
          "Davies-Bouldin index (k≥2 only)\nlower = better",
          "Davies-Bouldin index", "#C44E52")

    # (D) compactness vs separation (both defined at k=1)
    comp = pd.to_numeric(df.get("compactness"), errors="coerce")
    sep = pd.to_numeric(df.get("separation"), errors="coerce")
    m = comp.notna() & sep.notna()
    a = ax[1, 1]
    if m.any():
        a.scatter(comp[m], sep[m], s=18, alpha=0.5, color="#8172B3", edgecolors="none")
        a.axhline(float(match_thr), color="red", ls="--", lw=0.8,
                  label=f"match threshold {match_thr:g} Å")
        a.set_ylim(0, 30)
        a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "no data", ha="center", va="center",
               transform=a.transAxes, color="0.5")
    a.set_title("Site compactness vs separation\ntight (low x) and far-apart (high y) = good")
    a.set_xlabel("Within-site spread (Å, median distance to centre)")
    a.set_ylabel("Nearest inter-site distance (Å; higher = better)")
    a.grid(alpha=0.25); a.set_axisbelow(True)

    # (E) bootstrap stability
    _hist(ax[2, 0], _col("boot_stability"),
          f"Bootstrap stability ({stability_boot} resamples)\nreproducibility of the partition",
          "Mean best-cluster Jaccard overlap", "#CCB974",
          ref=[(0.75, "stable ≥0.75", "#2ca02c"),
               (0.5, "dissolved <0.5", "#d62728")])

    # (F) precision@1 ranking ablation
    a = ax[2, 1]
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
    # Spearman enrichment, Cochran-Q omnibus and the Holm-adjusted pairwise McNemar
    # family are moved off the panel → cluster_quality_metrics_stats.txt.
    st = _ablation_paired_stats(df)
    a.set_title("Precision@1 — rank-1 pocket within "
                f"{match_thr:g} Å of crystal\n(default = consensus)")
    a.set_xlabel("Precision@1 (fraction of complexes)")
    a.grid(axis="x", alpha=0.25); a.set_axisbelow(True)

    _label_panels(ax)
    fig.suptitle("Cluster-quality assessment — internal validity, stability, and "
                 f"ranking accuracy  (n={len(df)} complexes)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = out_dir / "cluster_quality_metrics.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    if st or ranking_rho is not None:
        txt = _format_ablation_stats(st, ranking_rho, match_thr, len(df))
        txt_path = out_dir / "cluster_quality_metrics_stats.txt"
        txt_path.write_text(txt + "\n")
        print("\n" + txt)
        print(f"  Stats TXT:  {txt_path}")
    return p


# ════════════════════════════════════════════════════════════════════════
# Driver
# ════════════════════════════════════════════════════════════════════════

def _fig_descriptor_quality(df_complex, features_csv, out_dir,
                            per_pose_csv=None, dd_variant="diffdock",
                            eq_variant="equibind_unguided_smina",
                            ad_variant="autodock"):
    """Which ligand types dock well? Pose accuracy (oracle RMSD to crystal) vs the
    ligand's Lipinski Ro5 compliance, physicochemical descriptors, and PCA space.

    Every success measure is reported for TWO categories: ``near-native`` (a pose
    RMSD <= 2 A, any validity) and ``PB-valid & <= 2 A`` (the near-native pose is
    also PoseBusters-valid) — the gap is the fraction of near-native poses lost to
    physical invalidity. Both oracles are recomputed from the raw per-pose CSV (over
    the chosen AutoDock/DiffDock/EquiBind variants) so the comparison holds
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
    tool_method = {"autodock": ad_variant, "diffdock": dd_variant, "equibind": eq_variant}

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
    # Ro5 pass vs fail are INDEPENDENT groups -> Mann-Whitney U + Cliff's delta.
    try:
        a = pd.to_numeric(grp[0][1]["oracle_all"], errors="coerce").to_numpy()
        b = pd.to_numeric(grp[1][1]["oracle_all"], errors="coerce").to_numpy()
        mw = su.mannwhitney_cliffs(a, b)
        if mw is not None and min(mw["n_a"], mw["n_b"]) >= 3:
            lo, hi = mw["delta_ci"]
            ax.text(0.5, 0.98,
                    f"Mann-Whitney {su.fmt_p(mw['p'])} {su.p_stars(mw['p'])} · "
                    f"Cliff's δ={mw['cliffs_delta']:.2f} [{lo:.2f}, {hi:.2f}]",
                    transform=ax.transAxes, ha="center", va="top", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7",
                              alpha=0.85))
            _reach_stats_sink(out_dir, "descriptor_accuracy_vs_ro5", {
                "measure": "oracle_all_rmsd_A", "design": "independent (Ro5 pass vs fail)",
                "test": "Mann-Whitney U + Cliff's delta",
                "n_pass": mw["n_a"], "n_fail": mw["n_b"],
                "median_pass": mw["medians"][0], "median_fail": mw["medians"][1],
                "U": mw["U"], "p": mw["p"], "cliffs_delta": mw["cliffs_delta"],
                "delta_ci": [lo, hi]})
        else:
            ax.text(0.5, 0.98, "n too small — exploratory", transform=ax.transAxes,
                    ha="center", va="top", fontsize=8, color="0.4")
    except Exception as e:                                     # pragma: no cover
        print(f"  [stats] Ro5 Mann-Whitney skipped: {e}")
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
    rhos = []          # [feature, rho, p_raw, n]
    for f in feats:
        v = pd.to_numeric(d[f], errors="coerce")
        m = v.notna() & d["best_oracle_rmsd"].notna()
        if m.sum() > 10:
            try:
                rho, prho, nn = su.spearman(v[m].to_numpy(),
                                            d["best_oracle_rmsd"][m].to_numpy())
            except Exception:
                rho, prho, nn = spearmanr(v[m], d["best_oracle_rmsd"][m])[0], np.nan, int(m.sum())
            rhos.append([f, float(rho), float(prho), int(nn)])
    if rhos:
        # BH-FDR across the descriptor family; stars from the q-values.
        try:
            qs = su.bh_fdr([r[2] for r in rhos])
        except Exception as e:                                # pragma: no cover
            print(f"  [stats] descriptor BH-FDR skipped: {e}")
            qs = [np.nan] * len(rhos)
        for r, q in zip(rhos, qs):
            r.append(float(q))                                # r[4] = q_bh
        rhos.sort(key=lambda x: x[1])
        vals = [r[1] for r in rhos]
        fig, ax = plt.subplots(figsize=(6.8, 5.4))
        ax.barh(range(len(vals)), vals, color=["#C44E52" if r > 0 else "#4C72B0" for r in vals])
        ax.set_yticks(range(len(vals))); ax.set_yticklabels([nice(r[0]) for r in rhos], fontsize=8)
        ax.axvline(0, color="k", lw=0.8)
        try:
            for i, r in enumerate(rhos):
                star = su.p_stars(r[4]) if len(r) > 4 else ""
                if star and star != "ns":
                    ax.text(r[1] + (0.012 if r[1] >= 0 else -0.012), i, star,
                            va="center", ha="left" if r[1] >= 0 else "right",
                            fontsize=9, fontweight="bold")
        except Exception as e:                                # pragma: no cover
            print(f"  [stats] descriptor stars skipped: {e}")
        ax.set_title("Which properties track pose error?\n"
                     "(Spearman ρ vs RMSD; red = worse, blue = better; * = BH-FDR q<0.05)")
        ax.set_xlabel("Spearman correlation with best oracle RMSD")
        ax.grid(alpha=0.2, axis="x"); ax.set_axisbelow(True)
        _save(fig, "descriptor_property_correlations.png")
        try:
            _reach_stats_sink(out_dir, "descriptor_property_correlations", {
                "measure": "best_oracle_rmsd", "n_complexes": int(len(d)),
                "multiplicity": f"BH-FDR across {len(rhos)} descriptor rows",
                "rows": [{"feature": r[0], "rho": r[1], "p_raw": r[2], "n": r[3],
                          "q_bh": r[4] if len(r) > 4 else None} for r in rhos]})
        except Exception as e:                                # pragma: no cover
            print(f"  [stats] descriptor correlation sidecar skipped: {e}")

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
_CACHE_SCHEMA = 6   # 6: load_heavy_atom_mol now strips explicit Hs (Phase 0.3); schema-5 caches carry H-inclusive centroids


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


def _analysis_signature(args, eq_variant: str, dd_variant: str,
                        ad_variant: str = "autodock") -> dict:
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
        "eq_variant": eq_variant, "dd_variant": dd_variant, "ad_variant": ad_variant,
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
        # Method exclusion changes which poses enter the analysis, so it MUST be
        # part of the key. Without it a cached run silently returns the unfiltered
        # analysis and the exclusion looks like it did nothing.
        "exclude_methods": mf.resolve_patterns(args),
        # D5 crystal-site convention. Keyed ONLY when it departs from the reference
        # so every existing reference cache keeps its exact signature; an any-copy
        # cache can therefore never be mistaken for a reference one or vice versa.
        **({"crystal_copies": args.crystal_copies}
           if getattr(args, "crystal_copies", "reference") != "reference" else {}),
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


# ════════════════════════════════════════════════════════════════════════
# Rank-1 pose quality cross-tab — near-native (RMSD < 2 Å) × crystal cluster
# ════════════════════════════════════════════════════════════════════════
# Produced in BOTH pose-set flavours on EVERY run, independent of --pb-valid-only,
# so a single invocation always yields the all-poses AND the PB-valid table:
#   • all_poses — rank-1 = each tool's native top pick (Vina mode 1 / DiffDock
#     confidence-1 / EquiBind best-affinity), clustered over ALL poses; INCLUDES
#     rank-1 poses that fail PoseBusters.
#   • pb_valid  — rank-1 = each tool's best PoseBusters-valid pose (survivors
#     re-ranked), clustered over the PB-valid poses only.
# Each flavour uses its OWN clustering (the crystal-closest cluster is defined on
# that pose set), so the two tables reproduce what a plain vs --pb-valid-only run
# each give. The flavour matching the main run reuses its analysis; the other is
# computed with a lite clustering pass (bootstrap + placement-mode analysis off).

def _compute_ok_lite(csv, eq_variant, ids, pb_valid_only, dd_variant, args,
                     cents_seed=None, ad_variant="autodock"):
    """Per-complex clustering pass for the rank-1 cross-tab with the expensive
    bootstrap-stability + placement-mode passes OFF (n_boot=0, do_placement=False):
    the cross-tab only needs each pose's rank, RMSD-to-crystal and crystal-cluster
    membership. ``cents_seed`` reuses centroids the main run already extracted; only
    poses missing from it are (re)loaded, so when the main run is all-poses the
    PB-valid subset pass needs no new centroid extraction. Same pose loading and
    per-complex outlier removal as the main analysis loop, so counts match a real
    run of the corresponding mode."""
    df = load_poses(csv, eq_variant, ids, pb_valid_only, dd_variant, ad_variant)
    if df.empty:
        return []
    complexes = sorted(df["protein"].unique())
    if args.limit:
        complexes = complexes[:args.limit]
        df = df[df["protein"].isin(complexes)].copy()
    cents = dict(cents_seed) if cents_seed else {}
    missing = sorted(set(df["pose_file"].unique()) - set(cents))
    if missing:
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
            for path, cx, cy, cz, na in ex.map(_extract_centroid, missing, chunksize=64):
                if cx is not None:
                    cents[path] = (cx, cy, cz)
    keep_idx = []                                   # same per-complex outlier removal as main
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
    results = []
    for cid in complexes:
        sub = df[df["protein"] == cid]
        if sub.empty:
            continue
        crystal = crystal_centroid(Path(args.benchmark_dir), cid,
                                   getattr(args, "crystal_copies", "reference"))
        stem = f"{cid}_protein"
        fp_out = Path(args.fpocket_dir) / f"{stem}_out"
        pr_csv = Path(args.p2rank_dir) / f"{stem}.pdb_predictions.csv"
        fp_pockets = parse_fpocket(fp_out)[:args.top_n_pockets] if fp_out.exists() else []
        pr_pockets = parse_p2rank(pr_csv)[:args.top_n_pockets] if pr_csv.exists() else []
        res = analyze_complex(cid, sub, cents, crystal, fp_pockets, pr_pockets,
                              args.match_thr, False, args.mode_rmsd_thr,
                              args.rmsd_pose_cap, args.site_cluster,
                              args.pocket_radius, args.rank_by, 0)
        results.append(res)
    return [r for r in results if not r.get("skipped")]


def _rank1_quality_counts(ok):
    """Per-tool contingency of the rank-1 pose over (near-native RMSD < 2 Å) ×
    (in the crystal-closest cluster), read from the rank-1 rows of ``ok`` (the
    ``_per_rank`` payload). ``ic``/``rm`` may be NaN (no crystal / uncomputable
    RMSD); a NaN counts as a failure of that axis."""
    agg = {t: dict(n=0, both=0, near_only=0, cluster_only=0, neither=0, nan_rmsd=0)
           for t in TOOLS}
    for r in ok:
        if not r.get("has_crystal"):
            continue
        for (_cid, t, rp, _cd, rm, ic, _cl) in r.get("_per_rank", []):
            if rp != 1 or t not in agg:
                continue
            a = agg[t]
            a["n"] += 1
            has_rmsd = (rm == rm)
            near = has_rmsd and (rm < 2.0)
            if not has_rmsd:
                a["nan_rmsd"] += 1
            incl = bool(ic) if not (isinstance(ic, float) and ic != ic) else False
            if near and incl:
                a["both"] += 1
            elif near and not incl:
                a["near_only"] += 1
            elif (not near) and incl:
                a["cluster_only"] += 1
            else:
                a["neither"] += 1
    return agg


_RANK1_QUALITY_MODES = [
    ("all_poses", False,
     "all rank-1 poses — native top pick per tool (incl. PB-invalid); all-poses clustering"),
    ("pb_valid", True,
     "PB-valid rank-1 poses — best pose passing PoseBusters (re-ranked); PB-valid clustering"),
]


def _fig_rank1_quality_crosstab(df_ct, out_dir, eq_variant=None, dd_variant="diffdock"):
    """Two stacked tables (all_poses / pb_valid) of the per-tool rank-1 quality
    contingency for the notebook + thesis. Purely descriptive counts."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if df_ct.empty:
        return None
    cols = [("tool_display", "tool"), ("n", "N"),
            ("near_and_in_cluster", "<2Å & in\ncluster"),
            ("fail_near_andor_cluster", "fail ≥1\n(and/or)"),
            ("not_near_native", "not <2Å"),
            ("not_in_crystal_cluster", "not in\ncluster"),
            ("fail_both", "complete\nmiss")]
    modes = [("all_poses", "All rank-1 poses (native top pick — incl. PB-invalid)"),
             ("pb_valid", "PB-valid rank-1 poses (best PoseBusters-valid pose)")]
    fig, axes = plt.subplots(len(modes), 1, figsize=(9.2, 4.6), constrained_layout=True)
    for ax, (key, sub_title) in zip(np.atleast_1d(axes), modes):
        ax.axis("off")
        ax.set_title(sub_title, fontsize=9.5, fontweight="bold", loc="left")
        d = df_ct[df_ct["pose_set"] == key]
        if d.empty:
            ax.text(0.5, 0.5, "(no data)", ha="center")
            continue
        cell = [[str(r[c]) for c, _ in cols] for _, r in d.iterrows()]
        tbl = ax.table(cellText=cell, colLabels=[h for _, h in cols],
                       cellLoc="center", loc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8.5)
        tbl.scale(1, 1.6)
        for (row, _c), cellobj in tbl.get_celld().items():
            if row == 0:
                cellobj.set_facecolor("#e8eef4")
                cellobj.set_text_props(fontweight="bold")
    fig.suptitle("Rank-1 pose quality — near-native (RMSD < 2 Å) × crystal cluster, per tool\n"
                 f"Benchmark; {_equibind_rank_note(eq_variant)}",
                 fontsize=11, fontweight="bold")
    p = out_dir / "rank1_quality_crosstab.png"
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


def _rank1_quality_crosstabs(csv, eq_variant, ids, dd_variant, args, out_dir,
                             ok_main, cents_seed=None, render_png=True,
                             ad_variant="autodock"):
    """Write the per-tool rank-1 quality cross-tab in BOTH pose-set flavours
    (all_poses + pb_valid): a CSV, a human-readable TXT and (unless render_png is
    False) a PNG. The flavour whose ``--pb-valid-only`` state matches the main run
    reuses ``ok_main``; the other is computed with a lite clustering pass."""
    records, txt_blocks = [], []
    for key, pb_only, desc in _RANK1_QUALITY_MODES:
        if bool(pb_only) == bool(args.pb_valid_only) and ok_main:
            ok_m = ok_main
        else:
            ok_m = _compute_ok_lite(csv, eq_variant, ids, pb_only, dd_variant, args,
                                    cents_seed=cents_seed, ad_variant=ad_variant)
        agg = _rank1_quality_counts(ok_m)
        head = (f"{'tool':22s} {'N':>4} | {'<2Å&inCl':>8} {'<2Åonly':>7} {'inClonly':>8} "
                f"{'neither':>7} | {'not<2Å':>6} {'notInCl':>7} {'fail≥1':>6} {'failboth':>8}")
        lines = [f"[{key}]  {desc}", head, "-" * len(head)]
        for t in TOOLS:
            a = agg[t]
            n = a["n"]
            if n == 0:
                continue
            both = a["both"]
            not_near = a["cluster_only"] + a["neither"]
            not_cl = a["near_only"] + a["neither"]
            fail_any = n - both
            records.append(dict(
                pose_set=key, tool=t, tool_display=_TOOL_DISPLAY[t], n=n,
                near_and_in_cluster=both, near_not_in_cluster=a["near_only"],
                far_but_in_cluster=a["cluster_only"], neither=a["neither"],
                not_near_native=not_near, not_in_crystal_cluster=not_cl,
                fail_near_andor_cluster=fail_any, fail_both=a["neither"],
                pct_fully_correct=round(100 * both / n, 1),
                pct_fail_any=round(100 * fail_any / n, 1),
                nan_rmsd=a["nan_rmsd"]))
            lines.append(f"{_TOOL_DISPLAY[t]:22s} {n:>4} | {both:>8} {a['near_only']:>7} "
                         f"{a['cluster_only']:>8} {a['neither']:>7} | {not_near:>6} "
                         f"{not_cl:>7} {fail_any:>6} {a['neither']:>8}")
        txt_blocks.append("\n".join(lines))

    df_ct = pd.DataFrame(records)
    csv_path = out_dir / "rank1_quality_crosstab.csv"
    df_ct.to_csv(csv_path, index=False)

    header = [
        "Rank-1 pose quality — near-native (RMSD < 2 Å) × crystal-closest cluster, per tool",
        f"Benchmark; {_equibind_rank_note(eq_variant)}; DiffDock variant = {dd_variant}; "
        f"AutoDock variant = {ad_variant}.",
        *_crystal_site_lines(),
        "Both flavours are written on EVERY run, independent of --pb-valid-only:",
        "  • all_poses — rank-1 = each tool's native top pick (AutoDock mode 1, i.e. the",
        "    optimizer's re-ranked mode 1 for the gnina/smina variants / DiffDock",
        "    confidence-1 / EquiBind best-affinity); clustering over ALL poses; includes",
        "    rank-1 poses that FAIL PoseBusters.",
        "  • pb_valid  — rank-1 = each tool's best PoseBusters-valid pose (valid poses",
        "    re-ranked); clustering over the PB-valid poses only.",
        "=" * 78,
    ]
    key_lines = [
        "",
        "Column key:",
        "  N        = complexes with a rank-1 pose for that tool",
        "  <2Å&inCl = rank-1 near-native AND in the crystal cluster (fully correct)",
        "  <2Åonly  = near-native but NOT in the crystal cluster (rare edge case)",
        "  inClonly = in the crystal cluster but RMSD ≥ 2 Å (right pocket, wrong pose)",
        "  neither  = both fail (complete miss)",
        "  not<2Å   = RMSD ≥ 2 Å            (= inClonly + neither)",
        "  notInCl  = outside the crystal cluster (= <2Åonly + neither)",
        "  fail≥1   = not < 2 Å AND/OR not in cluster (= N − <2Å&inCl)  ← the 'and/or' count",
        "  failboth = neither near-native nor in cluster",
    ]
    txt_path = out_dir / "rank1_quality_crosstab.txt"
    txt_path.write_text("\n".join(header) + "\n\n" + "\n\n".join(txt_blocks) + "\n"
                        + "\n".join(key_lines) + "\n")

    png_path = _fig_rank1_quality_crosstab(df_ct, out_dir, eq_variant, dd_variant) \
        if render_png else None
    print(f"  Rank-1 quality cross-tab (all_poses + pb_valid) → {csv_path.name}, "
          f"{txt_path.name}" + (f", {png_path.name}" if png_path else ""))
    return csv_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")
    ap.add_argument("--benchmark-dir", default="Data/PoseBuster Benchmark Set")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--p2rank-dir", default="pocket_results/p2rank_results")
    ap.add_argument("--ids-file", default=None)
    # ── AutoDock tool slot: choose the variant by component, or by full name ──
    ap.add_argument("--autodock-variant", default=None,
                    help="Full AutoDock method string ('autodock' | 'autodock_gnina' "
                         "| 'autodock_gnina_refinement' | 'autodock_vinardo' | ...) or "
                         "a bare optimizer keyword. Overrides --autodock-scoring/"
                         "--autodock-optimize. Default: compose from the components.")
    ap.add_argument("--autodock-scoring", choices=_AD_SCORINGS, default="vina",
                    help="AutoDock scoring function: 'vina' -> autodock* | 'vinardo' "
                         "-> autodock_vinardo* (default: vina).")
    ap.add_argument("--autodock-optimize", choices=_AD_OPTIMIZERS, default="none",
                    help="Post-docking optimizer applied to the Vina poses: 'none' "
                         "(raw Vina modes) | 'gnina' (minimise + CNN re-rank) | "
                         "'gnina_refinement' (CNN gradients drive the geometry) | "
                         "'smina' (default: none).")
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
    ap.add_argument("--allow-mixed-convention", action="store_true",
                    help="Proceed when --crystal-copies and the per-pose table's reference_convention "
                         "disagree (default: abort, because centroid distances and the RMSD column "
                         "would then refer to different copies).")
    ap.add_argument("--crystal-copies", choices=_CRYSTAL_COPIES_CHOICES, default="reference",
                    help="Crystal-site convention (plan v2 D5). 'reference' (default) "
                         "measures every crystal distance to the single <ID>_ligand.sdf "
                         "instance — today's outputs byte for byte. 'any' loads every "
                         "record of <ID>_ligands.sdf (record 0 = reference) on the same "
                         "heavy-atom basis as the poses; every crystal distance becomes "
                         "the minimum over copies, a cluster is correct through ANY copy "
                         "(the cluster nearest to a copy counts when within --match-thr "
                         "of it), the crystal-closest cluster is the one nearest to the "
                         "nearest copy, and a cluster nearest to two copies counts once. "
                         "Part of the cache signature; stated "
                         "in every stats *.txt header and in summary.json under 'any'.")
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
    mf.add_method_filter_args(ap)
    args = ap.parse_args(argv)
    global _CRYSTAL_COPIES, _TABLE_CONVENTION
    _CRYSTAL_COPIES = args.crystal_copies      # read by the stats-sidecar headers
    if _CRYSTAL_COPIES == "any":
        print(_CRYSTAL_ANY_RULE)

    csv = Path(args.per_pose_csv)
    if not csv.exists():
        print(f"per-pose CSV not found: {csv}")
        return 1
    # The RMSD-derived quantities of this report read the hub's ``rmsd`` column as is,
    # so the crystal-site convention chosen here must match the convention that
    # table was scored under; otherwise centroid distances and RMSDs refer to
    # different copies. Tables written before hub schema 6 carry no column and are
    # single-instance by construction.
    _hdr = pd.read_csv(csv, nrows=0).columns
    if "reference_convention" in _hdr:
        _vals = pd.read_csv(csv, usecols=["reference_convention"], low_memory=False)["reference_convention"].dropna().unique().tolist()
        _TABLE_CONVENTION = _vals[0] if len(_vals) == 1 else ",".join(map(str, sorted(_vals)))
    else:
        _TABLE_CONVENTION = None
    _table_is_nearest = _TABLE_CONVENTION == "nearest"
    if (_CRYSTAL_COPIES == "any") != _table_is_nearest:
        msg = (f"MIXED CONVENTION: --crystal-copies {_CRYSTAL_COPIES} but the per-pose table's "
               f"reference_convention is {_TABLE_CONVENTION or 'absent (single instance)'}; centroid "
               "distances and the RMSD column would refer to different copies. Use the matching "
               "table or pass --allow-mixed-convention.")
        if not args.allow_mixed_convention:
            print("ERROR: " + msg)
            return 2
        print("WARNING: " + msg)
    ids = _load_ids(Path(args.ids_file)) if args.ids_file else None
    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve which docking-mode variant fills each tool slot (exits with a helpful
    # list of what's available if the selection can't be matched in the CSV).
    cat = _variant_catalog(csv)
    eq_variant = resolve_equibind_variant(
        cat, args.equibind_variant, args.equibind_pocket,
        args.equibind_refine, args.equibind_clamp)
    dd_variant = resolve_diffdock_variant(cat, args.diffdock_variant, args.diffdock_refine)
    ad_variant = resolve_autodock_variant(cat, args.autodock_variant,
                                          args.autodock_scoring, args.autodock_optimize)
    # Figures label every slot by the variant actually clustered, so a plot never
    # says "AutoDock Vina" while showing gnina-reranked poses — and never names one
    # slot's refiner while leaving the other two bare, which is what made AutoDock
    # look like the only optimised arm on the cluster axes.
    _TOOL_DISPLAY["autodock"] = _autodock_display(ad_variant)
    _TOOL_DISPLAY["diffdock"] = _diffdock_display(dd_variant)
    _TOOL_DISPLAY["equibind"] = _equibind_display(eq_variant)

    # ── reuse the cached per-complex analysis when the inputs are unchanged ──
    sig = _analysis_signature(args, eq_variant, dd_variant, ad_variant)
    cache_path = out_dir / "analysis_cache.pkl"
    main_cents: Dict[str, Tuple[float, float, float]] = {}   # reused by the rank-1 cross-tab
    ok = None if args.force else _load_analysis_cache(cache_path, sig)
    if ok is not None:
        print(f"Reusing cached analysis: {len(ok)} complexes (inputs unchanged; "
              f"skipping re-analysis — use --force to recompute). "
              f"autodock={ad_variant} | equibind={eq_variant} | diffdock={dd_variant}")
    else:
        df = load_poses(csv, eq_variant, ids, args.pb_valid_only, dd_variant, ad_variant)
        df = mf.apply_method_filter(df, "method", args, label="pose-cluster",
                                    out_dir=out_dir)
        if df.empty:
            print("No poses after filtering (check --equibind-* / --diffdock-* / --ids-file).")
            return 1
        complexes = sorted(df["protein"].unique())
        if args.limit:
            complexes = complexes[:args.limit]
            df = df[df["protein"].isin(complexes)].copy()
        print(f"Complexes: {len(complexes)} | poses: {len(df)} | "
              f"tools: {sorted(df['tool'].unique())} | autodock={ad_variant} | "
              f"equibind={eq_variant} | diffdock={dd_variant}")

        # ── parallel centroid extraction ─────────────────────────────────
        files = sorted(df["pose_file"].unique())
        cents: Dict[str, Tuple[float, float, float]] = {}
        print(f"Extracting centroids for {len(files)} pose SDFs ({args.workers} workers)...")
        with ProcessPoolExecutor(max_workers=max(1, args.workers)) as ex:
            for path, cx, cy, cz, na in ex.map(_extract_centroid, files, chunksize=64):
                if cx is not None:
                    cents[path] = (cx, cy, cz)
        print(f"  loaded {len(cents)}/{len(files)} centroids")
        main_cents = cents          # seed the rank-1 cross-tab so it skips re-extraction

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
            crystal = crystal_centroid(Path(args.benchmark_dir), cid, args.crystal_copies)
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
                "dist_to_crystal": round(_cdist(p["center"], cr), 3) if cr is not None else np.nan,
                "medoid_dist_to_crystal": (round(_cdist(p.get("medoid_center", p["center"]), cr), 3)
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

    # ── rank-1 pose quality cross-tab — BOTH pose-set flavours, every run ──
    # (all_poses + pb_valid, near-native × crystal-cluster), independent of the
    # --pb-valid-only flag; wrapped so an add-on failure never breaks the report.
    try:
        _rank1_quality_crosstabs(csv, eq_variant, ids, dd_variant, args, out_dir,
                                 ok_main=ok, cents_seed=main_cents,
                                 render_png=not args.no_plot, ad_variant=ad_variant)
    except Exception as e:                                   # pragma: no cover
        print(f"  [rank1-quality] cross-tab skipped: {e}")

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
        modal_frac = pd.to_numeric(dc.get(f"{t}_top5_modal_frac"), errors="coerce").dropna() \
            if f"{t}_top5_modal_frac" in dc else pd.Series(dtype=float)
        nc = pd.to_numeric(dc.get(f"{t}_top5_n_clusters"), errors="coerce").dropna() \
            if f"{t}_top5_n_clusters" in dc else pd.Series(dtype=float)
        br = pd.to_numeric(dc.get(f"{t}_top5_best_rank"), errors="coerce").dropna() \
            if f"{t}_top5_best_rank" in dc else pd.Series(dtype=float)
        if len(modal_frac):
            print(f"  {t:<10}{modal_frac.mean():>16.2f}{nc.mean():>16.2f}"
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
            rr.append(p["rank"]); ddc.append(_cdist(p["center"], cr))
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
        "autodock_variant": ad_variant,
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
        # drop the underscore-prefixed formatter-only twin so the JSON schema is unchanged
        "ranking_ablation_significance": ({k: v for k, v in ablation_sig.items()
                                           if not str(k).startswith("_")}
                                          if ablation_sig is not None else None),
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
    if _CRYSTAL_COPIES == "any":
        # D5: name the convention in the machine-readable summary too (absent under
        # the reference convention so that summary.json stays byte-identical there).
        summary["crystal_copies"] = "any"
        summary["crystal_copies_rule"] = _CRYSTAL_ANY_RULE
        summary["per_pose_table_convention"] = _TABLE_CONVENTION
        summary["n_multi_copy_complexes"] = (
            int((pd.to_numeric(df_complex.get("crystal_n_copies"), errors="coerce") > 1).sum())
            if "crystal_n_copies" in df_complex else None)
        summary["n_complexes_no_correct_cluster"] = (
            int((dcx["correct_cluster_is_hit"].map(_TRUTH_MAP) == False).sum())   # noqa: E712
            if "correct_cluster_is_hit" in dcx else None)
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
                                        dd_variant=dd_variant, eq_variant=eq_variant,
                                        ad_variant=ad_variant)]
        for f in figs:
            if not f:
                continue
            for pth in (f if isinstance(f, (list, tuple)) else [f]):
                print(f"  Figure: {pth}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
