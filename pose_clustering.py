# %% [markdown]
# # Pose Clustering v3 — Centroid, RMSD, and Weighted Hybrid
#
# Three distance modes for clustering docked poses:
#
# | Mode          | Distance metric                                      |
# |---------------|------------------------------------------------------|
# | `centroid`    | Euclidean distance between heavy-atom centroids      |
# | `rmsd`        | Heavy-atom RMSD (best-fit via RDKit, fallback atom)  |
# | `hybrid`      | `w * norm_rmsd + (1-w) * norm_centroid`              |
#
# Each mode runs K-Medoids (optimal k via silhouette) and Ward
# hierarchical clustering independently, per protein-ligand pair.
#
# ## How to Use
#
# ### Required Arguments
# - `--results-csv`   Path to the PoseBusters filtered results CSV
#                     (e.g. `posebusters_filtered_results.csv`)
#
# ### Optional Arguments
# - `--mode`          Distance mode: `centroid`, `rmsd`, or `hybrid` (default: `hybrid`)
# - `--hybrid-weight` RMSD weight for hybrid mode, 0–1 (default: `0.5`).
#                     Formula: `w × norm(RMSD) + (1-w) × norm(Centroid)`.
#                     Ignored unless `--mode hybrid`.
# - `--output-dir`    Output root directory (default: `<results-csv-parent>/clustering_v3_results`).
#                     A sub-directory named after the mode is created automatically.
#
# ### Command Line
# ```bash
# python pose_clustering.py \
#     --results-csv /path/to/posebusters_filtered_results.csv \
#     --mode hybrid \
#     --hybrid-weight 0.7 \
#     --output-dir /path/to/output
# ```
#
# ### Jupyter Notebook
# ```python
# %run pose_clustering.py \
#     --results-csv /path/to/posebusters_filtered_results.csv \
#     --mode hybrid \
#     --hybrid-weight 0.7 \
#     --output-dir /path/to/output
# ```
#
# ### Outputs
# - Per-combo silhouette search plots (K-Medoids & Hierarchical)
# - 3D centroid scatter plots per protein-ligand combination
# - Cluster stats, pocket overlap tables, and conformity scores (printed)
# - All figures saved as PNG to `<output-dir>/<mode>/`

# %% [markdown]
# # Configuration

# %%
# ==========================================================================
# USER CONFIGURATION
# ==========================================================================
from __future__ import annotations
import argparse

_parser = argparse.ArgumentParser(
    description="Pose Clustering v3 — Centroid, RMSD, and Weighted Hybrid",
    add_help=False,
)
_parser.add_argument("--results-csv", type=str, required=True,
                     help="Path to posebusters_filtered_results.csv")
_parser.add_argument("--mode", type=str, default="hybrid",
                     choices=["centroid", "rmsd", "hybrid"],
                     help="Distance mode: centroid, rmsd, or hybrid (default: hybrid)")
_parser.add_argument("--hybrid-weight", type=float, default=0.5, dest="hybrid_weight",
                     help="RMSD weight for hybrid mode: dist = w*norm_rmsd + (1-w)*norm_centroid (default: 0.5)")
_parser.add_argument("--output-dir", type=str, default=None, dest="output_dir",
                     help="Output root directory (default: <results-csv-parent>/clustering_v3_results)")

_args, _unknown = _parser.parse_known_args()

POSEBUSTERS_CSV = _args.results_csv
DISTANCE_MODE = _args.mode
HYBRID_WEIGHT_RMSD = _args.hybrid_weight

# Derive output directory from argument or default
if _args.output_dir:
    OUTPUT_DIR = _args.output_dir
else:
    from pathlib import Path as _Path
    OUTPUT_DIR = str(_Path(POSEBUSTERS_CSV).parent / "clustering_v3_results")

# Which docking methods to keep (must match values in the "docking_method" column)
METHODS_TO_INCLUDE = [
    "autodock",
    "diffdock",
    "equibind_guided",
]

# ---------- Clustering parameters ----------
K_RANGE = range(2, 11)       # range of k to evaluate
RANDOM_STATE = 42
MAX_CENTROID_DIST = 100       # Å — outlier filter

# %% [markdown]
# # Imports & Validation

# %%


import math
import re
import warnings
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    raise ImportError("RDKit is required.  Install with: conda install -c conda-forge rdkit")

try:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
    from scipy.spatial.distance import pdist, squareform
    SKLEARN_AVAILABLE = True
except ImportError:
    raise ImportError("scikit-learn and scipy are required.")

import kmedoids as _kmedoids_pkg

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import seaborn as sns

# ============================================================================
# Derived constants
# ============================================================================
METHODS_LIST = list(METHODS_TO_INCLUDE)

output_root = Path(OUTPUT_DIR)
output_dir = output_root / DISTANCE_MODE
output_dir.mkdir(parents=True, exist_ok=True)

_PALETTE = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6", "#f39c12", "#1abc9c"]
_MARKERS = ["o", "^", "s", "D", "v", "P"]
METHOD_COLORS  = {m: _PALETTE[i % len(_PALETTE)] for i, m in enumerate(METHODS_LIST)}
METHOD_MARKERS = {m: _MARKERS[i % len(_MARKERS)] for i, m in enumerate(METHODS_LIST)}

print("=" * 100)
print(f"POSE CLUSTERING v3  —  mode = {DISTANCE_MODE.upper()}")
if DISTANCE_MODE == "hybrid":
    print(f"  Hybrid weight:  RMSD × {HYBRID_WEIGHT_RMSD:.2f}  +  Centroid × {1 - HYBRID_WEIGHT_RMSD:.2f}")
print(f"  Output:  {output_dir}")
print("=" * 100)

# %% [markdown]
# # Load & Filter Data from CSV

# %%
# ============================================================================
# LOAD CSV & FILTER  (identical logic to v2)
# ============================================================================

print(f"\nLoading PoseBusters results from:\n  {POSEBUSTERS_CSV}")

df_raw = pd.read_csv(POSEBUSTERS_CSV)
print(f"  Total rows in CSV: {len(df_raw)}")

df_raw = df_raw[df_raw["docking_method"].isin(METHODS_TO_INCLUDE)].copy()
print(f"  Rows after filtering to {METHODS_TO_INCLUDE}: {len(df_raw)}")

if df_raw.empty:
    raise RuntimeError(
        f"No rows for docking_method in {METHODS_TO_INCLUDE}.\n"
        f"Available: {sorted(pd.read_csv(POSEBUSTERS_CSV)['docking_method'].unique())}"
    )

# ── Identify boolean PoseBusters test columns ────────────────────────────
_exclude_cols = {
    'mol_true_loaded', 'mol_cond_loaded',
    'number_short_outlier_bonds', 'number_long_outlier_bonds',
    'number_outlier_angles', 'number_clashes',
    'number_non-aromatic_rings_pass', 'number_aromatic_rings_pass',
    'number_non-aromatic_rings_checked', 'number_aromatic_rings_checked',
    'number_double_bonds_checked', 'number_double_bonds_pass',
    'number_valid_bonds', 'number_valid_angles', 'number_valid_noncov_pairs',
    'number_noncov_pairs', 'number_bonds', 'number_angles',
    'num_h_added',
}
_meta_cols = {
    'file_path', 'filepath', 'file', 'path', 'sdf_file', 'sdf_path',
    'method', 'docking_method', 'protein', 'ligand', 'pose_rank', 'rank',
    'molecule', 'mol_name', 'name', 'complex', 'protein_path', 'ligand_path',
    'mol_pred', 'mol_true', 'mol_cond', 'pose_file', 'pose_name',
    'file_format', 'protein_file_used', 'posebusters_mode',
    'diffdock_confidence',
}
_numeric_cols = {
    'mol_pred_energy', 'ensemble_avg_energy', 'energy_ratio',
    'shortest_bond_relative_length', 'longest_bond_relative_length',
    'most_extreme_relative_angle',
    'shortest_noncovalent_relative_distance',
    'aromatic_ring_maximum_distance_from_plane',
    'non-aromatic_ring_maximum_distance_from_plane',
    'double_bond_maximum_distance_from_plane',
}

test_cols = []
for col in df_raw.columns:
    cl = col.lower().strip()
    if cl in _meta_cols or col in _exclude_cols:
        continue
    if cl.startswith("number_") or cl.startswith("num_"):
        continue
    if cl in _numeric_cols:
        continue
    uv = df_raw[col].dropna().unique()
    if set(uv).issubset({True, False, 1, 0, 1.0, 0.0, "True", "False", "true", "false"}):
        test_cols.append(col)

print(f"  Boolean test columns detected ({len(test_cols)}): {test_cols}")

for tc in test_cols:
    if df_raw[tc].dtype == object:
        df_raw[tc] = df_raw[tc].map({"True": True, "true": True, "False": False, "false": False})
    df_raw[tc] = df_raw[tc].astype(bool)

df_raw["all_passed"] = df_raw[test_cols].all(axis=1)
df_passed = df_raw[df_raw["all_passed"]].copy()
print(f"  Poses passing ALL PoseBusters tests: {len(df_passed)} / {len(df_raw)}")

for m in METHODS_TO_INCLUDE:
    nm = (df_passed["docking_method"] == m).sum()
    nt = (df_raw["docking_method"] == m).sum()
    print(f"    {m}: {nm} passed / {nt} total")

# %% [markdown]
# # Feature Extraction — Centroids & Heavy-Atom Coords

# %%
# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def _extract_rank(pose_name: str) -> int:
    """Extract numeric rank from a pose_name string."""
    if pd.isna(pose_name):
        return 1
    s = str(pose_name)
    for pat in (r"rank(\d+)", r"pose(\d+)", r"model(\d+)"):
        m = re.search(pat, s)
        if m:
            return int(m.group(1))
    return 1


def load_heavy_atom_mol(sdf_path: str) -> Optional[Chem.Mol]:
    """Load an SDF and return a heavy-atom-only mol (or None)."""
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
    """Return (N, 3) heavy-atom coordinate array."""
    conf = mol.GetConformer()
    return np.array([
        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
        for i in range(mol.GetNumAtoms())
    ])


def centroid_from_mol(mol: Chem.Mol) -> np.ndarray:
    """Return (3,) centroid vector from heavy atoms."""
    return heavy_atom_coords(mol).mean(axis=0)


def compute_heavy_atom_rmsd(mol_a: Chem.Mol, mol_b: Chem.Mol) -> float:
    """Best-fit heavy-atom RMSD via RDKit; fallback to centroid distance."""
    try:
        if mol_a.GetNumAtoms() == mol_b.GetNumAtoms():
            return rdMolAlign.GetBestRMS(mol_a, mol_b)
        else:
            # Atom count mismatch — fall back to centroid-centroid distance
            return float(np.linalg.norm(
                centroid_from_mol(mol_a) - centroid_from_mol(mol_b)
            ))
    except Exception:
        try:
            return float(np.linalg.norm(
                centroid_from_mol(mol_a) - centroid_from_mol(mol_b)
            ))
        except Exception:
            return np.nan


# ============================================================================
# EXTRACT FEATURES PER POSE
# ============================================================================

print("\n" + "=" * 100)
print("Extracting heavy-atom features from PoseBusters-PASSED poses ...")
print("=" * 100)

all_poses: List[dict] = []
loaded_mols: Dict[str, Chem.Mol] = {}   # keyed by file_path
n_ok = {m: 0 for m in METHODS_TO_INCLUDE}
n_fail = {m: 0 for m in METHODS_TO_INCLUDE}

for _, row in df_passed.iterrows():
    method = row["docking_method"]
    protein = re.sub(r"_cleaned$", "", str(row["protein"]))
    ligand = row["ligand"]
    pose_file = str(row["pose_file"])
    pose_name = row.get("pose_name", "")
    rank = _extract_rank(pose_name)

    if not Path(pose_file).exists():
        n_fail[method] += 1
        continue

    mol = load_heavy_atom_mol(pose_file)
    if mol is None:
        n_fail[method] += 1
        continue

    centroid = centroid_from_mol(mol)
    n_ok[method] += 1
    loaded_mols[pose_file] = mol
    all_poses.append({
        "method": method,
        "protein": protein,
        "ligand": ligand,
        "pose_rank": rank,
        "file_path": pose_file,
        "centroid_x": centroid[0],
        "centroid_y": centroid[1],
        "centroid_z": centroid[2],
    })

df_poses = pd.DataFrame(all_poses)

print(f"\n  Feature extraction summary:")
for m in METHODS_TO_INCLUDE:
    print(f"    {m}: {n_ok[m]} ok / {n_fail[m]} failed")
print(f"  Total poses loaded: {len(df_poses)}")
if not df_poses.empty:
    print(f"  Unique protein-ligand combos: {df_poses.groupby(['protein', 'ligand']).ngroups}")

# ── Outlier removal ──
n_before = len(df_poses)
keep = pd.Series(True, index=df_poses.index)
for (protein, ligand), grp in df_poses.groupby(["protein", "ligand"]):
    med = grp[["centroid_x", "centroid_y", "centroid_z"]].median().values
    d = np.linalg.norm(grp[["centroid_x", "centroid_y", "centroid_z"]].values - med, axis=1)
    bad = d > MAX_CENTROID_DIST
    if bad.any():
        keep.loc[grp.index[bad]] = False
        for i, idx in enumerate(grp.index[bad]):
            r = df_poses.loc[idx]
            print(f"    OUTLIER: {r['method']} rank {r['pose_rank']} "
                  f"({r['protein']}+{r['ligand']}), dist={d[np.where(bad)[0][i]]:.0f} Å")
df_poses = df_poses[keep].reset_index(drop=True)
n_removed = n_before - len(df_poses)
if n_removed:
    print(f"  Removed {n_removed} outliers (>{MAX_CENTROID_DIST} Å from median)")

# %% [markdown]
# # Distance Matrix Construction

# %%
# ============================================================================
# BUILD DISTANCE MATRICES  (per protein-ligand combination)
# ============================================================================
#
# For each protein-ligand group the script computes:
#   centroid_dist[i,j] = Euclidean distance between centroid_i and centroid_j
#   rmsd_dist[i,j]     = Heavy-atom RMSD between pose_i and pose_j
#   hybrid_dist[i,j]   = w * norm(rmsd) + (1-w) * norm(centroid)
#
# Only the matrix required by the chosen DISTANCE_MODE is actually used for
# clustering, but all three are computed so they can be compared.
# ============================================================================

import sys
from functools import lru_cache


def _pairwise_centroid_dm(centroids: np.ndarray) -> np.ndarray:
    """Euclidean distance matrix from (N,3) centroid array."""
    return squareform(pdist(centroids, metric="euclidean"))


def _pairwise_rmsd_dm(file_paths: List[str], mol_cache: Dict[str, Chem.Mol]) -> np.ndarray:
    """Full N×N heavy-atom RMSD distance matrix."""
    n = len(file_paths)
    dm = np.zeros((n, n))
    for i in range(n):
        mol_i = mol_cache.get(file_paths[i])
        if mol_i is None:
            dm[i, :] = np.nan
            continue
        for j in range(i + 1, n):
            mol_j = mol_cache.get(file_paths[j])
            if mol_j is None:
                dm[i, j] = dm[j, i] = np.nan
                continue
            r = compute_heavy_atom_rmsd(mol_i, mol_j)
            dm[i, j] = dm[j, i] = r
    return dm


def _normalise(dm: np.ndarray) -> np.ndarray:
    """Min-max normalise a distance matrix to [0, 1]."""
    finite = dm[np.isfinite(dm)]
    if len(finite) == 0:
        return np.zeros_like(dm)
    lo, hi = finite.min(), finite.max()
    if hi - lo < 1e-12:
        return np.zeros_like(dm)
    out = (dm - lo) / (hi - lo)
    out[~np.isfinite(out)] = 1.0  # NaN → max distance
    return out


def build_distance_matrix(
    group_df: pd.DataFrame,
    mol_cache: Dict[str, Chem.Mol],
    mode: str = "centroid",
    w_rmsd: float = 0.5,
) -> np.ndarray:
    """Return the distance matrix for the requested mode.

    Parameters
    ----------
    group_df : DataFrame with centroid_x/y/z and file_path columns.
    mol_cache : dict mapping file_path → RDKit Mol.
    mode : "centroid", "rmsd", or "hybrid".
    w_rmsd : weight for RMSD in hybrid mode (0–1).

    Returns
    -------
    dm : (N, N) numpy distance matrix.
    """
    centroids = group_df[["centroid_x", "centroid_y", "centroid_z"]].values
    file_paths = group_df["file_path"].tolist()

    if mode == "centroid":
        return _pairwise_centroid_dm(centroids)

    if mode == "rmsd":
        return _pairwise_rmsd_dm(file_paths, mol_cache)

    if mode == "hybrid":
        dm_c = _pairwise_centroid_dm(centroids)
        dm_r = _pairwise_rmsd_dm(file_paths, mol_cache)
        nc = _normalise(dm_c)
        nr = _normalise(dm_r)
        return w_rmsd * nr + (1.0 - w_rmsd) * nc

    raise ValueError(f"Unknown distance mode: {mode!r}")


print(f"\nDistance mode:  {DISTANCE_MODE}")
if DISTANCE_MODE == "hybrid":
    print(f"  Formula:  {HYBRID_WEIGHT_RMSD:.2f} × norm(RMSD) + {1 - HYBRID_WEIGHT_RMSD:.2f} × norm(Centroid)")

# %% [markdown]
# # Clustering Engines

# %%
# ============================================================================
# CLUSTERING ENGINES
# ============================================================================


def _kmedoids_with_dm(dm: np.ndarray, k: int, random_state: int = 42):
    """Run K-Medoids (FasterPAM) on a precomputed distance matrix."""
    result = _kmedoids_pkg.fasterpam(dm, k, random_state=random_state)
    labels = np.array(result.labels)
    inertia = result.loss
    return labels, inertia


def find_optimal_k_kmedoids(dm: np.ndarray, k_range, random_state=42):
    """Sweep k and pick the one with the best silhouette score."""
    n = dm.shape[0]
    results = {"k": [], "silhouette": [], "inertia": []}
    for k in k_range:
        if k >= n:
            break
        labels, inertia = _kmedoids_with_dm(dm, k, random_state)
        sil = silhouette_score(dm, labels, metric="precomputed") if len(set(labels)) > 1 else -1
        results["k"].append(k)
        results["silhouette"].append(sil)
        results["inertia"].append(inertia)
    if not results["k"]:
        return {"optimal_k": 1, "best_silhouette": -1, "k": [], "silhouette": [], "inertia": []}
    best_idx = int(np.argmax(results["silhouette"]))
    results["optimal_k"] = results["k"][best_idx]
    results["best_silhouette"] = results["silhouette"][best_idx]
    return results


def find_optimal_k_ward(centroids: np.ndarray, k_range):
    """Ward clustering uses Euclidean coordinates (always centroid-based internally)."""
    results = {"k": [], "silhouette": [], "inertia": []}
    for k in k_range:
        if k >= len(centroids):
            break
        model = AgglomerativeClustering(n_clusters=k, linkage="ward")
        labels = model.fit_predict(centroids)
        sil = silhouette_score(centroids, labels) if len(set(labels)) > 1 else -1
        # pseudo-inertia: sum of within-cluster distances to centroid
        inertia = 0.0
        for c in set(labels):
            members = centroids[labels == c]
            center = members.mean(axis=0)
            inertia += np.sum(np.linalg.norm(members - center, axis=1) ** 2)
        results["k"].append(k)
        results["silhouette"].append(sil)
        results["inertia"].append(inertia)
    if not results["k"]:
        return {"optimal_k": 1, "best_silhouette": -1, "k": [], "silhouette": [], "inertia": []}
    best_idx = int(np.argmax(results["silhouette"]))
    results["optimal_k"] = results["k"][best_idx]
    results["best_silhouette"] = results["silhouette"][best_idx]
    return results


def ward_on_dm(dm: np.ndarray, k: int) -> np.ndarray:
    """Ward-like agglomerative clustering on a precomputed distance matrix.

    Ward linkage proper requires Euclidean coordinates, so when the distance
    matrix comes from RMSD or Hybrid mode we use *average* linkage instead
    (which accepts precomputed distances) as the best hierarchical analogue.
    """
    condensed = squareform(dm, checks=False)
    Z = linkage(condensed, method="average")
    return fcluster(Z, t=k, criterion="maxclust") - 1   # 0-indexed


def find_optimal_k_hierarchical_dm(dm: np.ndarray, k_range):
    """Hierarchical clustering on precomputed distance matrix (average linkage)."""
    n = dm.shape[0]
    results = {"k": [], "silhouette": [], "inertia": []}
    for k in k_range:
        if k >= n:
            break
        labels = ward_on_dm(dm, k)
        sil = silhouette_score(dm, labels, metric="precomputed") if len(set(labels)) > 1 else -1
        results["k"].append(k)
        results["silhouette"].append(sil)
        results["inertia"].append(0.0)
    if not results["k"]:
        return {"optimal_k": 1, "best_silhouette": -1, "k": [], "silhouette": [], "inertia": []}
    best_idx = int(np.argmax(results["silhouette"]))
    results["optimal_k"] = results["k"][best_idx]
    results["best_silhouette"] = results["silhouette"][best_idx]
    return results


print("Clustering engines loaded.")

# %% [markdown]
# # Shared Stats & Overlap Helpers

# %%
# ============================================================================
# SHARED HELPERS
# ============================================================================

def compute_cluster_stats(df: pd.DataFrame, cluster_col: str) -> pd.DataFrame:
    """Compute cluster centers, sizes, and radii."""
    rows = []
    for (protein, ligand, cid), grp in df.groupby(["protein", "ligand", cluster_col]):
        coords = grp[["centroid_x", "centroid_y", "centroid_z"]].values
        center = coords.mean(axis=0)
        radius = float(np.max(np.linalg.norm(coords - center, axis=1))) if len(coords) > 1 else 0.0
        mc = grp["method"].value_counts().to_dict()
        row = {
            "protein": protein, "ligand": ligand, "cluster_id": cid,
            "cluster_method": cluster_col.replace("cluster_", ""),
            "n_poses": len(grp),
            "center_x": center[0], "center_y": center[1], "center_z": center[2],
            "radius_A": round(radius, 2),
        }
        for m in METHODS_LIST:
            row[m] = mc.get(m, 0)
        row["n_methods"] = len(mc)
        row["methods_present"] = ", ".join(sorted(mc.keys()))
        rows.append(row)
    return pd.DataFrame(rows)


def identify_pocket_overlap(df_stats: pd.DataFrame) -> pd.DataFrame:
    out = []
    for _, row in df_stats.iterrows():
        present = [m for m in METHODS_LIST if row.get(m, 0) > 0]
        if len(present) >= 2:
            ptype = "Overlapping"
        elif len(present) == 1:
            ptype = f"Unique ({present[0]})"
        else:
            ptype = "Empty"
        out.append({**row.to_dict(), "pocket_type": ptype, "n_methods_in_pocket": len(present)})
    return pd.DataFrame(out)


def score_conformity(group_df: pd.DataFrame, cluster_col: str) -> dict:
    """Composite conformity score (higher = more agreement across methods)."""
    nc = group_df[cluster_col].nunique()
    nm = group_df["method"].nunique()
    np_ = len(group_df)

    overlap_poses = 0
    overlap_clusters = 0
    for _, cg in group_df.groupby(cluster_col):
        if cg["method"].nunique() >= 2:
            overlap_clusters += 1
            overlap_poses += len(cg)
    frac = overlap_poses / np_ if np_ else 0

    method_cents = {}
    for method, mg in group_df.groupby("method"):
        method_cents[method] = mg[["centroid_x", "centroid_y", "centroid_z"]].values.mean(axis=0)
    imd = [np.linalg.norm(c1 - c2) for (_, c1), (_, c2) in combinations(method_cents.items(), 2)]
    mean_imd = float(np.mean(imd)) if imd else 0.0

    all_c = group_df[["centroid_x", "centroid_y", "centroid_z"]].values
    spread = float(np.mean(np.linalg.norm(all_c - all_c.mean(axis=0), axis=1)))

    cscore = (
        frac * 0.4
        + max(0, 1 - mean_imd / 50) * 0.3
        + max(0, 1 - spread / 30) * 0.2
        + max(0, 1 - nc / 10) * 0.1
    )
    return {
        "n_clusters": nc, "n_methods": nm, "n_poses": np_,
        "overlapping_clusters": overlap_clusters,
        "overlap_fraction": round(frac, 3),
        "mean_inter_method_dist": round(mean_imd, 2),
        "spatial_spread": round(spread, 2),
        "conformity_score": round(cscore, 4),
    }


print("Shared helper functions loaded.")

# %% [markdown]
# # Run Clustering

# %%
# ============================================================================
# MAIN CLUSTERING LOOP
# ============================================================================

if df_poses.empty:
    raise RuntimeError("No poses available after filtering.")

print("\n" + "=" * 100)
print(f"CLUSTERING  —  mode = {DISTANCE_MODE.upper()}")
print("=" * 100)

df_clust = df_poses.copy()
df_clust["cluster_kmedoids"] = -1
df_clust["cluster_hierarchical"] = -1

km_info: Dict[str, dict] = {}
hi_info: Dict[str, dict] = {}

combo_groups = list(df_clust.groupby(["protein", "ligand"]))
n_combos = len(combo_groups)

for gi, ((protein, ligand), grp) in enumerate(combo_groups, 1):
    combo_key = f"{protein}__{ligand}"
    idx = grp.index
    n_poses = len(grp)

    print(f"\n  [{gi}/{n_combos}] {combo_key}  ({n_poses} poses)")

    # ── Trivial case ──
    if n_poses < 3:
        df_clust.loc[idx, "cluster_kmedoids"] = 0
        df_clust.loc[idx, "cluster_hierarchical"] = 0
        km_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        hi_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        print(f"    → <3 poses, assigning to cluster 0")
        continue

    # ── Build distance matrix ──
    dm = build_distance_matrix(grp, loaded_mols, mode=DISTANCE_MODE, w_rmsd=HYBRID_WEIGHT_RMSD)

    # Handle NaN in dm (replace with max finite value)
    finite_vals = dm[np.isfinite(dm)]
    if len(finite_vals) == 0:
        print(f"    ⚠ Distance matrix is all NaN — skipping")
        df_clust.loc[idx, "cluster_kmedoids"] = 0
        df_clust.loc[idx, "cluster_hierarchical"] = 0
        km_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        hi_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        continue
    dm[~np.isfinite(dm)] = finite_vals.max()

    max_k = min(max(K_RANGE), n_poses - 1)
    actual_range = range(2, max_k + 1)
    if len(actual_range) == 0:
        df_clust.loc[idx, "cluster_kmedoids"] = 0
        df_clust.loc[idx, "cluster_hierarchical"] = 0
        km_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        hi_info[combo_key] = {"n_poses": n_poses, "optimal_k": 1, "best_silhouette": -1, "search_results": None}
        continue

    # ── K-Medoids (precomputed distance matrix) ──
    km_res = find_optimal_k_kmedoids(dm, actual_range, RANDOM_STATE)
    k_opt_km = km_res["optimal_k"]
    labels_km, _ = _kmedoids_with_dm(dm, k_opt_km, RANDOM_STATE)
    df_clust.loc[idx, "cluster_kmedoids"] = labels_km
    km_info[combo_key] = {
        "n_poses": n_poses, "optimal_k": k_opt_km,
        "best_silhouette": km_res["best_silhouette"],
        "search_results": km_res,
    }
    print(f"    K-Medoids: k={k_opt_km}  sil={km_res['best_silhouette']:.3f}")

    # ── Hierarchical ──
    if DISTANCE_MODE == "centroid":
        # Pure centroid mode → use proper Ward linkage on coordinates
        centroids = grp[["centroid_x", "centroid_y", "centroid_z"]].values
        hi_res = find_optimal_k_ward(centroids, actual_range)
        k_opt_hi = hi_res["optimal_k"]
        model_hi = AgglomerativeClustering(n_clusters=k_opt_hi, linkage="ward")
        labels_hi = model_hi.fit_predict(centroids)
        # Store linkage matrix for dendrogram
        Z_hi = linkage(centroids, method="ward")
    else:
        # RMSD or Hybrid → average-linkage on precomputed DM
        hi_res = find_optimal_k_hierarchical_dm(dm, actual_range)
        k_opt_hi = hi_res["optimal_k"]
        labels_hi = ward_on_dm(dm, k_opt_hi)
        # Store linkage matrix for dendrogram
        condensed = squareform(dm, checks=False)
        Z_hi = linkage(condensed, method="average")

    df_clust.loc[idx, "cluster_hierarchical"] = labels_hi
    hi_info[combo_key] = {
        "n_poses": n_poses, "optimal_k": k_opt_hi,
        "best_silhouette": hi_res["best_silhouette"],
        "search_results": hi_res,
        "linkage_matrix": Z_hi,
        "pose_labels": grp["method"].values.tolist(),
    }
    print(f"    Hierarchical: k={k_opt_hi}  sil={hi_res['best_silhouette']:.3f}")

df_clust["cluster_kmedoids"] = df_clust["cluster_kmedoids"].astype(int)
df_clust["cluster_hierarchical"] = df_clust["cluster_hierarchical"].astype(int)

print("\n" + "=" * 100)
print("CLUSTERING COMPLETE")
print("=" * 100)

# %% [markdown]
# # Results Summary

# %%
# ============================================================================
# RESULTS TABLES
# ============================================================================

print(f"\n{'='*100}")
print(f"K-MEDOIDS RESULTS  (mode={DISTANCE_MODE})")
print(f"{'='*100}")
print(f"  {'Protein-Ligand':<55} {'k':>5} {'Silhouette':>12} {'Poses':>8}")
print("  " + "-" * 85)
for ck, info in sorted(km_info.items()):
    print(f"  {ck:<55} {info['optimal_k']:>5} {info['best_silhouette']:>12.3f} {info['n_poses']:>8}")

df_km_stats = compute_cluster_stats(df_clust, "cluster_kmedoids")
df_km_overlap = identify_pocket_overlap(df_km_stats)

print(f"\n{'='*100}")
print("K-MEDOIDS — Binding Pocket Clusters")
print(f"{'='*100}")
km_display_cols = ["protein", "ligand", "cluster_id", "n_poses", "radius_A"] + \
                  METHODS_LIST + ["pocket_type", "n_methods_in_pocket"]
print(df_km_overlap[km_display_cols].sort_values(["protein", "ligand", "cluster_id"]).to_string(index=False))

# ── Pocket overlap summary ──
print(f"\n{'='*100}")
print("K-MEDOIDS — Overlapping vs Unique")
print(f"{'='*100}")
agg_dict = {"n_poses": ["count", "sum"]}
for m in METHODS_LIST:
    agg_dict[m] = "sum"
ps = df_km_overlap.groupby("pocket_type").agg(agg_dict)
ps.columns = ["N_Clusters", "Total_Poses"] + [f"{m}_Poses" for m in METHODS_LIST]
ps = ps.reset_index()
print(ps.to_string(index=False))

# ── Hierarchical ──
print(f"\n{'='*100}")
print(f"HIERARCHICAL RESULTS  (mode={DISTANCE_MODE})")
print(f"{'='*100}")
print(f"  {'Protein-Ligand':<55} {'k':>5} {'Silhouette':>12} {'Poses':>8}")
print("  " + "-" * 85)
for ck, info in sorted(hi_info.items()):
    print(f"  {ck:<55} {info['optimal_k']:>5} {info['best_silhouette']:>12.3f} {info['n_poses']:>8}")

df_hi_stats = compute_cluster_stats(df_clust, "cluster_hierarchical")
df_hi_overlap = identify_pocket_overlap(df_hi_stats)

print(f"\n{'='*100}")
print("HIERARCHICAL — Binding Pocket Clusters")
print(f"{'='*100}")
hi_display_cols = ["protein", "ligand", "cluster_id", "n_poses", "radius_A"] + \
                  METHODS_LIST + ["pocket_type", "n_methods_in_pocket"]
print(df_hi_overlap[hi_display_cols].sort_values(["protein", "ligand", "cluster_id"]).to_string(index=False))

# %% [markdown]
# # Conformity Scoring

# %%
# ============================================================================
# CONFORMITY / DIVERSITY SCORING
# ============================================================================

print(f"\n{'='*100}")
print("CONFORMITY SCORING")
print(f"{'='*100}")

score_rows = []
for (protein, ligand), grp in df_clust.groupby(["protein", "ligand"]):
    s_km = score_conformity(grp, "cluster_kmedoids")
    s_hi = score_conformity(grp, "cluster_hierarchical")
    # Average of both clusterers
    avg_score = round((s_km["conformity_score"] + s_hi["conformity_score"]) / 2, 4)
    score_rows.append({
        "protein": protein, "ligand": ligand,
        "n_poses": s_km["n_poses"], "n_methods": s_km["n_methods"],
        "km_k": km_info.get(f"{protein}__{ligand}", {}).get("optimal_k", 1),
        "km_sil": km_info.get(f"{protein}__{ligand}", {}).get("best_silhouette", -1),
        "km_conformity": s_km["conformity_score"],
        "hi_k": hi_info.get(f"{protein}__{ligand}", {}).get("optimal_k", 1),
        "hi_sil": hi_info.get(f"{protein}__{ligand}", {}).get("best_silhouette", -1),
        "hi_conformity": s_hi["conformity_score"],
        "avg_conformity": avg_score,
        "overlap_frac": s_km["overlap_fraction"],
        "mean_inter_method_dist": s_km["mean_inter_method_dist"],
        "spatial_spread": s_km["spatial_spread"],
    })

df_scores = pd.DataFrame(score_rows).sort_values("avg_conformity", ascending=False)
print(df_scores.to_string(index=False))

# %% [markdown]
# # Silhouette Search Plots

# %%
# ============================================================================
# SILHOUETTE SCORE SEARCH PLOTS
# ============================================================================

combos_with_search = [(k, v) for k, v in km_info.items() if v.get("search_results")]
n_plots = len(combos_with_search)

if n_plots > 0:
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 4 * n_rows), squeeze=False)
    axes_flat = axes.flatten()

    for i, (ck, info) in enumerate(combos_with_search):
        ax = axes_flat[i]
        sr = info["search_results"]
        ax.plot(sr["k"], sr["silhouette"], "o-", color="#3498db")
        ax.axvline(sr["optimal_k"], ls="--", color="#e74c3c", alpha=0.7)
        ax.set_title(ck.replace("__", "\n"), fontsize=9)
        ax.set_xlabel("k"); ax.set_ylabel("Silhouette")
        ax.grid(alpha=0.3)

    for j in range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.suptitle(f"K-Medoids Silhouette Search  ({DISTANCE_MODE} mode)",
                 fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    plt.savefig(output_dir / "silhouette_search_kmedoids.png", dpi=150, bbox_inches="tight")
    plt.show()

# %% [markdown]
# # Hierarchical Dendrograms

# %%
# ============================================================================
# HIERARCHICAL CLUSTERING DENDROGRAMS
# ============================================================================

combos_with_dendro = [(k, v) for k, v in hi_info.items() if v.get("linkage_matrix") is not None]

# Maximum width cap so dendrograms stay readable on screen / paper
_DENDRO_MAX_WIDTH = 18   # inches

for ck, info in combos_with_dendro:
    Z_hi = info["linkage_matrix"]
    k_opt = info["optimal_k"]
    pose_labels = info.get("pose_labels", None)
    n_poses = info["n_poses"]

    # Width: scale with leaf count but clamp to max
    fig_w = min(_DENDRO_MAX_WIDTH, max(8, n_poses * 0.25))
    fig_h = 5.5
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Color threshold: cut at the level that gives k_opt clusters
    if k_opt > 1 and len(Z_hi) >= k_opt:
        color_threshold = Z_hi[-(k_opt - 1), 2] if k_opt <= len(Z_hi) else None
    else:
        color_threshold = None

    # Leaf font size & line width: shrink for large trees so nothing overlaps
    leaf_fs = max(2, min(8, int(140 / max(n_poses, 1))))
    line_w  = max(0.3, min(1.5, 60 / max(n_poses, 1)))

    dn = dendrogram(
        Z_hi,
        ax=ax,
        labels=pose_labels,
        color_threshold=color_threshold,
        above_threshold_color="grey",
        leaf_rotation=90,
        leaf_font_size=leaf_fs,
    )

    # Thin out the dendrogram lines for large trees
    for coll in ax.collections:
        coll.set_linewidth(line_w)
    for line in ax.lines:
        line.set_linewidth(line_w)

    ax.set_title(f"Hierarchical Dendrogram — {ck.replace('__', ' / ')}\n"
                 f"({DISTANCE_MODE} mode, k={k_opt}, {n_poses} poses)",
                 fontsize=11, fontweight="bold", pad=12)
    linkage_type = "Ward" if DISTANCE_MODE == "centroid" else "Average"
    ax.set_ylabel(f"Distance ({linkage_type} linkage)")
    ax.set_xlabel("Pose (method)")

    # Draw a horizontal line at the cut threshold
    if color_threshold is not None:
        ax.axhline(y=color_threshold, ls="--", color="#e74c3c", alpha=0.6,
                   label=f"Cut for k={k_opt}")
        ax.legend(loc="upper right", fontsize=9)

    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    combo_safe = ck.replace("/", "_")
    fname = f"dendrogram_{combo_safe}_{DISTANCE_MODE}.png"
    fig.savefig(output_dir / fname, dpi=150, bbox_inches="tight")
    plt.show()
    plt.close(fig)

# %% [markdown]
# # 3D Scatter Visualisation

# %%
# ============================================================================
# 3D CLUSTER SCATTER PLOTS  (K-Medoids & Hierarchical side by side)
# ============================================================================

combo_list = list(df_clust.groupby(["protein", "ligand"]))
n_combos = len(combo_list)

legend_elements = [
    plt.Line2D([0], [0], marker=METHOD_MARKERS[m], color="w",
               markerfacecolor="gray", markersize=10, label=m)
    for m in METHODS_LIST
] + [
    plt.Line2D([0], [0], marker="X", color="w", markerfacecolor="black",
               markersize=12, label="Cluster Center"),
]

for clust_col, clust_label in [("cluster_kmedoids", "K-Medoids"), ("cluster_hierarchical", "Hierarchical")]:
    for ci, ((protein, ligand), grp) in enumerate(combo_list):
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(111, projection="3d")
        ax.set_title(f"{clust_label} — {protein} / {ligand}\n"
                     f"({DISTANCE_MODE} mode)  Color=Cluster, Shape=Method",
                     fontsize=11, fontweight="bold", pad=15)

        clusters = sorted(grp[clust_col].unique())
        cmap = cm.get_cmap("tab10", max(len(clusters), 1))

        for _, row in grp.iterrows():
            c_idx = clusters.index(row[clust_col])
            ax.scatter(
                row["centroid_x"], row["centroid_y"], row["centroid_z"],
                c=[cmap(c_idx)],
                marker=METHOD_MARKERS.get(row["method"], "o"),
                s=60, alpha=0.8, edgecolors="k", linewidths=0.3,
            )

        # Cluster centers
        for cid in clusters:
            mask = grp[clust_col] == cid
            center = grp.loc[mask, ["centroid_x", "centroid_y", "centroid_z"]].mean()
            c_idx = clusters.index(cid)
            ax.scatter(center["centroid_x"], center["centroid_y"], center["centroid_z"],
                       marker="X", s=200, c=[cmap(c_idx)], edgecolors="k", linewidths=1.5)

        ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
        ax.legend(handles=legend_elements, fontsize=8, loc="upper left",
                  bbox_to_anchor=(0.0, 1.0))
        fig.tight_layout()
        combo_safe = f"{protein}__{ligand}".replace("/", "_")
        fname = f"3d_{clust_label.lower().replace('-', '')}_{combo_safe}_{DISTANCE_MODE}.png"
        fig.savefig(output_dir / fname, dpi=150, bbox_inches="tight")
        plt.show()
        plt.close(fig)

# %% [markdown]
# # Heatmaps

# %%
# ============================================================================
# POSE-COUNT HEATMAPS
# ============================================================================

for df_ov, label in [(df_km_overlap, "K-Medoids"), (df_hi_overlap, "Hierarchical")]:
    df_h = df_ov.copy()
    df_h["pocket_label"] = (
        df_h["protein"].str[:20] + " / " +
        df_h["ligand"].str[:15] + " / C" +
        df_h["cluster_id"].astype(str)
    )
    heat = df_h.set_index("pocket_label")[METHODS_LIST]

    if len(heat) > 0:
        fig, ax = plt.subplots(figsize=(max(6, len(METHODS_LIST) * 1.5), max(4, len(heat) * 0.4)))
        im = ax.imshow(heat.values, aspect="auto", cmap="YlOrRd")
        ax.set_xticks(range(len(METHODS_LIST)))
        ax.set_xticklabels(METHODS_LIST, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(heat)))
        ax.set_yticklabels(heat.index, fontsize=7)
        for r in range(heat.shape[0]):
            for c in range(heat.shape[1]):
                ax.text(c, r, str(int(heat.values[r, c])), ha="center", va="center", fontsize=8)
        ax.set_title(f"{label} — Pose Counts per Cluster ({DISTANCE_MODE} mode)",
                     fontsize=12, fontweight="bold", pad=12)
        plt.colorbar(im, ax=ax, shrink=0.7)
        fig.tight_layout(rect=[0, 0, 1, 0.97])
        fname = f"heatmap_{label.lower().replace('-', '')}_{DISTANCE_MODE}.png"
        plt.savefig(output_dir / fname, dpi=150, bbox_inches="tight")
        plt.show()

# %% [markdown]
# # Cross-Method RMSD per Cluster

# %%
# ============================================================================
# CROSS-METHOD CLOSEST-POSE RMSD PER CLUSTER
# ============================================================================

print(f"\n{'='*100}")
print("CROSS-METHOD RMSD PER CLUSTER")
print(f"{'='*100}")

_ABBR = {m: m[:2].upper() for m in METHODS_LIST}

rmsd_rows = []

for clust_col, clust_label in [("cluster_kmedoids", "K-Medoids"), ("cluster_hierarchical", "Hierarchical")]:
    for (protein, ligand), grp in df_clust.groupby(["protein", "ligand"]):
        for cid in sorted(grp[clust_col].unique()):
            cg = grp[grp[clust_col] == cid]
            methods_present = sorted(cg["method"].unique())
            if len(methods_present) < 2:
                continue

            for m1, m2 in combinations(methods_present, 2):
                files_1 = cg[cg["method"] == m1]["file_path"].tolist()
                files_2 = cg[cg["method"] == m2]["file_path"].tolist()

                best_rmsd = np.inf
                for f1 in files_1:
                    mol1 = loaded_mols.get(f1)
                    if mol1 is None:
                        continue
                    for f2 in files_2:
                        mol2 = loaded_mols.get(f2)
                        if mol2 is None:
                            continue
                        r = compute_heavy_atom_rmsd(mol1, mol2)
                        if np.isfinite(r) and r < best_rmsd:
                            best_rmsd = r
                if np.isfinite(best_rmsd):
                    rmsd_rows.append({
                        "cluster_method": clust_label,
                        "protein": protein, "ligand": ligand,
                        "cluster_id": cid,
                        "method_1": _ABBR[m1], "method_2": _ABBR[m2],
                        "closest_rmsd_A": round(best_rmsd, 2),
                    })

df_xrmsd = pd.DataFrame(rmsd_rows)
if df_xrmsd.empty:
    print("  No multi-method clusters found.")
else:
    print(df_xrmsd.to_string(index=False))

# %% [markdown]
# # Cross-Method Cluster RMSD Analysis (Balanced)

# %%
# ============================================================================
# CROSS-METHOD CLUSTER RMSD ANALYSIS  (balanced pose selection + heatmaps)
# ============================================================================
#
# For each cluster that contains poses from multiple docking methods:
#   1. Determine the maximum balanced count N (= min poses across methods).
#   2. From each method, pick the N poses closest to the cluster centroid.
#   3. Compute all pairwise cross-method RMSDs among the selected poses.
#   4. Record per-cluster statistics and produce heatmaps.
# ============================================================================

def _pairwise_rmsd_coords(c1, c2):
    """Compute RMSD between two coordinate arrays (atom-by-atom or centroid fallback)."""
    if c1 is None or c2 is None:
        return np.nan
    if len(c1) == len(c2):
        d = c1 - c2
        return float(np.sqrt(np.mean(np.sum(d ** 2, axis=1))))
    return float(np.linalg.norm(c1.mean(axis=0) - c2.mean(axis=0)))


def rmsd_overview_across_clusters(df_cluster_rmsd_in, methods_list, clust_label,
                                  output_dir_path):
    """
    Generate a comprehensive RMSD overview across clusters and protein-ligand
    pairs: aggregates, rankings, box-plot, heatmap, histogram.

    Parameters
    ----------
    df_cluster_rmsd_in : pd.DataFrame
        Per-cluster cross-method RMSD table.
    methods_list : list[str]
        List of docking method names.
    clust_label : str
        Label for the clustering algorithm (e.g. "K-Medoids").
    output_dir_path : Path
        Directory where figures and CSVs are saved.

    Returns
    -------
    df_overview : pd.DataFrame
        Per protein-ligand aggregate RMSD table.
    """

    _TOOL_ABR = {m: m[:2].upper() for m in methods_list}

    # ------------------------------------------------------------------
    # 1. Aggregate per protein-ligand pair
    # ------------------------------------------------------------------
    agg = (
        df_cluster_rmsd_in
        .groupby(["protein", "ligand"])
        .agg(
            n_clusters=("cluster", "nunique"),
            n_multi_method_clusters=("cluster", "count"),
            mean_rmsd_across_clusters=("mean_rmsd", "mean"),
            median_rmsd_across_clusters=("mean_rmsd", "median"),
            min_cluster_rmsd=("mean_rmsd", "min"),
            max_cluster_rmsd=("mean_rmsd", "max"),
            std_rmsd_across_clusters=("mean_rmsd", "std"),
            global_min_rmsd=("min_rmsd", "min"),
            global_max_rmsd=("max_rmsd", "max"),
            avg_n_methods=("n_methods", "mean"),
        )
        .reset_index()
    )
    for c in agg.columns:
        if agg[c].dtype == float:
            agg[c] = agg[c].round(3)
    agg = agg.sort_values("mean_rmsd_across_clusters", ascending=True)

    # ------------------------------------------------------------------
    # 2. Print ranked table
    # ------------------------------------------------------------------
    print("\n" + "=" * 130)
    print(f"CROSS-CLUSTER RMSD OVERVIEW — PROTEIN-LIGAND RANKING  ({clust_label})")
    print("=" * 130)
    print(agg.to_string(index=False))

    q1 = agg["mean_rmsd_across_clusters"].quantile(0.25)
    q3 = agg["mean_rmsd_across_clusters"].quantile(0.75)
    print(f"\n  RMSD quartile thresholds: Q1 = {q1:.2f} Å, Q3 = {q3:.2f} Å")

    top_agree = agg[agg["mean_rmsd_across_clusters"] <= q1]
    top_diverse = agg[agg["mean_rmsd_across_clusters"] >= q3]

    if not top_agree.empty:
        print(f"\n  HIGH AGREEMENT (mean RMSD <= Q1 = {q1:.2f} Å):")
        for _, r in top_agree.iterrows():
            print(f"    {r['protein']:<35} + {r['ligand']:<20}  "
                  f"mean={r['mean_rmsd_across_clusters']:.3f} Å  "
                  f"({int(r['n_multi_method_clusters'])} multi-method clusters)")

    if not top_diverse.empty:
        print(f"\n  HIGH DIVERSITY (mean RMSD >= Q3 = {q3:.2f} Å):")
        for _, r in top_diverse.iterrows():
            print(f"    {r['protein']:<35} + {r['ligand']:<20}  "
                  f"mean={r['mean_rmsd_across_clusters']:.3f} Å  "
                  f"({int(r['n_multi_method_clusters'])} multi-method clusters)")

    # ------------------------------------------------------------------
    # 3. Box-plot: cluster RMSD distribution per protein-ligand pair
    # ------------------------------------------------------------------
    combo_order = agg.apply(
        lambda r: f"{r['protein'][:20]}\n{r['ligand'][:15]}", axis=1
    ).values

    df_plot = df_cluster_rmsd_in.copy()
    df_plot["combo"] = (
        df_plot["protein"].str[:20] + "\n" + df_plot["ligand"].str[:15]
    )
    df_plot["combo"] = pd.Categorical(df_plot["combo"], categories=combo_order,
                                       ordered=True)

    n_combos_ov = len(combo_order)
    fig_box, ax_box = plt.subplots(
        figsize=(max(10, n_combos_ov * 1.0 + 2), 6)
    )
    sns.boxplot(
        data=df_plot, x="combo", y="mean_rmsd", ax=ax_box,
        color="#74b9ff", linewidth=1.2, fliersize=4,
    )
    sns.stripplot(
        data=df_plot, x="combo", y="mean_rmsd", ax=ax_box,
        color="#2d3436", alpha=0.5, size=4, jitter=0.15,
    )
    ax_box.axhline(y=2.0, color="green", linestyle="--", alpha=0.6,
                    label="2.0 Å (good agreement)")
    ax_box.axhline(y=5.0, color="red", linestyle="--", alpha=0.6,
                    label="5.0 Å (high diversity)")
    ax_box.set_xlabel("Protein-Ligand Pair", fontsize=11)
    ax_box.set_ylabel("Mean Cross-Method RMSD (Å)", fontsize=11)
    ax_box.set_title(
        f"RMSD Distribution Across Clusters per Protein-Ligand Pair ({clust_label})",
        fontsize=13, fontweight="bold",
    )
    ax_box.tick_params(axis="x", rotation=45, labelsize=8)
    ax_box.legend(fontsize=9, loc="upper left")
    ax_box.grid(axis="y", alpha=0.3)
    fig_box.tight_layout()
    tag = clust_label.lower().replace("-", "")
    fig_box.savefig(output_dir_path / f"rmsd_overview_boxplot_{tag}.png",
                    dpi=200, bbox_inches="tight")
    plt.show()

    # ------------------------------------------------------------------
    # 4. Heatmap: protein-ligand × cluster mean RMSD
    # ------------------------------------------------------------------
    df_cluster_rmsd_in = df_cluster_rmsd_in.copy()
    df_cluster_rmsd_in["combo_label"] = (
        df_cluster_rmsd_in["protein"].str[:22] + " / " +
        df_cluster_rmsd_in["ligand"].str[:15]
    )
    combo_mean = (
        df_cluster_rmsd_in.groupby("combo_label")["mean_rmsd"]
        .mean()
        .sort_values()
    )
    sorted_combos = combo_mean.index.tolist()

    heatmap_rows = []
    max_clusters = 0
    for combo_lbl in sorted_combos:
        sub = df_cluster_rmsd_in[
            df_cluster_rmsd_in["combo_label"] == combo_lbl
        ].sort_values("cluster")
        vals = sub["mean_rmsd"].values
        max_clusters = max(max_clusters, len(vals))
        heatmap_rows.append(vals)

    heat_matrix_ov = np.full((len(sorted_combos), max_clusters), np.nan)
    for i, vals in enumerate(heatmap_rows):
        heat_matrix_ov[i, :len(vals)] = vals

    cluster_col_labels = [f"C{j}" for j in range(max_clusters)]
    cmap_rg_ov = sns.color_palette("RdYlGn_r", as_cmap=True)

    fig_hm, ax_hm = plt.subplots(
        figsize=(max(6, max_clusters * 1.5 + 3),
                 max(5, len(sorted_combos) * 0.55 + 2))
    )
    im_hm = ax_hm.imshow(
        heat_matrix_ov, cmap=cmap_rg_ov, aspect="auto",
        vmin=0, vmax=max(np.nanmax(heat_matrix_ov), 5),
    )
    ax_hm.set_xticks(range(max_clusters))
    ax_hm.set_xticklabels(cluster_col_labels, fontsize=10, fontweight="bold")
    ax_hm.set_yticks(range(len(sorted_combos)))
    ax_hm.set_yticklabels(sorted_combos, fontsize=8)
    for i in range(heat_matrix_ov.shape[0]):
        for j in range(heat_matrix_ov.shape[1]):
            v = heat_matrix_ov[i, j]
            if np.isnan(v):
                ax_hm.text(j, i, "—", ha="center", va="center",
                           fontsize=9, color="lightgray")
            else:
                color = "white" if v > np.nanmax(heat_matrix_ov) * 0.55 else "black"
                ax_hm.text(j, i, f"{v:.1f}", ha="center", va="center",
                           fontsize=9, fontweight="bold", color=color)
    cbar_hm = plt.colorbar(im_hm, ax=ax_hm, shrink=0.8, pad=0.02)
    cbar_hm.set_label("Mean Cross-Method RMSD (Å)", fontsize=10)
    ax_hm.set_xlabel("Cluster (sorted per pair)", fontsize=11)
    ax_hm.set_title(
        f"Cross-Cluster RMSD Overview ({clust_label})\n"
        f"Rows sorted by overall mean RMSD  |  Green = agreement, Red = diversity",
        fontsize=13, fontweight="bold", pad=12,
    )
    fig_hm.tight_layout()
    fig_hm.savefig(output_dir_path / f"rmsd_overview_heatmap_{tag}.png",
                   dpi=200, bbox_inches="tight")
    plt.show()

    # ------------------------------------------------------------------
    # 5. Global RMSD histogram
    # ------------------------------------------------------------------
    fig_hist, ax_hist = plt.subplots(figsize=(9, 5))
    ax_hist.hist(
        df_cluster_rmsd_in["mean_rmsd"], bins=30, color="#636e72",
        edgecolor="white", alpha=0.85, label="All clusters",
    )
    ax_hist.axvline(
        df_cluster_rmsd_in["mean_rmsd"].mean(), color="#d63031",
        linestyle="--", linewidth=2,
        label=f"Mean = {df_cluster_rmsd_in['mean_rmsd'].mean():.2f} Å",
    )
    ax_hist.axvline(
        df_cluster_rmsd_in["mean_rmsd"].median(), color="#0984e3",
        linestyle="--", linewidth=2,
        label=f"Median = {df_cluster_rmsd_in['mean_rmsd'].median():.2f} Å",
    )
    ax_hist.set_xlabel("Mean Cross-Method RMSD (Å)", fontsize=11)
    ax_hist.set_ylabel("Number of Clusters", fontsize=11)
    ax_hist.set_title(
        f"Global RMSD Distribution Across All Multi-Method Clusters ({clust_label})",
        fontsize=13, fontweight="bold",
    )
    ax_hist.legend(fontsize=10)
    ax_hist.grid(axis="y", alpha=0.3)
    fig_hist.tight_layout()
    fig_hist.savefig(output_dir_path / f"rmsd_overview_histogram_{tag}.png",
                     dpi=200, bbox_inches="tight")
    plt.show()

    # ------------------------------------------------------------------
    # 6. Global summary statistics
    # ------------------------------------------------------------------
    print("\n" + "=" * 130)
    print(f"CROSS-CLUSTER RMSD OVERVIEW — GLOBAL SUMMARY  ({clust_label})")
    print("=" * 130)
    print(f"  Protein-ligand pairs analysed       : {len(agg)}")
    print(f"  Multi-method clusters analysed       : {len(df_cluster_rmsd_in)}")
    print(f"  Overall mean RMSD across clusters    : "
          f"{df_cluster_rmsd_in['mean_rmsd'].mean():.3f} Å")
    print(f"  Overall median RMSD across clusters  : "
          f"{df_cluster_rmsd_in['mean_rmsd'].median():.3f} Å")
    print(f"  Overall std RMSD across clusters     : "
          f"{df_cluster_rmsd_in['mean_rmsd'].std():.3f} Å")
    print(f"  Min cluster mean RMSD                : "
          f"{df_cluster_rmsd_in['mean_rmsd'].min():.3f} Å")
    print(f"  Max cluster mean RMSD                : "
          f"{df_cluster_rmsd_in['mean_rmsd'].max():.3f} Å")
    print(f"  Clusters with RMSD <= 2 Å            : "
          f"{(df_cluster_rmsd_in['mean_rmsd'] <= 2.0).sum()}")
    print(f"  Clusters with RMSD 2-5 Å             : "
          f"{((df_cluster_rmsd_in['mean_rmsd'] > 2.0) & (df_cluster_rmsd_in['mean_rmsd'] <= 5.0)).sum()}")
    print(f"  Clusters with RMSD > 5 Å             : "
          f"{(df_cluster_rmsd_in['mean_rmsd'] > 5.0).sum()}")

    pair_mean_cols_ov = [c for c in df_cluster_rmsd_in.columns
                         if c.startswith("mean_") and " vs " in c]
    if pair_mean_cols_ov:
        print(f"\n  Method-pair global mean RMSD:")
        for pc in pair_mean_cols_ov:
            label_pc = pc.replace("mean_", "")
            val = df_cluster_rmsd_in[pc].mean()
            print(f"    {label_pc:<20}: {val:.3f} Å")

    agg.to_csv(output_dir_path / f"rmsd_overview_per_protlig_{tag}.csv", index=False)
    print(f"\n  Saved: {output_dir_path / f'rmsd_overview_per_protlig_{tag}.csv'}")
    print(f"  Saved: {output_dir_path / f'rmsd_overview_boxplot_{tag}.png'}")
    print(f"  Saved: {output_dir_path / f'rmsd_overview_heatmap_{tag}.png'}")
    print(f"  Saved: {output_dir_path / f'rmsd_overview_histogram_{tag}.png'}")

    return agg


# ---------- Run cross-method cluster RMSD for both clustering methods ----------

for _clust_col, _clust_label in [("cluster_kmedoids", "K-Medoids"),
                                  ("cluster_hierarchical", "Hierarchical")]:

    print("\n" + "=" * 120)
    print(f"CROSS-METHOD CLUSTER RMSD ANALYSIS  (clustering: {_clust_label})")
    print("=" * 120)

    _combo_list = sorted(df_clust.groupby(["protein", "ligand"]).groups.keys())
    print(f"Protein-ligand combinations: {len(_combo_list)}")

    cluster_rmsd_rows = []

    for protein, ligand in _combo_list:
        grp = df_clust[(df_clust["protein"] == protein) & (df_clust["ligand"] == ligand)]
        cluster_ids = sorted(grp[_clust_col].unique())

        for cid in cluster_ids:
            c_grp = grp[grp[_clust_col] == cid]

            # --- Load coordinates per method ---
            method_coords = {}  # method -> list of (dist_to_centroid, coords)
            cluster_centroid = c_grp[["centroid_x", "centroid_y", "centroid_z"]].mean().values

            for meth in METHODS_LIST:
                m_sub = c_grp[c_grp["method"] == meth].sort_values("pose_rank")
                entries = []
                for _, rw in m_sub.iterrows():
                    mol = loaded_mols.get(rw["file_path"])
                    if mol is not None:
                        coords = heavy_atom_coords(mol)
                        heavy_centroid = coords.mean(axis=0)
                        dist = float(np.linalg.norm(heavy_centroid - cluster_centroid))
                        entries.append((dist, coords, rw["file_path"], rw["pose_rank"]))
                if entries:
                    entries.sort(key=lambda x: x[0])
                    method_coords[meth] = entries

            active_methods = sorted(method_coords.keys())
            if len(active_methods) < 2:
                continue

            # Balanced N = min available poses across all present methods
            n_balanced = min(len(method_coords[m]) for m in active_methods)
            if n_balanced == 0:
                continue

            # Keep only the N closest poses per method
            selected = {m: method_coords[m][:n_balanced] for m in active_methods}

            # --- Cross-method pairwise RMSD ---
            all_cross_rmsds = []
            pair_rmsds = {}  # (m1, m2) -> list of rmsd values

            for m1, m2 in combinations(active_methods, 2):
                pair_key = f"{_ABBR[m1]} vs {_ABBR[m2]}"
                rmsds = []
                for _, c1, _, _ in selected[m1]:
                    for _, c2, _, _ in selected[m2]:
                        r = _pairwise_rmsd_coords(c1, c2)
                        if not np.isnan(r):
                            rmsds.append(r)
                pair_rmsds[pair_key] = rmsds
                all_cross_rmsds.extend(rmsds)

            if not all_cross_rmsds:
                continue

            # Per-tool pose counts (all poses in cluster, not just selected)
            tool_pose_counts = {}
            for meth in METHODS_LIST:
                tool_pose_counts[meth] = int((c_grp["method"] == meth).sum())

            row_data = {
                "protein": protein,
                "ligand": ligand,
                "cluster": cid,
                "n_methods": len(active_methods),
                "methods": "+".join(_ABBR[m] for m in active_methods),
                "n_poses_per_method": n_balanced,
                "mean_rmsd": round(np.mean(all_cross_rmsds), 3),
                "median_rmsd": round(np.median(all_cross_rmsds), 3),
                "min_rmsd": round(np.min(all_cross_rmsds), 3),
                "max_rmsd": round(np.max(all_cross_rmsds), 3),
                "std_rmsd": round(np.std(all_cross_rmsds), 3),
            }
            # Per-tool pose counts
            for meth in METHODS_LIST:
                row_data[f"n_{_ABBR[meth]}"] = tool_pose_counts[meth]
            # Per-pair breakdown
            for pk, rv in pair_rmsds.items():
                if rv:
                    row_data[f"mean_{pk}"] = round(np.mean(rv), 3)
                    row_data[f"min_{pk}"] = round(np.min(rv), 3)
                else:
                    row_data[f"mean_{pk}"] = np.nan
                    row_data[f"min_{pk}"] = np.nan

            cluster_rmsd_rows.append(row_data)

    df_cluster_rmsd = pd.DataFrame(cluster_rmsd_rows)

    if df_cluster_rmsd.empty:
        print("\nNo multi-method clusters found — skipping heatmaps.")
    else:
        # Print table
        display_cols = ["protein", "ligand", "cluster", "n_methods", "methods",
                        "n_poses_per_method", "mean_rmsd", "median_rmsd",
                        "min_rmsd", "max_rmsd", "std_rmsd"]
        print(f"\nPer-cluster cross-method RMSD ({len(df_cluster_rmsd)} clusters):\n")
        print(df_cluster_rmsd[display_cols].to_string(index=False))

        # ================================================================
        # HEATMAP 1: Per-cluster mean cross-method RMSD + tool presence
        # ================================================================
        df_cluster_rmsd["label"] = (
            df_cluster_rmsd["protein"].str[:22] + " / " +
            df_cluster_rmsd["ligand"].str[:16] + " / C" +
            df_cluster_rmsd["cluster"].astype(str)
        )

        df_hm1 = df_cluster_rmsd.sort_values("mean_rmsd", ascending=True).copy()

        tool_count_cols = [f"n_{_ABBR[m]}" for m in METHODS_LIST]
        tool_col_labels = [_ABBR[m] for m in METHODS_LIST]
        n_tool_cols = len(tool_count_cols)

        tool_vals = df_hm1[tool_count_cols].values.astype(float)
        rmsd_vals = df_hm1[["mean_rmsd"]].values
        n_rows_h1 = len(df_hm1)

        cmap_rg = sns.color_palette("RdYlGn_r", as_cmap=True)

        fig_h1, (ax_tools, ax_rmsd) = plt.subplots(
            1, 2,
            figsize=(max(10, n_tool_cols * 1.8 + 8), max(4, n_rows_h1 * 0.45 + 2)),
            gridspec_kw={"width_ratios": [n_tool_cols, 1], "wspace": 0.05},
            sharey=True,
        )

        # Left panel: tool presence / pose-count matrix
        cmap_tools = sns.color_palette("YlOrRd", as_cmap=True)
        tool_max = max(tool_vals.max(), 1)
        im_tools = ax_tools.imshow(tool_vals, cmap=cmap_tools, aspect="auto",
                                    vmin=0, vmax=tool_max)
        ax_tools.set_xticks(range(n_tool_cols))
        ax_tools.set_xticklabels(tool_col_labels, fontsize=11, fontweight="bold")
        ax_tools.set_yticks(range(n_rows_h1))
        ax_tools.set_yticklabels(df_hm1["label"].values, fontsize=8)
        for i in range(n_rows_h1):
            for j in range(n_tool_cols):
                v = int(tool_vals[i, j])
                txt = str(v) if v > 0 else "—"
                color = "white" if v > tool_max * 0.6 else ("black" if v > 0 else "gray")
                ax_tools.text(j, i, txt, ha="center", va="center",
                              fontsize=9, fontweight="bold", color=color)
        cbar_t = plt.colorbar(im_tools, ax=ax_tools, shrink=0.6, pad=0.02,
                               orientation="horizontal", location="bottom")
        cbar_t.set_label("Poses in Cluster", fontsize=9)
        ax_tools.set_title("Docking Tool Poses", fontsize=12, fontweight="bold")

        # Right panel: mean cross-method RMSD
        im_rmsd = ax_rmsd.imshow(rmsd_vals, cmap=cmap_rg, aspect="auto",
                                  vmin=0, vmax=max(df_hm1["mean_rmsd"].max(), 5))
        ax_rmsd.set_xticks([0])
        ax_rmsd.set_xticklabels(["Mean\nRMSD (Å)"], fontsize=10)
        for i, row in enumerate(df_hm1.itertuples()):
            color = "white" if row.mean_rmsd > rmsd_vals.max() * 0.55 else "black"
            ax_rmsd.text(0, i, f"{row.mean_rmsd:.2f}", ha="center", va="center",
                         fontsize=10, fontweight="bold", color=color)
        cbar_r = plt.colorbar(im_rmsd, ax=ax_rmsd, shrink=0.6, pad=0.02,
                               orientation="horizontal", location="bottom")
        cbar_r.set_label("RMSD (Å)", fontsize=9)
        ax_rmsd.set_title("RMSD", fontsize=12, fontweight="bold")

        _tag = _clust_label.lower().replace("-", "")
        fig_h1.suptitle(
            f"Cross-Method RMSD per Cluster ({_clust_label})\n"
            f"Green = low RMSD (agreement)  |  Red = high RMSD (diversity)",
            fontsize=13, fontweight="bold", y=1.02,
        )
        fig_h1.tight_layout()
        fig_h1.savefig(output_dir / f"cross_method_rmsd_per_cluster_{_tag}.png",
                       dpi=200, bbox_inches="tight")
        plt.show()

        # ================================================================
        # HEATMAP 2: Protein-Ligand × Method-Pair mean RMSD
        # ================================================================
        pair_mean_cols = [c for c in df_cluster_rmsd.columns if c.startswith("mean_") and " vs " in c]

        if pair_mean_cols:
            agg_dict_h2 = {pc: "mean" for pc in pair_mean_cols}
            agg_dict_h2["mean_rmsd"] = "mean"
            agg_dict_h2["cluster"] = "count"
            df_pl_rmsd = (
                df_cluster_rmsd.groupby(["protein", "ligand"])
                .agg(agg_dict_h2)
                .rename(columns={"cluster": "n_clusters"})
                .reset_index()
            )
            df_pl_rmsd["combo"] = (
                df_pl_rmsd["protein"].str[:22] + "\n" + df_pl_rmsd["ligand"].str[:16]
            )

            pair_labels = [c.replace("mean_", "") for c in pair_mean_cols]
            heat_matrix = df_pl_rmsd[pair_mean_cols].values
            row_labels = df_pl_rmsd["combo"].values

            fig_h2, ax_h2 = plt.subplots(
                figsize=(max(7, len(pair_labels) * 2.5 + 2),
                         max(4, len(row_labels) * 0.7 + 2))
            )
            im2 = ax_h2.imshow(heat_matrix, cmap=cmap_rg, aspect="auto",
                                vmin=0, vmax=max(np.nanmax(heat_matrix), 5))
            ax_h2.set_xticks(range(len(pair_labels)))
            ax_h2.set_xticklabels(pair_labels, fontsize=11, fontweight="bold")
            ax_h2.set_yticks(range(len(row_labels)))
            ax_h2.set_yticklabels(row_labels, fontsize=9)
            for i in range(heat_matrix.shape[0]):
                for j in range(heat_matrix.shape[1]):
                    v = heat_matrix[i, j]
                    if np.isnan(v):
                        ax_h2.text(j, i, "—", ha="center", va="center", fontsize=10, color="gray")
                    else:
                        color = "white" if v > np.nanmax(heat_matrix) * 0.55 else "black"
                        ax_h2.text(j, i, f"{v:.2f}", ha="center", va="center",
                                   fontsize=10, fontweight="bold", color=color)
            cbar2 = plt.colorbar(im2, ax=ax_h2, shrink=0.8, pad=0.02)
            cbar2.set_label("Mean RMSD (Å)", fontsize=10)
            ax_h2.set_title(
                f"Mean Cross-Method RMSD per Protein-Ligand Pair ({_clust_label})\n"
                f"Averaged across multi-method clusters",
                fontsize=13, fontweight="bold", pad=12,
            )
            fig_h2.tight_layout()
            fig_h2.savefig(output_dir / f"cross_method_rmsd_per_protlig_{_tag}.png",
                           dpi=200, bbox_inches="tight")
            plt.show()
        else:
            print("  (Only one method pair found — skipping pair-level heatmap)")

        # ================================================================
        # HEATMAP 3: Full cluster × (tool presence + method-pair RMSD)
        # ================================================================
        if pair_mean_cols and len(df_cluster_rmsd) > 1:
            df_hm3 = df_cluster_rmsd.sort_values(["protein", "ligand", "cluster"]).copy()
            labels_h3 = df_hm3["label"].values
            n_rows_h3 = len(labels_h3)

            tool_matrix_h3 = df_hm3[tool_count_cols].values.astype(float)
            rmsd_matrix_h3 = df_hm3[pair_mean_cols].values
            n_pair_cols = len(pair_labels)

            fig_h3, (ax3_tools, ax3_rmsd) = plt.subplots(
                1, 2,
                figsize=(max(10, (n_tool_cols + n_pair_cols) * 2 + 4),
                         max(5, n_rows_h3 * 0.45 + 2)),
                gridspec_kw={"width_ratios": [n_tool_cols, n_pair_cols], "wspace": 0.08},
                sharey=True,
            )

            tool_max_h3 = max(tool_matrix_h3.max(), 1)
            im3_t = ax3_tools.imshow(tool_matrix_h3, cmap=cmap_tools, aspect="auto",
                                      vmin=0, vmax=tool_max_h3)
            ax3_tools.set_xticks(range(n_tool_cols))
            ax3_tools.set_xticklabels(tool_col_labels, fontsize=10, fontweight="bold")
            ax3_tools.set_yticks(range(n_rows_h3))
            ax3_tools.set_yticklabels(labels_h3, fontsize=8)
            for i in range(n_rows_h3):
                for j in range(n_tool_cols):
                    v = int(tool_matrix_h3[i, j])
                    txt = str(v) if v > 0 else "—"
                    color = "white" if v > tool_max_h3 * 0.6 else ("black" if v > 0 else "gray")
                    ax3_tools.text(j, i, txt, ha="center", va="center",
                                   fontsize=9, fontweight="bold", color=color)
            cbar3_t = plt.colorbar(im3_t, ax=ax3_tools, shrink=0.5, pad=0.02,
                                    orientation="horizontal", location="bottom")
            cbar3_t.set_label("Poses", fontsize=9)
            ax3_tools.set_title("Docking Tools", fontsize=12, fontweight="bold")

            im3_r = ax3_rmsd.imshow(rmsd_matrix_h3, cmap=cmap_rg, aspect="auto",
                                    vmin=0, vmax=max(np.nanmax(rmsd_matrix_h3), 5))
            ax3_rmsd.set_xticks(range(n_pair_cols))
            ax3_rmsd.set_xticklabels(pair_labels, fontsize=10, fontweight="bold")
            for i in range(n_rows_h3):
                for j in range(n_pair_cols):
                    v = rmsd_matrix_h3[i, j]
                    if np.isnan(v):
                        ax3_rmsd.text(j, i, "—", ha="center", va="center", fontsize=9, color="gray")
                    else:
                        color = "white" if v > np.nanmax(rmsd_matrix_h3) * 0.55 else "black"
                        ax3_rmsd.text(j, i, f"{v:.2f}", ha="center", va="center",
                                      fontsize=9, fontweight="bold", color=color)
            cbar3_r = plt.colorbar(im3_r, ax=ax3_rmsd, shrink=0.5, pad=0.02,
                                    orientation="horizontal", location="bottom")
            cbar3_r.set_label("Mean RMSD (Å)", fontsize=9)
            ax3_rmsd.set_title("Cross-Method RMSD", fontsize=12, fontweight="bold")

            fig_h3.suptitle(
                f"Cluster Detail: Tool Presence + Cross-Method RMSD ({_clust_label})\n"
                f"Green = high pose agreement  |  Red = high pose diversity",
                fontsize=13, fontweight="bold", y=1.02,
            )
            fig_h3.tight_layout()
            fig_h3.savefig(output_dir / f"cross_method_rmsd_cluster_x_pair_{_tag}.png",
                           dpi=200, bbox_inches="tight")
            plt.show()

        # ================================================================
        # SUMMARY STATISTICS
        # ================================================================
        print("\n" + "=" * 120)
        print(f"CROSS-METHOD CLUSTER RMSD — SUMMARY  ({_clust_label})")
        print("=" * 120)
        print(f"  Clustering method : {_clust_label}")
        print(f"  Multi-method clusters analysed : {len(df_cluster_rmsd)}")
        print(f"  Protein-ligand combinations    : {df_cluster_rmsd.groupby(['protein','ligand']).ngroups}")
        print(f"\n  Overall cross-method RMSD across all clusters:")
        print(f"    Mean   : {df_cluster_rmsd['mean_rmsd'].mean():.3f} Å")
        print(f"    Median : {df_cluster_rmsd['median_rmsd'].median():.3f} Å")
        print(f"    Min    : {df_cluster_rmsd['min_rmsd'].min():.3f} Å")
        print(f"    Max    : {df_cluster_rmsd['max_rmsd'].max():.3f} Å")

        low_rmsd = df_cluster_rmsd[df_cluster_rmsd["mean_rmsd"] <= 2.0]
        high_rmsd = df_cluster_rmsd[df_cluster_rmsd["mean_rmsd"] > 5.0]
        print(f"\n  Clusters with mean RMSD <= 2.0 Å (high agreement) : {len(low_rmsd)}")
        print(f"  Clusters with mean RMSD >  5.0 Å (high diversity)  : {len(high_rmsd)}")

        if not low_rmsd.empty:
            print(f"\n  Best-agreement clusters (lowest RMSD):")
            for _, r in low_rmsd.nsmallest(5, "mean_rmsd").iterrows():
                print(f"    {r['protein']} + {r['ligand']} C{r['cluster']}  "
                      f"({r['methods']}, N={r['n_poses_per_method']})  "
                      f"mean={r['mean_rmsd']:.3f} Å")

        if not high_rmsd.empty:
            print(f"\n  Most-diverse clusters (highest RMSD):")
            for _, r in high_rmsd.nlargest(5, "mean_rmsd").iterrows():
                print(f"    {r['protein']} + {r['ligand']} C{r['cluster']}  "
                      f"({r['methods']}, N={r['n_poses_per_method']})  "
                      f"mean={r['mean_rmsd']:.3f} Å")

        # Save per-cluster CSV
        df_cluster_rmsd.to_csv(output_dir / f"cross_method_cluster_rmsd_{_tag}.csv", index=False)
        print(f"\n  Saved: {output_dir / f'cross_method_cluster_rmsd_{_tag}.csv'}")

        # ================================================================
        # RMSD OVERVIEW ACROSS CLUSTERS
        # ================================================================
        rmsd_overview_across_clusters(
            df_cluster_rmsd_in=df_cluster_rmsd,
            methods_list=METHODS_LIST,
            clust_label=_clust_label,
            output_dir_path=output_dir,
        )

print("\n" + "=" * 100)
print("CROSS-METHOD CLUSTER RMSD ANALYSIS COMPLETE")
print("=" * 100)

# %% [markdown]
# # K-Medoids vs Hierarchical Comparison

# %%
# ============================================================================
# K-MEDOIDS vs HIERARCHICAL AGREEMENT
# ============================================================================

print(f"\n{'='*100}")
print(f"K-MEDOIDS vs HIERARCHICAL AGREEMENT  (mode={DISTANCE_MODE})")
print(f"{'='*100}")

from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

agree_rows = []
for (protein, ligand), grp in df_clust.groupby(["protein", "ligand"]):
    lab_km = grp["cluster_kmedoids"].values
    lab_hi = grp["cluster_hierarchical"].values
    ari = adjusted_rand_score(lab_km, lab_hi) if len(set(lab_km)) > 1 and len(set(lab_hi)) > 1 else np.nan
    nmi = normalized_mutual_info_score(lab_km, lab_hi) if len(set(lab_km)) > 1 and len(set(lab_hi)) > 1 else np.nan
    agree_rows.append({
        "protein": protein, "ligand": ligand,
        "k_kmedoids": grp["cluster_kmedoids"].nunique(),
        "k_hierarchical": grp["cluster_hierarchical"].nunique(),
        "ARI": round(ari, 3) if np.isfinite(ari) else "N/A",
        "NMI": round(nmi, 3) if np.isfinite(nmi) else "N/A",
        "n_poses": len(grp),
    })

df_agree = pd.DataFrame(agree_rows)
print(df_agree.to_string(index=False))

# %% [markdown]
# # Save All Results

# %%
# ============================================================================
# SAVE ALL OUTPUT
# ============================================================================

# Poses with cluster labels
csv_poses = output_dir / f"clustering_{DISTANCE_MODE}_all_poses.csv"
df_clust.to_csv(csv_poses, index=False)
print(f"  Saved: {csv_poses}")

# K-Medoids pocket stats
csv_km = output_dir / f"clustering_{DISTANCE_MODE}_kmedoids_pockets.csv"
df_km_overlap.to_csv(csv_km, index=False)
print(f"  Saved: {csv_km}")

# Hierarchical pocket stats
csv_hi = output_dir / f"clustering_{DISTANCE_MODE}_hierarchical_pockets.csv"
df_hi_overlap.to_csv(csv_hi, index=False)
print(f"  Saved: {csv_hi}")

# Conformity scores
csv_scores = output_dir / f"clustering_{DISTANCE_MODE}_conformity_scores.csv"
df_scores.to_csv(csv_scores, index=False)
print(f"  Saved: {csv_scores}")

# Cross-method RMSD
if not df_xrmsd.empty:
    csv_xrmsd = output_dir / f"clustering_{DISTANCE_MODE}_cross_method_rmsd.csv"
    df_xrmsd.to_csv(csv_xrmsd, index=False)
    print(f"  Saved: {csv_xrmsd}")

# Agreement
csv_agree = output_dir / f"clustering_{DISTANCE_MODE}_agreement.csv"
df_agree.to_csv(csv_agree, index=False)
print(f"  Saved: {csv_agree}")

# ── Summary JSON ──
import json
summary = {
    "distance_mode": DISTANCE_MODE,
    "hybrid_weight_rmsd": HYBRID_WEIGHT_RMSD if DISTANCE_MODE == "hybrid" else None,
    "methods": METHODS_LIST,
    "n_combos": n_combos,
    "total_poses": len(df_clust),
    "kmedoids": {ck: {"k": v["optimal_k"], "silhouette": round(v["best_silhouette"], 3)}
                 for ck, v in km_info.items()},
    "hierarchical": {ck: {"k": v["optimal_k"], "silhouette": round(v["best_silhouette"], 3)}
                     for ck, v in hi_info.items()},
}
json_path = output_dir / f"clustering_{DISTANCE_MODE}_summary.json"
with open(json_path, "w") as f:
    json.dump(summary, f, indent=2)
print(f"  Saved: {json_path}")

print(f"\n{'='*100}")
print("ALL DONE")
print(f"{'='*100}")
