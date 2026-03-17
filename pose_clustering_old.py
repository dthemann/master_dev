# %% [markdown]
# # Setting Filepaths

# %%
# Paths to docked pose folders (using Linux paths)

wd = "/home/manndo/master_dev/posebuster_proved/"
autodock_poses_dir = wd + "/autodock"
diffdock_poses_dir = wd + "/diffdock"
equibind_guided_poses_dir = wd + "/equibind_guided"
equibind_exclusion_poses_dir = wd + "/equibind_exclusion"

docking_directories = {
    "autodock":            autodock_poses_dir,
    "diffdock":            diffdock_poses_dir,
    "equibind_guided":     equibind_guided_poses_dir
}

# Optional: protein files for intermolecular distance checks
proteins_dir = wd + "/Orai"

# %% [markdown]
# # Setup Clustering Data 

# %%
"""
Binding Pocket Clustering Analysis — Shared Setup
===================================================
Collect all docked poses and extract ligand centroids.
Only includes poses that PASSED ALL PoseBusters tests.
Shared functions used by both K-Medoids and Ward clustering cells.
"""

import numpy as np
import pandas as pd
import math
from pathlib import Path
from collections import defaultdict
from itertools import combinations
import re
import warnings
warnings.filterwarnings('ignore')

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("ERROR: RDKit is required. Install with: conda install -c conda-forge rdkit")

try:
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
    from sklearn.preprocessing import StandardScaler
    from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
    from scipy.spatial.distance import pdist, squareform
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False
    print("ERROR: scikit-learn and scipy are required. Install with: conda install scikit-learn scipy")

try:
    from sklearn_extra.cluster import KMedoids
    KMEDOIDS_AVAILABLE = True
except ImportError:
    KMEDOIDS_AVAILABLE = False
    print("WARNING: sklearn_extra not found. Install with: pip install scikit-learn-extra")
    print("         Will use a simple K-Medoids fallback implementation.")

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches

# ============================================================================
# DERIVE METHODS LIST FROM THE DICT
# ============================================================================
METHODS_LIST = list(docking_directories.keys())

# ============================================================================
# CONFIGURATION
# ============================================================================
output_dir = Path(wd + "/posebusters_results")
converted_dir = output_dir / "converted_pdbqt"
output_dir.mkdir(exist_ok=True)

# PoseBusters-proved pose directories  (wd already IS posebuster_proved)
posebuster_proved_dir = Path(wd)

K_RANGE = range(2, 11)
RANDOM_STATE = 42

# ============================================================================
# LOAD POSEBUSTERS-PROVED MANIFEST
# ============================================================================

manifest_path = posebuster_proved_dir / "posebuster_proved_manifest.csv"
pb_results_path = Path(wd).parent / "posebusters_results" / "posebusters_filtered_results.csv"

# Build a set of file paths that passed all PoseBusters tests
proved_files = set()

if manifest_path.exists():
    df_manifest = pd.read_csv(manifest_path)
    if "pose_file" in df_manifest.columns:
        proved_files = set(df_manifest["pose_file"].dropna().astype(str))
    elif df_manifest.columns[0]:
        # First column is likely the file path
        proved_files = set(df_manifest.iloc[:, 0].dropna().astype(str))
    print(f"Loaded {len(proved_files)} PoseBusters-proved file paths from manifest.")
elif pb_results_path.exists():
    # Fallback: recompute from the full results CSV
    df_pb = pd.read_csv(pb_results_path)
    
    # Identify boolean test columns (same logic as PoseBusters.ipynb)
    exclude_cols = {
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
    metadata_cols = {
        'file_path', 'filepath', 'file', 'path', 'sdf_file', 'sdf_path',
        'method', 'docking_method', 'protein', 'ligand', 'pose_rank', 'rank',
        'molecule', 'mol_name', 'name', 'complex', 'protein_path', 'ligand_path',
        'mol_pred', 'mol_true', 'mol_cond',
    }
    
    test_cols = []
    for col in df_pb.columns:
        col_lower = col.lower().strip()
        if col_lower in metadata_cols or col in exclude_cols:
            continue
        if col_lower.startswith('number_') or col_lower.startswith('num_'):
            continue
        unique_vals = df_pb[col].dropna().unique()
        if set(unique_vals).issubset({True, False, 1, 0, 1.0, 0.0, 'True', 'False', 'true', 'false'}):
            test_cols.append(col)
    
    # Convert to bool and filter
    for tc in test_cols:
        if df_pb[tc].dtype == object:
            df_pb[tc] = df_pb[tc].map({'True': True, 'true': True, 'False': False, 'false': False})
        df_pb[tc] = df_pb[tc].astype(bool)
    
    df_pb["all_passed"] = df_pb[test_cols].all(axis=1)
    df_passed = df_pb[df_pb["all_passed"]]
    
    if "pose_file" in df_passed.columns:
        proved_files = set(df_passed["pose_file"].dropna().astype(str))
    print(f"Computed {len(proved_files)} PoseBusters-proved file paths from results CSV.")
else:
    print("WARNING: No PoseBusters manifest or results found!")
    print(f"  Looked for: {manifest_path}")
    print(f"  Looked for: {pb_results_path}")
    print("  Clustering will proceed with NO filtering (all poses).")

# Also build set of proved filenames (basename only) for directory-based matching
proved_filenames = {Path(f).name for f in proved_files}

print(f"  Unique proved filenames: {len(proved_filenames)}")


# ============================================================================
# SHARED FUNCTIONS
# ============================================================================

def get_ligand_centroid_from_sdf(sdf_file: str):
    """Extract heavy-atom centroid from an SDF file. Returns (x, y, z) or None."""
    try:
        mol = Chem.MolFromMolFile(sdf_file, removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(sdf_file, removeHs=False, sanitize=False)
            if mol is not None:
                mol = Chem.RemoveHs(mol)
        if mol is None or mol.GetNumConformers() == 0:
            return None
        conf = mol.GetConformer()
        coords = np.array([
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])
        return coords.mean(axis=0)
    except Exception:
        return None


def is_posebuster_proved(file_path: str) -> bool:
    """Check whether a pose file passed all PoseBusters tests."""
    if not proved_files and not proved_filenames:
        return True  # No filter available — include everything
    
    # Check exact path match
    if file_path in proved_files:
        return True
    
    # Check by resolved absolute path
    abs_path = str(Path(file_path).resolve())
    if abs_path in proved_files:
        return True
    
    # Check by filename match (handles copied files in posebuster_proved/)
    fname = Path(file_path).name
    if fname in proved_filenames:
        return True
    
    return False


def _parse_pose_info(sdf_file, base_dir, method_name):
    """
    Extract (protein, ligand, rank) from an SDF file path.
    
    Handles three patterns:
      1. Subdir pattern: SDFs inside directories named {A}__{B}
         → ligand = A, protein = B  (DiffDock / EquiBind convention)
      2. Flat pattern: SDF files with {A}__{B}_vina_out_model{N} naming
         → protein = A, ligand = B  (AutoDock converted files)
      3. Prefixed flat: {protein}__{ligand}__{original_name}.sdf
         → protein = parts[0], ligand = parts[1]  (renamed by pose_busters.py)
    """
    protein = ligand = None
    rank = 1

    try:
        rel_path = sdf_file.relative_to(base_dir)
        parts_list = list(rel_path.parts)
    except ValueError:
        parts_list = [sdf_file.name]

    # --- Strategy 1: parent directory named {A}__{B} ---
    if len(parts_list) > 1:
        parent_name = parts_list[0]
        # Skip special directories
        if parent_name in ("prepared_proteins", "converted_pdbqt"):
            return None
        if "_pymol" in parent_name:
            return None
        dir_parts = parent_name.split("__")
        if len(dir_parts) == 2:
            ligand, protein = dir_parts[0], dir_parts[1]
            # Strip trailing _cleaned from protein name
            protein = re.sub(r'_cleaned$', '', protein)

    # --- Strategy 2: flat file with _vina_out_model in stem ---
    if protein is None or ligand is None:
        stem = sdf_file.stem
        if "_vina_out_model" in stem:
            base = stem.split("_vina_out_model")[0]
            model_str = stem.split("_vina_out_model")[-1]
            fp = base.split("__")
            if len(fp) == 2:
                protein, ligand = fp[0], fp[1]
                rank = int(model_str) if model_str.isdigit() else 0

    # --- Strategy 3: flat file with {protein}__{ligand}__{rest} naming ---
    # Handles DiffDock/EquiBind files renamed by pose_busters.py
    if protein is None or ligand is None:
        stem = sdf_file.stem
        parts = stem.split("__")
        if len(parts) >= 3:
            protein, ligand = parts[0], parts[1]
            # Strip trailing _cleaned from protein name
            protein = re.sub(r'_cleaned$', '', protein)
            # Extract rank from the remainder
            rest = "__".join(parts[2:])
            rank_match = re.search(r'rank(\d+)', rest)
            if rank_match:
                rank = int(rank_match.group(1))
            else:
                model_match = re.search(r'model(\d+)', rest)
                if model_match:
                    rank = int(model_match.group(1))

    # --- Extract rank from filename if not already set ---
    if protein is not None and rank == 1:
        rank_match = re.search(r'rank(\d+)', sdf_file.stem)
        if rank_match:
            rank = int(rank_match.group(1))

    if protein is None or ligand is None:
        return None

    return protein, ligand, rank


def collect_proved_pose_centroids():
    """
    Collect centroids ONLY for poses that passed all PoseBusters tests.
    Iterates over docking_directories dict to collect poses from each method.
    """
    all_poses = []
    n_skipped = {m: 0 for m in docking_directories}
    n_included = {m: 0 for m in docking_directories}

    for method_name, dir_path in docking_directories.items():
        dir_path = Path(dir_path)
        if not dir_path.exists():
            print(f"  WARNING: {method_name} directory not found: {dir_path}")
            continue

        # Collect all SDF files recursively
        sdf_files = sorted(dir_path.rglob("*.sdf"))
        if not sdf_files:
            print(f"  WARNING: No SDF files found for {method_name} in: {dir_path}")
            continue

        for sdf_file in sdf_files:
            parsed = _parse_pose_info(sdf_file, dir_path, method_name)
            if parsed is None:
                continue
            protein, ligand, rank = parsed

            if not is_posebuster_proved(str(sdf_file)):
                n_skipped[method_name] += 1
                continue

            centroid = get_ligand_centroid_from_sdf(str(sdf_file))
            if centroid is not None:
                n_included[method_name] += 1
                all_poses.append({
                    "method": method_name, "protein": protein,
                    "ligand": ligand, "pose_rank": rank,
                    "file_path": str(sdf_file),
                    "centroid_x": centroid[0], "centroid_y": centroid[1],
                    "centroid_z": centroid[2],
                })

    print(f"\n  PoseBusters filtering summary:")
    for method_name in docking_directories:
        total = n_included[method_name] + n_skipped[method_name]
        print(f"    {method_name}: {n_included[method_name]} included / "
              f"{n_skipped[method_name]} skipped (of {total} total)")

    return pd.DataFrame(all_poses)


def compute_cluster_stats(df, cluster_col):
    """Compute cluster centers, sizes, and radii for a given cluster column."""
    stats = []
    for (protein, ligand, cluster_id), group in df.groupby(["protein", "ligand", cluster_col]):
        coords = group[["centroid_x", "centroid_y", "centroid_z"]].values
        center = coords.mean(axis=0)
        radius = np.max(np.linalg.norm(coords - center, axis=1)) if len(coords) > 1 else 0.0
        method_counts = group["method"].value_counts().to_dict()
        row_data = {
            "protein": protein, "ligand": ligand, "cluster_id": cluster_id,
            "cluster_method": cluster_col.replace("cluster_", ""),
            "n_poses": len(group),
            "center_x": center[0], "center_y": center[1], "center_z": center[2],
            "radius_A": round(radius, 2),
        }
        # Add per-method counts dynamically
        for m in METHODS_LIST:
            row_data[m] = method_counts.get(m, 0)
        row_data["n_methods"] = len(method_counts)
        row_data["methods_present"] = ", ".join(sorted(method_counts.keys()))
        stats.append(row_data)
    return pd.DataFrame(stats)


def identify_pocket_overlap(df_cluster_stats):
    """Identify overlapping (multi-method) and unique (single-method) pockets."""
    overlap_data = []
    for _, row in df_cluster_stats.iterrows():
        methods_in_cluster = [m for m in METHODS_LIST if row.get(m, 0) > 0]

        if len(methods_in_cluster) >= 2:
            pocket_type = "Overlapping"
        elif len(methods_in_cluster) == 1:
            pocket_type = f"Unique ({methods_in_cluster[0]})"
        else:
            pocket_type = "Empty"

        overlap_data.append({
            **row.to_dict(),
            "pocket_type": pocket_type,
            "n_methods_in_pocket": len(methods_in_cluster),
        })
    return pd.DataFrame(overlap_data)


# ============================================================================
# COLLECT POSEBUSTERS-PROVED POSE CENTROIDS ONLY
# ============================================================================

print("=" * 100)
print("Collecting ligand centroids from PoseBusters-PROVED poses only...")
print("=" * 100)

df_poses = collect_proved_pose_centroids()

if df_poses.empty:
    print("ERROR: No PoseBusters-proved poses found!")
else:
    print(f"\n  Total proved poses collected: {len(df_poses)}")
    print(f"\n  Proved poses by docking method:")
    for method, count in df_poses["method"].value_counts().items():
        print(f"    {method}: {count}")
    print(f"\n  Unique protein-ligand combinations: "
          f"{df_poses.groupby(['protein', 'ligand']).ngroups}")
    
    
# ====================================================================
# REMOVE CENTROID OUTLIERS
# ====================================================================
# DiffDock low-confidence poses can be placed thousands of Å away from
# the protein. Filter per protein-ligand group: any pose whose centroid
# is > MAX_CENTROID_DIST Å from the group median is dropped.

MAX_CENTROID_DIST = 100  # Å
n_before = len(df_poses)
keep_mask = pd.Series(True, index=df_poses.index)

for (protein, ligand), group in df_poses.groupby(["protein", "ligand"]):
    median_xyz = group[["centroid_x", "centroid_y", "centroid_z"]].median().values
    dists = np.linalg.norm(
        group[["centroid_x", "centroid_y", "centroid_z"]].values - median_xyz,
        axis=1
    )
    outliers = dists > MAX_CENTROID_DIST
    if outliers.any():
        keep_mask.loc[group.index[outliers]] = False
        for i, idx_out in enumerate(group.index[outliers]):
            row = df_poses.loc[idx_out]
            print(f"    OUTLIER REMOVED: {row['method']} rank {row['pose_rank']} "
                  f"({row['protein']} + {row['ligand']}), "
                  f"dist from median = {dists[np.where(outliers)[0][i]]:.0f} Å")
df_poses = df_poses[keep_mask].reset_index(drop=True)
n_removed = n_before - len(df_poses)

if n_removed > 0:
    print(f"\n  Outlier filter (>{MAX_CENTROID_DIST} Å from median): "
          f"removed {n_removed} poses")

# %%
"""
CELL: K-Medoids Clustering of Binding Pockets
================================================
Cluster all docked poses using K-Medoids.
- Automatically selects optimal k via silhouette score
- Shows poses per docking method in each pocket cluster
- Identifies overlapping vs unique binding pockets
"""

# ============================================================================
# K-MEDOIDS FALLBACK IMPLEMENTATION
# ============================================================================

class SimpleKMedoids:
    """Simple K-Medoids (PAM-like) implementation as fallback."""
    def __init__(self, n_clusters=3, max_iter=300, random_state=42):
        self.n_clusters = n_clusters
        self.max_iter = max_iter
        self.random_state = random_state
        self.labels_ = None
        self.cluster_centers_ = None
        self.medoid_indices_ = None
        self.inertia_ = None

    def fit(self, X):
        np.random.seed(self.random_state)
        n_samples = X.shape[0]
        dist_matrix = squareform(pdist(X, metric='euclidean'))
        medoid_indices = np.random.choice(n_samples, self.n_clusters, replace=False)

        for _ in range(self.max_iter):
            labels = np.argmin(dist_matrix[:, medoid_indices], axis=1)
            new_medoid_indices = medoid_indices.copy()
            for k in range(self.n_clusters):
                cluster_indices = np.where(labels == k)[0]
                if len(cluster_indices) == 0:
                    continue
                cluster_dists = dist_matrix[np.ix_(cluster_indices, cluster_indices)]
                best_idx = cluster_indices[np.argmin(cluster_dists.sum(axis=1))]
                new_medoid_indices[k] = best_idx
            if np.array_equal(medoid_indices, new_medoid_indices):
                break
            medoid_indices = new_medoid_indices

        self.medoid_indices_ = medoid_indices
        self.labels_ = np.argmin(dist_matrix[:, medoid_indices], axis=1)
        self.cluster_centers_ = X[medoid_indices]
        self.inertia_ = sum(
            dist_matrix[i, medoid_indices[self.labels_[i]]] for i in range(n_samples)
        )
        return self

    def fit_predict(self, X):
        self.fit(X)
        return self.labels_


# ============================================================================
# K-MEDOIDS: FIND OPTIMAL K
# ============================================================================

def find_optimal_k_kmedoids(coords, k_range, random_state=42):
    """Find optimal k for K-Medoids using silhouette score."""
    results = {"k": [], "silhouette": [], "inertia": []}
    for k in k_range:
        if k >= len(coords):
            break
        if KMEDOIDS_AVAILABLE:
            model = KMedoids(n_clusters=k, metric='euclidean',
                             random_state=random_state, method='pam')
        else:
            model = SimpleKMedoids(n_clusters=k, random_state=random_state)
        labels = model.fit_predict(coords)
        inertia = model.inertia_
        sil = silhouette_score(coords, labels) if len(set(labels)) > 1 else -1
        results["k"].append(k)
        results["silhouette"].append(sil)
        results["inertia"].append(inertia)
    best_idx = np.argmax(results["silhouette"])
    results["optimal_k"] = results["k"][best_idx]
    results["best_silhouette"] = results["silhouette"][best_idx]
    return results


# ============================================================================
# DYNAMIC COLORS / MARKERS — derived from docking_directories
# ============================================================================
_PALETTE = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6", "#f39c12", "#1abc9c"]
_MARKERS = ["o", "^", "s", "D", "v", "P"]

METHOD_COLORS  = {m: _PALETTE[i % len(_PALETTE)] for i, m in enumerate(METHODS_LIST)}
METHOD_MARKERS = {m: _MARKERS[i % len(_MARKERS)] for i, m in enumerate(METHODS_LIST)}


# ============================================================================
# K-MEDOIDS: CLUSTER PER PROTEIN-LIGAND COMBINATION
# ============================================================================

if df_poses.empty:
    print("ERROR: No poses available. Run the data collection cell first.")
else:
    print("=" * 100)
    print("K-MEDOIDS CLUSTERING OF BINDING POCKETS")
    print("=" * 100)

    df_km = df_poses.copy()
    df_km["cluster_kmedoids"] = -1
    km_cluster_info = {}

    for (protein, ligand), group in df_km.groupby(["protein", "ligand"]):
        combo_key = f"{protein}__{ligand}"
        idx = group.index
        coords = group[["centroid_x", "centroid_y", "centroid_z"]].values
        n_poses = len(coords)

        if n_poses < 3:
            df_km.loc[idx, "cluster_kmedoids"] = 0
            km_cluster_info[combo_key] = {
                "n_poses": n_poses, "optimal_k": 1,
                "best_silhouette": -1, "search_results": None,
            }
            continue

        max_k = min(max(K_RANGE), n_poses - 1)
        actual_k_range = range(2, max_k + 1)

        if len(actual_k_range) == 0:
            df_km.loc[idx, "cluster_kmedoids"] = 0
            km_cluster_info[combo_key] = {
                "n_poses": n_poses, "optimal_k": 1,
                "best_silhouette": -1, "search_results": None,
            }
            continue

        # Find optimal k
        km_results = find_optimal_k_kmedoids(coords, actual_k_range,
                                              random_state=RANDOM_STATE)
        k_opt = km_results["optimal_k"]

        # Final fit with optimal k
        if KMEDOIDS_AVAILABLE:
            km_model = KMedoids(n_clusters=k_opt, metric='euclidean',
                                random_state=RANDOM_STATE, method='pam')
        else:
            km_model = SimpleKMedoids(n_clusters=k_opt, random_state=RANDOM_STATE)

        labels = km_model.fit_predict(coords)
        df_km.loc[idx, "cluster_kmedoids"] = labels

        km_cluster_info[combo_key] = {
            "n_poses": n_poses, "optimal_k": k_opt,
            "best_silhouette": km_results["best_silhouette"],
            "search_results": km_results,
        }

    df_km["cluster_kmedoids"] = df_km["cluster_kmedoids"].astype(int)

    # ========================================================================
    # RESULTS
    # ========================================================================

    # --- Optimal k per combination ---
    print(f"\n  {'Protein-Ligand':<55} {'Optimal k':>10} {'Silhouette':>12} {'Poses':>8}")
    print("  " + "-" * 90)
    for combo_key, info in sorted(km_cluster_info.items()):
        print(f"  {combo_key:<55} {info['optimal_k']:>10} "
              f"{info['best_silhouette']:>12.3f} {info['n_poses']:>8}")

    # --- Cluster statistics table ---
    df_km_stats = compute_cluster_stats(df_km, "cluster_kmedoids")
    df_km_overlap = identify_pocket_overlap(df_km_stats)

    print("\n" + "=" * 100)
    print("K-MEDOIDS — Binding Pocket Clusters (poses per method)")
    print("=" * 100)

    km_display_cols = ["protein", "ligand", "cluster_id", "n_poses", "radius_A"] + \
                      METHODS_LIST + ["pocket_type", "n_methods_in_pocket"]
    km_display = df_km_overlap[km_display_cols].sort_values(
        ["protein", "ligand", "cluster_id"])
    print(km_display.to_string(index=False))

    # --- Overlapping vs unique ---
    print("\n" + "=" * 100)
    print("K-MEDOIDS — Overlapping vs Unique Binding Pockets")
    print("=" * 100)

    agg_dict = {"n_poses": ["count", "sum"]}
    for m in METHODS_LIST:
        agg_dict[m] = "sum"
    pocket_summary_km = df_km_overlap.groupby("pocket_type").agg(agg_dict)
    pocket_summary_km.columns = (
        ["N_Clusters", "Total_Poses"] +
        [f"{m}_Poses" for m in METHODS_LIST]
    )
    pocket_summary_km = pocket_summary_km.reset_index()
    print(pocket_summary_km.to_string(index=False))

    # --- Aggregated per protein-ligand ---
    print("\n" + "=" * 100)
    print("K-MEDOIDS — Aggregated per Protein-Ligand")
    print("=" * 100)

    def _agg_func_km(g):
        result = {
            "n_clusters": len(g),
            "overlapping_pockets": (g["pocket_type"] == "Overlapping").sum(),
            "unique_pockets": (g["n_methods_in_pocket"] == 1).sum(),
            "total_poses": g["n_poses"].sum(),
        }
        for m in METHODS_LIST:
            result[f"{m}_in_overlapping"] = g.loc[
                g["pocket_type"] == "Overlapping", m
            ].sum()
        return pd.Series(result)

    agg_km = df_km_overlap.groupby(["protein", "ligand"]).apply(_agg_func_km).reset_index()
    int_cols = [c for c in agg_km.columns if c not in ["protein", "ligand"]]
    for c in int_cols:
        agg_km[c] = agg_km[c].astype(int)
    print(agg_km.to_string(index=False))

    # --- Silhouette score plots ---
    combos_with_search = [(k, v) for k, v in km_cluster_info.items()
                          if v["search_results"] is not None]

    if combos_with_search:
        n_plots = len(combos_with_search)
        n_cols_sil = 4
        n_rows_sil = (n_plots + n_cols_sil - 1) // n_cols_sil
        fig, axes = plt.subplots(n_rows_sil, n_cols_sil,
                                  figsize=(6 * n_cols_sil, 4 * n_rows_sil), squeeze=False)
        fig.suptitle("K-Medoids: Silhouette Score vs k (optimal k selection)",
                     fontsize=14, fontweight='bold')

        for idx, (combo_key, info) in enumerate(combos_with_search):
            row_i, col_i = divmod(idx, n_cols_sil)
            ax = axes[row_i, col_i]
            sr = info["search_results"]
            ax.plot(sr["k"], sr["silhouette"], 'o-', color='#2c3e50', linewidth=2, markersize=8)
            ax.axvline(x=sr["optimal_k"], color='red', linestyle='--', alpha=0.7,
                       label=f'Optimal k={sr["optimal_k"]}')
            best_i = sr["k"].index(sr["optimal_k"])
            ax.plot(sr["optimal_k"], sr["silhouette"][best_i], '*', color='red',
                    markersize=20, zorder=5)
            ax.set_xlabel("k")
            ax.set_ylabel("Silhouette Score")
            ax.set_title(combo_key.replace("__", "\n"), fontsize=10)
            ax.legend(fontsize=9)
            ax.grid(alpha=0.3)
            ax.set_xticks(sr["k"])

        # Hide unused axes
        for idx in range(n_plots, n_rows_sil * n_cols_sil):
            row_i, col_i = divmod(idx, n_cols_sil)
            axes[row_i, col_i].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_dir / "kmedoids_silhouette_scores.png", dpi=150, bbox_inches='tight')
        plt.show()
        

    # --- 3D scatter ---
    combo_groups = list(df_km.groupby(["protein", "ligand"]))
    n_combos = len(combo_groups)
    n_cols_3d = 4
    n_rows_3d = (n_combos + n_cols_3d - 1) // n_cols_3d

    fig3d = plt.figure(figsize=(8 * n_cols_3d, 7 * n_rows_3d))
    fig3d.suptitle("K-Medoids — 3D Binding Pocket Clusters\n"
                   "Color = Cluster, Shape = Docking Method",
                   fontsize=14, fontweight='bold', y=1.02)

    for idx, ((protein, ligand), group) in enumerate(combo_groups):
        row_i, col_i = divmod(idx, n_cols_3d)
        ax = fig3d.add_subplot(n_rows_3d, n_cols_3d, idx + 1, projection='3d')
        n_cl = group["cluster_kmedoids"].nunique()
        cmap = cm.get_cmap('tab10', max(n_cl, 1))

        for _, row in group.iterrows():
            ax.scatter(row["centroid_x"], row["centroid_y"], row["centroid_z"],
                       c=[cmap(int(row["cluster_kmedoids"]))],
                       marker=METHOD_MARKERS.get(row["method"], "o"),
                       s=80, alpha=0.8, edgecolors='black', linewidth=0.5)

        for c_id in range(n_cl):
            mask = group["cluster_kmedoids"] == c_id
            if mask.sum() > 0:
                center = group.loc[mask, ["centroid_x", "centroid_y", "centroid_z"]].mean()
                ax.scatter(*center, c='black', marker='X', s=200,
                           edgecolors='white', linewidth=2, zorder=10)
                ax.text(center["centroid_x"], center["centroid_y"],
                        center["centroid_z"], f'  C{c_id}', fontsize=10, fontweight='bold')

        ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
        ax.set_title(f"{protein}\n{ligand}\n({n_cl} clusters)", fontsize=11)
        ax.view_init(elev=20, azim=45)

    # Dynamic legend from METHODS_LIST
    legend_elements = [
        plt.Line2D([0], [0], marker=METHOD_MARKERS[m], color='w',
                   markerfacecolor='gray', markersize=10, label=m)
        for m in METHODS_LIST
    ] + [
        plt.Line2D([0], [0], marker='X', color='w', markerfacecolor='black',
                   markersize=12, label='Cluster Center'),
    ]
    fig3d.legend(handles=legend_elements, loc='lower center',
                 ncol=min(len(legend_elements), 6),
                 fontsize=10, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout()
    plt.savefig(output_dir / "kmedoids_3d_clusters.png", dpi=150, bbox_inches='tight')
    plt.show()

    # --- Heatmap ---
    df_heat = df_km_overlap.copy()
    df_heat["pocket_label"] = (df_heat["protein"].str[:20] + " / " +
                                df_heat["ligand"].str[:15] + " / C" +
                                df_heat["cluster_id"].astype(str))
    heat_data = df_heat.set_index("pocket_label")[METHODS_LIST]

    if len(heat_data) > 0:
        n_methods = len(METHODS_LIST)
        fig_h, ax_h = plt.subplots(figsize=(max(6, 2 * n_methods + 2),
                                            max(3, len(heat_data) * 0.5 + 2)))
        im = ax_h.imshow(heat_data.values, cmap='YlOrRd', aspect='auto')
        ax_h.set_xticks(range(n_methods))
        ax_h.set_xticklabels(METHODS_LIST, fontsize=10, rotation=30, ha='right')
        ax_h.set_yticks(range(len(heat_data)))
        ax_h.set_yticklabels(heat_data.index, fontsize=8)
        for i in range(len(heat_data)):
            for j in range(n_methods):
                val = int(heat_data.values[i, j])
                color = 'white' if val > heat_data.values.max() * 0.6 else 'black'
                ax_h.text(j, i, str(val), ha='center', va='center',
                          fontsize=10, fontweight='bold', color=color)
        plt.colorbar(im, ax=ax_h, shrink=0.8).set_label("Number of Poses")
        ax_h.set_title("K-Medoids: Poses per Method per Pocket Cluster",
                        fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_dir / "kmedoids_heatmap.png", dpi=150, bbox_inches='tight')
        plt.show()

    # --- Save ---
    df_km.to_csv(output_dir / "clustering_kmedoids_all_poses.csv", index=False)
    df_km_overlap.to_csv(output_dir / "clustering_kmedoids_pockets.csv", index=False)
    print(f"\nSaved: {output_dir / 'clustering_kmedoids_all_poses.csv'}")
    print(f"Saved: {output_dir / 'clustering_kmedoids_pockets.csv'}")
    print("=" * 100)

# %% [markdown]
# # Hierarchical Clustering

# %%
"""
CELL: Hierarchical Clustering (Ward Criterium) of Binding Pockets
===================================================================
Cluster all docked poses using Agglomerative Clustering with Ward linkage.
- Automatically selects optimal k via silhouette score
- Shows poses per docking method in each pocket cluster
- Identifies overlapping vs unique binding pockets
- Includes dendrograms
"""

# Short abbreviations for dendrogram labels (dynamic from METHODS_LIST)
_ABBREV = {m: m[:2].upper() for m in METHODS_LIST}

# ============================================================================
# WARD: FIND OPTIMAL K
# ============================================================================

def find_optimal_k_ward(coords, k_range):
    """Find optimal k for Ward clustering using silhouette score."""
    results = {"k": [], "silhouette": [], "inertia": []}
    for k in k_range:
        if k >= len(coords):
            break
        model = AgglomerativeClustering(n_clusters=k, linkage='ward')
        labels = model.fit_predict(coords)
        # Within-cluster sum of squared distances
        inertia = 0
        for c in range(k):
            mask = labels == c
            if mask.sum() > 0:
                center = coords[mask].mean(axis=0)
                inertia += np.sum((coords[mask] - center) ** 2)
        sil = silhouette_score(coords, labels) if len(set(labels)) > 1 else -1
        results["k"].append(k)
        results["silhouette"].append(sil)
        results["inertia"].append(inertia)
    best_idx = np.argmax(results["silhouette"])
    results["optimal_k"] = results["k"][best_idx]
    results["best_silhouette"] = results["silhouette"][best_idx]
    return results


# ============================================================================
# WARD: CLUSTER PER PROTEIN-LIGAND COMBINATION
# ============================================================================

if df_poses.empty:
    print("ERROR: No poses available. Run the data collection cell first.")
else:
    print("=" * 100)
    print("HIERARCHICAL CLUSTERING (WARD CRITERIUM) OF BINDING POCKETS")
    print("=" * 100)

    df_ward = df_poses.copy()
    df_ward["cluster_ward"] = -1
    ward_cluster_info = {}

    for (protein, ligand), group in df_ward.groupby(["protein", "ligand"]):
        combo_key = f"{protein}__{ligand}"
        idx = group.index
        coords = group[["centroid_x", "centroid_y", "centroid_z"]].values
        n_poses = len(coords)

        if n_poses < 3:
            df_ward.loc[idx, "cluster_ward"] = 0
            ward_cluster_info[combo_key] = {
                "n_poses": n_poses, "optimal_k": 1,
                "best_silhouette": -1, "search_results": None,
            }
            continue

        max_k = min(max(K_RANGE), n_poses - 1)
        actual_k_range = range(2, max_k + 1)

        if len(actual_k_range) == 0:
            df_ward.loc[idx, "cluster_ward"] = 0
            ward_cluster_info[combo_key] = {
                "n_poses": n_poses, "optimal_k": 1,
                "best_silhouette": -1, "search_results": None,
            }
            continue

        # Find optimal k
        ward_results = find_optimal_k_ward(coords, actual_k_range)
        k_opt = ward_results["optimal_k"]

        # Final fit with optimal k
        ward_model = AgglomerativeClustering(n_clusters=k_opt, linkage='ward')
        labels = ward_model.fit_predict(coords)
        df_ward.loc[idx, "cluster_ward"] = labels

        ward_cluster_info[combo_key] = {
            "n_poses": n_poses, "optimal_k": k_opt,
            "best_silhouette": ward_results["best_silhouette"],
            "search_results": ward_results,
        }

    df_ward["cluster_ward"] = df_ward["cluster_ward"].astype(int)

    # ========================================================================
    # RESULTS
    # ========================================================================

    # --- Optimal k per combination ---
    print(f"\n  {'Protein-Ligand':<55} {'Optimal k':>10} {'Silhouette':>12} {'Poses':>8}")
    print("  " + "-" * 90)
    for combo_key, info in sorted(ward_cluster_info.items()):
        print(f"  {combo_key:<55} {info['optimal_k']:>10} "
              f"{info['best_silhouette']:>12.3f} {info['n_poses']:>8}")

    # --- Cluster statistics table ---
    df_ward_stats = compute_cluster_stats(df_ward, "cluster_ward")
    df_ward_overlap = identify_pocket_overlap(df_ward_stats)

    print("\n" + "=" * 100)
    print("WARD — Binding Pocket Clusters (poses per method)")
    print("=" * 100)

    ward_display_cols = ["protein", "ligand", "cluster_id", "n_poses", "radius_A"] + \
                        METHODS_LIST + ["pocket_type", "n_methods_in_pocket"]
    ward_display = df_ward_overlap[ward_display_cols].sort_values(
        ["protein", "ligand", "cluster_id"])
    print(ward_display.to_string(index=False))

    # --- Overlapping vs unique ---
    print("\n" + "=" * 100)
    print("WARD — Overlapping vs Unique Binding Pockets")
    print("=" * 100)

    agg_dict = {"n_poses": ["count", "sum"]}
    for m in METHODS_LIST:
        agg_dict[m] = "sum"
    pocket_summary_ward = df_ward_overlap.groupby("pocket_type").agg(agg_dict)
    pocket_summary_ward.columns = (
        ["N_Clusters", "Total_Poses"] +
        [f"{m}_Poses" for m in METHODS_LIST]
    )
    pocket_summary_ward = pocket_summary_ward.reset_index()
    print(pocket_summary_ward.to_string(index=False))

    # --- Aggregated per protein-ligand ---
    print("\n" + "=" * 100)
    print("WARD — Aggregated per Protein-Ligand")
    print("=" * 100)

    def _agg_func_ward(g):
        result = {
            "n_clusters": len(g),
            "overlapping_pockets": (g["pocket_type"] == "Overlapping").sum(),
            "unique_pockets": (g["n_methods_in_pocket"] == 1).sum(),
            "total_poses": g["n_poses"].sum(),
        }
        for m in METHODS_LIST:
            result[f"{m}_in_overlapping"] = g.loc[
                g["pocket_type"] == "Overlapping", m
            ].sum()
        return pd.Series(result)

    agg_ward = df_ward_overlap.groupby(["protein", "ligand"]).apply(
        _agg_func_ward).reset_index()
    int_cols = [c for c in agg_ward.columns if c not in ["protein", "ligand"]]
    for c in int_cols:
        agg_ward[c] = agg_ward[c].astype(int)
    print(agg_ward.to_string(index=False))

    # --- Silhouette score plots ---
    combos_with_search = [(k, v) for k, v in ward_cluster_info.items()
                          if v["search_results"] is not None]

    if combos_with_search:
        n_plots = len(combos_with_search)
        n_cols_sil = min(4, n_plots)
        n_rows_sil = math.ceil(n_plots / n_cols_sil)
        fig, axes = plt.subplots(n_rows_sil, n_cols_sil,
                                 figsize=(6 * n_cols_sil, 4 * n_rows_sil), squeeze=False)
        fig.suptitle("Ward: Silhouette Score vs k (optimal k selection)",
                     fontsize=14, fontweight='bold')

        for idx, (combo_key, info) in enumerate(combos_with_search):
            row_i, col_i = divmod(idx, n_cols_sil)
            ax = axes[row_i, col_i]
            sr = info["search_results"]
            ax.plot(sr["k"], sr["silhouette"], 'o-', color='#8e44ad', linewidth=2, markersize=8)
            ax.axvline(x=sr["optimal_k"], color='red', linestyle='--', alpha=0.7,
                       label=f'Optimal k={sr["optimal_k"]}')
            best_i = sr["k"].index(sr["optimal_k"])
            ax.plot(sr["optimal_k"], sr["silhouette"][best_i], '*', color='red',
                    markersize=20, zorder=5)
            ax.set_xlabel("k")
            ax.set_ylabel("Silhouette Score")
            ax.set_title(combo_key.replace("__", "\n"), fontsize=10)
        # Hide unused axes
        for idx in range(n_plots, n_rows_sil * n_cols_sil):
            row_i, col_i = divmod(idx, n_cols_sil)
            axes[row_i, col_i].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_dir / "ward_silhouette_scores.png", dpi=150, bbox_inches='tight')

    # --- Dendrograms ---
    combo_groups = list(df_ward.groupby(["protein", "ligand"]))
    n_combos = len(combo_groups)
    n_cols = min(4, n_combos)
    n_rows = math.ceil(n_combos / n_cols)

    fig_dend, axes_d = plt.subplots(n_rows, n_cols,
                                    figsize=(8 * n_cols, 6 * n_rows), squeeze=False)
    fig_dend.suptitle("Hierarchical Clustering Dendrograms (Ward Criterium)",
                      fontsize=16, fontweight='bold')

    for idx, ((protein, ligand), group) in enumerate(combo_groups):
        row_i, col_i = divmod(idx, n_cols)
        ax = axes_d[row_i, col_i]
        coords = group[["centroid_x", "centroid_y", "centroid_z"]].values

        if len(coords) < 2:
            ax.text(0.5, 0.5, "Only 1 pose", transform=ax.transAxes,
                    ha='center', va='center')
            ax.set_title(f"{protein}\n{ligand}")
            continue

        Z = linkage(coords, method='ward')

        labels_dend = []
        for _, row in group.iterrows():
            abbrev = _ABBREV.get(row["method"], "??")
            labels_dend.append(f"{abbrev}-R{row['pose_rank']}")

        combo_key = f"{protein}__{ligand}"
        info = ward_cluster_info.get(combo_key, {})
        k_opt = info.get("optimal_k", 2)

        color_threshold = None
        if k_opt > 1 and len(Z) >= k_opt:
            color_threshold = (Z[-k_opt + 1, 2] + Z[-k_opt, 2]) / 2

        dendrogram(Z, labels=labels_dend, ax=ax, leaf_rotation=90, leaf_font_size=8,
                   color_threshold=color_threshold)
        ax.set_title(f"{protein}\n{ligand}\n(optimal k={k_opt})", fontsize=11)
        ax.set_ylabel("Distance (Ward)")
        if color_threshold is not None:
            ax.axhline(y=color_threshold, color='red', linestyle='--', alpha=0.5,
                       label=f'k={k_opt} cut')
            ax.legend(fontsize=8)

    # Hide unused axes
    for idx in range(n_combos, n_rows * n_cols):
        row_i, col_i = divmod(idx, n_cols)
        axes_d[row_i, col_i].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_dir / "ward_dendrograms.png", dpi=150, bbox_inches='tight')
    plt.show()
    n_cols_3d = min(4, n_combos)
    n_rows_3d = math.ceil(n_combos / n_cols_3d)
    fig3d = plt.figure(figsize=(8 * n_cols_3d, 7 * n_rows_3d))
    fig3d.suptitle("Ward — 3D Binding Pocket Clusters\n"
                   "Color = Cluster, Shape = Docking Method",
                   fontsize=14, fontweight='bold', y=1.02)

    for idx, ((protein, ligand), group) in enumerate(combo_groups):
        ax = fig3d.add_subplot(n_rows_3d, n_cols_3d, idx + 1, projection='3d')
        n_cl = group["cluster_ward"].nunique()
        cmap = cm.get_cmap('tab10', max(n_cl, 1))

        for _, row in group.iterrows():
            ax.scatter(row["centroid_x"], row["centroid_y"], row["centroid_z"],
                       c=[cmap(int(row["cluster_ward"]))],
                       marker=METHOD_MARKERS.get(row["method"], "o"),
                       s=80, alpha=0.8, edgecolors='black', linewidth=0.5)

        for c_id in range(n_cl):
            mask = group["cluster_ward"] == c_id
            if mask.sum() > 0:
                center = group.loc[mask, ["centroid_x", "centroid_y", "centroid_z"]].mean()
                ax.scatter(*center, c='black', marker='X', s=200,
                           edgecolors='white', linewidth=2, zorder=10)
                ax.text(center["centroid_x"], center["centroid_y"],
                        center["centroid_z"], f'  C{c_id}', fontsize=10, fontweight='bold')

        ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
        ax.set_title(f"{protein}\n{ligand}\n({n_cl} clusters)", fontsize=11)
        ax.view_init(elev=20, azim=45)

    # Dynamic legend
    legend_elements = [
        plt.Line2D([0], [0], marker=METHOD_MARKERS[m], color='w',
                   markerfacecolor='gray', markersize=10, label=m)
        for m in METHODS_LIST
    ] + [
        plt.Line2D([0], [0], marker='X', color='w', markerfacecolor='black',
                   markersize=12, label='Cluster Center'),
    ]
    fig3d.legend(handles=legend_elements, loc='lower center',
                 ncol=min(len(legend_elements), 6),
                 fontsize=10, bbox_to_anchor=(0.5, -0.05))
    plt.tight_layout()
    plt.savefig(output_dir / "ward_3d_clusters.png", dpi=150, bbox_inches='tight')
    plt.show()

    # --- Heatmap ---
    df_heat = df_ward_overlap.copy()
    df_heat["pocket_label"] = (df_heat["protein"].str[:20] + " / " +
                                df_heat["ligand"].str[:15] + " / C" +
                                df_heat["cluster_id"].astype(str))
    heat_data = df_heat.set_index("pocket_label")[METHODS_LIST]
    n_methods = len(METHODS_LIST)

    if len(heat_data) > 0:
        fig_h, ax_h = plt.subplots(figsize=(max(6, 2 * n_methods + 2),
                                            max(3, len(heat_data) * 0.5 + 2)))
        im = ax_h.imshow(heat_data.values, cmap='YlOrRd', aspect='auto')
        ax_h.set_xticks(range(n_methods))
        ax_h.set_xticklabels(METHODS_LIST, fontsize=10, rotation=30, ha='right')
        ax_h.set_yticks(range(len(heat_data)))
        ax_h.set_yticklabels(heat_data.index, fontsize=8)
        for i in range(len(heat_data)):
            for j in range(n_methods):
                val = int(heat_data.values[i, j])
                color = 'white' if val > heat_data.values.max() * 0.6 else 'black'
                ax_h.text(j, i, str(val), ha='center', va='center',
                          fontsize=10, fontweight='bold', color=color)
        plt.colorbar(im, ax=ax_h, shrink=0.8).set_label("Number of Poses")
        ax_h.set_title("Ward: Poses per Method per Pocket Cluster",
                        fontsize=13, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_dir / "ward_heatmap.png", dpi=150, bbox_inches='tight')
        plt.show()

    # --- Compare with K-Medoids if available ---
    if 'df_km' in dir() and not df_km.empty:
        print("\n" + "=" * 100)
        print("COMPARISON: K-Medoids vs Ward Clustering")
        print("=" * 100)

        # Merge cluster labels
        df_combined = df_km[["file_path", "cluster_kmedoids"]].merge(
            df_ward[["file_path", "cluster_ward"]], on="file_path"
        )
        df_combined = df_poses.merge(df_combined, on="file_path")

        comparison = []
        for (protein, ligand), group in df_combined.groupby(["protein", "ligand"]):
            lkm = group["cluster_kmedoids"].values
            lw = group["cluster_ward"].values
            if len(set(lkm)) < 2 or len(set(lw)) < 2:
                ari = 1.0 if np.array_equal(lkm, lw) else 0.0
                nmi = 1.0 if np.array_equal(lkm, lw) else 0.0
            else:
                ari = adjusted_rand_score(lkm, lw)
                nmi = normalized_mutual_info_score(lkm, lw)
            comparison.append({
                "protein": protein, "ligand": ligand, "n_poses": len(group),
                "k_kmedoids": len(set(lkm)), "k_ward": len(set(lw)),
                "ARI": round(ari, 3), "NMI": round(nmi, 3),
                "same_k": len(set(lkm)) == len(set(lw)),
            })

        df_comparison = pd.DataFrame(comparison)

        print("\n  ARI = Adjusted Rand Index (1.0 = identical, 0.0 = random)")
        print("  NMI = Normalized Mutual Information (1.0 = identical, 0.0 = independent)\n")
        print(df_comparison.to_string(index=False))
        print(f"\n  Mean ARI: {df_comparison['ARI'].mean():.3f}")
        print(f"  Mean NMI: {df_comparison['NMI'].mean():.3f}")
        n_same = df_comparison['same_k'].sum()
        print(f"  Same k:   {n_same}/{len(df_comparison)} "
              f"({100*n_same/len(df_comparison):.1f}%)")

        # ARI/NMI bar chart
        fig_c, ax_c = plt.subplots(figsize=(10, 5))
        x = np.arange(len(df_comparison))
        w = 0.35
        ax_c.bar(x - w/2, df_comparison["ARI"], w, label="ARI", color="#3498db", alpha=0.85)
        ax_c.bar(x + w/2, df_comparison["NMI"], w, label="NMI", color="#e67e22", alpha=0.85)
        ax_c.axhline(y=1.0, color='green', linestyle='--', alpha=0.3, label='Perfect')
        ax_c.axhline(y=0.0, color='red', linestyle='--', alpha=0.3, label='Random')
        combo_labels = [f"{r['protein'][:20]}\n{r['ligand'][:15]}"
                        for _, r in df_comparison.iterrows()]
        ax_c.set_xticks(x)
        ax_c.set_xticklabels(combo_labels, fontsize=8)
        ax_c.set_ylabel("Score")
        ax_c.set_title("K-Medoids vs Ward: Clustering Agreement (ARI & NMI)",
                        fontsize=14, fontweight='bold')
        ax_c.legend()
        ax_c.grid(axis='y', alpha=0.3)
        ax_c.set_ylim(-0.1, 1.15)
        plt.tight_layout()
        plt.savefig(output_dir / "clustering_agreement_ari_nmi.png", dpi=150, bbox_inches='tight')
        plt.show()
        plt.tight_layout()
        df_comparison.to_csv(output_dir / "clustering_comparison.csv", index=False)

    # --- Save ---
    df_ward.to_csv(output_dir / "clustering_ward_all_poses.csv", index=False)
    df_ward_overlap.to_csv(output_dir / "clustering_ward_pockets.csv", index=False)
    print(f"\nSaved: {output_dir / 'clustering_ward_all_poses.csv'}")
    print(f"Saved: {output_dir / 'clustering_ward_pockets.csv'}")
    print("=" * 100)

# %%
"""
CELL: Conformity vs Diversity Analysis
========================================
Select 3 protein-ligand combinations with highly CONFORMING clusters
and 3 with highly DIVERSE clusters across docking methods.

Visualize:
1. 3D scatter plots of poses and clusters
2. Heatmap of inter-pocket distances
3. RMSD comparison across docking methods (balanced sampling)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations, product
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from sklearn.metrics import silhouette_score
from scipy.spatial.distance import pdist, squareform, cdist
from scipy.cluster.hierarchy import linkage, fcluster
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

# ============================================================================
# CONFIGURATION
# ============================================================================
N_BALANCED_POSES = 5  # number of poses per method for balanced RMSD comparison
# METHOD_COLORS, METHOD_MARKERS, METHODS_LIST already defined from Cell 3

# Short abbreviation map for display
TOOL_ABBREV = {m: m[:2].upper() for m in METHODS_LIST}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_mol_from_sdf(sdf_path):
    """Load molecule from SDF, return mol with conformer or None."""
    try:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=False)
        if mol is not None and mol.GetNumConformers() > 0:
            return mol
    except Exception:
        pass
    return None


def get_heavy_atom_coords(sdf_path):
    """Get heavy-atom coordinates from SDF file."""
    try:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=False)
            if mol is not None:
                mol = Chem.RemoveHs(mol)
        if mol is None or mol.GetNumConformers() == 0:
            return None
        conf = mol.GetConformer()
        return np.array([
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])
    except Exception:
        return None


def compute_centroid_rmsd(coords1, coords2):
    """Compute RMSD between two sets of coordinates using centroid alignment.
    If atom counts differ, use centroid distance as proxy."""
    if coords1 is None or coords2 is None:
        return np.nan
    if len(coords1) == len(coords2):
        # Kabsch-like: center both, compute RMSD
        c1 = coords1 - coords1.mean(axis=0)
        c2 = coords2 - coords2.mean(axis=0)
        diff = c1 - c2
        return np.sqrt(np.mean(np.sum(diff ** 2, axis=1)))
    else:
        # Different atom counts: use centroid distance
        return np.linalg.norm(coords1.mean(axis=0) - coords2.mean(axis=0))


def compute_rdkit_rmsd(mol1, mol2):
    """Compute RMSD using RDKit's GetBestRMS (handles atom mapping)."""
    try:
        if mol1 is None or mol2 is None:
            return np.nan
        # Remove hydrogens for comparison
        m1 = Chem.RemoveHs(mol1)
        m2 = Chem.RemoveHs(mol2)
        if m1.GetNumAtoms() != m2.GetNumAtoms():
            # Fall back to centroid distance
            c1 = m1.GetConformer()
            c2 = m2.GetConformer()
            coords1 = np.array([[c1.GetAtomPosition(i).x, c1.GetAtomPosition(i).y, c1.GetAtomPosition(i).z]
                                 for i in range(m1.GetNumAtoms())])
            coords2 = np.array([[c2.GetAtomPosition(i).x, c2.GetAtomPosition(i).y, c2.GetAtomPosition(i).z]
                                 for i in range(m2.GetNumAtoms())])
            return np.linalg.norm(coords1.mean(axis=0) - coords2.mean(axis=0))
        rmsd = rdMolAlign.GetBestRMS(m1, m2)
        return rmsd
    except Exception:
        return np.nan


def score_conformity(group_df, cluster_col):
    """
    Score how conforming a protein-ligand combination is across methods.
    High conformity = all methods place poses in the same cluster(s),
    low spread, high overlap.
    
    Returns a dict with conformity metrics.
    """
    n_clusters = group_df[cluster_col].nunique()
    n_methods = group_df["method"].nunique()
    n_poses = len(group_df)
    
    # Count how many clusters contain poses from ≥ 2 methods
    overlapping_clusters = 0
    total_overlap_poses = 0
    for c_id, c_group in group_df.groupby(cluster_col):
        methods_in_cluster = c_group["method"].nunique()
        if methods_in_cluster >= 2:
            overlapping_clusters += 1
            total_overlap_poses += len(c_group)
    
    # Fraction of poses in overlapping clusters
    overlap_fraction = total_overlap_poses / n_poses if n_poses > 0 else 0
    
    # Compute mean pairwise centroid distance between methods
    method_centroids = {}
    for method, m_group in group_df.groupby("method"):
        coords = m_group[["centroid_x", "centroid_y", "centroid_z"]].values
        method_centroids[method] = coords.mean(axis=0)
    
    inter_method_dist = []
    for (m1, c1), (m2, c2) in combinations(method_centroids.items(), 2):
        inter_method_dist.append(np.linalg.norm(c1 - c2))
    mean_inter_method_dist = np.mean(inter_method_dist) if inter_method_dist else 0
    
    # Compute spatial spread (mean distance of all poses from global centroid)
    all_coords = group_df[["centroid_x", "centroid_y", "centroid_z"]].values
    global_centroid = all_coords.mean(axis=0)
    spread = np.mean(np.linalg.norm(all_coords - global_centroid, axis=1))
    
    # Conformity score: high overlap + low inter-method distance + low spread + few clusters
    # Normalize to 0-1 range approximately
    conformity_score = (
        overlap_fraction * 0.4 +
        max(0, 1 - mean_inter_method_dist / 50) * 0.3 +
        max(0, 1 - spread / 30) * 0.2 +
        max(0, 1 - n_clusters / 10) * 0.1
    )
    
    return {
        "n_clusters": n_clusters,
        "n_methods": n_methods,
        "n_poses": n_poses,
        "overlapping_clusters": overlapping_clusters,
        "overlap_fraction": round(overlap_fraction, 3),
        "mean_inter_method_dist": round(mean_inter_method_dist, 2),
        "spatial_spread": round(spread, 2),
        "conformity_score": round(conformity_score, 4),
    }


# ============================================================================
# STEP 1: SCORE ALL PROTEIN-LIGAND COMBINATIONS
# ============================================================================

# Check which data is available (prefer Ward, fall back to K-Medoids)
if 'df_ward' in dir() and not df_ward.empty:
    df_analysis = df_ward.copy()
    cluster_col = "cluster_ward"
    clustering_label = "Ward"
    print("Using Ward clustering results for selection.")
elif 'df_km' in dir() and not df_km.empty:
    df_analysis = df_km.copy()
    cluster_col = "cluster_kmedoids"
    clustering_label = "K-Medoids"
    print("Using K-Medoids clustering results for selection.")
else:
    raise RuntimeError("No clustering results found. Run K-Medoids or Ward cell first.")

# Only consider combinations with poses from at least 2 methods
combo_scores = []
for (protein, ligand), group in df_analysis.groupby(["protein", "ligand"]):
    n_methods = group["method"].nunique()
    if n_methods < 2:
        continue
    
    scores = score_conformity(group, cluster_col)
    scores["protein"] = protein
    scores["ligand"] = ligand
    combo_scores.append(scores)

df_scores = pd.DataFrame(combo_scores).sort_values("conformity_score", ascending=False)

print("\n" + "=" * 120)
print("CONFORMITY SCORES FOR ALL PROTEIN-LIGAND COMBINATIONS")
print("(Higher score = more conforming clusters across methods)")
print("=" * 120)
print(df_scores[[
    "protein", "ligand", "n_poses", "n_methods", "n_clusters",
    "overlapping_clusters", "overlap_fraction",
    "mean_inter_method_dist", "spatial_spread", "conformity_score"
]].to_string(index=False))

# ============================================================================
# STEP 2: SELECT TOP 3 CONFORMING AND TOP 3 DIVERSE
# ============================================================================

n_available = len(df_scores)
n_conform = min(3, n_available // 2) if n_available >= 2 else min(1, n_available)
n_diverse = min(3, n_available - n_conform)

top_conform = df_scores.head(n_conform)
top_diverse = df_scores.tail(n_diverse)

# Make sure no overlap
diverse_keys = set(zip(top_diverse["protein"], top_diverse["ligand"]))
conform_keys = set(zip(top_conform["protein"], top_conform["ligand"]))
if diverse_keys & conform_keys:
    # Remove overlapping from diverse, take next available
    remaining = df_scores[~df_scores.apply(
        lambda r: (r["protein"], r["ligand"]) in conform_keys, axis=1
    )]
    top_diverse = remaining.tail(n_diverse)

print("\n" + "=" * 120)
print(f"SELECTED — TOP {n_conform} HIGHLY CONFORMING combinations:")
print("=" * 120)
for _, row in top_conform.iterrows():
    print(f"  {row['protein']} + {row['ligand']:>30}  "
          f"(score={row['conformity_score']:.4f}, "
          f"overlap={row['overlap_fraction']:.1%}, "
          f"inter-method dist={row['mean_inter_method_dist']:.1f} Å, "
          f"clusters={row['n_clusters']})")

print(f"\nSELECTED — TOP {n_diverse} HIGHLY DIVERSE combinations:")
print("=" * 120)
for _, row in top_diverse.iterrows():
    print(f"  {row['protein']} + {row['ligand']:>30}  "
          f"(score={row['conformity_score']:.4f}, "
          f"overlap={row['overlap_fraction']:.1%}, "
          f"inter-method dist={row['mean_inter_method_dist']:.1f} Å, "
          f"clusters={row['n_clusters']})")

selected_conform = list(zip(top_conform["protein"], top_conform["ligand"]))
selected_diverse = list(zip(top_diverse["protein"], top_diverse["ligand"]))
all_selected = selected_conform + selected_diverse
category_map = {}
for p, l in selected_conform:
    category_map[(p, l)] = "CONFORMING"
for p, l in selected_diverse:
    category_map[(p, l)] = "DIVERSE"


# ============================================================================
# STEP 3: 3D SCATTER PLOTS
# ============================================================================

print("\n" + "=" * 120)
print("GENERATING 3D SCATTER PLOTS...")
print("=" * 120)

n_selected = len(all_selected)
fig_3d = plt.figure(figsize=(8 * min(n_selected, 3), 7 * ((n_selected + 2) // 3)))
fig_3d.suptitle(
    f"3D Binding Pocket Clusters ({clustering_label})\n"
    "Top row: CONFORMING | Bottom row: DIVERSE\n"
    "Color = Cluster, Shape = Docking Method",
    fontsize=16, fontweight='bold', y=1.03
)

for plot_idx, (protein, ligand) in enumerate(all_selected):
    group = df_analysis[
        (df_analysis["protein"] == protein) & (df_analysis["ligand"] == ligand)
    ]
    category = category_map[(protein, ligand)]
    
    ax = fig_3d.add_subplot(
        (n_selected + 2) // 3, min(n_selected, 3), plot_idx + 1, projection='3d'
    )
    
    n_cl = group[cluster_col].nunique()
    cmap_clusters = cm.get_cmap('tab10', max(n_cl, 1))
    
    # Plot each pose
    for _, row in group.iterrows():
        c_id = int(row[cluster_col])
        method = row["method"]
        ax.scatter(
            row["centroid_x"], row["centroid_y"], row["centroid_z"],
            c=[cmap_clusters(c_id)],
            marker=METHOD_MARKERS.get(method, "o"),
            s=100, alpha=0.85, edgecolors='black', linewidth=0.5
        )
    
    # Plot cluster centers
    for c_id in group[cluster_col].unique():
        mask = group[cluster_col] == c_id
        if mask.sum() > 0:
            center = group.loc[mask, ["centroid_x", "centroid_y", "centroid_z"]].mean()
            ax.scatter(
                center["centroid_x"], center["centroid_y"], center["centroid_z"],
                c='black', marker='X', s=250, edgecolors='white', linewidth=2, zorder=10
            )
            ax.text(
                center["centroid_x"], center["centroid_y"], center["centroid_z"],
                f'  C{c_id}', fontsize=9, fontweight='bold', color='black'
            )
    
    # Title with category badge
    badge_color = '#27ae60' if category == "CONFORMING" else '#c0392b'
    ax.set_title(
        f"[{category}]\n{protein}\n{ligand}\n({n_cl} clusters)",
        fontsize=10, fontweight='bold', color=badge_color
    )
    ax.set_xlabel("X (Å)", fontsize=8)
    ax.set_ylabel("Y (Å)", fontsize=8)
    ax.set_zlabel("Z (Å)", fontsize=8)
    ax.view_init(elev=25, azim=45)

# Dynamic legend
legend_elements = [
    plt.Line2D([0], [0], marker=METHOD_MARKERS[m], color='w',
               markerfacecolor='gray', markersize=10, label=m)
    for m in METHODS_LIST
] + [
    plt.Line2D([0], [0], marker='X', color='w', markerfacecolor='black',
               markersize=12, label='Cluster Center'),
]
fig_3d.legend(handles=legend_elements, loc='lower center',
              ncol=min(len(legend_elements), 6),
              fontsize=10, bbox_to_anchor=(0.5, -0.04))
plt.tight_layout()
plt.savefig(output_dir / "selected_3d_clusters.png", dpi=150, bbox_inches='tight')
plt.show()
print(f"Saved: {output_dir / 'selected_3d_clusters.png'}")


# ============================================================================
# STEP 4: INTER-POCKET DISTANCE HEATMAPS
# ============================================================================

print("\n" + "=" * 120)
print("GENERATING INTER-POCKET DISTANCE HEATMAPS...")
print("=" * 120)

n_rows = (n_selected + 2) // 3
n_cols = min(n_selected, 3)
fig_dist, axes_dist = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 6 * n_rows),
                                    squeeze=False)
fig_dist.suptitle("Inter-Pocket Centroid Distance (Å)", fontsize=16, fontweight='bold')

for plot_idx, (protein, ligand) in enumerate(all_selected):
    row_i, col_i = divmod(plot_idx, n_cols)
    ax = axes_dist[row_i][col_i]
    category = category_map[(protein, ligand)]
    
    group = df_analysis[
        (df_analysis["protein"] == protein) & (df_analysis["ligand"] == ligand)
    ]
    
    # Compute cluster centers
    cluster_centers = {}
    for c_id, c_group in group.groupby(cluster_col):
        coords = c_group[["centroid_x", "centroid_y", "centroid_z"]].values
        center = coords.mean(axis=0)
        methods_present = sorted(c_group["method"].unique())
        methods_str = "/".join([TOOL_ABBREV.get(m, "?") for m in methods_present])
        cluster_centers[c_id] = {
            "center": center,
            "label": f"C{c_id}\n({methods_str})\nn={len(c_group)}",
            "n_poses": len(c_group),
        }
    
    cluster_ids = sorted(cluster_centers.keys())
    n_cl = len(cluster_ids)
    
    if n_cl < 2:
        ax.text(0.5, 0.5, f"Only {n_cl} cluster", transform=ax.transAxes,
                ha='center', va='center', fontsize=14)
        badge_color = '#27ae60' if category == "CONFORMING" else '#c0392b'
        ax.set_title(f"[{category}] {protein}\n{ligand}", fontsize=10,
                     fontweight='bold', color=badge_color)
        continue
    
    # Distance matrix between cluster centers
    centers_array = np.array([cluster_centers[c]["center"] for c in cluster_ids])
    dist_matrix = squareform(pdist(centers_array, metric='euclidean'))
    
    # Plot heatmap
    im = ax.imshow(dist_matrix, cmap='YlOrRd', aspect='auto')
    
    labels = [cluster_centers[c]["label"] for c in cluster_ids]
    ax.set_xticks(range(n_cl))
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha='right')
    ax.set_yticks(range(n_cl))
    ax.set_yticklabels(labels, fontsize=7)
    
    # Annotate distances
    for i in range(n_cl):
        for j in range(n_cl):
            color = 'white' if dist_matrix[i, j] > dist_matrix.max() * 0.6 else 'black'
            ax.text(j, i, f"{dist_matrix[i, j]:.1f}", ha='center', va='center',
                    fontsize=9, fontweight='bold', color=color)
    
    plt.colorbar(im, ax=ax, shrink=0.7, label="Distance (Å)")
    badge_color = '#27ae60' if category == "CONFORMING" else '#c0392b'
    ax.set_title(f"[{category}] {protein}\n{ligand}", fontsize=10,
                 fontweight='bold', color=badge_color)

# Hide unused axes
for idx in range(n_selected, n_rows * n_cols):
    row_i, col_i = divmod(idx, n_cols)
    axes_dist[row_i][col_i].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / "selected_pocket_distances.png", dpi=150, bbox_inches='tight')
plt.show()
print(f"Saved: {output_dir / 'selected_pocket_distances.png'}")


# ============================================================================
# STEP 5: BALANCED RMSD COMPARISON ACROSS DOCKING METHODS
# ============================================================================

print("\n" + "=" * 120)
print("COMPUTING BALANCED RMSD COMPARISONS...")
print(f"(Using up to {N_BALANCED_POSES} poses per method per pocket)")
print("=" * 120)

all_rmsd_records = []

for (protein, ligand) in all_selected:
    category = category_map[(protein, ligand)]
    group = df_analysis[
        (df_analysis["protein"] == protein) & (df_analysis["ligand"] == ligand)
    ]
    
    # For each cluster (pocket), sample balanced poses
    for c_id, c_group in group.groupby(cluster_col):
        method_poses = {}
        for method in METHODS_LIST:
            m_group = c_group[c_group["method"] == method]
            if len(m_group) == 0:
                continue
            # Sample up to N_BALANCED_POSES (prefer top-ranked)
            m_sorted = m_group.sort_values("pose_rank")
            sampled = m_sorted.head(N_BALANCED_POSES)
            
            # Load coordinates for each pose
            pose_coords = []
            pose_files = []
            for _, row in sampled.iterrows():
                coords = get_heavy_atom_coords(row["file_path"])
                if coords is not None:
                    pose_coords.append(coords)
                    pose_files.append(row["file_path"])
            
            if pose_coords:
                method_poses[method] = {
                    "coords": pose_coords,
                    "files": pose_files,
                    "n_used": len(pose_coords),
                }
        
        # Compute pairwise RMSD between all methods in this pocket
        available_methods = sorted(method_poses.keys())
        
        if len(available_methods) < 2:
            continue
        
        # Balance: use the minimum number of poses across available methods
        min_n = min(method_poses[m]["n_used"] for m in available_methods)
        
        for m1, m2 in combinations(available_methods, 2):
            coords1_list = method_poses[m1]["coords"][:min_n]
            coords2_list = method_poses[m2]["coords"][:min_n]
            
            # Compute all pairwise RMSDs between poses from m1 and m2
            pairwise_rmsds = []
            for c1 in coords1_list:
                for c2 in coords2_list:
                    r = compute_centroid_rmsd(c1, c2)
                    if not np.isnan(r):
                        pairwise_rmsds.append(r)
            
            if pairwise_rmsds:
                all_rmsd_records.append({
                    "protein": protein,
                    "ligand": ligand,
                    "category": category,
                    "cluster_id": c_id,
                    "method_1": m1,
                    "method_2": m2,
                    "n_poses_m1": len(coords1_list),
                    "n_poses_m2": len(coords2_list),
                    "mean_rmsd": round(np.mean(pairwise_rmsds), 2),
                    "median_rmsd": round(np.median(pairwise_rmsds), 2),
                    "min_rmsd": round(np.min(pairwise_rmsds), 2),
                    "max_rmsd": round(np.max(pairwise_rmsds), 2),
                    "std_rmsd": round(np.std(pairwise_rmsds), 2),
                })
        
        # Also compute intra-method RMSD
        for method in available_methods:
            coords_list = method_poses[method]["coords"][:min_n]
            if len(coords_list) < 2:
                continue
            intra_rmsds = []
            for c1, c2 in combinations(coords_list, 2):
                r = compute_centroid_rmsd(c1, c2)
                if not np.isnan(r):
                    intra_rmsds.append(r)
            if intra_rmsds:
                all_rmsd_records.append({
                    "protein": protein,
                    "ligand": ligand,
                    "category": category,
                    "cluster_id": c_id,
                    "method_1": method,
                    "method_2": method,
                    "n_poses_m1": len(coords_list),
                    "n_poses_m2": len(coords_list),
                    "mean_rmsd": round(np.mean(intra_rmsds), 2),
                    "median_rmsd": round(np.median(intra_rmsds), 2),
                    "min_rmsd": round(np.min(intra_rmsds), 2),
                    "max_rmsd": round(np.max(intra_rmsds), 2),
                    "std_rmsd": round(np.std(intra_rmsds), 2),
                })

df_rmsd = pd.DataFrame(all_rmsd_records)

if not df_rmsd.empty:
    print("\n" + "=" * 120)
    print("RMSD COMPARISON TABLE (balanced poses per method)")
    print("=" * 120)
    print(df_rmsd.to_string(index=False))
    
    # ========================================================================
    # STEP 6: RMSD HEATMAPS PER SELECTED COMBINATION
    # ========================================================================
    
    print("\n" + "=" * 120)
    print("GENERATING RMSD HEATMAPS...")
    print("=" * 120)
    
    fig_rmsd, axes_rmsd = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 6 * n_rows),
                                        squeeze=False)
    fig_rmsd.suptitle(
        "Mean RMSD Between Docking Methods per Pocket Cluster (Å)\n"
        f"(Balanced: up to {N_BALANCED_POSES} poses per method)",
        fontsize=16, fontweight='bold'
    )
    
    n_m = len(METHODS_LIST)
    
    for plot_idx, (protein, ligand) in enumerate(all_selected):
        row_i, col_i = divmod(plot_idx, n_cols)
        ax = axes_rmsd[row_i][col_i]
        category = category_map[(protein, ligand)]
        
        df_sub = df_rmsd[
            (df_rmsd["protein"] == protein) & (df_rmsd["ligand"] == ligand)
        ]
        
        if df_sub.empty:
            ax.text(0.5, 0.5, "No RMSD data", transform=ax.transAxes,
                    ha='center', va='center', fontsize=14)
            badge_color = '#27ae60' if category == "CONFORMING" else '#c0392b'
            ax.set_title(f"[{category}] {protein}\n{ligand}", fontsize=10,
                         fontweight='bold', color=badge_color)
            continue
        
        # Get unique clusters
        clusters = sorted(df_sub["cluster_id"].unique())
        
        # Overall aggregated RMSD matrix
        rmsd_matrix = np.full((n_m, n_m), np.nan)
        
        for i, m1 in enumerate(METHODS_LIST):
            for j, m2 in enumerate(METHODS_LIST):
                # Get all RMSD values for this method pair across all clusters
                mask = (
                    ((df_sub["method_1"] == m1) & (df_sub["method_2"] == m2)) |
                    ((df_sub["method_1"] == m2) & (df_sub["method_2"] == m1))
                )
                sub = df_sub[mask]
                if not sub.empty:
                    rmsd_matrix[i, j] = sub["mean_rmsd"].mean()
        
        # Plot
        im = ax.imshow(rmsd_matrix, cmap='RdYlGn_r', aspect='auto',
                       vmin=0, vmax=max(np.nanmax(rmsd_matrix), 5))
        
        method_labels = [TOOL_ABBREV.get(m, m[:6]) for m in METHODS_LIST]
        ax.set_xticks(range(n_m))
        ax.set_xticklabels(method_labels, fontsize=9, rotation=45)
        ax.set_yticks(range(n_m))
        ax.set_yticklabels(method_labels, fontsize=9)
        
        # Annotate
        for i in range(n_m):
            for j in range(n_m):
                val = rmsd_matrix[i, j]
                if np.isnan(val):
                    ax.text(j, i, "N/A", ha='center', va='center',
                            fontsize=10, fontweight='bold', color='gray')
                else:
                    color = 'white' if val > np.nanmax(rmsd_matrix) * 0.6 else 'black'
                    label_str = f"{val:.1f}" if i != j else f"{val:.1f}\n(intra)"
                    ax.text(j, i, label_str, ha='center', va='center',
                            fontsize=10, fontweight='bold', color=color)
        
        plt.colorbar(im, ax=ax, shrink=0.7, label="Mean RMSD (Å)")
        badge_color = '#27ae60' if category == "CONFORMING" else '#c0392b'
        n_cl = len(clusters)
        ax.set_title(
            f"[{category}] {protein}\n{ligand}\n({n_cl} pocket{'s' if n_cl > 1 else ''})",
            fontsize=10, fontweight='bold', color=badge_color
        )
    
    # Hide unused axes
    for idx in range(n_selected, n_rows * n_cols):
        row_i, col_i = divmod(idx, n_cols)
        axes_rmsd[row_i][col_i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_dir / "selected_rmsd_heatmaps.png", dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Saved: {output_dir / 'selected_rmsd_heatmaps.png'}")
    
    # ========================================================================
    # STEP 7: RMSD BAR CHART COMPARISON (Conforming vs Diverse)
    # ========================================================================
    
    # Aggregate RMSD by category and method pair
    df_rmsd_inter = df_rmsd[df_rmsd["method_1"] != df_rmsd["method_2"]].copy()
    df_rmsd_inter["method_pair"] = df_rmsd_inter.apply(
        lambda r: " vs ".join(sorted([
            TOOL_ABBREV.get(r["method_1"], r["method_1"][:6]),
            TOOL_ABBREV.get(r["method_2"], r["method_2"][:6])
        ])),
        axis=1
    )
    
    agg = df_rmsd_inter.groupby(["category", "method_pair"]).agg({
        "mean_rmsd": ["mean", "std"]
    }).reset_index()
    agg.columns = ["category", "method_pair", "avg_rmsd", "std_rmsd"]
    
    if not agg.empty:
        fig_bar, ax_bar = plt.subplots(figsize=(12, 6))
        
        pairs = sorted(agg["method_pair"].unique())
        x = np.arange(len(pairs))
        width = 0.35
        
        for i, cat in enumerate(["CONFORMING", "DIVERSE"]):
            cat_data = agg[agg["category"] == cat]
            vals = []
            errs = []
            for p in pairs:
                sub = cat_data[cat_data["method_pair"] == p]
                if not sub.empty:
                    vals.append(sub["avg_rmsd"].values[0])
                    errs.append(sub["std_rmsd"].values[0] if not np.isnan(sub["std_rmsd"].values[0]) else 0)
                else:
                    vals.append(0)
                    errs.append(0)
            
            color = '#27ae60' if cat == "CONFORMING" else '#c0392b'
            bars = ax_bar.bar(x + (i - 0.5) * width, vals, width, yerr=errs,
                             label=cat, color=color, alpha=0.8, capsize=5)
            
            # Add value labels
            for bar, val in zip(bars, vals):
                if val > 0:
                    ax_bar.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.2,
                               f'{val:.1f}', ha='center', va='bottom', fontsize=9,
                               fontweight='bold')
        
        ax_bar.set_xticks(x)
        ax_bar.set_xticklabels(pairs, fontsize=11)
        ax_bar.set_ylabel("Mean RMSD (Å)", fontsize=12)
        ax_bar.set_title(
            "Inter-Method RMSD: Conforming vs Diverse Protein-Ligand Combinations",
            fontsize=14, fontweight='bold'
        )
        ax_bar.legend(fontsize=11)
        ax_bar.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(output_dir / "selected_rmsd_conform_vs_diverse.png", dpi=150,
                    bbox_inches='tight')
        plt.show()
        print(f"Saved: {output_dir / 'selected_rmsd_conform_vs_diverse.png'}")
    
    # Save all RMSD data
    df_rmsd.to_csv(output_dir / "selected_rmsd_comparison.csv", index=False)
    df_scores.to_csv(output_dir / "conformity_scores_all.csv", index=False)
    print(f"\nSaved: {output_dir / 'selected_rmsd_comparison.csv'}")
    print(f"Saved: {output_dir / 'conformity_scores_all.csv'}")

else:
    print("\nWARNING: No RMSD data could be computed. Check that pose SDF files are accessible.")


# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 120)
print("SUMMARY")
print("=" * 120)
print(f"\n  Clustering method used: {clustering_label}")
print(f"  Total protein-ligand combinations scored: {len(df_scores)}")
print(f"  Selected conforming:  {len(selected_conform)}")
print(f"  Selected diverse:     {len(selected_diverse)}")
if not df_rmsd.empty:
    conform_rmsd = df_rmsd[
        (df_rmsd["category"] == "CONFORMING") & (df_rmsd["method_1"] != df_rmsd["method_2"])
    ]["mean_rmsd"]
    diverse_rmsd = df_rmsd[
        (df_rmsd["category"] == "DIVERSE") & (df_rmsd["method_1"] != df_rmsd["method_2"])
    ]["mean_rmsd"]
    print(f"\n  CONFORMING — inter-method RMSD: "
          f"mean={conform_rmsd.mean():.2f} Å, "
          f"median={conform_rmsd.median():.2f} Å" if not conform_rmsd.empty else "")
    print(f"  DIVERSE    — inter-method RMSD: "
          f"mean={diverse_rmsd.mean():.2f} Å, "
          f"median={diverse_rmsd.median():.2f} Å" if not diverse_rmsd.empty else "")
print("\n" + "=" * 120)

# %% [markdown]
# # Cluster Conformity

# %%
"""
CELL: Cluster Composition & Closest-Pose Inter-Method RMSD
============================================================
For EACH protein-ligand combination (hierarchical / Ward clustering):
  1. Number of clusters found
  2. Per cluster: how many poses each docking method contributed
  3. Balanced closest-pose RMSD between different docking methods
     → trim each method to the same N poses (= min across present
       methods, top-ranked first), then report the MINIMUM pairwise
       RMSD between every pair of methods (= closest pose agreement)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ============================================================================
# CONFIGURATION  (derived from docking_directories)
# ============================================================================
TOOL_ABBREV   = {m: m[:2].upper() for m in METHODS_LIST}

# ============================================================================
# HELPERS
# ============================================================================

def _load_heavy_coords(sdf_path):
    """Return heavy-atom coordinate array (N x 3) or None."""
    try:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=False)
            if mol is not None:
                mol = Chem.RemoveHs(mol)
        if mol is None or mol.GetNumConformers() == 0:
            return None
        conf = mol.GetConformer()
        return np.array([
            [conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])
    except Exception:
        return None


def _rmsd(c1, c2):
    """RMSD between two coord arrays; centroid distance if atom counts differ."""
    if c1 is None or c2 is None:
        return np.nan
    if len(c1) == len(c2):
        d = c1 - c2
        return float(np.sqrt(np.mean(np.sum(d ** 2, axis=1))))
    return float(np.linalg.norm(c1.mean(axis=0) - c2.mean(axis=0)))

# ============================================================================
# SELECT DATA SOURCE  (prefer Ward, fall back to K-Medoids)
# ============================================================================

if 'df_ward' in dir() and not df_ward.empty:
    df_src = df_ward.copy()
    clust_col  = "cluster_ward"
    clust_name = "Ward"
elif 'df_km' in dir() and not df_km.empty:
    df_src = df_km.copy()
    clust_col  = "cluster_kmedoids"
    clust_name = "K-Medoids"
else:
    raise RuntimeError("No clustering results found — run clustering cells first.")

print(f"Using {clust_name} clustering  ({clust_col})")
print(f"Total poses in dataset: {len(df_src)}")

combo_list = sorted(df_src.groupby(["protein", "ligand"]).groups.keys())
print(f"Protein-ligand combinations: {len(combo_list)}")

# ============================================================================
# STEP 1 — Cluster overview & pose composition table
# ============================================================================

print("\n" + "=" * 120)
print("CLUSTER OVERVIEW PER PROTEIN-LIGAND COMBINATION")
print("=" * 120)

overview_rows = []

for protein, ligand in combo_list:
    grp = df_src[(df_src["protein"] == protein) & (df_src["ligand"] == ligand)]
    cluster_ids = sorted(grp[clust_col].unique())
    n_clusters = len(cluster_ids)

    for cid in cluster_ids:
        c_grp = grp[grp[clust_col] == cid]
        method_counts = {m: int((c_grp["method"] == m).sum()) for m in METHODS_LIST}
        overview_rows.append({
            "protein": protein,
            "ligand": ligand,
            "n_clusters": n_clusters,
            "cluster": cid,
            "total_poses": len(c_grp),
            **{TOOL_ABBREV[m]: method_counts[m] for m in METHODS_LIST},
            "methods_present": ", ".join(
                TOOL_ABBREV[m] for m in METHODS_LIST if method_counts[m] > 0
            ),
        })

df_overview = pd.DataFrame(overview_rows)
print(df_overview.to_string(index=False))

# ============================================================================
# STEP 2 — Balanced closest-pose RMSD between docking methods
# ============================================================================
# For every cluster that contains poses from ≥ 2 methods:
#   • balance to the same N poses per method (min count, top-rank first)
#   • compute ALL pairwise RMSD between methods
#   • report the CLOSEST (minimum) RMSD per method-pair

print("\n" + "=" * 120)
print("CLOSEST-POSE RMSD BETWEEN DOCKING METHODS  (balanced, per cluster)")
print("=" * 120)

rmsd_records  = []   # summary: one row per method-pair per cluster
detail_records = []  # every pairwise RMSD value

for protein, ligand in combo_list:
    grp = df_src[(df_src["protein"] == protein) & (df_src["ligand"] == ligand)]
    cluster_ids = sorted(grp[clust_col].unique())

    for cid in cluster_ids:
        c_grp = grp[grp[clust_col] == cid]

        # Load coordinates per method (top-ranked first)
        method_data = {}
        for meth in METHODS_LIST:
            m_sub = c_grp[c_grp["method"] == meth].sort_values("pose_rank")
            coords = []
            for _, rw in m_sub.iterrows():
                c = _load_heavy_coords(rw["file_path"])
                if c is not None:
                    coords.append(c)
            if coords:
                method_data[meth] = coords

        active_methods = sorted(method_data.keys())
        if len(active_methods) < 2:
            continue  # need ≥ 2 methods for inter-method RMSD

        # Balance: use the minimum pose count across present methods
        min_n = min(len(method_data[m]) for m in active_methods)
        if min_n == 0:
            continue
        for m in active_methods:
            method_data[m] = method_data[m][:min_n]

        # Compute inter-method pairwise RMSD → track closest
        for t1, t2 in combinations(active_methods, 2):
            pair_label = f"{TOOL_ABBREV[t1]} vs {TOOL_ABBREV[t2]}"
            all_rmsds = []
            for i, c1 in enumerate(method_data[t1]):
                for j, c2 in enumerate(method_data[t2]):
                    r = _rmsd(c1, c2)
                    if not np.isnan(r):
                        all_rmsds.append(r)
                        detail_records.append({
                            "protein": protein, "ligand": ligand,
                            "cluster": cid, "pair": pair_label,
                            "pose_i": f"{TOOL_ABBREV[t1]}-{i+1}",
                            "pose_j": f"{TOOL_ABBREV[t2]}-{j+1}",
                            "rmsd": round(r, 3),
                        })
            if all_rmsds:
                rmsd_records.append({
                    "protein": protein, "ligand": ligand,
                    "cluster": cid,
                    "pair": pair_label,
                    "n_poses_per_method": min_n,
                    "n_comparisons": len(all_rmsds),
                    "closest_rmsd": round(np.min(all_rmsds), 3),
                    "mean_rmsd": round(np.mean(all_rmsds), 3),
                    "median_rmsd": round(np.median(all_rmsds), 3),
                    "max_rmsd": round(np.max(all_rmsds), 3),
                    "std_rmsd": round(np.std(all_rmsds), 3),
                })

df_rmsd = pd.DataFrame(rmsd_records)
df_detail = pd.DataFrame(detail_records)

if df_rmsd.empty:
    print("\nNo clusters with ≥ 2 docking methods found — cannot compute inter-method RMSD.")
else:
    print(f"\nTotal method-pair RMSD rows: {len(df_rmsd)}")
    print(f"Total pairwise comparisons : {len(df_detail)}\n")

    display_cols = [
        "protein", "ligand", "cluster", "pair",
        "n_poses_per_method", "closest_rmsd",
        "mean_rmsd", "median_rmsd", "max_rmsd", "std_rmsd",
    ]
    print(df_rmsd[display_cols].to_string(index=False))

    # ====================================================================
    # STEP 3 — Per-combo figure: one subplot per cluster
    #          Box-plot of all pairwise RMSD + highlight closest
    # ====================================================================

    print("\n" + "=" * 120)
    print("GENERATING PER-COMBO CLUSTER RMSD PLOTS …")
    print("=" * 120)

    pair_color_map = {}
    base_colors = ["#e74c3c", "#3498db", "#2ecc71", "#9b59b6", "#f39c12", "#1abc9c"]
    all_pairs = sorted(df_detail["pair"].unique())
    for idx, p in enumerate(all_pairs):
        pair_color_map[p] = base_colors[idx % len(base_colors)]

    for protein, ligand in combo_list:
        sub = df_detail[
            (df_detail["protein"] == protein) & (df_detail["ligand"] == ligand)
        ]
        if sub.empty:
            continue

        cluster_ids = sorted(sub["cluster"].unique())
        n_clusters = len(cluster_ids)

        # balanced-N lookup
        combo_rmsd = df_rmsd[
            (df_rmsd["protein"] == protein) & (df_rmsd["ligand"] == ligand)
        ]
        bal_n_map = combo_rmsd.groupby("cluster")["n_poses_per_method"].first().to_dict()

        # closest RMSD lookup  {(cluster, pair): value}
        closest_map = {}
        for _, rw in combo_rmsd.iterrows():
            closest_map[(rw["cluster"], rw["pair"])] = rw["closest_rmsd"]

        n_cols = min(n_clusters, 4)
        n_rows = (n_clusters + n_cols - 1) // n_cols

        fig, axes = plt.subplots(
            n_rows, n_cols,
            figsize=(5.5 * n_cols, 5 * n_rows),
            squeeze=False,
        )
        fig.suptitle(
            f"{protein}  +  {ligand}\n"
            f"Inter-method pairwise RMSD per cluster  "
            f"(balanced, {clust_name})",
            fontsize=13, fontweight="bold", y=1.03,
        )

        for ci, cid in enumerate(cluster_ids):
            ri, co = divmod(ci, n_cols)
            ax = axes[ri][co]
            c_sub = sub[sub["cluster"] == cid]
            pairs_here = sorted(c_sub["pair"].unique())

            if not pairs_here:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                        ha="center", va="center", fontsize=12)
                ax.set_title(f"Cluster {cid}", fontsize=10)
                continue

            data_for_box = [c_sub[c_sub["pair"] == p]["rmsd"].values for p in pairs_here]
            colors = [pair_color_map.get(p, "gray") for p in pairs_here]

            bp = ax.boxplot(
                data_for_box, patch_artist=True, widths=0.5,
                medianprops=dict(color="black", linewidth=1.5),
                showfliers=False,
            )
            for patch, c in zip(bp["boxes"], colors):
                patch.set_facecolor(c)
                patch.set_alpha(0.5)

            # Overlay individual RMSD values
            rng = np.random.default_rng(42)
            for pi, (p, vals) in enumerate(zip(pairs_here, data_for_box)):
                jitter = rng.uniform(-0.15, 0.15, size=len(vals))
                ax.scatter(
                    pi + 1 + jitter, vals,
                    color=pair_color_map.get(p, "gray"),
                    edgecolors="black", linewidth=0.4,
                    s=26, alpha=0.75, zorder=5,
                )
                # Mark closest RMSD with a star
                cl_val = closest_map.get((cid, p))
                if cl_val is not None:
                    ax.scatter(
                        pi + 1, cl_val,
                        marker="*", s=220, color="gold",
                        edgecolors="black", linewidth=0.8, zorder=10,
                        label="closest" if pi == 0 else "",
                    )

            bal_n = bal_n_map.get(cid, "?")
            ax.set_xticks(range(1, len(pairs_here) + 1))
            ax.set_xticklabels(pairs_here, fontsize=8, rotation=35, ha="right")
            ax.set_ylabel("RMSD (Å)", fontsize=9)
            ax.set_title(
                f"Cluster {cid}  ({bal_n} poses/method)",
                fontsize=10, fontweight="bold",
            )
            ax.grid(axis="y", alpha=0.3)

        # Hide unused axes
        for ci in range(n_clusters, n_rows * n_cols):
            r, c = divmod(ci, n_cols)
            axes[r][c].set_visible(False)

        # Legend
        legend_handles = [
            mpatches.Patch(facecolor=pair_color_map[p], edgecolor="black",
                           alpha=0.65, label=p)
            for p in all_pairs if p in sub["pair"].values
        ]
        legend_handles.append(
            plt.Line2D([0], [0], marker="*", color="w", markerfacecolor="gold",
                       markeredgecolor="black", markersize=14,
                       label="Closest RMSD")
        )
        fig.legend(
            handles=legend_handles,
            loc="lower center", ncol=min(len(legend_handles), 6),
            fontsize=9, bbox_to_anchor=(0.5, -0.05),
        )

        plt.tight_layout()
        safe = f"{protein}__{ligand}".replace(" ", "_")
        fig_path = output_dir / f"cluster_rmsd_{safe}.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.show()
        print(f"  Saved: {fig_path}")

    # ====================================================================
    # STEP 4 — Heatmap: closest inter-method RMSD per combo (avg clusters)
    # ====================================================================

    print("\n" + "=" * 120)
    print("CLOSEST INTER-METHOD RMSD HEATMAP  (averaged across clusters)")
    print("=" * 120)

    combo_labels = [f"{p}\n{l}" for p, l in combo_list]
    inter_pairs  = sorted(df_rmsd["pair"].unique())

    heat_data = pd.DataFrame(np.nan, index=combo_labels, columns=inter_pairs)
    agg_heat = (
        df_rmsd.groupby(["protein", "ligand", "pair"])["closest_rmsd"]
        .mean().reset_index()
    )
    for _, rw in agg_heat.iterrows():
        lbl = f"{rw['protein']}\n{rw['ligand']}"
        if lbl in heat_data.index:
            heat_data.loc[lbl, rw["pair"]] = rw["closest_rmsd"]

    fig_hm, ax_h = plt.subplots(
        figsize=(max(6, len(inter_pairs) * 2.5),
                 max(4, len(combo_labels) * 0.75))
    )
    mat_h = heat_data.values.astype(float)
    im_h = ax_h.imshow(mat_h, cmap="RdYlGn_r", aspect="auto")

    ax_h.set_xticks(range(len(inter_pairs)))
    ax_h.set_xticklabels(inter_pairs, fontsize=10, rotation=30, ha="right")
    ax_h.set_yticks(range(len(combo_labels)))
    ax_h.set_yticklabels(combo_labels, fontsize=8)

    vmax_h = np.nanmax(mat_h) if not np.all(np.isnan(mat_h)) else 1.0
    for i in range(len(combo_labels)):
        for j in range(len(inter_pairs)):
            v = mat_h[i, j]
            if np.isnan(v):
                ax_h.text(j, i, "N/A", ha="center", va="center",
                          fontsize=10, color="gray")
            else:
                clr = "white" if v > vmax_h * 0.6 else "black"
                ax_h.text(j, i, f"{v:.1f}", ha="center", va="center",
                          fontsize=10, fontweight="bold", color=clr)

    plt.colorbar(im_h, ax=ax_h, shrink=0.75, label="Closest RMSD (Å)")
    ax_h.set_title(
        "Closest Inter-Method RMSD per Protein-Ligand Combination\n"
        f"(balanced poses, {clust_name} clusters, avg. over clusters)",
        fontsize=13, fontweight="bold",
    )
    plt.tight_layout()
    plt.savefig(output_dir / "closest_rmsd_heatmap.png",
                dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Saved: {output_dir / 'closest_rmsd_heatmap.png'}")

    # ====================================================================
    # STEP 5 — Save to CSV
    # ====================================================================
    df_overview.to_csv(output_dir / "cluster_overview.csv", index=False)
    df_rmsd.to_csv(output_dir / "closest_rmsd_per_cluster.csv", index=False)
    df_detail.to_csv(output_dir / "pairwise_rmsd_detail.csv", index=False)
    print(f"\nSaved: {output_dir / 'cluster_overview.csv'}")
    print(f"Saved: {output_dir / 'closest_rmsd_per_cluster.csv'}")
    print(f"Saved: {output_dir / 'pairwise_rmsd_detail.csv'}")

    # ====================================================================
    # SUMMARY
    # ====================================================================
    print("\n" + "=" * 120)
    print("SUMMARY")
    print("=" * 120)
    n_combos = len(combo_list)
    print(f"  Clustering method   : {clust_name}")
    print(f"  Protein-ligand combos: {n_combos}")
    print(f"  Total clusters       : {len(df_overview)}")

    agg = df_rmsd.groupby("pair")["closest_rmsd"].agg(
        ["mean", "median", "min", "max", "std"]).round(2)
    print(f"\n  Closest RMSD across all combos & clusters:")
    for pair, row in agg.iterrows():
        print(f"    {pair:>12s}  —  mean {row['mean']:.2f} Å,  "
              f"median {row['median']:.2f} Å,  "
              f"min {row['min']:.2f} Å,  max {row['max']:.2f} Å")
    print("=" * 120)

# %%
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import cm
import seaborn as sns

"""
CELL: Comprehensive K-Medoids vs Ward Cluster Comparison
=========================================================
1. Overview table: # clusters per protein-ligand pair for each method
2. Side-by-side cluster composition (poses per tool per cluster per method)
3. Cross-tabulation: which poses land in the SAME vs DIFFERENT cluster
4. Heatmaps of cluster composition & cross-method agreement
5. For multi-tool clusters: RMSD of top-N ranked poses (equal counts)
"""

# ============================================================================
# CONFIGURATION  (derived from docking_directories)
# ============================================================================
TOP_N_POSES = 3  # Number of top-ranked poses per tool for RMSD comparison

TOOL_ABBREV = {m: m[:2].upper() for m in METHODS_LIST}
TOOL_COLORS = {m: METHOD_COLORS[m] for m in METHODS_LIST}
tool_labels = [TOOL_ABBREV[m] for m in METHODS_LIST]

# ============================================================================
# VERIFY PREREQUISITES
# ============================================================================
if 'df_km' not in dir() or df_km.empty:
    raise RuntimeError("K-Medoids results (df_km) not found. Run cell 5 first.")
if 'df_ward' not in dir() or df_ward.empty:
    raise RuntimeError("Ward results (df_ward) not found. Run cell 7 first.")

# ============================================================================
# HELPER: RMSD between two SDF molecules
# ============================================================================

def _load_mol(path):
    try:
        mol = Chem.MolFromMolFile(str(path), removeHs=False, sanitize=False)
        if mol and mol.GetNumConformers() > 0:
            return mol
    except Exception:
        pass
    return None

def _rmsd_pair(mol_a, mol_b):
    """Best-fit RMSD between two molecules (heavy atoms)."""
    if mol_a is None or mol_b is None:
        return np.nan
    try:
        ha = Chem.RemoveHs(mol_a)
        hb = Chem.RemoveHs(mol_b)
    except Exception:
        return np.nan
    if ha.GetNumAtoms() != hb.GetNumAtoms():
        # centroid distance proxy
        try:
            ca = ha.GetConformer(); cb = hb.GetConformer()
            xa = np.array([[ca.GetAtomPosition(i).x, ca.GetAtomPosition(i).y,
                            ca.GetAtomPosition(i).z] for i in range(ha.GetNumAtoms())])
            xb = np.array([[cb.GetAtomPosition(i).x, cb.GetAtomPosition(i).y,
                            cb.GetAtomPosition(i).z] for i in range(hb.GetNumAtoms())])
            return float(np.linalg.norm(xa.mean(0) - xb.mean(0)))
        except Exception:
            return np.nan
    try:
        return float(rdMolAlign.GetBestRMS(ha, hb))
    except Exception:
        return np.nan


# ====================================================================
# STEP 1 — OVERVIEW TABLE: clusters per protein-ligand pair
# ====================================================================

print("=" * 130)
print("STEP 1 — CLUSTER COUNT OVERVIEW PER PROTEIN-LIGAND PAIR")
print("          (K-Medoids vs Hierarchical / Ward)")
print("=" * 130)

overview_rows = []
for (protein, ligand), grp_km in df_km.groupby(["protein", "ligand"]):
    grp_ward = df_ward[(df_ward["protein"] == protein) & (df_ward["ligand"] == ligand)]
    n_km   = grp_km["cluster_kmedoids"].nunique()
    n_ward = grp_ward["cluster_ward"].nunique() if not grp_ward.empty else 0
    n_poses = len(grp_km)    # same poses in both frames
    overview_rows.append({
        "protein": protein, "ligand": ligand,
        "n_poses": n_poses,
        "k_kmedoids": n_km, "k_ward": n_ward,
        "same_k": "Yes" if n_km == n_ward else "No",
    })

df_overview = pd.DataFrame(overview_rows)
print(df_overview.to_string(index=False))

# ====================================================================
# STEP 2 — DETAILED CLUSTER COMPOSITION TABLE
#           (poses per docking tool in every cluster, both methods)
# ====================================================================

print("\n" + "=" * 130)
print("STEP 2 — CLUSTER COMPOSITION: POSES PER DOCKING TOOL IN EACH CLUSTER")
print("=" * 130)

comp_rows = []
for (protein, ligand), grp_km in df_km.groupby(["protein", "ligand"]):
    grp_ward = df_ward[(df_ward["protein"] == protein) & (df_ward["ligand"] == ligand)]

    for c_id in sorted(grp_km["cluster_kmedoids"].unique()):
        sub = grp_km[grp_km["cluster_kmedoids"] == c_id]
        mc = sub["method"].value_counts()
        row_data = {
            "protein": protein, "ligand": ligand,
            "clustering": "K-Medoids", "cluster_id": c_id,
            "n_poses": len(sub),
        }
        for m in METHODS_LIST:
            row_data[TOOL_ABBREV[m]] = mc.get(m, 0)
        row_data["n_tools"] = sum(1 for m in METHODS_LIST if mc.get(m, 0) > 0)
        comp_rows.append(row_data)

    for c_id in sorted(grp_ward["cluster_ward"].unique()):
        sub = grp_ward[grp_ward["cluster_ward"] == c_id]
        mc = sub["method"].value_counts()
        row_data = {
            "protein": protein, "ligand": ligand,
            "clustering": "Ward", "cluster_id": c_id,
            "n_poses": len(sub),
        }
        for m in METHODS_LIST:
            row_data[TOOL_ABBREV[m]] = mc.get(m, 0)
        row_data["n_tools"] = sum(1 for m in METHODS_LIST if mc.get(m, 0) > 0)
        comp_rows.append(row_data)

df_comp = pd.DataFrame(comp_rows)
print(df_comp.to_string(index=False))


# ====================================================================
# STEP 3 — CROSS-METHOD AGREEMENT
#           For each pose: is it in the same or different cluster
#           under K-Medoids vs Ward?
# ====================================================================

print("\n" + "=" * 130)
print("STEP 3 — POSE-LEVEL AGREEMENT: SAME vs DIFFERENT CLUSTER ACROSS METHODS")
print("=" * 130)

# Merge K-Medoids and Ward labels per pose (file_path is unique key)
df_merged = df_km[["file_path", "protein", "ligand", "method", "pose_rank",
                    "cluster_kmedoids"]].merge(
    df_ward[["file_path", "cluster_ward"]], on="file_path", how="inner"
)

# For every pair of poses within the same protein-ligand group, check
# whether they are co-clustered under BOTH methods
cross_rows = []
for (protein, ligand), grp in df_merged.groupby(["protein", "ligand"]):
    n = len(grp)
    if n < 2:
        continue
    km_labels = grp["cluster_kmedoids"].values
    w_labels  = grp["cluster_ward"].values
    methods   = grp["method"].values

    same_both  = 0  # same cluster in BOTH methods
    same_km    = 0  # same in K-Medoids only
    same_ward  = 0  # same in Ward only
    diff_both  = 0  # different cluster in BOTH methods
    total_pairs = 0

    for i in range(n):
        for j in range(i + 1, n):
            s_km = km_labels[i] == km_labels[j]
            s_w  = w_labels[i]  == w_labels[j]
            if s_km and s_w:
                same_both += 1
            elif s_km and not s_w:
                same_km += 1
            elif not s_km and s_w:
                same_ward += 1
            else:
                diff_both += 1
            total_pairs += 1

    cross_rows.append({
        "protein": protein, "ligand": ligand,
        "total_pairs": total_pairs,
        "same_both": same_both,
        "same_KM_only": same_km,
        "same_Ward_only": same_ward,
        "diff_both": diff_both,
        "agreement_pct": round(100 * (same_both + diff_both) / total_pairs, 1)
            if total_pairs else 0,
    })

df_cross = pd.DataFrame(cross_rows)
print(df_cross.to_string(index=False))

# ---- Per-tool breakdown within each cluster -----
print("\n" + "=" * 130)
print("STEP 3b — PER-TOOL SAME/DIFFERENT CLUSTER COUNTS")
print("=" * 130)

tool_cross_rows = []
for (protein, ligand), grp in df_merged.groupby(["protein", "ligand"]):
    for tool in METHODS_LIST:
        t_grp = grp[grp["method"] == tool]
        if len(t_grp) < 2:
            continue
        km_l = t_grp["cluster_kmedoids"].values
        wd_l = t_grp["cluster_ward"].values
        same_b = diff_b = 0
        for i in range(len(t_grp)):
            for j in range(i + 1, len(t_grp)):
                if (km_l[i] == km_l[j]) == (wd_l[i] == wd_l[j]):
                    same_b += 1
                else:
                    diff_b += 1
        tool_cross_rows.append({
            "protein": protein, "ligand": ligand,
            "tool": TOOL_ABBREV[tool],
            "pose_pairs": same_b + diff_b,
            "agree": same_b, "disagree": diff_b,
            "agree_pct": round(100 * same_b / (same_b + diff_b), 1) if (same_b + diff_b) else 0,
        })

df_tool_cross = pd.DataFrame(tool_cross_rows)
if not df_tool_cross.empty:
    print(df_tool_cross.to_string(index=False))
else:
    print("  (Not enough poses per tool to compare)")


# ====================================================================
# STEP 4 — HEATMAPS
# ====================================================================

print("\n" + "=" * 130)
print("STEP 4 — GENERATING HEATMAPS")
print("=" * 130)

combo_list = list(df_overview[["protein", "ligand"]].itertuples(index=False, name=None))
n_combos = len(combo_list)

# ---- 4a  Cluster count comparison bar chart ----
fig_ov, ax_ov = plt.subplots(figsize=(max(8, n_combos * 1.2), 5))
x = np.arange(n_combos)
w = 0.35
ax_ov.bar(x - w/2, df_overview["k_kmedoids"], w, label="K-Medoids", color="#3498db", alpha=0.85)
ax_ov.bar(x + w/2, df_overview["k_ward"],     w, label="Ward",      color="#8e44ad", alpha=0.85)
labels_ov = [f"{p[:18]}\n{l[:15]}" for p, l in combo_list]
ax_ov.set_xticks(x)
ax_ov.set_xticklabels(labels_ov, fontsize=8, rotation=45, ha="right")
ax_ov.set_ylabel("Number of Clusters (k)")
ax_ov.set_title("Optimal k: K-Medoids vs Ward per Protein-Ligand Pair",
                fontsize=14, fontweight="bold")
ax_ov.legend()
ax_ov.grid(axis="y", alpha=0.3)
for i in range(n_combos):
    ax_ov.text(i - w/2, df_overview["k_kmedoids"].iloc[i] + 0.1,
               str(df_overview["k_kmedoids"].iloc[i]), ha="center", fontsize=9, fontweight="bold")
    ax_ov.text(i + w/2, df_overview["k_ward"].iloc[i] + 0.1,
               str(df_overview["k_ward"].iloc[i]), ha="center", fontsize=9, fontweight="bold")
plt.tight_layout()
plt.savefig(output_dir / "comparison_cluster_counts.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"  Saved: {output_dir / 'comparison_cluster_counts.png'}")


# ---- 4b  Composition heatmaps (side-by-side per protein-ligand pair) ----
n_cols_side = 2  # KM | Ward
n_rows_side = n_combos
n_methods = len(METHODS_LIST)

fig_side, axes_side = plt.subplots(n_rows_side, n_cols_side,
                                    figsize=(max(10, 3 * n_methods + 4),
                                             max(3, 2.5 * n_rows_side)),
                                    squeeze=False)
fig_side.suptitle("Cluster Composition: Poses per Docking Tool\n"
                  "Left = K-Medoids · Right = Ward",
                  fontsize=15, fontweight="bold")

for row_idx, (protein, ligand) in enumerate(combo_list):
    for col_idx, (clust_method, clust_col, df_src) in enumerate([
        ("K-Medoids", "cluster_kmedoids", df_km),
        ("Ward",      "cluster_ward",     df_ward),
    ]):
        ax = axes_side[row_idx][col_idx]
        grp = df_src[(df_src["protein"] == protein) & (df_src["ligand"] == ligand)]
        cluster_ids = sorted(grp[clust_col].unique())

        mat = np.zeros((len(cluster_ids), n_methods), dtype=int)
        for ci, cid in enumerate(cluster_ids):
            sub = grp[grp[clust_col] == cid]
            mc = sub["method"].value_counts()
            for mi, m in enumerate(METHODS_LIST):
                mat[ci, mi] = mc.get(m, 0)

        im = ax.imshow(mat, cmap="YlOrRd", aspect="auto",
                       vmin=0, vmax=max(mat.max(), 1))
        ax.set_xticks(range(n_methods))
        ax.set_xticklabels([TOOL_ABBREV[m] for m in METHODS_LIST], fontsize=9)
        ax.set_yticks(range(len(cluster_ids)))
        ax.set_yticklabels([f"C{c}" for c in cluster_ids], fontsize=9)
        for ci in range(len(cluster_ids)):
            for cj in range(n_methods):
                v = int(mat[ci, cj])
                clr = "white" if v > mat.max() * 0.6 else "black"
                ax.text(cj, ci, str(v), ha="center", va="center",
                        fontsize=10, fontweight="bold", color=clr)

        title_str = f"{clust_method}\n{protein[:22]} / {ligand[:18]}"
        ax.set_title(title_str, fontsize=10, fontweight="bold")

        if col_idx == 1:
            plt.colorbar(im, ax=ax, shrink=0.7, label="# Poses")

for idx in range(n_combos, n_rows_side):
    for c in range(n_cols_side):
        axes_side[idx][c].set_visible(False)

plt.tight_layout()
plt.savefig(output_dir / "comparison_composition_heatmaps.png",
            dpi=150, bbox_inches="tight")
plt.show()
print(f"  Saved: {output_dir / 'comparison_composition_heatmaps.png'}")


# ---- 4c  Agreement heatmap per protein-ligand pair ----
if not df_cross.empty:
    fig_ag, ax_ag = plt.subplots(figsize=(10, max(3, n_combos * 0.7 + 2)))
    agree_data = df_cross.set_index(
        df_cross["protein"].str[:20] + " / " + df_cross["ligand"].str[:15]
    )[["same_both", "same_KM_only", "same_Ward_only", "diff_both"]]

    im_ag = ax_ag.imshow(agree_data.values, cmap="Blues", aspect="auto")
    ax_ag.set_xticks(range(4))
    ax_ag.set_xticklabels(["Same Both", "Same KM Only", "Same Ward Only", "Diff Both"],
                          fontsize=10)
    ax_ag.set_yticks(range(len(agree_data)))
    ax_ag.set_yticklabels(agree_data.index, fontsize=9)
    for i in range(len(agree_data)):
        for j in range(4):
            v = int(agree_data.values[i, j])
            clr = "white" if v > agree_data.values.max() * 0.6 else "black"
            ax_ag.text(j, i, str(v), ha="center", va="center",
                       fontsize=10, fontweight="bold", color=clr)
    plt.colorbar(im_ag, ax=ax_ag, shrink=0.7, label="# Pose Pairs")
    ax_ag.set_title("Pose-Pair Agreement: K-Medoids vs Ward\n"
                    "(How many pose pairs land in same/different clusters)",
                    fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_dir / "comparison_agreement_heatmap.png",
                dpi=150, bbox_inches="tight")
    plt.show()
    print(f"  Saved: {output_dir / 'comparison_agreement_heatmap.png'}")


# ====================================================================
# STEP 5 — MULTI-TOOL CLUSTERS: RMSD FOR TOP-N RANKED POSES
# ====================================================================

print("\n" + "=" * 130)
print(f"STEP 5 — RMSD FOR TOP-{TOP_N_POSES} RANKED POSES IN MULTI-TOOL CLUSTERS")
print("          (equal number of poses per tool; clusters with ≥2 tools)")
print("=" * 130)

rmsd_records = []

for clust_method, clust_col, df_src in [
    ("K-Medoids", "cluster_kmedoids", df_km),
    ("Ward",      "cluster_ward",     df_ward),
]:
    for (protein, ligand), grp in df_src.groupby(["protein", "ligand"]):
        for c_id, c_grp in grp.groupby(clust_col):
            tools_in = [t for t in METHODS_LIST if (c_grp["method"] == t).any()]
            if len(tools_in) < 2:
                continue

            # For each tool, take top-N ranked poses
            tool_mols = {}
            for tool in tools_in:
                t_sub = c_grp[c_grp["method"] == tool].sort_values("pose_rank").head(TOP_N_POSES)
                mols = []
                for _, row in t_sub.iterrows():
                    m = _load_mol(row["file_path"])
                    if m is not None:
                        mols.append((row["pose_rank"], m, row["file_path"]))
                if mols:
                    tool_mols[tool] = mols

            active_tools = [t for t in tools_in if t in tool_mols]
            if len(active_tools) < 2:
                continue

            # Use equal number of poses = min across tools
            n_eq = min(len(tool_mols[t]) for t in active_tools)

            # Inter-tool RMSD (all pairs of tools)
            for t1, t2 in combinations(active_tools, 2):
                pair_rmsds = []
                for _, m1, _ in tool_mols[t1][:n_eq]:
                    for _, m2, _ in tool_mols[t2][:n_eq]:
                        r = _rmsd_pair(m1, m2)
                        if not np.isnan(r):
                            pair_rmsds.append(r)
                if pair_rmsds:
                    rmsd_records.append({
                        "clustering": clust_method,
                        "protein": protein, "ligand": ligand,
                        "cluster_id": c_id,
                        "tool_1": TOOL_ABBREV[t1], "tool_2": TOOL_ABBREV[t2],
                        "n_poses_each": n_eq,
                        "type": "inter",
                        "mean_rmsd": round(np.mean(pair_rmsds), 2),
                        "min_rmsd":  round(np.min(pair_rmsds), 2),
                        "max_rmsd":  round(np.max(pair_rmsds), 2),
                    })

            # Intra-tool RMSD (self-consistency)
            for tool in active_tools:
                if n_eq < 2:
                    continue
                pairs = []
                sel = tool_mols[tool][:n_eq]
                for a in range(len(sel)):
                    for b in range(a + 1, len(sel)):
                        r = _rmsd_pair(sel[a][1], sel[b][1])
                        if not np.isnan(r):
                            pairs.append(r)
                if pairs:
                    rmsd_records.append({
                        "clustering": clust_method,
                        "protein": protein, "ligand": ligand,
                        "cluster_id": c_id,
                        "tool_1": TOOL_ABBREV[tool], "tool_2": TOOL_ABBREV[tool],
                        "n_poses_each": n_eq,
                        "type": "intra",
                        "mean_rmsd": round(np.mean(pairs), 2),
                        "min_rmsd":  round(np.min(pairs), 2),
                        "max_rmsd":  round(np.max(pairs), 2),
                    })

df_rmsd_cmp = pd.DataFrame(rmsd_records)

if not df_rmsd_cmp.empty:
    print("\n  Full RMSD table:")
    print(df_rmsd_cmp.to_string(index=False))

    # ---- summary by clustering method ----
    print("\n" + "=" * 130)
    print("RMSD SUMMARY BY CLUSTERING METHOD & TOOL PAIR")
    print("=" * 130)
    rmsd_summary = df_rmsd_cmp.groupby(["clustering", "type", "tool_1", "tool_2"]).agg(
        n_clusters = ("mean_rmsd", "count"),
        avg_rmsd   = ("mean_rmsd", "mean"),
        med_rmsd   = ("mean_rmsd", "median"),
        min_rmsd   = ("min_rmsd",  "min"),
        max_rmsd   = ("max_rmsd",  "max"),
    ).round(2).reset_index()
    print(rmsd_summary.to_string(index=False))

    # ---- 5a  RMSD heatmap per clustering method (aggregated) ----
    fig_r, axes_r = plt.subplots(1, 2, figsize=(8 * n_methods, 6))
    fig_r.suptitle(f"Mean Inter-Tool RMSD (top-{TOP_N_POSES} poses, equal counts)\n"
                   "Diagonal = intra-tool",
                   fontsize=14, fontweight="bold")

    for ax_i, clust_name in enumerate(["K-Medoids", "Ward"]):
        ax = axes_r[ax_i]
        sub = df_rmsd_cmp[df_rmsd_cmp["clustering"] == clust_name]
        mat = np.full((n_methods, n_methods), np.nan)
        for i, t1 in enumerate(tool_labels):
            for j, t2 in enumerate(tool_labels):
                mask = ((sub["tool_1"] == t1) & (sub["tool_2"] == t2)) | \
                       ((sub["tool_1"] == t2) & (sub["tool_2"] == t1))
                vals = sub.loc[mask, "mean_rmsd"]
                if not vals.empty:
                    mat[i, j] = vals.mean()
        vmax_r = np.nanmax(mat) if not np.all(np.isnan(mat)) else 10
        im_r = ax.imshow(mat, cmap="RdYlGn_r", aspect="auto", vmin=0,
                         vmax=max(vmax_r, 3))
        ax.set_xticks(range(n_methods))
        ax.set_xticklabels([TOOL_ABBREV[m] for m in METHODS_LIST], fontsize=10)
        ax.set_yticks(range(n_methods))
        ax.set_yticklabels([TOOL_ABBREV[m] for m in METHODS_LIST], fontsize=10)
        for i in range(n_methods):
            for j in range(n_methods):
                v = mat[i, j]
                if np.isnan(v):
                    ax.text(j, i, "N/A", ha="center", va="center",
                            fontsize=11, color="gray")
                else:
                    lbl = f"{v:.1f}" if i != j else f"{v:.1f}\n(intra)"
                    clr = "white" if v > vmax_r * 0.6 else "black"
                    ax.text(j, i, lbl, ha="center", va="center",
                            fontsize=11, fontweight="bold", color=clr)
        plt.colorbar(im_r, ax=ax, shrink=0.8, label="Mean RMSD (Å)")
        ax.set_title(clust_name, fontsize=13, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_dir / "comparison_rmsd_heatmap.png", dpi=150, bbox_inches="tight")
    plt.show()
    print(f"  Saved: {output_dir / 'comparison_rmsd_heatmap.png'}")


    # ---- 5b  Per protein-ligand RMSD heatmaps (multi-tool clusters) ----
    inter_only = df_rmsd_cmp[df_rmsd_cmp["type"] == "inter"]
    if not inter_only.empty:
        pl_combos = list(inter_only.groupby(["protein", "ligand"]).groups.keys())
        n_pl = len(pl_combos)
        n_cols_pl = min(2, n_pl)
        n_rows_pl = (n_pl + n_cols_pl - 1) // n_cols_pl

        fig_pl, axes_pl = plt.subplots(n_rows_pl, n_cols_pl * 2,
                                        figsize=(7 * n_cols_pl * 2,
                                                 5 * n_rows_pl),
                                        squeeze=False)
        fig_pl.suptitle(f"Per Protein-Ligand: Inter-Tool RMSD in Multi-Tool Clusters\n"
                        f"(top-{TOP_N_POSES} poses, equal counts · Left = K-Medoids, Right = Ward)",
                        fontsize=14, fontweight="bold")

        for pl_idx, (prot, lig) in enumerate(pl_combos):
            row_pl = pl_idx // n_cols_pl
            col_pl_base = (pl_idx % n_cols_pl) * 2  # 2 sub-cols per combo

            for ci, cname in enumerate(["K-Medoids", "Ward"]):
                ax_c = axes_pl[row_pl][col_pl_base + ci]
                sub_c = inter_only[(inter_only["protein"] == prot) &
                                   (inter_only["ligand"] == lig) &
                                   (inter_only["clustering"] == cname)]
                mat_c = np.full((n_methods, n_methods), np.nan)
                np.fill_diagonal(mat_c, 0)
                for _, rw in sub_c.iterrows():
                    i = tool_labels.index(rw["tool_1"])
                    j = tool_labels.index(rw["tool_2"])
                    v = rw["mean_rmsd"]
                    mat_c[i, j] = v
                    mat_c[j, i] = v

                vmax_c = np.nanmax(mat_c) if not np.all(np.isnan(mat_c)) else 5
                im_c = ax_c.imshow(mat_c, cmap="RdYlGn_r", aspect="auto",
                                   vmin=0, vmax=max(vmax_c, 3))
                ax_c.set_xticks(range(n_methods))
                ax_c.set_xticklabels(tool_labels, fontsize=9)
                ax_c.set_yticks(range(n_methods))
                ax_c.set_yticklabels(tool_labels, fontsize=9)
                for ii in range(n_methods):
                    for jj in range(n_methods):
                        vv = mat_c[ii, jj]
                        if np.isnan(vv):
                            ax_c.text(jj, ii, "-", ha="center", va="center",
                                      fontsize=10, color="gray")
                        else:
                            clr = "white" if vv > vmax_c * 0.5 else "black"
                            ax_c.text(jj, ii, f"{vv:.1f}", ha="center", va="center",
                                      fontsize=10, fontweight="bold", color=clr)
                plt.colorbar(im_c, ax=ax_c, shrink=0.7, label="RMSD (Å)")
                ax_c.set_title(f"{cname}\n{prot[:22]} / {lig[:16]}",
                               fontsize=10, fontweight="bold")

        # hide unused
        for extra in range(n_pl, n_rows_pl * n_cols_pl):
            rr = extra // n_cols_pl
            cc = (extra % n_cols_pl) * 2
            axes_pl[rr][cc].set_visible(False)
            axes_pl[rr][cc + 1].set_visible(False)

        plt.tight_layout()
        plt.savefig(output_dir / "comparison_per_pl_rmsd_heatmaps.png",
                    dpi=150, bbox_inches="tight")
        plt.show()
        print(f"  Saved: {output_dir / 'comparison_per_pl_rmsd_heatmaps.png'}")

else:
    print("  No multi-tool clusters found for RMSD calculation.")


# ====================================================================
# STEP 6 — SAVE ALL COMPARISON DATA
# ====================================================================

df_overview.to_csv(output_dir / "comparison_cluster_overview.csv", index=False)
df_comp.to_csv(output_dir / "comparison_cluster_composition.csv", index=False)
df_cross.to_csv(output_dir / "comparison_agreement.csv", index=False)
if not df_rmsd_cmp.empty:
    df_rmsd_cmp.to_csv(output_dir / "comparison_rmsd.csv", index=False)

print("\n" + "=" * 130)
print("ALL COMPARISON DATA SAVED")
print("=" * 130)
print(f"  {output_dir / 'comparison_cluster_overview.csv'}")
print(f"  {output_dir / 'comparison_cluster_composition.csv'}")
print(f"  {output_dir / 'comparison_agreement.csv'}")
if not df_rmsd_cmp.empty:
    print(f"  {output_dir / 'comparison_rmsd.csv'}")

print("\n" + "=" * 130)
print("SUMMARY")
print("=" * 130)
print(f"  Protein-ligand pairs analyzed:  {n_combos}")
print(f"  K-Medoids clusters (total):     {df_comp[df_comp['clustering'] == 'K-Medoids']['cluster_id'].count()}")
print(f"  Ward clusters (total):          {df_comp[df_comp['clustering'] == 'Ward']['cluster_id'].count()}")
same_k = (df_overview["same_k"] == "Yes").sum()
print(f"  Pairs with identical k:         {same_k}/{n_combos}")
if not df_cross.empty:
    avg_agree = df_cross["agreement_pct"].mean()
    print(f"  Mean pose-pair agreement:       {avg_agree:.1f}%")
if not df_rmsd_cmp.empty:
    inter_mean = df_rmsd_cmp[df_rmsd_cmp["type"] == "inter"]["mean_rmsd"].mean()
    print(f"  Mean inter-tool RMSD (all):     {inter_mean:.2f} Å")
print("=" * 130)

# %%
"""
RMSD Method Audit — Atom-by-Atom vs Centroid Fallback vs RDKit Best-Fit
=========================================================================
For every multi-method cluster, load each pose pair and report:
  • atom-by-atom RMSD  (equal heavy-atom counts → direct RMSD)
  • centroid distance   (unequal atom counts → fallback)
  • RDKit GetBestRMS    (equal atom counts + successful alignment)
Broken down by protein, ligand, cluster, and docking method pair.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import rdMolAlign

# ============================================================================
# HELPERS
# ============================================================================

def _load_mol_and_coords(sdf_path):
    """Return (mol_with_Hs_removed, heavy_atom_coords_Nx3) or (None, None)."""
    try:
        mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True, sanitize=False)
        if mol is None:
            mol = Chem.MolFromMolFile(str(sdf_path), removeHs=False, sanitize=False)
            if mol is not None:
                mol = Chem.RemoveHs(mol)
        if mol is None or mol.GetNumConformers() == 0:
            return None, None
        conf = mol.GetConformer()
        coords = np.array([
            [conf.GetAtomPosition(i).x,
             conf.GetAtomPosition(i).y,
             conf.GetAtomPosition(i).z]
            for i in range(mol.GetNumAtoms())
        ])
        return mol, coords
    except Exception:
        return None, None


def classify_rmsd_pair(mol1, coords1, mol2, coords2):
    """
    Classify which RMSD method applies and compute the value.

    Returns dict with keys:
        rmsd_method: 'atom_by_atom' | 'centroid_fallback' | 'failed'
        rdkit_bestfit: True/False  (whether GetBestRMS succeeded)
        rmsd_value: float or NaN
        rdkit_rmsd: float or NaN  (GetBestRMS value, if available)
        n_atoms_1, n_atoms_2: heavy-atom counts
    """
    result = {
        "rmsd_method": "failed",
        "rdkit_bestfit": False,
        "rmsd_value": np.nan,
        "rdkit_rmsd": np.nan,
        "n_atoms_1": 0,
        "n_atoms_2": 0,
    }

    if coords1 is None or coords2 is None:
        return result

    result["n_atoms_1"] = len(coords1)
    result["n_atoms_2"] = len(coords2)

    if len(coords1) == len(coords2):
        # Atom-by-atom RMSD (direct, no alignment)
        d = coords1 - coords2
        result["rmsd_value"] = float(np.sqrt(np.mean(np.sum(d ** 2, axis=1))))
        result["rmsd_method"] = "atom_by_atom"

        # Also try RDKit best-fit (Kabsch alignment + atom mapping)
        if mol1 is not None and mol2 is not None:
            try:
                result["rdkit_rmsd"] = float(rdMolAlign.GetBestRMS(mol1, mol2))
                result["rdkit_bestfit"] = True
            except Exception:
                pass  # alignment failed (e.g., atom-type mismatch)
    else:
        # Centroid distance fallback
        result["rmsd_value"] = float(
            np.linalg.norm(coords1.mean(axis=0) - coords2.mean(axis=0))
        )
        result["rmsd_method"] = "centroid_fallback"

    return result


# ============================================================================
# SELECT DATA SOURCE  (same logic as earlier cells)
# ============================================================================

if 'df_ward' in dir() and not df_ward.empty:
    df_src = df_ward.copy()
    clust_col  = "cluster_ward"
    clust_name = "Ward"
elif 'df_km' in dir() and not df_km.empty:
    df_src = df_km.copy()
    clust_col  = "cluster_kmedoids"
    clust_name = "K-Medoids"
else:
    raise RuntimeError("No clustering results found — run clustering cells first.")

# Use the global METHODS_LIST and TOOL_ABBREV (from Cell 2)
# Already defined: METHODS_LIST = list(docking_directories.keys())
#                  TOOL_ABBREV  = {m: m[:2].upper() for m in METHODS_LIST}

combo_list = sorted(df_src.groupby(["protein", "ligand"]).groups.keys())

print(f"Using {clust_name} clustering")
print(f"Poses: {len(df_src)} | Protein-ligand combos: {len(combo_list)}")
print(f"Methods: {METHODS_LIST}")

# ============================================================================
# AUDIT EVERY INTER- & INTRA-METHOD POSE PAIR
# ============================================================================

print("\n" + "=" * 100)
print("SCANNING ALL POSE PAIRS PER CLUSTER...")
print("=" * 100)

pair_records = []

for protein, ligand in combo_list:
    grp = df_src[(df_src["protein"] == protein) & (df_src["ligand"] == ligand)]
    cluster_ids = sorted(grp[clust_col].unique())

    for cid in cluster_ids:
        c_grp = grp[grp[clust_col] == cid]

        # Load molecules + coords once per cluster
        pose_cache = {}  # idx → (method, rank, mol, coords, path)
        for idx, row in c_grp.iterrows():
            mol, coords = _load_mol_and_coords(row["file_path"])
            pose_cache[idx] = (row["method"], row["pose_rank"], mol, coords, row["file_path"])

        loaded_idxs = [i for i, v in pose_cache.items() if v[3] is not None]

        # All unique pairs within this cluster
        for i, j in combinations(loaded_idxs, 2):
            m1, r1, mol1, c1, p1 = pose_cache[i]
            m2, r2, mol2, c2, p2 = pose_cache[j]

            pair_type = "inter" if m1 != m2 else "intra"
            res = classify_rmsd_pair(mol1, c1, mol2, c2)

            pair_records.append({
                "protein": protein,
                "ligand": ligand,
                "cluster": cid,
                "method_1": TOOL_ABBREV.get(m1, m1),
                "method_2": TOOL_ABBREV.get(m2, m2),
                "pair_type": pair_type,
                "n_atoms_1": res["n_atoms_1"],
                "n_atoms_2": res["n_atoms_2"],
                "rmsd_method": res["rmsd_method"],
                "rmsd_value": res["rmsd_value"],
                "rdkit_bestfit": res["rdkit_bestfit"],
                "rdkit_rmsd": res["rdkit_rmsd"],
            })

df_audit = pd.DataFrame(pair_records)

print(f"\nTotal pose pairs evaluated: {len(df_audit)}")

# ============================================================================
# SUMMARY 1 — Global counts by RMSD method
# ============================================================================

print("\n" + "=" * 100)
print("RMSD METHOD BREAKDOWN (all pairs)")
print("=" * 100)

method_counts = df_audit["rmsd_method"].value_counts()
for rm, cnt in method_counts.items():
    pct = 100 * cnt / len(df_audit)
    print(f"  {rm:22s}: {cnt:6d}  ({pct:5.1f}%)")

rdkit_ok = df_audit["rdkit_bestfit"].sum()
print(f"\n  RDKit GetBestRMS succeeded: {rdkit_ok} / {len(df_audit)} "
      f"({100*rdkit_ok/max(len(df_audit),1):.1f}%)")

# ============================================================================
# SUMMARY 2 — Breakdown by docking-method pair
# ============================================================================

print("\n" + "=" * 100)
print("RMSD METHOD BY DOCKING-TOOL PAIR")
print("=" * 100)

pair_summary = (
    df_audit
    .groupby(["method_1", "method_2", "pair_type", "rmsd_method"])
    .agg(
        n_pairs   = ("rmsd_value", "count"),
        mean_rmsd = ("rmsd_value", "mean"),
        rdkit_ok  = ("rdkit_bestfit", "sum"),
        mean_rdkit = ("rdkit_rmsd", "mean"),
    )
    .reset_index()
)
pair_summary["mean_rmsd"]  = pair_summary["mean_rmsd"].round(2)
pair_summary["mean_rdkit"] = pair_summary["mean_rdkit"].round(2)
pair_summary["rdkit_ok"]   = pair_summary["rdkit_ok"].astype(int)
print(pair_summary.to_string(index=False))

# ============================================================================
# SUMMARY 3 — Per protein-ligand-cluster breakdown
# ============================================================================

print("\n" + "=" * 100)
print("PER PROTEIN / LIGAND / CLUSTER  (multi-method clusters only)")
print("=" * 100)

cluster_summary = (
    df_audit[df_audit["pair_type"] == "inter"]
    .groupby(["protein", "ligand", "cluster"])
    .agg(
        total_pairs     = ("rmsd_value", "count"),
        atom_by_atom    = ("rmsd_method", lambda s: (s == "atom_by_atom").sum()),
        centroid_fb     = ("rmsd_method", lambda s: (s == "centroid_fallback").sum()),
        failed          = ("rmsd_method", lambda s: (s == "failed").sum()),
        rdkit_bestfit   = ("rdkit_bestfit", "sum"),
        mean_direct     = ("rmsd_value", "mean"),
        mean_rdkit      = ("rdkit_rmsd", "mean"),
    )
    .reset_index()
)
cluster_summary["mean_direct"] = cluster_summary["mean_direct"].round(2)
cluster_summary["mean_rdkit"]  = cluster_summary["mean_rdkit"].round(2)
cluster_summary["rdkit_bestfit"] = cluster_summary["rdkit_bestfit"].astype(int)

print(cluster_summary.to_string(index=False))

# ============================================================================
# VISUALIZATION — Stacked bar: RMSD method distribution per tool pair
# ============================================================================

import matplotlib.pyplot as plt

fig, axes = plt.subplots(1, 2, figsize=(16, 5))

for ax, pt in zip(axes, ["inter", "intra"]):
    sub = df_audit[df_audit["pair_type"] == pt]
    if sub.empty:
        ax.set_title(f"No {pt}-method pairs")
        continue

    # Create label for each tool pair
    sub = sub.copy()
    sub["pair_label"] = sub["method_1"] + " vs " + sub["method_2"]
    if pt == "intra":
        sub["pair_label"] = sub["method_1"] + " (self)"

    labels = sorted(sub["pair_label"].unique())
    atom_counts = []
    centroid_counts = []
    failed_counts = []

    for lbl in labels:
        s = sub[sub["pair_label"] == lbl]
        atom_counts.append((s["rmsd_method"] == "atom_by_atom").sum())
        centroid_counts.append((s["rmsd_method"] == "centroid_fallback").sum())
        failed_counts.append((s["rmsd_method"] == "failed").sum())

    x = np.arange(len(labels))
    w = 0.6

    ax.bar(x, atom_counts, w, label="Atom-by-atom RMSD", color="#2ecc71")
    ax.bar(x, centroid_counts, w, bottom=atom_counts, label="Centroid fallback", color="#e67e22")
    ax.bar(x, failed_counts, w,
           bottom=[a+c for a, c in zip(atom_counts, centroid_counts)],
           label="Failed", color="#e74c3c", alpha=0.6)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=9)
    ax.set_ylabel("Number of pose pairs")
    ax.set_title(f"{'Inter' if pt == 'inter' else 'Intra'}-method pairs: RMSD method used")
    ax.legend(fontsize=8)

plt.tight_layout()
output_path = output_dir / "rmsd_method_audit.png"
plt.savefig(output_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nSaved: {output_path}")


