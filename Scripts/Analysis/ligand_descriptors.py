"""Shared per-ligand descriptor loader for the docking pose-cluster reports.

Reuses the cached RDKit feature table built by ``PoseBusters_DataSet_Analysis.ipynb``
(``PoseBusters_Benchmark_Analysis/ligand_protein_features.csv``) so every report
describes ligands identically — Lipinski Rule-of-Five compliance, the
physicochemical descriptors, and a PCA of ligand chemical space. Keyed by ``entry``
= ``<PDBID>_<CCD>``, which matches the ``protein`` / complex-id column used by the
docking pose tables.

Used by ``pose_cluster_crystal_pocket_report.py`` (descriptor → RMSD-to-crystal)
and ``orai_pose_cluster_report.py`` (descriptor → cross-tool agreement / validity).
"""
from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

DEFAULT_FEATURES_CSV = "PoseBusters_Benchmark_Analysis/ligand_protein_features.csv"

# Continuous / count descriptors used for correlation + PCA (binary flags dropped),
# mirroring LIG_NUM / PCA_FEATURES in the dataset notebook.
PCA_FEATURES: List[str] = [
    "heavy_atoms", "mw", "rot_bonds", "hbd", "hba", "tpsa", "logp",
    "n_rings", "n_arom_rings", "fsp3", "qed", "n_heteroatoms", "n_halogens", "n_stereo",
]

# spelled-out axis labels (no bare abbreviations) — matches the figure conventions
NICE = {
    "heavy_atoms": "Heavy atoms", "mw": "Molecular weight (Da)",
    "rot_bonds": "Rotatable bonds", "hbd": "H-bond donors", "hba": "H-bond acceptors",
    "tpsa": "Topological polar surface area (Å²)", "logp": "Calculated logP (lipophilicity)",
    "n_rings": "Rings", "n_arom_rings": "Aromatic rings", "fsp3": "Fraction sp³ carbons",
    "qed": "Quantitative estimate of drug-likeness (QED)", "n_heteroatoms": "Heteroatoms",
    "n_halogens": "Halogen atoms", "n_stereo": "Stereocentres",
    "ro5_violations": "Lipinski Rule-of-Five violations",
    "startconf_rmsd": "Start-conformer to crystal RMSD (Å)",
}


def load_descriptors(csv: str = DEFAULT_FEATURES_CSV, key: str = "entry") -> Optional[pd.DataFrame]:
    """Per-ligand descriptor table indexed by complex id, with Lipinski Ro5 added.

    Returns None if the cache CSV is missing (so callers can degrade gracefully).
    """
    p = Path(csv)
    if not p.exists():
        return None
    df = pd.read_csv(p)
    if key not in df.columns:
        return None
    # Lipinski Rule-of-Five — identical definition to the dataset notebook.
    df["ro5_violations"] = (
        (df["mw"] > 500).astype(int) + (df["logp"] > 5).astype(int)
        + (df["hbd"] > 5).astype(int) + (df["hba"] > 10).astype(int)
    )
    df["ro5_pass"] = df["ro5_violations"] <= 1        # standard "<=1 violation" rule
    return df.set_index(key)


def add_pca(df: pd.DataFrame, features: Optional[List[str]] = None,
            var_target: float = 0.80) -> Tuple[pd.DataFrame, pd.DataFrame, np.ndarray, int]:
    """Append PC1..PCk score columns (standardised PCA of `features`) to `df`.

    Returns (df_with_PCs, loadings, explained_variance_ratio, n_pcs_for_var_target).
    Rows with any missing feature are dropped from the fit but kept in df (PC = NaN).
    """
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    feats = [f for f in (features or PCA_FEATURES) if f in df.columns]
    X = df[feats].apply(pd.to_numeric, errors="coerce").dropna()
    if len(X) < 3 or len(feats) < 2:
        return df.copy(), pd.DataFrame(), np.array([]), 0
    Xz = StandardScaler().fit_transform(X)
    pca = PCA().fit(Xz)
    scores = pca.transform(Xz)
    cols = [f"PC{i + 1}" for i in range(scores.shape[1])]
    sc = pd.DataFrame(scores, index=X.index, columns=cols)
    loadings = pd.DataFrame(pca.components_.T * np.sqrt(pca.explained_variance_),
                            index=feats, columns=cols)
    evr = pca.explained_variance_ratio_
    n_keep = int(np.searchsorted(np.cumsum(evr), var_target) + 1)
    return df.join(sc), loadings, evr, n_keep


def nice(col: str) -> str:
    return NICE.get(col, col)
