"""Compare docked poses (AutoDock Vina, DiffDock, EquiBind) against the
crystallographic reference poses of the PoseBusters Benchmark Set.

Implements the metrics outlined in the project notes:

    * Symmetry-corrected heavy-atom RMSD vs. crystal ligand
      (no superposition — protein frame is shared across all three methods).
    * Centroid (center-of-mass) distance to crystal ligand.
    * UFF strain energy of the predicted ligand (predicted − ETKDG/UFF minimum).
    * Steric clash count vs. protein heavy atoms (vdW-based, scale 0.75).
    * Native contact recovery (residue-level Jaccard within 4 Å).
    * Protein–ligand interaction-fingerprint recovery via ProLIF
      (Tanimoto similarity vs. crystal IFP; HBond / Hydrophobic /
      π-stacking / ionic / halogen-bond bits).
    * PB-Valid flag re-used from posebusters_filtered_results.csv.

Aggregations reported per docking method:
    * Top-1 success at RMSD ≤ 2 / 5 Å, centroid ≤ 4 Å, PB-Valid + RMSD ≤ 2 Å.
    * Best-of-N (oracle) success at the same thresholds.
    * Median / mean RMSD, RMSD CDF plot.
    * Per receptor-ligand pair best/top-1 RMSD CSV.

----------------------------------------------------------------------
Quick usage
----------------------------------------------------------------------
    # Run with defaults (uses posebusters_filtered_results.csv as the pose index)
    python Scripts/Analysis/posebusters_pose_comparison.py

    # Run on a small subset for testing
    python Scripts/Analysis/posebusters_pose_comparison.py --limit-pairs 5

    # Custom paths
    python Scripts/Analysis/posebusters_pose_comparison.py \
        --pb-csv  posebusters_results/benchmark/dock/posebusters_filtered_results.csv \
        --benchmark-dir "Data/PoseBuster Benchmark Set" \
        --out-dir posebusters_results/benchmark/dock/pose_comparison_report \
        --workers 8

From a Jupyter cell:
    !python Scripts/Analysis/posebusters_pose_comparison.py
"""

from __future__ import annotations
import argparse
import math
import multiprocessing as mp
import os
import re
import subprocess
import tempfile
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdForceFieldHelpers
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

RDLogger.DisableLog("rdApp.*")

try:
    import prolif as plf  # type: ignore
    _HAS_PROLIF = True
except Exception:
    _HAS_PROLIF = False

# ProLIF interaction set (avoid VdWContact + Metal* which require vdW radii
# for every protein element and crash on cofactor metals like Co/Mn).
PLIF_INTERACTIONS = [
    "Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking",
    "Anionic", "Cationic", "CationPi", "PiCation",
    "XBDonor", "XBAcceptor",
]

# ───────────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────────

# vdW radii (Å) — used for clash detection.
_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "NA": 2.27, "MG": 1.73, "K": 2.75, "CA": 2.31, "MN": 2.05,
    "FE": 2.00, "CO": 2.00, "NI": 1.97, "CU": 1.96, "ZN": 1.39,
}
_DEFAULT_VDW = 1.70
CLASH_SCALE = 0.75       # clash if d < scale * (r_lig + r_prot)
CONTACT_CUTOFF = 4.0     # Å — residue counted as contact if any atom within cutoff
RMSD_THRESHOLDS = (1.0, 2.0, 5.0)
CENTROID_THRESHOLD = 4.0

PB_CRITICAL_CHECKS = [
    "mol_pred_loaded", "sanitization", "all_atoms_connected",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "protein-ligand_maximum_distance",
    "minimum_distance_to_protein", "volume_overlap_with_protein",
]

TOOL_LABEL = {
    "autodock": "AutoDock Vina",
    "diffdock": "DiffDock",
    "equibind_guided": "EquiBind (guided)",
}

# ───────────────────────────────────────────────────────────────────
# Helpers
# ───────────────────────────────────────────────────────────────────


def _to_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(bool)
    return (
        s.astype(str).str.strip().str.lower()
        .map({"true": True, "1": True, "1.0": True,
              "false": False, "0": False, "0.0": False, "nan": False, "": False})
        .fillna(False).astype(bool)
    )


def parse_rank(method: str, pose_file: str) -> int:
    """Best-effort rank extraction from output filenames (1 = top)."""
    name = Path(pose_file).name
    if method == "autodock":
        m = re.search(r"_model(\d+)\.sdf$", name)
        return int(m.group(1)) if m else 99
    if method == "diffdock":
        m = re.search(r"rank(\d+)", name)
        return int(m.group(1)) if m else 99
    # equibind: no canonical ranking; prefer "unguided_001" then alphabetical
    m = re.search(r"unguided_(\d+)", name)
    if m:
        return int(m.group(1))
    return 50  # all guided poses share rank ~50 → only oracle metrics meaningful


def load_protein_heavy_atoms(pdb_path: Path) -> tuple[np.ndarray, list[str], list[tuple[str, int]]]:
    """Lightweight PDB parser → (coords [N×3], elements, residue_id list).

    Skips waters and hydrogens. Residue id is (chain, resseq).
    """
    coords, elems, resids = [], [], []
    with open(pdb_path) as fh:
        for ln in fh:
            if not (ln.startswith("ATOM  ") or ln.startswith("HETATM")):
                continue
            resname = ln[17:20].strip()
            if resname in {"HOH", "WAT", "DOD"}:
                continue
            elem = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
            if elem == "H":
                continue
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError:
                continue
            chain = ln[21].strip() or "A"
            try:
                resseq = int(ln[22:26])
            except ValueError:
                continue
            coords.append((x, y, z))
            elems.append(elem)
            resids.append((chain, resseq))
    return np.asarray(coords, dtype=np.float32), elems, resids


def load_first_mol(sdf_path: Path) -> Chem.Mol | None:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    for m in suppl:
        if m is not None:
            try:
                Chem.SanitizeMol(m)
            except Exception:
                try:
                    m.UpdatePropertyCache(strict=False)
                    Chem.GetSymmSSSR(m)
                except Exception:
                    return None
            return m
    return None


def reassign_template(pose: Chem.Mol, template: Chem.Mol) -> Chem.Mol | None:
    """Reassign bond orders from crystal template (handles PDBQT round-trips)."""
    try:
        return AssignBondOrdersFromTemplate(template, pose)
    except Exception:
        return None


def symmetry_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float:
    """Heavy-atom RMSD over all symmetry-equivalent matchings, NO superposition.

    Protein frame is shared across all docking methods, so atom positions are
    directly comparable. We enumerate substructure matches of `ref` in `pose`
    (which captures topological symmetry) and take the minimum coordinate RMSD.
    """
    pose = Chem.RemoveHs(pose)
    ref = Chem.RemoveHs(ref)
    pose_conf = pose.GetConformer()
    ref_conf = ref.GetConformer()
    ref_coords = np.array([list(ref_conf.GetAtomPosition(i))
                           for i in range(ref.GetNumAtoms())])

    matches = pose.GetSubstructMatches(ref, uniquify=False, useChirality=False)
    if not matches:
        # Fall back to identity mapping if substructure search fails.
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return float("nan")
        matches = [tuple(range(ref.GetNumAtoms()))]

    best = math.inf
    for mp_ in matches:
        pose_coords = np.array([list(pose_conf.GetAtomPosition(j)) for j in mp_])
        diff = pose_coords - ref_coords
        rmsd = float(np.sqrt((diff * diff).sum() / len(ref_coords)))
        if rmsd < best:
            best = rmsd
    return best if best != math.inf else float("nan")


def centroid_distance(pose: Chem.Mol, ref: Chem.Mol) -> float:
    pose = Chem.RemoveHs(pose); ref = Chem.RemoveHs(ref)
    pc = pose.GetConformer(); rc = ref.GetConformer()
    a = np.mean([list(pc.GetAtomPosition(i)) for i in range(pose.GetNumAtoms())], axis=0)
    b = np.mean([list(rc.GetAtomPosition(i)) for i in range(ref.GetNumAtoms())], axis=0)
    return float(np.linalg.norm(a - b))


def uff_strain(pose: Chem.Mol) -> float:
    """Strain = E(pose) − E(UFF-minimized pose). Returns NaN on failure."""
    try:
        m = Chem.AddHs(pose, addCoords=True)
        ff = AllChem.UFFGetMoleculeForceField(m)
        if ff is None:
            return float("nan")
        e_now = ff.CalcEnergy()
        m_min = Chem.Mol(m)
        ff_min = AllChem.UFFGetMoleculeForceField(m_min)
        ff_min.Minimize(maxIts=200)
        e_min = ff_min.CalcEnergy()
        return float(e_now - e_min)
    except Exception:
        return float("nan")


def _load_prolif_protein(pdb_path: Path) -> "plf.Molecule | None":
    """Add hydrogens via OpenBabel, load via RDKit, wrap as ProLIF Molecule.

    Returns None if obabel is missing or the PDB cannot be parsed.
    """
    if not _HAS_PROLIF:
        return None
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf:
            tmp = tf.name
        subprocess.run(["obabel", str(pdb_path), "-O", tmp, "-h"],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=120)
        rdmol = Chem.MolFromPDBFile(tmp, removeHs=False, sanitize=False,
                                    proximityBonding=True)
        try:
            os.unlink(tmp)
        except OSError:
            pass
        if rdmol is None:
            return None
        return plf.Molecule(rdmol)
    except Exception:
        return None


def _ligand_for_plif(mol: Chem.Mol) -> "plf.Molecule | None":
    if not _HAS_PROLIF or mol is None:
        return None
    try:
        m = Chem.AddHs(mol, addCoords=True)
        return plf.Molecule.from_rdkit(m)
    except Exception:
        return None


def compute_plif_recovery(crystal_mol: Chem.Mol,
                          pose_mols: list[Chem.Mol],
                          prot_mol: "plf.Molecule | None") -> list[float]:
    """Tanimoto similarity of pose interaction fingerprint vs crystal IFP.

    All ligands are scored in a single ProLIF run so bit indices are aligned.
    Returns a list of floats (NaN where ProLIF could not score the pose).
    """
    n = len(pose_mols)
    if not _HAS_PROLIF or prot_mol is None:
        return [float("nan")] * n

    crystal_lig = _ligand_for_plif(crystal_mol)
    if crystal_lig is None:
        return [float("nan")] * n

    pose_ligs = [_ligand_for_plif(m) for m in pose_mols]
    valid_idx = [i for i, x in enumerate(pose_ligs) if x is not None]
    valid_ligs = [crystal_lig] + [pose_ligs[i] for i in valid_idx]

    try:
        fp = plf.Fingerprint(interactions=PLIF_INTERACTIONS)
        fp.run_from_iterable(valid_ligs, prot_mol, n_jobs=1, progress=False)
        bvs = fp.to_bitvectors()
    except Exception:
        return [float("nan")] * n

    if not bvs:
        return [float("nan")] * n
    crystal_bv = bvs[0]
    out = [float("nan")] * n
    for k, i in enumerate(valid_idx, start=1):
        try:
            out[i] = float(DataStructs.TanimotoSimilarity(crystal_bv, bvs[k]))
        except Exception:
            out[i] = float("nan")
    return out


def clash_and_contacts(pose: Chem.Mol,
                       prot_xyz: np.ndarray,
                       prot_elem: list[str],
                       prot_resid: list[tuple[str, int]]) -> tuple[int, set]:
    """Return (n_clashes, set_of_contact_residue_ids)."""
    pose = Chem.RemoveHs(pose)
    pc = pose.GetConformer()
    lig_xyz = np.array([list(pc.GetAtomPosition(i)) for i in range(pose.GetNumAtoms())])
    lig_elem = [a.GetSymbol().upper() for a in pose.GetAtoms()]

    # Prefilter: bounding-box pruning to avoid O(N×M) when protein is huge.
    if len(prot_xyz) == 0:
        return 0, set()
    lo = lig_xyz.min(0) - 6.0
    hi = lig_xyz.max(0) + 6.0
    mask = ((prot_xyz >= lo).all(1) & (prot_xyz <= hi).all(1))
    p_xyz = prot_xyz[mask]
    p_elem = [e for e, m in zip(prot_elem, mask) if m]
    p_res = [r for r, m in zip(prot_resid, mask) if m]
    if len(p_xyz) == 0:
        return 0, set()

    # Pairwise distances (small after pruning).
    diff = lig_xyz[:, None, :] - p_xyz[None, :, :]
    d2 = (diff * diff).sum(-1)
    d = np.sqrt(d2)

    n_clash = 0
    for i, le in enumerate(lig_elem):
        rl = _VDW.get(le, _DEFAULT_VDW)
        for j, pe in enumerate(p_elem):
            rp = _VDW.get(pe, _DEFAULT_VDW)
            if d[i, j] < CLASH_SCALE * (rl + rp):
                n_clash += 1

    # Native contact residues.
    contact_mask = (d < CONTACT_CUTOFF).any(0)
    contacts = {p_res[j] for j in np.where(contact_mask)[0]}
    return n_clash, contacts


# ───────────────────────────────────────────────────────────────────
# Per-pair worker
# ───────────────────────────────────────────────────────────────────

@dataclass
class PoseRecord:
    method: str
    protein: str
    ligand: str
    pose_file: str
    pose_name: str
    rank: int
    rmsd: float
    centroid_dist: float
    strain_energy: float
    n_clashes: int
    contact_recovery: float
    plif_recovery: float
    pb_valid: bool


def process_pair(args) -> list[dict]:
    pair_key, group_records, benchmark_dir, root = args
    protein, ligand = pair_key
    out: list[dict] = []

    pdb_id = protein  # naming convention: protein == "<PDB>_<CCD>"
    cdir = Path(benchmark_dir) / pdb_id
    crystal_sdf = cdir / f"{pdb_id}_ligand.sdf"
    protein_pdb = cdir / f"{pdb_id}_protein.pdb"
    if not crystal_sdf.exists() or not protein_pdb.exists():
        return out

    crystal = load_first_mol(crystal_sdf)
    if crystal is None:
        return out
    crystal_h = Chem.RemoveHs(crystal)

    prot_xyz, prot_elem, prot_resid = load_protein_heavy_atoms(protein_pdb)

    # Native contact set from crystal pose.
    _, native_contacts = clash_and_contacts(crystal_h, prot_xyz, prot_elem, prot_resid)
    n_native = len(native_contacts) or 1

    # ProLIF protein (loaded once per pair). NaN PLIF if ProLIF/obabel missing.
    prot_plif = _load_prolif_protein(protein_pdb)

    # First pass: load every pose mol so we can batch the PLIF call.
    loaded: list[tuple[dict, Chem.Mol, Path]] = []
    for rec in group_records:
        pose_path = Path(rec["pose_file"])
        if not pose_path.is_absolute():
            pose_path = root / pose_path
        if not pose_path.exists():
            continue
        pose = load_first_mol(pose_path)
        if pose is None:
            continue
        pose = reassign_template(pose, crystal_h) or pose
        loaded.append((rec, pose, pose_path))

    # Batch PLIF: single fingerprint run with crystal + all poses.
    plif_vals = compute_plif_recovery(crystal_h,
                                      [p for _, p, _ in loaded],
                                      prot_plif)

    for (rec, pose, pose_path), plif_val in zip(loaded, plif_vals):
        try:
            rmsd = symmetry_rmsd(pose, crystal_h)
            cdist = centroid_distance(pose, crystal_h)
        except Exception:
            rmsd, cdist = float("nan"), float("nan")

        strain = uff_strain(pose)
        try:
            n_clash, contacts = clash_and_contacts(pose, prot_xyz, prot_elem, prot_resid)
            recov = len(native_contacts & contacts) / n_native
        except Exception:
            n_clash, recov = 0, float("nan")

        out.append(asdict(PoseRecord(
            method=rec["docking_method"],
            protein=protein,
            ligand=ligand,
            pose_file=str(pose_path),
            pose_name=rec.get("pose_name", pose_path.stem),
            rank=parse_rank(rec["docking_method"], str(pose_path)),
            rmsd=rmsd,
            centroid_dist=cdist,
            strain_energy=strain,
            n_clashes=int(n_clash),
            contact_recovery=float(recov),
            plif_recovery=float(plif_val),
            pb_valid=bool(rec["pb_valid"]),
        )))
    return out


# ───────────────────────────────────────────────────────────────────
# Aggregation & plots
# ───────────────────────────────────────────────────────────────────


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Per-method top-1 and oracle (best-of-N) success rates."""
    rows = []
    for method, sub in df.groupby("method"):
        # top-1 = pose with smallest rank per pair (lower is better)
        top1 = sub.sort_values("rank").groupby(["protein", "ligand"]).head(1)
        # oracle = pose with smallest rmsd per pair (drop pairs where all RMSDs are NaN)
        sub_valid = sub.dropna(subset=["rmsd"])
        if len(sub_valid):
            oracle_idx = sub_valid.groupby(["protein", "ligand"])["rmsd"].idxmin()
            oracle = sub_valid.loc[oracle_idx.dropna()]
        else:
            oracle = sub_valid

        n_pairs = top1[["protein", "ligand"]].drop_duplicates().shape[0]

        def _rate(s, mask):
            return float(mask.mean()) if len(s) else float("nan")

        row = {"method": method, "n_pairs": n_pairs, "n_poses": len(sub),
               "median_rmsd_top1": float(top1["rmsd"].median(skipna=True)),
               "median_rmsd_oracle": float(oracle["rmsd"].median(skipna=True))}
        for thr in RMSD_THRESHOLDS:
            row[f"top1_rmsd_le_{thr}A_%"] = 100 * _rate(top1, top1["rmsd"] <= thr)
            row[f"oracle_rmsd_le_{thr}A_%"] = 100 * _rate(oracle, oracle["rmsd"] <= thr)
        row[f"top1_centroid_le_{CENTROID_THRESHOLD}A_%"] = \
            100 * _rate(top1, top1["centroid_dist"] <= CENTROID_THRESHOLD)
        row["top1_pb_valid_and_rmsd2_%"] = \
            100 * _rate(top1, (top1["rmsd"] <= 2.0) & top1["pb_valid"])
        row["oracle_pb_valid_and_rmsd2_%"] = \
            100 * _rate(oracle, (oracle["rmsd"] <= 2.0) & oracle["pb_valid"])
        row["mean_strain"] = float(sub["strain_energy"].mean(skipna=True))
        row["mean_n_clashes"] = float(sub["n_clashes"].mean(skipna=True))
        row["mean_contact_recovery"] = float(sub["contact_recovery"].mean(skipna=True))
        row["mean_plif_recovery"] = float(sub["plif_recovery"].mean(skipna=True))
        rows.append(row)
    return pd.DataFrame(rows).set_index("method").round(2)


def plot_rmsd_cdf(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for method, sub in df.groupby("method"):
        top1 = sub.sort_values("rank").groupby(["protein", "ligand"]).head(1)
        vals = top1["rmsd"].dropna().sort_values().values
        if len(vals) == 0:
            continue
        ys = np.arange(1, len(vals) + 1) / len(vals) * 100
        ax.plot(vals, ys, label=f"{TOOL_LABEL.get(method, method)} (n={len(vals)})")
    ax.axvline(2.0, color="grey", ls="--", alpha=0.6, label="2 Å")
    ax.set_xlim(0, 15)
    ax.set_xlabel("Top-1 heavy-atom RMSD vs crystal (Å)")
    ax.set_ylabel("Cumulative % of pairs")
    ax.set_title("Top-1 RMSD CDF — predicted vs crystal pose")
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_success_bars(summary: pd.DataFrame, out: Path) -> None:
    metrics = ["top1_rmsd_le_2.0A_%", "oracle_rmsd_le_2.0A_%",
               "top1_pb_valid_and_rmsd2_%", "oracle_pb_valid_and_rmsd2_%",
               f"top1_centroid_le_{CENTROID_THRESHOLD}A_%"]
    labels = ["Top-1 RMSD ≤ 2 Å", "Oracle RMSD ≤ 2 Å",
              "Top-1 PB-Valid + RMSD ≤ 2 Å", "Oracle PB-Valid + RMSD ≤ 2 Å",
              f"Top-1 centroid ≤ {CENTROID_THRESHOLD} Å"]
    methods = list(summary.index)
    x = np.arange(len(metrics))
    w = 0.8 / max(len(methods), 1)
    fig, ax = plt.subplots(figsize=(11, 5))
    for i, m in enumerate(methods):
        vals = [summary.loc[m, k] for k in metrics]
        ax.bar(x + (i - (len(methods) - 1) / 2) * w, vals, w,
               label=TOOL_LABEL.get(m, m), edgecolor="black")
        for j, v in enumerate(vals):
            ax.text(x[j] + (i - (len(methods) - 1) / 2) * w, v + 0.5,
                    f"{v:.1f}", ha="center", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("% of receptor-ligand pairs")
    ax.set_title("Pose-quality success rates per docking tool")
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_rmsd_box(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    methods = sorted(df["method"].unique())
    data = []
    for m in methods:
        top1 = (df[df.method == m].sort_values("rank")
                .groupby(["protein", "ligand"]).head(1))
        data.append(top1["rmsd"].dropna().values)
    ax.boxplot(data, tick_labels=[TOOL_LABEL.get(m, m) for m in methods], showmeans=True)
    ax.axhline(2.0, color="grey", ls="--", alpha=0.6)
    ax.set_ylabel("Top-1 RMSD (Å)")
    ax.set_title("Distribution of top-1 RMSD per method")
    ax.grid(axis="y", alpha=0.3); fig.tight_layout()
    fig.savefig(out, dpi=160); plt.close(fig)


def plot_strain_vs_rmsd(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    palette = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind_guided": "#2ca02c"}
    for method, sub in df.groupby("method"):
        top1 = sub.sort_values("rank").groupby(["protein", "ligand"]).head(1)
        ax.scatter(top1["rmsd"], top1["strain_energy"],
                   s=10, alpha=0.5, label=TOOL_LABEL.get(method, method),
                   color=palette.get(method))
    ax.set_xlim(0, 15); ax.set_ylim(0, 200)
    ax.axvline(2.0, color="grey", ls="--", alpha=0.4)
    ax.set_xlabel("Top-1 RMSD (Å)")
    ax.set_ylabel("UFF strain energy (kcal/mol)")
    ax.set_title("Pose strain vs RMSD (top-1)")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    fig.savefig(out, dpi=160); plt.close(fig)


# Metric labels and clipping bounds for the comparison plots.
_METRIC_INFO = [
    ("rmsd",            "RMSD (Å)",                 (0, 15)),
    ("centroid_dist",   "Centroid distance (Å)",    (0, 15)),
    ("contact_recovery", "Residue contact recovery", (0, 1)),
    ("plif_recovery",   "PLIF Tanimoto vs crystal", (0, 1)),
    ("strain_energy",   "UFF strain (kcal/mol)",    (0, 200)),
    ("n_clashes",       "Heavy-atom clashes",       (0, 60)),
]


def _top1_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    return (df.sort_values("rank")
              .groupby(["method", "protein", "ligand"]).head(1))


def plot_metric_distributions(df: pd.DataFrame, out: Path) -> None:
    """Per-tool violin/strip distributions for every comparison metric (top-1)."""
    top1 = _top1_per_pair(df)
    methods = sorted(top1["method"].unique())
    palette = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind_guided": "#2ca02c"}
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, (key, label, (lo, hi)) in zip(axes.flat, _METRIC_INFO):
        data, kept_methods = [], []
        for m in methods:
            v = top1.loc[top1.method == m, key].dropna().values
            if len(v) == 0:
                continue
            data.append(np.clip(v, lo, hi))
            kept_methods.append(m)
        if not data:
            ax.text(0.5, 0.5, f"No data for\n{label}",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=10, color="grey")
            ax.set_title(label); ax.set_xticks([]); ax.set_yticks([])
            continue
        parts = ax.violinplot(data, showmedians=True, widths=0.8)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(palette.get(kept_methods[i], "grey"))
            body.set_alpha(0.6)
        ax.set_xticks(range(1, len(kept_methods) + 1))
        ax.set_xticklabels([TOOL_LABEL.get(m, m) for m in kept_methods],
                           rotation=12, ha="right", fontsize=9)
        ax.set_ylabel(label)
        ax.set_title(label)
        ax.grid(axis="y", alpha=0.3)
        if key == "rmsd":
            ax.axhline(2.0, color="red", ls="--", alpha=0.5)
    fig.suptitle("Top-1 metric distributions per docking tool", fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_top_vs_bottom(df: pd.DataFrame, out: Path, n: int = 10) -> None:
    """Per tool: compare metrics for the N best- and N worst-RMSD top-1 poses."""
    top1 = _top1_per_pair(df).dropna(subset=["rmsd"])
    methods = sorted(top1["method"].unique())
    metrics = [k for k, _, _ in _METRIC_INFO if k != "rmsd"]
    metric_labels = {k: lab for k, lab, _ in _METRIC_INFO}
    metric_clip = {k: rng for k, _, rng in _METRIC_INFO}

    fig, axes = plt.subplots(len(methods), len(metrics),
                             figsize=(3.0 * len(metrics), 3.0 * len(methods)),
                             sharey=False)
    if len(methods) == 1:
        axes = np.array([axes])

    rows_summary = []
    for r, m in enumerate(methods):
        sub = top1[top1.method == m].sort_values("rmsd")
        best = sub.head(n)
        worst = sub.tail(n)
        rows_summary.append({
            "method": m, "n": n,
            "best_mean_rmsd": float(best["rmsd"].mean()),
            "worst_mean_rmsd": float(worst["rmsd"].mean()),
        })
        for c, key in enumerate(metrics):
            ax = axes[r, c]
            lo, hi = metric_clip[key]
            b = np.clip(best[key].dropna().values, lo, hi)
            w = np.clip(worst[key].dropna().values, lo, hi)
            if len(b) == 0 and len(w) == 0:
                ax.text(0.5, 0.5, "No data", ha="center", va="center",
                        transform=ax.transAxes, fontsize=10, color="grey")
                ax.set_xticks([]); ax.set_yticks([])
                if r == 0:
                    ax.set_title(metric_labels[key], fontsize=10)
                if c == 0:
                    ax.set_ylabel(TOOL_LABEL.get(m, m), fontsize=10)
                continue
            positions, vals, colors = [], [], []
            if len(b):
                positions.append(1); vals.append(b); colors.append("#2ca02c")
            if len(w):
                positions.append(2); vals.append(w); colors.append("#d62728")
            bp = ax.boxplot(vals, positions=positions, widths=0.6,
                            patch_artist=True, showmeans=True)
            for patch, col in zip(bp["boxes"], colors):
                patch.set_facecolor(col); patch.set_alpha(0.6)
            ax.set_xticks([1, 2]); ax.set_xticklabels([f"Best {n}", f"Worst {n}"], fontsize=8)
            if r == 0:
                ax.set_title(metric_labels[key], fontsize=10)
            if c == 0:
                ax.set_ylabel(TOOL_LABEL.get(m, m), fontsize=10)
            ax.grid(axis="y", alpha=0.3)
    fig.suptitle(f"Best vs worst {n} top-1 poses per tool (by RMSD)", fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
    return pd.DataFrame(rows_summary)


def plot_cross_tool_pairwise(df: pd.DataFrame, out: Path) -> None:
    """Per-pair top-1 RMSD scatter between every tool combination + delta hist."""
    top1 = _top1_per_pair(df).dropna(subset=["rmsd"])
    pivot = top1.pivot_table(index=["protein", "ligand"],
                             columns="method", values="rmsd")
    methods = list(pivot.columns)
    if len(methods) < 2:
        return
    n = len(methods)
    fig, axes = plt.subplots(n, n, figsize=(3.2 * n, 3.2 * n))
    if n == 1:
        axes = np.array([[axes]])
    for i, mi in enumerate(methods):
        for j, mj in enumerate(methods):
            ax = axes[i, j]
            if i == j:
                vals = pivot[mi].dropna()
                ax.hist(np.clip(vals, 0, 15), bins=30, color="#888", edgecolor="black")
                ax.axvline(2.0, color="red", ls="--", alpha=0.6)
                ax.set_title(TOOL_LABEL.get(mi, mi), fontsize=10)
                ax.set_xlim(0, 15)
                ax.set_xlabel("Top-1 RMSD (Å)")
            elif i < j:
                # Delta histogram in upper triangle: row method − col method.
                pair = pivot[[mi, mj]].dropna()
                if len(pair) == 0:
                    ax.set_axis_off(); continue
                d = (pair[mi] - pair[mj]).values
                ax.hist(np.clip(d, -10, 10), bins=30, color="#4c72b0", edgecolor="black")
                ax.axvline(0.0, color="red", ls="--", alpha=0.6)
                wins_i = int((pair[mi] < pair[mj]).sum())
                ax.set_title(f"{TOOL_LABEL.get(mi,mi)[:8]} − {TOOL_LABEL.get(mj,mj)[:8]}\n"
                             f"row better in {wins_i}/{len(pair)}", fontsize=8)
            else:
                # Scatter in lower triangle.
                pair = pivot[[mj, mi]].dropna()
                if len(pair) == 0:
                    ax.set_axis_off(); continue
                ax.scatter(pair[mj], pair[mi], s=10, alpha=0.5)
                lim = 15
                ax.plot([0, lim], [0, lim], color="red", ls="--", alpha=0.5)
                ax.set_xlim(0, lim); ax.set_ylim(0, lim)
                ax.set_xlabel(TOOL_LABEL.get(mj, mj), fontsize=9)
                ax.set_ylabel(TOOL_LABEL.get(mi, mi), fontsize=9)
            ax.grid(alpha=0.3)
    fig.suptitle("Cross-tool comparison — per-pair top-1 RMSD", fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_cross_tool_metric_means(summary: pd.DataFrame, out: Path) -> None:
    """Side-by-side bar chart of mean physical/interaction metrics per tool."""
    cols = ["mean_contact_recovery", "mean_plif_recovery",
            "mean_n_clashes", "mean_strain"]
    labels = ["Mean contact recovery", "Mean PLIF Tanimoto",
              "Mean clashes / pose", "Mean UFF strain (kcal/mol, log)"]
    methods = list(summary.index)
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.2))
    palette = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind_guided": "#2ca02c"}
    for ax, col, lab in zip(axes, cols, labels):
        if col not in summary.columns:
            ax.text(0.5, 0.5, f"{lab}\nnot available", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="grey")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(lab, fontsize=10)
            continue
        raw_vals = [summary.loc[m, col] for m in methods]
        # Skip column if every value is NaN (e.g. ProLIF unavailable).
        if all(pd.isna(v) for v in raw_vals):
            ax.text(0.5, 0.5, f"{lab}\nall NaN", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="grey")
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(lab, fontsize=10)
            continue
        vals = [0.0 if pd.isna(v) else v for v in raw_vals]
        colors = [palette.get(m, "grey") for m in methods]
        if col == "mean_strain":
            vals_plot = [max(v, 1e-3) for v in vals]
            ax.set_yscale("log")
        else:
            vals_plot = vals
        ax.bar([TOOL_LABEL.get(m, m) for m in methods], vals_plot,
               color=colors, edgecolor="black")
        for i, v in enumerate(raw_vals):
            txt = "NaN" if pd.isna(v) else f"{v:.2g}"
            ax.text(i, vals_plot[i], txt, ha="center",
                    va="bottom", fontsize=9)
        ax.set_title(lab, fontsize=10)
        ax.tick_params(axis="x", rotation=15)
        ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Cross-tool comparison of physical & interaction metrics", fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────────

def _build_pose_index(pb_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(pb_csv, low_memory=False)
    needed = {"docking_method", "protein", "ligand", "pose_file", "pose_name"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"PB CSV missing columns: {missing}")
    bool_df = pd.DataFrame({c: _to_bool(df[c])
                            for c in PB_CRITICAL_CHECKS if c in df.columns})
    df["pb_valid"] = bool_df.all(axis=1) if not bool_df.empty else False
    df["docking_method"] = df["docking_method"].astype(str).str.lower()
    return df[["docking_method", "protein", "ligand", "pose_file",
               "pose_name", "pb_valid"]]


def main() -> None:
    # How to run from a notebook cell:
    #   !python Scripts/Analysis/posebusters_pose_comparison.py
    #
    # Optional flags:
    #   --pb-csv         PoseBusters filtered-results CSV (defines pose index).
    #   --benchmark-dir  Folder containing per-PDB crystal SDFs / protein PDBs.
    #   --out-dir        Output folder for CSVs + plots.
    #   --workers        Number of worker processes (default = cpu_count).
    #   --limit-pairs    Process only the first N receptor-ligand pairs (debug).
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pb-csv", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "posebusters_filtered_results.csv"))
    ap.add_argument("--benchmark-dir", type=Path,
                    default=Path("Data/PoseBuster Benchmark Set"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "pose_comparison_report"))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--limit-pairs", type=int, default=0,
                    help="Process only the first N pairs (debug).")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    root = Path.cwd()

    print(f"Loading pose index from: {args.pb_csv}")
    idx = _build_pose_index(args.pb_csv)
    print(f"  poses: {len(idx):,}  pairs: "
          f"{idx.groupby(['protein','ligand']).ngroups:,}")

    # Group per receptor-ligand pair so each worker loads protein & crystal once.
    grouped: dict[tuple[str, str], list[dict]] = {}
    for rec in idx.to_dict("records"):
        grouped.setdefault((rec["protein"], rec["ligand"]), []).append(rec)
    pairs = list(grouped.items())
    if args.limit_pairs:
        pairs = pairs[:args.limit_pairs]
    print(f"  pairs to process: {len(pairs):,} (workers={args.workers})")

    work = [((p, l), recs, str(args.benchmark_dir.resolve()), root)
            for (p, l), recs in pairs]

    out_records: list[dict] = []
    if args.workers > 1:
        with mp.Pool(args.workers) as pool:
            for i, batch in enumerate(pool.imap_unordered(process_pair, work, chunksize=4), 1):
                out_records.extend(batch)
                if i % 25 == 0 or i == len(work):
                    print(f"  [{i}/{len(work)}] pairs scored "
                          f"({len(out_records):,} poses)")
    else:
        for i, w in enumerate(work, 1):
            out_records.extend(process_pair(w))
            if i % 25 == 0 or i == len(work):
                print(f"  [{i}/{len(work)}] pairs scored")

    if not out_records:
        print("No pose records produced. Check input paths.")
        return

    df = pd.DataFrame(out_records)
    per_pose_csv = args.out_dir / "per_pose_metrics.csv"
    df.to_csv(per_pose_csv, index=False)
    print(f"\nWrote per-pose metrics → {per_pose_csv}  ({len(df):,} rows)")

    summary = aggregate(df)
    summary_csv = args.out_dir / "summary_per_method.csv"
    summary.to_csv(summary_csv)
    print(f"Wrote per-method summary → {summary_csv}\n")
    print(summary.to_string())

    # Per-pair top-1 RMSD pivot for inspection.
    top1 = (df.sort_values("rank")
              .groupby(["method", "protein", "ligand"]).head(1))
    pivot = top1.pivot_table(index=["protein", "ligand"],
                             columns="method", values="rmsd")
    pivot.to_csv(args.out_dir / "top1_rmsd_per_pair.csv")

    # Plots.
    plot_rmsd_cdf(df, args.out_dir / "01_rmsd_cdf.png")
    plot_success_bars(summary, args.out_dir / "02_success_bars.png")
    plot_rmsd_box(df, args.out_dir / "03_rmsd_boxplot.png")
    plot_strain_vs_rmsd(df, args.out_dir / "04_strain_vs_rmsd.png")
    plot_metric_distributions(df, args.out_dir / "05_metric_distributions.png")
    extremes = plot_top_vs_bottom(df, args.out_dir / "06_top_vs_bottom_extremes.png", n=10)
    if extremes is not None:
        extremes.to_csv(args.out_dir / "top_vs_bottom_summary.csv", index=False)
    plot_cross_tool_pairwise(df, args.out_dir / "07_cross_tool_pairwise.png")
    plot_cross_tool_metric_means(summary, args.out_dir / "08_cross_tool_metric_means.png")
    print(f"\nFigures + CSVs → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
