# %%
from __future__ import annotations

import csv
import inspect
import json
import os
import shutil
import subprocess
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import yaml
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

# %% [markdown]
# # Config

# %%
# Paths
workspace_root = Path("/workspace")
receptors_dir = workspace_root / "master_dev/Orai"

# List of ligand directories — each gets its own output folder, error log, and summary
drugs_dirs = [
    Path("/root/workspace/storage/ligands_sdf_large_approvedA100"),
    Path("/root/workspace/storage/ligands_sdf_small_symmetric_approved/")
]

# DiffDock paths
DIFFDOCK_DIR = Path("/workspace/DiffDock")
DIFFDOCK_CONDA_ENV = "diffdock"

# Python executable — MUST be the diffdock conda env, not the notebook kernel
DIFFDOCK_PYTHON = os.path.expanduser("~/anaconda3/envs/diffdock/bin/python")

# Number of samples (poses) DiffDock generates per complex
NUM_SAMPLES: int = 10

# Device for DiffDock (cpu or cuda:0)
DIFFDOCK_DEVICE = "cuda:0"  # Change to "cpu" if you don't have a compatible GPU

# Per-complex time limit in seconds (kills subprocess if exceeded)
TIMEOUT_PER_COMPLEX: int = 300  # 15 minutes per complex; set to 0 for no limit

# Number of denoising steps DiffDock uses (default 20; fewer = faster, lower quality)
INFERENCE_STEPS: int = 20  # Try 10 for faster runs, 20 for best quality

# Base output directory — each drugs_dir gets a subfolder named after it
DIFFDOCK_OUTPUT_BASE = workspace_root / "diffdock_results"
DIFFDOCK_OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

# Overwrite settings
OVERWRITE_EXISTING = False  # Set to True to re-dock existing combinations
OVERWRITE_ERROR_LOG = True  # Set to True to retry previously failed ligands

# Error log filename (one per output directory)
ERROR_LOG_FILENAME = "failed_docking.json"

# ── Performance settings ──
# Use CSV batch mode: submits all complexes in one DiffDock call (much faster)
USE_BATCH_MODE = True

# Batch size: how many complexes per DiffDock subprocess call in batch mode
# Sweet spot: 50 balances model-load overhead vs preprocessing time
# Set to 0 for all-in-one (risky with large runs — long preprocessing stalls GPU)
BATCH_SIZE: int = 15

# Number of parallel workers for single-complex fallback mode (USE_BATCH_MODE=False)
# Only used when USE_BATCH_MODE is False. Set to 1 for sequential.
# On GPU, keep at 1 (GPU can't run multiple DiffDock in parallel).
# On CPU, try 2-4 depending on cores.
NUM_WORKERS: int = 1

# Validate directories
for drugs_dir in drugs_dirs:
    if not drugs_dir.exists():
        raise FileNotFoundError(f"Ligand directory missing: {drugs_dir}")
if not receptors_dir.exists():
    raise FileNotFoundError(f"Receptor directory missing: {receptors_dir}")

print("=" * 80)
print("DiffDock Batch Docking Configuration")
print("=" * 80)
print(f"Number of samples per complex: {NUM_SAMPLES}")
print(f"Inference steps: {INFERENCE_STEPS}")
print(f"Timeout per complex: {TIMEOUT_PER_COMPLEX}s ({TIMEOUT_PER_COMPLEX/60:.0f} min)" if TIMEOUT_PER_COMPLEX else "Timeout per complex: unlimited")
print(f"DiffDock directory: {DIFFDOCK_DIR}")
print(f"DiffDock conda env: {DIFFDOCK_CONDA_ENV}")
print(f"DiffDock device: {DIFFDOCK_DEVICE}")
print(f"Output base directory: {DIFFDOCK_OUTPUT_BASE}")
print(f"Overwrite existing: {OVERWRITE_EXISTING}")
print(f"Batch mode: {USE_BATCH_MODE} (batch size: {'all' if BATCH_SIZE == 0 else BATCH_SIZE})")
if not USE_BATCH_MODE:
    print(f"Parallel workers: {NUM_WORKERS}")
print(f"\nLigand directories ({len(drugs_dirs)}):")
for d in drugs_dirs:
    n_files = len(list(d.glob("*.sdf")) + list(d.glob("*.mol2")) + list(d.glob("*.pdb")))
    output_dir = DIFFDOCK_OUTPUT_BASE / d.name
    print(f"  • {d.name:30s} → output: {output_dir.name}/  ({n_files} ligand files)")
print()

# Validate DiffDock installation
if not DIFFDOCK_DIR.exists():
    print(f"⚠️  WARNING: DiffDock directory not found at {DIFFDOCK_DIR}")
    print("   Please update DIFFDOCK_DIR to point to your DiffDock installation")
else:
    print(f"✓ DiffDock directory found")

# %% [markdown]
# # Helper Functions & Data Classes

# %%
@dataclass
class DockingResult:
    """Result of docking a protein-ligand combination."""
    protein_name: str
    ligand_name: str
    protein_path: Path
    ligand_path: Path
    output_dir: Path
    status: str  # "success", "failed", "skipped"
    num_poses: int = 0
    pose_files: List[Path] = field(default_factory=list)
    error_message: str = ""
    elapsed_time: float = 0.0
    
    def to_dict(self) -> dict:
        return {
            "protein_name": self.protein_name,
            "ligand_name": self.ligand_name,
            "protein_path": str(self.protein_path),
            "ligand_path": str(self.ligand_path),
            "output_dir": str(self.output_dir),
            "status": self.status,
            "num_poses": self.num_poses,
            "pose_files": [str(p) for p in self.pose_files],
            "error_message": self.error_message,
            "elapsed_time": self.elapsed_time,
        }


def get_file_stem(path: Path) -> str:
    """Get clean filename stem without extension."""
    return path.stem.replace("_ligand", "").replace("_protein", "")


def collect_files(root: Path, extensions: List[str]) -> List[Path]:
    """Return files with given extensions located directly inside root (no recursion)."""
    files = []
    for ext in extensions:
        files.extend(sorted(p for p in root.glob(f"*{ext}") if p.is_file()))
    return sorted(set(files))


def extract_confidence_from_filename(filename: str) -> float:
    """Extract confidence score from DiffDock output filename (e.g., rank1_confidence-0.85.sdf)."""
    try:
        if "confidence" in filename.lower():
            parts = filename.lower().split("confidence")
            if len(parts) > 1:
                conf_str = parts[1].replace("-", "").replace("_", "").replace(".sdf", "")
                return float(conf_str)
    except:
        pass
    return 0.0


def extract_rank_from_filename(filename: str) -> int:
    """Extract rank from DiffDock output filename (e.g., rank1_confidence-0.85.sdf)."""
    try:
        if "rank" in filename.lower():
            parts = filename.lower().split("rank")
            if len(parts) > 1:
                rank_str = ""
                for c in parts[1]:
                    if c.isdigit():
                        rank_str += c
                    else:
                        break
                if rank_str:
                    return int(rank_str)
    except:
        pass
    return 0


def prepare_protein_for_diffdock(input_pdb: Path, output_dir: Path) -> Path:
    """
    Prepare a protein PDB file for DiffDock by fixing common issues.
    """
    if input_pdb.is_dir():
        pdb_files = sorted(input_pdb.glob("*.pdb"))
        if not pdb_files:
            raise FileNotFoundError(f"No .pdb files found in directory: {input_pdb}")
        print(f"Preparing {len(pdb_files)} PDB files from {input_pdb.name}/")
        for pdb in pdb_files:
            prepared = _prepare_single_protein(pdb, output_dir)
            print(f"  ✓ {prepared.name}")
        return output_dir
    else:
        return _prepare_single_protein(input_pdb, output_dir)


def _prepare_single_protein(input_pdb: Path, output_dir: Path) -> Path:
    """Prepare a single PDB file for DiffDock."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pdb = output_dir / f"{input_pdb.stem}_prepared.pdb"
    
    residue_mapping = {
        'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
        'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
    }
    
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                res_name = line[17:20].strip()
                if res_name in residue_mapping:
                    new_res = residue_mapping[res_name]
                    line = line[:17] + f"{new_res:>3}" + line[20:]
                f_out.write(line)
            elif line.startswith(('END', 'TER')):
                f_out.write(line)
    
    return output_pdb


def prepare_ligand_for_diffdock(input_sdf: Path, output_dir: Path) -> Path:
    """Prepare a ligand SDF file for DiffDock by fixing valence/charge issues."""
    output_dir.mkdir(parents=True, exist_ok=True)
    output_sdf = output_dir / f"{input_sdf.stem}_prepared.sdf"
    
    default_valence = {
        'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 3,
        'F': 1, 'Cl': 1, 'Br': 1, 'I': 1, 'B': 3,
    }
    
    try:
        with open(input_sdf, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 5:
            return input_sdf
        
        counts_line = lines[3]
        try:
            num_atoms = int(counts_line[:3])
            num_bonds = int(counts_line[3:6])
        except ValueError:
            return input_sdf
        
        if len(lines) < 4 + num_atoms + num_bonds:
            return input_sdf
        
        atom_symbols = []
        explicit_valences = []
        
        for i in range(4, 4 + num_atoms):
            line = lines[i]
            parts = line.split()
            if len(parts) >= 10:
                symbol = parts[3]
                atom_symbols.append(symbol)
                try:
                    val = int(parts[9])
                except (ValueError, IndexError):
                    val = 0
                explicit_valences.append(val)
            else:
                atom_symbols.append('?')
                explicit_valences.append(0)
        
        bond_counts = [0] * num_atoms
        for i in range(4 + num_atoms, 4 + num_atoms + num_bonds):
            line = lines[i]
            parts = line.split()
            if len(parts) >= 3:
                try:
                    a1 = int(parts[0]) - 1
                    a2 = int(parts[1]) - 1
                    bond_order = int(parts[2])
                except ValueError:
                    continue
                if 0 <= a1 < num_atoms and 0 <= a2 < num_atoms:
                    bond_counts[a1] += bond_order
                    bond_counts[a2] += bond_order
        
        charges_to_add = {}
        for idx in range(num_atoms):
            symbol = atom_symbols[idx]
            if symbol not in default_valence:
                continue
            dv = default_valence[symbol]
            bc = bond_counts[idx]
            ev = explicit_valences[idx]
            if bc > dv and ev > 0:
                charge = bc - dv
                charges_to_add[idx + 1] = charge
        
        if not charges_to_add:
            return input_sdf
        
        charge_atoms = list(charges_to_add.items())
        chg_line = f"M  CHG  {len(charge_atoms)}"
        for atom_idx, charge in charge_atoms:
            chg_line += f"  {atom_idx:3d}  {charge:3d}"
        chg_line += "\n"
        
        with open(output_sdf, 'w') as f:
            for line in lines:
                if line.strip() == "M  END":
                    f.write(chg_line)
                f.write(line)
        
        fixed_atoms = ", ".join(f"atom #{k} ({atom_symbols[k-1]}) -> +{v}" for k, v in charges_to_add.items())
        print(f"  ⚠ Fixed missing charges in {input_sdf.name}: {fixed_atoms}")
        
        return output_sdf
    
    except Exception as e:
        print(f"  ⚠ Could not parse {input_sdf.name} for charge fixing ({e}), using original")
        return input_sdf


def load_error_log(output_dir: Path) -> Dict[str, dict]:
    """Load the error log for a given output directory."""
    error_log_path = output_dir / ERROR_LOG_FILENAME
    if error_log_path.exists():
        try:
            with open(error_log_path, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, Exception) as e:
            print(f"  ⚠ Could not read error log {error_log_path}: {e}")
    return {}


def save_error_log(output_dir: Path, error_log: Dict[str, dict]):
    """Save the error log for a given output directory."""
    error_log_path = output_dir / ERROR_LOG_FILENAME
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(error_log_path, 'w') as f:
        json.dump(error_log, f, indent=2)


def record_failure(output_dir: Path, combo_name: str, result: "DockingResult"):
    """Record a single failure to the error log."""
    error_log = load_error_log(output_dir)
    error_log[combo_name] = {
        "error": result.error_message,
        "timestamp": datetime.now().isoformat(),
        "protein": result.protein_name,
        "ligand": result.ligand_name,
        "elapsed_time": round(result.elapsed_time, 2),
    }
    save_error_log(output_dir, error_log)


def is_known_failure(output_dir: Path, combo_name: str, error_log: Optional[Dict] = None) -> bool:
    """Check if a combination is in the error log."""
    if error_log is None:
        error_log = load_error_log(output_dir)
    return combo_name in error_log

# %% [markdown]
# # Docking Log & Molecular Properties

# %%
DOCKING_LOG_FILENAME = "docking_log.csv"

DOCKING_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name",
    "status", "error_reason", "num_poses", "elapsed_time_s",
    # Ligand properties
    "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula",
    # Protein properties
    "prot_num_residues", "prot_num_atoms", "prot_num_chains",
    # Config & hardware
    "num_samples", "inference_steps", "device", "gpu_model",
    # Best pose confidence
    "best_confidence",
]


def get_gpu_model() -> str:
    """Return the GPU model name via nvidia-smi, or 'N/A' if unavailable."""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=gpu_name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=5,
        )
        names = [line.strip() for line in result.stdout.strip().splitlines() if line.strip()]
        return names[0] if names else "N/A"
    except Exception:
        return "N/A"


def get_ligand_properties(sdf_path: Path) -> Dict:
    """Extract molecular properties from a ligand SDF file using RDKit."""
    props = {
        "lig_molecular_weight": None, "lig_heavy_atoms": None,
        "lig_total_atoms": None, "lig_rotatable_bonds": None,
        "lig_num_rings": None, "lig_aromatic_rings": None,
        "lig_hbd": None, "lig_hba": None, "lig_tpsa": None,
        "lig_logp": None, "lig_formula": None,
    }
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        mol = next(iter(supplier), None)
        if mol is None:
            return props
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

        props["lig_molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        props["lig_heavy_atoms"] = mol.GetNumHeavyAtoms()
        props["lig_total_atoms"] = mol.GetNumAtoms()
        props["lig_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ri = mol.GetRingInfo()
        props["lig_num_rings"] = ri.NumRings()
        props["lig_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        props["lig_hbd"] = Lipinski.NumHDonors(mol)
        props["lig_hba"] = Lipinski.NumHAcceptors(mol)
        props["lig_tpsa"] = round(Descriptors.TPSA(mol), 2)
        props["lig_logp"] = round(Descriptors.MolLogP(mol), 2)
        props["lig_formula"] = rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        print(f"  ⚠ Could not compute ligand properties for {sdf_path.name}: {e}")
    return props


def get_protein_properties(pdb_path: Path) -> Dict:
    """Extract basic protein properties by parsing PDB ATOM records."""
    props = {"prot_num_residues": None, "prot_num_atoms": None, "prot_num_chains": None}
    try:
        residues = set()
        chains = set()
        atom_count = 0
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1
                    chain_id = line[21]
                    res_seq = line[22:27].strip()
                    residues.add((chain_id, res_seq))
                    chains.add(chain_id)
        props["prot_num_residues"] = len(residues)
        props["prot_num_atoms"] = atom_count
        props["prot_num_chains"] = len(chains)
    except Exception as e:
        print(f"  ⚠ Could not parse protein {pdb_path.name}: {e}")
    return props


def ligand_props_valid(props: Dict) -> bool:
    """Return True if at least one key ligand property was computed successfully."""
    return props.get("lig_heavy_atoms") is not None


def precompute_properties(
    proteins: List[Path],
    ligands: List[Path],
) -> Tuple[Dict[Path, Dict], Dict[Path, Dict]]:
    """Pre-compute molecular properties for all proteins and ligands before docking.

    Returns (lig_props_cache, prot_props_cache) keyed by file path.
    """
    lig_props_cache: Dict[Path, Dict] = {}
    prot_props_cache: Dict[Path, Dict] = {}

    print("Pre-computing ligand properties...")
    for lig in ligands:
        lig_props_cache[lig] = get_ligand_properties(lig)
    print(f"  {len(lig_props_cache)} ligands analysed")

    print("Pre-computing protein properties...")
    for prot in proteins:
        prot_props_cache[prot] = get_protein_properties(prot)
    print(f"  {len(prot_props_cache)} proteins analysed")

    return lig_props_cache, prot_props_cache


# Detect GPU once at definition time — reused for every log row
_GPU_MODEL = get_gpu_model()
print(f"GPU detected: {_GPU_MODEL}")


def _build_log_row(
    result: "DockingResult",
    lig_props: Dict,
    prot_props: Dict,
    samples: int,
    steps: int,
    device: str,
) -> Dict:
    """Build a single row dict for the docking log CSV."""
    best_conf = 0.0
    if result.pose_files:
        confs = [extract_confidence_from_filename(p.name) for p in result.pose_files]
        best_conf = max(confs) if confs else 0.0

    row = {
        "timestamp": datetime.now().isoformat(),
        "combo_name": f"{result.ligand_name}__{result.protein_name}",
        "protein_name": result.protein_name,
        "ligand_name": result.ligand_name,
        "status": result.status,
        "error_reason": result.error_message if result.status == "failed" else "",
        "num_poses": result.num_poses,
        "elapsed_time_s": round(result.elapsed_time, 2),
        "best_confidence": round(best_conf, 4),
        "num_samples": samples,
        "inference_steps": steps,
        "device": device,
        "gpu_model": _GPU_MODEL,
    }
    row.update(lig_props)
    row.update(prot_props)
    return row


def init_docking_log(output_dir: Path) -> Tuple[Path, set]:
    """Initialise the docking log CSV and return (log_path, existing_combo_names).

    Creates the file with a header row if it does not yet exist.
    """
    log_path = output_dir / DOCKING_LOG_FILENAME
    existing_combos: set = set()

    if log_path.exists():
        try:
            df = pd.read_csv(log_path)
            existing_combos = set(df["combo_name"])
        except Exception:
            pass
    else:
        pd.DataFrame(columns=DOCKING_LOG_COLUMNS).to_csv(log_path, index=False)

    return log_path, existing_combos


def append_log_rows(
    results: List["DockingResult"],
    log_path: Path,
    existing_combos: set,
    lig_props_cache: Dict[Path, Dict],
    prot_props_cache: Dict[Path, Dict],
    samples: int = NUM_SAMPLES,
    steps: int = INFERENCE_STEPS,
    device: str = DIFFDOCK_DEVICE,
) -> None:
    """Append log rows for a batch of results to the CSV, skipping duplicates."""
    rows = []
    for r in results:
        combo = f"{r.ligand_name}__{r.protein_name}"
        if combo in existing_combos:
            continue
        lig_props = lig_props_cache.get(r.ligand_path, get_ligand_properties(r.ligand_path))
        prot_props = prot_props_cache.get(r.protein_path, get_protein_properties(r.protein_path))
        rows.append(_build_log_row(r, lig_props, prot_props, samples, steps, device))
        existing_combos.add(combo)

    if not rows:
        return

    df = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)

    df.to_csv(log_path, mode='a', header=False, index=False)
    print("Docking log functions defined.")

    print(f"  📝 Appended {len(rows)} entries to {log_path.name} ({len(existing_combos)} total)")


# %% [markdown]
# # Summary & Reporting

# %%
def generate_summary(results: List[DockingResult]) -> Dict:
    """Generate summary of docking results."""
    summary = {
        "timestamp": datetime.now().isoformat(),
        "configuration": {
            "num_samples": NUM_SAMPLES,
            "device": DIFFDOCK_DEVICE,
            "diffdock_dir": str(DIFFDOCK_DIR),
        },
        "overall": {
            "total_combinations": len(results),
            "successful": sum(1 for r in results if r.status == "success"),
            "failed": sum(1 for r in results if r.status == "failed"),
            "skipped": sum(1 for r in results if r.status == "skipped"),
            "total_poses": sum(r.num_poses for r in results),
            "total_time_seconds": sum(r.elapsed_time for r in results),
        },
        "by_protein": defaultdict(lambda: {"combinations": 0, "poses": 0, "success": 0}),
        "by_ligand": defaultdict(lambda: {"combinations": 0, "poses": 0, "success": 0}),
        "combinations": [],
    }
    
    for r in results:
        summary["by_protein"][r.protein_name]["combinations"] += 1
        summary["by_protein"][r.protein_name]["poses"] += r.num_poses
        if r.status == "success":
            summary["by_protein"][r.protein_name]["success"] += 1
            
        summary["by_ligand"][r.ligand_name]["combinations"] += 1
        summary["by_ligand"][r.ligand_name]["poses"] += r.num_poses
        if r.status == "success":
            summary["by_ligand"][r.ligand_name]["success"] += 1
        
        summary["combinations"].append({
            "protein": r.protein_name,
            "ligand": r.ligand_name,
            "status": r.status,
            "num_poses": r.num_poses,
            "time_seconds": round(r.elapsed_time, 2),
            "error": r.error_message if r.error_message else None,
        })
    
    summary["by_protein"] = dict(summary["by_protein"])
    summary["by_ligand"] = dict(summary["by_ligand"])
    
    return summary


def print_summary(summary: Dict):
    """Print formatted docking summary."""
    print("\n" + "=" * 80)
    print("DIFFDOCK DOCKING SUMMARY")
    print("=" * 80)
    
    overall = summary["overall"]
    print(f"\nOverall Statistics:")
    print(f"  Total combinations: {overall['total_combinations']}")
    print(f"  Successful: {overall['successful']}")
    print(f"  Failed: {overall['failed']}")
    print(f"  Skipped: {overall['skipped']}")
    print(f"  Total poses generated: {overall['total_poses']}")
    print(f"  Total time: {overall['total_time_seconds']:.1f}s ({overall['total_time_seconds']/60:.1f} min)")
    
    print("\n" + "-" * 80)
    print("Results by Protein:")
    print("-" * 80)
    for protein, stats in sorted(summary["by_protein"].items()):
        print(f"  {protein:40s} | Combos: {stats['combinations']:3d} | Poses: {stats['poses']:4d}")
    
    print("\n" + "-" * 80)
    print("Results by Ligand:")
    print("-" * 80)
    for ligand, stats in sorted(summary["by_ligand"].items()):
        print(f"  {ligand:40s} | Combos: {stats['combinations']:3d} | Poses: {stats['poses']:4d}")

# %% [markdown]
# # DiffDock Inference Engine

# %%
def _build_diffdock_config(samples: int, steps: int, output_dir: Path) -> Path:
    """Create a runtime config YAML that merges DiffDock defaults with our settings."""
    default_yaml = DIFFDOCK_DIR / "default_inference_args.yaml"
    with open(default_yaml) as f:
        config = yaml.safe_load(f)

    config["samples_per_complex"] = samples
    config["inference_steps"] = steps
    config["actual_steps"] = steps - 1

    output_dir.mkdir(parents=True, exist_ok=True)
    custom_yaml = output_dir / "custom_inference_args.yaml"
    with open(custom_yaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return custom_yaml


def _collect_ranked_poses(output_dir: Path) -> List[Path]:
    """Collect only the actual ranked pose SDF files from DiffDock output."""
    output_sdfs = list(output_dir.rglob("rank*_confidence*.sdf"))
    output_sdfs = [s for s in output_sdfs if s.stat().st_size > 0]
    return sorted(output_sdfs, key=lambda p: extract_rank_from_filename(p.name))


def _write_batch_csv(
    combos: List[Tuple[str, Path, Path]],
    csv_path: Path,
) -> Path:
    """Write a protein_ligand_csv for DiffDock batch inference.

    combos: list of (complex_name, prepared_protein_path, prepared_ligand_path)
    """
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['complex_name', 'protein_path', 'ligand_description', 'protein_sequence'])
        for name, prot, lig in combos:
            writer.writerow([name, str(prot.absolute()), str(lig.absolute()), ''])
    return csv_path


def _run_diffdock_batch_subprocess(
    csv_path: Path,
    out_dir: Path,
    config_yaml: Path,
    samples: int,
    device: str,
    timeout: int,
) -> subprocess.CompletedProcess:
    """Run DiffDock in CSV batch mode — streams stdout/stderr to log files."""
    cmd = [
        DIFFDOCK_PYTHON, "-m", "inference",
        "--config", str(config_yaml),
        "--protein_ligand_csv", str(csv_path),
        "--out_dir", str(out_dir),
        "--samples_per_complex", str(samples),
        "--no_final_step_noise",
    ]

    env = os.environ.copy()
    if device == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""

    effective_timeout = timeout if timeout > 0 else None

    log_dir = out_dir / "_batch_logs"
    log_dir.mkdir(exist_ok=True)
    stdout_log = log_dir / f"batch_{csv_path.stem}_stdout.log"
    stderr_log = log_dir / f"batch_{csv_path.stem}_stderr.log"

    with open(stdout_log, "w") as fout, open(stderr_log, "w") as ferr:
        result = subprocess.run(
            cmd,
            stdout=fout,
            stderr=ferr,
            text=True,
            timeout=effective_timeout,
            cwd=str(DIFFDOCK_DIR),
            env=env,
        )

    result.stdout = stdout_log.read_text() if stdout_log.stat().st_size < 5_000_000 else f"[see {stdout_log}]"
    result.stderr = stderr_log.read_text() if stderr_log.stat().st_size < 5_000_000 else f"[see {stderr_log}]"

    return result


def create_protein_ligand_csv(
    proteins: List[Path],
    ligands: List[Path],
    output_csv: Path,
    prepare: bool = True,
) -> List[Tuple[str, Path, Path]]:
    """Create a CSV file for DiffDock batch inference."""
    combinations = []

    prepared_protein_cache: Dict[Path, Path] = {}
    prepared_ligand_cache: Dict[Path, Path] = {}
    if prepare:
        prepared_prot_dir = output_csv.parent / "prepared_proteins"
        for protein in proteins:
            if protein not in prepared_protein_cache:
                prepared_protein_cache[protein] = prepare_protein_for_diffdock(protein, prepared_prot_dir)
                print(f"  Prepared protein: {prepared_protein_cache[protein].name}")

        prepared_lig_dir = output_csv.parent / "prepared_ligands"
        for ligand in ligands:
            if ligand not in prepared_ligand_cache:
                if ligand.suffix.lower() == ".sdf":
                    prepared_ligand_cache[ligand] = prepare_ligand_for_diffdock(ligand, prepared_lig_dir)
                else:
                    prepared_ligand_cache[ligand] = ligand

    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['complex_name', 'protein_path', 'ligand_description', 'protein_sequence'])

        for protein in proteins:
            protein_name = get_file_stem(protein)
            protein_to_use = prepared_protein_cache.get(protein, protein)
            for ligand in ligands:
                ligand_name = get_file_stem(ligand)
                ligand_to_use = prepared_ligand_cache.get(ligand, ligand)
                complex_name = f"{ligand_name}__{protein_name}"

                writer.writerow([
                    complex_name,
                    str(protein_to_use.absolute()),
                    str(ligand_to_use.absolute()),
                    ''
                ])
                combinations.append((complex_name, protein, ligand))

    print(f"Created CSV with {len(combinations)} combinations: {output_csv}")
    return combinations


def _run_diffdock_subprocess(
    protein_path: Path,
    ligand_path: Path,
    output_dir: Path,
    output_base_dir: Path,
    samples: int,
    steps: int,
    device: str,
    timeout: int,
) -> Tuple[subprocess.CompletedProcess, List[Path]]:
    """Run DiffDock subprocess for a single complex."""
    config_yaml = _build_diffdock_config(samples, steps, output_base_dir)

    cmd = [
        DIFFDOCK_PYTHON, "-m", "inference",
        "--config", str(config_yaml),
        "--protein_path", str(protein_path),
        "--ligand_description", str(ligand_path),
        "--out_dir", str(output_dir),
        "--samples_per_complex", str(samples),
        "--no_final_step_noise",
    ]

    env = os.environ.copy()
    if device == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""

    effective_timeout = timeout if timeout > 0 else None

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=effective_timeout,
        cwd=str(DIFFDOCK_DIR),
        env=env,
    )

    output_sdfs = _collect_ranked_poses(output_dir)
    return result, output_sdfs


def run_diffdock_single(
    protein_path: Path,
    ligand_path: Path,
    output_dir: Path,
    output_base_dir: Path,
    samples: int = NUM_SAMPLES,
    steps: int = INFERENCE_STEPS,
    device: str = DIFFDOCK_DEVICE,
    timeout: int = TIMEOUT_PER_COMPLEX,
) -> Tuple[bool, List[Path], str]:
    """Run DiffDock inference for a single protein-ligand pair (with CUDA OOM fallback)."""
    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result, output_sdfs = _run_diffdock_subprocess(
            protein_path, ligand_path, output_dir, output_base_dir,
            samples, steps, device, timeout,
        )

        if output_sdfs:
            return True, output_sdfs, ""

        stderr = result.stderr or ""
        if ("OutOfMemoryError" in stderr or "CUDA out of memory" in stderr) and device != "cpu":
            print(f"  ⚠ CUDA OOM — retrying on CPU...")
            if output_dir.exists():
                shutil.rmtree(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            result, output_sdfs = _run_diffdock_subprocess(
                protein_path, ligand_path, output_dir, output_base_dir,
                samples, steps, "cpu", timeout,
            )

            if output_sdfs:
                return True, output_sdfs, ""

            return False, [], f"CPU fallback also failed. stderr: {result.stderr[-500:] if result.stderr else ''}"

        return False, [], f"No output SDF files found. stderr: {stderr[-500:]}"

    except subprocess.TimeoutExpired:
        return False, [], f"Docking timeout ({timeout}s)"
    except Exception as e:
        return False, [], str(e)


def _prepare_all_inputs(
    proteins: List[Path],
    ligands: List[Path],
    output_base_dir: Path,
) -> Tuple[Dict[Path, Path], Dict[Path, Path]]:
    """Prepare all proteins and ligands upfront (cached). Returns two caches."""
    prepared_prot_dir = output_base_dir / "prepared_proteins"
    prepared_lig_dir = output_base_dir / "prepared_ligands"

    prot_cache: Dict[Path, Path] = {}
    for protein in proteins:
        if protein not in prot_cache:
            prot_cache[protein] = prepare_protein_for_diffdock(protein, prepared_prot_dir)
    print(f"  Prepared {len(prot_cache)} proteins (cached in {prepared_prot_dir.name}/)")

    lig_cache: Dict[Path, Path] = {}
    for ligand in ligands:
        if ligand not in lig_cache:
            if ligand.suffix.lower() == ".sdf":
                lig_cache[ligand] = prepare_ligand_for_diffdock(ligand, prepared_lig_dir)
            else:
                lig_cache[ligand] = ligand
    print(f"  Prepared {len(lig_cache)} ligands (cached in {prepared_lig_dir.name}/)")

    return prot_cache, lig_cache

# %% [markdown]
# # Batch Docking Pipeline

# %%
def dock_single_combination(
    protein_path: Path,
    ligand_path: Path,
    output_base_dir: Path,
    prot_cache: Optional[Dict[Path, Path]] = None,
    lig_cache: Optional[Dict[Path, Path]] = None,
    samples: int = NUM_SAMPLES,
    steps: int = INFERENCE_STEPS,
    device: str = DIFFDOCK_DEVICE,
    timeout: int = TIMEOUT_PER_COMPLEX,
    error_log: Optional[Dict[str, dict]] = None,
) -> DockingResult:
    """Dock a single protein-ligand combination (uses caches if provided)."""
    protein_name = get_file_stem(protein_path)
    ligand_name = get_file_stem(ligand_path)
    combo_name = f"{ligand_name}__{protein_name}"
    output_dir = output_base_dir / combo_name
    
    result = DockingResult(
        protein_name=protein_name,
        ligand_name=ligand_name,
        protein_path=protein_path,
        ligand_path=ligand_path,
        output_dir=output_dir,
        status="pending",
    )
    
    if not OVERWRITE_ERROR_LOG and error_log and combo_name in error_log:
        prev = error_log[combo_name]
        result.status = "skipped"
        result.error_message = f"Previously failed: {prev.get('error', 'unknown')}"
        return result
    
    if not OVERWRITE_EXISTING and output_dir.exists():
        existing_sdfs = _collect_ranked_poses(output_dir)
        if existing_sdfs:
            result.status = "skipped"
            result.num_poses = len(existing_sdfs)
            result.pose_files = existing_sdfs
            return result
    
    start_time = time.time()
    
    if prot_cache and protein_path in prot_cache:
        prepared_protein = prot_cache[protein_path]
    else:
        prepared_dir = output_base_dir / "prepared_proteins"
        prepared_protein = prepare_protein_for_diffdock(protein_path, prepared_dir)
    
    if lig_cache and ligand_path in lig_cache:
        prepared_ligand = lig_cache[ligand_path]
    else:
        prepared_ligand_dir = output_base_dir / "prepared_ligands"
        if ligand_path.suffix.lower() == ".sdf":
            prepared_ligand = prepare_ligand_for_diffdock(ligand_path, prepared_ligand_dir)
        else:
            prepared_ligand = ligand_path
    
    success, pose_files, error = run_diffdock_single(
        protein_path=prepared_protein,
        ligand_path=prepared_ligand,
        output_dir=output_dir,
        output_base_dir=output_base_dir,
        samples=samples,
        steps=steps,
        device=device,
        timeout=timeout,
    )
    
    result.elapsed_time = time.time() - start_time
    
    if success and pose_files:
        result.status = "success"
        result.num_poses = len(pose_files)
        result.pose_files = pose_files
    else:
        result.status = "failed"
        result.error_message = error
        record_failure(output_base_dir, combo_name, result)
    
    return result


def _run_batch_csv_mode(
    combos_to_dock: List[Tuple[str, Path, Path, Path, Path]],
    output_dir: Path,
    prot_cache: Dict[Path, Path],
    lig_cache: Dict[Path, Path],
    samples: int,
    steps: int,
    device: str,
    timeout: int,
    error_log: Dict[str, dict],
    log_path: Optional[Path] = None,
    existing_combos: Optional[set] = None,
    lig_props_cache: Optional[Dict[Path, Dict]] = None,
    prot_props_cache: Optional[Dict[Path, Dict]] = None,
) -> List[DockingResult]:
    """
    Run DiffDock in CSV batch mode — one subprocess handles many complexes.
    Writes docking log incrementally after each batch.
    
    combos_to_dock: list of (combo_name, protein_path, ligand_path, orig_prot, orig_lig)
    """
    results = []
    
    if not combos_to_dock:
        return results
    
    batch_size = BATCH_SIZE if BATCH_SIZE > 0 else len(combos_to_dock)
    batches = [combos_to_dock[i:i + batch_size] for i in range(0, len(combos_to_dock), batch_size)]
    
    config_yaml = _build_diffdock_config(samples, steps, output_dir)
    
    combos_completed = 0
    total_combos = len(combos_to_dock)
    
    progress_file = output_dir / "_batch_logs" / "_progress.txt"
    progress_file.parent.mkdir(exist_ok=True)
    
    for batch_idx, batch in enumerate(batches):
        batch_start_idx = batch_idx * batch_size
        batch_end_idx = min(batch_start_idx + len(batch), total_combos)
        
        print(f"\n  ── Batch {batch_idx + 1}/{len(batches)} ({len(batch)} complexes) "
              f"[{batch_start_idx + 1}–{batch_end_idx} of {total_combos}] ──")
        
        for j, (combo_name, _, _, orig_prot, orig_lig) in enumerate(batch):
            prot_name = get_file_stem(orig_prot)
            lig_name = get_file_stem(orig_lig)
            print(f"    {batch_start_idx + j + 1:>4}/{total_combos}  {lig_name} + {prot_name}")
        
        progress_file.write_text(
            f"{combos_completed}/{total_combos}|batch {batch_idx + 1}/{len(batches)}|"
            + "|".join(f"{get_file_stem(op)}+{get_file_stem(ol)}" for _, _, _, op, ol in batch)
        )
        
        batch_start = time.time()
        batch_results = []
        
        csv_rows = []
        for combo_name, prot_path, lig_path, orig_prot, orig_lig in batch:
            prepared_prot = prot_cache.get(orig_prot, prot_path)
            prepared_lig = lig_cache.get(orig_lig, lig_path)
            csv_rows.append((combo_name, prepared_prot, prepared_lig))
        
        csv_path = output_dir / f"_batch_{batch_idx}.csv"
        _write_batch_csv(csv_rows, csv_path)
        
        batch_timeout = timeout * len(batch) if timeout > 0 else 0
        
        try:
            proc = _run_diffdock_batch_subprocess(
                csv_path=csv_path,
                out_dir=output_dir,
                config_yaml=config_yaml,
                samples=samples,
                device=device,
                timeout=batch_timeout,
            )
            
            batch_elapsed = time.time() - batch_start
            stderr = proc.stderr or ""
            
            batch_log_dir = output_dir / "_batch_logs"
            batch_log_dir.mkdir(exist_ok=True)
            if stderr:
                (batch_log_dir / f"batch_{batch_idx}_stderr.log").write_text(stderr)
            if proc.stdout:
                (batch_log_dir / f"batch_{batch_idx}_stdout.log").write_text(proc.stdout)
            
            if proc.returncode != 0:
                print(f"  ⚠ Batch subprocess exit code {proc.returncode}")
                error_lines = [l for l in stderr.splitlines() if any(k in l.lower() for k in ['error', 'exception', 'traceback', 'failed', 'cannot', 'no module'])]
                if error_lines:
                    for line in error_lines[:5]:
                        print(f"    {line.strip()[:200]}")
            
            if ("OutOfMemoryError" in stderr or "CUDA out of memory" in stderr) and device != "cpu":
                # Halve the batch and retry on GPU; fall back to CPU only at size 1
                sub_batches = [batch]
                retry_size = max(1, len(batch) // 2)
                while retry_size >= 1:
                    sub_batches = [batch[i:i + retry_size] for i in range(0, len(batch), retry_size)]
                    print(f"  ⚠ CUDA OOM on batch {batch_idx + 1} — retrying with sub-batch size {retry_size} on GPU...")
                    oom_again = False
                    for sb_idx, sub_batch in enumerate(sub_batches):
                        sub_csv_rows = [(cn, prot_cache.get(op, pp), lig_cache.get(ol, lp))
                                        for cn, pp, lp, op, ol in sub_batch]
                        sub_csv_path = output_dir / f"_batch_{batch_idx}_sub{sb_idx}.csv"
                        _write_batch_csv(sub_csv_rows, sub_csv_path)
                        sub_timeout = timeout * len(sub_batch) if timeout > 0 else 0
                        try:
                            proc = _run_diffdock_batch_subprocess(
                                csv_path=sub_csv_path,
                                out_dir=output_dir,
                                config_yaml=config_yaml,
                                samples=samples,
                                device=device,
                                timeout=sub_timeout,
                            )
                            sub_stderr = proc.stderr or ""
                            if "OutOfMemoryError" in sub_stderr or "CUDA out of memory" in sub_stderr:
                                oom_again = True
                                break
                        finally:
                            if sub_csv_path.exists():
                                sub_csv_path.unlink()
                    if not oom_again:
                        break
                    if retry_size == 1:
                        # GPU OOM even at size 1 — last resort: CPU
                        print(f"  ⚠ CUDA OOM even at batch size 1 — falling back to CPU...")
                        for sb_idx, sub_batch in enumerate(sub_batches):
                            sub_csv_rows = [(cn, prot_cache.get(op, pp), lig_cache.get(ol, lp))
                                            for cn, pp, lp, op, ol in sub_batch]
                            sub_csv_path = output_dir / f"_batch_{batch_idx}_cpu{sb_idx}.csv"
                            _write_batch_csv(sub_csv_rows, sub_csv_path)
                            sub_timeout = timeout * len(sub_batch) if timeout > 0 else 0
                            try:
                                proc = _run_diffdock_batch_subprocess(
                                    csv_path=sub_csv_path,
                                    out_dir=output_dir,
                                    config_yaml=config_yaml,
                                    samples=samples,
                                    device="cpu",
                                    timeout=sub_timeout,
                                )
                            finally:
                                if sub_csv_path.exists():
                                    sub_csv_path.unlink()
                        break
                    retry_size = max(1, retry_size // 2)
                batch_elapsed = time.time() - batch_start
                stderr = proc.stderr or ""
            
            for combo_name, prot_path, lig_path, orig_prot, orig_lig in batch:
                combo_dir = output_dir / combo_name
                pose_files = _collect_ranked_poses(combo_dir) if combo_dir.exists() else []
                
                r = DockingResult(
                    protein_name=get_file_stem(orig_prot),
                    ligand_name=get_file_stem(orig_lig),
                    protein_path=orig_prot,
                    ligand_path=orig_lig,
                    output_dir=combo_dir,
                    status="success" if pose_files else "failed",
                    num_poses=len(pose_files),
                    pose_files=pose_files,
                    error_message="" if pose_files else f"No poses from batch. stderr tail: {stderr[-300:]}",
                    elapsed_time=batch_elapsed / len(batch),
                )
                
                if r.status == "failed":
                    record_failure(output_dir, combo_name, r)
                    error_log[combo_name] = {
                        "error": r.error_message,
                        "timestamp": datetime.now().isoformat(),
                        "protein": r.protein_name,
                        "ligand": r.ligand_name,
                        "elapsed_time": round(r.elapsed_time, 2),
                    }
                
                batch_results.append(r)
                status_icon = "✓" if r.status == "success" else "✗"
                print(f"    {status_icon} {combo_name}: {r.num_poses} poses")
            
            combos_completed += len(batch)
            print(f"  Batch {batch_idx + 1} done in {batch_elapsed:.1f}s "
                  f"({batch_elapsed/len(batch):.1f}s/complex) "
                  f"— {combos_completed}/{total_combos} total")
        
        except subprocess.TimeoutExpired:
            print(f"  ✗ Batch {batch_idx + 1} timed out ({batch_timeout}s)")
            for combo_name, prot_path, lig_path, orig_prot, orig_lig in batch:
                r = DockingResult(
                    protein_name=get_file_stem(orig_prot),
                    ligand_name=get_file_stem(orig_lig),
                    protein_path=orig_prot,
                    ligand_path=orig_lig,
                    output_dir=output_dir / combo_name,
                    status="failed",
                    error_message=f"Batch timeout ({batch_timeout}s)",
                )
                record_failure(output_dir, combo_name, r)
                batch_results.append(r)
        
        except Exception as e:
            print(f"  ✗ Batch {batch_idx + 1} error: {e}")
            for combo_name, prot_path, lig_path, orig_prot, orig_lig in batch:
                r = DockingResult(
                    protein_name=get_file_stem(orig_prot),
                    ligand_name=get_file_stem(orig_lig),
                    protein_path=orig_prot,
                    ligand_path=orig_lig,
                    output_dir=output_dir / combo_name,
                    status="failed",
                    error_message=str(e),
                )
                record_failure(output_dir, combo_name, r)
                batch_results.append(r)
        
        finally:
            if csv_path.exists():
                csv_path.unlink()
        
        # Write log incrementally after each batch
        results.extend(batch_results)
        if log_path is not None and existing_combos is not None:
            append_log_rows(
                batch_results, log_path, existing_combos,
                lig_props_cache or {}, prot_props_cache or {},
                samples, steps, device,
            )
    
    if progress_file.exists():
        progress_file.unlink()
    
    return results


def run_diffdock_sequential(
    proteins: List[Path],
    ligands: List[Path],
    output_dir: Path,
    samples: int = NUM_SAMPLES,
    steps: int = INFERENCE_STEPS,
    device: str = DIFFDOCK_DEVICE,
    timeout: int = TIMEOUT_PER_COMPLEX,
) -> List[DockingResult]:
    """Run DiffDock for all protein-ligand combinations.
    
    Supports three modes:
      - USE_BATCH_MODE=True:  CSV batch (fastest — one DiffDock call per batch)
      - NUM_WORKERS>1:        parallel single-complex subprocesses (CPU only)
      - else:                 sequential single-complex subprocesses
    
    Molecular properties are pre-computed before docking starts and the
    docking log is written incrementally after each batch / result.
    """
    total = len(proteins) * len(ligands)
    results = []
    
    error_log = load_error_log(output_dir)
    if OVERWRITE_ERROR_LOG:
        error_log = {}
        save_error_log(output_dir, error_log)
        print("Error log cleared (OVERWRITE_ERROR_LOG=True)")
    
    num_prev_failed = len(error_log)
    
    print("=" * 80)
    print("DiffDock Docking")
    print("=" * 80)
    mode_str = "CSV batch" if USE_BATCH_MODE else (f"parallel ({NUM_WORKERS} workers)" if NUM_WORKERS > 1 else "sequential")
    print(f"Mode: {mode_str}")
    print(f"Output directory: {output_dir}")
    print(f"Proteins: {len(proteins)}  |  Ligands: {len(ligands)}  |  Total: {total}")
    print(f"Previously failed (will skip): {num_prev_failed}")
    print(f"Samples: {samples}  |  Steps: {steps}  |  Timeout: {timeout}s")
    
    if num_prev_failed > 0:
        print(f"\nPreviously failed (from {ERROR_LOG_FILENAME}):")
        for i, (name, info) in enumerate(sorted(error_log.items())):
            if i < 10:
                print(f"  ✗ {name}: {info.get('error', 'unknown')[:80]}")
        if num_prev_failed > 10:
            print(f"  ... and {num_prev_failed - 10} more")

    # ── Pre-compute molecular properties before docking ──
    lig_props_cache, prot_props_cache = precompute_properties(proteins, ligands)

    # ── Exclude ligands that RDKit cannot parse (all properties None) ──
    invalid_ligs = [lig for lig in ligands if not ligand_props_valid(lig_props_cache.get(lig, {}))]
    if invalid_ligs:
        print(f"\n  ⚠ Excluding {len(invalid_ligs)} unreadable ligand(s):")
        for lig in invalid_ligs:
            print(f"    ✗ {lig.name}")
        ligands = [lig for lig in ligands if ligand_props_valid(lig_props_cache.get(lig, {}))]
        total = len(proteins) * len(ligands)
        print(f"  Remaining: {len(ligands)} ligands × {len(proteins)} proteins = {total} combinations")

    # ── Initialise docking log (create CSV header if needed) ──
    log_path, existing_combos = init_docking_log(output_dir)
    print(f"Docking log: {log_path} ({len(existing_combos)} existing entries)")

    # ── Prepare all inputs upfront (shared cache) ──
    print("\nPreparing inputs...")
    prot_cache, lig_cache = _prepare_all_inputs(proteins, ligands, output_dir)
    
    # ── Classify combos: skip vs dock ──
    combos_to_dock = []
    skipped_results = []
    
    for protein in proteins:
        for ligand in ligands:
            prot_stem = get_file_stem(protein)
            lig_stem = get_file_stem(ligand)
            combo_name = f"{lig_stem}__{prot_stem}"
            combo_dir = output_dir / combo_name
            
            if not OVERWRITE_ERROR_LOG and combo_name in error_log:
                skipped_results.append(DockingResult(
                    protein_name=prot_stem, ligand_name=lig_stem,
                    protein_path=protein, ligand_path=ligand,
                    output_dir=combo_dir, status="skipped",
                    error_message=f"Previously failed: {error_log[combo_name].get('error', 'unknown')[:100]}",
                ))
                continue
            
            if not OVERWRITE_EXISTING and combo_dir.exists():
                existing = _collect_ranked_poses(combo_dir)
                if existing:
                    skipped_results.append(DockingResult(
                        protein_name=prot_stem, ligand_name=lig_stem,
                        protein_path=protein, ligand_path=ligand,
                        output_dir=combo_dir, status="skipped",
                        num_poses=len(existing), pose_files=existing,
                    ))
                    continue
            
            prepared_prot = prot_cache.get(protein, protein)
            prepared_lig = lig_cache.get(ligand, ligand)
            combos_to_dock.append((combo_name, prepared_prot, prepared_lig, protein, ligand))
    
    results.extend(skipped_results)
    n_skip_fail = sum(1 for r in skipped_results if r.error_message.startswith("Previously"))
    n_skip_done = len(skipped_results) - n_skip_fail
    print(f"\nSkipping {n_skip_done} already docked + {n_skip_fail} previously failed = {len(skipped_results)} total")
    print(f"Will dock: {len(combos_to_dock)} combinations")

    # Log skipped results immediately so they appear in the CSV
    if skipped_results:
        append_log_rows(
            skipped_results, log_path, existing_combos,
            lig_props_cache, prot_props_cache,
            samples, steps, device,
        )
    
    if not combos_to_dock:
        print("Nothing to dock!")
        save_error_log(output_dir, error_log)
        return results
    
    # ── Dock ──
    start_all = time.time()
    
    if USE_BATCH_MODE:
        batch_results = _run_batch_csv_mode(
            combos_to_dock=combos_to_dock,
            output_dir=output_dir,
            prot_cache=prot_cache,
            lig_cache=lig_cache,
            samples=samples,
            steps=steps,
            device=device,
            timeout=timeout,
            error_log=error_log,
            log_path=log_path,
            existing_combos=existing_combos,
            lig_props_cache=lig_props_cache,
            prot_props_cache=prot_props_cache,
        )
        results.extend(batch_results)
    
    elif NUM_WORKERS > 1:
        print(f"\nRunning {len(combos_to_dock)} complexes with {NUM_WORKERS} parallel workers...")
        
        def _dock_one(args):
            combo_name, prep_prot, prep_lig, orig_prot, orig_lig = args
            combo_dir = output_dir / combo_name
            combo_dir.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            success, poses, err = run_diffdock_single(
                prep_prot, prep_lig, combo_dir, output_dir,
                samples, steps, device, timeout,
            )
            elapsed = time.time() - t0
            r = DockingResult(
                protein_name=get_file_stem(orig_prot),
                ligand_name=get_file_stem(orig_lig),
                protein_path=orig_prot, ligand_path=orig_lig,
                output_dir=combo_dir,
                status="success" if success else "failed",
                num_poses=len(poses), pose_files=poses,
                error_message=err, elapsed_time=elapsed,
            )
            return combo_name, r
        
        with ThreadPoolExecutor(max_workers=NUM_WORKERS) as executor:
            futures = {executor.submit(_dock_one, c): c[0] for c in combos_to_dock}
            batch_buffer = []
            for i, future in enumerate(as_completed(futures), 1):
                combo_name, r = future.result()
                results.append(r)
                batch_buffer.append(r)
                if r.status == "failed":
                    record_failure(output_dir, combo_name, r)
                    error_log[combo_name] = {
                        "error": r.error_message,
                        "timestamp": datetime.now().isoformat(),
                        "protein": r.protein_name,
                        "ligand": r.ligand_name,
                        "elapsed_time": round(r.elapsed_time, 2),
                    }
                icon = "✓" if r.status == "success" else "✗"
                print(f"  [{i}/{len(combos_to_dock)}] {icon} {combo_name}: {r.num_poses} poses ({r.elapsed_time:.1f}s)")
                # Flush to log every BATCH_SIZE results (or at the end)
                if len(batch_buffer) >= BATCH_SIZE or i == len(combos_to_dock):
                    append_log_rows(
                        batch_buffer, log_path, existing_combos,
                        lig_props_cache, prot_props_cache,
                        samples, steps, device,
                    )
                    batch_buffer = []
    
    else:
        for i, (combo_name, _, _, orig_prot, orig_lig) in enumerate(combos_to_dock, 1):
            print(f"\n[{i}/{len(combos_to_dock)}] {get_file_stem(orig_lig)} + {get_file_stem(orig_prot)}")
            
            r = dock_single_combination(
                protein_path=orig_prot,
                ligand_path=orig_lig,
                output_base_dir=output_dir,
                prot_cache=prot_cache,
                lig_cache=lig_cache,
                samples=samples, steps=steps,
                device=device, timeout=timeout,
                error_log=error_log,
            )
            
            if r.status == "failed":
                error_log[combo_name] = {
                    "error": r.error_message,
                    "timestamp": datetime.now().isoformat(),
                    "protein": r.protein_name,
                    "ligand": r.ligand_name,
                    "elapsed_time": round(r.elapsed_time, 2),
                }
            
            results.append(r)
            icon = {"success": "✓", "failed": "✗", "skipped": "⊘"}.get(r.status, "?")
            print(f"  {icon} {r.status} | Poses: {r.num_poses} | Time: {r.elapsed_time:.1f}s")
            if r.error_message:
                print(f"    Error: {r.error_message[:200]}")
            
            # Write to log after each individual docking
            append_log_rows(
                [r], log_path, existing_combos,
                lig_props_cache, prot_props_cache,
                samples, steps, device,
            )
    
    total_time = time.time() - start_all
    n_success = sum(1 for r in results if r.status == "success")
    n_failed = sum(1 for r in results if r.status == "failed")

    print(f"\nDocking complete in {total_time:.1f}s ({total_time/60:.1f} min)")
    print("Batch docking functions defined.")

    print(f"  Success: {n_success}  |  Failed: {n_failed}  |  Skipped: {len(skipped_results)}")

    

    save_error_log(output_dir, error_log)    
    return results

    print(f"Error log saved: {output_dir / ERROR_LOG_FILENAME} ({len(error_log)} entries)")
    print(f"Docking log: {log_path} ({len(existing_combos)} entries)")

# %% [markdown]
# # Pre-flight Check

# %%
errors = []

# Check DiffDock directory + config
if not DIFFDOCK_DIR.exists():
    errors.append(f"✗ DiffDock dir missing: {DIFFDOCK_DIR}")
default_yaml = DIFFDOCK_DIR / "default_inference_args.yaml"
if not default_yaml.exists():
    errors.append(f"✗ DiffDock config missing: {default_yaml}")

# Check Python executable
if not Path(DIFFDOCK_PYTHON).exists():
    errors.append(f"✗ DiffDock Python not found: {DIFFDOCK_PYTHON}")

# Check that key functions are the correct (latest) versions
src = inspect.getsource(_run_diffdock_subprocess)
if "subprocess.run" not in src:
    errors.append("✗ _run_diffdock_subprocess is incomplete (missing subprocess.run) — re-run the correct cell")
if "sys.executable" in src:
    errors.append("⚠ _run_diffdock_subprocess uses sys.executable instead of DIFFDOCK_PYTHON")

# Check inputs exist
for drugs_dir in drugs_dirs:
    proteins = collect_files(receptors_dir, [".pdb"])
    ligands = collect_files(drugs_dir, [".sdf", ".mol2"])
    print(f"  {drugs_dir.name}: {len(proteins)} proteins, {len(ligands)} ligands")
    if not proteins:
        errors.append(f"✗ No .pdb files in {receptors_dir}")
    if not ligands:
        errors.append(f"✗ No .sdf/.mol2 files in {drugs_dir}")

if errors:
    print("\n⚠ ISSUES FOUND:")
    for e in errors:
        print(f"  {e}")
    print("\nFix these before running docking.")
else:
    print("\n✓ All checks passed — ready to dock.")

# %%


# %%
def _monitor_progress(output_dir: Path, stop_event: threading.Event, poll_interval: float = 10.0):
    """Background thread that periodically prints docking progress."""
    log_dir = output_dir / "_batch_logs"
    progress_file = log_dir / "_progress.txt"
    start_time = time.time()
    
    while not stop_event.is_set():
        stop_event.wait(poll_interval)
        if stop_event.is_set():
            break
        
        elapsed = time.time() - start_time
        
        # Count completed poses (directories with rank1 SDF files)
        rank1_files = list(output_dir.glob("*/rank1_*.sdf"))
        n_done = len(rank1_files)
        
        # Count empty dirs (in-progress or failed)
        all_combo_dirs = [d for d in output_dir.iterdir() 
                         if d.is_dir() and not d.name.startswith("_") 
                         and d.name not in ("prepared_proteins", "prepared_ligands")]
        n_empty = sum(1 for d in all_combo_dirs if not list(d.glob("rank1_*.sdf")))
        
        # Read progress file written by batch mode
        current_batch_info = ""
        current_combos = ""
        if progress_file.exists():
            try:
                parts = progress_file.read_text().strip().split("|")
                if len(parts) >= 2:
                    current_batch_info = f"[{parts[0]} combos] {parts[1]}"
                if len(parts) >= 3:
                    current_combos = ", ".join(parts[2:6])  # show up to 4 current combos
                    if len(parts) > 6:
                        current_combos += f" (+{len(parts) - 6} more)"
            except Exception:
                pass
        
        # GPU utilization
        try:
            gpu_info = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            gpu_line = gpu_info.stdout.strip()
        except Exception:
            gpu_line = "N/A"
        
        # Read last line of stdout log for DiffDock's own progress
        diffdock_status = ""
        log_files = sorted(log_dir.glob("batch_*_stdout.log")) if log_dir.exists() else []
        if log_files:
            try:
                with open(log_files[-1], 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        stripped = line.strip()
                        if stripped:
                            diffdock_status = stripped[:120]
                            break
            except Exception:
                pass
        
        # stderr last meaningful line (skip noise)
        STDERR_SKIP = {"FutureWarning", "DeprecationWarning", "KeyboardInterrupt", "UserWarning"}
        stderr_status = ""
        stderr_files = sorted(log_dir.glob("batch_*_stderr.log")) if log_dir.exists() else []
        if stderr_files:
            try:
                with open(stderr_files[-1], 'r') as f:
                    lines = f.readlines()
                    for line in reversed(lines):
                        stripped = line.strip()
                        if stripped and not any(skip in stripped for skip in STDERR_SKIP):
                            stderr_status = stripped[:120]
                            break
            except Exception:
                pass
        
        mins, secs = divmod(int(elapsed), 60)
        hrs, mins = divmod(mins, 60)
        time_str = f"{hrs}h{mins:02d}m{secs:02d}s" if hrs else f"{mins}m{secs:02d}s"
        
        line1 = f"  ⏱ {time_str} | ✓ {n_done} done | ⏳ {n_empty} pending | GPU: {gpu_line}"
        if current_batch_info:
            line1 += f" | {current_batch_info}"
        print(line1)
        if current_combos:
            print(f"    Docking: {current_combos}")
        if diffdock_status:
            print(f"    DiffDock: {diffdock_status}")
        if stderr_status:
            print(f"    Status:   {stderr_status}")


all_diffdock_results = []

for drugs_dir in drugs_dirs:
    DIFFDOCK_OUTPUT_DIR = DIFFDOCK_OUTPUT_BASE / drugs_dir.name
    DIFFDOCK_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    print(f"\n{'=' * 80}")
    print(f"Processing ligand directory: {drugs_dir.name}")
    print(f"Output: {DIFFDOCK_OUTPUT_DIR}")
    print(f"{'=' * 80}")
    
    # Collect input files
    proteins = collect_files(receptors_dir, [".pdb"])
    ligands = collect_files(drugs_dir, [".sdf", ".mol2"])
    
    print(f"Found {len(proteins)} proteins and {len(ligands)} ligands")
    
    if not proteins or not ligands:
        print("⚠ No proteins or ligands found, skipping.")
        continue
    
    # Start progress monitor
    stop_monitor = threading.Event()
    monitor_thread = threading.Thread(
        target=_monitor_progress,
        args=(DIFFDOCK_OUTPUT_DIR, stop_monitor, 15),
        daemon=True,
    )
    monitor_thread.start()
    
    try:
        # Run docking
        results = run_diffdock_sequential(
            proteins=proteins,
            ligands=ligands,
            output_dir=DIFFDOCK_OUTPUT_DIR,
        )
    finally:
        stop_monitor.set()
        monitor_thread.join(timeout=5)
    
    all_diffdock_results.extend(results)
    
    # Generate and print summary
    summary = generate_summary(results)
    print_summary(summary)
    
    # Save summary JSON
    summary_path = DIFFDOCK_OUTPUT_DIR / "docking_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSummary saved: {summary_path}")
    
    # Docking log was written incrementally during docking — no final write needed

diffdock_results = all_diffdock_results
print(f"\n{'=' * 80}")
print(f"ALL DONE — {len(diffdock_results)} total results across {len(drugs_dirs)} ligand directories")
print(f"{'=' * 80}")

# %% [markdown]
# # Results

# %%
# Convert results to DataFrame
results_data = []
for r in diffdock_results: 
    for i, pose_file in enumerate(r.pose_files):
        confidence = extract_confidence_from_filename(pose_file.name)
        rank = extract_rank_from_filename(pose_file.name)
        results_data.append({
            "protein": r.protein_name,
            "ligand": r.ligand_name,
            "rank": rank if rank > 0 else i + 1,
            "confidence": confidence,
            "sdf_path": str(pose_file),
            "status": r.status,
        })

results_df = pd.DataFrame(results_data)

print("=" * 80)
print("DOCKING RESULTS DATAFRAME")
print("=" * 80)
print(f"\nTotal poses: {len(results_df)}")

if not results_df.empty:
    print("\nPoses per protein-ligand combination:")
    pose_counts = results_df.groupby(["protein", "ligand"]).agg({
        "rank": "count",
        "confidence": "mean",
    }).reset_index()
    pose_counts.columns = ["protein", "ligand", "num_poses", "avg_confidence"]
    display(pose_counts)
    
    print("\nConfidence score distribution:")
    if results_df['confidence'].sum() > 0:
        print(f"  Mean confidence: {results_df['confidence'].mean():.3f}")
        print(f"  Max confidence: {results_df['confidence'].max():.3f}")
        print(f"  Min confidence: {results_df['confidence'].min():.3f}")

# Save to CSV
csv_path = DIFFDOCK_OUTPUT_DIR / "diffdock_poses.csv"
results_df.to_csv(csv_path, index=False)
print(f"\nResults saved to: {csv_path}")

# %%

# ============================================================================
# LIST ALL GENERATED POSES
# ============================================================================

def list_pose_files(output_dir: Path) -> Dict[str, List[Path]]:
    """List all generated pose files organized by combination."""
    poses_by_combo = {}
    
    if not output_dir.exists():
        return poses_by_combo
    
    for combo_dir in sorted(output_dir.iterdir()):
        if not combo_dir.is_dir():
            continue
        
        # Only count actual ranked poses (rank*_confidence*.sdf)
        pose_files = sorted(
            _collect_ranked_poses(combo_dir),
            key=lambda p: extract_rank_from_filename(p.name),
        )
        if pose_files:
            poses_by_combo[combo_dir.name] = pose_files
    
    return poses_by_combo


pose_files = list_pose_files(DIFFDOCK_OUTPUT_DIR)

print("=" * 80)
print("GENERATED POSE FILES")
print("=" * 80)
print(f"\nOutput directory: {DIFFDOCK_OUTPUT_DIR}")
print(f"Total combinations with poses: {len(pose_files)}")
print()

total_poses = 0
for combo_name, files in pose_files.items():
    print(f"\n{combo_name}/")
    for f in files[:5]:  # Show first 5 poses per combination
        conf = extract_confidence_from_filename(f.name)
        print(f"  └── {f.name} (confidence: {conf:.3f})")
    if len(files) > 5:
        print(f"  └── ... and {len(files) - 5} more")
    total_poses += len(files)

print("\n" + "-" * 80)
print(f"TOTAL POSES GENERATED: {total_poses}")
print("-" * 80)


# %% [markdown]
# # Docking Log

# %%
log_path = DIFFDOCK_OUTPUT_DIR / DOCKING_LOG_FILENAME

if log_path.exists():
    docking_log_df = pd.read_csv(log_path)
    
    print("=" * 80)
    print("DOCKING LOG")
    print("=" * 80)
    print(f"\nLog file: {log_path}")
    print(f"Total entries: {len(docking_log_df)}")
    
    # Status breakdown
    status_counts = docking_log_df["status"].value_counts()
    print(f"\nStatus breakdown:")
    for status, count in status_counts.items():
        print(f"  {status}: {count}")
    
    # Show failed entries with reasons
    failed = docking_log_df[docking_log_df["status"] == "failed"]
    if not failed.empty:
        print(f"\n{'─' * 80}")
        print(f"FAILED DOCKINGS ({len(failed)}):")
        print(f"{'─' * 80}")
        for _, row in failed.iterrows():
            reason = row["error_reason"][:120] if pd.notna(row["error_reason"]) else "unknown"
            print(f"  ✗ {row['combo_name']}: {reason}")
    
    # Ligand property summary
    lig_cols = ["ligand_name", "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
                "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
                "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula"]
    lig_summary = docking_log_df[lig_cols].drop_duplicates(subset=["ligand_name"]).sort_values("ligand_name")
    print(f"\n{'─' * 80}")
    print("LIGAND PROPERTIES:")
    print(f"{'─' * 80}")
    display(lig_summary.reset_index(drop=True))
    
    # Protein property summary
    prot_cols = ["protein_name", "prot_num_residues", "prot_num_atoms", "prot_num_chains"]
    prot_summary = docking_log_df[prot_cols].drop_duplicates(subset=["protein_name"]).sort_values("protein_name")
    print(f"\n{'─' * 80}")
    print("PROTEIN PROPERTIES:")
    print(f"{'─' * 80}")
    display(prot_summary.reset_index(drop=True))
    
    # Timing summary
    docked = docking_log_df[docking_log_df["status"].isin(["success", "failed"])]
    if not docked.empty:
        print(f"\n{'─' * 80}")
        print("TIMING:")
        print(f"{'─' * 80}")
        print(f"  Mean time per complex: {docked['elapsed_time_s'].mean():.1f}s")
        print(f"  Max time: {docked['elapsed_time_s'].max():.1f}s")
        print(f"  Total: {docked['elapsed_time_s'].sum():.1f}s ({docked['elapsed_time_s'].sum()/60:.1f} min)")
else:
    print(f"No docking log found at {log_path}. Run docking first.")


