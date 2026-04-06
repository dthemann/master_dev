# %%
from __future__ import annotations

import csv
import itertools
import json
import os
import re
import subprocess
import threading
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from IPython.display import display
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

# %%
# # Verify we are running inside the vina-gpu conda environment
# import sys, subprocess
# 
# env_name = sys.executable
# print(f"Python executable: {env_name}")
# assert "vina-gpu" in env_name, (
#     f"This notebook requires the 'vina-gpu' conda environment.\n"
#     f"Select it via: Kernel → Change Kernel → vina-gpu\n"
#     f"Current interpreter: {env_name}"
# )
# 
# # Verify the Vina-GPU binary is available
# vina_gpu_bin = "/workspace/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1/AutoDock-Vina-GPU-2-1"
# result = subprocess.run([vina_gpu_bin, "--version"], capture_output=True, text=True, timeout=10)
# print(f"Vina-GPU binary: {vina_gpu_bin}")
# print(f"Version output: {result.stdout.strip() or result.stderr.strip()}")
# 
# # Verify GPU is available
# gpu_check = subprocess.run(
#     ["nvidia-smi", "--query-gpu=gpu_name,memory.total", "--format=csv,noheader"],
#     capture_output=True, text=True, timeout=10,
# )
# print(f"GPU: {gpu_check.stdout.strip()}")
# print("\n✓ vina-gpu environment verified")

# %%
# Batch clean directory
from Scripts.clean_pdb_like_pymol import batch_clean_pdbs
batch_clean_pdbs('Orai/', pattern='*.pdb', fix_histidines=True)

# %% [markdown]
# # Config

# %%
from Scripts.prep_docking import run_workflow

# ── Paths ──
workspace_root = Path("/workspace")
wd = workspace_root / "master_dev"
wdlogs = wd / "logs"

orai_folder = wd / "Orai"
receptors_dir = orai_folder  # alias for consistency with DiffDock naming

# List of ligand directories — each gets its own output folder, error log, and summary
drugs_dirs: List[Path] = [
    Path("/workspace/storage/ligands_sdf_large_approved"),
    Path("/workspace/storage/ligands_sdf_small_symmetric_approved"),  # not present on this pod
]

# ── Vina-GPU executable ──
VINA_GPU_DIR: Path = Path("/workspace/Vina-GPU-2.1/AutoDock-Vina-GPU-2.1")
VINA_BIN: Path = VINA_GPU_DIR / "AutoDock-Vina-GPU-2-1"
OPENCL_BINARY_PATH: Path = VINA_GPU_DIR  # directory with Kernel1_Opt.bin / Kernel2_Opt.bin

# ── Vina-GPU docking parameters ──
THREAD: int = 8000           # number of GPU computing lanes (higher = more GPU usage)
SEARCH_DEPTH: int = 32       # Monte-Carlo search depth (replaces CPU exhaustiveness)
NUM_MODES: int = 10          # Number of binding modes to output per ligand-receptor pair
ENERGY_RANGE: int = 3
SEED: Optional[int] = None
TIMEOUT_PER_COMPLEX: int = 600  # seconds; set to 0 for no limit

# ── Output ──
VINA_OUTPUT_BASE = workspace_root / "vina_results"
VINA_OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

# ── Overwrite settings ──
OVERWRITE_EXISTING: bool = False  # Set to True to re-dock existing combinations
OVERWRITE_ERROR_LOG: bool = True  # Set to True to retry previously failed ligands

# ── Log filenames ──
ERROR_LOG_FILENAME = "failed_docking.json"
DOCKING_LOG_FILENAME = "docking_log.csv"

# ── Preparation folders ──
output_folder_mgltools = Path("docking_ready_mgltools")
output_folder_obabel = Path("docking_read_obabel")
output_folder_meeko = Path("docking_ready_meeko")

output_folder_mgltools.mkdir(parents=True, exist_ok=True)
output_folder_meeko.mkdir(parents=True, exist_ok=True)
output_folder_obabel.mkdir(parents=True, exist_ok=True)
wdlogs.mkdir(parents=True, exist_ok=True)

# ── Validate directories ──
for drugs_dir in drugs_dirs:
    if not drugs_dir.exists():
        raise FileNotFoundError(f"Ligand directory missing: {drugs_dir}")
if not receptors_dir.exists():
    raise FileNotFoundError(f"Receptor directory missing: {receptors_dir}")

print("=" * 80)
print("AutoDock Vina-GPU Batch Docking Configuration")
print("=" * 80)
print(f"Vina-GPU executable: {VINA_BIN}")
print(f"OpenCL binary path: {OPENCL_BINARY_PATH}")
print(f"GPU threads: {THREAD}")
print(f"Search depth: {SEARCH_DEPTH}")
print(f"Num modes: {NUM_MODES}")
print(f"Energy range: {ENERGY_RANGE} kcal/mol")
print(f"Timeout per complex: {TIMEOUT_PER_COMPLEX}s ({TIMEOUT_PER_COMPLEX/60:.0f} min)" if TIMEOUT_PER_COMPLEX else "Timeout per complex: unlimited")
print(f"Output base directory: {VINA_OUTPUT_BASE}")
print(f"Overwrite existing: {OVERWRITE_EXISTING}")
print(f"\nLigand directories ({len(drugs_dirs)}):")
for d in drugs_dirs:
    n_files = len(list(d.glob("*.pdbqt")))
    output_dir = VINA_OUTPUT_BASE / d.name
    print(f"  • {d.name:30s} → output: {output_dir.name}/  ({n_files} ligand files)")
print()

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
    best_affinity: Optional[float] = None
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
            "best_affinity": self.best_affinity,
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


def parse_vina_affinities(log_text: str) -> List[dict]:
    """Parse Vina stdout/log text to extract affinity scores for each mode."""
    modes: List[dict] = []
    mode_pattern = r'^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
    for m in re.finditer(mode_pattern, log_text, re.MULTILINE):
        modes.append({
            'mode': int(m.group(1)),
            'affinity': float(m.group(2)),
            'rmsd_lb': float(m.group(3)),
            'rmsd_ub': float(m.group(4)),
        })
    return modes


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


def record_failure(output_dir: Path, combo_name: str, result: DockingResult):
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
    "search_depth", "thread", "num_modes", "energy_range", "gpu_model",
    # Best pose affinity
    "best_affinity_kcal",
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


def get_ligand_properties_from_sdf(sdf_path: Path) -> Dict:
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


def get_ligand_properties_from_pdbqt(pdbqt_path: Path) -> Dict:
    """Extract basic ligand properties by parsing a PDBQT file (fallback when no SDF)."""
    props = {
        "lig_molecular_weight": None, "lig_heavy_atoms": None,
        "lig_total_atoms": None, "lig_rotatable_bonds": None,
        "lig_num_rings": None, "lig_aromatic_rings": None,
        "lig_hbd": None, "lig_hba": None, "lig_tpsa": None,
        "lig_logp": None, "lig_formula": None,
    }
    try:
        atom_count = 0
        heavy_count = 0
        rotatable = 0
        with open(pdbqt_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_count += 1
                    element = line[77:79].strip() if len(line) > 77 else ""
                    if element and element not in ("H", "HD"):
                        heavy_count += 1
                elif line.startswith("BRANCH"):
                    rotatable += 1
        props["lig_total_atoms"] = atom_count
        props["lig_heavy_atoms"] = heavy_count if heavy_count > 0 else None
        props["lig_rotatable_bonds"] = rotatable // 2 if rotatable > 0 else None
    except Exception as e:
        print(f"  ⚠ Could not parse PDBQT ligand {pdbqt_path.name}: {e}")
    return props


def get_ligand_properties(lig_path: Path) -> Dict:
    """Get ligand properties — prefer SDF, fall back to PDBQT parsing."""
    if lig_path.suffix.lower() in (".sdf", ".mol2"):
        return get_ligand_properties_from_sdf(lig_path)
    return get_ligand_properties_from_pdbqt(lig_path)


def get_protein_properties(pdb_path: Path) -> Dict:
    """Extract basic protein properties by parsing PDB/PDBQT ATOM records."""
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
    return props.get("lig_heavy_atoms") is not None or props.get("lig_total_atoms") is not None


def precompute_properties(
    proteins: List[Path],
    ligands: List[Path],
) -> Tuple[Dict[Path, Dict], Dict[Path, Dict]]:
    """Pre-compute molecular properties for all proteins and ligands before docking."""
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


_GPU_MODEL = get_gpu_model()
print(f"GPU detected: {_GPU_MODEL}")


def _build_log_row(
    result: DockingResult,
    lig_props: Dict,
    prot_props: Dict,
) -> Dict:
    """Build a single row dict for the docking log CSV."""
    row = {
        "timestamp": datetime.now().isoformat(),
        "combo_name": f"{result.ligand_name}__{result.protein_name}",
        "protein_name": result.protein_name,
        "ligand_name": result.ligand_name,
        "status": result.status,
        "error_reason": result.error_message if result.status == "failed" else "",
        "num_poses": result.num_poses,
        "elapsed_time_s": round(result.elapsed_time, 2),
        "best_affinity_kcal": result.best_affinity,
        "search_depth": SEARCH_DEPTH,
        "thread": THREAD,
        "num_modes": NUM_MODES,
        "energy_range": ENERGY_RANGE,
        "gpu_model": _GPU_MODEL,
    }
    row.update(lig_props)
    row.update(prot_props)
    return row


def init_docking_log(output_dir: Path) -> Tuple[Path, set]:
    """Initialise the docking log CSV and return (log_path, existing_combo_names)."""
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
    results: List[DockingResult],
    log_path: Path,
    existing_combos: set,
    lig_props_cache: Dict[Path, Dict],
    prot_props_cache: Dict[Path, Dict],
) -> None:
    """Append log rows for a batch of results to the CSV, skipping duplicates."""
    rows = []
    for r in results:
        combo = f"{r.ligand_name}__{r.protein_name}"
        if combo in existing_combos:
            continue
        lig_props = lig_props_cache.get(r.ligand_path, get_ligand_properties(r.ligand_path))
        prot_props = prot_props_cache.get(r.protein_path, get_protein_properties(r.protein_path))
        rows.append(_build_log_row(r, lig_props, prot_props))
        existing_combos.add(combo)

    if not rows:
        return

    df = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
    df.to_csv(log_path, mode='a', header=False, index=False)
    print(f"  📝 Appended {len(rows)} entries to {log_path.name} ({len(existing_combos)} total)")


print("Docking log functions defined.")

# %% [markdown]
# # Summary & Reporting

# %%
def generate_summary(results: List[DockingResult]) -> Dict:
    """Generate summary of docking results."""
    summary = {
        "timestamp": datetime.now().isoformat(),
        "configuration": {
            "search_depth": SEARCH_DEPTH,
            "thread": THREAD,
            "num_modes": NUM_MODES,
            "energy_range": ENERGY_RANGE,
            "vina_bin": str(VINA_BIN),
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
            "best_affinity": r.best_affinity,
            "error": r.error_message if r.error_message else None,
        })

    summary["by_protein"] = dict(summary["by_protein"])
    summary["by_ligand"] = dict(summary["by_ligand"])
    return summary


def print_summary(summary: Dict):
    """Print formatted docking summary."""
    print("\n" + "=" * 80)
    print("AUTODOCK VINA-GPU DOCKING SUMMARY")
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
# # Autodock Vina-GPU Function

# %%
def _normalize_key(path: Path) -> str:
    """Extract a clean stem from a .pdbqt or .box.txt filename."""
    name = path.name
    if name.endswith('.box.txt'):
        return name[:-len('.box.txt')]
    if name.endswith('.pdbqt'):
        return name[:-len('.pdbqt')]
    return path.stem


def _collect_receptors(
    base_dir: Path,
    prepared_manifest: dict,
    workflow_data: dict | None,
) -> list[dict]:
    """Gather receptor PDBQT + box.txt pairs from manifest, workflow, or directory scan."""
    groups: defaultdict[str, dict] = defaultdict(dict)

    for entry in prepared_manifest.get('proteins', []):
        entry_path = Path(entry.get('path', '')).resolve()
        if not entry_path.exists():
            continue
        key = _normalize_key(entry_path)
        if entry_path.suffix.lower() == '.pdbqt':
            groups[key]['pdbqt'] = entry_path
        elif entry_path.name.lower().endswith('.box.txt'):
            groups[key]['box'] = entry_path

    if workflow_data:
        meeko_outputs = workflow_data.get('meeko_outputs') or {}
        meeko_box_files = workflow_data.get('meeko_box_files') or {}
        for original, pdbqt in meeko_outputs.items():
            box = meeko_box_files.get(original)
            if not pdbqt or not box:
                continue
            pdbqt_path = Path(pdbqt).resolve()
            box_path = Path(box).resolve()
            if not (pdbqt_path.exists() and box_path.exists()):
                continue
            key = _normalize_key(pdbqt_path)
            groups[key].setdefault('pdbqt', pdbqt_path)
            groups[key].setdefault('box', box_path)

    receptors = []
    for key, data in groups.items():
        pdbqt = data.get('pdbqt')
        box = data.get('box')
        if pdbqt and box:
            receptors.append({'name': key, 'pdbqt': pdbqt, 'box': box})

    if not receptors:
        for pdbqt_path in sorted(base_dir.glob('*.pdbqt')):
            box_path = pdbqt_path.with_suffix('.box.txt')
            if box_path.exists():
                receptors.append({
                    'name': _normalize_key(pdbqt_path),
                    'pdbqt': pdbqt_path.resolve(),
                    'box': box_path.resolve(),
                })

    receptors.sort(key=lambda r: r['name'])
    return receptors


def _collect_ligands(
    base_dir: Path,
    prepared_manifest: dict,
    workflow_data: dict | None,
) -> list[dict]:
    """Gather ligand PDBQT files from manifest, workflow, or directory scan."""
    ligands: dict[str, Path] = {}

    for entry in prepared_manifest.get('ligands', []):
        entry_path = Path(entry.get('path', '')).resolve()
        if entry_path.exists() and entry_path.suffix.lower() == '.pdbqt':
            ligands.setdefault(_normalize_key(entry_path), entry_path)

    if workflow_data:
        meeko_ligs = workflow_data.get('meeko_ligand_outputs') or {}
        for original, pdbqt in meeko_ligs.items():
            if not pdbqt:
                continue
            pdbqt_path = Path(pdbqt).resolve()
            if pdbqt_path.exists():
                ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path)

    ligand_dir = base_dir / 'ligands'
    if not ligands and ligand_dir.exists():
        for pdbqt_path in sorted(ligand_dir.glob('*.pdbqt')):
            ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path.resolve())

    if not ligands:
        for pdbqt_path in sorted(base_dir.glob('*.pdbqt')):
            if pdbqt_path.with_suffix('.box.txt').exists():
                continue
            ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path.resolve())

    return [{'name': n, 'pdbqt': p} for n, p in sorted(ligands.items())]


def run_autodock_vina_gpu(
    base_dir: Path | str,
    log_dir: Path | str,
    prepared_manifest: dict,
    protein_workflow_data: dict | None = None,
    ligand_workflow_data: dict | None = None,
    *,
    vina_bin: Path | str = VINA_BIN,
    opencl_binary_path: Path | str = OPENCL_BINARY_PATH,
    dock_subdir: str = 'docking',
    thread: int = 8000,
    search_depth: int = 32,
    num_modes: int | None = 10,
    energy_range: int | None = 3,
    seed: int | None = None,
    overwrite: bool = False,
    timeout: int = 0,
    output_base_dir: Path | None = None,
    lig_props_cache: Dict[Path, Dict] | None = None,
    prot_props_cache: Dict[Path, Dict] | None = None,
) -> Tuple[pd.DataFrame, List[DockingResult]]:
    """Run AutoDock Vina-GPU for all receptor-ligand combinations.

    Parameters
    ----------
    base_dir : Path or str
        Root directory containing prepared receptor/ligand PDBQT files.
    log_dir : Path or str
        Directory where Vina log files are written.
    prepared_manifest : dict
        Manifest with 'proteins' and 'ligands' entries (paths + metadata).
    protein_workflow_data : dict, optional
        Workflow output dict for receptors (e.g. Meeko outputs).
    ligand_workflow_data : dict, optional
        Workflow output dict for ligands (e.g. Meeko ligand outputs).
    vina_bin : Path or str
        Path to the Vina-GPU executable.
    opencl_binary_path : Path or str
        Path to directory containing precompiled OpenCL kernels.
    dock_subdir : str
        Subdirectory under *base_dir* for docking output files.
    thread : int
        Number of GPU computing lanes.
    search_depth : int
        Monte-Carlo search depth (replaces CPU exhaustiveness).
    num_modes : int or None
        Maximum number of binding modes to generate.
    energy_range : int or None
        Maximum energy difference between best and worst mode (kcal/mol).
    seed : int or None
        Random seed for reproducibility.
    overwrite : bool
        Re-run jobs whose output already exists.
    timeout : int
        Per-complex timeout in seconds (0 = no limit).
    output_base_dir : Path, optional
        If provided, docking log CSV and error log are written here.
    lig_props_cache : dict, optional
        Pre-computed ligand properties keyed by Path.
    prot_props_cache : dict, optional
        Pre-computed protein properties keyed by Path.

    Returns
    -------
    (pd.DataFrame, List[DockingResult])
        Results table and list of DockingResult objects.
    """
    base_dir = Path(base_dir).resolve()
    log_dir = Path(log_dir).resolve()
    vina_bin = Path(vina_bin).expanduser().resolve()
    opencl_binary_path = Path(opencl_binary_path).expanduser().resolve()
    dock_dir = base_dir / dock_subdir

    dock_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    pd.set_option('display.max_colwidth', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 0)

    # ── Collect inputs ──────────────────────────────────────────────────
    receptors = _collect_receptors(base_dir, prepared_manifest, protein_workflow_data)
    ligands = _collect_ligands(base_dir, prepared_manifest, ligand_workflow_data)
    print(f'Docking receptors: {len(receptors)} | ligands: {len(ligands)}')

    vina_ok = vina_bin.is_file() and os.access(vina_bin, os.X_OK)
    if not vina_ok:
        print(f"WARNING: Vina-GPU executable not found or not executable ('{vina_bin}').")
    if not receptors:
        print('WARNING: No receptor PDBQT + box.txt pairs found; aborting docking step.')
    if not ligands:
        print('WARNING: No ligand PDBQT files found; aborting docking step.')

    raw_results: list[dict] = []
    docking_results: List[DockingResult] = []
    if not vina_ok or not receptors or not ligands:
        return pd.DataFrame(raw_results), docking_results

    # ── Error log ───────────────────────────────────────────────────────
    effective_output_dir = output_base_dir or dock_dir
    error_log = load_error_log(effective_output_dir)
    if OVERWRITE_ERROR_LOG:
        error_log = {}
        save_error_log(effective_output_dir, error_log)

    # ── Docking log ─────────────────────────────────────────────────────
    csv_log_path, existing_combos = init_docking_log(effective_output_dir)
    print(f"Docking log: {csv_log_path} ({len(existing_combos)} existing entries)")

    # ── Build & run ─────────────────────────────────────────────────────
    def _build_cmd(receptor_pdbqt, ligand_pdbqt, box_file, out_path):
        cmd = [
            str(vina_bin),
            '--receptor', str(receptor_pdbqt),
            '--ligand', str(ligand_pdbqt),
            '--config', str(box_file),
            '--thread', str(thread),
            '--search_depth', str(search_depth),
            '--opencl_binary_path', str(opencl_binary_path),
            '--out', str(out_path),
        ]
        if num_modes is not None:
            cmd += ['--num_modes', str(num_modes)]
        if energy_range is not None:
            cmd += ['--energy_range', str(energy_range)]
        if seed is not None:
            cmd += ['--seed', str(seed)]
        return cmd

    def _run(cmd, log_path, complex_timeout=0):
        effective_timeout = complex_timeout if complex_timeout > 0 else None
        try:
            with open(log_path, 'w') as log_f:
                proc = subprocess.Popen(
                    cmd, stdout=log_f, stderr=subprocess.STDOUT,
                )
                try:
                    retcode = proc.wait(timeout=effective_timeout)
                except subprocess.TimeoutExpired:
                    proc.kill()
                    proc.wait()
                    return False, f"Timeout ({complex_timeout}s)"
            output = log_path.read_text()
            if retcode != 0:
                return False, output.strip()
            return True, output
        except Exception as exc:
            return False, str(exc)

    jobs = list(itertools.product(receptors, ligands))
    print(f"Total docking jobs: {len(jobs)}")
    start_all = time.time()

    for i, (rec, lig) in enumerate(jobs, 1):
        rec_path = rec['pdbqt']
        lig_path = lig['pdbqt']
        box_path = rec['box']
        rec_name = _normalize_key(rec_path)
        lig_name = _normalize_key(lig_path)
        combo_name = f"{lig_name}__{rec_name}"
        stem = f"{rec_path.stem}__{lig_path.stem}_vina_out"
        out_path = dock_dir / f"{stem}.pdbqt"
        vina_log_path = log_dir / f"{stem}.log"

        # Skip previously failed
        if not OVERWRITE_ERROR_LOG and combo_name in error_log:
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="skipped",
                error_message=f"Previously failed: {error_log[combo_name].get('error', 'unknown')[:100]}",
            )
            docking_results.append(r)
            raw_results.append({
                'receptor': str(rec_path), 'ligand': str(lig_path),
                'output': str(out_path), 'log': str(vina_log_path),
                'status': 'skip-failed',
            })
            continue

        # Skip existing
        if out_path.exists() and not overwrite:
            # Parse existing log to get pose count + affinity
            n_poses = 0
            best_aff = None
            if vina_log_path.exists():
                modes = parse_vina_affinities(vina_log_path.read_text())
                n_poses = len(modes)
                if modes:
                    best_aff = modes[0]['affinity']

            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="skipped",
                num_poses=n_poses, best_affinity=best_aff,
                pose_files=[out_path],
            )
            docking_results.append(r)
            raw_results.append({
                'receptor': str(rec_path), 'ligand': str(lig_path),
                'output': str(out_path), 'log': str(vina_log_path),
                'status': 'skip-exists',
            })
            continue

        # Build and execute command
        cmd = _build_cmd(rec_path, lig_path, box_path, out_path)
        t0 = time.time()
        ok, output = _run(cmd, vina_log_path, timeout)
        elapsed = time.time() - t0

        status_str = 'ok' if ok and out_path.exists() and out_path.stat().st_size > 0 else 'fail'
        if status_str != 'ok' and out_path.exists():
            try:
                out_path.unlink()
            except Exception:
                pass

        # Parse affinities from output
        modes = parse_vina_affinities(output or '')
        n_poses = len(modes)
        best_aff = modes[0]['affinity'] if modes else None
        status_str = 'ok' if ok and out_path.exists() and out_path.stat().st_size > 0 else 'fail'
        r = DockingResult(
            protein_name=rec_name, ligand_name=lig_name,
            protein_path=rec_path, ligand_path=lig_path,
            output_dir=dock_dir,
            status="success" if status_str == 'ok' else "failed",
            num_poses=n_poses,
            best_affinity=best_aff,
            pose_files=[out_path] if status_str == 'ok' else [],
            error_message="" if status_str == 'ok' else (output or "Unknown error")[:300],
            elapsed_time=elapsed,
        )

        if r.status == "failed":
            record_failure(effective_output_dir, combo_name, r)
            error_log[combo_name] = {
                "error": r.error_message,
                "timestamp": datetime.now().isoformat(),
                "protein": r.protein_name,
                "ligand": r.ligand_name,
                "elapsed_time": round(r.elapsed_time, 2),
            }

        docking_results.append(r)

        raw_results.append({
            'receptor': str(rec_path), 'ligand': str(lig_path),
            'output': str(out_path), 'log': str(vina_log_path),
            'status': status_str,
        })

        icon = "✓" if status_str == 'ok' else "✗"
        aff_str = f"{best_aff:.2f} kcal/mol" if best_aff is not None else "N/A"
        print(f"  [{i}/{len(jobs)}] {icon} {rec_name} x {lig_name} -> {status_str} "
              f"| {n_poses} poses | best: {aff_str} | {elapsed:.1f}s")

        # Append to docking log incrementally
        append_log_rows(
            [r], csv_log_path, existing_combos,
            lig_props_cache or {}, prot_props_cache or {},
        )

    total_time = time.time() - start_all
    n_success = sum(1 for r in docking_results if r.status == "success")
    n_failed = sum(1 for r in docking_results if r.status == "failed")
    n_skipped = sum(1 for r in docking_results if r.status == "skipped")

    print(f"\nDocking complete in {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"  Success: {n_success}  |  Failed: {n_failed}  |  Skipped: {n_skipped}")

    save_error_log(effective_output_dir, error_log)

    results_df = pd.DataFrame(raw_results)
    display(results_df)
    try:
        print(results_df['status'].value_counts().to_dict())
    except Exception:
        pass
    return results_df, docking_results

# %% [markdown]
# # MGL-TOOLS

# %% [markdown]
# ## Prepare Ligands

# %%
# # Prepare ligands from first drugs_dir (adjust if needed for each dir)
# ligands_folder = drugs_dirs[0]  # Use first ligand directory for preparation
# 
# ligand_outputs = run_workflow(
#     input_dir=ligands_folder,
#     contains="ligands",
#     output_dir=output_folder_mgltools,
#     skip_pdb_validation=True, 
#     process_postfixes=False,
#     repair_terminals=False,
#     converter="mgltools",
#     convert_ligands=True,  # Enable ligand SDF → PDBQT conversion
#     log_file=Path("logs/ligand_conversion_overview_mgltools.txt"),
# )

# %% [markdown]
# ## Prepare Proteins

# %%
# protein_outputs = run_workflow(
#     input_dir=orai_folder,
#     contains="proteins",
#     output_dir=output_folder_mgltools,
#     skip_pdb_validation=True, 
#     process_postfixes=False,
#     repair_terminals=False,
#     converter="mgltools",
#     convert_proteins=True,  # Enable protein conversion
#     log_file=Path("logs/protein_conversion_overview_mgltools.txt"),
# )

# %% [markdown]
# # MEEKO

# %% [markdown]
# ## Prepare Ligands

# %%
# Check if ligands have already been prepared (PDBQT files exist in output dir)
existing_lig_pdbqts = sorted(output_folder_meeko.glob("*_ligand*.pdbqt")) + sorted(
    (output_folder_meeko / "ligands").glob("*.pdbqt")
) if (output_folder_meeko / "ligands").exists() else sorted(output_folder_meeko.glob("*_ligand*.pdbqt"))

# Also check for any .pdbqt that doesn't have a matching .box.txt (those are ligands, not proteins)
if not existing_lig_pdbqts:
    existing_lig_pdbqts = sorted(
        p for p in output_folder_meeko.glob("*.pdbqt")
        if not p.with_suffix(".box.txt").exists()
        and not Path(str(p).replace(".pdbqt", ".box.txt")).exists()
    )

if existing_lig_pdbqts:
    print(f"✓ Found {len(existing_lig_pdbqts)} existing ligand PDBQT files in {output_folder_meeko} — skipping preparation.")
    # Reconstruct ligand_outputs dict matching run_workflow structure
    ligand_outputs = {
        "meeko_ligand_outputs": {str(p): str(p) for p in existing_lig_pdbqts},
        "converted_ligands": {},
        "ligand_files": set(),
        "protein_files": set(),
        "protein_results": {},
        "fixed_proteins": [],
        "log_entries": [],
    }
    for p in existing_lig_pdbqts:
        print(f"  • {p.name}")
else:
    print("No existing ligand PDBQTs found — running Meeko preparation...")
    ligand_outputs = run_workflow(
        input_dir=drugs_dirs[0],  # Use first ligand directory for preparation
        contains="ligands",
        output_dir=output_folder_meeko,
        process_postfixes=False,
        convert_ligands_with_meeko=True,
        log_file=Path("logs/ligand_conversion_overview_ligands.txt"),
    )

# %% [markdown]
# ## Prepare Proteins

# %%
# Check if proteins have already been prepared (PDBQT + box.txt pairs exist)
existing_prot_pdbqts = sorted(
    p for p in output_folder_meeko.glob("*.pdbqt")
    if p.with_suffix(".box.txt").exists() or Path(str(p).replace(".pdbqt", ".box.txt")).exists()
)

if existing_prot_pdbqts:
    print(f"✓ Found {len(existing_prot_pdbqts)} existing protein PDBQT+box pairs in {output_folder_meeko} — skipping preparation.")
    # Reconstruct protein_outputs dict matching run_workflow structure
    protein_outputs = {
        "meeko_outputs": {str(p): str(p) for p in existing_prot_pdbqts},
        "meeko_box_files": {str(p): str(p.with_suffix(".box.txt")) for p in existing_prot_pdbqts},
        "converter_protein_outputs": {str(p): str(p) for p in existing_prot_pdbqts},
        "converter_box_files": {str(p): str(p.with_suffix(".box.txt")) for p in existing_prot_pdbqts},
        "protein_results": {},
        "fixed_proteins": [str(p) for p in existing_prot_pdbqts],
        "converted_ligands": {},
        "ligand_files": set(),
        "protein_files": set(),
        "log_entries": [],
    }
    for p in existing_prot_pdbqts:
        box = p.with_suffix(".box.txt")
        print(f"  • {p.name}  +  {box.name}")
else:
    print("No existing protein PDBQTs found — running preparation...")
    protein_outputs = run_workflow(
        input_dir=orai_folder,
        contains="proteins",
        output_dir=output_folder_meeko,
        skip_pdb_validation=False,
        process_postfixes=False,
        repair_terminals=False,
        converter="openbabel",
        convert_proteins=True,
        log_file=Path("logs/protein_conversion_overview_obabel.txt"),
    )

# %%
from pathlib import Path
from pprint import pprint

def build_prepared_manifest(output_dir: Path, protein_data: dict | None, ligand_data: dict | None):
    output_dir = Path(output_dir).resolve()
    manifest = {"proteins": [], "ligands": []}

    def record(section: str, source: str, file_path: str | None):
        if not file_path:
            return
        path_obj = Path(file_path).resolve()
        try:
            path_obj.relative_to(output_dir)
        except ValueError:
            return
        manifest[section].append({
            "name": path_obj.name,
            "path": str(path_obj),
            "filetype": path_obj.suffix.lstrip(".").lower() or "unknown",
            "source": source,
        })

    if protein_data:
        for original, prepared in (protein_data.get("protein_results") or {}).items():
            record("proteins", original, prepared)
        for original, pdbqt in (protein_data.get("meeko_outputs") or {}).items():
            record("proteins", f"{original} [pdbqt]", pdbqt)
        for original, box_file in (protein_data.get("meeko_box_files") or {}).items():
            record("proteins", f"{original} [box]", box_file)

    if ligand_data:
        for original, formats in (ligand_data.get("converted_ligands") or {}).items():
            for fmt, path in formats.items():
                record("ligands", f"{original} [{fmt}]", path)
        for original, pdbqt in (ligand_data.get("meeko_ligand_outputs") or {}).items():
            record("ligands", f"{original} [pdbqt]", pdbqt)

    return manifest

# %%
# prepared_manifest = build_prepared_manifest(output_folder_mgltools, protein_outputs, ligand_outputs)
# print(f"Proteins prepared: {len(prepared_manifest['proteins'])}")
# print(f"Ligands prepared: {len(prepared_manifest['ligands'])}")
# pprint(prepared_manifest)

# %%
prepared_manifest = build_prepared_manifest(output_folder_meeko, protein_outputs, ligand_outputs)
print(f"Proteins prepared: {len(prepared_manifest['proteins'])}")
print(f"Ligands prepared: {len(prepared_manifest['ligands'])}")
pprint(prepared_manifest)

# %% [markdown]
# # Pre-flight Check

# %%
errors = []

# Check Vina-GPU executable
if not VINA_BIN.exists():
    errors.append(f"✗ Vina-GPU executable not found: {VINA_BIN}")
elif not os.access(VINA_BIN, os.X_OK):
    errors.append(f"✗ Vina-GPU executable not executable: {VINA_BIN}")

# Check OpenCL kernel files
kernel1 = OPENCL_BINARY_PATH / "Kernel1_Opt.bin"
kernel2 = OPENCL_BINARY_PATH / "Kernel2_Opt.bin"
if not kernel1.exists() or not kernel2.exists():
    errors.append(f"✗ OpenCL kernel binaries missing in: {OPENCL_BINARY_PATH}")

# Check inputs exist for each ligand dir
for drugs_dir in drugs_dirs:
    ligand_files = collect_files(output_folder_meeko, [".pdbqt"])
    print(f"  {drugs_dir.name}: {len(ligand_files)} ligand files")
    if not ligand_files:
        errors.append(f"✗ No ligand files in {drugs_dir}")

recep_files = collect_files(receptors_dir, [".pdb", ".pdbqt"])
print(f"  Receptors: {len(recep_files)} files in {receptors_dir}")
if not recep_files:
    errors.append(f"✗ No receptor files in {receptors_dir}")

# Check GPU availability
try:
    gpu_result = subprocess.run(
        ["nvidia-smi", "--query-gpu=gpu_name", "--format=csv,noheader"],
        capture_output=True, text=True, timeout=5,
    )
    gpu_name = gpu_result.stdout.strip()
    if not gpu_name:
        errors.append("✗ No GPU detected via nvidia-smi")
    else:
        print(f"  GPU: {gpu_name}")
except Exception:
    errors.append("✗ nvidia-smi not available — GPU required for Vina-GPU")

if errors:
    print("\n⚠ ISSUES FOUND:")
    for e in errors:
        print(f"  {e}")
    print("\nFix these before running docking.")
else:
    print("\n✓ All checks passed — ready to dock with Vina-GPU.")

# %%
def _monitor_progress(output_dir: Path, stop_event: threading.Event, poll_interval: float = 10.0):
    """Background thread that periodically prints docking progress."""
    start_time = time.time()

    while not stop_event.is_set():
        stop_event.wait(poll_interval)
        if stop_event.is_set():
            break

        elapsed = time.time() - start_time

        # Count completed output files
        n_done = len(list(output_dir.rglob("*_vina_out.pdbqt")))

        # GPU utilization
        try:
            gpu_info = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total",
                 "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            gpu_line = gpu_info.stdout.strip()
        except Exception:
            gpu_line = "N/A"

        mins, secs = divmod(int(elapsed), 60)
        hrs, mins = divmod(mins, 60)
        time_str = f"{hrs}h{mins:02d}m{secs:02d}s" if hrs else f"{mins}m{secs:02d}s"

        print(f"  ⏱ {time_str} | ✓ {n_done} done | GPU: {gpu_line}")

# %% [markdown]
# # AutoDock Vina-GPU Docking

# %% [markdown]
# ## MGL-Tools Docking (multi-directory loop)

# %%
# def _monitor_progress(output_dir: Path, stop_event: threading.Event, poll_interval: float = 10.0):
#     """Background thread that periodically prints docking progress."""
#     start_time = time.time()
# 
#     while not stop_event.is_set():
#         stop_event.wait(poll_interval)
#         if stop_event.is_set():
#             break
# 
#         elapsed = time.time() - start_time
# 
#         # Count completed output files
#         n_done = len(list(output_dir.rglob("*_vina_out.pdbqt")))
# 
#         # GPU utilization
#         try:
#             gpu_info = subprocess.run(
#                 ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total",
#                  "--format=csv,noheader,nounits"],
#                 capture_output=True, text=True, timeout=5
#             )
#             gpu_line = gpu_info.stdout.strip()
#         except Exception:
#             gpu_line = "N/A"
# 
#         mins, secs = divmod(int(elapsed), 60)
#         hrs, mins = divmod(mins, 60)
#         time_str = f"{hrs}h{mins:02d}m{secs:02d}s" if hrs else f"{mins}m{secs:02d}s"
# 
#         print(f"  ⏱ {time_str} | ✓ {n_done} done | GPU: {gpu_line}")
# 
# 
# all_vina_results: List[DockingResult] = []
# 
# for drugs_dir in drugs_dirs:
#     VINA_OUTPUT_DIR = VINA_OUTPUT_BASE / drugs_dir.name
#     VINA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
# 
#     print(f"\n{'=' * 80}")
#     print(f"Processing ligand directory (Vina-GPU): {drugs_dir.name}")
#     print(f"Output: {VINA_OUTPUT_DIR}")
#     print(f"{'=' * 80}")
# 
#     # Collect ligand files for property pre-computation
#     lig_files = collect_files(drugs_dir, [".pdbqt", ".sdf", ".mol2"])
#     rec_files = collect_files(receptors_dir, [".pdb", ".pdbqt"])
# 
#     print(f"Found {len(rec_files)} receptors and {len(lig_files)} ligands")
# 
#     if not rec_files or not lig_files:
#         print("⚠ No receptors or ligands found, skipping.")
#         continue
# 
#     # Pre-compute molecular properties
#     lig_props_cache, prot_props_cache = precompute_properties(rec_files, lig_files)
# 
#     # Start progress monitor
#     stop_monitor = threading.Event()
#     monitor_thread = threading.Thread(
#         target=_monitor_progress,
#         args=(VINA_OUTPUT_DIR, stop_monitor, 15),
#         daemon=True,
#     )
#     monitor_thread.start()
# 
#     try:
#         results_df, docking_results = run_autodock_vina_gpu(
#             base_dir=output_folder_mgltools,
#             log_dir=wdlogs,
#             prepared_manifest=prepared_manifest,
#             protein_workflow_data=protein_outputs,
#             ligand_workflow_data=ligand_outputs,
#             thread=THREAD,
#             search_depth=SEARCH_DEPTH,
#             num_modes=NUM_MODES,
#             energy_range=ENERGY_RANGE,
#             overwrite=OVERWRITE_EXISTING,
#             timeout=TIMEOUT_PER_COMPLEX,
#             output_base_dir=VINA_OUTPUT_DIR,
#             lig_props_cache=lig_props_cache,
#             prot_props_cache=prot_props_cache,
#         )
#     finally:
#         stop_monitor.set()
#         monitor_thread.join(timeout=5)
# 
#     all_vina_results.extend(docking_results)
# 
#     # Generate and print summary
#     summary = generate_summary(docking_results)
#     print_summary(summary)
# 
#     # Save summary JSON
#     summary_path = VINA_OUTPUT_DIR / "docking_summary.json"
#     with open(summary_path, 'w') as f:
#         json.dump(summary, f, indent=2, default=str)
#     print(f"\nSummary saved: {summary_path}")
# 
# print(f"\n{'=' * 80}")
# print(f"ALL DONE (Vina-GPU) — {len(all_vina_results)} total results across {len(drugs_dirs)} ligand directories")
# print(f"{'=' * 80}")

# %% [markdown]
# ## Meeko Docking (multi-directory loop)

# %%
all_vina_results_meeko: List[DockingResult] = []

# The actual docking output goes to base_dir/docking, not VINA_OUTPUT_DIR
dock_output_dir = (output_folder_meeko / "docking").resolve()

for drugs_dir in drugs_dirs:
    VINA_OUTPUT_DIR = VINA_OUTPUT_BASE / (drugs_dir.name + "_meeko")
    VINA_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n{'=' * 80}")
    print(f"Processing ligand directory (Meeko / Vina-GPU): {drugs_dir.name}")
    print(f"Output: {VINA_OUTPUT_DIR}")
    print(f"{'=' * 80}")

    lig_files = collect_files(output_folder_meeko, [".pdbqt"])
    rec_files = collect_files(receptors_dir, [".pdb", ".pdbqt"])

    if not rec_files or not lig_files:
        print("⚠ No receptors or ligands found, skipping.")
        continue

    lig_props_cache, prot_props_cache = precompute_properties(rec_files, lig_files)

    stop_monitor = threading.Event()
    monitor_thread = threading.Thread(
        target=_monitor_progress,
        args=(dock_output_dir, stop_monitor, 15),  # Watch actual dock output dir
        daemon=True,
    )
    monitor_thread.start()

    try:
        results_df_meeko, docking_results_meeko = run_autodock_vina_gpu(
            base_dir=output_folder_meeko,
            log_dir=wdlogs,
            prepared_manifest=prepared_manifest,
            protein_workflow_data=protein_outputs,
            ligand_workflow_data=ligand_outputs,
            thread=THREAD,
            search_depth=SEARCH_DEPTH,
            num_modes=NUM_MODES,
            energy_range=ENERGY_RANGE,
            overwrite=OVERWRITE_EXISTING,
            timeout=TIMEOUT_PER_COMPLEX,
            output_base_dir=VINA_OUTPUT_DIR,
            lig_props_cache=lig_props_cache,
            prot_props_cache=prot_props_cache,
        )
    finally:
        stop_monitor.set()
        monitor_thread.join(timeout=5)

    all_vina_results_meeko.extend(docking_results_meeko)

    summary = generate_summary(docking_results_meeko)
    print_summary(summary)

    summary_path = VINA_OUTPUT_DIR / "docking_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\nSummary saved: {summary_path}")

print(f"\n{'=' * 80}")
print(f"ALL DONE (Meeko / Vina-GPU) — {len(all_vina_results_meeko)} total results across {len(drugs_dirs)} ligand directories")
print(f"{'=' * 80}")

# %% [markdown]
# # Results

# %%
# Combine all results for display — handle case when MGL-Tools docking is skipped
all_vina_results = all_vina_results if 'all_vina_results' in dir() else []
all_combined_results = all_vina_results + all_vina_results_meeko

results_data = []
for r in all_combined_results:
    results_data.append({
        "protein": r.protein_name,
        "ligand": r.ligand_name,
        "status": r.status,
        "num_poses": r.num_poses,
        "best_affinity_kcal": r.best_affinity,
        "elapsed_time_s": round(r.elapsed_time, 2),
        "output_dir": str(r.output_dir),
    })

results_df_all = pd.DataFrame(results_data)

print("=" * 80)
print("VINA-GPU DOCKING RESULTS DATAFRAME")
print("=" * 80)
print(f"\nTotal combinations: {len(results_df_all)}")

if not results_df_all.empty:
    print("\nStatus breakdown:")
    print(results_df_all['status'].value_counts().to_string())
    print(f"\nAffinity statistics (successful dockings):")
    success = results_df_all[results_df_all['status'] == 'success']
    if not success.empty and success['best_affinity_kcal'].notna().any():
        print(f"  Mean best affinity: {success['best_affinity_kcal'].mean():.2f} kcal/mol")
        print(f"  Min (best): {success['best_affinity_kcal'].min():.2f} kcal/mol")
        print(f"  Max (worst): {success['best_affinity_kcal'].max():.2f} kcal/mol")
    print(f"\nTiming statistics:")
    docked = results_df_all[results_df_all['status'].isin(['success', 'failed'])]
    if not docked.empty:
        print(f"  Mean time per complex: {docked['elapsed_time_s'].mean():.1f}s")
        print(f"  Total time: {docked['elapsed_time_s'].sum():.1f}s ({docked['elapsed_time_s'].sum()/60:.1f} min)")

display(results_df_all)

# Save combined results CSV
combined_csv = VINA_OUTPUT_BASE / "vina_gpu_combined_results.csv"
results_df_all.to_csv(combined_csv, index=False)
print(f"\nCombined results saved: {combined_csv}")

# %% [markdown]
# # Docking Log

# %%
# Display the docking log for the last output directory
for drugs_dir in drugs_dirs:
    for suffix in [drugs_dir.name, drugs_dir.name + "_meeko"]:
        log_dir_check = VINA_OUTPUT_BASE / suffix
        log_path = log_dir_check / DOCKING_LOG_FILENAME

        if not log_path.exists():
            continue

        docking_log_df = pd.read_csv(log_path)

        print("=" * 80)
        print(f"DOCKING LOG: {log_dir_check.name}")
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
        lig_cols = [c for c in ["ligand_name", "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
                    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
                    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula"] if c in docking_log_df.columns]
        if lig_cols:
            lig_summary = docking_log_df[lig_cols].drop_duplicates(subset=["ligand_name"]).sort_values("ligand_name")
            print(f"\n{'─' * 80}")
            print("LIGAND PROPERTIES:")
            print(f"{'─' * 80}")
            display(lig_summary.reset_index(drop=True))

        # Protein property summary
        prot_cols = [c for c in ["protein_name", "prot_num_residues", "prot_num_atoms", "prot_num_chains"] if c in docking_log_df.columns]
        if prot_cols:
            prot_summary = docking_log_df[prot_cols].drop_duplicates(subset=["protein_name"]).sort_values("protein_name")
            print(f"\n{'─' * 80}")
            print("PROTEIN PROPERTIES:")
            print(f"{'─' * 80}")
            display(prot_summary.reset_index(drop=True))

        # Timing summary
        docked = docking_log_df[docking_log_df["status"].isin(["success", "failed"])]
        if not docked.empty and "elapsed_time_s" in docked.columns:
            print(f"\n{'─' * 80}")
            print("TIMING:")
            print(f"{'─' * 80}")
            print(f"  Mean time per complex: {docked['elapsed_time_s'].mean():.1f}s")
            print(f"  Max time: {docked['elapsed_time_s'].max():.1f}s")
            print(f"  Total: {docked['elapsed_time_s'].sum():.1f}s ({docked['elapsed_time_s'].sum()/60:.1f} min)")

        # Affinity summary
        if "best_affinity_kcal" in docking_log_df.columns:
            success_log = docking_log_df[docking_log_df["status"] == "success"]
            if not success_log.empty and success_log["best_affinity_kcal"].notna().any():
                print(f"\n{'─' * 80}")
                print("AFFINITY:")
                print(f"{'─' * 80}")
                print(f"  Mean best affinity: {success_log['best_affinity_kcal'].mean():.2f} kcal/mol")
                print(f"  Best (most negative): {success_log['best_affinity_kcal'].min():.2f} kcal/mol")
                print(f"  Worst: {success_log['best_affinity_kcal'].max():.2f} kcal/mol")

        print()

# %%
import importlib
import Scripts.compare_meeko_mgltools_docking

importlib.reload(Scripts.compare_meeko_mgltools_docking)

cdt = Scripts.compare_meeko_mgltools_docking.CompareDockingTools(
    meeko_dir=output_folder_meeko / "docking",  # Make sure this path points to where the log files are
    mgltools_dir=output_folder_mgltools / "docking",
    output_dir=wd / "docking_comparison",
)
cdt.analyze()

# %%
"""
Compare docked poses from MGLTools vs Meeko using RMSD calculation.
Calculates RMSD between best poses for each protein-ligand pair.
"""
import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")

def parse_pdbqt_poses(pdbqt_path):
    """Parse a multi-model PDBQT file and extract coordinates for each pose."""
    poses = []
    
    if not pdbqt_path.exists():
        return poses
    
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    
    # Split by MODEL/ENDMDL
    models = re.split(r'MODEL\s+\d+', content)
    
    for model in models:
        if not model.strip():
            continue
        
        atoms = []
        atom_names = []
        
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atom_name = line[12:16].strip()
                    atoms.append([x, y, z])
                    atom_names.append(atom_name)
                except (ValueError, IndexError):
                    continue
        
        if atoms:
            poses.append({
                'coords': np.array(atoms),
                'atom_names': atom_names,
                'num_atoms': len(atoms)
            })
    
    return poses


def calculate_rmsd(coords1, coords2):
    """Calculate RMSD between two coordinate arrays."""
    if coords1.shape != coords2.shape:
        # Try to match by minimum atoms
        min_atoms = min(len(coords1), len(coords2))
        coords1 = coords1[:min_atoms]
        coords2 = coords2[:min_atoms]
    
    diff = coords1 - coords2
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd


def calculate_rmsd_with_alignment(coords1, coords2):
    """Calculate RMSD with Kabsch alignment (optimal superposition)."""
    if coords1.shape != coords2.shape:
        min_atoms = min(len(coords1), len(coords2))
        coords1 = coords1[:min_atoms]
        coords2 = coords2[:min_atoms]
    
    # Center both structures
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2
    
    # Kabsch algorithm for optimal rotation
    H = coords1_centered.T @ coords2_centered
    U, S, Vt = np.linalg.svd(H)
    
    # Handle reflection case
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    
    R = Vt.T @ D @ U.T
    coords2_aligned = coords2_centered @ R
    
    # Calculate RMSD
    diff = coords1_centered - coords2_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd


def get_normalized_protein_name(name):
    """Normalize protein name by removing _retry1 suffix."""
    return name.replace('_retry1', '')


# Collect all docked pose files
meeko_files = {f.stem.replace('_vina_out', ''): f 
               for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_files = {f.stem.replace('_vina_out', ''): f 
                  for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}

print(f"Found {len(meeko_files)} Meeko docked poses")
print(f"Found {len(mgltools_files)} MGLTools docked poses")
print()

# Create mapping with normalized protein names
meeko_normalized = {}
for key, path in meeko_files.items():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = path

mgltools_normalized = {}
for key, path in mgltools_files.items():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = path

# Find common pairs
common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())
print(f"Common protein-ligand pairs for RMSD comparison: {len(common_pairs)}")
print()

# Calculate RMSD for each pair
results = []

print("=" * 100)
print("RMSD COMPARISON: MEEKO vs MGLTools DOCKED POSES")
print("=" * 100)
print()
print(f"{'Protein':<30} {'Ligand':<25} {'RMSD (Å)':>12} {'Aligned RMSD':>14} {'Atoms':>8}")
print("-" * 100)

for protein, ligand in sorted(common_pairs):
    meeko_path = meeko_normalized[(protein, ligand)]
    mgltools_path = mgltools_normalized[(protein, ligand)]
    
    # Parse poses (get best pose = first one)
    meeko_poses = parse_pdbqt_poses(meeko_path)
    mgltools_poses = parse_pdbqt_poses(mgltools_path)
    
    if not meeko_poses or not mgltools_poses:
        print(f"{protein:<30} {ligand:<25} {'N/A':>12} {'N/A':>14} {'N/A':>8}")
        continue
    
    # Compare best poses (mode 1)
    meeko_best = meeko_poses[0]
    mgltools_best = mgltools_poses[0]
    
    # Calculate RMSD (direct and aligned)
    try:
        rmsd_direct = calculate_rmsd(meeko_best['coords'], mgltools_best['coords'])
        rmsd_aligned = calculate_rmsd_with_alignment(meeko_best['coords'], mgltools_best['coords'])
        num_atoms = min(meeko_best['num_atoms'], mgltools_best['num_atoms'])
        
        print(f"{protein:<30} {ligand:<25} {rmsd_direct:>12.3f} {rmsd_aligned:>14.3f} {num_atoms:>8}")
        
        results.append({
            'protein': protein,
            'ligand': ligand,
            'rmsd_direct': rmsd_direct,
            'rmsd_aligned': rmsd_aligned,
            'meeko_atoms': meeko_best['num_atoms'],
            'mgltools_atoms': mgltools_best['num_atoms'],
            'atoms_compared': num_atoms
        })
    except Exception as e:
        print(f"{protein:<30} {ligand:<25} {'Error':>12} {str(e)[:14]:>14}")

print()

# Summary statistics
if results:
    df_rmsd = pd.DataFrame(results)
    
    print("=" * 100)
    print("RMSD SUMMARY STATISTICS")
    print("=" * 100)
    print()
    print(f"Total pairs compared: {len(df_rmsd)}")
    print()
    print("Direct RMSD (without alignment):")
    print(f"  Mean:   {df_rmsd['rmsd_direct'].mean():.3f} Å")
    print(f"  Median: {df_rmsd['rmsd_direct'].median():.3f} Å")
    print(f"  Min:    {df_rmsd['rmsd_direct'].min():.3f} Å")
    print(f"  Max:    {df_rmsd['rmsd_direct'].max():.3f} Å")
    print(f"  Std:    {df_rmsd['rmsd_direct'].std():.3f} Å")
    print()
    print("Aligned RMSD (with Kabsch superposition):")
    print(f"  Mean:   {df_rmsd['rmsd_aligned'].mean():.3f} Å")
    print(f"  Median: {df_rmsd['rmsd_aligned'].median():.3f} Å")
    print(f"  Min:    {df_rmsd['rmsd_aligned'].min():.3f} Å")
    print(f"  Max:    {df_rmsd['rmsd_aligned'].max():.3f} Å")
    print(f"  Std:    {df_rmsd['rmsd_aligned'].std():.3f} Å")
    print()
    
    # Classification based on RMSD thresholds
    print("=" * 100)
    print("POSE SIMILARITY CLASSIFICATION (based on aligned RMSD)")
    print("=" * 100)
    print()
    very_similar = (df_rmsd['rmsd_aligned'] < 2.0).sum()
    similar = ((df_rmsd['rmsd_aligned'] >= 2.0) & (df_rmsd['rmsd_aligned'] < 4.0)).sum()
    different = (df_rmsd['rmsd_aligned'] >= 4.0).sum()
    
    print(f"Very similar poses (RMSD < 2.0 Å):     {very_similar} ({100*very_similar/len(df_rmsd):.1f}%)")
    print(f"Similar poses (2.0 Å ≤ RMSD < 4.0 Å): {similar} ({100*similar/len(df_rmsd):.1f}%)")
    print(f"Different poses (RMSD ≥ 4.0 Å):       {different} ({100*different/len(df_rmsd):.1f}%)")
    print()
    
    # Group by ligand
    print("=" * 100)
    print("RMSD BY LIGAND")
    print("=" * 100)
    print()
    for ligand in df_rmsd['ligand'].unique():
        ligand_data = df_rmsd[df_rmsd['ligand'] == ligand]
        print(f"Ligand: {ligand}")
        print(f"  Mean aligned RMSD: {ligand_data['rmsd_aligned'].mean():.3f} Å")
        print(f"  Range: {ligand_data['rmsd_aligned'].min():.3f} - {ligand_data['rmsd_aligned'].max():.3f} Å")
        print()
    
    # Save results
    output_csv = wdlogs / "meeko_vs_mgltools_rmsd_comparison.csv"
    df_rmsd.to_csv(output_csv, index=False)
    print(f"RMSD comparison saved to: {output_csv}")
    
    # Display DataFrame
    display(df_rmsd)

# %%
"""
Comprehensive overview of docked positions for each protein-ligand pair
from both MGLTools and Meeko conversion methods.
"""
import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")

def parse_vina_log(log_path):
    """Parse Vina log file to get affinity scores for all modes."""
    results = {'modes': []}
    if not log_path.exists():
        return results
    
    with open(log_path, 'r') as f:
        content = f.read()
    
    mode_pattern = r'^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
    modes = re.findall(mode_pattern, content, re.MULTILINE)
    
    for mode in modes:
        mode_num, affinity, rmsd_lb, rmsd_ub = mode
        results['modes'].append({
            'mode': int(mode_num),
            'affinity': float(affinity),
            'rmsd_lb': float(rmsd_lb),
            'rmsd_ub': float(rmsd_ub)
        })
    return results

def parse_pdbqt_poses(pdbqt_path):
    """Parse PDBQT file to get coordinates and centroid for each pose."""
    poses = []
    if not pdbqt_path.exists():
        return poses
    
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    
    models = re.split(r'MODEL\s+\d+', content)
    
    for model in models:
        if not model.strip():
            continue
        
        atoms = []
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        
        if atoms:
            coords = np.array(atoms)
            centroid = np.mean(coords, axis=0)
            poses.append({
                'coords': coords,
                'centroid': centroid,
                'num_atoms': len(atoms),
                'min_coords': coords.min(axis=0),
                'max_coords': coords.max(axis=0)
            })
    
    return poses

def get_normalized_protein_name(name):
    """Normalize protein name by removing _retry1 suffix."""
    return name.replace('_retry1', '')

# Collect files
meeko_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
meeko_logs = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.log')}
mgltools_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_logs = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.log')}

# Create normalized mappings
meeko_normalized = {}
for key in meeko_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = key

mgltools_normalized = {}
for key in mgltools_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = key

# Find common pairs
common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())

print("=" * 120)
print("COMPREHENSIVE OVERVIEW: DOCKED POSITIONS FOR EACH PROTEIN-LIGAND PAIR")
print("=" * 120)
print()

all_data = []

for protein, ligand in sorted(common_pairs):
    print("=" * 120)
    print(f"PROTEIN: {protein}  |  LIGAND: {ligand}")
    print("=" * 120)
    
    # Get file keys
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    # Parse poses
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    
    # Parse logs for affinities
    meeko_log = parse_vina_log(meeko_logs.get(meeko_key, Path('')))
    mgltools_log = parse_vina_log(mgltools_logs.get(mgltools_key, Path('')))
    
    print()
    print(f"{'Mode':<6} {'Method':<10} {'Affinity':>10} {'Centroid X':>12} {'Centroid Y':>12} {'Centroid Z':>12} {'Atoms':>8}")
    print("-" * 80)
    
    # Display up to 9 modes for each method
    max_modes = min(9, max(len(meeko_poses), len(mgltools_poses)))
    
    for mode_idx in range(max_modes):
        # MGLTools pose
        if mode_idx < len(mgltools_poses):
            mgl_pose = mgltools_poses[mode_idx]
            mgl_affinity = mgltools_log['modes'][mode_idx]['affinity'] if mode_idx < len(mgltools_log['modes']) else None
            aff_str = f"{mgl_affinity:.3f}" if mgl_affinity else "N/A"
            print(f"{mode_idx+1:<6} {'MGLTools':<10} {aff_str:>10} {mgl_pose['centroid'][0]:>12.3f} {mgl_pose['centroid'][1]:>12.3f} {mgl_pose['centroid'][2]:>12.3f} {mgl_pose['num_atoms']:>8}")
            
            all_data.append({
                'protein': protein,
                'ligand': ligand,
                'mode': mode_idx + 1,
                'method': 'MGLTools',
                'affinity': mgl_affinity,
                'centroid_x': mgl_pose['centroid'][0],
                'centroid_y': mgl_pose['centroid'][1],
                'centroid_z': mgl_pose['centroid'][2],
                'min_x': mgl_pose['min_coords'][0],
                'min_y': mgl_pose['min_coords'][1],
                'min_z': mgl_pose['min_coords'][2],
                'max_x': mgl_pose['max_coords'][0],
                'max_y': mgl_pose['max_coords'][1],
                'max_z': mgl_pose['max_coords'][2],
                'num_atoms': mgl_pose['num_atoms']
            })
        
        # Meeko pose
        if mode_idx < len(meeko_poses):
            meeko_pose = meeko_poses[mode_idx]
            meeko_affinity = meeko_log['modes'][mode_idx]['affinity'] if mode_idx < len(meeko_log['modes']) else None
            aff_str = f"{meeko_affinity:.3f}" if meeko_affinity else "N/A"
            print(f"{mode_idx+1:<6} {'Meeko':<10} {aff_str:>10} {meeko_pose['centroid'][0]:>12.3f} {meeko_pose['centroid'][1]:>12.3f} {meeko_pose['centroid'][2]:>12.3f} {meeko_pose['num_atoms']:>8}")
            
            all_data.append({
                'protein': protein,
                'ligand': ligand,
                'mode': mode_idx + 1,
                'method': 'Meeko',
                'affinity': meeko_affinity,
                'centroid_x': meeko_pose['centroid'][0],
                'centroid_y': meeko_pose['centroid'][1],
                'centroid_z': meeko_pose['centroid'][2],
                'min_x': meeko_pose['min_coords'][0],
                'min_y': meeko_pose['min_coords'][1],
                'min_z': meeko_pose['min_coords'][2],
                'max_x': meeko_pose['max_coords'][0],
                'max_y': meeko_pose['max_coords'][1],
                'max_z': meeko_pose['max_coords'][2],
                'num_atoms': meeko_pose['num_atoms']
            })
        
        # Calculate distance between centroids for this mode
        if mode_idx < len(mgltools_poses) and mode_idx < len(meeko_poses):
            dist = np.linalg.norm(mgltools_poses[mode_idx]['centroid'] - meeko_poses[mode_idx]['centroid'])
            print(f"       {'Δ Distance':<10} {dist:>10.3f} Å between centroids")
        print()
    
    # Summary for this pair
    if meeko_poses and mgltools_poses:
        best_meeko = meeko_log['modes'][0]['affinity'] if meeko_log['modes'] else None
        best_mgltools = mgltools_log['modes'][0]['affinity'] if mgltools_log['modes'] else None
        
        print(f"  SUMMARY:")
        print(f"    Best MGLTools affinity: {best_mgltools:.3f} kcal/mol" if best_mgltools else "    Best MGLTools affinity: N/A")
        print(f"    Best Meeko affinity:    {best_meeko:.3f} kcal/mol" if best_meeko else "    Best Meeko affinity: N/A")
        
        if best_meeko and best_mgltools:
            delta = best_meeko - best_mgltools
            better = "Meeko" if delta < 0 else "MGLTools"
            print(f"    Δ Affinity: {delta:+.3f} kcal/mol ({better} is better)")
        
        # Distance between best poses
        dist_best = np.linalg.norm(mgltools_poses[0]['centroid'] - meeko_poses[0]['centroid'])
        print(f"    Distance between best pose centroids: {dist_best:.3f} Å")
    print()

# Create DataFrame and display
df_overview = pd.DataFrame(all_data)

print("=" * 120)
print("SUMMARY TABLE: BEST POSES (MODE 1) COMPARISON")
print("=" * 120)
print()

# Pivot for easy comparison
df_best = df_overview[df_overview['mode'] == 1].copy()
df_pivot = df_best.pivot_table(
    index=['protein', 'ligand'],
    columns='method',
    values=['affinity', 'centroid_x', 'centroid_y', 'centroid_z'],
    aggfunc='first'
)

# Calculate centroid distance
comparison_summary = []
for (protein, ligand) in df_best.groupby(['protein', 'ligand']).groups.keys():
    mgl_row = df_best[(df_best['protein'] == protein) & (df_best['ligand'] == ligand) & (df_best['method'] == 'MGLTools')]
    meeko_row = df_best[(df_best['protein'] == protein) & (df_best['ligand'] == ligand) & (df_best['method'] == 'Meeko')]
    
    if len(mgl_row) > 0 and len(meeko_row) > 0:
        mgl = mgl_row.iloc[0]
        meeko = meeko_row.iloc[0]
        
        dist = np.sqrt((mgl['centroid_x'] - meeko['centroid_x'])**2 + 
                      (mgl['centroid_y'] - meeko['centroid_y'])**2 + 
                      (mgl['centroid_z'] - meeko['centroid_z'])**2)
        
        comparison_summary.append({
            'protein': protein,
            'ligand': ligand,
            'mgltools_affinity': mgl['affinity'],
            'meeko_affinity': meeko['affinity'],
            'affinity_delta': meeko['affinity'] - mgl['affinity'],
            'mgltools_centroid': f"({mgl['centroid_x']:.1f}, {mgl['centroid_y']:.1f}, {mgl['centroid_z']:.1f})",
            'meeko_centroid': f"({meeko['centroid_x']:.1f}, {meeko['centroid_y']:.1f}, {meeko['centroid_z']:.1f})",
            'centroid_distance': dist
        })

df_summary = pd.DataFrame(comparison_summary)
display(df_summary)

# Save detailed overview
output_csv = wdlogs / "docked_poses_overview.csv"
df_overview.to_csv(output_csv, index=False)
print(f"\nDetailed overview saved to: {output_csv}")

# %%
"""
Binding Pocket Analysis: Compare poses in the same binding pockets
Only compare MGLTools vs Meeko poses when they dock to similar locations.
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")

# Pocket clustering threshold (Angstroms) - poses within this distance are considered same pocket
POCKET_THRESHOLD = 10.0

def parse_vina_log(log_path):
    """Parse Vina log file to get affinity scores for all modes."""
    results = {'modes': []}
    if not log_path.exists():
        return results
    with open(log_path, 'r') as f:
        content = f.read()
    mode_pattern = r'^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
    modes = re.findall(mode_pattern, content, re.MULTILINE)
    for mode in modes:
        mode_num, affinity, rmsd_lb, rmsd_ub = mode
        results['modes'].append({
            'mode': int(mode_num),
            'affinity': float(affinity),
            'rmsd_lb': float(rmsd_lb),
            'rmsd_ub': float(rmsd_ub)
        })
    return results

def parse_pdbqt_poses(pdbqt_path):
    """Parse PDBQT file to get coordinates and centroid for each pose."""
    poses = []
    if not pdbqt_path.exists():
        return poses
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    models = re.split(r'MODEL\s+\d+', content)
    for model in models:
        if not model.strip():
            continue
        atoms = []
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        if atoms:
            coords = np.array(atoms)
            centroid = np.mean(coords, axis=0)
            poses.append({
                'coords': coords,
                'centroid': centroid,
                'num_atoms': len(atoms)
            })
    return poses

def get_normalized_protein_name(name):
    return name.replace('_retry1', '')

def cluster_poses_by_pocket(all_centroids, threshold=POCKET_THRESHOLD):
    """Cluster poses into binding pockets based on centroid distance."""
    if len(all_centroids) <= 1:
        return [0] * len(all_centroids)
    
    centroids_array = np.array(all_centroids)
    distances = pdist(centroids_array)
    
    if len(distances) == 0:
        return [0] * len(all_centroids)
    
    linkage_matrix = linkage(distances, method='average')
    clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return clusters.tolist()

# Collect files
meeko_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
meeko_logs = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.log')}
mgltools_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_logs = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.log')}

# Create normalized mappings
meeko_normalized = {}
for key in meeko_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = key

mgltools_normalized = {}
for key in mgltools_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = key

# Find common pairs
common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())

print("=" * 120)
print("POSE COUNT OVERVIEW: MGLTools vs Meeko")
print("=" * 120)
print()
print(f"{'Protein':<25} {'Ligand':<25} {'MGLTools':>10} {'Meeko':>10} {'Total':>10}")
print("-" * 90)

pose_counts = []
for protein, ligand in sorted(common_pairs):
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    
    n_meeko = len(meeko_poses)
    n_mgltools = len(mgltools_poses)
    
    print(f"{protein:<25} {ligand:<25} {n_mgltools:>10} {n_meeko:>10} {n_mgltools + n_meeko:>10}")
    
    pose_counts.append({
        'protein': protein,
        'ligand': ligand,
        'mgltools_poses': n_mgltools,
        'meeko_poses': n_meeko,
        'total_poses': n_mgltools + n_meeko
    })

df_counts = pd.DataFrame(pose_counts)
print("-" * 90)
print(f"{'TOTAL':<25} {'':<25} {df_counts['mgltools_poses'].sum():>10} {df_counts['meeko_poses'].sum():>10} {df_counts['total_poses'].sum():>10}")
print()

# Binding pocket analysis
print("=" * 120)
print(f"BINDING POCKET ANALYSIS (clustering threshold: {POCKET_THRESHOLD} Å)")
print("=" * 120)
print()

pocket_analysis = []
detailed_pocket_data = []

for protein, ligand in sorted(common_pairs):
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    meeko_log = parse_vina_log(meeko_logs.get(meeko_key, Path('')))
    mgltools_log = parse_vina_log(mgltools_logs.get(mgltools_key, Path('')))
    
    # Combine all centroids for clustering
    all_centroids = []
    pose_info = []
    
    for i, pose in enumerate(mgltools_poses):
        all_centroids.append(pose['centroid'])
        affinity = mgltools_log['modes'][i]['affinity'] if i < len(mgltools_log['modes']) else None
        pose_info.append({'method': 'MGLTools', 'mode': i+1, 'affinity': affinity, 'centroid': pose['centroid']})
    
    for i, pose in enumerate(meeko_poses):
        all_centroids.append(pose['centroid'])
        affinity = meeko_log['modes'][i]['affinity'] if i < len(meeko_log['modes']) else None
        pose_info.append({'method': 'Meeko', 'mode': i+1, 'affinity': affinity, 'centroid': pose['centroid']})
    
    if not all_centroids:
        continue
    
    # Cluster into pockets
    clusters = cluster_poses_by_pocket(all_centroids)
    
    # Assign clusters to poses
    for idx, cluster_id in enumerate(clusters):
        pose_info[idx]['pocket'] = cluster_id
    
    # Analyze pockets
    unique_pockets = sorted(set(clusters))
    
    print(f"\n{'='*100}")
    print(f"PROTEIN: {protein}  |  LIGAND: {ligand}")
    print(f"{'='*100}")
    print(f"Identified {len(unique_pockets)} distinct binding pocket(s)")
    print()
    
    for pocket_id in unique_pockets:
        pocket_poses = [p for p in pose_info if p['pocket'] == pocket_id]
        mgl_in_pocket = [p for p in pocket_poses if p['method'] == 'MGLTools']
        meeko_in_pocket = [p for p in pocket_poses if p['method'] == 'Meeko']
        
        # Calculate pocket centroid
        pocket_centroids = [p['centroid'] for p in pocket_poses]
        pocket_center = np.mean(pocket_centroids, axis=0)
        
        print(f"  POCKET {pocket_id}: Center ({pocket_center[0]:.1f}, {pocket_center[1]:.1f}, {pocket_center[2]:.1f})")
        print(f"    MGLTools poses: {len(mgl_in_pocket)}, Meeko poses: {len(meeko_in_pocket)}")
        
        # Show best affinity from each method in this pocket
        if mgl_in_pocket:
            best_mgl = min(mgl_in_pocket, key=lambda x: x['affinity'] if x['affinity'] else float('inf'))
            print(f"    Best MGLTools: Mode {best_mgl['mode']}, Affinity: {best_mgl['affinity']:.3f} kcal/mol" if best_mgl['affinity'] else "")
        if meeko_in_pocket:
            best_meeko = min(meeko_in_pocket, key=lambda x: x['affinity'] if x['affinity'] else float('inf'))
            print(f"    Best Meeko:    Mode {best_meeko['mode']}, Affinity: {best_meeko['affinity']:.3f} kcal/mol" if best_meeko['affinity'] else "")
        
        # Compare if both methods have poses in this pocket
        if mgl_in_pocket and meeko_in_pocket:
            best_mgl_aff = min(p['affinity'] for p in mgl_in_pocket if p['affinity'])
            best_meeko_aff = min(p['affinity'] for p in meeko_in_pocket if p['affinity'])
            delta = best_meeko_aff - best_mgl_aff
            better = "Meeko" if delta < 0 else "MGLTools"
            print(f"    → SAME POCKET COMPARISON: Δ = {delta:+.3f} kcal/mol ({better} better)")
            
            detailed_pocket_data.append({
                'protein': protein,
                'ligand': ligand,
                'pocket_id': pocket_id,
                'pocket_center_x': pocket_center[0],
                'pocket_center_y': pocket_center[1],
                'pocket_center_z': pocket_center[2],
                'mgltools_poses': len(mgl_in_pocket),
                'meeko_poses': len(meeko_in_pocket),
                'mgltools_best_affinity': best_mgl_aff,
                'meeko_best_affinity': best_meeko_aff,
                'affinity_delta': delta,
                'better_method': better
            })
        else:
            # Only one method found this pocket
            method_only = "MGLTools" if mgl_in_pocket else "Meeko"
            print(f"    → UNIQUE POCKET: Only found by {method_only}")
            
            best_aff = None
            if mgl_in_pocket:
                best_aff = min(p['affinity'] for p in mgl_in_pocket if p['affinity'])
            elif meeko_in_pocket:
                best_aff = min(p['affinity'] for p in meeko_in_pocket if p['affinity'])
            
            detailed_pocket_data.append({
                'protein': protein,
                'ligand': ligand,
                'pocket_id': pocket_id,
                'pocket_center_x': pocket_center[0],
                'pocket_center_y': pocket_center[1],
                'pocket_center_z': pocket_center[2],
                'mgltools_poses': len(mgl_in_pocket),
                'meeko_poses': len(meeko_in_pocket),
                'mgltools_best_affinity': min((p['affinity'] for p in mgl_in_pocket if p['affinity']), default=None),
                'meeko_best_affinity': min((p['affinity'] for p in meeko_in_pocket if p['affinity']), default=None),
                'affinity_delta': None,
                'better_method': f"Only {method_only}"
            })
        print()
    
    pocket_analysis.append({
        'protein': protein,
        'ligand': ligand,
        'num_pockets': len(unique_pockets),
        'shared_pockets': sum(1 for p in unique_pockets if 
                             any(x['pocket']==p and x['method']=='MGLTools' for x in pose_info) and
                             any(x['pocket']==p and x['method']=='Meeko' for x in pose_info)),
        'mgltools_only_pockets': sum(1 for p in unique_pockets if 
                                     any(x['pocket']==p and x['method']=='MGLTools' for x in pose_info) and
                                     not any(x['pocket']==p and x['method']=='Meeko' for x in pose_info)),
        'meeko_only_pockets': sum(1 for p in unique_pockets if 
                                  any(x['pocket']==p and x['method']=='Meeko' for x in pose_info) and
                                  not any(x['pocket']==p and x['method']=='MGLTools' for x in pose_info))
    })

# Summary tables
print("\n" + "=" * 120)
print("SUMMARY: POCKET DISTRIBUTION BY PROTEIN-LIGAND PAIR")
print("=" * 120)
df_pocket_summary = pd.DataFrame(pocket_analysis)
display(df_pocket_summary)

print("\n" + "=" * 120)
print("DETAILED POCKET COMPARISON (where both methods found the same pocket)")
print("=" * 120)
df_pocket_detail = pd.DataFrame(detailed_pocket_data)
df_shared = df_pocket_detail[df_pocket_detail['affinity_delta'].notna()]
if len(df_shared) > 0:
    display(df_shared)
    
    print("\n" + "=" * 80)
    print("SAME-POCKET COMPARISON STATISTICS")
    print("=" * 80)
    meeko_wins = (df_shared['better_method'] == 'Meeko').sum()
    mgl_wins = (df_shared['better_method'] == 'MGLTools').sum()
    print(f"Shared pockets analyzed: {len(df_shared)}")
    print(f"Meeko better: {meeko_wins} ({100*meeko_wins/len(df_shared):.1f}%)")
    print(f"MGLTools better: {mgl_wins} ({100*mgl_wins/len(df_shared):.1f}%)")
    print(f"Mean Δ affinity: {df_shared['affinity_delta'].mean():.3f} kcal/mol")
else:
    print("No shared pockets found between methods.")

# Save results
df_pocket_detail.to_csv(wdlogs / "binding_pocket_analysis.csv", index=False)
print(f"\nBinding pocket analysis saved to: binding_pocket_analysis.csv")

# %%
"""
Visualization: Compare MGLTools vs Meeko Docked Poses by Binding Pocket
- Number of poses per pocket for each method
- RMSD comparison for first n-ranked poses (n = min poses from both methods)
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")
POCKET_THRESHOLD = 10.0

def parse_pdbqt_poses(pdbqt_path):
    """Parse PDBQT file to get coordinates and centroid for each pose."""
    poses = []
    if not pdbqt_path.exists():
        return poses
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    models = re.split(r'MODEL\s+\d+', content)
    for model in models:
        if not model.strip():
            continue
        atoms = []
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        if atoms:
            coords = np.array(atoms)
            centroid = np.mean(coords, axis=0)
            poses.append({'coords': coords, 'centroid': centroid, 'num_atoms': len(atoms)})
    return poses

def parse_vina_log(log_path):
    """Parse Vina log file to get affinity scores."""
    results = {'modes': []}
    if not log_path.exists():
        return results
    with open(log_path, 'r') as f:
        content = f.read()
    mode_pattern = r'^\s*(\d+)\s+([\-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
    modes = re.findall(mode_pattern, content, re.MULTILINE)
    for mode in modes:
        mode_num, affinity, rmsd_lb, rmsd_ub = mode
        results['modes'].append({
            'mode': int(mode_num),
            'affinity': float(affinity),
            'rmsd_lb': float(rmsd_lb),
            'rmsd_ub': float(rmsd_ub)
        })
    return results

def get_normalized_protein_name(name):
    return name.replace('_retry1', '')

def cluster_poses_by_pocket(all_centroids, threshold=POCKET_THRESHOLD):
    if len(all_centroids) <= 1:
        return [0] * len(all_centroids)
    centroids_array = np.array(all_centroids)
    distances = pdist(centroids_array)
    if len(distances) == 0:
        return [0] * len(all_centroids)
    linkage_matrix = linkage(distances, method='average')
    clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return clusters.tolist()

def calculate_rmsd_aligned(coords1, coords2):
    """Calculate RMSD with Kabsch alignment."""
    if coords1.shape != coords2.shape:
        min_atoms = min(len(coords1), len(coords2))
        coords1 = coords1[:min_atoms]
        coords2 = coords2[:min_atoms]
    centroid1 = np.mean(coords1, axis=0)
    centroid2 = np.mean(coords2, axis=0)
    coords1_centered = coords1 - centroid1
    coords2_centered = coords2 - centroid2
    H = coords1_centered.T @ coords2_centered
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    coords2_aligned = coords2_centered @ R
    diff = coords1_centered - coords2_aligned
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd

# Collect files
meeko_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
meeko_logs = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.log')}
mgltools_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_logs = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.log')}

# Create normalized mappings
meeko_normalized = {}
for key in meeko_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = key

mgltools_normalized = {}
for key in mgltools_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = key

common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())

# Analyze each protein-ligand pair
pocket_pose_counts = []  # For bar chart: poses per pocket
rmsd_data = []           # For RMSD comparison

for protein, ligand in sorted(common_pairs):
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    meeko_log = parse_vina_log(meeko_logs.get(meeko_key, Path('')))
    mgltools_log = parse_vina_log(mgltools_logs.get(mgltools_key, Path('')))
    
    if not meeko_poses or not mgltools_poses:
        continue
    
    # Combine all centroids for clustering
    all_centroids = []
    pose_info = []
    
    for i, pose in enumerate(mgltools_poses):
        all_centroids.append(pose['centroid'])
        affinity = mgltools_log['modes'][i]['affinity'] if i < len(mgltools_log['modes']) else None
        pose_info.append({'method': 'MGLTools', 'mode': i+1, 'affinity': affinity, 
                          'centroid': pose['centroid'], 'coords': pose['coords']})
    
    for i, pose in enumerate(meeko_poses):
        all_centroids.append(pose['centroid'])
        affinity = meeko_log['modes'][i]['affinity'] if i < len(meeko_log['modes']) else None
        pose_info.append({'method': 'Meeko', 'mode': i+1, 'affinity': affinity,
                          'centroid': pose['centroid'], 'coords': pose['coords']})
    
    # Cluster into pockets
    clusters = cluster_poses_by_pocket(all_centroids)
    for idx, cluster_id in enumerate(clusters):
        pose_info[idx]['pocket'] = cluster_id
    
    # Analyze each pocket
    unique_pockets = sorted(set(clusters))
    
    for pocket_id in unique_pockets:
        mgl_pocket_poses = sorted([p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'MGLTools'],
                                   key=lambda x: x['mode'])
        meeko_pocket_poses = sorted([p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'Meeko'],
                                     key=lambda x: x['mode'])
        
        n_mgl = len(mgl_pocket_poses)
        n_meeko = len(meeko_pocket_poses)
        
        pocket_pose_counts.append({
            'protein': protein,
            'ligand': ligand,
            'pocket': pocket_id,
            'pair_pocket': f"{protein}_{ligand}_P{pocket_id}",
            'mgltools_poses': n_mgl,
            'meeko_poses': n_meeko
        })
        
        # Calculate RMSD for first n-ranked poses (n = min of both)
        if n_mgl > 0 and n_meeko > 0:
            n_compare = min(n_mgl, n_meeko)
            for rank in range(n_compare):
                mgl_pose = mgl_pocket_poses[rank]
                meeko_pose = meeko_pocket_poses[rank]
                try:
                    rmsd = calculate_rmsd_aligned(mgl_pose['coords'], meeko_pose['coords'])
                    rmsd_data.append({
                        'protein': protein,
                        'ligand': ligand,
                        'pocket': pocket_id,
                        'pair_pocket': f"{protein}_{ligand}_P{pocket_id}",
                        'rank': rank + 1,
                        'rmsd': rmsd,
                        'mgltools_affinity': mgl_pose['affinity'],
                        'meeko_affinity': meeko_pose['affinity']
                    })
                except Exception:
                    pass

df_pocket_counts = pd.DataFrame(pocket_pose_counts)
df_rmsd = pd.DataFrame(rmsd_data)

print(f"Total pockets analyzed: {len(df_pocket_counts)}")
print(f"Total RMSD comparisons: {len(df_rmsd)}")

# ============ FIGURE 1: Poses per Pocket (Grouped Bar Chart) ============
fig1, ax1 = plt.subplots(figsize=(14, 6))

x = np.arange(len(df_pocket_counts))
width = 0.35

bars1 = ax1.bar(x - width/2, df_pocket_counts['mgltools_poses'], width, 
                label='MGLTools', color='#2ecc71', edgecolor='black')
bars2 = ax1.bar(x + width/2, df_pocket_counts['meeko_poses'], width, 
                label='Meeko', color='#3498db', edgecolor='black')

ax1.set_xlabel('Protein-Ligand Pair / Pocket', fontsize=12)
ax1.set_ylabel('Number of Poses', fontsize=12)
ax1.set_title('Number of Docked Poses per Binding Pocket: MGLTools vs Meeko', fontsize=14, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(df_pocket_counts['pair_pocket'], rotation=45, ha='right', fontsize=8)
ax1.legend(fontsize=11)
ax1.grid(axis='y', alpha=0.3)

for bar in bars1:
    height = bar.get_height()
    if height > 0:
        ax1.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                     xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=7)
for bar in bars2:
    height = bar.get_height()
    if height > 0:
        ax1.annotate(f'{int(height)}', xy=(bar.get_x() + bar.get_width()/2, height),
                     xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=7)

plt.tight_layout()
plt.show()

# ============ FIGURE 2: RMSD by Rank for Each Pocket ============
if len(df_rmsd) > 0:
    unique_pair_pockets = df_rmsd['pair_pocket'].unique()
    n_pockets = len(unique_pair_pockets)
    
    if n_pockets <= 6:
        ncols = min(3, n_pockets)
        nrows = (n_pockets + ncols - 1) // ncols
    else:
        ncols = 4
        nrows = (n_pockets + ncols - 1) // ncols
    
    fig2, axes2 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), squeeze=False)
    axes2_flat = axes2.flatten()
    
    for idx, pair_pocket in enumerate(sorted(unique_pair_pockets)):
        ax = axes2_flat[idx]
        pocket_df = df_rmsd[df_rmsd['pair_pocket'] == pair_pocket].sort_values('rank')
        
        ranks = pocket_df['rank'].values
        rmsds = pocket_df['rmsd'].values
        
        bars = ax.bar(ranks, rmsds, color='#e74c3c', edgecolor='black', alpha=0.8)
        ax.axhline(y=2.0, color='green', linestyle='--', linewidth=1.5, label='2.0 Å threshold')
        ax.axhline(y=4.0, color='orange', linestyle='--', linewidth=1.5, label='4.0 Å threshold')
        
        ax.set_xlabel('Pose Rank', fontsize=10)
        ax.set_ylabel('RMSD (Å)', fontsize=10)
        ax.set_title(pair_pocket, fontsize=10, fontweight='bold')
        ax.set_xticks(ranks)
        ax.grid(axis='y', alpha=0.3)
        
        for bar, rmsd in zip(bars, rmsds):
            ax.annotate(f'{rmsd:.2f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                        xytext=(0, 3), textcoords='offset points', ha='center', va='bottom', fontsize=8)
    
    # Hide unused subplots
    for idx in range(n_pockets, len(axes2_flat)):
        axes2_flat[idx].set_visible(False)
    
    axes2_flat[0].legend(loc='upper right', fontsize=8)
    fig2.suptitle('RMSD between MGLTools and Meeko Poses by Rank (Aligned)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.show()

# ============ FIGURE 3: Summary Box Plot of RMSD by Pocket ============
if len(df_rmsd) > 0:
    fig3, ax3 = plt.subplots(figsize=(12, 5))
    
    pockets_order = sorted(df_rmsd['pair_pocket'].unique())
    rmsd_by_pocket = [df_rmsd[df_rmsd['pair_pocket'] == p]['rmsd'].values for p in pockets_order]
    
    bp = ax3.boxplot(rmsd_by_pocket, labels=pockets_order, patch_artist=True)
    
    for patch in bp['boxes']:
        patch.set_facecolor('#9b59b6')
        patch.set_alpha(0.7)
    
    ax3.axhline(y=2.0, color='green', linestyle='--', linewidth=2, label='2.0 Å (very similar)')
    ax3.axhline(y=4.0, color='orange', linestyle='--', linewidth=2, label='4.0 Å (similar)')
    
    ax3.set_xlabel('Protein-Ligand Pair / Pocket', fontsize=12)
    ax3.set_ylabel('RMSD (Å)', fontsize=12)
    ax3.set_title('RMSD Distribution: MGLTools vs Meeko Poses per Pocket', fontsize=14, fontweight='bold')
    ax3.set_xticklabels(pockets_order, rotation=45, ha='right', fontsize=9)
    ax3.legend(fontsize=10)
    ax3.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.show()

# ============ FIGURE 4: Heatmap of RMSD by Protein-Ligand and Rank ============
if len(df_rmsd) > 0:
    pivot_df = df_rmsd.pivot_table(index='pair_pocket', columns='rank', values='rmsd', aggfunc='mean')
    
    fig4, ax4 = plt.subplots(figsize=(10, max(6, len(pivot_df)*0.4)))
    
    im = ax4.imshow(pivot_df.values, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=max(6, pivot_df.values.max()))
    
    ax4.set_xticks(np.arange(len(pivot_df.columns)))
    ax4.set_yticks(np.arange(len(pivot_df.index)))
    ax4.set_xticklabels([f'Rank {c}' for c in pivot_df.columns], fontsize=10)
    ax4.set_yticklabels(pivot_df.index, fontsize=9)
    
    for i in range(len(pivot_df.index)):
        for j in range(len(pivot_df.columns)):
            val = pivot_df.values[i, j]
            if not np.isnan(val):
                ax4.text(j, i, f'{val:.2f}', ha='center', va='center', 
                         color='white' if val > 3 else 'black', fontsize=9, fontweight='bold')
    
    cbar = plt.colorbar(im, ax=ax4)
    cbar.set_label('RMSD (Å)', fontsize=11)
    
    ax4.set_xlabel('Pose Rank', fontsize=12)
    ax4.set_ylabel('Protein-Ligand Pair / Pocket', fontsize=12)
    ax4.set_title('RMSD Heatmap: MGLTools vs Meeko by Rank and Pocket', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.show()

# Summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)
print(f"\nPose counts per pocket:")
display(df_pocket_counts)

if len(df_rmsd) > 0:
    print(f"\nRMSD Statistics (all comparisons):")
    print(f"  Mean RMSD: {df_rmsd['rmsd'].mean():.3f} Å")
    print(f"  Median RMSD: {df_rmsd['rmsd'].median():.3f} Å")
    print(f"  Min RMSD: {df_rmsd['rmsd'].min():.3f} Å")
    print(f"  Max RMSD: {df_rmsd['rmsd'].max():.3f} Å")
    print(f"\nRMSD by Rank:")
    print(df_rmsd.groupby('rank')['rmsd'].agg(['mean', 'std', 'count']).round(3))
    print(f"\nDetailed RMSD data:")
    display(df_rmsd)

# %%
"""
Pocket Overview: Compare binding pocket distribution between MGLTools and Meeko
- Total pockets per file type
- Overlapping vs unique pockets
- Distance between unique pockets
- Pose counts per pocket type
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from collections import defaultdict

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")
POCKET_THRESHOLD = 10.0

def parse_pdbqt_poses(pdbqt_path):
    poses = []
    if not pdbqt_path.exists():
        return poses
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    models = re.split(r'MODEL\s+\d+', content)
    for model in models:
        if not model.strip():
            continue
        atoms = []
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        if atoms:
            coords = np.array(atoms)
            centroid = np.mean(coords, axis=0)
            poses.append({'coords': coords, 'centroid': centroid, 'num_atoms': len(atoms)})
    return poses

def get_normalized_protein_name(name):
    return name.replace('_retry1', '')

def cluster_poses_by_pocket(all_centroids, threshold=POCKET_THRESHOLD):
    if len(all_centroids) <= 1:
        return [0] * len(all_centroids)
    centroids_array = np.array(all_centroids)
    distances = pdist(centroids_array)
    if len(distances) == 0:
        return [0] * len(all_centroids)
    linkage_matrix = linkage(distances, method='average')
    clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return clusters.tolist()

# Collect files
meeko_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}

# Create normalized mappings
meeko_normalized = {}
for key in meeko_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = key

mgltools_normalized = {}
for key in mgltools_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = key

common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())

# Analyze pockets for each protein-ligand pair
all_pocket_data = []
pair_summaries = []

for protein, ligand in sorted(common_pairs):
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    
    if not meeko_poses or not mgltools_poses:
        continue
    
    # Combine all centroids for clustering
    all_centroids = []
    pose_info = []
    
    for i, pose in enumerate(mgltools_poses):
        all_centroids.append(pose['centroid'])
        pose_info.append({'method': 'MGLTools', 'mode': i+1, 'centroid': pose['centroid']})
    
    for i, pose in enumerate(meeko_poses):
        all_centroids.append(pose['centroid'])
        pose_info.append({'method': 'Meeko', 'mode': i+1, 'centroid': pose['centroid']})
    
    # Cluster into pockets
    clusters = cluster_poses_by_pocket(all_centroids)
    for idx, cluster_id in enumerate(clusters):
        pose_info[idx]['pocket'] = cluster_id
    
    unique_pockets = sorted(set(clusters))
    
    mgl_only_pockets = 0
    meeko_only_pockets = 0
    shared_pockets = 0
    
    for pocket_id in unique_pockets:
        mgl_in_pocket = [p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'MGLTools']
        meeko_in_pocket = [p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'Meeko']
        
        pocket_centroids = [p['centroid'] for p in pose_info if p['pocket'] == pocket_id]
        pocket_center = np.mean(pocket_centroids, axis=0)
        
        has_mgl = len(mgl_in_pocket) > 0
        has_meeko = len(meeko_in_pocket) > 0
        
        if has_mgl and has_meeko:
            pocket_type = 'Shared'
            shared_pockets += 1
        elif has_mgl:
            pocket_type = 'MGLTools Only'
            mgl_only_pockets += 1
        else:
            pocket_type = 'Meeko Only'
            meeko_only_pockets += 1
        
        all_pocket_data.append({
            'protein': protein,
            'ligand': ligand,
            'pair': f"{protein}_{ligand}",
            'pocket_id': pocket_id,
            'pocket_type': pocket_type,
            'mgltools_poses': len(mgl_in_pocket),
            'meeko_poses': len(meeko_in_pocket),
            'total_poses': len(mgl_in_pocket) + len(meeko_in_pocket),
            'center_x': pocket_center[0],
            'center_y': pocket_center[1],
            'center_z': pocket_center[2]
        })
    
    pair_summaries.append({
        'protein': protein,
        'ligand': ligand,
        'pair': f"{protein}_{ligand}",
        'total_pockets': len(unique_pockets),
        'shared_pockets': shared_pockets,
        'mgltools_only': mgl_only_pockets,
        'meeko_only': meeko_only_pockets
    })

df_pockets = pd.DataFrame(all_pocket_data)
df_pairs = pd.DataFrame(pair_summaries)

# ============ FIGURE 1: Overview - Pockets per File Type ============
fig1, axes1 = plt.subplots(1, 3, figsize=(15, 5))

# 1a: Total pockets contributed by each method
total_shared = df_pairs['shared_pockets'].sum()
total_mgl_only = df_pairs['mgltools_only'].sum()
total_meeko_only = df_pairs['meeko_only'].sum()

# Pockets that MGLTools finds (shared + mgl_only)
mgl_total = total_shared + total_mgl_only
meeko_total = total_shared + total_meeko_only

ax1a = axes1[0]
methods = ['MGLTools', 'Meeko']
pocket_counts = [mgl_total, meeko_total]
colors = ['#2ecc71', '#3498db']
bars = ax1a.bar(methods, pocket_counts, color=colors, edgecolor='black', linewidth=1.5)
ax1a.set_ylabel('Number of Pockets Found', fontsize=12)
ax1a.set_title('Total Pockets per File Type', fontsize=13, fontweight='bold')
for bar, count in zip(bars, pocket_counts):
    ax1a.annotate(f'{count}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                  xytext=(0, 5), textcoords='offset points', ha='center', fontsize=14, fontweight='bold')
ax1a.grid(axis='y', alpha=0.3)

# 1b: Pocket breakdown (Shared vs Unique)
ax1b = axes1[1]
categories = ['Shared\n(Both)', 'MGLTools\nOnly', 'Meeko\nOnly']
counts = [total_shared, total_mgl_only, total_meeko_only]
colors_breakdown = ['#9b59b6', '#2ecc71', '#3498db']
bars = ax1b.bar(categories, counts, color=colors_breakdown, edgecolor='black', linewidth=1.5)
ax1b.set_ylabel('Number of Pockets', fontsize=12)
ax1b.set_title('Pocket Classification', fontsize=13, fontweight='bold')
for bar, count in zip(bars, counts):
    ax1b.annotate(f'{count}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                  xytext=(0, 5), textcoords='offset points', ha='center', fontsize=14, fontweight='bold')
ax1b.grid(axis='y', alpha=0.3)

# 1c: Pie chart of pocket distribution
ax1c = axes1[2]
sizes = [total_shared, total_mgl_only, total_meeko_only]
labels = [f'Shared\n({total_shared})', f'MGLTools Only\n({total_mgl_only})', f'Meeko Only\n({total_meeko_only})']
explode = (0.05, 0.05, 0.05)
ax1c.pie(sizes, explode=explode, labels=labels, colors=colors_breakdown, autopct='%1.1f%%',
         shadow=True, startangle=90, textprops={'fontsize': 11})
ax1c.set_title('Pocket Distribution', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()

# ============ FIGURE 2: Pockets per Protein-Ligand Pair (Stacked) ============
fig2, ax2 = plt.subplots(figsize=(14, 6))

x = np.arange(len(df_pairs))
width = 0.6

bars_shared = ax2.bar(x, df_pairs['shared_pockets'], width, label='Shared', color='#9b59b6', edgecolor='black')
bars_mgl = ax2.bar(x, df_pairs['mgltools_only'], width, bottom=df_pairs['shared_pockets'], 
                   label='MGLTools Only', color='#2ecc71', edgecolor='black')
bars_meeko = ax2.bar(x, df_pairs['meeko_only'], width, 
                     bottom=df_pairs['shared_pockets'] + df_pairs['mgltools_only'],
                     label='Meeko Only', color='#3498db', edgecolor='black')

ax2.set_xlabel('Protein-Ligand Pair', fontsize=12)
ax2.set_ylabel('Number of Pockets', fontsize=12)
ax2.set_title('Pocket Distribution per Protein-Ligand Pair', fontsize=14, fontweight='bold')
ax2.set_xticks(x)
ax2.set_xticklabels(df_pairs['pair'], rotation=45, ha='right', fontsize=9)
ax2.legend(fontsize=11)
ax2.grid(axis='y', alpha=0.3)

# Add total count on top
for i, (s, m, e) in enumerate(zip(df_pairs['shared_pockets'], df_pairs['mgltools_only'], df_pairs['meeko_only'])):
    total = s + m + e
    ax2.annotate(f'{total}', xy=(i, total), xytext=(0, 3), textcoords='offset points', 
                 ha='center', fontsize=9, fontweight='bold')

plt.tight_layout()
plt.show()

# ============ FIGURE 3: Poses per Pocket Type ============
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5))

# 3a: Total poses by pocket type
ax3a = axes3[0]
pose_by_type = df_pockets.groupby('pocket_type').agg({
    'mgltools_poses': 'sum',
    'meeko_poses': 'sum',
    'total_poses': 'sum'
}).reindex(['Shared', 'MGLTools Only', 'Meeko Only'])

x_types = np.arange(len(pose_by_type))
width = 0.35

bars_mgl = ax3a.bar(x_types - width/2, pose_by_type['mgltools_poses'], width, 
                    label='MGLTools Poses', color='#2ecc71', edgecolor='black')
bars_meeko = ax3a.bar(x_types + width/2, pose_by_type['meeko_poses'], width,
                      label='Meeko Poses', color='#3498db', edgecolor='black')

ax3a.set_xlabel('Pocket Type', fontsize=12)
ax3a.set_ylabel('Total Number of Poses', fontsize=12)
ax3a.set_title('Docked Poses by Pocket Type', fontsize=13, fontweight='bold')
ax3a.set_xticks(x_types)
ax3a.set_xticklabels(['Shared', 'MGLTools\nOnly', 'Meeko\nOnly'], fontsize=10)
ax3a.legend(fontsize=10)
ax3a.grid(axis='y', alpha=0.3)

for bar in bars_mgl:
    h = bar.get_height()
    if h > 0:
        ax3a.annotate(f'{int(h)}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')
for bar in bars_meeko:
    h = bar.get_height()
    if h > 0:
        ax3a.annotate(f'{int(h)}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')

# 3b: Average poses per pocket by type
ax3b = axes3[1]
avg_by_type = df_pockets.groupby('pocket_type').agg({
    'mgltools_poses': 'mean',
    'meeko_poses': 'mean'
}).reindex(['Shared', 'MGLTools Only', 'Meeko Only'])

bars_mgl_avg = ax3b.bar(x_types - width/2, avg_by_type['mgltools_poses'], width,
                        label='MGLTools Poses', color='#2ecc71', edgecolor='black')
bars_meeko_avg = ax3b.bar(x_types + width/2, avg_by_type['meeko_poses'], width,
                          label='Meeko Poses', color='#3498db', edgecolor='black')

ax3b.set_xlabel('Pocket Type', fontsize=12)
ax3b.set_ylabel('Average Poses per Pocket', fontsize=12)
ax3b.set_title('Average Poses per Pocket by Type', fontsize=13, fontweight='bold')
ax3b.set_xticks(x_types)
ax3b.set_xticklabels(['Shared', 'MGLTools\nOnly', 'Meeko\nOnly'], fontsize=10)
ax3b.legend(fontsize=10)
ax3b.grid(axis='y', alpha=0.3)

for bar in bars_mgl_avg:
    h = bar.get_height()
    if h > 0:
        ax3b.annotate(f'{h:.1f}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10)
for bar in bars_meeko_avg:
    h = bar.get_height()
    if h > 0:
        ax3b.annotate(f'{h:.1f}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10)

plt.tight_layout()
plt.show()

# ============ FIGURE 4: Distance Between Unique Pockets ============
# Calculate distances between unique pockets within each protein-ligand pair
distance_data = []

for pair in df_pockets['pair'].unique():
    pair_pockets = df_pockets[df_pockets['pair'] == pair]
    
    mgl_only = pair_pockets[pair_pockets['pocket_type'] == 'MGLTools Only']
    meeko_only = pair_pockets[pair_pockets['pocket_type'] == 'Meeko Only']
    shared = pair_pockets[pair_pockets['pocket_type'] == 'Shared']
    
    # Distance between MGLTools-only and Meeko-only pockets
    for _, mgl_pocket in mgl_only.iterrows():
        mgl_center = np.array([mgl_pocket['center_x'], mgl_pocket['center_y'], mgl_pocket['center_z']])
        
        for _, meeko_pocket in meeko_only.iterrows():
            meeko_center = np.array([meeko_pocket['center_x'], meeko_pocket['center_y'], meeko_pocket['center_z']])
            dist = np.linalg.norm(mgl_center - meeko_center)
            distance_data.append({
                'pair': pair,
                'comparison': 'MGLTools-only vs Meeko-only',
                'distance': dist
            })
        
        # Distance from MGLTools-only to nearest shared pocket
        for _, shared_pocket in shared.iterrows():
            shared_center = np.array([shared_pocket['center_x'], shared_pocket['center_y'], shared_pocket['center_z']])
            dist = np.linalg.norm(mgl_center - shared_center)
            distance_data.append({
                'pair': pair,
                'comparison': 'MGLTools-only vs Shared',
                'distance': dist
            })
    
    # Distance from Meeko-only to nearest shared pocket
    for _, meeko_pocket in meeko_only.iterrows():
        meeko_center = np.array([meeko_pocket['center_x'], meeko_pocket['center_y'], meeko_pocket['center_z']])
        
        for _, shared_pocket in shared.iterrows():
            shared_center = np.array([shared_pocket['center_x'], shared_pocket['center_y'], shared_pocket['center_z']])
            dist = np.linalg.norm(meeko_center - shared_center)
            distance_data.append({
                'pair': pair,
                'comparison': 'Meeko-only vs Shared',
                'distance': dist
            })

df_distances = pd.DataFrame(distance_data)

if len(df_distances) > 0:
    fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))
    
    # 4a: Box plot of distances by comparison type
    ax4a = axes4[0]
    comparison_types = df_distances['comparison'].unique()
    distance_groups = [df_distances[df_distances['comparison'] == c]['distance'].values for c in comparison_types]
    
    bp = ax4a.boxplot(distance_groups, labels=[c.replace(' vs ', '\nvs\n') for c in comparison_types], 
                      patch_artist=True)
    colors_box = ['#e74c3c', '#f39c12', '#1abc9c']
    for patch, color in zip(bp['boxes'], colors_box[:len(bp['boxes'])]):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax4a.axhline(y=POCKET_THRESHOLD, color='red', linestyle='--', linewidth=2, 
                 label=f'Clustering threshold ({POCKET_THRESHOLD} Å)')
    ax4a.set_ylabel('Distance (Å)', fontsize=12)
    ax4a.set_title('Distance Between Unique Pockets', fontsize=13, fontweight='bold')
    ax4a.legend(fontsize=9)
    ax4a.grid(axis='y', alpha=0.3)
    
    # 4b: Histogram of all unique pocket distances
    ax4b = axes4[1]
    all_distances = df_distances['distance'].values
    ax4b.hist(all_distances, bins=20, color='#e74c3c', edgecolor='black', alpha=0.7)
    ax4b.axvline(x=POCKET_THRESHOLD, color='red', linestyle='--', linewidth=2,
                 label=f'Clustering threshold ({POCKET_THRESHOLD} Å)')
    ax4b.axvline(x=np.mean(all_distances), color='blue', linestyle='-', linewidth=2,
                 label=f'Mean ({np.mean(all_distances):.1f} Å)')
    ax4b.set_xlabel('Distance (Å)', fontsize=12)
    ax4b.set_ylabel('Frequency', fontsize=12)
    ax4b.set_title('Distribution of Distances Between Unique Pockets', fontsize=13, fontweight='bold')
    ax4b.legend(fontsize=10)
    ax4b.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.show()
else:
    print("No unique pockets to compare distances.")

# ============ FIGURE 5: Venn-style visualization per pair ============
fig5, ax5 = plt.subplots(figsize=(12, 6))

# Create grouped data for Venn-like comparison
pairs = df_pairs['pair'].values
shared_counts = df_pairs['shared_pockets'].values
mgl_counts = df_pairs['mgltools_only'].values
meeko_counts = df_pairs['meeko_only'].values

x = np.arange(len(pairs))
width = 0.25

ax5.bar(x - width, mgl_counts, width, label='MGLTools Only', color='#2ecc71', edgecolor='black')
ax5.bar(x, shared_counts, width, label='Shared (Overlapping)', color='#9b59b6', edgecolor='black')
ax5.bar(x + width, meeko_counts, width, label='Meeko Only', color='#3498db', edgecolor='black')

ax5.set_xlabel('Protein-Ligand Pair', fontsize=12)
ax5.set_ylabel('Number of Pockets', fontsize=12)
ax5.set_title('Overlapping vs Unique Pockets per Pair', fontsize=14, fontweight='bold')
ax5.set_xticks(x)
ax5.set_xticklabels(pairs, rotation=45, ha='right', fontsize=9)
ax5.legend(fontsize=11)
ax5.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()

# ============ SUMMARY TABLES ============
print("=" * 100)
print("SUMMARY: POCKET OVERVIEW BY FILE TYPE")
print("=" * 100)
print()
print(f"Total unique pockets identified: {len(df_pockets)}")
print(f"  - Shared pockets (found by both): {total_shared} ({100*total_shared/len(df_pockets):.1f}%)")
print(f"  - MGLTools-only pockets: {total_mgl_only} ({100*total_mgl_only/len(df_pockets):.1f}%)")
print(f"  - Meeko-only pockets: {total_meeko_only} ({100*total_meeko_only/len(df_pockets):.1f}%)")
print()
print(f"Pockets found by MGLTools: {mgl_total}")
print(f"Pockets found by Meeko: {meeko_total}")
print()

print("Per Protein-Ligand Pair Summary:")
display(df_pairs)

print("\nDetailed Pocket Data:")
display(df_pockets)

if len(df_distances) > 0:
    print("\nDistance Statistics Between Unique Pockets:")
    print(df_distances.groupby('comparison')['distance'].agg(['mean', 'std', 'min', 'max', 'count']).round(2))

# %% [markdown]
# # Clustering

# %%
"""
Comprehensive Pocket Clustering Analysis: MGLTools vs Meeko
- Cluster pockets and show differences
- Poses per file type in shared vs unique pockets
- Pocket separation distances in the environment
- Ratio of poses in shared vs unique pockets per file format
"""
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.cluster.hierarchy import fcluster, linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict

# Directories
MEEKO_DOCKING_DIR = Path(wd / "docking_ready_meeko/docking")
MGLTOOLS_DOCKING_DIR = Path(wd / "docking_ready_mgltools/docking")
POCKET_THRESHOLD = 10.0

def parse_pdbqt_poses(pdbqt_path):
    poses = []
    if not pdbqt_path.exists():
        return poses
    with open(pdbqt_path, 'r') as f:
        content = f.read()
    models = re.split(r'MODEL\s+\d+', content)
    for model in models:
        if not model.strip():
            continue
        atoms = []
        for line in model.split('\n'):
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append([x, y, z])
                except (ValueError, IndexError):
                    continue
        if atoms:
            coords = np.array(atoms)
            centroid = np.mean(coords, axis=0)
            poses.append({'coords': coords, 'centroid': centroid, 'num_atoms': len(atoms)})
    return poses

def get_normalized_protein_name(name):
    return name.replace('_retry1', '')

def cluster_poses_by_pocket(all_centroids, threshold=POCKET_THRESHOLD):
    if len(all_centroids) <= 1:
        return [1] * len(all_centroids)
    centroids_array = np.array(all_centroids)
    distances = pdist(centroids_array)
    if len(distances) == 0:
        return [1] * len(all_centroids)
    linkage_matrix = linkage(distances, method='average')
    clusters = fcluster(linkage_matrix, t=threshold, criterion='distance')
    return clusters.tolist()

# Collect files
meeko_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MEEKO_DOCKING_DIR.glob('*_vina_out.pdbqt')}
mgltools_pdbqt = {f.stem.replace('_vina_out', ''): f for f in MGLTOOLS_DOCKING_DIR.glob('*_vina_out.pdbqt')}

# Create normalized mappings
meeko_normalized = {}
for key in meeko_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        meeko_normalized[(protein_norm, ligand)] = key

mgltools_normalized = {}
for key in mgltools_pdbqt.keys():
    parts = key.split('__')
    if len(parts) == 2:
        protein_norm = get_normalized_protein_name(parts[0])
        ligand = parts[1]
        mgltools_normalized[(protein_norm, ligand)] = key

common_pairs = set(meeko_normalized.keys()) & set(mgltools_normalized.keys())

# Analyze all pockets across all pairs
all_pocket_data = []
all_pose_data = []
pair_pocket_summaries = []

for protein, ligand in sorted(common_pairs):
    meeko_key = meeko_normalized[(protein, ligand)]
    mgltools_key = mgltools_normalized[(protein, ligand)]
    
    meeko_poses = parse_pdbqt_poses(meeko_pdbqt[meeko_key])
    mgltools_poses = parse_pdbqt_poses(mgltools_pdbqt[mgltools_key])
    
    if not meeko_poses or not mgltools_poses:
        continue
    
    # Combine all centroids for clustering
    all_centroids = []
    pose_info = []
    
    for i, pose in enumerate(mgltools_poses):
        all_centroids.append(pose['centroid'])
        pose_info.append({
            'method': 'MGLTools', 
            'mode': i+1, 
            'centroid': pose['centroid'],
            'protein': protein,
            'ligand': ligand
        })
    
    for i, pose in enumerate(meeko_poses):
        all_centroids.append(pose['centroid'])
        pose_info.append({
            'method': 'Meeko', 
            'mode': i+1, 
            'centroid': pose['centroid'],
            'protein': protein,
            'ligand': ligand
        })
    
    # Cluster into pockets
    clusters = cluster_poses_by_pocket(all_centroids)
    for idx, cluster_id in enumerate(clusters):
        pose_info[idx]['pocket'] = cluster_id
    
    unique_pockets = sorted(set(clusters))
    
    # Analyze each pocket
    shared_poses_mgl = 0
    shared_poses_meeko = 0
    unique_poses_mgl = 0
    unique_poses_meeko = 0
    
    for pocket_id in unique_pockets:
        mgl_in_pocket = [p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'MGLTools']
        meeko_in_pocket = [p for p in pose_info if p['pocket'] == pocket_id and p['method'] == 'Meeko']
        
        pocket_centroids = [p['centroid'] for p in pose_info if p['pocket'] == pocket_id]
        pocket_center = np.mean(pocket_centroids, axis=0)
        
        has_mgl = len(mgl_in_pocket) > 0
        has_meeko = len(meeko_in_pocket) > 0
        
        if has_mgl and has_meeko:
            pocket_type = 'Shared'
            shared_poses_mgl += len(mgl_in_pocket)
            shared_poses_meeko += len(meeko_in_pocket)
        elif has_mgl:
            pocket_type = 'MGLTools Only'
            unique_poses_mgl += len(mgl_in_pocket)
        else:
            pocket_type = 'Meeko Only'
            unique_poses_meeko += len(meeko_in_pocket)
        
        all_pocket_data.append({
            'protein': protein,
            'ligand': ligand,
            'pair': f"{protein}_{ligand}",
            'pocket_id': pocket_id,
            'pocket_type': pocket_type,
            'mgltools_poses': len(mgl_in_pocket),
           'meeko_poses': len(meeko_in_pocket),
            'total_poses': len(mgl_in_pocket) + len(meeko_in_pocket),
            'center_x': pocket_center[0],
            'center_y': pocket_center[1],
            'center_z': pocket_center[2]
        })
        
        # Store individual pose data
        for p in mgl_in_pocket + meeko_in_pocket:
            all_pose_data.append({
                'protein': protein,
                'ligand': ligand,
                'pair': f"{protein}_{ligand}",
                'pocket_id': pocket_id,
                'pocket_type': pocket_type,
                'method': p['method'],
                'mode': p['mode'],
                'centroid_x': p['centroid'][0],
                'centroid_y': p['centroid'][1],
                'centroid_z': p['centroid'][2]
            })
    
    pair_pocket_summaries.append({
        'protein': protein,
        'ligand': ligand,
        'pair': f"{protein}_{ligand}",
        'total_pockets': len(unique_pockets),
        'shared_pockets': sum(1 for p in all_pocket_data if p['pair'] == f"{protein}_{ligand}" and p['pocket_type'] == 'Shared'),
        'mgltools_only_pockets': sum(1 for p in all_pocket_data if p['pair'] == f"{protein}_{ligand}" and p['pocket_type'] == 'MGLTools Only'),
        'meeko_only_pockets': sum(1 for p in all_pocket_data if p['pair'] == f"{protein}_{ligand}" and p['pocket_type'] == 'Meeko Only'),
        'shared_poses_mgl': shared_poses_mgl,
        'shared_poses_meeko': shared_poses_meeko,
        'unique_poses_mgl': unique_poses_mgl,
        'unique_poses_meeko': unique_poses_meeko
    })

df_pockets = pd.DataFrame(all_pocket_data)
df_poses = pd.DataFrame(all_pose_data)
df_summaries = pd.DataFrame(pair_pocket_summaries)

# ============ GLOBAL STATISTICS ============
total_shared_pockets = len(df_pockets[df_pockets['pocket_type'] == 'Shared'])
total_mgl_only_pockets = len(df_pockets[df_pockets['pocket_type'] == 'MGLTools Only'])
total_meeko_only_pockets = len(df_pockets[df_pockets['pocket_type'] == 'Meeko Only'])
total_pockets = len(df_pockets)

# Poses in shared vs unique pockets
shared_pocket_df = df_pockets[df_pockets['pocket_type'] == 'Shared']
mgl_only_pocket_df = df_pockets[df_pockets['pocket_type'] == 'MGLTools Only']
meeko_only_pocket_df = df_pockets[df_pockets['pocket_type'] == 'Meeko Only']

mgl_poses_in_shared = shared_pocket_df['mgltools_poses'].sum()
meeko_poses_in_shared = shared_pocket_df['meeko_poses'].sum()
mgl_poses_in_unique = mgl_only_pocket_df['mgltools_poses'].sum()
meeko_poses_in_unique = meeko_only_pocket_df['meeko_poses'].sum()

total_mgl_poses = mgl_poses_in_shared + mgl_poses_in_unique
total_meeko_poses = meeko_poses_in_shared + meeko_poses_in_unique

print("=" * 100)
print("COMPREHENSIVE POCKET CLUSTERING ANALYSIS: MGLTools vs Meeko")
print("=" * 100)
print()
print(f"Clustering threshold: {POCKET_THRESHOLD} Å")
print(f"Protein-ligand pairs analyzed: {len(common_pairs)}")
print()

print("=" * 80)
print("POCKET OVERVIEW")
print("=" * 80)
print(f"Total pockets identified: {total_pockets}")
print(f"  - Shared pockets (both methods): {total_shared_pockets} ({100*total_shared_pockets/total_pockets:.1f}%)")
print(f"  - MGLTools-only pockets: {total_mgl_only_pockets} ({100*total_mgl_only_pockets/total_pockets:.1f}%)")
print(f"  - Meeko-only pockets: {total_meeko_only_pockets} ({100*total_meeko_only_pockets/total_pockets:.1f}%)")
print()

print("=" * 80)
print("POSES IN SHARED vs UNIQUE POCKETS")
print("=" * 80)
print()
print("MGLTools File Format:")
print(f"  - Poses in SHARED pockets: {mgl_poses_in_shared}")
print(f"  - Poses in UNIQUE pockets: {mgl_poses_in_unique}")
print(f"  - Total MGLTools poses: {total_mgl_poses}")
print(f"  - Ratio (Shared/Unique): {mgl_poses_in_shared/max(1, mgl_poses_in_unique):.2f}")
print(f"  - % in Shared pockets: {100*mgl_poses_in_shared/max(1, total_mgl_poses):.1f}%")
print()
print("Meeko File Format:")
print(f"  - Poses in SHARED pockets: {meeko_poses_in_shared}")
print(f"  - Poses in UNIQUE pockets: {meeko_poses_in_unique}")
print(f"  - Total Meeko poses: {total_meeko_poses}")
print(f"  - Ratio (Shared/Unique): {meeko_poses_in_shared/max(1, meeko_poses_in_unique):.2f}")
print(f"  - % in Shared pockets: {100*meeko_poses_in_shared/max(1, total_meeko_poses):.1f}%")
print()

# ============ FIGURE 1: Pocket Clustering Overview ============
fig1, axes1 = plt.subplots(2, 2, figsize=(14, 12))

# 1a: Pocket type distribution
ax1a = axes1[0, 0]
pocket_types = ['Shared\n(Both)', 'MGLTools\nOnly', 'Meeko\nOnly']
pocket_counts = [total_shared_pockets, total_mgl_only_pockets, total_meeko_only_pockets]
colors = ['#9b59b6', '#2ecc71', '#3498db']
bars = ax1a.bar(pocket_types, pocket_counts, color=colors, edgecolor='black', linewidth=1.5)
ax1a.set_ylabel('Number of Pockets', fontsize=12)
ax1a.set_title('Pocket Distribution by Type', fontsize=13, fontweight='bold')
for bar, count in zip(bars, pocket_counts):
    ax1a.annotate(f'{count}\n({100*count/total_pockets:.1f}%)', 
                  xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                  xytext=(0, 5), textcoords='offset points', ha='center', fontsize=11, fontweight='bold')
ax1a.grid(axis='y', alpha=0.3)

# 1b: Poses per pocket type and file format
ax1b = axes1[0, 1]
x = np.arange(3)
width = 0.35

pose_counts_mgl = [mgl_poses_in_shared, mgl_poses_in_unique, 0]
pose_counts_meeko = [meeko_poses_in_shared, 0, meeko_poses_in_unique]

bars_mgl = ax1b.bar(x - width/2, pose_counts_mgl, width, label='MGLTools', color='#2ecc71', edgecolor='black')
bars_meeko = ax1b.bar(x + width/2, pose_counts_meeko, width, label='Meeko', color='#3498db', edgecolor='black')

ax1b.set_xlabel('Pocket Type', fontsize=12)
ax1b.set_ylabel('Number of Poses', fontsize=12)
ax1b.set_title('Poses per File Type in Each Pocket Category', fontsize=13, fontweight='bold')
ax1b.set_xticks(x)
ax1b.set_xticklabels(['Shared', 'MGLTools Only', 'Meeko Only'])
ax1b.legend(fontsize=11)
ax1b.grid(axis='y', alpha=0.3)

for bar in bars_mgl:
    h = bar.get_height()
    if h > 0:
        ax1b.annotate(f'{int(h)}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')
for bar in bars_meeko:
    h = bar.get_height()
    if h > 0:
        ax1b.annotate(f'{int(h)}', xy=(bar.get_x() + bar.get_width()/2, h),
                      xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')

# 1c: Ratio of poses in shared vs unique pockets
ax1c = axes1[1, 0]
methods = ['MGLTools', 'Meeko']
shared_ratios = [mgl_poses_in_shared/max(1, total_mgl_poses)*100, 
                 meeko_poses_in_shared/max(1, total_meeko_poses)*100]
unique_ratios = [mgl_poses_in_unique/max(1, total_mgl_poses)*100,
                 meeko_poses_in_unique/max(1, total_meeko_poses)*100]

x = np.arange(len(methods))
width = 0.35

bars_shared = ax1c.bar(x - width/2, shared_ratios, width, label='In Shared Pockets', color='#9b59b6', edgecolor='black')
bars_unique = ax1c.bar(x + width/2, unique_ratios, width, label='In Unique Pockets', color='#e74c3c', edgecolor='black')

ax1c.set_ylabel('% of Total Poses', fontsize=12)
ax1c.set_title('Pose Distribution: Shared vs Unique Pockets', fontsize=13, fontweight='bold')
ax1c.set_xticks(x)
ax1c.set_xticklabels(methods, fontsize=11)
ax1c.legend(fontsize=10)
ax1c.grid(axis='y', alpha=0.3)
ax1c.set_ylim(0, 110)

for bar in bars_shared:
    h = bar.get_height()
    ax1c.annotate(f'{h:.1f}%', xy=(bar.get_x() + bar.get_width()/2, h),
                  xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')
for bar in bars_unique:
    h = bar.get_height()
    ax1c.annotate(f'{h:.1f}%', xy=(bar.get_x() + bar.get_width()/2, h),
                  xytext=(0, 3), textcoords='offset points', ha='center', fontsize=10, fontweight='bold')

# 1d: Pie chart - Overall pose distribution
ax1d = axes1[1, 1]
pose_categories = [mgl_poses_in_shared, meeko_poses_in_shared, mgl_poses_in_unique, meeko_poses_in_unique]
category_labels = [f'MGLTools\nin Shared\n({mgl_poses_in_shared})',
                   f'Meeko\nin Shared\n({meeko_poses_in_shared})',
                   f'MGLTools\nUnique\n({mgl_poses_in_unique})',
                   f'Meeko\nUnique\n({meeko_poses_in_unique})']
colors_pie = ['#27ae60', '#2980b9', '#16a085', '#1abc9c']
explode = (0.02, 0.02, 0.05, 0.05)

ax1d.pie(pose_categories, explode=explode, labels=category_labels, colors=colors_pie,
         autopct='%1.1f%%', shadow=True, startangle=90, textprops={'fontsize': 9})
ax1d.set_title('Overall Pose Distribution', fontsize=13, fontweight='bold')

plt.tight_layout()
plt.show()

# ============ FIGURE 2: Pocket Separation Distances ============
# Calculate distances between all pocket centers
pocket_centers = df_pockets[['center_x', 'center_y', 'center_z']].values
pocket_types_arr = df_pockets['pocket_type'].values
pocket_pairs_arr = df_pockets['pair'].values

if len(pocket_centers) > 1:
    distance_matrix = squareform(pdist(pocket_centers))
    
    # Distances within same pair vs across pairs
    within_pair_distances = []
    across_pair_distances = []
    
    for i in range(len(pocket_centers)):
        for j in range(i+1, len(pocket_centers)):
            dist = distance_matrix[i, j]
            if pocket_pairs_arr[i] == pocket_pairs_arr[j]:
                within_pair_distances.append({
                    'distance': dist,
                    'type1': pocket_types_arr[i],
                    'type2': pocket_types_arr[j],
                    'pair': pocket_pairs_arr[i]
                })
            else:
                across_pair_distances.append({
                    'distance': dist,
                    'type1': pocket_types_arr[i],
                    'type2': pocket_types_arr[j]
                })
    
    df_within = pd.DataFrame(within_pair_distances)
    df_across = pd.DataFrame(across_pair_distances)
    
    fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))
    
    # 2a: Within-pair pocket distances
    ax2a = axes2[0]
    if len(df_within) > 0:
        ax2a.hist(df_within['distance'], bins=20, color='#3498db', edgecolor='black', alpha=0.7)
        ax2a.axvline(x=POCKET_THRESHOLD, color='red', linestyle='--', linewidth=2,
                     label=f'Clustering threshold ({POCKET_THRESHOLD} Å)')
        ax2a.axvline(x=df_within['distance'].mean(), color='green', linestyle='-', linewidth=2,
                     label=f'Mean ({df_within["distance"].mean():.1f} Å)')
    ax2a.set_xlabel('Distance (Å)', fontsize=12)
    ax2a.set_ylabel('Frequency', fontsize=12)
    ax2a.set_title('Pocket Separation (Same Protein-Ligand Pair)', fontsize=12, fontweight='bold')
    ax2a.legend(fontsize=9)
    ax2a.grid(axis='y', alpha=0.3)
    
    # 2b: Comparison by pocket type combination
    ax2b = axes2[1]
    if len(df_within) > 0:
        df_within['type_combo'] = df_within.apply(
            lambda x: f"{min(x['type1'], x['type2'])} - {max(x['type1'], x['type2'])}", axis=1)
        type_combos = df_within['type_combo'].unique()
        distances_by_combo = [df_within[df_within['type_combo'] == tc]['distance'].values for tc in type_combos]
        
        bp = ax2b.boxplot(distances_by_combo, labels=[tc.replace(' - ', '\n-\n') for tc in type_combos],
                          patch_artist=True)
        colors_box = plt.cm.Set2(np.linspace(0, 1, len(bp['boxes'])))
        for patch, color in zip(bp['boxes'], colors_box):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax2b.axhline(y=POCKET_THRESHOLD, color='red', linestyle='--', linewidth=2)
    ax2b.set_ylabel('Distance (Å)', fontsize=12)
    ax2b.set_title('Pocket Distance by Type Combination', fontsize=12, fontweight='bold')
    ax2b.grid(axis='y', alpha=0.3)
    
    # 2c: 3D scatter of pocket centers (if matplotlib supports it)
    ax2c = axes2[2]
    ax2c = fig2.add_subplot(1, 3, 3, projection='3d')
    
    colors_3d = {'Shared': '#9b59b6', 'MGLTools Only': '#2ecc71', 'Meeko Only': '#3498db'}
    for pocket_type in ['Shared', 'MGLTools Only', 'Meeko Only']:
        mask = df_pockets['pocket_type'] == pocket_type
        ax2c.scatter(df_pockets.loc[mask, 'center_x'],
                     df_pockets.loc[mask, 'center_y'],
                     df_pockets.loc[mask, 'center_z'],
                     c=colors_3d[pocket_type], label=pocket_type, s=50, alpha=0.7)
    
    ax2c.set_xlabel('X (Å)')
    ax2c.set_ylabel('Y (Å)')
    ax2c.set_zlabel('Z (Å)')
    ax2c.set_title('Pocket Centers in 3D Space', fontsize=12, fontweight='bold')
    ax2c.legend(fontsize=9)
    
    plt.tight_layout()
    plt.show()

# ============ FIGURE 3: Per-Pair Cluster Analysis ============
fig3, axes3 = plt.subplots(2, 2, figsize=(16, 12))

# 3a: Poses per pocket per pair (stacked by type)
ax3a = axes3[0, 0]
pairs = df_summaries['pair'].values
x = np.arange(len(pairs))
width = 0.6

# Stack: shared MGLTools, shared Meeko, unique MGLTools, unique Meeko
ax3a.bar(x, df_summaries['shared_poses_mgl'], width, label='MGLTools (Shared)', color='#27ae60', edgecolor='black')
ax3a.bar(x, df_summaries['shared_poses_meeko'], width, bottom=df_summaries['shared_poses_mgl'],
         label='Meeko (Shared)', color='#2980b9', edgecolor='black')
ax3a.bar(x, df_summaries['unique_poses_mgl'], width, 
         bottom=df_summaries['shared_poses_mgl'] + df_summaries['shared_poses_meeko'],
         label='MGLTools (Unique)', color='#16a085', edgecolor='black')
ax3a.bar(x, df_summaries['unique_poses_meeko'], width,
         bottom=df_summaries['shared_poses_mgl'] + df_summaries['shared_poses_meeko'] + df_summaries['unique_poses_mgl'],
         label='Meeko (Unique)', color='#1abc9c', edgecolor='black')

ax3a.set_xlabel('Protein-Ligand Pair', fontsize=11)
ax3a.set_ylabel('Number of Poses', fontsize=11)
ax3a.set_title('Pose Distribution per Pair (Shared vs Unique)', fontsize=12, fontweight='bold')
ax3a.set_xticks(x)
ax3a.set_xticklabels(pairs, rotation=45, ha='right', fontsize=8)
ax3a.legend(fontsize=9, loc='upper right')
ax3a.grid(axis='y', alpha=0.3)

# 3b: Number of pockets per pair by type
ax3b = axes3[0, 1]
width = 0.25
ax3b.bar(x - width, df_summaries['mgltools_only_pockets'], width, label='MGLTools Only', color='#2ecc71', edgecolor='black')
ax3b.bar(x, df_summaries['shared_pockets'], width, label='Shared', color='#9b59b6', edgecolor='black')
ax3b.bar(x + width, df_summaries['meeko_only_pockets'], width, label='Meeko Only', color='#3498db', edgecolor='black')

ax3b.set_xlabel('Protein-Ligand Pair', fontsize=11)
ax3b.set_ylabel('Number of Pockets', fontsize=11)
ax3b.set_title('Pockets per Pair by Type', fontsize=12, fontweight='bold')
ax3b.set_xticks(x)
ax3b.set_xticklabels(pairs, rotation=45, ha='right', fontsize=8)
ax3b.legend(fontsize=10)
ax3b.grid(axis='y', alpha=0.3)

# 3c: Ratio heatmap - shared vs unique poses per pair
ax3c = axes3[1, 0]
ratio_data = []
for _, row in df_summaries.iterrows():
    total_mgl = row['shared_poses_mgl'] + row['unique_poses_mgl']
    total_meeko = row['shared_poses_meeko'] + row['unique_poses_meeko']
    ratio_data.append({
        'pair': row['pair'],
        'mgl_shared_pct': 100 * row['shared_poses_mgl'] / max(1, total_mgl),
        'meeko_shared_pct': 100 * row['shared_poses_meeko'] / max(1, total_meeko),
        'mgl_unique_pct': 100 * row['unique_poses_mgl'] / max(1, total_mgl),
        'meeko_unique_pct': 100 * row['unique_poses_meeko'] / max(1, total_meeko)
    })

df_ratios = pd.DataFrame(ratio_data)
heatmap_data = df_ratios[['mgl_shared_pct', 'mgl_unique_pct', 'meeko_shared_pct', 'meeko_unique_pct']].values

im = ax3c.imshow(heatmap_data, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
ax3c.set_xticks(np.arange(4))
ax3c.set_xticklabels(['MGLTools\nShared', 'MGLTools\nUnique', 'Meeko\nShared', 'Meeko\nUnique'], fontsize=9)
ax3c.set_yticks(np.arange(len(df_ratios)))
ax3c.set_yticklabels(df_ratios['pair'], fontsize=8)

for i in range(len(df_ratios)):
    for j in range(4):
        val = heatmap_data[i, j]
        ax3c.text(j, i, f'{val:.0f}%', ha='center', va='center',
                  color='white' if val > 60 or val < 40 else 'black', fontsize=9, fontweight='bold')

cbar = plt.colorbar(im, ax=ax3c)
cbar.set_label('% of Poses', fontsize=10)
ax3c.set_title('Pose % in Shared vs Unique Pockets by Pair', fontsize=12, fontweight='bold')

# 3d: Summary table visualization
ax3d = axes3[1, 1]
ax3d.axis('off')

summary_text = f"""
SUMMARY: POSES IN SHARED vs UNIQUE POCKETS

┌─────────────────┬────────────────┬────────────────┐
│                 │   MGLTools     │     Meeko      │
├─────────────────┼────────────────┼────────────────┤
│ Total Poses     │ {total_mgl_poses:>10}     │ {total_meeko_poses:>10}     │
│ In Shared       │ {mgl_poses_in_shared:>10}     │ {meeko_poses_in_shared:>10}     │
│ In Unique       │ {mgl_poses_in_unique:>10}     │ {meeko_poses_in_unique:>10}     │
│ % Shared        │ {100*mgl_poses_in_shared/max(1,total_mgl_poses):>9.1f}%    │ {100*meeko_poses_in_shared/max(1,total_meeko_poses):>9.1f}%    │
│ % Unique        │ {100*mgl_poses_in_unique/max(1,total_mgl_poses):>9.1f}%    │ {100*meeko_poses_in_unique/max(1,total_meeko_poses):>9.1f}%    │
│ Ratio Sh/Uniq   │ {mgl_poses_in_shared/max(1,mgl_poses_in_unique):>10.2f}    │ {meeko_poses_in_shared/max(1,meeko_poses_in_unique):>10.2f}    │
└─────────────────┴────────────────┴────────────────┘

POCKET STATISTICS:
• Total Pockets: {total_pockets}
• Shared: {total_shared_pockets} ({100*total_shared_pockets/total_pockets:.1f}%)
• MGLTools Only: {total_mgl_only_pockets} ({100*total_mgl_only_pockets/total_pockets:.1f}%)
• Meeko Only: {total_meeko_only_pockets} ({100*total_meeko_only_pockets/total_pockets:.1f}%)
"""

ax3d.text(0.1, 0.95, summary_text, transform=ax3d.transAxes, fontsize=11,
          verticalalignment='top', fontfamily='monospace',
          bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.show()

# ============ FIGURE 4: Detailed Pose Count per Cluster ============
fig4, ax4 = plt.subplots(figsize=(16, 6))

# Sort pockets by total poses
df_pockets_sorted = df_pockets.sort_values('total_poses', ascending=False).reset_index(drop=True)

x = np.arange(len(df_pockets_sorted))
width = 0.35

# Color by pocket type
type_colors = {'Shared': '#9b59b6', 'MGLTools Only': '#2ecc71', 'Meeko Only': '#3498db'}
bar_colors_mgl = [type_colors[t] for t in df_pockets_sorted['pocket_type']]
bar_colors_meeko = [type_colors[t] for t in df_pockets_sorted['pocket_type']]

bars_mgl = ax4.bar(x - width/2, df_pockets_sorted['mgltools_poses'], width, 
                   label='MGLTools Poses', color='#2ecc71', edgecolor='black', alpha=0.8)
bars_meeko = ax4.bar(x + width/2, df_pockets_sorted['meeko_poses'], width,
                     label='Meeko Poses', color='#3498db', edgecolor='black', alpha=0.8)

# Mark pocket types with different edge colors
for i, (bar_m, bar_e, ptype) in enumerate(zip(bars_mgl, bars_meeko, df_pockets_sorted['pocket_type'])):
    if ptype == 'Shared':
        bar_m.set_edgecolor('#9b59b6')
        bar_e.set_edgecolor('#9b59b6')
        bar_m.set_linewidth(2)
        bar_e.set_linewidth(2)

ax4.set_xlabel('Pocket (sorted by total poses)', fontsize=12)
ax4.set_ylabel('Number of Poses', fontsize=12)
ax4.set_title('Poses per Cluster/Pocket (All Pairs Combined)', fontsize=14, fontweight='bold')
ax4.set_xticks(x[::max(1, len(x)//20)])  # Show every nth tick
ax4.set_xticklabels([f"P{i+1}" for i in x[::max(1, len(x)//20)]], fontsize=8)
ax4.legend(fontsize=11)
ax4.grid(axis='y', alpha=0.3)

# Add legend for pocket types
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#9b59b6', edgecolor='black', label='Shared Pocket'),
                   Patch(facecolor='#2ecc71', edgecolor='black', label='MGLTools Only'),
                   Patch(facecolor='#3498db', edgecolor='black', label='Meeko Only')]
ax4.legend(handles=legend_elements + ax4.get_legend_handles_labels()[0][:2], 
           loc='upper right', fontsize=9)

plt.tight_layout()
plt.show()

# ============ DETAILED TABLES ============
print("\n" + "=" * 100)
print("DETAILED POCKET DATA")
print("=" * 100)
display(df_pockets[['pair', 'pocket_id', 'pocket_type', 'mgltools_poses', 'meeko_poses', 'total_poses']])

print("\n" + "=" * 100)
print("PER-PAIR SUMMARY")
print("=" * 100)
display(df_summaries)

print("\n" + "=" * 100)
print("POSE RATIOS (Shared/Unique) BY PAIR")
print("=" * 100)
display(df_ratios.round(1))


