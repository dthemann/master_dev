#!/usr/bin/env python3
"""
DiffDock batch docking CLI.

Usage:
    python run_diffdock.py                           # uses default config next to this script
    python run_diffdock.py -c my_config.yaml         # custom config
    python run_diffdock.py --dry-run                  # validate setup without docking
    python run_diffdock.py --optimize smina           # re-optimise poses (none|smina|gnina|all)

Optionally re-optimises each DiffDock pose against the receptor with smina
and/or gnina (local minimisation, same pocket); set via the 'optimization'
config key or the --optimize flag. Optimised poses and scores are written
alongside the originals (optimized_<tool>/ + optimization_log.csv).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
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
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

_SCRIPT_DIR = Path(__file__).resolve().parent


# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

def load_config(config_path: Path) -> dict:
    """Load and validate a YAML configuration file."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    required = ["receptors_dir", "ligand_dirs", "output_dir", "diffdock_dir", "diffdock_python"]
    for key in required:
        if key not in cfg or cfg[key] is None:
            raise ValueError(f"Missing required config key: '{key}'")

    # Normalise device
    device = cfg.get("device", "cuda:0").lower().strip()
    cfg["device"] = device

    return cfg


def print_config(cfg: dict) -> None:
    """Print a formatted summary of the configuration."""
    print("=" * 80)
    print("DiffDock Batch Docking Configuration")
    print("=" * 80)
    print(f"DiffDock dir:        {cfg['diffdock_dir']}")
    print(f"DiffDock python:     {cfg['diffdock_python']}")
    print(f"Device:              {cfg['device']}")
    print(f"Samples per complex: {cfg.get('num_samples', 10)}")
    print(f"Inference steps:     {cfg.get('inference_steps', 20)}")

    opt = str(cfg.get("optimization", "none")).strip().lower()
    print(f"Pose optimization:   {opt}")
    if opt not in ("none", "off", ""):
        print(f"  search mode:       {cfg.get('optimize_search', 'minimize')}")
        print(f"  autobox add:       {cfg.get('optimize_autobox_add', 4.0)} Å")
        if opt in ("gnina", "all"):
            print(f"  gnina GPU:         {'on' if cfg.get('gnina_use_gpu', False) else 'off (--no_gpu)'}")
    tpc = cfg.get("timeout_per_complex", 0)
    print(f"Timeout per complex: {tpc}s ({tpc / 60:.0f} min)" if tpc else "Timeout per complex: unlimited")

    batch_mode = cfg.get("batch_mode", True)
    batch_size = cfg.get("batch_size", 15)
    print(f"Batch mode:          {'ON' if batch_mode else 'OFF'}")
    print(f"Batch size:          {batch_size or 'all at once'}")
    if batch_mode:
        print(f"Group by receptor:   {'YES' if cfg.get('group_by_receptor', False) else 'no'}")
    if not batch_mode:
        print(f"Workers:             {cfg.get('max_workers', 1)}")

    print(f"Overwrite existing:  {cfg.get('overwrite_existing', False)}")
    print(f"Overwrite error log: {cfg.get('overwrite_error_log', True)}")
    print(f"Output dir:          {cfg['output_dir']}")
    print(f"Log dir:             {cfg.get('log_dir', 'logs')}")
    print(f"Receptor dir:        {cfg['receptors_dir']}")
    print(f"Ligand dirs ({len(cfg['ligand_dirs'])}):")
    for d in cfg["ligand_dirs"]:
        print(f"  • {d}")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# Data classes & helpers
# ═══════════════════════════════════════════════════════════════════════════════

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
    return path.stem.replace("_ligand", "").replace("_protein", "")


def collect_files(root: Path, extensions: List[str]) -> List[Path]:
    files = []
    for ext in extensions:
        files.extend(sorted(p for p in root.glob(f"*{ext}") if p.is_file()))
    return sorted(set(files))


def extract_confidence_from_filename(filename: str) -> float:
    # DiffDock encodes the (usually negative) confidence directly in the name,
    # e.g. rank1_confidence-1.65.sdf → -1.65. Keep the sign.
    m = re.search(r"confidence(-?\d+(?:\.\d+)?)", filename.lower())
    return float(m.group(1)) if m else 0.0


def extract_rank_from_filename(filename: str) -> int:
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
    except Exception:
        pass
    return 0


# ═══════════════════════════════════════════════════════════════════════════════
# Protein & ligand preparation
# ═══════════════════════════════════════════════════════════════════════════════

_RESIDUE_MAPPING = {
    'HSD': 'HIS', 'HSE': 'HIS', 'HSP': 'HIS',
    'HIE': 'HIS', 'HID': 'HIS', 'HIP': 'HIS',
}


def prepare_protein_for_diffdock(input_pdb: Path, output_dir: Path) -> Path:
    """Prepare a protein PDB for DiffDock (remap non-standard residues)."""
    if input_pdb.is_dir():
        pdb_files = sorted(input_pdb.glob("*.pdb"))
        if not pdb_files:
            raise FileNotFoundError(f"No .pdb files in directory: {input_pdb}")
        for pdb in pdb_files:
            _prepare_single_protein(pdb, output_dir)
        return output_dir
    return _prepare_single_protein(input_pdb, output_dir)


def _prepare_single_protein(input_pdb: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    output_pdb = output_dir / f"{input_pdb.stem}_prepared.pdb"

    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM'):
                res_name = line[17:20].strip()
                if res_name in _RESIDUE_MAPPING:
                    new_res = _RESIDUE_MAPPING[res_name]
                    line = line[:17] + f"{new_res:>3}" + line[20:]
                f_out.write(line)
            elif line.startswith(('END', 'TER')):
                f_out.write(line)

    return output_pdb


def prepare_ligand_for_diffdock(input_sdf: Path, output_dir: Path) -> Path:
    """Prepare a ligand SDF for DiffDock by fixing valence/charge issues."""
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
            # Only infer a charge when the SDF does NOT already declare an
            # explicit valence (ev == 0). If ev > 0, the writer already told
            # the parser the intended valence (e.g. P with valence 5 in a
            # phosphonate) — adding a charge would corrupt the molecule.
            if bc > dv and ev == 0:
                charge = bc - dv
                charges_to_add[idx + 1] = charge

        if not charges_to_add:
            return input_sdf

        # MDL V2000 M  CHG line uses fixed 4-char columns:
        # "M  CHGnnn aaa ccc aaa ccc..." where nnn=count(>3), each atom/charge=>4
        charge_atoms = list(charges_to_add.items())
        chg_line = f"M  CHG{len(charge_atoms):>3}"
        for atom_idx, charge in charge_atoms:
            chg_line += f"{atom_idx:>4}{charge:>4}"
        chg_line += "\n"

        with open(output_sdf, 'w') as f:
            for line in lines:
                if line.strip() == "M  END":
                    f.write(chg_line)
                f.write(line)

        return output_sdf

    except Exception:
        return input_sdf


# ═══════════════════════════════════════════════════════════════════════════════
# Error log
# ═══════════════════════════════════════════════════════════════════════════════

ERROR_LOG_FILENAME = "failed_docking.json"


def load_error_log(output_dir: Path) -> Dict[str, dict]:
    error_log_path = output_dir / ERROR_LOG_FILENAME
    if error_log_path.exists():
        try:
            with open(error_log_path, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, Exception):
            pass
    return {}


def save_error_log(output_dir: Path, error_log: Dict[str, dict]):
    error_log_path = output_dir / ERROR_LOG_FILENAME
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(error_log_path, 'w') as f:
        json.dump(error_log, f, indent=2)


def record_failure(output_dir: Path, combo_name: str, result: DockingResult):
    error_log = load_error_log(output_dir)
    error_log[combo_name] = {
        "error": result.error_message,
        "timestamp": datetime.now().isoformat(),
        "protein": result.protein_name,
        "ligand": result.ligand_name,
        "elapsed_time": round(result.elapsed_time, 2),
    }
    save_error_log(output_dir, error_log)


# ═══════════════════════════════════════════════════════════════════════════════
# Docking log & molecular properties
# ═══════════════════════════════════════════════════════════════════════════════

DOCKING_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name",
    "status", "error_reason", "num_poses", "elapsed_time_s",
    "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula",
    "prot_num_residues", "prot_num_atoms", "prot_num_chains",
    "num_samples", "inference_steps", "device", "gpu_model",
    "best_confidence",
]


def get_gpu_model() -> str:
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
    except Exception:
        pass
    return props


def get_protein_properties(pdb_path: Path) -> Dict:
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
    except Exception:
        pass
    return props


def ligand_props_valid(props: Dict) -> bool:
    return props.get("lig_heavy_atoms") is not None


def precompute_properties(
    proteins: List[Path],
    ligands: List[Path],
) -> Tuple[Dict[Path, Dict], Dict[Path, Dict]]:
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


def _build_log_row(
    result: DockingResult,
    lig_props: Dict,
    prot_props: Dict,
    cfg: dict,
    gpu_model: str,
) -> Dict:
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
        "num_samples": cfg.get("num_samples", 10),
        "inference_steps": cfg.get("inference_steps", 20),
        "device": cfg.get("device", "cuda:0"),
        "gpu_model": gpu_model,
    }
    row.update(lig_props)
    row.update(prot_props)
    return row


def init_docking_log(output_dir: Path) -> Tuple[Path, set]:
    log_path = output_dir / "docking_log.csv"
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
    cfg: dict,
    gpu_model: str,
    overwrite: bool = False,
) -> None:
    rows = []
    replaced = 0
    for r in results:
        combo = f"{r.ligand_name}__{r.protein_name}"
        if combo in existing_combos and not overwrite:
            continue
        lig_props = lig_props_cache.get(r.ligand_path, get_ligand_properties(r.ligand_path))
        prot_props = prot_props_cache.get(r.protein_path, get_protein_properties(r.protein_path))
        rows.append(_build_log_row(r, lig_props, prot_props, cfg, gpu_model))
        if combo in existing_combos:
            replaced += 1
        existing_combos.add(combo)

    if not rows:
        return

    if replaced > 0 and log_path.exists():
        df_existing = pd.read_csv(log_path)
        replace_combos = {row["combo_name"] for row in rows}
        df_existing = df_existing[~df_existing["combo_name"].isin(replace_combos)]
        df_new = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
        frames = [df for df in (df_existing, df_new) if not df.empty]
        df_all = pd.concat(frames, ignore_index=True) if frames else df_new
        df_all.to_csv(log_path, index=False)
        print(f"  📝 Replaced {replaced} + appended {len(rows) - replaced} entries in {log_path.name}")
    else:
        df = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
        df.to_csv(log_path, mode='a', header=False, index=False)
        print(f"  📝 Appended {len(rows)} entries to {log_path.name}")


# ═══════════════════════════════════════════════════════════════════════════════
# Summary & reporting
# ═══════════════════════════════════════════════════════════════════════════════

def generate_summary(results: List[DockingResult], cfg: dict) -> Dict:
    summary = {
        "timestamp": datetime.now().isoformat(),
        "configuration": {
            "num_samples": cfg.get("num_samples", 10),
            "inference_steps": cfg.get("inference_steps", 20),
            "device": cfg.get("device", "cuda:0"),
            "diffdock_dir": str(cfg["diffdock_dir"]),
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


# ═══════════════════════════════════════════════════════════════════════════════
# DiffDock inference engine
# ═══════════════════════════════════════════════════════════════════════════════

def _build_diffdock_config(
    samples: int,
    steps: int,
    output_dir: Path,
    diffdock_dir: Path,
) -> Path:
    default_yaml = diffdock_dir / "default_inference_args.yaml"
    with open(default_yaml) as f:
        config = yaml.safe_load(f)

    config["samples_per_complex"] = samples
    config["inference_steps"] = steps
    config["actual_steps"] = steps - 1

    output_dir.mkdir(parents=True, exist_ok=True)
    custom_yaml = (output_dir / "custom_inference_args.yaml").resolve()
    with open(custom_yaml, "w") as f:
        yaml.dump(config, f, default_flow_style=False)

    return custom_yaml


def _collect_ranked_poses(output_dir: Path) -> List[Path]:
    # rglob (not glob) because single-complex mode nests poses under a
    # DiffDock-named subdir (complex_0/). Exclude optimized_<tool>/ subfolders,
    # whose re-scored poses (rank1_confidence-1.65_smina.sdf) also match the
    # pattern — otherwise a re-run would double-count and re-optimize them.
    output_sdfs = [
        s for s in output_dir.rglob("rank*_confidence*.sdf")
        if s.stat().st_size > 0 and not s.parent.name.startswith("optimized_")
    ]
    return sorted(output_sdfs, key=lambda p: extract_rank_from_filename(p.name))


def _write_batch_csv(
    combos: List[Tuple[str, Path, Path]],
    csv_path: Path,
) -> Path:
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
    diffdock_python: str,
    diffdock_dir: Path,
    inference_seed: Optional[int] = None,
) -> subprocess.CompletedProcess:
    cmd = [
        diffdock_python, "-m", "inference",
        "--config", str(Path(config_yaml).resolve()),
        "--protein_ligand_csv", str(Path(csv_path).resolve()),
        "--out_dir", str(Path(out_dir).resolve()),
        "--samples_per_complex", str(samples),
        "--no_final_step_noise",
    ]

    env = os.environ.copy()
    if device == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""
    # DiffDock ships no seed option, so its sampling differs on every run and a
    # re-run silently re-samples the whole dataset. inference.py reads this and
    # seeds python/numpy/torch when it is set.
    if inference_seed is not None:
        env["DIFFDOCK_SEED"] = str(int(inference_seed))

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
            cwd=str(diffdock_dir),
            env=env,
        )

    result.stdout = stdout_log.read_text() if stdout_log.stat().st_size < 5_000_000 else f"[see {stdout_log}]"
    result.stderr = stderr_log.read_text() if stderr_log.stat().st_size < 5_000_000 else f"[see {stderr_log}]"

    return result


def _run_diffdock_subprocess(
    protein_path: Path,
    ligand_path: Path,
    output_dir: Path,
    output_base_dir: Path,
    samples: int,
    steps: int,
    device: str,
    timeout: int,
    diffdock_python: str,
    diffdock_dir: Path,
) -> Tuple[subprocess.CompletedProcess, List[Path]]:
    config_yaml = _build_diffdock_config(samples, steps, output_base_dir, diffdock_dir)

    cmd = [
        diffdock_python, "-m", "inference",
        "--config", str(config_yaml.resolve()),
        "--protein_path", str(protein_path.resolve()),
        "--ligand_description", str(ligand_path.resolve()),
        "--out_dir", str(output_dir.resolve()),
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
        cwd=str(diffdock_dir),
        env=env,
    )

    output_sdfs = _collect_ranked_poses(output_dir)
    return result, output_sdfs


def run_diffdock_single(
    protein_path: Path,
    ligand_path: Path,
    output_dir: Path,
    output_base_dir: Path,
    cfg: dict,
) -> Tuple[bool, List[Path], str]:
    """Run DiffDock for a single protein-ligand pair (with CUDA OOM fallback)."""
    samples = cfg.get("num_samples", 10)
    steps = cfg.get("inference_steps", 20)
    device = cfg.get("device", "cuda:0")
    timeout = cfg.get("timeout_per_complex", 300)
    diffdock_python = cfg["diffdock_python"]
    diffdock_dir = Path(cfg["diffdock_dir"])

    output_dir.mkdir(parents=True, exist_ok=True)

    try:
        result, output_sdfs = _run_diffdock_subprocess(
            protein_path, ligand_path, output_dir, output_base_dir,
            samples, steps, device, timeout,
            diffdock_python, diffdock_dir,
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
                diffdock_python, diffdock_dir,
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


# ═══════════════════════════════════════════════════════════════════════════════
# Batch docking pipeline
# ═══════════════════════════════════════════════════════════════════════════════

def _run_batch_csv_mode(
    combos_to_dock: List[Tuple[str, Path, Path, Path, Path]],
    output_dir: Path,
    prot_cache: Dict[Path, Path],
    lig_cache: Dict[Path, Path],
    cfg: dict,
    gpu_model: str,
    error_log: Dict[str, dict],
    log_path: Optional[Path] = None,
    existing_combos: Optional[set] = None,
    lig_props_cache: Optional[Dict[Path, Dict]] = None,
    prot_props_cache: Optional[Dict[Path, Dict]] = None,
) -> List[DockingResult]:
    results = []

    if not combos_to_dock:
        return results

    samples = cfg.get("num_samples", 10)
    steps = cfg.get("inference_steps", 20)
    device = cfg.get("device", "cuda:0")
    timeout = cfg.get("timeout_per_complex", 300)
    batch_size = cfg.get("batch_size", 15) or len(combos_to_dock)
    group_by_receptor = cfg.get("group_by_receptor", False)
    diffdock_python = cfg["diffdock_python"]
    diffdock_dir = Path(cfg["diffdock_dir"])

    if group_by_receptor:
        receptor_groups: Dict[Path, list] = defaultdict(list)
        for combo in combos_to_dock:
            receptor_groups[combo[3]].append(combo)  # combo[3] == orig_prot
        batches = []
        for group in receptor_groups.values():
            batches.extend(
                group[i:i + batch_size] for i in range(0, len(group), batch_size)
            )
    else:
        batches = [combos_to_dock[i:i + batch_size] for i in range(0, len(combos_to_dock), batch_size)]
    config_yaml = _build_diffdock_config(samples, steps, output_dir, diffdock_dir)

    combos_completed = 0
    total_combos = len(combos_to_dock)

    progress_file = output_dir / "_batch_logs" / "_progress.txt"
    progress_file.parent.mkdir(exist_ok=True)

    for batch_idx, batch in enumerate(batches):
        batch_start_idx = sum(len(batches[k]) for k in range(batch_idx))
        batch_end_idx = min(batch_start_idx + len(batch), total_combos)

        print(f"\n  ── Batch {batch_idx + 1}/{len(batches)} ({len(batch)} complexes) "
              f"[{batch_start_idx + 1}–{batch_end_idx} of {total_combos}] ──")

        for j, (combo_name, _, _, orig_prot, orig_lig) in enumerate(batch):
            print(f"    {batch_start_idx + j + 1:>4}/{total_combos}  "
                  f"{get_file_stem(orig_lig)} + {get_file_stem(orig_prot)}")

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
                diffdock_python=diffdock_python,
                diffdock_dir=diffdock_dir,
                inference_seed=cfg.get("inference_seed"),
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
                error_lines = [
                    l for l in stderr.splitlines()
                    if any(k in l.lower() for k in
                           ['error', 'exception', 'traceback', 'failed', 'cannot', 'no module'])
                ]
                for line in error_lines[:5]:
                    print(f"    {line.strip()[:200]}")

            # OOM retry logic
            if ("OutOfMemoryError" in stderr or "CUDA out of memory" in stderr) and device != "cpu":
                retry_size = max(1, len(batch) // 2)
                while retry_size >= 1:
                    sub_batches = [batch[i:i + retry_size] for i in range(0, len(batch), retry_size)]
                    print(f"  ⚠ CUDA OOM — retrying with sub-batch size {retry_size}...")
                    oom_again = False
                    for sb_idx, sub_batch in enumerate(sub_batches):
                        sub_csv_rows = [
                            (cn, prot_cache.get(op, pp), lig_cache.get(ol, lp))
                            for cn, pp, lp, op, ol in sub_batch
                        ]
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
                                diffdock_python=diffdock_python,
                                diffdock_dir=diffdock_dir,
                                inference_seed=cfg.get("inference_seed"),
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
                        print(f"  ⚠ CUDA OOM even at batch size 1 — falling back to CPU...")
                        for sb_idx, sub_batch in enumerate(sub_batches):
                            sub_csv_rows = [
                                (cn, prot_cache.get(op, pp), lig_cache.get(ol, lp))
                                for cn, pp, lp, op, ol in sub_batch
                            ]
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
                                    diffdock_python=diffdock_python,
                                    diffdock_dir=diffdock_dir,
                                    inference_seed=cfg.get("inference_seed"),
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

        results.extend(batch_results)
        if log_path is not None and existing_combos is not None:
            append_log_rows(
                batch_results, log_path, existing_combos,
                lig_props_cache or {}, prot_props_cache or {},
                cfg, gpu_model,
            )

    if progress_file.exists():
        progress_file.unlink()

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Post-pose optimization (smina / gnina)
# ═══════════════════════════════════════════════════════════════════════════════
#
# DiffDock predicts ligand coordinates directly; it never minimises against a
# force field or excluded-volume term. An optional local re-optimisation relaxes
# each pose against the receptor (boxed around the pose, so it stays in the same
# pocket) and attaches a physics/CNN score:
#
#   smina  → --minimize / --local_only, writes <minimizedAffinity> (kcal/mol).
#   gnina  → same CLI + a CNN rescorer, writes <CNNscore>/<CNNaffinity> too.
#
# Both share the AutoDock-Vina interface, so one wrapper drives either. The
# original DiffDock SDFs are never modified — optimised poses are written to a
# per-tool subfolder (optimized_smina/ , optimized_gnina/) and scored in
# optimization_log.csv.

OPTIMIZER_TOOLS = ("smina", "gnina")

OPTIMIZATION_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name", "tool",
    "pose_rank", "pose_file", "diffdock_confidence",
    "minimized_affinity", "cnn_score", "cnn_affinity",
    "optimized_file", "status", "message", "elapsed_time_s",
]

# SDF data tags written by smina/gnina, mapped to our log column names.
_SDF_TAG_RE = re.compile(r">\s*<([^>]+)>\s*\n([^\n]*)")
_SCORE_TAGS = {
    "minimizedAffinity": "minimized_affinity",
    "CNNscore": "cnn_score",
    "CNNaffinity": "cnn_affinity",
}


def resolve_optimizers(cfg: dict) -> List[str]:
    """Map the ``optimization`` config value to a list of tools to run."""
    mode = str(cfg.get("optimization", "none")).strip().lower()
    if mode in ("none", "off", "false", ""):
        return []
    if mode == "all":
        return list(OPTIMIZER_TOOLS)
    if mode in OPTIMIZER_TOOLS:
        return [mode]
    raise ValueError(
        f"optimization must be one of none/smina/gnina/all (got {mode!r})"
    )


def _resolve_optimizer_exe(tool: str, cfg: dict) -> str:
    """Locate a tool's executable: explicit config path → PATH → ''."""
    explicit = cfg.get(f"{tool}_executable")
    if explicit:
        path = os.path.expanduser(str(explicit))
        if Path(path).exists():
            return path
        return shutil.which(path) or ""
    return shutil.which(tool) or ""


def _gnina_subprocess_env(cfg: dict) -> dict:
    """Environment for the gnina subprocess, with CUDA/cuDNN libs on the path.

    The prebuilt gnina CUDA binary is dynamically linked against the CUDA
    runtime + cuDNN 9 shared objects (libcudart.so.12, libcudnn.so.9, libcublas,
    …). Those are not installed system-wide but ARE bundled inside the diffdock
    env's torch install, so we add that env's ``nvidia/*/lib`` dirs to
    LD_LIBRARY_PATH. gnina resolves these at load time regardless of --no_gpu,
    so this is needed for CPU runs too. Honours an explicit ``gnina_lib_dirs``
    config list (prepended first).
    """
    env = os.environ.copy()
    lib_dirs: List[str] = []

    for d in cfg.get("gnina_lib_dirs") or []:
        p = os.path.expanduser(str(d))
        if Path(p).is_dir():
            lib_dirs.append(p)

    diffdock_python = cfg.get("diffdock_python")
    if diffdock_python:
        env_root = Path(os.path.expanduser(diffdock_python)).resolve().parent.parent
        lib_dirs.extend(
            str(p) for p in sorted(env_root.glob("lib/python*/site-packages/nvidia/*/lib"))
            if p.is_dir()
        )

    if lib_dirs:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(lib_dirs + ([existing] if existing else []))
    return env


def _probe_optimizer(tool: str, exe: str, cfg: dict) -> Tuple[bool, str]:
    """Run ``<exe> --version`` to confirm the binary actually loads.

    Catches missing shared libraries (e.g. libcudnn.so.9 for the gnina CUDA
    build) before a real docking run — used by --dry-run.
    """
    env = _gnina_subprocess_env(cfg) if tool == "gnina" else None
    try:
        proc = subprocess.run([exe, "--version"], capture_output=True, text=True,
                              timeout=30, env=env)
    except Exception as e:
        return False, str(e)
    lines = (proc.stdout or proc.stderr or "").strip().splitlines()
    if proc.returncode != 0:
        return False, lines[-1] if lines else f"rc={proc.returncode}"
    return True, lines[0] if lines else "ok"


def _parse_sdf_scores(sdf_path: Path) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    try:
        text = Path(sdf_path).read_text()
    except Exception:
        return scores
    for m in _SDF_TAG_RE.finditer(text):
        tag = m.group(1).strip()
        col = _SCORE_TAGS.get(tag)
        if col and col not in scores:
            try:
                scores[col] = float(m.group(2).strip())
            except ValueError:
                pass
    return scores


def run_optimizer_tool(
    tool: str,
    exe: str,
    receptor: Path,
    pose_sdf: Path,
    out_sdf: Path,
    cfg: dict,
) -> Tuple[bool, Dict[str, float], str]:
    """Locally re-optimise a single pose with smina or gnina.

    Returns ``(success, scores, message)``. On failure the caller keeps the
    original DiffDock pose, so a pose is never lost.
    """
    search = str(cfg.get("optimize_search", "minimize")).strip().lower()
    autobox_add = cfg.get("optimize_autobox_add", 4.0)
    cpu = cfg.get("optimize_cpu", 1)
    seed = cfg.get("optimize_seed", 0)
    timeout_s = cfg.get("optimize_timeout", 300)

    out_sdf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(pose_sdf),
        "--autobox_ligand", str(pose_sdf),
        "--autobox_add", str(autobox_add),
        "--out", str(out_sdf),
        "--cpu", str(cpu),
        "--seed", str(seed),
        "--num_modes", "1",
    ]
    # --local_only and --minimize both keep the ligand near its input pose and
    # never run a global search, so the optimised pose stays in the same pocket.
    cmd.append("--minimize" if search == "minimize" else "--local_only")
    # gnina runs a CNN by default; keep it off the GPU unless asked (the GPU is
    # usually busy with DiffDock). smina has no --no_gpu flag.
    if tool == "gnina" and not cfg.get("gnina_use_gpu", False):
        cmd.append("--no_gpu")

    # gnina needs the diffdock env's bundled CUDA/cuDNN libs on LD_LIBRARY_PATH.
    proc_env = _gnina_subprocess_env(cfg) if tool == "gnina" else None

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s, env=proc_env)
    except subprocess.TimeoutExpired:
        return False, {}, f"{tool} timed out after {timeout_s}s"
    except FileNotFoundError:
        return False, {}, f"{tool} executable not found: {exe}"
    except Exception as e:  # pragma: no cover - defensive
        return False, {}, f"{tool} invocation failed: {e}"

    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        return False, {}, f"{tool} rc={proc.returncode}: {err[-1] if err else 'unknown error'}"

    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        return False, {}, f"{tool} produced no output pose"

    return True, _parse_sdf_scores(out_sdf), f"{tool} {search} ok"


def _load_optimization_keys(log_path: Path) -> set:
    keys: set = set()
    if log_path.exists():
        try:
            df = pd.read_csv(log_path)
            for _, row in df.iterrows():
                keys.add((row["combo_name"], row["tool"], row["pose_file"]))
        except Exception:
            pass
    return keys


def _write_optimization_log(log_path: Path, rows: List[dict], overwrite: bool) -> None:
    if not rows:
        return
    df_new = pd.DataFrame(rows, columns=OPTIMIZATION_LOG_COLUMNS)
    if log_path.exists():
        try:
            df_old = pd.read_csv(log_path)
        except Exception:
            df_old = pd.DataFrame(columns=OPTIMIZATION_LOG_COLUMNS)
        if overwrite and not df_old.empty:
            new_keys = set(zip(df_new["combo_name"], df_new["tool"], df_new["pose_file"]))
            mask = [
                (c, t, p) not in new_keys
                for c, t, p in zip(df_old["combo_name"], df_old["tool"], df_old["pose_file"])
            ]
            df_old = df_old[mask]
        frames = [df for df in (df_old, df_new) if not df.empty]
        df_all = pd.concat(frames, ignore_index=True) if frames else df_new
    else:
        df_all = df_new
    log_path.parent.mkdir(parents=True, exist_ok=True)
    df_all.to_csv(log_path, index=False)


def optimize_results(
    results: List[DockingResult],
    prot_cache: Dict[Path, Path],
    cfg: dict,
    output_dir: Path,
    overwrite: bool = False,
) -> None:
    """Run the configured optimizer(s) over every successful pose in-place.

    Operates on both freshly docked and previously docked (skipped) complexes,
    so enabling ``optimization`` and re-running optimises existing poses without
    re-docking. Idempotent: existing optimised SDFs are reused unless ``overwrite``.
    """
    tools = resolve_optimizers(cfg)
    if not tools:
        return

    exes: Dict[str, str] = {}
    for tool in tools:
        exe = _resolve_optimizer_exe(tool, cfg)
        if exe:
            exes[tool] = exe
        else:
            print(f"  ⚠ {tool} not found on PATH — set '{tool}_executable' in the "
                  f"config or install it into the env. Skipping {tool}.")
    tools = [t for t in tools if t in exes]
    if not tools:
        print("  ⚠ No optimizer executables available — skipping optimization.")
        return

    targets = [r for r in results if r.status in ("success", "skipped") and r.pose_files]

    print("\n" + "=" * 80)
    print(f"Pose optimization: {', '.join(tools)}")
    print("=" * 80)
    print(f"Optimizing {len(targets)} complexes × {len(tools)} tool(s) "
          f"(search={cfg.get('optimize_search', 'minimize')})")

    log_path = output_dir / "optimization_log.csv"
    seen = _load_optimization_keys(log_path)

    rows: List[dict] = []
    n_ok = n_fail = n_reused = 0

    for idx, r in enumerate(targets, 1):
        receptor = Path(prot_cache.get(r.protein_path, r.protein_path))
        combo_name = f"{r.ligand_name}__{r.protein_name}"

        for tool in tools:
            tool_dir = r.output_dir / f"optimized_{tool}"
            for pose in r.pose_files:
                out_sdf = tool_dir / f"{pose.stem}_{tool}.sdf"
                key = (combo_name, tool, pose.name)

                if key in seen and not overwrite:
                    continue

                reuse = out_sdf.exists() and out_sdf.stat().st_size > 0 and not overwrite
                if reuse:
                    ok, scores, msg = True, _parse_sdf_scores(out_sdf), "reused existing"
                    elapsed = 0.0
                    n_reused += 1
                else:
                    t0 = time.time()
                    ok, scores, msg = run_optimizer_tool(
                        tool, exes[tool], receptor, pose, out_sdf, cfg,
                    )
                    elapsed = time.time() - t0
                    if ok:
                        n_ok += 1
                    else:
                        n_fail += 1

                rows.append({
                    "timestamp": datetime.now().isoformat(),
                    "combo_name": combo_name,
                    "protein_name": r.protein_name,
                    "ligand_name": r.ligand_name,
                    "tool": tool,
                    "pose_rank": extract_rank_from_filename(pose.name),
                    "pose_file": pose.name,
                    "diffdock_confidence": round(extract_confidence_from_filename(pose.name), 4),
                    "minimized_affinity": scores.get("minimized_affinity"),
                    "cnn_score": scores.get("cnn_score"),
                    "cnn_affinity": scores.get("cnn_affinity"),
                    "optimized_file": str(out_sdf) if ok else "",
                    "status": "success" if ok else "failed",
                    "message": msg,
                    "elapsed_time_s": round(elapsed, 2),
                })
                seen.add(key)

        print(f"  [{idx}/{len(targets)}] {combo_name}: "
              f"{len(r.pose_files)} pose(s) × {len(tools)} tool(s)")

    _write_optimization_log(log_path, rows, overwrite=overwrite)
    print(f"\nOptimization complete: {n_ok} optimized, {n_reused} reused, {n_fail} failed")

    # Per-tool wall-clock time spent this run (reused poses contribute 0s and
    # are excluded from the pose count / average).
    time_by_tool: Dict[str, float] = defaultdict(float)
    runs_by_tool: Dict[str, int] = defaultdict(int)
    for row in rows:
        time_by_tool[row["tool"]] += row["elapsed_time_s"]
        if row["message"] != "reused existing":
            runs_by_tool[row["tool"]] += 1
    if rows:
        print("Time by tool (this run, excludes reused):")
        for tool in tools:
            secs = time_by_tool.get(tool, 0.0)
            n = runs_by_tool.get(tool, 0)
            avg = secs / n if n else 0.0
            print(f"  {tool}: {secs:.1f}s total ({secs / 60:.1f} min) "
                  f"over {n} pose(s) ({avg:.1f}s/pose)")

    print(f"Optimization log: {log_path}")


# ═══════════════════════════════════════════════════════════════════════════════
# Main docking orchestrator
# ═══════════════════════════════════════════════════════════════════════════════

def run_diffdock(
    proteins: List[Path],
    ligands: List[Path],
    output_dir: Path,
    cfg: dict,
    gpu_model: str,
) -> List[DockingResult]:
    """Run DiffDock for all protein-ligand combinations.

    Supports three modes (driven by cfg):
      - batch_mode=True:   CSV batch (fastest — one DiffDock call per batch)
      - max_workers>1:     parallel single-complex subprocesses
      - else:              sequential single-complex subprocesses
    """
    samples = cfg.get("num_samples", 10)
    steps = cfg.get("inference_steps", 20)
    device = cfg.get("device", "cuda:0")
    timeout = cfg.get("timeout_per_complex", 300)
    batch_mode = cfg.get("batch_mode", True)
    batch_size = cfg.get("batch_size", 15)
    max_workers = cfg.get("max_workers", 1)
    overwrite_existing = cfg.get("overwrite_existing", False)
    overwrite_error_log = cfg.get("overwrite_error_log", True)

    total = len(proteins) * len(ligands)
    results: List[DockingResult] = []

    error_log = load_error_log(output_dir)
    if overwrite_error_log:
        error_log = {}
        save_error_log(output_dir, error_log)
        print("Error log cleared (overwrite_error_log=true)")

    print("=" * 80)
    print("DiffDock Docking")
    print("=" * 80)
    if batch_mode:
        gbr = cfg.get("group_by_receptor", False)
        mode_str = "CSV batch (group-by-receptor)" if gbr else "CSV batch"
    else:
        mode_str = f"parallel ({max_workers} workers)" if max_workers > 1 else "sequential"
    print(f"Mode: {mode_str}")
    print(f"Output directory: {output_dir}")
    print(f"Proteins: {len(proteins)}  |  Ligands: {len(ligands)}  |  Total: {total}")
    print(f"Previously failed (will skip): {len(error_log)}")
    print(f"Samples: {samples}  |  Steps: {steps}  |  Timeout: {timeout}s")

    if error_log:
        print(f"\nPreviously failed (from {ERROR_LOG_FILENAME}):")
        for i, (name, info) in enumerate(sorted(error_log.items())):
            if i < 10:
                print(f"  ✗ {name}: {info.get('error', 'unknown')[:80]}")
        if len(error_log) > 10:
            print(f"  ... and {len(error_log) - 10} more")

    # Pre-compute molecular properties
    lig_props_cache, prot_props_cache = precompute_properties(proteins, ligands)

    # Exclude ligands RDKit cannot parse
    invalid_ligs = [lig for lig in ligands if not ligand_props_valid(lig_props_cache.get(lig, {}))]
    if invalid_ligs:
        print(f"\n  ⚠ Excluding {len(invalid_ligs)} unreadable ligand(s):")
        for lig in invalid_ligs:
            print(f"    ✗ {lig.name}")
        ligands = [lig for lig in ligands if ligand_props_valid(lig_props_cache.get(lig, {}))]
        total = len(proteins) * len(ligands)
        print(f"  Remaining: {len(ligands)} ligands × {len(proteins)} proteins = {total} combinations")

    # Initialise docking log
    log_path, existing_combos = init_docking_log(output_dir)
    print(f"Docking log: {log_path} ({len(existing_combos)} existing entries)")

    # Prepare all inputs
    print("\nPreparing inputs...")
    prot_cache, lig_cache = _prepare_all_inputs(proteins, ligands, output_dir)

    # Classify combos: skip vs dock
    combos_to_dock = []
    skipped_results = []

    for protein in proteins:
        for ligand in ligands:
            prot_stem = get_file_stem(protein)
            lig_stem = get_file_stem(ligand)
            combo_name = f"{lig_stem}__{prot_stem}"
            combo_dir = output_dir / combo_name

            if not overwrite_error_log and combo_name in error_log:
                skipped_results.append(DockingResult(
                    protein_name=prot_stem, ligand_name=lig_stem,
                    protein_path=protein, ligand_path=ligand,
                    output_dir=combo_dir, status="skipped",
                    error_message=f"Previously failed: {error_log[combo_name].get('error', 'unknown')[:100]}",
                ))
                continue

            if not overwrite_existing and combo_dir.exists():
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

    if skipped_results:
        append_log_rows(
            skipped_results, log_path, existing_combos,
            lig_props_cache, prot_props_cache,
            cfg, gpu_model,
        )

    if not combos_to_dock:
        print("Nothing to dock!")
        # Still run the configured post-pose optimisation over already-docked
        # poses (skipped results carry their pose_files). Without this, enabling
        # optimization for an already-docked complex would silently do nothing,
        # because the optimize_results() call below is only reached when there is
        # something new to dock.
        optimize_results(results, prot_cache, cfg, output_dir, overwrite=overwrite_existing)
        save_error_log(output_dir, error_log)
        return results

    start_all = time.time()

    if batch_mode:
        batch_results = _run_batch_csv_mode(
            combos_to_dock=combos_to_dock,
            output_dir=output_dir,
            prot_cache=prot_cache,
            lig_cache=lig_cache,
            cfg=cfg,
            gpu_model=gpu_model,
            error_log=error_log,
            log_path=log_path,
            existing_combos=existing_combos,
            lig_props_cache=lig_props_cache,
            prot_props_cache=prot_props_cache,
        )
        results.extend(batch_results)

    elif max_workers > 1:
        print(f"\nRunning {len(combos_to_dock)} complexes with {max_workers} parallel workers...")

        def _dock_one(args):
            combo_name, prep_prot, prep_lig, orig_prot, orig_lig = args
            combo_dir = output_dir / combo_name
            combo_dir.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            success, poses, err = run_diffdock_single(
                prep_prot, prep_lig, combo_dir, output_dir, cfg,
            )
            elapsed = time.time() - t0
            return combo_name, DockingResult(
                protein_name=get_file_stem(orig_prot),
                ligand_name=get_file_stem(orig_lig),
                protein_path=orig_prot, ligand_path=orig_lig,
                output_dir=combo_dir,
                status="success" if success else "failed",
                num_poses=len(poses), pose_files=poses,
                error_message=err, elapsed_time=elapsed,
            )

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
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
                print(f"  [{i}/{len(combos_to_dock)}] {icon} {combo_name}: "
                      f"{r.num_poses} poses ({r.elapsed_time:.1f}s)")
                if len(batch_buffer) >= batch_size or i == len(combos_to_dock):
                    append_log_rows(
                        batch_buffer, log_path, existing_combos,
                        lig_props_cache, prot_props_cache,
                        cfg, gpu_model,
                    )
                    batch_buffer = []

    else:
        for i, (combo_name, prep_prot, prep_lig, orig_prot, orig_lig) in enumerate(combos_to_dock, 1):
            print(f"\n[{i}/{len(combos_to_dock)}] {get_file_stem(orig_lig)} + {get_file_stem(orig_prot)}")

            combo_dir = output_dir / combo_name
            combo_dir.mkdir(parents=True, exist_ok=True)
            t0 = time.time()
            success, poses, err = run_diffdock_single(
                prep_prot, prep_lig, combo_dir, output_dir, cfg,
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

            if r.status == "failed":
                record_failure(output_dir, combo_name, r)
                error_log[combo_name] = {
                    "error": r.error_message,
                    "timestamp": datetime.now().isoformat(),
                    "protein": r.protein_name,
                    "ligand": r.ligand_name,
                    "elapsed_time": round(r.elapsed_time, 2),
                }

            results.append(r)
            icon = {"success": "✓", "failed": "✗", "skipped": "⊘"}.get(r.status, "?")
            print(f"  {icon} {r.status} | Poses: {r.num_poses} | Time: {elapsed:.1f}s")
            if r.error_message:
                print(f"    Error: {r.error_message[:200]}")

            append_log_rows(
                [r], log_path, existing_combos,
                lig_props_cache, prot_props_cache,
                cfg, gpu_model,
            )

    total_time = time.time() - start_all
    n_success = sum(1 for r in results if r.status == "success")
    n_failed = sum(1 for r in results if r.status == "failed")

    print(f"\nDocking complete in {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"  Success: {n_success}  |  Failed: {n_failed}  |  Skipped: {len(skipped_results)}")

    # Optional post-pose optimization (smina / gnina) over all poses produced.
    optimize_results(results, prot_cache, cfg, output_dir, overwrite=overwrite_existing)

    save_error_log(output_dir, error_log)
    print(f"Error log saved: {output_dir / ERROR_LOG_FILENAME} ({len(error_log)} entries)")
    print(f"Docking log: {log_path} ({len(existing_combos)} entries)")

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# Progress monitor
# ═══════════════════════════════════════════════════════════════════════════════

def _monitor_progress(
    output_dir: Path,
    stop_event: threading.Event,
    poll_interval: float = 10.0,
):
    log_dir = output_dir / "_batch_logs"
    progress_file = log_dir / "_progress.txt"
    start_time = time.time()

    while not stop_event.is_set():
        stop_event.wait(poll_interval)
        if stop_event.is_set():
            break

        elapsed = time.time() - start_time

        # Count completed poses
        rank1_files = list(output_dir.glob("*/rank1_*.sdf"))
        n_done = len(rank1_files)

        # Count empty dirs (in-progress or failed)
        all_combo_dirs = [
            d for d in output_dir.iterdir()
            if d.is_dir() and not d.name.startswith("_")
            and d.name not in ("prepared_proteins", "prepared_ligands")
        ]
        n_empty = sum(1 for d in all_combo_dirs if not list(d.glob("rank1_*.sdf")))

        # Read progress file
        current_batch_info = ""
        current_combos = ""
        if progress_file.exists():
            try:
                parts = progress_file.read_text().strip().split("|")
                if len(parts) >= 2:
                    current_batch_info = f"[{parts[0]} combos] {parts[1]}"
                if len(parts) >= 3:
                    current_combos = ", ".join(parts[2:6])
                    if len(parts) > 6:
                        current_combos += f" (+{len(parts) - 6} more)"
            except Exception:
                pass

        # GPU utilization
        try:
            gpu_info = subprocess.run(
                ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used,memory.total",
                 "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5,
            )
            gpu_line = gpu_info.stdout.strip()
        except Exception:
            gpu_line = "N/A"

        # DiffDock stdout tail
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


# ═══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="DiffDock batch docking pipeline",
    )
    parser.add_argument(
        "-c", "--config",
        default=str(_SCRIPT_DIR / "diffdock_docking_config.yaml"),
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate configuration and inputs without running docking",
    )
    parser.add_argument(
        "--optimize", choices=["none", "smina", "gnina", "all"], default=None,
        help="Post-pose optimization tool(s); overrides config 'optimization'",
    )
    args = parser.parse_args()

    # Load config
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)
    cfg = load_config(config_path)

    # CLI override for the optimization mode, then validate it up front.
    if args.optimize is not None:
        cfg["optimization"] = args.optimize
    try:
        resolve_optimizers(cfg)
    except ValueError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # Resolve paths relative to CWD
    receptors_dir = Path(cfg["receptors_dir"])
    output_base = Path(cfg["output_dir"])
    log_dir = Path(cfg.get("log_dir", "logs"))
    ligand_dirs = [Path(d) for d in cfg["ligand_dirs"]]
    diffdock_dir = Path(cfg["diffdock_dir"])
    diffdock_python = os.path.expanduser(cfg["diffdock_python"])
    cfg["diffdock_python"] = diffdock_python

    output_base.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Validate directories
    errors = []
    for d in ligand_dirs:
        if not d.exists():
            errors.append(f"Ligand directory missing: {d}")
    if not receptors_dir.exists():
        errors.append(f"Receptor directory missing: {receptors_dir}")
    if not diffdock_dir.exists():
        errors.append(f"DiffDock directory missing: {diffdock_dir}")
    if not Path(diffdock_python).exists():
        errors.append(f"DiffDock Python not found: {diffdock_python}")

    default_yaml = diffdock_dir / "default_inference_args.yaml"
    if diffdock_dir.exists() and not default_yaml.exists():
        errors.append(f"DiffDock config missing: {default_yaml}")

    if errors:
        print("\n⚠ ISSUES FOUND:")
        for e in errors:
            print(f"  ✗ {e}")
        sys.exit(1)

    print_config(cfg)

    gpu_model = get_gpu_model()
    print(f"GPU: {gpu_model}")

    if args.dry_run:
        # Validate inputs
        for d in ligand_dirs:
            proteins = collect_files(receptors_dir, [".pdb"])
            ligands = collect_files(d, [".sdf", ".mol2"])
            print(f"  {d.name}: {len(proteins)} proteins, {len(ligands)} ligands")
        # Check optimizer executables resolve AND actually load (probe --version)
        for tool in resolve_optimizers(cfg):
            exe = _resolve_optimizer_exe(tool, cfg)
            if not exe:
                print(f"  optimizer {tool}: NOT FOUND (install into the env or set {tool}_executable)")
                continue
            ok, msg = _probe_optimizer(tool, exe, cfg)
            print(f"  optimizer {tool}: {exe}")
            print(f"      {'✓ ' + msg if ok else '✗ FAILS TO RUN: ' + msg}")
        print("\n[DRY RUN] Configuration validated. Exiting before docking.")
        sys.exit(0)

    # Dock per ligand directory
    all_results: List[DockingResult] = []

    for drugs_dir in ligand_dirs:
        result_subfolder = drugs_dir.name
        diffdock_output_dir = output_base / result_subfolder
        diffdock_output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'=' * 80}")
        print(f"Processing ligand directory: {drugs_dir.name}")
        print(f"Output: {diffdock_output_dir}")
        print(f"{'=' * 80}")

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
            args=(diffdock_output_dir, stop_monitor, 15),
            daemon=True,
        )
        monitor_thread.start()

        try:
            results = run_diffdock(
                proteins=proteins,
                ligands=ligands,
                output_dir=diffdock_output_dir,
                cfg=cfg,
                gpu_model=gpu_model,
            )
        finally:
            stop_monitor.set()
            monitor_thread.join(timeout=5)

        all_results.extend(results)

        summary = generate_summary(results, cfg)
        print_summary(summary)

        summary_path = diffdock_output_dir / "docking_summary.json"
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2, default=str)
        print(f"\nSummary saved: {summary_path}")

    print(f"\n{'=' * 80}")
    print(f"ALL DONE — {len(all_results)} total results across {len(ligand_dirs)} ligand directories")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
