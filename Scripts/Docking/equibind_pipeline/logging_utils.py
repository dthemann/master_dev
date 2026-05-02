"""Docking log (CSV) and error log (JSON) helpers."""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

from .config import CFG
from .properties import get_gpu_model

DOCKING_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name",
    "status", "error_reason", "elapsed_time_s",
    "total_poses", "fpocket_poses", "p2rank_poses", "unguided_poses", "failed_poses",
    "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula",
    "prot_num_residues", "prot_num_atoms", "prot_num_chains",
    "n_top_pockets", "poses_per_pocket", "n_unguided_poses",
    "force_pocket", "uff_minimize", "device", "gpu_model",
]

# Detect GPU model once on import.
_GPU_MODEL = get_gpu_model()
print(f"GPU detected: {_GPU_MODEL}")


def _build_log_row(combo_name: str, protein_name: str, ligand_name: str,
                   status: str, error_reason: str, elapsed_time: float,
                   total_poses: int, fpocket_poses: int, p2rank_poses: int,
                   unguided_poses: int, failed_poses: int,
                   lig_props: Dict, prot_props: Dict) -> Dict:
    row = {
        "timestamp": datetime.now().isoformat(),
        "combo_name": combo_name,
        "protein_name": protein_name,
        "ligand_name": ligand_name,
        "status": status,
        "error_reason": error_reason,
        "elapsed_time_s": round(elapsed_time, 2),
        "total_poses": total_poses,
        "fpocket_poses": fpocket_poses,
        "p2rank_poses": p2rank_poses,
        "unguided_poses": unguided_poses,
        "failed_poses": failed_poses,
        "n_top_pockets": CFG.n_top_pockets,
        "poses_per_pocket": CFG.poses_per_pocket,
        "n_unguided_poses": CFG.n_unguided_poses,
        "force_pocket": CFG.force_pocket,
        "uff_minimize": CFG.uff_minimize,
        "device": CFG.equibind_device,
        "gpu_model": _GPU_MODEL,
    }
    row.update(lig_props)
    row.update(prot_props)
    return row


def init_docking_log(output_dir: Path) -> Tuple[Path, set]:
    """Create / load the docking log; return (path, set of existing combos)."""
    log_path = output_dir / CFG.docking_log_filename
    existing: set = set()
    if log_path.exists():
        try:
            df = pd.read_csv(log_path)
            existing = set(df["combo_name"])
        except Exception:
            pass
    else:
        pd.DataFrame(columns=DOCKING_LOG_COLUMNS).to_csv(log_path, index=False)
    return log_path, existing


def append_log_row(log_path: Path, existing_combos: set,
                   combo_name: str, protein_name: str, ligand_name: str,
                   status: str, error_reason: str, elapsed_time: float,
                   combo_results: Iterable, lig_props: Dict, prot_props: Dict) -> None:
    """Append a single combination's row to the docking log."""
    if combo_name in existing_combos:
        return

    fpocket_poses = sum(1 for r in combo_results if r.mode == "fpocket" and r.success)
    p2rank_poses = sum(1 for r in combo_results if r.mode == "p2rank" and r.success)
    unguided_poses = sum(1 for r in combo_results if r.mode == "unguided" and r.success)
    failed_poses = sum(1 for r in combo_results if not r.success)
    total_poses = fpocket_poses + p2rank_poses + unguided_poses

    row = _build_log_row(
        combo_name, protein_name, ligand_name, status, error_reason, elapsed_time,
        total_poses, fpocket_poses, p2rank_poses, unguided_poses, failed_poses,
        lig_props, prot_props,
    )
    pd.DataFrame([row], columns=DOCKING_LOG_COLUMNS).to_csv(
        log_path, mode="a", header=False, index=False)
    existing_combos.add(combo_name)


# ── Error log (JSON) ────────────────────────────────────────────────────

def load_error_log(output_dir: Path) -> Dict[str, dict]:
    p = output_dir / CFG.error_log_filename
    if p.exists():
        try:
            with open(p, "r") as f:
                return json.load(f)
        except Exception as e:
            print(f"  \u26a0 Could not read error log {p}: {e}")
    return {}


def save_error_log(output_dir: Path, error_log: Dict[str, dict]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / CFG.error_log_filename, "w") as f:
        json.dump(error_log, f, indent=2)


def record_failure(output_dir: Path, combo_name: str, error_message: str,
                   protein_name: str, ligand_name: str, elapsed_time: float) -> None:
    log = load_error_log(output_dir)
    log[combo_name] = {
        "error": error_message,
        "timestamp": datetime.now().isoformat(),
        "protein": protein_name,
        "ligand": ligand_name,
        "elapsed_time": round(elapsed_time, 2),
    }
    save_error_log(output_dir, log)
