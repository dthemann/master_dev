"""Centralised configuration for the EquiBind docking pipeline.

All tunables live here. Override via environment variables (``EQ_*``) or by
mutating ``CFG`` before calling :func:`equibind_pipeline.pipeline.run_pipeline`.
"""
from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import List


def _env_int(name: str, default: int) -> int:
    v = os.environ.get(name)
    return int(v) if v else default


def _env_float(name: str, default: float) -> float:
    v = os.environ.get(name)
    return float(v) if v else default


def _env_bool(name: str, default: bool) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v.strip().lower() in ("1", "true", "yes", "on")


def _env_path(name: str, default: Path) -> Path:
    v = os.environ.get(name)
    return Path(v) if v else default


@dataclass
class Config:
    # ── Workspace / inputs ────────────────────────────────────────────────
    workspace_root: Path = field(default_factory=Path.cwd)
    drugs_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_DRUGS_DIR", Path("storage/ligands_sdf_small_symmetric_approved")))
    receptors_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_RECEPTORS_DIR", Path.cwd() / "Orai"))
    receptor_name_filter: str = os.environ.get("EQ_RECEPTOR_FILTER", "_cleaned")

    # ── Pocket result dirs ────────────────────────────────────────────────
    fpocket_results_folder: Path = field(default_factory=lambda: _env_path(
        "EQ_FPOCKET_DIR", Path.cwd() / "pocket_results" / "fpocket_results"))
    p2rank_folder: Path = field(default_factory=lambda: _env_path(
        "EQ_P2RANK_DIR", Path.cwd() / "pocket_results" / "p2rank_results"))

    # ── EquiBind ──────────────────────────────────────────────────────────
    equibind_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_EQUIBIND_DIR", Path.home() / "tools" / "EquiBind"))
    equibind_device: str = os.environ.get("EQ_DEVICE", "cuda")
    gpu_batch_size: int = _env_int("EQ_GPU_BATCH_SIZE", 8)

    # ── Output ────────────────────────────────────────────────────────────
    pocket_guided_output_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_OUTPUT_DIR", Path.cwd() / "small_equibind_pocket_guided"))

    # ── Pocket-guided settings ────────────────────────────────────────────
    n_top_pockets: int = _env_int("EQ_N_TOP_POCKETS", 3)
    poses_per_pocket: int = _env_int("EQ_POSES_PER_POCKET", 5)
    n_unguided_poses: int = _env_int("EQ_N_UNGUIDED_POSES", 5)
    pocket_match_threshold: float = _env_float("EQ_POCKET_MATCH_THRESHOLD", 8.0)

    # ── Pocket enforcement ────────────────────────────────────────────────
    force_pocket: bool = _env_bool("EQ_FORCE_POCKET", True)
    clamp_pose_to_pocket: bool = _env_bool("EQ_CLAMP_POSE", True)

    # ── Protein cropping ──────────────────────────────────────────────────
    use_protein_cropping: bool = _env_bool("EQ_USE_CROPPING", True)
    pocket_crop_radius: float = _env_float("EQ_CROP_RADIUS", 8.0)
    pocket_crop_buffer: float = _env_float("EQ_CROP_BUFFER", 3.0)
    include_full_residues: bool = True

    # ── RDKit conformer generation ────────────────────────────────────────
    rdkit_seeds: List[int] = field(default_factory=lambda: [
        42, 123, 456, 789, 1001, 2022, 3141, 5926, 8675, 9999])
    conformers_per_seed: int = _env_int("EQ_CONFORMERS_PER_SEED", 1)

    # ── UFF minimization ──────────────────────────────────────────────────
    uff_minimize: bool = _env_bool("EQ_UFF_MINIMIZE", True)
    uff_max_iters: int = _env_int("EQ_UFF_MAX_ITERS", 200)
    uff_energy_tol: float = _env_float("EQ_UFF_ENERGY_TOL", 1e-4)
    uff_force_tol: float = _env_float("EQ_UFF_FORCE_TOL", 1e-3)
    uff_proximity_radius: float = _env_float("EQ_UFF_PROXIMITY_RADIUS", 6.0)
    uff_add_hydrogens: bool = _env_bool("EQ_UFF_ADD_H", True)
    uff_vdw_thresh: float = _env_float("EQ_UFF_VDW_THRESH", 0.1)
    uff_strip_nonstandard: bool = _env_bool("EQ_UFF_STRIP_NONSTANDARD", True)

    # ── Parallelism ───────────────────────────────────────────────────────
    n_parallel_workers: int = _env_int(
        "EQ_WORKERS", min(os.cpu_count() or 22, 26))
    parallel_conformer_gen: bool = _env_bool("EQ_PARALLEL_CONFORMER_GEN", True)

    # ── Skip / overwrite ──────────────────────────────────────────────────
    skip_existing: bool = _env_bool("EQ_SKIP_EXISTING", True)

    # ── Filenames ─────────────────────────────────────────────────────────
    docking_log_filename: str = "docking_log.csv"
    error_log_filename: str = "failed_docking.json"

    # ── Derived / discovered ──────────────────────────────────────────────
    reduce_executable: str = field(default_factory=lambda: shutil.which("reduce") or "")

    def validate(self) -> None:
        if not self.drugs_dir.exists():
            raise FileNotFoundError(f"Ligand directory missing: {self.drugs_dir}")
        if not self.receptors_dir.exists():
            raise FileNotFoundError(f"Receptor directory missing: {self.receptors_dir}")

    def ensure_output_dirs(self) -> None:
        self.pocket_guided_output_dir.mkdir(parents=True, exist_ok=True)

    def banner(self) -> None:
        print("=" * 80)
        print("EquiBind Batch Docking Configuration")
        print("=" * 80)
        print(f"  Workspace:        {self.workspace_root}")
        print(f"  Ligand dir:       {self.drugs_dir}")
        print(f"  Receptor dir:     {self.receptors_dir}")
        print(f"  EquiBind dir:     {self.equibind_dir}")
        print(f"  EquiBind device:  {self.equibind_device}")
        print(f"  Output:           {self.pocket_guided_output_dir}")
        print(f"  Pocket results:   fpocket={self.fpocket_results_folder}")
        print(f"                    p2rank={self.p2rank_folder}")
        print(f"  Top pockets:      {self.n_top_pockets} per method")
        print(f"  Poses/pocket:     {self.poses_per_pocket}")
        print(f"  Unguided poses:   {self.n_unguided_poses}")
        print(f"  Workers:          {self.n_parallel_workers}")
        print(f"  GPU batch size:   {self.gpu_batch_size}")
        print(f"  UFF minimize:     {self.uff_minimize}")
        print(f"  Skip existing:    {self.skip_existing}")
        print(f"  Force pocket:     {self.force_pocket}")
        print(f"  Protein cropping: {self.use_protein_cropping} "
              f"(radius={self.pocket_crop_radius}\u00c5 + buffer={self.pocket_crop_buffer}\u00c5)")
        print("=" * 80)


# Singleton — import this everywhere
CFG = Config()
