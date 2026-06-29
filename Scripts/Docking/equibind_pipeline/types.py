"""Dataclasses for jobs flowing through Phase 1 -> 2 -> 3 of the pipeline."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Tuple

import numpy as np

from .pockets import PocketInfo


@dataclass
class TimingRecord:
    protein: str = ""
    ligand: str = ""
    phase1_prep_s: float = 0.0
    phase2_dock_s: float = 0.0
    phase3_post_s: float = 0.0
    total_s: float = 0.0
    n_poses_attempted: int = 0
    n_poses_success: int = 0


@dataclass
class GuidedPoseResult:
    protein_name: str
    ligand_name: str
    mode: str
    pocket_id: Optional[int]
    pocket_unique_id: Optional[str]
    pose_num: int
    pocket_center: Optional[Tuple[float, float, float]]
    pose_centroid: Optional[Tuple[float, float, float]]
    sdf_path: Optional[Path]
    success: bool
    error: str = ""
    uff_energy_before: Optional[float] = None
    uff_energy_after: Optional[float] = None
    uff_minimized: bool = False
    # Variant provenance (None == single-variant / legacy run):
    clamp_variant: Optional[str] = None     # "clampON" | "clampOFF"
    refine_variant: Optional[str] = None     # "raw" | "smina" | "gnina"
    refine_affinity: Optional[float] = None  # smina/gnina minimized affinity (kcal/mol)
    cnn_score: Optional[float] = None        # gnina CNN pose score (0..1); None for smina/raw
    cnn_affinity: Optional[float] = None     # gnina CNN predicted affinity (pK); None otherwise
    prep_time_s: float = 0.0        # Phase 1 input prep (shared across a pose's variants)
    dock_time_s: float = 0.0        # Phase 2 GPU inference (shared across a pose's variants)
    post_time_s: float = 0.0        # Phase 3 corrections + clamp + UFF (per clamp variant)
    refine_time_s: float = 0.0      # smina/gnina re-search cost (per refine variant; ~0 for raw)

    def to_dict(self) -> dict:
        return {
            "protein_name": self.protein_name, "ligand_name": self.ligand_name,
            "mode": self.mode, "pocket_id": self.pocket_id,
            "pocket_unique_id": self.pocket_unique_id, "pose_num": self.pose_num,
            "pocket_center": list(self.pocket_center) if self.pocket_center else None,
            "pose_centroid": list(self.pose_centroid) if self.pose_centroid else None,
            "sdf_path": str(self.sdf_path) if self.sdf_path else None,
            "success": self.success, "error": self.error,
            "uff_minimized": self.uff_minimized,
            "uff_energy_before": self.uff_energy_before,
            "uff_energy_after": self.uff_energy_after,
            "clamp_variant": self.clamp_variant,
            "refine_variant": self.refine_variant,
            "refine_affinity": self.refine_affinity,
            "cnn_score": self.cnn_score,
            "cnn_affinity": self.cnn_affinity,
            "prep_time_s": round(self.prep_time_s, 4),
            "dock_time_s": round(self.dock_time_s, 4),
            "post_time_s": round(self.post_time_s, 4),
            "refine_time_s": round(self.refine_time_s, 4),
        }


@dataclass
class DockJob:
    """A single docking job that flows through all 3 phases."""
    job_id: str
    combo_name: str
    protein_pdb: Path
    protein_key: str
    prepared_protein: Path
    ligand_sdf: Path
    ligand_mol: Optional[Any]
    lig_graph: Optional[Any]
    lig_coord: Optional[Any]
    geometry_graph: Optional[Any]
    seed: int
    pocket: Optional[PocketInfo]
    pocket_center: Optional[Tuple[float, float, float]]
    crop_offset: Optional[np.ndarray]
    combo_dir: Path
    final_sdf: Path
    pose_num: int
    protein_name: str
    ligand_name: str
    mode: str
    prep_time: float = 0.0
    # filled after Phase 2
    predicted_coords: Optional[Any] = None
    dock_time: float = 0.0
    gpu_success: bool = False
    error: str = ""


@dataclass
class P3Intermediate:
    """State passed from Stage 1 -> Stage 3 (kept in main process; not pickled)."""
    out_sdf: Path
    pose_dir: Path
    minimized_sdf: Path
    final_sdf: Path
    protein_pdb: Path
    protein_name: str
    ligand_name: str
    mode: str
    pocket_id: Optional[int]
    pocket_unique_id: Optional[str]
    pose_num: int
    pocket_center: Optional[Tuple[float, float, float]]
    prep_time: float
    dock_time: float
    pre_uff_time: float
    # "" for a single-variant run; "__clampON"/"__clampOFF" when clamp_mode=both.
    clamp_suffix: str = ""
    clamp_variant: Optional[str] = None
