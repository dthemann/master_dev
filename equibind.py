# %%
# EquiBind Batch Docking Pipeline with Multiple Pose Generation
#
# This notebook docks all protein-ligand combinations using EquiBind,
# generating multiple diverse poses per combination using RDKit conformer generation.
#
# Strategy: Generate multiple 3D conformers with RDKit's EmbedMultipleConfs(),
# then run EquiBind on each conformer to produce diverse docked poses.
#
# Uses conda environment 'equibind' to run multiligand_inference.py:
#   conda run -n equibind python ~/docking_tools/EquiBind/multiligand_inference.py \
#     -o ./equibind_out -r protein.pdb -l ligand.sdf --device cpu

# %% [markdown]
# # Imports

# %%
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import os
import shutil
import subprocess
import time
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from enum import IntEnum
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdForceFieldHelpers
from rdkit.Chem import rdMolTransforms
from rdkit.ForceField import rdForceField

# Suppress RDKit deprecation warnings (e.g. GetValence)
RDLogger.DisableLog('rdApp.*')

# %% [markdown]
# # Config

# %%

# ── Configuration ──────────────────────────────────────────────────────────

# Paths
workspace_root = Path.cwd()
drugs_dir = workspace_root / "Drugs"
receptors_dir = workspace_root / "Orai"

# Pocket result directories (from previous runs)
fpocket_results_folder = workspace_root / "fpocket_results"
p2rank_folder = workspace_root / "p2rank_results"

# EquiBind paths
EQUIBIND_DIR = Path("/home/manndo/docking_tools/EquiBind")
EQUIBIND_MULTILIGAND_SCRIPT = EQUIBIND_DIR / "multiligand_inference.py"
EQUIBIND_DEVICE = "cuda"  # "cpu" or "cuda"

# Output directories
POCKET_GUIDED_OUTPUT_DIR = workspace_root / "equibind_pocket_guided"
POCKET_GUIDED_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EQUIBIND_BATCH_DIR = workspace_root / "equibind_batches"
EQUIBIND_OUTPUT_DIR = workspace_root / "equibind_docked_poses"
EQUIBIND_BATCH_DIR.mkdir(parents=True, exist_ok=True)
EQUIBIND_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

LOG_DIR = POCKET_GUIDED_OUTPUT_DIR / "wdlogs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

# Pocket-guided docking settings
N_TOP_POCKETS: int = 5          # Top-ranked pockets to use from each method
POSES_PER_POCKET: int = 3       # Poses to generate per pocket (using different conformers)
N_UNGUIDED_POSES: int = 10      # Natural EquiBind poses per protein-ligand pair
POCKET_MATCH_THRESHOLD: float = 8.0  # Distance (Å) to consider a pose "inside" a pocket

# Pocket enforcement settings
FORCE_POCKET: bool = True        # If True, reject poses whose centroid drifts outside the pocket
MAX_POCKET_TRIALS: int = 200      # Max conformer trials per desired pose when FORCE_POCKET is True
CLAMP_POSE_TO_POCKET: bool = True  # If True, clamp drifted poses back to pocket (instead of rejecting)

# Protein cropping settings (constrains EquiBind's search space to the pocket region)
USE_PROTEIN_CROPPING: bool = True   # If True, dock against cropped protein (pocket only)
POCKET_CROP_RADIUS: float = 6.0    # Residues within this radius of pocket center are kept
POCKET_CROP_BUFFER: float = 3.0     # Additional buffer for context

# Batch docking settings
NUM_POSES: int = 30             # Unique poses per protein-ligand combination
NUM_CONFORMERS: int = 10        # RDKit conformers to generate (>= NUM_POSES for diversity)
MAX_ATTEMPTS: int = NUM_POSES * 30  # Fallback if conformers fail
POSE_RMSD_THRESHOLD: float = 1.0    # Min RMSD (Å) between poses to consider distinct

# RDKit conformer generation parameters
RDKIT_SEEDS: list = [42, 123, 456, 789, 1001, 2022, 3141, 5926, 8675, 9999]
CONFORMERS_PER_SEED: int = 5
RDKIT_RANDOM_SEED: int = 42
RDKIT_NUM_THREADS: int = 0              # 0 = use all available threads
RDKIT_PRUNE_RMS_THRESH: float = 0.2     # Prune similar conformers during generation

# UFF post-docking minimization settings
UFF_MINIMIZE: bool = True            # Enable UFF minimization of docked poses
UFF_MAX_ITERS: int = 500             # Maximum minimization iterations
UFF_ENERGY_TOL: float = 1e-4         # Energy convergence tolerance (kcal/mol)
UFF_FORCE_TOL: float = 1e-3          # Force convergence tolerance
UFF_PROTEIN_CONSTRAINT_WT: float = 100.0  # Constraint weight for protein atoms (kcal/mol/Å²)
UFF_PROXIMITY_RADIUS: float = 8.0    # Only include protein residues within this Å of ligand
UFF_ADD_HYDROGENS: bool = True       # Add hydrogens before minimization for better UFF terms
UFF_VDW_THRESH: float = 0.1          # Non-bonded interaction cutoff threshold

# Overwrite / skip settings
SKIP_EXISTING: bool = False
OVERWRITE_EXISTING: bool = True  # Set to True to re-dock existing combinations

# Check for reduce executable (for adding hydrogens to proteins)
REDUCE_EXECUTABLE = shutil.which("reduce")

# ── Validate ───────────────────────────────────────────────────────────────
if not drugs_dir.exists():
    raise FileNotFoundError(f"Ligand directory missing: {drugs_dir}")
if not receptors_dir.exists():
    raise FileNotFoundError(f"Receptor directory missing: {receptors_dir}")

print("=" * 80)
print("EquiBind Batch Docking Configuration")
print("=" * 80)
print(f"Number of poses per combination: {NUM_POSES}")
print(f"RDKit conformers to generate:    {NUM_CONFORMERS}")
print(f"Maximum attempts per combination:{MAX_ATTEMPTS}")
print(f"Pose RMSD threshold:             {POSE_RMSD_THRESHOLD} Å")
print(f"RDKit prune RMS threshold:       {RDKIT_PRUNE_RMS_THRESH} Å")
print(f"EquiBind script:  {EQUIBIND_MULTILIGAND_SCRIPT}")
print(f"EquiBind device:  {EQUIBIND_DEVICE}")
print(f"Batch directory:  {EQUIBIND_BATCH_DIR}")
print(f"Output directory: {EQUIBIND_OUTPUT_DIR}")
print()
print(f"UFF post-docking minimization:   {'ON' if UFF_MINIMIZE else 'OFF'}")
if UFF_MINIMIZE:
    print(f"  Max iterations:                {UFF_MAX_ITERS}")
    print(f"  Protein constraint weight:     {UFF_PROTEIN_CONSTRAINT_WT} kcal/mol/Å²")
    print(f"  Proximity radius:              {UFF_PROXIMITY_RADIUS} Å")
print()

if not EQUIBIND_DIR.exists():
    print(f"⚠️  WARNING: EquiBind directory not found at {EQUIBIND_DIR}")
    print("   Please update EQUIBIND_DIR to point to your EquiBind installation")
elif not EQUIBIND_MULTILIGAND_SCRIPT.exists():
    print(f"⚠️  WARNING: multiligand_inference.py not found at {EQUIBIND_MULTILIGAND_SCRIPT}")
else:
    print(f"✓ EquiBind script found")

# %% [markdown]
# # Signal Monitor

# %%
# ── Signal Monitor ──────────────────────────────────────────────────────────
# Three verbosity levels control what gets printed during docking.
#
#   MonitorLevel.ALL      → Everything: EquiBind invocations, returned poses,
#                           pocket distance checks, conformer generation, etc.
#   MonitorLevel.WARNING  → Rejected poses, dock failures, pocket drift,
#                           plus everything from CRITICAL.
#   MonitorLevel.CRITICAL → Only fatal errors (timeouts, missing output, total
#                           combo failures).
# ────────────────────────────────────────────────────────────────────────────


class MonitorLevel(IntEnum):
    """Signal monitoring verbosity levels."""
    ALL = 1          # Full trace
    WARNING = 2      # Warnings + critical only
    CRITICAL = 3     # Fatal errors only

# ── Set the active monitoring level here ──────────────────────────────────
# Change this to MonitorLevel.WARNING or MonitorLevel.CRITICAL to reduce output.
MONITOR_LEVEL = MonitorLevel.ALL


class DockingMonitor:
    """Structured signal monitor for the EquiBind docking pipeline.

    Every message is tagged with a level; only messages at or above the
    configured threshold are printed.  Colour prefixes make it easy to
    scan terminal output.
    """

    _PREFIXES = {
        MonitorLevel.ALL:      "    ·",
        MonitorLevel.WARNING:  " ⚠️ ",
        MonitorLevel.CRITICAL: " 🔴",
    }

    def __init__(self, level: MonitorLevel = MonitorLevel.ALL):
        self.level = level
        self.call_count = 0
        self.success_count = 0
        self.fail_count = 0
        self.rejected_count = 0
        self.in_pocket_count = 0
        self.outside_pocket_count = 0

    # ── core dispatcher ──────────────────────────────────────────────────
    def _emit(self, lvl: MonitorLevel, msg: str) -> None:
        if lvl >= self.level:
            prefix = self._PREFIXES.get(lvl, "")
            print(f"{prefix} {msg}")

    # ── convenience loggers ──────────────────────────────────────────────
    def info(self, msg: str) -> None:
        """Verbose trace (ALL level)."""
        self._emit(MonitorLevel.ALL, msg)

    def warning(self, msg: str) -> None:
        """Non-fatal issue (WARNING level)."""
        self._emit(MonitorLevel.WARNING, msg)

    def critical(self, msg: str) -> None:
        """Fatal / blocking error (CRITICAL level)."""
        self._emit(MonitorLevel.CRITICAL, msg)

    # ── EquiBind call / return ───────────────────────────────────────────
    def equibind_call(self, protein: str, ligand: str, output_dir: str,
                      seed: int, device: str) -> None:
        """Log an outgoing EquiBind invocation."""
        self.call_count += 1
        self.info(
            f"[CALL #{self.call_count}] EquiBind ← "
            f"protein={protein}, ligand={ligand}, "
            f"seed={seed}, device={device}, out={output_dir}"
        )

    def equibind_return(self, success: bool, output_sdf: Optional[str],
                        error: str = "") -> None:
        """Log what EquiBind returned."""
        if success:
            self.success_count += 1
            self.info(f"[RECV] EquiBind → OK  output={output_sdf}")
        else:
            self.fail_count += 1
            self.warning(f"[RECV] EquiBind → FAIL  error={error[:200]}")

    # ── Pocket proximity ─────────────────────────────────────────────────
    def pose_accepted_in_pocket(self, pose_id: str, pocket_id: str,
                                distance: float, threshold: float) -> None:
        """Log a pose that landed inside the target pocket."""
        self.in_pocket_count += 1
        self.info(
            f"[POCKET ✓] {pose_id} INSIDE {pocket_id}  "
            f"(dist={distance:.1f}Å ≤ {threshold:.1f}Å)"
        )

    def pose_rejected_from_pocket(self, pose_id: str, pocket_id: str,
                                  distance: float, threshold: float) -> None:
        """Log a pose that drifted outside the target pocket."""
        self.rejected_count += 1
        self.warning(
            f"[POCKET ✗] {pose_id} OUTSIDE {pocket_id}  "
            f"(dist={distance:.1f}Å > {threshold:.1f}Å) → REJECTED"
        )

    def pose_outside_all_pockets(self, pose_id: str,
                                 nearest_pocket: str,
                                 nearest_dist: float,
                                 threshold: float) -> None:
        """Log an unguided pose that is not near any known pocket."""
        self.outside_pocket_count += 1
        self.warning(
            f"[POCKET ✗] {pose_id} not in any pocket  "
            f"(nearest={nearest_pocket}, dist={nearest_dist:.1f}Å > {threshold:.1f}Å)"
        )

    def pose_inside_pocket(self, pose_id: str, pocket_id: str,
                           distance: float, threshold: float,
                           pocket_source: str) -> None:
        """Log an unguided pose that falls inside a known pocket."""
        self.in_pocket_count += 1
        self.info(
            f"[POCKET ✓] {pose_id} in {pocket_source} pocket {pocket_id}  "
            f"(dist={distance:.1f}Å ≤ {threshold:.1f}Å)"
        )

    # ── Section headers ──────────────────────────────────────────────────
    def header(self, msg: str) -> None:
        """Always printed regardless of level."""
        print(f"\n{'─' * 70}")
        print(f"  {msg}")
        print(f"{'─' * 70}")

    def section(self, msg: str) -> None:
        """Printed at ALL level."""
        self._emit(MonitorLevel.ALL, f"── {msg} ──")

    # ── Summary ──────────────────────────────────────────────────────────
    def print_summary(self) -> None:
        """Print accumulated counters (always shown)."""
        print(f"\n{'═' * 70}")
        print("  SIGNAL MONITOR SUMMARY")
        print(f"{'═' * 70}")
        print(f"  EquiBind calls:        {self.call_count}")
        print(f"  Successful returns:    {self.success_count}")
        print(f"  Failed returns:        {self.fail_count}")
        print(f"  Poses inside pocket:   {self.in_pocket_count}")
        print(f"  Poses rejected/outside:{self.rejected_count + self.outside_pocket_count}")
        print(f"    ├ guided rejected:   {self.rejected_count}")
        print(f"    └ unguided outside:  {self.outside_pocket_count}")
        print(f"{'═' * 70}")


# Instantiate the global monitor using the configured level
monitor = DockingMonitor(level=MONITOR_LEVEL)
print(f"DockingMonitor initialised  —  level = {MONITOR_LEVEL.name}")
print(f"  ALL      = full trace of every call, return, and pocket check")
print(f"  WARNING  = only rejected poses, failures, and pocket drift")
print(f"  CRITICAL = only fatal errors")

# %% [markdown]
# # Validate GPU Support

# %%
import torch
import subprocess

# ============================================================================
# VALIDATE GPU SUPPORT
# ============================================================================


print("=" * 80)
print("GPU SUPPORT VALIDATION")
print("=" * 80)

# Check PyTorch CUDA availability
print("\nPyTorch CUDA Status:")
print(f"  PyTorch version: {torch.__version__}")
print(f"  CUDA available: {torch.cuda.is_available()}")

if torch.cuda.is_available():
    print(f"  CUDA version: {torch.version.cuda}")
    print(f"  Number of GPUs: {torch.cuda.device_count()}")
    print(f"  Current GPU: {torch.cuda.current_device()}")
    print(f"  GPU Name: {torch.cuda.get_device_name(0)}")
    print(f"  GPU Memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.2f} GB")
    
    # Test GPU with a simple tensor operation
    try:
        x = torch.randn(100, 100).cuda()
        y = torch.randn(100, 100).cuda()
        z = torch.matmul(x, y)
        print(f"  ✓ GPU tensor operations working")
    except Exception as e:
        print(f"  ✗ GPU tensor test failed: {e}")
else:
    print("  ⚠️  No CUDA-capable GPU detected")
    print("  EquiBind will use CPU (slower)")

# Check nvidia-smi
print("\nNVIDIA GPU Information (nvidia-smi):")
try:
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        print(result.stdout)
    else:
        print("  ⚠️  nvidia-smi not available or failed")
except FileNotFoundError:
    print("  ⚠️  nvidia-smi command not found")
except Exception as e:
    print(f"  ⚠️  Error running nvidia-smi: {e}")

# Recommendation
print("\n" + "=" * 80)
if torch.cuda.is_available():
    print("✓ GPU ACCELERATION AVAILABLE")
    print(f"  Recommended setting: EQUIBIND_DEVICE = 'cuda'")
    print(f"  Current setting: EQUIBIND_DEVICE = '{EQUIBIND_DEVICE}'")
else:
    print("⚠️  GPU NOT AVAILABLE - Using CPU")
    print(f"  Recommended setting: EQUIBIND_DEVICE = 'cpu'")
    print(f"  Current setting: EQUIBIND_DEVICE = '{EQUIBIND_DEVICE}'")
print("=" * 80)

# %% [markdown]
# # Pocket Cropping

# %%
# ============================================================================
# PROTEIN CROPPING FOR POCKET-CONSTRAINED DOCKING
# ============================================================================
# This cell implements protein cropping to constrain EquiBind's search space
# to a specific binding pocket. By extracting only residues within a defined
# radius of the pocket center, EquiBind physically cannot predict binding
# poses outside the pocket region.
#
# Key functions:
#   - crop_protein_to_pocket(): Extract residues near pocket center
#   - get_cropped_protein_for_pocket(): Cache management for cropped proteins
#   - Offset tracking to transform coordinates back to original frame
#
# Configuration is in cell 5:
#   USE_PROTEIN_CROPPING, POCKET_CROP_RADIUS, POCKET_CROP_BUFFER
# ============================================================================

from typing import Set, Dict, Tuple, Optional, List
from pathlib import Path
import numpy as np

# Include all atoms of a residue if any atom is in range
INCLUDE_FULL_RESIDUES: bool = True

# Cache: (protein_path, pocket_unique_id) → (cropped_pdb_path, offset_vector)
_cropped_protein_cache: Dict[Tuple[str, str], Tuple[Path, np.ndarray]] = {}


def _parse_pdb_atoms(pdb_path: Path) -> List[dict]:
    """Parse ATOM/HETATM records from a PDB file.
    
    Returns list of dicts with keys:
        line, record_type, atom_num, atom_name, res_name, chain, res_num, x, y, z
    """
    atoms = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    atoms.append({
                        'line': line,
                        'record_type': line[0:6].strip(),
                        'atom_num': int(line[6:11]),
                        'atom_name': line[12:16].strip(),
                        'res_name': line[17:20].strip(),
                        'chain': line[21:22].strip() or 'A',
                        'res_num': int(line[22:26]),
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                    })
                except (ValueError, IndexError):
                    continue
    return atoms


def _get_residues_near_center(
    atoms: List[dict],
    center: Tuple[float, float, float],
    radius: float,
    include_full_residues: bool = True,
) -> Set[Tuple[str, int]]:
    """Find residues with at least one atom within radius of center.
    
    Returns set of (chain, res_num) tuples.
    """
    center_np = np.array(center)
    near_residues: Set[Tuple[str, int]] = set()
    
    for atom in atoms:
        coord = np.array([atom['x'], atom['y'], atom['z']])
        dist = np.linalg.norm(coord - center_np)
        if dist <= radius:
            near_residues.add((atom['chain'], atom['res_num']))
    
    return near_residues


def crop_protein_to_pocket(
    protein_pdb: Path,
    pocket_center: Tuple[float, float, float],
    output_pdb: Path,
    crop_radius: float = POCKET_CROP_RADIUS,
    buffer: float = POCKET_CROP_BUFFER,
    include_full_residues: bool = INCLUDE_FULL_RESIDUES,
) -> Tuple[bool, np.ndarray, int, int]:
    """Crop a protein PDB to residues near a pocket center.
    
    Args:
        protein_pdb: Input full protein PDB
        pocket_center: (x, y, z) pocket center coordinates
        output_pdb: Output cropped PDB path
        crop_radius: Radius around pocket center to include
        buffer: Additional buffer for context residues
        include_full_residues: If True, include all atoms of selected residues
        
    Returns:
        (success, offset_vector, n_residues_kept, n_atoms_kept)
        
        offset_vector is the translation applied to center the pocket at origin.
        To convert cropped coordinates back to original frame: coord + offset_vector
    """
    effective_radius = crop_radius + buffer
    
    # Parse full protein
    atoms = _parse_pdb_atoms(protein_pdb)
    if not atoms:
        return False, np.zeros(3), 0, 0
    
    # Find residues within range
    near_residues = _get_residues_near_center(
        atoms, pocket_center, effective_radius, include_full_residues
    )
    
    if not near_residues:
        monitor.warning(f"No residues found within {effective_radius}Å of pocket center")
        return False, np.zeros(3), 0, 0
    
    # Filter atoms to keep
    if include_full_residues:
        kept_atoms = [a for a in atoms if (a['chain'], a['res_num']) in near_residues]
    else:
        # Only keep atoms actually within radius
        center_np = np.array(pocket_center)
        kept_atoms = []
        for a in atoms:
            coord = np.array([a['x'], a['y'], a['z']])
            if np.linalg.norm(coord - center_np) <= effective_radius:
                kept_atoms.append(a)
    
    if not kept_atoms:
        return False, np.zeros(3), 0, 0
    
    # Compute offset: we'll translate so pocket center is at origin
    # This helps EquiBind by centering the binding region
    offset = np.array(pocket_center)
    
    # Write cropped PDB with translated coordinates
    output_pdb.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_pdb, 'w', encoding='utf-8') as f:
        f.write(f"REMARK   CROPPED PROTEIN - pocket center at origin\n")
        f.write(f"REMARK   Original pocket center: {pocket_center[0]:.3f} {pocket_center[1]:.3f} {pocket_center[2]:.3f}\n")
        f.write(f"REMARK   Crop radius: {effective_radius:.1f} A\n")
        f.write(f"REMARK   Residues kept: {len(near_residues)}\n")
        f.write(f"REMARK   To restore original coords: add offset ({offset[0]:.3f}, {offset[1]:.3f}, {offset[2]:.3f})\n")
        
        atom_num = 1
        for a in kept_atoms:
            # Translate to center pocket at origin
            new_x = a['x'] - offset[0]
            new_y = a['y'] - offset[1]
            new_z = a['z'] - offset[2]
            
            # Reconstruct PDB line with new coordinates
            line = a['line']
            new_line = (
                f"{line[0:6]}"  # record type
                f"{atom_num:5d}"  # atom number (renumbered)
                f"{line[11:30]}"  # atom name, res name, chain, res num
                f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}"  # coordinates
                f"{line[54:]}"  # rest of line (occupancy, B-factor, element)
            )
            f.write(new_line)
            atom_num += 1
        
        f.write("END\n")
    
    n_residues = len(near_residues)
    n_atoms = len(kept_atoms)
    
    return True, offset, n_residues, n_atoms


def get_cropped_protein_for_pocket(
    protein_pdb: Path,
    pocket: 'PocketInfo',  # Forward reference - defined in cell 12
    prep_dir: Path,
    crop_radius: float = POCKET_CROP_RADIUS,
) -> Tuple[Optional[Path], np.ndarray]:
    """Get or create a cropped protein PDB for a specific pocket.
    
    Uses caching to avoid re-cropping for the same protein/pocket combination.
    
    Args:
        protein_pdb: Full protein PDB path
        pocket: PocketInfo object with center coordinates
        prep_dir: Directory to store cropped PDB files
        crop_radius: Radius around pocket center
        
    Returns:
        (cropped_pdb_path, offset_vector) or (None, zeros) on failure
        
        The offset_vector should be added to docked coordinates to get
        original-frame coordinates.
    """
    cache_key = (str(protein_pdb), pocket.unique_id)
    
    if cache_key in _cropped_protein_cache:
        return _cropped_protein_cache[cache_key]
    
    # Create cropped protein
    cropped_pdb = prep_dir / f"{protein_pdb.stem}_crop_{pocket.unique_id}.pdb"
    
    success, offset, n_res, n_atoms = crop_protein_to_pocket(
        protein_pdb=protein_pdb,
        pocket_center=pocket.center,
        output_pdb=cropped_pdb,
        crop_radius=crop_radius,
    )
    
    if success:
        monitor.info(f"Cropped protein for {pocket.unique_id}: {n_res} residues, {n_atoms} atoms "
                     f"(radius={crop_radius + POCKET_CROP_BUFFER:.1f}Å)")
        _cropped_protein_cache[cache_key] = (cropped_pdb, offset)
        return cropped_pdb, offset
    else:
        monitor.warning(f"Failed to crop protein for pocket {pocket.unique_id}")
        return None, np.zeros(3)


def translate_pose_back_to_original_frame(
    sdf_path: Path,
    offset: np.ndarray,
    output_path: Path,
) -> bool:
    """Translate a docked pose back to the original protein coordinate frame.
    
    After docking against a cropped protein (centered at origin), this function
    shifts the ligand coordinates back to the original frame.
    
    Args:
        sdf_path: Input SDF from docking against cropped protein
        offset: The offset vector returned by crop_protein_to_pocket
        output_path: Output SDF in original coordinate frame
        
    Returns:
        True on success, False on failure
    """
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            return False
        
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            # Add offset to restore original coordinates
            conf.SetAtomPosition(i, (
                pos.x + offset[0],
                pos.y + offset[1],
                pos.z + offset[2],
            ))
        
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        
        return output_path.exists() and output_path.stat().st_size > 0
    except Exception as e:
        monitor.warning(f"Failed to translate pose back to original frame: {e}")
        return False


def clamp_pose_to_pocket(
    sdf_path: Path,
    pocket_center: Tuple[float, float, float],
    max_distance: float,
    output_path: Path,
) -> Tuple[bool, float, float]:
    """Clamp a docked pose so its centroid is within max_distance of pocket center.
    
    If the pose centroid is farther than max_distance from the pocket center,
    the entire pose is translated to bring the centroid exactly to max_distance
    from the pocket center (along the line from pocket to current centroid).
    
    This ensures ALL poses stay within the defined pocket region.
    
    Args:
        sdf_path: Input SDF file
        pocket_center: Target pocket center (x, y, z)
        max_distance: Maximum allowed distance from pocket center
        output_path: Output SDF with clamped coordinates
        
    Returns:
        (success, original_distance, new_distance)
    """
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            return False, 0.0, 0.0
        
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        
        # Compute current centroid
        positions = conf.GetPositions()
        centroid = positions.mean(axis=0)
        pocket_np = np.array(pocket_center)
        
        # Compute distance from pocket center
        displacement = centroid - pocket_np
        original_distance = np.linalg.norm(displacement)
        
        if original_distance <= max_distance:
            # Already within bounds, just copy
            shutil.copy2(sdf_path, output_path)
            return True, original_distance, original_distance
        
        # Compute shift needed: move centroid to exactly max_distance from pocket
        # Unit vector from pocket to centroid
        direction = displacement / original_distance
        target_centroid = pocket_np + direction * max_distance
        shift = target_centroid - centroid
        
        # Apply shift to all atoms
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (
                pos.x + shift[0],
                pos.y + shift[1],
                pos.z + shift[2],
            ))
        
        # Write clamped SDF
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        
        new_distance = max_distance  # By construction
        
        return output_path.exists() and output_path.stat().st_size > 0, original_distance, new_distance
        
    except Exception as e:
        monitor.warning(f"Failed to clamp pose to pocket: {e}")
        return False, 0.0, 0.0


def clear_cropped_protein_cache():
    """Clear the cropped protein cache (useful for memory management)."""
    global _cropped_protein_cache
    _cropped_protein_cache.clear()
    monitor.info("Cropped protein cache cleared")


# ── Statistics helper ─────────────────────────────────────────────────────

def estimate_search_space_reduction(
    protein_pdb: Path,
    pocket_center: Tuple[float, float, float],
    crop_radius: float = POCKET_CROP_RADIUS,
) -> Tuple[int, int, float]:
    """Estimate the search space reduction from cropping.
    
    Returns:
        (original_atoms, cropped_atoms, reduction_percent)
    """
    atoms = _parse_pdb_atoms(protein_pdb)
    original_count = len(atoms)
    
    effective_radius = crop_radius + POCKET_CROP_BUFFER
    near_residues = _get_residues_near_center(atoms, pocket_center, effective_radius)
    cropped_atoms = [a for a in atoms if (a['chain'], a['res_num']) in near_residues]
    cropped_count = len(cropped_atoms)
    
    reduction = 100.0 * (1.0 - cropped_count / original_count) if original_count > 0 else 0.0
    
    return original_count, cropped_count, reduction


print("✓ Protein cropping functions loaded")
print(f"  Default crop radius: {POCKET_CROP_RADIUS}Å (+ {POCKET_CROP_BUFFER}Å buffer)")
print(f"  Include full residues: {INCLUDE_FULL_RESIDUES}")

# %%
# ============================================================================
# UFF POST-DOCKING MINIMIZATION
# ============================================================================
# After EquiBind produces a docked pose, the ligand geometry may contain
# steric clashes or non-ideal bond lengths/angles. This section applies
# RDKit's Universal Force Field (UFF) to minimize the ligand in the
# context of nearby protein atoms.
#
# Strategy:
#   1. Load the docked ligand and protein structures.
#   2. Extract protein residues within UFF_PROXIMITY_RADIUS of the ligand.
#   3. Combine ligand + nearby protein into a single RDKit molecule.
#   4. Build UFF force field; constrain protein atoms (high force constant).
#   5. Minimize — only ligand atoms are free to move.
#   6. Extract minimized ligand coordinates and write to output SDF.
#
# Configuration (in Config cell above):
#   UFF_MINIMIZE, UFF_MAX_ITERS, UFF_ENERGY_TOL, UFF_PROTEIN_CONSTRAINT_WT,
#   UFF_PROXIMITY_RADIUS, UFF_ADD_HYDROGENS
# ============================================================================

from copy import deepcopy


def _extract_nearby_protein_atoms(
    protein_pdb: Path,
    ligand_coords: np.ndarray,
    radius: float = UFF_PROXIMITY_RADIUS,
) -> Optional[Chem.Mol]:
    """Extract protein residues with any atom within `radius` of any ligand atom.

    Returns an RDKit Mol of the nearby protein subset, or None on failure.
    """
    try:
        protein_mol = Chem.MolFromPDBFile(str(protein_pdb), removeHs=False, sanitize=False)
        if protein_mol is None:
            return None

        # Try sanitization; fall back to partial sanitize on failure
        try:
            Chem.SanitizeMol(protein_mol)
        except Exception:
            try:
                Chem.SanitizeMol(
                    protein_mol,
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                )
            except Exception:
                pass

        prot_conf = protein_mol.GetConformer()
        prot_positions = prot_conf.GetPositions()

        # Find protein atoms near any ligand atom
        near_atom_indices = set()
        for lig_xyz in ligand_coords:
            dists = np.linalg.norm(prot_positions - lig_xyz, axis=1)
            near_atom_indices.update(np.where(dists <= radius)[0])

        if not near_atom_indices:
            return None

        # Expand to full residues so UFF has complete residue topology
        atom_info = protein_mol.GetAtomWithIdx(0).GetPDBResidueInfo()
        if atom_info is not None:
            # Collect residue ids for nearby atoms
            nearby_residues = set()
            for idx in near_atom_indices:
                info = protein_mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
                if info:
                    nearby_residues.add((info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode()))
            # Include all atoms from nearby residues
            full_atom_indices = set()
            for i in range(protein_mol.GetNumAtoms()):
                info = protein_mol.GetAtomWithIdx(i).GetPDBResidueInfo()
                if info and (info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode()) in nearby_residues:
                    full_atom_indices.add(i)
            near_atom_indices = full_atom_indices

        # Create editable molecule with only nearby atoms
        emol = Chem.RWMol(protein_mol)
        atoms_to_remove = sorted(
            set(range(protein_mol.GetNumAtoms())) - near_atom_indices, reverse=True
        )
        for idx in atoms_to_remove:
            emol.RemoveAtom(idx)

        nearby_mol = emol.GetMol()
        return nearby_mol

    except Exception as e:
        monitor.warning(f"Failed to extract nearby protein atoms: {e}")
        return None


def _uff_minimize_pose(
    docked_sdf: Path,
    protein_pdb: Path,
    output_sdf: Path,
    max_iters: int = UFF_MAX_ITERS,
    energy_tol: float = UFF_ENERGY_TOL,
    force_tol: float = UFF_FORCE_TOL,
    constraint_weight: float = UFF_PROTEIN_CONSTRAINT_WT,
    proximity_radius: float = UFF_PROXIMITY_RADIUS,
    add_hs: bool = UFF_ADD_HYDROGENS,
) -> Tuple[bool, float, float, str]:
    """Minimize a docked ligand using UFF with protein environment constraints.

    The protein atoms are held in place with harmonic position constraints
    while the ligand is free to relax.  This resolves steric clashes and
    improves bond geometries without letting the ligand drift away.

    Args:
        docked_sdf:       Path to the docked ligand SDF.
        protein_pdb:      Path to the protein PDB (full or cropped).
        output_sdf:       Where to write the minimized ligand SDF.
        max_iters:        Maximum optimization steps.
        energy_tol:       Energy convergence criterion (kcal/mol).
        force_tol:        Force convergence criterion.
        constraint_weight: Harmonic constraint weight for protein atoms.
        proximity_radius:  Include protein residues within this Å of ligand.
        add_hs:           Add hydrogens before FF setup for better UFF terms.

    Returns:
        (success, energy_before, energy_after, error_message)
    """
    try:
        # ── Load docked ligand ──
        suppl = Chem.SDMolSupplier(str(docked_sdf), removeHs=False, sanitize=False)
        lig_mol = next(iter(suppl), None)
        if lig_mol is None:
            return False, 0.0, 0.0, "Could not load docked ligand SDF"

        try:
            Chem.SanitizeMol(lig_mol)
        except Exception:
            # Attempt partial sanitization for difficult molecules
            try:
                Chem.SanitizeMol(
                    lig_mol,
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                )
            except Exception as e:
                return False, 0.0, 0.0, f"Ligand sanitization failed: {e}"

        lig_conf = lig_mol.GetConformer()
        lig_coords = lig_conf.GetPositions()
        n_lig_atoms = lig_mol.GetNumAtoms()

        # ── Extract nearby protein atoms ──
        prot_nearby = _extract_nearby_protein_atoms(protein_pdb, lig_coords, proximity_radius)
        if prot_nearby is None or prot_nearby.GetNumAtoms() == 0:
            # No protein context — minimize ligand alone
            monitor.info("UFF: No nearby protein atoms found, minimizing ligand in vacuum")
            if add_hs:
                lig_mol = Chem.AddHs(lig_mol, addCoords=True)
                n_lig_atoms_with_h = lig_mol.GetNumAtoms()
            else:
                n_lig_atoms_with_h = n_lig_atoms

            if not rdForceFieldHelpers.UFFHasAllMoleculeParams(lig_mol):
                return False, 0.0, 0.0, "UFF missing parameters for ligand atoms"

            ff = rdForceFieldHelpers.UFFGetMoleculeForceField(lig_mol)
            if ff is None:
                return False, 0.0, 0.0, "Could not create UFF force field for ligand"

            e_before = ff.CalcEnergy()
            converged = ff.Minimize(maxIts=max_iters, energyTol=energy_tol, forceTol=force_tol)
            e_after = ff.CalcEnergy()

            # Remove added Hs to match original atom count
            if add_hs:
                lig_mol = Chem.RemoveHs(lig_mol)

            writer = Chem.SDWriter(str(output_sdf))
            writer.write(lig_mol)
            writer.close()
            return True, e_before, e_after, ""

        n_prot_atoms = prot_nearby.GetNumAtoms()

        # ── Optionally add hydrogens ──
        if add_hs:
            try:
                lig_mol = Chem.AddHs(lig_mol, addCoords=True)
            except Exception:
                pass  # proceed without added Hs
            try:
                prot_nearby = Chem.AddHs(prot_nearby, addCoords=True)
            except Exception:
                pass

        n_lig_with_h = lig_mol.GetNumAtoms()
        n_prot_with_h = prot_nearby.GetNumAtoms()

        # ── Combine ligand + protein into one molecule ──
        combined = Chem.CombineMols(lig_mol, prot_nearby)

        # Verify the combined mol has a valid conformer
        if combined.GetNumConformers() == 0:
            return False, 0.0, 0.0, "Combined molecule has no conformer"

        # ── Check UFF parameter coverage ──
        if not rdForceFieldHelpers.UFFHasAllMoleculeParams(combined):
            # Fallback: try ligand-only minimization
            monitor.warning("UFF missing parameters for some atoms in combined system, "
                            "falling back to ligand-only minimization")
            if rdForceFieldHelpers.UFFHasAllMoleculeParams(lig_mol):
                ff = rdForceFieldHelpers.UFFGetMoleculeForceField(lig_mol)
                if ff is not None:
                    e_before = ff.CalcEnergy()
                    ff.Minimize(maxIts=max_iters, energyTol=energy_tol, forceTol=force_tol)
                    e_after = ff.CalcEnergy()
                    if add_hs:
                        lig_mol = Chem.RemoveHs(lig_mol)
                    writer = Chem.SDWriter(str(output_sdf))
                    writer.write(lig_mol)
                    writer.close()
                    return True, e_before, e_after, "ligand-only (protein params missing)"
            return False, 0.0, 0.0, "UFF parameters missing for combined system"

        # ── Build UFF force field ──
        ff = rdForceFieldHelpers.UFFGetMoleculeForceField(
            combined,
            vdwThresh=UFF_VDW_THRESH,
        )
        if ff is None:
            return False, 0.0, 0.0, "Could not create UFF force field for combined system"

        # ── Constrain protein atoms (indices n_lig_with_h .. end) ──
        for i in range(n_lig_with_h, n_lig_with_h + n_prot_with_h):
            ff.AddFixedPoint(i)

        # ── Minimize ──
        e_before = ff.CalcEnergy()
        converged = ff.Minimize(
            maxIts=max_iters,
            energyTol=energy_tol,
            forceTol=force_tol,
        )
        e_after = ff.CalcEnergy()

        convergence_msg = "converged" if converged == 0 else f"not converged (code {converged})"

        # ── Extract minimized ligand coordinates ──
        combined_conf = combined.GetConformer()
        out_lig = deepcopy(lig_mol)
        out_conf = out_lig.GetConformer()
        for i in range(n_lig_with_h):
            pos = combined_conf.GetAtomPosition(i)
            out_conf.SetAtomPosition(i, pos)

        # Remove added Hs to match original atom count
        if add_hs:
            out_lig = Chem.RemoveHs(out_lig)

        # ── Write minimized ligand ──
        output_sdf.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(out_lig)
        writer.close()

        monitor.info(f"UFF minimization: E {e_before:.1f} → {e_after:.1f} kcal/mol "
                     f"(ΔE={e_after - e_before:.1f}), {convergence_msg}, "
                     f"{n_lig_with_h} ligand + {n_prot_with_h} protein atoms")

        return True, e_before, e_after, convergence_msg

    except Exception as e:
        monitor.warning(f"UFF minimization failed: {e}")
        return False, 0.0, 0.0, str(e)


print("✓ UFF post-docking minimization functions loaded")
print(f"  UFF minimization enabled: {UFF_MINIMIZE}")
print(f"  Max iterations: {UFF_MAX_ITERS}")
print(f"  Protein constraint weight: {UFF_PROTEIN_CONSTRAINT_WT} kcal/mol/Å²")
print(f"  Proximity radius: {UFF_PROXIMITY_RADIUS} Å")
print(f"  Add hydrogens: {UFF_ADD_HYDROGENS}")

# %%
import sys
import dgl
import dgl.function as fn
from copy import deepcopy
from types import SimpleNamespace

import yaml
from rdkit.Geometry import Point3D

# ── DGL compatibility shim ────────────────────────────────────────────────
# Newer DGL versions removed dgl.function.copy_edge (renamed to copy_e).
# EquiBind's model code still uses copy_edge, so we patch it back in.
if not hasattr(fn, 'copy_edge'):
    fn.copy_edge = fn.copy_e
    print("  ✓ Patched dgl.function.copy_edge → copy_e (DGL compat shim)")

# ── 0. Import EquiBind in-process ─────────────────────────────────────────
#    Add EquiBind to sys.path so we can import its modules directly,
#    avoiding subprocess + conda activation overhead per call.

# %% [markdown]
# # Pocket-Guided Docking

# %%
# ============================================================================
# POCKET-GUIDED EQUIBIND DOCKING (Multi-Pose per Pocket)
# ============================================================================
# Runs EquiBind docking guided by fpocket and p2rank pocket predictions,
# plus unguided (natural) EquiBind poses. Then produces statistics on how
# many unguided poses overlap with fpocket / p2rank pockets.
#
# Generates N poses per pocket using different RDKit conformers, each
# translated to the pocket center as starting geometry.
#
# REFACTORED: EquiBind is now loaded IN-PROCESS (no subprocess/conda overhead).
#   - Model loaded once globally.
#   - Receptor graph loaded once per protein and cached.
#   - Each docking call is a direct Python function call (~100x faster).
#
# All runtime output is governed by the DockingMonitor (see Signal Monitor
# cell above). Set MONITOR_LEVEL in Config to control verbosity:
#   MonitorLevel.ALL      → full trace
#   MonitorLevel.WARNING  → rejected poses & failures only
#   MonitorLevel.CRITICAL → fatal errors only
# ============================================================================

import sys
import dgl
from copy import deepcopy
from types import SimpleNamespace

import yaml
from rdkit.Geometry import Point3D

# ── 0. Import EquiBind in-process ─────────────────────────────────────────
#    Add EquiBind to sys.path so we can import its modules directly,
#    avoiding subprocess + conda activation overhead per call.

_equibind_dir = str(EQUIBIND_DIR)
if _equibind_dir not in sys.path:
    sys.path.insert(0, _equibind_dir)

from models.equibind import EquiBind as EquiBindModel
from commons.process_mols import (
    get_receptor_inference,
    get_rec_graph,
    get_lig_graph_revised,
    get_geometry_graph,
)
from commons.geometry_utils import (
    rigid_transform_Kabsch_3D,
    get_torsions,
    get_dihedral_vonMises,
    apply_changes,
)
from commons.utils import seed_all


# ── 0a. Load EquiBind model ONCE ──────────────────────────────────────────

def _load_equibind_model(equibind_dir: Path, device: str = "cuda"):
    """Load EquiBind model and train args. Called once at startup."""
    ckpt_path = equibind_dir / "runs" / "flexible_self_docking" / "best_checkpoint.pt"
    train_args_path = ckpt_path.parent / "train_arguments.yaml"

    with open(train_args_path) as f:
        train_args = yaml.safe_load(f)

    # Critical: set noise_initial=0 for inference
    train_args["model_parameters"]["noise_initial"] = 0

    args = SimpleNamespace(**train_args)
    args.checkpoint = str(ckpt_path)
    args.use_rdkit_coords = args.dataset_params.get("use_rdkit_coords", True)
    args.device = device

    dev = torch.device("cuda:0" if torch.cuda.is_available() and device == "cuda" else "cpu")
    checkpoint = torch.load(str(ckpt_path), map_location=dev)

    model = EquiBindModel(
        device=dev,
        lig_input_edge_feats_dim=15,
        rec_input_edge_feats_dim=27,
        **args.model_parameters,
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(dev)
    model.eval()

    seed_all(args.seed)

    monitor.info(f"EquiBind model loaded from {ckpt_path}")
    monitor.info(f"  Device: {dev}")
    return model, args, dev


def _load_receptor_graph(protein_pdb: Path, args: SimpleNamespace):
    """Build the EquiBind receptor graph for a protein. Cached per protein."""
    dp = args.dataset_params
    rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(str(protein_pdb))
    rec_graph = get_rec_graph(
        rec, rec_coords, c_alpha_coords, n_coords, c_coords,
        use_rec_atoms=dp["use_rec_atoms"],
        rec_radius=dp["rec_graph_radius"],
        surface_max_neighbors=dp["surface_max_neighbors"],
        surface_graph_cutoff=dp["surface_graph_cutoff"],
        surface_mesh_cutoff=dp["surface_mesh_cutoff"],
        c_alpha_max_neighbors=dp["c_alpha_max_neighbors"],
    )
    return rec_graph


# Load model once now
_eb_model, _eb_args, _eb_device = _load_equibind_model(EQUIBIND_DIR, EQUIBIND_DEVICE)

# Cache: protein_path_str → rec_graph
_rec_graph_cache: Dict[str, object] = {}


# ── 0b. In-process single-ligand docking ──────────────────────────────────

def _run_corrections_inproc(lig, lig_coord, predicted_coords):
    """Post-process: torsion fitting + Kabsch alignment (identical to EquiBind's run_corrections)."""
    input_coords = lig_coord.detach().cpu()
    prediction = predicted_coords.detach().cpu()

    lig_input = deepcopy(lig)
    conf = lig_input.GetConformer()
    for i in range(lig_input.GetNumAtoms()):
        x, y, z = input_coords.numpy()[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    lig_equibind = deepcopy(lig)
    conf = lig_equibind.GetConformer()
    for i in range(lig_equibind.GetNumAtoms()):
        x, y, z = prediction.numpy()[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    coords_pred = lig_equibind.GetConformer().GetPositions()
    Z_pt_cloud = coords_pred
    rotable_bonds = get_torsions([lig_input])
    new_dihedrals = np.zeros(len(rotable_bonds))
    for idx_t, r in enumerate(rotable_bonds):
        new_dihedrals[idx_t] = get_dihedral_vonMises(lig_input, lig_input.GetConformer(), r, Z_pt_cloud)
    optimized_mol = apply_changes(lig_input, new_dihedrals, rotable_bonds)
    optimized_conf = optimized_mol.GetConformer()
    coords_pred_optimized = optimized_conf.GetPositions()
    R, t = rigid_transform_Kabsch_3D(coords_pred_optimized.T, coords_pred.T)
    coords_pred_optimized = (R @ coords_pred_optimized.T).T + t.squeeze()
    for i in range(optimized_mol.GetNumAtoms()):
        x, y, z = coords_pred_optimized[i]
        optimized_conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    return optimized_mol


def _dock_single_ligand_inproc(
    mol,
    rec_graph,
    model,
    args: SimpleNamespace,
    device: torch.device,
) -> Optional[Chem.Mol]:
    """Dock a single RDKit Mol in-process. Returns optimized Mol or None."""
    dp = args.dataset_params
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else "unnamed"

    try:
        lig_graph = get_lig_graph_revised(
            mol, name,
            max_neighbors=dp["lig_max_neighbors"],
            use_rdkit_coords=args.use_rdkit_coords,
            radius=dp["lig_graph_radius"],
        )
    except Exception as e:
        monitor.warning(f"Ligand graph build failed for {name}: {e}")
        return None

    # Ensure new_x is set (required by model forward)
    lig_graph.ndata["new_x"] = lig_graph.ndata["x"]

    geometry_graph = get_geometry_graph(mol) if dp.get("geometry_regularization", True) else None

    # Save input coords for post-processing
    lig_coord = lig_graph.ndata["new_x"].clone()

    # Move graphs to device
    lig_graph = lig_graph.to(device)
    rec_graph_dev = rec_graph.to(device)
    geom_graph_dev = geometry_graph.to(device) if geometry_graph is not None else None

    try:
        with torch.no_grad():
            predictions = model(lig_graph, rec_graph_dev, geom_graph_dev)
        predicted_coords = predictions[0][0]  # (n_atoms, 3)
    except Exception as e:
        monitor.warning(f"EquiBind forward failed for {name}: {e}")
        return None

    try:
        optimized_mol = _run_corrections_inproc(mol, lig_coord, predicted_coords)
        return optimized_mol
    except Exception as e:
        monitor.warning(f"Post-processing failed for {name}: {e}")
        # Fallback: return mol with raw predicted coords
        out_mol = deepcopy(mol)
        conf = out_mol.GetConformer()
        coords_np = predicted_coords.detach().cpu().numpy()
        for i in range(out_mol.GetNumAtoms()):
            conf.SetAtomPosition(i, Point3D(*coords_np[i].tolist()))
        return out_mol


# ── 1. Parse fpocket results ──────────────────────────────────────────────

@dataclass
class PocketInfo:
    """Generic pocket descriptor from fpocket or p2rank."""
    source: str            # "fpocket" or "p2rank"
    pocket_id: int         # rank / pocket number (original)
    unique_id: str         # globally unique ID (e.g. "fpocket_cleaned_p001")
    score: float
    center: Tuple[float, float, float]
    radius: float          # rough extent
    protein_name: str
    extra: dict = field(default_factory=dict)

    def to_dict(self):
        return {
            "source": self.source,
            "pocket_id": self.pocket_id,
            "unique_id": self.unique_id,
            "score": self.score,
            "center": list(self.center),
            "radius": self.radius,
            "protein_name": self.protein_name,
            "extra": self.extra,
        }


def _parse_fpocket_pocket_centroid(pocket_pdb: Path) -> Tuple[Tuple[float, float, float], float]:
    """Compute centroid and approximate radius from fpocket pocket_*_atm.pdb."""
    coords = []
    with open(pocket_pdb, "r") as fh:
        for line in fh:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except (ValueError, IndexError):
                    continue
    if not coords:
        return (0.0, 0.0, 0.0), 0.0
    arr = np.array(coords)
    centroid = tuple(arr.mean(axis=0))
    dists = np.linalg.norm(arr - arr.mean(axis=0), axis=1)
    radius = float(dists.max())
    return centroid, radius


def parse_fpocket_results(fpocket_dir: Path, protein_name: str, n_top: int = N_TOP_POCKETS) -> List[PocketInfo]:
    """Parse top-n fpocket pockets for a given protein.
    
    Uses unique IDs that include the directory variant (e.g. 'cleaned' vs raw)
    to avoid collisions when multiple _out directories match.
    """
    out_dirs = sorted(fpocket_dir.glob(f"{protein_name}*_out"))
    pockets: List[PocketInfo] = []
    for out_dir in out_dirs:
        pockets_dir = out_dir / "pockets"
        if not pockets_dir.exists():
            continue

        # Derive a short variant tag from the directory name
        dir_suffix = out_dir.name[len(protein_name):]
        variant = dir_suffix.replace("_out", "").strip("_") or "raw"

        # Parse info file for scores
        info_file = list(out_dir.glob("*_info.txt"))
        scores: Dict[int, float] = {}
        if info_file:
            with open(info_file[0]) as fh:
                current_pocket = None
                for line in fh:
                    line_s = line.strip()
                    if line_s.startswith("Pocket ") and ":" in line_s:
                        try:
                            current_pocket = int(line_s.split()[1])
                        except (ValueError, IndexError):
                            current_pocket = None
                    elif current_pocket is not None and line_s.startswith("Score"):
                        try:
                            scores[current_pocket] = float(line_s.split(":")[-1].strip())
                        except ValueError:
                            pass

        # Iterate over pocket PDB files
        pocket_pdbs = sorted(pockets_dir.glob("pocket*_atm.pdb"))
        for ppdb in pocket_pdbs:
            try:
                pocket_num = int(ppdb.stem.replace("pocket", "").replace("_atm", ""))
            except ValueError:
                continue
            centroid, radius = _parse_fpocket_pocket_centroid(ppdb)
            unique_id = f"fpocket_{variant}_p{pocket_num:03d}"
            pockets.append(PocketInfo(
                source="fpocket",
                pocket_id=pocket_num,
                unique_id=unique_id,
                score=scores.get(pocket_num, 0.0),
                center=centroid,
                radius=radius,
                protein_name=protein_name,
                extra={"pdb_file": str(ppdb), "out_dir": str(out_dir), "variant": variant},
            ))

    pockets.sort(key=lambda p: p.score, reverse=True)
    return pockets[:n_top]


# ── 2. Parse p2rank results ──────────────────────────────────────────────

def parse_p2rank_results(p2rank_dir: Path, protein_name: str, n_top: int = N_TOP_POCKETS) -> List[PocketInfo]:
    """Parse top-n p2rank pockets for a given protein from predictions CSV."""
    pockets: List[PocketInfo] = []
    pred_files = sorted(p2rank_dir.glob(f"{protein_name}*_predictions.csv"))
    for pred_file in pred_files:
        fname_stem = pred_file.stem
        base = fname_stem.replace("_predictions", "")
        suffix_part = base[len(protein_name):].replace(".pdb", "").strip("_")
        variant = suffix_part or "raw"

        with open(pred_file) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                try:
                    def _g(key):
                        for k, v in row.items():
                            if k.strip() == key:
                                return v.strip()
                        return None

                    rank = int(_g("rank"))
                    score = float(_g("score"))
                    prob = float(_g("probability"))
                    cx = float(_g("center_x"))
                    cy = float(_g("center_y"))
                    cz = float(_g("center_z"))
                    name = (_g("name") or "").strip()
                    sas = float(_g("sas_points") or 0)
                except (TypeError, ValueError):
                    continue

                radius = max(5.0, np.sqrt(sas) * 0.5)
                unique_id = f"p2rank_{variant}_p{rank:03d}"
                pockets.append(PocketInfo(
                    source="p2rank",
                    pocket_id=rank,
                    unique_id=unique_id,
                    score=score,
                    center=(cx, cy, cz),
                    radius=radius,
                    protein_name=protein_name,
                    extra={
                        "probability": prob,
                        "sas_points": sas,
                        "pred_file": str(pred_file),
                        "name": name,
                        "variant": variant,
                    },
                ))

    pockets.sort(key=lambda p: p.score, reverse=True)
    return pockets[:n_top]


# ── 3. Pocket-guided EquiBind docking (IN-PROCESS) ───────────────────────

def _translate_sdf_to_pocket(sdf_path: Path, pocket_center: Tuple[float, float, float],
                              output_path: Path) -> bool:
    """Translate a ligand SDF so its centroid is at the pocket center."""
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            mol = Chem.MolFromMol2File(str(sdf_path), removeHs=False)
        if mol is None:
            return False

        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        positions = conf.GetPositions()
        centroid = positions.mean(axis=0)
        shift = np.array(pocket_center) - centroid

        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (pos.x + shift[0], pos.y + shift[1], pos.z + shift[2]))

        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        return output_path.exists() and output_path.stat().st_size > 0
    except Exception as e:
        monitor.warning(f"Could not translate ligand to pocket: {e}")
        return False


def _run_equibind(protein_pdb: Path, ligand_file: Path, output_dir: Path,
                  seed: int = 42, device: str = "cpu") -> Tuple[bool, Optional[Path], str]:
    """Run EquiBind IN-PROCESS (no subprocess) for a single ligand.

    Uses the globally loaded model and caches receptor graphs per protein.
    The model, args, and device are stored in module-level variables
    _eb_model, _eb_args, _eb_device set at cell startup.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Monitor: log the outgoing call ──
    monitor.equibind_call(
        protein=protein_pdb.name,
        ligand=ligand_file.name,
        output_dir=str(output_dir),
        seed=seed,
        device=device,
    )

    # ── Get or cache receptor graph ──
    protein_key = str(protein_pdb)
    if protein_key not in _rec_graph_cache:
        monitor.info(f"Building receptor graph for {protein_pdb.name} (cached for reuse)")
        try:
            _rec_graph_cache[protein_key] = _load_receptor_graph(protein_pdb, _eb_args)
        except Exception as e:
            err_msg = f"Receptor graph build failed: {e}"
            monitor.equibind_return(False, None, err_msg)
            return False, None, err_msg

    rec_graph = _rec_graph_cache[protein_key]

    # ── Load ligand mol ──
    try:
        suffix = ligand_file.suffix.lower()
        if suffix == ".sdf":
            suppl = Chem.SDMolSupplier(str(ligand_file), sanitize=False, removeHs=False)
            mol = next(iter(suppl), None)
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(ligand_file), removeHs=False)
        elif suffix == ".pdb":
            mol = Chem.MolFromPDBFile(str(ligand_file), removeHs=False)
        else:
            mol = None

        if mol is None:
            err_msg = f"Could not load ligand from {ligand_file.name}"
            monitor.equibind_return(False, None, err_msg)
            return False, None, err_msg

        Chem.SanitizeMol(mol)
        if not mol.HasProp("_Name"):
            mol.SetProp("_Name", ligand_file.stem)
    except Exception as e:
        err_msg = f"Ligand loading failed: {e}"
        monitor.equibind_return(False, None, err_msg)
        return False, None, err_msg

    # ── Seed for this call ──
    seed_all(seed)

    # ── Run in-process docking ──
    docked_mol = _dock_single_ligand_inproc(mol, rec_graph, _eb_model, _eb_args, _eb_device)

    if docked_mol is None:
        err_msg = "EquiBind inference returned None"
        monitor.equibind_return(False, None, err_msg)
        return False, None, err_msg

    # ── Write output SDF ──
    output_sdf = output_dir / "output.sdf"
    try:
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(docked_mol)
        writer.close()
    except Exception as e:
        err_msg = f"Failed to write output SDF: {e}"
        monitor.equibind_return(False, None, err_msg)
        return False, None, err_msg

    if output_sdf.exists() and output_sdf.stat().st_size > 0:
        monitor.equibind_return(True, str(output_sdf))
        return True, output_sdf, ""
    else:
        err_msg = "Output SDF empty or missing"
        monitor.equibind_return(False, None, err_msg)
        return False, None, err_msg


def _sdf_centroid(sdf_path: Path) -> Optional[Tuple[float, float, float]]:
    """Compute centroid from SDF file."""
    try:
        coords = []
        with open(sdf_path) as f:
            lines = f.readlines()
        if len(lines) < 5:
            return None
        parts = lines[3].split()
        n_atoms = int(parts[0])
        for i in range(4, min(4 + n_atoms, len(lines))):
            p = lines[i].split()
            if len(p) >= 3:
                coords.append([float(p[0]), float(p[1]), float(p[2])])
        if not coords:
            return None
        return tuple(np.array(coords).mean(axis=0))
    except Exception:
        return None


def _pocket_distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """Euclidean distance between two 3D points."""
    return float(np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2))


def _generate_conformer_sdfs(ligand_file: Path, output_dir: Path,
                              seeds: list = None, per_seed: int = None) -> List[Tuple[Path, int]]:
    """Generate diverse RDKit conformers, return list of (sdf_path, seed)."""
    if seeds is None:
        seeds = RDKIT_SEEDS
    if per_seed is None:
        per_seed = CONFORMERS_PER_SEED

    output_dir.mkdir(parents=True, exist_ok=True)

    suffix = ligand_file.suffix.lower()
    try:
        if suffix == ".sdf":
            suppl = Chem.SDMolSupplier(str(ligand_file), removeHs=False)
            mol = next(iter(suppl), None)
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(ligand_file), removeHs=False)
        elif suffix == ".pdb":
            mol = Chem.MolFromPDBFile(str(ligand_file), removeHs=False)
        else:
            mol = None
    except Exception:
        mol = None

    if mol is None:
        monitor.warning(f"Could not load molecule from {ligand_file.name}")
        return []

    mol = Chem.AddHs(mol)
    conformer_files: List[Tuple[Path, int]] = []

    for seed in seeds:
        params = AllChem.ETKDGv3()
        params.randomSeed = seed
        params.numThreads = 0
        params.useRandomCoords = True
        params.maxIterations = 500
        try:
            cids = AllChem.EmbedMultipleConfs(mol, numConfs=per_seed, params=params)
        except Exception:
            cids = []
        for ci, cid in enumerate(cids):
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
            except Exception:
                pass
            cpath = output_dir / f"seed{seed}_c{ci}.sdf"
            try:
                w = Chem.SDWriter(str(cpath))
                w.write(mol, confId=cid)
                w.close()
                if cpath.exists() and cpath.stat().st_size > 0:
                    conformer_files.append((cpath, seed))
            except Exception:
                pass

    monitor.info(f"RDKit generated {len(conformer_files)} conformers for {ligand_file.name}")
    return conformer_files


@dataclass
class GuidedPoseResult:
    """One docked pose from pocket-guided or unguided EquiBind."""
    protein_name: str
    ligand_name: str
    mode: str              # "fpocket", "p2rank", "unguided"
    pocket_id: Optional[int]
    pocket_unique_id: Optional[str]
    pose_num: int          # pose number within this pocket (1..POSES_PER_POCKET)
    pocket_center: Optional[Tuple[float, float, float]]
    pose_centroid: Optional[Tuple[float, float, float]]
    sdf_path: Optional[Path]
    success: bool
    error: str = ""
    uff_energy_before: Optional[float] = None
    uff_energy_after: Optional[float] = None
    uff_minimized: bool = False

    def to_dict(self):
        return {
            "protein_name": self.protein_name,
            "ligand_name": self.ligand_name,
            "mode": self.mode,
            "pocket_id": self.pocket_id,
            "pocket_unique_id": self.pocket_unique_id,
            "pose_num": self.pose_num,
            "pocket_center": list(self.pocket_center) if self.pocket_center else None,
            "pose_centroid": list(self.pose_centroid) if self.pose_centroid else None,
            "sdf_path": str(self.sdf_path) if self.sdf_path else None,
            "success": self.success,
            "error": self.error,
            "uff_minimized": self.uff_minimized,
            "uff_energy_before": self.uff_energy_before,
            "uff_energy_after": self.uff_energy_after,
        }


def _ensure_ligand_sdf(ligand_file: Path, prep_dir: Path) -> Path:
    """Ensure ligand is available as SDF; convert if needed."""
    ligand_name = ligand_file.stem
    lig_sdf = prep_dir / f"{ligand_name}.sdf"
    if lig_sdf.exists():
        return lig_sdf

    suffix = ligand_file.suffix.lower()
    if suffix == ".sdf":
        shutil.copy2(ligand_file, lig_sdf)
    elif suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(ligand_file), removeHs=False)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            w = Chem.SDWriter(str(lig_sdf)); w.write(mol); w.close()
        else:
            shutil.copy2(ligand_file, lig_sdf)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(ligand_file), removeHs=False)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            w = Chem.SDWriter(str(lig_sdf)); w.write(mol); w.close()
        else:
            shutil.copy2(ligand_file, lig_sdf)
    monitor.info(f"Prepared ligand SDF: {lig_sdf.name}")
    return lig_sdf


def _dock_pocket_multi_pose(
    pocket: PocketInfo,
    prepared_protein: Path,
    lig_sdf: Path,
    conformer_files: List[Tuple[Path, int]],
    combo_dir: Path,
    prep_dir: Path,
    protein_name: str,
    ligand_name: str,
    poses_per_pocket: int,
    skip_existing: bool,
    force_pocket: bool = FORCE_POCKET,
    max_pocket_trials: int = MAX_POCKET_TRIALS,
    pocket_threshold: float = POCKET_MATCH_THRESHOLD,
    use_cropping: bool = USE_PROTEIN_CROPPING,
) -> List[GuidedPoseResult]:
    """Dock N poses into a single pocket using different conformers.

    If use_cropping is True, the protein is cropped to only include residues
    near the pocket center, physically constraining EquiBind's search space.
    
    If force_pocket is True, each docked pose is checked: if the resulting
    centroid is farther than pocket_threshold from the pocket center, the
    pose is rejected and the next conformer is tried, up to
    max_pocket_trials attempts per desired pose.
    """
    mode = pocket.source  # "fpocket" or "p2rank"
    uid = pocket.unique_id
    results: List[GuidedPoseResult] = []
    
    # ── Protein cropping for pocket-constrained docking ──
    cropped_protein = None
    crop_offset = np.zeros(3)
    if use_cropping:
        cropped_protein, crop_offset = get_cropped_protein_for_pocket(
            protein_pdb=prepared_protein,
            pocket=pocket,
            prep_dir=prep_dir,
            crop_radius=POCKET_CROP_RADIUS,
        )
        if cropped_protein is not None:
            # Build receptor graph for cropped protein (cached separately)
            cropped_key = str(cropped_protein)
            if cropped_key not in _rec_graph_cache:
                monitor.info(f"Building receptor graph for cropped protein {cropped_protein.name}")
                _rec_graph_cache[cropped_key] = _load_receptor_graph(cropped_protein, _eb_args)
            monitor.info(f"Using cropped protein for {uid} (pocket at origin)")
        else:
            monitor.warning(f"Protein cropping failed for {uid}, falling back to full protein")

    # Check how many poses already exist for this pocket
    existing_count = 0
    for pose_num in range(1, poses_per_pocket + 1):
        final_sdf = combo_dir / f"{uid}_pose{pose_num:02d}.sdf"
        if final_sdf.exists() and skip_existing:
            centroid = _sdf_centroid(final_sdf)
            results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode=mode, pocket_id=pocket.pocket_id,
                pocket_unique_id=uid, pose_num=pose_num,
                pocket_center=pocket.center, pose_centroid=centroid,
                sdf_path=final_sdf, success=True,
            ))
            existing_count += 1

    if existing_count >= poses_per_pocket:
        monitor.info(f"{uid}: all {poses_per_pocket} poses already cached, skipping")
        return results

    # Generate poses using conformers translated to pocket center
    poses_done = existing_count
    conf_idx = 0

    while poses_done < poses_per_pocket and conf_idx < len(conformer_files):
        pose_num = poses_done + 1
        final_sdf = combo_dir / f"{uid}_pose{pose_num:02d}.sdf"
        if final_sdf.exists() and skip_existing:
            conf_idx += 1
            continue

        # Try up to max_pocket_trials conformers for this pose slot
        trials_for_this_pose = 0
        accepted = False

        while trials_for_this_pose < max_pocket_trials and conf_idx < len(conformer_files):
            conf_path, seed = conformer_files[conf_idx]
            conf_idx += 1
            trials_for_this_pose += 1

            # Translate this conformer to the target location:
            # - If using cropped protein: translate to ORIGIN (pocket is at origin in cropped frame)
            # - Otherwise: translate to pocket center in original frame
            translated = prep_dir / f"lig_{uid}_pose{pose_num:02d}_t{trials_for_this_pose:02d}.sdf"
            if cropped_protein is not None:
                # Cropped protein has pocket at origin, so translate ligand to origin
                target_center = (0.0, 0.0, 0.0)
            else:
                target_center = pocket.center
                
            ok = _translate_sdf_to_pocket(conf_path, target_center, translated)
            if not ok:
                ok = _translate_sdf_to_pocket(lig_sdf, target_center, translated)
            if not ok:
                monitor.warning(f"{uid} pose {pose_num} trial {trials_for_this_pose}: "
                                f"ligand translation failed, skipping conformer")
                continue

            # Use cropped protein if available, otherwise full protein
            dock_protein = cropped_protein if cropped_protein is not None else prepared_protein
            
            pose_dir = combo_dir / f"{uid}_run{pose_num:02d}_t{trials_for_this_pose:02d}"
            success, out_sdf, err = _run_equibind(dock_protein, translated, pose_dir,
                                                   seed=seed, device=EQUIBIND_DEVICE)

            if not success or out_sdf is None:
                # Clean up failed run directory
                try:
                    shutil.rmtree(pose_dir, ignore_errors=True)
                except Exception:
                    pass
                monitor.warning(f"{uid} pose {pose_num} trial {trials_for_this_pose}: DOCK FAIL")
                continue

            # If using cropped protein, translate pose back to original coordinate frame
            if cropped_protein is not None and np.any(crop_offset != 0):
                restored_sdf = pose_dir / "output_restored.sdf"
                if translate_pose_back_to_original_frame(out_sdf, crop_offset, restored_sdf):
                    out_sdf = restored_sdf
                else:
                    monitor.warning(f"{uid} pose {pose_num}: failed to restore coordinates")

            # Compute centroid of the docked pose (now in original frame)
            centroid = _sdf_centroid(out_sdf)

            # ── Pocket enforcement with optional clamping ──
            if force_pocket and centroid is not None:
                dist = _pocket_distance(centroid, pocket.center)
                if dist > pocket_threshold:
                    pose_label = f"{uid}_pose{pose_num:02d}"
                    
                    # Option 1: Clamp pose back to pocket (if enabled)
                    if CLAMP_POSE_TO_POCKET:
                        clamped_sdf = pose_dir / "output_clamped.sdf"
                        clamp_ok, orig_dist, new_dist = clamp_pose_to_pocket(
                            sdf_path=out_sdf,
                            pocket_center=pocket.center,
                            max_distance=pocket_threshold,
                            output_path=clamped_sdf,
                        )
                        if clamp_ok:
                            out_sdf = clamped_sdf
                            centroid = _sdf_centroid(out_sdf)
                            monitor.info(f"[CLAMP] {pose_label} clamped from {orig_dist:.1f}Å → {new_dist:.1f}Å")
                        else:
                            monitor.warning(f"{pose_label}: clamping failed, rejecting pose")
                            try:
                                shutil.rmtree(pose_dir, ignore_errors=True)
                            except Exception:
                                pass
                            continue  # try next conformer
                    else:
                        # Option 2: Reject pose (original behavior)
                        monitor.pose_rejected_from_pocket(
                            pose_id=pose_label, pocket_id=uid,
                            distance=dist, threshold=pocket_threshold,
                        )
                        # Clean up rejected pose directory
                        try:
                            shutil.rmtree(pose_dir, ignore_errors=True)
                        except Exception:
                            pass
                        continue  # try next conformer

            # ── UFF post-docking minimization ──
            _uff_minimized = False
            _uff_e_before = None
            _uff_e_after = None
            if UFF_MINIMIZE:
                minimized_sdf = pose_dir / "output_uff_minimized.sdf"
                uff_ok, e_before, e_after, uff_msg = _uff_minimize_pose(
                    docked_sdf=out_sdf,
                    protein_pdb=prepared_protein,
                    output_sdf=minimized_sdf,
                )
                if uff_ok and minimized_sdf.exists():
                    out_sdf = minimized_sdf
                    _uff_minimized = True
                    _uff_e_before = e_before
                    _uff_e_after = e_after
                    monitor.info(f"{uid} pose {pose_num}: UFF minimized "
                                 f"(E {e_before:.1f}→{e_after:.1f}, {uff_msg})")
                else:
                    monitor.warning(f"{uid} pose {pose_num}: UFF minimization failed "
                                    f"({uff_msg}), using raw docked pose")

            # Pose accepted
            shutil.copy2(out_sdf, final_sdf)
            # Clean up run directory after successful copy
            try:
                shutil.rmtree(pose_dir, ignore_errors=True)
            except Exception:
                pass
            centroid = _sdf_centroid(final_sdf)  # re-read from final location
            results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode=mode, pocket_id=pocket.pocket_id,
                pocket_unique_id=uid, pose_num=pose_num,
                pocket_center=pocket.center, pose_centroid=centroid,
                sdf_path=final_sdf, success=True,
                uff_minimized=_uff_minimized,
                uff_energy_before=_uff_e_before,
                uff_energy_after=_uff_e_after,
            ))
            poses_done += 1
            if force_pocket and centroid:
                dist = _pocket_distance(centroid, pocket.center)
                pose_label = f"{uid}_pose{pose_num:02d}"
                monitor.pose_accepted_in_pocket(
                    pose_id=pose_label, pocket_id=uid,
                    distance=dist, threshold=pocket_threshold,
                )
                monitor.info(f"{uid} pose {pose_num}/{poses_per_pocket}: OK "
                             f"(trial {trials_for_this_pose}/{max_pocket_trials})")
            else:
                monitor.info(f"{uid} pose {pose_num}/{poses_per_pocket}: OK")
            accepted = True
            break  # move to next pose slot

        if not accepted:
            error_msg = (f"No pose stayed in pocket after {trials_for_this_pose} trials"
                         if force_pocket else "All conformers exhausted")
            results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode=mode, pocket_id=pocket.pocket_id,
                pocket_unique_id=uid, pose_num=pose_num,
                pocket_center=pocket.center, pose_centroid=None,
                sdf_path=None, success=False,
                error=error_msg,
            ))
            monitor.critical(f"{uid} pose {pose_num}/{poses_per_pocket}: EXHAUSTED "
                             f"({trials_for_this_pose} trials, none in pocket)")

    return results


def dock_guided_by_pockets(
    protein_pdb: Path,
    ligand_file: Path,
    fpocket_pockets: List[PocketInfo],
    p2rank_pockets: List[PocketInfo],
    poses_per_pocket: int = POSES_PER_POCKET,
    n_unguided: int = N_UNGUIDED_POSES,
    skip_existing: bool = SKIP_EXISTING,
    force_pocket: bool = FORCE_POCKET,
    max_pocket_trials: int = MAX_POCKET_TRIALS,
) -> List[GuidedPoseResult]:
    """
    Dock one protein-ligand pair:
      1. For each fpocket pocket: generate N poses (different conformers)
      2. For each p2rank pocket: same
      3. Generate n_unguided natural poses (different conformers, no translation)

    If force_pocket=True, docked poses whose centroid drifts beyond
    POCKET_MATCH_THRESHOLD from the pocket center are rejected, and
    up to max_pocket_trials conformers are tried per pose slot.

    Returns list of GuidedPoseResult.
    """
    protein_name = protein_pdb.stem
    ligand_name = ligand_file.stem
    combo = f"{ligand_name}__{protein_name}"

    combo_dir = POCKET_GUIDED_OUTPUT_DIR / combo
    combo_dir.mkdir(parents=True, exist_ok=True)

    results_file = combo_dir / "guided_results.json"

    # Check for cached results
    if skip_existing and results_file.exists():
        try:
            cached = json.loads(results_file.read_text())
            expected_total = (len(fpocket_pockets) + len(p2rank_pockets)) * poses_per_pocket + n_unguided
            cached_results = []
            for r in cached:
                cached_results.append(GuidedPoseResult(
                    protein_name=r["protein_name"],
                    ligand_name=r["ligand_name"],
                    mode=r["mode"],
                    pocket_id=r.get("pocket_id"),
                    pocket_unique_id=r.get("pocket_unique_id"),
                    pose_num=r.get("pose_num", 1),
                    pocket_center=tuple(r["pocket_center"]) if r.get("pocket_center") else None,
                    pose_centroid=tuple(r["pose_centroid"]) if r.get("pose_centroid") else None,
                    sdf_path=Path(r["sdf_path"]) if r.get("sdf_path") else None,
                    success=r["success"],
                    error=r.get("error", ""),
                ))
            if len(cached_results) >= expected_total:
                monitor.info(f"[SKIP] {combo}: loaded {len(cached_results)} cached results")
                return cached_results
        except Exception:
            pass

    # Prepare protein (copy + reduce)
    prep_dir = combo_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    prepared_protein = prep_dir / f"{protein_name}_protein.pdb"
    if not prepared_protein.exists():
        reduce_exe = shutil.which("reduce")
        if reduce_exe:
            try:
                res = subprocess.run([reduce_exe, "-Quiet", str(protein_pdb)],
                                     capture_output=True, text=True, timeout=60)
                if res.returncode == 0 and res.stdout.strip():
                    prepared_protein.write_text(res.stdout)
                    monitor.info(f"Protein prepared with reduce: {prepared_protein.name}")
                else:
                    shutil.copy2(protein_pdb, prepared_protein)
                    monitor.info(f"Protein copied (reduce returned empty): {prepared_protein.name}")
            except Exception:
                shutil.copy2(protein_pdb, prepared_protein)
        else:
            shutil.copy2(protein_pdb, prepared_protein)
            monitor.info(f"Protein copied (reduce not available): {prepared_protein.name}")

    # Pre-cache the receptor graph for this protein (done once, reused for all pockets/ligands)
    protein_key = str(prepared_protein)
    if protein_key not in _rec_graph_cache:
        monitor.info(f"Pre-caching receptor graph for {prepared_protein.name} ...")
        _rec_graph_cache[protein_key] = _load_receptor_graph(prepared_protein, _eb_args)
        monitor.info(f"Receptor graph cached for {prepared_protein.name}")

    # Ensure ligand SDF
    lig_sdf = _ensure_ligand_sdf(ligand_file, prep_dir)

    # Generate conformers (shared across all pockets for this combo)
    conf_dir = prep_dir / "conformers"
    conformer_files = _generate_conformer_sdfs(ligand_file, conf_dir, RDKIT_SEEDS, CONFORMERS_PER_SEED)
    if not conformer_files:
        conformer_files = [(lig_sdf, 42)]
        monitor.warning(f"No RDKit conformers generated — using original ligand SDF as fallback")
    monitor.info(f"Using {len(conformer_files)} conformers for multi-pose docking")
    if force_pocket:
        monitor.info(f"Pocket enforcement ON: max {max_pocket_trials} trials/pose, "
                     f"threshold {POCKET_MATCH_THRESHOLD}Å")
    if USE_PROTEIN_CROPPING:
        monitor.info(f"Protein cropping ON: {POCKET_CROP_RADIUS}Å radius + {POCKET_CROP_BUFFER}Å buffer")

    all_results: List[GuidedPoseResult] = []

    # ── A. Guided by fpocket pockets ──
    monitor.section(f"fpocket-guided  ({len(fpocket_pockets)} pockets × {poses_per_pocket} poses "
                    f"= {len(fpocket_pockets) * poses_per_pocket} total)")
    for pocket in fpocket_pockets:
        monitor.info(f"Pocket {pocket.unique_id}  score={pocket.score:.2f}  "
                     f"center=({pocket.center[0]:.1f}, {pocket.center[1]:.1f}, {pocket.center[2]:.1f})  "
                     f"radius={pocket.radius:.1f}Å")
        pocket_results = _dock_pocket_multi_pose(
            pocket=pocket,
            prepared_protein=prepared_protein,
            lig_sdf=lig_sdf,
            conformer_files=conformer_files,
            combo_dir=combo_dir,
            prep_dir=prep_dir,
            protein_name=protein_name,
            ligand_name=ligand_name,
            poses_per_pocket=poses_per_pocket,
            skip_existing=skip_existing,
            force_pocket=force_pocket,
            max_pocket_trials=max_pocket_trials,
        )
        ok_count = sum(1 for r in pocket_results if r.success)
        fail_count = sum(1 for r in pocket_results if not r.success)
        monitor.info(f"  → {ok_count}/{poses_per_pocket} poses OK, {fail_count} failed")
        if fail_count > 0:
            monitor.warning(f"  {pocket.unique_id}: {fail_count} pose(s) could not be placed in pocket")
        all_results.extend(pocket_results)

    # ── B. Guided by p2rank pockets ──
    monitor.section(f"p2rank-guided  ({len(p2rank_pockets)} pockets × {poses_per_pocket} poses "
                    f"= {len(p2rank_pockets) * poses_per_pocket} total)")
    for pocket in p2rank_pockets:
        monitor.info(f"Pocket {pocket.unique_id}  score={pocket.score:.2f}  "
                     f"center=({pocket.center[0]:.1f}, {pocket.center[1]:.1f}, {pocket.center[2]:.1f})  "
                     f"radius={pocket.radius:.1f}Å")
        pocket_results = _dock_pocket_multi_pose(
            pocket=pocket,
            prepared_protein=prepared_protein,
            lig_sdf=lig_sdf,
            conformer_files=conformer_files,
            combo_dir=combo_dir,
            prep_dir=prep_dir,
            protein_name=protein_name,
            ligand_name=ligand_name,
            poses_per_pocket=poses_per_pocket,
            skip_existing=skip_existing,
            force_pocket=force_pocket,
            max_pocket_trials=max_pocket_trials,
        )
        ok_count = sum(1 for r in pocket_results if r.success)
        fail_count = sum(1 for r in pocket_results if not r.success)
        monitor.info(f"  → {ok_count}/{poses_per_pocket} poses OK, {fail_count} failed")
        if fail_count > 0:
            monitor.warning(f"  {pocket.unique_id}: {fail_count} pose(s) could not be placed in pocket")
        all_results.extend(pocket_results)

    # ── C. Unguided (natural) EquiBind poses ──
    monitor.section(f"Unguided EquiBind poses  ({n_unguided} requested)")

    unguided_count = 0
    for conf_path, seed in conformer_files:
        if unguided_count >= n_unguided:
            break

        pose_num = unguided_count + 1
        final_sdf = combo_dir / f"unguided_{pose_num:03d}.sdf"

        if final_sdf.exists() and skip_existing:
            centroid = _sdf_centroid(final_sdf)
            all_results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode="unguided", pocket_id=None,
                pocket_unique_id=None, pose_num=pose_num,
                pocket_center=None, pose_centroid=centroid,
                sdf_path=final_sdf, success=True,
            ))
            unguided_count += 1
            monitor.info(f"Unguided pose {pose_num}/{n_unguided}: cached")
            continue

        pose_dir = combo_dir / f"unguided_run_{pose_num:03d}"
        success, out_sdf, err = _run_equibind(prepared_protein, conf_path, pose_dir,
                                               seed=seed, device=EQUIBIND_DEVICE)
        if success and out_sdf:
            # ── UFF post-docking minimization for unguided poses ──
            _uff_minimized = False
            _uff_e_before = None
            _uff_e_after = None
            if UFF_MINIMIZE:
                minimized_sdf = pose_dir / "output_uff_minimized.sdf"
                uff_ok, e_before, e_after, uff_msg = _uff_minimize_pose(
                    docked_sdf=out_sdf,
                    protein_pdb=prepared_protein,
                    output_sdf=minimized_sdf,
                )
                if uff_ok and minimized_sdf.exists():
                    out_sdf = minimized_sdf
                    _uff_minimized = True
                    _uff_e_before = e_before
                    _uff_e_after = e_after
                    monitor.info(f"Unguided pose {pose_num}: UFF minimized "
                                 f"(E {e_before:.1f}→{e_after:.1f}, {uff_msg})")
                else:
                    monitor.warning(f"Unguided pose {pose_num}: UFF minimization failed "
                                    f"({uff_msg}), using raw docked pose")

            shutil.copy2(out_sdf, final_sdf)
            centroid = _sdf_centroid(final_sdf)
            all_results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode="unguided", pocket_id=None,
                pocket_unique_id=None, pose_num=pose_num,
                pocket_center=None, pose_centroid=centroid,
                sdf_path=final_sdf, success=True,
                uff_minimized=_uff_minimized,
                uff_energy_before=_uff_e_before,
                uff_energy_after=_uff_e_after,
            ))
            unguided_count += 1
            monitor.info(f"Unguided pose {pose_num}/{n_unguided}: OK  "
                         f"centroid=({centroid[0]:.1f}, {centroid[1]:.1f}, {centroid[2]:.1f})"
                         if centroid else f"Unguided pose {pose_num}/{n_unguided}: OK")
        else:
            all_results.append(GuidedPoseResult(
                protein_name=protein_name, ligand_name=ligand_name,
                mode="unguided", pocket_id=None,
                pocket_unique_id=None, pose_num=pose_num,
                pocket_center=None, pose_centroid=None,
                sdf_path=None, success=False, error=err[:200],
            ))
            monitor.warning(f"Unguided pose {pose_num}/{n_unguided}: FAIL  error={err[:120]}")

        try:
            shutil.rmtree(pose_dir, ignore_errors=True)
        except Exception:
            pass

    # ── D. Check which unguided poses fall inside known pockets ──
    monitor.section("Unguided pose ↔ pocket proximity check")
    all_fp_pockets = fpocket_pockets
    all_pr_pockets = p2rank_pockets
    all_known_pockets = all_fp_pockets + all_pr_pockets

    for r in all_results:
        if r.mode != "unguided" or not r.success or r.pose_centroid is None:
            continue
        pose_label = f"unguided_{r.pose_num:03d} ({ligand_name})"
        in_any = False
        nearest_id = None
        nearest_dist = float("inf")

        for pocket in all_known_pockets:
            dist = _pocket_distance(r.pose_centroid, pocket.center)
            if dist < nearest_dist:
                nearest_dist = dist
                nearest_id = pocket.unique_id
            if dist <= POCKET_MATCH_THRESHOLD:
                in_any = True
                monitor.pose_inside_pocket(
                    pose_id=pose_label, pocket_id=pocket.unique_id,
                    distance=dist, threshold=POCKET_MATCH_THRESHOLD,
                    pocket_source=pocket.source,
                )

        if not in_any and nearest_id is not None:
            monitor.pose_outside_all_pockets(
                pose_id=pose_label, nearest_pocket=nearest_id,
                nearest_dist=nearest_dist, threshold=POCKET_MATCH_THRESHOLD,
            )

    # Summary for this combo
    ok_fp = sum(1 for r in all_results if r.mode == "fpocket" and r.success)
    ok_pr = sum(1 for r in all_results if r.mode == "p2rank" and r.success)
    ok_ug = sum(1 for r in all_results if r.mode == "unguided" and r.success)
    fail_total = sum(1 for r in all_results if not r.success)
    monitor.header(f"COMBO DONE  {combo}")
    print(f"  fpocket-guided: {ok_fp}  |  p2rank-guided: {ok_pr}  |  "
          f"unguided: {ok_ug}  |  failed: {fail_total}  |  "
          f"total OK: {ok_fp + ok_pr + ok_ug}")
    if fail_total > 0:
        monitor.warning(f"{combo}: {fail_total} pose(s) failed across all modes")

    # Save results
    with open(results_file, "w") as f:
        json.dump([r.to_dict() for r in all_results], f, indent=2)

    return all_results


# ── 4. Batch run across all protein-ligand pairs ─────────────────────────

def _match_protein_name_to_pockets(protein_pdb: Path) -> str:
    return protein_pdb.stem


def collect_input_files(directory: Path, extensions: List[str]) -> List[Path]:
    files = []
    for ext in extensions:
        files.extend(sorted(p for p in directory.glob(f"*{ext}") if p.is_file()))
    return sorted(set(files))


# Collect inputs — use only .pdb ligands to avoid duplicate docking of the same molecule
ligand_files = collect_input_files(drugs_dir, [".pdb"])
receptor_files = collect_input_files(receptors_dir, [".pdb"])

n_fp_total = N_TOP_POCKETS * POSES_PER_POCKET
n_pr_total = N_TOP_POCKETS * POSES_PER_POCKET
n_per_combo = n_fp_total + n_pr_total + N_UNGUIDED_POSES
total_combos = len(receptor_files) * len(ligand_files)

monitor.header("POCKET-GUIDED EQUIBIND DOCKING (Multi-Pose per Pocket)")
print(f"  Monitor level: {MONITOR_LEVEL.name}")
print(f"  EquiBind mode: IN-PROCESS (no subprocess overhead)")
print(f"  Proteins: {len(receptor_files)}")
for p in receptor_files:
    print(f"    - {p.name}")
print(f"  Ligands: {len(ligand_files)}")
for l in ligand_files:
    print(f"    - {l.name}")
print(f"\n  Pockets configuration:")
print(f"    Top fpocket pockets per protein: {N_TOP_POCKETS}")
print(f"    Top p2rank pockets per protein:  {N_TOP_POCKETS}")
print(f"    Poses per pocket:                {POSES_PER_POCKET}")
print(f"    fpocket poses per combo:         {n_fp_total}  ({N_TOP_POCKETS} × {POSES_PER_POCKET})")
print(f"    p2rank poses per combo:          {n_pr_total}  ({N_TOP_POCKETS} × {POSES_PER_POCKET})")
print(f"    Unguided poses per combo:        {N_UNGUIDED_POSES}")
print(f"    TOTAL poses per combo:           {n_per_combo}")
print(f"    Total combinations:              {total_combos}")
print(f"    Expected grand total poses:      {total_combos * n_per_combo}")
print(f"\n  Pocket enforcement:                {'ON' if FORCE_POCKET else 'OFF'}")
print(f"    Max trials per pose:             {MAX_POCKET_TRIALS}")
print(f"    Pocket match threshold:          {POCKET_MATCH_THRESHOLD} Å")
print(f"  Output directory:                  {POCKET_GUIDED_OUTPUT_DIR}")
print(f"  Logs directory:                    {LOG_DIR}")
print()

# Collect pockets for each protein
protein_fpocket: Dict[str, List[PocketInfo]] = {}
protein_p2rank: Dict[str, List[PocketInfo]] = {}

for prot in receptor_files:
    pname = _match_protein_name_to_pockets(prot)
    fp_pockets = parse_fpocket_results(fpocket_results_folder, pname, N_TOP_POCKETS)
    pr_pockets = parse_p2rank_results(p2rank_folder, pname, N_TOP_POCKETS)
    protein_fpocket[pname] = fp_pockets
    protein_p2rank[pname] = pr_pockets
    monitor.info(f"{pname}: fpocket={len(fp_pockets)} pockets {[p.unique_id for p in fp_pockets]}")
    monitor.info(f"{pname}: p2rank ={len(pr_pockets)} pockets {[p.unique_id for p in pr_pockets]}")

print()

# Run docking
all_guided_results: List[GuidedPoseResult] = []
start_time = time.time()

for idx, (protein, ligand) in enumerate(itertools.product(receptor_files, ligand_files), 1):
    pname = _match_protein_name_to_pockets(protein)
    monitor.header(f"[{idx}/{total_combos}] {protein.name} × {ligand.name}")

    fp_pockets = protein_fpocket.get(pname, [])
    pr_pockets = protein_p2rank.get(pname, [])

    if not fp_pockets and not pr_pockets:
        monitor.warning(f"No pockets found for protein {pname} — only unguided poses will be generated")

    results = dock_guided_by_pockets(
        protein_pdb=protein,
        ligand_file=ligand,
        fpocket_pockets=fp_pockets,
        p2rank_pockets=pr_pockets,
        poses_per_pocket=POSES_PER_POCKET,
        n_unguided=N_UNGUIDED_POSES,
        skip_existing=SKIP_EXISTING,
        force_pocket=FORCE_POCKET,
        max_pocket_trials=MAX_POCKET_TRIALS,
    )
    all_guided_results.extend(results)

elapsed = time.time() - start_time


# ── 5. Statistics: how many unguided poses fall in fpocket / p2rank pockets ──

def _distance(a, b):
    return np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


monitor.header("STATISTICS: UNGUIDED POSES vs FPOCKET / P2RANK POCKETS")

# Build lookup of all pockets for each protein
all_pockets_by_protein: Dict[str, Dict[str, List[PocketInfo]]] = defaultdict(lambda: {"fpocket": [], "p2rank": []})
for pname, pockets in protein_fpocket.items():
    all_pockets_by_protein[pname]["fpocket"] = pockets
for pname, pockets in protein_p2rank.items():
    all_pockets_by_protein[pname]["p2rank"] = pockets

# Analyze unguided poses
unguided_results = [r for r in all_guided_results if r.mode == "unguided" and r.success and r.pose_centroid]
fpocket_guided = [r for r in all_guided_results if r.mode == "fpocket" and r.success]
p2rank_guided = [r for r in all_guided_results if r.mode == "p2rank" and r.success]

stats_rows = []
total_unguided = 0
total_in_fpocket = 0
total_in_p2rank = 0
total_in_either = 0
total_in_neither = 0

# Group unguided by protein
unguided_by_protein: Dict[str, List[GuidedPoseResult]] = defaultdict(list)
for r in unguided_results:
    unguided_by_protein[r.protein_name].append(r)

for pname, poses in sorted(unguided_by_protein.items()):
    fp_pockets = all_pockets_by_protein[pname]["fpocket"]
    pr_pockets = all_pockets_by_protein[pname]["p2rank"]

    for pose in poses:
        total_unguided += 1
        in_fp = False
        in_pr = False
        nearest_fp_dist = float("inf")
        nearest_pr_dist = float("inf")
        nearest_fp_id = None
        nearest_pr_id = None

        for pocket in fp_pockets:
            dist = _distance(pose.pose_centroid, pocket.center)
            if dist < nearest_fp_dist:
                nearest_fp_dist = dist
                nearest_fp_id = pocket.unique_id
            if dist <= POCKET_MATCH_THRESHOLD:
                in_fp = True

        for pocket in pr_pockets:
            dist = _distance(pose.pose_centroid, pocket.center)
            if dist < nearest_pr_dist:
                nearest_pr_dist = dist
                nearest_pr_id = pocket.unique_id
            if dist <= POCKET_MATCH_THRESHOLD:
                in_pr = True

        if in_fp:
            total_in_fpocket += 1
        if in_pr:
            total_in_p2rank += 1
        if in_fp or in_pr:
            total_in_either += 1
        if not in_fp and not in_pr:
            total_in_neither += 1

        stats_rows.append({
            "protein": pname,
            "ligand": pose.ligand_name,
            "pose_sdf": str(pose.sdf_path.name) if pose.sdf_path else "",
            "pose_centroid_x": round(pose.pose_centroid[0], 2),
            "pose_centroid_y": round(pose.pose_centroid[1], 2),
            "pose_centroid_z": round(pose.pose_centroid[2], 2),
            "in_fpocket": in_fp,
            "nearest_fpocket_id": nearest_fp_id,
            "nearest_fpocket_dist_A": round(nearest_fp_dist, 2) if nearest_fp_dist < float("inf") else None,
            "in_p2rank": in_pr,
            "nearest_p2rank_id": nearest_pr_id,
            "nearest_p2rank_dist_A": round(nearest_pr_dist, 2) if nearest_pr_dist < float("inf") else None,
            "in_any_pocket": in_fp or in_pr,
        })

import pandas as pd

stats_df = pd.DataFrame(stats_rows)
pd.set_option("display.max_colwidth", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 0)

print(f"\n  Total unguided poses analyzed: {total_unguided}")
print(f"  Poses within an fpocket pocket (<{POCKET_MATCH_THRESHOLD}Å): {total_in_fpocket} ({100*total_in_fpocket/max(total_unguided,1):.1f}%)")
print(f"  Poses within a p2rank pocket  (<{POCKET_MATCH_THRESHOLD}Å): {total_in_p2rank} ({100*total_in_p2rank/max(total_unguided,1):.1f}%)")
print(f"  Poses in EITHER fpocket or p2rank:  {total_in_either} ({100*total_in_either/max(total_unguided,1):.1f}%)")
print(f"  Poses in NEITHER pocket:            {total_in_neither} ({100*total_in_neither/max(total_unguided,1):.1f}%)")

print(f"\n  Total guided docking runs:")
print(f"    fpocket-guided: {len(fpocket_guided)} successful poses")
print(f"    p2rank-guided:  {len(p2rank_guided)} successful poses")
print(f"    unguided:       {len(unguided_results)} successful poses")
print(f"\n  Total elapsed time: {elapsed:.1f}s ({elapsed/60:.1f} min)")

# Per-protein summary
print("\n" + "-" * 80)
print("PER-PROTEIN SUMMARY")
print("-" * 80)
if not stats_df.empty:
    summary = stats_df.groupby("protein").agg(
        total_unguided=("in_any_pocket", "count"),
        in_fpocket=("in_fpocket", "sum"),
        in_p2rank=("in_p2rank", "sum"),
        in_any=("in_any_pocket", "sum"),
    ).reset_index()
    summary["pct_in_fpocket"] = (100 * summary["in_fpocket"] / summary["total_unguided"]).round(1)
    summary["pct_in_p2rank"] = (100 * summary["in_p2rank"] / summary["total_unguided"]).round(1)
    summary["pct_in_any"] = (100 * summary["in_any"] / summary["total_unguided"]).round(1)
    display(summary)

print("\nDetailed per-pose results:")
if not stats_df.empty:
    display(stats_df)

# Save stats
stats_csv = POCKET_GUIDED_OUTPUT_DIR / "unguided_vs_pockets_statistics.csv"
stats_df.to_csv(stats_csv, index=False)
print(f"\nStatistics saved to: {stats_csv}")

# Save full summary JSON
full_summary = {
    "timestamp": datetime.now().isoformat(),
    "config": {
        "n_top_pockets": N_TOP_POCKETS,
        "poses_per_pocket": POSES_PER_POCKET,
        "n_unguided_poses": N_UNGUIDED_POSES,
        "pocket_match_threshold_A": POCKET_MATCH_THRESHOLD,
        "force_pocket": FORCE_POCKET,
        "max_pocket_trials": MAX_POCKET_TRIALS,
        "monitor_level": MONITOR_LEVEL.name,
        "uff_minimize": UFF_MINIMIZE,
        "uff_max_iters": UFF_MAX_ITERS,
        "uff_protein_constraint_wt": UFF_PROTEIN_CONSTRAINT_WT,
        "uff_proximity_radius_A": UFF_PROXIMITY_RADIUS,
    },
    "totals": {
        "unguided_poses": total_unguided,
        "in_fpocket": total_in_fpocket,
        "in_p2rank": total_in_p2rank,
        "in_either": total_in_either,
        "in_neither": total_in_neither,
        "fpocket_guided_success": len(fpocket_guided),
        "p2rank_guided_success": len(p2rank_guided),
        "uff_minimized_count": sum(1 for r in all_guided_results if r.uff_minimized),
        "elapsed_s": round(elapsed, 1),
    },
    "monitor_counters": {
        "equibind_calls": monitor.call_count,
        "equibind_success": monitor.success_count,
        "equibind_fail": monitor.fail_count,
        "poses_in_pocket": monitor.in_pocket_count,
        "poses_rejected": monitor.rejected_count,
        "poses_outside_all": monitor.outside_pocket_count,
    },
    "all_results": [r.to_dict() for r in all_guided_results],
}
summary_json = POCKET_GUIDED_OUTPUT_DIR / "pocket_guided_summary.json"
with open(summary_json, "w") as f:
    json.dump(full_summary, f, indent=2)
print(f"\nFull summary saved to: {summary_json}")

# ── Print signal monitor summary ──
monitor.print_summary()
