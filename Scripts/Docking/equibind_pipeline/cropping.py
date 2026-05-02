"""Protein cropping for pocket-constrained EquiBind docking."""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from rdkit import Chem

from .config import CFG
from .monitor import monitor
from .pdb_utils import parse_pdb_atoms, residues_near_center
from .pockets import PocketInfo


# Cache: (protein_path, pocket_unique_id) -> (cropped_pdb_path, offset_vector)
_cropped_protein_cache: Dict[Tuple[str, str], Tuple[Path, np.ndarray]] = {}


def crop_protein_to_pocket(
    protein_pdb: Path,
    pocket_center: Tuple[float, float, float],
    output_pdb: Path,
    crop_radius: Optional[float] = None,
    buffer: Optional[float] = None,
) -> Tuple[bool, np.ndarray, int, int]:
    """Crop a PDB to residues near a pocket center, translated so center=origin.

    Returns ``(success, offset_vector, n_residues_kept, n_atoms_kept)``.
    To restore original coords from cropped: ``coord + offset_vector``.
    """
    crop_radius = CFG.pocket_crop_radius if crop_radius is None else crop_radius
    buffer = CFG.pocket_crop_buffer if buffer is None else buffer
    effective_radius = crop_radius + buffer

    atoms = parse_pdb_atoms(protein_pdb)
    if not atoms:
        return False, np.zeros(3), 0, 0

    near_residues = residues_near_center(atoms, pocket_center, effective_radius)
    if not near_residues:
        monitor.warning(f"No residues found within {effective_radius}\u00c5 of pocket center")
        return False, np.zeros(3), 0, 0

    if CFG.include_full_residues:
        kept_atoms = [a for a in atoms if (a["chain"], a["res_num"]) in near_residues]
    else:
        center_np = np.asarray(pocket_center, dtype=np.float64)
        kept_atoms = [
            a for a in atoms
            if np.linalg.norm(np.array([a["x"], a["y"], a["z"]]) - center_np) <= effective_radius
        ]

    if not kept_atoms:
        return False, np.zeros(3), 0, 0

    # ── Deduplicate atom names within each residue (EquiBind fails on dupes) ─
    deduped: list = []
    seen_names: Dict[Tuple[str, int], set] = {}
    n_dupes = 0
    for a in kept_atoms:
        key = (a["chain"], a["res_num"])
        s = seen_names.setdefault(key, set())
        if a["atom_name"] in s:
            n_dupes += 1
            continue
        s.add(a["atom_name"])
        deduped.append(a)
    if n_dupes:
        monitor.info(f"Dropped {n_dupes} duplicate atom name(s) from cropped PDB")
    kept_atoms = deduped

    offset = np.asarray(pocket_center, dtype=np.float64)
    output_pdb.parent.mkdir(parents=True, exist_ok=True)

    with open(output_pdb, "w", encoding="utf-8") as f:
        f.write("REMARK   CROPPED PROTEIN - pocket center at origin\n")
        f.write(f"REMARK   Original pocket center: {pocket_center[0]:.3f} "
                f"{pocket_center[1]:.3f} {pocket_center[2]:.3f}\n")
        f.write(f"REMARK   Crop radius: {effective_radius:.1f} A\n")
        f.write(f"REMARK   Residues kept: {len(near_residues)}\n")
        f.write(f"REMARK   To restore original coords: add offset "
                f"({offset[0]:.3f}, {offset[1]:.3f}, {offset[2]:.3f})\n")
        for atom_num, a in enumerate(kept_atoms, start=1):
            line = a["line"]
            new_x = a["x"] - offset[0]
            new_y = a["y"] - offset[1]
            new_z = a["z"] - offset[2]
            f.write(
                f"{line[0:6]}{atom_num:5d}{line[11:30]}"
                f"{new_x:8.3f}{new_y:8.3f}{new_z:8.3f}{line[54:]}"
            )
        f.write("END\n")

    return True, offset, len(near_residues), len(kept_atoms)


def get_cropped_protein_for_pocket(
    protein_pdb: Path, pocket: PocketInfo, prep_dir: Path,
    crop_radius: Optional[float] = None,
) -> Tuple[Optional[Path], np.ndarray]:
    """Cached cropped-protein producer."""
    crop_radius = CFG.pocket_crop_radius if crop_radius is None else crop_radius
    cache_key = (str(protein_pdb), pocket.unique_id)
    if cache_key in _cropped_protein_cache:
        return _cropped_protein_cache[cache_key]

    cropped_pdb = prep_dir / f"{protein_pdb.stem}_crop_{pocket.unique_id}.pdb"
    success, offset, n_res, n_atoms = crop_protein_to_pocket(
        protein_pdb=protein_pdb, pocket_center=pocket.center,
        output_pdb=cropped_pdb, crop_radius=crop_radius,
    )
    if not success:
        monitor.warning(f"Failed to crop protein for pocket {pocket.unique_id}")
        return None, np.zeros(3)

    monitor.info(
        f"Cropped protein for {pocket.unique_id}: {n_res} residues, {n_atoms} atoms "
        f"(radius={crop_radius + CFG.pocket_crop_buffer:.1f}\u00c5)"
    )
    _cropped_protein_cache[cache_key] = (cropped_pdb, offset)
    return cropped_pdb, offset


def translate_pose_back_to_original_frame(
    sdf_path: Path, offset: np.ndarray, output_path: Path,
) -> bool:
    """Add ``offset`` to all coordinates in the docked SDF."""
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            return False
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        ox, oy, oz = float(offset[0]), float(offset[1]), float(offset[2])
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x + ox, p.y + oy, p.z + oz))
        w = Chem.SDWriter(str(output_path))
        w.write(mol); w.close()
        return output_path.exists() and output_path.stat().st_size > 0
    except Exception as e:
        monitor.warning(f"Failed to translate pose back to original frame: {e}")
        return False


def clamp_pose_to_pocket(
    sdf_path: Path, pocket_center: Tuple[float, float, float],
    max_distance: float, output_path: Path,
) -> Tuple[bool, float, float]:
    """Translate a pose so its centroid is within ``max_distance`` of ``pocket_center``.

    Returns ``(success, original_distance, new_distance)``.
    """
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            return False, 0.0, 0.0
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        positions = conf.GetPositions()
        centroid = positions.mean(axis=0)
        pocket_np = np.asarray(pocket_center, dtype=np.float64)
        displacement = centroid - pocket_np
        original_distance = float(np.linalg.norm(displacement))

        if original_distance <= max_distance:
            shutil.copy2(sdf_path, output_path)
            return True, original_distance, original_distance

        direction = displacement / original_distance
        target_centroid = pocket_np + direction * max_distance
        shift = target_centroid - centroid
        sx, sy, sz = float(shift[0]), float(shift[1]), float(shift[2])
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x + sx, p.y + sy, p.z + sz))
        w = Chem.SDWriter(str(output_path)); w.write(mol); w.close()
        return (output_path.exists() and output_path.stat().st_size > 0,
                original_distance, max_distance)
    except Exception as e:
        monitor.warning(f"Failed to clamp pose to pocket: {e}")
        return False, 0.0, 0.0


def clear_cropped_protein_cache() -> None:
    _cropped_protein_cache.clear()
