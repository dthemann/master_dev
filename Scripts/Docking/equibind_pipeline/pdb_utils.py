"""Fast PDB parsing helpers shared by cropping and other modules."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Set, Tuple

import numpy as np


def parse_pdb_atoms(pdb_path: Path) -> List[dict]:
    """Parse ATOM/HETATM records from a PDB file.

    Returns list of dicts with keys: line, record_type, atom_num, atom_name,
    res_name, chain, res_num, x, y, z.
    """
    atoms: List[dict] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                atoms.append({
                    "line": line,
                    "record_type": line[0:6].strip(),
                    "atom_num": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "chain": line[21:22].strip() or "A",
                    "res_num": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                })
            except (ValueError, IndexError):
                continue
    return atoms


def atom_coords(atoms: Iterable[dict]) -> np.ndarray:
    """Return an (N,3) numpy array of atom coordinates."""
    return np.fromiter(
        (c for a in atoms for c in (a["x"], a["y"], a["z"])),
        dtype=np.float64,
    ).reshape(-1, 3)


def residues_near_center(
    atoms: List[dict],
    center: Tuple[float, float, float],
    radius: float,
) -> Set[Tuple[str, int]]:
    """Return {(chain, res_num)} where any atom is within ``radius`` of ``center``.

    Vectorized via numpy distance computation.
    """
    if not atoms:
        return set()
    coords = atom_coords(atoms)
    diffs = coords - np.asarray(center, dtype=np.float64)
    dists2 = np.einsum("ij,ij->i", diffs, diffs)
    mask = dists2 <= radius * radius
    return {(atoms[i]["chain"], atoms[i]["res_num"]) for i in np.nonzero(mask)[0]}
