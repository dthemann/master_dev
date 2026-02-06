#!/usr/bin/env python3
"""
Manual GPF creation - Python 3

Adds a small CLI so you can create a .gpf for AutoGrid4 easily, before using
AutoDock Vina with '--scoring ad4'.
"""

from __future__ import annotations

import argparse
import math
import os
from typing import Iterable, Sequence, Tuple


def _as_tuple3(values: Iterable[float]) -> Tuple[float, float, float]:
    vals = list(values)
    if len(vals) != 3:
        raise ValueError("Expected exactly 3 values (x, y, z)")
    return float(vals[0]), float(vals[1]), float(vals[2])


def create_gpf_from_pdb(
    receptor_pdbqt: str,
    output_gpf: str,
    grid_center: Sequence[float],
    grid_points: Tuple[int, int, int] = (100, 100, 100),
    spacing: float = 0.375,
    # Canonical AutoDock4 ligand atom types (subset, ordered to match typical parameter libraries)
    # Excludes raw 'H' (AD4 uses 'HD' for polar hydrogens) and exotic types unless explicitly needed.
    ligand_types: Sequence[str] = (
        "HD", "C", "A", "N", "NA", "OA", "F", "P", "SA", "S", "Cl", "Br", "I"
    ),
    parameter_file: str | None = None,
) -> str:
    """Create a GPF file for AutoGrid4.

    Parameters
    ----------
    receptor_pdbqt : str
        Path to receptor PDBQT file.
    output_gpf : str
        Path to write .gpf file.
    grid_center : Sequence[float]
        Center of the grid in Angstroms (x, y, z).
    grid_points : Tuple[int, int, int]
        Number of grid points along x, y, z (npts). Default (100,100,100).
    spacing : float
        Grid spacing in Angstroms. Default 0.375.
    ligand_types : Sequence[str]
        AD4 atom types expected for ligand. Used to emit map lines.

    Returns
    -------
    str
        The path to the created .gpf file.
    """

    # Normalize inputs
    center = _as_tuple3(grid_center)
    npts = tuple(int(v) for v in grid_points)
    if len(npts) != 3:
        raise ValueError("grid_points must have length 3")

    # Use basename (without directory) for map/grid file stems to avoid duplication and ensure autogrid output paths are clean.
    receptor_root = os.path.splitext(os.path.basename(receptor_pdbqt))[0]

    # Sanitize ligand_types: remove duplicates, strip whitespace, drop unsupported raw 'H'
    lt_clean = []
    for t in ligand_types:
        t = str(t).strip()
        if not t or t == "H":  # AD4 uses HD for polar H; skip bare H to avoid unknown atom type errors
            continue
        if t not in lt_clean:
            lt_clean.append(t)
    ligand_types = tuple(lt_clean)

    # One map line per ligand type (must match count exactly to avoid "Too many map keywords" error)
    map_lines = [f"map {receptor_root}.{t}.map" for t in ligand_types]

    # Build GPF body cleanly using a list, then join
    gpf_lines = [
        "# Grid Parameter File",
        f"receptor {os.path.basename(receptor_pdbqt)}",  # basename so autogrid finds file when run in GPF dir
        f"gridfld {receptor_root}.maps.fld",
    ]
    if parameter_file:
        gpf_lines.append(f"parameter_file {parameter_file}")
    gpf_lines.extend([
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        f"spacing {spacing}",
        f"gridcenter {center[0]} {center[1]} {center[2]}",
        f"ligand_types {' '.join(ligand_types)}",
        "smooth 0.5",
        "",
    ])
    gpf_lines.extend(map_lines)
    gpf_lines.extend([
        f"elecmap {receptor_root}.e.map",
        "dielectric -0.1465",
        f"dsolvmap {receptor_root}.d.map",
        "",
    ])
    gpf_content = "\n".join(gpf_lines)

    os.makedirs(os.path.dirname(output_gpf) or ".", exist_ok=True)
    with open(output_gpf, "w") as f:
        f.write(gpf_content)

    print(f"GPF file created: {output_gpf}")
    return output_gpf


def _compute_npts_from_size(size_xyz: Sequence[float], spacing: float) -> Tuple[int, int, int]:
    """Compute integer npts from box size in Angstroms and spacing.

    The AD4 convention is: physical_length ≈ (npts - 1) * spacing.
    We round to the nearest integer to cover the requested size.
    """
    sx, sy, sz = _as_tuple3(size_xyz)
    # Use npts = round(size/spacing) + 1 to better span desired size
    nx = max(2, int(round(sx / spacing)) + 1)
    ny = max(2, int(round(sy / spacing)) + 1)
    nz = max(2, int(round(sz / spacing)) + 1)
    return nx, ny, nz


def main():
    parser = argparse.ArgumentParser(description="Create an AutoGrid4 .gpf file")
    parser.add_argument("--receptor", required=True, help="Path to receptor .pdbqt")
    parser.add_argument("--out", required=False, help="Output .gpf path (default: <receptor>.gpf)")
    parser.add_argument(
        "--center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        required=True,
        help="Grid center (Å): x y z",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--npts",
        nargs=3,
        type=int,
        metavar=("NX", "NY", "NZ"),
        help="Grid points (npts) along x y z",
    )
    group.add_argument(
        "--size",
        nargs=3,
        type=float,
        metavar=("SX", "SY", "SZ"),
        help="Physical grid size (Å) along x y z; converted to npts using spacing",
    )
    parser.add_argument("--spacing", type=float, default=0.375, help="Grid spacing (Å)")
    parser.add_argument(
        "--ligand-types",
        nargs="+",
        default=[
            # Canonical AutoDock4 ordering expected by many tools (matches internal parameter tables)
            # 12 heavy atom types + polar hydrogens last:
            "C", "A", "N", "NA", "OA", "S", "SA", "P", "F", "Cl", "Br", "I", "HD"
        ],
        help=(
            "AD4 ligand atom types to include maps for. Defaults to canonical ordering: "
            "C A N NA OA S SA P F Cl Br I HD"
        ),
    )
    parser.add_argument(
        "--parameter-file",
        type=str,
        default=None,
        help="Path to AD4 parameters file (e.g., AD4_parameters.dat) to include as 'parameter_file' directive.",
    )

    args = parser.parse_args()

    receptor = args.receptor
    if not os.path.isfile(receptor):
        raise SystemExit(f"ERROR: receptor file not found: {receptor}")

    out_gpf = args.out or receptor.replace(".pdbqt", ".gpf")

    if args.npts is not None:
        npts = tuple(int(v) for v in args.npts)
    elif args.size is not None:
        npts = _compute_npts_from_size(args.size, args.spacing)
    else:
        npts = (100, 100, 100)

    create_gpf_from_pdb(
        receptor_pdbqt=receptor,
        output_gpf=out_gpf,
        grid_center=args.center,
        grid_points=npts,
        spacing=args.spacing,
        ligand_types=args.ligand_types,
        parameter_file=args.parameter_file,
    )


if __name__ == "__main__":
    main()
