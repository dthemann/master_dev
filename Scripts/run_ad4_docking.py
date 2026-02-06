#!/usr/bin/env python3
"""
End-to-end AD4-scoring docking with AutoDock Vina, generating maps via AutoGrid4.

Steps:
 1) Create a .gpf (grid parameter file) using Scripts/create_gpf.py
 2) Run autogrid4 to generate AD4 maps (.map, .fld, .xyz)
 3) Run vina with --scoring ad4

Requirements:
 - autogrid4 available on PATH
 - vina (>=1.2) available at VINA_BIN (default: ~/AutoDock-Vina/build/linux/release/vina)
 - receptor and ligand in PDBQT format

Example:
    ./run_ad4_docking.py \
        --receptor prepared/Orai1WT-MDSnap-Fr300.pdbqt \
        --ligand prepared/ligands/2abp-nh2-OPT.pdbqt \
        --center 0 10 -6 \
        --size 30 30 30 \
        --spacing 0.375 \
        --out prepared/docking/Orai1WT-MDSnap-Fr300__2abp-nh2-OPT_vina_out.pdbqt
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from typing import Sequence, Tuple

# Path to the Vina binary to use (expand ~)
VINA_BIN = os.path.expanduser("~/AutoDock-Vina/build/linux/release/vina")

# Import the GPF creator from sibling module
sys.path.append(os.path.dirname(__file__))
from create_gpf import create_gpf_from_pdb  # noqa: E402

def run_cmd(cmd: Sequence[str], cwd: str | None = None) -> int:
        print("$", " ".join(shlex.quote(c) for c in cmd))
        try:
                proc = subprocess.run(cmd, cwd=cwd, check=False)
                return proc.returncode
        except FileNotFoundError:
                print(f"ERROR: Command not found: {cmd[0]}")
                return 127


def compute_npts_from_size(size_xyz: Sequence[float], spacing: float) -> Tuple[int, int, int]:
        sx, sy, sz = (float(size_xyz[0]), float(size_xyz[1]), float(size_xyz[2]))
        # npts ~ round(size/spacing)+1 to span requested size per AD4 convention
        import math

        nx = max(2, int(round(sx / spacing)) + 1)
        ny = max(2, int(round(sy / spacing)) + 1)
        nz = max(2, int(round(sz / spacing)) + 1)
        return nx, ny, nz


def run_autogrid(gpf_path: str) -> int:
        if not os.path.isfile(gpf_path):
                print(f"ERROR: GPF not found: {gpf_path}")
                return 2
        # Log file alongside GPF
        base = os.path.splitext(gpf_path)[0]
        glg = base + ".glg"
        cmd = ["autogrid4", "-p", os.path.basename(gpf_path), "-l", os.path.basename(glg)]
        # Run in the GPF directory so output maps are written next to receptor/GPF
        return run_cmd(cmd, cwd=os.path.dirname(gpf_path) or None)


def run_vina_ad4(
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        out_pdbqt: str,
        exhaustiveness: int = 16,
        num_modes: int = 9,
        seed: int | None = None,
) -> int:
        """Run AutoDock Vina in AD4 scoring mode using precomputed maps.

        Vina 1.2 requires precomputed AutoGrid maps when using --scoring ad4; a receptor
        model cannot be supplied directly (error: "No receptor allowed"). We therefore
        point Vina at the maps prefix produced by autogrid4 (basename of receptor).
        """
        # Ensure the configured Vina binary exists and is executable
        if not (os.path.isfile(VINA_BIN) and os.access(VINA_BIN, os.X_OK)):
                print(f"ERROR: vina not found or not executable at: {VINA_BIN}")
                return 127

        maps_prefix = os.path.splitext(os.path.basename(receptor_pdbqt))[0]
        maps_dir = os.path.dirname(receptor_pdbqt)
        prefix_path = os.path.join(maps_dir, maps_prefix)  # pass full path prefix so vina finds maps

        cmd = [
                VINA_BIN,
                "--ligand",
                ligand_pdbqt,
                "--maps",
                prefix_path,
                "--scoring",
                "ad4",
                "--exhaustiveness",
                str(exhaustiveness),
                "--num_modes",
                str(num_modes),
                "--out",
                out_pdbqt,
        ]
        if seed is not None:
                cmd += ["--seed", str(seed)]
        return run_cmd(cmd)


def main():
        ap = argparse.ArgumentParser(description="Create AD4 maps with AutoGrid4 and dock with Vina --scoring ad4")
        ap.add_argument("--receptor", required=True, help="Path to receptor .pdbqt")
        ap.add_argument("--ligand", required=True, help="Path to ligand .pdbqt")
        ap.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"), required=True, help="Grid center (Å)")
        size_grp = ap.add_mutually_exclusive_group(required=True)
        size_grp.add_argument("--size", nargs=3, type=float, metavar=("SX", "SY", "SZ"), help="Box size (Å)")
        size_grp.add_argument("--npts", nargs=3, type=int, metavar=("NX", "NY", "NZ"), help="Grid points (npts)")
        ap.add_argument("--spacing", type=float, default=0.375, help="Grid spacing (Å), used when --size is given")
        ap.add_argument(
            "--ligand-types",
            nargs="+",
            default=["C", "A", "N", "NA", "OA", "S", "SA", "P", "F", "Cl", "Br", "I", "B"],
            help="AD4 ligand atom types to generate maps for (includes B; exclude HD unless polar H map needed)",
        )
        ap.add_argument(
            "--parameter-file",
            type=str,
            default=os.path.expanduser("~/AutoDock-Vina/data/AD4_parameters.dat"),
            help="Path to AD4 parameters file; included in generated GPF if exists.",
        )
        ap.add_argument("--maps-only", action="store_true", help="Stop after generating maps (skip docking)")
        ap.add_argument("--out", required=True, help="Output docked pose .pdbqt path for Vina")
        ap.add_argument("--exhaustiveness", type=int, default=16)
        ap.add_argument("--num-modes", type=int, default=9)
        ap.add_argument("--seed", type=int)

        args = ap.parse_args()

        receptor = args.receptor
        ligand = args.ligand
        if not os.path.isfile(receptor):
                print(f"ERROR: receptor not found: {receptor}")
                return 2
        if not os.path.isfile(ligand):
                print(f"ERROR: ligand not found: {ligand}")
                return 2

        # Paths and grid setup
        gpf_path = receptor.replace(".pdbqt", ".gpf")

        if args.npts is not None:
                npts = tuple(int(v) for v in args.npts)
                size_for_vina = tuple(float(v) for v in args.size) if args.size else (
                        (npts[0] - 1) * args.spacing,
                        (npts[1] - 1) * args.spacing,
                        (npts[2] - 1) * args.spacing,
                )
        else:
                # Convert size to npts for GPF
                npts = compute_npts_from_size(args.size, args.spacing)
                size_for_vina = tuple(float(v) for v in args.size)

        # 1) Create GPF
        # Ensure parameter file presence (optional)
        parameter_file = args.parameter_file if (args.parameter_file and os.path.isfile(args.parameter_file)) else None
        create_gpf_from_pdb(
            receptor_pdbqt=receptor,
            output_gpf=gpf_path,
            grid_center=args.center,
            grid_points=npts,
            spacing=args.spacing,
            ligand_types=args.ligand_types,
            parameter_file=parameter_file,
        )

        # 2) Run AutoGrid4
        rc = run_autogrid(gpf_path)
        if rc != 0:
                print(f"ERROR: autogrid4 failed with exit code {rc}")
                return rc

        if args.maps_only:
                print("Maps generated; skipping docking as requested (--maps-only)")
                return 0

        # 3) Dock with Vina using AD4 scoring
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        rc = run_vina_ad4(
                receptor_pdbqt=receptor,
                ligand_pdbqt=ligand,
                out_pdbqt=args.out,
                exhaustiveness=args.exhaustiveness,
                num_modes=args.num_modes,
                seed=args.seed,
        )
        if rc != 0:
                print(f"ERROR: vina failed with exit code {rc}")
                return rc

        print(f"Docking completed: {args.out}")
        return 0


if __name__ == "__main__":
        raise SystemExit(main())
