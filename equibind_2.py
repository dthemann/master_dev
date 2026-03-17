#!/usr/bin/env python3
"""
EquiBind batch docking script.

Docks every ligand SDF in Drugs/ against every receptor PDB in Orai/
using EquiBind's multiligand_inference.py, then collects the results
into per-combo directories with individual pose SDF files plus a
summary CSV and JSON.
"""

import json
import csv
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# ── Configuration ──────────────────────────────────────────────────────────

# Paths
workspace_root = Path.cwd()
drugs_dir = workspace_root / "Drugs"
receptors_dir = workspace_root / "Orai"
output_root = workspace_root / "equibind_docked_poses_2"

# EquiBind paths
EQUIBIND_DIR = Path("/home/manndo/docking_tools/EquiBind")
EQUIBIND_MULTILIGAND_SCRIPT = EQUIBIND_DIR / "multiligand_inference.py"

# Conda environment name for EquiBind
CONDA_ENV = "equibind"

# Whether to skip combos whose output directory already exists
SKIP_EXISTING = True
# Device for EquiBind (cuda or cpu)
DEVICE = "cpu"


# ── Helpers ────────────────────────────────────────────────────────────────

def gather_ligands(drugs_dir: Path) -> list[Path]:
    """Return sorted list of ligand .sdf files in the drugs directory."""
    return sorted(drugs_dir.glob("*.sdf"))


def gather_receptors(receptors_dir: Path) -> list[Path]:
    """Return sorted list of receptor .pdb files in the receptors directory."""
    return sorted(receptors_dir.glob("*.pdb"))


def combo_name(ligand_path: Path, receptor_path: Path) -> str:
    """Build a combo directory name: <ligand_stem>__<receptor_stem>."""
    return f"{ligand_path.stem}__{receptor_path.stem}"


def run_equibind(
    ligand_sdf: Path,
    receptor_pdb: Path,
    output_dir: Path,
    device: str = "cuda",
) -> dict:
    """
    Run EquiBind multiligand_inference on one (ligand, receptor) pair.

    Returns a result dict with keys:
        protein, ligand, combo, output_dir, success, output_sdf, error, elapsed_s
    """
    combo = combo_name(ligand_sdf, receptor_pdb)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "conda", "run", "-n", CONDA_ENV, "--no-capture-output",
        "python", str(EQUIBIND_MULTILIGAND_SCRIPT),
        "-l", str(ligand_sdf),
        "-r", str(receptor_pdb),
        "-o", str(output_dir),
        "--device", device,
    ]

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(EQUIBIND_DIR),
            timeout=600,  # generous 10-min timeout per combo
        )
        elapsed = round(time.time() - t0, 1)

        output_sdf = output_dir / "output.sdf"
        if proc.returncode != 0:
            return dict(
                protein=receptor_pdb.name,
                ligand=ligand_sdf.name,
                combo=combo,
                output_dir=str(output_dir),
                success=False,
                output_sdf=None,
                error=proc.stderr.strip()[-500:] if proc.stderr else f"exit code {proc.returncode}",
                elapsed_s=elapsed,
            )
        if not output_sdf.exists() or output_sdf.stat().st_size == 0:
            return dict(
                protein=receptor_pdb.name,
                ligand=ligand_sdf.name,
                combo=combo,
                output_dir=str(output_dir),
                success=False,
                output_sdf=None,
                error="No output.sdf produced (missing or empty)",
                elapsed_s=elapsed,
            )
        return dict(
            protein=receptor_pdb.name,
            ligand=ligand_sdf.name,
            combo=combo,
            output_dir=str(output_dir),
            success=True,
            output_sdf=str(output_sdf),
            error="",
            elapsed_s=elapsed,
        )
    except subprocess.TimeoutExpired:
        elapsed = round(time.time() - t0, 1)
        return dict(
            protein=receptor_pdb.name,
            ligand=ligand_sdf.name,
            combo=combo,
            output_dir=str(output_dir),
            success=False,
            output_sdf=None,
            error="Timed out after 600 s",
            elapsed_s=elapsed,
        )
    except Exception as exc:
        elapsed = round(time.time() - t0, 1)
        return dict(
            protein=receptor_pdb.name,
            ligand=ligand_sdf.name,
            combo=combo,
            output_dir=str(output_dir),
            success=False,
            output_sdf=None,
            error=str(exc),
            elapsed_s=elapsed,
        )


def sdf_centroid(sdf_path: Path) -> tuple[float, float, float] | None:
    """Compute the centroid of an SDF file by parsing ATOM coords."""
    try:
        from rdkit import Chem
        mol = Chem.SDMolSupplier(str(sdf_path), removeHs=False)[0]
        if mol is None:
            return None
        pos = mol.GetConformer().GetPositions()
        centroid = pos.mean(axis=0)
        return tuple(centroid.tolist())
    except Exception:
        return None


def split_output_sdf(output_sdf: Path, combo_dir: Path) -> list[dict]:
    """
    Split an output.sdf (may contain multiple molecules) into
    individual pose_XX.sdf files. Returns a list of pose metadata dicts.
    """
    try:
        from rdkit import Chem
    except ImportError:
        # Fallback: just copy the whole SDF as pose_01.sdf
        import shutil
        pose_path = combo_dir / "pose_01.sdf"
        shutil.copy2(output_sdf, pose_path)
        return [dict(pose_id=1, sdf_path=str(pose_path), centroid=None, conformer_id=1)]

    try:
        supplier = Chem.SDMolSupplier(str(output_sdf), removeHs=False)
    except OSError:
        return []
    poses = []
    for idx, mol in enumerate(supplier, start=1):
        if mol is None:
            continue
        pose_path = combo_dir / f"pose_{idx:02d}.sdf"
        writer = Chem.SDWriter(str(pose_path))
        writer.write(mol)
        writer.close()
        centroid = mol.GetConformer().GetPositions().mean(axis=0).tolist()
        poses.append(dict(
            pose_id=idx,
            sdf_path=str(pose_path),
            centroid=centroid,
            conformer_id=idx,
        ))
    return poses


def build_metadata(
    ligand_path: Path,
    receptor_path: Path,
    combo_dir: Path,
    poses: list[dict],
    status: str = "complete",
) -> dict:
    """Build per-combo docking_metadata.json content."""
    return dict(
        protein_name=receptor_path.stem,
        ligand_name=ligand_path.stem,
        protein_path=str(receptor_path),
        ligand_path=str(ligand_path),
        output_dir=str(combo_dir),
        status=status,
        poses=poses,
    )


# ── Main ───────────────────────────────────────────────────────────────────

def main() -> None:
    ligands = gather_ligands(drugs_dir)
    receptors = gather_receptors(receptors_dir)

    if not ligands:
        print(f"No .sdf files found in {drugs_dir}")
        sys.exit(1)
    if not receptors:
        print(f"No .pdb files found in {receptors_dir}")
        sys.exit(1)

    print(f"Ligands  ({len(ligands)}): {[l.name for l in ligands]}")
    print(f"Receptors ({len(receptors)}): {[r.name for r in receptors]}")
    total = len(ligands) * len(receptors)
    print(f"Total combos: {total}\n")

    output_root.mkdir(parents=True, exist_ok=True)
    summary: list[dict] = []
    csv_rows: list[dict] = []
    done = 0

    for receptor in receptors:
        for ligand in ligands:
            done += 1
            cname = combo_name(ligand, receptor)
            combo_dir = output_root / cname
            tag = f"[{done}/{total}]"

            # ── Skip if already complete ──────────────────────────────
            if SKIP_EXISTING and combo_dir.exists() and (combo_dir / "output.sdf").exists():
                print(f"{tag} SKIP (exists): {cname}")
                continue

            print(f"{tag} Running: {cname} ...")
            result = run_equibind(ligand, receptor, combo_dir, device=DEVICE)
            summary.append(result)

            if result["success"]:
                print(f"     ✓ done in {result['elapsed_s']}s")
                # Split output.sdf into individual pose files
                poses = split_output_sdf(Path(result["output_sdf"]), combo_dir)
                meta = build_metadata(ligand, receptor, combo_dir, poses)
                with open(combo_dir / "docking_metadata.json", "w") as f:
                    json.dump(meta, f, indent=2)

                # Accumulate CSV rows
                for p in poses:
                    cx, cy, cz = (p["centroid"] if p["centroid"] else [None, None, None])
                    csv_rows.append(dict(
                        protein=receptor.stem,
                        ligand=ligand.stem,
                        pose_id=p["pose_id"],
                        conformer_id=p["conformer_id"],
                        sdf_path=p["sdf_path"],
                        centroid_x=cx,
                        centroid_y=cy,
                        centroid_z=cz,
                        status="complete",
                        num_conformers_generated=len(poses),
                    ))
            else:
                print(f"     ✗ FAILED: {result['error'][:120]}")

    # ── Write summary JSON ─────────────────────────────────────────────
    summary_path = output_root / "docking_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSummary written to {summary_path}")

    # ── Write / append poses CSV ───────────────────────────────────────
    csv_path = output_root / "equibind_poses.csv"
    fieldnames = [
        "protein", "ligand", "pose_id", "conformer_id", "sdf_path",
        "centroid_x", "centroid_y", "centroid_z", "status",
        "num_conformers_generated",
    ]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Poses CSV written to {csv_path}")

    # ── Print stats ────────────────────────────────────────────────────
    ok = sum(1 for r in summary if r["success"])
    fail = sum(1 for r in summary if not r["success"])
    print(f"\nResults: {ok} succeeded, {fail} failed out of {len(summary)} run")


if __name__ == "__main__":
    main()
