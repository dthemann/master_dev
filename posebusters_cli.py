#!/usr/bin/env python3
"""
PoseBusters CLI — Validate docked poses from AutoDock, DiffDock, or EquiBind.

Usage:
    python posebusters_cli.py \
        --receptors /path/to/proteins/ \
        --ligands   /path/to/ligands/ \
        --poses     /path/to/docked_poses/ \
        --tool      autodock|diffdock|equibind \
        --output    /path/to/output/ \
        [--workers N] [--overwrite] [--save-interval N]
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from multiprocessing import Pool, cpu_count
from pathlib import Path

import pandas as pd


# ============================================================================
# CLI ARGUMENT PARSING
# ============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate docked poses with PoseBusters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--receptors", required=True, type=Path,
        help="Directory containing receptor PDB files.",
    )
    parser.add_argument(
        "--ligands", required=True, type=Path,
        help="Directory containing reference ligand files (also used as bond-order templates for AutoDock PDBQT conversion).",
    )
    parser.add_argument(
        "--ligand-format", type=str, default="sdf", dest="ligand_format",
        help="File extension of reference ligand files, e.g. sdf, mol2, pdb, pdbqt (default: sdf).",
    )
    parser.add_argument(
        "--poses", required=True, type=Path,
        help="Directory containing docked pose files.",
    )
    parser.add_argument(
        "--tool", required=True, choices=["autodock", "diffdock", "equibind"],
        help="Docking tool used to generate the poses.",
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="Output directory for PoseBusters results and logs.",
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Number of parallel workers (default: all CPU cores).",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Discard previous results and re-validate all poses.",
    )
    parser.add_argument(
        "--save-interval", type=int, default=100, dest="save_interval",
        help="Save checkpoint every N poses (default: 100).",
    )
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    for name, path in [("--receptors", args.receptors), ("--ligands", args.ligands), ("--poses", args.poses)]:
        if not path.is_dir():
            sys.exit(f"ERROR: {name} directory does not exist: {path}")
    args.output.mkdir(parents=True, exist_ok=True)


# ============================================================================
# PROTEIN FILE RESOLUTION
# ============================================================================

def discover_proteins(receptors_dir: Path) -> list[Path]:
    return sorted(receptors_dir.glob("**/*.pdb"))


def find_protein_file(protein_name: str, receptors_dir: Path, all_pdbs: list[Path]) -> str | None:
    """Multi-strategy protein PDB lookup."""
    # Strategy 1: Exact match candidates
    for suffix in ["", "_protein", "_clean", "_receptor"]:
        candidate = receptors_dir / f"{protein_name}{suffix}.pdb"
        if candidate.exists():
            return str(candidate)

    # Strategy 2: Recursive exact match
    for suffix in ["", "_protein", "_clean", "_receptor"]:
        matches = list(receptors_dir.glob(f"**/{protein_name}{suffix}.pdb"))
        if matches:
            return str(matches[0])

    # Strategy 3: Prefix match (shortest name wins)
    matches = list(receptors_dir.glob(f"**/{protein_name}*.pdb"))
    if matches:
        matches.sort(key=lambda p: len(p.name))
        return str(matches[0])

    # Strategy 4: Substring (case-insensitive)
    for pdb_file in all_pdbs:
        if protein_name.lower() in pdb_file.stem.lower():
            return str(pdb_file)

    # Strategy 5: Normalized matching (- → _)
    norm = protein_name.replace("-", "_").lower()
    for pdb_file in all_pdbs:
        stem_norm = pdb_file.stem.replace("-", "_").lower()
        if norm == stem_norm or norm in stem_norm:
            return str(pdb_file)

    return None


def build_protein_cache(
    protein_names: set[str], receptors_dir: Path, all_pdbs: list[Path],
) -> dict[str, str | None]:
    cache: dict[str, str | None] = {}
    for name in sorted(protein_names):
        pfile = find_protein_file(name, receptors_dir, all_pdbs)
        cache[name] = pfile
        status = f"✓ {Path(pfile).name}" if pfile else "✗ NOT FOUND"
        print(f"  {name}: {status}")
    n_found = sum(1 for v in cache.values() if v is not None)
    n_total = len(cache)
    print(f"\n  Resolved: {n_found}/{n_total} proteins")
    if n_found < n_total:
        print(f"  WARNING: {n_total - n_found} proteins have no PDB — poses will use 'mol' mode (no intermolecular checks)")
    return cache


# ============================================================================
# POSE COLLECTION — tool-specific
# ============================================================================

def collect_autodock_poses(
    poses_dir: Path, ligands_dir: Path, converted_dir: Path,
    ligand_format: str = "sdf",
) -> list[dict]:
    """Discover AutoDock Vina PDBQT results and convert to SDF."""
    poses: list[dict] = []
    pdbqt_files = sorted(poses_dir.glob("**/*_vina_out.pdbqt"))
    if not pdbqt_files:
        print("  WARNING: No *_vina_out.pdbqt files found.")
        return poses

    print(f"  Found {len(pdbqt_files)} PDBQT result files")

    for pdbqt_file in pdbqt_files:
        stem = pdbqt_file.stem.replace("_vina_out", "")
        # Remove scoring suffix like _vina, _vinardo, _ad4
        stem = re.sub(r'_(vina|vinardo|ad4)$', '', stem)
        parts = stem.split("__")
        if len(parts) != 2:
            print(f"  WARNING: Cannot parse protein/ligand from {pdbqt_file.name}, skipping.")
            continue
        protein, ligand = parts

        print(f"  Converting {pdbqt_file.name} → SDF...")
        sdf_files = convert_pdbqt_to_sdf_files(str(pdbqt_file), converted_dir, ligands_dir, ligand_format)

        for sdf_file in sdf_files:
            sdf_path = Path(sdf_file)
            model_num = sdf_path.stem.split("_model")[-1] if "_model" in sdf_path.stem else "1"
            poses.append({
                "method": "autodock",
                "protein": protein,
                "ligand": ligand,
                "pose_file": sdf_file,
                "pose_name": f"{protein}__{ligand}_pose{model_num}",
                "file_format": "sdf",
            })
    return poses


def collect_diffdock_poses(poses_dir: Path) -> list[dict]:
    """Discover DiffDock SDF results from subdirectories."""
    poses: list[dict] = []
    skip_names = {"prepared_proteins", "converted_pdbqt", "prepared_ligands"}

    for subdir in sorted(poses_dir.iterdir()):
        if not subdir.is_dir():
            continue
        if "_ligand__" in subdir.name or subdir.name in skip_names:
            continue
        parts = subdir.name.split("__")
        if len(parts) != 2:
            continue
        ligand, protein = parts

        for sdf_file in sorted(subdir.glob("**/*.sdf")):
            confidence = None
            if "confidence" in sdf_file.stem.lower():
                match = re.search(r'confidence[_-]?([\d.]+)', sdf_file.stem, re.IGNORECASE)
                if match:
                    try:
                        confidence = float(match.group(1))
                    except ValueError:
                        pass
            try:
                rel = sdf_file.relative_to(subdir)
                pose_suffix = str(rel)
            except ValueError:
                pose_suffix = sdf_file.name

            info: dict = {
                "method": "diffdock",
                "protein": protein,
                "ligand": ligand,
                "pose_file": str(sdf_file),
                "pose_name": f"{ligand}__{protein}/{pose_suffix}",
                "file_format": "sdf",
            }
            if confidence is not None:
                info["confidence"] = confidence
            poses.append(info)

    if not poses:
        print("  WARNING: No DiffDock SDF poses found.")
    return poses


def collect_equibind_poses(poses_dir: Path) -> list[dict]:
    """Discover EquiBind SDF results from subdirectories."""
    poses: list[dict] = []

    for subdir in sorted(poses_dir.iterdir()):
        if not subdir.is_dir() or "_pymol" in subdir.name:
            continue
        dir_name = subdir.name
        clean = dir_name.replace("_spatial_sites", "") if dir_name.endswith("_spatial_sites") else dir_name
        parts = clean.split("__")
        if len(parts) != 2:
            continue
        ligand, protein = parts

        for sdf_file in sorted(subdir.glob("**/*.sdf")):
            try:
                rel = sdf_file.relative_to(subdir)
                pose_suffix = str(rel)
            except ValueError:
                pose_suffix = sdf_file.name

            poses.append({
                "method": "equibind",
                "protein": protein,
                "ligand": ligand,
                "pose_file": str(sdf_file),
                "pose_name": f"{ligand}__{protein}/{pose_suffix}",
                "file_format": "sdf",
            })

    if not poses:
        print("  WARNING: No EquiBind SDF poses found.")
    return poses


# ============================================================================
# PDBQT → SDF CONVERSION (AutoDock only)
# ============================================================================

def split_pdbqt_models(pdbqt_file: str) -> list[str]:
    """Split a multi-model PDBQT file into separate model blocks."""
    with open(pdbqt_file, 'r') as f:
        content = f.read()

    models: list[str] = []
    current: list[str] = []
    in_model = False

    for line in content.split('\n'):
        if line.startswith('MODEL'):
            in_model = True
            current = [line]
        elif line.startswith('ENDMDL'):
            current.append(line)
            models.append('\n'.join(current))
            current = []
            in_model = False
        elif in_model:
            current.append(line)

    if not models and content.strip():
        models = [content]
    return models


_RDKIT_READERS = {
    "sdf": lambda p: __import__("rdkit.Chem", fromlist=["Chem"]).MolFromMolFile(str(p), removeHs=True, sanitize=True),
    "mol2": lambda p: __import__("rdkit.Chem", fromlist=["Chem"]).MolFromMol2File(str(p), removeHs=True, sanitize=True),
    "pdb": lambda p: __import__("rdkit.Chem", fromlist=["Chem"]).MolFromPDBFile(str(p), removeHs=True, sanitize=True),
    "pdbqt": lambda p: __import__("rdkit.Chem", fromlist=["Chem"]).MolFromPDBFile(str(p), removeHs=True, sanitize=False),
}


def _find_template_mol(ligand_name: str, ligands_dir: Path, ligand_format: str = "sdf"):
    """Search the ligands directory for a bond-order template in the given format."""
    ext = ligand_format.lower().lstrip(".")
    reader = _RDKIT_READERS.get(ext)
    if reader is None:
        print(f"    WARNING: No RDKit reader for '.{ext}' — skipping template search.")
        return None, None

    for pattern in [f"{ligand_name}.{ext}", f"{ligand_name}_*.{ext}", f"*{ligand_name}*.{ext}"]:
        for match in sorted(ligands_dir.glob(f"**/{pattern}")):
            try:
                mol = reader(match)
                if mol is not None:
                    return mol, match.name
            except Exception:
                pass
    return None, None


def convert_pdbqt_to_sdf_files(
    pdbqt_file: str, output_dir: Path, ligands_dir: Path,
    ligand_format: str = "sdf",
) -> list[str]:
    """Convert a multi-pose PDBQT to individual SDF files."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pdbqt_path = Path(pdbqt_file)
    base_name = pdbqt_path.stem
    models = split_pdbqt_models(pdbqt_file)
    if not models:
        print(f"    Warning: No models in {pdbqt_file}")
        return []

    # Resolve template ligand
    stem_clean = base_name.replace("_vina_out", "")
    stem_clean = re.sub(r'_(vina|vinardo|ad4)$', '', stem_clean)
    parts = stem_clean.split("__")
    template_mol = None
    if len(parts) == 2:
        _, ligand_name = parts
        template_mol, tpl_name = _find_template_mol(ligand_name, ligands_dir, ligand_format)
        if template_mol is not None:
            print(f"    Using bond-order template: {tpl_name}")

    converted: list[str] = []

    for i, model_content in enumerate(models, start=1):
        output_sdf = output_dir / f"{base_name}_model{i}.sdf"
        if output_sdf.exists():
            converted.append(str(output_sdf))
            continue

        temp_pdbqt = output_dir / f"{base_name}_model{i}.pdbqt"
        temp_pdb = output_dir / f"{base_name}_model{i}.pdb"

        with open(temp_pdbqt, 'w') as f:
            f.write(model_content)

        mol_final = None

        # Strategy 1: obabel PDBQT→PDB then RDKit with template
        if template_mol is not None:
            try:
                result = subprocess.run(
                    ['obabel', str(temp_pdbqt), '-O', str(temp_pdb)],
                    capture_output=True, text=True,
                )
                if result.returncode == 0 and temp_pdb.exists():
                    raw_mol = Chem.MolFromPDBFile(str(temp_pdb), removeHs=True, sanitize=False)
                    if raw_mol is not None:
                        try:
                            mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
                            Chem.SanitizeMol(mol_final)
                        except Exception as e:
                            print(f"    Template assignment failed model {i}: {e}")
                            mol_final = None
                if temp_pdb.exists():
                    temp_pdb.unlink()
            except FileNotFoundError:
                pass

        # Strategy 2: direct obabel PDBQT→SDF
        if mol_final is None:
            try:
                result = subprocess.run(
                    ['obabel', str(temp_pdbqt), '-O', str(output_sdf)],
                    capture_output=True, text=True,
                )
                if result.returncode == 0 and output_sdf.exists():
                    if template_mol is not None:
                        try:
                            raw_mol = Chem.MolFromMolFile(str(output_sdf), removeHs=True, sanitize=False)
                            if raw_mol is not None:
                                mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
                                Chem.SanitizeMol(mol_final)
                        except Exception:
                            converted.append(str(output_sdf))
                            if temp_pdbqt.exists():
                                temp_pdbqt.unlink()
                            continue
                    else:
                        converted.append(str(output_sdf))
                        if temp_pdbqt.exists():
                            temp_pdbqt.unlink()
                        continue
            except FileNotFoundError:
                try:
                    raw_mol = Chem.MolFromPDBFile(str(temp_pdbqt), removeHs=True, sanitize=False)
                    if raw_mol is not None and template_mol is not None:
                        mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
                        Chem.SanitizeMol(mol_final)
                    elif raw_mol is not None:
                        try:
                            Chem.SanitizeMol(raw_mol)
                        except Exception:
                            pass
                        mol_final = raw_mol
                except Exception as e:
                    print(f"    RDKit fallback failed model {i}: {e}")

        if mol_final is not None:
            writer = Chem.SDWriter(str(output_sdf))
            writer.write(mol_final)
            writer.close()
            if output_sdf.exists():
                converted.append(str(output_sdf))
        elif not output_sdf.exists():
            print(f"    Warning: Could not convert model {i} from {pdbqt_path.name}")

        if temp_pdbqt.exists():
            temp_pdbqt.unlink()

    return converted


# ============================================================================
# POSEBUSTERS PARALLEL VALIDATION
# ============================================================================

# Worker globals (initialised once per child process)
_worker_buster_dock = None
_worker_buster_mol = None
_worker_protein_cache: dict = {}
_worker_config_mode: str = "dock"


def _init_worker(config_mode: str, protein_cache: dict) -> None:
    global _worker_buster_dock, _worker_buster_mol, _worker_protein_cache, _worker_config_mode
    from posebusters import PoseBusters

    _worker_config_mode = config_mode
    _worker_protein_cache = protein_cache or {}
    if config_mode == "dock":
        _worker_buster_dock = PoseBusters(config="dock")
    _worker_buster_mol = PoseBusters(config="mol")


def _process_single_pose(pose_info: dict):
    """Validate one pose; returns (DataFrame | None, error_msg | None)."""
    global _worker_buster_dock, _worker_buster_mol, _worker_protein_cache, _worker_config_mode

    pose_file = pose_info["pose_file"]
    protein_name = pose_info["protein"]
    config_mode = _worker_config_mode

    if not os.path.exists(pose_file):
        return None, f"File not found: {pose_file}"

    try:
        protein_file = None
        used_mode = config_mode

        if config_mode == "dock":
            protein_file = _worker_protein_cache.get(protein_name)
            if protein_file is not None and Path(protein_file).exists():
                df = _worker_buster_dock.bust(pose_file, None, protein_file, full_report=True)
                used_mode = "dock"
            else:
                df = _worker_buster_mol.bust(pose_file, None, None, full_report=True)
                used_mode = "mol (fallback)"
        else:
            df = _worker_buster_mol.bust(pose_file, None, None, full_report=True)
            used_mode = "mol"

        df["docking_method"] = pose_info["method"]
        df["protein"] = protein_name
        df["ligand"] = pose_info["ligand"]
        df["pose_file"] = pose_info["pose_file"]
        df["pose_name"] = pose_info["pose_name"]
        df["file_format"] = pose_info.get("file_format", "sdf")
        df["protein_file_used"] = protein_file if protein_file else "none"
        df["posebusters_mode"] = used_mode
        if "confidence" in pose_info:
            df["diffdock_confidence"] = pose_info["confidence"]

        return df, None
    except Exception as e:
        return None, f"Error processing {Path(pose_file).name}: {e}"


def _save_checkpoint(
    result_frames: list[pd.DataFrame],
    output_file: Path,
    n_new: int,
    n_prev: int,
    n_pending: int,
    final: bool = False,
) -> None:
    df_out = pd.concat(result_frames, ignore_index=True)
    df_out.to_csv(output_file, index=False)
    tag = "FINAL SAVE" if final else "CHECKPOINT"
    print(f"    [{tag}] {n_new}/{n_pending} new  |  total rows: {len(df_out)}  →  {output_file.name}")


def validate_poses(
    poses_list: list[dict],
    protein_cache: dict[str, str | None],
    output_file: Path,
    save_interval: int = 100,
    overwrite: bool = False,
    num_workers: int | None = None,
) -> pd.DataFrame:
    """Run PoseBusters validation with parallel processing and checkpointing."""
    if not poses_list:
        print("No poses to validate!")
        return pd.DataFrame()

    # Resume / overwrite
    previous_results = pd.DataFrame()
    already_done: set[str] = set()

    if output_file.exists():
        if overwrite:
            print(f"  OVERWRITE — deleting {output_file.name}")
            output_file.unlink()
        else:
            previous_results = pd.read_csv(output_file)
            if "pose_file" in previous_results.columns:
                already_done = set(previous_results["pose_file"].dropna().astype(str))
            print(f"  Loaded {len(previous_results)} previous results ({len(already_done)} unique pose files)")

    pending = [p for p in poses_list if str(p["pose_file"]) not in already_done]
    n_skipped = len(poses_list) - len(pending)
    if n_skipped:
        print(f"  Skipping {n_skipped} already-inspected poses")
    if not pending:
        print("  All poses already inspected.")
        return previous_results

    n_workers = num_workers if num_workers is not None else cpu_count()
    n_workers = min(n_workers, len(pending))
    print(f"\n  Workers: {n_workers}  ({cpu_count()} CPUs available)")
    print(f"  Remaining: {len(pending)}  (save every {save_interval})\n")

    config_mode = "dock"
    all_results: list[pd.DataFrame] = []
    if not previous_results.empty:
        all_results.append(previous_results)

    n_new = n_dock = n_mol_fallback = n_errors = 0

    with Pool(processes=n_workers, initializer=_init_worker, initargs=(config_mode, protein_cache)) as pool:
        for result_df, error_msg in pool.imap_unordered(_process_single_pose, pending):
            if error_msg is not None:
                n_errors += 1
                print(f"  {error_msg}")
                continue
            if result_df is not None:
                mode_val = result_df["posebusters_mode"].iloc[0] if "posebusters_mode" in result_df.columns else ""
                if mode_val == "dock":
                    n_dock += 1
                elif "fallback" in str(mode_val):
                    n_mol_fallback += 1
                all_results.append(result_df)
                n_new += 1

            if n_new % 10 == 0 or n_new == 1:
                print(f"  Processed {n_new}/{len(pending)}  (overall {len(already_done) + n_new}/{len(poses_list)})")
            if n_new > 0 and n_new % save_interval == 0:
                _save_checkpoint(all_results, output_file, n_new, len(already_done), len(pending))

    if all_results:
        _save_checkpoint(all_results, output_file, n_new, len(already_done), len(pending), final=True)

    print(f"\n  Complete:")
    print(f"    Newly processed:  {n_new}")
    print(f"    Previously done:  {n_skipped}")
    print(f"    Total in CSV:     {len(already_done) + n_new}")
    print(f"    Dock mode:        {n_dock}")
    print(f"    Mol fallback:     {n_mol_fallback}")
    print(f"    Errors:           {n_errors}")

    return pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    args = parse_args()
    validate_args(args)

    converted_dir = args.output / "converted_pdbqt"
    converted_dir.mkdir(parents=True, exist_ok=True)

    results_csv = args.output / "posebusters_results.csv"

    print("=" * 80)
    print("PoseBusters CLI — Docked Pose Validation")
    print("=" * 80)
    print(f"  Receptors:     {args.receptors}")
    print(f"  Ligands:       {args.ligands}")
    print(f"  Poses:         {args.poses}")
    print(f"  Tool:          {args.tool}")
    print(f"  Output:        {args.output}")
    print(f"  Ligand format: {args.ligand_format}")
    print(f"  Workers:       {args.workers or 'all CPUs'}")
    print(f"  Overwrite:     {args.overwrite}")
    print(f"  Save interval: {args.save_interval}")

    # ── Step 1: Discover proteins ───────────────────────────────────────
    print(f"\n{'─' * 80}")
    print("Step 1: Discovering receptor PDB files")
    print(f"{'─' * 80}")
    all_pdbs = discover_proteins(args.receptors)
    print(f"  Found {len(all_pdbs)} PDB files in {args.receptors}")
    for pdb in all_pdbs:
        print(f"    {pdb.relative_to(args.receptors)}")

    # ── Step 2: Collect poses ───────────────────────────────────────────
    print(f"\n{'─' * 80}")
    print(f"Step 2: Collecting {args.tool} poses from {args.poses}")
    print(f"{'─' * 80}")

    if args.tool == "autodock":
        poses = collect_autodock_poses(args.poses, args.ligands, converted_dir, args.ligand_format)
    elif args.tool == "diffdock":
        poses = collect_diffdock_poses(args.poses)
    elif args.tool == "equibind":
        poses = collect_equibind_poses(args.poses)
    else:
        sys.exit(f"Unknown tool: {args.tool}")

    if not poses:
        sys.exit("ERROR: No poses found. Check --poses directory and --tool choice.")

    # Summarise
    proteins_in_poses = {p["protein"] for p in poses}
    ligands_in_poses = {p["ligand"] for p in poses}
    print(f"\n  Total poses:    {len(poses)}")
    print(f"  Unique proteins: {len(proteins_in_poses)}")
    print(f"  Unique ligands:  {len(ligands_in_poses)}")

    # ── Step 3: Resolve protein files ───────────────────────────────────
    print(f"\n{'─' * 80}")
    print("Step 3: Resolving protein PDB files")
    print(f"{'─' * 80}")
    protein_cache = build_protein_cache(proteins_in_poses, args.receptors, all_pdbs)

    # Filter out poses whose protein PDB could not be resolved
    missing_proteins = {name for name, pfile in protein_cache.items() if pfile is None}
    if missing_proteins:
        before = len(poses)
        poses = [p for p in poses if p["protein"] not in missing_proteins]
        skipped = before - len(poses)
        print(f"\n  Skipped {skipped} pose(s) for {len(missing_proteins)} unresolved protein(s) — no PDB available")

    # ── Step 4: Run PoseBusters ─────────────────────────────────────────
    print(f"\n{'─' * 80}")
    print("Step 4: Running PoseBusters validation")
    print(f"{'─' * 80}")

    t0 = time.time()
    results_df = validate_poses(
        poses,
        protein_cache=protein_cache,
        output_file=results_csv,
        save_interval=args.save_interval,
        overwrite=args.overwrite,
        num_workers=args.workers,
    )
    elapsed = time.time() - t0

    # ── Summary ─────────────────────────────────────────────────────────
    print(f"\n{'=' * 80}")
    print("SUMMARY")
    print(f"{'=' * 80}")
    if not results_df.empty:
        print(f"  Results CSV:  {results_csv}")
        print(f"  Total rows:   {len(results_df)}")
        print(f"  Elapsed:      {elapsed:.1f}s ({elapsed / 60:.1f} min)")

        if "posebusters_mode" in results_df.columns:
            print(f"\n  Mode breakdown:")
            for mode, count in results_df["posebusters_mode"].value_counts().items():
                print(f"    {mode}: {count}")

        # Quick pass/fail summary
        metadata_cols = {
            'docking_method', 'protein', 'ligand', 'pose_file', 'pose_name',
            'file_format', 'protein_file_used', 'posebusters_mode', 'diffdock_confidence',
            'all_passed',
        }
        exclude_prefixes = ('number_', 'num_')
        test_cols = []
        for col in results_df.columns:
            if col in metadata_cols or col.lower().startswith(exclude_prefixes):
                continue
            uv = set(results_df[col].dropna().unique())
            if uv.issubset({True, False, 1, 0, 1.0, 0.0, 'True', 'False', 'true', 'false'}):
                test_cols.append(col)

        if test_cols:
            for tc in test_cols:
                if results_df[tc].dtype == object:
                    results_df[tc] = results_df[tc].map({'True': True, 'true': True, 'False': False, 'false': False})
                results_df[tc] = results_df[tc].astype(bool)

            # ── Critical tests (physical validity) ──
            CRITICAL_TESTS = [
                "mol_pred_loaded",
                "sanitization",
                "all_atoms_connected",
                "bond_lengths",
                "bond_angles",
                "internal_steric_clash",
                "minimum_distance_to_protein",
                "volume_overlap_with_protein",
            ]
            critical_available = [t for t in CRITICAL_TESTS if t in test_cols]
            if critical_available:
                results_df["critical_passed"] = results_df[critical_available].all(axis=1)
                n_crit = results_df["critical_passed"].sum()
                print(f"\n  Critical tests ({len(critical_available)} of {len(CRITICAL_TESTS)}):")
                for ct in critical_available:
                    rate = results_df[ct].mean() * 100
                    n_ok = results_df[ct].sum()
                    print(f"    {ct:40s}  {n_ok:>5}/{len(results_df)}  ({rate:5.1f}%)")
                print(f"  {'─' * 60}")
                print(f"  Passed ALL critical: {n_crit}/{len(results_df)} ({100 * n_crit / len(results_df):.1f}%)")
            else:
                results_df["critical_passed"] = True

            # ── All tests ──
            results_df["all_passed"] = results_df[test_cols].all(axis=1)
            n_pass = results_df["all_passed"].sum()
            print(f"\n  All tests ({len(test_cols)}):")
            for tc in sorted(test_cols):
                rate = results_df[tc].mean() * 100
                n_ok = results_df[tc].sum()
                print(f"    {tc:40s}  {n_ok:>5}/{len(results_df)}  ({rate:5.1f}%)")
            print(f"  {'─' * 60}")
            print(f"  Passed ALL {len(test_cols)} tests: {n_pass}/{len(results_df)} ({100 * n_pass / len(results_df):.1f}%)")

            # Save updated CSV with both columns
            results_df.to_csv(results_csv, index=False)

            # Bottleneck tests
            print(f"\n  Bottleneck tests (< 99% pass rate):")
            for tc in test_cols:
                rate = results_df[tc].mean() * 100
                if rate < 99.0:
                    print(f"    {tc}: {rate:.1f}%")
    else:
        print("  No results produced.")

    print(f"\n{'=' * 80}")
    print("Done!")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
