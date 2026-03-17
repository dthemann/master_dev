# %% [markdown]
# # Setting Directories

# %%
import sys
import argparse

# ============================================================================
# CLI ARGUMENTS
# ============================================================================
# Usage:
#   python pose_busters.py \
#       --docking-dir autodock:/path/to/autodock_poses \
#       --docking-dir diffdock:/path/to/diffdock_poses \
#       [--proteins-dir /path/to/proteins]  (required when --mode dock) \
#       [--wd /path/to/workspace] [--filter _cleaned] [--overwrite] \
#       [--mode dock] [--output-dir /path/to/output]
#
# --docking-dir   One or more method:path pairs (repeatable). Example:
#                   --docking-dir autodock:/data/autodock_poses
#                   --docking-dir diffdock:/data/diffdock_results
# --proteins-dir  Path to the directory containing protein PDB files.
#                 Required when --mode dock.
# --wd            Working directory (default: /home/manndo/master_dev)
# --filter        Optional substring filter for protein names
# --overwrite     Discard previous results and re-run all poses
# --mode          PoseBusters mode: "dock" (with protein) or "mol" (ligand only)
# --output-dir    Output directory for results (default: <wd>/posebusters_results/<mode>)
# ============================================================================
_parser = argparse.ArgumentParser(
    description="PoseBusters validation of docked poses",
    add_help=False,  # avoid conflict when run as notebook cell
)
_parser.add_argument("--docking-dir", type=str, action="append", dest="docking_dirs",
                     metavar="METHOD:PATH",
                     help=("Docking method and its poses directory as method:path. "
                           "Can be specified multiple times. Example: "
                           "--docking-dir autodock:/data/autodock --docking-dir diffdock:/data/diffdock"))
_parser.add_argument("--proteins-dir", type=str, default=None, dest="proteins_dir",
                     help="Path to directory containing protein PDB files (required for --mode dock)")
_parser.add_argument("--wd", type=str, default="/home/manndo/master_dev",
                     help="Working directory (used for output paths and ligand lookups)")
_parser.add_argument("--filter", type=str, default=None, dest="protein_filter",
                     help="Only include poses whose protein name contains this substring (e.g. _cleaned)")
_parser.add_argument("--overwrite", action="store_true", default=False,
                     help="Discard previous results and re-run all poses from scratch")
_parser.add_argument("--mode", type=str, default="dock", choices=["dock", "mol"],
                     help="PoseBusters config mode: 'dock' (with protein) or 'mol' (ligand only)")
_parser.add_argument("--output-dir", type=str, default=None, dest="output_dir",
                     help="Output directory for results (default: <wd>/posebusters_results/<mode>)")

# parse_known_args tolerates extra flags / notebook runner args
_args, _unknown = _parser.parse_known_args()

wd = _args.wd
PROTEIN_FILTER = _args.protein_filter
OVERWRITE = _args.overwrite
CONFIG_MODE = _args.mode

# ============================================================================
# PARSE --docking-dir arguments into docking_directories dict
# ============================================================================
if not _args.docking_dirs:
    _parser.error("At least one --docking-dir METHOD:PATH is required.")

docking_directories = {}
for entry in _args.docking_dirs:
    if ":" not in entry:
        _parser.error(
            f"Invalid --docking-dir format: '{entry}'. "
            f"Expected METHOD:PATH (e.g. autodock:/data/autodock_poses)"
        )
    method, path = entry.split(":", 1)
    method = method.strip()
    path = path.strip()
    if not method:
        _parser.error(f"Empty method name in --docking-dir '{entry}'")
    if not path:
        _parser.error(f"Empty path in --docking-dir '{entry}'")
    docking_directories[method] = path

# ============================================================================
# PROTEINS DIRECTORY (required for mode "dock")
# ============================================================================
proteins_dir = _args.proteins_dir
if CONFIG_MODE == "dock" and not proteins_dir:
    _parser.error("--proteins-dir is required when --mode dock is used.")

# ============================================================================
# PRINT CONFIGURATION
# ============================================================================
print(f"Working directory: {wd}")
print(f"PoseBusters mode: {CONFIG_MODE}")
print(f"\nDocking directories ({len(docking_directories)}):")
for method_name, method_path in docking_directories.items():
    print(f"  {method_name}: {method_path}")
if proteins_dir:
    print(f"\nProteins directory: {proteins_dir}")
if PROTEIN_FILTER:
    print(f"Protein filter active: only names containing '{PROTEIN_FILTER}'")
else:
    print("No --filter provided — all poses included (no protein-name filtering).")
if OVERWRITE:
    print("OVERWRITE mode is ON — previous results will be discarded.")

# %%
"""
Overview: Pose Counts by Docking Method and Protein-Ligand Combination
This cell analyzes all docking results directories and creates a summary table.
It iterates over `docking_directories` and filters each directory for poses
based on the docking method key.
"""
import pandas as pd
from pathlib import Path
import re


def count_autodock_poses(poses_dir):
    """Count AutoDock Vina poses from PDBQT files"""
    data = []
    base_path = Path(poses_dir)

    for pdbqt_file in base_path.glob("*_vina_out.pdbqt"):
        stem = pdbqt_file.stem.replace("_vina_out", "")
        parts = stem.split("__")

        if len(parts) == 2:
            protein, ligand = parts

            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue

            # Count MODEL entries in PDBQT file
            with open(pdbqt_file, 'r') as f:
                content = f.read()
                pose_count = content.count('MODEL')
                if pose_count == 0:
                    pose_count = 1  # Single pose file without MODEL markers

            data.append({
                "Protein": protein,
                "Ligand": ligand,
                "Poses": pose_count
            })

    return pd.DataFrame(data)


def count_diffdock_poses(poses_dir):
    """Count DiffDock poses from SDF files in subdirectories"""
    data = []
    base_path = Path(poses_dir)

    for subdir in base_path.iterdir():
        if not subdir.is_dir():
            continue
        # Skip helper directories
        if "_ligand__" in subdir.name or subdir.name in ["prepared_proteins", "converted_pdbqt", "prepared_ligands", "Orai"]:
            continue

        parts = subdir.name.split("__")
        if len(parts) == 2:
            ligand, protein = parts

            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue

            # Count all SDF files in subdirectory (including complex_* folders)
            pose_count = len(list(subdir.glob("**/*.sdf")))

            if pose_count > 0:
                data.append({
                    "Protein": protein,
                    "Ligand": ligand,
                    "Poses": pose_count
                })

    return pd.DataFrame(data)


def count_equibind_poses(poses_dir):
    """Count EquiBind poses from SDF files, supporting spatial_sites directory structure.

    Handles directories named like:
      {ligand}__{protein}_spatial_sites/site_XX/pose_XX.sdf
    as well as flat directories:
      {ligand}__{protein}/*.sdf
    """
    data = []
    base_path = Path(poses_dir)

    for subdir in base_path.iterdir():
        if not subdir.is_dir():
            continue
        # Skip pymol export/fixed variants for main count
        if "_pymol" in subdir.name:
            continue

        dir_name = subdir.name

        # Handle spatial_sites suffix: strip it to extract the protein name
        is_spatial = dir_name.endswith("_spatial_sites")
        if is_spatial:
            dir_name_clean = dir_name.replace("_spatial_sites", "")
        else:
            dir_name_clean = dir_name

        parts = dir_name_clean.split("__")
        if len(parts) == 2:
            ligand, protein = parts

            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue

            # Count SDF files recursively (covers site_XX subdirs and flat layouts)
            pose_count = len(list(subdir.glob("**/*.sdf")))

            if pose_count > 0:
                data.append({
                    "Protein": protein,
                    "Ligand": ligand,
                    "Poses": pose_count
                })

    return pd.DataFrame(data)


# Map docking method keys to their counting functions
POSE_COUNTERS = {
    "autodock": count_autodock_poses,
    "diffdock": count_diffdock_poses,
    "equibind": count_equibind_poses,
    "equibind_guided": count_equibind_poses,
    "equibind_exclusion": count_equibind_poses,
    "equibind_docked_poses": count_equibind_poses,
    "equibind_docked_poses_2": count_equibind_poses,
}

# Collect pose counts from all methods defined in docking_directories
print("\n" + "=" * 80)
print("DOCKING POSES OVERVIEW")
print("=" * 80)

method_dfs = {}
for method, poses_dir in docking_directories.items():
    counter_fn = POSE_COUNTERS.get(method)
    if counter_fn is None:
        print(f"  WARNING: No pose counter registered for method '{method}', skipping.")
        continue

    if not Path(poses_dir).exists():
        print(f"  WARNING: Directory not found for '{method}': {poses_dir}, skipping.")
        continue

    df = counter_fn(poses_dir)
    col_name = method.replace(" ", "_").title()
    df = df.rename(columns={"Poses": col_name})
    method_dfs[col_name] = df

# Merge all method DataFrames into one table
df_combined = pd.DataFrame()
for col_name, df in method_dfs.items():
    if df_combined.empty:
        df_combined = df
    else:
        df_combined = df_combined.merge(df, on=["Protein", "Ligand"], how="outer")

if not df_combined.empty:
    method_cols = [c for c in df_combined.columns if c not in ("Protein", "Ligand")]
    df_combined = df_combined.fillna(0)
    for col in method_cols:
        df_combined[col] = df_combined[col].astype(int)
    df_combined["Total"] = df_combined[method_cols].sum(axis=1)
    df_combined = df_combined.sort_values(["Protein", "Ligand"])

    # Compact summary
    print(f"\n  {len(df_combined)} protein-ligand combinations | "
          f"{df_combined['Protein'].nunique()} proteins | "
          f"{df_combined['Ligand'].nunique()} ligands")
    for col in method_cols:
        print(f"  {col}: {df_combined[col].sum():,} poses")
    print(f"  Combined Total: {df_combined['Total'].sum():,} poses")

    # Save to CSV
    output_file = Path(wd + "/posebusters_results") / "pose_counts_overview.csv"
    output_file.parent.mkdir(exist_ok=True)
    df_combined.to_csv(output_file, index=False)
    print(f"  Saved to: {output_file}")

else:
    print("\nNo docking results found!")

# %%
"""
Filter Poses: Keep Only Protein-Ligand Combinations Docked in ALL Methods
==========================================================================
This cell filters the results to include only protein-ligand combinations
that have poses from every docking method defined in `docking_directories`.
"""

# ============================================================================
# POSE COLLECTORS — return list of row dicts for each method
# ============================================================================

def collect_autodock_rows(poses_dir):
    """Collect per-file rows from AutoDock Vina PDBQT outputs."""
    rows = []
    base_path = Path(poses_dir)
    for pdbqt_file in base_path.glob("*_vina_out.pdbqt"):
        stem = pdbqt_file.stem.replace("_vina_out", "")
        parts = stem.split("__")
        if len(parts) == 2:
            protein, ligand = parts
            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue
            with open(pdbqt_file, 'r') as f:
                content = f.read()
                pose_count = content.count('MODEL')
                if pose_count == 0:
                    pose_count = 1
            rows.append({
                "docking_tool": "autodock",
                "protein": protein,
                "ligand": ligand,
                "file_path": str(pdbqt_file),
                "pose_count": pose_count
            })
    return rows


def collect_diffdock_rows(poses_dir):
    """Collect per-directory rows from DiffDock SDF outputs."""
    rows = []
    base_path = Path(poses_dir)
    for subdir in base_path.iterdir():
        if not subdir.is_dir():
            continue
        if "_ligand__" in subdir.name or subdir.name in [
            "prepared_proteins", "converted_pdbqt", "prepared_ligands", "Orai"
        ]:
            continue
        parts = subdir.name.split("__")
        if len(parts) == 2:
            ligand, protein = parts
            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue
            pose_count = len(list(subdir.glob("**/*.sdf")))
            if pose_count > 0:
                rows.append({
                    "docking_tool": "diffdock",
                    "protein": protein,
                    "ligand": ligand,
                    "file_path": str(subdir),
                    "pose_count": pose_count
                })
    return rows


def collect_equibind_rows(poses_dir):
    """Collect per-directory rows from EquiBind SDF outputs."""
    rows = []
    base_path = Path(poses_dir)
    for subdir in base_path.iterdir():
        if not subdir.is_dir() or "_pymol" in subdir.name:
            continue
        dir_name = subdir.name
        if dir_name.endswith("_spatial_sites"):
            dir_name_clean = dir_name.replace("_spatial_sites", "")
        else:
            dir_name_clean = dir_name
        parts = dir_name_clean.split("__")
        if len(parts) == 2:
            ligand, protein = parts
            # Apply optional protein-name filter
            if PROTEIN_FILTER and PROTEIN_FILTER not in protein:
                continue
            pose_count = len(list(subdir.glob("**/*.sdf")))
            if pose_count > 0:
                rows.append({
                    "docking_tool": dir_name.split("__")[0] if "__" in dir_name else "equibind",
                    "protein": protein,
                    "ligand": ligand,
                    "file_path": str(subdir),
                    "pose_count": pose_count
                })
    return rows


# Registry: maps docking_directories keys → collector functions
ROW_COLLECTORS = {
    "autodock": collect_autodock_rows,
    "diffdock": collect_diffdock_rows,
    "equibind": collect_equibind_rows,
    "equibind_guided": collect_equibind_rows,
    "equibind_exclusion": collect_equibind_rows,
    "equibind_docked_poses": collect_equibind_rows,
    "equibind_docked_poses_2": collect_equibind_rows,
}

# ============================================================================
# BUILD filtered_poses_df driven by docking_directories
# ============================================================================
rows = []
for method_key, poses_dir in docking_directories.items():
    collector = ROW_COLLECTORS.get(method_key)
    if collector is None:
        print(f"  WARNING: No row collector for '{method_key}', skipping.")
        continue
    if not Path(poses_dir).exists():
        print(f"  WARNING: Directory not found for '{method_key}': {poses_dir}, skipping.")
        continue
    method_rows = collector(poses_dir)
    # Tag every row with the dictionary key as docking_tool
    for r in method_rows:
        r["docking_tool"] = method_key
    rows.extend(method_rows)

filtered_poses_df = pd.DataFrame(rows)

# ============================================================================
# FILTER: keep only protein-ligand combos present in ALL methods
# ============================================================================
if not filtered_poses_df.empty:
    combos_by_method = {}
    for method in filtered_poses_df["docking_tool"].unique():
        method_df = filtered_poses_df[filtered_poses_df["docking_tool"] == method]
        combos = set(zip(method_df["protein"], method_df["ligand"]))
        combos_by_method[method] = combos

    # Intersection across all methods
    all_combo_sets = list(combos_by_method.values())
    common_combos = all_combo_sets[0]
    for s in all_combo_sets[1:]:
        common_combos = common_combos.intersection(s)

    if common_combos:
        # Filter to common combos only
        filtered_poses_df = filtered_poses_df[
            filtered_poses_df.apply(
                lambda row: (row["protein"], row["ligand"]) in common_combos,
                axis=1
            )
        ].copy()

        print("\n" + "=" * 80)
        print("FILTERED SUBSET (combinations present in all methods)")
        print("=" * 80)
        print(f"  Common protein-ligand combinations: {len(common_combos)}")

        methods_in_df = sorted(filtered_poses_df["docking_tool"].unique())
        grand_total_all = 0
        for method in methods_in_df:
            method_df = filtered_poses_df[filtered_poses_df["docking_tool"] == method]
            total = method_df['pose_count'].sum()
            grand_total_all += total
            print(f"  {method}: {len(method_df)} combos, {total:,} poses")
        print(f"  Total: {grand_total_all:,} poses across {len(methods_in_df)} method(s)")
    else:
        print("\nWARNING: No protein-ligand combinations found in all methods!")
        print("Cannot filter further. Keeping current filtered_poses_df.")
else:
    print("filtered_poses_df is empty. Check docking_directories paths.")

# %%
"""
PoseBusters Analysis - Validate docked poses from the FILTERED SUBSET only
Checks: molecule validity, stereochemistry, bond geometry, ring flatness, etc.

NOTE: This cell uses the filtered_poses_df from the previous cell to only validate
poses that match the defined prefix/postfix filter criteria.
"""
import os
import subprocess
import re
import pandas as pd
from pathlib import Path
from posebusters import PoseBusters

# ============================================================================
# DISCOVER AVAILABLE PROTEIN FILES
# ============================================================================
proteins_base = Path(proteins_dir)
all_protein_pdbs = sorted(proteins_base.glob("**/*.pdb"))

print(f"\nProtein PDB files: {len(all_protein_pdbs)} found in {proteins_base}")


def find_protein_file(protein_name: str) -> str | None:
    """
    Find the PDB file for a given protein name by searching the proteins directory.
    """
    # Strategy 1: Exact match
    exact_candidates = [
        proteins_base / f"{protein_name}.pdb",
        proteins_base / f"{protein_name}_protein.pdb",
        proteins_base / f"{protein_name}_clean.pdb",
        proteins_base / f"{protein_name}_receptor.pdb",
    ]
    for candidate in exact_candidates:
        if candidate.exists():
            return str(candidate)

    # Strategy 2: Recursive exact match
    for candidate_name in [
        f"{protein_name}.pdb",
        f"{protein_name}_protein.pdb",
        f"{protein_name}_clean.pdb",
        f"{protein_name}_receptor.pdb",
    ]:
        matches = list(proteins_base.glob(f"**/{candidate_name}"))
        if matches:
            return str(matches[0])

    # Strategy 3: Recursive glob with protein name as prefix
    matches = list(proteins_base.glob(f"**/{protein_name}*.pdb"))
    if matches:
        matches.sort(key=lambda p: len(p.name))
        return str(matches[0])

    # Strategy 4: Fuzzy - check if protein_name is contained in any PDB filename
    for pdb_file in all_protein_pdbs:
        if protein_name.lower() in pdb_file.stem.lower():
            return str(pdb_file)

    # Strategy 5: Normalized matching
    protein_name_normalized = protein_name.replace("-", "_").lower()
    for pdb_file in all_protein_pdbs:
        pdb_stem_normalized = pdb_file.stem.replace("-", "_").lower()
        if protein_name_normalized == pdb_stem_normalized:
            return str(pdb_file)
        if protein_name_normalized in pdb_stem_normalized:
            return str(pdb_file)

    return None


# Pre-build a protein file lookup cache for all unique protein names
_protein_file_cache = {}

if 'filtered_poses_df' in dir() and not filtered_poses_df.empty:
    unique_proteins = filtered_poses_df["protein"].unique()
    for pname in sorted(unique_proteins):
        _protein_file_cache[pname] = find_protein_file(pname)

    n_found = sum(1 for v in _protein_file_cache.values() if v is not None)
    n_missing = sum(1 for v in _protein_file_cache.values() if v is None)
    print(f"  Protein resolution: {n_found}/{len(unique_proteins)} matched")
    if n_missing > 0:
        missing = [k for k, v in _protein_file_cache.items() if v is None]
        print(f"  WARNING: {n_missing} proteins not found: {missing}")
        print(f"  These poses will fall back to 'mol' mode (no intermolecular checks)")

# ============================================================================
# PDBQT CONVERSION FUNCTIONS
# ============================================================================

def split_pdbqt_models(pdbqt_file: str) -> list[str]:
    """Split a multi-model PDBQT file into separate model blocks."""
    with open(pdbqt_file, 'r') as f:
        content = f.read()

    models = []
    current_model = []
    in_model = False

    for line in content.split('\n'):
        if line.startswith('MODEL'):
            in_model = True
            current_model = [line]
        elif line.startswith('ENDMDL'):
            current_model.append(line)
            models.append('\n'.join(current_model))
            current_model = []
            in_model = False
        elif in_model:
            current_model.append(line)

    if not models and content.strip():
        models = [content]

    return models

def convert_pdbqt_to_sdf_files(pdbqt_file: str, output_dir: Path) -> list[str]:
    """
    Convert a multi-pose PDBQT file to multiple SDF files.
    Uses RDKit's AssignBondOrdersFromTemplate to fix bond orders.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolops

    pdbqt_path = Path(pdbqt_file)
    base_name = pdbqt_path.stem

    models = split_pdbqt_models(pdbqt_file)
    if not models:
        print(f"  Warning: No models found in {pdbqt_file}")
        return []

    # Try to find the original ligand SDF to use as bond-order template
    stem_clean = base_name.replace("_vina_out", "")
    parts = stem_clean.split("__")
    template_mol = None

    if len(parts) == 2:
        protein_name, ligand_name = parts
        ligand_search_dirs = [
            Path(wd) / "ligands",
            Path(wd) / "Ligands",
            Path(wd) / "docking_ready_mgltools",
            Path(wd) / "docking_ready_mgltools" / "ligands",
            Path(wd),
        ]

        for search_dir in ligand_search_dirs:
            if not search_dir.exists():
                continue
            for pattern in [
                f"{ligand_name}.sdf",
                f"{ligand_name}_*.sdf",
                f"*{ligand_name}*.sdf",
            ]:
                matches = list(search_dir.glob(pattern))
                if matches:
                    try:
                        template_mol = Chem.MolFromMolFile(str(matches[0]), removeHs=True, sanitize=True)
                        if template_mol is not None:
                            print(f"    Using bond-order template: {matches[0].name}")
                            break
                    except:
                        pass
            if template_mol is not None:
                break

    converted_files = []

    for i, model_content in enumerate(models, start=1):
        temp_pdbqt = output_dir / f"{base_name}_model{i}.pdbqt"
        temp_pdb = output_dir / f"{base_name}_model{i}.pdb"
        output_sdf = output_dir / f"{base_name}_model{i}.sdf"

        if output_sdf.exists():
            converted_files.append(str(output_sdf))
            continue

        with open(temp_pdbqt, 'w') as f:
            f.write(model_content)

        mol_final = None

        # Strategy 1: obabel PDBQT→PDB, then RDKit with template
        if template_mol is not None:
            try:
                result = subprocess.run(
                    ['obabel', str(temp_pdbqt), '-O', str(temp_pdb)],
                    capture_output=True, text=True
                )
                if result.returncode == 0 and temp_pdb.exists():
                    raw_mol = Chem.MolFromPDBFile(str(temp_pdb), removeHs=True, sanitize=False)
                    if raw_mol is not None:
                        try:
                            mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
                            Chem.SanitizeMol(mol_final)
                        except Exception as e:
                            print(f"    Template assignment failed for model {i}: {e}")
                            mol_final = None
                if temp_pdb.exists():
                    temp_pdb.unlink()
            except FileNotFoundError:
                pass

        # Strategy 2: Direct obabel PDBQT→SDF (fallback)
        if mol_final is None:
            try:
                result = subprocess.run(
                    ['obabel', str(temp_pdbqt), '-O', str(output_sdf)],
                    capture_output=True, text=True
                )
                if result.returncode == 0 and output_sdf.exists():
                    if template_mol is not None:
                        try:
                            raw_mol = Chem.MolFromMolFile(str(output_sdf), removeHs=True, sanitize=False)
                            if raw_mol is not None:
                                mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw_mol)
                                Chem.SanitizeMol(mol_final)
                        except:
                            converted_files.append(str(output_sdf))
                            if temp_pdbqt.exists():
                                temp_pdbqt.unlink()
                            continue
                    else:
                        converted_files.append(str(output_sdf))
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
                        except:
                            pass
                        mol_final = raw_mol
                except Exception as e:
                    print(f"    RDKit fallback failed for model {i}: {e}")

        if mol_final is not None:
            writer = Chem.SDWriter(str(output_sdf))
            writer.write(mol_final)
            writer.close()
            if output_sdf.exists():
                converted_files.append(str(output_sdf))
        elif not output_sdf.exists():
            print(f"    Warning: Could not convert model {i} from {pdbqt_path.name}")

        if temp_pdbqt.exists():
            temp_pdbqt.unlink()

    return converted_files


# ============================================================================
# POSE COLLECTION — method-specific logic for gathering individual pose files
# ============================================================================

def collect_autodock_pose_files(row, converted_dir):
    """Collect individual pose files from an AutoDock Vina PDBQT result."""
    poses = []
    file_path = Path(row["file_path"])
    protein = row["protein"]
    ligand = row["ligand"]
    method_key = row["docking_tool"]

    if file_path.exists() and file_path.suffix == ".pdbqt":
        print(f"  Converting {file_path.name} to SDF...")
        sdf_files = convert_pdbqt_to_sdf_files(str(file_path), converted_dir)
        for sdf_file in sdf_files:
            sdf_path = Path(sdf_file)
            model_num = sdf_path.stem.split("_model")[-1] if "_model" in sdf_path.stem else "1"
            poses.append({
                "method": method_key,
                "protein": protein,
                "ligand": ligand,
                "pose_file": sdf_file,
                "pose_name": f"{protein}__{ligand}_pose{model_num}",
                "original_pdbqt": str(file_path),
                "file_format": "sdf"
            })
    return poses


def collect_diffdock_pose_files(row, converted_dir):
    """Collect individual pose files from a DiffDock result directory."""
    poses = []
    file_path = Path(row["file_path"])
    protein = row["protein"]
    ligand = row["ligand"]
    method_key = row["docking_tool"]

    if file_path.is_dir():
        for sdf_file in sorted(file_path.glob("**/*.sdf")):
            confidence = None
            if "confidence" in sdf_file.stem.lower():
                match = re.search(r'confidence[_-]?([\d.]+)', sdf_file.stem, re.IGNORECASE)
                if match:
                    try:
                        confidence = float(match.group(1))
                    except ValueError:
                        pass
            try:
                relative_path = sdf_file.relative_to(file_path)
                pose_name_suffix = str(relative_path)
            except ValueError:
                pose_name_suffix = sdf_file.name

            pose_info = {
                "method": method_key,
                "protein": protein,
                "ligand": ligand,
                "pose_file": str(sdf_file),
                "pose_name": f"{ligand}__{protein}/{pose_name_suffix}",
                "file_format": "sdf"
            }
            if confidence is not None:
                pose_info["confidence"] = confidence
            poses.append(pose_info)
    return poses


def collect_equibind_pose_files(row, converted_dir):
    """Collect individual pose files from an EquiBind result directory."""
    poses = []
    file_path = Path(row["file_path"])
    protein = row["protein"]
    ligand = row["ligand"]
    method_key = row["docking_tool"]

    if file_path.is_dir():
        for sdf_file in sorted(file_path.glob("**/*.sdf")):
            try:
                relative_path = sdf_file.relative_to(file_path)
                pose_name_suffix = str(relative_path)
            except ValueError:
                pose_name_suffix = sdf_file.name

            poses.append({
                "method": method_key,
                "protein": protein,
                "ligand": ligand,
                "pose_file": str(sdf_file),
                "pose_name": f"{ligand}__{protein}/{pose_name_suffix}",
                "file_format": "sdf"
            })
    return poses


# Registry: maps docking_directories keys → pose file collectors
POSE_FILE_COLLECTORS = {
    "autodock": collect_autodock_pose_files,
    "diffdock": collect_diffdock_pose_files,
    "equibind": collect_equibind_pose_files,
    "equibind_guided": collect_equibind_pose_files,
    "equibind_exclusion": collect_equibind_pose_files,
    "equibind_docked_poses": collect_equibind_pose_files,
    "equibind_docked_poses_2": collect_equibind_pose_files,
}


def collect_poses_from_filtered_subset(filtered_df: pd.DataFrame) -> list[dict]:
    """
    Collect all individual pose files from the filtered subset DataFrame.
    Dispatches to the correct collector based on the docking_tool column,
    which matches the keys in docking_directories.

    Args:
        filtered_df: DataFrame with columns: docking_tool, protein, ligand, file_path, pose_count

    Returns:
        List of dicts with pose information ready for PoseBusters analysis
    """
    all_poses = []

    for _, row in filtered_df.iterrows():
        method_key = row["docking_tool"]
        collector = POSE_FILE_COLLECTORS.get(method_key)
        if collector is None:
            print(f"  WARNING: No pose file collector for method '{method_key}', skipping.")
            continue
        all_poses.extend(collector(row, converted_dir))

    return all_poses

# %% [markdown]
# # PoseBuster Config

# %%
# ============================================================================
# CONFIGURATION
# ============================================================================
# CONFIG_MODE is set via --mode CLI flag (see top of script)
# OVERWRITE is set via --overwrite CLI flag (see top of script)

# Output directory for results
if _args.output_dir:
    output_dir = Path(_args.output_dir)
else:
    output_base_dir = Path(wd) / "posebusters_results"
    output_base_dir.mkdir(exist_ok=True)
    output_dir = output_base_dir / CONFIG_MODE

output_dir.mkdir(parents=True, exist_ok=True)
print(f"Output directory: {output_dir}")

# Temporary directory for converted PDBQT files
converted_dir = output_dir / "converted_pdbqt"
converted_dir.mkdir(exist_ok=True)

# Incremental save interval (save results every N poses)
SAVE_INTERVAL = 100

# %% [markdown]
# # PoseBuster Function

# %%
def analyze_poses_with_posebusters(
    poses_list: list[dict],
    config: str = "mol",
    output_file: Path | None = None,
    save_interval: int = 100,
    overwrite: bool = False,
) -> pd.DataFrame:
    """
    Run PoseBusters validation on a list of poses with incremental checkpointing.

    - Saves results to `output_file` every `save_interval` poses.
    - On restart, loads the existing CSV and skips poses already processed
      (matched by `pose_file`), unless `overwrite=True`.
    - If `overwrite=True`, deletes any existing results and starts fresh.

    Args:
        poses_list:    List of dicts with 'pose_file', 'method', 'protein', 'ligand' keys
        config:        PoseBusters config mode ("mol" or "dock")
        output_file:   Path for the results CSV (used for checkpointing)
        save_interval: Save accumulated results every N poses (default 100)
        overwrite:     If True, ignore/delete previous results and re-run everything

    Returns:
        DataFrame with PoseBusters results merged with pose metadata
    """
    if not poses_list:
        print("No poses to analyze!")
        return pd.DataFrame()

    # ------------------------------------------------------------------
    # Handle existing results: resume or overwrite
    # ------------------------------------------------------------------
    previous_results = pd.DataFrame()
    already_done: set[str] = set()

    if output_file is not None and output_file.exists():
        if overwrite:
            print(f"  OVERWRITE enabled — deleting previous results: {output_file.name}")
            output_file.unlink()
        else:
            previous_results = pd.read_csv(output_file)
            if "pose_file" in previous_results.columns:
                already_done = set(previous_results["pose_file"].dropna().astype(str))
            print(f"  Loaded {len(previous_results)} previous results from {output_file.name}")
            print(f"  {len(already_done)} unique pose files already inspected — will skip them")

    # Filter out already-processed poses
    pending_poses = [
        p for p in poses_list if str(p["pose_file"]) not in already_done
    ]

    n_skipped = len(poses_list) - len(pending_poses)
    if n_skipped > 0:
        print(f"  Skipping {n_skipped} already-inspected poses")
    if not pending_poses:
        print("  All poses already inspected — nothing to do.")
        return previous_results

    # ------------------------------------------------------------------
    # Initialise PoseBusters
    # ------------------------------------------------------------------
    if config == "dock":
        print("Initializing PoseBusters in 'dock' mode (with protein)...")
        buster_dock = PoseBusters(config="dock")
    else:
        print("Initializing PoseBusters in 'mol' mode (ligand only)...")
    buster_mol = PoseBusters(config="mol")  # always available as fallback

    all_results: list[pd.DataFrame] = []
    if not previous_results.empty:
        all_results.append(previous_results)

    total = len(pending_poses)
    n_new = 0      # newly processed in this run
    n_dock = 0
    n_mol_fallback = 0
    n_errors = 0

    print(f"\n  {total} poses remaining to process (save every {save_interval})...\n")

    for idx, pose_info in enumerate(pending_poses, 1):
        pose_file = pose_info["pose_file"]
        protein_name = pose_info["protein"]

        # Progress indicator
        if idx % 10 == 0 or idx == 1:
            print(f"  Processing pose {idx}/{total}  "
                  f"(overall {len(already_done) + n_new + 1}/{len(poses_list)})...")

        if not os.path.exists(pose_file):
            print(f"  Warning: File not found: {pose_file}")
            continue

        try:
            protein_file = None
            used_mode = config

            if config == "dock":
                protein_file = _protein_file_cache.get(protein_name)
                if protein_file is None:
                    protein_file = find_protein_file(protein_name)

                if protein_file is not None and Path(protein_file).exists():
                    df = buster_dock.bust(
                        pose_file, None, protein_file, full_report=True
                    )
                    n_dock += 1
                    used_mode = "dock"
                else:
                    print(f"  Warning: No protein PDB for '{protein_name}' "
                          f"— using 'mol' mode for {Path(pose_file).name}")
                    df = buster_mol.bust(pose_file, None, None, full_report=True)
                    n_mol_fallback += 1
                    used_mode = "mol (fallback)"
            else:
                df = buster_mol.bust(pose_file, None, None, full_report=True)
                used_mode = "mol"

            # Attach metadata
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

            all_results.append(df)
            n_new += 1

        except Exception as e:
            n_errors += 1
            print(f"  Error processing {Path(pose_file).name}: {e}")

        # ----------------------------------------------------------
        # Checkpoint: save every save_interval new results
        # ----------------------------------------------------------
        if output_file is not None and n_new > 0 and n_new % save_interval == 0:
            _save_checkpoint(all_results, output_file, n_new, len(already_done), total)

    # Final save
    if output_file is not None and all_results:
        _save_checkpoint(all_results, output_file, n_new, len(already_done), total, final=True)

    # Summary
    print(f"\n  Processing complete:")
    print(f"    Newly processed:  {n_new}")
    print(f"    Previously done:  {n_skipped}")
    print(f"    Total in CSV:     {len(already_done) + n_new}")
    if config == "dock":
        print(f"    Dock mode (with protein): {n_dock}")
        print(f"    Mol mode (fallback):      {n_mol_fallback}")
    print(f"    Errors: {n_errors}")

    if all_results:
        return pd.concat(all_results, ignore_index=True)
    else:
        return pd.DataFrame()


def _save_checkpoint(
    result_frames: list[pd.DataFrame],
    output_file: Path,
    n_new: int,
    n_prev: int,
    n_pending: int,
    final: bool = False,
):
    """Concatenate accumulated frames and write to CSV."""
    df_out = pd.concat(result_frames, ignore_index=True)
    df_out.to_csv(output_file, index=False)
    tag = "FINAL SAVE" if final else "CHECKPOINT"
    print(f"    [{tag}] {n_new}/{n_pending} new poses saved  |  "
          f"total rows in CSV: {len(df_out)}  →  {output_file.name}")

# %%
# ============================================================================
# MAIN ANALYSIS - USING FILTERED SUBSET ONLY
# ============================================================================
results_csv_path = output_dir / "posebusters_filtered_results.csv"

# Check if filtered_poses_df exists from previous cell
if 'filtered_poses_df' not in dir() or filtered_poses_df.empty:
    print("\nERROR: No filtered poses available!")
    print("Please run the filtering cell first to create filtered_poses_df")
    results_df = pd.DataFrame()
else:

    # Collect all poses from the filtered subset
    print("\nCollecting individual pose files...")
    all_poses = collect_poses_from_filtered_subset(filtered_poses_df)
    print(f"  {len(all_poses)} poses to analyse")

    # Run PoseBusters analysis (with checkpointing + resume)
    print(f"\nRunning PoseBusters validation (mode='{CONFIG_MODE}', overwrite={OVERWRITE})...")

    results_df = analyze_poses_with_posebusters(
        all_poses,
        config=CONFIG_MODE,
        output_file=results_csv_path,
        save_interval=SAVE_INTERVAL,
        overwrite=OVERWRITE,
    )

    # Final summary
    if not results_df.empty:
        print(f"\n  Results: {len(results_df)} rows saved to {results_csv_path}")
    else:
        print("\n  No results to save!")

print("\n" + "=" * 80)
print("Done!")
print("=" * 80)
print(f"\nResults saved to: {output_dir}")
print("To analyse & visualise results, run:")
print(f"  python pose_busters_analysis.py --results-dir {output_dir}")


