#!/usr/bin/env python3
"""
Automated Docking Preparation Script
Converting XYZ ligands and validating Orai1 protein structures
For use with AutoDock Vina, Equibind, and Diffdock
"""

import argparse
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# Default number of parallel workers (capped to avoid overwhelming the system)
DEFAULT_MAX_WORKERS: int = min(os.cpu_count() or 4, 8)


def discover_input_files(folder: Path, contents: str, recursive: bool = False) -> List[str]:
    """Return sorted list of files from *folder* matching the requested contents."""

    if not folder.exists() or not folder.is_dir():
        raise FileNotFoundError(f"Input directory not found: {folder}")

    contents = contents.lower()
    patterns = {
        "ligands": ("*.xyz", "*.sdf", "*.mol2", "*.pdb"),
        "proteins": ("*.pdb",),
    }

    if contents not in patterns:
        raise ValueError(f"Unknown contents type '{contents}'. Use 'ligands' or 'proteins'.")

    searcher = folder.rglob if recursive else folder.glob

    def iter_matches() -> Iterable[Path]:
        for pattern in patterns[contents]:
            for match in searcher(pattern):
                if match.is_file():
                    yield match.resolve()

    files = sorted({str(path) for path in iter_matches()})

    if not files:
        printable = folder if folder.exists() else f"{folder} (missing)"
        raise FileNotFoundError(
            f"No {contents} files found in {printable}. Ensure the folder contains the correct formats."
        )

    return files


DEFAULT_LIGAND_FILES = [
    "2abp-nh2-OPT.xyz",
    "gsk7975a-deprot-OPT.xyz",
    "Synta-66-OPT-Singlet.xyz",
]

DEFAULT_PROTEIN_FILES = [
    "Orai1WT-START-Fr0.pdb",
    "Orai1WT-MDSnap-Fr300.pdb",
    "Orai1WT-MDSnap-Fr400.pdb",
    "Orai1WT-MDSnap-Fr499.pdb",
]


def _printer(enabled: bool):
    return print if enabled else (lambda *args, **kwargs: None)


def _run_parallel(
    func,
    items,
    *,
    max_workers: int = DEFAULT_MAX_WORKERS,
    desc: str = "items",
    verbose: bool = True,
):
    """Execute *func* over *items* in parallel using threads.

    *func* receives a single item and must return a (key, value) tuple.
    Returns a dict of {key: value} preserving all results.
    Exceptions per-item are caught and stored as (key, exception).
    """
    log = _printer(verbose)
    results = {}
    workers = min(max_workers, len(items)) if items else 1
    log(f"⚡ Processing {len(items)} {desc} with {workers} parallel workers")

    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_item = {executor.submit(func, item): item for item in items}
        for future in as_completed(future_to_item):
            item = future_to_item[future]
            try:
                key, value = future.result()
                results[key] = value
            except Exception as exc:
                log(f"✗ Parallel task failed for {item}: {exc}")
                results[item] = exc
    return results


def _strip_extension_from_stem(stem: str, ext: str) -> str:
    """Remove target extension from stem if it already ends with it.
    
    This prevents double extensions like 'file.pdbqt.pdbqt' when the input
    file already has the target extension embedded in its name.
    """
    if not ext:
        return stem
    ext_lower = ext.lower().lstrip(".")
    while stem.lower().endswith("." + ext_lower):
        stem = stem[: -(len(ext_lower) + 1)]
    return stem


def _build_output_path(
    source: Path,
    *,
    output_dir: Optional[Path],
    base_postfix: str = "",
    default_suffix: str = "",
    apply_process_suffix: bool = True,
    extension: Optional[str] = None,
    tool_postfix: str = "",
    tool_prefix: str = "",
) -> Path:
    """Construct an output path honoring destination directory and prefixes/postfixes."""

    dest = output_dir if output_dir else source.parent
    dest.mkdir(parents=True, exist_ok=True)

    suffix = ""
    if apply_process_suffix and default_suffix:
        suffix += default_suffix
    if tool_postfix:
        suffix += tool_postfix
    if base_postfix:
        suffix += base_postfix

    ext = extension if extension is not None else source.suffix
    if ext and not ext.startswith("."):
        ext = "." + ext

    # Start with the stem (filename without its own suffix)
    stem = source.stem
    
    # Remove target extension from stem if already present (prevents double extensions)
    stem = _strip_extension_from_stem(stem, ext)

    # Add prefix to the beginning of the stem
    if tool_prefix:
        stem = tool_prefix + stem

    return dest / f"{stem}{suffix}{ext}"


def _calculate_box_from_pdb(
    pdb_path: Path,
    padding: float = 2.0,
) -> Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]]:
    """Compute bounding-box center and size (with padding) from ATOM records."""

    min_x = min_y = min_z = float("inf")
    max_x = max_y = max_z = float("-inf")
    atom_count = 0

    try:
        with pdb_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                min_x = min(min_x, x)
                min_y = min(min_y, y)
                min_z = min(min_z, z)
                max_x = max(max_x, x)
                max_y = max(max_y, y)
                max_z = max(max_z, z)
                atom_count += 1
    except FileNotFoundError:
        return None

    if atom_count == 0:
        return None

    span_x = max_x - min_x
    span_y = max_y - min_y
    span_z = max_z - min_z

    center = (
        (max_x + min_x) / 2.0,
        (max_y + min_y) / 2.0,
        (max_z + min_z) / 2.0,
    )

    size = (
        span_x + padding,
        span_y + padding,
        span_z + padding,
    )

    return center, size


def _write_box_file(
    box_path: Path,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
) -> None:
    """Write a minimal Vina-style box.txt file."""

    cx, cy, cz = center
    sx, sy, sz = size
    content = (
        f"center_x = {cx:.3f}\n"
        f"center_y = {cy:.3f}\n"
        f"center_z = {cz:.3f}\n"
        f"size_x = {sx:.1f}\n"
        f"size_y = {sy:.1f}\n"
        f"size_z = {sz:.1f}\n"
    )
    box_path.write_text(content)


def fix_pdb_terminal_residues(
    input_pdb: str,
    output_pdb: Optional[str] = None,
    *,
    problematic_atoms: Optional[Sequence[str]] = None,
    verbose: bool = True,
    output_dir: Optional[Path] = None,
    custom_postfix: str = "",
    process_postfixes: bool = True,
) -> str:
    """Remove problematic terminal atoms and return repaired file path."""

    log = _printer(verbose)
    input_path = Path(input_pdb)
    if output_pdb:
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)
    else:
        output_path = _build_output_path(
            input_path,
            output_dir=output_dir,
            base_postfix=custom_postfix,
            default_suffix="_REPAIRED",
            apply_process_suffix=process_postfixes,
            extension=input_path.suffix,
        )

    atoms_to_remove = set(
        problematic_atoms
        or ["CAY", "HY1", "HY2", "HY3", "CY", "OY", "CG2", "CG1"]
    )

    with input_path.open("r") as handle:
        lines = handle.readlines()

    log(f"Processing {input_path} ({len(lines)} total lines)...")
    fixed_lines: List[str] = []
    kept = removed = 0

    for line in lines:
        if not line.startswith("ATOM"):
            fixed_lines.append(line)
            continue

        try:
            atom_name = line[12:16].strip()
            res_num = int(line[22:26].strip())
        except ValueError:
            fixed_lines.append(line)
            continue

        if res_num == 1 and atom_name in atoms_to_remove:
            removed += 1
            continue

        fixed_lines.append(line)
        kept += 1

    with output_path.open("w") as handle:
        handle.writelines(fixed_lines)

    log(
        f"Terminal repair results — kept: {kept} atoms, removed: {removed}. Output: {output_path}"
    )

    return output_path.as_posix()


def _extract_problematic_residues(error_text: str) -> Optional[List[str]]:
    """
    Extract problematic residues from Meeko error messages.
    Parses strings like "Input residues {'M:69': 'HSP', 'M:113': 'HSE', ...}"
    Returns list of residue specs like ['M:69', 'M:113', ...]
    """
    
    # Pattern to match the dictionary of problematic residues
    match = re.search(r"Input residues\s+\{([^}]+)\}", error_text)
    if not match:
        return None
    
    residue_dict_str = "{" + match.group(1) + "}"
    try:
        # Evaluate the dictionary string
        residue_dict = eval(residue_dict_str)
        return sorted(list(residue_dict.keys()))
    except Exception:
        return None


def _remove_residues_from_pdb(
    input_pdb: str,
    residue_specs: List[str],
    output_pdb: str,
    verbose: bool = True,
) -> bool:
    """
    Remove specified residues from a PDB file.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file
    residue_specs : List[str]
        List of residue specs to remove (e.g., ['M:69', 'M:113'])
    output_pdb : str
        Path to output PDB file
    verbose : bool
        Print progress messages
    
    Returns:
    --------
    bool : True if successful, False otherwise
    """
    
    log = _printer(verbose)
    
    # Parse residue specs into (chain, residue_number)
    residues_to_remove = set()
    for spec in residue_specs:
        parts = spec.split(':')
        if len(parts) == 2:
            chain = parts[0]
            try:
                res_num = int(parts[1])
                residues_to_remove.add((chain, res_num))
            except ValueError:
                log(f"⚠ Could not parse residue spec: {spec}")
                continue
    
    log(f"Removing {len(residues_to_remove)} problematic residue(s): {residue_specs}")
    
    try:
        input_path = Path(input_pdb)
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with input_path.open("r") as infile:
            lines = infile.readlines()
        
        kept_lines = []
        removed_count = 0
        
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                # Keep non-atom records
                kept_lines.append(line)
                continue
            
            try:
                chain = line[21].strip() if len(line) > 21 else ""
                res_num = int(line[22:26].strip()) if line[22:26].strip() else 0
                
                if (chain, res_num) in residues_to_remove:
                    removed_count += 1
                    continue
            except (ValueError, IndexError):
                pass
            
            kept_lines.append(line)
        
        with output_path.open("w") as outfile:
            outfile.writelines(kept_lines)
        
        log(f"✓ Removed {removed_count} atoms from problematic residues")
        return True
        
    except Exception as exc:
        log(f"✗ Failed to remove residues: {exc}")
        return False


def convert_protein_with_meeko(
    protein_pdb: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    padding: float = 2.0,
    verbose: bool = True,
    add_tool_postfix: bool = False,
    use_converter_prefix: bool = False,
) -> Tuple[str, Optional[str]]:
    """
    Convert a protein PDB to PDBQT using Meeko with automatic retry and fallback strategies.
    
    Retry strategy:
    1. First attempt: Standard Meeko conversion
    2. Second attempt: Retry with extended options (box_enveloping, allow_bad_res)
    3. Third attempt: Remove problematic residues and retry
    
    File extensions indicate conversion mode:
    - .pdbqt: Standard successful conversion
    - _retry1.pdbqt: First retry with extended options
    - _retry2.pdbqt: Second retry after removing problematic residues
    """

    log = _printer(verbose)
    input_path = Path(protein_pdb)
    
    # Build output directory
    dest = output_dir if output_dir else input_path.parent
    dest.mkdir(parents=True, exist_ok=True)
    
    # Build the base filename (stem) without any extensions
    # Strip .pdbqt if already present to prevent double extensions
    stem = input_path.stem
    stem = _strip_extension_from_stem(stem, ".pdbqt")
    
    # Build prefix/postfix/suffix
    prefix = "mko_" if use_converter_prefix else ""
    postfix = "_mko" if add_tool_postfix else ""
    suffix = "_PDBQT" if process_postfixes else ""
    
    # Build base name without extension
    base_name = f"{prefix}{stem}{suffix}{postfix}{custom_postfix}"
    
    # ATTEMPT 1: Standard Meeko conversion
    log(f"Running Meeko receptor prep for {input_path.name} → {base_name}.pdbqt")
    output_path_no_ext = dest / base_name
    expected_output = dest / f"{base_name}.pdbqt"
    
    cmd = [
        "mk_prepare_receptor.py",
        "-i",
        input_path.as_posix(),
        "-o",
        output_path_no_ext.as_posix(),
        "-p",
        "-v",
        "-g",
        "--box_enveloping",
        input_path.as_posix(),
        "--padding",
        str(padding),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0 and expected_output.exists():
        log(f"✓ Meeko conversion successful")
        if verbose and result.stdout:
            log(result.stdout.strip())
        box_result: Optional[str] = None
        box_path = expected_output.with_suffix(".box.txt")
        box_params = _calculate_box_from_pdb(input_path, padding=padding)
        if box_params:
            try:
                _write_box_file(box_path, *box_params)
                box_result = box_path.as_posix()
                log(f"✓ Box file generated: {box_path}")
            except Exception as exc:
                log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
        return expected_output.as_posix(), box_result
    
    # ATTEMPT 2: Retry with extended options (box_enveloping, allow_bad_res)
    log(f"⚠ Meeko attempt 1 failed, retrying with extended options...")
    if verbose and result.stderr:
        log(f"Error: {result.stderr.strip()[:200]}")
    
    base_name_retry1 = f"{base_name}_retry1"
    output_path_no_ext_retry1 = dest / base_name_retry1
    expected_output_retry1 = dest / f"{base_name_retry1}.pdbqt"
    
    cmd_retry1 = [
        "mk_prepare_receptor.py",
        "-i",
        input_path.as_posix(),
        "-o",
        output_path_no_ext_retry1.as_posix(),
        "-p",
        "-v",
        "-g",
        "--box_enveloping",
        input_path.as_posix(),
        "--padding",
        str(padding),
        "--allow_bad_res",
    ]
    
    result_retry1 = subprocess.run(cmd_retry1, capture_output=True, text=True)
    
    if result_retry1.returncode == 0 and expected_output_retry1.exists():
        log(f"✓ Meeko retry 1 successful (extended options)")
        if verbose and result_retry1.stdout:
            log(result_retry1.stdout.strip())
        box_result: Optional[str] = None
        box_path = expected_output_retry1.with_suffix(".box.txt")
        box_params = _calculate_box_from_pdb(input_path, padding=padding)
        if box_params:
            try:
                _write_box_file(box_path, *box_params)
                box_result = box_path.as_posix()
                log(f"✓ Box file generated: {box_path}")
            except Exception as exc:
                log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
        return expected_output_retry1.as_posix(), box_result
    
    # ATTEMPT 3: Remove problematic residues and retry
    log(f"⚠ Meeko attempt 2 failed, attempting to remove problematic residues...")
    if verbose and result_retry1.stderr:
        log(f"Error: {result_retry1.stderr.strip()[:200]}")
    
    problematic_residues = _extract_problematic_residues(result_retry1.stderr)
    if not problematic_residues:
        problematic_residues = _extract_problematic_residues(result.stderr)
    
    if problematic_residues:
        log(f"Identified problematic residues: {problematic_residues}")
        
        # Create a cleaned PDB without problematic residues
        cleaned_pdb = dest / f"{input_path.stem}_cleaned.pdb"
        if _remove_residues_from_pdb(
            input_path.as_posix(),
            problematic_residues,
            cleaned_pdb.as_posix(),
            verbose=verbose,
        ):
            # Retry with cleaned PDB
            base_name_retry2 = f"{base_name}_retry2"
            output_path_no_ext_retry2 = dest / base_name_retry2
            expected_output_retry2 = dest / f"{base_name_retry2}.pdbqt"
            
            cmd_retry2 = [
                "mk_prepare_receptor.py",
                "-i",
                cleaned_pdb.as_posix(),
                "-o",
                output_path_no_ext_retry2.as_posix(),
                "-p",
                "-v",
                "-g",
                "--box_enveloping",
                cleaned_pdb.as_posix(),
                "--padding",
                str(padding),
            ]
            
            result_retry2 = subprocess.run(cmd_retry2, capture_output=True, text=True)
            
            if result_retry2.returncode == 0 and expected_output_retry2.exists():
                log(f"✓ Meeko retry 2 successful (after removing problematic residues)")
                if verbose and result_retry2.stdout:
                    log(result_retry2.stdout.strip())
                box_result: Optional[str] = None
                box_path = expected_output_retry2.with_suffix(".box.txt")
                box_params = _calculate_box_from_pdb(cleaned_pdb, padding=padding)
                if box_params:
                    try:
                        _write_box_file(box_path, *box_params)
                        box_result = box_path.as_posix()
                        log(f"✓ Box file generated: {box_path}")
                    except Exception as exc:
                        log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
                return expected_output_retry2.as_posix(), box_result
    
    # All attempts failed
    error_msg = f"Meeko receptor preparation failed for {protein_pdb} after all retry attempts"
    log(f"✗ {error_msg}")
    raise RuntimeError(error_msg)


# Meeko assigns extended atom types (CG0, CG1, CG2, G0, G1, G2, Si, B) that
# are not recognized by standard AutoDock Vina / Vina-CUDA.  Map them back
# to the closest AD4 type so that docking does not fail.
_MEEKO_AD4_TYPE_MAP: Dict[str, str] = {
    "CG0": "C",
    "CG1": "C",
    "CG2": "C",
    "G0":  "C",
    "G1":  "C",
    "G2":  "C",
    "Si":  "S",
    "B":   "C",
}


def sanitize_pdbqt_atom_types(
    pdbqt_path: Path,
    type_map: Optional[Dict[str, str]] = None,
    verbose: bool = True,
) -> int:
    """Remap non-standard atom types in a PDBQT file in-place.

    Returns the number of atoms that were remapped.
    """
    log = _printer(verbose)
    if type_map is None:
        type_map = _MEEKO_AD4_TYPE_MAP

    text = pdbqt_path.read_text()
    lines = text.splitlines(keepends=True)
    n_fixed = 0

    for i, line in enumerate(lines):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        # AD4 atom type starts at column 78 (0-indexed 77) to end of line
        if len(line.rstrip()) < 78:
            continue
        atype = line[77:].strip()
        if atype in type_map:
            new_type = type_map[atype]
            # Right-pad to preserve fixed-width format
            padded = new_type.ljust(len(line.rstrip()) - 77)
            lines[i] = line[:77] + padded + "\n"
            n_fixed += 1

    if n_fixed:
        pdbqt_path.write_text("".join(lines))
        log(f"  Sanitized {n_fixed} non-standard atom type(s) in {pdbqt_path.name}")

    return n_fixed


def convert_ligand_with_meeko(
    ligand_path: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    verbose: bool = True,
    add_tool_postfix: bool = False,
    use_converter_prefix: bool = False,
) -> str:
    """Convert a ligand structure (PDB/SDF/MOL2) to PDBQT using Meeko."""

    log = _printer(verbose)
    input_path = Path(ligand_path)
    target_extension = ".pdbqt" if input_path.suffix.lower() != ".pdbqt" else input_path.suffix
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_LIGAND_PDBQT",
        apply_process_suffix=process_postfixes,
        extension=target_extension,
        tool_postfix="_mko" if add_tool_postfix else "",
        tool_prefix="mko_" if use_converter_prefix else "",
    )

    log(f"Running Meeko ligand prep for {input_path.name} → {output_path.name}")
    cmd = [
        "mk_prepare_ligand.py",
        "-i",
        input_path.as_posix(),
        "-o",
        output_path.as_posix(),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"Meeko ligand preparation failed for {ligand_path}: {exc}"
        raise RuntimeError(message) from exc

    # Remap non-standard Meeko atom types so Vina-CUDA can parse the file
    sanitize_pdbqt_atom_types(output_path, verbose=verbose)

    return output_path.as_posix()


def convert_protein_with_mgltools(
    protein_pdb: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    padding: float = 2.0,
    verbose: bool = True,
    add_tool_postfix: bool = False,
    use_converter_prefix: bool = False,
) -> Tuple[str, Optional[str]]:
    """Convert a protein PDB to PDBQT using MGL Tools (prepare_receptor4.py)."""

    log = _printer(verbose)
    input_path = Path(protein_pdb)
    target_extension = ".pdbqt" if input_path.suffix.lower() != ".pdbqt" else input_path.suffix
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_PDBQT",
        apply_process_suffix=process_postfixes,
        extension=target_extension,
        tool_postfix="_mgl" if add_tool_postfix else "",
        tool_prefix="mgl_" if use_converter_prefix else "",
    )

    log(f"Running MGL Tools receptor prep for {input_path.name} → {output_path.name}")
    
    # Use ADFRsuite prepare_receptor binary (works better than calling Python scripts directly)
    prepare_receptor_bin = Path(MGLTOOLS_PATH) / "bin" / "prepare_receptor" if MGLTOOLS_PATH else None
    
    if prepare_receptor_bin and prepare_receptor_bin.exists():
        # Use the ADFRsuite binary directly (recommended)
        cmd = [
            str(prepare_receptor_bin),
            "-r",
            input_path.as_posix(),
            "-o",
            output_path.as_posix(),
            "-A",
            "hydrogens",
            "-U",
            "nphs_lps_waters_nonstdres",
        ]
    else:
        # Fallback to Python script (requires MGLTools Python environment)
        cmd = [
            "prepare_receptor4.py",
            "-r",
            input_path.as_posix(),
            "-o",
            output_path.as_posix(),
            "-A",
            "hydrogens",
            "-U",
            "nphs_lps_waters_nonstdres",
        ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"MGL Tools receptor preparation failed for {protein_pdb}: {exc}"
        raise RuntimeError(message) from exc

    box_result: Optional[str] = None
    box_path = output_path.with_suffix(".box.txt")
    box_params = _calculate_box_from_pdb(input_path, padding=padding)
    if box_params:
        try:
            _write_box_file(box_path, *box_params)
            box_result = box_path.as_posix()
            log(f"✓ Box file generated: {box_path}")
        except Exception as exc:
            log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
    else:
        log(f"⚠ Could not determine box parameters for {input_path.name}; box file skipped")

    return output_path.as_posix(), box_result


def convert_ligand_with_mgltools(
    ligand_path: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    verbose: bool = True,
    add_tool_postfix: bool = False,
    use_converter_prefix: bool = False,
) -> str:
    """Convert a ligand structure (PDB/SDF/MOL2) to PDBQT using MGL Tools (prepare_ligand4.py)."""

    log = _printer(verbose)
    input_path = Path(ligand_path)
    target_extension = ".pdbqt" if input_path.suffix.lower() != ".pdbqt" else input_path.suffix
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_LIGAND_PDBQT",
        apply_process_suffix=process_postfixes,
        extension=target_extension,
        tool_postfix="_mgl" if add_tool_postfix else "",
        tool_prefix="mgl_" if use_converter_prefix else "",
    )

    # MGL Tools prepare_ligand4.py requires MOL2 or PDB input
    # If input is SDF, convert to MOL2 first using OpenBabel
    intermediate_path = input_path
    if input_path.suffix.lower() == ".sdf":
        mol2_path = input_path.with_suffix(".mol2")
        if not mol2_path.exists():
            log(f"Converting SDF to MOL2 for MGL Tools compatibility...")
            try:
                subprocess.run(
                    ["obabel", input_path.as_posix(), "-O", mol2_path.as_posix()],
                    capture_output=True,
                    text=True,
                    check=True,
                )
            except (subprocess.CalledProcessError, FileNotFoundError) as exc:
                message = f"OpenBabel SDF→MOL2 conversion failed: {exc}"
                raise RuntimeError(message) from exc
        intermediate_path = mol2_path

    log(f"Running MGL Tools ligand prep for {intermediate_path.name} → {output_path.name}")
    
    # Use ADFRsuite prepare_ligand binary (works better than calling Python scripts directly)
    prepare_ligand_bin = Path(MGLTOOLS_PATH) / "bin" / "prepare_ligand" if MGLTOOLS_PATH else None
    
    if prepare_ligand_bin and prepare_ligand_bin.exists():
        # Use the ADFRsuite binary directly (recommended)
        cmd = [
            str(prepare_ligand_bin),
            "-l",
            intermediate_path.as_posix(),
            "-o",
            output_path.as_posix(),
            "-A",
            "hydrogens",
        ]
    else:
        # Fallback to Python script (requires MGLTools Python environment)
        cmd = [
            "prepare_ligand4.py",
            "-l",
            intermediate_path.as_posix(),
            "-o",
            output_path.as_posix(),
            "-A",
            "hydrogens",
        ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"MGL Tools ligand preparation failed for {ligand_path}: {exc}"
        raise RuntimeError(message) from exc

    return output_path.as_posix()


def convert_protein_with_pymol(
    protein_pdb: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    padding: float = 2.0,
    add_hydrogens: bool = True,
    verbose: bool = True,
    add_tool_postfix: bool = False,
) -> Tuple[str, Optional[str]]:
    """Convert a protein PDB using PyMOL for cleaning and preparation."""

    log = _printer(verbose)
    input_path = Path(protein_pdb)
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_PYMOL",
        apply_process_suffix=process_postfixes,
        extension=".pdb",
        tool_postfix="_pym" if add_tool_postfix else "",
    )

    log(f"Running PyMOL protein prep for {input_path.name} → {output_path.name}")

    # Build PyMOL script for protein preparation
    pymol_commands = [
        f"load {input_path.as_posix()}, protein",
        "remove solvent",
        "remove resn HOH",
        "remove hydrogens" if not add_hydrogens else "",
        "h_add" if add_hydrogens else "",
        f"save {output_path.as_posix()}, protein",
        "quit",
    ]
    pymol_script = "; ".join(cmd for cmd in pymol_commands if cmd)

    cmd = [
        "pymol",
        "-cq",
        "-d",
        pymol_script,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"PyMOL protein preparation failed for {protein_pdb}: {exc}"
        raise RuntimeError(message) from exc

    box_result: Optional[str] = None
    box_path = output_path.with_suffix(".box.txt")
    box_params = _calculate_box_from_pdb(input_path, padding=padding)
    if box_params:
        try:
            _write_box_file(box_path, *box_params)
            box_result = box_path.as_posix()
            log(f"✓ Box file generated: {box_path}")
        except Exception as exc:
            log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
    else:
        log(f"⚠ Could not determine box parameters for {input_path.name}; box file skipped")

    return output_path.as_posix(), box_result


def convert_ligand_with_pymol(
    ligand_path: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    output_format: str = "mol2",
    add_hydrogens: bool = True,
    verbose: bool = True,
    add_tool_postfix: bool = False,
) -> str:
    """Convert a ligand structure using PyMOL for cleaning and format conversion."""

    log = _printer(verbose)
    input_path = Path(ligand_path)
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_PYMOL",
        apply_process_suffix=process_postfixes,
        extension=output_format,
        tool_postfix="_pym" if add_tool_postfix else "",
    )

    log(f"Running PyMOL ligand prep for {input_path.name} → {output_path.name}")

    # Build PyMOL script for ligand preparation
    pymol_commands = [
        f"load {input_path.as_posix()}, ligand",
        "h_add" if add_hydrogens else "",
        f"save {output_path.as_posix()}, ligand",
        "quit",
    ]
    pymol_script = "; ".join(cmd for cmd in pymol_commands if cmd)

    cmd = [
        "pymol",
        "-cq",
        "-d",
        pymol_script,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"PyMOL ligand preparation failed for {ligand_path}: {exc}"
        raise RuntimeError(message) from exc

    return output_path.as_posix()


def convert_protein_with_openbabel(
    protein_pdb: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    padding: float = 2.0,
    add_hydrogens: bool = True,
    verbose: bool = True,
    add_tool_postfix: bool = False,
) -> Tuple[str, Optional[str]]:
    """Convert a protein PDB to PDBQT using OpenBabel."""

    log = _printer(verbose)
    input_path = Path(protein_pdb)
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_PDBQT",
        apply_process_suffix=process_postfixes,
        extension=".pdbqt",
        tool_postfix="_obabel" if add_tool_postfix else "",
    )

    log(f"Running OpenBabel protein conversion for {input_path.name} → {output_path.name}")
    cmd = [
        "obabel",
        input_path.as_posix(),
        "-O",
        output_path.as_posix(),
        "-xh",  # output polar hydrogens only (PDBQT-specific)
        "-xr",  # rigid receptor mode for PDBQT
    ]
    if add_hydrogens:
        cmd.append("-h")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"OpenBabel protein conversion failed for {protein_pdb}: {exc}"
        raise RuntimeError(message) from exc

    box_result: Optional[str] = None
    box_path = output_path.with_suffix(".box.txt")
    box_params = _calculate_box_from_pdb(input_path, padding=padding)
    if box_params:
        try:
            _write_box_file(box_path, *box_params)
            box_result = box_path.as_posix()
            log(f"✓ Box file generated: {box_path}")
        except Exception as exc:
            log(f"⚠ Failed to write box file for {input_path.name}: {exc}")
    else:
        log(f"⚠ Could not determine box parameters for {input_path.name}; box file skipped")

    return output_path.as_posix(), box_result


def convert_ligand_with_openbabel(
    ligand_path: str,
    *,
    output_dir: Optional[Path],
    custom_postfix: str,
    process_postfixes: bool,
    verbose: bool = True,
    add_tool_postfix: bool = False,
) -> str:
    """Convert a ligand structure to PDBQT using OpenBabel."""

    log = _printer(verbose)
    input_path = Path(ligand_path)
    output_path = _build_output_path(
        input_path,
        output_dir=output_dir,
        base_postfix=custom_postfix,
        default_suffix="_LIGAND_PDBQT",
        apply_process_suffix=process_postfixes,
        extension=".pdbqt",
        tool_postfix="_obabel" if add_tool_postfix else "",
    )

    log(f"Running OpenBabel ligand conversion for {input_path.name} → {output_path.name}")
    cmd = [
        "obabel",
        input_path.as_posix(),
        "-O",
        output_path.as_posix(),
        # "-xh",  # output polar hydrogens only (PDBQT-specific)
        "-xr",  # rigid receptor mode for PDBQT
        "-h",   # add hydrogens
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stdout:
            log(result.stdout.strip())
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"OpenBabel ligand conversion failed for {ligand_path}: {exc}"
        raise RuntimeError(message) from exc

    return output_path.as_posix()


def convert_with_openbabel_generic(
    input_path: str,
    output_extension: str,
    *,
    output_dir: Optional[Path] = None,
    custom_postfix: str = "",
    add_hydrogens: bool = False,
    verbose: bool = True,
) -> str:
    """Generic OpenBabel conversion between any supported formats."""
    
    log = _printer(verbose)
    source = Path(input_path)
    
    dest = output_dir if output_dir else source.parent
    dest.mkdir(parents=True, exist_ok=True)
    
    ext = output_extension if output_extension.startswith(".") else f".{output_extension}"
    stem = _strip_extension_from_stem(source.stem, ext)
    output_path = dest / f"{stem}{custom_postfix}{ext}"
    
    log(f"OpenBabel: {source.name} → {output_path.name}")
    
    cmd = ["obabel", source.as_posix(), "-O", output_path.as_posix()]
    if add_hydrogens:
        cmd.append("-h")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"OpenBabel conversion failed for {input_path}: {exc}"
        raise RuntimeError(message) from exc
    
    return output_path.as_posix()


def convert_with_pymol_generic(
    input_path: str,
    output_extension: str,
    *,
    output_dir: Optional[Path] = None,
    custom_postfix: str = "",
    add_hydrogens: bool = False,
    verbose: bool = True,
) -> str:
    """Generic PyMOL conversion between formats."""
    
    log = _printer(verbose)
    source = Path(input_path)
    
    dest = output_dir if output_dir else source.parent
    dest.mkdir(parents=True, exist_ok=True)
    
    ext = output_extension if output_extension.startswith(".") else f".{output_extension}"
    stem = _strip_extension_from_stem(source.stem, ext)
    output_path = dest / f"{stem}{custom_postfix}{ext}"
    
    log(f"PyMOL: {source.name} → {output_path.name}")
    
    pymol_commands = [
        f"load {source.as_posix()}, mol",
        "h_add" if add_hydrogens else "",
        f"save {output_path.as_posix()}, mol",
        "quit",
    ]
    pymol_script = "; ".join(cmd for cmd in pymol_commands if cmd)
    
    cmd = ["pymol", "-cq", "-d", pymol_script]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if verbose and result.stderr:
            log(result.stderr.strip())
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        message = f"PyMOL conversion failed for {input_path}: {exc}"
        raise RuntimeError(message) from exc
    
    return output_path.as_posix()

# Global variable for MGLTools path - can be set by users
MGLTOOLS_PATH: Optional[Path] = "/opt/ADFRsuite-1.0"

def set_mgltools_path(path: str) -> None:
    """Set the path to MGLTools installation directory."""
    global MGLTOOLS_PATH
    MGLTOOLS_PATH = Path(path)
    if not MGLTOOLS_PATH.exists():
        raise FileNotFoundError(f"MGLTools path not found: {MGLTOOLS_PATH}")
    print(f"✓ MGLTools path set to: {MGLTOOLS_PATH}")


def get_mgltools_python() -> str:
    """Get the path to MGLTools Python executable."""
    if MGLTOOLS_PATH is None:
        # Try common installation paths
        common_paths = [
            Path("/usr/local/MGLTools"),
            Path("/opt/mgltools"),
            Path.home() / "MGLTools",
            Path.home() / "mgltools",
            Path("/home/manndo/MGLTools"),
        ]
        for path in common_paths:
            if path.exists():
                set_mgltools_path(str(path))
                break
    
    if MGLTOOLS_PATH is None:
        raise RuntimeError(
            "MGLTools path not set. Use set_mgltools_path('/path/to/MGLTools') "
            "or ensure MGLTools is installed in a standard location."
        )
    
    # Look for pythonsh or python executable
    pythonsh = MGLTOOLS_PATH / "bin" / "pythonsh"
    if pythonsh.exists():
        return str(pythonsh)
    
    python_exe = MGLTOOLS_PATH / "bin" / "python"
    if python_exe.exists():
        return str(python_exe)
    
    # Try MGLToolsPckgs path structure
    pythonsh_alt = MGLTOOLS_PATH / "MGLToolsPckgs" / "bin" / "pythonsh"
    if pythonsh_alt.exists():
        return str(pythonsh_alt)
    
    raise FileNotFoundError(f"Could not find Python executable in MGLTools: {MGLTOOLS_PATH}")


def get_mgltools_script(script_name: str) -> str:
    """Get the path to an MGLTools script."""
    if MGLTOOLS_PATH is None:
        get_mgltools_python()  # This will set MGLTOOLS_PATH or raise error
    
    # Common script locations
    script_paths = [
        MGLTOOLS_PATH / "MGLToolsPckgs" / "AutoDockTools" / "Utilities24" / script_name,
        MGLTOOLS_PATH / "Utilities24" / script_name,
        MGLTOOLS_PATH / "bin" / script_name,
        MGLTOOLS_PATH / "MGLToolsPckgs" / "bin" / script_name,
    ]
    
    for script_path in script_paths:
        if script_path.exists():
            return str(script_path)
    
    raise FileNotFoundError(f"Could not find MGLTools script: {script_name}")


def convert_with_mgltools_generic(
    input_path: str,
    output_extension: str,
    *,
    output_dir: Optional[Path] = None,
    custom_postfix: str = "",
    add_hydrogens: bool = True,
    verbose: bool = True,
    mgltools_path: Optional[str] = None,
) -> str:
    """
    Generic MGLTools conversion between formats.
    
    MGLTools primarily handles PDB/MOL2 to PDBQT conversions.
    For other format conversions, it uses OpenBabel as a fallback.
    
    Parameters:
    -----------
    input_path : str
        Path to input file
    output_extension : str
        Target format extension (pdb, mol2, pdbqt, sdf)
    output_dir : Path
        Output directory
    custom_postfix : str
        Custom suffix for output filename
    add_hydrogens : bool
        Add hydrogens during conversion
    verbose : bool
        Print progress messages
    mgltools_path : str
        Path to MGLTools installation (optional, uses global if not set)
    """
    
    log = _printer(verbose)
    source = Path(input_path)
    
    # Set MGLTools path if provided
    if mgltools_path:
        set_mgltools_path(mgltools_path)
    
    dest = output_dir if output_dir else source.parent
    dest.mkdir(parents=True, exist_ok=True)
    
    ext = output_extension if output_extension.startswith(".") else f".{output_extension}"
    stem = _strip_extension_from_stem(source.stem, ext)
    output_path = dest / f"{stem}{custom_postfix}{ext}"
    
    log(f"MGLTools: {source.name} → {output_path.name}")
    
    input_ext = source.suffix.lower()
    target_ext = ext.lower()
    
    # MGLTools is best for PDBQT conversions
    if target_ext == ".pdbqt":
        try:
            pythonsh = get_mgltools_python()
            
            if input_ext in [".pdb", ".mol2"]:
                # Use prepare_ligand4.py or prepare_receptor4.py
                # Determine if it's likely a ligand or receptor based on file size
                file_size = source.stat().st_size
                if file_size < 50000:  # Likely a ligand (small file)
                    script = get_mgltools_script("prepare_ligand4.py")
                    cmd = [pythonsh, script, "-l", source.as_posix(), "-o", output_path.as_posix()]
                else:  # Likely a receptor (large file)
                    script = get_mgltools_script("prepare_receptor4.py")
                    cmd = [pythonsh, script, "-r", source.as_posix(), "-o", output_path.as_posix()]
                
                if add_hydrogens:
                    cmd.extend(["-A", "hydrogens"])
                
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                if verbose and result.stdout:
                    log(result.stdout.strip())
                return output_path.as_posix()
            
        except (FileNotFoundError, RuntimeError) as exc:
            log(f"⚠ MGLTools PDBQT conversion failed, falling back to OpenBabel: {exc}")
    
    # For non-PDBQT conversions or if MGLTools fails, use OpenBabel
    # MGLTools doesn't natively support SDF format well
    return convert_with_openbabel_generic(
        input_path,
        output_extension,
        output_dir=output_dir,
        custom_postfix=custom_postfix,
        add_hydrogens=add_hydrogens,
        verbose=verbose,
    )


def run_conversion_chain(
    input_file: str,
    converter: str,
    output_base_dir: Path,
    *,
    chain: Optional[List[str]] = None,
    verbose: bool = True,
    mgltools_path: Optional[str] = None,
) -> Dict[str, str]:
    """
    Run a chain of format conversions using the specified converter tool.
    
    Parameters:
    -----------
    input_file : str
        Path to input file (PDB, XYZ, etc.)
    converter : str
        Conversion tool to use: 'openbabel', 'pymol', 'mgltools', or 'meeko'
    output_base_dir : Path
        Base directory for storing converted files
    chain : List[str]
        List of format extensions for conversion chain (e.g., ['sdf', 'mol2', 'pdb'])
    verbose : bool
        Print progress messages
    mgltools_path : str
        Path to MGLTools installation directory (required for mgltools converter)
    
    Returns:
    --------
    Dict[str, str]
        Dictionary mapping step names to output file paths
    """
    
    log = _printer(verbose)
    
    if chain is None:
        chain = ["sdf", "mol2", "pdb"]
    
    # Set MGLTools path if provided
    if mgltools_path and converter == "mgltools":
        set_mgltools_path(mgltools_path)
    
    source = Path(input_file)
    tool_dir = output_base_dir / converter
    tool_dir.mkdir(parents=True, exist_ok=True)
    
    results = {"input": input_file}
    current_file = input_file
    
    log(f"\n{'='*50}")
    log(f"Conversion chain using {converter.upper()}")
    log(f"Input: {source.name}")
    log(f"Chain: {source.suffix} → " + " → ".join(f".{ext}" for ext in chain))
    log(f"{'='*50}")
    
    for i, target_ext in enumerate(chain):
        step_name = f"step_{i+1}_{target_ext}"
        step_dir = tool_dir / f"step_{i+1}_{target_ext}"
        step_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            if converter == "openbabel":
                output_file = convert_with_openbabel_generic(
                    current_file,
                    target_ext,
                    output_dir=step_dir,
                    verbose=verbose,
                )
            elif converter == "pymol":
                output_file = convert_with_pymol_generic(
                    current_file,
                    target_ext,
                    output_dir=step_dir,
                    verbose=verbose,
                )
            elif converter == "mgltools":
                # Use MGLTools generic conversion (falls back to OpenBabel for unsupported formats)
                output_file = convert_with_mgltools_generic(
                    current_file,
                    target_ext,
                    output_dir=step_dir,
                    verbose=verbose,
                )
            elif converter == "meeko":
                # Meeko is specifically for PDBQT, use OpenBabel for intermediate steps
                output_file = convert_with_openbabel_generic(
                    current_file,
                    target_ext,
                    output_dir=step_dir,
                    verbose=verbose,
                )
            else:
                raise ValueError(f"Unknown converter: {converter}")
            
            results[step_name] = output_file
            current_file = output_file
            log(f"✓ Step {i+1}: Created {Path(output_file).name}")
            
        except RuntimeError as exc:
            log(f"✗ Step {i+1} failed: {exc}")
            results[step_name] = None
            break
    
    results["final_output"] = current_file
    return results


def compare_pdb_files(
    original_pdb: str,
    converted_pdb: str,
    *,
    verbose: bool = True,
) -> Dict[str, object]:
    """
    Compare two PDB files and report differences.
    
    Parameters:
    -----------
    original_pdb : str
        Path to original PDB file
    converted_pdb : str
        Path to converted PDB file
    verbose : bool
        Print detailed comparison
    
    Returns:
    --------
    Dict with comparison results
    """
    
    log = _printer(verbose)
    
    def parse_pdb_atoms(pdb_path: str) -> List[Dict]:
        """Parse ATOM/HETATM records from PDB file."""
        atoms = []
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        try:
                            atoms.append({
                                'record': line[:6].strip(),
                                'serial': int(line[6:11].strip()) if line[6:11].strip() else 0,
                                'name': line[12:16].strip(),
                                'resname': line[17:20].strip(),
                                'chain': line[21].strip() if len(line) > 21 else '',
                                'resnum': int(line[22:26].strip()) if line[22:26].strip() else 0,
                                'x': float(line[30:38].strip()) if line[30:38].strip() else 0.0,
                                'y': float(line[38:46].strip()) if line[38:46].strip() else 0.0,
                                'z': float(line[46:54].strip()) if line[46:54].strip() else 0.0,
                                'element': line[76:78].strip() if len(line) > 76 else '',
                            })
                        except (ValueError, IndexError):
                            continue
        except FileNotFoundError:
            pass
        return atoms
    
    orig_atoms = parse_pdb_atoms(original_pdb)
    conv_atoms = parse_pdb_atoms(converted_pdb)
    
    # Calculate statistics
    orig_count = len(orig_atoms)
    conv_count = len(conv_atoms)
    atom_diff = conv_count - orig_count
    
    # Get unique elements
    orig_elements = set(a['element'] for a in orig_atoms if a['element'])
    conv_elements = set(a['element'] for a in conv_atoms if a['element'])
    
    # Get unique residues
    orig_residues = set((a['resname'], a['resnum'], a['chain']) for a in orig_atoms)
    conv_residues = set((a['resname'], a['resnum'], a['chain']) for a in conv_atoms)
    
    # Calculate coordinate differences for matching atoms
    coord_diffs = []
    orig_by_name = {(a['name'], a['resnum'], a['chain']): a for a in orig_atoms}
    for conv_a in conv_atoms:
        key = (conv_a['name'], conv_a['resnum'], conv_a['chain'])
        if key in orig_by_name:
            orig_a = orig_by_name[key]
            diff = (
                (conv_a['x'] - orig_a['x'])**2 +
                (conv_a['y'] - orig_a['y'])**2 +
                (conv_a['z'] - orig_a['z'])**2
            ) ** 0.5
            coord_diffs.append(diff)
    
    avg_coord_diff = sum(coord_diffs) / len(coord_diffs) if coord_diffs else 0.0
    max_coord_diff = max(coord_diffs) if coord_diffs else 0.0
    
    comparison = {
        'original_file': original_pdb,
        'converted_file': converted_pdb,
        'original_atom_count': orig_count,
        'converted_atom_count': conv_count,
        'atom_difference': atom_diff,
        'original_elements': sorted(orig_elements),
        'converted_elements': sorted(conv_elements),
        'elements_added': sorted(conv_elements - orig_elements),
        'elements_removed': sorted(orig_elements - conv_elements),
        'original_residues': len(orig_residues),
        'converted_residues': len(conv_residues),
        'residues_added': len(conv_residues - orig_residues),
        'residues_removed': len(orig_residues - conv_residues),
        'avg_coordinate_difference': avg_coord_diff,
        'max_coordinate_difference': max_coord_diff,
        'matching_atoms': len(coord_diffs),
    }
    
    if verbose:
        log(f"\n{'='*60}")
        log(f"PDB COMPARISON: {Path(original_pdb).name} vs {Path(converted_pdb).name}")
        log(f"{'='*60}")
        log(f"Original atoms:  {orig_count:>6}")
        log(f"Converted atoms: {conv_count:>6}")
        log(f"Difference:      {atom_diff:>+6}")
        log(f"")
        log(f"Original elements:  {', '.join(comparison['original_elements'])}")
        log(f"Converted elements: {', '.join(comparison['converted_elements'])}")
        if comparison['elements_added']:
            log(f"Elements added:     {', '.join(comparison['elements_added'])}")
        if comparison['elements_removed']:
            log(f"Elements removed:   {', '.join(comparison['elements_removed'])}")
        log(f"")
        log(f"Original residues:  {comparison['original_residues']}")
        log(f"Converted residues: {comparison['converted_residues']}")
        log(f"")
        log(f"Matching atoms:           {comparison['matching_atoms']}")
        log(f"Avg coordinate diff (Å):  {avg_coord_diff:.4f}")
        log(f"Max coordinate diff (Å):  {max_coord_diff:.4f}")
    
    return comparison


def run_conversion_workflow(
    input_files: List[str],
    converters: List[str],
    output_base_dir: Path,
    *,
    chain: Optional[List[str]] = None,
    verbose: bool = True,
    mgltools_path: Optional[str] = None,
    max_workers: int = DEFAULT_MAX_WORKERS,
) -> Dict[str, Dict]:
    """
    Run conversion chains for multiple files using multiple converters in parallel.
    
    Parameters:
    -----------
    input_files : List[str]
        List of input file paths
    converters : List[str]
        List of converter tools to use
    output_base_dir : Path
        Base directory for output
    chain : List[str]
        Conversion chain (default: ['sdf', 'mol2', 'pdb'])
    verbose : bool
        Print progress
    mgltools_path : str
        Path to MGLTools installation directory (required for mgltools converter)
    max_workers : int
        Maximum number of parallel threads
    
    Returns:
    --------
    Dict with all conversion results and comparisons
    """
    
    log = _printer(verbose)
    
    if chain is None:
        chain = ["sdf", "mol2", "pdb"]
    
    # Set MGLTools path if provided and mgltools is in converters
    if mgltools_path and "mgltools" in converters:
        set_mgltools_path(mgltools_path)
    
    output_base_dir = Path(output_base_dir)
    output_base_dir.mkdir(parents=True, exist_ok=True)
    
    all_results = {}
    comparisons = {}
    
    log("\n" + "=" * 70)
    log("MULTI-TOOL CONVERSION WORKFLOW")
    log("=" * 70)
    log(f"Input files: {len(input_files)}")
    log(f"Converters: {', '.join(converters)}")
    log(f"Chain: " + " → ".join(chain))
    if mgltools_path:
        log(f"MGLTools path: {mgltools_path}")
    
    # Build all (file, converter) tasks for parallel execution
    tasks: List[Tuple[str, str]] = []
    for input_file in input_files:
        for conv in converters:
            tasks.append((input_file, conv))

    def _run_one_chain(task):
        input_file, conv = task
        file_name = Path(input_file).stem
        try:
            result = run_conversion_chain(
                input_file,
                conv,
                output_base_dir / file_name,
                chain=chain,
                verbose=verbose,
                mgltools_path=mgltools_path,
            )
            comparison = None
            final_output = result.get("final_output")
            if (final_output and
                final_output.endswith('.pdb') and
                input_file.endswith('.pdb') and
                Path(final_output).exists()):
                comparison = compare_pdb_files(
                    input_file,
                    final_output,
                    verbose=verbose,
                )
            return (file_name, conv), (result, comparison)
        except Exception as exc:
            log(f"✗ {conv} failed for {Path(input_file).name}: {exc}")
            return (file_name, conv), ({"error": str(exc)}, None)

    log(f"⚡ Running {len(tasks)} conversion chains in parallel...")
    workers = min(max_workers, len(tasks)) if tasks else 1
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(_run_one_chain, t): t for t in tasks}
        for future in as_completed(futures):
            (file_name, conv), (result, comparison) = future.result()
            all_results.setdefault(file_name, {})[conv] = result
            if comparison:
                comparisons.setdefault(file_name, {})[conv] = comparison
    
    return {
        "results": all_results,
        "comparisons": comparisons,
    }


class XYZConverter:
    """Convert XYZ files to PDB, SDF, and MOL2 formats"""
    
    def __init__(
        self,
        verbose: bool = True,
        output_dir: Optional[Path] = None,
        custom_postfix: str = "",
    ):
        self.verbose = verbose
        self.converted_files = {}
        self.output_dir = Path(output_dir) if output_dir else None
        self.custom_postfix = custom_postfix
        self._check_openbabel()
    
    def _check_openbabel(self):
        """Check if OpenBabel is installed"""
        try:
            subprocess.run(['obabel', '-V'], capture_output=True, check=True)
            if self.verbose:
                print("✓ OpenBabel is installed")
        except (FileNotFoundError, subprocess.CalledProcessError) as exc:
            message = (
                "OpenBabel not found. Install with:\n"
                "  conda install -c conda-forge openbabel\n"
                "  OR pip install openbabel-wheel"
            )
            raise RuntimeError(message) from exc
    
    def convert_xyz(self, xyz_file: str, formats: List[str] = None) -> Dict[str, str]:
        """
        Convert XYZ file to specified formats
        
        Parameters:
        -----------
        xyz_file : str
            Path to input XYZ file
        formats : List[str]
            Output formats (pdb, sdf, mol2)
        
        Returns:
        --------
        Dict[str, str]
            Dictionary mapping format to output file path
        """
        
        if formats is None:
            formats = ['pdb', 'sdf', 'mol2']
        
        if not os.path.exists(xyz_file):
            print(f"✗ File not found: {xyz_file}")
            return {}
        
        source_path = Path(xyz_file)
        base = source_path.stem
        results = {}
        
        if self.verbose:
            print(f"\n{'='*60}")
            print(f"Converting: {xyz_file}")
            print(f"{'='*60}")
        
        for fmt in formats:
            output_path = _build_output_path(
                source_path,
                output_dir=self.output_dir,
                base_postfix=self.custom_postfix,
                default_suffix="",
                apply_process_suffix=False,
                extension=fmt,
            )

            try:
                result = subprocess.run(
                    ["obabel", xyz_file, "-O", output_path.as_posix(), "-h"],
                    capture_output=True,
                    text=True,
                    timeout=30,
                )
                
                if result.returncode == 0 and output_path.exists():
                    size = output_path.stat().st_size
                    if self.verbose:
                        print(f"✓ {fmt.upper():4} ({size:8} bytes) → {output_path}")
                    results[fmt] = output_path.as_posix()
                else:
                    if self.verbose:
                        print(f"✗ {fmt.upper()} conversion failed")
                        if result.stderr:
                            print(f"  Error: {result.stderr[:100]}")
                    
            except subprocess.TimeoutExpired:
                print(f"✗ Conversion timeout for {fmt.upper()}")
            except Exception as e:
                print(f"✗ Error converting to {fmt.upper()}: {e}")
        
        self.converted_files[xyz_file] = results
        return results
    
    def batch_convert(
        self,
        xyz_files: List[str],
        formats: List[str] = None,
        max_workers: int = DEFAULT_MAX_WORKERS,
    ) -> Dict:
        """Convert multiple XYZ files in parallel."""
        print(f"\nPROCESSING {len(xyz_files)} LIGAND(S)")
        print("="*60)

        existing = [f for f in xyz_files if os.path.exists(f)]
        for f in xyz_files:
            if not os.path.exists(f):
                print(f"⚠ Skipping: {f} (not found)")

        if not existing:
            return {}

        def _convert_one(xyz_file):
            results = self.convert_xyz(xyz_file, formats)
            return xyz_file, results

        all_results = _run_parallel(
            _convert_one,
            existing,
            max_workers=max_workers,
            desc="ligand XYZ files",
            verbose=self.verbose,
        )
        # Filter out any exceptions
        return {k: v for k, v in all_results.items() if not isinstance(v, Exception)}


class PDBValidator:
    """Validate and fix PDB files using PDBFixer"""
    
    def __init__(
        self,
        verbose: bool = True,
        repair_terminals: bool = False,
        output_dir: Optional[Path] = None,
        custom_postfix: str = "",
        process_postfixes: bool = True,
        require_pdbfixer: bool = False,
    ):
        self.verbose = verbose
        self.fixed_files = {}
        self.repair_terminals = repair_terminals
        self.output_dir = Path(output_dir) if output_dir else None
        self.custom_postfix = custom_postfix
        self.process_postfixes = process_postfixes
        self.pdbfixer_available = self._check_pdbfixer(require=require_pdbfixer)
    
    def _check_pdbfixer(self, require: bool = False):
        """Check if PDBFixer is installed"""
        try:
            from pdbfixer import PDBFixer
            if self.verbose:
                print("✓ PDBFixer is installed")
            return True
        except ImportError as exc:
            if require:
                message = (
                    "PDBFixer not found. Install with:\n"
                    "  pip install pdbfixer\n"
                    "  OR conda install -c conda-forge pdbfixer"
                )
                raise RuntimeError(message) from exc
            else:
                if self.verbose:
                    print("⚠ PDBFixer not available - will skip validation")
                return False
    
    def validate_pdb(
        self,
        pdb_file: str,
        *,
        add_hydrogens: bool = False,
        remove_heterogens: bool = True,
        repair_terminals: Optional[bool] = None,
    ) -> Optional[str]:
        """
        Validate and fix PDB file
        
        Parameters:
        -----------
        pdb_file : str
            Path to input PDB file
        add_hydrogens : bool
            Add missing hydrogens
        remove_heterogens : bool
            Remove water and heterogens
        
        Returns:
        --------
        str : Path to fixed PDB file, or None if error
        """
        
        if not os.path.exists(pdb_file):
            print(f"✗ File not found: {pdb_file}")
            return None
        
        # If PDBFixer is not available, return original file
        if not self.pdbfixer_available:
            if self.verbose:
                print(f"⚠ Skipping validation for {pdb_file} (PDBFixer not available)")
            return pdb_file
        
        try:
            from pdbfixer import PDBFixer
            from openmm.app import PDBFile
            
            if self.verbose:
                print(f"\n{'='*60}")
                print(f"Validating: {pdb_file}")
                print(f"{'='*60}")
            
            fixer = PDBFixer(filename=pdb_file)
            
            issues_found = []
            
            # Check for missing residues
            fixer.findMissingResidues()
            if fixer.missingResidues:
                fixer.addMissingAtoms()
                issues_found.append(f"Missing residues: {len(fixer.missingResidues)}")
                if self.verbose:
                    print(f"• Fixed {len(fixer.missingResidues)} missing residues")
            else:
                if self.verbose:
                    print("✓ No missing residues")
            
            # Check for non-standard residues
            fixer.findNonstandardResidues()
            if fixer.nonstandardResidues:
                fixer.replaceNonstandardResidues()
                issues_found.append(f"Non-standard residues: {len(fixer.nonstandardResidues)}")
                if self.verbose:
                    print(f"• Replaced {len(fixer.nonstandardResidues)} non-standard residues")
            else:
                if self.verbose:
                    print("✓ All residues are standard")
            
            # Remove heterogens and water
            if remove_heterogens:
                fixer.removeHeterogens(keepWater=False)
                if self.verbose:
                    print("• Removed heterogens (including waters)")
            
            # Add hydrogens if requested
            if add_hydrogens:
                fixer.addMissingHydrogens(7.4)
                if self.verbose:
                    print("• Added missing hydrogens (pH 7.4)")
            
            # Save fixed structure
            base_path = Path(pdb_file)
            output_path = _build_output_path(
                base_path,
                output_dir=self.output_dir,
                base_postfix=self.custom_postfix,
                default_suffix="_FIXED",
                apply_process_suffix=self.process_postfixes,
                extension=base_path.suffix,
            )

            with open(output_path, 'w') as f:
                PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
            if self.verbose:
                print(f"\n✓ Fixed PDB saved: {output_path}")

            if repair_terminals is None:
                repair_terminals = self.repair_terminals

            final_output = output_path.as_posix()
            if repair_terminals:
                repaired_path = _build_output_path(
                    base_path,
                    output_dir=self.output_dir,
                    base_postfix=self.custom_postfix,
                    default_suffix="_REPAIRED",
                    apply_process_suffix=self.process_postfixes,
                    extension=base_path.suffix,
                )
                final_output = fix_pdb_terminal_residues(
                    output_path,
                    output_pdb=repaired_path.as_posix(),
                    verbose=self.verbose,
                    output_dir=self.output_dir,
                    custom_postfix=self.custom_postfix,
                    process_postfixes=self.process_postfixes,
                )
            
            self.fixed_files[pdb_file] = final_output
            return final_output
            
        except Exception as e:
            print(f"✗ Error validating {pdb_file}:")
            print(f"  {str(e)[:200]}")
            return None
    
    def batch_validate(
        self,
        pdb_files: List[str],
        *,
        add_hydrogens: bool = False,
        repair_terminals: Optional[bool] = None,
        max_workers: int = DEFAULT_MAX_WORKERS,
    ) -> Dict[str, Optional[str]]:
        """Validate multiple PDB files in parallel."""
        print(f"\nPROCESSING {len(pdb_files)} PROTEIN STRUCTURE(S)")
        print("="*60)

        all_results: Dict[str, Optional[str]] = {}
        existing = [f for f in pdb_files if os.path.exists(f)]
        for f in pdb_files:
            if not os.path.exists(f):
                print(f"⚠ Skipping: {f} (not found)")
                all_results[f] = None

        if not existing:
            return all_results

        def _validate_one(pdb_file):
            result = self.validate_pdb(
                pdb_file,
                add_hydrogens=add_hydrogens,
                repair_terminals=repair_terminals,
            )
            return pdb_file, result

        parallel_results = _run_parallel(
            _validate_one,
            existing,
            max_workers=max_workers,
            desc="protein PDB files",
            verbose=self.verbose,
        )
        for k, v in parallel_results.items():
            all_results[k] = None if isinstance(v, Exception) else v
        return all_results


class DockingPreparation:
    """Prepare files for specific docking software"""
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
    
    def prepare_for_vina(self, protein_pdb: str, ligand_sdf: str) -> Dict:
        """
        Prepare files for AutoDock Vina
        Requires conversion to PDBQT format
        """
        
        print(f"\n{'='*60}")
        print("AUTODOCK VINA PREPARATION")
        print(f"{'='*60}")
        
        print(f"\nProtein: {protein_pdb}")
        print(f"Ligand:  {ligand_sdf}")
        
        print("\nRequirements:")
        print("  • Convert both to PDBQT format")
        print("  • Define docking box coordinates")
        print("  • Specify search parameters")

        print("\nConversion Commands:")
        print("\n1. Using Meeko (Recommended):")
        print("   pip install meeko rdkit")
        print(f"   mk_prepare_ligand.py -i {ligand_sdf} -o ligand.pdbqt")
        print(f"   mk_prepare_receptor.py -r {protein_pdb} -o protein.pdbqt")

        print("\n2. Using AutoDock Tools (Legacy):")
        print(f"   prepare_receptor -r {protein_pdb} -o protein.pdbqt")
        print(f"   obabel {ligand_sdf} -xpdbqt -o ligand.pdbqt")
        
        return {
            'software': 'AutoDock Vina',
            'protein': protein_pdb,
            'ligand': ligand_sdf,
            'required_formats': ['PDBQT', 'PDBQT']
        }
    
    def prepare_for_equibind(self, protein_pdb: str, ligand_sdf: str) -> Dict:
        """
        Prepare files for Equibind
        Accepts PDB and SDF directly
        """
        
        print(f"\n{'='*60}")
        print("EQUIBIND PREPARATION")
        print(f"{'='*60}")
        
        print(f"\nProtein: {protein_pdb}")
        print(f"Ligand:  {ligand_sdf}")
        
        print("\nStatus: ✓ Ready for Equibind")
        print("\nNotes:")
        print("  • Equibind accepts PDB format directly")
        print("  • Ligand should be in SDF with 3D coordinates")
        print("  • Ensure protein is cleaned (use PDBFixer)")
        print("  • Equibind can handle standard PDB variations")
        
        return {
            'software': 'Equibind',
            'protein': protein_pdb,
            'ligand': ligand_sdf,
            'required_formats': ['PDB', 'SDF'],
            'ready': True
        }
    
    def prepare_for_diffdock(self, protein_pdb: str, ligand_file: str) -> Dict:
        """
        Prepare files for DiffDock
        Accepts multiple ligand formats
        """
        
        print(f"\n{'='*60}")
        print("DIFFDOCK PREPARATION")
        print(f"{'='*60}")
        
        print(f"\nProtein: {protein_pdb}")
        print(f"Ligand:  {ligand_file}")
        
        # Determine ligand format
        if ligand_file.endswith('.sdf'):
            fmt = 'SDF'
        elif ligand_file.endswith('.mol2'):
            fmt = 'MOL2'
        else:
            fmt = 'Unknown'
        
        print(f"\nLigand Format: {fmt}")
        
        print("\nStatus: ✓ Ready for DiffDock")
        print("\nNotes:")
        print("  • DiffDock accepts PDB format directly")
        print("  • Ligand can be SDF, MOL2, or SMILES")
        print("  • Ensure 3D coordinates for structure files")
        print("  • SMILES input will auto-generate 3D conformers")
        print("  • Clean protein with PDBFixer for best results")
        
        return {
            'software': 'DiffDock',
            'protein': protein_pdb,
            'ligand': ligand_file,
            'ligand_format': fmt,
            'ready': True
        }
    
    def generate_summary(self, protein_file: str, 
                        converted_ligands: Dict, 
                        docking_tools: List[str]) -> str:
        """Generate preparation summary"""
        
        summary = f"\n{'='*60}\n"
        summary += "DOCKING PREPARATION SUMMARY\n"
        summary += f"{'='*60}\n\n"
        
        summary += f"PROTEIN: {protein_file}\n"
        summary += f"LIGANDS: {len(converted_ligands)} ligand(s)\n\n"
        
        summary += "CONVERSION STATUS:\n"
        for ligand, formats in converted_ligands.items():
            for fmt, file in formats.items():
                if file:
                    size = os.path.getsize(file) if os.path.exists(file) else 0
                    summary += f"  ✓ {ligand} → {fmt.upper()} ({size:,} bytes)\n"
        
        summary += f"\nREADY FOR DOCKING:\n"
        for tool in docking_tools:
            summary += f"  ✓ {tool}\n"
        
        summary += f"\n{'='*60}\n"
        
        return summary


def resolve_input_files(
    *,
    ligand_files: Optional[Iterable[str]] = None,
    protein_files: Optional[Iterable[str]] = None,
    input_dir: Optional[Path] = None,
    contains: Optional[str] = None,
    recursive: bool = False,
) -> Tuple[List[str], List[str], List[str]]:
    """Normalize user inputs into ligand and protein file lists.

    Returns (xyz_files, pdb_files, ready_ligands) where *ready_ligands*
    are SDF/MOL2/PDB files that do NOT need OpenBabel XYZ conversion.
    """

    _XYZ_EXTS = {".xyz"}
    _READY_EXTS = {".sdf", ".mol2", ".pdb"}

    xyz_files: List[str] = []
    ready_ligands: List[str] = []
    pdb_files = list(protein_files or [])

    def _classify_ligands(paths: Iterable[str]) -> None:
        for p in paths:
            ext = Path(p).suffix.lower()
            if ext in _XYZ_EXTS:
                xyz_files.append(p)
            elif ext in _READY_EXTS:
                ready_ligands.append(p)
            else:
                xyz_files.append(p)  # fallback: try OpenBabel

    if ligand_files:
        _classify_ligands(ligand_files)

    if input_dir:
        if not contains:
            raise ValueError("'contains' must be provided when specifying input_dir")
        discovered = discover_input_files(input_dir, contains, recursive)
        if contains.lower() == "ligands":
            _classify_ligands(discovered)
        else:
            pdb_files.extend(discovered)

    if not xyz_files and not ready_ligands and not pdb_files:
        xyz_files = list(DEFAULT_LIGAND_FILES)
        pdb_files = list(DEFAULT_PROTEIN_FILES)

    return xyz_files, pdb_files, ready_ligands


# Valid converter options
CONVERTER_OPTIONS = ("openbabel", "meeko", "mgltools", "pymol")


def run_workflow(
    *,
    ligand_files: Optional[Iterable[str]] = None,
    protein_files: Optional[Iterable[str]] = None,
    input_dir: Optional[Path] = None,
    contains: Optional[str] = None,
    recursive: bool = False,
    ligand_formats: Optional[List[str]] = None,
    add_hydrogens: bool = False,
    verbose: bool = True,
    repair_terminals: bool = False,
    output_dir: Optional[Path] = None,
    custom_postfix: str = "",
    process_postfixes: bool = True,
    convert_with_meeko: bool = False,
    convert_ligands_with_meeko: bool = False,
    converter: Optional[str] = None,
    convert_proteins: bool = False,
    convert_ligands: bool = False,
    add_tool_postfix: bool = False,
    use_converter_prefix: bool = False,
    log_file: Optional[Path] = None,
    skip_pdb_validation: bool = False,
    max_workers: int = DEFAULT_MAX_WORKERS,
) -> Dict[str, object]:
    """Execute the preparation workflow and return collected outputs.

    Parameters
    ----------
    max_workers : int
        Maximum number of parallel threads for conversion tasks.
    """

    log = _printer(verbose)

    log_entries: List[Dict[str, object]] = []

    def add_entry(tool: str, stage: str, input_file: str, output_file: Optional[str], status: str, message: str = ""):
        log_entries.append(
            {
                "tool": tool,
                "stage": stage,
                "input": input_file,
                "output": output_file,
                "status": status,
                "message": message,
            }
        )

    log("\n" + "=" * 60)
    log("ORAI1 PROTEIN-LIGAND DOCKING PREPARATION WORKFLOW")
    log("=" * 60)

    output_dir = Path(output_dir) if output_dir else None

    xyz_files, pdb_files, ready_ligands = resolve_input_files(
        ligand_files=ligand_files,
        protein_files=protein_files,
        input_dir=input_dir,
        contains=contains,
        recursive=recursive,
    )

    log("\nSTEP 1: CONVERTING LIGAND FILES (XYZ → PDB/SDF/MOL2)")
    converted_ligands: Dict[str, Dict[str, str]] = {}
    if xyz_files:
        xyz_converter = XYZConverter(
            verbose=verbose,
            output_dir=output_dir,
            custom_postfix=custom_postfix,
        )
        xyz_converter.batch_convert(list(xyz_files), formats=ligand_formats, max_workers=max_workers)
        converted_ligands = xyz_converter.converted_files
        # Log ligand conversions
        for original_xyz, fmt_map in converted_ligands.items():
            for fmt, out_path in fmt_map.items():
                status = "success" if out_path else "failed"
                message = "" if out_path else "No output produced"
                add_entry("openbabel", f"ligand_convert_{fmt}", original_xyz, out_path, status, message)

    # SDF/MOL2/PDB ligands are already in a dockable format — register them
    # directly so downstream steps (Meeko, etc.) can pick them up without an
    # intermediate OpenBabel round-trip that may corrupt query features.
    if ready_ligands:
        log(f"• {len(ready_ligands)} ligand(s) already in SDF/MOL2/PDB format — using as-is")
        import shutil
        for lig_path in ready_ligands:
            lig = Path(lig_path)
            ext = lig.suffix.lower().lstrip(".")
            if output_dir:
                dest = Path(output_dir) / lig.name
                if not dest.exists():
                    shutil.copy2(lig, dest)
                converted_ligands[lig_path] = {ext: dest.as_posix()}
            else:
                converted_ligands[lig_path] = {ext: lig_path}
            add_entry("passthrough", f"ligand_ready_{ext}", lig_path,
                      converted_ligands[lig_path][ext], "success", "Already in target format")

    if not xyz_files and not ready_ligands:
        log("• No ligand files supplied; skipping conversion")

    log("\n\nSTEP 2: VALIDATING PROTEIN STRUCTURES")
    protein_results: Dict[str, Optional[str]] = {}
    if pdb_files:
        if skip_pdb_validation:
            log("• PDB validation skipped (using original files)")
            # Use original PDB files without validation
            protein_results = {pdb: pdb for pdb in pdb_files}
            for pdb_in in pdb_files:
                add_entry("pdbfixer", "protein_validate", pdb_in, pdb_in, "skipped", "Validation skipped by user")
        else:
            validator = PDBValidator(
                verbose=verbose,
                repair_terminals=repair_terminals,
                output_dir=output_dir,
                custom_postfix=custom_postfix,
                process_postfixes=process_postfixes,
            )
            protein_results = validator.batch_validate(
                list(pdb_files),
                add_hydrogens=add_hydrogens,
                repair_terminals=repair_terminals,
                max_workers=max_workers,
            )
            for pdb_in, pdb_out in protein_results.items():
                status = "success" if pdb_out else "failed"
                message = "" if pdb_out else "Validation/repair failed"
                add_entry("pdbfixer", "protein_validate", pdb_in, pdb_out, status, message)
    else:
        log("• No protein files supplied; skipping validation")

    # Validate converter option
    if converter and converter not in CONVERTER_OPTIONS:
        raise ValueError(
            f"Invalid converter '{converter}'. Choose from: {', '.join(CONVERTER_OPTIONS)}"
        )

    log("\n\nSTEP 3: PREPARING FOR MOLECULAR DOCKING")
    prep = DockingPreparation(verbose=verbose)
    fixed_proteins = [v for v in protein_results.values() if v]

    selected_protein = None
    selected_ligand = None
    selected_ligand_formats: Optional[Dict[str, str]] = None
    docking_tools: List[str] = []

    meeko_outputs: Dict[str, Optional[str]] = {}
    meeko_box_files: Dict[str, Optional[str]] = {}
    meeko_ligand_outputs: Dict[str, Optional[str]] = {}

    if fixed_proteins and converted_ligands:
        selected_protein = fixed_proteins[0]
        ligand_sdf = None
        for ligand, fmt_map in converted_ligands.items():
            if 'sdf' in fmt_map:
                selected_ligand = ligand
                ligand_sdf = fmt_map['sdf']
                selected_ligand_formats = fmt_map
                break

        if ligand_sdf:
            prep.prepare_for_equibind(selected_protein, ligand_sdf)
            docking_tools.append('Equibind')
            prep.prepare_for_diffdock(selected_protein, ligand_sdf)
            docking_tools.append('DiffDock')

            if selected_ligand_formats and 'mol2' in selected_ligand_formats:
                prep.prepare_for_vina(selected_protein, ligand_sdf)
                docking_tools.append('AutoDock Vina')

    if convert_with_meeko and fixed_proteins:
        def _meeko_protein(protein):
            pdbqt_path, box_path = convert_protein_with_meeko(
                protein,
                output_dir=output_dir,
                custom_postfix=custom_postfix,
                process_postfixes=process_postfixes,
                verbose=verbose,
                add_tool_postfix=add_tool_postfix,
                use_converter_prefix=use_converter_prefix,
            )
            return protein, (pdbqt_path, box_path)

        log("⚡ Parallelizing Meeko protein conversions...")
        workers = min(max_workers, len(fixed_proteins))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_meeko_protein, p): p for p in fixed_proteins}
            for future in as_completed(futures):
                protein = futures[future]
                try:
                    _, (pdbqt_path, box_path) = future.result()
                    meeko_outputs[protein] = pdbqt_path
                    meeko_box_files[protein] = box_path
                    log(f"✓ Meeko PDBQT generated: {pdbqt_path}")
                    add_entry("meeko", "protein_pdbqt", protein, pdbqt_path, "success")
                except RuntimeError as exc:
                    meeko_outputs[protein] = None
                    meeko_box_files[protein] = None
                    log(f"✗ Meeko conversion failed for {protein}: {exc}")
                    add_entry("meeko", "protein_pdbqt", protein, None, "failed", str(exc))

    if convert_ligands_with_meeko and converted_ligands:
        preferred_formats = ("sdf", "mol2", "pdb")
        # Resolve source files first
        ligand_tasks: List[Tuple[str, str]] = []  # (original_key, source_path)
        for original_xyz, format_map in converted_ligands.items():
            ligand_source = None
            for fmt in preferred_formats:
                candidate = format_map.get(fmt)
                if candidate:
                    ligand_source = candidate
                    break
            if not ligand_source:
                meeko_ligand_outputs[original_xyz] = None
                log(f"⚠ No suitable format found for Meeko ligand prep: {original_xyz}")
            else:
                ligand_tasks.append((original_xyz, ligand_source))

        def _meeko_ligand(task):
            original_xyz, ligand_source = task
            ligand_pdbqt = convert_ligand_with_meeko(
                ligand_source,
                output_dir=output_dir,
                custom_postfix=custom_postfix,
                process_postfixes=process_postfixes,
                verbose=verbose,
                add_tool_postfix=add_tool_postfix,
                use_converter_prefix=use_converter_prefix,
            )
            return original_xyz, (ligand_source, ligand_pdbqt)

        if ligand_tasks:
            log("⚡ Parallelizing Meeko ligand conversions...")
            workers = min(max_workers, len(ligand_tasks))
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {executor.submit(_meeko_ligand, t): t for t in ligand_tasks}
                for future in as_completed(futures):
                    original_xyz, ligand_source = futures[future]
                    try:
                        _, (src, ligand_pdbqt) = future.result()
                        meeko_ligand_outputs[original_xyz] = ligand_pdbqt
                        log(f"✓ Meeko ligand PDBQT generated: {ligand_pdbqt}")
                        add_entry("meeko", "ligand_pdbqt", src, ligand_pdbqt, "success")
                    except RuntimeError as exc:
                        meeko_ligand_outputs[original_xyz] = None
                        log(f"✗ Meeko ligand conversion failed for {ligand_source}: {exc}")
                        add_entry("meeko", "ligand_pdbqt", ligand_source, None, "failed", str(exc))

    # Converter-based protein conversion
    converter_protein_outputs: Dict[str, Optional[str]] = {}
    converter_box_files: Dict[str, Optional[str]] = {}
    if converter and convert_proteins and fixed_proteins:
        log(f"\n\nSTEP 4: CONVERTING PROTEINS USING {converter.upper()}")

        def _convert_protein(protein):
            if converter == "openbabel":
                return protein, convert_protein_with_openbabel(
                    protein, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, add_hydrogens=add_hydrogens,
                    verbose=verbose, add_tool_postfix=add_tool_postfix,
                )
            elif converter == "meeko":
                return protein, convert_protein_with_meeko(
                    protein, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix, use_converter_prefix=use_converter_prefix,
                )
            elif converter == "mgltools":
                return protein, convert_protein_with_mgltools(
                    protein, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix, use_converter_prefix=use_converter_prefix,
                )
            elif converter == "pymol":
                return protein, convert_protein_with_pymol(
                    protein, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, add_hydrogens=add_hydrogens,
                    verbose=verbose, add_tool_postfix=add_tool_postfix,
                )
            else:
                raise ValueError(f"Unknown converter: {converter}")

        log(f"⚡ Parallelizing {converter.upper()} protein conversions...")
        workers = min(max_workers, len(fixed_proteins))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_convert_protein, p): p for p in fixed_proteins}
            for future in as_completed(futures):
                protein = futures[future]
                try:
                    _, (pdbqt_path, box_path) = future.result()
                    converter_protein_outputs[protein] = pdbqt_path
                    converter_box_files[protein] = box_path
                    log(f"✓ {converter.upper()} protein output: {pdbqt_path}")
                    add_entry(converter, "protein_convert", protein, pdbqt_path, "success")
                except RuntimeError as exc:
                    converter_protein_outputs[protein] = None
                    converter_box_files[protein] = None
                    log(f"✗ {converter.upper()} protein conversion failed for {protein}: {exc}")
                    add_entry(converter, "protein_convert", protein, None, "failed", str(exc))

    # Converter-based ligand conversion
    converter_ligand_outputs: Dict[str, Optional[str]] = {}
    if converter and convert_ligands and converted_ligands:
        log(f"\n\nSTEP 5: CONVERTING LIGANDS USING {converter.upper()}")
        preferred_formats = ("sdf", "mol2", "pdb")
        # Resolve source files first
        ligand_conv_tasks: List[Tuple[str, str]] = []
        for original_xyz, format_map in converted_ligands.items():
            ligand_source = None
            for fmt in preferred_formats:
                candidate = format_map.get(fmt)
                if candidate:
                    ligand_source = candidate
                    break
            if not ligand_source:
                converter_ligand_outputs[original_xyz] = None
                log(f"⚠ No suitable format found for {converter.upper()} ligand prep: {original_xyz}")
            else:
                ligand_conv_tasks.append((original_xyz, ligand_source))

        def _convert_ligand(task):
            original_xyz, ligand_source = task
            if converter == "openbabel":
                out = convert_ligand_with_openbabel(
                    ligand_source, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix,
                )
            elif converter == "meeko":
                out = convert_ligand_with_meeko(
                    ligand_source, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix, use_converter_prefix=use_converter_prefix,
                )
            elif converter == "mgltools":
                out = convert_ligand_with_mgltools(
                    ligand_source, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix, use_converter_prefix=use_converter_prefix,
                )
            elif converter == "pymol":
                out = convert_ligand_with_pymol(
                    ligand_source, output_dir=output_dir, custom_postfix=custom_postfix,
                    process_postfixes=process_postfixes, verbose=verbose,
                    add_tool_postfix=add_tool_postfix,
                )
            else:
                raise ValueError(f"Unknown converter: {converter}")
            return original_xyz, (ligand_source, out)

        if ligand_conv_tasks:
            log(f"⚡ Parallelizing {converter.upper()} ligand conversions...")
            workers = min(max_workers, len(ligand_conv_tasks))
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {executor.submit(_convert_ligand, t): t for t in ligand_conv_tasks}
                for future in as_completed(futures):
                    original_xyz, ligand_source = futures[future]
                    try:
                        _, (src, ligand_output) = future.result()
                        converter_ligand_outputs[original_xyz] = ligand_output
                        log(f"✓ {converter.upper()} ligand output: {ligand_output}")
                        add_entry(converter, "ligand_convert", src, ligand_output, "success")
                    except RuntimeError as exc:
                        converter_ligand_outputs[original_xyz] = None
                        log(f"✗ {converter.upper()} ligand conversion failed for {ligand_source}: {exc}")
                        add_entry(converter, "ligand_convert", ligand_source, None, "failed", str(exc))

    log("\n\n" + "=" * 60)
    log("WORKFLOW COMPLETE")
    log("=" * 60)

    log("\nNext Steps:")
    if converted_ligands:
        log("1. Review converted ligand files in the output directory")
    if protein_results:
        log("2. Use fixed PDB files for downstream docking")
    log("3. For AutoDock Vina, convert to PDBQT using Meeko when ligands and proteins are available")
    log("4. For Equibind & Diffdock, use cleaned PDB/SDF files directly")
    log("\nRefer to docking_guide.md for detailed instructions")

    summary = None
    if selected_protein and docking_tools:
        summary = prep.generate_summary(selected_protein, converted_ligands, docking_tools)

    overview_table = None
    if log_entries:
        headers = ["Tool", "Stage", "Input", "Output", "Status", "Message"]
        rows = []
        for entry in log_entries:
            rows.append([
                entry.get("tool", ""),
                entry.get("stage", ""),
                Path(entry.get("input", "")).name if entry.get("input") else "",
                Path(entry.get("output", "")).name if entry.get("output") else "",
                entry.get("status", ""),
                entry.get("message", ""),
            ])

        # compute column widths
        widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(str(cell)))

        def fmt_row(parts):
            return " | ".join(str(part).ljust(widths[i]) for i, part in enumerate(parts))

        sep = "-+-".join("-" * w for w in widths)
        lines = [fmt_row(headers), sep]
        lines.extend(fmt_row(r) for r in rows)
        overview_table = "\n".join(lines)

        log("\nOVERVIEW TABLE")
        log(overview_table)

        if log_file:
            log_path = Path(log_file)
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_path.write_text(overview_table)
            log(f"\nLog written to {log_path}")

    return {
        "ligand_files": xyz_files,
        "protein_files": pdb_files,
        "converted_ligands": converted_ligands,
        "protein_results": protein_results,
        "fixed_proteins": fixed_proteins,
        "selected_protein": selected_protein,
        "selected_ligand": selected_ligand,
        "selected_ligand_formats": selected_ligand_formats,
        "docking_tools": docking_tools,
        "summary": summary,
        "meeko_outputs": meeko_outputs,
        "meeko_box_files": meeko_box_files,
        "meeko_ligand_outputs": meeko_ligand_outputs,
        "converter": converter,
        "converter_protein_outputs": converter_protein_outputs,
        "converter_box_files": converter_box_files,
        "converter_ligand_outputs": converter_ligand_outputs,
        "log_entries": log_entries,
        "overview_table": overview_table,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Automate ligand conversion and protein validation for docking workflows."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        help="Path to a directory containing ligand or protein files to process",
    )
    parser.add_argument(
        "--contains",
        choices=["ligands", "proteins"],
        help="Specify whether --input-dir stores ligand (.xyz) or protein (.pdb) files",
    )
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recursively search subdirectories within --input-dir",
    )
    parser.add_argument(
        "--ligand-formats",
        nargs="+",
        help="Override ligand conversion formats (default: pdb sdf mol2)",
    )
    parser.add_argument(
        "--add-hydrogens",
        action="store_true",
        help="Add missing hydrogens during protein validation",
    )
    parser.add_argument(
        "--repair-terminals",
        action="store_true",
        help="Run terminal residue repair routine on validated proteins",
    )
    parser.add_argument(
        "--skip-pdb-validation",
        action="store_true",
        help="Skip PDB validation with PDBFixer (use original protein files)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Destination folder for converted ligands and validated proteins",
    )
    parser.add_argument(
        "--postfix",
        default="",
        help="Custom postfix appended to generated files (e.g., _CLEAN)",
    )
    parser.add_argument(
        "--process-postfixes",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Include automatic process postfixes such as _FIXED/_REPAIRED",
    )
    parser.add_argument(
        "--meeko",
        action="store_true",
        help="Run mk_prepare_receptor.py on validated proteins to produce PDBQT files",
    )
    parser.add_argument(
        "--meeko-ligands",
        action="store_true",
        help="Run mk_prepare_ligand.py on converted ligands to produce PDBQT files",
    )
    parser.add_argument(
        "--converter",
        choices=CONVERTER_OPTIONS,
        help=f"Conversion tool to use: {', '.join(CONVERTER_OPTIONS)}",
    )
    parser.add_argument(
        "--convert-proteins",
        action="store_true",
        help="Convert proteins using the specified --converter tool",
    )
    parser.add_argument(
        "--convert-ligands",
        action="store_true",
        help="Convert ligands using the specified --converter tool",
    )
    parser.add_argument(
        "--add-tool-postfix",
        action="store_true",
        help="Add tool-specific postfix to converted filenames (e.g., _mgl, _mko, _pym, _obabel)",
    )
    parser.add_argument(
        "--use-converter-prefix",
        action="store_true",
        help="Add tool-specific prefix to converted filenames (e.g., mgl_file.pdbqt, mko_file.pdbqt)",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        help="Path to write conversion overview log file",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=DEFAULT_MAX_WORKERS,
        help=f"Maximum number of parallel threads for conversion tasks (default: {DEFAULT_MAX_WORKERS})",
    )
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.input_dir and not args.contains:
        parser.error("--contains is required when --input-dir is provided")

    if args.converter and not (args.convert_proteins or args.convert_ligands):
        parser.error("--converter requires --convert-proteins and/or --convert-ligands")

    try:
        run_workflow(
            input_dir=args.input_dir,
            contains=args.contains,
            recursive=args.recursive,
            ligand_formats=args.ligand_formats,
            add_hydrogens=args.add_hydrogens,
            verbose=True,
            repair_terminals=args.repair_terminals,
            output_dir=args.output_dir,
            custom_postfix=args.postfix,
            process_postfixes=args.process_postfixes,
            convert_with_meeko=args.meeko,
            convert_ligands_with_meeko=args.meeko_ligands,
            converter=args.converter,
            convert_proteins=args.convert_proteins,
            convert_ligands=args.convert_ligands,
            add_tool_postfix=args.add_tool_postfix,
            use_converter_prefix=args.use_converter_prefix,
            log_file=args.log_file,
            skip_pdb_validation=args.skip_pdb_validation,
            max_workers=args.max_workers,
        )
    except Exception as exc:  # pragma: no cover - CLI convenience
        print(f"✗ {exc}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
