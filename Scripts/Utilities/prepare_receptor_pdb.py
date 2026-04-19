#!/usr/bin/env python3
"""
Convert CHARMM-style PDB files to standard PDB format compatible with
AutoDock Vina receptor preparation tools (mk_prepare_receptor.py and
prepare_receptor4.py from MGLTools).

Handles:
1. Segment ID (MONA/MONB/...) → unique chain ID (A/B/...) mapping
2. CHARMM histidine names (HSP/HSD/HSE) → standard HIS
3. Removal of CHARMM capping groups (ACE-like CAY/CY/OY, CT3-like NT/CAT)
4. Hydrogen removal (recommended; receptor prep tools re-add them)
5. Element symbol population (columns 77-78)
6. Occupancy set to 1.00 (MGLTools silently drops atoms with occupancy < 1.0)
7. TER records between chains
8. Backbone-first atom ordering within residues
9. Alternate location cleanup (keep only ' ' or 'A')
"""

from pathlib import Path
from typing import Optional, List, Dict, Tuple
import re
import string


# ---------------------------------------------------------------------------
# CHARMM capping-group atom names
# ---------------------------------------------------------------------------
# N-terminal acetyl cap (ACE / CT1 patch): CAY, CY, OY + hydrogens
NTER_CAP_ATOMS = {'CAY', 'CY', 'OY', 'HY1', 'HY2', 'HY3'}
# C-terminal N-methylamide cap (CT3 patch): NT, CAT, HNT + hydrogens
CTER_CAP_ATOMS = {'NT', 'CAT', 'HNT', 'HT1', 'HT2', 'HT3'}
ALL_CAP_ATOMS = NTER_CAP_ATOMS | CTER_CAP_ATOMS

# Standard amino acid 3-letter codes (used to decide ATOM vs HETATM)
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
    'TYR', 'VAL',
    # Also accept protonation-state variants before renaming
    'HSP', 'HSD', 'HSE', 'HIP', 'HIE', 'HID',
}


def parse_pdb_line(line: str) -> Optional[dict]:
    """Parse a PDB ATOM/HETATM line into components (fixed-width columns)."""
    if not (line.startswith('ATOM') or line.startswith('HETATM')):
        return None

    return {
        'record': line[0:6].strip(),
        'serial': line[6:11].strip(),
        'name': line[12:16].strip(),
        'altLoc': line[16:17].strip(),
        'resName': line[17:20].strip(),
        'chainID': line[21:22].strip(),
        'resSeq': line[22:26].strip(),
        'iCode': line[26:27].strip(),
        'x': line[30:38].strip(),
        'y': line[38:46].strip(),
        'z': line[46:54].strip(),
        'occupancy': line[54:60].strip() if len(line) > 54 else '1.00',
        'tempFactor': line[60:66].strip() if len(line) > 60 else '0.00',
        'segment': line[72:76].strip() if len(line) > 72 else '',
        'element': line[76:78].strip() if len(line) > 76 else '',
        'charge': line[78:80].strip() if len(line) > 78 else '',
    }


def get_element_from_atom_name(atom_name: str, res_name: str = '') -> str:
    """Infer element symbol from atom name following PDB conventions.
    
    Uses residue context to disambiguate (e.g. CA in ALA = carbon-alpha,
    CA as a residue itself = calcium).
    """
    name = atom_name.strip()

    # Two-letter elements that could be ambiguous
    if res_name in STANDARD_RESIDUES:
        # In a standard amino acid context, these are always single-letter
        if name == 'CA':
            return 'C'  # carbon-alpha, not calcium
        if name == 'CD' or name.startswith('CD'):
            return 'C'  # delta carbon, not cadmium
    else:
        # Non-standard residue: check for known metal ion residue names
        if res_name in ('CA', 'CAL'):
            return 'CA'
        if res_name in ('ZN', 'ZN2'):
            return 'ZN'
        if res_name in ('MG', 'MG2'):
            return 'MG'
        if res_name in ('FE', 'FE2', 'FE3'):
            return 'FE'

    # Explicit two-letter elements (order matters: check before single-letter)
    if name[:2] in ('CL', 'BR', 'MG', 'ZN', 'FE', 'MN', 'CU', 'CO', 'NI', 'SE'):
        return name[:2].upper()

    # Hydrogen: any atom name starting with H, or digit+H (e.g. 1HB)
    if name[0] == 'H' or (len(name) > 1 and name[0].isdigit() and name[1] == 'H'):
        return 'H'

    # Standard single-letter elements
    first = name[0].upper()
    if first in ('C', 'N', 'O', 'S', 'P', 'F', 'I'):
        return first

    return first


def get_atom_order_priority(atom_name: str) -> int:
    """Return sort priority for atom ordering (lower = earlier).
    
    Standard PDB order: N, CA, C, O, then sidechain heavy atoms,
    then hydrogens (attached to each heavy atom in order).
    """
    name = atom_name.strip()

    backbone_order = {'N': 0, 'CA': 1, 'C': 2, 'O': 3, 'OXT': 4}
    if name in backbone_order:
        return backbone_order[name]

    # Side chain heavy atoms by Greek letter depth
    if not name.startswith('H') and not (len(name) > 1 and name[0].isdigit() and name[1] == 'H'):
        if name.startswith('C'):
            return 10
        elif name.startswith('N'):
            return 15
        elif name.startswith('O'):
            return 20
        elif name.startswith('S'):
            return 25
        else:
            return 50

    # Hydrogens come last, sub-sorted by parent atom
    return 100


def _build_segment_to_chain_map(segments: list) -> Dict[str, str]:
    """Map unique segment IDs to chain letters A-Z, a-z (up to 52 chains)."""
    chain_letters = list(string.ascii_uppercase + string.ascii_lowercase)
    mapping = {}
    for seg in segments:
        if seg not in mapping:
            if len(mapping) >= len(chain_letters):
                raise ValueError(f"Too many segments ({len(segments)}) to map to chain IDs (max 52)")
            mapping[seg] = chain_letters[len(mapping)]
    return mapping


def _format_atom_name(name: str, element: str) -> str:
    """Format atom name to PDB standard 4-character field (columns 13-16).
    
    PDB convention:
    - 1-char element names: name starts at column 14 (right-justified in 2-char element field)
      e.g. ' CA ', ' N  ', ' OG1'
    - 2-char element names: name starts at column 13
      e.g. 'FE  ', 'CL  ', 'BR1 '
    - 4-char atom names always start at column 13
      e.g. 'HD11', '1HG2'
    """
    if len(name) >= 4:
        return f"{name:<4}"
    if len(element) == 2:
        return f"{name:<4}"
    # 1-char element: pad left with one space
    return f" {name:<3}"


def clean_pdb_like_pymol(
    input_pdb: str,
    output_pdb: str,
    remove_hydrogens: bool = True,
    fix_histidines: bool = True,
    remove_caps: bool = True,
    assign_chains_from_segments: bool = True,
    set_occupancy: float = 1.00,
    keep_altloc: str = 'A',
    verbose: bool = True
) -> str:
    """
    Convert a CHARMM-style PDB to standard PDB format suitable for
    AutoDock Vina receptor preparation (mk_prepare_receptor.py / prepare_receptor4.py).

    Parameters
    ----------
    input_pdb : str
        Path to input PDB file.
    output_pdb : str
        Path to output PDB file.
    remove_hydrogens : bool
        Remove all hydrogen atoms (default True; receptor prep tools re-add them).
    fix_histidines : bool
        Convert HSP/HSE/HSD/HIP/HIE/HID → HIS.
    remove_caps : bool
        Remove CHARMM N-terminal (CAY/CY/OY) and C-terminal (NT/CAT/HNT) capping atoms.
    assign_chains_from_segments : bool
        Map CHARMM segment IDs (MONA, MONB, ...) to unique chain letters.
        Only applied when multiple segments share the same chain ID.
    set_occupancy : float
        Force occupancy to this value (MGLTools drops atoms with occupancy < 1.0).
    keep_altloc : str
        Which alternate location indicator to keep ('A' or ' '). Others are dropped.
    verbose : bool
        Print progress messages.

    Returns
    -------
    str
        Path to the output file.
    """
    input_path = Path(input_pdb)
    output_path = Path(output_pdb)

    if verbose:
        print(f"Converting PDB: {input_path.name} → {output_path.name}")

    with open(input_path, 'r') as f:
        lines = f.readlines()

    # --- First pass: collect all atoms and discover segments ---
    header_lines = []
    raw_atoms: List[dict] = []
    footer_lines = []
    in_atoms = False
    segment_order: List[str] = []  # preserve segment encounter order
    seen_segments: set = set()

    stats = {
        'total_input': 0,
        'dropped_hydrogen': 0,
        'dropped_cap': 0,
        'dropped_altloc': 0,
        'his_renamed': 0,
    }

    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            in_atoms = True
            atom = parse_pdb_line(line)
            if atom is None:
                continue
            stats['total_input'] += 1

            seg = atom['segment']
            if seg and seg not in seen_segments:
                segment_order.append(seg)
                seen_segments.add(seg)

            raw_atoms.append(atom)
        elif not in_atoms:
            header_lines.append(line)
        else:
            # Collect non-ATOM lines after atoms (CONECT, MASTER, END, etc.)
            if not line.startswith('END'):
                footer_lines.append(line)

    # --- Build segment → chain mapping if needed ---
    segment_chain_map: Dict[str, str] = {}
    if assign_chains_from_segments and len(segment_order) > 1:
        segment_chain_map = _build_segment_to_chain_map(segment_order)
        if verbose:
            print(f"  Segment → Chain mapping: {segment_chain_map}")
    elif assign_chains_from_segments and len(segment_order) == 1:
        # Single segment: keep original chain ID or assign 'A'
        seg = segment_order[0]
        segment_chain_map = {seg: 'A'}

    # --- Second pass: filter and transform atoms ---
    atoms_by_chain: Dict[str, List[dict]] = {}

    for atom in raw_atoms:
        name = atom['name']

        # 1. Remove capping group atoms
        if remove_caps and name in ALL_CAP_ATOMS:
            stats['dropped_cap'] += 1
            continue

        # 2. Remove hydrogens
        if remove_hydrogens:
            if name.startswith('H') or (len(name) > 1 and name[0].isdigit() and name[1] == 'H'):
                stats['dropped_hydrogen'] += 1
                continue

        # 3. Handle alternate locations: keep only ' ' or keep_altloc
        alt = atom['altLoc']
        if alt and alt != ' ' and alt != keep_altloc:
            stats['dropped_altloc'] += 1
            continue
        # Clear altloc indicator so output is clean
        atom['altLoc'] = ' '

        # 4. Fix histidine names
        if fix_histidines and atom['resName'] in ('HSP', 'HSE', 'HSD', 'HIP', 'HIE', 'HID'):
            stats['his_renamed'] += 1
            atom['resName'] = 'HIS'

        # 5. Assign chain ID from segment
        if segment_chain_map and atom['segment'] in segment_chain_map:
            atom['chainID'] = segment_chain_map[atom['segment']]
        elif not atom['chainID'].strip():
            atom['chainID'] = 'A'

        # 6. Populate element column
        atom['element'] = get_element_from_atom_name(name, atom['resName'])

        # 7. Force occupancy
        atom['occupancy'] = f"{set_occupancy:.2f}"

        # Group by chain for TER record insertion
        chain = atom['chainID']
        if chain not in atoms_by_chain:
            atoms_by_chain[chain] = []
        atoms_by_chain[chain].append(atom)

    # --- Sort atoms within each residue ---
    for chain in atoms_by_chain:
        # Group by residue within each chain
        residues: Dict[Tuple, List[dict]] = {}
        for atom in atoms_by_chain[chain]:
            res_key = (atom['resSeq'], atom['iCode'])
            if res_key not in residues:
                residues[res_key] = []
            residues[res_key].append(atom)

        # Sort within each residue
        for res_key in residues:
            residues[res_key].sort(key=lambda a: (get_atom_order_priority(a['name']), a['name']))

        # Flatten back in residue-number order
        sorted_atoms = []
        for res_key in sorted(residues.keys(), key=lambda k: (int(k[0]), k[1])):
            sorted_atoms.extend(residues[res_key])
        atoms_by_chain[chain] = sorted_atoms

    # --- Write output ---
    with open(output_path, 'w') as f:
        # Header (CRYST1, REMARK, etc.) – pad/truncate to exactly 80 chars per PDB v3.3
        for line in header_lines:
            content = line.rstrip('\n\r')[:80]
            f.write(f"{content:<80}\n")

        atom_serial = 1
        chain_order = sorted(atoms_by_chain.keys())

        for chain_idx, chain in enumerate(chain_order):
            for atom in atoms_by_chain[chain]:
                record = 'ATOM' if atom['resName'] in STANDARD_RESIDUES else 'HETATM'
                element = atom['element']
                formatted_name = _format_atom_name(atom['name'], element)

                line = (
                    f"{record:<6}"
                    f"{atom_serial:>5} "
                    f"{formatted_name}"
                    f"{atom['altLoc']:1}"
                    f"{atom['resName']:>3} "
                    f"{chain:1}"
                    f"{atom['resSeq']:>4}"
                    f"{atom['iCode']:1}   "
                    f"{float(atom['x']):>8.3f}"
                    f"{float(atom['y']):>8.3f}"
                    f"{float(atom['z']):>8.3f}"
                    f"{float(atom['occupancy']):>6.2f}"
                    f"{float(atom['tempFactor']):>6.2f}"
                    f"          "  # columns 67-76 (blank; no segment ID in standard PDB)
                    f"{element:>2}"
                    f"  "  # charge columns 79-80
                    "\n"
                )
                f.write(line)
                atom_serial += 1

            # TER record after each chain
            if atoms_by_chain[chain]:
                last = atoms_by_chain[chain][-1]
                ter = (
                    f"TER   {atom_serial:>5}      "
                    f"{last['resName']:>3} "
                    f"{chain:1}"
                    f"{last['resSeq']:>4}"
                    f"{last['iCode']:1}"
                )
                f.write(f"{ter:<80}\n")
                atom_serial += 1

        f.write(f"{'END':<80}\n")

    if verbose:
        total_output = sum(len(v) for v in atoms_by_chain.values())
        print(f"  Input atoms:      {stats['total_input']}")
        print(f"  Output atoms:     {total_output}")
        print(f"  Chains written:   {len(chain_order)} ({', '.join(chain_order)})")
        print(f"  Hydrogens removed:{stats['dropped_hydrogen']}")
        print(f"  Cap atoms removed:{stats['dropped_cap']}")
        print(f"  Altloc dropped:   {stats['dropped_altloc']}")
        print(f"  HIS renamed:      {stats['his_renamed']}")
        print(f"  ✓ Saved: {output_path}")

    return str(output_path)


def batch_clean_pdbs(
    input_dir: str,
    output_dir: Optional[str] = None,
    pattern: str = "*.pdb",
    postfix: str = "_cleaned",
    remove_hydrogens: bool = True,
    fix_histidines: bool = True,
    remove_caps: bool = True,
    assign_chains_from_segments: bool = True,
    set_occupancy: float = 1.00,
    keep_altloc: str = 'A',
    verbose: bool = True
) -> List[str]:
    """
    Convert multiple CHARMM PDB files to standard PDB format.

    Parameters
    ----------
    input_dir : str
        Directory containing input PDB files.
    output_dir : str, optional
        Directory for output files. If None and postfix is set, files are
        written to the input directory with the postfix appended.
        If an output_dir is given, files keep their original name (no postfix).
    pattern : str
        Glob pattern for PDB files.
    postfix : str
        Suffix appended to the stem of each output filename when writing to
        the same directory as the input (default: '_cleaned'). Ignored when
        output_dir is specified.
    remove_hydrogens : bool
        Remove all hydrogen atoms.
    fix_histidines : bool
        Convert HSP/HSE/HSD → HIS.
    remove_caps : bool
        Remove CHARMM capping group atoms.
    assign_chains_from_segments : bool
        Map segment IDs to unique chain letters.
    set_occupancy : float
        Force occupancy to this value.
    keep_altloc : str
        Which alternate location to keep.
    verbose : bool
        Print progress messages.

    Returns
    -------
    List[str]
        Paths to cleaned files.
    """
    input_path = Path(input_dir)
    
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    else:
        output_path = input_path
    
    pdb_files = list(input_path.glob(pattern))
    
    if verbose:
        print(f"\nCleaning {len(pdb_files)} PDB files from {input_dir}")
        print("=" * 60)
    
    cleaned_files = []
    skipped = 0
    for pdb_file in pdb_files:
        # Skip files that already have the postfix
        if postfix and postfix in pdb_file.stem:
            continue

        # When writing to a separate output dir, keep original name;
        # otherwise append the postfix to the stem.
        if output_dir:
            output_file = output_path / pdb_file.name
        else:
            output_file = output_path / f"{pdb_file.stem}{postfix}.pdb"

        # Skip if a cleaned version already exists
        if output_file.exists():
            if verbose:
                print(f"⏭ Skipping {pdb_file.name}: cleaned version already exists ({output_file.name})")
            skipped += 1
            cleaned_files.append(str(output_file))
            continue

        try:
            cleaned = clean_pdb_like_pymol(
                str(pdb_file),
                str(output_file),
                remove_hydrogens=remove_hydrogens,
                fix_histidines=fix_histidines,
                remove_caps=remove_caps,
                assign_chains_from_segments=assign_chains_from_segments,
                set_occupancy=set_occupancy,
                keep_altloc=keep_altloc,
                verbose=verbose
            )
            cleaned_files.append(cleaned)
        except Exception as e:
            if verbose:
                print(f"✗ Error cleaning {pdb_file.name}: {e}")
    
    if verbose:
        print(f"\n✓ Successfully cleaned {len(cleaned_files)}/{len(pdb_files)} files ({skipped} skipped, already existed)")
    
    return cleaned_files


def main():
    """Command-line interface."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Convert CHARMM-style PDB files to standard PDB format "
                    "compatible with AutoDock Vina receptor preparation tools."
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input PDB file or directory'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output PDB file or directory'
    )
    parser.add_argument(
        '--keep-hydrogens',
        action='store_true',
        help='Keep hydrogen atoms (default: remove them)'
    )
    parser.add_argument(
        '--no-fix-histidines',
        action='store_true',
        help='Do NOT convert HSP/HSE/HSD to HIS'
    )
    parser.add_argument(
        '--keep-caps',
        action='store_true',
        help='Keep CHARMM capping group atoms (CAY/CY/OY, NT/CAT, etc.)'
    )
    parser.add_argument(
        '--no-chain-from-segment',
        action='store_true',
        help='Do NOT assign chain IDs from segment IDs'
    )
    parser.add_argument(
        '--occupancy',
        type=float,
        default=1.00,
        help='Set occupancy value for all atoms (default: 1.00)'
    )
    parser.add_argument(
        '--keep-altloc',
        default='A',
        help='Which altloc indicator to keep (default: A)'
    )
    parser.add_argument(
        '--postfix',
        default='_cleaned',
        help='Suffix appended to output filenames in same-directory mode (default: _cleaned)'
    )
    parser.add_argument(
        '--pattern',
        default='*.pdb',
        help='File pattern for batch processing (default: *.pdb)'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        default=True,
        help='Verbose output'
    )

    args = parser.parse_args()
    input_path = Path(args.input)

    kwargs = dict(
        remove_hydrogens=not args.keep_hydrogens,
        fix_histidines=not args.no_fix_histidines,
        remove_caps=not args.keep_caps,
        assign_chains_from_segments=not args.no_chain_from_segment,
        set_occupancy=args.occupancy,
        keep_altloc=args.keep_altloc,
        verbose=args.verbose,
    )

    if input_path.is_file():
        output = args.output if args.output else str(
            input_path.parent / f"{input_path.stem}{args.postfix}.pdb"
        )
        clean_pdb_like_pymol(str(input_path), output, **kwargs)
    elif input_path.is_dir():
        batch_clean_pdbs(
            str(input_path),
            args.output,
            pattern=args.pattern,
            postfix=args.postfix,
            **kwargs,
        )
    else:
        print(f"Error: {input_path} is not a valid file or directory")
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
