#!/usr/bin/env python3
"""
Clean PDB files similar to how PyMOL does it.

PyMOL's key transformations:
1. Reorders atoms: backbone atoms first (N, CA, C, O), then sidechain
2. Adds proper element symbols in the last column
3. Removes/fixes problematic hydrogens
4. Standardizes formatting
"""

from pathlib import Path
from typing import Optional, List, Tuple
import re


def parse_pdb_line(line: str) -> dict:
    """Parse a PDB ATOM/HETATM line into components."""
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
        'occupancy': line[54:60].strip() if len(line) > 54 else '0.00',
        'tempFactor': line[60:66].strip() if len(line) > 60 else '0.00',
        'segment': line[72:76].strip() if len(line) > 72 else '',
        'element': line[76:78].strip() if len(line) > 76 else '',
        'charge': line[78:80].strip() if len(line) > 78 else '',
    }


def get_element_from_atom_name(atom_name: str) -> str:
    """Infer element symbol from atom name."""
    # Remove leading/trailing spaces and numbers
    name = atom_name.strip()
    
    # Common patterns
    if name.startswith('H'):
        return 'H'
    elif name.startswith('C'):
        return 'C'
    elif name.startswith('N'):
        return 'N'
    elif name.startswith('O'):
        return 'O'
    elif name.startswith('S'):
        return 'S'
    elif name.startswith('P'):
        return 'P'
    elif name.startswith('F'):
        return 'F'
    elif name.startswith('CL'):
        return 'CL'
    elif name.startswith('BR'):
        return 'BR'
    elif name.startswith('I'):
        return 'I'
    elif name.startswith('MG'):
        return 'MG'
    elif name.startswith('CA') and len(name) > 2:  # CA as calcium, not C-alpha
        return 'CA'
    elif name.startswith('ZN'):
        return 'ZN'
    elif name.startswith('FE'):
        return 'FE'
    
    # Default: first character
    return name[0].upper()


def get_atom_order_priority(atom_name: str) -> int:
    """Return sort priority for atom ordering (lower = earlier)."""
    name = atom_name.strip()
    
    # Backbone atoms come first
    if name == 'N':
        return 0
    elif name == 'CA':
        return 1
    elif name == 'C':
        return 2
    elif name == 'O':
        return 3
    elif name == 'OXT':
        return 4
    
    # N-terminal cap
    elif name in ['CAY', 'CY', 'OY']:
        if name == 'CAY':
            return 5
        elif name == 'CY':
            return 6
        else:  # OY
            return 7
    
    # Side chain heavy atoms
    elif name.startswith('C'):
        return 10
    elif name.startswith('N'):
        return 15
    elif name.startswith('O'):
        return 20
    elif name.startswith('S'):
        return 25
    
    # Hydrogens come last
    elif name.startswith('H'):
        # Order hydrogens by their parent atom name
        parent = name[1:] if len(name) > 1 else ''
        if parent == 'N':
            return 100
        elif parent == 'A':  # HA
            return 101
        else:
            return 150
    
    # Everything else
    return 50


def clean_pdb_like_pymol(
    input_pdb: str,
    output_pdb: str,
    remove_hydrogens: bool = False,
    fix_histidines: bool = True,
    verbose: bool = True
) -> str:
    """
    Clean PDB file to match PyMOL's output format.
    
    Parameters:
    -----------
    input_pdb : str
        Path to input PDB file
    output_pdb : str
        Path to output PDB file
    remove_hydrogens : bool
        Remove all hydrogen atoms
    fix_histidines : bool
        Convert HSP/HSE/HSD to HIS
    verbose : bool
        Print progress messages
    
    Returns:
    --------
    str : Path to output file
    """
    input_path = Path(input_pdb)
    output_path = Path(output_pdb)
    
    if verbose:
        print(f"Cleaning PDB: {input_path.name} → {output_path.name}")
    
    # Read input file
    with open(input_path, 'r') as f:
        lines = f.readlines()
    
    # Process atoms by residue
    header_lines = []
    atoms_by_residue = {}
    footer_lines = []
    in_atoms = False
    
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            in_atoms = True
            atom_data = parse_pdb_line(line)
            if atom_data:
                # Skip hydrogens if requested
                if remove_hydrogens and atom_data['name'].startswith('H'):
                    continue
                
                # Fix histidine names
                if fix_histidines and atom_data['resName'] in ['HSP', 'HSE', 'HSD']:
                    atom_data['resName'] = 'HIS'
                
                # Add element symbol if missing
                if not atom_data['element']:
                    atom_data['element'] = get_element_from_atom_name(atom_data['name'])
                
                # Group by residue
                res_key = (atom_data['chainID'], atom_data['resSeq'], atom_data['iCode'])
                if res_key not in atoms_by_residue:
                    atoms_by_residue[res_key] = []
                atoms_by_residue[res_key].append(atom_data)
        elif not in_atoms:
            header_lines.append(line)
        else:
            footer_lines.append(line)
    
    # Sort atoms within each residue (backbone first, then sidechain, then hydrogens)
    for res_key in atoms_by_residue:
        atoms_by_residue[res_key].sort(key=lambda x: (
            get_atom_order_priority(x['name']),
            x['name']
        ))
    
    # Write output file
    with open(output_path, 'w') as f:
        # Write header
        for line in header_lines:
            f.write(line)
        
        # Write atoms in residue order
        atom_serial = 1
        for res_key in sorted(atoms_by_residue.keys()):
            for atom in atoms_by_residue[res_key]:
                # Convert HETATM to ATOM for standard residues
                record = 'ATOM' if atom['resName'] not in ['HOH', 'WAT'] else 'HETATM'
                
                # Format PDB line with proper element column
                line = (
                    f"{record:<6}"
                    f"{atom_serial:>5} "
                    f"{atom['name']:>4}"
                    f"{atom['altLoc']:1}"
                    f"{atom['resName']:>3} "
                    f"{atom['chainID']:1}"
                    f"{atom['resSeq']:>4}"
                    f"{atom['iCode']:1}   "
                    f"{float(atom['x']):>8.3f}"
                    f"{float(atom['y']):>8.3f}"
                    f"{float(atom['z']):>8.3f}"
                    f"{float(atom['occupancy']):>6.2f}"
                    f"{float(atom['tempFactor']):>6.2f}"
                    f"      "
                    f"{atom['segment']:>4}"
                    f"{atom['element']:>2}"
                    f"{atom['charge']:>2}"
                    "\n"
                )
                f.write(line)
                atom_serial += 1
        
        # Write footer
        for line in footer_lines:
            f.write(line)
        
        # Add END if not present
        if not footer_lines or not footer_lines[-1].startswith('END'):
            f.write('END\n')
    
    if verbose:
        print(f"✓ Cleaned PDB saved: {output_path}")
        print(f"  Atoms processed: {atom_serial - 1}")
    
    return str(output_path)


def batch_clean_pdbs(
    input_dir: str,
    output_dir: Optional[str] = None,
    pattern: str = "*.pdb",
    remove_hydrogens: bool = False,
    fix_histidines: bool = True,
    verbose: bool = True
) -> List[str]:
    """
    Clean multiple PDB files.
    
    Parameters:
    -----------
    input_dir : str
        Directory containing input PDB files
    output_dir : str, optional
        Directory for output files (default: same as input with _cleaned suffix)
    pattern : str
        Glob pattern for PDB files
    remove_hydrogens : bool
        Remove all hydrogen atoms
    fix_histidines : bool
        Convert HSP/HSE/HSD to HIS
    verbose : bool
        Print progress messages
    
    Returns:
    --------
    List[str] : Paths to cleaned files
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
        # Skip files that are already cleaned versions
        if "_cleaned" in pdb_file.stem:
            continue

        # Add _cleaned suffix if output is same directory
        if output_dir:
            output_file = output_path / pdb_file.name
        else:
            output_file = output_path / f"{pdb_file.stem}_cleaned.pdb"

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
        description="Clean PDB files similar to PyMOL's output format"
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
        '--remove-hydrogens',
        action='store_true',
        help='Remove all hydrogen atoms'
    )
    parser.add_argument(
        '--fix-histidines',
        action='store_true',
        default=True,
        help='Convert HSP/HSE/HSD to HIS (default: True)'
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
    
    if input_path.is_file():
        # Single file
        output = args.output if args.output else str(input_path.parent / f"{input_path.stem}_cleaned.pdb")
        clean_pdb_like_pymol(
            str(input_path),
            output,
            remove_hydrogens=args.remove_hydrogens,
            fix_histidines=args.fix_histidines,
            verbose=args.verbose
        )
    elif input_path.is_dir():
        # Batch process directory
        batch_clean_pdbs(
            str(input_path),
            args.output,
            pattern=args.pattern,
            remove_hydrogens=args.remove_hydrogens,
            fix_histidines=args.fix_histidines,
            verbose=args.verbose
        )
    else:
        print(f"Error: {input_path} is not a valid file or directory")
        return 1
    
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
