#!/usr/bin/env python3
"""
Filter prepared proteins and ligands for Autodock Vina docking.

Prioritizes protein versions based on:
- For Orais: specific priority list (_retry1 > _clean_retry1 > _pm_saved_retry1, etc.)
- For others: alphabetically sorted

Allows specifying how many versions to keep for each protein-ligand combination.
"""

import os
from pathlib import Path
from collections import defaultdict
import re


# Priority order for Orai proteins
ORAI_PRIORITY = [
    "_retry1.pdbqt",
    "_clean_retry1.pdbqt",
    "_pm_saved_retry1.pdbqt",
    "_clean_pm_saved.pdbqt",
    "_amber_clean_pm_saved.pdbqt",
]


def extract_base_protein_name(filename):
    """Extract the base protein name without version suffix."""
    # Remove .pdbqt extension
    name = filename.replace(".pdbqt", "")
    
    # Remove common version suffixes
    # Match patterns like _retry1, _clean, _pm_saved, _amber, etc.
    base = re.sub(r'(_retry\d+|_clean|_pm_saved|_amber|_roundtrip|_pymol_fixed|_protein|_docking_ready.*|_\(\d+\))+$', '', name)
    
    return base


def is_orai_protein(filename):
    """Check if filename contains Orai protein."""
    return "Orai" in filename


def get_priority_score(filename, is_orai):
    """
    Calculate priority score for filtering.
    Lower score = higher priority (kept first).
    """
    if is_orai:
        # Check against Orai priority list
        for idx, priority_suffix in enumerate(ORAI_PRIORITY):
            if filename.endswith(priority_suffix):
                return idx
        # If not in priority list, give it a high score (low priority)
        return len(ORAI_PRIORITY)
    else:
        # For non-Orai: alphabetically (use filename as-is)
        return (len(ORAI_PRIORITY) + 1, filename)


def filter_proteins_and_ligands(receptor_dir, ligand_dir, max_versions=1, output_dir=None):
    """
    Filter prepared proteins and ligands.
    
    Parameters:
    -----------
    receptor_dir : str
        Path to prepared receptors directory
    ligand_dir : str
        Path to prepared ligands directory
    max_versions : int
        Maximum number of versions to keep for each protein-ligand combination
    output_dir : str, optional
        Directory to save filtered file lists. If None, only returns data.
    
    Returns:
    --------
    dict : Dictionary with filtered proteins and ligands
    """
    
    receptor_dir = Path(receptor_dir)
    ligand_dir = Path(ligand_dir)
    
    # Get all receptor and ligand files
    receptors = sorted([f.name for f in receptor_dir.glob("*.pdbqt")])
    ligands = sorted([f.name for f in ligand_dir.glob("*.pdbqt")])
    
    print(f"Found {len(receptors)} receptor files and {len(ligands)} ligand files")
    print(f"\nGrouping by protein base name (keeping up to {max_versions} version(s) per protein)...\n")
    
    # Group proteins by base name
    protein_groups = defaultdict(list)
    for protein in receptors:
        base_name = extract_base_protein_name(protein)
        is_orai = is_orai_protein(protein)
        priority_score = get_priority_score(protein, is_orai)
        protein_groups[base_name].append((protein, priority_score, is_orai))
    
    # Sort each group by priority and select top versions
    filtered_proteins = []
    protein_selection_info = []
    
    for base_name in sorted(protein_groups.keys()):
        proteins_for_base = protein_groups[base_name]
        is_orai = proteins_for_base[0][2]
        
        if is_orai:
            # Sort by priority score (lower = higher priority)
            sorted_proteins = sorted(proteins_for_base, key=lambda x: x[1])
        else:
            # Sort alphabetically
            sorted_proteins = sorted(proteins_for_base, key=lambda x: x[0])
        
        # Select top max_versions
        selected = sorted_proteins[:max_versions]
        filtered_proteins.extend([p[0] for p in selected])
        
        # Store selection info
        selection_info = {
            'base_name': base_name,
            'is_orai': is_orai,
            'total_available': len(proteins_for_base),
            'selected_versions': [p[0] for p in selected],
            'not_selected': [p[0] for p in sorted_proteins[max_versions:]]
        }
        protein_selection_info.append(selection_info)
    
    # Display summary
    for info in protein_selection_info:
        status = "ORAI" if info['is_orai'] else "OTHER"
        print(f"[{status}] {info['base_name']}")
        print(f"  Total available: {info['total_available']}")
        print(f"  Keeping {len(info['selected_versions'])} version(s):")
        for ver in info['selected_versions']:
            print(f"    ✓ {ver}")
        if info['not_selected']:
            print(f"  Not using ({len(info['not_selected'])}):")
            for ver in info['not_selected'][:3]:  # Show first 3
                print(f"    ✗ {ver}")
            if len(info['not_selected']) > 3:
                print(f"    ... and {len(info['not_selected']) - 3} more")
        print()
    
    # Match ligands (keep all available ligands as they have fewer versions)
    filtered_ligands = ligands
    
    result = {
        'filtered_proteins': filtered_proteins,
        'filtered_ligands': filtered_ligands,
        'selection_info': protein_selection_info,
        'total_protein_combinations': len(filtered_proteins) * len(filtered_ligands)
    }
    
    print(f"\n{'='*60}")
    print(f"FILTERING SUMMARY")
    print(f"{'='*60}")
    print(f"Filtered proteins: {len(filtered_proteins)}")
    print(f"Available ligands: {len(filtered_ligands)}")
    print(f"Total protein-ligand combinations: {len(filtered_proteins) * len(filtered_ligands)}")
    print(f"{'='*60}\n")
    
    # Save results if output_dir specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save protein list
        protein_list_file = output_dir / f"filtered_proteins_v{max_versions}.txt"
        with open(protein_list_file, 'w') as f:
            for protein in filtered_proteins:
                f.write(f"{protein}\n")
        print(f"✓ Saved protein list to: {protein_list_file}")
        
        # Save ligand list
        ligand_list_file = output_dir / f"filtered_ligands_v{max_versions}.txt"
        with open(ligand_list_file, 'w') as f:
            for ligand in filtered_ligands:
                f.write(f"{ligand}\n")
        print(f"✓ Saved ligand list to: {ligand_list_file}")
        
        # Save detailed report
        report_file = output_dir / f"filtering_report_v{max_versions}.txt"
        with open(report_file, 'w') as f:
            f.write(f"PROTEIN FILTERING REPORT (max_versions={max_versions})\n")
            f.write(f"{'='*60}\n\n")
            for info in protein_selection_info:
                f.write(f"{'ORAI' if info['is_orai'] else 'OTHER'} - {info['base_name']}\n")
                f.write(f"  Total available: {info['total_available']}\n")
                f.write(f"  Selected versions:\n")
                for ver in info['selected_versions']:
                    f.write(f"    - {ver}\n")
                f.write(f"\n")
        print(f"✓ Saved detailed report to: {report_file}")
        print()
    
    return result


def main():
    """Main execution function."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Filter prepared proteins and ligands for Autodock Vina docking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Keep only the best version for each protein
  python filter_prepared_proteins.py --max-versions 1
  
  # Keep top 3 versions for each protein
  python filter_prepared_proteins.py --max-versions 3
  
  # Keep all versions
  python filter_prepared_proteins.py --max-versions 999
        """
    )
    
    parser.add_argument(
        '--receptor-dir',
        type=str,
        default='/home/manndo/MasterProject/prepared/receptors',
        help='Path to prepared receptors directory'
    )
    parser.add_argument(
        '--ligand-dir',
        type=str,
        default='/home/manndo/MasterProject/prepared/ligands',
        help='Path to prepared ligands directory'
    )
    parser.add_argument(
        '--max-versions',
        type=int,
        default=1,
        help='Maximum number of versions to keep for each protein (default: 1)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='/home/manndo/MasterProject/docking_ready',
        help='Directory to save filtered file lists (default: docking_ready)'
    )
    parser.add_argument(
        '--no-save',
        action='store_true',
        help='Do not save output files'
    )
    
    args = parser.parse_args()
    
    output_dir = None if args.no_save else args.output_dir
    
    result = filter_proteins_and_ligands(
        receptor_dir=args.receptor_dir,
        ligand_dir=args.ligand_dir,
        max_versions=args.max_versions,
        output_dir=output_dir
    )
    
    return result


if __name__ == "__main__":
    main()
