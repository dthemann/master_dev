"""
Filter prepared proteins to select the simplest versions based on priority.

For Orai proteins: Follow specific priority order
For other proteins: Alphabetically

Allows user to specify how many versions of each protein to include in the final set.
"""

import os
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple


class ProteinFilter:
    # Priority suffixes for Orai proteins
    ORAI_PRIORITY = [
        "_retry1.pdbqt",
        "_clean_retry1.pdbqt",
        "_pm_saved_retry1.pdbqt",
        "_clean_pm_saved.pdbqt",
        "_amber_clean_pm_saved.pdbqt",
    ]
    
    def __init__(self, receptor_folder: str):
        self.receptor_folder = receptor_folder
        self.proteins = self._load_proteins()
    
    def _load_proteins(self) -> Dict[str, List[str]]:
        """Load all protein PDBQT files and group by base protein name."""
        proteins = defaultdict(list)
        
        if not os.path.exists(self.receptor_folder):
            print(f"Error: Receptor folder '{self.receptor_folder}' not found")
            return proteins
        
        for filename in os.listdir(self.receptor_folder):
            if not filename.endswith(".pdbqt"):
                continue
            
            # Skip ligand files
            if "ligand" in filename.lower():
                continue
            
            # Extract base protein name (everything before the first underscore or variation suffix)
            base_name = self._extract_base_name(filename)
            if base_name:
                proteins[base_name].append(filename)
        
        return proteins
    
    def _extract_base_name(self, filename: str) -> str:
        """Extract the base protein name from the filename."""
        # Remove .pdbqt extension
        name = filename[:-6]  # Remove .pdbqt
        
        # Handle docking_ready prefix
        if name.startswith("docking_ready_"):
            name = name[14:]  # Remove "docking_ready_"
        
        # Handle cs_ prefix
        if name.startswith("cs_"):
            name = name[3:]  # Remove "cs_"
        
        # Extract base protein name
        # Match patterns like "Orai1WT-MDSnap-Fr300", "Orai1WT-START-Fr0"
        # These should match the part BEFORE any underscore suffixes
        match = re.match(r"(Orai1WT-[A-Za-z]+-Fr\d+)", name)
        if match:
            return match.group(1)
        
        # For non-Orai proteins, return the first part before underscore
        # or the whole name if no underscore
        base = name.split("_")[0] if "_" in name else name
        # Remove any trailing numbers in parentheses
        base = re.sub(r'\(\d+\)$', '', base)
        return base if base else None
    
    def _is_orai_protein(self, base_name: str) -> bool:
        """Check if the protein is an Orai protein."""
        return "Orai" in base_name
    
    def _get_priority_score(self, filename: str, base_name: str) -> Tuple[int, str]:
        """
        Get priority score for sorting.
        Returns (priority_score, filename) where lower score = higher priority.
        """
        if self._is_orai_protein(base_name):
            # Check against Orai priority list
            for priority, suffix in enumerate(self.ORAI_PRIORITY):
                if filename.endswith(suffix):
                    return (priority, filename)
            
            # Files not matching priority list get lowest priority (highest number)
            return (len(self.ORAI_PRIORITY), filename)
        else:
            # For non-Orai proteins, use alphabetical order
            return (0, filename)
    
    def filter_proteins(self, num_versions: int = 1) -> Dict[str, List[str]]:
        """
        Filter proteins to keep only the top N versions based on priority.
        
        Args:
            num_versions: Number of versions to keep for each protein (default: 1)
        
        Returns:
            Dictionary with base protein names and their selected filenames
        """
        filtered = {}
        
        for base_name in sorted(self.proteins.keys()):
            filenames = self.proteins[base_name]
            is_orai = self._is_orai_protein(base_name)
            
            if is_orai:
                # Sort by Orai priority
                sorted_files = sorted(
                    filenames,
                    key=lambda f: self._get_priority_score(f, base_name)
                )
            else:
                # Sort alphabetically
                sorted_files = sorted(filenames)
            
            # Keep top N versions
            filtered[base_name] = sorted_files[:num_versions]
        
        return filtered
    
    def print_summary(self, filtered_proteins: Dict[str, List[str]]):
        """Print a summary of filtered proteins."""
        print("\n" + "="*80)
        print("FILTERED PROTEIN SUMMARY")
        print("="*80)
        
        orai_count = 0
        other_count = 0
        
        for base_name in sorted(filtered_proteins.keys()):
            files = filtered_proteins[base_name]
            is_orai = self._is_orai_protein(base_name)
            
            if is_orai:
                orai_count += len(files)
                print(f"\n[ORAI] {base_name}")
            else:
                other_count += len(files)
                print(f"\n[OTHER] {base_name}")
            
            for i, filename in enumerate(files, 1):
                print(f"  {i}. {filename}")
        
        print("\n" + "="*80)
        print(f"Total Orai proteins selected: {orai_count}")
        print(f"Total other proteins selected: {other_count}")
        print(f"Total protein files selected: {orai_count + other_count}")
        print("="*80 + "\n")
    
    def save_filtered_list(self, filtered_proteins: Dict[str, List[str]], 
                          output_file: str = "filtered_proteins.txt"):
        """Save the filtered protein list to a file."""
        with open(output_file, 'w') as f:
            for base_name in sorted(filtered_proteins.keys()):
                for filename in filtered_proteins[base_name]:
                    f.write(f"{filename}\n")
        print(f"Filtered protein list saved to: {output_file}\n")


def main():
    receptor_folder = "/home/manndo/MasterProject/prepared/receptors"
    
    print("Protein Filter for Docking Preparation")
    print("-" * 80)
    
    filter_obj = ProteinFilter(receptor_folder)
    
    print(f"\nTotal unique protein bases found: {len(filter_obj.proteins)}")
    print(f"Total protein PDBQT files: {sum(len(v) for v in filter_obj.proteins.values())}\n")
    
    # Interactive menu
    while True:
        try:
            num_versions = int(input(
                "How many versions of each protein should make it into the final set? "
                "(e.g., 1, 2, 3): "
            ))
            if num_versions < 1:
                print("Please enter a number >= 1")
                continue
            break
        except ValueError:
            print("Please enter a valid integer")
    
    # Filter proteins
    filtered = filter_obj.filter_proteins(num_versions=num_versions)
    
    # Print summary
    filter_obj.print_summary(filtered)
    
    # Save to file
    output_file = f"filtered_proteins_{num_versions}_version{'s' if num_versions > 1 else ''}.txt"
    filter_obj.save_filtered_list(filtered, output_file)
    
    # Show Orai priority details
    print("\nORAI PROTEIN PRIORITY ORDER:")
    print("-" * 80)
    for i, suffix in enumerate(filter_obj.ORAI_PRIORITY, 1):
        print(f"{i}. {suffix}")
    print("\nOther proteins are sorted alphabetically")
    print("-" * 80 + "\n")
    
    # Option to save specific versions
    print("\nOPTIONS:")
    print("1. Generate filtered lists for multiple version counts")
    print("2. View detailed information for a specific protein")
    print("3. Exit")
    
    choice = input("\nSelect an option (1-3): ").strip()
    
    if choice == "1":
        multi_versions = input(
            "Enter version counts separated by commas (e.g., 1,2,3): "
        )
        try:
            versions = [int(v.strip()) for v in multi_versions.split(",")]
            for v in versions:
                if v >= 1:
                    filtered = filter_obj.filter_proteins(num_versions=v)
                    output_file = f"filtered_proteins_{v}_version{'s' if v > 1 else ''}.txt"
                    filter_obj.save_filtered_list(filtered, output_file)
            print("\nAll filtered lists generated successfully!")
        except ValueError:
            print("Invalid input format")
    
    elif choice == "2":
        protein_name = input("Enter protein base name (e.g., Orai1WT-MDSnap-Fr300): ").strip()
        if protein_name in filter_obj.proteins:
            print(f"\nAll versions of {protein_name}:")
            is_orai = filter_obj._is_orai_protein(protein_name)
            if is_orai:
                print("(Orai priority order shown with [P] marker)")
                sorted_files = sorted(
                    filter_obj.proteins[protein_name],
                    key=lambda f: filter_obj._get_priority_score(f, protein_name)
                )
                for i, f in enumerate(sorted_files, 1):
                    match = next(
                        (p for p, suffix in enumerate(filter_obj.ORAI_PRIORITY, 1)
                         if f.endswith(suffix)),
                        None
                    )
                    priority_marker = f" [P{match}]" if match else " [P-]"
                    print(f"{i}. {f}{priority_marker}")
            else:
                print("(Alphabetically sorted)")
                for i, f in enumerate(sorted(filter_obj.proteins[protein_name]), 1):
                    print(f"{i}. {f}")
        else:
            print(f"Protein '{protein_name}' not found")


if __name__ == "__main__":
    main()
