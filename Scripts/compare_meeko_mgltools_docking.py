#!/usr/bin/env python3
"""
Compare docking results between Meeko and MGLTools conversion methods.
Analyzes the successfully docked poses from both workflows.
"""

import os
import re
import pandas as pd
from pathlib import Path

# Define working directory paths
wd = {
    "meeko": "/home/manndo/MasterProject/docking_ready_meeko/docking",
    "mgltools": "/home/manndo/MasterProject/docking_ready_mgltools/docking",
    "output": "/home/manndo/MasterProject"
}

MEEKO_DOCKING_DIR = wd["meeko"]
MGLTOOLS_DOCKING_DIR = wd["mgltools"]

def parse_vina_log(log_path):
    """Parse a Vina log file and extract docking results."""
    results = {
        'success': False,
        'modes': [],
        'best_affinity': None,
        'receptor': None,
        'ligand': None,
        'grid_center': None,
        'grid_size': None,
        'exhaustiveness': None
    }
    
    if not os.path.exists(log_path):
        return results
    
    try:
        with open(log_path, 'r') as f:
            content = f.read()
        
        # Extract receptor and ligand
        receptor_match = re.search(r'Rigid receptor: (.+)', content)
        ligand_match = re.search(r'Ligand: (.+)', content)
        if receptor_match:
            results['receptor'] = receptor_match.group(1).strip()
        if ligand_match:
            results['ligand'] = ligand_match.group(1).strip()
        
        # Extract grid info
        grid_center_match = re.search(r'Grid center: X ([\d.-]+) Y ([\d.-]+) Z ([\d.-]+)', content)
        grid_size_match = re.search(r'Grid size\s*: X ([\d.]+) Y ([\d.]+) Z ([\d.]+)', content)
        exhaustiveness_match = re.search(r'Exhaustiveness: (\d+)', content)
        
        if grid_center_match:
            results['grid_center'] = (float(grid_center_match.group(1)), 
                                      float(grid_center_match.group(2)), 
                                      float(grid_center_match.group(3)))
        if grid_size_match:
            results['grid_size'] = (float(grid_size_match.group(1)), 
                                    float(grid_size_match.group(2)), 
                                    float(grid_size_match.group(3)))
        if exhaustiveness_match:
            results['exhaustiveness'] = int(exhaustiveness_match.group(1))
        
        # Extract docking modes/poses
        mode_pattern = r'^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
        modes = re.findall(mode_pattern, content, re.MULTILINE)
        
        if modes:
            results['success'] = True
            for mode in modes:
                mode_num, affinity, rmsd_lb, rmsd_ub = mode
                results['modes'].append({
                    'mode': int(mode_num),
                    'affinity': float(affinity),
                    'rmsd_lb': float(rmsd_lb),
                    'rmsd_ub': float(rmsd_ub)
                })
            results['best_affinity'] = results['modes'][0]['affinity'] if results['modes'] else None
    
    except Exception as e:
        print(f"Error parsing {log_path}: {e}")
    
    return results


def collect_docking_results(docking_dir, method_name):
    """Collect all docking results from a directory."""
    results = []
    
    log_files = [f for f in os.listdir(docking_dir) if f.endswith('_vina_out.log')]
    
    for log_file in log_files:
        log_path = os.path.join(docking_dir, log_file)
        
        # Parse the filename to extract protein and ligand names
        # Format: Protein__Ligand_vina_out.log
        base_name = log_file.replace('_vina_out.log', '')
        parts = base_name.split('__')
        
        if len(parts) == 2:
            protein_name = parts[0]
            ligand_name = parts[1]
        else:
            protein_name = base_name
            ligand_name = "unknown"
        
        # Parse the log file
        parsed = parse_vina_log(log_path)
        
        # Create result entry
        result = {
            'method': method_name,
            'protein': protein_name,
            'ligand': ligand_name,
            'log_file': log_file,
            'success': parsed['success'],
            'best_affinity': parsed['best_affinity'],
            'num_modes': len(parsed['modes']),
            'exhaustiveness': parsed['exhaustiveness']
        }
        
        # Add all modes' affinities
        for i, mode in enumerate(parsed['modes'][:9], 1):
            result[f'mode_{i}_affinity'] = mode['affinity']
            result[f'mode_{i}_rmsd_lb'] = mode['rmsd_lb']
            result[f'mode_{i}_rmsd_ub'] = mode['rmsd_ub']
        
        results.append(result)
    
    return results


class CompareDockingTools:
    """Compare docking results between Meeko and MGLTools conversion methods."""
    
    def __init__(self, meeko_dir=None, mgltools_dir=None, output_dir=None, verbose=True):
        """Initialize comparison with directory paths.
        
        Parameters:
        -----------
        meeko_dir : str or Path, optional
            Path to Meeko docking directory (default: from wd dict)
        mgltools_dir : str or Path, optional
            Path to MGLTools docking directory (default: from wd dict)
        output_dir : str or Path, optional
            Path for output files (default: from wd dict)
        verbose : bool, optional
            Print detailed output (default: True)
        """
        self.meeko_dir = str(meeko_dir) if meeko_dir else MEEKO_DOCKING_DIR
        self.mgltools_dir = str(mgltools_dir) if mgltools_dir else MGLTOOLS_DOCKING_DIR
        self.output_dir = str(output_dir) if output_dir else wd.get("output", "/home/manndo/MasterProject")
        self.verbose = verbose
        
        # Results storage
        self.df_meeko = None
        self.df_mgltools = None
        self.df_comparison = None
        self.comparison_data = []
        
    def collect_results(self):
        """Collect docking results from both directories."""
        if self.verbose:
            print("=" * 80)
            print("COMPARISON OF DOCKING RESULTS: MEEKO vs MGLTools CONVERSION")
            print("=" * 80)
            print()
            print(f"Meeko directory: {self.meeko_dir}")
            print(f"MGLTools directory: {self.mgltools_dir}")
            print(f"Output directory: {self.output_dir}")
            print()
        
        # Check if directories exist
        if not os.path.exists(self.meeko_dir):
            if self.verbose:
                print(f"⚠ Warning: Meeko directory not found: {self.meeko_dir}")
            meeko_results = []
        else:
            if self.verbose:
                print("Collecting Meeko docking results...")
            meeko_results = collect_docking_results(self.meeko_dir, "Meeko")
            if self.verbose:
                print(f"  Found {len(meeko_results)} docking runs")
        
        if not os.path.exists(self.mgltools_dir):
            if self.verbose:
                print(f"⚠ Warning: MGLTools directory not found: {self.mgltools_dir}")
            mgltools_results = []
        else:
            if self.verbose:
                print("Collecting MGLTools docking results...")
            mgltools_results = collect_docking_results(self.mgltools_dir, "MGLTools")
            if self.verbose:
                print(f"  Found {len(mgltools_results)} docking runs")
        
        if self.verbose:
            print()
        
        self.df_meeko = pd.DataFrame(meeko_results)
        self.df_mgltools = pd.DataFrame(mgltools_results)
        
        return self.df_meeko, self.df_mgltools
    
    def analyze(self):
        """Run complete comparison analysis."""
        if self.df_meeko is None or self.df_mgltools is None:
            self.collect_results()
        
        self._print_summary_statistics()
        self._compare_pairs()
        self._compare_affinities()
        self._save_results()
        
        return self.df_comparison
    
    def _print_summary_statistics(self):
        """Print summary statistics for both methods."""
        if not self.verbose:
            return
            
        print("=" * 80)
        print("SUMMARY STATISTICS")
        print("=" * 80)
        print()
        
        print("MEEKO Conversion Method:")
        print(f"  Total docking runs: {len(self.df_meeko)}")
        if len(self.df_meeko) > 0 and 'success' in self.df_meeko.columns:
            print(f"  Successful dockings: {self.df_meeko['success'].sum()}")
            print(f"  Failed dockings: {(~self.df_meeko['success']).sum()}")
            if self.df_meeko['success'].sum() > 0:
                print(f"  Best affinity (overall): {self.df_meeko['best_affinity'].min():.3f} kcal/mol")
                print(f"  Worst affinity (overall): {self.df_meeko['best_affinity'].max():.3f} kcal/mol")
                print(f"  Mean best affinity: {self.df_meeko['best_affinity'].mean():.3f} kcal/mol")
        else:
            print("  No docking results found")
        print()
        
        print("MGLTools Conversion Method:")
        print(f"  Total docking runs: {len(self.df_mgltools)}")
        if len(self.df_mgltools) > 0 and 'success' in self.df_mgltools.columns:
            print(f"  Successful dockings: {self.df_mgltools['success'].sum()}")
            print(f"  Failed dockings: {(~self.df_mgltools['success']).sum()}")
            if self.df_mgltools['success'].sum() > 0:
                print(f"  Best affinity (overall): {self.df_mgltools['best_affinity'].min():.3f} kcal/mol")
                print(f"  Worst affinity (overall): {self.df_mgltools['best_affinity'].max():.3f} kcal/mol")
                print(f"  Mean best affinity: {self.df_mgltools['best_affinity'].mean():.3f} kcal/mol")
        else:
            print("  No docking results found")
        print()
    
    def _compare_pairs(self):
        """Compare protein-ligand pairs between methods."""
        if not self.verbose:
            return
        
        # Check if we have data to compare
        if len(self.df_meeko) == 0 and len(self.df_mgltools) == 0:
            print("⚠ No docking results to compare")
            return set()
            
        print("=" * 80)
        print("DETAILED COMPARISON BY PROTEIN-LIGAND PAIR")
        print("=" * 80)
        print()
        
        # Normalize protein names
        if len(self.df_meeko) > 0 and 'protein' in self.df_meeko.columns:
            self.df_meeko['protein_normalized'] = self.df_meeko['protein'].str.replace('_retry1', '', regex=False)
        if len(self.df_mgltools) > 0 and 'protein' in self.df_mgltools.columns:
            self.df_mgltools['protein_normalized'] = self.df_mgltools['protein']
        
        # Get unique pairs
        meeko_pairs = set(zip(self.df_meeko['protein_normalized'], self.df_meeko['ligand']))
        mgltools_pairs = set(zip(self.df_mgltools['protein_normalized'], self.df_mgltools['ligand']))
        
        common_pairs = meeko_pairs & mgltools_pairs
        meeko_only = meeko_pairs - mgltools_pairs
        mgltools_only = mgltools_pairs - meeko_pairs
        
        print(f"Common protein-ligand pairs: {len(common_pairs)}")
        print(f"Meeko-only pairs: {len(meeko_only)}")
        print(f"MGLTools-only pairs: {len(mgltools_only)}")
        print()
        
        if meeko_only:
            print("Pairs only in Meeko results:")
            for pair in sorted(meeko_only):
                print(f"  {pair[0]} + {pair[1]}")
            print()
        
        if mgltools_only:
            print("Pairs only in MGLTools results:")
            for pair in sorted(mgltools_only):
                print(f"  {pair[0]} + {pair[1]}")
            print()
        
        return common_pairs
    
    def _compare_affinities(self):
        """Compare binding affinities for common pairs."""
        # Check if we have data to compare
        if len(self.df_meeko) == 0 or len(self.df_mgltools) == 0:
            if self.verbose:
                print("⚠ Insufficient data for affinity comparison")
            return
        
        # Normalize protein names if not already done
        if 'protein_normalized' not in self.df_meeko.columns and 'protein' in self.df_meeko.columns:
            self.df_meeko['protein_normalized'] = self.df_meeko['protein'].str.replace('_retry1', '', regex=False)
        if 'protein_normalized' not in self.df_mgltools.columns and 'protein' in self.df_mgltools.columns:
            self.df_mgltools['protein_normalized'] = self.df_mgltools['protein']
        
        meeko_pairs = set(zip(self.df_meeko['protein_normalized'], self.df_meeko['ligand']))
        mgltools_pairs = set(zip(self.df_mgltools['protein_normalized'], self.df_mgltools['ligand']))
        common_pairs = meeko_pairs & mgltools_pairs
        
        if not self.verbose:
            # Silent mode - just collect data
            for protein, ligand in sorted(common_pairs):
                meeko_row = self.df_meeko[(self.df_meeko['protein_normalized'] == protein) & (self.df_meeko['ligand'] == ligand)]
                mgltools_row = self.df_mgltools[(self.df_mgltools['protein_normalized'] == protein) & (self.df_mgltools['ligand'] == ligand)]
                
                meeko_affinity = meeko_row['best_affinity'].values[0] if len(meeko_row) > 0 else None
                mgltools_affinity = mgltools_row['best_affinity'].values[0] if len(mgltools_row) > 0 else None
                
                if meeko_affinity is not None and mgltools_affinity is not None:
                    delta = meeko_affinity - mgltools_affinity
                    better = "Meeko" if meeko_affinity < mgltools_affinity else ("MGLTools" if mgltools_affinity < meeko_affinity else "Equal")
                    
                    self.comparison_data.append({
                        'protein': protein,
                        'ligand': ligand,
                        'meeko_affinity': meeko_affinity,
                        'mgltools_affinity': mgltools_affinity,
                        'delta': delta,
                        'better_method': better
                    })
        else:
            print("=" * 80)
            print("AFFINITY COMPARISON FOR COMMON PAIRS")
            print("=" * 80)
            print()
            
            print(f"{'Protein':<30} {'Ligand':<25} {'Meeko':>10} {'MGLTools':>10} {'Δ (M-MGL)':>10} {'Better':>10}")
            print("-" * 105)
            
            for protein, ligand in sorted(common_pairs):
                meeko_row = self.df_meeko[(self.df_meeko['protein_normalized'] == protein) & (self.df_meeko['ligand'] == ligand)]
                mgltools_row = self.df_mgltools[(self.df_mgltools['protein_normalized'] == protein) & (self.df_mgltools['ligand'] == ligand)]
                
                meeko_affinity = meeko_row['best_affinity'].values[0] if len(meeko_row) > 0 else None
                mgltools_affinity = mgltools_row['best_affinity'].values[0] if len(mgltools_row) > 0 else None
                
                if meeko_affinity is not None and mgltools_affinity is not None:
                    delta = meeko_affinity - mgltools_affinity
                    better = "Meeko" if meeko_affinity < mgltools_affinity else ("MGLTools" if mgltools_affinity < meeko_affinity else "Equal")
                    
                    print(f"{protein:<30} {ligand:<25} {meeko_affinity:>10.3f} {mgltools_affinity:>10.3f} {delta:>+10.3f} {better:>10}")
                    
                    self.comparison_data.append({
                        'protein': protein,
                        'ligand': ligand,
                        'meeko_affinity': meeko_affinity,
                        'mgltools_affinity': mgltools_affinity,
                        'delta': delta,
                        'better_method': better
                    })
            print()
        
        if self.comparison_data:
            self.df_comparison = pd.DataFrame(self.comparison_data)
            
            if self.verbose:
                self._print_comparison_summary()
    
    def _print_comparison_summary(self):
        """Print summary of affinity comparison."""
        if self.df_comparison is None or not self.verbose:
            return
            
        print("=" * 80)
        print("SUMMARY OF AFFINITY COMPARISON")
        print("=" * 80)
        print()
        
        meeko_wins = (self.df_comparison['better_method'] == 'Meeko').sum()
        mgltools_wins = (self.df_comparison['better_method'] == 'MGLTools').sum()
        ties = (self.df_comparison['better_method'] == 'Equal').sum()
        
        print(f"Meeko produced better (lower) affinity: {meeko_wins} times")
        print(f"MGLTools produced better (lower) affinity: {mgltools_wins} times")
        print(f"Equal affinity: {ties} times")
        print()
        
        print(f"Average delta (Meeko - MGLTools): {self.df_comparison['delta'].mean():.3f} kcal/mol")
        print(f"  - Negative delta means Meeko is better on average")
        print(f"  - Positive delta means MGLTools is better on average")
        print()
        
        print(f"Maximum difference: {self.df_comparison['delta'].abs().max():.3f} kcal/mol")
        print(f"Minimum difference: {self.df_comparison['delta'].abs().min():.3f} kcal/mol")
        print()
    
    def _save_results(self):
        """Save comparison results to CSV files."""
        os.makedirs(self.output_dir, exist_ok=True)
        
        if self.df_comparison is not None:
            output_csv = os.path.join(self.output_dir, "meeko_vs_mgltools_comparison.csv")
            self.df_comparison.to_csv(output_csv, index=False)
            if self.verbose:
                print(f"Comparison data saved to: {output_csv}")
                print()
        
        # Save detailed results
        all_results = pd.concat([self.df_meeko, self.df_mgltools], ignore_index=True)
        detailed_csv = os.path.join(self.output_dir, "meeko_vs_mgltools_detailed.csv")
        all_results.to_csv(detailed_csv, index=False)
        if self.verbose:
            print(f"Detailed results saved to: {detailed_csv}")


def main(meeko_dir=None, mgltools_dir=None, output_dir=None):
    """Compare docking results between Meeko and MGLTools.
    
    Parameters:
    -----------
    meeko_dir : str, optional
        Path to Meeko docking directory (default: from wd dict)
    mgltools_dir : str, optional
        Path to MGLTools docking directory (default: from wd dict)
    output_dir : str, optional
        Path for output files (default: from wd dict)
    """
    # Use the CompareDockingTools class
    comparator = CompareDockingTools(
        meeko_dir=meeko_dir,
        mgltools_dir=mgltools_dir,
        output_dir=output_dir,
        verbose=True
    )
    
    # Run complete analysis
    comparator.analyze()


if __name__ == "__main__":
    main()
