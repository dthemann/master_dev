#!/usr/bin/env python3
"""
Compare docking results between Meeko and MGLTools conversion methods.
Analyzes the successfully docked poses from both workflows.
"""

import os
import re
import pandas as pd
from pathlib import Path

# Define paths
MEEKO_DOCKING_DIR = "/home/manndo/MasterProject/docking_ready_meeko/docking"
MGLTOOLS_DOCKING_DIR = "/home/manndo/MasterProject/docking_ready_mgltools/docking"

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


def main():
    print("=" * 80)
    print("COMPARISON OF DOCKING RESULTS: MEEKO vs MGLTools CONVERSION")
    print("=" * 80)
    print()
    
    # Collect results from both methods
    print("Collecting Meeko docking results...")
    meeko_results = collect_docking_results(MEEKO_DOCKING_DIR, "Meeko")
    print(f"  Found {len(meeko_results)} docking runs")
    
    print("Collecting MGLTools docking results...")
    mgltools_results = collect_docking_results(MGLTOOLS_DOCKING_DIR, "MGLTools")
    print(f"  Found {len(mgltools_results)} docking runs")
    print()
    
    # Create DataFrames
    df_meeko = pd.DataFrame(meeko_results)
    df_mgltools = pd.DataFrame(mgltools_results)
    
    # ============================================================
    # SUMMARY STATISTICS
    # ============================================================
    print("=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print()
    
    print("MEEKO Conversion Method:")
    print(f"  Total docking runs: {len(df_meeko)}")
    print(f"  Successful dockings: {df_meeko['success'].sum()}")
    print(f"  Failed dockings: {(~df_meeko['success']).sum()}")
    if df_meeko['success'].sum() > 0:
        print(f"  Best affinity (overall): {df_meeko['best_affinity'].min():.3f} kcal/mol")
        print(f"  Worst affinity (overall): {df_meeko['best_affinity'].max():.3f} kcal/mol")
        print(f"  Mean best affinity: {df_meeko['best_affinity'].mean():.3f} kcal/mol")
    print()
    
    print("MGLTools Conversion Method:")
    print(f"  Total docking runs: {len(df_mgltools)}")
    print(f"  Successful dockings: {df_mgltools['success'].sum()}")
    print(f"  Failed dockings: {(~df_mgltools['success']).sum()}")
    if df_mgltools['success'].sum() > 0:
        print(f"  Best affinity (overall): {df_mgltools['best_affinity'].min():.3f} kcal/mol")
        print(f"  Worst affinity (overall): {df_mgltools['best_affinity'].max():.3f} kcal/mol")
        print(f"  Mean best affinity: {df_mgltools['best_affinity'].mean():.3f} kcal/mol")
    print()
    
    # ============================================================
    # DETAILED COMPARISON BY PROTEIN-LIGAND PAIR
    # ============================================================
    print("=" * 80)
    print("DETAILED COMPARISON BY PROTEIN-LIGAND PAIR")
    print("=" * 80)
    print()
    
    # Normalize protein names for comparison (remove _retry1 suffix from meeko)
    df_meeko['protein_normalized'] = df_meeko['protein'].str.replace('_retry1', '', regex=False)
    df_mgltools['protein_normalized'] = df_mgltools['protein']
    
    # Get unique protein-ligand combinations
    meeko_pairs = set(zip(df_meeko['protein_normalized'], df_meeko['ligand']))
    mgltools_pairs = set(zip(df_mgltools['protein_normalized'], df_mgltools['ligand']))
    
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
    
    # ============================================================
    # AFFINITY COMPARISON FOR COMMON PAIRS
    # ============================================================
    print("=" * 80)
    print("AFFINITY COMPARISON FOR COMMON PAIRS")
    print("=" * 80)
    print()
    
    comparison_data = []
    
    print(f"{'Protein':<30} {'Ligand':<25} {'Meeko':>10} {'MGLTools':>10} {'Δ (M-MGL)':>10} {'Better':>10}")
    print("-" * 105)
    
    for protein, ligand in sorted(common_pairs):
        meeko_row = df_meeko[(df_meeko['protein_normalized'] == protein) & (df_meeko['ligand'] == ligand)]
        mgltools_row = df_mgltools[(df_mgltools['protein_normalized'] == protein) & (df_mgltools['ligand'] == ligand)]
        
        meeko_affinity = meeko_row['best_affinity'].values[0] if len(meeko_row) > 0 else None
        mgltools_affinity = mgltools_row['best_affinity'].values[0] if len(mgltools_row) > 0 else None
        
        if meeko_affinity is not None and mgltools_affinity is not None:
            delta = meeko_affinity - mgltools_affinity
            better = "Meeko" if meeko_affinity < mgltools_affinity else ("MGLTools" if mgltools_affinity < meeko_affinity else "Equal")
            
            print(f"{protein:<30} {ligand:<25} {meeko_affinity:>10.3f} {mgltools_affinity:>10.3f} {delta:>+10.3f} {better:>10}")
            
            comparison_data.append({
                'protein': protein,
                'ligand': ligand,
                'meeko_affinity': meeko_affinity,
                'mgltools_affinity': mgltools_affinity,
                'delta': delta,
                'better_method': better
            })
    
    print()
    
    # ============================================================
    # SUMMARY OF COMPARISON
    # ============================================================
    if comparison_data:
        df_comparison = pd.DataFrame(comparison_data)
        
        print("=" * 80)
        print("SUMMARY OF AFFINITY COMPARISON")
        print("=" * 80)
        print()
        
        meeko_wins = (df_comparison['better_method'] == 'Meeko').sum()
        mgltools_wins = (df_comparison['better_method'] == 'MGLTools').sum()
        ties = (df_comparison['better_method'] == 'Equal').sum()
        
        print(f"Meeko produced better (lower) affinity: {meeko_wins} times")
        print(f"MGLTools produced better (lower) affinity: {mgltools_wins} times")
        print(f"Equal affinity: {ties} times")
        print()
        
        print(f"Average delta (Meeko - MGLTools): {df_comparison['delta'].mean():.3f} kcal/mol")
        print(f"  - Negative delta means Meeko is better on average")
        print(f"  - Positive delta means MGLTools is better on average")
        print()
        
        print(f"Maximum difference: {df_comparison['delta'].abs().max():.3f} kcal/mol")
        print(f"Minimum difference: {df_comparison['delta'].abs().min():.3f} kcal/mol")
        print()
        
        # ============================================================
        # COMPARISON BY LIGAND
        # ============================================================
        print("=" * 80)
        print("COMPARISON BY LIGAND")
        print("=" * 80)
        print()
        
        for ligand in df_comparison['ligand'].unique():
            ligand_data = df_comparison[df_comparison['ligand'] == ligand]
            print(f"Ligand: {ligand}")
            print(f"  Mean Meeko affinity: {ligand_data['meeko_affinity'].mean():.3f} kcal/mol")
            print(f"  Mean MGLTools affinity: {ligand_data['mgltools_affinity'].mean():.3f} kcal/mol")
            print(f"  Mean delta: {ligand_data['delta'].mean():.3f} kcal/mol")
            meeko_better = (ligand_data['better_method'] == 'Meeko').sum()
            mgl_better = (ligand_data['better_method'] == 'MGLTools').sum()
            print(f"  Meeko better: {meeko_better}, MGLTools better: {mgl_better}")
            print()
        
        # ============================================================
        # COMPARISON BY PROTEIN
        # ============================================================
        print("=" * 80)
        print("COMPARISON BY PROTEIN")
        print("=" * 80)
        print()
        
        for protein in df_comparison['protein'].unique():
            protein_data = df_comparison[df_comparison['protein'] == protein]
            print(f"Protein: {protein}")
            print(f"  Mean Meeko affinity: {protein_data['meeko_affinity'].mean():.3f} kcal/mol")
            print(f"  Mean MGLTools affinity: {protein_data['mgltools_affinity'].mean():.3f} kcal/mol")
            print(f"  Mean delta: {protein_data['delta'].mean():.3f} kcal/mol")
            meeko_better = (protein_data['better_method'] == 'Meeko').sum()
            mgl_better = (protein_data['better_method'] == 'MGLTools').sum()
            print(f"  Meeko better: {meeko_better}, MGLTools better: {mgl_better}")
            print()
        
        # Save comparison to CSV
        output_csv = "/home/manndo/MasterProject/meeko_vs_mgltools_comparison.csv"
        df_comparison.to_csv(output_csv, index=False)
        print(f"Comparison data saved to: {output_csv}")
        print()
        
        # Also save detailed results
        all_results = pd.concat([df_meeko, df_mgltools], ignore_index=True)
        detailed_csv = "/home/manndo/MasterProject/meeko_vs_mgltools_detailed.csv"
        all_results.to_csv(detailed_csv, index=False)
        print(f"Detailed results saved to: {detailed_csv}")


if __name__ == "__main__":
    main()
