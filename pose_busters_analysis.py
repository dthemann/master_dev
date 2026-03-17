#!/usr/bin/env python3
"""
PoseBusters Results Analysis & Visualization
=============================================
Standalone script to analyze PoseBusters validation results.
Produces survival statistics, copies proved poses, and generates plots.

Usage:
    python pose_busters_analysis.py --results-dir /path/to/posebusters_results/dock
    python pose_busters_analysis.py --results-dir /path/to/posebusters_results/mol --output-dir /path/to/output

Arguments:
    --results-dir   Directory containing posebusters_filtered_results.csv
    --output-dir    Where to save plots and analysis outputs (default: same as --results-dir)
    --wd            Working directory for locating source pose files (for copying proved poses)
"""
import argparse
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

# ============================================================================
# CLI ARGUMENTS
# ============================================================================
_parser = argparse.ArgumentParser(
    description="Analyze and visualize PoseBusters validation results",
)
_parser.add_argument("--results-dir", type=str, required=True,
                     help="Directory containing posebusters_filtered_results.csv")
_parser.add_argument("--output-dir", type=str, default=None,
                     help="Output directory for plots and analysis (default: same as --results-dir)")
_parser.add_argument("--wd", type=str, default=None,
                     help="Working directory for locating source pose files (for copy step)")

_args, _unknown = _parser.parse_known_args()

results_dir = Path(_args.results_dir)
output_dir = Path(_args.output_dir) if _args.output_dir else results_dir
wd = _args.wd

output_dir.mkdir(parents=True, exist_ok=True)

print(f"Results directory: {results_dir}")
print(f"Output directory:  {output_dir}")
if wd:
    print(f"Working directory: {wd}")

# ============================================================================
# CONSTANTS
# ============================================================================
METADATA_COLS = {
    'file_path', 'filepath', 'file', 'path', 'sdf_file', 'sdf_path',
    'method', 'docking_method', 'protein', 'ligand', 'pose_rank', 'rank',
    'molecule', 'mol_name', 'name', 'complex', 'protein_path', 'ligand_path',
    'mol_pred', 'mol_true', 'mol_cond',
}

EXCLUDE_COLS = {
    'mol_true_loaded', 'mol_cond_loaded',
    'number_short_outlier_bonds', 'number_long_outlier_bonds',
    'number_outlier_angles', 'number_clashes',
    'number_non-aromatic_rings_pass', 'number_aromatic_rings_pass',
    'number_non-aromatic_rings_checked', 'number_aromatic_rings_checked',
    'number_double_bonds_checked', 'number_double_bonds_pass',
    'number_valid_bonds', 'number_valid_angles', 'number_valid_noncov_pairs',
    'number_noncov_pairs', 'number_bonds', 'number_angles',
    'num_h_added',
    'not_too_far_away_organic_cofactors',
    'not_too_far_away_inorganic_cofactors',
    'not_too_far_away_waters',
}

TEST_DISPLAY = {
    'mol_pred_loaded': 'Molecule Loaded',
    'sanitization': 'Sanitization',
    'inchi_convertible': 'InChI Convertible',
    'all_atoms_connected': 'All Atoms Connected',
    'no_radicals': 'No Radicals',
    'bond_lengths': 'Bond Lengths',
    'bond_angles': 'Bond Angles',
    'internal_steric_clash': 'No Steric Clash',
    'aromatic_ring_flatness': 'Aromatic Flatness',
    'non-aromatic_ring_non-flatness': 'Non-Arom. Ring Shape',
    'double_bond_flatness': 'Double Bond Flatness',
    'internal_energy': 'Internal Energy',
    'passes_valence_checks': 'Valence Checks',
    'passes_kekulization': 'Kekulization',
    'no_radicals_before_sanitization': 'No Pre-Sanit. Radicals',
}

COLOR_PALETTE = [
    '#3498db', '#e74c3c', '#2ecc71', '#9b59b6', '#f39c12',
    '#1abc9c', '#e67e22', '#34495e', '#d35400', '#8e44ad',
]


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
def detect_test_columns(df):
    """Identify boolean PoseBusters test columns in a DataFrame."""
    test_cols = []
    for col in df.columns:
        cl = col.lower().strip()
        if cl in METADATA_COLS or col in EXCLUDE_COLS or cl in EXCLUDE_COLS:
            continue
        if cl.startswith('number_') or cl.startswith('num_'):
            continue
        uv = set(df[col].dropna().unique())
        if uv.issubset({True, False, 1, 0, 1.0, 0.0, 'True', 'False', 'true', 'false'}):
            test_cols.append(col)
    if not test_cols:
        test_cols = [c for c in df.columns if df[c].dtype == bool]
    return test_cols


def coerce_bool_columns(df, test_cols):
    """Convert test columns to proper bool dtype."""
    for tc in test_cols:
        if df[tc].dtype == object:
            df[tc] = df[tc].map({'True': True, 'true': True, 'False': False, 'false': False})
        df[tc] = df[tc].astype(bool)
    return df


def find_results_csv(directory):
    """Find the PoseBusters results CSV in a directory."""
    candidates = [
        directory / "posebusters_filtered_results.csv",
    ]
    for c in candidates:
        if c.exists():
            return c
    # Fallback: search for any results-like CSV
    csvs = list(directory.glob("*.csv"))
    for c in csvs:
        if 'result' in c.name.lower() or 'posebust' in c.name.lower():
            return c
    if csvs:
        return csvs[0]
    return None


# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================
print("\n" + "=" * 100)
print("POSEBUSTERS RESULTS ANALYSIS")
print("=" * 100)

results_csv = find_results_csv(results_dir)
if results_csv is None:
    print(f"\nERROR: No PoseBusters results CSV found in {results_dir}")
    sys.exit(1)

print(f"\nLoading: {results_csv}")
df = pd.read_csv(results_csv)
print(f"Total entries: {len(df)}")
print(f"Columns: {list(df.columns)}")

# ============================================================================
# STEP 2: IDENTIFY TEST COLUMNS
# ============================================================================
test_cols = detect_test_columns(df)
if not test_cols:
    print("\nERROR: Could not auto-detect test columns.")
    for col in df.columns:
        print(f"  {col}: {df[col].dtype}, unique: {df[col].nunique()}")
    sys.exit(1)

df = coerce_bool_columns(df, test_cols)
df["all_passed"] = df[test_cols].all(axis=1)

print(f"\nIdentified {len(test_cols)} PoseBusters test columns:")
for tc in test_cols:
    n_pass = df[tc].sum()
    print(f"  {tc}: {n_pass}/{len(df)} passed ({100*n_pass/len(df):.1f}%)")

excluded_found = [col for col in df.columns if col in EXCLUDE_COLS or col.lower() in EXCLUDE_COLS
                  or col.lower().startswith('number_') or col.lower().startswith('num_')]
if excluded_found:
    print(f"\nExcluded {len(excluded_found)} non-test columns:")
    for ec in excluded_found:
        print(f"  - {ec}")

# ============================================================================
# STEP 3: QUICK DIAGNOSTIC — bottleneck tests
# ============================================================================
print(f"\n{'=' * 100}")
print("QUICK DIAGNOSTIC: Bottleneck Tests")
print(f"{'=' * 100}")

bottleneck_tests = ['no_radicals', 'non-aromatic_ring_non-flatness', 'internal_steric_clash']
for test in bottleneck_tests:
    if test in df.columns:
        print(f"\n{test} pass rate by method:")
        print(df.groupby('docking_method')[test].mean().round(3) * 100)

# ============================================================================
# STEP 4: SURVIVAL SUMMARY
# ============================================================================
df_passed = df[df["all_passed"]].copy()

print(f"\n{'=' * 100}")
print(f"POSES THAT PASSED ALL {len(test_cols)} POSEBUSTERS TESTS: "
      f"{len(df_passed)} / {len(df)} ({100*len(df_passed)/len(df):.1f}%)")
print(f"{'=' * 100}")

print(f"\nBottleneck tests (lowest pass rates):")
pass_rates = {tc: df[tc].sum() / len(df) * 100 for tc in test_cols}
for tc, rate in sorted(pass_rates.items(), key=lambda x: x[1]):
    if rate < 99.0:
        print(f"  {tc}: {rate:.1f}%")

# Detect method column
method_col = None
for candidate in ['docking_method', 'method', 'Method', 'DOCKING_METHOD']:
    if candidate in df.columns:
        method_col = candidate
        break

if method_col:
    methods = sorted(df[method_col].unique())
    METHOD_COLORS = {m: COLOR_PALETTE[i % len(COLOR_PALETTE)] for i, m in enumerate(methods)}

    print(f"\nSurvival by method:")
    for m in methods:
        g = df[df[method_col] == m]
        n_pass = g['all_passed'].sum()
        print(f"  {m}: {n_pass}/{len(g)} passed ({n_pass/len(g)*100:.1f}%)")
else:
    methods = ['all']
    METHOD_COLORS = {'all': COLOR_PALETTE[0]}
    print("WARNING: No method column found. All poses treated as one group.")

# ============================================================================
# STEP 5: COPY PROVED POSES (only if --wd provided)
# ============================================================================
file_col = None
for candidate in ['pose_file', 'file_path', 'filepath', 'sdf_file', 'sdf_path', 'path',
                   'file', 'mol_pred', 'File', 'FILE_PATH']:
    if candidate in df_passed.columns:
        file_col = candidate
        break
if file_col is None:
    for col in df_passed.columns:
        if 'path' in col.lower() or 'file' in col.lower():
            file_col = col
            break

if wd and file_col and method_col and not df_passed.empty:
    wd_path = Path(wd)

    # Infer mode from results directory name or default to "results"
    mode_name = results_dir.name if results_dir.name in ("dock", "mol") else "results"
    copy_base = wd_path / "posebuster_proved" / mode_name

    # Detect method for passed poses
    if method_col:
        df_passed["_method"] = df_passed[method_col]
    else:
        df_passed["_method"] = "unknown"

    # Build folder mapping from unique methods
    unique_methods = df_passed["_method"].unique()
    folders = {m: copy_base / m for m in unique_methods}
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)

    copied_count = {m: 0 for m in unique_methods}
    skipped_not_found = 0

    print(f"\n{'=' * 100}")
    print("COPYING POSEBUSTERS-PROVED POSES...")
    print(f"{'=' * 100}")

    for _, row in df_passed.iterrows():
        src_path = Path(str(row[file_col]))
        method = row["_method"]

        protein_name = str(row["protein"]) if "protein" in row.index and pd.notna(row.get("protein")) else ""
        ligand_name = str(row["ligand"]) if "ligand" in row.index and pd.notna(row.get("ligand")) else ""

        if not src_path.is_absolute():
            src_path = wd_path / src_path

        if not src_path.exists():
            alt_paths = [
                results_dir / src_path.name,
                results_dir / "converted_pdbqt" / src_path.name,
                wd_path / "posebusters_results" / "converted_pdbqt" / src_path.name,
            ]
            found = False
            for alt in alt_paths:
                if alt.exists():
                    src_path = alt
                    found = True
                    break
            if not found:
                skipped_not_found += 1
                continue

        if method not in folders:
            continue

        dest_dir = folders[method]

        if protein_name and ligand_name and "__" not in src_path.name:
            dest_name = f"{protein_name}__{ligand_name}__{src_path.name}"
        else:
            dest_name = src_path.name

        dest_file = dest_dir / dest_name
        if dest_file.exists():
            stem = dest_file.stem
            suffix = dest_file.suffix
            counter = 1
            while dest_file.exists():
                dest_file = dest_dir / f"{stem}_dup{counter}{suffix}"
                counter += 1

        shutil.copy2(str(src_path), str(dest_file))
        copied_count[method] = copied_count.get(method, 0) + 1

    print(f"\n  Output directory: {copy_base}")
    print(f"\n  Files copied per method:")
    total_copied = 0
    for m, count in copied_count.items():
        if count > 0:
            print(f"    {m:>25}: {count:>4} files  ->  {folders[m]}")
            total_copied += count
    print(f"\n  Total files copied: {total_copied}")
    if skipped_not_found > 0:
        print(f"  WARNING: {skipped_not_found} files skipped (source not found)")

    # Save manifest
    manifest_cols = [file_col, "_method"]
    for extra_col in ["protein", "ligand"]:
        if extra_col in df_passed.columns:
            manifest_cols.append(extra_col)
    manifest_cols += test_cols
    manifest = df_passed[manifest_cols].copy()
    manifest.rename(columns={"_method": "docking_method"}, inplace=True)
    manifest_path = copy_base / "posebuster_proved_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"  Manifest saved: {manifest_path}")
elif not wd:
    print("\n  Skipping pose copy step (no --wd provided).")


# ============================================================================
# GRAPHICAL REPRESENTATIONS
# ============================================================================
print(f"\n{'=' * 100}")
print("GENERATING PLOTS")
print(f"{'=' * 100}")

plot_dir = output_dir
if not method_col:
    print("  WARNING: No method column — skipping per-method plots.")
    sys.exit(0)

# ============================================================================
# FIGURE 1: Overall pass rate heatmap (methods x tests)
# ============================================================================
variable_tests = [t for t in test_cols if df[t].mean() < 1.0]
trivial_tests = [t for t in test_cols if t not in variable_tests]

pass_rate_df = df.groupby(method_col)[test_cols].mean() * 100
pass_rates_var = pass_rate_df[variable_tests] if variable_tests else pass_rate_df

fig1, ax1 = plt.subplots(figsize=(max(10, len(test_cols) * 0.8), max(4, len(methods) * 0.8)))
cmap = LinearSegmentedColormap.from_list('ryg', ['#e74c3c', '#f39c12', '#27ae60'])

im = ax1.imshow(pass_rates_var.values, cmap=cmap, aspect='auto', vmin=50, vmax=100)
ax1.set_yticks(range(len(pass_rates_var.index)))
ax1.set_yticklabels(pass_rates_var.index, fontsize=11, fontweight='bold')
ax1.set_xticks(range(len(pass_rates_var.columns)))
xlabels = [TEST_DISPLAY.get(c, c.replace('_', ' ').title()) for c in pass_rates_var.columns]
ax1.set_xticklabels(xlabels, rotation=45, ha='right', fontsize=10)

for i in range(pass_rates_var.shape[0]):
    for j in range(pass_rates_var.shape[1]):
        val = pass_rates_var.values[i, j]
        color = 'white' if val < 75 else 'black'
        ax1.text(j, i, f'{val:.1f}%', ha='center', va='center', fontsize=10,
                fontweight='bold', color=color)

cbar = plt.colorbar(im, ax=ax1, shrink=0.8, pad=0.02)
cbar.set_label('Pass Rate (%)', fontsize=11)

if trivial_tests:
    trivial_str = ', '.join(TEST_DISPLAY.get(t, t) for t in trivial_tests)
    ax1.set_xlabel(f'Tests at 100% (not shown): {trivial_str}', fontsize=8, style='italic')

ax1.set_title('PoseBusters Test Pass Rates by Docking Method', fontsize=14, fontweight='bold', pad=12)
fig1.tight_layout()
fig1.savefig(plot_dir / "pb_passrate_heatmap.png", dpi=200, bbox_inches='tight')
plt.show()

# ============================================================================
# FIGURE 2: Stacked bar — pass / fail counts per method
# ============================================================================
fig2, ax2 = plt.subplots(figsize=(max(7, len(methods) * 2), 5))

counts = df.groupby(method_col)['all_passed'].value_counts().unstack(fill_value=0)
if True not in counts.columns:
    counts[True] = 0
if False not in counts.columns:
    counts[False] = 0
counts = counts.rename(columns={True: 'Passed All', False: 'Failed >=1'})
counts = counts[['Passed All', 'Failed >=1']]

bars_pass = ax2.bar(range(len(counts)), counts['Passed All'],
                    color=[METHOD_COLORS.get(m, '#888') for m in counts.index],
                    edgecolor='white', linewidth=1.5, label='Passed All Tests')
bars_fail = ax2.bar(range(len(counts)), counts['Failed >=1'],
                    bottom=counts['Passed All'],
                    color=[METHOD_COLORS.get(m, '#888') for m in counts.index],
                    alpha=0.3, edgecolor='white', linewidth=1.5, hatch='///',
                    label='Failed >=1 Test')

for i, (m, row) in enumerate(counts.iterrows()):
    total = row.sum()
    pct = row['Passed All'] / total * 100
    ax2.text(i, total + 5, f'{int(row["Passed All"])}/{int(total)}\n({pct:.1f}%)',
             ha='center', va='bottom', fontsize=11, fontweight='bold')

ax2.set_xticks(range(len(counts)))
ax2.set_xticklabels(counts.index, fontsize=12, fontweight='bold', rotation=45, ha='right')
ax2.set_ylabel('Number of Poses', fontsize=12)
ax2.set_title('PoseBusters Validation Summary by Docking Method', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
ax2.set_ylim(0, counts.sum(axis=1).max() * 1.25)
ax2.grid(axis='y', alpha=0.3)
fig2.tight_layout()
fig2.savefig(plot_dir / "pb_pass_fail_bars.png", dpi=200, bbox_inches='tight')
plt.show()

# ============================================================================
# FIGURE 3: Per-test failure rate comparison (grouped bar chart)
# ============================================================================
if variable_tests:
    fail_rates = (1 - df.groupby(method_col)[variable_tests].mean()) * 100

    fig3, ax3 = plt.subplots(figsize=(max(10, len(variable_tests) * 1.5), 5))
    x = np.arange(len(variable_tests))
    width = 0.8 / len(methods)

    for i, m in enumerate(methods):
        vals = fail_rates.loc[m].values
        ax3.bar(x + i * width - 0.4 + width / 2, vals, width,
                label=m, color=METHOD_COLORS.get(m, '#888'), edgecolor='white')
        for j, v in enumerate(vals):
            if v > 2:
                ax3.text(x[j] + i * width - 0.4 + width / 2, v + 0.3,
                         f'{v:.1f}%', ha='center', va='bottom', fontsize=7,
                         fontweight='bold', rotation=90)

    ax3.set_xticks(x)
    xlabels3 = [TEST_DISPLAY.get(c, c.replace('_', ' ').title()) for c in variable_tests]
    ax3.set_xticklabels(xlabels3, rotation=40, ha='right', fontsize=10)
    ax3.set_ylabel('Failure Rate (%)', fontsize=12)
    ax3.set_title('PoseBusters Failure Rates by Test and Method', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=10, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
    ax3.grid(axis='y', alpha=0.3)
    ax3.set_ylim(0, fail_rates.values.max() * 1.35 + 1)
    fig3.tight_layout()
    fig3.savefig(plot_dir / "pb_failure_rates.png", dpi=200, bbox_inches='tight')
    plt.show()

# ============================================================================
# FIGURE 4: Per protein-ligand pass rate (faceted by method)
# ============================================================================
if 'protein' in df.columns and 'ligand' in df.columns:
    # Compute pass rate per combo per method
    combo_method_pass = df.groupby([method_col, 'protein', 'ligand'])['all_passed'].mean() * 100
    combo_method_pass = combo_method_pass.reset_index()
    combo_method_pass.columns = ['method', 'protein', 'ligand', 'pass_rate']
    combo_method_pass['combo'] = combo_method_pass['protein'] + ' | ' + combo_method_pass['ligand']

    # Overall average pass rate per combo (for sorting)
    combo_avg = combo_method_pass.groupby('combo')['pass_rate'].mean().sort_values(ascending=True)
    sorted_combos = list(combo_avg.index)

    n_combos = len(sorted_combos)
    n_methods_fig4 = len(methods)
    bar_height = 0.8 / n_methods_fig4
    y_pos = np.arange(n_combos)

    fig4, ax4 = plt.subplots(figsize=(10, max(6, n_combos * 0.4)))

    for i, m in enumerate(methods):
        sub = combo_method_pass[combo_method_pass['method'] == m]
        # Map each combo to its sorted position
        vals = []
        for combo in sorted_combos:
            match = sub[sub['combo'] == combo]
            vals.append(match['pass_rate'].values[0] if len(match) > 0 else 0)
        ax4.barh(y_pos + i * bar_height - 0.4 + bar_height / 2, vals, bar_height,
                 label=m, color=METHOD_COLORS.get(m, '#888'), edgecolor='white', linewidth=0.5)

    ax4.set_yticks(y_pos)
    ax4.set_yticklabels(sorted_combos, fontsize=8)
    ax4.set_xlabel('Pass Rate (%)', fontsize=12)
    ax4.set_title('PoseBusters Pass Rate by Protein-Ligand Combination', fontsize=14, fontweight='bold')
    ax4.set_xlim(0, 110)
    ax4.axvline(75, color='gray', linestyle='--', alpha=0.5)
    ax4.grid(axis='x', alpha=0.3)
    ax4.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
    fig4.tight_layout()
    fig4.savefig(plot_dir / "pb_pass_by_combo.png", dpi=200, bbox_inches='tight')
    plt.show()

# ============================================================================
# FIGURE 5: Number of tests failed distribution (violin + strip)
# ============================================================================
df['n_failed'] = df[test_cols].apply(lambda row: (~row).sum(), axis=1)

fig5, ax5 = plt.subplots(figsize=(max(8, len(methods) * 2.5), 5))

for i, m in enumerate(methods):
    vals = df[df[method_col] == m]['n_failed']
    parts = ax5.violinplot([vals], positions=[i], showmedians=True, widths=0.7)
    for pc in parts['bodies']:
        pc.set_facecolor(METHOD_COLORS.get(m, '#888'))
        pc.set_alpha(0.4)
    for key in ['cbars', 'cmins', 'cmaxes', 'cmedians']:
        parts[key].set_color(METHOD_COLORS.get(m, '#888'))
    jitter = np.random.normal(0, 0.08, len(vals))
    ax5.scatter(np.full(len(vals), i) + jitter, vals, s=10, alpha=0.3,
                color=METHOD_COLORS.get(m, '#888'))

ax5.set_xticks(range(len(methods)))
ax5.set_xticklabels(methods, fontsize=12, fontweight='bold')
ax5.set_ylabel('Number of Tests Failed', fontsize=12)
ax5.set_title('Distribution of Failed Tests per Pose', fontsize=14, fontweight='bold')
ax5.set_ylim(-0.5, df['n_failed'].max() + 1)
ax5.grid(axis='y', alpha=0.3)

for i, m in enumerate(methods):
    med = df[df[method_col] == m]['n_failed'].median()
    ax5.text(i, med + 0.3, f'median={med:.0f}', ha='center', fontsize=9, fontweight='bold')

fig5.tight_layout()
fig5.savefig(plot_dir / "pb_nfailed_violin.png", dpi=200, bbox_inches='tight')
plt.show()

# ============================================================================
# FIGURE 6: % of poses that FAIL per method per test (all tests)
# ============================================================================
fail_rates_all = (1 - df.groupby(method_col)[test_cols].mean()) * 100

fig6, ax6 = plt.subplots(figsize=(max(14, len(test_cols) * 1.2), 6))
x6 = np.arange(len(test_cols))
w6 = 0.8 / len(methods)

for i, m in enumerate(methods):
    vals = fail_rates_all.loc[m].values
    ax6.bar(x6 + i * w6 - 0.4 + w6 / 2, vals, w6,
            label=m, color=METHOD_COLORS[m], edgecolor='white', linewidth=0.5)

ax6.set_xticks(x6)
xlabels6 = [TEST_DISPLAY.get(c, c.replace('_', ' ').title()) for c in test_cols]
ax6.set_xticklabels(xlabels6, rotation=50, ha='right', fontsize=9)
ax6.set_ylabel('Failure Rate (%)', fontsize=12)
ax6.set_title('Percentage of Poses Failing Each PoseBusters Test - All Methods',
              fontsize=14, fontweight='bold')
ax6.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
ax6.grid(axis='y', alpha=0.3)
ax6.set_ylim(0, fail_rates_all.values.max() * 1.15 + 1)
fig6.tight_layout()
fig6.savefig(plot_dir / "pb_fail_pct_per_method_all_tests.png", dpi=200, bbox_inches='tight')
plt.show()

# ============================================================================
# FIGURE 7: Only tests where at least one method has failures
# ============================================================================
failing_tests = [t for t in test_cols if fail_rates_all[t].max() > 0]

if failing_tests:
    fail_rates_filt = fail_rates_all[failing_tests]

    fig7, ax7 = plt.subplots(figsize=(max(10, len(failing_tests) * 1.5), 6))
    x7 = np.arange(len(failing_tests))
    w7 = 0.8 / len(methods)

    for i, m in enumerate(methods):
        vals = fail_rates_filt.loc[m].values
        ax7.bar(x7 + i * w7 - 0.4 + w7 / 2, vals, w7,
                label=m, color=METHOD_COLORS[m], edgecolor='white', linewidth=0.5)
        for j, v in enumerate(vals):
            if v > 1:
                ax7.text(x7[j] + i * w7 - 0.4 + w7 / 2, v + 0.3,
                         f'{v:.1f}%', ha='center', va='bottom', fontsize=7,
                         fontweight='bold', rotation=90)

    ax7.set_xticks(x7)
    xlabels7 = [TEST_DISPLAY.get(c, c.replace('_', ' ').title()) for c in failing_tests]
    ax7.set_xticklabels(xlabels7, rotation=45, ha='right', fontsize=10)
    ax7.set_ylabel('Failure Rate (%)', fontsize=12)
    ax7.set_title('PoseBusters Tests With Failures - Method Comparison',
                  fontsize=14, fontweight='bold')
    ax7.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1), borderaxespad=0)
    ax7.grid(axis='y', alpha=0.3)
    ax7.set_ylim(0, fail_rates_filt.values.max() * 1.35 + 1)
    fig7.tight_layout()
    fig7.savefig(plot_dir / "pb_fail_pct_failing_tests_only.png", dpi=200, bbox_inches='tight')
    plt.show()

    print(f"\n  Tests with failures ({len(failing_tests)} / {len(test_cols)}):")

    # Build tabulated output
    col_width = max(len(m) for m in methods) + 2
    test_width = max(len(TEST_DISPLAY.get(ft, ft)) for ft in failing_tests) + 2
    header = f"    {'Test':<{test_width}}" + "".join(f"{m:>{col_width}}" for m in methods)
    print(header)
    print("    " + "-" * (test_width + col_width * len(methods)))
    for ft in failing_tests:
        display_name = TEST_DISPLAY.get(ft, ft)
        row = f"    {display_name:<{test_width}}"
        for m in methods:
            row += f"{fail_rates_all.loc[m, ft]:>{col_width - 1}.1f}%"
        print(row)
else:
    print("  All poses pass all tests - no failing-tests chart needed.")

# ============================================================================
# SUMMARY
# ============================================================================
print(f"\n{'=' * 100}")
print("ANALYSIS COMPLETE")
print(f"{'=' * 100}")
print(f"\nAll figures saved to: {plot_dir}")
print("  pb_passrate_heatmap.png           - Test pass rate heatmap (methods x tests)")
print("  pb_pass_fail_bars.png             - Pass/fail stacked bars per method")
print("  pb_failure_rates.png              - Per-test failure rates (grouped bars)")
print("  pb_pass_by_combo.png              - Pass rate by protein-ligand combination")
print("  pb_nfailed_violin.png             - Distribution of #tests failed per pose")
print("  pb_fail_pct_per_method_all_tests.png  - % failure per method per test (all)")
print("  pb_fail_pct_failing_tests_only.png    - % failure for tests with failures only")
