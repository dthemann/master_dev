#!/usr/bin/env python3
"""
Compare AutoDock Vina docking results produced with two receptor preparation
pipelines: Meeko vs. MGL Tools.

Sections
--------
1. Affinity comparison   – best_affinity_kcal per (protein, ligand) pair
2. Pose RMSD             – heavy-atom RMSD between matched best poses
3. PoseBusters checks    – per-test and overall pass-rate comparison

All figures are saved to an output directory and also shown interactively.
"""
from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# ── Optional: RDKit for RMSD calculation ─────────────────────────────────────
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign, rdmolops
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================================
# CONFIGURATION — edit these paths if your layout differs
# ============================================================================

MEEKO_DOCK_LOG   = Path("Dockings/vina_results/Test_Set_meeko/docking_log_Test_Set_meeko.csv")
MGL_DOCK_LOG     = Path("Dockings/vina_results/Test_Set_mgl_tools/docking_log_Test_Set_mgl_tools.csv")

MEEKO_POSES_DIR  = Path("Dockings/vina_results/Test_Set_meeko/docking")
MGL_POSES_DIR    = Path("Dockings/vina_results/Test_Set_mgl_tools/docking")

MEEKO_PB_CSV     = Path("posebusters_results/autodock_test_set_meeko/posebusters_results.csv")
MGL_PB_CSV       = Path("posebusters_results/autodock_test_set_mgl_tools/posebusters_results.csv")

OUTPUT_DIR       = Path("comparison_meeko_vs_mgltools")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Regex to strip the converter suffix from a protein name so that both
# pipelines map to the same "base" protein identity.
_TOOL_SUFFIX_RE = re.compile(
    r"_(meeko(_meeko)?(_retry\d+)?|mgl_tools(_mgl_tools)?)$"
)


def base_protein(name: str) -> str:
    """Return the protein name without the converter suffix."""
    return _TOOL_SUFFIX_RE.sub("", name)


# ============================================================================
# 1. LOAD & MERGE DOCKING LOGS
# ============================================================================

def load_docking_logs() -> pd.DataFrame:
    meeko = pd.read_csv(MEEKO_DOCK_LOG)
    mgl = pd.read_csv(MGL_DOCK_LOG)

    meeko["pipeline"] = "meeko"
    mgl["pipeline"] = "mgl_tools"

    meeko["base_protein"] = meeko["protein_name"].apply(base_protein)
    mgl["base_protein"] = mgl["protein_name"].apply(base_protein)

    combined = pd.concat([meeko, mgl], ignore_index=True)
    return combined


# ============================================================================
# 2. AFFINITY COMPARISON
# ============================================================================

def compare_affinities(df: pd.DataFrame) -> pd.DataFrame:
    """Pivot so each row is (base_protein, ligand) with columns for each pipeline's affinity."""
    success = df[df["status"] == "success"].copy()
    pivot = success.pivot_table(
        index=["base_protein", "ligand_name"],
        columns="pipeline",
        values="best_affinity_kcal",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None

    # Only keep rows present in both pipelines
    paired = pivot.dropna(subset=["meeko", "mgl_tools"])
    return paired


def plot_affinity_scatter(paired: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(8, 8))

    proteins = paired["base_protein"].unique()
    cmap = plt.colormaps.get_cmap("tab10").resampled(max(len(proteins), 1))
    for i, prot in enumerate(sorted(proteins)):
        sub = paired[paired["base_protein"] == prot]
        ax.scatter(sub["mgl_tools"], sub["meeko"], label=prot, color=cmap(i), s=40, alpha=0.7)

    lo = min(paired[["meeko", "mgl_tools"]].min())
    hi = max(paired[["meeko", "mgl_tools"]].max())
    margin = 0.5
    ax.plot([lo - margin, hi + margin], [lo - margin, hi + margin], "k--", alpha=0.4, lw=1)

    ax.set_xlabel("MGL Tools – Best Affinity (kcal/mol)", fontsize=12)
    ax.set_ylabel("Meeko – Best Affinity (kcal/mol)", fontsize=12)
    ax.set_title("Best Docking Affinity: Meeko vs. MGL Tools", fontsize=14)
    ax.legend(fontsize=8, title="Protein", loc="lower right")
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "affinity_scatter.png", dpi=200)
    plt.close(fig)


def plot_affinity_difference(paired: pd.DataFrame) -> None:
    paired = paired.copy()
    paired["diff"] = paired["meeko"] - paired["mgl_tools"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Histogram of differences
    ax = axes[0]
    ax.hist(paired["diff"], bins=30, edgecolor="black", alpha=0.7, color="steelblue")
    ax.axvline(0, color="red", ls="--", lw=1.2)
    mean_diff = paired["diff"].mean()
    ax.axvline(mean_diff, color="orange", ls="-", lw=1.5, label=f"Mean = {mean_diff:.3f}")
    ax.set_xlabel("ΔAffinity (Meeko − MGL Tools) kcal/mol")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Affinity Differences")
    ax.legend()

    # Bland–Altman style plot
    ax = axes[1]
    paired["mean_aff"] = (paired["meeko"] + paired["mgl_tools"]) / 2
    ax.scatter(paired["mean_aff"], paired["diff"], s=25, alpha=0.6, color="steelblue")
    ax.axhline(0, color="red", ls="--", lw=1)
    ax.axhline(mean_diff, color="orange", ls="-", lw=1.2)
    sd = paired["diff"].std()
    ax.axhline(mean_diff + 1.96 * sd, color="grey", ls=":", lw=1)
    ax.axhline(mean_diff - 1.96 * sd, color="grey", ls=":", lw=1)
    ax.set_xlabel("Mean Affinity (kcal/mol)")
    ax.set_ylabel("ΔAffinity (Meeko − MGL Tools)")
    ax.set_title("Bland–Altman Plot")

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "affinity_difference.png", dpi=200)
    plt.close(fig)


# ============================================================================
# 3. POSE RMSD (heavy-atom RMSD via PDBQT coordinate comparison)
# ============================================================================

def _extract_heavy_coords_from_pdbqt(path: Path, model: int = 1) -> np.ndarray:
    """Extract heavy-atom (x, y, z) from a specific MODEL in a multi-model PDBQT."""
    coords: List[Tuple[float, float, float]] = []
    current_model = 0
    in_target = False

    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                current_model = int(line.split()[1])
                in_target = (current_model == model)
                continue
            if line.startswith("ENDMDL") and in_target:
                break
            if not in_target:
                continue
            if line.startswith(("ATOM", "HETATM")):
                element = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0]
                if element.upper() == "H":
                    continue
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))

    # Fallback: single-model file
    if not coords and model == 1:
        with open(path) as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    element = line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0]
                    if element.upper() == "H":
                        continue
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))

    return np.array(coords)


def compute_rmsd(coords_a: np.ndarray, coords_b: np.ndarray) -> float | None:
    """Compute RMSD between two coordinate arrays of equal length."""
    if coords_a.shape != coords_b.shape or len(coords_a) == 0:
        return None
    diff = coords_a - coords_b
    return float(np.sqrt((diff ** 2).sum(axis=1).mean()))


def compute_pose_rmsds(paired_aff: pd.DataFrame) -> pd.DataFrame:
    """For each matched (protein, ligand), compute heavy-atom RMSD between best poses (model 1)."""
    rows = []
    for _, row in paired_aff.iterrows():
        bp = row["base_protein"]
        lig = row["ligand_name"]

        # Find matching PDBQT files
        meeko_pattern = f"*{lig}*_vina_out.pdbqt"
        mgl_pattern = f"*{lig}*_vina_out.pdbqt"

        meeko_hits = sorted(MEEKO_POSES_DIR.glob(meeko_pattern))
        mgl_hits = sorted(MGL_POSES_DIR.glob(mgl_pattern))

        # Filter to the correct base protein
        meeko_file = None
        mgl_file = None
        for f in meeko_hits:
            if base_protein(f.stem.split("__")[0]) == bp:
                meeko_file = f
                break
        for f in mgl_hits:
            if base_protein(f.stem.split("__")[0]) == bp:
                mgl_file = f
                break

        if meeko_file is None or mgl_file is None:
            continue

        coords_m = _extract_heavy_coords_from_pdbqt(meeko_file, model=1)
        coords_g = _extract_heavy_coords_from_pdbqt(mgl_file, model=1)

        rmsd = compute_rmsd(coords_m, coords_g)
        rows.append({
            "base_protein": bp,
            "ligand_name": lig,
            "rmsd": rmsd,
            "meeko_affinity": row["meeko"],
            "mgl_affinity": row["mgl_tools"],
            "n_atoms_meeko": len(coords_m),
            "n_atoms_mgl": len(coords_g),
        })

    return pd.DataFrame(rows)


def plot_rmsd(rmsd_df: pd.DataFrame) -> None:
    valid = rmsd_df.dropna(subset=["rmsd"])
    if valid.empty:
        print("  No valid RMSD values to plot.")
        return

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Histogram
    ax = axes[0]
    ax.hist(valid["rmsd"], bins=30, edgecolor="black", alpha=0.7, color="coral")
    ax.axvline(valid["rmsd"].median(), color="navy", ls="--", lw=1.5, label=f'Median = {valid["rmsd"].median():.2f} Å')
    ax.set_xlabel("RMSD (Å)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of Pose RMSD\n(Meeko vs. MGL Tools, best pose)")
    ax.legend()

    # RMSD vs affinity difference
    ax = axes[1]
    valid = valid.copy()
    valid["aff_diff"] = valid["meeko_affinity"] - valid["mgl_affinity"]
    ax.scatter(valid["aff_diff"], valid["rmsd"], s=30, alpha=0.6, color="teal")
    ax.set_xlabel("ΔAffinity (Meeko − MGL Tools) kcal/mol")
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Pose RMSD vs. Affinity Difference")
    ax.axvline(0, color="red", ls="--", lw=0.8, alpha=0.5)

    # RMSD per protein (box plot)
    ax = axes[2]
    proteins = sorted(valid["base_protein"].unique())
    data_per_prot = [valid[valid["base_protein"] == p]["rmsd"].values for p in proteins]
    short_names = [p.replace("Orai1WT-", "").replace("_cleaned", "") for p in proteins]
    bp = ax.boxplot(data_per_prot, tick_labels=short_names, patch_artist=True)
    colors = plt.cm.Set2(np.linspace(0, 1, len(proteins)))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
    ax.set_ylabel("RMSD (Å)")
    ax.set_title("Pose RMSD by Protein")
    ax.tick_params(axis="x", rotation=30)

    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pose_rmsd.png", dpi=200)
    plt.close(fig)


# ============================================================================
# 4. POSEBUSTERS COMPARISON
# ============================================================================

# Boolean test columns of interest (physical validity checks)
PB_TEST_COLS = [
    "mol_pred_loaded",
    "sanitization",
    "all_atoms_connected",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
]


def load_posebusters() -> pd.DataFrame:
    meeko = pd.read_csv(MEEKO_PB_CSV)
    mgl = pd.read_csv(MGL_PB_CSV)

    meeko["pipeline"] = "meeko"
    mgl["pipeline"] = "mgl_tools"

    meeko["base_protein"] = meeko["protein"].apply(base_protein)
    mgl["base_protein"] = mgl["protein"].apply(base_protein)

    combined = pd.concat([meeko, mgl], ignore_index=True)

    # Normalise boolean columns
    for col in PB_TEST_COLS + ["all_passed", "critical_passed"]:
        if col in combined.columns:
            combined[col] = combined[col].map(
                {True: True, False: False, "True": True, "False": False, 1: True, 0: False}
            ).astype(bool)

    return combined


def plot_pb_overall_pass_rate(pb: pd.DataFrame) -> None:
    """Bar chart: overall pass rate (all_passed & critical_passed) per pipeline."""
    metrics = []
    for col in ["all_passed", "critical_passed"]:
        if col not in pb.columns:
            continue
        for pipe in ["meeko", "mgl_tools"]:
            sub = pb[pb["pipeline"] == pipe]
            rate = sub[col].mean() * 100
            metrics.append({"metric": col, "pipeline": pipe, "pass_rate": rate, "n": len(sub)})
    if not metrics:
        return

    mdf = pd.DataFrame(metrics)

    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(mdf["metric"].nunique())
    width = 0.35
    metric_names = mdf["metric"].unique()

    for i, pipe in enumerate(["meeko", "mgl_tools"]):
        vals = [mdf[(mdf["metric"] == m) & (mdf["pipeline"] == pipe)]["pass_rate"].values[0] for m in metric_names]
        ns =   [mdf[(mdf["metric"] == m) & (mdf["pipeline"] == pipe)]["n"].values[0] for m in metric_names]
        bars = ax.bar(x + i * width, vals, width, label=f"{pipe} (n={ns[0]})")
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1, f"{v:.1f}%", ha="center", fontsize=9)

    ax.set_xticks(x + width / 2)
    ax.set_xticklabels([m.replace("_", " ").title() for m in metric_names])
    ax.set_ylabel("Pass Rate (%)")
    ax.set_title("PoseBusters Overall Pass Rates")
    ax.set_ylim(0, 110)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pb_overall_pass_rate.png", dpi=200)
    plt.close(fig)


def plot_pb_per_test(pb: pd.DataFrame) -> None:
    """Grouped bar chart: pass rate per individual PB test, by pipeline."""
    available = [c for c in PB_TEST_COLS if c in pb.columns]
    if not available:
        print("  No PB test columns found.")
        return

    rates: Dict[str, Dict[str, float]] = {}
    for col in available:
        rates[col] = {}
        for pipe in ["meeko", "mgl_tools"]:
            sub = pb[pb["pipeline"] == pipe]
            rates[col][pipe] = sub[col].mean() * 100

    rdf = pd.DataFrame(rates).T
    rdf.index = [c.replace("_", " ").replace("-", " ").title() for c in rdf.index]

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(len(rdf))
    width = 0.35

    ax.barh(x - width / 2, rdf["meeko"], width, label="Meeko", color="steelblue")
    ax.barh(x + width / 2, rdf["mgl_tools"], width, label="MGL Tools", color="coral")

    ax.set_yticks(x)
    ax.set_yticklabels(rdf.index)
    ax.set_xlabel("Pass Rate (%)")
    ax.set_title("PoseBusters Per-Test Pass Rates: Meeko vs. MGL Tools")
    ax.set_xlim(0, 105)
    ax.legend(loc="lower right")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pb_per_test.png", dpi=200)
    plt.close(fig)


def plot_pb_per_protein(pb: pd.DataFrame) -> None:
    """Pass rate per protein, grouped by pipeline."""
    if "all_passed" not in pb.columns:
        return

    rates = pb.groupby(["base_protein", "pipeline"])["all_passed"].mean().unstack(fill_value=0) * 100
    short = [n.replace("Orai1WT-", "").replace("_cleaned", "") for n in rates.index]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(rates))
    width = 0.35

    if "meeko" in rates.columns:
        ax.bar(x - width / 2, rates["meeko"], width, label="Meeko", color="steelblue")
    if "mgl_tools" in rates.columns:
        ax.bar(x + width / 2, rates["mgl_tools"], width, label="MGL Tools", color="coral")

    ax.set_xticks(x)
    ax.set_xticklabels(short, rotation=30, ha="right")
    ax.set_ylabel("All-Passed Rate (%)")
    ax.set_title("PoseBusters All-Passed Rate by Protein")
    ax.set_ylim(0, 110)
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pb_per_protein.png", dpi=200)
    plt.close(fig)


def plot_pb_per_ligand_comparison(pb: pd.DataFrame) -> None:
    """Scatter: per-ligand all_passed rate in Meeko vs. MGL Tools."""
    if "all_passed" not in pb.columns:
        return

    lig_rates = pb.groupby(["ligand", "pipeline"])["all_passed"].mean().unstack(fill_value=np.nan)
    paired = lig_rates.dropna()
    if paired.empty or "meeko" not in paired.columns or "mgl_tools" not in paired.columns:
        print("  Not enough paired ligand data for scatter.")
        return

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(paired["mgl_tools"] * 100, paired["meeko"] * 100, s=30, alpha=0.6, color="teal")
    ax.plot([0, 100], [0, 100], "k--", alpha=0.4)
    ax.set_xlabel("MGL Tools – All-Passed Rate (%)")
    ax.set_ylabel("Meeko – All-Passed Rate (%)")
    ax.set_title("Per-Ligand PoseBusters Pass Rate")
    ax.set_xlim(-5, 105)
    ax.set_ylim(-5, 105)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "pb_per_ligand_scatter.png", dpi=200)
    plt.close(fig)


# ============================================================================
# 5. SUMMARY TABLE
# ============================================================================

def print_summary(paired_aff: pd.DataFrame, rmsd_df: pd.DataFrame, pb: pd.DataFrame) -> None:
    print("\n" + "=" * 80)
    print("COMPARISON SUMMARY: Meeko vs. MGL Tools")
    print("=" * 80)

    # Affinity
    diff = paired_aff["meeko"] - paired_aff["mgl_tools"]
    print(f"\n  Paired affinity comparisons:  {len(paired_aff)}")
    print(f"  Mean ΔAffinity (Meeko − MGL): {diff.mean():.3f} ± {diff.std():.3f} kcal/mol")
    print(f"  Median ΔAffinity:             {diff.median():.3f} kcal/mol")
    corr = paired_aff[["meeko", "mgl_tools"]].corr().iloc[0, 1]
    print(f"  Pearson correlation:           {corr:.4f}")

    # RMSD
    valid_rmsd = rmsd_df.dropna(subset=["rmsd"])
    print(f"\n  Pose RMSD comparisons:  {len(valid_rmsd)}")
    if not valid_rmsd.empty:
        print(f"  Mean RMSD:              {valid_rmsd['rmsd'].mean():.3f} Å")
        print(f"  Median RMSD:            {valid_rmsd['rmsd'].median():.3f} Å")
        print(f"  Max RMSD:               {valid_rmsd['rmsd'].max():.3f} Å")
        pct_low = (valid_rmsd["rmsd"] < 2.0).mean() * 100
        print(f"  Poses with RMSD < 2 Å:  {pct_low:.1f}%")
    mismatch = rmsd_df[rmsd_df["n_atoms_meeko"] != rmsd_df["n_atoms_mgl"]]
    if len(mismatch) > 0:
        print(f"  Atom-count mismatches (RMSD=N/A): {len(mismatch)}")

    # PoseBusters
    for col in ["all_passed", "critical_passed"]:
        if col not in pb.columns:
            continue
        print(f"\n  PoseBusters '{col}':")
        for pipe in ["meeko", "mgl_tools"]:
            sub = pb[pb["pipeline"] == pipe]
            rate = sub[col].mean() * 100
            print(f"    {pipe:12s}:  {rate:.1f}%  (n={len(sub)})")

    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 80)
    print("MEEKO vs. MGL TOOLS — DOCKING COMPARISON")
    print("=" * 80)

    # ── 1. Docking affinities ───────────────────────────────────────────
    print("\n── Loading docking logs ──")
    dock_df = load_docking_logs()
    paired_aff = compare_affinities(dock_df)
    print(f"  Paired (protein, ligand) comparisons: {len(paired_aff)}")

    print("\n── Plotting affinity comparison ──")
    plot_affinity_scatter(paired_aff)
    plot_affinity_difference(paired_aff)

    # ── 2. Pose RMSD ───────────────────────────────────────────────────
    print("\n── Computing pose RMSDs ──")
    rmsd_df = compute_pose_rmsds(paired_aff)
    valid = rmsd_df.dropna(subset=["rmsd"])
    print(f"  Valid RMSD comparisons: {len(valid)} / {len(rmsd_df)}")
    if not valid.empty:
        print(f"  Mean RMSD: {valid['rmsd'].mean():.3f} Å   Median: {valid['rmsd'].median():.3f} Å")
    plot_rmsd(rmsd_df)

    # ── 3. PoseBusters ─────────────────────────────────────────────────
    print("\n── Loading PoseBusters results ──")
    pb = load_posebusters()
    print(f"  Total PB rows: {len(pb)}  (meeko={len(pb[pb['pipeline']=='meeko'])}, mgl_tools={len(pb[pb['pipeline']=='mgl_tools'])})")

    print("\n── Plotting PoseBusters comparison ──")
    plot_pb_overall_pass_rate(pb)
    plot_pb_per_test(pb)
    plot_pb_per_protein(pb)
    plot_pb_per_ligand_comparison(pb)

    # ── 4. Summary ─────────────────────────────────────────────────────
    print_summary(paired_aff, rmsd_df, pb)

    # Save RMSD data
    rmsd_df.to_csv(OUTPUT_DIR / "pose_rmsd_data.csv", index=False)
    paired_aff.to_csv(OUTPUT_DIR / "paired_affinities.csv", index=False)
    print(f"Results saved to {OUTPUT_DIR.resolve()}/")


if __name__ == "__main__":
    main()
