"""Visualize PoseBusters validity per docking tool & receptor-ligand pair.

Reads posebusters_filtered_results.csv and produces:
    1. Bar chart: total / valid pose counts per docking tool.
    2. Heatmap: valid-pose counts per receptor-ligand x tool.
    3. Stacked bar: top-N receptor-ligand pairs comparing valid poses per tool.
    4. CSV summaries written next to the plots.

A pose is considered "physically valid" when ALL of PoseBusters' most
critical structural checks pass (the canonical PB-Valid criterion).

Quick usage:
    # Run with defaults
    python Scripts/Analysis/posebusters_validity_report.py

    # Run with custom input/output locations
    python Scripts/Analysis/posebusters_validity_report.py \
        --csv posebusters_results/benchmark/dock/posebusters_filtered_results.csv \
        --out-dir posebusters_results/benchmark/dock/validity_report
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# PoseBusters "most critical" structural-validity checks (PB-Valid set).
# Boolean columns; True == pass.
CRITICAL_CHECKS: list[str] = [
    "mol_pred_loaded",
    "sanitization",
    "all_atoms_connected",
    "bond_lengths",
    "bond_angles",
    "internal_steric_clash",
    "aromatic_ring_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_flatness",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "volume_overlap_with_protein",
]

TOOL_ORDER = ["autodock", "diffdock", "equibind_guided"]
TOOL_LABEL = {
    "autodock": "AutoDock Vina",
    "diffdock": "DiffDock",
    "equibind_guided": "EquiBind (guided)",
}


def _to_bool(series: pd.Series) -> pd.Series:
    """Coerce mixed bool/float/string PB outputs to boolean (NaN -> False)."""
    if series.dtype == bool:
        return series
    if pd.api.types.is_numeric_dtype(series):
        return series.fillna(0).astype(bool)
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .map({"true": True, "1": True, "1.0": True,
              "false": False, "0": False, "0.0": False, "nan": False, "": False})
        .fillna(False)
        .astype(bool)
    )


def load_and_score(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path, low_memory=False)
    missing = [c for c in CRITICAL_CHECKS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing PoseBusters columns in CSV: {missing}")

    bool_df = pd.DataFrame({c: _to_bool(df[c]) for c in CRITICAL_CHECKS})
    df["pb_valid"] = bool_df.all(axis=1)
    df["pair"] = df["protein"].astype(str) + " / " + df["ligand"].astype(str)
    df["docking_method"] = df["docking_method"].astype(str).str.lower()
    return df


def per_tool_summary(df: pd.DataFrame) -> pd.DataFrame:
    grp = (
        df.groupby("docking_method")
        .agg(total_poses=("pb_valid", "size"),
             valid_poses=("pb_valid", "sum"))
        .reindex(TOOL_ORDER)
        .dropna(how="all")
    )
    grp["valid_fraction"] = grp["valid_poses"] / grp["total_poses"]
    return grp


def per_pair_tool_matrix(df: pd.DataFrame, value: str) -> pd.DataFrame:
    """Pivot table: rows=pair, cols=tool. value in {'valid', 'total', 'fraction'}."""
    if value == "valid":
        mat = df.pivot_table(index="pair", columns="docking_method",
                             values="pb_valid", aggfunc="sum", fill_value=0)
    elif value == "total":
        mat = df.pivot_table(index="pair", columns="docking_method",
                             values="pb_valid", aggfunc="size", fill_value=0)
    elif value == "fraction":
        valid = df.pivot_table(index="pair", columns="docking_method",
                               values="pb_valid", aggfunc="sum", fill_value=0)
        total = df.pivot_table(index="pair", columns="docking_method",
                               values="pb_valid", aggfunc="size", fill_value=0)
        mat = (valid / total.replace(0, np.nan)).fillna(0)
    else:
        raise ValueError(value)

    cols = [c for c in TOOL_ORDER if c in mat.columns]
    return mat[cols]


# ───────────────────────────── plots ─────────────────────────────

def plot_per_tool(summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(summary))
    w = 0.38
    ax.bar(x - w / 2, summary["total_poses"], w, label="Total poses",
           color="#bdbdbd", edgecolor="black")
    ax.bar(x + w / 2, summary["valid_poses"], w, label="PB-Valid poses",
           color="#2ca02c", edgecolor="black")

    for i, (tot, val) in enumerate(zip(summary["total_poses"], summary["valid_poses"])):
        pct = (val / tot * 100) if tot else 0
        ax.text(i + w / 2, val, f"{int(val)}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=9)
        ax.text(i - w / 2, tot, f"{int(tot)}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([TOOL_LABEL.get(t, t) for t in summary.index])
    ax.set_ylabel("Number of poses")
    ax.set_title("PoseBusters Benchmark — Pose validity per docking tool\n"
                 "(valid = passes all critical PB checks)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_heatmap(valid_mat: pd.DataFrame, out: Path, top_n: int = 60) -> None:
    mat = valid_mat.copy()
    mat["__sum"] = mat.sum(axis=1)
    mat = mat.sort_values("__sum", ascending=False).head(top_n).drop(columns="__sum")
    mat.columns = [TOOL_LABEL.get(c, c) for c in mat.columns]

    fig, ax = plt.subplots(figsize=(8, max(8, 0.22 * len(mat))))
    sns.heatmap(mat, annot=True, fmt=".0f", cmap="YlGnBu", cbar_kws={"label": "Valid poses"},
                linewidths=0.3, linecolor="white", ax=ax)
    ax.set_xlabel("Docking tool")
    ax.set_ylabel("Receptor / Ligand")
    ax.set_title(f"Valid poses per receptor-ligand pair (top {len(mat)})")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_grouped_bars(valid_mat: pd.DataFrame, out: Path, top_n: int = 30) -> None:
    mat = valid_mat.copy()
    mat["__sum"] = mat.sum(axis=1)
    mat = mat.sort_values("__sum", ascending=False).head(top_n).drop(columns="__sum")

    pairs = mat.index.tolist()
    tools = list(mat.columns)
    x = np.arange(len(pairs))
    w = 0.8 / max(len(tools), 1)
    colors = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind_guided": "#2ca02c"}

    fig, ax = plt.subplots(figsize=(max(10, 0.45 * len(pairs)), 6))
    for i, t in enumerate(tools):
        ax.bar(x + (i - (len(tools) - 1) / 2) * w, mat[t].values, w,
               label=TOOL_LABEL.get(t, t), color=colors.get(t, None), edgecolor="black")

    ax.set_xticks(x)
    ax.set_xticklabels(pairs, rotation=90, fontsize=8)
    ax.set_ylabel("Valid poses")
    ax.set_title(f"Valid poses per docking tool — top {len(pairs)} receptor-ligand pairs")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_validity_distribution(df: pd.DataFrame, out: Path) -> None:
    """Distribution of per-pair valid-pose counts by tool (boxplot)."""
    valid = df.pivot_table(index="pair", columns="docking_method",
                           values="pb_valid", aggfunc="sum", fill_value=0)
    cols = [c for c in TOOL_ORDER if c in valid.columns]
    valid = valid[cols]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.boxplot([valid[c].values for c in cols],
               labels=[TOOL_LABEL.get(c, c) for c in cols], showmeans=True)
    ax.set_ylabel("Valid poses per receptor-ligand pair")
    ax.set_title("Distribution of valid poses across receptor-ligand pairs")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_check_passrate(df: pd.DataFrame, out: Path) -> None:
    """Per-check pass rate per tool (which PB criterion is failing most)."""
    rows = []
    for tool, sub in df.groupby("docking_method"):
        for chk in CRITICAL_CHECKS:
            rows.append({"tool": tool, "check": chk,
                         "pass_rate": _to_bool(sub[chk]).mean()})
    pr = pd.DataFrame(rows).pivot(index="check", columns="tool", values="pass_rate")
    pr = pr[[c for c in TOOL_ORDER if c in pr.columns]]
    pr.columns = [TOOL_LABEL.get(c, c) for c in pr.columns]

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(pr * 100, annot=True, fmt=".1f", cmap="RdYlGn", vmin=0, vmax=100,
                cbar_kws={"label": "Pass rate (%)"}, ax=ax)
    ax.set_title("Per-check pass rate (%) per docking tool")
    ax.set_xlabel("")
    ax.set_ylabel("PoseBusters critical check")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ───────────────────────────── main ─────────────────────────────

def main() -> None:
    # How to run this file from a notebook cell:
    #   !python Scripts/Analysis/posebusters_validity_report.py
    #
    # Helpful optional flags:
    #   --csv             Path to PoseBusters filtered-results CSV.
    #   --out-dir         Folder where plots and summary CSVs are saved.
    #   --top-n-heatmap   Number of receptor-ligand pairs shown in heatmap.
    #   --top-n-bars      Number of receptor-ligand pairs shown in grouped bars.
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "posebusters_filtered_results.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("posebusters_results/benchmark/dock/validity_report"))
    ap.add_argument("--top-n-heatmap", type=int, default=60)
    ap.add_argument("--top-n-bars", type=int, default=30)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_and_score(args.csv)

    summary = per_tool_summary(df)
    valid_mat = per_pair_tool_matrix(df, "valid")
    total_mat = per_pair_tool_matrix(df, "total")
    frac_mat = per_pair_tool_matrix(df, "fraction")

    summary.to_csv(args.out_dir / "summary_per_tool.csv")
    valid_mat.to_csv(args.out_dir / "valid_poses_per_pair_tool.csv")
    total_mat.to_csv(args.out_dir / "total_poses_per_pair_tool.csv")
    frac_mat.round(3).to_csv(args.out_dir / "valid_fraction_per_pair_tool.csv")

    plot_per_tool(summary, args.out_dir / "01_per_tool_validity.png")
    plot_heatmap(valid_mat, args.out_dir / "02_valid_heatmap_top.png",
                 top_n=args.top_n_heatmap)
    plot_grouped_bars(valid_mat, args.out_dir / "03_valid_grouped_bars_top.png",
                      top_n=args.top_n_bars)
    plot_validity_distribution(df, args.out_dir / "04_valid_per_pair_distribution.png")
    plot_check_passrate(df, args.out_dir / "05_per_check_passrate.png")

    print(f"Total poses scored        : {len(df):,}")
    print(f"Total PB-valid poses      : {int(df['pb_valid'].sum()):,} "
          f"({df['pb_valid'].mean() * 100:.2f}%)")
    print(f"Receptor-ligand pairs     : {df['pair'].nunique()}")
    print()
    print("── Per-tool summary ─────────────────────────────")
    print(summary.to_string(float_format=lambda x: f"{x:.3f}"))
    print()
    print(f"Wrote CSVs and 5 figures to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
