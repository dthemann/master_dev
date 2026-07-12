"""What receptor / pocket properties separate well- vs badly-docked complexes, per tool.

Companion to ligand_docking_difficulty.py — same best-N vs worst-N split (by each tool's
oracle best RMSD-to-crystal per complex), but comparing RECEPTOR and binding-POCKET
descriptors instead of ligand ones. This is the receptor-side explanation of the
rankings — and, in particular, of a tool like DiffDock whose difficulty is *not*
explained by ligand chemistry.

Receptor/pocket features come from benchmark_receptor_properties.csv (written by
PoseBusters_DataSet_Analysis.ipynb): protein size (heavy atoms, residues, chains), the
number of hetero/cofactor atoms near the structure, how many candidate pockets fpocket /
p2rank found (site ambiguity), and the top/best pocket volume, druggability and scores.

Effect size = Cliff's delta (worst-N vs best-N): +1 → the property is always larger in
the harder-to-dock complexes, -1 → always smaller, ~0 → no separation.

Usage:
    python Scripts/Analysis/receptor_docking_difficulty.py            # defaults (benchmark)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

TOOLS = ["autodock", "diffdock", "equibind"]
FEATURES = ["prot_heavy_atoms", "n_residues", "n_chains", "prot_hetatm",
            "n_cofactor_types", "n_fpocket", "n_p2rank", "top1_volume", "max_volume",
            "top1_druggability", "max_druggability", "top1_fpocket_score",
            "top_p2rank_score"]
NICE = {"prot_heavy_atoms": "Protein heavy atoms", "n_residues": "Residues",
        "n_chains": "Chains", "prot_hetatm": "HETATM atoms (cofactors/ions)",
        "n_cofactor_types": "Cofactor types", "n_fpocket": "# fpocket pockets",
        "n_p2rank": "# p2rank pockets", "top1_volume": "Top-pocket volume",
        "max_volume": "Max pocket volume", "top1_druggability": "Top-pocket druggability",
        "max_druggability": "Max druggability", "top1_fpocket_score": "Top fpocket score",
        "top_p2rank_score": "Top p2rank score"}


def cliffs_delta(worst, best) -> float:
    a = np.asarray(worst, float); b = np.asarray(best, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = int((a[:, None] > b[None, :]).sum())
    lt = int((a[:, None] < b[None, :]).sum())
    return (gt - lt) / (a.size * b.size)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-metrics", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "pose_comparison_report/per_pose_metrics.csv"))
    ap.add_argument("--features", type=Path,
                    default=Path("posebusters_results/benchmark/receptor_comparison/"
                                 "benchmark_receptor_properties.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/receptor_difficulty"))
    ap.add_argument("--n", type=int, default=60)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    n = args.n

    met = pd.read_csv(args.per_pose_metrics, low_memory=False)
    feat = pd.read_csv(args.features)
    met["tool"] = met["method"].astype(str).str.replace(r"_.*$", "", regex=True)
    oracle = (met.groupby(["tool", "protein"])["rmsd"].min().reset_index()
              .rename(columns={"protein": "entry", "rmsd": "oracle_rmsd"})
              .merge(feat, on="entry", how="left"))
    feats = [f for f in FEATURES if f in oracle.columns]

    delta_tbl, med_rows = {}, []
    for tool in TOOLS:
        sub = oracle[oracle["tool"] == tool].dropna(subset=["oracle_rmsd"]).sort_values("oracle_rmsd")
        if sub.empty:
            continue
        best, worst = sub.head(n), sub.tail(n)
        for p in feats:
            d = cliffs_delta(worst[p], best[p])
            delta_tbl.setdefault(NICE.get(p, p), {})[tool] = d
            med_rows.append(dict(tool=tool, property=NICE.get(p, p),
                                 well_median=best[p].median(), badly_median=worst[p].median(),
                                 cliffs_delta=d))

    cols = [t for t in TOOLS if any(r["tool"] == t for r in med_rows)]
    dtbl = pd.DataFrame(delta_tbl).T.reindex(columns=cols)
    dtbl["mean"] = dtbl.mean(axis=1)
    dtbl = dtbl.sort_values("mean", key=lambda s: s.abs(), ascending=False)
    pd.DataFrame(med_rows).to_csv(args.out_dir / "receptor_property_best_vs_worst.csv", index=False)

    print(f"=== Cliff's delta (worst-{n} vs best-{n} docked); + = larger in HARD complexes ===")
    print(dtbl.round(2).to_string())
    print("\n=== receptor/pocket medians, well vs badly docked ===")
    for tool in cols:
        t = pd.DataFrame([r for r in med_rows if r["tool"] == tool]).set_index("property")
        t = t.reindex([p for p in dtbl.index if p in t.index])
        print(f"\n-- {tool} --")
        print(t[["well_median", "badly_median", "cliffs_delta"]].round(2).to_string())

    fig, ax = plt.subplots(figsize=(1.6 * len(cols) + 4.0, 0.42 * len(dtbl) + 1.6))
    sns.heatmap(dtbl[cols], annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-0.7, vmax=0.7, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Cliff's δ  (+ = larger in hard-to-dock complexes)"}, ax=ax)
    ax.set_xlabel("Docking tool"); ax.set_ylabel("Receptor / pocket property")
    ax.set_title(f"What receptor/pocket properties separate the best-{n} from worst-{n}\n"
                 "docked complexes, per tool (oracle RMSD-to-crystal)")
    fig.tight_layout()
    fig.savefig(args.out_dir / "receptor_difficulty_delta_heatmap.png", dpi=160)
    plt.close(fig)
    print(f"\nWrote CSV + heatmap to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
