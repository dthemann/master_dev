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
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import stats_utils as su

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
    ap.add_argument("--diffdock-variant", default="all",
                    choices=("all", "raw", "smina", "gnina"),
                    help="Restrict DiffDock to one optimizer variant before the per-tool "
                         "min-RMSD oracle (default 'all' pools raw+smina+gnina).")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    n = args.n

    met = pd.read_csv(args.per_pose_metrics, low_memory=False)
    feat = pd.read_csv(args.features)
    if args.diffdock_variant != "all":
        _target = {"raw": "diffdock", "smina": "diffdock_smina",
                   "gnina": "diffdock_gnina"}[args.diffdock_variant]
        _dd = met["method"].astype(str).str.startswith("diffdock")
        met = met[(~_dd) | (met["method"].astype(str) == _target)].copy()
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

    # ---- NEW: per-cell Mann–Whitney (worst vs best) + BH across the heatmap,
    #      Cliff's-δ bootstrap CI, and a stronger full-sample Spearman.
    #      All wrapped so a stats failure degrades to the test-free figure. ----
    star_map = {}  # (nice_property, tool) -> stars string for BH-significant cells
    try:
        cell_records = []
        for tool in cols:
            sub = (oracle[oracle["tool"] == tool].dropna(subset=["oracle_rmsd"])
                   .sort_values("oracle_rmsd"))
            best_g, worst_g = sub.head(n), sub.tail(n)
            for p in feats:
                w = worst_g[p].to_numpy(float); b = best_g[p].to_numpy(float)
                mw = su.mannwhitney_cliffs(w, b)
                if mw is None:
                    continue
                d, lo, hi = su.cliffs_delta_ci(w, b)
                cell_records.append(dict(
                    tool=tool, property=NICE.get(p, p), descriptor=p,
                    cliffs_delta=d, delta_ci_lo=lo, delta_ci_hi=hi,
                    mw_U=mw["U"], p_raw=mw["p"],
                    n_worst=mw["n_a"], n_best=mw["n_b"]))
        if cell_records:
            qs = su.bh_fdr([r["p_raw"] for r in cell_records])
            for r, q in zip(cell_records, qs):
                r["p_bh"] = float(q); r["star"] = su.p_stars(q)
                if q == q and q < 0.05:
                    star_map[(r["property"], r["tool"])] = su.p_stars(q)
        with open(args.out_dir / "receptor_difficulty_delta_heatmap_stats.json", "w") as fh:
            json.dump({"n_split": n, "unit": "one oracle-best-RMSD value per complex",
                       "family": "mannwhitney_cliffs, BH across all heatmap cells",
                       "cells": cell_records}, fh, indent=2)

        # Full-sample Spearman: descriptor vs continuous oracle RMSD over ALL
        # complexes (not just the N=60 tails), BH across the descriptors/tool.
        spear_rows = []
        for tool in cols:
            sub = oracle[oracle["tool"] == tool].dropna(subset=["oracle_rmsd"])
            recs = []
            for p in feats:
                rho, pv, nn = su.spearman(sub[p].to_numpy(float),
                                          sub["oracle_rmsd"].to_numpy(float))
                recs.append(dict(tool=tool, property=NICE.get(p, p), descriptor=p,
                                 spearman_rho=rho, spearman_p=pv, n=nn))
            for r, q in zip(recs, su.bh_fdr([r["spearman_p"] for r in recs])):
                r["spearman_p_bh"] = float(q); r["star"] = su.p_stars(q)
            spear_rows.extend(recs)
        pd.DataFrame(spear_rows).to_csv(args.out_dir / "difficulty_stats.csv", index=False)
        print("\n=== full-sample Spearman ρ (descriptor vs oracle RMSD, all complexes; BH) ===")
        print(pd.DataFrame(spear_rows).pivot(index="property", columns="tool",
              values="spearman_rho").reindex(dtbl.index).round(2).to_string())
    except Exception as e:
        print(f"[warn] difficulty stats skipped ({e}); drawing test-free heatmap")
        star_map = {}

    fig, ax = plt.subplots(figsize=(1.6 * len(cols) + 4.0, 0.42 * len(dtbl) + 1.6))
    sns.heatmap(dtbl[cols], annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-0.7, vmax=0.7, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Cliff's δ  (+ = larger in hard-to-dock complexes)"}, ax=ax)
    for i, prop in enumerate(dtbl.index):          # overlay BH-significant stars
        for j, tool in enumerate(cols):
            s = star_map.get((prop, tool))
            if s:
                ax.text(j + 0.5, i + 0.16, s, ha="center", va="center",
                        fontsize=9, fontweight="bold", color="black")
    ax.set_xlabel("Docking tool"); ax.set_ylabel("Receptor / pocket property")
    ax.set_title(f"What receptor/pocket properties separate the best-{n} from worst-{n}\n"
                 "docked complexes, per tool (oracle RMSD-to-crystal)\n"
                 "* / ** / *** = Mann–Whitney worst-vs-best, BH-corrected q<.05/.01/.001",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out_dir / "receptor_difficulty_delta_heatmap.png", dpi=160)
    plt.close(fig)
    print(f"\nWrote CSV + heatmap to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
