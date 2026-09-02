"""What kinds of ligands dock well vs badly, per docking tool.

Joins per-pose docking accuracy (per_pose_metrics.csv, written by
posebusters_pose_comparison.py) to the benchmark ligand descriptor table
(ligand_protein_features.csv, from PoseBusters_DataSet_Analysis.ipynb) and asks, for
each docking tool, which ligand physicochemical properties separate the complexes it
docks well from the ones it docks badly.

Success metric = the *oracle* (best) RMSD-to-crystal a tool reaches for a complex — the
minimum RMSD over all that tool's poses and variants (e.g. DiffDock pools
diffdock/_smina/_gnina; EquiBind pools every pocket×refine×clamp variant). Per tool we
rank complexes by oracle RMSD and take the N best-docked (lowest) and N worst-docked
(highest), then compare their ligand properties.

Effect size = Cliff's delta (worst-N vs best-N): +1 → the property is always larger in
the harder-to-dock ligands, -1 → always smaller, ~0 → no separation. |δ| ≳ 0.33 is a
medium and ≳ 0.47 a large effect.

Outputs (to --out-dir):
    ligand_property_best_vs_worst.csv   — per tool × property: well/badly medians + δ
    ligand_difficulty_delta_heatmap.png — property × tool Cliff's-δ heatmap
and a printed report (per-tool success rates + the δ table).

Usage:
    python Scripts/Analysis/ligand_docking_difficulty.py            # defaults (benchmark)
    python Scripts/Analysis/ligand_docking_difficulty.py --n 40
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
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)

TOOLS = ["autodock", "diffdock", "equibind"]
CONT = ["mw", "heavy_atoms", "rot_bonds", "hbd", "hba", "tpsa", "logp", "n_rings",
        "n_arom_rings", "fsp3", "qed", "n_stereo", "n_heteroatoms", "n_halogens"]
BOOL = ["has_halogen", "has_sulfur", "has_phosphorus", "has_metal"]
NICE = {"mw": "Mol. weight", "heavy_atoms": "Heavy atoms", "rot_bonds": "Rotatable bonds",
        "hbd": "H-bond donors", "hba": "H-bond acceptors", "tpsa": "TPSA", "logp": "logP",
        "n_rings": "Rings", "n_arom_rings": "Aromatic rings", "fsp3": "Fsp3",
        "qed": "QED (drug-likeness)", "n_stereo": "Stereocentres",
        "n_heteroatoms": "Heteroatoms", "n_halogens": "Halogen atoms"}


def cliffs_delta(worst, best) -> float:
    """Cliff's delta of *worst* vs *best* (fraction of worst>best minus worst<best)."""
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
                    default=Path("PoseBusters_Benchmark_Analysis/ligand_protein_features.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/ligand_difficulty"))
    ap.add_argument("--n", type=int, default=60,
                    help="Size of the best/worst docked groups compared (default 60).")
    ap.add_argument("--diffdock-variant", default="all",
                    choices=("all", "raw", "smina", "gnina"),
                    help="Restrict DiffDock to one optimizer variant before the per-tool "
                         "min-RMSD oracle (default 'all' pools raw+smina+gnina).")
    mf.add_method_filter_args(ap)
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
    # Before the fold below, which collapses every variant key onto its bare
    # engine name and makes the arms indistinguishable.
    met = mf.apply_method_filter(met, "method", args, label="ligand-difficulty",
                                 out_dir=args.out_dir)
    met["tool"] = met["method"].astype(str).str.replace(r"_.*$", "", regex=True)
    oracle = (met.groupby(["tool", "protein"])["rmsd"].min().reset_index()
              .rename(columns={"protein": "entry", "rmsd": "oracle_rmsd"})
              .merge(feat, on="entry", how="left"))

    delta_tbl, med_rows, succ = {}, [], {}
    for tool in TOOLS:
        sub = oracle[oracle["tool"] == tool].dropna(subset=["oracle_rmsd"]).sort_values("oracle_rmsd")
        if sub.empty:
            continue
        best, worst = sub.head(n), sub.tail(n)          # low RMSD = docks well
        succ[tool] = dict(n_complexes=len(sub),
                          best_rmsd_med=best["oracle_rmsd"].median(),
                          worst_rmsd_med=worst["oracle_rmsd"].median(),
                          pct_success_le2=100 * (sub["oracle_rmsd"] <= 2).mean())
        for p in CONT:
            d = cliffs_delta(worst[p], best[p])
            delta_tbl.setdefault(NICE.get(p, p), {})[tool] = d
            med_rows.append(dict(tool=tool, property=NICE.get(p, p),
                                 well_median=best[p].median(), badly_median=worst[p].median(),
                                 cliffs_delta=d))
        for p in BOOL:
            med_rows.append(dict(tool=tool, property=p,
                                 well_median=100 * best[p].mean(),
                                 badly_median=100 * worst[p].mean(),
                                 cliffs_delta=(worst[p].mean() - best[p].mean())))

    dtbl = pd.DataFrame(delta_tbl).T.reindex(columns=[t for t in TOOLS if t in succ])
    dtbl["mean"] = dtbl.mean(axis=1)
    dtbl = dtbl.sort_values("mean", key=lambda s: s.abs(), ascending=False)
    med = pd.DataFrame(med_rows)
    med.to_csv(args.out_dir / "ligand_property_best_vs_worst.csv", index=False)

    print(f"=== per-tool docking success (oracle best RMSD; success = ≤ 2 Å) — best/worst {n} ===")
    print(pd.DataFrame(succ).T.round(2).to_string())
    print(f"\n=== Cliff's delta (worst-{n} vs best-{n}); + = larger in HARD-to-dock ligands ===")
    print(dtbl.round(2).to_string())
    print("\n=== element flags: prevalence % (well vs badly docked) ===")
    b = med[med.property.isin(BOOL)].pivot_table(index="property",
            columns="tool", values=["well_median", "badly_median"])
    print(b.round(0).to_string())

    cols = [t for t in TOOLS if t in succ]

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
            for p in CONT:
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
        with open(args.out_dir / "ligand_difficulty_delta_heatmap_stats.json", "w") as fh:
            json.dump({"n_split": n, "unit": "one oracle-best-RMSD value per complex",
                       "family": "mannwhitney_cliffs, BH across all heatmap cells",
                       "cells": cell_records}, fh, indent=2)

        # Full-sample Spearman: descriptor vs continuous oracle RMSD over ALL
        # complexes (not just the N=60 tails), BH across the 14 descriptors/tool.
        spear_rows = []
        for tool in cols:
            sub = oracle[oracle["tool"] == tool].dropna(subset=["oracle_rmsd"])
            recs = []
            for p in CONT:
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

    fig, ax = plt.subplots(figsize=(1.6 * len(cols) + 3.5, 0.42 * len(dtbl) + 1.6))
    sns.heatmap(dtbl[cols], annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-0.7, vmax=0.7, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Cliff's δ  (+ = larger in hard-to-dock ligands)"}, ax=ax)
    for i, prop in enumerate(dtbl.index):          # overlay BH-significant stars
        for j, tool in enumerate(cols):
            s = star_map.get((prop, tool))
            if s:
                ax.text(j + 0.5, i + 0.16, s, ha="center", va="center",
                        fontsize=9, fontweight="bold", color="black")
    ax.set_xlabel("Docking tool"); ax.set_ylabel("Ligand property")
    ax.set_title(f"What ligand properties separate the best-{n} from worst-{n}\n"
                 "docked complexes, per tool (oracle RMSD-to-crystal)\n"
                 "* / ** / *** = Mann–Whitney worst-vs-best, BH-corrected q<.05/.01/.001",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(args.out_dir / "ligand_difficulty_delta_heatmap.png", dpi=160)
    plt.close(fig)
    print(f"\nWrote CSV + heatmap to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
