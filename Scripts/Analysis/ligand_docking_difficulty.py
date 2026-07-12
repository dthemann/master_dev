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
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    n = args.n

    met = pd.read_csv(args.per_pose_metrics, low_memory=False)
    feat = pd.read_csv(args.features)
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
    fig, ax = plt.subplots(figsize=(1.6 * len(cols) + 3.5, 0.42 * len(dtbl) + 1.6))
    sns.heatmap(dtbl[cols], annot=True, fmt=".2f", cmap="RdBu_r", center=0,
                vmin=-0.7, vmax=0.7, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "Cliff's δ  (+ = larger in hard-to-dock ligands)"}, ax=ax)
    ax.set_xlabel("Docking tool"); ax.set_ylabel("Ligand property")
    ax.set_title(f"What ligand properties separate the best-{n} from worst-{n}\n"
                 "docked complexes, per tool (oracle RMSD-to-crystal)")
    fig.tight_layout()
    fig.savefig(args.out_dir / "ligand_difficulty_delta_heatmap.png", dpi=160)
    plt.close(fig)
    print(f"\nWrote CSV + heatmap to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
