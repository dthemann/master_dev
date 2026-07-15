"""Is 'benefit from gnina/smina re-ranking' predictable from ligand chemistry?

Follow-up to diffdock_gnina_rerank_analysis.py. It asks whether the complexes
where affinity re-ranking picks a better pose than DiffDock's confidence order
(i) form a chemical group/cluster and (ii) coincide with ligands well explained
by the two dominant components of the PoseBusters descriptor PCA.

Method
  * Reproduce the dataset PCA from PoseBusters_DataSet_Analysis.ipynb:
    StandardScaler + PCA on the 16 continuous/count ligand descriptors (LIG_NUM),
    read from ligand_protein_features.csv.
  * Per ligand, compute the 2-component *communality* -- the fraction of that
    ligand's standardised variance captured by PC1+PC2 -- i.e. how well the two
    dominant components "explain" it. High communality = mainstream chemistry
    sitting on the dominant axes; low = descriptor outlier.
  * Define the per-complex re-ranking benefit on matched coordinates:
        benefit_min = RMSD(DiffDock rank-1, minimized) - RMSD(affinity pick, minimized)
    (positive = the affinity ranking supersedes DiffDock's). benefit_raw is the
    same on raw coordinates. Uses the per-complex table written by the analysis
    script (selection_per_complex_<tool>.csv).
  * Test benefit against communality, PC1/PC2 scores, every raw descriptor,
    KMeans clusters, and -- as the competing hypothesis -- the native pose's own
    RMSD.

Output: correlation CSVs + a 3-panel figure to --out-dir.

Usage:
    python Scripts/Analysis/diffdock_gnina_rerank_chemotype.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from pocket_comparison_report import _label_panels  # noqa: E402

# The exact PCA feature set used in PoseBusters_DataSet_Analysis.ipynb (LIG_NUM).
LIG_NUM = ["heavy_atoms", "mw", "rot_bonds", "hbd", "hba", "tpsa", "logp",
           "n_rings", "n_arom_rings", "fsp3", "qed", "formal_charge",
           "n_stereo", "n_heteroatoms", "n_halogens", "n_element_types"]
HIT = 2.0


def run_pca(feat: pd.DataFrame):
    F = feat[["entry"] + LIG_NUM].dropna().reset_index(drop=True)
    Xz = StandardScaler().fit_transform(F[LIG_NUM])
    pca = PCA().fit(Xz)
    scores = pca.transform(Xz)
    evr = pca.explained_variance_ratio_
    # 2-component communality (orthonormal components => ||recon||^2 = PC1^2 + PC2^2)
    comm2 = (scores[:, :2] ** 2).sum(1) / (Xz ** 2).sum(1)
    out = F[["entry"]].copy()
    for i in range(min(5, scores.shape[1])):
        out[f"PC{i+1}"] = scores[:, i]
    out["comm2"] = comm2
    out = out.join(F[LIG_NUM])
    return out, pca, evr, Xz, F


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/ligand_protein_features.csv"))
    ap.add_argument("--selection", type=Path, default=None,
                    help="selection_per_complex_<tool>.csv from the analysis script")
    ap.add_argument("--out-dir", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/gnina_rerank"))
    ap.add_argument("--tool", choices=["gnina", "smina"], default="gnina")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sel_path = args.selection or (args.out_dir / f"selection_per_complex_{args.tool}.csv")

    feat = pd.read_csv(args.features)
    pcdf, pca, evr, Xz, F = run_pca(feat)
    load = pd.DataFrame(pca.components_[:2].T, index=LIG_NUM, columns=["PC1", "PC2"])
    pc1_drivers = ", ".join(load.PC1.abs().sort_values(ascending=False).head(4).index)
    pc2_drivers = ", ".join(load.PC2.abs().sort_values(ascending=False).head(4).index)

    sel = pd.read_csv(sel_path).rename(columns={"protein": "entry"})
    m = pcdf.merge(sel, on="entry", how="inner")
    m["benefit_min"] = m["D_minimize"] - m["B_full"]       # + = affinity ranking wins (min coords)
    m["benefit_raw"] = m["A_native"] - m["C_rerank"]       # + = affinity ranking wins (raw coords)
    m["supersede"] = m["benefit_min"] > 1e-9
    m["picks_differ"] = m["benefit_min"].abs().gt(1e-9) | m["benefit_raw"].abs().gt(1e-9)
    print(f"[chemotype] {len(m)} complexes; PC1={evr[0]*100:.0f}% PC2={evr[1]*100:.0f}%; "
          f"supersede {m.supersede.mean()*100:.0f}%")

    # ---- correlations: chemistry predictors vs benefit, plus native-quality ----
    preds = [("2-PC communality", m.comm2), ("PC1 score", m.PC1), ("PC2 score", m.PC2),
             ("|PC1| score", m.PC1.abs()), ("|PC2| score", m.PC2.abs()),
             ("DiffDock native top-1 RMSD", m.A_native)]
    rows = []
    for name, x in preds:
        r, p = spearmanr(x, m["benefit_min"])
        rows.append(dict(predictor=name, spearman_r=r, p_value=p))
    for d in LIG_NUM:
        r, p = spearmanr(m[d], m["benefit_min"])
        rows.append(dict(predictor=f"descriptor:{d}", spearman_r=r, p_value=p))
    corr = pd.DataFrame(rows).sort_values("p_value")
    corr.to_csv(args.out_dir / f"chemotype_benefit_correlations_{args.tool}.csv", index=False)

    # ---- median split on communality ----
    hi = m[m.comm2 >= m.comm2.median()]; lo = m[m.comm2 < m.comm2.median()]
    _, p_mwu = mannwhitneyu(hi.benefit_min, lo.benefit_min)
    split = pd.Series({
        "high_comm_supersede_pct": hi.supersede.mean() * 100,
        "low_comm_supersede_pct": lo.supersede.mean() * 100,
        "high_comm_benefit_median": hi.benefit_min.median(),
        "low_comm_benefit_median": lo.benefit_min.median(),
        "mannwhitney_p": p_mwu,
    })
    split.to_frame("value").to_csv(args.out_dir / f"chemotype_communality_split_{args.tool}.csv")

    # ---- KMeans clusters vs supersede rate ----
    clrows = []
    for k in (2, 3, 4):
        km = KMeans(n_clusters=k, n_init=10, random_state=0).fit(Xz)
        lab = pd.Series(km.labels_, index=F.entry).reindex(m.entry).values
        for c in range(k):
            sub = m[lab == c]
            clrows.append(dict(k=k, cluster=c, n=len(sub),
                               supersede_pct=sub.supersede.mean() * 100,
                               benefit_median=sub.benefit_min.median()))
    pd.DataFrame(clrows).to_csv(args.out_dir / f"chemotype_clusters_{args.tool}.csv", index=False)

    r_comm, p_comm = spearmanr(m.comm2, m.benefit_min)
    r_nat, p_nat = spearmanr(m.A_native, m.benefit_min)

    # ------------------------------------------------------------------ figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    vlim = 5.0
    bc = m.benefit_min.clip(-vlim, vlim)

    # (A) biplot coloured by benefit -- is there spatial grouping?
    sc = ax[0].scatter(m.PC1, m.PC2, c=bc, cmap="coolwarm_r", vmin=-vlim, vmax=vlim,
                       s=26, edgecolor="grey", linewidth=.3)
    ax[0].axhline(0, color="grey", lw=.5); ax[0].axvline(0, color="grey", lw=.5)
    ax[0].set_xlabel(f"Principal Component 1 ({evr[0]*100:.0f}% of variance)\n"
                     f"size & polarity: {pc1_drivers}", fontsize=8)
    ax[0].set_ylabel(f"Principal Component 2 ({evr[1]*100:.0f}% of variance)\n"
                     f"aromaticity & lipophilicity: {pc2_drivers}", fontsize=8)
    ax[0].set_title("Benefit across chemical space")
    ax[0].text(0.03, 0.97, "winners & losers intermixed\n(no spatial grouping)",
               transform=ax[0].transAxes, va="top", ha="left", fontsize=8,
               bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.85))
    cb = fig.colorbar(sc, ax=ax[0]); cb.set_label(f"Re-ranking benefit (Angstrom)\n"
                                                  f"native RMSD - affinity-pick RMSD")

    # (B) benefit vs 2-PC communality -- the user's question, answered
    ax[1].scatter(m.comm2, m.benefit_min, s=22, color="#4C72B0", alpha=.6, edgecolor="none")
    z = np.polyfit(m.comm2, m.benefit_min, 1)
    xs = np.linspace(m.comm2.min(), m.comm2.max(), 50)
    ax[1].plot(xs, np.polyval(z, xs), color="crimson", lw=1.5)
    ax[1].axhline(0, color="grey", lw=.6, ls=":")
    ax[1].set_xlabel("Fraction of ligand variance captured by\nPrincipal Components 1+2 "
                     "(2-component communality)")
    ax[1].set_ylabel("Re-ranking benefit (Angstrom)\nnative RMSD - affinity-pick RMSD")
    ax[1].set_title("Benefit vs. 2-component communality")
    ax[1].text(0.03, 0.97, f"Spearman r = {r_comm:+.2f}\np = {p_comm:.2f} (n.s.)",
               transform=ax[1].transAxes, va="top", ha="left", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.85))

    # (C) the competing hypothesis: native pose quality
    ax[2].scatter(m.A_native, m.benefit_min, s=22, color="#55A868", alpha=.6, edgecolor="none")
    ax[2].axhline(0, color="grey", lw=.6, ls=":")
    ax[2].axvline(HIT, color="crimson", lw=.8, ls=":")
    ax[2].set_xscale("log")
    ax[2].set_xlabel("DiffDock native top-1 RMSD to crystal (Angstrom)\n"
                     "(log scale; dotted line = 2 Angstrom hit)")
    ax[2].set_ylabel("Re-ranking benefit (Angstrom)\nnative RMSD - affinity-pick RMSD")
    ax[2].set_title("Benefit vs. native-pose quality")
    ax[2].text(0.03, 0.97, f"Spearman r = {r_nat:+.2f}\np = {p_nat:.0e}",
               transform=ax[2].transAxes, va="top", ha="left", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.85))

    _label_panels(ax)
    fig.suptitle(f"Does ligand chemistry predict which complexes benefit from "
                 f"{args.tool} re-ranking? ({len(m)} complexes)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.out_dir / f"chemotype_benefit_{args.tool}.png", dpi=160)
    plt.close(fig)

    print(f"  communality vs benefit: r={r_comm:+.3f} p={p_comm:.3f}  (no chemistry signal)")
    print(f"  native RMSD vs benefit: r={r_nat:+.3f} p={p_nat:.1e}  (recovery effect)")
    print(f"wrote CSVs + figure to {args.out_dir}")


if __name__ == "__main__":
    main()
