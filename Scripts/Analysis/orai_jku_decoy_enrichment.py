#!/usr/bin/env python
"""Property-matched decoy enrichment for Orai x JKU interaction proxies.

Orai1 has no crystal reference, so a docked JKU pose cannot be scored by RMSD.
This script asks, per proxy, whether the three experimental CRAC modulators
(Synta-66, GSK-7975A, 2-APB) engage the pharmacologically-defined Orai1 anchors
*more than random benchmark ligands docked onto the same receptor* -- i.e. it
uses Orai x Benchmark as a same-receptor negative control.

Because the actives are polarity/lipophilicity extremes within the benchmark
distribution, a raw actives-vs-all-decoys contrast is confounded with
physicochemistry. So for each active we draw a DUD-E-style property-matched
decoy pool (nearest neighbours in standardised descriptor space, topological
dissimilarity Tanimoto < 0.3 already guaranteed) and express the active as a
percentile / z / Cliff's delta against that matched null. Everything is
aggregated to ONE value per ligand before any comparison; with only three
actives we report per-active percentiles descriptively and attach no p-value to
"actives beat decoys".

Both pandamap pose-summaries are already pb_valid + no_tm filtered. We restrict
to the tool variants shared exactly by both datasets (autodock, diffdock_smina);
EquiBind is excluded because JKU uses the gnina refiner while the benchmark map
only carries the smina EquiBind variant, and JKU EquiBind is 2 poses (2-APB only).
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- config -----
ROOT = Path("/home/manndo/master_dev")
MATCHED_VARIANTS = ["autodock", "diffdock_smina"]

# descriptors used for the DUD-E match -- identical source (pandamap report
# ligand_descriptors) for BOTH datasets, so no cross-source computation bias.
MATCH_DESCS = ["mw", "logp", "hbd", "hba", "n_aromatic_rings", "n_halogen"]

ACTIVE_SDF = {  # docked protonation states, for secondary (rdkit) descriptors
    "2abp-nh2-OPT":          ROOT / "Drugs/JKU/2abp-nh2-OPT.sdf",
    "Synta-66-OPT-Singlet":  ROOT / "Drugs/JKU/Synta-66-OPT-Singlet.sdf",
    "gsk7975a-prot-OPT":     ROOT / "Drugs/JKU/gsk7975a-prot-OPT.sdf",
}

# directional / specific interaction types (quality, not bulk hydrophobic count)
DIRECTIONAL = {"hydrogen_bonds", "halogen_bonds", "pi_pi_stacking",
               "cation_pi", "pi_cation", "salt_bridge"}

# Orai1 anchor residue sets (human numbering; ring of six across chains A-F).
ANCHORS = {
    "filter_E106":        {106},                                   # selectivity filter glutamate ring
    "vestibule_polar":    {108, 113, 114},                         # Q108/H113/D114 mutagenesis-backed polar anchors
    "synta_vestibule":    {107, 108, 109, 110, 111, 112, 113, 115, 202},  # R2/loop1+loop3 Synta66 pocket
    "gate_V102_F99":      {99, 102},                               # hydrophobic gate
}
PHARM_ANCHOR = set().union(*ANCHORS.values()) | {199, 201}          # + loop3 F199/P201
FUNNEL_R1 = {76, 80, 83, 87, 91}                                    # basic cytosolic funnel = OFF-target (decoy-enriched)

# proxy definitions: (key, label, direction) direction=+1 on-target (high=good), -1 off-target
PROXIES = [
    ("frac_dir_anchor",   "directional contact to a pharmacological anchor (frac poses)", +1),
    ("frac_synta_vest",   "contacts Synta66 vestibule R2/loop3 (frac poses)",             +1),
    ("frac_vest_polar",   "contacts Q108/H113/D114 polar anchors (frac poses)",           +1),
    ("frac_gate",         "contacts V102/F99 gate (frac poses)",                          +1),
    ("frac_filter",       "contacts E106 filter ring (frac poses)",                       +1),
    ("mean_directional",  "directional interactions per pose (mean)",                     +1),
    ("mean_dir_fraction", "directional fraction of all contacts (mean)",                  +1),
    ("hbond_per_polar",   "H-bonds per ligand polar atom (mean H-bonds / (HBD+HBA))",     +1),
    ("frac_funnel",       "contacts R83/K87/R91 basic funnel [OFF-target] (frac poses)",  -1),
    ("mean_total",        "total interactions per pose (mean) [quantity baseline]",        0),
    ("mean_nres",         "distinct residues contacted per pose (mean) [quantity baseline]",0),
]


# ------------------------------------------------------------- descriptors ---
def load_match_descriptors(jku_report: Path, bench_report: Path):
    """Identical-source descriptors (pandamap report) for both datasets."""
    keep = ["ligand"] + MATCH_DESCS
    j = pd.read_csv(jku_report).drop_duplicates("ligand")[keep].reset_index(drop=True)
    b = pd.read_csv(bench_report).drop_duplicates("ligand")[keep].reset_index(drop=True)
    return j, b


def rdkit_secondary(sdf_paths: dict) -> pd.DataFrame:
    """tpsa / rot_bonds / qed for the actives (secondary balance check only)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, QED, rdMolDescriptors
    except Exception as e:  # pragma: no cover
        print(f"  [warn] rdkit unavailable ({e}); skipping secondary descriptors")
        return pd.DataFrame(columns=["ligand", "tpsa", "rot_bonds", "qed"])
    rows = []
    for lig, path in sdf_paths.items():
        mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)), None)
        if mol is None:
            print(f"  [warn] could not parse {path}")
            continue
        rows.append(dict(ligand=lig,
                         tpsa=round(Descriptors.TPSA(mol), 2),
                         rot_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
                         qed=round(QED.qed(mol), 3)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------- proxies ------
def pose_proxy_table(pose_summary: Path, interactions: Path) -> pd.DataFrame:
    """One row per pose with every proxy ingredient, matched-variant only."""
    ps = pd.read_csv(pose_summary)
    ps = ps[ps["method"].isin(MATCHED_VARIANTS)].copy()
    ps["directional"] = ps[list(DIRECTIONAL)].sum(axis=1)
    ps["dir_fraction"] = np.where(ps["total_interactions"] > 0,
                                  ps["directional"] / ps["total_interactions"], 0.0)
    key = ["method", "protein", "ligand", "pose_name"]
    pose = ps[key + ["ligand", "hydrogen_bonds", "directional", "dir_fraction",
                     "total_interactions", "n_residues"]].copy()
    pose = pose.loc[:, ~pose.columns.duplicated()]

    ic = pd.read_csv(interactions,
                     usecols=["interaction_type", "resnum", "method",
                              "protein", "ligand", "pose_name"])
    ic = ic[ic["method"].isin(MATCHED_VARIANTS)].copy()
    ic["pid"] = ic[key].astype(str).agg("|".join, axis=1)
    pose["pid"] = pose[key].astype(str).agg("|".join, axis=1)

    # per-pose contact flags
    grp = ic.groupby("pid")
    ressets = grp["resnum"].agg(lambda s: set(s.tolist()))
    dir_anchor = ic[ic["interaction_type"].isin(DIRECTIONAL) &
                    ic["resnum"].isin(PHARM_ANCHOR)].groupby("pid").size()

    flags = pd.DataFrame(index=ressets.index)
    flags["frac_filter"]      = ressets.apply(lambda s: len(s & ANCHORS["filter_E106"]) > 0)
    flags["frac_vest_polar"]  = ressets.apply(lambda s: len(s & ANCHORS["vestibule_polar"]) > 0)
    flags["frac_synta_vest"]  = ressets.apply(lambda s: len(s & ANCHORS["synta_vestibule"]) > 0)
    flags["frac_gate"]        = ressets.apply(lambda s: len(s & ANCHORS["gate_V102_F99"]) > 0)
    flags["frac_funnel"]      = ressets.apply(lambda s: len(s & FUNNEL_R1) > 0)
    flags["frac_dir_anchor"]  = flags.index.isin(dir_anchor.index)
    flags = flags.reset_index()

    pose = pose.merge(flags, on="pid", how="left")
    boolcols = ["frac_filter", "frac_vest_polar", "frac_synta_vest",
                "frac_gate", "frac_funnel", "frac_dir_anchor"]
    pose[boolcols] = pose[boolcols].where(pose[boolcols].notna(), False).astype(bool)
    return pose


def aggregate_per_ligand(pose: pd.DataFrame, descr: pd.DataFrame,
                         min_poses: int) -> pd.DataFrame:
    g = pose.groupby("ligand")
    out = pd.DataFrame({
        "n_poses":          g.size(),
        "frac_dir_anchor":  g["frac_dir_anchor"].mean(),
        "frac_synta_vest":  g["frac_synta_vest"].mean(),
        "frac_vest_polar":  g["frac_vest_polar"].mean(),
        "frac_gate":        g["frac_gate"].mean(),
        "frac_filter":      g["frac_filter"].mean(),
        "frac_funnel":      g["frac_funnel"].mean(),
        "mean_directional": g["directional"].mean(),
        "mean_dir_fraction":g["dir_fraction"].mean(),
        "mean_hbonds":      g["hydrogen_bonds"].mean(),
        "mean_total":       g["total_interactions"].mean(),
        "mean_nres":        g["n_residues"].mean(),
    }).reset_index()
    out = out.merge(descr, on="ligand", how="left")
    out["hbond_per_polar"] = out["mean_hbonds"] / (out["hbd"] + out["hba"]).replace(0, np.nan)
    return out[out["n_poses"] >= min_poses].reset_index(drop=True)


# ------------------------------------------------------------- matching -------
def match_decoys(active_row, bench, mu, sd, tanimoto, n_decoys, max_tan=0.3):
    az = ((active_row[MATCH_DESCS].astype(float) - mu) / sd).values
    cand = bench.copy()
    cand["max_tanimoto"] = cand["ligand"].map(tanimoto).fillna(0.0)
    cand = cand[cand["max_tanimoto"] < max_tan].copy()
    bz = (cand[MATCH_DESCS].astype(float) - mu) / sd
    cand["dist"] = np.sqrt(((bz.values - az) ** 2).sum(axis=1))
    cand = cand.sort_values("dist").head(n_decoys).reset_index(drop=True)
    return cand


def cliffs_delta(active_val, decoy_vals):
    gt = np.sum(active_val > decoy_vals)
    lt = np.sum(active_val < decoy_vals)
    return (gt - lt) / len(decoy_vals)


def enrich(active_val, decoy_vals):
    from scipy.stats import percentileofscore
    decoy_vals = np.asarray(decoy_vals, float)
    decoy_vals = decoy_vals[~np.isnan(decoy_vals)]
    if len(decoy_vals) == 0 or np.isnan(active_val):
        return dict(percentile=np.nan, z=np.nan, cliff=np.nan,
                    decoy_median=np.nan, n_decoys=0)
    sd = decoy_vals.std(ddof=1)
    return dict(
        percentile=round(percentileofscore(decoy_vals, active_val, kind="mean"), 1),
        z=round((active_val - decoy_vals.mean()) / sd, 2) if sd > 0 else np.nan,
        cliff=round(cliffs_delta(active_val, decoy_vals), 3),
        decoy_median=round(float(np.median(decoy_vals)), 4),
        n_decoys=int(len(decoy_vals)),
    )


# --------------------------------------------------------------- figures ------
def make_figures(enr_df, per_lig, actives, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [p[0] for p in PROXIES]
    labels = {p[0]: p[1] for p in PROXIES}
    piv = enr_df.pivot(index="proxy", columns="active", values="percentile").reindex(order)
    piv = piv[actives]

    fig, ax = plt.subplots(figsize=(8.2, 0.55 * len(order) + 1.8))
    im = ax.imshow(piv.values, cmap="RdBu_r", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(len(actives)))
    ax.set_xticklabels([a.split("-OPT")[0].replace("-Singlet", "") for a in actives], fontsize=9)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([labels[k] for k in order], fontsize=8)
    for i in range(len(order)):
        for j in range(len(actives)):
            v = piv.values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        color="white" if (v < 25 or v > 75) else "black", fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("percentile of active within property-matched decoy null", fontsize=8)
    ax.set_title("Orai x JKU interaction proxies vs same-receptor matched decoys\n"
                 "(50 = indistinguishable from decoys; high = on-target enrichment)",
                 fontsize=9.5)
    fig.tight_layout()
    fig.savefig(out / "enrichment_percentile_heatmap.png", dpi=150)
    plt.close(fig)

    # strip plots: decoy null distribution + active markers, key proxies
    key = ["frac_dir_anchor", "frac_synta_vest", "frac_gate", "frac_filter",
           "mean_dir_fraction", "frac_funnel"]
    fig, axes = plt.subplots(2, 3, figsize=(12, 6.5))
    rng = np.random.default_rng(0)
    for ax, k in zip(axes.ravel(), key):
        for ai, a in enumerate(actives):
            pool = per_lig[per_lig["active_pool"] == a]
            dv = pool[k].dropna().values
            x = ai + (rng.random(len(dv)) - 0.5) * 0.5
            ax.scatter(x, dv, s=8, alpha=0.35, color="0.6")
            av = per_lig[(per_lig["ligand"] == a) & (per_lig["is_active"])][k]
            if len(av):
                ax.scatter([ai], av.values[:1], s=90, color="crimson",
                           marker="D", zorder=5, edgecolor="k", linewidth=0.5)
        ax.set_xticks(range(len(actives)))
        ax.set_xticklabels([a.split("-OPT")[0].replace("-Singlet", "") for a in actives],
                           fontsize=8, rotation=15)
        ax.set_title(dict((p[0], p[1]) for p in PROXIES)[k], fontsize=8)
    fig.suptitle("Active (red) vs its property-matched decoy pool (grey), per proxy", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out / "enrichment_stripplots.png", dpi=150)
    plt.close(fig)


# ----------------------------------------------------------------- main -------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jku-pose-summary", default=ROOT / "pandamap_results/orai_jku/pandamap_pose_summary.csv")
    ap.add_argument("--jku-interactions", default=ROOT / "pandamap_results/orai_jku/pandamap_interactions.csv")
    ap.add_argument("--bench-pose-summary", default=ROOT / "pandamap_results/orai_benchmark/pandamap_pose_summary.csv")
    ap.add_argument("--bench-interactions", default=ROOT / "pandamap_results/orai_benchmark/pandamap_interactions.csv")
    ap.add_argument("--jku-descriptors", default=ROOT / "pandamap_results/orai_jku/report/ligand_descriptors.csv")
    ap.add_argument("--bench-descriptors", default=ROOT / "pandamap_results/orai_benchmark/report/ligand_descriptors.csv")
    ap.add_argument("--bench-features", default=ROOT / "PoseBusters_Benchmark_Analysis/ligand_protein_features.csv")
    ap.add_argument("--chem-similarity", default=ROOT / "posebusters_results/orai_benchmark/region_consensus/chem_similarity.csv")
    ap.add_argument("--out", default=ROOT / "pandamap_results/orai_jku/decoy_enrichment")
    ap.add_argument("--n-decoys", type=int, default=30)
    ap.add_argument("--min-poses", type=int, default=5)
    args = ap.parse_args()

    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    print(f"[1] descriptors  (match source = pandamap report, identical for both)")
    jku_d, bench_d = load_match_descriptors(Path(args.jku_descriptors), Path(args.bench_descriptors))
    sec = rdkit_secondary(ACTIVE_SDF)
    print("    active match descriptors:\n", jku_d.to_string(index=False))

    # descriptor percentiles of the actives within the benchmark pool (context)
    print("\n    active descriptor percentile within 308-benchmark pool:")
    for _, r in jku_d.iterrows():
        pcts = {d: round(float((bench_d[d] < r[d]).mean() * 100), 0) for d in MATCH_DESCS}
        print(f"      {r['ligand']:<22} " + "  ".join(f"{d}={pcts[d]:.0f}" for d in MATCH_DESCS))

    print(f"\n[2] pose proxy tables (variants={MATCHED_VARIANTS})")
    jku_pose = pose_proxy_table(Path(args.jku_pose_summary), Path(args.jku_interactions))
    bench_pose = pose_proxy_table(Path(args.bench_pose_summary), Path(args.bench_interactions))
    print(f"    JKU poses {len(jku_pose)} ({jku_pose['ligand'].nunique()} lig); "
          f"BENCH poses {len(bench_pose)} ({bench_pose['ligand'].nunique()} lig)")

    jku_lig = aggregate_per_ligand(jku_pose, jku_d, args.min_poses)
    bench_lig = aggregate_per_ligand(bench_pose, bench_d, args.min_poses)
    jku_lig["is_active"] = True; bench_lig["is_active"] = False
    actives = [a for a in ACTIVE_SDF if a in set(jku_lig["ligand"])]
    print(f"    per-ligand: {len(jku_lig)} actives, {len(bench_lig)} candidate decoys (>= {args.min_poses} poses)")

    # matching
    print(f"\n[3] property-matched decoy pools (n={args.n_decoys} each, Tanimoto<0.3)")
    tanimoto = pd.read_csv(args.chem_similarity).set_index("ligand")["max_tanimoto"].to_dict()
    mu = bench_d[MATCH_DESCS].astype(float).mean()
    sd = bench_d[MATCH_DESCS].astype(float).std(ddof=0).replace(0, 1.0)
    bench_feat = pd.read_csv(args.bench_features).rename(columns={"entry": "ligand"})

    pools, balance_rows, per_lig_rows, enr_rows = {}, [], [], []
    for a in actives:
        arow = jku_lig[jku_lig["ligand"] == a].iloc[0]
        pool_desc = match_decoys(arow, bench_d, mu, sd, tanimoto, args.n_decoys)
        pool_ligs = pool_desc["ligand"].tolist()
        pools[a] = pool_ligs
        pool_lig = bench_lig[bench_lig["ligand"].isin(pool_ligs)].copy()
        pool_lig["active_pool"] = a
        per_lig_rows.append(pool_lig)

        # balance table (match descriptors + secondary tpsa/rot_bonds/qed)
        a_sec = sec[sec["ligand"] == a]
        feat_pool = bench_feat[bench_feat["ligand"].isin(pool_ligs)]
        for d in MATCH_DESCS:
            dv = pool_desc[d].astype(float).values
            z = (float(arow[d]) - dv.mean()) / dv.std(ddof=1) if dv.std(ddof=1) > 0 else np.nan
            balance_rows.append(dict(active=a, descriptor=d, source="match",
                                     active_value=round(float(arow[d]), 3),
                                     decoy_median=round(float(np.median(dv)), 3),
                                     decoy_iqr=round(float(np.subtract(*np.percentile(dv, [75, 25]))), 3),
                                     std_diff=round(z, 2) if not np.isnan(z) else np.nan))
        for d, col in [("tpsa", "tpsa"), ("rot_bonds", "rot_bonds"), ("qed", "qed")]:
            if len(a_sec) and col in bench_feat.columns:
                dv = feat_pool[col].astype(float).values
                av = float(a_sec[d].values[0])
                z = (av - dv.mean()) / dv.std(ddof=1) if dv.std(ddof=1) > 0 else np.nan
                balance_rows.append(dict(active=a, descriptor=d, source="secondary(cross-src)",
                                         active_value=round(av, 3),
                                         decoy_median=round(float(np.median(dv)), 3),
                                         decoy_iqr=round(float(np.subtract(*np.percentile(dv, [75, 25]))), 3),
                                         std_diff=round(z, 2) if not np.isnan(z) else np.nan))

        # enrichment per proxy
        for key_, label, direction in PROXIES:
            e = enrich(float(arow[key_]), pool_lig[key_].values)
            enr_rows.append(dict(active=a, proxy=key_, label=label, direction=direction,
                                 active_value=round(float(arow[key_]), 4), **e))

    per_lig = pd.concat(per_lig_rows, ignore_index=True)
    # attach actives to their own pool rows for plotting
    act_plot = []
    for a in actives:
        r = jku_lig[jku_lig["ligand"] == a].copy(); r["active_pool"] = a
        act_plot.append(r)
    per_lig_plot = pd.concat([per_lig] + act_plot, ignore_index=True)
    enr_df = pd.DataFrame(enr_rows)
    bal_df = pd.DataFrame(balance_rows)

    # ---- write outputs
    print(f"\n[4] writing outputs -> {out}")
    pd.DataFrame([(a, d) for a, ds in pools.items() for d in ds],
                 columns=["active", "decoy_ligand"]).to_csv(out / "decoy_pool.csv", index=False)
    bal_df.to_csv(out / "property_balance.csv", index=False)
    pd.concat([jku_lig, bench_lig], ignore_index=True).to_csv(out / "per_ligand_proxies.csv", index=False)
    enr_df.to_csv(out / "per_active_enrichment.csv", index=False)
    make_figures(enr_df, per_lig_plot, actives, out)

    summary = dict(
        matched_variants=MATCHED_VARIANTS, n_decoys_per_active=args.n_decoys,
        min_poses=args.min_poses, actives=actives,
        n_candidate_decoys=int(len(bench_lig)),
        note=("Percentiles are per-active vs a property-matched decoy null; with n=3 actives "
              "these are descriptive, not a powered actives-vs-decoys test. EquiBind excluded "
              "(variant mismatch: JKU=gnina, benchmark map=smina). 2-APB is a boron compound "
              "and the polar/low-logP outlier -- interpret separately."),
        max_abs_std_diff_balance=round(float(bal_df["std_diff"].abs().max()), 2),
    )
    json.dump(summary, open(out / "summary.json", "w"), indent=2)

    # ---- console report
    print("\n" + "=" * 78)
    print("PROPERTY BALANCE  (matched descriptors; |std_diff| small => active inside pool)")
    for a in actives:
        m = bal_df[(bal_df["active"] == a) & (bal_df["source"] == "match")]
        s = bal_df[(bal_df["active"] == a) & (bal_df["source"] != "match")]
        mworst = m.reindex(m["std_diff"].abs().sort_values(ascending=False).index).iloc[0]
        print(f"  {a:<22} matched max|std_diff|={m['std_diff'].abs().max():.2f} "
              f"(worst {mworst['descriptor']}: act={mworst['active_value']} "
              f"decoy_med={mworst['decoy_median']})  |  unmatched(reported): "
              + ", ".join(f"{r['descriptor']} Δ={r['std_diff']:+.1f}σ" for _, r in s.iterrows()))
    print("\nPER-ACTIVE ENRICHMENT  (percentile within matched decoy null)")
    for key_, label, direction in PROXIES:
        tag = {1: "on ", -1: "OFF", 0: "qty"}[direction]
        cells = []
        for a in actives:
            r = enr_df[(enr_df["active"] == a) & (enr_df["proxy"] == key_)].iloc[0]
            cells.append(f"{a.split('-OPT')[0].replace('-Singlet',''):>10}:{r['percentile']:>5.0f}")
        print(f"  [{tag}] {label[:46]:<46} " + " ".join(cells))
    print("=" * 78)
    print(f"done. figures + CSVs in {out}")


if __name__ == "__main__":
    sys.exit(main())
