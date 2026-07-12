"""Oracle vs top-n success, how success/validity decline with pose rank, and how they
depend on ligand attributes — for each tool AND its top variant.

A pose is a docking "hit" when RMSD-to-crystal ≤ 2 Å, and "valid" when it passes every
PoseBusters check. Ranked entities analysed (each needs a real pose ordering):

  * AutoDock Vina                — score rank 1..30.
  * DiffDock (raw)               — confidence rank 1..30.
  * DiffDock top variant         — the diffdock optimizer variant with the highest
                                   PB-valid fraction (usually diffdock_gnina); the
                                   smina/gnina poses inherit the parent confidence rank
                                   (needs posebusters_pose_comparison.py --force
                                   --diffdock-variant all so the ranks are present).
  * EquiBind top variant         — EquiBind has no native ranking, so its smina-refined
                                   variants are ranked by smina docking score
                                   (smina_affinity, best = rank 1), the same idea as
                                   AutoDock's score rank. Raw/gnina EquiBind poses are
                                   unranked and cannot appear here.

Reports (to --out-dir):
  topk_success_per_tool.csv / rank_success_curves.png
        top-1/5/10/oracle RMSD-hit success; success@k curve + per-rank hit probability.
  validity_by_rank.csv / rank_validity_curves.png                              (NEW)
        per-rank PB-validity P(rank-r pose is valid) — the "validity declines with rank"
        view — plus cumulative top-k "≥1 valid pose" per complex.
  success_by_attribute_split.csv / rank_success_by_flexibility.png
        success vs ligand attributes (median split; rotatable-bond tertiles).

Usage:
    python Scripts/Analysis/rank_success_analysis.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HIT = 2.0
K_ATTR_TOOLS = ("autodock", "diffdock")     # the two with a clean confidence/score rank
ATTRS = {"mw": "MW", "rot_bonds": "Rotatable bonds", "tpsa": "TPSA", "qed": "QED"}


def pretty(m: str) -> str:
    if m == "autodock":
        return "AutoDock Vina"
    if m == "diffdock":
        return "DiffDock (raw)"
    if m == "diffdock_smina":
        return "DiffDock (smina-opt)"
    if m == "diffdock_gnina":
        return "DiffDock (gnina-opt)"
    if m.startswith("equibind"):
        toks = [t for t in m.split("_")[1:]]
        return "EquiBind (" + ", ".join(toks) + ")"
    return m


def resolve_rank(sub: pd.DataFrame, source: str, K: int) -> pd.DataFrame:
    """Return *sub* restricted to ranks 1..K. source='rank' uses the metrics rank
    column; source='smina_affinity' derives a score rank per complex (best = 1)."""
    sub = sub.copy()
    if source == "smina_affinity":
        sub = sub.dropna(subset=["smina_affinity"])
        sub["rank"] = (sub.groupby("protein")["smina_affinity"]
                       .rank(method="first", ascending=True).astype(int))
    sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
    return sub[sub["rank"].between(1, K)]


def pick_entities(met: pd.DataFrame) -> list[dict]:
    """AutoDock, raw DiffDock, the best-validity DiffDock variant, and the best-validity
    *rankable* (smina-refined) EquiBind variant — those actually present with a rank."""
    ents = []
    if (met["method"] == "autodock").any():
        ents.append(dict(method="autodock", source="rank"))
    if (met["method"] == "diffdock").any():
        ents.append(dict(method="diffdock", source="rank"))
    dd = met[met["method"].astype(str).str.match(r"diffdock_(smina|gnina)$")]
    dd = dd[dd["rank"].between(1, 999) & (dd["rank"] != 999)]      # ranked variants only
    if len(dd):
        top = dd.groupby("method")["pb_valid"].mean().idxmax()
        ents.append(dict(method=top, source="rank"))
    eqs = met[met["method"].astype(str).str.startswith("equibind")
              & met["method"].astype(str).str.contains("smina")
              & met["smina_affinity"].notna()]
    if len(eqs):
        top = eqs.groupby("method")["pb_valid"].mean().idxmax()
        ents.append(dict(method=top, source="smina_affinity"))
    return ents


def cum_event_at_k(sub: pd.DataFrame, mask: pd.Series, K: int) -> np.ndarray:
    """% of complexes with the event (mask True) somewhere in ranks 1..k, for k=1..K."""
    complexes = sub["protein"].unique()
    first = (sub[mask].groupby("protein")["rank"].min()
             .reindex(complexes).fillna(np.inf).to_numpy(float))
    ks = np.arange(1, K + 1)
    return (first[:, None] <= ks[None, :]).mean(axis=0) * 100 if len(complexes) else np.zeros(K)


def per_rank_prob(sub: pd.DataFrame, col: str, K: int) -> np.ndarray:
    """P(the rank-r pose has *col* True), r=1..K (over complexes with a rank-r pose)."""
    out = []
    for r in range(1, K + 1):
        rr = sub[sub["rank"] == r].groupby("protein")[col].max()
        out.append(float(rr.mean()) * 100 if len(rr) else np.nan)
    return np.array(out)


def success_at_k(sub: pd.DataFrame, K: int) -> np.ndarray:
    return cum_event_at_k(sub, sub["rmsd"] <= HIT, K)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-metrics", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "pose_comparison_report/per_pose_metrics.csv"))
    ap.add_argument("--features", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/ligand_protein_features.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/rank_success"))
    ap.add_argument("--max-rank", type=int, default=30)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    K = args.max_rank
    ks = np.arange(1, K + 1)

    met = pd.read_csv(args.per_pose_metrics, low_memory=False)
    met["pb_valid"] = met["pb_valid"].astype(bool)
    feat = pd.read_csv(args.features)[["entry"] + list(ATTRS)]

    ents = pick_entities(met)
    data = {}
    for e in ents:
        sub = resolve_rank(met[met["method"] == e["method"]], e["source"], K)
        if len(sub):
            data[pretty(e["method"])] = sub.merge(feat, left_on="protein",
                                                  right_on="entry", how="left")
    print("Ranked entities:", ", ".join(data))

    # ── RMSD success + per-rank hit (all ranked entities) ──
    rows = []
    for lbl, sub in data.items():
        s = success_at_k(sub, K)
        rows.append(dict(entity=lbl, top1=s[0], top5=s[min(4, K - 1)],
                         top10=s[min(9, K - 1)], oracle=s[-1]))
    summary = pd.DataFrame(rows).set_index("entity").round(1)
    summary.to_csv(args.out_dir / "topk_success_per_tool.csv")
    print("\n=== oracle vs top-n RMSD success (% of complexes with a pose ≤ 2 Å) ===")
    print(summary.to_string())

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for lbl, sub in data.items():
        a1.plot(ks, success_at_k(sub, K), marker="o", ms=3, label=lbl)
        a2.plot(ks, per_rank_prob(sub, "_hit", K) if "_hit" in sub else
                per_rank_prob(sub.assign(_hit=sub["rmsd"] <= HIT), "_hit", K),
                marker="o", ms=3, label=lbl)
    a1.set_title("Top-k success — a hit among ranks 1..k")
    a1.set_ylabel("Complexes with a pose ≤ 2 Å (%)")
    a2.set_title("Per-rank hit probability — is the rank-r pose correct?")
    a2.set_ylabel("P(rank-r pose ≤ 2 Å) (%)")
    for ax in (a1, a2):
        ax.set_xlabel("Pose rank"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("Oracle vs top-n: success rises with k; per-pose reliability falls with rank",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out_dir / "rank_success_curves.png", dpi=160)
    plt.close(fig)

    # ── PB-validity vs rank (NEW) ──
    vrows = []
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(13, 5))
    for lbl, sub in data.items():
        pr = per_rank_prob(sub, "pb_valid", K)
        cum = cum_event_at_k(sub, sub["pb_valid"], K)
        b1.plot(ks, pr, marker="o", ms=3, label=lbl)
        b2.plot(ks, cum, marker="o", ms=3, label=lbl)
        vrows.append(dict(entity=lbl, valid_rank1=pr[0],
                          valid_rank_last=pr[-1],
                          cum_valid_top1=cum[0], cum_valid_oracle=cum[-1]))
    b1.set_title("Per-rank PB-validity — is the rank-r pose valid?")
    b1.set_ylabel("P(rank-r pose is PB-valid) (%)")
    b2.set_title("Cumulative — a PB-valid pose among ranks 1..k")
    b2.set_ylabel("Complexes with a valid pose in top-k (%)")
    for ax in (b1, b2):
        ax.set_xlabel("Pose rank"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.suptitle("PoseBusters validity vs pose rank", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out_dir / "rank_validity_curves.png", dpi=160)
    plt.close(fig)
    vtbl = pd.DataFrame(vrows).set_index("entity").round(1)
    vtbl.to_csv(args.out_dir / "validity_by_rank.csv")
    print("\n=== PB-validity by rank (rank-1 vs last, and cumulative top-1 vs oracle) ===")
    print(vtbl.to_string())

    # ── attribute dependence (clean-rank tools only) ──
    split_rows = []
    for tool in K_ATTR_TOOLS:
        lbl = pretty(tool)
        if lbl not in data:
            continue
        sub = data[lbl]
        for a in ATTRS:
            per_c = sub.groupby("protein")[a].first(); med = per_c.median()
            for half, sel in (("low", per_c <= med), ("high", per_c > med)):
                cx = sel[sel].index
                s = success_at_k(sub[sub["protein"].isin(cx)], K)
                split_rows.append(dict(tool=lbl, attribute=ATTRS[a], half=half,
                                       n=len(cx), top1=s[0], oracle=s[-1]))
    split = pd.DataFrame(split_rows).round(1)
    split.to_csv(args.out_dir / "success_by_attribute_split.csv", index=False)

    fig, axes = plt.subplots(1, len(K_ATTR_TOOLS), figsize=(6.4 * len(K_ATTR_TOOLS), 5),
                             squeeze=False)
    for ax, tool in zip(axes[0], K_ATTR_TOOLS):
        lbl = pretty(tool)
        if lbl not in data:
            ax.axis("off"); continue
        sub = data[lbl]
        rb = sub.groupby("protein")["rot_bonds"].first()
        q = rb.quantile([1 / 3, 2 / 3]).to_numpy()
        bins = {f"≤{q[0]:.0f} rot-bonds (rigid)": rb <= q[0],
                f"{q[0]:.0f}–{q[1]:.0f}": (rb > q[0]) & (rb <= q[1]),
                f">{q[1]:.0f} (flexible)": rb > q[1]}
        for blbl, sel in bins.items():
            s = success_at_k(sub[sub["protein"].isin(sel[sel].index)], K)
            ax.plot(ks, s, marker="o", ms=3, label=f"{blbl}  (n={int(sel.sum())})")
        ax.set_title(f"{lbl} — top-k success by ligand flexibility")
        ax.set_xlabel("Pose rank"); ax.set_ylabel("Complexes with a pose ≤ 2 Å (%)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "rank_success_by_flexibility.png", dpi=160)
    plt.close(fig)
    print(f"\nWrote CSVs + figures to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
