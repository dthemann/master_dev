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
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import stats_utils as su

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
        top = ("diffdock_gnina" if (dd["method"] == "diffdock_gnina").any()
               else dd.groupby("method")["pb_valid"].mean().idxmax())
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


# ── statistics helpers (unit of analysis = one value per complex) ──
def first_event_rank(sub: pd.DataFrame, mask: pd.Series) -> pd.Series:
    """Per-complex first rank at which *mask* is True (inf if never), indexed by
    protein. This collapses the many correlated poses of a complex to one number
    so the paired McNemar tests below see complexes, not poses."""
    complexes = sub["protein"].unique()
    return (sub[mask].groupby("protein")["rank"].min()
            .reindex(complexes).fillna(np.inf))


def cum_counts(sub: pd.DataFrame, mask: pd.Series, K: int):
    """(successes[k], n_complexes) for the cumulative 'event in ranks 1..k'."""
    first = first_event_rank(sub, mask).to_numpy(float)
    n = len(first)
    ks = np.arange(1, K + 1)
    succ = (first[:, None] <= ks[None, :]).sum(axis=0).astype(int)
    return succ, n


def per_rank_counts(sub: pd.DataFrame, col: str, K: int):
    """(successes[r], n_complexes_with_a_rank_r_pose) for the per-rank event."""
    succ, tot = [], []
    for r in range(1, K + 1):
        rr = sub[sub["rank"] == r].groupby("protein")[col].max()
        tot.append(int(len(rr)))
        succ.append(int(rr.sum()) if len(rr) else 0)
    return np.array(succ), np.array(tot)


def wilson_bands(succ, tot):
    """Per-point Wilson 95% CI (percent) for rate = succ/tot. Returns
    (rate%, lo%, hi%) arrays, NaN where tot==0."""
    rate, lo, hi = [], [], []
    for k, n in zip(np.asarray(succ), np.asarray(tot)):
        if n:
            l, h = su.wilson_ci(int(k), int(n))
            rate.append(100.0 * k / n); lo.append(100.0 * l); hi.append(100.0 * h)
        else:
            rate.append(np.nan); lo.append(np.nan); hi.append(np.nan)
    return np.array(rate), np.array(lo), np.array(hi)


def mcnemar_across_tools(data: dict, depth: int) -> list[dict]:
    """Pairwise exact McNemar of the per-complex boolean 'a hit within top-`depth`'
    across every pair of ranked entities (paired on the complexes both cover).
    Holm-corrects the pairwise family at this depth."""
    labels = list(data)
    firsts = {lbl: first_event_rank(sub, sub["_hit"]) for lbl, sub in data.items()}
    out = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            la, lb = labels[i], labels[j]
            fa, fb = firsts[la], firsts[lb]
            common = fa.index.intersection(fb.index)
            if len(common) < 2:
                continue
            av = (fa.loc[common].to_numpy() <= depth).astype(int)
            bv = (fb.loc[common].to_numpy() <= depth).astype(int)
            n10, n01, p = su.mcnemar_exact(av, bv)
            out.append({"depth": depth, "a": la, "b": lb, "n": int(len(common)),
                        "a_hits": int(av.sum()), "b_hits": int(bv.sum()),
                        "a_only": n10, "b_only": n01, "p_raw": float(p)})
    if out:
        for d, pa in zip(out, su.holm([o["p_raw"] for o in out])):
            d["p_holm"] = float(pa); d["star"] = su.p_stars(pa)
    return out


def abbrev(lbl: str) -> str:
    """Compact tag for on-figure annotation of pairwise tests."""
    if lbl.startswith("AutoDock"):
        return "AD"
    if lbl.startswith("DiffDock"):
        tag = "DD"
        if "gnina" in lbl:
            return tag + "-g"
        if "smina" in lbl:
            return tag + "-s"
        return tag
    if lbl.startswith("EquiBind"):
        return "EB"
    return lbl[:4]


def _annotate_mcnemar(ax, mcn: dict, depth_x: dict) -> None:
    """Mark the compared depths with dashed verticals and a compact corner box of
    pairwise McNemar stars (Holm-adjusted). Full numbers live in the sidecar."""
    if not mcn:
        return
    for x in depth_x.values():
        ax.axvline(x, color="0.5", ls="--", lw=0.8, alpha=0.5, zorder=0)
    lines = []
    for tag, x in depth_x.items():
        rows = mcn.get(tag, [])
        if not rows:
            continue
        lines.append(f"top-{x} McNemar:")
        for r in rows:
            lines.append(f"  {abbrev(r['a'])}–{abbrev(r['b'])} "
                         f"{r.get('star', 'ns')} (p={su.fmt_p(r.get('p_holm'))})")
    if not lines:
        return
    ax.text(0.985, 0.02, "\n".join(lines), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.2, family="monospace",
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))


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
            sub = sub.merge(feat, left_on="protein", right_on="entry", how="left")
            sub["_hit"] = sub["rmsd"] <= HIT
            data[pretty(e["method"])] = sub
    print("Ranked entities:", ", ".join(data))

    stats: dict = {"hit_threshold_A": HIT, "max_rank": K, "entities": list(data),
                   "unit_of_analysis": "one value per complex (paired across tools)"}

    # ── RMSD success + per-rank hit (all ranked entities) ──
    DEPTHS = {"top1": 0, "top5": min(4, K - 1), "top10": min(9, K - 1), "oracle": K - 1}
    rows = []
    stats["success_curves"] = {"success_at_k": {}, "per_rank_hit": {}, "mcnemar_top_k": {}}
    for lbl, sub in data.items():
        succ_c, n_c = cum_counts(sub, sub["_hit"], K)
        rate_c, lo_c, hi_c = wilson_bands(succ_c, np.full(K, n_c))
        row = dict(entity=lbl, n=n_c)
        for name, idx in DEPTHS.items():
            row[name] = round(float(rate_c[idx]), 1)
            row[f"{name}_lo"] = round(float(lo_c[idx]), 1)
            row[f"{name}_hi"] = round(float(hi_c[idx]), 1)
        rows.append(row)
        stats["success_curves"]["success_at_k"][lbl] = {
            "k": ks.tolist(), "n_complexes": int(n_c),
            "rate": rate_c.round(3).tolist(), "lo": lo_c.round(3).tolist(),
            "hi": hi_c.round(3).tolist()}
    summary = pd.DataFrame(rows).set_index("entity")
    summary.to_csv(args.out_dir / "topk_success_per_tool.csv")
    print("\n=== oracle vs top-n RMSD success (% of complexes with a pose ≤ 2 Å; Wilson 95% CI) ===")
    print(summary.round(1).to_string())

    # paired McNemar comparing tools at top-1 and top-5
    mcn = {}
    try:
        for tag, depth in (("top1", 1), ("top5", 5)):
            mcn[tag] = mcnemar_across_tools(data, depth)
            stats["success_curves"]["mcnemar_top_k"][tag] = mcn[tag]
    except Exception as ex:                                # noqa: BLE001
        print(f"[warn] McNemar across tools failed, figure drawn without stars: {ex}")
        mcn = {}

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5))
    for lbl, sub in data.items():
        succ_c, n_c = cum_counts(sub, sub["_hit"], K)
        rate_c, lo_c, hi_c = wilson_bands(succ_c, np.full(K, n_c))
        line, = a1.plot(ks, rate_c, marker="o", ms=3, label=lbl)
        succ_r, tot_r = per_rank_counts(sub, "_hit", K)
        rate_r, lo_r, hi_r = wilson_bands(succ_r, tot_r)
        line2, = a2.plot(ks, rate_r, marker="o", ms=3, label=lbl)
        try:
            a1.fill_between(ks, lo_c, hi_c, color=line.get_color(),
                            alpha=0.15, linewidth=0)
            a2.fill_between(ks, lo_r, hi_r, color=line2.get_color(),
                            alpha=0.15, linewidth=0)
        except Exception as ex:                            # noqa: BLE001
            print(f"[warn] Wilson band draw failed for {lbl}: {ex}")
        stats["success_curves"]["per_rank_hit"][lbl] = {
            "rate": rate_r.round(3).tolist(), "lo": lo_r.round(3).tolist(),
            "hi": hi_r.round(3).tolist(), "n_per_rank": tot_r.tolist()}
    a1.set_title("Top-k success — a hit among ranks 1..k")
    a1.set_ylabel("Complexes with a pose ≤ 2 Å (%)")
    a2.set_title("Per-rank hit probability — is the rank-r pose correct?")
    a2.set_ylabel("P(rank-r pose ≤ 2 Å) (%)")
    for ax in (a1, a2):
        ax.set_xlabel("Pose rank"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    _annotate_mcnemar(a1, mcn, {"top1": 1, "top5": 5})
    fig.suptitle("Oracle vs top-n: success rises with k; per-pose reliability falls with rank"
                 "\nbands = Wilson 95% CI · brackets at k=1,5 = paired McNemar (Holm)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.out_dir / "rank_success_curves.png", dpi=160)
    plt.close(fig)

    # ── PB-validity vs rank (NEW) ──
    vrows = []
    stats["validity_curves"] = {"per_rank_valid": {}, "cum_valid_at_k": {},
                                "mcnemar_valid_top_k": {}}
    fig, (b1, b2) = plt.subplots(1, 2, figsize=(13, 5))
    for lbl, sub in data.items():
        succ_r, tot_r = per_rank_counts(sub, "pb_valid", K)
        pr, plo, phi = wilson_bands(succ_r, tot_r)
        succ_c, n_c = cum_counts(sub, sub["pb_valid"], K)
        cum, clo, chi = wilson_bands(succ_c, np.full(K, n_c))
        lb, = b1.plot(ks, pr, marker="o", ms=3, label=lbl)
        lc, = b2.plot(ks, cum, marker="o", ms=3, label=lbl)
        try:
            b1.fill_between(ks, plo, phi, color=lb.get_color(), alpha=0.15, linewidth=0)
            b2.fill_between(ks, clo, chi, color=lc.get_color(), alpha=0.15, linewidth=0)
        except Exception as ex:                            # noqa: BLE001
            print(f"[warn] validity Wilson band draw failed for {lbl}: {ex}")
        vrows.append(dict(entity=lbl, valid_rank1=round(float(pr[0]), 1),
                          valid_rank1_lo=round(float(plo[0]), 1),
                          valid_rank1_hi=round(float(phi[0]), 1),
                          valid_rank_last=round(float(pr[-1]), 1),
                          cum_valid_top1=round(float(cum[0]), 1),
                          cum_valid_top1_lo=round(float(clo[0]), 1),
                          cum_valid_top1_hi=round(float(chi[0]), 1),
                          cum_valid_oracle=round(float(cum[-1]), 1),
                          cum_valid_oracle_lo=round(float(clo[-1]), 1),
                          cum_valid_oracle_hi=round(float(chi[-1]), 1)))
        stats["validity_curves"]["per_rank_valid"][lbl] = {
            "rate": pr.round(3).tolist(), "lo": plo.round(3).tolist(),
            "hi": phi.round(3).tolist(), "n_per_rank": tot_r.tolist()}
        stats["validity_curves"]["cum_valid_at_k"][lbl] = {
            "rate": cum.round(3).tolist(), "lo": clo.round(3).tolist(),
            "hi": chi.round(3).tolist(), "n_complexes": int(n_c)}
    # paired McNemar of 'a valid pose within top-k' across tools
    vmcn = {}
    try:
        data_valid = {lbl: sub.assign(_hit=sub["pb_valid"].astype(bool))
                      for lbl, sub in data.items()}
        for tag, depth in (("top1", 1), ("top5", 5)):
            vmcn[tag] = mcnemar_across_tools(data_valid, depth)
            stats["validity_curves"]["mcnemar_valid_top_k"][tag] = vmcn[tag]
    except Exception as ex:                                # noqa: BLE001
        print(f"[warn] validity McNemar failed, figure drawn without stars: {ex}")
        vmcn = {}
    b1.set_title("Per-rank PB-validity — is the rank-r pose valid?")
    b1.set_ylabel("P(rank-r pose is PB-valid) (%)")
    b2.set_title("Cumulative — a PB-valid pose among ranks 1..k")
    b2.set_ylabel("Complexes with a valid pose in top-k (%)")
    for ax in (b1, b2):
        ax.set_xlabel("Pose rank"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    _annotate_mcnemar(b2, vmcn, {"top1": 1, "top5": 5})
    fig.suptitle("PoseBusters validity vs pose rank"
                 "\nbands = Wilson 95% CI · brackets at k=1,5 = paired McNemar (Holm)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.out_dir / "rank_validity_curves.png", dpi=160)
    plt.close(fig)
    vtbl = pd.DataFrame(vrows).set_index("entity")
    vtbl.to_csv(args.out_dir / "validity_by_rank.csv")
    print("\n=== PB-validity by rank (rank-1 vs last, and cumulative top-1 vs oracle; Wilson 95% CI) ===")
    print(vtbl.round(1).to_string())

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

    # ordered flexibility tertiles are DISJOINT complex subsets -> Cochran-Armitage
    # trend (unpaired) at fixed depths; rigid<mid<flexible scored 0,1,2.
    stats["flexibility_trend"] = {}
    CA_DEPTHS = {"top1": 1, "top5": 5}
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
        blabels = list(bins)
        tert_n = [int(sel.sum()) for sel in bins.values()]
        for blbl, sel in bins.items():
            s = success_at_k(sub[sub["protein"].isin(sel[sel].index)], K)
            ax.plot(ks, s, marker="o", ms=3, label=f"{blbl}  (n={int(sel.sum())})")
        ax.set_title(f"{lbl} — top-k success by ligand flexibility")
        ax.set_xlabel("Pose rank"); ax.set_ylabel("Complexes with a pose ≤ 2 Å (%)")
        ax.grid(alpha=0.3); ax.legend(fontsize=8)

        # Cochran-Armitage trend at each depth (ordered rigid->flexible)
        rec = {"tertiles": blabels, "tertile_n": tert_n, "trend": {}}
        lines = ["Trend vs flexibility (Cochran-Armitage):"]
        small = sum(tert_n) < 6 or min(tert_n) < 1
        for tag, depth in CA_DEPTHS.items():
            succ = []
            for sel in bins.values():
                sc, _ = cum_counts(sub[sub["protein"].isin(sel[sel].index)],
                                   sub["_hit"], depth)
                succ.append(int(sc[depth - 1]) if len(sc) >= depth else 0)
            rec["trend"][tag] = {"successes": succ, "totals": tert_n}
            if small:
                lines.append(f"  {tag}: n too small — exploratory")
                continue
            try:
                z, p, sign = su.cochran_armitage(succ, tert_n)
                arrow = "↑" if sign > 0 else "↓" if sign < 0 else "~"
                rec["trend"][tag].update(z=z, p=p, sign=int(sign),
                                         star=su.p_stars(p))
                lines.append(f"  {tag}: z={z:+.2f} {arrow}flex "
                             f"{su.p_stars(p)} (p={su.fmt_p(p)})")
            except Exception as ex:                        # noqa: BLE001
                print(f"[warn] Cochran-Armitage failed for {lbl}/{tag}: {ex}")
                lines.append(f"  {tag}: n/a")
        stats["flexibility_trend"][lbl] = rec
        try:
            ax.text(0.985, 0.02, "\n".join(lines), transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=6.6, family="monospace",
                    bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
        except Exception as ex:                            # noqa: BLE001
            print(f"[warn] flexibility annotation failed for {lbl}: {ex}")
    fig.suptitle("Top-k success by ligand flexibility; corner box tests a monotone "
                 "success-vs-flexibility trend across ordered tertiles", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(args.out_dir / "rank_success_by_flexibility.png", dpi=160)
    plt.close(fig)

    # ── numeric sidecar (all test results, recoverable) ──
    try:
        with open(args.out_dir / "rank_success_stats.json", "w") as fh:
            json.dump(stats, fh, indent=2)
    except Exception as ex:                                # noqa: BLE001
        print(f"[warn] could not write rank_success_stats.json: {ex}")
    print(f"\nWrote CSVs + figures + rank_success_stats.json to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
