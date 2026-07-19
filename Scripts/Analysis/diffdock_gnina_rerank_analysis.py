"""How gnina (or smina) minimization re-ranks DiffDock poses, and whether that
re-ranking picks a better top-1 pose than DiffDock's own confidence order.

DiffDock emits up to 30 poses per complex ordered by a learned *confidence*
score (rank 1 = most confident). The optional refinement step minimizes every
pose with gnina/smina and records a physics/CNN ``minimized_affinity``
(more negative = stronger predicted binding). Sorting a complex's poses by that
affinity yields a *different* ordering. This script quantifies the difference and
tests whether it helps.

Two questions, two data sources:

  A. RANKING STABILITY  (needs only the optimizer logs; every complex with a
     successful minimization counts, no crystal required)
       For each complex we compare DiffDock's confidence rank against the rank
       obtained by sorting the same poses on ``minimized_affinity``:
         * per-complex Spearman / Kendall rank correlation,
         * the fate of DiffDock's rank-1 pose under the affinity order,
         * the DiffDock origin of the affinity order's new rank-1 pose,
         * per-pose rank displacement.

  B. SELECTION BENEFIT  (needs per-pose crystal RMSD from per_pose_metrics.csv;
     restricted to the benchmark complexes that have a crystal reference)
       A 2x2 decomposition of the top-1 pose the pipeline would hand a user,
       scored by crystal RMSD (hit = RMSD <= --hit):

                        select by DiffDock rank      select by affinity
         raw coords     A native                     C rerank-only
         min. coords    D minimize-only              B full pipeline

       A vs C isolates the *re-ranking* effect (same coordinates, different
       pick); A vs D isolates the *coordinate minimization* effect on the native
       pick; A vs B is the full refinement pipeline. Oracle (best-of-pool) RMSD
       for each coordinate source is reported as the achievable ceiling.

Outputs (CSVs + PNGs) go to --out-dir.

Usage:
    python Scripts/Analysis/diffdock_gnina_rerank_analysis.py            # gnina (default)
    python Scripts/Analysis/diffdock_gnina_rerank_analysis.py --tool smina
"""
from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, kendalltau

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import stats_utils as su  # noqa: E402

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
# Shared (A),(B),(C)... panel labeller for multi-panel figures.
from pocket_comparison_report import _label_panels  # noqa: E402

# Colour-blind-safe palette shared with the sibling reports.
C_DD      = "#4C72B0"   # DiffDock native / raw
C_RERANK  = "#DD8452"   # affinity re-ranking
C_MINIM   = "#8172B3"   # coordinate minimization
C_FULL    = "#55A868"   # full refinement pipeline
C_ORACLE  = "#937860"   # oracle ceiling


# --------------------------------------------------------------------------- #
# Part A -- ranking stability from the optimizer logs                         #
# --------------------------------------------------------------------------- #
def load_optimizer_logs(root: Path, tool: str) -> pd.DataFrame:
    """Concatenate every per-complex optimization_log.csv, keep the successful
    minimizations for ``tool``, and drop rows without an affinity."""
    logs = sorted(glob.glob(str(root / "*" / "optimization_log.csv")))
    if not logs:
        raise SystemExit(f"no optimization_log.csv under {root}")
    df = pd.concat((pd.read_csv(f) for f in logs), ignore_index=True)
    df = df[(df["tool"] == tool) & (df["status"] == "success")].copy()
    df["minimized_affinity"] = pd.to_numeric(df["minimized_affinity"], errors="coerce")
    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="coerce")
    df["diffdock_confidence"] = pd.to_numeric(df["diffdock_confidence"], errors="coerce")
    df = df.dropna(subset=["minimized_affinity", "pose_rank"])
    # A complex could in principle log a pose_rank twice; keep the first.
    df = df.drop_duplicates(["protein_name", "pose_rank"])
    return df


def add_ranks(df: pd.DataFrame) -> pd.DataFrame:
    """Add DiffDock rank (from pose_rank) and the affinity rank (1 = most
    negative minimized_affinity) within each complex."""
    out = []
    for _, sub in df.groupby("protein_name", sort=False):
        sub = sub.copy()
        sub["dd_rank"] = sub["pose_rank"].rank(method="first").astype(int)
        sub["aff_rank"] = sub["minimized_affinity"].rank(method="first").astype(int)
        sub["n_poses"] = len(sub)
        sub["rank_change"] = sub["aff_rank"] - sub["dd_rank"]   # + demoted, - promoted
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def ranking_stability(df: pd.DataFrame, tool: str, out_dir: Path) -> pd.DataFrame:
    rows = []
    for combo, sub in df.groupby("protein_name", sort=False):
        if len(sub) < 3:
            continue
        # Correlate the raw scores, not the densified integer ranks: tied
        # minimized_affinity values (e.g. degenerate 0.0 minimizations) must be
        # averaged by scipy rather than broken in DiffDock order, which would
        # fabricate rank agreement. Both scores increase = worse, so a positive
        # coefficient still means the two orderings agree.
        rho, _ = spearmanr(sub["pose_rank"], sub["minimized_affinity"])
        tau, _ = kendalltau(sub["pose_rank"], sub["minimized_affinity"])
        aff_top_from = int(sub.loc[sub["aff_rank"] == 1, "dd_rank"].iloc[0])
        dd_top_to = int(sub.loc[sub["dd_rank"] == 1, "aff_rank"].iloc[0])
        rows.append(dict(protein=combo, n_poses=len(sub), spearman_rho=rho,
                         kendall_tau=tau, aff_top1_from_dd_rank=aff_top_from,
                         dd_top1_to_aff_rank=dd_top_to))
    per = pd.DataFrame(rows)
    per.to_csv(out_dir / f"ranking_stability_per_complex_{tool}.csv", index=False)

    dd1 = df[df["dd_rank"] == 1]
    aff1 = df[df["aff_rank"] == 1]
    absmove = df["rank_change"].abs()
    summ = pd.Series({
        "tool": tool,
        "n_complexes": df["protein_name"].nunique(),
        "n_poses": len(df),
        "spearman_rho_mean": per["spearman_rho"].mean(),
        "spearman_rho_median": per["spearman_rho"].median(),
        "pct_complexes_rho_below_0.2": (per["spearman_rho"] < 0.2).mean() * 100,
        "pct_complexes_rho_at_or_below_0": (per["spearman_rho"] <= 0).mean() * 100,
        "pct_complexes_rho_at_or_above_0.7": (per["spearman_rho"] >= 0.7).mean() * 100,
        "dd_top1_stays_aff_top1_pct": (dd1["aff_rank"] == 1).mean() * 100,
        "dd_top1_in_aff_top3_pct": (dd1["aff_rank"] <= 3).mean() * 100,
        "dd_top1_median_aff_rank": dd1["aff_rank"].median(),
        "aff_top1_was_dd_top1_pct": (aff1["dd_rank"] == 1).mean() * 100,
        "aff_top1_was_dd_top5_pct": (aff1["dd_rank"] <= 5).mean() * 100,
        "aff_top1_from_dd_bottom_half_pct": (aff1["dd_rank"] > aff1["n_poses"] / 2).mean() * 100,
        "aff_top1_median_dd_rank": aff1["dd_rank"].median(),
        "pose_median_abs_rank_move": absmove.median(),
        "pose_pct_exact_same_rank": (df["rank_change"] == 0).mean() * 100,
        "pose_pct_move_ge_5": (absmove >= 5).mean() * 100,
    })
    summ.to_frame("value").to_csv(out_dir / f"ranking_stability_summary_{tool}.csv")
    return per


def fig_ranking(df: pd.DataFrame, per: pd.DataFrame, tool: str, out_dir: Path) -> None:
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))

    # (A) rank-migration heatmap
    sub = df[(df["dd_rank"] <= 30) & (df["aff_rank"] <= 30)]
    H, _, _ = np.histogram2d(sub["dd_rank"], sub["aff_rank"], bins=np.arange(0.5, 31.5))
    im = ax[0].imshow(H.T, origin="lower", extent=[0.5, 30.5, 0.5, 30.5], aspect="auto",
                      cmap="viridis", norm=LogNorm(vmin=1, vmax=H.max()))
    ax[0].plot([0.5, 30.5], [0.5, 30.5], color="white", lw=1, ls="--", alpha=.75,
               label="rank unchanged")
    ax[0].set_xlabel("DiffDock confidence rank")
    ax[0].set_ylabel(f"Rank by {tool} minimized affinity")
    ax[0].set_title("Pose rank migration")
    ax[0].legend(loc="lower right", fontsize=8, framealpha=.8)
    cb = fig.colorbar(im, ax=ax[0]); cb.set_label("Number of poses")

    # (B) per-complex Spearman distribution
    ax[1].hist(per["spearman_rho"], bins=np.arange(-1, 1.05, 0.1),
               color=C_DD, edgecolor="white")
    ax[1].axvline(per["spearman_rho"].median(), color="crimson", lw=2,
                  label=f"median = {per['spearman_rho'].median():+.2f}")
    ax[1].axvline(0, color="grey", lw=1, ls=":")
    ax[1].set_xlabel("Spearman rank correlation per complex\n"
                     f"(DiffDock rank vs {tool} affinity rank)")
    ax[1].set_ylabel("Number of complexes")
    ax[1].set_title("Rank-order agreement")
    ax[1].legend()

    # (C) origin / fate of the top pose
    dd1 = df[df["dd_rank"] == 1]
    aff1 = df[df["aff_rank"] == 1]
    bins = np.arange(0.5, 31.5)
    ax[2].hist(dd1["aff_rank"], bins=bins, alpha=0.6, color=C_MINIM,
               label="affinity rank of DiffDock's top pose")
    ax[2].hist(aff1["dd_rank"], bins=bins, alpha=0.6, color=C_FULL,
               label=f"DiffDock rank of {tool}'s top pose")
    ax[2].set_xlabel("Rank position")
    ax[2].set_ylabel("Number of complexes")
    ax[2].set_title("Where the top pose comes from / goes to")
    ax[2].legend(fontsize=8)

    _label_panels(ax)
    fig.suptitle(f"DiffDock pose ranking before vs. after {tool} minimization "
                 f"({df['protein_name'].nunique()} complexes, {len(df):,} poses)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_dir / f"ranking_stability_{tool}.png", dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Part B -- does re-ranking pick a better pose (crystal RMSD)                  #
# --------------------------------------------------------------------------- #
def load_crystal_rmsd(per_pose: Path, tool: str) -> pd.DataFrame:
    """Per-pose crystal RMSD for raw DiffDock and the tool-minimized poses,
    keyed on (protein, rank). Drops the duplicate bare ``rank1.sdf`` row.

    LEFT-joins the minimized RMSD onto the full raw-DiffDock pose set so the raw
    coordinates (the native pick and the raw oracle) never depend on whether the
    minimizer produced/scored that pose; rmsd_min/pb_valid are NaN where absent."""
    d = pd.read_csv(per_pose, low_memory=False)
    d = d[~d["pose_name"].str.endswith("/rank1.sdf")]        # confidence-less duplicate
    raw = d[d["method"] == "diffdock"][["protein", "rank", "rmsd"]].rename(
        columns={"rmsd": "rmsd_raw"})
    opt = d[d["method"] == f"diffdock_{tool}"][["protein", "rank", "rmsd", "pb_valid"]].rename(
        columns={"rmsd": "rmsd_min"})
    m = raw.merge(opt, on=["protein", "rank"], how="left")
    return m


def selection_benefit(ranked: pd.DataFrame, rmsd: pd.DataFrame, tool: str,
                      hit: float, out_dir: Path) -> pd.DataFrame:
    """Join the affinity ranking to crystal RMSD and evaluate the 2x2 of
    (selection rule) x (coordinate source) plus oracle ceilings, per complex."""
    aff = ranked[["protein_name", "pose_rank", "minimized_affinity"]].rename(
        columns={"protein_name": "protein", "pose_rank": "rank"})
    j = rmsd.merge(aff, on=["protein", "rank"], how="left")

    rows = []
    for prot, sub in j.groupby("protein", sort=False):
        native = sub[sub["rank"] == 1]                       # DiffDock's confidence pick
        # affinity pick: best minimized_affinity among poses with a scored
        # minimized structure (lowest rank breaks ties = higher DiffDock confidence)
        scor = sub[sub["minimized_affinity"].notna() & sub["rmsd_min"].notna()].sort_values("rank")
        if native.empty or scor.empty:
            continue
        dd1 = native.iloc[0]
        aff1 = scor.loc[scor["minimized_affinity"].idxmin()]
        rows.append(dict(
            protein=prot, n_poses=len(sub),
            A_native=dd1["rmsd_raw"],               # DiffDock rank1, raw coords
            C_rerank=aff1["rmsd_raw"],              # affinity pick, raw coords
            D_minimize=dd1["rmsd_min"],             # DiffDock rank1, minimized coords (may be NaN)
            B_full=aff1["rmsd_min"],                # affinity pick, minimized coords
            oracle_raw=sub["rmsd_raw"].min(),               # over the full raw pool
            oracle_min=sub["rmsd_min"].min(),               # over scored minimized poses
            B_full_pb_valid=bool(aff1["pb_valid"]),
        ))
    per = pd.DataFrame(rows)
    per.to_csv(out_dir / f"selection_per_complex_{tool}.csv", index=False)

    strat = {"A_native": "DiffDock top-1 (native)",
             "C_rerank": "Affinity top-1, raw coords (re-rank only)",
             "D_minimize": "DiffDock top-1, minimized coords (minimize only)",
             "B_full": f"Affinity top-1, minimized coords (full {tool} pipeline)",
             "oracle_raw": "Oracle, raw coords (ceiling)",
             "oracle_min": "Oracle, minimized coords (ceiling)"}
    summ = pd.DataFrame([
        dict(strategy=k, label=v,
             success_pct=(per[k] <= hit).mean() * 100,
             median_rmsd=per[k].median(),
             mean_rmsd=per[k].mean())
        for k, v in strat.items()
    ])
    # paired comparison of the full pipeline vs the native pick
    d = per["B_full"] - per["A_native"]
    paired = pd.Series({
        "tool": tool, "hit_threshold_A": hit, "n_complexes": len(per),
        "native_success_pct": (per["A_native"] <= hit).mean() * 100,
        "full_pipeline_success_pct": (per["B_full"] <= hit).mean() * 100,
        "rerank_only_success_pct": (per["C_rerank"] <= hit).mean() * 100,
        "minimize_only_success_pct": (per["D_minimize"] <= hit).mean() * 100,
        "full_better_than_native_pct": (d < -1e-6).mean() * 100,
        "full_worse_than_native_pct": (d > 1e-6).mean() * 100,
        "median_rmsd_delta_full_minus_native": d.median(),
        # complexes the affinity pick rescues / breaks at the hit threshold
        "rescued_by_full_pct": ((per["A_native"] > hit) & (per["B_full"] <= hit)).mean() * 100,
        "broken_by_full_pct": ((per["A_native"] <= hit) & (per["B_full"] > hit)).mean() * 100,
    })
    summ.to_csv(out_dir / f"selection_strategy_summary_{tool}.csv", index=False)
    paired.to_frame("value").to_csv(out_dir / f"selection_paired_summary_{tool}.csv")
    return per


# Strategies compared on panel A (paired by complex). The four selectable rules
# feed the omnibus + pairwise McNemar; the two oracle ceilings only get Wilson CIs.
_CORE_STRATS = ["A_native", "C_rerank", "D_minimize", "B_full"]
_ALL_STRATS = _CORE_STRATS + ["oracle_raw", "oracle_min"]
_MIN_UNITS = 5   # below this, inference is not reported (figure still drawn)


def _find_pair(pairwise, x, y):
    """Locate the (orientation-agnostic) McNemar record for the {x, y} pair."""
    for pr in pairwise:
        if {pr["a"], pr["b"]} == {x, y}:
            return pr
    return None


def selection_stats(per: pd.DataFrame, tool: str, hit: float):
    """Real paired tests behind the selection figure. Returns a JSON-able dict
    (or None if too few complexes / on failure). Never raises."""
    try:
        n = len(per)
        # per-strategy Wilson CIs on the success proportion (success = <= hit)
        wilson = {}
        for k in _ALL_STRATS:
            kk = int((per[k] <= hit).sum())
            lo, hi = su.wilson_ci(kk, n)
            wilson[k] = {"k": kk, "n": n, "rate": kk / n if n else float("nan"),
                         "lo": lo, "hi": hi}
        rep = {"tool": tool, "hit": hit, "n": n,
               "panelA_strategy_success": {"wilson": wilson}}
        if n < _MIN_UNITS:
            rep["note"] = f"n={n} too small -- exploratory, no omnibus/pairwise test"
            return rep

        # (A) paired proportions across the four selectable strategies.
        succ = {k: (per[k] <= hit).astype(int) for k in _CORE_STRATS}
        pp = su.paired_proportions(succ, labels=_CORE_STRATS)
        fvn = _find_pair(pp["pairwise"], "A_native", "B_full")
        rep["panelA_strategy_success"].update({
            "omnibus": pp["omnibus"],
            "rates": {k: {"rate": v[0], "lo": v[1], "hi": v[2]}
                      for k, v in pp["rates"].items()},
            "pairwise": pp["pairwise"],
            "full_vs_native": fvn,
        })

        # (B) native vs full-pipeline top-1 RMSD, paired (drop incomplete pairs).
        pair = per[["A_native", "B_full"]].dropna()
        nat, full = pair["A_native"].to_numpy(), pair["B_full"].to_numpy()
        rb, pw, npair = su.wilcoxon_rankbiserial(nat, full)
        est, mlo, mhi = su.median_diff_ci(nat, full, paired=True)
        rep["panelB_rmsd_native_vs_full"] = {
            "wilcoxon": {"rank_biserial": rb, "p": pw, "star": su.p_stars(pw),
                         "n": npair,
                         "interpretation": "rank_biserial>0 & sig => full pipeline "
                                           "lower RMSD (closer to crystal)"},
            "median_diff_native_minus_full": {"est": est, "lo": mlo, "hi": mhi},
        }

        # (C) effect decomposition: each refinement step vs the native pick.
        decomp = {}
        for k, tag in [("C_rerank", "rerank_only_vs_native"),
                       ("D_minimize", "minimize_only_vs_native"),
                       ("B_full", "full_vs_native")]:
            decomp[tag] = _find_pair(pp["pairwise"], "A_native", k)
        rep["panelC_decomposition"] = decomp
        return rep
    except Exception as exc:   # degrade to the test-free figure
        print(f"[selection][warn] stats failed ({exc}); drawing without tests")
        return None


def fig_selection(per: pd.DataFrame, tool: str, hit: float, out_dir: Path,
                  stats=None) -> None:
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))

    # (A) success-rate bars for the 2x2 + oracle ceilings
    order = ["A_native", "C_rerank", "D_minimize", "B_full", "oracle_raw", "oracle_min"]
    names = ["DiffDock\ntop-1\n(native)", "Affinity pick\nraw coords\n(re-rank)",
             "DiffDock top-1\nmin. coords\n(minimize)",
             f"Affinity pick\nmin. coords\n(full {tool})",
             "Oracle\nraw", "Oracle\nmin."]
    cols = [C_DD, C_RERANK, C_MINIM, C_FULL, C_ORACLE, C_ORACLE]
    succ = [(per[k] <= hit).mean() * 100 for k in order]
    bars = ax[0].bar(range(len(order)), succ, color=cols, edgecolor="white")
    for b in bars[-2:]:
        b.set_alpha(0.5); b.set_hatch("//")
    for i, s in enumerate(succ):
        ax[0].text(i, s + 1, f"{s:.0f}", ha="center", va="bottom", fontsize=8)
    ax[0].set_xticks(range(len(order)))
    ax[0].set_xticklabels(names, fontsize=7)
    ax[0].set_ylabel(f"Top-1 poses within {hit:g} Angstrom of crystal (% of complexes)")
    titleA = f"Selection success ({len(per)} complexes)"
    ax[0].set_ylim(0, 100)

    # --- Wilson 95% CIs on every bar + omnibus + full-vs-native McNemar star ---
    if stats:
        try:
            wil = stats["panelA_strategy_success"].get("wilson", {})
            yerr_lo, yerr_hi = [], []
            for k in order:
                w = wil.get(k)
                if w:
                    yerr_lo.append((w["rate"] - w["lo"]) * 100)
                    yerr_hi.append((w["hi"] - w["rate"]) * 100)
                else:
                    yerr_lo.append(0.0); yerr_hi.append(0.0)
            ax[0].errorbar(range(len(order)), succ, yerr=[yerr_lo, yerr_hi],
                           fmt="none", ecolor="0.3", elinewidth=1, capsize=3, zorder=5)
            omni = stats["panelA_strategy_success"].get("omnibus")
            if omni and omni.get("p") == omni.get("p"):   # not NaN
                titleA += (f"\nCochran Q p={su.fmt_p(omni['p'])} "
                           f"{su.p_stars(omni['p'])} (paired, {len(_CORE_STRATS)} rules)")
            fvn = stats["panelA_strategy_success"].get("full_vs_native")
            if fvn:
                ia, ib = order.index("A_native"), order.index("B_full")
                ytop = max(succ[ia], succ[ib]) + 7
                ytop = min(ytop, 95)
                ax[0].plot([ia, ia, ib, ib],
                           [ytop - 2, ytop, ytop, ytop - 2], color="0.25", lw=1)
                ax[0].text((ia + ib) / 2, ytop + 0.5,
                           f"full vs native  McNemar p={su.fmt_p(fvn['p_holm'])} "
                           f"{su.p_stars(fvn['p_holm'])}",
                           ha="center", va="bottom", fontsize=7.5, color="0.2")
        except Exception as exc:
            print(f"[selection][warn] panel A annotation failed ({exc})")
    ax[0].set_title(titleA)

    # (B) paired per-complex RMSD, native vs full pipeline
    lo, hi = 0.2, max(per["A_native"].max(), per["B_full"].max()) * 1.1
    better = per["B_full"] < per["A_native"]
    ax[1].scatter(per.loc[better, "A_native"], per.loc[better, "B_full"], s=14,
                  color=C_FULL, alpha=.6, label="full pipeline closer")
    ax[1].scatter(per.loc[~better, "A_native"], per.loc[~better, "B_full"], s=14,
                  color=C_DD, alpha=.6, label="native closer / tie")
    ax[1].plot([lo, hi], [lo, hi], color="grey", lw=1, ls="--")
    ax[1].axhline(hit, color="crimson", lw=.8, ls=":")
    ax[1].axvline(hit, color="crimson", lw=.8, ls=":")
    ax[1].set_xscale("log"); ax[1].set_yscale("log")
    ax[1].set_xlim(lo, hi); ax[1].set_ylim(lo, hi)
    ax[1].set_xlabel("DiffDock native top-1 RMSD to crystal (Angstrom)")
    ax[1].set_ylabel(f"Full {tool} pipeline top-1 RMSD to crystal (Angstrom)")
    ax[1].set_title("Per-complex top-1 accuracy\n(below diagonal = re-ranking helped)")
    ax[1].legend(fontsize=8, loc="upper left")
    if stats and stats.get("panelB_rmsd_native_vs_full"):
        try:
            b = stats["panelB_rmsd_native_vs_full"]
            w = b["wilcoxon"]; md = b["median_diff_native_minus_full"]
            ax[1].text(0.97, 0.03,
                       f"Wilcoxon p={su.fmt_p(w['p'])} {su.p_stars(w['p'])}\n"
                       f"rank-biserial={w['rank_biserial']:+.2f} (n={w['n']})\n"
                       f"HL median (native-full)={md['est']:+.2f}\n"
                       f"[{md['lo']:+.2f}, {md['hi']:+.2f}] Angstrom",
                       transform=ax[1].transAxes, va="bottom", ha="right", fontsize=7.5,
                       bbox=dict(boxstyle="round", fc="white", ec="grey", alpha=.85))
        except Exception as exc:
            print(f"[selection][warn] panel B annotation failed ({exc})")

    # (C) effect decomposition relative to the native pick
    base = (per["A_native"] <= hit).mean() * 100
    eff = {
        "Re-ranking\nonly": (per["C_rerank"] <= hit).mean() * 100 - base,
        "Minimization\nonly": (per["D_minimize"] <= hit).mean() * 100 - base,
        f"Full {tool}\npipeline": (per["B_full"] <= hit).mean() * 100 - base,
    }
    cols3 = [C_RERANK, C_MINIM, C_FULL]
    bars = ax[2].bar(range(len(eff)), list(eff.values()), color=cols3, edgecolor="white")
    ax[2].axhline(0, color="black", lw=.8)
    # McNemar star for each step vs the native pick (paired, Holm across the family)
    decomp_star = ["", "", ""]
    if stats and stats.get("panelC_decomposition"):
        try:
            dc = stats["panelC_decomposition"]
            for i, tag in enumerate(["rerank_only_vs_native", "minimize_only_vs_native",
                                     "full_vs_native"]):
                rec = dc.get(tag)
                if rec is not None:
                    decomp_star[i] = su.p_stars(rec["p_holm"])
        except Exception as exc:
            print(f"[selection][warn] panel C annotation failed ({exc})")
    for i, v in enumerate(eff.values()):
        lbl = f"{v:+.1f}" + (f" {decomp_star[i]}" if decomp_star[i] else "")
        ax[2].text(i, v + (0.3 if v >= 0 else -0.3), lbl, ha="center",
                   va="bottom" if v >= 0 else "top", fontsize=9)
    ax[2].set_xticks(range(len(eff)))
    ax[2].set_xticklabels(list(eff.keys()), fontsize=8)
    ax[2].set_ylabel(f"Change in success rate vs DiffDock native\n"
                     f"(percentage points, hit = {hit:g} Angstrom)")
    titleC = "Effect decomposition"
    if any(decomp_star):
        titleC += "\n(paired McNemar vs native, Holm)"
    ax[2].set_title(titleC)

    _label_panels(ax)
    fig.suptitle(f"Does {tool} re-ranking pick a better DiffDock pose? "
                 f"Top-1 accuracy vs crystal", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_dir / f"selection_benefit_{tool}.png", dpi=160)
    plt.close(fig)


# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dockings-root", type=Path,
                    default=Path("Dockings/Benchmark_DiffDock"))
    ap.add_argument("--per-pose-metrics", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "pose_comparison_report/per_pose_metrics.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("PoseBusters_Benchmark_Analysis/gnina_rerank"))
    ap.add_argument("--tool", choices=["gnina", "smina"], default="gnina")
    ap.add_argument("--hit", type=float, default=2.0,
                    help="RMSD (Angstrom) success threshold")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    logs = load_optimizer_logs(args.dockings_root, args.tool)
    ranked = add_ranks(logs)
    print(f"[ranking] {ranked['protein_name'].nunique()} complexes, {len(ranked):,} "
          f"{args.tool} poses")
    per = ranking_stability(ranked, args.tool, args.out_dir)
    fig_ranking(ranked, per, args.tool, args.out_dir)

    if args.per_pose_metrics.exists():
        rmsd = load_crystal_rmsd(args.per_pose_metrics, args.tool)
        sel = selection_benefit(ranked, rmsd, args.tool, args.hit, args.out_dir)
        print(f"[selection] {len(sel)} complexes with a crystal reference")
        stats_rep = selection_stats(sel, args.tool, args.hit)
        if stats_rep is not None:
            (args.out_dir / f"selection_benefit_{args.tool}_stats.json").write_text(
                json.dumps(stats_rep, indent=2, default=str))
            omni = stats_rep.get("panelA_strategy_success", {}).get("omnibus")
            if omni:
                print(f"[selection] Cochran Q p={su.fmt_p(omni['p'])} "
                      f"(paired across {len(_CORE_STRATS)} strategies)")
        fig_selection(sel, args.tool, args.hit, args.out_dir, stats=stats_rep)
    else:
        print(f"[selection] skipped -- {args.per_pose_metrics} not found")

    print(f"wrote CSVs + figures to {args.out_dir}")


if __name__ == "__main__":
    main()
