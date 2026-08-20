#!/usr/bin/env python
"""Compare cross-tool agreement: Orai × Experimental (JKU) vs Orai × Benchmark ligands.

Companion to ``orai_pose_cluster_report.py``. That script runs once per dataset and
writes a single-dataset ``panels/orai_cross_tool_agreement.png`` (kept untouched by
this script). Here we read the ``per_pair.csv`` each run already produced — one row
per (Orai MD frame × ligand) pair carrying the inter-tool consensus-site distances
``au_di_dist`` / ``au_eq_dist`` / ``di_eq_dist`` and the pooled
``mean_inter_tool_dist`` — and overlay the two ligand categories.

    Experimental ligands  →  "Exp. Ligands"          (Orai × JKU)
    Benchmark ligands     →  "Benchmark ligands"

Two figures are written (originals are never overwritten):

  1. orai_cross_tool_agreement_compare.png
       The original distance-distribution graph, faceted into the two ligand
       categories so each keeps its three tool-pair curves.

  2. orai_tool_agreement_compare.png
       (A) the requested whisker/box plot of inter-tool consensus-site distance per
           tool-pair, Exp vs Benchmark (lower = tools land closer = agree); plus
       (B) the recommended "how often" view — agreement RATE (percent of pairs with
           the two tools ≤ thr apart) with Wilson 95% CIs. A rate + CI is the honest
           frequency answer: the strict all-tools-agree flag is near zero in both sets
           (tools scatter along the elongated Orai pore), so a box plot of that
           binary would be a flat line — the distance distribution and the rate are
           what actually carry signal.

Inferential comparisons (Mann-Whitney Exp-vs-Bench, two-proportion tests) are written
to a companion ``orai_tool_agreement_compare_stats.txt`` rather than drawn on the
panels — the panels carry descriptive content (boxes, rates, Wilson CIs) only.

Run in the analysis env (conda env ``vina``: scipy + matplotlib):

    python Scripts/Analysis/orai_cross_tool_agreement_compare.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
try:                                                # shared, unit-tested stats helpers
    import stats_utils as su
except Exception:                                   # pragma: no cover
    su = None

# The three inter-tool consensus-site distance columns, with a display label and
# colour. Colours are drawn from the SHARED docking-tool palette used across the
# thesis (AutoDock #1f77b4 blue / DiffDock #ff7f0e orange / EquiBind #2ca02c green),
# matching 09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png. The pooled
# per-pair mean is appended as a fourth "All (mean/pair)" group (neutral grey) for
# the box plot / rate bars.
TOOL_PAIRS = [
    ("au_di_dist", "AutoDock Vina\n↔ DiffDock", "#1f77b4"),   # AutoDock blue
    ("au_eq_dist", "AutoDock Vina\n↔ EquiBind", "#ff7f0e"),   # DiffDock orange
    ("di_eq_dist", "DiffDock\n↔ EquiBind", "#2ca02c"),        # EquiBind green
]
_MEAN_COL = "mean_inter_tool_dist"
_AGG_GREY = "#7f7f7f"       # neutral colour for the pooled "All pairs" group

# Ligand-category encoding, matching 09f: the two datasets are told apart by SOLID
# (Benchmark) vs HATCHED + lighter (Experimental) fill of the SAME hue — the hue
# carries the tool-pair / tool, NOT the dataset — instead of a separate colour.
CAT_STYLE = {
    "bench": {"hatch": "",    "alpha": 0.80},   # Benchmark  → solid
    "exp":   {"hatch": "///", "alpha": 0.45},   # Experimental / JKU → hatched, lighter
}
_CAT_LEGEND_GREY = "#7f7f7f"     # neutral swatch for the Bench-vs-Exp legend (09f style)
# Shared per-tool palette (for the cluster-quality overlay, x = single tool).
TOOL_COLORS = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind": "#2ca02c"}
# raw-vs-filtered pose population: same house convention as CAT_STYLE — the hue
# carries the TOOL, the population is told apart by SOLID (raw = full cloud) vs
# HATCHED + lighter (filtered = outside-pore ∩ PB-valid subset).
POSE_SET_STYLE = {
    "raw":      {"hatch": "",    "alpha": 0.80, "label": "raw cloud (all generated poses)"},
    "filtered": {"hatch": "///", "alpha": 0.45, "label": "filtered (outside-pore ∩ PB-valid)"},
}

# Ligand names dropped from EVERY analysis by default (case-insensitive substring).
# gsk7975a-deprot-OPT is the DEPROTONATED protomer of GSK-7975A; the protonated form
# (gsk7975a-prot-OPT) is the one carried through the agreement analysis, so the deprot
# entry double-counts the same compound. Excluded to keep one entry per compound and
# to keep the quality overlay's ligand set in step with the agreement figures. Extend
# via --exclude-ligand.
DEFAULT_EXCLUDE_LIGANDS = ["gsk7975a-deprot"]


# ── IO ───────────────────────────────────────────────────────────────────

def _num(d: pd.DataFrame, c: str) -> np.ndarray:
    return pd.to_numeric(d[c], errors="coerce").to_numpy() if c in d.columns else np.array([])


def drop_excluded(d: Optional[pd.DataFrame], patterns: List[str]) -> Optional[pd.DataFrame]:
    """Drop rows whose 'ligand' matches any case-insensitive substring in `patterns`."""
    if d is None or "ligand" not in getattr(d, "columns", []):
        return d
    subs = [s.strip().lower() for s in (patterns or []) if s and s.strip()]
    if not subs:
        return d
    lig = d["ligand"].astype(str).str.lower()
    mask = lig.apply(lambda x: any(s in x for s in subs))
    return d[~mask].copy()


def load_dataset(label: str, key: str, csv: Path,
                 exclude: Optional[List[str]] = None) -> dict:
    """Load one per_pair.csv and read its match threshold from the sibling stats
    sidecar (falls back to the CLI --thr)."""
    d = drop_excluded(pd.read_csv(csv), exclude if exclude is not None else DEFAULT_EXCLUDE_LIGANDS)
    thr = None
    for name in ("orai_cluster_stats.json", "summary.json"):
        sj = csv.parent / name
        if sj.exists():
            try:
                thr = float(json.loads(sj.read_text()).get("match_thr_A"))
                break
            except Exception:
                pass
    return {"label": label, "key": key, "csv": csv, "df": d, "thr": thr,
            "n_pairs": int(len(d)),
            "n_ligands": int(d["ligand"].nunique()) if "ligand" in d else 0,
            "n_frames": int(d["frame"].nunique()) if "frame" in d else 0}


_PAIR_DIST_COLS = ("au_di_dist", "au_eq_dist", "di_eq_dist")


def strict_all3(ds: dict) -> tuple:
    """Measure the strict all-three-tools-agree rate on the units that actually have
    all three tools. Returns (k, n, pct). Units with a missing pairwise distance are
    two-tool units and are excluded from n rather than counted as disagreements."""
    d = ds["df"]
    if not all(c in d.columns for c in _PAIR_DIST_COLS):
        return (0, 0, float("nan"))
    thr = ds.get("thr") or 5.0
    dist = d[list(_PAIR_DIST_COLS)].apply(pd.to_numeric, errors="coerce")
    three = dist.notna().all(axis=1)
    n = int(three.sum())
    if not n:
        return (0, 0, float("nan"))
    k = int((three & (dist <= thr).all(axis=1)).sum())
    return (k, n, 100.0 * k / n)


def strict_all3_phrase(datasets: list) -> str:
    """One-line human summary of the strict flag across the loaded datasets."""
    bits = []
    for ds in datasets:
        k, n, pct = strict_all3(ds)
        bits.append(f"{ds['label']} {k}/{n}" + ("" if not n else f" = {pct:.1f}%"))
    return "; ".join(bits)


def keep_ligands(ds: dict, keep: str) -> dict:
    """Restrict a loaded dataset's per_pair rows to ligands whose name contains any
    of the comma-separated (case-insensitive) substrings in `keep`; no-op if empty.
    Used to carve out subsets such as the pore blockers (GSK-7975A, Synta66)."""
    if not keep.strip():
        return ds
    subs = [s.strip().lower() for s in keep.split(",") if s.strip()]
    d = ds["df"]
    if "ligand" not in d.columns:
        raise SystemExit(f"--keep given but {ds['csv']} has no 'ligand' column")
    mask = d["ligand"].astype(str).str.lower().apply(lambda x: any(s in x for s in subs))
    kept = d[mask].copy()
    if kept.empty:
        raise SystemExit(f"--keep '{keep}' matched no ligand in {ds['csv']}; "
                         f"available: {sorted(d['ligand'].astype(str).unique())}")
    ds = dict(ds)
    ds["df"] = kept
    ds["n_pairs"] = int(len(kept))
    ds["n_ligands"] = int(kept["ligand"].nunique())
    ds["n_frames"] = int(kept["frame"].nunique()) if "frame" in kept else 0
    ds["kept_ligands"] = sorted(kept["ligand"].astype(str).unique())
    return ds


# ── Figure 1: faceted distance distribution (the "updated" original) ───────

def fig_distribution(datasets: List[dict], thr: float, out_dir: Path, suffix: str = ""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Shared bins across BOTH datasets and ALL tool-pairs so no curve is clipped —
    # inter-tool sites on the elongated Orai pore sit up to ~90 Å apart.
    allv = np.concatenate([_num(ds["df"], c)[~np.isnan(_num(ds["df"], c))]
                           for ds in datasets for c, *_ in TOOL_PAIRS] + [np.array([0.0])])
    top = max(40.0, float(np.ceil(np.nanmax(allv) / 3.0) * 3.0)) if allv.size else 40.0
    bins = np.arange(0.0, top + 3.0, 3.0)

    fig, axes = plt.subplots(len(datasets), 1, figsize=(7.8, 5.2 * len(datasets)),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    # Panel letters FIRST: _label_panels goes through ax.set_title(loc="left"), and
    # set_title rebuilds the shared title-offset transform from rcParams whenever it
    # is called without an explicit pad. Labelling after the titles were set would
    # therefore silently reset their pad to the 6 pt default and drop each headline
    # onto the subordinate line below it. Setting the centre title last wins, and the
    # pad it carries lifts the letter clear of that subordinate line too.
    _label_panels(axes)
    for ax, ds in zip(axes, datasets):
        for col, lab, colr in TOOL_PAIRS:
            v = _num(ds["df"], col)
            v = v[~np.isnan(v)]
            if not v.size:
                continue
            w = np.ones(v.size) / v.size          # fraction of THIS tool-pair's pairs
            med = float(np.median(v))
            ax.hist(v, bins=bins, weights=w, histtype="step", lw=2.2, color=colr)
            ax.axvline(med, color=colr, ls=":", lw=1.3, alpha=0.75)
        ax.axvline(thr, color="k", ls="--", lw=1.0)
        # Each curve is normalised by its OWN denominator: a (MD frame × ligand) unit
        # feeds a curve only when BOTH tools of that pair produced a consensus site
        # there. So the per-tool-pair n IS the headline, carried on the title line;
        # the panel's unit accounting drops to a smaller, greyed subordinate line.
        # Earlier titles headlined the panel row count, which no curve is normalised
        # by, and quoted that count without reconciling it against ligands × frames.
        per_curve = " · ".join(
            f"{short} {int(np.count_nonzero(~np.isnan(_num(ds['df'], col)))):,}"
            for col, short in zip([c for c, *_ in TOOL_PAIRS],
                                  ("AD↔DD", "AD↔EB", "DD↔EB")))
        # Reconcile the row count against the ligand × frame grid. Quoting the row
        # count on its own left a silent arithmetic gap (308 × 4 = 1,232 possible
        # units against 1,215 rows): the missing units are the frame–ligand
        # combinations for which no tool produced a docked pose to cluster, verified
        # against per_pose.csv, which spans exactly the same 1,215 combinations.
        possible = ds["n_ligands"] * ds["n_frames"]
        empty = possible - ds["n_pairs"]
        acct = (f"{'' if empty else 'all '}{ds['n_pairs']:,} of the {ds['n_ligands']:,} "
                f"ligands × {ds['n_frames']} frames = {possible:,} possible units have "
                f"docked poses" + (f", {empty:,} have none" if empty else ""))
        # The subordinate line is offset in POINTS from the axes top (not in axes
        # fractions) so its clearance from both the frame and the headline is the
        # same whatever the panel height; the title pad is set to clear it.
        ax.set_title(f"{ds['label']} — n per curve  {per_curve}",
                     fontsize=10, pad=25)
        ax.annotate(acct, xy=(0.5, 1.0), xycoords="axes fraction",
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=8.5, color="0.35")
        # The y-axis label carries the "defined" qualifier, which keeps it off the
        # title block and leaves the subordinate line free for the unit accounting.
        ax.set_ylabel("Fraction of units where both tools have a consensus site")
        ax.grid(alpha=0.25); ax.set_axisbelow(True)
    axes[-1].set_xlabel("Inter-tool consensus-site distance (Å)")
    # One shared tool-pair legend centred in the band between the suptitle and the
    # panels — per-panel n/median differ, so the shared legend stays generic (each
    # tool-pair's median is still marked by its dotted vertical line in each panel).
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=colr, lw=2.2, label=lab.replace(chr(10), " "))
               for _col, lab, colr in TOOL_PAIRS]
    handles.append(Line2D([0], [0], color="k", ls="--", lw=1.0, label=f"agree ≤ {thr:g} Å"))
    fig.suptitle("Cross-tool agreement on Orai — distance between tools' consensus "
                 "sites, by ligand category", fontsize=13)
    fig.legend(handles=handles, ncol=2, fontsize=11, frameon=True,
               loc="upper center", bbox_to_anchor=(0.5, 0.945))
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = out_dir / f"orai_cross_tool_agreement_compare{suffix}.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ── Figure 2: whisker plot + agreement-rate bars ──────────────────────────

def _agree_rate(v: np.ndarray, thr: float) -> Tuple[int, int]:
    """(#pairs with the two tools within thr, #pairs where the pair is defined)."""
    v = v[~np.isnan(v)]
    return int(np.sum(v <= thr)), int(v.size)


def fig_agreement(datasets: List[dict], thr: float, out_dir: Path, suffix: str = "",
                  title: Optional[str] = None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    groups = TOOL_PAIRS + [(_MEAN_COL, "All pairs\n(mean/pair)", _AGG_GREY)]
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(10, 12.4))

    # ── (A) whisker/box plot of inter-tool distance, Exp vs Bench per tool-pair ──
    # Hue = tool-pair; the two datasets are solid (Bench) vs hatched/lighter (Exp).
    width = 0.34
    offs = np.linspace(-(len(datasets) - 1), (len(datasets) - 1), len(datasets)) * (width / 2 + 0.03)
    for di, ds in enumerate(datasets):
        st = CAT_STYLE[ds["key"]]
        data, pos, colrs = [], [], []
        for gi, (col, _lab, gcolr) in enumerate(groups):
            v = _num(ds["df"], col); v = v[~np.isnan(v)]
            data.append(v if v.size else np.array([np.nan]))
            pos.append(gi + offs[di]); colrs.append(gcolr)
        bp = axA.boxplot(data, positions=pos, widths=width, patch_artist=True,
                         showfliers=False, medianprops=dict(color="black", lw=1.6),
                         whiskerprops=dict(color="#555555"),
                         capprops=dict(color="#555555"))
        for box, gcolr in zip(bp["boxes"], colrs):
            box.set(facecolor=gcolr, alpha=st["alpha"], edgecolor=gcolr,
                    hatch=st["hatch"], lw=1.4)
        # Overlay the individual pairs (small n on the Exp side) as jittered points.
        for gi, v in enumerate(data):
            vv = v[~np.isnan(v)]
            if vv.size and vv.size <= 60:
                jitter = (np.linspace(-1, 1, vv.size) if vv.size > 1 else np.zeros(1)) * width * 0.28
                axA.scatter(np.full(vv.size, pos[gi]) + jitter, vv, s=14,
                            color=colrs[gi], edgecolor="white", lw=0.4, zorder=3)
    axA.axhline(thr, color="k", ls="--", lw=1.0)
    axA.text(-0.45, thr + 1.0, f"agree ≤ {thr:g} Å", va="bottom", ha="left",
             fontsize=8, color="0.3")
    axA.set_xticks(range(len(groups)))
    axA.set_xticklabels([g[1] for g in groups], fontsize=9)
    axA.set_ylabel("Inter-tool consensus-site distance (Å)")
    axA.set_title("Distance between tools' consensus sites\n(lower = tools land closer / agree)",
                  fontsize=11)
    axA.grid(alpha=0.25, axis="y"); axA.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_CAT_LEGEND_GREY,
                             alpha=CAT_STYLE[ds["key"]]["alpha"],
                             hatch=CAT_STYLE[ds["key"]]["hatch"],
                             edgecolor=_CAT_LEGEND_GREY, label=ds["label"])
               for ds in datasets]
    # Legend is drawn once at the figure level (between title and axes), below —
    # both panels share the same Exp-vs-Benchmark categories, so per-panel legends
    # would just duplicate it.

    # ── (B) agreement RATE (% of pairs ≤ thr) with Wilson 95% CIs ──
    # The fourth "All pairs (pooled)" bar keeps the whole panel on ONE denominator
    # unit — pair-comparisons — so every bar (and the caption's "percentage of pairs")
    # means the same thing; see the __all_pairs__ branch below.
    rate_groups = TOOL_PAIRS + [("__all_pairs__", "All pairs\n(pooled)", _AGG_GREY)]
    bar_colors = [g[2] for g in rate_groups]
    xb = np.arange(len(rate_groups))
    bw = 0.8 / len(datasets)
    ymax = 1.0
    for di, ds in enumerate(datasets):
        st = CAT_STYLE[ds["key"]]
        rates, los, his, labels = [], [], [], []
        for col, _lab, _c in rate_groups:
            if col == "__all_pairs__":
                # Pooled over all three tool-pairs: every DEFINED pair-comparison is
                # one observation, so k and n are the SUMS of the per-pair counts. The
                # denominator is therefore pair-comparisons, matching the three per-pair
                # bars and the "percentage of pairs" axis/caption — NOT a per-(frame×
                # ligand)-row "≥1 pair agrees" OR, whose denominator would be receptor-
                # ligand rows (a different unit) and so overstate the rate.
                k = n = 0
                for c, *_ in TOOL_PAIRS:
                    ki, ni = _agree_rate(_num(ds["df"], c), thr)
                    k += ki; n += ni
            else:
                k, n = _agree_rate(_num(ds["df"], col), thr)
            r = k / n if n else np.nan
            lo, hi = su.wilson_ci(k, n) if (su and n) else (np.nan, np.nan)
            rates.append(r * 100 if n else 0.0)
            los.append((r - lo) * 100 if n else 0.0)
            his.append((hi - r) * 100 if n else 0.0)
            labels.append(f"{k}/{n}" if n else "n/a")
            if n:
                ymax = max(ymax, hi * 100)
        xpos = xb + (di - (len(datasets) - 1) / 2) * bw
        axB.bar(xpos, rates, bw * 0.92, color=bar_colors, alpha=st["alpha"],
                hatch=st["hatch"], edgecolor=bar_colors, lw=1.2)
        yerr = np.clip([los, his], 0.0, None)   # guard tiny FP negatives from Wilson CI
        axB.errorbar(xpos, rates, yerr=yerr, fmt="none", ecolor="0.25",
                     elinewidth=1.1, capsize=3)
        # k/n above each Wilson upper cap so a 0-of-small-n bar is honestly labelled.
        for x, r, hi, lab, bc in zip(xpos, rates, his, labels, bar_colors):
            axB.text(x, r + hi + 0.6, lab, ha="center", va="bottom", fontsize=11,
                     color=bc)
    axB.set_xticks(xb)
    axB.set_xticklabels([g[1] for g in rate_groups], fontsize=9)
    axB.set_ylabel(f"Pairs with the two tools ≤ {thr:g} Å apart (percent)")
    _bits = []
    for _ds in datasets:
        _k, _n, _pct = strict_all3(_ds)
        _short = str(_ds["label"]).split()[0].rstrip(".")
        _bits.append(f"{_short} {_k}/{_n}" + ("" if not _n else f" ({_pct:.1f}%)"))
    _s3 = ", ".join(_bits)
    axB.set_title("How often tools agree — rate with Wilson 95% CI\n"
                  f"strict all-3-tools-agree on 3-tool units: {_s3}",
                  fontsize=10)
    axB.set_ylim(0, ymax * 1.18 + 2)
    axB.grid(alpha=0.25, axis="y"); axB.set_axisbelow(True)

    _label_panels(np.array([axA, axB]))
    fig.suptitle(title or "Cross-tool agreement on Orai — Experimental vs Benchmark ligands",
                 fontsize=13, y=0.985)
    # One shared legend, centred in the band between the suptitle and the two panels.
    fig.legend(handles=handles, ncol=len(datasets), fontsize=12, frameon=True,
               loc="upper center", bbox_to_anchor=(0.5, 0.93))
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = out_dir / f"orai_tool_agreement_compare{suffix}.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ── companion stats (inferential tests, kept OFF the panels) ───────────────

def _frame_robust_median(df: pd.DataFrame, col: str) -> Tuple[float, int]:
    """Aggregate each ligand's MD frames to one value (median over its frames), then
    return (median across ligands, n_ligands). Removes the frame pseudoreplication the
    pooled per-(frame×ligand)-pair median carries — MD frames of one ligand are snapshots
    of the same receptor-ligand system, not independent observations."""
    if "ligand" not in df.columns or df.empty:
        return float("nan"), 0
    per_lig = (df.assign(_v=pd.to_numeric(df[col], errors="coerce"))
                 .groupby("ligand")["_v"].median().dropna())
    return (float(per_lig.median()) if len(per_lig) else float("nan"), int(len(per_lig)))


def write_stats(datasets: List[dict], thr: float, out_dir: Path, suffix: str = ""):
    if len(datasets) != 2:
        return None
    exp, bench = datasets
    L: List[str] = []
    L.append("Cross-tool agreement on Orai — Experimental (JKU) vs Benchmark ligands")
    L.append("=" * 74)
    L.append(f"agree threshold = {thr:g} Å")
    for ds in datasets:
        L.append(f"  {ds['label']:<28} n={ds['n_pairs']:>5} pairs, "
                 f"{ds['n_ligands']} ligands, {ds['n_frames']} MD frames  ({ds['csv']})")
    L.append("")
    L.append("NOTE: the strict all-tools-agree flag (all pairwise distances ≤ thr), measured")
    L.append("over the units that carry all three tools, is " + strict_all3_phrase(datasets) + ".")
    L.append("It is near zero because tools scatter along the elongated Orai pore, so agreement")
    L.append("is read from the CONTINUOUS distance and the per-tool-pair rate, not the binary.")
    L.append("NOTE: both pipelines use the SMINA-refined DiffDock variant, so DiffDock's")
    L.append("positions are comparable across datasets. The Experimental n is tiny")
    L.append("(exploratory).")
    L.append("")

    # ── unit of analysis + frame-robust median (pseudoreplication caveat) ──
    L.append("-- Unit of analysis (what every box / histogram / median pools over) --")
    L.append("Each per_pair.csv row is ONE (MD frame × ligand) pair, and each Orai MD frame")
    L.append("contributes one point PER ligand. All medians and box plots here (and in both")
    L.append("figures) are computed over that POOLED set of pairs — they are NOT collapsed to")
    L.append("one value per frame. MD frames of the same ligand are snapshots of the same")
    L.append("receptor-ligand system, so the pooled pairs are NOT independent: treating them")
    L.append("as independent observations is mild pseudoreplication and overstates the")
    L.append("effective n. A frame-robust summary aggregates each ligand's frames to one")
    L.append("value first (median over its frames), then takes the median across ligands.")
    L.append(f"  {'dataset':<28} {'pooled median':>14}   {'frame-robust median':>20}")
    for ds in datasets:
        pooled = float(np.nanmedian(_num(ds["df"], _MEAN_COL))) if ds["n_pairs"] else float("nan")
        fr, nlig = _frame_robust_median(ds["df"], _MEAN_COL)
        L.append(f"  {ds['label']:<28} {pooled:9.1f} Å   {fr:12.1f} Å   "
                 f"(n={ds['n_pairs']} pairs = {nlig} ligands × ≤{ds['n_frames']} frames)")
    L.append("  (column = mean_inter_tool_dist, Å; the same pooling applies to every tool-pair.)")
    L.append("")

    # ── distance: Mann-Whitney Exp vs Bench per group (Holm across groups) ──
    L.append("-- Inter-tool distance, Exp vs Benchmark (Mann-Whitney U, Holm-adjusted) --")
    groups = TOOL_PAIRS + [(_MEAN_COL, "All pairs (mean/pair)", None)]
    mw_rows, pvals = [], []
    for col, lab, _c in groups:
        a = _num(exp["df"], col); a = a[~np.isnan(a)]
        b = _num(bench["df"], col); b = b[~np.isnan(b)]
        res = su.mannwhitney_cliffs(a, b) if su else None
        if res:
            mw_rows.append((lab.replace("\n", " "), res)); pvals.append(res["p"])
        else:
            L.append(f"  {lab.replace(chr(10), ' '):<26} n/a (empty group)")
    adj = su.holm(pvals) if (su and pvals) else pvals
    for (lab, res), pa in zip(mw_rows, adj):
        me, mb = res["medians"]
        L.append(f"  {lab:<26} med Exp={me:5.1f} Å  Bench={mb:5.1f} Å  "
                 f"Cliff δ={res['cliffs_delta']:+.2f}  "
                 f"p_raw={su.fmt_p(res['p'])}  p_holm={su.fmt_p(pa)} {su.p_stars(pa)}")
    L.append("")

    # ── rate: two-proportion (Fisher) per tool-pair ──
    L.append(f"-- Agreement rate (pairs ≤ {thr:g} Å), Exp vs Benchmark (Fisher exact) --")
    for col, lab, _c in TOOL_PAIRS:
        ke, ne = _agree_rate(_num(exp["df"], col), thr)
        kb, nb = _agree_rate(_num(bench["df"], col), thr)
        tp = su.two_proportion_test(ke, ne, kb, nb) if (su and ne and nb) else None
        re = f"{ke}/{ne}={100*ke/ne:.0f}%" if ne else "n/a"
        rb = f"{kb}/{nb}={100*kb/nb:.0f}%" if nb else "n/a"
        line = f"  {lab.replace(chr(10),' '):<26} Exp {re:<12} Bench {rb:<12}"
        if tp:
            line += f"  Fisher p={su.fmt_p(tp['p'])} {su.p_stars(tp['p'])}"
        L.append(line)
    # pooled over all three pairings (the figure's fourth "All pairs (pooled)" bar):
    # sum the per-pair (k, n) so the denominator stays pair-comparisons, not rows.
    kep = nep = kbp = nbp = 0
    for col, _lab, _c in TOOL_PAIRS:
        ke, ne = _agree_rate(_num(exp["df"], col), thr);   kep += ke; nep += ne
        kb, nb = _agree_rate(_num(bench["df"], col), thr); kbp += kb; nbp += nb
    tp = su.two_proportion_test(kep, nep, kbp, nbp) if (su and nep and nbp) else None
    re = f"{kep}/{nep}={100*kep/nep:.0f}%" if nep else "n/a"
    rb = f"{kbp}/{nbp}={100*kbp/nbp:.0f}%" if nbp else "n/a"
    line = f"  {'All pairs (pooled)':<26} Exp {re:<12} Bench {rb:<12}"
    if tp:
        line += f"  Fisher p={su.fmt_p(tp['p'])} {su.p_stars(tp['p'])}"
    L.append(line)
    L.append("")
    L.append("Panels carry descriptive content only (boxes, rates, Wilson CIs); the")
    L.append("inferential tests above are intentionally kept off the figure.")

    p = out_dir / f"orai_tool_agreement_compare{suffix}_stats.txt"
    p.write_text("\n".join(L) + "\n")
    return p


# ── per-tool cluster QUALITY overlay (Exp vs Benchmark) ────────────────────
# Reads the cluster_quality_per_tool.csv each dataset's --cluster-quality run
# produced (one row per frame×ligand×tool, per pose_set). This is the headline
# "cognate (Exp) vs decoy-context (Benchmark) cluster-quality" comparison, plus
# the raw-vs-filtered tightening seen side by side. Inferential tests go to a
# companion *_stats.txt; the panels stay descriptive.

try:                                                # shared metric spec + labels
    from orai_pose_cluster_report import QUALITY_METRICS
except Exception:                                   # pragma: no cover
    QUALITY_METRICS = [
        ("silhouette",  "Silhouette",          "Silhouette (−1…1)",                       True),
        ("ch_score",    "Calinski–Harabasz",   "Calinski–Harabasz index",                 True),
        ("db_score",    "Davies–Bouldin",      "Davies–Bouldin index",                    False),
        ("compactness", "Compactness",         "Within-cluster spread (Å)",               False),
        ("separation",  "Separation",          "Nearest-cluster distance (Å)",            True),
        ("stability",   "Bootstrap stability", "Dominant-cluster Jaccard stability (0…1)", True),
    ]

_TOOL_PRETTY = {"autodock": "AutoDock Vina", "diffdock": "DiffDock", "equibind": "EquiBind"}

# Quality metrics whose panel gets a robust y-cap (heavy-tailed with a few extreme
# fliers that would otherwise flatten every box). Calinski–Harabasz blows up to ~1e5
# on the occasional degenerate EquiBind single-cluster case.
_ROBUST_YCAP_METRICS = {"ch_score"}


def _tool_name(t: str) -> str:
    return _TOOL_PRETTY.get(str(t), str(t))


def _tool_order(vals) -> List[str]:
    seen = {str(v) for v in vals}
    order = [t for t in ("autodock", "diffdock", "equibind") if t in seen]
    return order + [t for t in sorted(seen) if t not in order]


def load_quality(csv: Path, exclude: Optional[List[str]] = None) -> Optional[pd.DataFrame]:
    if not csv.exists():
        return None
    d = drop_excluded(pd.read_csv(csv), exclude if exclude is not None else DEFAULT_EXCLUDE_LIGANDS)
    return d if {"tool", "pose_set"}.issubset(d.columns) and len(d) else None


def _whisker_cap(bps, pad: float = 1.08) -> float:
    """Highest whisker end across a list of boxplot dicts (× pad); 0.0 if none.
    Caps a heavy-tailed panel just above the boxes+whiskers so extreme fliers
    (drawn only as jittered points) don't flatten the axis."""
    tops = [float(np.max(c.get_ydata())) for bp in bps for c in bp.get("caps", [])
            if len(c.get_ydata())]
    return max(tops) * pad if tops else 0.0


def fig_cluster_quality_compare(q_exp: pd.DataFrame, q_bench: pd.DataFrame,
                                exp_label: str, bench_label: str, out_dir: Path,
                                pose_set: str = "raw"):
    """Headline: per-tool intrinsic cluster quality, Exp vs Benchmark, one panel
    per metric (x = tool, two boxes = the two ligand categories)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dsets = [("exp", exp_label, q_exp), ("bench", bench_label, q_bench)]
    tools = _tool_order(pd.concat([q_exp["tool"], q_bench["tool"]]))
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    width = 0.36
    for ax, (key, short, ylab, hb) in zip(axes.ravel(), QUALITY_METRICS):
        capped = key in _ROBUST_YCAP_METRICS
        bps, panel_vals = [], []
        for di, (dkey, _lab, q) in enumerate(dsets):
            st = CAT_STYLE[dkey]
            qs = q[q["pose_set"] == pose_set]
            data, pos, ns, colrs = [], [], [], []
            for ti, t in enumerate(tools):
                v = pd.to_numeric(qs[qs.tool == t][key], errors="coerce").dropna().to_numpy()
                data.append(v if v.size else np.array([np.nan]))
                pos.append(ti + (di - 0.5) * (width + 0.04)); ns.append(v.size)
                colrs.append(TOOL_COLORS.get(str(t), "#888888"))
                panel_vals.extend(v.tolist())
            bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True,
                            showfliers=False, medianprops=dict(color="black"))
            bps.append(bp)
            for box, gcolr in zip(bp["boxes"], colrs):
                box.set(facecolor=gcolr, alpha=st["alpha"], edgecolor=gcolr,
                        hatch=st["hatch"], lw=1.3)
            for p_, v, gcolr in zip(pos, data, colrs):  # jitter (small Exp n)
                vv = v[~np.isnan(v)]
                if 0 < vv.size <= 60:
                    j = (np.linspace(-1, 1, vv.size) if vv.size > 1 else np.zeros(1)) * width * 0.28
                    ax.scatter(np.full(vv.size, p_) + j, vv, s=10, color=gcolr,
                               edgecolor="white", lw=0.3, zorder=3)
        # Robust y-cap for heavy-tailed metrics (Calinski–Harabasz): the boxes already
        # hide fliers (showfliers=False) but the jittered points autoscale the axis, so a
        # few EquiBind blow-ups flatten every box. Cap just above the tallest whisker and
        # note how many points fall off-scale (kept honest, not silently truncated).
        if capped:
            cap = _whisker_cap(bps)
            if cap:
                pv = np.asarray(panel_vals, float); pv = pv[~np.isnan(pv)]
                n_clip = int(np.sum(pv > cap))
                ax.set_ylim(0, cap)
                if n_clip:
                    ax.text(0.975, 0.96, f"{n_clip} outlier(s) off-scale (max "
                            f"{pv.max():,.0f})", transform=ax.transAxes, ha="right",
                            va="top", fontsize=8.5, color="0.35",
                            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", lw=0.6))
        ax.set_xticks(np.arange(len(tools)))
        ax.set_xticklabels([_tool_name(t) for t in tools], fontsize=8.5)
        ax.set_ylabel(ylab, fontsize=9)
        ax.set_title(f"{short}  ({'↑ better' if hb else '↓ better'})", fontsize=10)
        ax.grid(alpha=0.25, axis="y"); ax.set_axisbelow(True)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_CAT_LEGEND_GREY,
                             alpha=CAT_STYLE[k]["alpha"], hatch=CAT_STYLE[k]["hatch"],
                             edgecolor=_CAT_LEGEND_GREY, label=lab) for k, lab, _ in dsets]
    _label_panels(axes)
    fig.suptitle(f"Per-tool cluster QUALITY on Orai — Experimental vs Benchmark ligands "
                 f"('{pose_set}' pose cloud; reference-free)", fontsize=13, y=0.995)
    # One shared legend, centred in the band between the suptitle and the panels
    # (matching orai_tool_agreement_compare.png) rather than inside panel (A).
    fig.legend(handles=handles, ncol=len(handles), fontsize=12, frameon=True,
               loc="upper center", bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    p = out_dir / "orai_cluster_quality_compare.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def fig_cluster_quality_filtering_compare(q_exp: pd.DataFrame, q_bench: pd.DataFrame,
                                          exp_label: str, bench_label: str, out_dir: Path):
    """Raw cloud vs filtered (outside-pore ∩ PB-valid) tightening, side by side for the two
    ligand categories (compactness — lower = tighter)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dsets = [(exp_label, q_exp), (bench_label, q_bench)]
    if not all({"raw", "filtered"}.issubset(set(q["pose_set"])) for _, q in dsets):
        return None
    tools = _tool_order(pd.concat([q_exp["tool"], q_bench["tool"]]))
    fig, axes = plt.subplots(2, 1, figsize=(8.5, 12), sharey=True)
    width = 0.36
    for ax, (lab, q) in zip(axes, dsets):
        for si, ps in enumerate(("raw", "filtered")):
            st = POSE_SET_STYLE[ps]
            data, pos, colrs = [], [], []
            for ti, t in enumerate(tools):
                v = pd.to_numeric(q[(q.tool == t) & (q.pose_set == ps)]["compactness"],
                                  errors="coerce").dropna().to_numpy()
                data.append(v if v.size else np.array([np.nan]))
                pos.append(ti + (si - 0.5) * (width + 0.04))
                colrs.append(TOOL_COLORS.get(str(t), "#888888"))
            bp = ax.boxplot(data, positions=pos, widths=width, patch_artist=True,
                            showfliers=False, medianprops=dict(color="black"))
            for box, gcolr in zip(bp["boxes"], colrs):
                box.set(facecolor=gcolr, alpha=st["alpha"], hatch=st["hatch"],
                        edgecolor=gcolr, lw=1.3)
        ax.set_xticks(np.arange(len(tools)))
        ax.set_xticklabels([_tool_name(t) for t in tools], fontsize=9)
        ax.set_title(lab, fontsize=11)
        ax.grid(alpha=0.25, axis="y"); ax.set_axisbelow(True)
    for ax in axes:
        ax.set_ylabel("Cluster compactness — within-cluster spread (Å, lower = tighter)")
    # Neutral grey solid-vs-hatched legend (the hue carries the tool; solid/hatched
    # carries raw vs filtered), figure-level between title and panels — same as the
    # other compare figures.
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=_CAT_LEGEND_GREY,
                             alpha=POSE_SET_STYLE[s]["alpha"], hatch=POSE_SET_STYLE[s]["hatch"],
                             edgecolor=_CAT_LEGEND_GREY, label=POSE_SET_STYLE[s]["label"])
               for s in ("raw", "filtered")]
    _label_panels(np.asarray(axes))
    fig.suptitle("Does filtering to PB-valid poses outside the pore tighten each tool's clusters? "
                 "Experimental vs Benchmark", fontsize=13, y=0.985)
    fig.legend(handles=handles, ncol=len(handles), fontsize=12, frameon=True,
               loc="upper center", bbox_to_anchor=(0.5, 0.93))
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = out_dir / "orai_cluster_quality_filtering_compare.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def write_quality_stats(q_exp: pd.DataFrame, q_bench: pd.DataFrame,
                        exp_label: str, bench_label: str, out_dir: Path,
                        pose_set: str = "raw") -> Optional[Path]:
    """Exp-vs-Benchmark per (tool, metric): Mann-Whitney U + Cliff's δ, Holm across
    the (tool × metric) family. Unpaired (different ligands), so effect sizes carry
    more weight than p at the tiny Exp n."""
    if su is None:
        return None
    L: List[str] = []
    L.append("Per-tool cluster QUALITY — Experimental (JKU) vs Benchmark ligands on Orai")
    L.append("=" * 74)
    L.append(f"pose set = '{pose_set}';  unit = (frame,ligand,tool) row (UNPAIRED across")
    L.append("datasets — different ligand sets). Exp n is tiny → EXPLORATORY; read the")
    L.append("Cliff's δ effect size, not just the p-value. Both use SMINA-refined DiffDock.")
    for lab, q in ((exp_label, q_exp), (bench_label, q_bench)):
        qs = q[q["pose_set"] == pose_set]
        L.append(f"  {lab:<28} {len(qs)} rows, {qs['ligand'].nunique()} ligands, "
                 f"tools {sorted(qs['tool'].unique())}")
    L.append("")
    tools = _tool_order(pd.concat([q_exp["tool"], q_bench["tool"]]))
    rows, pvals = [], []
    for key, short, _yl, hb in QUALITY_METRICS:
        for t in tools:
            a = pd.to_numeric(q_exp[(q_exp.pose_set == pose_set) & (q_exp.tool == t)][key],
                              errors="coerce").dropna().to_numpy()
            b = pd.to_numeric(q_bench[(q_bench.pose_set == pose_set) & (q_bench.tool == t)][key],
                              errors="coerce").dropna().to_numpy()
            if a.size < 3 or b.size < 3:
                rows.append((short, hb, t, None, a.size, b.size)); continue
            res = su.mannwhitney_cliffs(a, b)
            rows.append((short, hb, t, res, a.size, b.size)); pvals.append(res["p"])
    adj = su.holm(pvals) if pvals else []
    ai = 0
    L.append(f"-- Mann-Whitney U + Cliff's δ, Exp vs Bench (Holm across tool×metric) --")
    for short, hb, t, res, na, nb in rows:
        head = f"  {short:<20} {_tool_name(t):<16} ({'higher' if hb else 'lower'}=better)"
        if res is None:
            L.append(f"{head}  n/a (n_exp={na}, n_bench={nb})"); continue
        me, mb = res["medians"]; pa = adj[ai]; ai += 1
        L.append(f"{head}  med Exp={me:6.2f}  Bench={mb:6.2f}  "
                 f"Cliff δ={res['cliffs_delta']:+.2f}  "
                 f"p_holm={su.fmt_p(pa)} {su.p_stars(pa)}  (n {na} vs {nb})")
    L.append("")
    L.append("Panels carry descriptive content only; these tests are kept off the figures.")
    p = out_dir / "orai_cluster_quality_compare_stats.txt"
    p.write_text("\n".join(L) + "\n")
    return p


# ── panel labels (local fallback if the shared helper is unavailable) ──────

def _label_panels(axes, fontsize=12):
    try:
        from pocket_comparison_report import _label_panels as shared
        return shared(axes, fontsize)
    except Exception:
        flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
        for i, ax in enumerate(flat):
            if ax is not None:
                ax.text(-0.06, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                        fontsize=fontsize, fontweight="bold", va="bottom", ha="right")


def _parse_dataset(spec: str) -> Tuple[str, str, Path]:
    """'key:Label=path.csv' → (label, key, Path). key is 'exp' or 'bench'."""
    key, rest = spec.split(":", 1)
    label, path = rest.split("=", 1)
    return label.strip(), key.strip(), Path(path.strip())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exp-csv",
                    default="posebusters_results/orai_jku/pose_clusters/per_pair.csv",
                    help="per_pair.csv for the Experimental (JKU) ligand run.")
    ap.add_argument("--bench-csv",
                    default="posebusters_results/orai_benchmark/pose_clusters/per_pair.csv",
                    help="per_pair.csv for the Benchmark ligand run.")
    ap.add_argument("--exp-label", default="Exp. Ligands")
    ap.add_argument("--bench-label", default="Benchmark ligands")
    ap.add_argument("--exp-keep", default="",
                    help="Comma-separated ligand-name substrings (case-insensitive); keep "
                         "only Experimental rows whose 'ligand' matches one. Empty = keep "
                         "all. E.g. 'gsk7975a,synta' for the pore blockers only.")
    ap.add_argument("--bench-keep", default="",
                    help="Same, applied to the Benchmark dataset.")
    ap.add_argument("--suffix", default="",
                    help="Appended to every output filename stem (e.g. '_poreblockers'), "
                         "so a filtered run does not overwrite the full-set figures.")
    ap.add_argument("--thr", type=float, default=5.0,
                    help="Agree threshold (Å); overridden per dataset by its stats "
                         "sidecar match_thr_A when present.")
    ap.add_argument("--out-dir",
                    default="posebusters_results/orai_jku/pose_clusters/panels",
                    help="Where the *_compare figures are written (originals kept).")
    ap.add_argument("--exp-quality-csv", default=None,
                    help="cluster_quality_per_tool.csv for the Experimental run "
                         "(default: sibling of --exp-csv). Enables the cluster-QUALITY "
                         "overlay; skipped if absent.")
    ap.add_argument("--bench-quality-csv", default=None,
                    help="cluster_quality_per_tool.csv for the Benchmark run "
                         "(default: sibling of --bench-csv).")
    ap.add_argument("--quality-pose-set", default="raw", choices=("raw", "filtered"),
                    help="Which pose population the quality-overlay boxplots use.")
    ap.add_argument("--exclude-ligand", default=",".join(DEFAULT_EXCLUDE_LIGANDS),
                    help="Comma-separated ligand-name substrings (case-insensitive) to drop "
                         "from BOTH datasets and BOTH the agreement and quality analyses. "
                         f"Default drops the GSK-7975A deprotonated protomer "
                         f"({', '.join(DEFAULT_EXCLUDE_LIGANDS)}); pass '' to keep everything.")
    args = ap.parse_args(argv)
    exclude = [s.strip() for s in args.exclude_ligand.split(",") if s.strip()]

    exp_p, bench_p = Path(args.exp_csv), Path(args.bench_csv)
    missing = [str(p) for p in (exp_p, bench_p) if not p.exists()]
    if missing:
        print("ERROR: missing per_pair.csv (run the pose-cluster cells first):")
        for m in missing:
            print("  -", m)
        return 1

    datasets = [load_dataset(args.exp_label, "exp", exp_p, exclude),
                load_dataset(args.bench_label, "bench", bench_p, exclude)]
    # Optional ligand subsetting (e.g. pore blockers only). Applied per dataset.
    subset = bool(args.exp_keep.strip() or args.bench_keep.strip())
    datasets[0] = keep_ligands(datasets[0], args.exp_keep)
    datasets[1] = keep_ligands(datasets[1], args.bench_keep)
    # Prefer each dataset's own recorded threshold; they match (5 Å) but guard anyway.
    thrs = {ds["thr"] for ds in datasets if ds["thr"]}
    thr = args.thr if len(thrs) != 1 else thrs.pop()

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    suffix = args.suffix
    if exclude:
        print(f"  [exclude] dropping ligands matching (case-insensitive): {', '.join(exclude)}")
    print(f"Comparing cross-tool agreement (thr={thr:g} Å){' [subset]' if subset else ''}:")
    for ds in datasets:
        kept = f"  [{', '.join(ds['kept_ligands'])}]" if ds.get("kept_ligands") else ""
        print(f"  {ds['label']:<28} n={ds['n_pairs']:>5} pairs, {ds['n_ligands']} ligands{kept}")

    f1 = fig_distribution(datasets, thr, out_dir, suffix)
    print(f"  Figure: {f1}")
    # For a ligand subset (e.g. pore blockers) the default "Experimental vs Benchmark"
    # title is misleading, so derive it from the actual dataset labels instead.
    agr_title = (f"Cross-tool agreement on Orai — {datasets[0]['label']} vs "
                 f"{datasets[1]['label']}") if subset else None
    f2 = fig_agreement(datasets, thr, out_dir, suffix, agr_title)
    print(f"  Figure: {f2}")
    if su is not None:
        fs = write_stats(datasets, thr, out_dir, suffix)
        if fs:
            print(f"  Stats:  {fs}")
    else:
        print("  [stats] stats_utils unavailable — companion stats skipped.")

    # The cluster-QUALITY overlay reads a per-tool aggregate (cluster_quality_per_tool.csv)
    # that is NOT resolved per ligand, so it cannot be restricted to a ligand subset —
    # skip it whenever a --*-keep filter is active to avoid a mislabelled figure.
    if subset:
        print("  [quality] cluster-quality overlay skipped for ligand-subset runs "
              "(per-tool table is not ligand-resolved).")
        return 0

    # ── per-tool cluster-QUALITY overlay (only if the --cluster-quality tables exist) ──
    exp_q_csv = Path(args.exp_quality_csv) if args.exp_quality_csv else \
        exp_p.parent / "cluster_quality_per_tool.csv"
    bench_q_csv = Path(args.bench_quality_csv) if args.bench_quality_csv else \
        bench_p.parent / "cluster_quality_per_tool.csv"
    q_exp, q_bench = load_quality(exp_q_csv, exclude), load_quality(bench_q_csv, exclude)
    if q_exp is not None and q_bench is not None:
        print(f"Cluster-quality overlay (pose set '{args.quality_pose_set}'):")
        fq = fig_cluster_quality_compare(q_exp, q_bench, args.exp_label, args.bench_label,
                                         out_dir, args.quality_pose_set)
        print(f"  Figure: {fq}")
        ff = fig_cluster_quality_filtering_compare(q_exp, q_bench, args.exp_label,
                                                   args.bench_label, out_dir)
        if ff:
            print(f"  Figure: {ff}")
        fqs = write_quality_stats(q_exp, q_bench, args.exp_label, args.bench_label,
                                  out_dir, args.quality_pose_set)
        if fqs:
            print(f"  Stats:  {fqs}")
    else:
        miss = [str(p) for p, qq in ((exp_q_csv, q_exp), (bench_q_csv, q_bench)) if qq is None]
        print("  [quality] cluster_quality_per_tool.csv not found for: "
              + ", ".join(miss) + " — run the report with --cluster-quality first; "
              "quality overlay skipped.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
