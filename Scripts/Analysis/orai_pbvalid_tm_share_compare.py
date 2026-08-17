#!/usr/bin/env python3
"""Cross-dataset comparison of PoseBusters validity and transmembrane survival.

This is an ADDITIVE companion to orai_transmembrane_exclusion.py. That script
draws the per-dataset ``fig_tool_tm_loss.png`` (absolute pose counts, one figure
per dataset), which cannot be read across datasets because Orai × Benchmark has
~100× more poses than Orai × Experimental. Here every quantity is a *share of
produced poses*, so the two ligand sets sit on one axis and are directly
comparable.

It reads the already-TM-classified per-pose tables written by the exclusion run
(``<dataset>/transmembrane_filter/tm_pose_classification.csv``) — no geometry is
recomputed, so the numbers agree by construction with tool_tm_loss_summary.csv /
group_toolchain_tm_comparison.csv.

Produces (in posebusters_results/orai_pbvalid_tm_share_compare/):
  fig_pbvalid_share_compare.png     (A) PB-valid poses as % of produced, per tool,
                                        Benchmark vs Experimental (grouped bars);
                                    (B) same but only PB-valid poses OUTSIDE the TM
                                        (usable yield after the TM filter).
  fig_tm_fate_share_compare.png     100%-normalised fate stack per tool × dataset
                                    (fig_tool_tm_loss re-expressed as shares) so the
                                    outside-TM / inside-TM / off-protein / invalid
                                    split is comparable between ligand sets.
  fig_tm_loss_relative_compare.png  the direct cross-dataset difference in the
                                    'lost to TM %' headline of fig_tool_tm_loss:
                                    (A) dumbbell of the TM-loss rate (inside-TM /
                                    PB-valid) per tool, Experimental vs Benchmark,
                                    Wilson CIs; (B) diverging bar of Δ = Benchmark −
                                    Experimental in PERCENTAGE POINTS (both terms are
                                    proportions, so the gap is pp, not a ratio).
  fig_pbvalid_outside_tm_whisker.png per-(frame × ligand) distribution of the
                                    PB-valid-and-outside-TM share — the honest
                                    (pseudoreplication-free) view of panel (B).
  pbvalid_tm_share_pooled.csv       pooled counts + shares + Wilson CIs.
  pbvalid_tm_share_per_unit.csv     per (dataset, tool, frame, ligand) shares.
  pbvalid_tm_share_compare_stats.txt exact numbers + statistical tests + caveats.

Inferential test results live only in the txt sidecar; the figures carry
descriptive content (shares, Wilson CIs, n) — matching the project convention.
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

# Reuse the exclusion module's palettes/labels so colours + names stay in lock-step
# with fig_tool_tm_loss / the survival box. Fall back to local copies if the heavy
# import is unavailable (e.g. run outside the vina env).
try:
    from orai_transmembrane_exclusion import (
        _FATE, CATEGORY_COLOUR, CATEGORY_LABELS, CATEGORY_ORDER, _label_panels,
        _toolchain_label,
    )
except Exception:                                   # pragma: no cover
    _FATE = [
        ("survive",  "PB-valid · outside TM (kept)",       "#55A868", ""),
        ("lost_tm",  "PB-valid · inside TM (lost)",        "#C44E52", "///"),
        ("unplaced", "PB-valid · off-protein/unevaluable", "#DD8452", ""),
        ("nonvalid", "not PB-valid",                       "#CFCFCF", ""),
    ]
    CATEGORY_LABELS = {"orai_jku": "Exp. Ligands",
                       "orai_benchmark": "Benchmark ligands × Orai"}
    CATEGORY_ORDER = ["orai_jku", "orai_benchmark"]
    CATEGORY_COLOUR = {"Exp. Ligands": "#8172B3",
                       "Benchmark ligands × Orai": "#64B5CD"}

    def _label_panels(axes, fontsize: int = 13) -> None:
        import string
        flat = list(axes) if isinstance(axes, (list, tuple)) else [axes]
        for i, ax in enumerate([a for a in flat if a is not None]):
            ax.set_title(f"({string.ascii_uppercase[i]})", loc="left",
                         fontweight="bold", fontsize=fontsize)

    def _toolchain_label(base, variant):            # mirror of the exclusion helper
        name = {"autodock": "AutoDock Vina", "diffdock": "DiffDock",
                "equibind": "EquiBind"}.get(base, str(base).title())
        if base == "autodock" or variant in ("none", ""):
            return name
        if variant == "gnina":
            return f"{name}*"
        if variant == "raw":
            return name
        return f"{name} ({variant})"

try:
    import stats_utils as su
    _HAVE_SU = True
except Exception:                                   # pragma: no cover
    su = None
    _HAVE_SU = False

# Per-tool colour palette for the whisker figure. Reuse posebusters_pose_comparison's
# TOOL_COLORS (AutoDock blue / DiffDock orange / EquiBind green) so the whisker matches
# the 09f_pbvalid_yield_dominant figure exactly; fall back to the current hex values
# (matplotlib tab palette) if the heavy import is unavailable.
try:
    from posebusters_pose_comparison import TOOL_COLORS as _PPC_TOOL_COLORS
except Exception:                                   # pragma: no cover
    _PPC_TOOL_COLORS = None
_TOOL_COLOUR_KEY = {"AutoDock": "autodock", "DiffDock": "diffdock",
                    "EquiBind": "equibind_unguided"}
_TOOL_COLOUR_FALLBACK = {"autodock": "#1f77b4", "diffdock": "#ff7f0e",
                         "equibind_unguided": "#2ca02c"}


def _tool_colour(short_name: str) -> str:
    """Tool colour keyed on the x-axis short name (AutoDock/DiffDock/EquiBind)."""
    key = _TOOL_COLOUR_KEY.get(short_name, "")
    if _PPC_TOOL_COLORS is not None:
        c = _PPC_TOOL_COLORS.get(key)
        if c:
            return c
    return _TOOL_COLOUR_FALLBACK.get(key, "#888888")


# Per-dataset FILL STYLE + short legend label for the whisker figure (fill colour encodes
# the tool, so the ligand SET is told apart by solid-vs-hatched, matching 09f's raw/opt
# convention). 'Benchmark ligands × Orai' is shortened to just 'Benchmark' in the legend.
_DATASET_STYLE = {
    "orai_benchmark": {"alpha": 0.80, "hatch": "",    "label": "Benchmark",
                       "full": "Orai × Benchmark ligands"},
    "orai_jku":       {"alpha": 0.45, "hatch": "///", "label": "Experimental",
                       "full": "Orai × Experimental ligands"},
}

# The reported toolchain per tool. Resolved from the DiffDock/EquiBind refine
# variants so this stays in lock-step with the parent exclusion run's selection
# (tool_tm_loss_summary.csv). Defaults = the current pipeline pick (smina-refined
# DiffDock, gnina-refined EquiBind; AutoDock has no refinement). main() overrides
# both globals from --diffdock-variant / --equibind-variant.
SELECTED_TOOLCHAINS = ["AutoDock Vina", "DiffDock (smina)", "EquiBind*"]
# Short display name per selected toolchain (x-tick on the grouped-bar figures).
TOOL_SHORT = {"AutoDock Vina": "AutoDock", "DiffDock (smina)": "DiffDock",
              "EquiBind*": "EquiBind"}


def _resolve_toolchains(dd_variant: str, eb_variant: str):
    """(selected-toolchain-labels, short-name map) for the given refine variants.
    'all'/'original' fall back to the reported pick (smina DiffDock, gnina EquiBind)
    since this comparison needs exactly one variant per tool."""
    dd = "smina" if str(dd_variant) in ("all", "original") else str(dd_variant)
    eb = "gnina" if str(eb_variant) in ("all", "original") else str(eb_variant)
    ad_l = "AutoDock Vina"
    dd_l = _toolchain_label("diffdock", dd)
    eb_l = _toolchain_label("equibind", eb)
    return [ad_l, dd_l, eb_l], {ad_l: "AutoDock", dd_l: "DiffDock", eb_l: "EquiBind"}

_KEPT = "kept"
_TM = "transmembrane"
_OFF = ("invalid_placement", "unevaluable")

# fate key -> label/colour/hatch (order = stack order, valid-and-kept at the bottom)
_FATE_ORDER = [k for k, *_ in _FATE]
_FATE_LABEL = {k: lab for k, lab, _c, _h in _FATE}
_FATE_COLOUR = {k: c for k, _lab, c, _h in _FATE}
_FATE_HATCH = {k: h for k, _lab, _c, h in _FATE}


# ---------------------------------------------------------------------------
def _classification_path(dataset: str, results_root: Path) -> Path:
    return (results_root / dataset / "transmembrane_filter" /
            "tm_pose_classification.csv")


def _coerce_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0"))


def _fate_of(status: pd.Series, pb: pd.Series) -> pd.Series:
    """survive / lost_tm / unplaced / nonvalid — identical semantics to _prep_fate."""
    out = np.where(~pb, "nonvalid",
          np.where(status == _KEPT, "survive",
          np.where(status == _TM, "lost_tm", "unplaced")))
    return pd.Series(out, index=status.index)


def load_dataset(dataset: str, results_root: Path, top_n: int = 0) -> pd.DataFrame:
    """Selected-toolchain per-pose rows tagged with pb_valid + fate for one dataset.

    ``top_n`` > 0 equalises the pose budget across tools by keeping only each tool's
    top-N ranked poses per (frame × ligand): pose_rank is the tool's native ranking
    (AutoDock = Vina mode rank, DiffDock = confidence rank, EquiBind = refinement-energy
    rank). Poses with an unrankable (null) pose_rank are dropped. This only removes
    poses where a tool produced more than N per unit (here: Orai × Experimental AutoDock
    & EquiBind at 30/unit); tools already ≤ N (DiffDock, all Benchmark) are untouched.
    """
    path = _classification_path(dataset, results_root)
    if not path.exists():
        raise FileNotFoundError(f"missing classification table: {path}")
    df = pd.read_csv(path, low_memory=False)
    df = df[df["toolchain"].isin(SELECTED_TOOLCHAINS)].copy()
    if df.empty:
        raise ValueError(f"no selected toolchains present in {path}")
    if top_n and "pose_rank" in df.columns:
        pr = pd.to_numeric(df["pose_rank"], errors="coerce")
        df = df[pr.le(top_n)].copy()
    df["pb"] = _coerce_bool(df["pb_valid"])
    df["fate"] = _fate_of(df["analysis_status"], df["pb"])
    df["dataset"] = dataset
    df["category"] = CATEGORY_LABELS.get(dataset, dataset)
    return df


# ---------------------------------------------------------------------------
def pooled_table(frames: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per (dataset, toolchain): produced/valid/kept/tm/off counts + shares
    + Wilson CIs on the two key shares (valid/produced and kept/produced)."""
    rows = []
    for ds in CATEGORY_ORDER:
        df = frames.get(ds)
        if df is None:
            continue
        for tc in SELECTED_TOOLCHAINS:
            s = df[df["toolchain"] == tc]
            produced = int(len(s))
            if produced == 0:
                continue
            valid = int(s["pb"].sum())
            kept = int(((s["fate"] == "survive")).sum())
            tm = int((s["fate"] == "lost_tm").sum())
            off = int((s["fate"] == "unplaced").sum())
            nonvalid = int((s["fate"] == "nonvalid").sum())
            vlo, vhi = su.wilson_ci(valid, produced) if _HAVE_SU else (np.nan, np.nan)
            klo, khi = su.wilson_ci(kept, produced) if _HAVE_SU else (np.nan, np.nan)
            rows.append(dict(
                dataset=ds, category=CATEGORY_LABELS.get(ds, ds), toolchain=tc,
                tool=TOOL_SHORT.get(tc, tc), produced=produced, pb_valid=valid,
                outside_tm=kept, inside_tm=tm, off_protein_uneval=off, non_valid=nonvalid,
                pbvalid_share=valid / produced, pbvalid_share_lo=vlo, pbvalid_share_hi=vhi,
                outside_tm_share=kept / produced, outside_tm_share_lo=klo, outside_tm_share_hi=khi,
            ))
    return pd.DataFrame(rows)


def per_unit_table(frames: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """One row per (dataset, toolchain, frame, ligand): produced/valid/kept counts
    and the two per-unit shares. The (frame × ligand) pair is the independent unit."""
    rows = []
    for ds in CATEGORY_ORDER:
        df = frames.get(ds)
        if df is None:
            continue
        for (tc, prot, lig), s in df.groupby(["toolchain", "protein", "ligand"]):
            produced = int(len(s))
            if produced == 0:
                continue
            valid = int(s["pb"].sum())
            kept = int((s["fate"] == "survive").sum())
            rows.append(dict(
                dataset=ds, category=CATEGORY_LABELS.get(ds, ds), toolchain=tc,
                tool=TOOL_SHORT.get(tc, tc), frame=prot, ligand=lig,
                produced=produced, pb_valid=valid, outside_tm=kept,
                pbvalid_share=valid / produced, outside_tm_share=kept / produced,
            ))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
def _bar_panel(ax, pooled: pd.DataFrame, value: str, lo: str, hi: str, ylabel: str):
    """Grouped bars: x = tool, one bar per dataset, y = share (%). Wilson CI whiskers."""
    tools = [TOOL_SHORT[t] for t in SELECTED_TOOLCHAINS]
    x = np.arange(len(tools), dtype=float)
    width = 0.38
    for j, ds in enumerate(CATEGORY_ORDER):
        cat = CATEGORY_LABELS.get(ds, ds)
        colour = CATEGORY_COLOUR.get(cat, "#888888")
        offs = (j - 0.5) * width
        ys, elo, ehi = [], [], []
        for tc in SELECTED_TOOLCHAINS:
            r = pooled[(pooled.dataset == ds) & (pooled.toolchain == tc)]
            if len(r):
                r = r.iloc[0]
                ys.append(100 * r[value]); elo.append(100 * (r[value] - r[lo]))
                ehi.append(100 * (r[hi] - r[value]))
            else:
                ys.append(np.nan); elo.append(0); ehi.append(0)
        bars = ax.bar(x + offs, ys, width, color=colour, edgecolor="black",
                      linewidth=0.6, label=cat, zorder=2)
        ax.errorbar(x + offs, ys, yerr=[elo, ehi], fmt="none", ecolor="#333333",
                    elinewidth=1.0, capsize=3, zorder=3)
        for xi, yi in zip(x + offs, ys):
            if not np.isnan(yi):
                ax.text(xi, yi + 1.5, f"{yi:.1f}", ha="center", va="bottom",
                        fontsize=8.5, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(tools)
    ax.set_ylim(0, 108); ax.set_ylabel(ylabel)
    ax.grid(axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)


def fig_share_compare(pooled: pd.DataFrame, out: Path):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6))
    _bar_panel(axes[0], pooled, "pbvalid_share", "pbvalid_share_lo", "pbvalid_share_hi",
               "PoseBusters-valid poses\n(% of produced poses)")
    axes[0].set_title("PoseBusters-valid share of produced poses", fontsize=11)
    _bar_panel(axes[1], pooled, "outside_tm_share", "outside_tm_share_lo",
               "outside_tm_share_hi",
               "PoseBusters-valid poses outside the\ntransmembrane (% of produced poses)")
    axes[1].set_title("PoseBusters-valid AND outside the transmembrane\n"
                      "(usable yield after the transmembrane filter)", fontsize=11)
    _label_panels([axes[0], axes[1]])
    handles = [Patch(facecolor=CATEGORY_COLOUR.get(CATEGORY_LABELS[ds], "#888"),
                     edgecolor="black", label=CATEGORY_LABELS[ds]) for ds in CATEGORY_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.005))
    fig.text(0.5, -0.03, "Error bars are Wilson 95% CIs on the pooled pose count "
             "(pose-level, pseudoreplicated) — for generalisation across the small "
             "Experimental set,\nread the per-(frame × ligand) whisker figure and the "
             "statistics sidecar, not the width of these bars.",
             ha="center", va="top", fontsize=8, color="#555555")
    fig.suptitle("Orai docking — PoseBusters validity of produced poses, per tool "
                 "(error bars: Wilson 95% CI on the pooled pose count)",
                 fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def fig_fate_stack(pooled: pd.DataFrame, out: Path):
    """100%-normalised fate stack per tool × dataset (fig_tool_tm_loss as shares)."""
    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    tools = SELECTED_TOOLCHAINS
    x = np.arange(len(tools), dtype=float)
    width = 0.38
    for j, ds in enumerate(CATEGORY_ORDER):
        offs = (j - 0.5) * width
        for tc, xi in zip(tools, x + offs):
            r = pooled[(pooled.dataset == ds) & (pooled.toolchain == tc)]
            if not len(r):
                continue
            r = r.iloc[0]
            produced = r["produced"]
            seg = {"survive": r["outside_tm"], "lost_tm": r["inside_tm"],
                   "unplaced": r["off_protein_uneval"], "nonvalid": r["non_valid"]}
            bottom = 0.0
            for k in _FATE_ORDER:
                frac = 100 * seg[k] / produced
                ax.bar(xi, frac, width, bottom=bottom, color=_FATE_COLOUR[k],
                       hatch=_FATE_HATCH[k] or None, edgecolor="white", linewidth=0.4,
                       zorder=2)
                if frac >= 6:
                    ax.text(xi, bottom + frac / 2, f"{frac:.0f}", ha="center",
                            va="center", fontsize=7.5, color="#222222")
                bottom += frac
            ax.text(xi, 101.5, CATEGORY_LABELS[ds].split()[0], ha="center",
                    va="bottom", fontsize=8, rotation=0, color="#555555")
            ax.text(xi, 104, f"n={produced}", ha="center", va="bottom", fontsize=7.5,
                    color="#555555")
    ax.set_xticks(x); ax.set_xticklabels([TOOL_SHORT[t] for t in tools])
    ax.set_ylim(0, 110); ax.set_ylabel("Share of produced poses (%)")
    ax.set_title("Fate of every produced pose, per tool — Experimental vs Benchmark ligands\n"
                 "left bar of each pair = Experimental, right = Benchmark   "
                 "(PB = PoseBusters-valid, TM = transmembrane)",
                 fontsize=11.5, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    handles = [Patch(facecolor=_FATE_COLOUR[k], hatch=_FATE_HATCH[k] or None,
                     edgecolor="white", label=_FATE_LABEL[k]) for k in _FATE_ORDER]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.18),
              ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def fig_whisker(per_unit: pd.DataFrame, out: Path):
    """Per-(frame × ligand) distribution of the PB-valid-and-outside-TM share.

    Styled to match 09f_pbvalid_yield_dominant: box FILL encodes the TOOL (AutoDock blue /
    DiffDock orange / EquiBind green), and the two ligand SETS within a tool differ only by
    fill STYLE (Benchmark solid, Experimental hatched). The SET legend is a neutral-grey
    solid-vs-hatched pair placed above the axes."""
    fig, ax = plt.subplots(figsize=(11.5, 6.4))
    tools = SELECTED_TOOLCHAINS
    x = np.arange(len(tools), dtype=float)
    width = 0.34
    rng = np.random.default_rng(0)
    for j, ds in enumerate(CATEGORY_ORDER):
        style = _DATASET_STYLE.get(ds, {"alpha": 0.55, "hatch": ""})
        offs = (j - 0.5) * width
        for tc, xi in zip(tools, x + offs):
            colour = _tool_colour(TOOL_SHORT.get(tc, tc))
            vals = 100 * per_unit[(per_unit.dataset == ds) &
                                  (per_unit.toolchain == tc)]["outside_tm_share"].values
            if len(vals) == 0:
                continue
            bp = ax.boxplot([vals], positions=[xi], widths=width * 0.9,
                            patch_artist=True, showfliers=False, whis=(5, 95),
                            medianprops=dict(color="black", linewidth=1.4),
                            manage_ticks=False)
            for box in bp["boxes"]:
                box.set(facecolor=colour, alpha=style["alpha"], edgecolor=colour,
                        linewidth=1.1)
                if style["hatch"]:
                    box.set_hatch(style["hatch"])
            jit = xi + (rng.random(len(vals)) - 0.5) * width * 0.55
            ax.scatter(jit, vals, s=8, color=colour, edgecolor="black",
                       linewidth=0.2, alpha=0.55, zorder=3)
            ax.text(xi, -6, f"n={len(vals)}", ha="center", va="top", fontsize=11,
                    color="#555555")
    ax.set_xticks(x); ax.set_xticklabels([TOOL_SHORT[t] for t in tools])
    ax.set_ylim(-10, 105)
    ax.set_ylabel("PoseBusters-valid poses outside the transmembrane\n"
                  "(% of that unit's produced poses)")
    ax.grid(axis="y", alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    # SET legend: fill colour already encodes the tool, so the ligand set is shown as a
    # neutral-grey solid-vs-hatched pair. Sits above the axes (between title and plot),
    # matching 09f_pbvalid_yield_dominant's legend placement.
    handles = [Patch(facecolor="0.5", alpha=_DATASET_STYLE[ds]["alpha"], edgecolor="0.5",
                     hatch=_DATASET_STYLE[ds]["hatch"] or None,
                     label=_DATASET_STYLE[ds]["label"]) for ds in CATEGORY_ORDER]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.01),
              ncol=len(handles), frameon=False, borderaxespad=0.0, fontsize=13)
    fig.tight_layout()
    fig.subplots_adjust(top=0.85)
    fig.suptitle("PoseBusters-valid poses surviving outside the transmembrane, "
                 "per (frame × ligand) unit\n"
                 "each point = one frame × ligand; box = median / IQR / 5–95% whiskers",
                 fontsize=12, fontweight="bold", y=0.995)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
def _loss_rate_rows(pooled: pd.DataFrame) -> List[dict]:
    """Per tool: TM-loss RATE = inside_tm / pb_valid for each dataset, with Wilson CI.
    This is exactly the 'lost to TM %' headline of fig_tool_tm_loss (a share of the
    PoseBusters-valid poses), so the two datasets sit on one comparable axis."""
    rows: List[dict] = []
    for tc in SELECTED_TOOLCHAINS:
        rec = {"toolchain": tc, "tool": TOOL_SHORT.get(tc, tc), "by_ds": {}}
        for ds in CATEGORY_ORDER:
            r = pooled[(pooled.dataset == ds) & (pooled.toolchain == tc)]
            if not len(r):
                continue
            r = r.iloc[0]
            valid, tm = int(r["pb_valid"]), int(r["inside_tm"])
            if valid == 0:
                continue
            lo, hi = su.wilson_ci(tm, valid) if _HAVE_SU else (np.nan, np.nan)
            rec["by_ds"][ds] = dict(rate=tm / valid, lo=lo, hi=hi, tm=tm, valid=valid)
        rows.append(rec)
    return rows


def fig_tm_loss_relative(pooled: pd.DataFrame, out: Path):
    """Cross-dataset difference in the TM-loss rate of PB-valid poses.

    (A) dumbbell — per tool, the loss rate (inside-TM / PB-valid) for each dataset as
        two connected dots on a shared % axis, with Wilson 95% CI whiskers, and the gap
        annotated in percentage points.
    (B) diverging bar — Δ = Benchmark − Experimental (percentage points) per tool.

    The difference is shown in percentage POINTS, not a relative-% ratio: both terms are
    already proportions, and the Experimental base rests on ~3 ligands (≤12 units) so a
    ratio would be unstable/misleading. The ratio is reported (caveated) in the sidecar.
    """
    rows = _loss_rate_rows(pooled)
    if not rows:
        print("  [warn] fig_tm_loss_relative: no loss-rate rows — skipping")
        return
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(9.2, 9.6))
    y = np.arange(len(rows))[::-1]                    # first tool (AutoDock) at the top
    ylim = (-0.6, len(rows) - 0.4)

    # Dot / bar FILL encodes the TOOL (AutoDock blue / DiffDock orange / EquiBind green,
    # posebusters_pose_comparison.TOOL_COLORS); the two ligand SETS are told apart by fill
    # STYLE (Benchmark solid, Experimental hatched) — the same palette + legend convention as
    # the 09f_pbvalid_yield_dominant figure. Tool identity is also read off the y-axis.
    # ── Panel A: dumbbell of loss% per tool; one dot per dataset (solid vs hatched) ──
    for yi, r in zip(y, rows):
        col = _tool_colour(r["tool"])
        e = r["by_ds"].get("orai_jku"); b = r["by_ds"].get("orai_benchmark")
        if e and b:
            xe, xb = 100 * e["rate"], 100 * b["rate"]
            axA.plot([xe, xb], [yi, yi], color="0.6", lw=2.2, zorder=1)
            axA.annotate(f"Δ {xb - xe:+.1f} pp", ((xe + xb) / 2, yi + 0.17),
                         ha="center", va="bottom", fontsize=8.5, color="#333333")
        present = [(ds, r["by_ds"][ds]) for ds in CATEGORY_ORDER if r["by_ds"].get(ds)]
        for ds, p in present:
            style = _DATASET_STYLE.get(ds, {"alpha": 0.6, "hatch": ""})
            x = 100 * p["rate"]
            axA.errorbar([x], [yi],
                         xerr=[[100 * (p["rate"] - p["lo"])], [100 * (p["hi"] - p["rate"])]],
                         fmt="none", ecolor=col, elinewidth=1.3, capsize=3, zorder=2)
            # Alpha baked into the FILL so the white edge/hatch stay opaque — the crisp
            # ring separates the two dots on rows where they nearly coincide (EquiBind,
            # ~1.5 pp apart). Fill hue still encodes the tool; solid vs hatched the set.
            axA.scatter([x], [yi], s=200, facecolor=to_rgba(col, style["alpha"]),
                        edgecolor="white", linewidth=1.3, hatch=style["hatch"] or None,
                        zorder=3)
        # value labels below the dots; push outward so near-coincident dots
        # (e.g. EquiBind, ~1.5 pp apart) don't overprint each other.
        if len(present) == 2:
            (dsL, pL), (dsR, pR) = sorted(present, key=lambda t: t[1]["rate"])
            axA.text(100 * pL["rate"], yi - 0.26, f"{100*pL['rate']:.1f}%", ha="right",
                     va="top", fontsize=8, fontweight="bold", color="#222222")
            axA.text(100 * pR["rate"], yi - 0.26, f"{100*pR['rate']:.1f}%", ha="left",
                     va="top", fontsize=8, fontweight="bold", color="#222222")
        else:
            for ds, p in present:
                axA.text(100 * p["rate"], yi - 0.26, f"{100*p['rate']:.1f}%", ha="center",
                         va="top", fontsize=8, fontweight="bold", color="#222222")
    axA.set_yticks(y); axA.set_yticklabels([r["tool"] for r in rows])
    axA.set_ylim(*ylim); axA.set_xlim(-4, 108)
    axA.set_xlabel("PoseBusters-valid poses lost to the transmembrane\n"
                   "(% of PoseBusters-valid poses)")
    axA.set_title("Transmembrane loss per tool — Benchmark vs Experimental", fontsize=11)
    axA.grid(axis="x", alpha=0.3, zorder=0); axA.set_axisbelow(True)

    # ── Panel B: diverging bar of Δ (Benchmark − Experimental) in percentage points ──
    #    bars keep the TOOL colour; direction (which set loses more) is read off the zero
    #    line and the arrowed axis label, not a separate red/green scale.
    deltas, bcolours = [], []
    for r in rows:
        e = r["by_ds"].get("orai_jku"); b = r["by_ds"].get("orai_benchmark")
        deltas.append(100 * (b["rate"] - e["rate"]) if (e and b) else np.nan)
        bcolours.append(_tool_colour(r["tool"]))
    axB.barh(y, deltas, color=bcolours, edgecolor="black", linewidth=0.6, height=0.55,
             zorder=2)
    axB.axvline(0, color="#333333", lw=1.0, zorder=1)
    for yi, d in zip(y, deltas):
        if np.isnan(d):
            continue
        axB.text(d + (1.3 if d >= 0 else -1.3), yi, f"{d:+.1f} pp", va="center",
                 ha="left" if d >= 0 else "right", fontsize=9.5, fontweight="bold")
    m = np.nanmax(np.abs(deltas)) if np.any(~np.isnan(deltas)) else 1.0
    # Stacked below A (no longer sharing A's y-axis alongside it), so B carries its own
    # tool labels rather than borrowing A's.
    axB.set_yticks(y); axB.set_yticklabels([r["tool"] for r in rows])
    axB.set_ylim(*ylim); axB.set_xlim(-1.6 * m - 4, 1.6 * m + 4)
    axB.set_xlabel("Δ transmembrane loss = Benchmark − Experimental (percentage points)\n"
                   "← Benchmark loses fewer          Benchmark loses more →", fontsize=9)
    axB.set_title("Relative difference across datasets", fontsize=11)
    axB.grid(axis="x", alpha=0.3, zorder=0); axB.set_axisbelow(True)

    _label_panels([axA, axB])
    # Explicit margins (not tight_layout) so the single-line suptitle + legend sit flush
    # above the axes with no dead band; hspace holds A's two-line xlabel above B's title.
    fig.subplots_adjust(left=0.14, right=0.96, top=0.88, bottom=0.08, hspace=0.42)
    # SET legend (matches 09f_pbvalid_yield_dominant): fill colour already encodes the tool
    # (read off the y-axis), so the ligand set is a neutral-grey solid-vs-hatched pair in a
    # single row between the suptitle and the panel titles — Benchmark first, as in 09f.
    ds_handles = [Patch(facecolor="0.5", alpha=_DATASET_STYLE[ds]["alpha"], edgecolor="0.5",
                        hatch=_DATASET_STYLE[ds]["hatch"] or None,
                        label=_DATASET_STYLE[ds]["full"])
                  for ds in ("orai_benchmark", "orai_jku")]
    fig.legend(handles=ds_handles, loc="upper center", ncol=2, framealpha=0.93, fontsize=9,
               bbox_to_anchor=(0.5, 0.935), borderaxespad=0.0, handlelength=1.6,
               columnspacing=1.6)
    fig.suptitle("Orai docking — cross-dataset difference in transmembrane loss of "
                 "PoseBusters-valid poses, per tool",
                 fontsize=12.5, fontweight="bold", y=0.985)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    # Caption/method note lives in a companion txt (project convention: keep the figure
    # descriptive, push the footnote off-panel). Written next to the png as <stem>.txt.
    note = (
        "fig_tm_loss_relative_compare.png — caption / method note\n"
        "Orai docking — cross-dataset difference in transmembrane loss of "
        "PoseBusters-valid poses, per tool.\n"
        "loss = poses inside the transmembrane / PoseBusters-valid poses; the gap is in "
        "percentage points (both terms are proportions).\n\n"
        "Dot / bar colour encodes the docking tool (AutoDock blue, DiffDock orange, "
        "EquiBind green); the two ligand sets are told apart by fill (Benchmark solid, "
        "Experimental hatched).\n"
        "Dots are pooled pose-level rates with Wilson 95% CIs; the difference is in "
        "percentage POINTS, not a relative-% ratio (the Experimental base rests on ~3 "
        "ligands / ≤12 units, so a ratio would be unstable — it is reported, "
        "caveated, in the statistics sidecar pbvalid_tm_share_compare_stats.txt).\n"
    )
    note_path = out.with_suffix(".txt")
    note_path.write_text(note)
    print(f"  wrote {out}")
    print(f"  wrote {note_path}")


# ---------------------------------------------------------------------------
def _fmt_pct(x):
    return f"{100 * x:5.1f}%"


def write_stats(pooled: pd.DataFrame, per_unit: pd.DataFrame, out: Path,
                data_source_note: str):
    L: List[str] = []
    A = L.append
    A("PoseBusters validity & transmembrane survival — Orai × Benchmark vs Orai × Experimental")
    A("#" * 88)
    A("")
    A("Metric definitions (denominator = poses the tool PRODUCED):")
    A("  PB-valid share            = PoseBusters-valid poses / produced poses")
    A("  outside-TM (usable) share = (PoseBusters-valid AND outside the transmembrane) / produced poses")
    A("Selected toolchain per tool (matches tool_tm_loss_summary.csv): "
      + ", ".join(SELECTED_TOOLCHAINS) + ".")
    A(data_source_note)
    A("")
    A("=" * 88)
    A("1. POOLED counts and shares (all produced poses)")
    A("=" * 88)
    for ds in CATEGORY_ORDER:
        A("")
        A(f"── {CATEGORY_LABELS.get(ds, ds)}  ({ds}) ──")
        A(f"  {'tool':10s} {'produced':>8s} {'pb_valid':>9s} {'valid%':>7s} {'valid% 95%CI':>15s}"
          f" {'outTM':>7s} {'outTM%':>7s} {'outTM% 95%CI':>15s}")
        for tc in SELECTED_TOOLCHAINS:
            r = pooled[(pooled.dataset == ds) & (pooled.toolchain == tc)]
            if not len(r):
                continue
            r = r.iloc[0]
            vci = f"[{100*r.pbvalid_share_lo:.1f}-{100*r.pbvalid_share_hi:.1f}]"
            kci = f"[{100*r.outside_tm_share_lo:.1f}-{100*r.outside_tm_share_hi:.1f}]"
            A(f"  {r.tool:10s} {r.produced:8d} {r.pb_valid:9d} {_fmt_pct(r.pbvalid_share):>7s}"
              f" {vci:>15s} {r.outside_tm:7d} {_fmt_pct(r.outside_tm_share):>7s} {kci:>15s}")
        A(f"  (fate split per tool: inside-TM / off-protein-or-unevaluable also in "
          f"pbvalid_tm_share_pooled.csv)")
    A("")
    A("=" * 88)
    A("2. Benchmark vs Experimental, per tool")
    A("=" * 88)
    A("")
    A("2a. Pooled two-proportion Fisher exact test (POSE-LEVEL — pseudoreplicated, see caveat).")
    A("    Δ = Benchmark rate − Experimental rate (percentage points). Poses within a")
    A("    (frame × ligand) unit are correlated, so pose-level p-values OVERSTATE")
    A("    significance; they are reported only as a descriptive pooled contrast.")
    for metric, num, den_valid in [("PB-valid share", "pb_valid", None),
                                   ("outside-TM share", "outside_tm", None)]:
        A(f"    · {metric}:")
        praw = []
        pairs = []
        for tc in SELECTED_TOOLCHAINS:
            rb = pooled[(pooled.dataset == "orai_benchmark") & (pooled.toolchain == tc)]
            re = pooled[(pooled.dataset == "orai_jku") & (pooled.toolchain == tc)]
            if not len(rb) or not len(re):
                continue
            rb, re = rb.iloc[0], re.iloc[0]
            kb, nb = int(rb[num]), int(rb["produced"])
            ke, ne = int(re[num]), int(re["produced"])
            if _HAVE_SU:
                res = su.two_proportion_test(kb, nb, ke, ne)
                p = res.get("p"); rd = res.get("risk_diff")   # bench_rate - exp_rate
            else:
                p, rd = float("nan"), float("nan")
            pairs.append((tc, kb, nb, ke, ne, rd, p)); praw.append(p)
        pholm = su.holm(praw) if (_HAVE_SU and praw) else praw
        for (tc, kb, nb, ke, ne, rd, p), ph in zip(pairs, pholm):
            star = su.p_stars(ph) if _HAVE_SU else ""
            pr = su.fmt_p(p) if _HAVE_SU else f"p={p}"
            phs = su.fmt_p(ph) if _HAVE_SU else f"p={ph}"
            A(f"        {TOOL_SHORT.get(tc, tc):9s} Bench {kb}/{nb} ({100*kb/nb:.1f}%) vs "
              f"Exp {ke}/{ne} ({100*ke/ne:.1f}%)  Δ={100*rd:+.1f}pp  Fisher {pr} (Holm {phs}) {star}")
    A("")
    A("2b. Per-(frame × ligand) unit Mann–Whitney U + Cliff's δ (pseudoreplication-free).")
    A("    δ>0 ⇒ Benchmark unit-shares higher than Experimental; the (frame × ligand) pair")
    A("    is the independent unit. Experimental has only 3 ligands × 4 frames → EXPLORATORY.")
    for metric, col in [("PB-valid share", "pbvalid_share"),
                        ("outside-TM share", "outside_tm_share")]:
        A(f"    · {metric}:")
        praw, pairs = [], []
        for tc in SELECTED_TOOLCHAINS:
            b = per_unit[(per_unit.dataset == "orai_benchmark") &
                         (per_unit.toolchain == tc)][col].values
            e = per_unit[(per_unit.dataset == "orai_jku") &
                         (per_unit.toolchain == tc)][col].values
            if len(b) == 0 or len(e) == 0:
                continue
            if _HAVE_SU:
                mw = su.mannwhitney_cliffs(b, e)     # Bench vs Exp
                U, p, delta = mw.get("U"), mw.get("p"), mw.get("cliffs_delta")
            else:
                U, p, delta = float("nan"), float("nan"), float("nan")
            pairs.append((tc, len(b), len(e), np.median(b), np.median(e), U, delta, p))
            praw.append(p)
        pholm = su.holm(praw) if (_HAVE_SU and praw) else praw
        for (tc, nb, ne, mb, me, U, delta, p), ph in zip(pairs, pholm):
            star = su.p_stars(ph) if _HAVE_SU else ""
            pr = su.fmt_p(p) if _HAVE_SU else f"p={p}"
            phs = su.fmt_p(ph) if _HAVE_SU else f"p={ph}"
            A(f"        {TOOL_SHORT.get(tc, tc):9s} median Bench {100*mb:.1f}% (n={nb}) vs "
              f"Exp {100*me:.1f}% (n={ne})  U={U:.1f} {pr} (Holm {phs}) {star} Cliff's δ={delta:+.2f}")
    A("")
    A("2c. TM-loss RATE of PB-valid poses (inside-TM / PB-valid) — the 'lost to TM %'")
    A("    headline of fig_tool_tm_loss, compared across datasets (feeds fig_tm_loss_relative).")
    A("    Δpp = Benchmark loss% − Experimental loss% (percentage POINTS). The ratio column")
    A("    is Benchmark/Experimental (×), reported for reference ONLY — with the tiny")
    A("    Experimental base it is unstable, so the figure shows Δpp, not the ratio.")
    praw, pairs = [], []
    for tc in SELECTED_TOOLCHAINS:
        rb = pooled[(pooled.dataset == "orai_benchmark") & (pooled.toolchain == tc)]
        re_ = pooled[(pooled.dataset == "orai_jku") & (pooled.toolchain == tc)]
        if not len(rb) or not len(re_):
            continue
        rb, re_ = rb.iloc[0], re_.iloc[0]
        kb, nb = int(rb["inside_tm"]), int(rb["pb_valid"])
        ke, ne = int(re_["inside_tm"]), int(re_["pb_valid"])
        if nb == 0 or ne == 0:
            continue
        p = su.two_proportion_test(kb, nb, ke, ne).get("p") if _HAVE_SU else float("nan")
        pairs.append((tc, kb, nb, ke, ne)); praw.append(p)
    pholm = su.holm(praw) if (_HAVE_SU and praw) else praw
    for (tc, kb, nb, ke, ne), ph, p in zip(pairs, pholm, praw):
        lb, le = kb / nb, ke / ne
        ratio = (lb / le) if le > 0 else float("inf")
        star = su.p_stars(ph) if _HAVE_SU else ""
        pr = su.fmt_p(p) if _HAVE_SU else f"p={p}"
        phs = su.fmt_p(ph) if _HAVE_SU else f"p={ph}"
        A(f"        {TOOL_SHORT.get(tc, tc):9s} Bench {kb}/{nb} ({100*lb:.1f}%) vs "
          f"Exp {ke}/{ne} ({100*le:.1f}%)  Δ={100*(lb-le):+.1f}pp  ratio={ratio:.2f}×  "
          f"Fisher {pr} (Holm {phs}) {star}")
    A("")
    A("=" * 88)
    A("3. Per-(frame × ligand) distribution — outside-TM share (feeds the whisker figure)")
    A("=" * 88)
    A(f"  {'dataset':26s} {'tool':10s} {'n_units':>7s} {'median':>7s} {'IQR':>17s} {'mean':>6s}")
    for ds in CATEGORY_ORDER:
        for tc in SELECTED_TOOLCHAINS:
            v = 100 * per_unit[(per_unit.dataset == ds) &
                               (per_unit.toolchain == tc)]["outside_tm_share"].values
            if len(v) == 0:
                continue
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            A(f"  {CATEGORY_LABELS.get(ds, ds):26s} {TOOL_SHORT.get(tc, tc):10s} {len(v):7d}"
              f" {med:6.1f}% [{q1:5.1f}-{q3:5.1f}]% {np.mean(v):5.1f}%")
    A("")
    A("=" * 88)
    A("CAVEATS")
    A("=" * 88)
    A("• Pseudoreplication: each tool emits many correlated poses per (frame × ligand).")
    A("  Section 2a pooled tests are on pose counts and overstate significance; Section 2b")
    A("  on per-unit shares is the honest unit-of-analysis test.")
    A("• Experimental (Orai × JKU) = 3 ligands × 4 frames (≤12 units/tool). Its quartiles")
    A("  are unstable — read the jittered points, not the box. Formal inference is exploratory.")
    A("• Frames are 4 snapshots of the SAME Orai receptor, so units within a ligand are also")
    A("  correlated; the per-unit test treats them as independent and is therefore anti-")
    A("  conservative for the ligand-generalisation question.")
    A("• EquiBind's poses are overwhelmingly placed in the transmembrane, so its PB-valid")
    A("  share is non-trivial but its outside-TM (usable) share collapses to ~0 in both sets;")
    A("  the whisker box for EquiBind is a near-degenerate distribution (a bar/point would")
    A("  carry the same information) — kept as a box only for visual parity across tools.")
    out.write_text("\n".join(L) + "\n")
    print(f"  wrote {out}")


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root", default="posebusters_results", type=Path,
                    help="root holding <dataset>/transmembrane_filter/tm_pose_classification.csv")
    ap.add_argument("--out-dir", default=None, type=Path,
                    help="output directory (default: <results-root>/orai_pbvalid_tm_share_compare)")
    ap.add_argument("--diffdock-variant", default="smina",
                    help="DiffDock refine variant to select (must match the parent "
                         "exclusion run; default smina → 'DiffDock (smina)').")
    ap.add_argument("--equibind-variant", default="gnina",
                    help="EquiBind refine variant to select (default gnina → 'EquiBind*').")
    ap.add_argument("--top-n-poses", type=int, default=0,
                    help="Keep only each tool's top-N ranked poses per (frame × ligand) "
                         "before computing shares (0 = no cap). Use 10 to equalise the "
                         "pose budget across tools (Orai × Experimental AutoDock/EquiBind "
                         "produce 30/unit; DiffDock and all Benchmark are already ≤10).")
    args = ap.parse_args(argv)

    global SELECTED_TOOLCHAINS, TOOL_SHORT
    SELECTED_TOOLCHAINS, TOOL_SHORT = _resolve_toolchains(args.diffdock_variant,
                                                          args.equibind_variant)

    root = args.results_root
    out_dir = args.out_dir or (root / "orai_pbvalid_tm_share_compare")
    out_dir.mkdir(parents=True, exist_ok=True)

    frames = {}
    src_notes = []
    for ds in CATEGORY_ORDER:
        p = _classification_path(ds, root)
        frames[ds] = load_dataset(ds, root, top_n=args.top_n_poses)
        mtime = datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        src_notes.append(f"    {ds}: {p}  (source mtime {mtime})")
    if args.top_n_poses:
        data_source_note_cap = (f"    POSE CAP: kept only each tool's top-{args.top_n_poses} "
                                f"ranked poses per (frame × ligand).")
    else:
        data_source_note_cap = ""
    data_source_note = "Source per-pose tables (TM-classified, no geometry recomputed):\n" + \
        "\n".join(src_notes + ([data_source_note_cap] if data_source_note_cap else []))

    pooled = pooled_table(frames)
    per_unit = per_unit_table(frames)
    pooled.to_csv(out_dir / "pbvalid_tm_share_pooled.csv", index=False)
    per_unit.to_csv(out_dir / "pbvalid_tm_share_per_unit.csv", index=False)
    print(f"  wrote {out_dir/'pbvalid_tm_share_pooled.csv'}")
    print(f"  wrote {out_dir/'pbvalid_tm_share_per_unit.csv'}")

    fig_share_compare(pooled, out_dir / "fig_pbvalid_share_compare.png")
    fig_fate_stack(pooled, out_dir / "fig_tm_fate_share_compare.png")
    fig_tm_loss_relative(pooled, out_dir / "fig_tm_loss_relative_compare.png")
    fig_whisker(per_unit, out_dir / "fig_pbvalid_outside_tm_whisker.png")
    write_stats(pooled, per_unit, out_dir / "pbvalid_tm_share_compare_stats.txt",
                data_source_note)
    print("done.")


if __name__ == "__main__":
    main()
