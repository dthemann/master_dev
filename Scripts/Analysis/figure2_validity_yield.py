#!/usr/bin/env python
"""Figure 2 — Post-hoc optimisation and PoseBusters validity.

Regenerates the thesis Figure 2 panel, a grouped box plot of per-complex
PoseBusters-valid yield for each engine family shown raw beside its optimised
variants.

Why this script exists
----------------------
`posebusters_validity_report.py` draws every variant present in the CSV. Figure 2
shows a curated nine, namely AutoDock Vina, DiffDock and unguided EquiBind, each
raw plus two post-hoc variants. Since that figure was first drawn the filtered
results CSV also gained AutoDock Vinardo and Uni-Dock2, so no combination of the
report script's existing flags reproduces the published panel set. This script
supplies the missing engine keep-set as an explicit `--series` list, so the panel
is a stated choice rather than whatever happens to be in the CSV.

Validity is not recomputed here. The CSV is scored through
`posebusters_validity_report.load_and_score`, which applies the canonical
PoseBusters check schema and performs the same EquiBind/AutoDock/DiffDock variant
splits as the report, so a series name here means exactly what it means there.

Statistics follow the thesis Methods. Each optimised variant is compared with its
family's raw arm by a paired Wilcoxon signed-rank test over the complexes both
arms populate, the effect is a Hodges-Lehmann estimate of the median paired
difference in percentage points, and p-values are Holm-corrected across the
comparisons drawn on the panel.

Usage
-----
    python Scripts/Analysis/figure2_validity_yield.py \
        --csv posebusters_results/benchmark_full_protein_vina_scoring/dock/posebusters_filtered_results.csv \
        --out thesis_latex/media/media/image2.png
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
import posebusters_validity_report as V  # noqa: E402

# Family colours follow the thesis house style, one hue per engine family. Raw
# arms are the same hue at reduced alpha and hatched, so "before optimisation"
# reads from the fill pattern and not from a second hue.
FAMILIES = [
    ("AutoDock Vina", "#5D9DC9", [
        ("autodock", "raw"),
        ("autodock_gnina", "gnina-opt"),
        ("autodock_gnina_refinement", "gnina CNN-refine"),
    ]),
    ("DiffDock", "#FFA251", [
        ("diffdock", "raw"),
        ("diffdock_smina", "smina-opt"),
        ("diffdock_gnina", "gnina-opt"),
    ]),
    ("EquiBind (blind / unguided)", "#67BA67", [
        ("equibind_unguided_raw", "raw"),
        ("equibind_unguided_smina", "smina-opt"),
        ("equibind_unguided_gnina", "gnina-opt"),
    ]),
]
RAW_ALPHA = 0.45


def hodges_lehmann(diffs: np.ndarray) -> float:
    """Median of all Walsh averages, the paired-difference estimator the thesis
    Methods specify. Uses the i <= j triangle so each pair counts once."""
    d = np.asarray(diffs, dtype=float)
    i, j = np.triu_indices(len(d), 0)
    return float(np.median((d[i] + d[j]) / 2.0))


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni step-down, NaN-safe and monotone."""
    m = len(pvals)
    order = sorted(range(m), key=lambda k: pvals[k])
    out = [float("nan")] * m
    running = 0.0
    for rank, k in enumerate(order):
        adj = min(1.0, (m - rank) * pvals[k])
        running = max(running, adj)
        out[k] = running
    return out


def stars(p: float) -> str:
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def per_complex_yield(df: pd.DataFrame, method: str) -> pd.Series:
    """Percent of a complex's produced poses that are PoseBusters-valid."""
    sub = df[df["docking_method"] == method]
    if sub.empty:
        raise SystemExit(f"series '{method}' is absent from the CSV")
    return sub.groupby("protein")["pb_valid"].mean() * 100.0


def build(df: pd.DataFrame, families: list) -> tuple[list[dict], list[dict]]:
    series, brackets = [], []
    pos = 0.0
    for fam, colour, members in families:
        raw_key = members[0][0]
        raw_y = per_complex_yield(df, raw_key)
        fam_start = pos
        for idx, (key, label) in enumerate(members):
            y = per_complex_yield(df, key)
            sub = df[df["docking_method"] == key]
            series.append(dict(
                key=key, label=label, family=fam, colour=colour, pos=pos,
                raw=(idx == 0), y=y.values,
                median=float(y.median()), mean=float(y.mean()),
                n_complexes=int(sub.loc[sub["pb_valid"], "protein"].nunique()),
                n_valid=int(sub["pb_valid"].sum()),
            ))
            if idx > 0:
                a, b = raw_y.align(y, join="inner")
                d = (b - a).dropna()
                brackets.append(dict(
                    family=fam, lo=fam_start, hi=pos, tier=idx,
                    delta=hodges_lehmann(d.values),
                    p=float(wilcoxon(b.loc[d.index], a.loc[d.index]).pvalue),
                    n=int(len(d)),
                ))
            pos += 1.0
        pos += 0.9  # gap between families
    for br, p_adj in zip(brackets, holm([b["p"] for b in brackets])):
        br["p_holm"] = p_adj
    return series, brackets


def draw(series: list[dict], brackets: list[dict], out: Path, n_complexes: int,
         title: str) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 6.3))

    for s in series:
        bp = ax.boxplot([s["y"]], positions=[s["pos"]], widths=0.62,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", linewidth=1.6),
                        whiskerprops=dict(color="black"),
                        capprops=dict(color="black"),
                        boxprops=dict(edgecolor="black", linewidth=1.0))
        box = bp["boxes"][0]
        box.set_facecolor(s["colour"])
        if s["raw"]:
            box.set_alpha(RAW_ALPHA)
            box.set_hatch("///")
        ax.plot([s["pos"]], [s["mean"]], marker="D", markersize=7,
                markerfacecolor="white", markeredgecolor="black", zorder=5)
        ax.annotate(f"{s['median']:.1f}%", (s["pos"] - 0.34, s["median"]),
                    ha="right", va="center", fontsize=9, fontweight="bold")

    top = 100.0
    for br in brackets:
        h = top + 12 + 16 * (br["tier"] - 1)
        ax.plot([br["lo"], br["lo"], br["hi"], br["hi"]],
                [h - 3, h, h, h - 3], lw=1.0, color="0.35")
        p = br["p_holm"]
        ptxt = "<1e-4" if p < 1e-4 else f"{p:.3g}"
        ax.text((br["lo"] + br["hi"]) / 2, h + 1.5,
                f"{stars(p)}  {br['delta']:+.1f} pp\n{ptxt}",
                ha="center", va="bottom", fontsize=9.5,
                color="#1a7a1a", fontweight="bold")

    ax.set_xticks([s["pos"] for s in series])
    ax.set_xticklabels([f"{s['label']}\n{s['n_complexes']}/{s['n_valid']}"
                        for s in series], fontsize=9.5)
    ax.set_ylabel("PoseBusters-valid poses produced\n"
                  "(% of a complex's poses, per complex)", fontsize=10)
    ax.set_ylim(-4, top + 12 + 16 * 2 + 14)
    ax.set_yticks(range(0, 101, 20))
    ax.grid(axis="y", alpha=0.3)
    ax.set_axisbelow(True)
    ax.set_title(f"{title}   (n = {n_complexes} complexes)",
                 fontsize=14, fontweight="bold", pad=58)

    # Family underlines beneath the tick labels, in axis-fraction space.
    for fam, _c, _m in FAMILIES:
        xs = [s["pos"] for s in series if s["family"] == fam]
        lo, hi = min(xs) - 0.45, max(xs) + 0.45
        ax.plot([lo, hi], [-0.155, -0.155], transform=ax.get_xaxis_transform(),
                clip_on=False, color="0.4", lw=1.2)
        ax.text((lo + hi) / 2, -0.20, fam, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=11.5, fontweight="bold")

    handles = [
        Patch(facecolor="0.75", edgecolor="black", hatch="///",
              label="raw / reference (before optimization)"),
        Patch(facecolor="0.75", edgecolor="black",
              label="post-hoc optimized (smina / gnina)"),
        Line2D([], [], marker="D", linestyle="none", markersize=8,
               markerfacecolor="white", markeredgecolor="black", label="mean"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.075),
              ncol=3, frameon=True, fontsize=10)

    fig.subplots_adjust(bottom=0.20, top=0.84)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def write_sidecar(series: list[dict], brackets: list[dict], out: Path,
                  csv_path: Path, n_complexes: int) -> None:
    L = [f"Figure 2 — Post-hoc optimisation and PoseBusters validity",
         "=" * 74,
         f"source CSV : {csv_path}",
         f"complexes  : {n_complexes}",
         "yield      : percent of a complex's produced poses that are PB-valid",
         "label      : complexes with at least one valid pose / total valid poses",
         "",
         f"{'series':32s} {'median':>8} {'mean':>8} {'cplx':>6} {'valid':>8}"]
    for s in series:
        L.append(f"{s['key']:32s} {s['median']:8.1f} {s['mean']:8.1f} "
                 f"{s['n_complexes']:6d} {s['n_valid']:8d}")
    L += ["", "Paired contrasts against each family's raw arm",
          "(Wilcoxon signed-rank, Hodges-Lehmann median paired difference, Holm-corrected)"]
    for br in brackets:
        L.append(f"  {br['family']:28s} tier {br['tier']}  "
                 f"{br['delta']:+6.1f} pp  p_raw={br['p']:.3g}  "
                 f"p_holm={br['p_holm']:.3g}  n={br['n']}")
    out.with_suffix(".txt").write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, required=True,
                    help="posebusters_filtered_results.csv for the run")
    ap.add_argument("--out", type=Path, required=True, help="output PNG path")
    ap.add_argument("--series", default=None,
                    help="engine keep-set override, semicolon-separated families as "
                         "'Label:colour:key=tick,key=tick,...'. Omit for the published nine.")
    ap.add_argument("--title", default="Post-Hoc Optimization and PoseBuster Validity")
    args = ap.parse_args()

    families = FAMILIES
    if args.series:
        families = []
        for fam_spec in args.series.split(";"):
            label, colour, members = fam_spec.split(":", 2)
            families.append((label, colour,
                             [tuple(m.split("=", 1)) for m in members.split(",")]))

    df = V.load_and_score(args.csv)
    n_complexes = int(df["protein"].nunique())
    series, brackets = build(df, families)
    draw(series, brackets, args.out, n_complexes, args.title)
    write_sidecar(series, brackets, args.out, args.csv, n_complexes)
    print(f"Wrote {args.out} and {args.out.with_suffix('.txt')}")
    for s in series:
        print(f"  {s['key']:32s} median={s['median']:5.1f}%  "
              f"{s['n_complexes']}/{s['n_valid']}")


if __name__ == "__main__":
    main()
