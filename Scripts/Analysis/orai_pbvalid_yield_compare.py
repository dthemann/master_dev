#!/usr/bin/env python
"""Orai1 PoseBusters-valid yield, control panel against experimental panel.

Produces thesis Figure 8 (`09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png`)
plus its per-complex table and its off-panel statistics.

PROVENANCE. Until 2026-09-02 this figure was the only float in the Results
chapter with no committed generator: it lived as an inline cell (id ``ee79ca6e``)
in ``Master_Docking_AD_Full_Protein.ipynb``. That made it the single item
blocking a script-driven rebuild of the chapter, so the cell was lifted here
verbatim and given a command line. The plotting body, the layout constants, the
colour grammar and the statistics are unchanged; only the hard-coded paths and
the three variant pins became arguments.

WHAT IT MEASURES. Yield is per complex, meaning one Orai frame crossed with one
ligand: the share of the poses a variant PRODUCED that pass all PoseBusters
checks. Each box is that per-complex yield distributed over complexes. The full
produced-pose CSV is used, transmembrane and off-pocket poses included, so the
denominator is everything the tool generated. It is a geometry-quality measure
and is deliberately independent of where the pose landed. The operational
placement endpoint is a different figure (`orai_pbvalid_tm_share_compare.py`).

DEFAULT TREE. ``posebusters_results/_orai_matched_root`` is the canonical Orai
source after the 2026-08-31 matched-EquiBind re-run. The plain
``posebusters_results/orai_{benchmark,jku}`` trees are one generation behind and
were what the notebook cell read. The shipped Figure 8 asset matches the matched
tree, so that is the default here.

Statistics: the two panels are DIFFERENT ligand sets, so the samples are
independent and the contrast is an unpaired Mann-Whitney U with Cliff's delta.
It is deliberately not the paired Wilcoxon that the benchmark raw-to-optimised
figure uses, because that one compares the same poses under two treatments.
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.patheffects as _pe  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

warnings.filterwarnings("ignore")

_ANALYSIS = str(Path(__file__).resolve().parent)
if _ANALYSIS not in sys.path:
    sys.path.insert(0, _ANALYSIS)

import posebusters_pose_comparison as ppc  # noqa: E402
import stats_utils as su  # noqa: E402
from pose_topn import top_n_allowlist  # noqa: E402

DEFAULT_ROOT = Path("posebusters_results/_orai_matched_root")

TOOL_OF = {"autodock": "autodock", "autodock_gnina": "autodock",
           "diffdock_smina": "diffdock", "diffdock_gnina": "diffdock",
           "equibind_unguided_gnina": "equibind", "equibind_unguided_smina": "equibind"}
TOOL_ORDER = ["autodock", "diffdock", "equibind"]
TOOL_NAME = {"autodock": "AutoDock Vina", "diffdock": "DiffDock", "equibind": "EquiBind"}
_VLABEL = {"autodock": "—", "autodock_gnina": "gnina",
           "diffdock_smina": "smina", "diffdock_gnina": "gnina",
           "equibind_unguided_gnina": "gnina", "equibind_unguided_smina": "smina"}

# Box fill colour encodes the TOOL, reusing posebusters_pose_comparison.TOOL_COLORS so
# this figure shares the palette of the benchmark yield boxplot exactly. The two ligand
# SETS are told apart by fill style alone, solid against hatched.
_TOOL_COLOUR_KEY = {"autodock": "autodock", "diffdock": "diffdock", "equibind": "equibind_unguided"}
_TOOL_COLOUR_FALLBACK = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind_unguided": "#2ca02c"}

EQ_PREF = ["equibind_unguided_gnina", "equibind_unguided_smina"]


def _tool_colour(tool: str) -> str:
    key = _TOOL_COLOUR_KEY.get(tool, "")
    return ppc.TOOL_COLORS.get(key) or _TOOL_COLOUR_FALLBACK.get(key, "#888888")


def _load_full(csv: Path, top_n: int, exclude_ligands: set[str]) -> pd.DataFrame:
    """Per-pose frame for one panel, capped at the reported ranking depth.

    ``top_n_allowlist`` applies each tool's OWN order: Vina rank for AutoDock,
    confidence for DiffDock, gnina energy for EquiBind. Appendix B records that
    the cap is applied on the gnina CNNaffinity order for AutoDock, which is why
    66 of the 120 retained experimental poses are Vina modes 11 to 30.
    """
    df = ppc._build_pose_index(Path(csv), split_equibind=True, select_diffdock=False)
    df = df[df["pose_file"].astype(str).isin(top_n_allowlist(str(csv), top_n))].copy()
    df = df.rename(columns={"docking_method": "method"})
    if exclude_ligands:
        df = df[~df["ligand"].astype(str).isin(exclude_ligands)].copy()
    return df


def _pick_equibind(dff: pd.DataFrame, override: str | None):
    """Resolve the EquiBind variant, preferring gnina unless its coverage lags.

    The coverage fallback is a legacy guard from when the control panel's gnina
    re-run was incomplete at one of four frames. On the matched tree both panels
    carry complete gnina coverage, so it never fires. It is kept so the script
    degrades visibly rather than silently on a partial tree.
    """
    cov = {k: dff[dff["method"] == k].groupby(["protein", "ligand"]).ngroups for k in EQ_PREF}
    if override and override != "auto":
        return override, cov
    gn, sm = cov.get("equibind_unguided_gnina", 0), cov.get("equibind_unguided_smina", 0)
    if gn >= max(1, sm) * 0.8:
        return ("equibind_unguided_gnina" if gn > 0 else None), cov
    return ("equibind_unguided_smina" if sm > 0 else None), cov


def build(args) -> int:
    csv_name = ("posebusters_filtered_results.no_tm.csv" if args.use_no_tm
                else "posebusters_filtered_results.csv")
    root = Path(args.results_root)
    datasets = {
        "orai_benchmark": {
            "label": "Orai × Benchmark ligands", "short": "Benchmark",
            "csv": Path(args.bench_csv) if args.bench_csv
            else root / "orai_benchmark" / "dock" / csv_name,
            "alpha": 0.80, "hatch": "",
            "equibind": args.equibind_variant_benchmark,
        },
        "orai_jku": {
            "label": "Orai × Experimental ligands (JKU)", "short": "Experimental",
            "csv": Path(args.exp_csv) if args.exp_csv
            else root / "orai_jku" / "dock" / csv_name,
            "alpha": 0.45, "hatch": "///",
            "equibind": args.equibind_variant_experimental,
        },
    }
    out_dir = Path(args.out_dir) if args.out_dir else root / "orai_pbvalid_yield_compare"
    out_dir.mkdir(parents=True, exist_ok=True)
    exclude = set(args.exclude_ligand or [])

    per_frames, sel = [], {}
    for ds, cfg in datasets.items():
        if not Path(cfg["csv"]).exists():
            sys.exit(f"{cfg['csv']} not found - run that panel's PoseBusters stage first.")
        dff = _load_full(cfg["csv"], args.top_n_poses, exclude)
        eq_key, eq_cov = _pick_equibind(dff, cfg["equibind"])
        keys = {"autodock": args.autodock_variant, "diffdock": args.diffdock_variant,
                "equibind": eq_key}
        sel[ds] = {"keys": keys, "eq_cov": eq_cov}
        specs = [(TOOL_OF[k], "dominant", k != "autodock", k) for k in keys.values() if k]
        per = ppc.pbvalid_yield_per_complex(dff, specs)
        per["dataset"] = ds
        per_frames.append(per)
        print(f"{ds}: EquiBind coverage {eq_cov} -> dominant EquiBind = {eq_key}")
    per_all = pd.concat(per_frames, ignore_index=True)
    per_all.to_csv(out_dir / "pbvalid_yield_dominant_per_complex.csv", index=False)

    # ---- layout: two boxes (Benchmark, Experimental) per tool group ----------
    intra, inter, box_w = 1.06, 1.05, 0.66
    ds_list = list(datasets)
    boxes, positions, spans, x = [], [], {}, 0.0
    for tool in TOOL_ORDER:
        xs = []
        for ds in ds_list:
            sub = per_all[(per_all["tool"] == tool) & (per_all["dataset"] == ds)]
            data = sub["pct_valid"].to_numpy(float)
            data = data[~np.isnan(data)]
            if data.size == 0:
                x += intra
                continue
            boxes.append(dict(tool=tool, dataset=ds, mkey=sel[ds]["keys"][tool], data=data, x=x,
                              n=int(data.size), median=float(np.median(data)),
                              color=_tool_colour(tool), alpha=datasets[ds]["alpha"],
                              hatch=datasets[ds]["hatch"]))
            positions.append(x)
            xs.append(x)
            x += intra
        if xs:
            spans[tool] = xs
        x += inter

    fig, ax = plt.subplots(figsize=(9.6, 6.4))
    bp = ax.boxplot([b["data"] for b in boxes], positions=positions, widths=box_w,
                    patch_artist=True, showmeans=True, showfliers=False,
                    medianprops=dict(color="#222222", lw=1.6),
                    meanprops=dict(marker="D", markerfacecolor="white",
                                   markeredgecolor="#222222", markersize=5))
    for patch, b in zip(bp["boxes"], boxes):
        patch.set_facecolor(b["color"])
        patch.set_alpha(b["alpha"])
        patch.set_edgecolor(b["color"])
        patch.set_linewidth(1.4)
        if b["hatch"]:
            patch.set_hatch(b["hatch"])

    for b in boxes:
        ax.text(b["x"] - box_w / 2, b["median"], f"{b['median']:.1f}%", ha="left", va="center",
                fontsize=8.6, fontweight="bold", color="#111111", zorder=6,
                path_effects=[_pe.withStroke(linewidth=3.2, foreground="white")])

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{datasets[b['dataset']]['short']}\n{_VLABEL[b['mkey']]} · n={b['n']}"
                        for b in boxes], fontsize=9.0)
    ax.tick_params(axis="x", length=0, pad=6)
    ax.tick_params(axis="y", labelsize=11.5)
    trans = ax.get_xaxis_transform()
    for tool, xs in spans.items():
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - box_w / 2, hi + box_w / 2], [-0.185, -0.185], transform=trans,
                color="0.45", lw=1.0, clip_on=False, zorder=1)
        ax.text(xc, -0.205, TOOL_NAME[tool], transform=trans, ha="center", va="top",
                fontsize=13.0, fontweight="bold", clip_on=False)

    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_ylabel("PoseBusters-valid poses produced\n"
                  "(% of a complex's produced poses, per complex)", fontsize=13.0)
    ax.set_ylim(-3, 103.0)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_bounds(0, 100)
    ax.spines["bottom"].set_bounds(min(positions) - 0.7, max(positions) + 0.7)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.grid(axis="y", alpha=0.3)

    handles = [Patch(facecolor="0.5", alpha=datasets[d]["alpha"], edgecolor="0.5",
                     hatch=datasets[d]["hatch"] or None, label=datasets[d]["label"])
               for d in ds_list]
    handles.append(plt.Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
                              markeredgecolor="#222222", markersize=8, label="mean"))
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.015),
              ncol=3, frameon=False, fontsize=10.0)
    fig.suptitle("Orai docking — PoseBusters-valid yield of produced poses",
                 fontsize=13, fontweight="bold", y=0.99)
    png = out_dir / "09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png"
    fig.savefig(png, dpi=160, bbox_inches="tight")
    plt.close(fig)

    # ---- off-panel statistics ----------------------------------------------
    lines = ["Orai PB-valid yield (dominant variant) - Benchmark vs Experimental (JKU) ligands",
             "Per tool: unpaired Mann-Whitney U on per-complex yield + Cliff's delta "
             "(d>0 => Benchmark higher)."]
    if _VLABEL.get(sel["orai_benchmark"]["keys"]["equibind"]) != \
            _VLABEL.get(sel["orai_jku"]["keys"]["equibind"]):
        lines.append("CAVEAT: the EquiBind refiner differs by set (see variant col) so that "
                     "row's contrast is confounded.")
    lines.append("")
    rows = []
    for tool in TOOL_ORDER:
        a = per_all[(per_all.tool == tool) & (per_all.dataset == "orai_benchmark")]["pct_valid"].to_numpy(float)
        b = per_all[(per_all.tool == tool) & (per_all.dataset == "orai_jku")]["pct_valid"].to_numpy(float)
        mw = su.mannwhitney_cliffs(a, b)
        if not mw:
            continue
        va = _VLABEL[sel["orai_benchmark"]["keys"][tool]]
        vb = _VLABEL[sel["orai_jku"]["keys"][tool]]
        rows.append(dict(tool=tool, variant_benchmark=va, variant_experimental=vb,
                         n_benchmark=mw["n_a"], n_experimental=mw["n_b"],
                         median_benchmark=round(mw["medians"][0], 1),
                         median_experimental=round(mw["medians"][1], 1),
                         cliffs_delta=round(mw["cliffs_delta"], 3), p=mw["p"],
                         p_stars=su.p_stars(mw["p"])))
        lines.append(f"{TOOL_NAME[tool]:14s} Benchmark({va}) med={mw['medians'][0]:5.1f}%  "
                     f"Experimental({vb}) med={mw['medians'][1]:5.1f}%  "
                     f"d={mw['cliffs_delta']:+.3f}  {su.fmt_p(mw['p'])} {su.p_stars(mw['p'])}")
    pd.DataFrame(rows).to_csv(out_dir / "pbvalid_yield_dominant_compare_stats.csv", index=False)
    (out_dir / "09f_pbvalid_yield_dominant_stats.txt").write_text("\n".join(lines) + "\n")

    print(f"\nWrote -> {png}")
    print(f"         {out_dir / 'pbvalid_yield_dominant_per_complex.csv'}")
    print(f"         {out_dir / 'pbvalid_yield_dominant_compare_stats.csv'} + ..._stats.txt")
    print("\n".join(lines))
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results-root", default=str(DEFAULT_ROOT),
                    help=f"root holding orai_benchmark/ and orai_jku/ (default: {DEFAULT_ROOT})")
    ap.add_argument("--bench-csv", default=None, help="override the control-panel PoseBusters CSV")
    ap.add_argument("--exp-csv", default=None, help="override the experimental-panel CSV")
    ap.add_argument("--out-dir", default=None,
                    help="default: <results-root>/orai_pbvalid_yield_compare")
    ap.add_argument("--use-no-tm", action="store_true",
                    help="read the transmembrane-excluded CSV instead of every produced pose. "
                         "OFF for the thesis figure, which measures geometry quality and must "
                         "not be conditioned on placement.")
    ap.add_argument("--top-n-poses", type=int, default=10,
                    help="ranking-depth cap per tool, in each tool's own order (default 10)")
    ap.add_argument("--autodock-variant", default="autodock_gnina",
                    help="AutoDock method key. The bare 'autodock' key is the RAW search and "
                         "silently drops the gnina rescoring the Orai chapter reports.")
    ap.add_argument("--diffdock-variant", default="diffdock_smina")
    ap.add_argument("--equibind-variant-benchmark", default="auto",
                    help="'auto' picks by coverage, preferring gnina")
    ap.add_argument("--equibind-variant-experimental", default="auto")
    ap.add_argument("--exclude-ligand", action="append",
                    default=None, metavar="NAME",
                    help="repeatable. Defaults to gsk7975a-deprot-OPT, the deprotonated "
                         "GSK-7975A duplicate, so the experimental panel is the three "
                         "intended compounds over four frames.")
    args = ap.parse_args(argv)
    if args.exclude_ligand is None:
        args.exclude_ligand = ["gsk7975a-deprot-OPT"]
    return build(args)


if __name__ == "__main__":
    sys.exit(main())
