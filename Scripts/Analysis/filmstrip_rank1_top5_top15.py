#!/usr/bin/env python
"""Depth filmstrip (fig 20d, _pbvalid variant) with buckets rank-1 / top-5 / top-15.

Companion to ``20d_form_vs_placement_by_family__depth_filmstrip_pbvalid.png`` (which
uses top-1/3/5). This regenerates the SAME mechanism-coloured form-vs-placement
filmstrip — PB-valid poses, RMSD ≤ 2 Å gate removed, poses outside the crystal-site
neighbourhood trimmed (docked centroid > 8 Å from the crystal site; these are NOT off
the receptor — nearly all remain in van der Waals contact with the protein), axes
clipped at 5 Å by ``AXIS_MAX``, where the sibling top-1/3/5 figure clips at 8 Å —
but with the ranking-depth columns set to rank-1 / top-5 / top-15, and dumps the
exact per-pose data behind every plotted point to CSV.

Variant selection reproduces the notebook (cell 21) EXACTLY: ``--collapse-plots-only``
implicitly turns on ``--best-equibind-only`` + ``--best-diffdock-only``, and
``--collapse-diffdock-variant diffdock_smina`` pins DiffDock. So the report collapses
each tool to one variant BEFORE plotting:
    AutoDock  → autodock_gnina           (relabelled 'autodock'; the arm the Results
                                          chapter reports, pinned by FORCED_AUTODOCK)
    DiffDock  → diffdock_smina           (relabelled 'diffdock' → shown 'DiffDock*')
    EquiBind  → equibind_unguided_gnina  (oracle-best EquiBind by PB-valid & RMSD≤2Å)
This must go through :func:`_select_best_diffdock` / :func:`_select_best_equibind`, not
a plain method filter: the collapse re-derives DiffDock's rank from the pose filename
(the cached smina rank is a flat 999), which the plain filter would get wrong.

Reuses the pipeline's own helpers so the methodology is byte-identical to the sibling
figure; only ``depths`` changes. Run in the ``vina`` conda env.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Whole-protein (blind) run. This is the campaign the Results chapter reports, and it
# is the only cache carrying optimizer / autodock_rank / optimized_rank, which
# ``P._read_cached_per_pose`` requires. The older boxed cache
# (posebusters_results/benchmark/dock/pose_comparison_report) lacks all three and also
# predates the AutoDock gnina arm entirely, so it cannot back this figure.
DEFAULT_REPORT = Path(
    "posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report")
DEPTHS = (1, 5, 15)
AXIS_MAX = 5.0                       # x/y axes clipped to this square (Å) for legibility
FORCED_DIFFDOCK = "diffdock_smina"   # matches --collapse-diffdock-variant in cell 21
# The whole-protein cache carries six autodock* variants and ``P._fam_key`` maps every
# one of them to the "autodock" family, so they would be pooled into a single scatter.
# The chapter reports AutoDock Vina + gnina, so that arm is pinned here and relabelled
# to "autodock", exactly as _select_best_diffdock relabels its chosen variant.
FORCED_AUTODOCK = "autodock_gnina"
FAM_LABEL = {"autodock": "AutoDock", "diffdock": "DiffDock", "equibind": "EquiBind"}

REPORT = DEFAULT_REPORT
CACHE = REPORT / "per_pose_metrics.csv"

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import posebusters_pose_comparison as P  # noqa: E402


def collapse_like_cell21(df: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, str]]:
    """Apply the exact best-variant collapse cell 21 does, and return the collapsed
    frame plus a family→source-variant map so the CSV can name the true variant behind
    the relabelled 'diffdock' family."""
    df = _drop_non_thesis_tools(df)
    df, best_ad = _select_autodock_arm(df, forced=FORCED_AUTODOCK)
    df, best_eq = P._select_best_equibind(df)
    df, best_dd = P._select_best_diffdock(df, forced=FORCED_DIFFDOCK)
    src = {"autodock": best_ad or "autodock",
           "diffdock": best_dd or FORCED_DIFFDOCK,
           "equibind": best_eq or "equibind"}
    print(f"collapse → AutoDock={src['autodock']} · DiffDock={src['diffdock']} "
          f"· EquiBind={src['equibind']}")
    return df, src


THESIS_TOOL_PREFIXES = ("autodock", "diffdock", "equibind")


def _drop_non_thesis_tools(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the three tool families this thesis reports.

    ``P._fam_key`` classifies by prefix and FALLS THROUGH to "equibind" for anything
    that is neither autodock* nor diffdock*, so any other engine in the cache is
    silently pooled into the EquiBind scatter. The whole-protein cache carries
    Uni-Dock2, which the thesis excludes: left in, it would inflate the EquiBind
    rank-1 cohort from 126 to 288 poses and make one complex contribute two "rank-1"
    points. That never reached the deployed figure, which shows 126, but the defect is
    real and the older boxed cache simply had no Uni-Dock rows to expose it.

    This is a WHITELIST on purpose. A blacklist of known-unwanted engines would let the
    next engine added to the cache re-enter EquiBind just as silently, which is exactly
    how Uni-Dock got in."""
    fam = df["method"].astype(str)
    keep = fam.str.startswith(THESIS_TOOL_PREFIXES)
    if not keep.all():
        dropped = sorted(fam[~keep].unique())
        print(f"dropped {int((~keep).sum()):,} non-thesis poses ({', '.join(dropped)})")
    return df[keep].copy()


def _select_autodock_arm(df: pd.DataFrame,
                         forced: str) -> tuple[pd.DataFrame, str | None]:
    """Keep exactly one autodock* variant and relabel it to 'autodock'.

    ``P._fam_key`` collapses every autodock* method into one family, so on a cache
    that carries several of them the scatter would silently pool all their poses.
    Mirrors ``P._select_best_diffdock``: drop the sibling variants, then rename the
    survivor so downstream family mapping and panel labelling are unchanged."""
    variants = sorted(m for m in df["method"].unique() if str(m).startswith("autodock"))
    if forced not in variants:
        # Fail loudly in every case. Returning the sole available variant would silently
        # rebuild the figure from the wrong arm — e.g. raw Vina on a cache that predates
        # the gnina arm — which is the exact defect this function exists to prevent.
        raise SystemExit(
            f"requested AutoDock arm {forced!r} is not in {CACHE}; available: {variants}")
    keep = ~df["method"].astype(str).str.startswith("autodock") | df["method"].eq(forced)
    out = df[keep].copy()
    out.loc[out["method"].eq(forced), "method"] = "autodock"
    return out, forced


def build_cohorts(df: pd.DataFrame) -> dict[int, pd.DataFrame]:
    """The exact per-depth cohorts the _pbvalid filmstrip plots: PB-valid poses
    within each method's top-d ranked (RMSD gate off), centroid-trimmed to the
    crystal-site neighbourhood,
    with form/inplace/placement/eff_rank added. Nested: top-1 ⊆ top-5 ⊆ top-15."""
    out = {}
    for d in DEPTHS:
        c = P._form_components(
            P._valid_topd_poses(df, d, None)).dropna(subset=["form", "inplace"])
        c = c[c["centroid_dist"].le(P.FAR_FROM_RECEPTOR_CENTROID_A)].copy()
        out[d] = c
    return out


def per_pose_table(cohorts: dict[int, pd.DataFrame], src: dict[str, str],
                   axis_max: float) -> pd.DataFrame:
    """One row per plotted pose (the top-15 superset), with the coordinates behind
    each scatter point and boolean membership flags for the three cumulative panels.
    top-1 ⊆ top-5 ⊆ top-15, so a rank-1 pose is drawn in all three columns — the flags
    reconstruct each panel exactly (panel = rows where in_top{d} is True)."""
    deep = cohorts[DEPTHS[-1]].copy()
    deep["family"] = deep["method"].map(P._fam_key)
    deep["variant"] = deep["family"].map(src)
    deep["mechanism"] = P._mechanism_region(deep).astype(str)
    deep["r_form_share"] = (deep["form"] ** 2 / deep["inplace"] ** 2).clip(0, 1)
    for d in DEPTHS:
        deep[f"in_top{d}"] = deep["eff_rank"].le(d)
    frame_col = f"in_frame_{axis_max:g}A"
    deep[frame_col] = deep["inplace"].le(axis_max) & deep["form"].le(axis_max)
    cols = {
        "family": "family", "variant": "variant", "protein": "protein",
        "ligand": "ligand", "eff_rank": "eff_rank",
        "inplace": "inplace_rmsd_A", "form": "form_bestfit_kabsch_rmsd_A",
        "placement": "placement_rmsd_A", "r_form_share": "r_form_share",
        "mechanism": "mechanism", "centroid_dist": "centroid_dist_A",
        "pb_valid": "pb_valid",
    }
    flags = [f"in_top{d}" for d in DEPTHS] + [frame_col]
    tbl = deep[list(cols) + flags].rename(columns=cols)
    tbl = tbl.sort_values(["family", "protein", "ligand", "eff_rank"]).reset_index(drop=True)
    for c in ("inplace_rmsd_A", "form_bestfit_kabsch_rmsd_A", "placement_rmsd_A",
              "r_form_share", "centroid_dist_A"):
        tbl[c] = tbl[c].round(4)
    return tbl


def cohort_summary(cohorts: dict[int, pd.DataFrame], src: dict[str, str]) -> pd.DataFrame:
    """Compact per-(family, depth) table mirroring the panel titles: pose count,
    distinct receptors/complexes, cohort centroid (median in-place, median form),
    drift from the top-1 centroid, and the placement/mixed/form-limited mechanism
    share — the numbers the panel headings report."""
    rows = []
    fams = ["autodock", "diffdock", "equibind"]
    cent1 = {}
    for d in DEPTHS:
        c = cohorts[d].copy()
        c["family"] = c["method"].map(P._fam_key)
        for fam in fams:
            s = c[c["family"] == fam]
            if s.empty:
                continue
            med_ip = float(s["inplace"].median())
            med_fm = float(s["form"].median())
            if d == DEPTHS[0]:
                cent1[fam] = (med_ip, med_fm)
            c1 = cent1.get(fam)
            mech = P._mechanism_region(s).value_counts(normalize=True).reindex(
                list(P._MECH_LABELS)).fillna(0.0) * 100
            rows.append({
                "family": FAM_LABEL[fam], "variant": src[fam],
                "depth": f"top-{d}", "n_poses": len(s),
                "n_receptors_protein": int(s["protein"].nunique()),
                "n_complexes": int(s.drop_duplicates(["protein", "ligand"]).shape[0]),
                "median_inplace_rmsd_A": round(med_ip, 3),
                "median_form_rmsd_A": round(med_fm, 3),
                "drift_inplace_vs_top1_A": round(med_ip - c1[0], 3) if c1 else 0.0,
                "drift_form_vs_top1_A": round(med_fm - c1[1], 3) if c1 else 0.0,
                "pct_placement_limited": round(float(mech["placement-limited"]), 1),
                "pct_mixed": round(float(mech["mixed"]), 1),
                "pct_form_limited": round(float(mech["form-limited"]), 1),
            })
    return pd.DataFrame(rows)


def _fmt_table(df: pd.DataFrame) -> str:
    """Fixed-width, index-free rendering with compact float/p-value formatting."""
    if df is None or df.empty:
        return "(no rows)"
    return df.to_string(index=False, float_format=lambda x: f"{x:.4g}")


def write_stats_txt(path: Path, summ: pd.DataFrame, stats: dict,
                    src: dict[str, str], depths=DEPTHS) -> None:
    """Write the numbers taken off the figure panels — the per-panel cohort data plus
    every companion statistical test — to one readable text sidecar. Section 1 is the
    exact per-panel data the compact titles no longer show; sections 2–5 are the
    honest-statistics tables (also dumped as the ``__stats_*.csv`` files)."""
    L = []
    L.append("Fig 20d — Form vs placement across ranking depth "
             f"(rank-1 / top-{'/'.join(str(d) for d in depths[1:])}, PB-valid)")
    L.append("Per-panel cohort data + companion statistical tests")
    L.append("Generated by Scripts/Analysis/filmstrip_rank1_top5_top15.py")
    L.append(f"Tool collapse: AutoDock={src['autodock']} · DiffDock={src['diffdock']} "
             f"· EquiBind={src['equibind']}")
    L.append("")

    L.append("=" * 78)
    L.append("1. PER-PANEL COHORT DATA  (the numbers removed from the figure panels)")
    L.append("=" * 78)
    L.append("One row per figure panel (tool family x ranking depth). n_poses = PB-valid")
    L.append("poses inside the crystal-site neighbourhood in the cohort (not off the")
    L.append("receptor); n_receptors_protein / n_complexes = distinct")
    L.append("receptors and protein-ligand complexes those poses span; median_inplace/form =")
    L.append("cohort centroid (the black X); drift_* = shift of that centroid from the top-1")
    L.append("cohort (the orange arrow); pct_placement/mixed/form = the mechanism composition")
    L.append("shares (the 'place/mix/form = ..%' line that used to sit under each panel).")
    L.append("")
    L.append(_fmt_table(summ))
    L.append("")

    if not stats:
        L.append("(SciPy unavailable — companion statistical tests were not computed.)")
        path.write_text("\n".join(L) + "\n", encoding="utf-8")
        return

    L.append("=" * 78)
    L.append("2. PER-TOOL, PER-DEPTH SUMMARY  (medians, IQR, coverage, mechanism %)")
    L.append("=" * 78)
    L.append("Cumulative in-place & form median/IQR, median r = form^2/in-place^2, mechanism")
    L.append("composition (%), and valid-complex COVERAGE of the benchmark universe (each")
    L.append("per-tool cell covers a DIFFERENT complex set, so cross-tool medians here are")
    L.append("not paired — see section 4 for the paired common-complex comparison).")
    L.append("")
    L.append(_fmt_table(stats.get("per_tool_depth")))
    L.append("")

    L.append("=" * 78)
    L.append("3. WITHIN-TOOL DEPTH TREND  (rank bands + Spearman + paired Friedman + drift)")
    L.append("=" * 78)
    L.append("Per (tool, metric): disjoint rank-band medians (r1 / r2-3 / r4-5), a pose-level")
    L.append("Spearman rank-vs-value trend, the within-complex PAIRED Friedman test (the")
    L.append("clustering-robust depth test, run only where >=5 complexes have all three")
    L.append("bands), and the top-1 -> top-deepest centroid drift with a receptor-bootstrap")
    L.append("95% CI (resamples receptors, not poses, to respect pseudoreplication).")
    L.append("")
    L.append(_fmt_table(stats.get("within_tool_trend")))
    L.append("")

    L.append("=" * 78)
    L.append("4. CROSS-TOOL COMPARISON ON COMMON COMPLEXES  (paired Friedman + Wilcoxon)")
    L.append("=" * 78)
    L.append("Per (depth, metric): the tools compared on the COMMON complex set only")
    L.append("(n_common), with an omnibus Friedman test and Wilcoxon signed-rank post-hoc")
    L.append("(Holm-adjusted p, with Cliff's delta effect size) for each tool pair.")
    L.append("")
    L.append(_fmt_table(stats.get("crosstool_paired")))
    L.append("")

    L.append("=" * 78)
    L.append("5. POSE DIVERSITY AT THE DEEPEST DEPTH  (mode-collapse vs genuine spread)")
    L.append("=" * 78)
    L.append("Within-complex spread of in-place RMSD across a tool's deepest cohort: a flat")
    L.append("depth profile is genuine ranking robustness only if the extra poses are")
    L.append("diverse; near-zero spread means near-duplicate (mode-collapsed) poses.")
    L.append("")
    L.append(_fmt_table(stats.get("pose_diversity")))
    L.append("")

    path.write_text("\n".join(L) + "\n", encoding="utf-8")


def main() -> None:
    if not CACHE.exists():
        sys.exit(f"per-pose cache not found: {CACHE}")
    df = P._read_cached_per_pose(CACHE)
    df, src = collapse_like_cell21(df)

    stem = "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15"
    png = REPORT / f"{stem}.png"
    caption = REPORT / f"{stem}__caption.txt"
    # Same code path as the sibling _pbvalid filmstrip; only `depths` differs.
    # caption_out moves the footnote off the figure into the sidecar .txt.
    # legend_loc="top" lifts the mechanism key below the title; title_mechanism_line
    # =False drops ONLY the second panel-title line ("place/mix/form = ..%") — the
    # mechanism breakdown moves to the __panel_data_and_stats.txt sidecar — while the
    # first line (family · top-d · n · rec.) stays on each panel; the font sizes
    # enlarge the panel headings and axis labelling.
    P.plot_form_vs_placement_depth_filmstrip(
        df, png, form_ok=P.FORM_OK_KABSCH_A, rmsd_gate=None, depths=DEPTHS,
        centroid_max=P.FAR_FROM_RECEPTOR_CENTROID_A, axis_max=AXIS_MAX,
        caption_out=caption, legend_loc="top", title_mechanism_line=False,
        title_fontsize=10, axis_fontsize=13, tick_fontsize=11)

    # Thesis Figure 5 layout. The two variants above move the footnote into a sidecar
    # and drop the per-panel mechanism line, but the figure embedded in the thesis
    # keeps BOTH on the image (footer block under the axes, "place/mix/form = ..%" as
    # the second title line) with the mechanism key inside the top-left panel. That is
    # simply the function's own defaults, so pass caption_out=None and leave
    # legend_loc / title_mechanism_line alone. Emitted separately so the sidecar
    # variants stay available and the thesis copy is reproducible from this script.
    stem_t = f"{stem}__thesis"
    png_t = REPORT / f"{stem_t}.png"
    P.plot_form_vs_placement_depth_filmstrip(
        df, png_t, form_ok=P.FORM_OK_KABSCH_A, rmsd_gate=None, depths=DEPTHS,
        centroid_max=P.FAR_FROM_RECEPTOR_CENTROID_A, axis_max=AXIS_MAX)

    cohorts = build_cohorts(df)
    tbl = per_pose_table(cohorts, src, axis_max=AXIS_MAX)
    summ = cohort_summary(cohorts, src)
    tbl.to_csv(REPORT / f"{stem}.csv", index=False)
    summ.to_csv(REPORT / f"{stem}__cohort_summary.csv", index=False)

    # Companion honest-statistics (the module the pipeline ships for its 1/3/5
    # filmstrip), recomputed for the 1/5/15 depths: valid-complex coverage,
    # within-tool depth trend (Spearman + paired Friedman + receptor-bootstrap drift),
    # cross-tool paired Friedman + Wilcoxon(Holm) + Cliff's δ, and pose diversity.
    stats = P.aggregate_filmstrip_statistics(
        df, depths=DEPTHS, centroid_max=P.FAR_FROM_RECEPTOR_CENTROID_A)
    if stats:
        for key, sdf in stats.items():
            sdf.to_csv(REPORT / f"{stem}__stats_{key}.csv", index=False)
        P.plot_filmstrip_statistics(stats, REPORT / f"{stem}__stats.png", depths=DEPTHS)
        print("wrote stats    → __stats.png + __stats_{per_tool_depth,within_tool_trend,"
              "crosstool_paired,pose_diversity}.csv")

    # Single readable sidecar: the per-panel data taken off the figure + every test.
    stats_txt = REPORT / f"{stem}__panel_data_and_stats.txt"
    write_stats_txt(stats_txt, summ, stats, src)

    # Additional graph: the SAME filmstrip with the per-panel numbers stripped entirely
    # from the headings. compact_titles=True reduces each panel title to just
    # "Family · top-d" — the n= (cohort pose count), the "· NN rec." (distinct receptors)
    # and the "(+NN poses vs top-1)" delta all move off the figure. Those numbers already
    # live in section 1 of the companion .txt (n_poses / n_receptors_protein / n_complexes),
    # which also carries every companion statistical test, so nothing is lost — the panel
    # just reads cleanly. Fonts enlarged relative to the sibling: panel headings, the
    # figure-level mechanism legend and the x/y axis labels/ticks all grow.
    stem_c = f"{stem}__compact_titles"
    png_c = REPORT / f"{stem_c}.png"
    caption_c = REPORT / f"{stem_c}__caption.txt"
    P.plot_form_vs_placement_depth_filmstrip(
        df, png_c, form_ok=P.FORM_OK_KABSCH_A, rmsd_gate=None, depths=DEPTHS,
        centroid_max=P.FAR_FROM_RECEPTOR_CENTROID_A, axis_max=AXIS_MAX,
        caption_out=caption_c, legend_loc="top", compact_titles=True,
        title_fontsize=16, axis_fontsize=16, tick_fontsize=13,
        legend_fontsize=15, legend_title_fontsize=16, suptitle_fontsize=16,
        shared_axis_labels=True)
    stats_txt_c = REPORT / f"{stem_c}__panel_data_and_stats.txt"
    write_stats_txt(stats_txt_c, summ, stats, src)

    print(f"\nwrote figure   → {png}")
    print(f"wrote caption  → {caption.name}")
    print(f"wrote txt      → {stats_txt.name}")
    print(f"wrote figure   → {png_c.name}  (compact titles, larger fonts)")
    print(f"wrote caption  → {caption_c.name}")
    print(f"wrote txt      → {stats_txt_c.name}")
    print(f"wrote per-pose → {stem}.csv  ({len(tbl):,} rows)")
    print(f"wrote summary  → {stem}__cohort_summary.csv")
    print("\nper-panel pose counts (family × depth):")
    for _, r in summ.iterrows():
        d = r["depth"].split("-")[1]
        assert (tbl["family"].eq(P._fam_key(r["family"].lower()))
                & tbl[f"in_top{d}"]).sum() == r["n_poses"]
        print(f"  {r['family']:9s} {r['depth']:7s} n={r['n_poses']:5d} "
              f"· {r['n_receptors_protein']} rec · "
              f"place/mix/form = {r['pct_placement_limited']:.0f}/"
              f"{r['pct_mixed']:.0f}/{r['pct_form_limited']:.0f}%")


if __name__ == "__main__":
    main()
