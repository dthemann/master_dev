#!/usr/bin/env python
"""Audit the pose basis of the Benchmark interaction-fingerprint analysis.

``run_pandamap.py`` selects poses on PoseBusters validity plus a rank-depth cap
(``poses_per_combo``). It applies NO accuracy gate — no ``rmsd <= 2 A`` and no
``bestfit_rmsd < 1 A``. That is deliberate and is stated in the thesis, but it
leaves four quantities undocumented that a reader of the interaction section
needs, and that ``pandamap_interaction_report.py`` never computes because it
does not join the per-pose RMSD table:

  A  COMPOSITION      what the profiled pool actually is, in RMSD terms.
  B  PLACEMENT SHARE  how much of the between-tool ordering is site-finding
                      rather than contact quality, by direct standardisation of
                      rank-1 F1 to a common RMSD-band distribution.
  C  SITE FILTER      the headline top-5-union Jaccard recomputed over poses
                      that reached the native site at all (``rmsd < 10 A``).
                      This is a SITE filter, not an accuracy gate: it removes
                      poses docked elsewhere on the receptor without selecting
                      on near-nativeness. It is a CONDITIONAL comparison, so the
                      per-tool count of complexes it drops is reported beside it.
  D  METAL AUDIT      PandaMap 4.1.0 emits ``metal_coordination`` position-
                      independently. This quantifies it and re-runs the headline
                      with the class removed, so the figure's robustness to the
                      defect is on record.

The numbers this writes are the ones quoted in the thesis appendix section
"Interaction Fingerprints". Re-run it whenever the arm is re-generated.

Run under the vina env (pandas/scipy), from the repo root:

    /home/manndo/anaconda3/envs/vina/bin/python \\
        Scripts/Analysis/interaction_pose_basis_audit.py \\
        --config Scripts/Analysis/pandamap_benchmark_full_protein_mgltools_config.yaml \\
        --exclude-methods unidock2,autodock_gnina

Writes ``<in-dir>/report/interaction_pose_basis_audit.{txt,csv}``.

Nearest-copy convention (opt-in, nearest-copy endpoint program): when
``crystal_interactions.csv`` carries ``copy_index`` (run_pandamap ``crystal_copies:
all``) AND the metrics table carries ``nearest_copy_index`` (hub
``--reference-convention nearest``), every top-5 union is scored against the union
of the fingerprints of the copies its poses are nearest to (plan v2 §1.8), joined on
POSE_KEY. With either column absent the audit is byte-identical to the historical
run (a per-copy crystal file against an instance table reduces to its reference copy).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
import method_filter as mf  # noqa: E402

# Pose key shared by pandamap_pose_summary.csv, recovery_detail_per_pose.csv and
# the pose_comparison per_pose_metrics.csv. The join is exact on all four fields;
# a partial join means the arms are mismatched and the audit aborts.
POSE_KEY = ["method", "protein", "ligand", "pose_name"]

# Charged classes PandaMap writes three times for one physical contact, and the
# class it emits position-independently. Both are appendix-documented defects.
CHARGED = ("ionic", "salt_bridge", "attractive_charge")
METAL = "metal_coordination"

# RMSD bands used for stratification and for the standardisation reference.
BANDS = [(0.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, np.inf)]

# Alternative schemes, to show the placement share is not an artefact of the banding.
ALT_BANDS = [
    ("3-band", [0, 2, 10, np.inf]),
    ("5-band (reported)", [0, 2, 5, 10, 20, np.inf]),
    ("8-band", [0, 1, 2, 3, 5, 8, 12, 20, np.inf]),
    ("14-band", [0, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15, 20, 30, np.inf]),
]

SITE_FILTER_A = 10.0   # "reached the native site" — see module docstring, section C
SITE_FILTER_SWEEP = (8.0, 10.0, 12.0, 15.0, 20.0)   # the threshold must not be load-bearing
NEAR_NATIVE_A = 2.0


def band_of(rmsd: float) -> str:
    for lo, hi in BANDS:
        if lo <= rmsd < hi:
            return f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}"
    return "nan"


def typed_fp(g: pd.DataFrame, *, collapse_charged: bool = False,
             drop_metal: bool = False) -> set:
    """Typed interaction fingerprint, mirroring pandamap_interaction_report.pose_fingerprint.

    Key is (interaction_type, resname, resnum, chain). The two optional repairs
    exist only to measure the defects' effect; the default reproduces the figure.
    """
    t = g["interaction_type"].astype(str)
    if collapse_charged:
        t = t.where(~t.isin(CHARGED), "charged")
    keep = np.ones(len(g), dtype=bool) if not drop_metal else (t != METAL).to_numpy()
    return set(zip(t[keep], g["resname"][keep], g["resnum"][keep], g["chain"][keep]))


def jaccard(a: set, b: set) -> float:
    u = len(a | b)
    return len(a & b) / u if u else np.nan


class CopySelector:
    """Nearest-copy crystal reference: per-copy crystal rows + per-pose nearest copy.

    Mirrors pandamap_interaction_report.CrystalCopySelector for the audit's two
    call sites. ``None`` (see :func:`build_copy_selector`) keeps the per-pair reference.
    """

    def __init__(self, crystal: pd.DataFrame, nearest: dict, ref_index: dict):
        ck = crystal["copy_index"].fillna(0).astype(int)
        self.copies = {(p, l, k): g for (p, l, k), g in
                       crystal.assign(_k=ck).groupby(["protein", "ligand", "_k"])}
        self.pairs = {(p, l) for (p, l, _k) in self.copies}
        self.nearest = nearest
        self.ref_index = ref_index
        self.fallback: set = set()

    def copy_of(self, m, p, l, pose_name) -> int:
        k = self.nearest.get((m, p, l, pose_name))
        if k is None:
            self.fallback.add((m, p, l, pose_name))
            return self.ref_index.get((p, l), 0)
        return k

    def union_fp(self, m, p, l, pose_names, **fp_kw) -> set:
        """Union of the (optionally repaired) fingerprints of every selected copy."""
        out: set = set()
        for k in sorted({self.copy_of(m, p, l, n) for n in pose_names}):
            g = self.copies.get((p, l, k))
            if g is not None:
                out |= typed_fp(g, **fp_kw)
        return out


def build_copy_selector(crystal: pd.DataFrame, metrics: pd.DataFrame) -> "CopySelector | None":
    """CopySelector when BOTH ``copy_index`` and ``nearest_copy_index`` exist, else None."""
    if "copy_index" not in crystal.columns or "nearest_copy_index" not in metrics.columns:
        return None
    mt = metrics.dropna(subset=["nearest_copy_index"]).drop_duplicates(POSE_KEY)
    nearest = {tuple(k): int(v) for *k, v in
               mt[POSE_KEY + ["nearest_copy_index"]].itertuples(index=False, name=None)}
    ref_index = {}
    if "ref_copy_index" in mt.columns:
        r = mt.dropna(subset=["ref_copy_index"]).drop_duplicates(["protein", "ligand"])
        ref_index = {(p, l): int(v) for p, l, v in
                     r[["protein", "ligand", "ref_copy_index"]].itertuples(index=False, name=None)}
    return CopySelector(crystal, nearest, ref_index)


def union_jaccard(inter: pd.DataFrame, crystal_fp: dict, depth: int,
                  copy_sel: "CopySelector | None" = None, **fp_kw) -> dict[str, dict]:
    """Mean Jaccard of each method's top-`depth` union against the crystal.

    Returns {method: {(protein, ligand): jaccard}}. Selection mirrors
    plot_fingerprint_similarity: rank-order, head(depth), union the rows.
    ``copy_sel`` (nearest-copy convention) replaces ``crystal_fp[(p, l)]`` by the union
    of the selected copies' fingerprints, repaired with the same ``fp_kw``.
    """
    keep = (inter[POSE_KEY + ["pose_rank"]].drop_duplicates()
            .sort_values("pose_rank")
            .groupby(["method", "protein", "ligand"], sort=False)
            .head(depth))
    sel = inter.merge(keep[POSE_KEY], on=POSE_KEY, how="inner")
    out: dict[str, dict] = {}
    for (m, p, l), g in sel.groupby(["method", "protein", "ligand"], sort=False):
        if copy_sel is not None:
            ref = copy_sel.union_fp(m, p, l, g["pose_name"].unique(), **fp_kw) \
                if (p, l) in copy_sel.pairs else None
        else:
            ref = crystal_fp.get((p, l))
        if ref:
            out.setdefault(m, {})[(p, l)] = jaccard(typed_fp(g, **fp_kw), ref)
    return out


def dd_vs_ad(per: dict, methods) -> tuple[set, float, set, float]:
    """DiffDock-vs-AutoDock Wilcoxon on both the all-tool cohort and the pairwise one.

    The thesis compares tools "on the complexes they share", so the all-tool cohort is the
    reportable one. The pairwise cohort is returned too because it is larger and the two
    p-values differ enough that quoting one beside the other's n misstates the test.
    """
    ad = next((m for m in methods if m.startswith("autodock")), None)
    dd = next((m for m in methods if m.startswith("diffdock")), None)
    common = set.intersection(*(set(per.get(m, {})) for m in methods)) if per else set()
    pair = (set(per.get(ad, {})) & set(per.get(dd, {}))) if ad and dd else set()

    def _p(cohort):
        if not (ad and dd) or len(cohort) <= 8:
            return float("nan")
        return float(wilcoxon([per[dd][x] for x in cohort], [per[ad][x] for x in cohort]).pvalue)

    return common, _p(common), pair, _p(pair)


def standardise(df: pd.DataFrame, value: str, reference: pd.Series) -> float:
    """Direct standardisation of `value` to the `reference` RMSD-band distribution."""
    means = df.groupby("band", observed=False)[value].mean()
    w = reference.reindex(means.index).fillna(0.0)
    ok = means.notna() & (w > 0)
    return float((means[ok] * w[ok]).sum() / w[ok].sum()) if ok.any() else np.nan


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, required=True,
                    help="pandamap config; supplies output_dir and pb_csv")
    ap.add_argument("--in-dir", type=Path, default=None,
                    help="override the config's output_dir")
    ap.add_argument("--metrics", type=Path, default=None,
                    help="per_pose_metrics.csv; default derived from the config's pb_csv")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="default <in-dir>/report")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="restrict to these '<PDBID>_<LIG>' complex ids, one per line "
                         "('#' comments ignored). Must match the restriction the pandamap "
                         "run and its report used, or the RMSD join guard fires.")
    mf.add_method_filter_args(ap)
    args = ap.parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text())
    root = Path(cfg.get("work_dir") or ".")
    in_dir = args.in_dir or (root / cfg["output_dir"])
    out_dir = args.out_dir or (in_dir / "report")
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.metrics is not None:
        metrics_csv = args.metrics
    else:
        # pb_csv lives at <dock>/posebusters_filtered_results.csv; the pose-comparison
        # report that carries rmsd/bestfit_rmsd sits beside it.
        metrics_csv = (root / cfg["pb_csv"]).parent / "pose_comparison_report" / "per_pose_metrics.csv"

    inter = pd.read_csv(in_dir / "pandamap_interactions.csv", low_memory=False)
    crystal = pd.read_csv(in_dir / "crystal_interactions.csv")
    # An INPUT written by pandamap_interaction_report.py, so it is read from the arm's own
    # report dir and not from --out-dir, which the caller may redirect anywhere.
    recov = pd.read_csv(in_dir / "report" / "recovery_detail_per_pose.csv")
    metrics = pd.read_csv(metrics_csv, low_memory=False)

    # Same semantics as run_pandamap._load_allowed_ids: the complex id is the 'protein' field.
    if args.ids_file is not None:
        allowed = {ln.strip() for ln in args.ids_file.read_text().splitlines()
                   if ln.strip() and not ln.lstrip().startswith("#")}
        before = len(inter)
        inter = inter[inter["protein"].isin(allowed)]
        recov = recov[recov["protein"].isin(allowed)]
        crystal = crystal[crystal["protein"].isin(allowed)]
        print(f"Restricted to {len(allowed)} ids from {args.ids_file.name}: "
              f"{before} → {len(inter)} interaction rows")

    patterns = mf.resolve_patterns(args)
    if patterns:
        inter = mf.apply_patterns(inter, "method", patterns, label="interactions")
        recov = mf.apply_patterns(recov, "method", patterns, label="recovery")

    # Nearest-copy crystal reference (opt-in; see module docstring). Decided BEFORE the
    # metrics frame is reduced to its RMSD columns.
    copy_sel = build_copy_selector(crystal, metrics)
    if copy_sel is None and "copy_index" in crystal.columns:
        crystal = crystal[crystal["copy_index"].fillna(0).astype(int) == 0].drop(columns=["copy_index"])
        print("crystal_interactions.csv carries copy_index but the metrics table has no "
              "nearest_copy_index — using the reference copy only (single-instance convention).")

    rm = metrics[POSE_KEY + ["rmsd", "bestfit_rmsd"]].drop_duplicates(POSE_KEY)
    inter = inter.merge(rm, on=POSE_KEY, how="left")
    recov = recov.merge(rm, on=POSE_KEY, how="left")

    poses = inter[POSE_KEY + ["pose_rank", "rmsd", "bestfit_rmsd"]].drop_duplicates(POSE_KEY)
    # Guard BOTH frames. Sections A and C run on `inter`, B and D on `recov`, and the two are
    # written by different steps — the interactions CSV is unrestricted while the recovery table
    # can be id-restricted. An unguarded gap there would silently bin poses into a "nan" band and
    # contaminate the standardisation reference.
    for label, frame in [("interactions", poses), ("recovery", recov)]:
        n_join = int(frame["rmsd"].notna().sum())
        if n_join != len(frame):
            raise SystemExit(
                f"{label}/metrics join is partial ({n_join}/{len(frame)}). The pandamap arm and the "
                f"pose-comparison report are out of step — regenerate one of them before auditing."
            )
    n_join = len(poses)
    poses = poses.assign(band=poses["rmsd"].map(band_of))
    recov = recov.assign(band=recov["rmsd"].map(band_of))
    methods = sorted(poses["method"].unique())

    crystal_fp = {(p, l): typed_fp(g) for (p, l), g in crystal.groupby(["protein", "ligand"])}
    L, kv = [], []

    def emit(line=""):
        L.append(line)

    def record(key, val):
        kv.append({"metric": key, "value": val})

    emit(f"Interaction-analysis pose-basis audit — {in_dir}")
    emit(f"methods: {', '.join(methods)}")
    if patterns:
        emit(f"excluded: {', '.join(patterns)}")
    emit(f"metrics : {metrics_csv}")
    emit(f"join    : {n_join}/{len(poses)} poses carry an RMSD")
    if copy_sel is not None:
        n_multi = len({(p, l) for (p, l, k) in copy_sel.copies if k != copy_sel.ref_index.get((p, l), 0)})
        emit(f"crystal : nearest-copy convention — each top-5 union is scored against the union of "
             f"the crystal copies its poses are nearest to ({n_multi} of {len(copy_sel.pairs)} pairs "
             f"carry >1 fingerprinted copy)")
        record("crystal_reference_convention", "nearest")
        record("crystal_pairs_multi_copy", n_multi)
    emit()

    # ── A. composition of the profiled pool ────────────────────────────────
    emit("A. COMPOSITION OF THE PROFILED POOL (no accuracy gate is applied)")
    emit(f"{'cohort':<34}{'n':>7}{'median RMSD':>13}{'<=2 A':>9}{'Kabsch<1':>10}{'>5 A':>8}")
    for label, sub in [("all profiled poses", poses),
                       ("rank-1 poses", poses[poses["pose_rank"] == 1])]:
        emit(f"{label:<34}{len(sub):>7}{sub['rmsd'].median():>12.2f}A"
             f"{(sub['rmsd'] <= NEAR_NATIVE_A).mean() * 100:>8.1f}%"
             f"{(sub['bestfit_rmsd'] < 1).mean() * 100:>9.1f}%"
             f"{(sub['rmsd'] > 5).mean() * 100:>7.1f}%")
    record("pool_median_rmsd_A", round(float(poses["rmsd"].median()), 2))
    record("pool_within_2A_pct", round(float((poses["rmsd"] <= NEAR_NATIVE_A).mean() * 100), 1))
    record("pool_n_poses", len(poses))
    emit()
    for m in methods:
        s = poses[poses["method"] == m]
        emit(f"  {m:<32}{len(s):>7}{s['rmsd'].median():>12.2f}A"
             f"{(s['rmsd'] <= NEAR_NATIVE_A).mean() * 100:>8.1f}%")
    emit()

    # ── B. how much of the ordering is placement ───────────────────────────
    emit("B. PLACEMENT SHARE OF THE RANK-1 CONTACT-RECOVERY ORDERING")
    emit("   Direct standardisation of rank-1 F1 to the pooled RMSD-band distribution.")
    r1 = recov[recov["pose_rank"] == 1]
    reference = r1["band"].value_counts(normalize=True)
    emit(f"{'method':<34}{'crude F1':>10}{'standardised':>14}{'n':>6}")
    crude, std = {}, {}
    for m in methods:
        s = r1[r1["method"] == m]
        crude[m] = float(s["f1"].mean())
        std[m] = standardise(s, "f1", reference)
        emit(f"  {m:<32}{crude[m]:>10.4f}{std[m]:>14.4f}{len(s):>6}")
    spread_c = max(crude.values()) - min(crude.values())
    spread_s = max(std.values()) - min(std.values())
    share = (1 - spread_s / spread_c) * 100 if spread_c else np.nan
    emit(f"  spread {spread_c:.4f} -> {spread_s:.4f}: {share:.1f}% of the ordering is placement")
    record("rank1_spread_crude", round(spread_c, 4))
    record("rank1_spread_standardised", round(spread_s, 4))
    record("rank1_placement_share_pct", round(float(share), 1))
    emit()
    # The placement share must not be an artefact of the band scheme, and the standardised
    # values of the two leading tools are close enough that their ORDER is scheme-dependent.
    # Both facts are quoted in the appendix, so both are computed rather than asserted.
    emit("   Sensitivity of that share to the band scheme:")
    for name, edges in ALT_BANDS:
        vals = {}
        for m in methods:
            s = r1[r1["method"] == m]
            bands = pd.cut(s["rmsd"], edges, right=False)
            ref = pd.cut(r1["rmsd"], edges, right=False).value_counts(normalize=True)
            means = s.groupby(bands, observed=False)["f1"].mean()
            w = ref.reindex(means.index).fillna(0.0)
            ok = means.notna() & (w > 0)
            vals[m] = float((means[ok] * w[ok]).sum() / w[ok].sum())
        sh = (1 - (max(vals.values()) - min(vals.values())) / spread_c) * 100
        order = " > ".join(m.split("_")[0][:8] for m in sorted(vals, key=vals.get, reverse=True))
        emit(f"     {name:<22}{sh:>6.1f}%   standardised order {order}")
        record(f"placement_share_{name.replace(' ', '_')}_pct", round(sh, 1))
    emit()
    emit("   Per-band rank-1 F1 (the tools at matched placement):")
    piv = r1.pivot_table(index="band", columns="method", values="f1", aggfunc="mean")
    cnt = r1.pivot_table(index="band", columns="method", values="f1", aggfunc="size")
    for b in [f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}" for lo, hi in BANDS]:
        if b not in piv.index:
            continue
        cells = "  ".join(f"{piv.loc[b, m]:.3f} (n={int(cnt.loc[b, m] or 0):>4})"
                          if pd.notna(piv.loc[b, m]) else f"{'-':>14}" for m in methods)
        emit(f"     {b:>8} A   {cells}")
    emit()

    # ── C. site filter ─────────────────────────────────────────────────────
    emit(f"C. SITE-FILTER SENSITIVITY OF THE TOP-5 UNION JACCARD (rmsd < {SITE_FILTER_A:g} A)")
    emit("   A site filter, not an accuracy gate. Conditional: the complexes each tool")
    emit("   never reached are dropped, so the dropped counts are reported beside it.")
    emit("   The Wilcoxon runs on the cohort every tool defines, which is the convention the")
    emit("   Results chapter uses. The pairwise AutoDock-DiffDock cohort is printed beside it,")
    emit("   because the two differ and quoting one against the other's n misreads the test.")
    emit(f"{'basis':<26}" + "".join(f"{m.split('_')[0][:8]:>18}" for m in methods)
         + f"{'common n':>10}{'p common':>11}{'pair n':>8}{'p pair':>11}")
    results = {}
    for label, sub in [("published (all poses)", inter),
                       (f"rmsd < {SITE_FILTER_A:g} A", inter[inter["rmsd"] < SITE_FILTER_A])]:
        per = union_jaccard(sub, crystal_fp, depth=5, copy_sel=copy_sel)
        results[label] = per
        common, p_com, pair, p_pair = dd_vs_ad(per, methods)
        cells = "".join(f"{np.mean(list(per.get(m, {}).values())):>10.4f} ({len(per.get(m, {})):>3})"
                        for m in methods)
        emit(f"{label:<26}{cells}{len(common):>10}{p_com:>11.2e}{len(pair):>8}{p_pair:>11.2e}")
        tag = label.split()[0]
        record(f"jaccard_{tag}_common_n", len(common))
        record(f"jaccard_{tag}_p_dd_vs_ad_common", float(p_com))
        record(f"jaccard_{tag}_pairwise_n", len(pair))
        record(f"jaccard_{tag}_p_dd_vs_ad_pairwise", float(p_pair))
        for m in methods:
            record(f"jaccard_{tag}_{m}", round(float(np.mean(list(per.get(m, {}).values()))), 4))
    base, filt = results["published (all poses)"], results[f"rmsd < {SITE_FILTER_A:g} A"]
    emit("   complexes dropped by the site filter (a tool never reached the site):")
    for m in methods:
        dropped = len(base.get(m, {})) - len(filt.get(m, {}))
        emit(f"     {m:<34}{dropped:>4}")
        record(f"site_filter_dropped_{m}", dropped)
    # The filter is a SITE filter, but it is not accuracy-neutral. State the enrichment rather
    # than let the text imply the retained pool is no more near-native than the whole pool.
    kept = poses[poses["rmsd"] < SITE_FILTER_A]
    enrich = float((kept["rmsd"] <= NEAR_NATIVE_A).mean() * 100)
    emit(f"   share of poses within {NEAR_NATIVE_A:g} A rises from "
         f"{(poses['rmsd'] <= NEAR_NATIVE_A).mean() * 100:.1f}% to {enrich:.1f}% under the filter")
    record("site_filter_within_2A_pct", round(enrich, 1))
    emit("   threshold sweep (the 10 A cut must not be load-bearing):")
    for thr in SITE_FILTER_SWEEP:
        per = union_jaccard(inter[inter["rmsd"] < thr], crystal_fp, depth=5, copy_sel=copy_sel)
        common, p_com, _, _ = dd_vs_ad(per, methods)
        cells = "  ".join(f"{np.mean(list(per.get(m, {}).values())):.4f}" for m in methods)
        emit(f"     < {thr:>4.0f} A   {cells}   common n={len(common):>4}   p={p_com:.2e}")
        record(f"site_sweep_{thr:g}A_p_common", float(p_com))
    emit()

    # ── D. metal-coordination position independence ────────────────────────
    emit("D. metal_coordination POSITION INDEPENDENCE (PandaMap 4.1.0 detection defect)")
    mc = inter[inter["interaction_type"] == METAL]
    if mc.empty:
        emit("   no metal-coordination contacts in this dataset — defect not applicable")
        record("metal_rows", 0)
    else:
        # "Every pose carries the same key-set" is only meaningful if every pose in the complex
        # actually records metal. A test whose denominator is the metal-bearing poses alone
        # cannot be falsified by a pose that recorded none, so coverage is checked separately.
        seen = poses.groupby(["protein", "ligand"])["pose_name"].nunique()
        same = uncovered = 0
        n_cplx = mc.groupby(["protein", "ligand"]).ngroups
        for (p, l), g in mc.groupby(["protein", "ligand"]):
            keysets = {frozenset(zip(gg["resname"], gg["resnum"], gg["chain"]))
                       for _, gg in g.groupby("pose_name")}
            same += len(keysets) == 1
            uncovered += g["pose_name"].nunique() < int(seen.get((p, l), 0))
        emit(f"   {len(mc)} rows over {mc['pose_name'].nunique()} poses in {n_cplx} complexes")
        emit(f"   complexes where EVERY pose carries an identical metal key-set: {same}/{n_cplx}")
        emit(f"   complexes where some profiled pose records NO metal contact: {uncovered}/{n_cplx}")
        record("metal_rows", len(mc))
        record("metal_complexes", n_cplx)
        record("metal_complexes_pose_invariant", same)
        record("metal_complexes_with_uncovered_pose", uncovered)

        has_metal = set(map(tuple, mc[POSE_KEY].drop_duplicates().to_numpy()))
        far = recov[recov["rmsd"] > 20].copy()
        far["has_metal"] = [tuple(x) in has_metal for x in far[POSE_KEY].to_numpy()]
        emit("   poses beyond 20 A, split on whether the pose carries a metal contact:")
        for flag, g in far.groupby("has_metal"):
            emit(f"     metal={str(flag):<5} n={len(g):>5}  mean F1 {g['f1'].mean():.4f}"
                 f"  share scoring exactly zero {(g['f1'] == 0).mean() * 100:>5.1f}%"
                 f"  recovering any native contact {int((g['f1'] > 0).sum()):>4}")
            record(f"far_field_metal_{flag}_mean_f1", round(float(g["f1"].mean()), 4))
            record(f"far_field_metal_{flag}_zero_pct", round(float((g["f1"] == 0).mean() * 100), 1))
            record(f"far_field_metal_{flag}_n", len(g))
            record(f"far_field_metal_{flag}_any_contact", int((g["f1"] > 0).sum()))
        emit()
        emit("   effect of the two detection defects on the published top-5 union Jaccard:")
        for label, kw in [("published", {}),
                          ("charged classes collapsed", {"collapse_charged": True}),
                          ("metal class dropped", {"drop_metal": True}),
                          ("both repaired", {"collapse_charged": True, "drop_metal": True})]:
            # The repairs apply to BOTH sides of the comparison. Repairing only the pose
            # fingerprint would inflate the union and understate every similarity.
            fp_ref = {(p, l): typed_fp(g, **kw) for (p, l), g in crystal.groupby(["protein", "ligand"])}
            per = union_jaccard(inter, fp_ref, depth=5, copy_sel=copy_sel, **kw)
            common, p_com, pair, p_pair = dd_vs_ad(per, methods)
            cells = "  ".join(f"{np.mean(list(per.get(m, {}).values())):.4f}" for m in methods)
            emit(f"     {label:<28}{cells}   p(DD vs AD) = {p_com:.4f} on the common "
                 f"n={len(common)} ({p_pair:.4f} pairwise n={len(pair)})")
            key = label.replace(" ", "_")
            record(f"defect_{key}_p_dd_vs_ad_common", round(float(p_com), 4))
            record(f"defect_{key}_p_dd_vs_ad_pairwise", round(float(p_pair), 4))
            for m in methods:
                record(f"defect_{label.replace(' ', '_')}_{m}",
                       round(float(np.mean(list(per.get(m, {}).values()))), 4))

    if copy_sel is not None and copy_sel.fallback:
        emit(f"   note: {len(copy_sel.fallback)} profiled poses carry no nearest_copy_index and were "
             "scored against the reference copy")
        record("crystal_poses_without_nearest_index", len(copy_sel.fallback))

    txt = out_dir / "interaction_pose_basis_audit.txt"
    csv = out_dir / "interaction_pose_basis_audit.csv"
    txt.write_text("\n".join(L) + "\n")
    pd.DataFrame(kv).to_csv(csv, index=False)
    print("\n".join(L))
    print(f"\nwrote {txt}\nwrote {csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
