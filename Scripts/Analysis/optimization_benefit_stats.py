#!/usr/bin/env python
"""Statistical tests for the optimization-benefit-by-rank analysis (figs 15c / 15d).

The plotted CSV (``optimization_benefit_by_rank.csv``) holds only marginal per-rank
percentages + n. Paired categorical inference lives in the *discordant* cells (how many
individual poses flip fail->pass vs pass->fail under optimization), which two marginal
percentages cannot reconstruct. So every test here goes back to the per-pose 0/1 outcomes
in ``per_pose_metrics.csv`` (the same rows fig 15c is aggregated from) and rebuilds the
raw / smina / gnina pose triples exactly as ``aggregate_optimization_benefit`` does.

Two pre-specified questions (Benchmark only — near-native needs the crystal):

(a) DOES OPTIMIZATION HELP?  raw vs smina vs gnina are the SAME poses measured three ways,
    so this is a paired-proportion problem. Per (tool, endpoint) at a fixed rank:
      * Cochran's Q omnibus across {raw, smina, gnina}
      * pairwise exact McNemar (raw-smina, raw-gnina, smina-gnina), Holm within the triple
      * Wilson CIs per rate + the paired risk difference (optimized - raw) in pp.
    Rank 1 (the operational top-1 pose) is the primary, clustering-proof test: at one fixed
    rank each complex contributes exactly one triple, so rows are independent across
    complexes. Ranks {5, 15, 30} are pre-specified supporting depths (NOT a 30-way sweep).
    Omnibus family is BH-FDR corrected. both_% is deterministic (near AND pbv) so it is
    reported descriptively, never tested as a third outcome.

(b) DOES THE BENEFIT DEPEND ON RANK?  the same 303 complexes recur at every rank, so this
    is an optimizer x rank interaction under clustering. Per (tool, endpoint, contrast):
      logit P(pass) = b0 + b1*opt + b2*rank + b3*(opt x rank)
    fitted on the pooled per-pose data; b3 is the single coefficient that tests
    rank-dependence. Inference is a CLUSTER (complex) bootstrap (resample whole complexes)
    because the vina env has no statsmodels/GEE sandwich. The raw rank slope b2 is reported
    on its own (raw-only fit): for DiffDock it is the confidence-rank quality gradient; for
    EquiBind (generation order) it is the pre-specified NEGATIVE CONTROL, expected ~flat.
    Group pp-gaps (early 1-10 vs late 21-30) accompany b3 for human-readable direction.

Runs in the `vina` conda env (scipy + sklearn; hand-rolled paired tests via stats_utils).
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)
from stats_utils import (bh_fdr, fmt_p, holm, p_stars, paired_proportions,  # noqa: E402
                         wilson_ci)

# ── how the raw/smina/gnina pose triples are keyed (mirrors _OPT_BENEFIT_TOOLS) ──
TOOLS = [
    ("DiffDock", "native",
     {"raw": "diffdock", "smina": "diffdock_smina", "gnina": "diffdock_gnina"}),
    ("EquiBind", "generation",
     {"raw": "equibind_unguided_raw", "smina": "equibind_unguided_smina",
      "gnina": "equibind_unguided_gnina"}),
]
ENDPOINTS = [("near_native", "Near-native (RMSD <= 2 A)"),
             ("pb_valid", "PoseBusters-valid")]
OPTIMIZERS = ["raw", "smina", "gnina"]
CONTRASTS = ["smina", "gnina"]          # each vs raw
PRIMARY_RANK = 1
SUPPORT_RANKS = [5, 15, 30]             # pre-specified, not a 30-way sweep
MAX_RANK = 30
EARLY, LATE = (1, 10), (21, 30)         # rank groups for the descriptive pp-gap
NEAR_THR = 2.0
_EQ_POSE_IDX_RE = re.compile(r"_(\d{2,3})__ref", re.IGNORECASE)


def _equibind_pose_index(name) -> int:
    m = _EQ_POSE_IDX_RE.search(str(name))
    return int(m.group(1)) if m else 10 ** 6


def _to_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(bool)
    return (s.astype(str).str.strip().str.lower()
            .map({"true": True, "1": True, "1.0": True, "false": False,
                  "0": False, "0.0": False, "nan": False, "": False})
            .fillna(False).astype(bool))


def build_per_pose(df: pd.DataFrame, max_rank: int = MAX_RANK) -> pd.DataFrame:
    """Long per-pose frame: tool, optimizer, complex, rank, near_native(0/1), pb_valid(0/1).

    Reproduces fig 15c's ranking: DiffDock by confidence rank (dedup the double-written
    rank-1), EquiBind by generation-order position parsed from the pose name.
    """
    df = df.copy()
    df["method"] = df["method"].astype(str)
    out = []
    for tool, kind, trio in TOOLS:
        for opt, mkey in trio.items():
            sub = df[df["method"] == mkey].copy()
            if sub.empty:
                continue
            if kind == "native":
                sub["_rk"] = pd.to_numeric(sub["rank"], errors="coerce")
                sub = sub.sort_values("_rk").drop_duplicates(["protein", "ligand", "_rk"])
            else:
                sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
                sub = sub.sort_values(["protein", "ligand", "_gi"], kind="mergesort")
                sub["_rk"] = sub.groupby(["protein", "ligand"]).cumcount() + 1
            sub = sub.dropna(subset=["_rk"])
            sub["_rk"] = sub["_rk"].astype(int)
            sub = sub[(sub["_rk"] >= 1) & (sub["_rk"] <= max_rank)]
            if sub.empty:
                continue
            near = (pd.to_numeric(sub["rmsd"], errors="coerce") <= NEAR_THR).astype(float)
            pbv = _to_bool(sub["pb_valid"]).astype(float)
            out.append(pd.DataFrame({
                "tool": tool, "optimizer": opt,
                "complex": (sub["protein"].astype(str) + "|" + sub["ligand"].astype(str)).values,
                "rank": sub["_rk"].values,
                "near_native": near.values, "pb_valid": pbv.values}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# ─────────────────────────────────────────────────────────────────────────────
# (a) does optimization help?  — Cochran's Q + pairwise exact McNemar at a rank
# ─────────────────────────────────────────────────────────────────────────────
def question_a(pp: pd.DataFrame) -> pd.DataFrame:
    ranks = [PRIMARY_RANK] + SUPPORT_RANKS
    rows = []
    for tool, _kind, _trio in TOOLS:
        for ep, _lab in ENDPOINTS:
            for rk in ranks:
                d = pp[(pp["tool"] == tool) & (pp["rank"] == rk)]
                wide = d.pivot_table(index="complex", columns="optimizer",
                                     values=ep, aggfunc="mean")
                wide = wide.reindex(columns=OPTIMIZERS).dropna()   # listwise-complete triples
                if len(wide) < 3:
                    continue
                res = paired_proportions({o: wide[o].to_numpy() for o in OPTIMIZERS},
                                         labels=OPTIMIZERS)
                rate = {o: res["rates"][o][0] for o in OPTIMIZERS}
                pw = {tuple(sorted((p["a"], p["b"]))): p for p in res["pairwise"]}
                base = {"tool": tool, "endpoint": ep, "rank": rk, "n": res["n"],
                        "primary": rk == PRIMARY_RANK,
                        "rate_raw": rate["raw"], "rate_smina": rate["smina"],
                        "rate_gnina": rate["gnina"],
                        "cochran_Q": res["omnibus"]["Q"], "cochran_p": res["omnibus"]["p"]}
                for opt in CONTRASTS:
                    p = pw[tuple(sorted(("raw", opt)))]
                    # discordant oriented as optimized-gains vs raw-gains
                    opt_gain = p["a_wins"] if p["a"] == opt else p["b_wins"]
                    raw_gain = p["b_wins"] if p["a"] == opt else p["a_wins"]
                    base[f"{opt}_minus_raw_pp"] = 100 * (rate[opt] - rate["raw"])
                    base[f"{opt}_opt_gain"] = opt_gain      # #poses raw-fail -> opt-pass
                    base[f"{opt}_raw_gain"] = raw_gain      # #poses raw-pass -> opt-fail
                    base[f"{opt}_mcnemar_p"] = p["p_raw"]
                    base[f"{opt}_mcnemar_holm"] = p["p_holm"]
                    base[f"{opt}_star"] = p["star"]
                rows.append(base)
    out = pd.DataFrame(rows)
    if not out.empty:                                    # BH-FDR across the omnibus family
        out["cochran_p_bh"] = bh_fdr(out["cochran_p"].to_numpy())
    return out


# ─────────────────────────────────────────────────────────────────────────────
# (b) does the benefit depend on rank?  — logistic opt x rank, cluster-bootstrapped
# ─────────────────────────────────────────────────────────────────────────────
def _fit_logit(X: np.ndarray, y: np.ndarray):
    from sklearn.linear_model import LogisticRegression
    if np.unique(y).size < 2:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = LogisticRegression(penalty=None, solver="lbfgs", max_iter=2000)
        m.fit(X, y)
    return m.coef_.ravel()


def _cluster_boot(rows_by_complex, complexes, build_fit, coef_idx,
                  n_boot=1000, seed=0):
    """Cluster (complex) bootstrap of one logistic coefficient.

    rows_by_complex : {complex -> row index array}; build_fit(idx) -> coef vector or None.
    Returns (point_est, lo, hi, boot_p) with a two-sided percentile-bootstrap p on 0.
    """
    full_idx = np.concatenate([rows_by_complex[c] for c in complexes])
    est_coef = build_fit(full_idx)
    if est_coef is None:
        return np.nan, np.nan, np.nan, np.nan
    est = float(est_coef[coef_idx])
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(complexes, complexes.size, replace=True)
        idx = np.concatenate([rows_by_complex[c] for c in pick])
        c = build_fit(idx)
        if c is not None:
            boots.append(c[coef_idx])
    boots = np.asarray(boots, float)
    if boots.size < max(50, n_boot // 4):
        return est, np.nan, np.nan, np.nan
    lo, hi = np.percentile(boots, [2.5, 97.5])
    n = boots.size
    p = 2 * min((np.sum(boots <= 0) + 1) / (n + 1),
                (np.sum(boots >= 0) + 1) / (n + 1))
    return est, float(lo), float(hi), float(min(p, 1.0))


def _pooled_rate(pp, tool, ep, opt, lo, hi):
    d = pp[(pp["tool"] == tool) & (pp["optimizer"] == opt)
           & (pp["rank"] >= lo) & (pp["rank"] <= hi)]
    return float(d[ep].mean()) if len(d) else np.nan


def question_b(pp: pd.DataFrame, n_boot=1000, seed=0) -> pd.DataFrame:
    from sklearn.linear_model import LogisticRegression  # noqa: F401 (import cost check)
    rows = []
    for tool, _kind, _trio in TOOLS:
        for ep, _lab in ENDPOINTS:
            # raw-only rank slope (confidence gradient / negative control)
            raw = pp[(pp["tool"] == tool) & (pp["optimizer"] == "raw")].reset_index(drop=True)
            rc = (raw["rank"].to_numpy(float) - 1.0)
            Xr = rc.reshape(-1, 1)
            yr = raw[ep].to_numpy(float)
            rbc = {c: np.where(raw["complex"].to_numpy() == c)[0]
                   for c in np.unique(raw["complex"].to_numpy())}
            comps_r = np.array(list(rbc.keys()))
            raw_slope, rs_lo, rs_hi, rs_p = _cluster_boot(
                rbc, comps_r,
                lambda idx: _fit_logit(Xr[idx], yr[idx]), 0, n_boot, seed)

            for opt in CONTRASTS:
                sub = pp[(pp["tool"] == tool)
                         & (pp["optimizer"].isin(["raw", opt]))].reset_index(drop=True)
                oi = (sub["optimizer"].to_numpy() == opt).astype(float)
                rkc = sub["rank"].to_numpy(float) - 1.0
                X = np.column_stack([oi, rkc, oi * rkc])   # [opt, rank, opt x rank]
                y = sub[ep].to_numpy(float)
                comp = sub["complex"].to_numpy()
                rbc2 = {c: np.where(comp == c)[0] for c in np.unique(comp)}
                comps = np.array(list(rbc2.keys()))
                b3, lo, hi, bp = _cluster_boot(
                    rbc2, comps,
                    lambda idx: _fit_logit(X[idx], y[idx]), 2, n_boot, seed)
                gap_early = 100 * (_pooled_rate(pp, tool, ep, opt, *EARLY)
                                   - _pooled_rate(pp, tool, ep, "raw", *EARLY))
                gap_late = 100 * (_pooled_rate(pp, tool, ep, opt, *LATE)
                                  - _pooled_rate(pp, tool, ep, "raw", *LATE))
                rows.append({
                    "tool": tool, "endpoint": ep, "contrast": f"{opt}_vs_raw",
                    "raw_rank_slope_logodds": raw_slope, "raw_slope_lo": rs_lo,
                    "raw_slope_hi": rs_hi, "raw_slope_p": rs_p,
                    "interaction_logodds_per_rank": b3, "inter_lo": lo, "inter_hi": hi,
                    "interaction_p": bp,
                    "gap_early_pp": gap_early, "gap_late_pp": gap_late,
                    "gap_change_late_minus_early_pp": gap_late - gap_early})
    out = pd.DataFrame(rows)
    if not out.empty:
        out["interaction_p_bh"] = bh_fdr(out["interaction_p"].to_numpy())
        out["star"] = [p_stars(p) for p in out["interaction_p_bh"]]
        out["raw_slope_star"] = [p_stars(p) for p in out["raw_slope_p"]]
    return out


# ─────────────────────────────────────────────────────────────────────────────
def _print_a(a: pd.DataFrame) -> None:
    print("\n" + "=" * 92)
    print("(a) DOES OPTIMIZATION HELP?  Cochran's Q + pairwise exact McNemar (Holm), paired 0/1")
    print("=" * 92)
    for tool, _k, _t in TOOLS:
        for ep, lab in ENDPOINTS:
            sub = a[(a["tool"] == tool) & (a["endpoint"] == ep)].sort_values("rank")
            if sub.empty:
                continue
            print(f"\n{tool} · {lab}")
            for _, r in sub.iterrows():
                tag = "  [PRIMARY]" if r["primary"] else ""
                print(f"  rank {int(r['rank']):>2} (n={int(r['n'])}){tag}:  "
                      f"raw {100*r['rate_raw']:5.1f}%  smina {100*r['rate_smina']:5.1f}%  "
                      f"gnina {100*r['rate_gnina']:5.1f}%   "
                      f"Cochran Q={r['cochran_Q']:.1f} p={fmt_p(r['cochran_p'])} "
                      f"(BH {fmt_p(r['cochran_p_bh'])})")
                for opt in CONTRASTS:
                    print(f"        raw->{opt:5s}: {r[f'{opt}_minus_raw_pp']:+5.1f} pp  "
                          f"(flips {int(r[f'{opt}_opt_gain'])} gain / {int(r[f'{opt}_raw_gain'])} lose)  "
                          f"McNemar p={fmt_p(r[f'{opt}_mcnemar_p'])} "
                          f"Holm={fmt_p(r[f'{opt}_mcnemar_holm'])} {r[f'{opt}_star']}")


def _print_b(b: pd.DataFrame) -> None:
    print("\n" + "=" * 92)
    print("(b) DOES THE BENEFIT DEPEND ON RANK?  logistic opt x rank, cluster(complex) bootstrap")
    print("=" * 92)
    for tool, _k, _t in TOOLS:
        for ep, lab in ENDPOINTS:
            sub = b[(b["tool"] == tool) & (b["endpoint"] == ep)]
            if sub.empty:
                continue
            r0 = sub.iloc[0]
            print(f"\n{tool} · {lab}")
            ctl = " (NEGATIVE CONTROL — generation order)" if tool == "EquiBind" else ""
            print(f"  raw rank slope = {r0['raw_rank_slope_logodds']:+.4f} log-odds/rank "
                  f"[{r0['raw_slope_lo']:+.4f}, {r0['raw_slope_hi']:+.4f}] "
                  f"p={fmt_p(r0['raw_slope_p'])} {r0['raw_slope_star']}{ctl}")
            for _, r in sub.iterrows():
                print(f"  {r['contrast']:11s}: opt x rank b3 = "
                      f"{r['interaction_logodds_per_rank']:+.4f}/rank "
                      f"[{r['inter_lo']:+.4f}, {r['inter_hi']:+.4f}] "
                      f"p={fmt_p(r['interaction_p'])} (BH {fmt_p(r['interaction_p_bh'])}) {r['star']}"
                      f"   |  benefit early(1-10)={r['gap_early_pp']:+.1f}pp "
                      f"late(21-30)={r['gap_late_pp']:+.1f}pp "
                      f"Δ={r['gap_change_late_minus_early_pp']:+.1f}pp")


def regenerate_annotated_figures(df: pd.DataFrame, a: pd.DataFrame, b: pd.DataFrame,
                                 out_dir: Path) -> None:
    """Overwrite 15c / 15d with the stats annotated on the panels. Reuses the pipeline's
    own plot functions (imported lazily) so the plotting stays single-sourced; only the
    annotation layer is added, driven by the ``stats={"a":..., "b":...}`` payload."""
    try:
        from posebusters_pose_comparison import (aggregate_optimization_benefit,
                                                 plot_optimization_benefit_by_group,
                                                 plot_optimization_benefit_by_rank)
    except Exception as exc:                          # pragma: no cover
        print(f"  (skipped figure annotation — could not import pipeline: {exc})")
        return
    rec = aggregate_optimization_benefit(df, max_rank=MAX_RANK)
    if rec is None or rec.empty:
        print("  (skipped figure annotation — no optimizer variants to plot)")
        return
    stats = {"a": a, "b": b}
    plot_optimization_benefit_by_rank(
        rec, out_dir / "15c_optimization_benefit_by_rank.png", stats=stats)
    plot_optimization_benefit_by_group(
        rec, out_dir / "15d_optimization_benefit_by_group.png", stats=stats)
    print(f"  annotated figures -> 15c / 15d_optimization_benefit_by_rank/group.png")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Whole-protein (blind) campaign — the run the Results chapter reports. The old
    # crystal-boxed dir (posebusters_results/benchmark/dock/pose_comparison_report) was
    # the default here and is a different campaign entirely, so a bare invocation used to
    # annotate figures 15c/15d that no longer belong to the reported numbers.
    default_dir = Path("posebusters_results/benchmark_full_protein_vina_scoring"
                       "/dock/pose_comparison_report")
    ap.add_argument("--report-dir", type=Path, default=default_dir,
                    help="dir holding per_pose_metrics.csv (default: %(default)s)")
    ap.add_argument("--out-dir", type=Path, default=None,
                    help="where to write results (default: --report-dir)")
    ap.add_argument("--n-boot", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--no-annotate", dest="annotate", action="store_false",
                    help="skip overwriting figs 15c/15d with the annotated versions")
    ap.add_argument("--annotate-only", action="store_true",
                    help="re-draw 15c/15d from the existing *_stats CSVs; skip the "
                         "(slow) recomputation — use to restore annotations a plain "
                         "report run wiped")
    mf.add_method_filter_args(ap)
    args = ap.parse_args()
    out_dir = args.out_dir or args.report_dir
    src = args.report_dir / "per_pose_metrics.csv"
    if not src.exists():
        ap.error(f"missing {src}")

    df = pd.read_csv(src, low_memory=False)
    # Feeds BOTH the --annotate-only redraw and the main path, so one call here
    # covers every consumer.
    df = mf.apply_method_filter(df, "method", args, label="optimization-benefit",
                                out_dir=out_dir)

    if args.annotate_only:                              # reuse prior results, just re-draw
        ap_csv = out_dir / "optimization_benefit_stats_pairwise.csv"
        bp_csv = out_dir / "optimization_benefit_stats_interaction.csv"
        if not (ap_csv.exists() and bp_csv.exists()):
            ap.error("--annotate-only needs existing *_stats CSVs; run without it first")
        a = pd.read_csv(ap_csv)
        b = pd.read_csv(bp_csv)
        if "primary" in a.columns:
            a["primary"] = a["primary"].astype(bool)
        regenerate_annotated_figures(df, a, b, out_dir)
        print("re-annotated 15c/15d from existing stats CSVs (no recompute)")
        return

    pp = build_per_pose(df)
    if pp.empty:
        ap.error("no optimizer-variant poses found in per_pose_metrics.csv")

    a = question_a(pp)
    b = question_b(pp, n_boot=args.n_boot, seed=args.seed)
    _print_a(a)
    _print_b(b)

    a.to_csv(out_dir / "optimization_benefit_stats_pairwise.csv", index=False)
    b.to_csv(out_dir / "optimization_benefit_stats_interaction.csv", index=False)
    sidecar = {
        "source": str(src), "n_boot": args.n_boot, "seed": args.seed,
        "primary_rank": PRIMARY_RANK, "support_ranks": SUPPORT_RANKS,
        "question_a_pairwise": a.to_dict(orient="records"),
        "question_b_interaction": b.to_dict(orient="records"),
    }
    with open(out_dir / "optimization_benefit_stats.json", "w") as fh:
        json.dump(sidecar, fh, indent=2, default=float)
    print(f"\nwrote optimization_benefit_stats_pairwise.csv / _interaction.csv / .json "
          f"-> {out_dir}")
    if args.annotate:
        regenerate_annotated_figures(df, a, b, out_dir)


if __name__ == "__main__":
    main()
