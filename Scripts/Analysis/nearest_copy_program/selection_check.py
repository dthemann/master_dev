#!/usr/bin/env python
"""Phase 4 selection re-check for the nearest-copy program (plan v2 D12, Phase 4).

Reads a hub report directory (default: the program's _nearest one) and recomputes, from its
per_pose_metrics.csv and its all-variants oracle sidecar, the three variant selections the thesis
relies on:
  * AutoDock rung: PB-valid & rmsd <= 2 & bestfit_rmsd <= 1 (triple gate) and the double gate at
    top-15 for every autodock_mgltools* arm, plus rank-1 / top-15 / top-30 double-gate recovery
  * DiffDock refiner: diffdock_smina vs diffdock_gnina, triple and double gate at top-15, and the
    hub's BEST_VARIANT_METRIC (oracle_pb_valid_and_rmsd2_%) from oracle_summary_all_variants.csv
  * EquiBind arm: the nine equibind_* keys on the same metric (rank-1 = lowest gnina / smina affinity,
    the rank column is a 999 sentinel)
Writes selection_check.csv and selection_check.txt into --out-dir (must be NEW or --overwrite) and
prints the verdict: whether autodock_mgltools_exh128_gnina, diffdock_smina and equibind_unguided_gnina
remain selected. Read-only on every input.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np, pandas as pd

REPO = Path("/home/manndo/master_dev")
DEFAULT_REPORT = REPO / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest"

def eff_rank(df: pd.DataFrame) -> pd.Series:
    r = pd.to_numeric(df["rank"], errors="coerce").astype(float)
    eq = df["method"].astype(str).str.startswith("equibind")
    if eq.any():
        e = df[eq].copy()
        g = pd.to_numeric(e.get("gnina_affinity"), errors="coerce")
        s = pd.to_numeric(e.get("smina_affinity"), errors="coerce")
        key = g.where(e["method"].str.endswith("_gnina"), s)
        e["_k"] = key
        order = e.sort_values(["method", "protein", "_k", "pose_file"], kind="mergesort")
        r.loc[order.index] = order.groupby(["method", "protein"]).cumcount().values + 1
    return r

def gates(df: pd.DataFrame, depth: int, kind: str) -> pd.Series:
    sub = df[df["eff_rank"] <= depth]
    valid = sub["pb_valid"].astype(str).str.lower().eq("true")
    near = pd.to_numeric(sub["rmsd"], errors="coerce") <= 2.0
    if kind == "double":
        ok = valid & near
    elif kind == "triple":
        ok = valid & near & (pd.to_numeric(sub["bestfit_rmsd"], errors="coerce") <= 1.0)
    elif kind == "near":
        ok = near
    else:
        raise ValueError(kind)
    return sub.assign(ok=ok).groupby("protein")["ok"].any()

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    if a.out_dir.exists() and any(a.out_dir.iterdir()) and not a.overwrite:
        sys.exit(f"REFUSING: {a.out_dir} exists and is not empty")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    ppm = a.report_dir / "per_pose_metrics.csv"
    df = pd.read_csv(ppm, low_memory=False)
    conv = df["reference_convention"].iloc[0] if "reference_convention" in df.columns else "instance (legacy table)"
    df["eff_rank"] = eff_rank(df)
    n_complex = df["protein"].nunique()
    rows = []
    fams = {
        "autodock": sorted(m for m in df["method"].unique() if str(m).startswith("autodock_mgltools")),
        "diffdock": sorted(m for m in df["method"].unique() if str(m).startswith("diffdock")),
        "equibind": sorted(m for m in df["method"].unique() if str(m).startswith("equibind")),
    }
    for fam, arms in fams.items():
        for m in arms:
            sub = df[df["method"] == m]
            rec = {"family": fam, "method": m, "n_complexes": int(sub["protein"].nunique())}
            for d in (1, 15, 30):
                rec[f"double_d{d}"] = int(gates(sub, d, "double").sum())
                rec[f"near_d{d}"] = int(gates(sub, d, "near").sum())
            rec["triple_d15"] = int(gates(sub, 15, "triple").sum())
            rec["triple_d30"] = int(gates(sub, 30, "triple").sum())
            rows.append(rec)
    out = pd.DataFrame(rows)
    # hub's own selection metric from the all-variants oracle sidecar, if present
    osum = a.report_dir / "oracle_summary_all_variants.csv"
    if osum.exists():
        o = pd.read_csv(osum)
        col = "oracle_pb_valid_and_rmsd2_%" if "oracle_pb_valid_and_rmsd2_%" in o.columns else None
        if col:
            o = o.rename(columns={col: "hub_oracle_pb_valid_and_rmsd2_pct"})[["method", "hub_oracle_pb_valid_and_rmsd2_pct"]]
            out = out.merge(o, on="method", how="left")
    out.to_csv(a.out_dir / "selection_check.csv", index=False)

    def winner(fam, col):
        """All arms attaining the family maximum (ties are reported, never broken silently)."""
        f = out[out["family"] == fam]
        if f.empty or col not in f or f[col].isna().all():
            return None, None
        vals = f.set_index("method")[col]
        best = vals.max()
        tops = sorted(vals[vals == best].index.tolist())
        return tops, vals.to_dict()

    def holds(tops, pinned):
        return tops is not None and pinned in tops
    lines = [f"Selection re-check on {ppm} ({conv}; {n_complex} complexes, {len(df):,} poses)", ""]
    verdicts = {}
    ad_w, ad_all = winner("autodock", "triple_d15")
    lines.append(f"AutoDock rung by top-15 TRIPLE gate: {ad_all}  -> {ad_w}")
    verdicts["autodock_triple_d15"] = holds(ad_w, "autodock_mgltools_exh128_gnina")
    ad_w2, _ = winner("autodock", "double_d15")
    lines.append(f"AutoDock rung by top-15 DOUBLE gate: -> {ad_w2}")
    dd_w, dd_all = winner("diffdock", "triple_d15")
    lines.append(f"DiffDock refiner by top-15 TRIPLE gate: {dd_all} -> {dd_w}")
    verdicts["diffdock_triple_d15"] = holds(dd_w, "diffdock_smina")
    if "hub_oracle_pb_valid_and_rmsd2_pct" in out.columns:
        dd_h, dd_hall = winner("diffdock", "hub_oracle_pb_valid_and_rmsd2_pct")
        lines.append(f"DiffDock refiner by hub BEST_VARIANT_METRIC (oracle_pb_valid_and_rmsd2_%): {dd_hall} -> {dd_h}")
        verdicts["diffdock_hub_metric"] = holds(dd_h, "diffdock_smina")
        eb_h, eb_hall = winner("equibind", "hub_oracle_pb_valid_and_rmsd2_pct")
        lines.append(f"EquiBind arm by hub BEST_VARIANT_METRIC: {eb_hall} -> {eb_h}")
        verdicts["equibind_hub_metric"] = holds(eb_h, "equibind_unguided_gnina")
    eb_w, eb_all = winner("equibind", "triple_d15")
    lines.append(f"EquiBind arm by top-15 TRIPLE gate: {eb_all} -> {eb_w}")
    verdicts["equibind_triple_d15"] = holds(eb_w, "equibind_unguided_gnina")
    lines.append("")
    lines.append("VERDICT: " + ("all three thesis selections HOLD (ties, if any, include the pinned arm and are listed above)" if all(verdicts.values()) else f"SELECTION CHANGES: {verdicts}"))
    (a.out_dir / "selection_check.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
