#!/usr/bin/env python
"""Per-family winner's-curse estimate for the variant selection (App. C, "One caveat applies to the
selected arm itself").

For each tool family the selected variant is the top-15 triple-gate winner (PB-valid AND rmsd <= 2 AND
bestfit_rmsd <= 1, existence over the top-15). The optimism of reporting the winner's recovery on the
complexes it was selected on is approximated by the standard winner's-curse formula

    bias ~= c_k * sigma_d / sqrt(2)

where k is the number of finalists in the family, c_k the expected maximum of k standard normals
(k = 2: 0.5642, k = 3: 0.8463, k = 4: 1.0294), and sigma_d the PAIRED standard error of the recovery
difference between the winner and its runner-up on the same complexes (sqrt(b + c - (b - c)^2 / n) / n,
with b and c the discordant counts), because the finalists are rescorings of the same poses. The same
estimate is given on the double gate (PB-valid AND rmsd <= 2) and at rank-1 for comparison. Read-only
on the per-pose table; the table's reference_convention column labels the block. Recipe recorded in the
project memory (2026-08-21) and reproduced on the canonical table before use.
"""
from __future__ import annotations
import argparse, json, math, re, sys
from pathlib import Path
import numpy as np, pandas as pd

C_K = {1: 0.0, 2: 0.5642, 3: 0.8463, 4: 1.0294, 5: 1.1630, 6: 1.2672, 7: 1.3522, 8: 1.4236, 9: 1.4850, 10: 1.5388}
FAMILIES = {
    "AutoDock (exh128)": ["autodock_mgltools_exh128", "autodock_mgltools_exh128_gnina"],
    "DiffDock": ["diffdock", "diffdock_smina", "diffdock_gnina"],
    "EquiBind (unguided)": ["equibind_unguided_raw", "equibind_unguided_smina", "equibind_unguided_gnina"],
}

def eff_rank(df):
    r = pd.to_numeric(df["rank"], errors="coerce").astype(float)
    eq = df["method"].astype(str).str.startswith("equibind")
    if eq.any():
        e = df[eq].copy()
        g = pd.to_numeric(e.get("gnina_affinity"), errors="coerce"); s = pd.to_numeric(e.get("smina_affinity"), errors="coerce")
        e["_k"] = g.where(e["method"].str.endswith("_gnina"), s)
        o = e.sort_values(["method", "protein", "_k", "pose_file"], kind="mergesort")
        r.loc[o.index] = o.groupby(["method", "protein"]).cumcount().values + 1
    return r

def recovered(df, arm, depth, gate, rmsd_col):
    sub = df[(df["method"] == arm) & (df["eff_rank"] <= depth)]
    valid = sub["pb_valid"].astype(str).str.lower().eq("true")
    near = pd.to_numeric(sub[rmsd_col], errors="coerce") <= 2.0
    ok = valid & near
    if gate == "triple":
        ok = ok & (pd.to_numeric(sub["bestfit_rmsd"], errors="coerce") <= 1.0)
    return sub.assign(ok=ok).groupby("protein")["ok"].any()

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--families", default=None, help="JSON dict name -> [method keys]; default the three thesis families")
    a = ap.parse_args()
    if a.out_dir.exists() and any(a.out_dir.iterdir()) and not a.overwrite:
        sys.exit(f"REFUSING: {a.out_dir} exists and is not empty")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    fams = json.loads(a.families) if a.families else dict(FAMILIES)
    df0 = pd.read_csv(a.per_pose_csv, usecols=["method"], low_memory=False)
    # the eight-rung ladder of App. C: exh18/32/64/92/128 raw plus gnina at 32/64/128; the base rung (exh32) is
    # keyed autodock_mgltools / autodock_mgltools_gnina in the per-pose table
    ladder = sorted(m for m in df0["method"].unique()
                    if re.fullmatch(r"autodock_mgltools(_exh\d+)?(_gnina)?", str(m)))
    if not a.families and len(ladder) >= 3:
        fams["AutoDock ladder (all exhaustiveness arms)"] = ladder   # the eight-configuration field of the Limitations
    df = pd.read_csv(a.per_pose_csv, low_memory=False)
    conv = df["reference_convention"].iloc[0] if "reference_convention" in df.columns else "instance"
    df["eff_rank"] = eff_rank(df)
    ids = sorted(df["protein"].unique()); n = len(ids)
    rows = []
    for fam, arms in fams.items():
        arms = [m for m in arms if m in set(df["method"])]
        for gate in ("triple", "double"):
            for depth in (1, 15):
                rec = {m: recovered(df, m, depth, gate, "rmsd").reindex(ids).fillna(False) for m in arms}
                counts = {m: int(v.sum()) for m, v in rec.items()}
                order = sorted(counts, key=lambda m: -counts[m])
                win, run = order[0], order[1] if len(order) > 1 else None
                if run is None:
                    continue
                b = int((rec[win] & ~rec[run]).sum()); c = int((~rec[win] & rec[run]).sum())
                sigma_d = math.sqrt(max(b + c - (b - c) ** 2 / n, 0.0)) / n * 100.0
                k = len(arms); ck = C_K.get(k, 1.0)
                rows.append(dict(convention=conv, family=fam, gate=gate, depth=depth, k=k, winner=win, runner_up=run,
                                 winner_count=counts[win], runner_count=counts[run], discordant_b=b, discordant_c=c,
                                 sigma_d_pp=round(sigma_d, 3), c_k=ck, bias_pp=round(ck * sigma_d / math.sqrt(2), 3),
                                 all_counts=json.dumps(counts)))
    out = pd.DataFrame(rows)
    out.to_csv(a.out_dir / "selection_bias.csv", index=False)
    head = out[(out.gate == "triple") & (out.depth == 15)]
    summ = {"convention": conv, "n_complexes": n, "formula": "c_k * sigma_d / sqrt(2), paired SE of winner vs runner-up",
            "top15_triple": {r.family: {"k": int(r.k), "winner": r.winner, "runner_up": r.runner_up, "bias_pp": r.bias_pp,
                                        "sigma_d_pp": r.sigma_d_pp, "discordant": [int(r.discordant_b), int(r.discordant_c)]}
                             for r in head.itertuples()}}
    (a.out_dir / "selection_bias_summary.json").write_text(json.dumps(summ, indent=2))
    print(out[["convention", "family", "gate", "depth", "k", "winner", "runner_up", "winner_count", "runner_count", "discordant_b", "discordant_c", "sigma_d_pp", "bias_pp"]].to_string(index=False))

if __name__ == "__main__":
    main()
