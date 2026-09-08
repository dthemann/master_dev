#!/usr/bin/env python
"""Recompute every appendix statistic the nearest-copy switch needs that no
rebuilt sidecar carries, from the rebuilt hub table
posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest/per_pose_metrics.csv
with Scripts/Analysis/stats_utils.py, exactly as thesis_assertions.py does for
the canonical values.  Writes computed_appendix_values.json beside this file.

Read-only on every input.  Nothing under thesis_latex/ or a canonical result
directory is touched.
"""
from __future__ import annotations
import json, math, sys, itertools
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path("/home/manndo/master_dev")
sys.path.insert(0, str(ROOT / "Scripts/Analysis"))
from stats_utils import (wilson_ci, newcombe_paired_diff_ci, mcnemar_exact,
                         mcnemar_power, tost_paired_proportions, cochran_q, holm)
from scipy.stats import norm

HUB = ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest"
CANON_HUB = ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report"
OUT = Path(__file__).resolve().parent / "computed_appendix_values.json"

df = pd.read_csv(HUB / "per_pose_metrics.csv", low_memory=False)
df["cid"] = df.protein.astype(str) + "_" + df.ligand.astype(str)
df["pb_valid_b"] = df.pb_valid.astype("boolean").fillna(False).astype(bool)
CIDS = sorted(df.cid.unique())
N = len(CIDS)
assert N == 303, N
R = {}

# ----------------------------------------------------------------------------- helpers
def eff_rank(x: pd.DataFrame, rankcol: str) -> pd.Series:
    """Effective rank: the named column, or for EquiBind (999 sentinel) the
    ascending order of its own refiner affinity within each complex."""
    if rankcol in ("gnina_affinity", "smina_affinity"):
        aff = pd.to_numeric(x[rankcol], errors="coerce")
        return aff.groupby(x["cid"]).rank(method="first", ascending=True)
    if rankcol in ("cnn_score_desc",):
        v = -pd.to_numeric(x["cnn_score"], errors="coerce")
        return v.groupby(x["cid"]).rank(method="first", ascending=True)
    if rankcol in ("minimized_affinity",):
        v = pd.to_numeric(x["minimized_affinity"], errors="coerce")
        return v.groupby(x["cid"]).rank(method="first", ascending=True)
    r = pd.to_numeric(x[rankcol], errors="coerce")
    if r.notna().sum() == 0 or (r.dropna() == 999).all():
        aff = pd.to_numeric(x["gnina_affinity"], errors="coerce")
        return aff.groupby(x["cid"]).rank(method="first", ascending=True)
    return r

def gate_vec(meth, depth, rankcol, gate="double", rmsd_col="rmsd",
             bestfit_col="bestfit_rmsd"):
    """Existence gate over the first `depth` ranks -> 0/1 numpy over CIDS."""
    x = df[df.method == meth].copy()
    x["_rk"] = eff_rank(x, rankcol)
    x = x[x._rk <= depth]
    near = x[rmsd_col] <= 2.0
    if gate == "near":
        ok = near
    elif gate == "double":
        ok = near & x.pb_valid_b
    elif gate == "triple":
        ok = near & x.pb_valid_b & (x[bestfit_col] <= 1.0)
    elif gate == "form":
        ok = (x[bestfit_col] <= 1.0) & (x[rmsd_col] < 1000)
    elif gate == "pb_double":            # PoseBusters' own in-place RMSD gate
        ok = (x["pb_rmsd"] <= 2.0) & x.pb_valid_b
    else:
        raise ValueError(gate)
    s = x.assign(ok=ok).groupby("cid").ok.any()
    return s.reindex(CIDS, fill_value=False).astype(int).values

def count(v): return int(v.sum())

def mcn(a, b):
    n10, n01, p = mcnemar_exact(a, b)
    return {"a_only": n10, "b_only": n01, "p": p}

def fmt(x, nd=1): return None if x is None else round(float(x), nd)

# ----------------------------------------------------------------------------- 1. ladder (App. C :173-248, Table 8)
LADDER = [("autodock_mgltools_exh18", "rank", "Exhaustiveness 18"),
          ("autodock_mgltools", "rank", "Exhaustiveness 32"),
          ("autodock_mgltools_gnina", "optimized_rank", "32 + gnina"),
          ("autodock_mgltools_exh64", "rank", "Exhaustiveness 64"),
          ("autodock_mgltools_exh64_gnina", "optimized_rank", "64 + gnina"),
          ("autodock_mgltools_exh92", "rank", "Exhaustiveness 92"),
          ("autodock_mgltools_exh128", "rank", "Exhaustiveness 128"),
          ("autodock_mgltools_exh128_gnina", "optimized_rank", "128 + gnina, selected")]
SEL = "autodock_mgltools_exh128_gnina"
lad = {}
for gate in ("double", "triple", "near"):
    lad[gate] = {}
    for meth, rc, label in LADDER:
        lad[gate][meth] = {d: gate_vec(meth, d, rc, gate) for d in range(1, 31)}

def ladder_contrasts(gate, depth, family="seven"):
    sel = lad[gate][SEL][depth]
    rows = []
    for meth, rc, label in LADDER:
        if meth == SEL: continue
        v = lad[gate][meth][depth]
        m = mcn(sel, v)
        rows.append({"arm": meth, "label": label, "k_arm": count(v), "k_selected": count(sel),
                     "selected_only": m["a_only"], "arm_only": m["b_only"], "p_raw": m["p"]})
    ps = holm([r["p_raw"] for r in rows])
    for r, ph in zip(rows, ps): r["p_holm_seven"] = float(ph)
    return rows

def all_pairs_holm(gate, depth):
    keys = [m for m, _, _ in LADDER]
    pairs, ps = [], []
    for a, b in itertools.combinations(keys, 2):
        m = mcn(lad[gate][a][depth], lad[gate][b][depth])
        pairs.append((a, b)); ps.append(m["p"])
    ph = holm(ps)
    return {f"{a}|{b}": {"p_raw": float(p), "p_holm_28": float(h)}
            for (a, b), p, h in zip(pairs, ps, ph)}

R["table_8"] = {meth: {"label": label, "rank_col": rc,
                       "double": {d: count(lad["double"][meth][d]) for d in (1, 5, 15, 30)},
                       "triple": {d: count(lad["triple"][meth][d]) for d in (1, 5, 15, 30)},
                       "near":   {d: count(lad["near"][meth][d]) for d in (1, 5, 15, 30)}}
                for meth, rc, label in LADDER}
R["table_8_holm_top15_double"] = ladder_contrasts("double", 15)
R["ladder_rank1_double_seven"] = ladder_contrasts("double", 1)
R["ladder_top15_triple_seven"] = ladder_contrasts("triple", 15)
R["ladder_rank1_triple_seven"] = ladder_contrasts("triple", 1)
R["ladder_top30_double_raw128_vs_selected"] = mcn(lad["double"]["autodock_mgltools_exh128"][30], lad["double"][SEL][30]) | {
    "raw128": count(lad["double"]["autodock_mgltools_exh128"][30]), "selected": count(lad["double"][SEL][30])}
p28d = all_pairs_holm("double", 15); p28t = all_pairs_holm("triple", 15)
R["ladder_28pair_holm_top15"] = {
    "double": {k: v for k, v in p28d.items() if SEL in k},
    "triple": {k: v for k, v in p28t.items() if SEL in k}}
# depth lead of the selected arm on the double gate
lead = []
for d in range(1, 31):
    s = count(lad["double"][SEL][d]); others = {m: count(lad["double"][m][d]) for m, _, _ in LADDER if m != SEL}
    mx = max(others.values())
    lead.append({"depth": d, "selected": s, "best_other": mx,
                 "best_other_arm": [m for m, v in others.items() if v == mx],
                 "state": "leads" if s > mx else ("ties" if s == mx else "trails")})
R["ladder_depth_lead_double"] = lead
R["ladder_depth_lead_summary"] = {
    "strict_lead_depths": [r["depth"] for r in lead if r["state"] == "leads"],
    "tie_depths": [r["depth"] for r in lead if r["state"] == "ties"],
    "trail_depths": [r["depth"] for r in lead if r["state"] == "trails"]}
# :184 rescoring at 128 gain/loss top-15 double; 32->128 rescored top-15
m = mcn(lad["double"]["autodock_mgltools_exh128"][15], lad["double"][SEL][15])
R["line184"] = {"raw128_top15": count(lad["double"]["autodock_mgltools_exh128"][15]),
                "sel_top15": count(lad["double"][SEL][15]),
                "gained_by_rescoring": m["b_only"], "lost_by_rescoring": m["a_only"],
                "gnina32_top15": count(lad["double"]["autodock_mgltools_gnina"][15])}
# parent-Vina-order carry of the rescored poses
x = df[df.method == SEL]
R["line184_parent_vina_order"] = {
    "gnina_poses_by_autodock_rank_d1": count(gate_vec(SEL, 1, "autodock_rank", "double")),
    "gnina_poses_by_autodock_rank_d15": count(gate_vec(SEL, 15, "autodock_rank", "double")),
    "raw128_d1": count(lad["double"]["autodock_mgltools_exh128"][1]),
    "raw128_d15": count(lad["double"]["autodock_mgltools_exh128"][15]),
    "note": "autodock_rank is the parent Vina rank carried on every rescored pose; "
            "'rank'=='optimized_rank' on this variant (CNNaffinity order)"}
# :173 three ranking heads on the rescored exh128 pool (double gate)
heads = {"cnn_affinity (optimized_rank)": "optimized_rank",
         "cnn_score": "cnn_score_desc", "minimized_vina_energy": "minimized_affinity"}
hv = {k: {d: gate_vec(SEL, d, rc, "double") for d in (1, 15, 30)} for k, rc in heads.items()}
mm = mcn(hv["cnn_affinity (optimized_rank)"][1], hv["cnn_score"][1])
R["line173_heads"] = {k: {d: count(v[d]) for d in (1, 15, 30)} for k, v in hv.items()} | {
    "rank1_discordant_affinity_only": mm["a_only"], "rank1_discordant_score_only": mm["b_only"],
    "rank1_mcnemar_p": mm["p"],
    "rank1_pct": {k: fmt(100 * count(v[1]) / N) for k, v in hv.items()},
    "same_set_top30": bool(all((hv[k][30] == hv["cnn_affinity (optimized_rank)"][30]).all() for k in hv))}
# :246 plateau, near-native without validity at rank-1
R["line246"] = {
    "raw_rank1_near": {m: count(lad["near"][m][1]) for m in ("autodock_mgltools_exh64", "autodock_mgltools_exh92", "autodock_mgltools_exh128")},
    "rescored_rank1_near": {m: count(lad["near"][m][1]) for m in ("autodock_mgltools_exh64_gnina", SEL)},
    "raw_rank1_near_all": {m: count(lad["near"][m][1]) for m, _, _ in LADDER}}
# :195 / :235 exh92->exh128 and exh64->exh92 at rank-1
for gate in ("near", "double", "triple"):
    for a, b in (("autodock_mgltools_exh92", "autodock_mgltools_exh128"), ("autodock_mgltools_exh64", "autodock_mgltools_exh92")):
        mmm = mcn(lad[gate][a][1], lad[gate][b][1])
        R.setdefault("ladder_steps_rank1", {})[f"{gate}:{a}->{b}"] = mmm | {"k_a": count(lad[gate][a][1]), "k_b": count(lad[gate][b][1])}
# first depth at which exh92->exh128 is individually resolvable (p<0.05), near and triple gates
for gate in ("near", "triple"):
    seq = []
    for d in range(1, 31):
        mmm = mcn(lad[gate]["autodock_mgltools_exh92"][d], lad[gate]["autodock_mgltools_exh128"][d])
        seq.append({"depth": d, "delta": count(lad[gate]["autodock_mgltools_exh128"][d]) - count(lad[gate]["autodock_mgltools_exh92"][d]), "p": mmm["p"]})
    R.setdefault("step_92_128_by_depth", {})[gate] = seq
# :244 step exh92->128 pool and rank-1, both gates
R["line244"] = {g: {"pool": count(lad[g]["autodock_mgltools_exh128"][30]) - count(lad[g]["autodock_mgltools_exh92"][30]),
                    "rank1": count(lad[g]["autodock_mgltools_exh128"][1]) - count(lad[g]["autodock_mgltools_exh92"][1])}
                for g in ("double", "triple", "near")}
# :248 winner's curse: standard estimator over the eight-arm family, both depths.
# Estimator (as used for the per-family figure in the Limitations): the expected
# maximum of eight independent binomial rates minus the true best, approximated by
# the normal order statistic E[max] = mu + sigma * e8 with e8 = 1.4236 (expected
# maximum of 8 standard normals) applied at the selected arm's rate.
e8 = 1.4236
for d in (1, 15):
    p = count(lad["double"][SEL][d]) / N
    R.setdefault("line248_winners_curse", {})[f"d{d}"] = {
        "selected_rate_pct": fmt(100 * p, 2),
        "sigma_pp": fmt(100 * math.sqrt(p * (1 - p) / N), 2),
        "bias_pp_normal_order_stat_8_arms": fmt(100 * e8 * math.sqrt(p * (1 - p) / N), 2),
        "note": "NO registered generator found for the eight-arm figure printed at :248 "
                "(grep for winner/curse over Scripts/Analysis finds none); this is a "
                "normal-order-statistic approximation, flagged UNRESOLVED"}

# ----------------------------------------------------------------------------- 2. EquiBind nine configurations (:267-269)
EQ = {"equibind_unguided_raw": "rank", "equibind_unguided_smina": "smina_affinity", "equibind_unguided_gnina": "gnina_affinity",
      "equibind_fpocket_raw": "rank", "equibind_fpocket_smina": "smina_affinity", "equibind_fpocket_gnina": "gnina_affinity",
      "equibind_p2rank_raw": "rank", "equibind_p2rank_smina": "smina_affinity", "equibind_p2rank_gnina": "gnina_affinity"}
eq = {}
for meth, rc in EQ.items():
    x = df[df.method == meth]
    eq[meth] = {"triple_top15": count(gate_vec(meth, 15, rc, "triple")),
                "double_top15": count(gate_vec(meth, 15, rc, "double")),
                "near_pool": count(gate_vec(meth, 30, rc, "near")),
                "near_pool_poses": int((x.rmsd <= 2).sum()),
                "form_pool": count(gate_vec(meth, 30, rc, "form")),
                "pb_valid_pool_complexes": int(x[x.pb_valid_b].cid.nunique()),
                "pb_valid_pct": fmt(100 * x.pb_valid_b.mean()),
                "n_poses": int(len(x))}
R["line267_equibind"] = eq

# ----------------------------------------------------------------------------- 3. Wilson table (:893-920), Table 22 (:1056)
tk = pd.read_csv(HUB / "topk_recovery_validity.csv")
wil = {}
for (var, mk), g in tk.groupby(["variant", "method_key"]):
    wil[var] = {}
    for k in (1, 5, 10, 15):
        r = g[g.k == k].iloc[0]
        lo, hi = wilson_ci(int(r.near_k), 303); lo2, hi2 = wilson_ci(int(r.valid_k), 303)
        wil[var][k] = {"near_k": int(r.near_k), "near_pct": fmt(100 * r.near_k / 303), "near_ci": [fmt(100 * lo), fmt(100 * hi)],
                       "valid_k": int(r.valid_k), "valid_pct": fmt(100 * r.valid_k / 303), "valid_ci": [fmt(100 * lo2), fmt(100 * hi2)],
                       "gap_k": int(r.near_k - r.valid_k), "gap_pct": fmt(100 * (r.near_k - r.valid_k) / 303),
                       "gap_ci": [fmt(100 * v) for v in wilson_ci(int(r.near_k - r.valid_k), 303)]}
R["wilson_table"] = wil

# ----------------------------------------------------------------------------- 4. Table 19 footnote: DiffDock re-rank (gnina) and gains
# gnina re-rank of unrefined DiffDock poses: pick per complex the diffdock_gnina pose
# with the lowest gnina minimized affinity (its 'gnina_affinity'/'minimized_affinity'),
# read the raw pose at the same confidence rank.
raw = df[df.method == "diffdock"][["cid", "rank", "rmsd"]]
def rerank_rate(refined_key, affcol):
    g = df[df.method == refined_key][["cid", "rank", affcol]].dropna(subset=[affcol])
    top = g.sort_values(affcol).groupby("cid").head(1)[["cid", "rank"]]
    j = top.merge(raw, on=["cid", "rank"], how="left")
    hit = (j.rmsd <= 2).groupby(j.cid).any().reindex(CIDS, fill_value=False)
    return {"k": int(hit.sum()), "pct": fmt(100 * hit.sum() / N), "n_joined": int(j.rmsd.notna().sum())}
R["table19_footnote_rerank"] = {
    "gnina_rerank_raw_coords_rank1": rerank_rate("diffdock_gnina", "gnina_affinity" if df[df.method == "diffdock_gnina"].gnina_affinity.notna().any() else "minimized_affinity"),
    "smina_rerank_raw_coords_rank1": rerank_rate("diffdock_smina", "smina_affinity" if df[df.method == "diffdock_smina"].smina_affinity.notna().any() else "minimized_affinity"),
    "confidence_order_rank1": {"k": count(gate_vec("diffdock", 1, "rank", "near")), "pct": fmt(100 * count(gate_vec("diffdock", 1, "rank", "near")) / N)},
    "note": "cross-check against PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv C_rerank (41.58 %)"}

# ----------------------------------------------------------------------------- 5. Table 22 (:1056-1062) from validity_gate_cost.csv
vgc = pd.read_csv(HUB / "validity_gate_cost.csv")
t22 = {}
for arm, g in vgc.groupby("arm"):
    t22[arm] = {int(r.k): {"cost_complexes": int(r.cost_complexes), "pct": fmt(100 * r.cost_complexes / 303),
                           "ci": [fmt(100 * v) for v in wilson_ci(int(r.cost_complexes), 303)]} for _, r in g.iterrows()}
R["table_22"] = t22

# ----------------------------------------------------------------------------- 6. Table 24 near-native half (:1172-1181) + :1188 counts
import yaml
spec = yaml.safe_load(open(Path(__file__).resolve().parent / "thesis_expected_values_nearest.yaml"))["table_24"]
checks = [c for grp in spec["groups"].values() for c in grp]
pb = pd.read_csv(ROOT / spec["pb_input"], low_memory=False, usecols=checks + ["pose_file"])
def tb(col): return col.astype("boolean").fillna(False).to_numpy()
pb["_valid"] = np.logical_and.reduce([tb(pb[c]) for c in checks])
j = df[["method", "pose_file", "rmsd", "rmsd_ref_instance"]].merge(pb, on="pose_file", how="left")
assert j["_valid"].isna().sum() == 0
def block(x):
    share = lambda cols: (~np.logical_and.reduce([tb(x[c]) for c in cols]))
    inter = share(spec["groups"]["intermolecular"])
    return {"poses": int(len(x)), "valid_pct": fmt(100 * tb(x["_valid"]).mean()),
            "chem_pct": fmt(100 * share(spec["groups"]["chemical"]).mean()),
            "intra_pct": fmt(100 * share(spec["groups"]["intramolecular"]).mean()),
            "inter_pct": fmt(100 * inter.mean()), "inter_fail_n": int(inter.sum()),
            "mindist_pct": fmt(100 * (~tb(x[spec["single_check"]])).mean()),
            "mindist_fail_n": int((~tb(x[spec["single_check"]])).sum())}
t24 = {}
for arm, meth in spec["arms"].items():
    d = j[j.method == meth]
    t24[arm] = {"all": block(d), "near_native_nearest": block(d[d.rmsd <= 2.0]),
                "near_native_instance": block(d[d.rmsd_ref_instance <= 2.0])}
R["table_24"] = t24

# ----------------------------------------------------------------------------- 7. H.9 (:1258-1260) recompute
A = {d: gate_vec(SEL, d, "optimized_rank", "double") for d in (1, 5, 15, 30)}
B = {d: gate_vec("diffdock_smina", d, "rank", "double") for d in (1, 5, 15, 30)}
z90 = norm.ppf(0.95)
h9 = {}
for d in (1, 5, 15, 30):
    a, b = A[d], B[d]
    ci = newcombe_paired_diff_ci(a, b); c90 = newcombe_paired_diff_ci(a, b, z=z90)
    _, _, p = mcnemar_exact(a, b)
    lo90, hi90 = -100 * c90["hi"], -100 * c90["lo"]
    row = {"autodock": count(a), "diffdock": count(b), "diff_pp": fmt(-100 * ci["diff"], 2),
           "ci95": [fmt(-100 * ci["hi"], 2), fmt(-100 * ci["lo"], 2)], "mcnemar_p": float(p),
           "ci90": [fmt(lo90, 2), fmt(hi90, 2)],
           "smallest_passing_margin_pp": math.ceil(max(abs(lo90), abs(hi90)) * 10) / 10,
           "both": ci["n11"], "autodock_only": ci["n10"], "diffdock_only": ci["n01"], "neither": ci["n00"],
           "discordant": ci["n_discordant"], "discordant_pct": fmt(100 * ci["n_discordant"] / N),
           "equiv": {m: bool(tost_paired_proportions(a, b, margin=m / 100)["equivalent"]) for m in (5, 10, 12, 15)}}
    pw = mcnemar_power(a, b)
    row |= {"mde_pp": fmt(100 * pw["mde"]), "observed_power": fmt(pw["observed_power"], 3),
            "n_needed": fmt(pw["n_needed"], 0), "n_needed_rounded_100": int(round(pw["n_needed"], -2)),
            "mde_at_428_pp": fmt(100 * pw["mde"] * math.sqrt(N / 428.0))}
    h9[d] = row
R["appendix_h9"] = h9

# ----------------------------------------------------------------------------- 8. Table 27 footnote (:1387)
T27 = {"AutoDock Vina + gnina": (SEL, "optimized_rank"), "DiffDock + smina": ("diffdock_smina", "rank"),
       "EquiBind + gnina": ("equibind_unguided_gnina", "gnina_affinity")}
fn = {}
for lab, (meth, rc) in T27.items():
    x = df[df.method == meth].copy(); x["_rk"] = eff_rank(x, rc)
    r1 = x[(x._rk == 1)]
    coh = r1[r1.pb_valid_b & (r1.centroid_dist <= 8.0)]
    own = gate_vec(meth, 1, rc, "double"); pbg = gate_vec(meth, 1, rc, "pb_double")
    mmm = mcn(own, pbg)
    fn[lab] = {"near_site_rank1_n": int(len(coh)), "median_pb_rmsd": fmt(coh.pb_rmsd.median(), 3), "median_rmsd": fmt(coh.rmsd.median(), 3),
               "rank1_recovery_own_gate": count(own), "rank1_recovery_pb_rmsd_gate": count(pbg),
               "rank1_recovery_pb_gate_pct": fmt(100 * count(pbg) / N), "decisions_moved_missed_to_recovered": mmm["b_only"],
               "decisions_moved_recovered_to_missed": mmm["a_only"]}
R["table27_footnote"] = fn | {"total_decisions": 3 * N,
                              "rank1_pbrmsd_gt_rmsd_rows": int((df.pb_rmsd > df.rmsd + 1e-3).sum())}

# ----------------------------------------------------------------------------- 9. Table 1 Form-loss rows (:666, B2)
T1 = {"AutoDock (raw)": "autodock_mgltools_exh128", "AutoDock + gnina": SEL, "DiffDock (raw)": "diffdock",
      "DiffDock + smina": "diffdock_smina", "DiffDock + gnina": "diffdock_gnina", "EquiBind (raw)": "equibind_unguided_raw",
      "EquiBind + smina": "equibind_unguided_smina", "EquiBind + gnina": "equibind_unguided_gnina",
      "EquiBind (raw) fpocket": "equibind_fpocket_raw", "EquiBind (raw) P2Rank": "equibind_p2rank_raw"}
fl = {}
for lab, meth in T1.items():
    x = df[df.method == meth]
    fn_ = x[(x.bestfit_rmsd <= 1) & (x.rmsd < 1000)].cid.nunique()
    fi = x[(x.bestfit_rmsd_ref_instance <= 1) & (x.rmsd_ref_instance < 1000)].cid.nunique()
    fl[lab] = {"form_complexes_nearest": int(fn_), "form_complexes_instance": int(fi), "delta": int(fn_ - fi),
               "form_poses_nearest": int(((x.bestfit_rmsd <= 1) & (x.rmsd < 1000)).sum()),
               "form_poses_instance": int(((x.bestfit_rmsd_ref_instance <= 1) & (x.rmsd_ref_instance < 1000)).sum())}
R["table1_form_loss"] = fl | {"rows_losing_a_complex": [k for k, v in fl.items() if v["delta"] < 0]}

# ----------------------------------------------------------------------------- 10. bands table (:992-999)
def bands(path):
    b = pd.read_csv(path)
    out = {}
    for tool in ("DiffDock", "EquiBind"):
        for opt in ("smina", "gnina"):
            for crit in ("near_%", "pbv_%", "both_%"):
                for lo, hi in ((1, 10), (11, 20), (21, 30)):
                    def pooled(o):
                        g = b[(b.tool == tool) & (b.optimizer == o) & (b["rank"] >= lo) & (b["rank"] <= hi)]
                        return (g[crit] * g.n).sum() / g.n.sum()
                    out[f"{tool}|{crit}|{lo}-{hi}|{opt}"] = fmt(pooled(opt) - pooled("raw"), 2)
    return out
R["bands_nearest"] = bands(HUB / "optimization_benefit_by_rank.csv")
R["bands_canonical_check"] = bands(CANON_HUB / "optimization_benefit_by_rank.csv")

# ----------------------------------------------------------------------------- 11. :1207, :871
R["line1207"] = {"exh32_raw_top15": count(lad["double"]["autodock_mgltools"][15]), "exh128_raw_top15": count(lad["double"]["autodock_mgltools_exh128"][15])}
R["line871"] = {"diffdock_smina_rank1_valid_pct": wil["DiffDock + smina"][1]["valid_pct"], "diffdock_gnina_rank1_valid_pct": wil["DiffDock (gnina-opt)"][1]["valid_pct"],
                "diffdock_smina_rank1_near_pct": wil["DiffDock + smina"][1]["near_pct"], "diffdock_gnina_rank1_near_pct": wil["DiffDock (gnina-opt)"][1]["near_pct"]}

# ----------------------------------------------------------------------------- 12. multi-copy facts (:666, :722)
mc = df.groupby("cid").n_copies.first()
R["copy_facts"] = {"multi_copy_303": int((mc > 1).sum()), "single_copy_303": int((mc == 1).sum()),
                   "poses_nearest_alternate": int((df.nearest_copy_is_ref == False).sum())}

json.dump(R, open(OUT, "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("wrote", OUT)
