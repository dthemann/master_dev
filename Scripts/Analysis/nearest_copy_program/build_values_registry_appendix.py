#!/usr/bin/env python
"""Build values_registry_appendix.json: one entry per affected location of
thesis_latex_nearest/body_appendix_short.tex (plan v2 Section 3.3 and every
Part A / Part B item of NEAREST_COPY_TEXT_CHANGES_2026-09-08.md that concerns
the appendix).  Every `new` value comes from the program's rebuilt `_nearest`
outputs or from a recompute on the rebuilt per_pose_metrics.csv with
stats_utils (computed_appendix_values.json, written by
compute_appendix_values.py).  Nothing is taken from the probe previews.

Entry format (shared with the main-text registry):
  id          stable key, "A<line>.<slug>"
  file        body_appendix_short.tex
  line        line number in thesis_latex_nearest/body_appendix_short.tex (== shipped copy)
  anchor      label / section / table the line belongs to
  plan_ref    plan v2 Section 3.3 row and / or change-list item
  kind        value | wording | table_cell | table_block | decision | verbatim
  current     what the line prints today (verbatim fragment or number)
  new         rebuilt value (number, list or text); "W" for wording-only; None if unresolved
  source      where `new` comes from (file + selector, or "recompute: ..." recipe)
  status      resolved | unchanged | wording | unresolved | decision
  note        anything the editor must know (traps, rounding, decisions taken)
"""
from __future__ import annotations
import json, math, re
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path("/home/manndo/master_dev")
PROG = ROOT / "Scripts/Analysis/nearest_copy_program"
HUB = ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest"
CANON_EXH = ROOT / "posebusters_results/autodock_exhaustiveness_returns"
EXH = ROOT / "posebusters_results/autodock_exhaustiveness_returns_nearest"
MET = ROOT / "posebusters_results/metal_stratum_nearest"
ITT = ROOT / "posebusters_results/itt_nearest/itt_summary.json"
REF = ROOT / "posebusters_results/reference_convention_nearest"
CLU = ROOT / "posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes_nearest"
PMR = ROOT / "pandamap_results/benchmark_matched_equibind_nearest/report"
SIG = PROG / "harness_signatures/harness_signature_A_oldvalues_newtree.csv"
TEX = ROOT / "thesis_latex_nearest/body_appendix_short.tex"
FILE = "body_appendix_short.tex"

C = json.load(open(PROG / "computed_appendix_values.json"))
sig = pd.read_csv(SIG)
sig_map = dict(zip(sig.key, sig.actual))
sig_exp = dict(zip(sig.key, sig.expected))
tex_lines = TEX.read_text().splitlines()

def L(n): return tex_lines[n - 1]
def S(p): return str(p.relative_to(ROOT)) if str(p).startswith(str(ROOT)) else str(p)

entries = []
def E(id_, line, anchor, plan_ref, kind, current, new, source, status="resolved", note=""):
    entries.append({"id": id_, "file": FILE, "line": line, "anchor": anchor, "plan_ref": plan_ref,
                    "kind": kind, "current": current, "new": new, "source": source,
                    "status": status, "note": note})

def r1(x): return round(float(x), 1)
def r2(x): return round(float(x), 2)
def pfmt(p):
    p = float(p)
    if p < 1e-3: return f"{p:.1e}"
    return f"{p:.4f}".rstrip("0")

# ============================================================================ A:165 metal stratum (D4, W30, A5)
mc = pd.read_csv(MET / "metal_stratum_contrasts.csv")
ms = json.load(open(MET / "metal_stratum_summary.json"))
def rec(conv, cut, strat, depth, arm):
    r = mc[(mc.convention == conv) & (mc.cutoff_A == cut) & (mc.stratum == strat) & (mc.depth == depth) & (mc.kind == "recovery") & (mc.arm == arm)].iloc[0]
    return int(r.k), r1(r.rate_pct), int(r.n)
def con(conv, cut, strat, depth):
    r = mc[(mc.convention == conv) & (mc.cutoff_A == cut) & (mc.stratum == strat) & (mc.depth == depth) & (mc.kind == "contrast") & (mc.arm == "autodock_mgltools_exh128_gnina") & (mc.arm_b == "diffdock_smina")].iloc[0]
    return {"autodock": int(r.k_a), "diffdock": int(r.k_b), "n": int(r.n), "diffdock_minus_autodock_pp": r1(r.diff_pp),
            "autodock_only": int(r.n_a_only), "diffdock_only": int(r.n_b_only), "p_mcnemar": round(float(r.p_mcnemar), 4)}
def fisher(conv, cut, depth):
    r = mc[(mc.convention == conv) & (mc.cutoff_A == cut) & (mc.kind == "fisher_discordant") & (mc.depth == depth) & (mc.arm_b == "diffdock_smina")].iloc[0]
    return round(float(r.p_fisher), 3)
def nq(conv, cut, strat, arm):
    r = mc[(mc.convention == conv) & (mc.cutoff_A == cut) & (mc.stratum == strat) & (mc.kind == "no_qualifying_pose_among_rank1_failures") & (mc.arm == arm)].iloc[0]
    return int(r.k), int(r.n), r1(r.rate_pct)
AD, DD = "autodock_mgltools_exh128_gnina", "diffdock_smina"
metal = {conv: {
    "counts": {"4A": [ms["counts_metal_adjacent"]["4A"], ms["counts_metal_free"]["4A"]], "5A": [78, 225], "6A": [ms["counts_metal_adjacent"]["6A"], ms["counts_metal_free"]["6A"]]},
    "rank1_rate_pct": {"metal_free": {"autodock": rec(conv, 5.0, "metal_free", 1, AD)[1], "diffdock": rec(conv, 5.0, "metal_free", 1, DD)[1]},
                       "metal_adjacent": {"autodock": rec(conv, 5.0, "metal_adjacent", 1, AD)[1], "diffdock": rec(conv, 5.0, "metal_adjacent", 1, DD)[1]}},
    "rank1_k": {"metal_free": {"autodock": rec(conv, 5.0, "metal_free", 1, AD)[0], "diffdock": rec(conv, 5.0, "metal_free", 1, DD)[0]},
                "metal_adjacent": {"autodock": rec(conv, 5.0, "metal_adjacent", 1, AD)[0], "diffdock": rec(conv, 5.0, "metal_adjacent", 1, DD)[0]}},
    "top15_contrast_5A": {"metal_free": con(conv, 5.0, "metal_free", 15), "metal_adjacent": con(conv, 5.0, "metal_adjacent", 15)},
    "fisher_between_strata_top15_5A": fisher(conv, 5.0, 15),
    "top15_contrast_4A_pp": [con(conv, 4.0, "metal_free", 15)["diffdock_minus_autodock_pp"], con(conv, 4.0, "metal_adjacent", 15)["diffdock_minus_autodock_pp"]],
    "top15_contrast_6A_pp": [con(conv, 6.0, "metal_free", 15)["diffdock_minus_autodock_pp"], con(conv, 6.0, "metal_adjacent", 15)["diffdock_minus_autodock_pp"]],
    "top15_contrast_4A_full": {"metal_free": con(conv, 4.0, "metal_free", 15), "metal_adjacent": con(conv, 4.0, "metal_adjacent", 15)},
    "top15_contrast_6A_full": {"metal_free": con(conv, 6.0, "metal_free", 15), "metal_adjacent": con(conv, 6.0, "metal_adjacent", 15)},
    "rank1_failures_no_qualifying_pose_pct": {"metal_free": {"autodock": nq(conv, 5.0, "metal_free", AD), "diffdock": nq(conv, 5.0, "metal_free", DD)},
                                              "metal_adjacent": {"autodock": nq(conv, 5.0, "metal_adjacent", AD), "diffdock": nq(conv, 5.0, "metal_adjacent", DD)}},
} for conv in ("nearest", "instance")}
split = ms["cofactor_ion_split"]
src165 = S(MET / "metal_stratum_contrasts.csv") + " (convention == nearest, cutoff_A == 5.0 unless stated) and " + S(MET / "metal_stratum_summary.json")
E("A165.stratum_sizes", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30 / D4", "value",
  "78 of 303; 4 A rule 73; 6 A rule 81", {"5A": 78, "4A": metal["nearest"]["counts"]["4A"][0], "6A": metal["nearest"]["counts"]["6A"][0], "metal_free_5A": 225},
  S(MET / "metal_stratum_summary.json") + " counts_metal_adjacent", "unchanged",
  "Registered rule (D4): residue-level, reference instance, <= 5 A. Membership is copy-invariant by construction; 73 / 78 / 81 reproduce the print.")
E("A165.cofactor_ion_split", 165, "App. C, Metal-adjacent stratum", "3.3:165 / A5 / D4", "value",
  "Fifteen of the 78 ... cofactor ... remaining 63 ... free ion", {"cofactor": split["5A"]["n_cofactor"], "ion": split["5A"]["n_ion"], "cofactor_resnames": split["5A"]["cofactor_resnames"]},
  S(MET / "metal_stratum_summary.json") + " cofactor_ion_split.5A", "resolved",
  "DECISION per task brief (D4): print the registered rule's split, 16 cofactor / 62 free ion, and FLAG the change from 15 / 63 in the report. The extra cofactor is a second SF4 (iron-sulfur cluster) complex that the 4 A rule excludes (4 A split 15 / 58). Sentence becomes 'Sixteen of the 78 ... and the remaining 62 involve a free ion.'")
E("A165.rank1_rates", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30", "value",
  "37.8% AutoDock / 39.6% DiffDock metal-free (225); 33.3% / 17.9% metal-adjacent (78)",
  {"metal_free": metal["nearest"]["rank1_rate_pct"]["metal_free"], "metal_adjacent": metal["nearest"]["rank1_rate_pct"]["metal_adjacent"], "k": metal["nearest"]["rank1_k"]},
  src165 + " kind == recovery, depth == 1", "resolved",
  f"Instance block reproduces the print exactly ({metal['instance']['rank1_rate_pct']}). Under nearest AutoDock still leads on the metal-adjacent stratum (43.6 vs 26.9) and the two are level on the metal-free one (51.1 vs 49.3), so the sentence 'AutoDock leads on the metal-adjacent stratum while the two are level on the metal-free one' holds.")
E("A165.top15_contrasts", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30", "value",
  "-6.7 pp metal-free (134 against 149 of 225, p = 0.142); -20.5 pp metal-adjacent (33 against 49 of 78, p = 0.020)",
  metal["nearest"]["top15_contrast_5A"], src165 + " kind == contrast, depth == 15", "resolved",
  "Print as DiffDock minus AutoDock: -8.4 points (141 against 160 of 225, exact McNemar p = 0.056) and -19.2 points (38 against 53 of 78, p = 0.028). Instance block: " + json.dumps(metal["instance"]["top15_contrast_5A"]))
E("A165.fisher_between_strata", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30", "value", "Fisher's exact test on the discordant pairs p = 0.256",
  metal["nearest"]["fisher_between_strata_top15_5A"], src165 + " kind == fisher_discordant, depth == 15", "resolved",
  f"Instance block gives {metal['instance']['fisher_between_strata_top15_5A']} (print 0.256). Nearest 0.440. 'reported as heterogeneity rather than as an established interaction' still holds.")
E("A165.failure_shares", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30", "value",
  "metal-free 53.6% AutoDock / 65.4% DiffDock; metal-adjacent 55.8% / 70.3%",
  {"metal_free": {"autodock_pct": metal["nearest"]["rank1_failures_no_qualifying_pose_pct"]["metal_free"]["autodock"][2], "diffdock_pct": metal["nearest"]["rank1_failures_no_qualifying_pose_pct"]["metal_free"]["diffdock"][2]},
   "metal_adjacent": {"autodock_pct": metal["nearest"]["rank1_failures_no_qualifying_pose_pct"]["metal_adjacent"]["autodock"][2], "diffdock_pct": metal["nearest"]["rank1_failures_no_qualifying_pose_pct"]["metal_adjacent"]["diffdock"][2]},
   "counts_k_of_n": metal["nearest"]["rank1_failures_no_qualifying_pose_pct"]},
  src165 + " kind == no_qualifying_pose_among_rank1_failures", "resolved",
  "58.2 / 70.2 on the metal-free stratum and 56.8 / 70.2 on the metal-adjacent one (64 of 110, 80 of 114, 25 of 44, 40 of 57). The reading 'AutoDock's misses about as often recoverable deeper on either stratum, DiffDock's more often absent' holds.")
E("A165.distance_rule_sensitivity", 165, "App. C, Metal-adjacent stratum", "3.3:165 / W30", "value",
  "4 A rule 73 / 6 A rule 81, contrasts -7.4 and -19.2 and -6.3 and -21.0",
  {"4A": {"n_adjacent": 73, "contrasts_pp": metal["nearest"]["top15_contrast_4A_pp"], "full": metal["nearest"]["top15_contrast_4A_full"]},
   "6A": {"n_adjacent": 81, "contrasts_pp": metal["nearest"]["top15_contrast_6A_pp"], "full": metal["nearest"]["top15_contrast_6A_full"]}},
  src165 + " cutoff_A in (4.0, 6.0), kind == contrast, depth == 15", "resolved",
  "Nearest: -9.1 and -17.8 (4 A) and -8.1 and -19.8 (6 A). 'insensitive in direction if not in size' holds.")
E("A165.wording_rule_and_crossover", 165, "App. C, Metal-adjacent stratum", "W30 / D4", "wording",
  "In 78 of the 303 analysed entries a metal ion or a metal-containing cofactor lies within 5 A of a crystal-ligand heavy atom",
  "W: residue-carrying-a-metal rule on the deposited reference instance, plus the D4 crossover sentence (7JHQ_VAJ counted metal-adjacent though scored against a metal-free copy; 7TB0_UD1 counted metal-free though scored against the one copy with a metal within 5 A). Replacement text is in NEAREST_COPY_TEXT_CHANGES W30.",
  "plan D4; metal_stratum_summary.json rule string", "wording",
  "The 7JHQ_VAJ / 7TB0_UD1 distances (21.8-39.5 A; K at 3.17 A) come from plan D4, not from a rebuilt sidecar; metal_stratum.csv carries per-complex distances if the editor wants to re-verify.")
E("A165.ligand_size_comparisons", 165, "App. C, Metal-adjacent stratum", "3.3:165", "value",
  "medians 25 / 22.5 (p 0.124, delta +0.12); 5 / 5 (p 0.596, -0.04); 358 / 361 Da (p 0.468, +0.06)", "unchanged",
  "stratum membership is copy-invariant (D4), the descriptor comparison does not touch a pose", "unchanged", "")

# ============================================================================ A:167 box
E("A167.box_copies_inside", 167, "App. C, AutoDock Vina box", "3.3:167 / W31", "wording", "No further padding is applied.",
  "No further padding is applied, and every deposited copy of the ligand lies inside the box.", "plan Phase 5.4 (receptor bounding box + 1 A contains every bound copy)", "wording", "")

# ============================================================================ A:173 ranking heads
h = C["line173_heads"]
E("A173.heads_rank1", 173, "App. C, CNNaffinity head", "3.3:173", "value",
  "111 (36.6%) CNNaffinity; 112 (37.0%) CNNscore; 90 (29.7%) minimised Vina energy",
  {"cnn_affinity": [h["cnn_affinity (optimized_rank)"]["1"], h["rank1_pct"]["cnn_affinity (optimized_rank)"]],
   "cnn_score": [h["cnn_score"]["1"], h["rank1_pct"]["cnn_score"]],
   "minimised_vina": [h["minimized_vina_energy"]["1"], h["rank1_pct"]["minimized_vina_energy"]]},
  "recompute: per_pose_metrics_nearest, method autodock_mgltools_exh128_gnina, per-complex order by optimized_rank / cnn_score desc / minimized_affinity asc, gate pb_valid & rmsd <= 2 (method reproduces the canonical 111 / 112 / 90 on the shipped table)", "resolved",
  "149 (49.2 %), 152 (50.2 %), 128 (42.2 %).")
E("A173.heads_discordance", 173, "App. C, CNNaffinity head", "3.3:173", "value", "31 complexes, 15 favouring CNNaffinity and 16 favouring CNNscore, exact McNemar p = 1.000",
  {"discordant": h["rank1_discordant_affinity_only"] + h["rank1_discordant_score_only"], "cnn_affinity_only": h["rank1_discordant_affinity_only"], "cnn_score_only": h["rank1_discordant_score_only"], "p": round(h["rank1_mcnemar_p"], 3)},
  "recompute (same recipe), stats_utils.mcnemar_exact", "resolved", "31 complexes, 14 favouring CNNaffinity and 17 favouring CNNscore, p = 0.720. 'not resolved' holds; 'does not report the arm at its most favourable rank-1 value' holds (CNNscore is higher).")
E("A173.heads_gap_and_depth", 173, "App. C, CNNaffinity head", "3.3:173", "value", "about seven percentage points; top-15 198 / 197 / 190; top-30 all 199",
  {"gap_pp_affinity_vs_vina": r1(h["rank1_pct"]["cnn_affinity (optimized_rank)"] - h["rank1_pct"]["minimized_vina_energy"]),
   "gap_pp_score_vs_vina": r1(h["rank1_pct"]["cnn_score"] - h["rank1_pct"]["minimized_vina_energy"]),
   "top15": [h["cnn_affinity (optimized_rank)"]["15"], h["cnn_score"]["15"], h["minimized_vina_energy"]["15"]],
   "top30": [h["cnn_affinity (optimized_rank)"]["30"], h["cnn_score"]["30"], h["minimized_vina_energy"]["30"]], "top30_same_set": h["same_set_top30"]},
  "recompute (same recipe)", "resolved", "Gaps 7.0 and 8.0 points, write 'by seven to eight percentage points'. Top-15 213 / 213 / 202, top-30 all 214 and the same set.")

# ============================================================================ A:182 selection paragraph
t8h = C["table_8_holm_top15_double"]
hol = {r["arm"]: r for r in t8h}
r1d = {r["arm"]: r for r in C["ladder_rank1_double_seven"]}
t15t = {r["arm"]: r for r in C["ladder_top15_triple_seven"]}
E("A182.top15_selected_and_holm_range", 182, "App. C, Selection paragraph", "3.3:182 / W32 / D14", "value",
  "198 of 303; Holm 2e-21 (exh 18) to 0.007 (own raw search)",
  {"selected_top15": 213, "holm_min": pfmt(hol["autodock_mgltools_exh18"]["p_holm_seven"]), "holm_min_arm": "Exhaustiveness 18",
   "holm_max": pfmt(hol["autodock_mgltools_exh128"]["p_holm_seven"]), "holm_max_arms": ["Exhaustiveness 128 (raw)", "64 + gnina"],
   "all_seven": {k: {"k_arm": v["k_arm"], "sel_only": v["selected_only"], "arm_only": v["arm_only"], "p_raw": pfmt(v["p_raw"]), "p_holm": pfmt(v["p_holm_seven"])} for k, v in hol.items()}},
  "recompute: gate pb_valid & rmsd <= 2, existence over top-15, exact McNemar vs autodock_mgltools_exh128_gnina, stats_utils.holm over the seven contrasts", "resolved",
  "7.3e-20 to 0.0052; the top Holm value is shared by the raw exh128 contrast and the 64 + gnina contrast (both 0.0052), so write 'to 0.005 against its own raw search and against exhaustiveness 64 rescored'.")
dl = C["ladder_depth_lead_summary"]
E("A182.depth_lead", 182, "App. C, Selection paragraph", "3.3:182 / W32", "value", "It leads at every depth from top-3 to top-22.",
  {"strict_lead_depths": f"{dl['strict_lead_depths'][0]}-{dl['strict_lead_depths'][-1]}", "tie_depths": dl["tie_depths"], "trail_depths": dl["trail_depths"]},
  "recompute: selected vs best of the seven other arms on the double gate at d = 1..30", "resolved",
  "Leads strictly at every depth from rank-1 to top-25, ties at top-26 and trails the raw exhaustiveness-128 arm from top-27. The claim is strengthened.")
g64 = r1d["autodock_mgltools_exh64_gnina"]
E("A182.rank1_seven_contrasts", 182, "App. C, Selection paragraph", "3.3:182 / W32 / D14", "value",
  "At rank-1 only two of the seven contrasts survive correction and the arm ties with exhaustiveness 64 rescored, at 111 against 110.",
  {"surviving": sum(v["p_holm_seven"] < 0.05 for v in r1d.values()), "failing": [k for k, v in r1d.items() if v["p_holm_seven"] >= 0.05],
   "vs_64_gnina": {"selected": g64["k_selected"], "arm": g64["k_arm"], "selected_only": g64["selected_only"], "arm_only": g64["arm_only"], "p_raw": round(g64["p_raw"], 3), "p_holm": round(g64["p_holm_seven"], 3)},
   "all_seven": {k: {"k_arm": v["k_arm"], "sel_only": v["selected_only"], "arm_only": v["arm_only"], "p_raw": pfmt(v["p_raw"]), "p_holm": pfmt(v["p_holm_seven"])} for k, v in r1d.items()}},
  "recompute (double gate, rank-1, seven-contrast Holm family)", "resolved",
  "Six of seven survive within the seven-contrast family (D14). The one that does not is 64 + gnina, 149 against 141, 16 against 8 discordant, exact p = 0.152, so the lead is nominal and unresolved. The raw-exh128 contrast has Holm p 0.037 here (survives), against Holm 0.060 in Table 19's sixteen-test family, which is the D14 contradiction the family naming resolves.")
t30 = C["ladder_top30_double_raw128_vs_selected"]
E("A182.top30_raw_ahead", 182, "App. C, Selection paragraph", "3.3:182", "value", "202 against 199", {"raw128": t30["raw128"], "selected": t30["selected"], "raw_only": t30["a_only"], "selected_only": t30["b_only"], "p": t30["p"]},
  "recompute (double gate, top-30)", "resolved", "217 against 214, 4 raw-only against 1 selected-only.")
E("A182.triple_endpoint", 182, "App. C, Selection paragraph", "3.3:182", "value",
  "160 complexes against 152 for its own raw search and 95 at exhaustiveness 18, Holm 6e-18 to 0.039",
  {"selected": t15t["autodock_mgltools_exh128"]["k_selected"], "raw128": t15t["autodock_mgltools_exh128"]["k_arm"], "exh18": t15t["autodock_mgltools_exh18"]["k_arm"],
   "holm_min": pfmt(min(v["p_holm_seven"] for v in t15t.values())), "holm_max": pfmt(max(v["p_holm_seven"] for v in t15t.values())),
   "all_seven": {k: {"k_arm": v["k_arm"], "sel_only": v["selected_only"], "arm_only": v["arm_only"], "p_raw": pfmt(v["p_raw"]), "p_holm": pfmt(v["p_holm_seven"])} for k, v in t15t.items()}},
  "recompute: triple gate (pb_valid & rmsd <= 2 & bestfit_rmsd <= 1), top-15, seven-contrast Holm; cross-check selection_check_nearest/selection_check.txt (175 / 166 / 110)", "resolved",
  "175 against 166 and 110 at exhaustiveness 18, Holm 6.5e-18 to 0.023 (the 0.023 is shared by raw exh128 and 64 + gnina).")
p28 = C["ladder_28pair_holm_top15"]
E("A182.twentyeight_pair_sentence", 182, "App. C, Selection paragraph", "3.3:182 / W32", "value",
  "the step over raw exhaustiveness 128 is the one contrast that would not survive a correction taken over all twenty-eight pairs of the ladder",
  {"triple_gate_holm28_vs_selected": {k.split("|")[0]: round(v["p_holm_28"], 4) for k, v in p28["triple"].items()},
   "triple_failing": [k.split("|")[0] for k, v in p28["triple"].items() if v["p_holm_28"] >= 0.05],
   "double_gate_holm28_vs_selected": {k.split("|")[0]: round(v["p_holm_28"], 4) for k, v in p28["double"].items()},
   "double_failing": [k.split("|")[0] for k, v in p28["double"].items() if v["p_holm_28"] >= 0.05]},
  "recompute: exact McNemar on all 28 arm pairs at top-15, stats_utils.holm over 28, read off the seven pairs involving the selected arm", "resolved",
  "On the three-part endpoint the sentence's gate, the steps over raw exhaustiveness 128 AND over exhaustiveness 64 rescored are the two contrasts that would not survive (Holm-28 p 0.107 each). On the double gate every one of the seven survives even the 28-pair correction (largest 0.035), so the sentence must stay on the triple gate as it does today. Holm step-down is less severe than the plan's 28 x p Bonferroni preview.")

# ============================================================================ A:184 components
l184 = C["line184"]; pv = C["line184_parent_vina_order"]
E("A184.components", 184, "App. C, two components", "3.3:184", "value",
  "32 -> 128 rescored top-15 154 -> 198; rescoring at 128 185 -> 198, gaining 17 and losing 4",
  {"gnina32_top15": l184["gnina32_top15"], "sel_top15": l184["sel_top15"], "raw128_top15": l184["raw128_top15"], "gained": l184["gained_by_rescoring"], "lost": l184["lost_by_rescoring"]},
  "recompute (double gate, top-15, existence)", "resolved", "170 to 213; 198 to 213, gaining 19 and losing 4.")
E("A184.parent_vina_order", 184, "App. C, two components", "3.3:184", "value",
  "Carried at the parent Vina order ... 93 complexes at rank-1 and 183 at top-15, below the raw search at 94 and 185",
  {"rescored_by_autodock_rank": [pv["gnina_poses_by_autodock_rank_d1"], pv["gnina_poses_by_autodock_rank_d15"]], "raw128": [pv["raw128_d1"], pv["raw128_d15"]]},
  "recompute: autodock_mgltools_exh128_gnina poses ordered by autodock_rank (parent Vina rank), double gate (reproduces the canonical 93 / 183 on the shipped table)", "resolved",
  "129 and 196, below the raw search at 130 and 198. Direction holds.")

# ============================================================================ captions :189 :195 :235 :241
st = C["ladder_steps_rank1"]; s92 = C["step_92_128_by_depth"]
def first_sig(seq):
    for r in seq:
        if r["p"] < 0.05: return r["depth"]
    return None
E("A189.fig18_caption", 189, "fig:appendix-exh-rescored", "3.3:189", "wording", "gains nothing on placement and gives a little back on the three-part gate",
  "unchanged (pool near-native 217 -> 217; pool triple 178 -> 175 at exh128)", "recompute (pool gates) and pose_validity_cascade.csv", "unchanged", "Regenerate image44 from autodock_exhaustiveness_returns_nearest/figures/exh_raw_vs_rescored.png.")
E("A195.fig19_caption", 195, "fig:appendix-exh-depth", "3.3:195 / W33", "value",
  "the top step is worth three complexes at rank-1 but does not become individually resolvable until about depth 22",
  {"near_gate_A": {"delta_rank1": st["near:autodock_mgltools_exh92->autodock_mgltools_exh128"]["k_b"] - st["near:autodock_mgltools_exh92->autodock_mgltools_exh128"]["k_a"], "gained": st["near:autodock_mgltools_exh92->autodock_mgltools_exh128"]["b_only"], "lost": st["near:autodock_mgltools_exh92->autodock_mgltools_exh128"]["a_only"], "p_rank1": round(st["near:autodock_mgltools_exh92->autodock_mgltools_exh128"]["p"], 4), "first_depth_p_lt_0.05": first_sig(s92["near"])},
   "triple_gate_B": {"delta_rank1": st["triple:autodock_mgltools_exh92->autodock_mgltools_exh128"]["k_b"] - st["triple:autodock_mgltools_exh92->autodock_mgltools_exh128"]["k_a"], "gained": st["triple:autodock_mgltools_exh92->autodock_mgltools_exh128"]["b_only"], "lost": st["triple:autodock_mgltools_exh92->autodock_mgltools_exh128"]["a_only"], "p_rank1": round(st["triple:autodock_mgltools_exh92->autodock_mgltools_exh128"]["p"], 4), "first_depth_p_lt_0.05": first_sig(s92["triple"])}},
  "recompute (exh92 vs exh128 exact McNemar by depth, both caption gates); cross-check " + S(EXH / "exh_yield_contrasts.csv") + " rmsd2 rank-1 exh92->exh128 (7, p 0.0391) and " + S(EXH / "exh_returns_report.txt") + " section 7", "resolved",
  "On the near-nativeness gate (A) the top step is worth seven complexes at rank-1 (8 gained, 1 lost) and is already resolvable there at an unadjusted p of 0.039. On the three-part gate (B) it is worth six at rank-1 (p 0.109) and first resolves at depth 2. W33 text uses the (A) reading; the editor may add the (B) clause.")
E("A235.fig20_caption", 235, "fig:appendix-exh-cost", "3.3:235", "wording", "Exhaustiveness 92 is the most expensive rung on the ladder and buys nothing at rank-1",
  "unchanged (exh64 -> exh92 at rank-1: near 124 -> 124, double 123 -> 123, triple 101 -> 100)", "recompute; " + S(EXH / "exh_returns_report.txt") + " section 5 (exh64->exh92 rank-1 +0 / -1)", "unchanged", "Regenerate image46.")
E("A241.fig21_caption", 241, "fig:appendix-exh-return", "3.3:241", "wording", "The top step carries no bar because its measured cost delta is negative", "unchanged (-0.110 wall-h is a cost fact)", S(EXH / "exh_returns_report.txt") + " section 2", "unchanged", "Regenerate image47.")

# ============================================================================ Table 8 :218-227
t8 = C["table_8"]
order = ["autodock_mgltools_exh18", "autodock_mgltools", "autodock_mgltools_gnina", "autodock_mgltools_exh64", "autodock_mgltools_exh64_gnina", "autodock_mgltools_exh92", "autodock_mgltools_exh128", "autodock_mgltools_exh128_gnina"]
cur_rows = {"autodock_mgltools_exh18": [72, 108, 121, 125, "2e-21", 2.14], "autodock_mgltools": [85, 130, 148, 154, "5e-12", 2.56], "autodock_mgltools_gnina": [100, 144, 154, 154, "2e-09", 2.56],
            "autodock_mgltools_exh64": [94, 145, 175, 182, "1.4e-04", 3.14], "autodock_mgltools_exh64_gnina": [110, 169, 182, 182, "0.0026", 3.14], "autodock_mgltools_exh92": [91, 150, 182, 194, "0.0031", 4.84],
            "autodock_mgltools_exh128": [94, 153, 185, 202, "0.0072", 4.73], "autodock_mgltools_exh128_gnina": [111, 184, 198, 199, "---", 4.73]}
search_h = {"autodock_mgltools_exh18": 2.14, "autodock_mgltools": 2.56, "autodock_mgltools_gnina": 2.56, "autodock_mgltools_exh64": 3.14, "autodock_mgltools_exh64_gnina": 3.14, "autodock_mgltools_exh92": 4.84, "autodock_mgltools_exh128": 4.73, "autodock_mgltools_exh128_gnina": 4.73}
for i, m in enumerate(order):
    row = t8[m]["double"]
    holmv = "---" if m == order[-1] else pfmt(hol[m]["p_holm_seven"])
    E(f"A{218+i}.table8.{m}", 218 + i, "tab:appendix-exhaustiveness-ladder", "3.3:218-227 / plan 1.5", "table_cell",
      cur_rows[m], [row["1"], row["5"], row["15"], row["30"], holmv, search_h[m]],
      "harness signature table_8[...] actual == recompute (double gate by rank / optimized_rank); Holm column recompute (exact McNemar vs selected at top-15, Holm over seven); Search h unchanged from " + S(EXH / "exh_cost_per_arm.csv"),
      "resolved", f"signature: {[sig_map.get(f'table_8[{m}].d{d}') for d in (1,5,15,30)]}. The task brief names exh_yield_contrasts.csv for the Holm column, but that file holds the raw-ladder STEP contrasts (exh18->32 etc.), not the arm-vs-selected family the footnote defines, so the column is recomputed.")
E("A227.table8_footnote", 227, "tab:appendix-exhaustiveness-ladder footnote", "3.3:218-227 / W34", "wording",
  "Complexes of 303 with at least one pose that is PoseBusters-valid and within 2 A among the first d ranks.",
  "... within 2 A of the nearest deposited copy of the ligand among the first d ranks. (rest unchanged: Search h, threads, 0.110 h inversion)", "W34", "wording", "")

# ============================================================================ A:230 D10 confirmatory cell
def exh_tables(root):
    mr = pd.read_csv(root / "exh_marginal_return.csv"); dr = pd.read_csv(root / "exh_diminishing_returns.csv"); yd = pd.read_csv(root / "exh_yield_by_depth.csv")
    tri = mr[(mr.gate == "triple") & (mr.depth == 0)]
    rates = {r.segment: {"rate_per_h": r2(r.rate_per_h) if pd.notna(r.rate_per_h) else None, "ci": [r2(r.rate_lo), r2(r.rate_hi)] if pd.notna(r.rate_lo) else None, "delta_yield": int(r.delta_yield), "delta_cost_h": round(float(r.delta_cost_h), 3), "wall_h_per_extra_complex": (round(float(r.wall_h_per_extra_complex), 4) if pd.notna(r.wall_h_per_extra_complex) else None), "status": r.rate_status} for _, r in tri.iterrows()}
    d = dr[(dr.gate == "triple") & (dr.depth == 0)].iloc[0]
    pool = {int(r.exhaustiveness): int(r.k_complexes) for _, r in yd[(yd.gate == "triple") & (yd.depth == 0)].iterrows()} if (yd.depth == 0).any() else None
    return {"per_hour": rates, "contrast": {"value": r2(d.contrast), "ci": [r2(d.lo), r2(d.hi)], "first": r2(d.rate_first_per_h), "last": r2(d.rate_last_per_h)}, "pool_triple": pool}
near_exh = exh_tables(EXH); can_exh = exh_tables(CANON_EXH)
casc = pd.read_csv(HUB / "pose_validity_cascade.csv"); cascC = pd.read_csv(ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/pose_validity_cascade.csv")
raw_keys = ["autodock_mgltools_exh18", "autodock_mgltools", "autodock_mgltools_exh64", "autodock_mgltools_exh92", "autodock_mgltools_exh128"]
pool_near = [int(casc[casc.method_key == k].triple_complexes.iloc[0]) for k in raw_keys]
pool_inst = [int(cascC[cascC.method_key == k].triple_complexes.iloc[0]) for k in raw_keys]
E("A230.confirmatory_cell_D10", 230, "App. C, per-hour returns", "3.3:230 / W35 / D10", "decision",
  "+71.6, +35.9 and +5.9 complexes per hour; contrast +65.6 (+32.1 to +137.2); rungs 98, 128, 149, 159, 164; first step adds thirty, last adds five",
  {"KEEP_instance": {"per_hour": [can_exh["per_hour"]["exh18->exh32"]["rate_per_h"], can_exh["per_hour"]["exh32->exh64"]["rate_per_h"], can_exh["per_hour"]["exh64->exh92"]["rate_per_h"]],
                     "contrast": can_exh["contrast"], "pool": pool_inst, "first_step": pool_inst[1] - pool_inst[0], "last_step": pool_inst[4] - pool_inst[3]},
   "ADD_nearest_sentence": {"per_hour": [near_exh["per_hour"]["exh18->exh32"]["rate_per_h"], near_exh["per_hour"]["exh32->exh64"]["rate_per_h"], near_exh["per_hour"]["exh64->exh92"]["rate_per_h"]],
                            "per_hour_ci": [near_exh["per_hour"]["exh18->exh32"]["ci"], near_exh["per_hour"]["exh32->exh64"]["ci"], near_exh["per_hour"]["exh64->exh92"]["ci"]],
                            "contrast": near_exh["contrast"], "pool": pool_near, "first_step": pool_near[1] - pool_near[0], "last_step": pool_near[4] - pool_near[3]}},
  "instance (kept, READ-ONLY): " + S(CANON_EXH / "exh_marginal_return.csv") + ", exh_diminishing_returns.csv (gate triple, depth 0 = oracle), canonical pose_validity_cascade.csv triple_complexes; nearest: same files under " + S(EXH) + " and " + S(HUB / "pose_validity_cascade.csv"),
  "decision", "D10 adopted as in the task brief: the whole paragraph stays on the single instance (98 / 128 / 149 / 159 / 164, +71.6 / +35.9 / +5.9, +65.6 [+32.1, +137.2]) and ONE sentence adds the nearest values: rungs 115, 143, 166, 172 and 178 (first step +28, last +6), per-hour returns +66.8, +39.3 and +3.5, contrast +62.4 (95 % +32.0 to +133.3). W35 text carries the sentence with [brackets] to fill from this entry.")

# ============================================================================ A:244 :246 :248
l244 = C["line244"]
E("A244.top_step", 244, "App. C, fourth step", "3.3:244-248", "value", "adds five complexes over the pool and three at rank-1; -0.110 hours; 18.7 percentage points",
  {"triple_pool_delta": l244["triple"]["pool"], "triple_rank1_delta": l244["triple"]["rank1"], "double_pool_delta": l244["double"]["pool"], "double_rank1_delta": l244["double"]["rank1"], "wall_h_delta": -0.110, "residual_arm_effect_pp": 18.7},
  "recompute (canonical gate of the sentence is the three-part one: instance 164-159 = 5 and 77-74 = 3 reproduce the print); " + S(EXH / "exh_returns_report.txt") + " section 2 (-0.110 wall-h, 18.7 pp unchanged)", "resolved",
  "Exhaustiveness 92 to 128 adds six complexes over the pool and six at rank-1 on the three-part gate. The -0.110 h and 18.7 pp are cost-side facts and do not move.")
mrn = near_exh["per_hour"]; mrc = can_exh["per_hour"]
E("A246.minutes_per_complex", 246, "App. C, what follows", "3.3:244-248", "value", "about ten minutes of wall time per additional complex against about one minute for the first",
  {"last_measurable_step_min": round(60 * mrn["exh64->exh92"]["wall_h_per_extra_complex"]), "first_step_min": round(60 * mrn["exh18->exh32"]["wall_h_per_extra_complex"]),
   "instance_check": {"last": round(60 * mrc["exh64->exh92"]["wall_h_per_extra_complex"]), "first": round(60 * mrc["exh18->exh32"]["wall_h_per_extra_complex"])}},
  S(EXH / "exh_marginal_return.csv") + " gate triple, depth 0, wall_h_per_extra_complex (report section 8 budget crossing: 1 / 2 / 17 min)", "resolved",
  "About seventeen minutes against about one minute (instance 10 / 1 reproduces the print). NOTE this sentence sits in the paragraph D10 keeps on the instance convention; if the paragraph stays on instance the 'ten minutes' stays, and the nearest figure belongs in the added sentence. Editor decides placement.")
l246 = C["line246"]
E("A246.plateau_argument", 246, "App. C, what follows", "3.3:246 / W36", "value",
  "raw rank-1 recovery tops out at 95 complexes anywhere above exhaustiveness 64, while rescoring lifts it to 113 and doubling the pool beneath the rescorer leaves it at exactly 113",
  {"raw_rank1_near_64_92_128": [l246["raw_rank1_near"]["autodock_mgltools_exh64"], l246["raw_rank1_near"]["autodock_mgltools_exh92"], l246["raw_rank1_near"]["autodock_mgltools_exh128"]],
   "rescored_rank1_near_64_128": [l246["rescored_rank1_near"]["autodock_mgltools_exh64_gnina"], l246["rescored_rank1_near"]["autodock_mgltools_exh128_gnina"]]},
  "recompute (rmsd <= 2, no validity, rank-1); cross-check " + S(EXH / "exh_raw_vs_rescored.csv") + " and report section 9 (124 -> 144, 131 -> 151)", "resolved",
  "124 / 124 / 131 raw and 144 / 151 rescored, so the ranking-ceiling argument collapses and is rewritten per W36 ('One thing follows, one follows in part, and a third does not').")
wc = C["line248_winners_curse"]
E("A248.winners_curse_eight_arms", 248, "App. C, caveat", "3.3:244-248", "value", "bias near two percentage points at rank-1 and near one and a half at top-15",
  None, "NO registered generator: grep for winner / curse / selection-bias over Scripts/Analysis and REGENERATE.md finds nothing that produces the eight-arm figure; the memory note 'variant-selection-bias-quantified' covers the per-family (three-arm) estimator only",
  "unresolved", f"A normal order-statistic approximation (E[max of 8] = 1.4236 sigma at the selected arm's rate) gives {wc['d1']['bias_pp_normal_order_stat_8_arms']} pp at rank-1 and {wc['d15']['bias_pp_normal_order_stat_8_arms']} pp at top-15, which does NOT reproduce the printed 2 / 1.5 on the instance rates either, so the printed estimator is unknown. Leave the two figures unchanged and flag, or locate the author's estimator.")

# ============================================================================ A:267-269 EquiBind configurations
eq = C["line267_equibind"]
cz = {r.method_key: r for _, r in casc.iterrows()}
E("A267.triple_top15", 267, "App. C, pocket-guided EquiBind", "3.3:267-269 / W37", "value", "2, 46 and 55 ... against 0, 7 and 14 under fpocket and 0, 5 and 18 under P2Rank, trails by 37",
  {"unguided_raw_smina_gnina": [int(cz["equibind_unguided_raw"].triple_complexes), eq["equibind_unguided_smina"]["triple_top15"], eq["equibind_unguided_gnina"]["triple_top15"]],
   "fpocket": [eq["equibind_fpocket_raw"]["triple_top15"], eq["equibind_fpocket_smina"]["triple_top15"], eq["equibind_fpocket_gnina"]["triple_top15"]],
   "p2rank": [eq["equibind_p2rank_raw"]["triple_top15"], eq["equibind_p2rank_smina"]["triple_top15"], eq["equibind_p2rank_gnina"]["triple_top15"]],
   "gap_best_guided_to_selected": eq["equibind_unguided_gnina"]["triple_top15"] - max(eq["equibind_fpocket_gnina"]["triple_top15"], eq["equibind_p2rank_gnina"]["triple_top15"])},
  S(HUB / "pose_validity_cascade.csv") + " triple_complexes (pool; the raw arms carry 2 / 0 / 0 over the whole pool so their top-15 value is bounded by it) and recompute top-15 by each refiner's own affinity (smina_affinity / gnina_affinity ascending); cross-check " + S(ROOT / "posebusters_results/selection_check_nearest/selection_check.txt"),
  "resolved", "2, 48 and 57 against 0, 10 and 20 under fpocket and 0, 9 and 23 under P2Rank, so the best guided configuration trails the selected unguided one by 34. The unguided-raw top-15 value of 2 is the pool value (both near-native complexes sit inside the 30-pose pool; the raw arm has no rank and the thesis orders it by unguided_NNN).")
E("A267.near_native_pool", 267, "App. C, pocket-guided EquiBind", "3.3:267-269 / W37", "value", "19, 64 and 83 unguided against 0, 14 and 24 and 0, 18 and 35, both raw guided arms recovering no near-native complex at all",
  {"unguided": [eq["equibind_unguided_raw"]["near_pool"] or int(cz["equibind_unguided_raw"].rmsd2_complexes), eq["equibind_unguided_smina"]["near_pool"], eq["equibind_unguided_gnina"]["near_pool"]],
   "fpocket": [int(cz["equibind_fpocket_raw"].rmsd2_complexes), eq["equibind_fpocket_smina"]["near_pool"], eq["equibind_fpocket_gnina"]["near_pool"]],
   "p2rank": [int(cz["equibind_p2rank_raw"].rmsd2_complexes), eq["equibind_p2rank_smina"]["near_pool"], eq["equibind_p2rank_gnina"]["near_pool"]],
   "p2rank_raw_poses": int(cz["equibind_p2rank_raw"].rmsd2_poses)},
  S(HUB / "pose_validity_cascade.csv") + " rmsd2_complexes", "resolved",
  "19, 67 and 86 against 0, 19 and 32 and 1, 28 and 43. The P2Rank raw arm now recovers one complex (two poses), so 'both raw guided arms recovering no near-native complex at all' becomes W37's 'the fpocket raw arm recovering no near-native complex and the P2Rank raw arm one'.")
E("A267.validity_pct", 267, "App. C, pocket-guided EquiBind", "3.3:267", "value", "2.9%, 51.3% and 62.0% unguided against 0.3%, 15.2% and 28.4% and 0.6%, 20.7% and 38.0%", "unchanged", S(HUB / "pose_validity_cascade.csv") + " pb_valid_%", "unchanged", "")
forms = {k: int(cz[k].kabsch1_complexes) for k in eq}
E("A267.form_range", 267, "App. C, pocket-guided EquiBind", "3.3:267", "value", "between 180 and 200 complexes across all nine configurations",
  {"min": min(forms.values()), "max": max(forms.values()), "per_configuration": forms}, S(HUB / "pose_validity_cascade.csv") + " kabsch1_complexes", "resolved", "Between 179 and 199 across the nine configurations.")
E("A269.pbvalid_complexes", 269, "App. C, pocket-guided EquiBind", "3.3:269", "value", "261 for fpocket with gnina and 282 for P2Rank with gnina against 266 unguided",
  [int(cz["equibind_fpocket_gnina"].pb_valid_complexes), int(cz["equibind_p2rank_gnina"].pb_valid_complexes), int(cz["equibind_unguided_gnina"].pb_valid_complexes)], S(HUB / "pose_validity_cascade.csv") + " pb_valid_complexes", "unchanged", "38 complexes / 570-pose shortfall / 8,520 total are production facts, unchanged.")

# ============================================================================ wording-only definitions
E("A276.accuracy_definition", 276, "App. C, refinement paragraph", "3.3:276 / W38", "wording", "symmetry-corrected heavy-atom RMSD to the crystal ligand", "... to the nearest deposited copy of the crystal ligand", "W38", "wording", "")
E("A487.load_all_clause", 487, "App. D, PoseBusters checks", "3.3:487 / W39", "wording", "It is kept off the validity axis here because near-nativeness is gated separately.",
  "+ 'That test loads every deposited instance of the ligand and reports the lowest RMSD, which is the copy rule the Methods adopt.'", "W39 (redock.yml:301-307 load_all: True; modules/rmsd.py:52-70)", "wording", "optional clause")

# ============================================================================ A:543 / :549 / :551 interaction audit
au = dict(pd.read_csv(PMR / "interaction_pose_basis_audit.csv").values)
nr = pd.read_csv(PMR / "native_recovery_by_rank.csv"); nr1 = nr[nr["rank"] == 1].set_index("method")
E("A543.profiled_pool", 543, "App. B, Interaction Fingerprints", "3.3:545-553 / W40a (anchor drifted to :543)", "value", "4,253 poses over 303 complexes with a median deviation of 7.57 A ... 21.3% of poses within 2 A",
  {"poses": int(au["pool_n_poses"]), "median_rmsd_A": float(au["pool_median_rmsd_A"]), "within_2A_pct": float(au["pool_within_2A_pct"])}, S(PMR / "interaction_pose_basis_audit.csv") + " pool_*", "resolved", "4,253 (unchanged), 5.01 A, 26.9 %.")
E("A543.standardised_f1", 543, "App. B, Interaction Fingerprints", "3.3:545-553 / W40a", "value", "F1 triple 0.565, 0.529 and 0.471 to 0.532, 0.527 and 0.513",
  {"crude": [round(nr1.loc["autodock_mgltools_exh128_gnina"].f1, 3), round(nr1.loc["diffdock_smina"].f1, 3), round(nr1.loc["equibind_unguided_gnina"].f1, 3)],
   "standardised": [0.667, 0.647, 0.628], "spread_crude": float(au["rank1_spread_crude"]), "spread_standardised": float(au["rank1_spread_standardised"])},
  S(PMR / "interaction_pose_basis_audit.txt") + " section B (crude / standardised columns) and " + S(PMR / "native_recovery_by_rank.csv") + " rank 1", "resolved",
  "0.732, 0.694 and 0.499 to 0.667, 0.647 and 0.628 (standardised values are printed to four places in the .txt, 0.6673 / 0.6469 / 0.6282).")
E("A543.placement_share", 543, "App. B, Interaction Fingerprints", "3.3:545-553 / W40a", "value", "Roughly four fifths ... 80 to 87% across band schemes of three to fourteen strata",
  {"reported_5band_pct": float(au["rank1_placement_share_pct"]), "range_pct": [float(au["placement_share_3-band_pct"]), float(au["placement_share_14-band_pct"])],
   "by_scheme": {k: float(v) for k, v in au.items() if k.startswith("placement_share")}}, S(PMR / "interaction_pose_basis_audit.csv") + " placement_share_*", "resolved", "83.2 %, running from 82 to 89 % (81.9 to 88.5).")
E("A543.standardised_order_swap", 543, "App. B, Interaction Fingerprints", "3.3:545-553", "wording", "The AutoDock Vina and DiffDock standardised values swap under finer banding and carry no ordering between them.",
  "W: under nearest the standardised order is autodock > diffdock > equibind at every band scheme (3, 5, 8, 14), so the swap sentence is FALSE and must be replaced, e.g. 'The standardised ordering AutoDock Vina above DiffDock above EquiBind holds at every banding, although the spread between the tools shrinks to 0.039.'",
  S(PMR / "interaction_pose_basis_audit.txt") + " section B 'standardised order' lines", "wording", "New defect the change list did not carry; recorded as a decision the editor must take.")
E("A543.site_filter", 543, "App. B, Interaction Fingerprints", "3.3:545-553", "value", "0.683 over 257 (AutoDock), 0.714 over 240 (DiffDock), 0.629 over 155 (EquiBind); p = 0.014 to 5.7e-5; holds 8 to 20 A; discards 44, 60 and 111; 21.3 to 38.1%",
  {"overlap": [float(au["jaccard_rmsd_autodock_mgltools_exh128_gnina"]), float(au["jaccard_rmsd_diffdock_smina"]), float(au["jaccard_rmsd_equibind_unguided_gnina"])],
   "complexes": [301 - int(au["site_filter_dropped_autodock_mgltools_exh128_gnina"]), 300 - int(au["site_filter_dropped_diffdock_smina"]), 266 - int(au["site_filter_dropped_equibind_unguided_gnina"])],
   "p_published_common": round(float(au["jaccard_published_p_dd_vs_ad_common"]), 4), "p_filtered_common": float(f"{float(au['jaccard_rmsd_p_dd_vs_ad_common']):.2g}"),
   "sweep_p_8_to_20": {k.replace("site_sweep_", "").replace("_p_common", ""): float(f"{float(v):.2g}") for k, v in au.items() if k.startswith("site_sweep")},
   "dropped": [int(au["site_filter_dropped_autodock_mgltools_exh128_gnina"]), int(au["site_filter_dropped_diffdock_smina"]), int(au["site_filter_dropped_equibind_unguided_gnina"])],
   "within_2A_from_to": [float(au["pool_within_2A_pct"]), float(au["site_filter_within_2A_pct"])]},
  S(PMR / "interaction_pose_basis_audit.txt") + " section C (published 0.5391 (301) / 0.5722 (300) / 0.4030 (266); rmsd < 10 A 0.6757 (274) / 0.7061 (262) / 0.6228 (166))", "resolved",
  "0.676 over 274, 0.706 over 262 and 0.623 over 166; DiffDock lead moves from p = 0.011 to p = 2.0e-4 on the shared cohort (common n 264 -> 144); holds at every threshold from 8 to 20 A (all p < 0.05); discards 27, 38 and 100; share within 2 A from 26.9 to 38.8 %. NOTE the 'p = 0.014' the current text quotes is the Results' published value; under nearest the published shared-cohort p is 0.0108, so the Results sentence and this one must agree on 0.011.")
E("A549.metal_far_field", 549, "App. B, Interaction Fingerprints, third defect", "3.3:545-553 (anchor :549)", "value",
  "6,296 records over 1,662 poses in 121 complexes; beyond 20 A no-metal poses mean F1 0.0003, three of 750 recover any contact; metal poses 0.1112, not one of 515 scores zero",
  {"records": int(au["metal_rows"]), "poses": 1662, "complexes": int(au["metal_complexes"]), "pose_invariant_complexes": int(au["metal_complexes_pose_invariant"]),
   "far_field_no_metal": {"n": int(au["far_field_metal_False_n"]), "mean_f1": float(au["far_field_metal_False_mean_f1"]), "any_contact": int(au["far_field_metal_False_any_contact"]), "zero_pct": float(au["far_field_metal_False_zero_pct"])},
   "far_field_metal": {"n": int(au["far_field_metal_True_n"]), "mean_f1": float(au["far_field_metal_True_mean_f1"]), "zero_pct": float(au["far_field_metal_True_zero_pct"])}},
  S(PMR / "interaction_pose_basis_audit.txt") + " section D", "resolved", "6,296 / 1,662 / 121 unchanged; beyond 20 A those carrying no metal contact recover a mean F1 of 0.0012 and only six of 335 recover any native contact, whereas those carrying one recover 0.1419 and not one of the 233 scores zero.")
E("A549.defect_effect_on_overlap", 549, "App. B, Interaction Fingerprints, third defect", "3.3:545-553 (anchor :549)", "value",
  "0.417, 0.450 and 0.380 to 0.407, 0.440 and 0.366; both repaired 0.408, 0.440 and 0.367; p 0.014 -> 0.017 (metal removed) -> 0.016 (both)",
  {"published": [float(au["defect_published_autodock_mgltools_exh128_gnina"]), float(au["defect_published_diffdock_smina"]), float(au["defect_published_equibind_unguided_gnina"])],
   "metal_dropped": [float(au["defect_metal_class_dropped_autodock_mgltools_exh128_gnina"]), float(au["defect_metal_class_dropped_diffdock_smina"]), float(au["defect_metal_class_dropped_equibind_unguided_gnina"])],
   "both_repaired": [float(au["defect_both_repaired_autodock_mgltools_exh128_gnina"]), float(au["defect_both_repaired_diffdock_smina"]), float(au["defect_both_repaired_equibind_unguided_gnina"])],
   "p_common": {"published": float(au["defect_published_p_dd_vs_ad_common"]), "metal_dropped": float(au["defect_metal_class_dropped_p_dd_vs_ad_common"]), "both_repaired": float(au["defect_both_repaired_p_dd_vs_ad_common"])}},
  S(PMR / "interaction_pose_basis_audit.csv") + " defect_*", "resolved", "0.539, 0.572 and 0.403 to 0.531, 0.564 and 0.389; both repaired 0.532, 0.566 and 0.391; p 0.011 -> 0.013 -> 0.012 on the shared cohort (n 264). Ordering and significance unaffected.")
E("A551.per_copy_fingerprint_sentence", 551, "App. B, Interaction Fingerprints", "W40b / plan 1.8", "wording", "(canonical receptor sentence)",
  "+ W40b sentence on the per-copy crystal fingerprint and nearest-copy scoring; state the union rule (a top-5 union is scored against the union of the copies its poses are nearest to)",
  S(PMR / "crystal_copy_selection.txt") + " and audit header ('each top-5 union is scored against the union of the crystal copies its poses are nearest to (138 of 303 pairs carry >1 fingerprinted copy)')", "wording",
  "Sidecar facts for the sentence: 138 of 303 pairs multi-copy; 60,539 of 241,913 poses nearest to an alternate copy; 9 profiled poses without an index fell back to the reference (3 pairs).")

# ============================================================================ A:664 :666 :722 :784
E("A664.footnote_file_roles", 664, "App. E, benchmark files footnote", "A6 / 3.3:664", "wording", "the single chosen crystal instance that serves as the ground-truth answer pose (*_ligand.sdf)",
  "A6 replacement (documented README roles); under Part B the clause 'and which this study scores against' becomes 'and which is the first record of the multi-instance file'", "Data/PoseBuster Benchmark Set/README.txt:51-63; reference_is_record0_in 303 of 303 (" + S(MET / "metal_stratum_summary.json") + ")", "wording", "")
rc = json.load(open(REF / "reference_convention_summary.json"))
fl = C["table1_form_loss"]
E("A666.convention_paragraph", 666, "App. E, leakage paragraph", "B2 / D3 / D11 / plan 1.3", "value",
  "the comparison is made against the single deposited reference instance rather than against the closest copy ... within that one instance and not equivalence between copies",
  {"multi_copy_308": rc["count_308"]["multi_copy"], "multi_copy_303": rc["multi_copy_303"], "median_alternate_separation_A_over_143": 36.0,
   "form_rows_losing_one_complex": fl["rows_losing_a_complex"], "n_form_rows_losing": len(fl["rows_losing_a_complex"]),
   "alternates_outside_reference_pocket": ["7A9E_R4W", "7TUO_KL9", "7VKZ_NOJ", "7Z1Q_NIO"], "complexes_gaining_recovery_through_them": 0},
  "143 / 138 from " + S(REF / "reference_convention_summary.json") + "; Form-loss rows recomputed on per_pose_metrics_nearest (bestfit_rmsd <= 1 & rmsd < 1000 against the *_ref_instance twins) for the ten Table 1 rows; 36.0 A median, the four alternates and the 0-gain fact are plan Data facts / D3 (verified twice per plan) with NO rebuilt sidecar",
  "resolved", "B2 text: 'one complex in each of five printed variants of Table 1' is confirmed (DiffDock raw 225 -> 224, DiffDock + smina 237 -> 236, EquiBind + gnina 200 -> 199, fpocket raw 180 -> 179, P2Rank raw 181 -> 180; EquiBind raw and + smina GAIN one each). The 36.0 A median and the four pocket-foreign alternates are not carried by any _nearest sidecar; flagged so the author can register a generator or accept the plan's verified fact.")
E("A722.multi_copy_share", 722, "App. E, fig:appendix-dataset-3 prose", "3.3:722 / W43", "value", "(no current text)", {"multi_copy_308": rc["count_308"]["multi_copy"], "share_pct": round(100 * rc["count_308"]["multi_copy"] / 308, 1), "multi_copy_303": rc["multi_copy_303"]},
  S(REF / "reference_convention_summary.json") + " count_308 / multi_copy_303", "resolved", "In 143 of the 308 entries, 46 %, ... 138 of those among the 303 analysed.")
E("A784.start_conformer_D13", 784, "App. F, start-conformer difficulty", "3.3:784 / W44 / D13", "wording", "1.56 A median; 65% within 2 A", "unchanged (kept on the single instance); add the D13 clause", "D13", "unchanged", "")

# ============================================================================ A:871 lead-in and Wilson table :893-920
w = C["wilson_table"]
E("A871.leadin_gates", 871, "App. H lead-in", "A4 / 3.3:871", "value", "at 34.0% and 34.3% rank-1 validity-aware recovery respectively",
  {"validity_aware_rank1": [w["DiffDock + smina"]["1"]["valid_pct"], w["DiffDock (gnina-opt)"]["1"]["valid_pct"]], "near_native_rank1": [w["DiffDock + smina"]["1"]["near_pct"], w["DiffDock (gnina-opt)"]["1"]["near_pct"]]},
  S(HUB / "topk_recovery_validity.csv") + " k == 1", "resolved", "A4 wording with 45.2 % against 45.5 % (near-native) and 43.6 % against 44.2 % (validity-aware).")
wil_order = [("AutoDock Vina (raw)", "AutoDock Vina", 893), ("AutoDock Vina + gnina", "AutoDock (gnina-opt)", 897), ("DiffDock (raw)", "DiffDock (raw)", 901), ("DiffDock + smina", "DiffDock + smina", 905),
             ("DiffDock + gnina", "DiffDock (gnina-opt)", 909), ("EquiBind (raw)", "EquiBind (raw)", 913), ("EquiBind + gnina", "EquiBind (gnina-opt)", 917)]
for label, key, ln in wil_order:
    for j, k in enumerate((1, 5, 10, 15)):
        cur = L(ln + j).strip()
        c = w[key][str(k)]
        E(f"A{ln+j}.table18.{key}.k{k}", ln + j, "tab:appendix-top-k-recovery", "3.3:893-920", "table_cell", cur,
          {"near": f"{c['near_pct']} [{c['near_ci'][0]}, {c['near_ci'][1]}]", "validity_aware": f"{c['valid_pct']} [{c['valid_ci'][0]}, {c['valid_ci'][1]}]", "near_k": c["near_k"], "valid_k": c["valid_k"]},
          S(HUB / "topk_recovery_validity.csv") + f" (variant '{key}', k {k}) + stats_utils.wilson_ci(k, 303); harness signature table_18 rows agree", "resolved", "")

# ============================================================================ Table 19 :947-964 and footnote
tks = pd.read_csv(HUB / "topk_recovery_stats.csv")
def star(p): return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
t19_rows = [("AutoDock Vina", "AutoDock (gnina-opt)", "AutoDock Vina + gnina", 947), ("DiffDock (raw)", "DiffDock + smina", "DiffDock + smina", 951), ("DiffDock (raw)", "DiffDock (gnina-opt)", "DiffDock + gnina", 955), ("EquiBind (raw)", "EquiBind (gnina-opt)", "EquiBind + gnina", 959)]
for base, comp, lab, ln in t19_rows:
    for j, k in enumerate((1, 5, 10, 15)):
        r = tks[(tks.baseline == base) & (tks.comparison == comp) & (tks.k == k) & (tks.family == "within_tool")].iloc[0]
        E(f"A{ln+j}.table19.{lab}.k{k}", ln + j, "tab:appendix-refinement-mcnemar", "3.3:947-964 / plan 1.5 star flips", "table_cell", L(ln + j).strip(),
          {"raw_pct": r1(r["baseline_rate_%"]), "opt_pct": r1(r["comparison_rate_%"]), "raw_only": int(r.baseline_only_wins), "opt_only": int(r.comparison_only_wins),
           "p_raw": pfmt(r.mcnemar_p), "p_holm": pfmt(r.p_holm), "sig": r.star},
          S(HUB / "topk_recovery_stats.csv") + f" (family within_tool, {base} -> {comp}, k {k}); harness signature table_19 rows agree", "resolved",
          "star flips vs print: " + ("AutoDock k=1 stays ns (Holm 0.060)" if lab.startswith("AutoDock") and k == 1 else "DiffDock k=5 loses its star (Holm 0.068)" if "DiffDock" in lab and k == 5 else "AutoDock k=5 / 10 / 15 become ***" if lab.startswith("AutoDock") and k > 1 else ""))
fnr = C["table19_footnote_rerank"]
rer_s = pd.read_csv(ROOT / "PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv").set_index("strategy")
rer_g = pd.read_csv("/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/gnina_rerank_nearest/selection_strategy_summary_gnina.csv").set_index("strategy")
E("A964.table19_footnote_rerank", 964, "tab:appendix-refinement-mcnemar footnote", "3.3:947-964", "value", "32.3% of complexes with gnina and 32.7% with smina, against the 33.7% of the confidence order",
  {"gnina_rerank_pct": r1(rer_g.loc["C_rerank"].success_pct), "smina_rerank_pct": r1(rer_s.loc["C_rerank"].success_pct), "confidence_pct": r1(rer_s.loc["A_native"].success_pct)},
  S(ROOT / "PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv") + " C_rerank / A_native; gnina value from the same registered script run with --tool gnina on the nearest table into the scratchpad (the program has NO gnina_rerank_nearest stage; REGENERATE.md registers only --tool smina)",
  "resolved", "41.6 % with gnina and 41.6 % with smina against the 43.6 % of the confidence order (both re-rank strategies select 126 of 303). Re-ordering still does not help. The gnina re-rank output should be added to the program as a `gnina_rerank_nearest` stage if the author wants a canonical sidecar.")
def tk(base, comp, k): return tks[(tks.baseline == base) & (tks.comparison == comp) & (tks.k == k)].iloc[0]
E("A964.table19_footnote_gains", 964, "tab:appendix-refinement-mcnemar footnote", "3.3:947-964", "value",
  "EquiBind starts at 1.3% and gains more than seventeen points at rank-1, AutoDock Vina starts at 31.4% and gains six, DiffDock starts at 33.7% and gains two; top-15 66.3% / 56.1% both / 27.4%; 24 and 25 recovered complexes",
  {"equibind_rank1": [r1(tk("EquiBind (raw)", "EquiBind (gnina-opt)", 1)["baseline_rate_%"]), r1(tk("EquiBind (raw)", "EquiBind (gnina-opt)", 1)["comparison_rate_%"])],
   "autodock_rank1": [r1(tk("AutoDock Vina", "AutoDock (gnina-opt)", 1)["baseline_rate_%"]), r1(tk("AutoDock Vina", "AutoDock (gnina-opt)", 1)["comparison_rate_%"])],
   "diffdock_rank1": [r1(tk("DiffDock (raw)", "DiffDock + smina", 1)["baseline_rate_%"]), r1(tk("DiffDock (raw)", "DiffDock + smina", 1)["comparison_rate_%"]), r1(tk("DiffDock (raw)", "DiffDock (gnina-opt)", 1)["comparison_rate_%"])],
   "top15_opt": {"autodock": r1(tk("AutoDock Vina", "AutoDock (gnina-opt)", 15)["comparison_rate_%"]), "diffdock_smina": r1(tk("DiffDock (raw)", "DiffDock + smina", 15)["comparison_rate_%"]), "diffdock_gnina": r1(tk("DiffDock (raw)", "DiffDock (gnina-opt)", 15)["comparison_rate_%"]), "equibind": r1(tk("EquiBind (raw)", "EquiBind (gnina-opt)", 15)["comparison_rate_%"])},
   "diffdock_top15_opt_only_wins": [int(tk("DiffDock (raw)", "DiffDock + smina", 15).comparison_only_wins), int(tk("DiffDock (raw)", "DiffDock (gnina-opt)", 15).comparison_only_wins)]},
  S(HUB / "topk_recovery_stats.csv"), "resolved",
  "EquiBind 1.3 -> 19.5 (gains eighteen), AutoDock 43.2 -> 49.8 (gains seven), DiffDock 43.6 -> 45.2 / 45.5 (gains two). Top-15 71.3 % for AutoDock against 60.1 % (smina) and 59.4 % (gnina) for DiffDock and 28.4 % for EquiBind, ordering of raw arms repeated. The two DiffDock refiners no longer reach 'the same' top-15 value: 60.1 against 59.4 by way of 27 and 25 recovered complexes, so the 'indistinguishable ... same 56.1 %' sentence is reworded to 'within one point of each other'.")

# ============================================================================ bands :992-999
def band_vals(path):
    b = pd.read_csv(path); out = {}
    for tool in ("DiffDock", "EquiBind"):
        for crit, lab in (("near_%", "near-native"), ("pbv_%", "PoseBusters-valid"), ("both_%", "near-native and PB-valid")):
            row = []
            for lo, hi in ((1, 10), (11, 20), (21, 30)):
                for opt in ("smina", "gnina"):
                    def pooled(o):
                        g = b[(b.tool == tool) & (b.optimizer == o) & (b["rank"] >= lo) & (b["rank"] <= hi)]
                        return (g[crit] * g.n).sum() / g.n.sum()
                    row.append(round(pooled(opt) - pooled("raw"), 3))
            out[f"{tool}|{lab}"] = row
    return out
bn = band_vals(HUB / "optimization_benefit_by_rank.csv"); bc = band_vals(ROOT / "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/optimization_benefit_by_rank.csv")
for j, key in enumerate(["DiffDock|near-native", "DiffDock|PoseBusters-valid", "DiffDock|near-native and PB-valid", "EquiBind|near-native", "EquiBind|PoseBusters-valid", "EquiBind|near-native and PB-valid"]):
    ln = 992 + j
    E(f"A{ln}.table20_bands.{key}", ln, "tab:appendix-refinement-bands", "3.3:992-999", "table_cell", L(ln).strip(),
      {"columns_1-10_smina_gnina_11-20_smina_gnina_21-30_smina_gnina": bn[key], "print_as": [(f"+{v:.0f}" if "valid" in key and "near" not in key else f"+{v:.1f}") for v in bn[key]], "canonical_recompute_check": bc[key]},
      S(HUB / "optimization_benefit_by_rank.csv") + " pose-count-weighted pooled rate per band (ranks 1-10 / 11-20 / 21-30), optimiser minus raw; recipe reproduces every printed canonical cell", "resolved", "")
E("A999.table20_footnote", 999, "tab:appendix-refinement-bands footnote", "3.3:992-999", "value",
  "near-native gnina bands span 0.1 to 0.7 (DiffDock) and 4.1 to 4.6 (EquiBind); smina bands 0.1 to 0.6 and 3.5 to 4.4; whole-run validity gains +60 / +63 and +48 / +59 against 84.5% / 87.1% from 24.4% and 51.3% / 62.0% from 2.9%",
  {"diffdock_near_gnina_span": [min(bn["DiffDock|near-native"][1::2]), max(bn["DiffDock|near-native"][1::2])], "equibind_near_gnina_span": [min(bn["EquiBind|near-native"][1::2]), max(bn["EquiBind|near-native"][1::2])],
   "diffdock_near_smina_span": [min(bn["DiffDock|near-native"][0::2]), max(bn["DiffDock|near-native"][0::2])], "equibind_near_smina_span": [min(bn["EquiBind|near-native"][0::2]), max(bn["EquiBind|near-native"][0::2])],
   "validity_gains_unchanged": "+60 / +63 and +48 / +59 (validity is convention-free)"},
  "same recipe", "resolved", "DiffDock gnina near bands 0.3 to 1.0, smina 0.4 to 0.8; EquiBind gnina 4.4 to 4.7, smina 3.6 to 4.5. The note to Table 19 quoting '0.1-0.7 and 4.1-4.6' must follow.")

# ============================================================================ cross-tool :1026-1031
ct_rows = [("AutoDock (gnina-opt)", "DiffDock (raw)", 1026), ("AutoDock (gnina-opt)", "DiffDock (gnina-opt)", 1027), ("DiffDock (gnina-opt)", "EquiBind (gnina-opt)", 1028), ("AutoDock (gnina-opt)", "EquiBind (gnina-opt)", 1029)]
for a, b, ln in ct_rows:
    cells = {}
    for k in (1, 5, 10, 15):
        r = tk(a, b, k)
        cells[f"k{k}"] = {"a_pct": r1(r["baseline_rate_%"]), "b_pct": r1(r["comparison_rate_%"]), "a_only": int(r.baseline_only_wins), "b_only": int(r.comparison_only_wins), "p_raw": pfmt(r.mcnemar_p), "p_holm": pfmt(r.p_holm), "sig": r.star}
    E(f"A{ln}.table21.{a}_vs_{b}", ln, "tab:appendix-cross-tool-mcnemar", "3.3:1026-1031", "table_cell", L(ln).strip(), cells, S(HUB / "topk_recovery_stats.csv") + " family between_tool", "resolved", "")
E("A1031.table21_footnote", 1031, "tab:appendix-cross-tool-mcnemar footnote", "3.3:1026-1031", "value",
  "vs DiffDock (raw) Holm 0.694 then 0.0001, 0.0002, 0.0003; vs DiffDock (gnina-opt) 0.709 then 0.0059, 0.0080, 0.029; rank-1 wins 62 / losses 51 (raw) and 60 / 55 (gnina-opt); EquiBind contrasts Holm < 3e-7",
  {"vs_diffdock_raw_holm": [pfmt(tk("AutoDock (gnina-opt)", "DiffDock (raw)", k).p_holm) for k in (1, 5, 10, 15)],
   "vs_diffdock_gnina_holm": [pfmt(tk("AutoDock (gnina-opt)", "DiffDock (gnina-opt)", k).p_holm) for k in (1, 5, 10, 15)],
   "rank1_wins_losses_raw": [int(tk("AutoDock (gnina-opt)", "DiffDock (raw)", 1).baseline_only_wins), int(tk("AutoDock (gnina-opt)", "DiffDock (raw)", 1).comparison_only_wins)],
   "rank1_wins_losses_gnina": [int(tk("AutoDock (gnina-opt)", "DiffDock (gnina-opt)", 1).baseline_only_wins), int(tk("AutoDock (gnina-opt)", "DiffDock (gnina-opt)", 1).comparison_only_wins)],
   "equibind_contrasts_max_holm": pfmt(max(tk(a, "EquiBind (gnina-opt)", k).p_holm for a in ("DiffDock (gnina-opt)", "AutoDock (gnina-opt)") for k in (1, 5, 10, 15)))},
  S(HUB / "topk_recovery_stats.csv"), "resolved", "Holm 0.242 then 1.2e-4, 6.5e-5, 3.8e-5 against DiffDock (raw); 0.294 then 0.0042, 0.0042, 0.0071 against DiffDock (gnina-opt). Rank-1 AutoDock wins 77 / loses 58 (raw) and 72 / 59 (gnina-opt), still a tie. EquiBind contrasts all below 1e-12 ('below 3e-7 throughout' remains true; tighten to 'below 1e-12').")

# ============================================================================ Table 22 :1056-1062
t22 = C["table_22"]
for lab, ln in (("AutoDock Vina + gnina", 1056), ("DiffDock (raw)", 1057), ("DiffDock + gnina", 1058), ("EquiBind (raw)", None), ("EquiBind + gnina", 1060)):
    pass
t22_rows = [("AutoDock Vina + gnina", "AutoDock (gnina-opt)", 1056), ("DiffDock (raw)", "DiffDock (raw)", 1057), ("DiffDock (gnina-opt)", "DiffDock (gnina-opt)", 1058), ("EquiBind (raw)", "EquiBind (raw)", 1059), ("EquiBind (gnina-opt)", "EquiBind (gnina-opt)", 1060)]
vgc_key = {"AutoDock (gnina-opt)": "AutoDock Vina + gnina", "DiffDock (raw)": "DiffDock (raw)", "DiffDock (gnina-opt)": "DiffDock + gnina", "EquiBind (gnina-opt)": "EquiBind + gnina"}
for lab, wkey, ln in t22_rows:
    cells = {f"k{k}": {"cost": w[wkey][str(k)]["gap_k"], "pct": w[wkey][str(k)]["gap_pct"], "ci": w[wkey][str(k)]["gap_ci"]} for k in (1, 5, 10, 15)}
    if wkey in vgc_key:
        for k in (1, 5, 15):
            assert t22[vgc_key[wkey]][str(k)]["cost_complexes"] == cells[f"k{k}"]["cost"], (lab, k)
    changed = lab == "DiffDock (raw)"
    E(f"A{ln}.table22.{lab}", ln, "tab:appendix-accurate-invalid", "3.3:1056-1062", "table_cell", L(ln).strip(), cells,
      S(HUB / "topk_recovery_validity.csv") + " near_k - valid_k with stats_utils.wilson_ci(gap, 303); k 1 / 5 / 15 cross-checked against " + S(HUB / "validity_gate_cost.csv") + "; harness signature table_22 agrees",
      "resolved" if changed else "unchanged", "DiffDock raw k=1 65 (21.5 [17.2, 26.4]), k=5 60 (19.8 [15.7, 24.7]); k=10 59 (19.5) and k=15 57 (18.8) unchanged" if changed else "")
E("A1062.table22_footnote", 1062, "tab:appendix-accurate-invalid footnote", "3.3:1056-1062", "value", "raw AutoDock Vina 0.3 [0.1, 1.8] at the first two depths and 0.0 [0.0, 1.3] at the last two",
  {f"k{k}": {"pct": w["AutoDock Vina"][str(k)]["gap_pct"], "ci": w["AutoDock Vina"][str(k)]["gap_ci"]} for k in (1, 5, 10, 15)}, S(HUB / "topk_recovery_validity.csv"), "unchanged", "footnote wording 'within 2 A of the crystal ligand' -> 'of the nearest deposited copy' (W)")

# ============================================================================ ITT :1074
it = json.load(open(ITT)); nn = it["conventions"]["nearest"]; ii = it["conventions"]["instance"]
E("A1074.itt", 1074, "App. I, Cohort and Exclusions", "3.3:1074 / plan 1.10", "value",
  "none of the five within 2 A at rank-1 and none within its top-15 pool, closest 6.4 A; 111 of 308 and 198 of 308 against 103 and 167; 2.6 pp; 10.2 -> 10.1; all 150 poses pass",
  {"rank1_recovered_dropped": nn["dropped_recovered"]["rank1"], "closest_rank1_A": nn["closest_rank1_A"], "closest_rank1_id": nn["closest_rank1_id"],
   "top15_recovered_dropped": nn["dropped_recovered"]["top15"], "closest_top15_A": nn["closest_top15_A"], "closest_top15_id": nn["closest_top15_id"],
   "all_poses_recovered_dropped": nn["dropped_recovered"]["all_poses"], "closest_all_A": nn["closest_all_A"], "closest_all_id": nn["closest_all_id"], "closest_all_rank": nn["closest_all_rank"],
   "counts_308": nn["counts_308"], "margin_303_pp": nn["margin_pp_303"], "margin_308_pp": nn["margin_pp_308"], "n_poses_scored": it["n_poses_scored"],
   "instance_check": {"counts_308": ii["counts_308"], "margins": ii["margin_pp_308"], "closest_top15_A": ii["closest_top15_A"]}},
  S(ITT), "resolved", "AutoDock Vina + gnina places none of the five within 2 A at rank-1 (closest 6.3 A) and none within its top-15 pool, where the closest pose sits at 2.1 A (7FRX_O88); over all thirty poses 8F4J_PHO reaches 1.5 A of an alternate copy at rank 28. Recovery stays at 149 of 308 rank-1 and 213 of 308 top-15 against 132 and 179; margin 5.6 -> 5.5 at rank-1 and 11.2 -> 11.0 at top-15. The validity clause ('all 150 poses pass') is NOT re-checked by the generator (posebusters_filtered_results.csv holds no rows for the five ids); it stands on the canonical statement.")

# ============================================================================ Table 24 :1172-1181 near-native half + :1188
t24 = C["table_24"]
t24_lines = {"AutoDock (raw)": 1172, "AutoDock + gnina": 1173, "DiffDock (raw)": 1174, "DiffDock + smina": 1175, "DiffDock + gnina": 1176, "EquiBind (raw)": 1177, "EquiBind + smina": 1178, "EquiBind + gnina": 1179}
for arm, ln in t24_lines.items():
    nnh = t24[arm]["near_native_nearest"]; ih = t24[arm]["near_native_instance"]
    E(f"A{ln}.table24.{arm}.near_native_half", ln, "tab:appendix-pb-decomposition-full", "3.3:1172-1181 / W46", "table_cell", L(ln).strip(),
      {"poses": nnh["poses"], "valid_pct": nnh["valid_pct"], "chem_pct": nnh["chem_pct"], "intra_pct": nnh["intra_pct"], "inter_pct": nnh["inter_pct"], "mindist_pct": nnh["mindist_pct"],
       "instance_check": [ih["poses"], ih["valid_pct"], ih["chem_pct"], ih["intra_pct"], ih["inter_pct"], ih["mindist_pct"]]},
      "recompute: posebusters_filtered_results.csv joined to per_pose_metrics_nearest on pose_file, rmsd <= 2, group shares per yaml table_24 (reproduces the printed near-native half on the *_ref_instance twin; harness signature table_24 agrees)", "resolved" if arm != "EquiBind (raw)" else "unchanged", "all-poses half unchanged")
E("A1181.table24_footnote", 1181, "tab:appendix-pb-decomposition-full footnote", "W46", "wording", "within 2 A of the crystal ligand by symmetry-corrected heavy-atom RMSD without superposition", "... of the nearest deposited copy of the crystal ligand ...", "W46", "wording", "")
E("A1188.coupling_strata", 1188, "App. H, Coupling", "3.3:1188-1196 / W47a", "value",
  "DiffDock 14.6% -> 5.4% over 1,715; EquiBind 37.0% -> 6.2% over 498; 98.9% against 62.0% narrows to roughly four points; EquiBind 31 of 498 against 5 of 319 for AutoDock",
  {"diffdock_smina": {"all_inter_pct": t24["DiffDock + smina"]["all"]["inter_pct"], "near_inter_pct": t24["DiffDock + smina"]["near_native_nearest"]["inter_pct"], "near_n": t24["DiffDock + smina"]["near_native_nearest"]["poses"], "near_inter_fail_n": t24["DiffDock + smina"]["near_native_nearest"]["inter_fail_n"]},
   "equibind_gnina": {"all_inter_pct": t24["EquiBind + gnina"]["all"]["inter_pct"], "near_inter_pct": t24["EquiBind + gnina"]["near_native_nearest"]["inter_pct"], "near_n": t24["EquiBind + gnina"]["near_native_nearest"]["poses"], "near_inter_fail_n": t24["EquiBind + gnina"]["near_native_nearest"]["inter_fail_n"]},
   "autodock_gnina": {"near_n": t24["AutoDock + gnina"]["near_native_nearest"]["poses"], "near_inter_fail_n": t24["AutoDock + gnina"]["near_native_nearest"]["inter_fail_n"]},
   "validity_spread_near_native": [t24["AutoDock + gnina"]["near_native_nearest"]["valid_pct"], t24["EquiBind + gnina"]["near_native_nearest"]["valid_pct"]]},
  "same recompute", "resolved", "DiffDock 14.6 % -> 4.6 % over its 2,169 near-native poses; EquiBind 37.0 % -> 6.0 % over 515. The 98.9 against 62.0 spread narrows to 98.0 against 93.8, still 'roughly four percentage points'. EquiBind fails the intermolecular group on 31 of its 515 near-native poses against 5 of 451 for AutoDock.")
vg = t22
E("A1196.validity_cost_headline", 1196, "App. H, what the requirement costs", "3.3:1188-1196 / W47b / A1b", "value",
  "at most four complexes of 303; two at rank-1, four at top-5 and three deeper for AutoDock Vina + gnina, three to four for DiffDock + smina, two to three for EquiBind + gnina; span 0.66 to 1.32 pp; raw DiffDock 47 at rank-1 and 57 to 59 deeper, span 15.5 to 19.5",
  {"autodock_gnina": [vg["AutoDock Vina + gnina"][str(k)]["cost_complexes"] for k in (1, 5, 15, 30)], "diffdock_smina": [vg["DiffDock + smina"][str(k)]["cost_complexes"] for k in (1, 5, 15, 30)],
   "equibind_gnina": [vg["EquiBind + gnina"][str(k)]["cost_complexes"] for k in (1, 5, 15, 30)], "max_three_variants": 5,
   "span_pp": [0.66, r2(100 * 5 / 303)], "diffdock_raw": [vg["DiffDock (raw)"][str(k)]["cost_complexes"] for k in (1, 5, 15, 30)], "diffdock_raw_span_pp": [r1(100 * 55 / 303), r1(100 * 65 / 303)]},
  S(HUB / "validity_gate_cost.csv"), "resolved", "At most five complexes; AutoDock 2 / 4 / 3 / 3 unchanged, DiffDock + smina 5 / 4 / 3 / 3, EquiBind 2 / 3 / 3 / 3; span 0.66 to 1.65 points. Raw DiffDock loses 65 at rank-1 and 55 to 60 at the deeper pools, a span of 18.2 to 21.5 points. A1b adds the pb_rmsd sentence (two decisions, one AutoDock Vina + gnina and one EquiBind + gnina, see A1387).")

# ============================================================================ :1203 :1207
E("A1203.cube_regime_clause", 1203, "App. H, Search Effort and Box Volume", "3.3:1203 / W48 / D11", "wording", "(25 A cube sentence)", "+ W48 sentence: only 8 of the 211 alternate copies have a centroid inside the cube; every box used here contains every copy",
  "plan D11 data facts (8 of 211 alternate centroids; 16 ids by any heavy atom); no rebuilt sidecar", "wording", "The 8 / 211 count is a plan fact without a `_nearest` generator; flagged.")
E("A1207.exh32_to_128_top15", 1207, "App. H, Search Effort and Box Volume", "3.3:1207", "value", "from 148 to 185 complexes of 303", C["line1207"], "recompute (Table 8 raw exh32 / exh128 top-15 double gate); harness signature table_8[autodock_mgltools].d15 = 161, [autodock_mgltools_exh128].d15 = 198", "resolved", "161 to 198.")

# ============================================================================ Table 25 :1247-1250
t25 = {arm: [int(sig_map[f"table_25[{arm}].{f}"]) for f in ("gained", "stay_invalid", "rescued", "net")] for arm in ("AutoDock Vina + gnina", "DiffDock + smina", "EquiBind + gnina")}
for ln, (label, idx) in zip((1247, 1248, 1249, 1250), (("gained", 0), ("stay_invalid", 1), ("rescued", 2), ("net", 3))):
    E(f"A{ln}.table25.{label}", ln, "tab:results-depth-gain", "3.3:1247-1250", "table_cell", L(ln).strip(),
      {"AutoDock Vina + gnina": t25["AutoDock Vina + gnina"][idx], "DiffDock + smina": t25["DiffDock + smina"][idx], "EquiBind + gnina": t25["EquiBind + gnina"][idx]},
      "harness signature table_25 actual (recompute on per_pose_metrics_nearest per yaml table_25 set algebra)", "resolved", "signs as printed: gained +, stay_invalid -, rescued +, net + with ***" + (" (stars: re-check the McNemar on the rebuilt table; nested depths, all three remain p < 0.001)" if label == "net" else ""))

# ============================================================================ H.9 :1258-1260
h9 = C["appendix_h9"]
E("A1258.h9_differences", 1258, "App. H.9, Resolution of the Two-Pipeline Contrast", "3.3:1258-1260 / W49", "value",
  "+2.6 [-4.2, +9.4] rank-1; +13.2 [+5.7, +20.4] top-5; +10.2 [+2.8, +17.5] top-15; +9.9 [+2.5, +17.1] top-30",
  {f"d{d}": {"diff_pp": r1(h9[str(d)]["diff_pp"]), "ci95": [r1(v) for v in h9[str(d)]["ci95"]]} for d in (1, 5, 15, 30)},
  "recompute: stats_utils.newcombe_paired_diff_ci on the double gate, AutoDock minus DiffDock; agrees with " + S(HUB / "bounded_claim.csv") + " (negated) and " + S(REF / "reference_convention_margins.csv"), "resolved",
  "+5.6 [-1.7, +12.8]; +13.2 [+5.7, +20.4] (identical); +11.2 [+3.9, +18.4]; +10.2 [+3.0, +17.3].")
E("A1260.h9_rank1_resolution", 1260, "App. H.9", "3.3:1258-1260 / W49", "value",
  "both 51, neither 140, discordant 112 (37.0%); MDD 9.8; power 0.12; roughly 4,200; 428 -> 8.2; p 7.6e-4, 0.009, 0.011",
  {"both": h9["1"]["both"], "neither": h9["1"]["neither"], "discordant": h9["1"]["discordant"], "discordant_pct": h9["1"]["discordant_pct"], "mde_pp": h9["1"]["mde_pp"], "observed_power": r2(h9["1"]["observed_power"]),
   "n_needed": h9["1"]["n_needed"], "n_needed_print": "roughly 1,100", "mde_at_428_pp": h9["1"]["mde_at_428_pp"], "mcnemar_p_5_15_30": [pfmt(h9[str(d)]["mcnemar_p"]) for d in (5, 15, 30)]},
  "recompute: stats_utils.mcnemar_power (MDD at 80 % power on the observed discordance; scaled by sqrt(303/428) for the 428-entry release), mcnemar_exact; agrees with bounded_claim.csv (power 0.322, mde 10.5, n_needed 1062)", "resolved",
  "Both 76, neither 98, discordant 129 (42.6 %); MDD 10.5 points; power 0.32; roughly 1,100 paired complexes; 428 entries leave the MDD near 8.8; top-5 / top-15 / top-30 separate at p = 7.6e-4, 0.0036 and 0.0075.")
E("A1260.h9_equivalence", 1260, "App. H.9", "3.3:1258-1260 / W49", "value",
  "only rank-1 is equivalent, within 10 points, smallest passing value 8.4; top-5 / 15 / 30 smallest passing 19.4, 16.4, 16.0; 90% interval -3.1 to +8.3, +6.9 to +19.3, +4.0 to +16.3, +3.7 to +16.0",
  {f"d{d}": {"ci90": [r1(v) for v in h9[str(d)]["ci90"]], "smallest_passing_margin": h9[str(d)]["smallest_passing_margin_pp"], "equivalent": h9[str(d)]["equiv"]} for d in (1, 5, 15, 30)},
  "recompute: stats_utils.tost_paired_proportions (90 % Newcombe interval) on the grid 5 / 10 / 12 / 15; margin = larger absolute 90 % bound rounded UP to one decimal (harness rule); harness signature prose[h9] agrees", "resolved",
  "Only rank-1 is equivalent, within 12 points and no longer within 10, smallest passing value 11.7; deeper depths 19.4, 17.3 and 16.2; 90 % intervals -0.5 to +11.7, +6.9 to +19.3, +5.1 to +17.2 and +4.2 to +16.2.")
E("A1265.threshold_leadin", 1265, "App. H, Threshold-Resolved Recovery", "3.3:1265 / W50", "wording", "for both the as-placed and the best-fit (Kabsch) criterion.", "..., each taken to the deposited copy of the ligand nearest to the pose.", "W50", "wording", "")

# ============================================================================ Table 26 :1303-1340
t26 = sig[sig.key.str.startswith("table_26") & sig.key.str.contains(r"\.t[\d.]+$")].copy()
t26["metric_tool"] = t26.key.str.extract(r"table_26\[(.*)\]\.")[0]
t26["depth"] = t26.key.str.extract(r"\.d(\d+)\.")[0].astype(int)
t26["thr"] = t26.key.str.extract(r"\.t([\d.]+)$")[0].astype(float)
line_map = {}
ln = 1303
for metric in ("RMSD", "Kabsch RMSD"):
    for tool in ("AutoDock Vina + gnina", "DiffDock + smina", "EquiBind + gnina"):
        for d in (1, 15, 30):
            line_map[(f"{metric}/{tool}", d)] = ln; ln += 2
for (mt, d), lnum in line_map.items():
    g = t26[(t26.metric_tool == mt) & (t26.depth == d)].sort_values("thr")
    cur = [float(v) for v in g.expected]; new = [float(v) for v in g.actual]
    moved = sum(abs(a - b) > 0.05 for a, b in zip(cur, new))
    E(f"A{lnum}.table26.{mt.replace(' ', '_')}.d{d}", lnum, "tab:results-near-native-form", "3.3:1303-1340 / plan 1.4", "table_cell",
      {"thresholds": list(g.thr), "pct": cur}, {"thresholds": list(g.thr), "pct": new},
      "harness signature table_26 actual (recompute: pb_valid & rmsd <= t / pb_valid & bestfit_rmsd <= t, existence over depth, per_pose_metrics_nearest); cross-check " + S(HUB / "topn_within_thresholds_pbvalid.csv") + " / topn_within_thresholds_kabsch_pbvalid.csv and the two *_report.txt",
      "resolved" if moved else "unchanged", f"{moved} of 10 cells move; the line spans two source lines ({lnum}-{lnum+1}) in the .tex")

# ============================================================================ Table 27 :1377-1387
pt = pd.read_csv(HUB / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__stats_per_tool_depth.csv")
cols = ["n_poses", "n_valid_complexes", "coverage_pct", "inplace_median", "inplace_q1", "inplace_q3", "form_median", "form_q1", "form_q3", "median_r", "pct_placement_limited", "pct_mixed", "pct_form_limited"]
ln = 1377
for tool, lab in (("AutoDock", "AutoDock Vina + gnina"), ("DiffDock", "DiffDock + smina"), ("EquiBind", "EquiBind + gnina")):
    for d in (1, 5, 15):
        r = pt[(pt.tool == tool) & (pt.depth == d)].iloc[0]
        new = [int(r.n_poses), int(r.n_valid_complexes), r1(r.coverage_pct), round(r.inplace_median, 3), round(r.inplace_q1, 3), round(r.inplace_q3, 3), round(r.form_median, 3), round(r.form_q1, 3), round(r.form_q3, 3), round(r.median_r, 3), r1(r.pct_placement_limited), r1(r.pct_mixed), r1(r.pct_form_limited)]
        E(f"A{ln}.table27.{tool}.d{d}", ln, "tab:results-placement-form", "3.3:1377-1387 / plan 1.6", "table_cell", L(ln).strip(), dict(zip(cols, new)),
          S(HUB / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__stats_per_tool_depth.csv") + f" (tool {tool}, depth {d}); harness signature table_27 agrees", "resolved",
          "EquiBind rank-1 in-place median prints 3.239 (the shipped 2.977 was already truncated to 2.9 in the body, A3)" if tool == "EquiBind" and d == 1 else "")
        ln += 1
f27 = C["table27_footnote"]
cap = (HUB / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__caption.txt").read_text()
E("A1387.table27_footnote", 1387, "tab:results-placement-form footnote", "3.3:1377-1387 / W51 / plan 1.6", "value",
  "1.490 against 1.543 (AutoDock Vina + gnina), 1.467 against 1.566 (DiffDock + smina); two of the 909 decisions (one AutoDock, one EquiBind, none DiffDock); 37.0% / 18.5% / unchanged 34.0%",
  {"autodock": {"median_pb_rmsd": f27["AutoDock Vina + gnina"]["median_pb_rmsd"], "median_own_rmsd": f27["AutoDock Vina + gnina"]["median_rmsd"], "n": f27["AutoDock Vina + gnina"]["near_site_rank1_n"]},
   "diffdock": {"median_pb_rmsd": f27["DiffDock + smina"]["median_pb_rmsd"], "median_own_rmsd": f27["DiffDock + smina"]["median_rmsd"], "n": f27["DiffDock + smina"]["near_site_rank1_n"]},
   "equibind": {"median_pb_rmsd": f27["EquiBind + gnina"]["median_pb_rmsd"], "median_own_rmsd": f27["EquiBind + gnina"]["median_rmsd"], "n": f27["EquiBind + gnina"]["near_site_rank1_n"]},
   "decisions_moved": {k: v["decisions_moved_missed_to_recovered"] + v["decisions_moved_recovered_to_missed"] for k, v in f27.items() if isinstance(v, dict)}, "total_decisions": f27["total_decisions"],
   "rank1_recovery_pb_gate_pct": {k: v["rank1_recovery_pb_gate_pct"] for k, v in f27.items() if isinstance(v, dict)}, "pb_rmsd_never_smaller_violations": f27["rank1_pbrmsd_gt_rmsd_rows"],
   "trim_note_from_caption": [l for l in cap.splitlines() if "dropped" in l or "beyond the frame" in l]},
  "recompute: rank-1 PB-valid poses with centroid_dist <= 8 A to the nearest copy, medians of pb_rmsd and rmsd (the 1.376 reproduces the table's in-place median); gate substitution pb_rmsd <= 2 vs rmsd <= 2 at rank-1 with pb_valid", "resolved",
  "1.376 against 1.491 A for AutoDock Vina + gnina and 1.532 against 1.585 A for DiffDock + smina (EquiBind 3.239 against 3.276). Still two of the 909 decisions, one AutoDock Vina + gnina and one EquiBind + gnina, none DiffDock, both from missed to recovered. Rank-1 recovery would then read 49.5 % and 19.1 % against an unchanged 43.6 %. pb_rmsd > rmsd on 0 rows ('never smaller' holds). Trim wording per W51.")

# ============================================================================ cluster statistics :1394 and Fig 38 :1399
hj = json.load(open(CLU / "crystal_cluster_homogeneity_stats.json"))
cq = (CLU / "cluster_quality_metrics_stats.txt").read_text()
omni_q = float(re.search(r"Cochran Q=([\d.]+)", cq).group(1))
A = hj["A_pose_counts"]; Rr = hj["A_reaches"]; B = hj["B_consensus"]; Cc = hj["C_tightness"]; D = hj["D_radius_trend"]
def ph(pairs, name): return [p for p in pairs if p["pair"] == name][0]
E("A1394.cluster_statistics", 1394, "App. H, Cluster-Composition and Quality Statistics", "3.3:1394-1399 / D5 / Phase 0.3 loader fix", "value",
  "Friedman chi2 33.6, W 0.055; Cochran Q 139.6; pose count DiffDock > AutoDock rb -0.28 (Holm p 2.7e-4, 248 pairs), > EquiBind +0.34 (p < 1e-4); AutoDock vs EquiBind +0.01, p 0.912; reach AutoDock vs DiffDock Holm p 0.470 on 38 / 31; both > EquiBind p < 1e-20; consensus 55 single-tool vs 51.8 expected, 105 two-tool vs 149.9; KW H 46.4, eps2 0.075; paired 110 complexes Friedman 31.1; JT z 6.07; ranking Cochran Q 174.0",
  {"pose_count_friedman": {"chi2": A["chi2"], "kendall_w": A["kendall_w"], "p": A["p"]},
   "reach_cochran_q": {"Q": Rr["Q"], "p": Rr["p"], "rates": Rr["reach_rate"]},
   "pose_count_pairs": {"AutoDock_vs_DiffDock": {"rank_biserial": round(ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs DiffDock*")["rank_biserial"], 3), "p_holm": pfmt(ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs DiffDock*")["p_holm"]), "n_pairs": ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs DiffDock*")["n_pairs"]},
                        "DiffDock_vs_EquiBind": {"rank_biserial": round(ph(A["posthoc_wilcoxon_holm"], "DiffDock* vs EquiBind*")["rank_biserial"], 3), "p_holm": pfmt(ph(A["posthoc_wilcoxon_holm"], "DiffDock* vs EquiBind*")["p_holm"])},
                        "AutoDock_vs_EquiBind": {"rank_biserial": round(ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs EquiBind*")["rank_biserial"], 3), "p_holm": pfmt(ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs EquiBind*")["p_holm"]), "n_pairs": ph(A["posthoc_wilcoxon_holm"], "AutoDock* vs EquiBind*")["n_pairs"]}},
   "reach_pairs": {"AutoDock_vs_DiffDock": {"p_holm": round(ph(Rr["posthoc_mcnemar_holm"], "AutoDock* vs DiffDock*")["p_holm"], 3), "disc": [ph(Rr["posthoc_mcnemar_holm"], "AutoDock* vs DiffDock*")["disc_i_only"], ph(Rr["posthoc_mcnemar_holm"], "AutoDock* vs DiffDock*")["disc_j_only"]]},
                   "AutoDock_vs_EquiBind_p_holm": pfmt(ph(Rr["posthoc_mcnemar_holm"], "AutoDock* vs EquiBind*")["p_holm"]), "DiffDock_vs_EquiBind_p_holm": pfmt(ph(Rr["posthoc_mcnemar_holm"], "DiffDock* vs EquiBind*")["p_holm"])},
   "consensus": {"observed_0_1_2_3": B["observed"], "expected_0_1_2_3": B["expected"], "single_tool_obs_exp": [B["observed"][1], B["expected"][1]], "two_tool_obs_exp": [B["observed"][2], B["expected"][2]], "three_tool_obs_exp": [B["observed"][3], B["expected"][3]]},
   "tightness": {"H": Cc["H"], "p": Cc["p"], "epsilon_sq": Cc["epsilon_sq"], "n": Cc["n"], "paired_subset": Cc["sensitivity_friedman_common_subset"]},
   "radius_trend": {"z": D["z"], "p": D["p"]}, "ranking_rule_cochran_q": omni_q},
  S(CLU / "crystal_cluster_homogeneity_stats.json") + " (+ .txt) and " + S(CLU / "cluster_quality_metrics_stats.txt") + " (any-copy rule, D5, per-pose table convention nearest, after the Phase 0.3 hydrogen-centroid loader fix)", "resolved",
  "Friedman chi2 36.6, W 0.06; Cochran Q 159.5. Pose count: DiffDock > AutoDock rb -0.31 (Holm p < 1e-4, 249 pairs), DiffDock > EquiBind +0.52 (p < 1e-4). CHANGE OF VERDICT: AutoDock and EquiBind now DO separate on pose count (rb +0.21, Holm p 0.004, 250 pairs), so 'AutoDock and EquiBind do not separate on pose count' must be rewritten. Reach: AutoDock vs DiffDock Holm p 0.161 on 31 / 20 discordant; both exceed EquiBind at p < 1e-23. Consensus: 37 single-tool against 33.7 expected and 120 two-tool against 147.4 (the all-three cell 137 against 119.9 is the enrichment). KW H 49.2, eps2 0.074; paired 119 complexes Friedman 37.9; JT z 6.31; ranking-rule Cochran Q 149.6.")
E("A1399.fig38_caption", 1399, "fig:results-r6 caption", "3.3:1394-1399 / W52 / D5", "wording", "4 A crystal-site recovery threshold", "... to any deposited copy of the ligand (D5 wording: for every copy, the cluster nearest to that copy counts when within 4 A; a cluster nearest to two copies counts once)", "W52; " + S(CLU / "cluster_quality_metrics_stats.txt") + " header", "wording", "Regenerate image8 from " + S(CLU / "cluster_quality_metrics.png"))

# ============================================================================ :1407-1415 interaction rank concordance
tau04 = (PMR / "04c_native_recovery_by_rank.txt").read_text(); tau10 = (PMR / "10b_contact_decomposition_by_rank.txt").read_text()
def grab(txt, tool, metric):
    blk = txt.split(tool, 1)[1]
    m = re.search(metric + r": median tau=([+-]?[\d.]+), rank-biserial=([+-]?[\d.]+), p=([^,]+), BH q=([^ ]+) ?(\S*)\s+\(n=(\d+)", blk)
    return {"tau": float(m.group(1)), "rank_biserial": float(m.group(2)), "p": m.group(3).strip(), "q": m.group(4).strip(), "n": int(m.group(6))}
E("A1407.kendall_tau_recovery", 1407, "App. H, Rank-Quality Concordance", "3.3:1407-1415", "value",
  "AutoDock tau -0.20 / -0.32 / -0.20 (P / R / F1) all q < 0.001; EquiBind recall and F1 -0.20 both significant, precision marginal q = 0.022; DiffDock tau ~ 0, none significant",
  {"AutoDock": {m: grab(tau04, "AutoDock Vina exh128 (gnina-opt):", m) for m in ("Precision", "Recall", "F1")},
   "DiffDock": {m: grab(tau04, "DiffDock (smina-opt):", m) for m in ("Precision", "Recall", "F1")},
   "EquiBind": {m: grab(tau04, "EquiBind (unguided, gnina-opt):", m) for m in ("Precision", "Recall", "F1")}},
  S(PMR / "04c_native_recovery_by_rank.txt") + " (per-copy crystal fingerprints, nearest-copy join)", "resolved",
  "AutoDock tau -0.32 / -0.40 / -0.40, all q < 0.001. EquiBind recall and F1 -0.20 (q < 0.001 and 7.7e-4), precision -0.11 (q 0.008). CHANGE OF VERDICT: DiffDock's confidence NOW shows a detected association on all three (tau -0.11 / -0.18 / -0.11, q 0.008 / 5.1e-4 / 0.002), so 'no detected association on any of the three ... absence of detected ordering' must be rewritten to a weak but detected ordering.")
cd = pd.read_csv(PMR / "10b_contact_decomposition_by_rank.csv").set_index(["method", "pose_rank"])
E("A1415.contact_decomposition", 1415, "App. H, Rank-Quality Concordance, Fig 39 prose", "3.3:1407-1415 / W53", "value",
  "AutoDock 15.3 spurious against DiffDock 14.1 at rank 1; AutoDock 17.6 by rank 5, DiffDock flat; EquiBind 16.7 matched against 18.9 missed, spurious 12.0 to 12.6; AutoDock tau -0.32 matched, +0.32 missed, +0.11 spurious all sig; EquiBind +-0.26; DiffDock flat",
  {"rank1_spurious": {"autodock": round(cd.loc[("autodock_mgltools_exh128_gnina", 1)].fp, 1), "diffdock": round(cd.loc[("diffdock_smina", 1)].fp, 1)},
   "autodock_rank5_spurious": round(cd.loc[("autodock_mgltools_exh128_gnina", 5)].fp, 1), "diffdock_spurious_r1_r5": [round(cd.loc[("diffdock_smina", k)].fp, 1) for k in (1, 5)],
   "equibind_rank1_matched_missed": [round(cd.loc[("equibind_unguided_gnina", 1)].tp, 1), round(cd.loc[("equibind_unguided_gnina", 1)].fn, 1)],
   "equibind_spurious_range": [round(min(cd.loc[("equibind_unguided_gnina", k)].fp for k in range(1, 6)), 1), round(max(cd.loc[("equibind_unguided_gnina", k)].fp for k in range(1, 6)), 1)],
   "tau": {"AutoDock": {m: grab(tau10, "AutoDock Vina exh128 (gnina-opt):", m) for m in ("matched", "missed", "spurious")},
           "DiffDock": {m: grab(tau10, "DiffDock (smina-opt):", m) for m in ("matched", "missed", "spurious")},
           "EquiBind": {m: grab(tau10, "EquiBind (unguided, gnina-opt):", m) for m in ("matched", "missed", "spurious")}}},
  S(PMR / "10b_contact_decomposition_by_rank.csv") + " and .txt; harness signature prose.fig39 agrees", "resolved",
  "At rank 1 AutoDock invents 9.1 spurious contacts against 7.8 for DiffDock; by rank 5 AutoDock's climbs to 13.5 while DiffDock's stays near 8.6. EquiBind recovers 17.7 matched against 17.8 missed at rank 1 and invents the fewest, 11.2 to 11.6. AutoDock tau -0.40 matched, +0.40 missed, +0.32 spurious, all significant; EquiBind +-0.20 on matched / missed. DiffDock matched -0.12 and missed +0.12 are now significant (q 0.002 / 0.003) with spurious flat, so 'DiffDock's counts are flat' becomes 'DiffDock's spurious count is flat and its matched and missed counts move weakly'. Regenerate image11.")

# ============================================================================ D9 sensitivity table (new material beside :1074)
tex_tab = (REF / "reference_convention_sensitivity.tex").read_text()
rcsv = pd.read_csv(REF / "reference_convention_sensitivity.csv")
rows = rcsv[rcsv.stratum != "all"][["arm_label", "depth", "stratum", "n_complexes", "recovered_nearest", "recovered_instance", "gain", "recovered_nearest_pct", "recovered_instance_pct", "form_complexes_nearest", "form_complexes_instance", "combined_complexes_nearest", "combined_complexes_instance"]]
E("A1074.D9_sensitivity_table_verbatim", 1074, "App. I, beside the intention-to-treat paragraph (label cohort-and-exclusions)", "D9 / B1b / B6 new material", "verbatim", "(no current text)", tex_tab,
  S(REF / "reference_convention_sensitivity.tex") + " (generated 2026-09-08T11:47:58 by reference_convention_table.py; paste verbatim; house-style check: the caption uses no mid-sentence colon or semicolon)", "resolved",
  "The single-instance convention is the sensitivity arm. Single-copy stratum rank-1: 88 / 80 / 53 of 165 = 53.3 / 48.5 / 32.1 % on the PB-valid and <= 2 A gate, identical under both conventions. The table's 'Recovery' block at d = 5 (200 / 184 AutoDock) differs from the M2 preview; the .tex is the source. The caption's \\AA{} macro must resolve in the thesis preamble (the body uses \\angstrom{}); adapt the macro, not the numbers.")
E("A1074.D9_single_multi_copy_rows", 1074, "App. I, D9 table companion rows", "D9", "table_block", "(no current text)",
  rows.to_dict(orient="records"), S(REF / "reference_convention_sensitivity.csv") + " (stratum in single_copy, multi_copy)", "resolved",
  "Multi-copy stratum (138): AutoDock rank-1 61 vs 23, DiffDock 52 vs 23, EquiBind 4 vs 2 under nearest vs instance; every gain sits in this stratum.")

json.dump({"meta": {"file": FILE, "copy": str(TEX), "built": pd.Timestamp.now().isoformat(timespec="seconds"),
                    "sources_root": str(ROOT), "computed_values": S(PROG / "computed_appendix_values.json"),
                    "signature": S(SIG), "rule": "every `new` is a rebuilt `_nearest` output or a stats_utils recompute on per_pose_metrics_nearest; probe previews are never used (D16)"},
           "n_entries": len(entries), "entries": entries},
          open(PROG / "values_registry_appendix.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
print("entries:", len(entries))
from collections import Counter
print(Counter(e["status"] for e in entries))
print([e["id"] for e in entries if e["status"] == "unresolved"])
