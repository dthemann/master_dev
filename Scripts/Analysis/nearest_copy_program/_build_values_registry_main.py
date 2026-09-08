"""Builds values_registry_main.json (plan Sections 3.1 / 3.2 + change-list Part A/B items for
Thesis_short.tex and body_main_short.tex). Every number is copied from a rebuilt _nearest sidecar or
recomputed from pose_comparison_report_nearest/per_pose_metrics.csv with stats_utils; the source of
each is recorded beside it. No thesis file is touched."""
import json
from pathlib import Path

HUB = "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest/"
CL = "posebusters_results/cluster_crystal_pocket_matched_equibind/autodock_mgltools_exh128_gnina__diffdock_smina_allposes_nearest/"
PM = "pandamap_results/benchmark_matched_equibind_nearest/report/"
EF = "posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged_nearest/"
EX = "posebusters_results/autodock_exhaustiveness_returns_nearest/exh_returns_report.txt"
MS = "posebusters_results/metal_stratum_nearest/"
RC = "posebusters_results/reference_convention_nearest/"
SC = "posebusters_results/selection_check_nearest/selection_check.txt"
SIG = "Scripts/Analysis/nearest_copy_program/harness_signatures/harness_signature_A_oldvalues_newtree.csv (column actual)"
PPM = HUB + "per_pose_metrics.csv"
RECOMP = ("recomputed 2026-09-08 from " + PPM + " with Scripts/Analysis/stats_utils.py; rank-1 = rank column "
          "(AutoDock gnina arm rank == optimized_rank, DiffDock confidence rank), EquiBind rank-1 = lowest gnina_affinity "
          "(rank is the 999 sentinel); gates: near-native rmsd <= 2, recovery = pb_valid AND rmsd <= 2, "
          "form = bestfit_rmsd <= 1 AND rmsd < 1000, triple = all three; existence over the top-d poses; n = 303 analysed ids")
TOPK = HUB + "topk_recovery_validity_gnina_arm.csv"
DEPTHS = HUB + "18_topn_within_thresholds_pbvalid_depths_report.txt"
KAB = HUB + "18_topn_within_thresholds_kabsch_pbvalid_depths_report.txt"
BC = HUB + "bounded_claim.csv (k = 1 row)"
F20 = HUB + "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__rank1_top5_top15__"
EFQ = EF + "effort_by_quality_stats.json"
EFS = EF + "effort_stats.json"
EFSUM = EF + "effort_summary.csv"
EFPC = EF + "per_complex_effort.csv"

R = {}

def add(key, current, new_values, source, wording_change, note=""):
    R[key] = {"current": current, "new_values": new_values, "source": source,
              "wording_change": bool(wording_change), "note": note}

# ------------------------------------------------------------------ headline block (shared)
HEAD = {
    "rank1_pct": {"AutoDock": 49.2, "DiffDock": 43.6, "EquiBind": 18.8},
    "rank1_k": {"AutoDock": 149, "DiffDock": 132, "EquiBind": 57},
    "top15_pct": {"AutoDock": 70.3, "DiffDock": 59.1, "EquiBind": 27.4},
    "top15_k": {"AutoDock": 213, "DiffDock": 179, "EquiBind": 83},
    "top30_pct": {"AutoDock": 70.6, "DiffDock": 60.4, "EquiBind": 27.4},
    "top30_k": {"AutoDock": 214, "DiffDock": 183, "EquiBind": 83},
    "rank1_margin_pp": 5.6, "rank1_newcombe95_lo": -1.7, "rank1_newcombe95_hi": 12.8,
    "rank1_mcnemar_p": 0.159, "rank1_discordant_AD_only": 73, "rank1_discordant_DD_only": 56,
    "rank1_failures_no_qualifying_pose_pct": {"AutoDock": 57.8, "DiffDock": 70.2, "EquiBind": 89.4},
    "six_in_ten_wording": "six to seven in ten (57.8 % and 70.2 % for the leading pair)",
    "equibind_success_complexes_top30": 83,
}
HEAD_SRC = (TOPK + " (pbvalid_near_native_k / pbvalid_near_native_pct, between_tools_mcnemar rows); " + RC +
            "reference_convention_margins.csv (Newcombe, McNemar); failure shares = (303 - rank1_k) vs (303 - top30_k) from the same csv")

# ================================================================== 3.1 Thesis_short.tex
add("T:61", "primary endpoint combining PoseBusters validity with a heavy-atom root-mean-square deviation of at most 2 A",
    {"optional_clause": "to the nearest deposited copy of the ligand"}, "plan D11 / change list B4a (wording only)", True,
    "W optional (plan 3.1 row 61). No number.")
add("T:63", "36.6 / 34.0 / 18.2 %; 2.6 points; -4.2 to +9.4; six in ten; 65.3 / 55.1 %; 80 of 303; 'corrected local geometry but recovered no binding mode ...'; 'most expensive per qualifying pose'",
    dict(HEAD, **{
        "abstract_swaps": {"36.6->": 49.2, "34.0->": 43.6, "18.2->": 18.8, "2.6->": 5.6, "-4.2->": -1.7, "+9.4->": 12.8, "65.3->": 70.3, "55.1->": 59.1, "80->": 83},
        "geometry_clause_check": {"diffdock_raw_rank1_near_native": 132, "diffdock_smina_rank1_near_native": 137, "verdict": "clause holds (refinement adds 5 rank-1 near-native complexes, no grossly misplaced mode recovered; oracle raw 161 -> minimised 186 over the pool)"},
        "cost_ordering_check": {"median_charged_s_per_qualifying_pose": {"AutoDock": 91.8, "DiffDock": 39.0, "EquiBind": 2.4}, "verdict": "AutoDock stays most expensive per qualifying pose; the clause holds, no W"},
    }),
    HEAD_SRC + "; geometry clause: " + HUB + "validity_gate_cost.csv near_native k=1 (DiffDock raw 132, DiffDock + smina 137) and PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv; cost: " + EFQ + " tiers.near2.per_method.median_s",
    True, "Wording: 'six in ten' -> 'six to seven in ten'. Cost clause unchanged (ordering holds).")
add("T:80", "Kurzfassung endpoint definition", {}, "plan 3.1 row 80", False, "no change")
add("T:89-90", "Sie korrigierte lokale Geometrie, nicht jedoch fehlende oder grob fehlplatzierte Bindungsmodi.",
    {"candidate_cut": "drop 'jedoch' (-7 chars) plus B4d second cut on :96-97 (-8) to offset the +12 chars of :94-95"},
    "change list B4d; clause verified as at T:63 (132 -> 137)", True,
    "Compensating cut for the full Kurzfassung page (memory: zero slack). Raster-check after editing.")
add("T:91/92/94/95/97/100", "36,6 / 34,0 / 18,2; 2,6; -4,2 bis +9,4; sechs von zehn; 65,3 / 55,1; 80",
    {"91": "49,2 / 43,6 / 18,8", "92": "5,6", "94": "-1,7 bis +12,8", "95": "sechs bis sieben von zehn", "97": "70,3 / 59,1", "100": "83",
     "length_delta_chars": "+12 on lines 94-95 (change list B4c)"}, HEAD_SRC, True,
    "German decimal commas. Same values as T:63.")
add("T:98-99", "In einheitlicher Belegungsrechnung war der AutoDock-Ablauf am teuersten pro qualifizierender Pose",
    {"median_charged_s_per_qualifying_pose": {"AutoDock": 91.8, "DiffDock": 39.0, "EquiBind": 2.4}, "verdict": "ordering holds, no wording change"},
    EFQ + " tiers.near2.per_method.median_s", False, "W-if row: the ordering did not flip.")

# ================================================================== 3.2 body_main_short.tex
add("M:88", "78 of the 303 analysed complexes place a metal within 5 A of the crystal ligand",
    {"n_metal_adjacent": 78, "n_metal_free": 225, "rule": "min heavy-atom distance from the reference instance (record 0 of _ligands.sdf) to any atom of a residue containing a metal element <= 5 A"},
    MS + "metal_stratum_summary.json (primary.n_metal_adjacent, rule)", True,
    "D4 wording (change list W1): 'a residue carrying a metal within 5 A of the deposited reference instance of the ligand'. 78 unchanged.")
add("M:109-115", "Methods definition (near-nativeness = RMSD to the crystallographic ligand); :113 'adds three recovered complexes'; :115 headline criterion",
    {"113_three_to_two": 2, "single_copy_complexes_303": 165, "multi_copy_303": 138, "multi_copy_308": 143},
    "A1 (Part A, valid today, verified on the canonical table: 7PGX_FMN AutoDock gnina, 6XBO_5MC EquiBind gnina, none DiffDock); counts: " + RC + "reference_convention_summary.json (single_copy_303, multi_copy_303, count_308.multi_copy)",
    True, "B1a/B1b/B1c text in the change list. The 'never smaller than a fully symmetrised RMSD' sentence must be re-verified on the nearest table by the editor (pb_rmsd <= rmsd on every row); the 'two' count is convention-invariant (plan 0.3).")
add("M:119", "clustering poses across tools and against the crystallographic pocket", {}, "change list W3 (D5)", True, "W only")
add("M:124-133", "The value is obtained from the PoseBusters implementation (:130)", {}, "change list W4 (D2)", True, "W only at :130; :133 decomposition caveat stays.")
add("M:140", "footnote: reference-dependent checks reduce DiffDock recovery in Table 2 by one at rank-1 and by two at top-15 and top-30; AutoDock and EquiBind unchanged",
    {}, "NO SIDECAR: no PoseBusters run with a co-crystallised mol_true exists under the nearest convention (run_posebusters.py sets no mol_true path); plan 3.2 row 140 = recompute + name the file passed",
    False, "UNRESOLVED. Either re-run the reference-dependent battery with _ligands.sdf as mol_true or mark the footnote as measured on the single reference instance and name _ligand.sdf.")
add("M:158", "4 A centroid threshold; cluster nearest the crystal not automatically recovered", {"pooled_oracle_any_copy_D5_pct": 99.7, "pooled_oracle_hub_D2_pct": 99.3},
    "change list W5; " + CL + "ensembles_summary.csv (AD+DD+EB oracle_hit_rate 0.9967)", True, "W (D5 wording as given in the task).")
add("M:165", "each tool is compared with the crystal-ligand fingerprint", {}, "change list W6 (union rule, plan 1.8)", True, "W only")
add("M:170", "within 2 A of the crystal ligand (cost definition)", {"qualifying_poses": {"AutoDock": 442, "DiffDock": 2068, "EquiBind": 483}}, "change list W7; " + EFSUM + " poses_near2", True, "W only; counts for context.")
add("M:191", "Variants chosen on PB-valid, within 2 A in place and within 1 A Kabsch", {}, "plan row 191-199", False, "definition sentence, no number")
add("M:193", "164 complexes against 160; two DiffDock refiners each recover 138; agree on 133, differ on five in each direction (p = 1.000); 87.1 vs 84.5 %; 62.0 vs 51.3 %; 55 versus 48",
    {"autodock_raw_vs_gnina_pooled_triple": [178, 175], "diffdock_smina_pooled_triple": 148, "diffdock_gnina_pooled_triple": 149, "agree_both": 145,
     "smina_only": 3, "gnina_only": 4, "mcnemar_p_exact": 1.0, "validity_pct_unchanged": [87.1, 84.5, 62.0, 51.3], "equibind_gnina_vs_smina_pooled_triple": [57, 50]},
    HUB + "pose_validity_cascade.csv (triple_complexes); agreement 145 / 3 / 4 and p: " + RECOMP, True,
    "W8a: 'each recover 138' must go (148 vs 149); 'differ on three and four'.")
add("M:195", "AutoDock exh128 + gnina reaches 160; DiffDock smina 133 vs gnina 132; EquiBind gnina 55 vs smina 46; no guided arm exceeds 18 at any depth; 198 complexes at top-15 rather than 160",
    {"autodock_top15_triple": 175, "diffdock_smina_top15_triple": 144, "diffdock_gnina_top15_triple": 143, "equibind_gnina_top15_triple": 57, "equibind_smina_top15_triple": 48,
     "guided_arm_ceiling_top15_triple": 23, "guided_arm_ceiling_arm": "equibind_p2rank_gnina (fpocket gnina 20)", "autodock_top15_double_gate": 213},
    SC + " (top-15 TRIPLE gate dictionaries); double gate 213: " + TOPK, True,
    "W8b. The '18' ceiling becomes 23 (guided arms saturate: top-15 == pooled for all six guided arms, " + RECOMP.split(';')[0] + ").")
add("M:197", "133 to 132; 834 to 828 poses; disagree on nine, four gnina-only; gnina higher near-native (176 vs 173); smina higher form (237 vs 235); footnote 32.7 / 33.7 / 32.3 %",
    {"top15_triple": [144, 143], "top15_triple_poses": [1031, 1022], "discordant_total": 7, "gnina_only": 3, "smina_only": 4, "mcnemar_p": 1.0,
     "pooled_near_native_complexes": [186, 186], "pooled_form_complexes": [236, 235],
     "footnote": {"smina_rerank_raw_coords_pct": 41.6, "diffdock_confidence_same_coords_pct": 43.6, "gnina_rerank_raw_coords_pct": 41.6,
                  "basis": "near-native RMSD <= 2, no validity gate, 303 complexes, raw coordinates (strategies C_rerank vs A_native)"}},
    RC + "reference_convention_sensitivity.csv (combined_complexes/poses at depth 15); discordant: " + RECOMP + "; pooled 186/186 and 236/235: " + SIG +
    " table_1 near_complexes / form_complexes; footnote: PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv (C_rerank 41.58, A_native 43.56) and Scripts/Analysis/nearest_copy_program/gnina_rerank_nearest/selection_strategy_summary_gnina.csv (C_rerank 41.58; run 2026-09-08 with --tool gnina on the nearest table, out-dir inside the program folder)",
    True, "W8c: 'level' not 'gnina higher'. Footnote: both re-rankings give 41.6 % against 43.6 %.")
add("M:199", "selected variants sentence", {}, "plan 1.5 (all three selections hold)", False, "no change")
add("M:205-214", "provenance comments: Form poses 4,012 / 3,903 / 3,873 with guard; 4,030 / 3,915 / 3,885 without; complexes 225 / 237 / 235; 4012/9013 = 44.5 %",
    {"form_poses_with_guard": [4016, 3887, 3854], "form_poses_without_guard": [4034, 3899, 3866], "form_complexes": [224, 236, 235], "pct_with_guard": [44.6, 43.2, 42.9],
     "columns_to_name": ["rmsd", "bestfit_rmsd", "pb_rmsd", "centroid_dist (nearest-copy values)", "rmsd_ref_instance", "bestfit_rmsd_ref_instance", "nearest_copy_index", "reference_convention"]},
    SIG + " (table_1 diffdock*.form_poses / form_pct); without guard: " + HUB + "pose_validity_cascade.csv kabsch1_poses", False,
    "Comment lines only. Form complex counts are no longer 'identical with or without the guard' claim-safe for DiffDock raw (224 with guard).")
add("M:240", "unidock2 row deliberately omitted", {}, "plan", False, "unchanged")
add("M:264/277-283", "Table 1 headers", {}, "change list W9 (optional footnote)", True, "optional footnote naming the nearest copy")

t1 = {
 "AutoDock (raw) exh128":   {"poses": 8976, "valid": [303, 8943, 99.6], "near": [217, 437, 4.9], "form": [225, 2829, 31.5], "combined": [178, 324, 3.6]},
 "AutoDock + gnina":        {"poses": 8976, "valid": [301, 8880, 98.9], "near": [217, 451, 5.0], "form": [227, 2823, 31.5], "combined": [175, 325, 3.6]},
 "DiffDock (raw)":          {"poses": 9013, "valid": [242, 2198, 24.4], "near": [161, 2118, 23.5], "form": [224, 4016, 44.6], "combined": [92, 869, 9.6]},
 "DiffDock + smina":        {"poses": 8993, "valid": [300, 7597, 84.5], "near": [186, 2169, 24.1], "form": [236, 3887, 43.2], "combined": [148, 1516, 16.9]},
 "DiffDock + gnina":        {"poses": 8993, "valid": [300, 7831, 87.1], "near": [186, 2171, 24.1], "form": [235, 3854, 42.9], "combined": [149, 1513, 16.8]},
 "EquiBind (raw) unguided": {"poses": 9090, "valid": [35, 263, 2.9],    "near": [19, 103, 1.1],   "form": [186, 2487, 27.4], "combined": [2, 21, 0.2]},
 "EquiBind + smina":        {"poses": 9090, "valid": [236, 4664, 51.3], "near": [67, 465, 5.1],   "form": [192, 2617, 28.8], "combined": [50, 253, 2.8]},
 "EquiBind + gnina":        {"poses": 9088, "valid": [266, 5633, 62.0], "near": [86, 515, 5.7],   "form": [199, 2661, 29.3], "combined": [57, 289, 3.2]},
 "EquiBind (raw) fpocket":  {"poses": 9090, "valid": [13, 27, 0.3],     "near": [0, 0, 0.0],      "form": [179, 2441, 26.9], "combined": [0, 0, 0.0]},
 "EquiBind (raw) P2Rank":   {"poses": 8520, "valid": [12, 47, 0.6],     "near": [1, 2, 0.0],      "form": [180, 2339, 27.5], "combined": [0, 0, 0.0]},
}
add("M:286-304", "Table 1 rows (ten arms; near / form / combined columns as printed 2026-09-08 morning)", {"rows_[complexes, poses, pct]": t1},
    SIG + " table_1[...] rows (near_*, form_*, triple_*; form uses the :205 guard bestfit<=1 AND rmsd<1000) cross-checked against " + HUB + "pose_validity_cascade.csv (rmsd2_* and triple_* identical; kabsch1_poses differ only by the guard)",
    False, "Validity columns unchanged. Form block LOSES a complex in five rows (DiffDock raw 225->224, DiffDock smina 237->236, EquiBind gnina 200->199, fpocket 180->179, P2Rank 181->180); EquiBind raw/smina GAIN one (185->186, 191->192). P2Rank raw near-native 0 -> 1 complex (2 poses, 0.0 %).")
add("M:324", "removes at most four of 303 complexes from the headline endpoint at any depth",
    {"max_validity_gate_cost_complexes": 5, "where": "DiffDock + smina rank-1 (137 near-native -> 132 gated)", "by_depth_1_5_15_30": {"AutoDock + gnina": [2, 4, 3, 3], "DiffDock + smina": [5, 4, 3, 3], "EquiBind + gnina": [2, 3, 3, 3]}, "span_pp": [0.66, 1.65]},
    HUB + "validity_gate_cost.csv", True, "W10: 'at most five'. Same at M:777.")
add("M:326", "recovered if any inspected pose is PoseBusters-valid and within 2 A", {}, "change list W11", True, "W only")
add("M:347-349", "Table 2 rows: 36.6 % (111) / 65.3 % (198) / 65.7 % (199) / +87 (+28.7 %) *** / +1 (+0.3 %); 34.0 % (103) / 55.1 % (167) / 55.8 % (169) / +64 (+21.1 %) *** / +2 (+0.7 %); 18.2 % (55) / 26.4 % (80) / 26.4 % (80) / +25 (+8.3 %) *** / +0",
    {"AutoDock": {"rank1": "49.2 % (149)", "top15": "70.3 % (213)", "top30": "70.6 % (214)", "gain_1_15": "+64 (+21.1 %) ***", "gain_15_30": "+1 (+0.3 %)"},
     "DiffDock": {"rank1": "43.6 % (132)", "top15": "59.1 % (179)", "top30": "60.4 % (183)", "gain_1_15": "+47 (+15.5 %) ***", "gain_15_30": "+4 (+1.3 %)"},
     "EquiBind": {"rank1": "18.8 % (57)", "top15": "27.4 % (83)", "top30": "27.4 % (83)", "gain_1_15": "+26 (+8.6 %) ***", "gain_15_30": "+0 (+0.0 %)"},
     "gain_15_30_one_sided_binomial_p": {"AutoDock b=1": 1.0, "DiffDock b=4,c=0": 0.125, "EquiBind b=0": 1.0}},
    TOPK + " (pbvalid_near_native_k / _pct, gain_from_prev_depth_k / _pp)", False,
    "Stars on the 1->15 step follow the canonical convention (Wilson interval excludes zero); the 15->30 DiffDock step b=4 c=0 p=0.125 carries no star.")
add("M:353-372", "footnote 'the two coincide at every depth shown here'; comments 111/198/199, raw 94/185/202, +88/+63/+26, 15->30 p = 1.0/0.5/1.0",
    {"printed_arm_1_15_30": [149, 213, 214], "raw_exh128_arm_1_15_30_gated": [130, 198, 217], "near_native_no_validity_gain_1_15": [65, 45, 27], "gated_gain_1_15": [64, 47, 26],
     "step_15_30_p": {"AutoDock": 1.0, "DiffDock": 0.125, "EquiBind": 1.0}, "coincide_with_table_26_at_2A": True},
    TOPK + "; raw arm: " + HUB + "validity_gate_cost.csv (AutoDock Vina (raw) gated k=1/15/30); coincidence with Table 26: " + DEPTHS + " (2 A rows 49.2/70.3/70.6, 43.6/59.1/60.4, 18.8/27.4/27.4)", False, "Comment lines; update counts.")
add("M:380", "as-placed crystal RMSD (i.e. near-nativeness) and best-fit Kabsch RMSD (i.e. form)", {}, "change list W12", True, "W only")
add("M:382", "Q p < 1e-7; 36.6 / 34.0 / 18.2 %; all p_holm < 3e-7; +2.6 (-4.2, +9.4) p 0.509; +13.2 at top-5 (+5.7, +20.4) p_holm 7.6e-4; +10.2 at top-15 (+2.8, +17.5) p_holm 0.009; +9.9 at top-30 (+2.5, +17.1) p_holm 0.011",
    {"cochran_Q_p_by_depth": {"rank1": "2.0e-17", "top15": "8.6e-29", "top30": "1.6e-29"}, "holds_p_lt_1e-7": True,
     "rank1_pct": [49.2, 43.6, 18.8], "vs_equibind_max_p_holm": "8.1e-13 (holds < 3e-7)",
     "rank1_AD_vs_DD": {"pp": 5.6, "newcombe95": [-1.7, 12.8], "mcnemar_p": 0.159},
     "top5": {"pp": 13.2, "newcombe95": [5.7, 20.4], "p_holm": "7.6e-4", "identical_under_both_conventions": True},
     "top15": {"pp": 11.2, "newcombe95": [3.9, 18.4], "p_holm": 0.0036},
     "top30": {"pp": 10.2, "newcombe95": [3.0, 17.3], "p_holm": 0.0075}},
    TOPK + " (between_tools_cochran_p, mcnemar_p_holm); intervals: " + RC + "reference_convention_margins.csv and " + BC, False,
    "'do not differ detectably at rank-1 ... but they do from top-5 onwards' holds.")
add("M:386", "112 discordant complexes; MDD 9.8 pp at 80 % power; observed +2.6",
    {"n_discordant": 129, "mde_pp_80pct_power": 10.5, "observed_power": 0.32, "observed_diff_pp": 5.6, "n_needed_for_observed_effect": 1062, "prose_rounding": "roughly 1,100"},
    BC + " (n_discordant 129, mde_pp 10.5, power 0.322, n_needed 1062); confirmed by stats_utils.mcnemar_power on the rank-1 vectors (" + RECOMP + ")", False, "")
add("M:388", "Q p < 1e-6; +28.7 (23.9-34.0); +21.1 (16.9-26.1); +8.3 (5.7-11.9); 36.6 -> 65.3, 87 vs 64; 31.0 -> 36.6 and 61.1 -> 65.3; 202 vs 199",
    {"within_tool_cochran_p": {"AutoDock": "1.6e-28", "DiffDock": "2.8e-21", "EquiBind": "5.1e-12"},
     "gain_1_15_pp_wilson95": {"AutoDock": [21.1, 16.9, 26.1], "DiffDock": [15.5, 11.9, 20.0], "EquiBind": [8.6, 5.9, 12.3]},
     "autodock_rank1_to_top15_pct": [49.2, 70.3], "added_complexes": {"AutoDock": 64, "DiffDock": 47},
     "raw_exh128_gated_rank1_pct": 42.9, "gnina_rank1_pct": 49.2, "raw_exh128_gated_top15_pct": 65.3, "gnina_top15_pct": 70.3,
     "whole_pool_raw_vs_rescored": [217, 214]},
    TOPK + " (within_tool_cochran_p, gain_wilson_lo/hi); raw arm 130/198 of 303: " + HUB + "validity_gate_cost.csv; pool 217 vs 214: " + HUB + "pose_validity_cascade.csv rmsd2_pbvalid_complexes", False,
    "'raising rank-1 from 42.9 % to 49.2 % and top-15 from 65.3 % to 70.3 %'.")
add("M:393", "caption: Best-of-top-N accuracy against the as-placed crystal RMSD", {}, "change list W13; regenerate media/media/image3.png", True, "W; Fig 3 caption at :401 holds (regenerate image4).")
add("M:396", "74.3 / 74.6 %; 95.4 / 97.0 %; Q p < 1e-9; 52.5 / 52.1 / 33.3 %; p_holm <= 4e-10; tie at p_holm = 1.000 with 47 / 46, 30 / 31, 25 / 26",
    {"best_of_top30_kabsch_1A_pct": [74.3, 74.6], "best_of_top30_kabsch_2A_pct": [95.4, 97.0], "cochran_p_by_depth": {"rank1": "1.2e-10", "top15": "2.0e-14", "top30": "1.3e-15"},
     "rank1_kabsch_1A_pct": [52.5, 52.1, 33.3], "vs_equibind_p_holm_max": "9.0e-9 (rank-1 8.997e-09; deeper <= 6.0e-11) -> write '<= 9e-9'",
     "AD_vs_DD_discordant": {"rank1": [48, 47], "top15": [30, 32], "top30": [25, 26]}, "AD_vs_DD_p_holm": {"rank1": 1.0, "top15": 0.899, "top30": 1.0}},
    KAB, True, "A2 (<= 9e-9) applies today; W14 tie wording ('not significant at any depth' or 'p_holm = 1.000, 0.899 and 1.000').")
add("M:401", "caption: Best-of-top-N accuracy against the best-fit (Kabsch) RMSD", {}, "plan", False, "holds; regenerate image4")
add("M:404", "33.3 % -> 54.1 %; 20.8-point gain; placement gains at most 0.7; form gains 2.3 / 2.3 / 3.0",
    {"equibind_kabsch_rank1_to_top30_pct": [33.3, 53.8], "depth_gain_pp": 20.5, "placement_top15_to_top30_max_pp": 1.3, "placement_by_tool_pp": [0.3, 1.3, 0.0],
     "form_1A_top15_to_top30_pp": {"AutoDock": 2.6, "DiffDock": 2.3, "EquiBind": 2.3}, "form_complexes_added": [8, 7, 7]},
    KAB + " (section B: top-15->top-30 +8/+7/+7 complexes = 2.6/2.3/2.3 pp; best-of-top-30 53.8 %); placement: " + DEPTHS + " (+1/+4/+0 complexes = 0.3/1.3/0.0 pp)", False, "")
add("M:406", "top-15 selection gains 87 and 64; six in ten; 98.9 vs 84.5 %; 'worth about three points of rank-1 recovery on the reported coordinates'",
    {"selection_gains_1_15": [64, 47], "rank1_failures_no_qualifying_pose_pct": [57.8, 70.2], "wording": "six to seven in ten", "validity_pct_unchanged": [98.9, 84.5],
     "refined_score_ordering_gain_pp": {"smina": 3.6, "gnina": 3.0, "basis": "B_full minus D_minimize (affinity re-rank on minimised coordinates vs confidence order on the same coordinates; RMSD <= 2, no validity gate)"}},
    TOPK + "; PoseBusters_Benchmark_Analysis/smina_rerank_nearest/selection_strategy_summary_smina.csv (48.84 - 45.21) and Scripts/Analysis/nearest_copy_program/gnina_rerank_nearest/selection_strategy_summary_gnina.csv (48.51 - 45.54)", True,
    "'about three points' -> smina 3.6 (gnina 3.0). Editor decides 'about three and a half' or keeps 'about three'.")
add("M:408", "footnote: 8 A from the crystal-ligand site; 6,169 PB-valid poses; 2,329 poses clipped beyond 5 A",
    {"near_site_pbvalid_top15_pool": 7565, "dropped_beyond_8A_top15": 4255, "clipped_beyond_5A_top15": 2854},
    F20 + "caption.txt (dropped 4255, beyond frame 2854) and " + F20 + "cohort_summary.csv (2832 + 2801 + 1932 = 7565)", True, "W15 wording 'nearest deposited copy'.")
add("M:410", "62.7 % remain; 1.5 -> 5.0 A; +3.55 (+3.00 to +3.99); Friedman p < 1e-15; placement-limited 48 -> 78 %; form-limited 27 -> 8 %; spread 5.72 A",
    {"rank1_coverage_pct": 82.2, "median_inplace_rank1_top15_A": [1.4, 4.9], "drift_A": 3.55, "drift_bootstrap95": [3.08, 3.95], "friedman_paired_p": "1.9e-38 (< 1e-15 holds)",
     "placement_limited_pct_rank1_top15": [47, 77], "form_limited_pct_rank1_top15": [28, 9], "within_complex_spread_A": 5.83},
    F20 + "stats_per_tool_depth.csv (coverage 82.18, medians 1.376 / 4.922, pct 46.59->77.30, 27.71->8.76), " + F20 + "stats_within_tool_trend.csv (drift 3.545, CI 3.079-3.949, p 1.9e-38), " + F20 + "stats_pose_diversity.csv (5.83)", False, "")
add("M:412", "+0.46 vs +3.55 and +1.09; 1.5-1.9 A in place and 0.9 in form; form shift +0.03 (-0.05 to +0.12); 137 complete complexes Friedman p = 0.006; spread 2.34 vs 5.72",
    {"drift_inplace_A": {"DiffDock": 0.53, "AutoDock": 3.55, "EquiBind": 0.96}, "diffdock_inplace_medians_1_5_15": [1.53, 1.68, 2.06], "range_wording": "1.5-2.1", "diffdock_form_median_A": 0.9,
     "form_shift_A": 0.05, "form_shift_bootstrap95": [-0.03, 0.13], "friedman_n_complexes": 203, "friedman_p": 0.00056, "spread_A": {"DiffDock": 2.81, "AutoDock": 5.83},
     "best_on_both_axes_top5_top15": "holds (Wilcoxon Holm p 2.0e-4 / 4.5e-9 in place, 7.3e-7 / 4.6e-8 form vs AutoDock; all < 0.001 vs EquiBind)"},
    F20 + "stats_within_tool_trend.csv, " + F20 + "stats_per_tool_depth.csv, " + F20 + "stats_crosstool_paired.csv, " + F20 + "stats_pose_diversity.csv", False,
    "'137 complete complexes' -> 203; 'Friedman p = 0.006' -> 5.6e-4; form shift +0.05 (-0.03 to +0.13) is still detected but negligible.")
add("M:414", "266 of 303; 26.4 % by top-15 and 26.4 % over all thirty; rank-1 median in-place 2.9 A",
    {"valid_complexes": 266, "top15_pct": 27.4, "top30_pct": 27.4, "rank1_near_site_median_inplace_A": 3.2},
    TOPK + "; " + F20 + "stats_per_tool_depth.csv (EquiBind depth 1 inplace_median 3.239); A3 today: canonical generator prints 2.977 -> 3.0", False, "A3 (2.9 -> 3.0) today, then 3.2 under nearest.")
add("M:419", "caption: Form versus in-place RMSD across ranking depth", {}, "plan (regenerate image5)", False, "caption holds")
add("M:422", "near-site interpretive paragraph ('few for DiffDock', 'form-limited only because ...')",
    {"pct_placement_limited_rank1_top5_top15": {"AutoDock": [46.6, 63.7, 77.3], "DiffDock": [41.5, 44.5, 47.1], "EquiBind": [60.6, 64.8, 70.5]},
     "pct_form_limited": {"AutoDock": [27.7, 16.7, 8.8], "DiffDock": [29.5, 24.3, 22.4], "EquiBind": [16.8, 15.0, 11.4]},
     "median_form_A_top15": {"AutoDock": 1.416, "DiffDock": 0.936, "EquiBind": 1.354}, "median_r_form_share_top15": {"AutoDock": 0.097, "DiffDock": 0.377, "EquiBind": 0.132}},
    F20 + "cohort_summary.csv and " + F20 + "stats_per_tool_depth.csv", True,
    "W16: re-read. The reading still holds in direction (DiffDock placement-limited share rises only 41.5 -> 47.1 %, AutoDock 46.6 -> 77.3 %; DiffDock has the lowest absolute form RMSD and the highest form-to-in-place ratio).")
add("M:427", "Sixteen of 303 complexes have no qualifying cluster", {"complexes_with_no_correct_cluster": 6},
    CL + "summary.json (n_complexes_no_correct_cluster 6; n_true_site_complexes 297)", True, "'Six of 303'.")
add("M:432", "caption: Top-N crystal-cluster co-recovery", {}, "regenerate media/media/image6.png", False, "caption holds")
add("M:435", "DiffDock gains 26 pp, AutoDock 24, EquiBind 4 from rank-1 to top-15; phi = +0.23 to +0.28 with Fisher p <= 2.8e-4",
    {"reach_pct_rank1_top15": {"AutoDock": [78, 90], "DiffDock": [73, 86], "EquiBind": [46, 51]}, "reach_k_rank1_top15": {"AutoDock": [237, 272], "DiffDock": [222, 261], "EquiBind": [138, 155]},
     "gain_pp_1_15": {"DiffDock": 12.9, "AutoDock": 11.6, "EquiBind": 5.6}, "top15_phi_range": [0.13, 0.22], "top15_fisher_p_max": 0.036,
     "intervals_overlap_AD_DD_every_depth": True},
    CL + "topN_crystal_cluster_matrix_stats.txt; exact reach counts from " + CL + "per_complex_summary.csv (<tool>_rank_in_correct_cluster <= N and correct_cluster_is_hit)", False,
    "phi now +0.13 to +0.22, Fisher p <= 0.036 (AutoDock + EquiBind pair drops to * at top-15).")
add("M:455-461", "Table 3: reach 59/79/81/83 [CI], 54/76/79/81, 43/46/47/47; phi +0.33/+0.23/+0.24/+0.23, +0.45/+0.27/+0.29/+0.28, +0.42/+0.30/+0.29/+0.26; all three 89/42, 114/82, 121/90, 125/96",
    {"AutoDock reach % [95% CI]": {"top1": [78, 73, 82], "top5": [86, 81, 89], "top10": [88, 84, 91], "top15": [90, 86, 93]},
     "DiffDock reach % [95% CI]": {"top1": [73, 68, 78], "top5": [83, 78, 86], "top10": [85, 80, 88], "top15": [86, 82, 90]},
     "EquiBind reach % [95% CI]": {"top1": [46, 40, 51], "top5": [49, 43, 54], "top10": [50, 45, 56], "top15": [51, 46, 57]},
     "AutoDock + DiffDock phi": [0.24, 0.21, 0.21, 0.21], "AutoDock + EquiBind phi": [0.19, 0.15, 0.15, 0.13], "DiffDock + EquiBind phi": [0.21, 0.22, 0.22, 0.22],
     "phi_fisher_p": {"AD+DD": ["<1e-4", "7.7e-4", "0.001", "0.001"], "AD+EB": ["7.7e-4", "0.009", "0.009", "0.036"], "DD+EB": ["4.0e-4", "1.2e-4", "1.0e-4", "1.9e-4"]},
     "All three reach observed / expected": {"top1": "103 / 79", "top5": "125 / 105", "top10": "132 / 113", "top15": "137 / 120"}, "permutation_p": "< 1e-4 at every depth",
     "footnote_check": "all three pairings p < 0.05 two-sided Fisher: holds (max 0.036)"},
    CL + "topN_crystal_cluster_matrix_stats.txt (matches " + SIG + " table_3 rows)", False, "Bold on phi cells: all still p < 0.05.")
add("M:469", "caption: Rates are taken over all 303 complexes, 18 of which contribute no pose to any panel",
    {"complexes_with_no_tool_pose_in_crystal_closest_cluster_top15": 9},
    CL + "crystal_cluster_homogeneity_stats.json (B_consensus observed [9, 37, 120, 137] for 0/1/2/3 tools); confirmed from per_complex_summary.csv", False, "'9 of which'; regenerate image7.")
add("M:474", "2.13 -> 3.72 A radius; partial rho -0.07 p 0.26; bootstrap Jaccard 0.78; oracle 95 %; about 56 %; 56.8 / 55.4 / 55.8 %; not separable; gains 7.3 / 5.9 / 6.3 on 49.5 %; Holm p 3.6e-4 / 0.004 / 0.008; 38-point gap",
    {"median_radius_A_1tool_3tools": [1.63, 3.62], "partial_rho": -0.033, "partial_p": 0.57, "bootstrap_jaccard_median": 0.786,
     "oracle_precision_at_1_pct": 98.0, "best_inexpensive_rule_pct": 71.3, "consensus_pct": 71.3, "medoid_pct": 69.0, "confidence_pct": 68.3, "legacy_size_pct": 59.7,
     "gains_on_legacy_pp": {"consensus": 11.6, "medoid": 9.2, "confidence": 8.6},
     "holm_mcnemar_vs_legacy": {"consensus": "4.8e-6", "confidence": "3.4e-4", "medoid": "3.6e-4"},
     "three_rules_separable": "consensus beats medoid 7/0, Holm p 0.047 (*); consensus vs confidence 0.326 ns; medoid vs confidence 0.868 ns",
     "oracle_gap_points": 26.7},
    CL + "crystal_cluster_homogeneity_stats.json (D_radius_trend medians, posecount_control); " + CL + "per_complex_summary.csv boot_stability median; " + CL + "cluster_quality_metrics_stats.txt and summary.json ranking_ablation(_significance)", True,
    "Wording: 'These three rules are not separable from one another' no longer holds strictly (consensus vs medoid Holm p 0.047). '38-point oracle gap' -> '27-point'. 'about 56 %' -> 'about 71 %'.")
add("M:476", "97.7 % at 0.34 A; AutoDock 88.1 % at 0.53; DiffDock 84.5 % at 0.56; together 97.7 % at 0.37; pooled centre 1.23 A",
    {"pooled_three_tool_oracle_pct": 99.7, "pooled_median_A": 0.239, "AD_pct_median": [92.4, 0.396], "DD_pct_median": [90.1, 0.386], "AD_DD_pct_median": [99.7, 0.268], "ensemble_cluster_centre_median_A": 1.088,
     "equibind_adds_nothing": "holds (AD+DD 0.9967 == AD+DD+EB 0.9967)", "convention_note": "cluster report any-copy rule (D5); hub per-pose D2 value would be 99.3 %"},
    CL + "ensembles_summary.csv (oracle_hit_rate, median_oracle_dist) and summary.json median_oracle_dist.ensemble_clusters", False, "Loader fix (heavy atoms only) is included in the rebuilt report.")
add("M:481", "Crystal-pose interactions are compared with those of PoseBusters-valid pipeline poses", {}, "change list W17", True, "W (per-copy fingerprint + union rule).")
add("M:486", "caption: Mean Jaccard similarity of interaction fingerprints, whole-protein search", {}, "regenerate media/media/image9.png", False, "caption holds")
add("M:489", "0.450 over 300 / 0.417 over 301 / 0.380 over 266; common 264: 0.474 / 0.431 / 0.381; DD > AD Wilcoxon p 0.014; vs EquiBind p 1.6e-3 and 5.2e-8; AutoDock 0.313 on 37 lacking EquiBind; pairwise 0.357-0.463; AD-DD 0.463 exceeds either tool's similarity to the crystal",
    {"mean_jaccard_top5_union_vs_crystal": {"DiffDock": [0.572, 300], "AutoDock": [0.539, 301], "EquiBind": [0.403, 266]},
     "common_cohort_264": {"DiffDock": 0.586, "AutoDock": 0.546, "EquiBind": 0.403},
     "wilcoxon_common": {"DD_vs_AD": 0.011, "AD_vs_EB": "1.1e-7", "DD_vs_EB": "7.2e-13"}, "wilcoxon_pairwise_300_DD_vs_AD": 0.026,
     "autodock_mean_on_37_lacking_equibind": 0.491, "pairwise_tool_similarity_range": [0.364, 0.463], "AD_DD_pair": 0.463,
     "pair_exceeds_crystal_similarity": False},
    PM + "interaction_pose_basis_audit.txt section C 'published (all poses)' (0.5391 (301) / 0.5722 (300) / 0.4030 (266), common n 264 p 1.08e-02, pair n 300 p 2.65e-02); remaining cells recomputed 2026-09-08 with Scripts/Analysis/interaction_pose_basis_audit.py functions (build_copy_selector, union_jaccard depth 5, typed_fp) on pandamap_results/benchmark_matched_equibind_nearest/{pandamap_interactions,crystal_interactions}.csv joined to " + PPM,
    True, "WORDING FLIP: the AutoDock--DiffDock pair (0.463) no longer exceeds either tool's similarity to the crystal (0.539 / 0.572); the sentence 'The two therefore resemble each other more than either resembles the reference' must be rewritten. Tool-tool Jaccard values are convention-invariant (pose fingerprints unchanged). No sidecar prints the matrix as numbers (05_fingerprint_similarity_top5.txt holds only the footer).")
add("M:491", "5 of 303 AutoDock complexes and 21 of 303 with DiffDock output lack a valid rank-1", {"n_k1": {"AutoDock": 298, "DiffDock": 282, "EquiBind": 266}}, PM + "native_recovery_by_rank.csv", False, "unchanged")
add("M:493", "F1 0.57 / 0.53 / 0.47; 280 shared: 0.5663 vs 0.5324 (p 0.25); zero-scored over 303: 0.556 / 0.492; by rank 5 AutoDock 0.42, DiffDock 0.51, EquiBind 0.44; DiffDock confidence has no detected association",
    {"rank1_f1": {"AutoDock": 0.73, "DiffDock": 0.69, "EquiBind": 0.50}, "shared_280": {"n": 280, "AutoDock": 0.7381, "DiffDock": 0.6986, "wilcoxon_p": 0.064},
     "zero_filled_303": {"AutoDock": 0.720, "DiffDock": 0.646, "EquiBind": 0.438}, "rank5_f1": {"AutoDock": 0.54, "DiffDock": 0.63, "EquiBind": 0.46},
     "kendall_tau_f1": {"AutoDock": "median tau -0.40, q < 1e-4", "DiffDock": "median tau -0.11, p 0.002, BH q 0.002 (NOW DETECTED)", "EquiBind": "median tau -0.20, q 7.7e-4"},
     "diffdock_top5_union_leads_without_rank1_lead": "holds (rank-1 AD 0.73 > DD 0.69; top-5 union DD 0.572 > AD 0.539)"},
    PM + "native_recovery_by_rank.csv; shared-cohort and zero-filled means recomputed from " + PM + "recovery_detail_per_pose.csv (pose_rank == 1, scipy wilcoxon); tau: " + PM + "04c_native_recovery_by_rank.txt", True,
    "WORDING FLIP: 'DiffDock confidence has no detected association with interaction recovery, so its top pose cannot be expected ... to outperform rank 5. This is absence of detected ordering' is now false (all three DiffDock tau tests significant after BH). 'nine BH-corrected tau tests' still nine.")
add("M:514-525", "Table 4 (k = 1..5 n / precision / recall / F1 per tool)",
    {"AutoDock": {"n": [298, 298, 299, 296, 299], "precision": [0.727, 0.690, 0.641, 0.615, 0.555], "recall": [0.741, 0.695, 0.633, 0.596, 0.532], "f1": [0.732, 0.689, 0.633, 0.601, 0.539]},
     "DiffDock": {"n": [282, 267, 277, 270, 275], "precision": [0.704, 0.682, 0.658, 0.645, 0.650], "recall": [0.692, 0.667, 0.631, 0.617, 0.622], "f1": [0.694, 0.671, 0.639, 0.626, 0.632]},
     "EquiBind": {"n": [266, 256, 247, 245, 243], "precision": [0.522, 0.518, 0.518, 0.508, 0.497], "recall": [0.486, 0.475, 0.475, 0.459, 0.444], "f1": [0.499, 0.491, 0.490, 0.476, 0.463]}},
    PM + "native_recovery_by_rank.csv (4-decimal values rounded to 3)", False, "n columns unchanged.")
add("M:530", "rank 1: 21.0 / 19.4 matched, 15.3 / 16.4 missed; rank 5 AutoDock 15.3 matched vs 21.0 missed, spurious 17.6; only pipeline with rising spurious; DiffDock near 18.6 matched / 17.4 missed; EquiBind misses more than it recovers at every rank",
    {"rank1_matched": {"AutoDock": 27.2, "DiffDock": 25.8}, "rank1_missed": {"AutoDock": 9.0, "DiffDock": 9.9},
     "autodock_rank5": {"matched": 19.5, "missed": 16.7, "spurious": 13.5, "spurious_rank1": 9.1},
     "diffdock_rank1_to_5": {"matched": [25.8, 23.0], "missed": [9.9, 12.8], "spurious": [7.8, 8.6]},
     "equibind_missed_gt_matched_every_rank": "holds (17.8 vs 17.7 at k=1 ... 19.0 vs 16.0 at k=5)",
     "spurious_trend_tau": {"AutoDock": "+0.32, q < 1e-4", "DiffDock": "+0.00 median, rank-biserial +0.17, p 0.019, BH q 0.021 (*)", "EquiBind": "ns"}},
    PM + "10b_contact_decomposition_by_rank.csv and .txt", True,
    "WORDING: 'AutoDock falls to 15.3 matched against 21.0 missed' is no longer true (19.5 matched still exceeds 16.7 missed); 'the only pipeline whose spurious count increases significantly with rank' no longer strictly holds (DiffDock spurious q 0.021).")
add("M:532", "no detected difference in nominated rank-1 contact fidelity; DiffDock's flatter pool remains unordered", {"shared_280_wilcoxon_p": 0.064, "diffdock_tau_f1_q": 0.002},
    "as M:493", True, "'DiffDock's flatter pool remains unordered' is contradicted by the BH-significant DiffDock tau; rewrite.")
add("M:697", "Calibration point estimates show the same ordering, although not a statistically separable difference",
    {"rank1_verdict": "DiffDock 43.6 % < AutoDock 49.2 %, p 0.159 (not separable); from top-5 AutoDock leads, resolved"}, "change list W18; " + TOPK, True, "W: 'at rank-1'.")
add("M:708", "within 2 A of the crystal ligand (cost definition)", {}, "change list W19", True, "W only")
add("M:715", "Fig 15 caption: 199, 169 and 80 qualifying complexes; 53 shared", {"qualifying_complexes": [214, 183, 83], "all_three_shared": 54},
    EFQ + " tiers.near2.per_method.n_complexes and n_paired", False, "regenerate image24")
add("M:737-741", "Table 6: AutoDock 137.2 / 75.6-167.4 / 199 / 1,752.8 / 118.9 / 311; DiffDock 49.4 / 14.1-127.3 / 169 / 2.93 / 49.9 / 1,622; EquiBind 2.4 / 1.0-7.3 / 80 / 165.3 / 1.4 / 466; Vina only 21.7 / 9.9-45.3 / 202 / 1,805.0 / 0.0 / 302",
    {"AutoDock": {"median_charged_s": 91.8, "iqr": [59.2, 153.7], "complexes": 214, "cpu_core_s_per_near2": 1233.3, "gpu_s_per_near2": 83.6, "poses_near2": 442},
     "DiffDock": {"median_charged_s": 39.0, "iqr": [11.7, 115.4], "complexes": 183, "cpu_core_s_per_near2": 2.30, "gpu_s_per_near2": 39.1, "poses_near2": 2068},
     "EquiBind": {"median_charged_s": 2.4, "iqr": [1.0, 7.3], "complexes": 83, "cpu_core_s_per_near2": 159.4, "gpu_s_per_near2": 1.4, "poses_near2": 483},
     "AutoDock, Vina search only": {"median_charged_s": 17.8, "iqr": [7.5, 37.4], "complexes": 217, "cpu_core_s_per_near2": 1264.8, "gpu_s_per_near2": 0.0, "poses_near2": 431}},
    "medians / IQR / complexes: " + EFQ + " tiers.near2.per_method (median_s, q1_s, q3_s, n_complexes); per-pose CPU/GPU and poses: " + EFSUM + " (cpu_core_s_per_near2, gpu_s_per_near2, poses_near2). Vina-only row: no sidecar; recomputed from " + EFPC + " (autodock cpu_core_s / 32 per complex) divided by the raw exh128 arm's PB-valid AND rmsd <= 2 pose count per complex from " + PPM + " (median / IQR over complexes with >= 1 such pose; CPU per pose = 151.4192 h * 3600 / 431). The same recipe reproduces the canonical row 21.7 / 9.9-45.3 / 202 / 1,805.0 / 302 exactly.",
    False, "Charged basis = cpu_core_s / 32 + gpu_s per complex.")
add("M:744", "most expensive route; Friedman chi2 85.2, p 3.2e-19, W 0.80 on 53; every pair passes Holm; ratios 0.14 and 0.02; three times the search alone; median successful complex yields one; 199 of 303 vs 169",
    {"ordering_holds_AutoDock_most_expensive": True, "shared_complexes": 54, "friedman_chi2": 87.1, "friedman_p": "1.2e-19", "kendall_W": 0.81, "every_pair_holm": "yes (max p_holm 3.9e-5)",
     "median_per_complex_cost_ratio": {"DiffDock_vs_AutoDock": 0.15, "EquiBind_vs_AutoDock": 0.02}, "gnina_multiplier_of_search": "15.00 / 4.73 = 3.2 (unchanged)",
     "median_qualifying_poses_per_succeeding_complex": 2, "complexes_yielding": {"AutoDock": 214, "DiffDock": 183}},
    EFQ + " tiers.near2 (n_paired 54, omnibus chi2 87.111 p 1.21e-19 kendall_w 0.8066, pairwise p_holm, median_ratios 6.741 and 57.27 -> 1/x); median poses per succeeding complex from " + EFPC + " (autodock near2 among near2 > 0, median 2.0)", True,
    "W20a: 'yields two'.")
add("M:746", "5.3, 5.3 and 137.2; 7.6, 9.2 and 49.4; 0.3, 0.4 and 2.4; 303 -> 301 -> 199; 303 -> 300 -> 169; 303 -> 266 -> 80; multipliers 25.9, 5.4, 6.0; fails in the other 223",
    {"median_charged_s_generated_valid_near2": {"AutoDock": [5.3, 5.3, 91.8], "DiffDock": [7.6, 9.2, 39.0], "EquiBind": [0.3, 0.4, 2.4]},
     "contributing_complexes": {"AutoDock": [303, 301, 214], "DiffDock": [303, 300, 183], "EquiBind": [303, 266, 83]},
     "valid_to_near2_multiplier": {"AutoDock": 17.3, "DiffDock": 4.3, "EquiBind": 6.4}, "equibind_fails_in_other": 220,
     "clause_check": "'a complex that succeeds contributes only one or two qualifying poses out of thirty' still holds (median 2)"},
    EFQ + " tiers.all / valid / near2 per_method median_s and n_complexes (multiplier = near2 median / valid median: 91.83/5.31, 39.04/9.17, 2.41/0.37)", False, "")
add("M:751", "Fig 16 caption", {}, "regenerate image25", False, "caption holds")
add("M:754", "most expensive on every axis; 151.42 CPU-core-h into 4.73 wall-h; 311, 1,622 and 466 qualifying poses; highest GPU cost per qualifying pose",
    {"qualifying_poses": [442, 2068, 483], "cpu_core_h": 151.42, "wall_h_search": 4.73, "gpu_s_per_near2": {"AutoDock": 83.6, "DiffDock": 39.1, "EquiBind": 1.4}, "most_expensive_every_axis": True},
    EFSUM, False, "")
add("M:756", "4.73 of 15.00 charged hours; 10.27 GPU; 21.39 CPU-core-h EquiBind; Cochran Q 61.1, p 5.5e-14",
    {"unchanged": [4.73, 15.00, 10.27, 21.39, 151.42], "cochran_Q_has_valid_pose": 61.1, "p": "5.5e-14"}, EFS + " validity.omnibus (Q 61.077, p 5.46e-14) - validity-based, convention-invariant", False, "no change")
add("M:758", "Vina median 21.7 s and no GPU time, still below DiffDock; CPU 1,805.0 s vs 1,752.8 s; adds nine qualifying poses",
    {"vina_only_median_charged_s": 17.8, "below_diffdock_39.0": True, "vina_cpu_core_s_per_near2": 1264.8, "rescored_cpu_core_s_per_near2": 1233.3, "rescoring_adds_qualifying_poses": 11, "wording": "eleven"},
    "as M:737-741 (Vina-only recompute; " + EFSUM + " cpu_core_s_per_near2 1233.28; 442 - 431 = 11 from " + HUB + "pose_validity_cascade.csv rmsd2_pbvalid_poses)", True, "W20b 'eleven'.")
add("M:760", "157.7 / 227.6 / 8.1 s per attempted docking; 271, 479 and 38 s per success; 186 vs 478 GPU-seconds",
    {"median_charged_s_per_attempt": [157.7, 227.6, 8.1], "charged_s_per_succeeding_complex": [252, 443, 37], "gpu_s_per_success": {"AutoDock": 173, "DiffDock": 442},
     "autodock_cheaper_than_diffdock_on_this_denominator": True, "cpu_more_expensive": True},
    EFS + " wall_distribution.medians_s (unchanged); per-success = totals from " + EFSUM + " divided by n_complexes with near2 > 0 (214 / 183 / 83): (151.4192/32 + 10.2699) h * 3600 / 214 = 252.4; (1.3221/32 + 22.4609) * 3600 / 183 = 442.7; (21.3922/32 + 0.1817) * 3600 / 83 = 36.9; GPU 10.2699*3600/214 = 172.8, 22.4609*3600/183 = 441.9 (recipe reproduces canonical 271 / 479 / 38 and 186 / 478)", False, "")
add("M:768", "led DiffDock on near-native recovery and native-pocket reach at every depth; reach never separated them; 4.2-point DiffDock lead to 9.4-point AutoDock lead",
    {"rank1_interval": [-1.7, 12.8], "reach_lead_every_depth": "holds (78/86/88/90 vs 73/83/85/86 %)", "reach_never_separated": "holds (paired exact McNemar AD vs DD reach top-1 50/35 p 0.128, top-5 37/27 p 0.26, top-10 33/24 p 0.29, top-15 31/20 p 0.161)"},
    RC + "reference_convention_margins.csv; reach: " + CL + "topN_crystal_cluster_matrix_stats.txt and crystal_cluster_homogeneity_stats.json A_reaches (top-15); other depths recomputed from " + CL + "per_complex_summary.csv with stats_utils.mcnemar_exact", True,
    "W21: interval 1.7 / 12.8; both reach clauses hold, so the bracketed clauses of W21 keep their current wording.")
add("M:777", "removes at most four of 303 complexes at any depth", {"max": 5}, HUB + "validity_gate_cost.csv", True, "W22 'five'")
add("M:779", "fixed-rank near-native rate barely changed", {"diffdock_smina_minus_raw_pp_by_rank_1_5_15_30": [1.65, 0.0, 1.33, 0.36]}, HUB + "optimization_benefit_stats_pairwise.csv (smina_minus_raw_pp)", False, "holds")
add("M:781", "133 to 132; margin survives only in the three-way conjunction; dropping either term leaves the arms exactly level; gnina level or ahead on every component; separates them at no depth; 83 vs 61 at top-15 survives Holm; selection margin nine",
    {"top15_triple": [144, 143], "top15_double_gate_smina_gnina": [179, 177], "top15_near_native_only": [182, 180], "pooled_near_native": [186, 186], "pooled_form": [236, 235],
     "smina_vs_gnina_gated_mcnemar_p_by_depth_1_5_15_30": [0.69, 1.0, 0.73, 1.0], "separates_at_no_depth": True,
     "equibind_near_native_top15": {"gnina": 86, "smina": 64, "mcnemar_p": "1.1e-4", "survives_holm": True},
     "equibind_selection_top15": {"gnina": 57, "smina": 48, "margin": 9, "mcnemar_p": 0.122}},
    SIG + " prose.refiners.* rows (diffdock_triple/double/near_top15, equibind_near_top15 86/64 p 0.0001131, equibind_selection_top15 57/48 p 0.1221); depth-wise smina vs gnina: " + RECOMP, True,
    "W23 wording. DECISION: the plan row 781 says the EquiBind selection margin becomes 'seven'; the rebuilt table gives 57 - 48 = 9, so 'nine complexes' is UNCHANGED (the plan's 'seven' is the DiffDock discordant count of :197).")
add("M:788", "about eight percentage points for EquiBind; six in ten; nine in ten",
    {"equibind_gain_1_15_pp": 8.6, "wording_options": ["about nine", "between eight and nine"], "leading_pair_failures_no_pose_pct": [57.8, 70.2], "equibind_pct": 89.4},
    TOPK, True, "W24a")
add("M:790", "indistinguishable on form; lead from top-5", {"form_p_holm_by_depth": [1.0, 0.899, 1.0]}, KAB, False, "holds")
add("M:792", "pocket in the pose cloud for 95 % of complexes; best inexpensive rule 57 %; AutoDock higher at top-15 but paired reach test did not separate",
    {"oracle_cluster_pct": 98.0, "best_rule_precision_at_1_pct": 71.3, "top15_reach_pct": [90, 86], "paired_reach_mcnemar_p_top15": 0.161},
    CL + "summary.json ranking_ablation (oracle 0.9802, ntools 0.7129); crystal_cluster_homogeneity_stats.json A_reaches", False, "'98 %' and 'about 71 %'.")
add("M:806", "Orai1 qualitative", {}, "plan", False, "hold")
add("M:817", "EquiBind shortest, then DiffDock, then AutoDock; 80 of 303 against 169; AutoDock most expensive on charged and hardware axes",
    {"ordering_holds": True, "medians_charged_s": [2.4, 39.0, 91.8], "equibind_success": 83, "diffdock_success": 183}, EFQ + " tiers.near2", False, "numbers 83 / 183")
add("M:824", "qualitative", {}, "plan", False, "hold")
add("M:834", "36.6 / 34.0 / 18.2 %; 65.3 / 55.1 / 26.4 %; +2.6 (-4.2, +9.4); +13.2 top-5 p_holm 7.6e-4; +10.2 top-15; 90 % compatibility interval within about eight points",
    {"rank1_pct": [49.2, 43.6, 18.8], "top15_pct": [70.3, 59.1, 27.4], "rank1_diff_pp": 5.6, "ci95": [-1.7, 12.8], "top5_pp": 13.2, "top5_p_holm": "7.6e-4", "top15_pp": 11.2,
     "newcombe90_pp": [-0.5, 11.7], "compatibility_wording": "within about twelve points", "tost_flags": {"10pp": False, "12pp": True, "15pp": True}},
    TOPK + "; " + BC + " (equiv_10pp False, equiv_12pp True); 90 % Newcombe recomputed with stats_utils.newcombe_paired_diff_ci(z = 1.645) on the rank-1 vectors (" + RECOMP + ")", True, "W25 'twelve points'.")
add("M:836", "did not separate (rank-1 native-interaction); DiffDock led only for the top-five union", {"rank1_f1_shared_280_wilcoxon_p": 0.064, "top5_union_DD_vs_AD_p_common_264": 0.011}, "as M:489 / M:493", False, "holds")
add("M:843", "0.1-0.7 pp for DiffDock, 4.1-4.6 EquiBind gnina, 3.5-4.4 smina; 31.0 % to 36.6 %; does not clear Holm; significance from top-5",
    {"band_gains_pp_1-10_11-20_21-30": {"DiffDock gnina": [0.6, 0.3, 1.0], "DiffDock smina": [0.5, 0.4, 0.8], "EquiBind gnina": [4.7, 4.5, 4.4], "EquiBind smina": [4.5, 3.6, 3.9]},
     "band_wording": "0.3-1.0 pp for DiffDock with gnina (0.4-0.8 with smina), 4.4-4.7 for EquiBind with gnina and 3.6-4.5 with smina",
     "autodock_rank1_gated_pct": [42.9, 49.2], "gated_discordant": [20, 39], "gated_p_raw": 0.018,
     "holm_table19_family_16_tests_near_native_no_validity": {"rates": [43.2, 49.8], "discordant": [19, 39], "p_raw": 0.012, "p_holm": 0.060, "verdict": "does not clear"},
     "holm_appC_seven_contrast_family": {"exh128_raw_vs_gnina_rank1_RMSD2": {"p_raw": 0.012, "p_holm": 0.107, "family": "the exhaustiveness report section 9 family (three rungs x two gates)"}},
     "gain_from_top5": "resolved at every depth from top-5 (Table 19 family Holm 4.3e-7 / 3.0e-4 / 4.0e-4 at k = 5 / 10 / 15)"},
    "bands recomputed from " + HUB + "optimization_benefit_by_rank.csv with the Table 25 recipe (pose-count-weighted pooled rate per band, opt minus raw; the recipe reproduces the printed canonical +0.6/+0.7/+0.1/+0.1/+0.5/+0.4 and 4.4/4.6/3.5/4.4/3.7/4.1); gated 42.9->49.2: " + HUB + "validity_gate_cost.csv (130 -> 149) and " + RECOMP + "; sixteen-test Holm: best-of-top-k near-native (no validity) for AutoDock raw->gnina, DiffDock raw->smina, raw->gnina, EquiBind raw->gnina at k = 1/5/10/15, stats_utils.mcnemar_exact + holm (" + RECOMP + "); seven-contrast family: " + EX + " section 9 (exh128 RMSD<=2 rank-1 131 -> 151, 39/19, p 0.012, Holm 0.107)",
    True, "D14: name the family. NOTE for the editor: under the REBUILT exhaustiveness report the rank-1 exh128 raw->gnina contrast has Holm 0.107 in that report's own six-contrast rescoring family and 0.060 in Table 19's sixteen-test family, so 'does not clear Holm' holds in BOTH families printed by the rebuilt sidecars; the plan's D14 preview (App. C seven-contrast Holm 0.037 survives) refers to the ladder-vs-selected family, which the appendix registry must confirm before the M:843 sentence claims a survival.")
add("M:848", "2.4 s / 49.4 s / 137.2 s; 169 / 199 / 80; twentyfold; quarter",
    {"median_charged_s_per_qualifying_pose": {"EquiBind": 2.4, "DiffDock": 39.0, "AutoDock": 91.8}, "succeeded_complexes": {"DiffDock": 183, "AutoDock": 214, "EquiBind": 83},
     "equibind_over_diffdock_fold": 16.2, "wording": "sixteenfold", "quarter_of_set": "83 / 303 = 27 % (holds)"},
    EFQ + " tiers.near2.per_method", True, "W27: 'twentyfold' -> 'sixteenfold'; ordering unchanged so the ordering sentences at :850 / :744 / :758 / :817 / Abstract / Kurzfassung keep their wording.")
add("M:850", "most expensive in charged, CPU-core and GPU seconds; 4.73 / 10.27 / 15.00 h; 68.5 %", {"ordering_holds": True, "unchanged": [4.73, 10.27, 15.00, 68.5]}, EFSUM + " (cpu_core_s_per_near2, gpu_s_per_near2)", False, "no change")
add("M:857", "six in ten; nine in ten; 87 and 64 additional complexes; 25 for EquiBind",
    {"wording": "six to seven in ten", "equibind_nine_in_ten": 89.4, "selection_gains_1_15": [64, 47, 26]}, TOPK, True, "W24b")
add("M:866", "winning margins eight / one / nine; 3.0 percentage points; winner's curse between 0.4 and 0.8",
    {"margins_top15_triple": {"AutoDock": 9, "DiffDock": 1, "EquiBind": 9}, "largest_margin_pp": 3.0,
     "winners_curse_pp": {"AutoDock k=2 (gnina vs raw exh128, discordant 10/1)": 0.44, "AutoDock k=3": 0.66, "DiffDock k=3 (smina vs gnina, 4/3)": 0.52, "EquiBind k=3 (gnina vs smina, 18/9)": 1.03},
     "range_wording_candidates": "between 0.4 and 1.0 (AutoDock k=2) or between 0.5 and 1.0 (k=3)"},
    SC + " (175 vs 166, 144 vs 143, 57 vs 48); winner's curse recomputed with the memory recipe (variant-selection-bias-quantified.md, 2026-08-24 basis: sigma_d = 100 * sqrt(n10 + n01) / 303 from the top-15 triple-gate winner-vs-runner-up discordance, bias = E[max of k standard normals] * sigma_d / sqrt(2)); discordances from " + RECOMP, True,
    "W28 'nine'. DECISION LEFT OPEN: the AutoDock family at exh128 has two arms (raw, gnina), the learned pipelines three; report 0.4-1.0 if k = 2 is accepted for AutoDock. No sidecar carries this estimator.")
add("M:868", "78 of 303; +6.7 pp on 225 metal-free vs +20.5 on 78 metal-adjacent",
    {"n_metal_adjacent": 78, "n_metal_free": 225, "top15_contrast_metal_free_pp": 8.4, "metal_free_counts": [160, 141], "metal_free_discordant": [54, 35], "metal_free_p": 0.056,
     "top15_contrast_metal_adjacent_pp": 19.2, "metal_adjacent_counts": [53, 38], "metal_adjacent_discordant": [28, 13], "metal_adjacent_p": 0.028,
     "between_strata_fisher_p": 0.44, "cofactor_ion_split_registered_rule": [16, 62], "convention_disclosure": {"multi_copy_308": 143, "multi_copy_303": 138, "median_nearest_other_copy_A_over_143": 36.0}},
    MS + "metal_stratum_contrasts.csv (convention nearest, cutoff 5.0, depth 15, contrast autodock vs diffdock rows; diff_pp is D - A); " + MS + "metal_stratum_summary.json (cofactor_ion_split.5A n_cofactor 16 n_ion 62); 143 / 138: " + RC + "reference_convention_summary.json; 36.0 A median: plan Data facts (not in a rebuilt sidecar)", True,
    "B3 text. FLAG (D4): the registered rule splits the 78 into 16 cofactor / 62 free-ion where App. C :165 prints 15 / 63. The 36.0 A median is the only number here without a rebuilt sidecar.")
add("M:870", "qualitative", {}, "plan", False, "hold")

out = Path("Scripts/Analysis/nearest_copy_program/values_registry_main.json")
out.write_text(json.dumps(R, indent=1, ensure_ascii=False))
print(out, len(R))
