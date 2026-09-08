"""Compare docked poses (AutoDock Vina, DiffDock, EquiBind) against the
crystallographic reference poses of the PoseBusters Benchmark Set.

Two distinct analyses are produced:

  PART A — Oracle (best-pose) comparison, all three tools
    For each method, the pose closest to the crystal is selected per
    receptor-ligand pair (oracle / best-of-N).  This shows each method's
    ceiling performance independent of any ranking.

  PART B — Ranking quality, AutoDock Vina & DiffDock only
    Compares each tool's top-ranked pose (rank 1, as scored by Vina energy /
    DiffDock confidence) against the oracle pose.  The gap between the two
    quantifies how well the ranking mechanism works.  EquiBind is excluded
    from Part B because it produces no canonical pose ranking.

EquiBind variant split (toggle with --no-split-equibind)
    EquiBind emits several flavours of pose, all tagged "equibind_guided" in the
    raw CSV. Every variant axis is unpacked into separate methods so each flavour
    gets its own oracle/cross-tool curves. The axes are read from the provenance
    columns the PoseBusters CSV now carries (pocket_source / refine_variant /
    clamp_variant, from each pose's SDF tags), with a pose-filename fallback:
      pocket   unguided (blind) | fpocket | p2rank  (fpocket & p2rank are
               reported as DISTINCT methods, not merged into "guided")
      refine   __refSMINA -> "smina";  __refGNINA -> "gnina";  __refRAW -> "raw" control
      clamp    __clampON / __clampOFF                    (centroid-clamp study)
    e.g. fpocket_raw_p001_pose01__refGNINA.sdf -> "equibind_fpocket_gnina". Axes
    not present collapse to "equibind_fpocket" / "equibind_p2rank" /
    "equibind_unguided". All EquiBind variants stay oracle-only (Part B is still
    AutoDock/DiffDock, which alone rank poses). EquiBind has no native ranking, so
    its ranked figures order poses by gnina AFFINITY (gnina_affinity, most-negative =
    rank 1) — read from each gnina pose's SDF tag. The per-pose provenance
    (pocket_source, clamp_variant, refine_variant, smina_affinity, gnina_affinity) is
    also carried into per_pose_metrics.csv. NOTE: refine_mode='on' rewrites the pose
    in place with no suffix/tag, so the refine split needs refine_mode='both'.

Best-variant filters (--best-equibind-only / --best-diffdock-only / --best-variants-only)
    Collapse the many EquiBind variants (and/or the three DiffDock optimizer
    variants) to the single best-performing one each (highest PB-Valid AND RMSD ≤
    2 Å = oracle_pb_valid_and_rmsd2_%). The individual family filters leave other
    methods untouched. The --best-variants-only umbrella restricts ordinary
    summaries and headline cross-tool plots to raw AutoDock plus the DiffDock* and
    EquiBind* winners. Explicit variant/optimization diagnostics retain the full
    frame, and per_pose_metrics.csv always carries every variant for drill-down.

Reference-ligand convention (--reference-convention {instance,nearest})
    Every crystal-referenced metric below is measured against ONE deposited copy
    of the ligand per pose. Under ``instance`` (the default, reproduces the frozen
    tables bit-for-bit) that copy is always the reference instance, record 0 of
    <ID>_ligands.sdf (== <ID>_ligand.sdf). Under ``nearest`` — the PoseBusters
    paper / ``check_rmsd`` convention — it is the copy of <ID>_ligands.sdf with the
    smallest symmetry-corrected in-place RMSD to the pose (NaN counts as +inf,
    record 0 breaks ties), and that SAME copy j* drives centroid, PoseBusters
    RMSD/Kabsch, best-fit, torsions, contact and PLIF recovery for the pose. The
    primary columns keep their names; ``*_ref_instance`` twins, ``n_copies``,
    ``ref_copy_index``, ``nearest_copy_index``, ``nearest_copy_is_ref`` and
    ``reference_convention`` are appended to per_pose_metrics.csv.

Metrics computed per pose:
    * Symmetry-corrected heavy-atom RMSD vs. the reference copy of the crystal
      ligand (no superposition; see the convention above)
    * Centroid distance to that copy (translation — how far it "moved")
    * Rigid-body rotation angle of the best-fit onto crystal (how it "turned")
    * Best-fit (superposed) RMSD — PoseBusters' kabsch RMSD; self-computed fallback
    * TFD + per-rotatable-bond torsion deviations & #flipped (how it "twisted")
    * UFF strain energy (predicted − ETKDG/UFF minimum)
    * Steric clash count vs. protein heavy atoms (vdW-based)
    * Native contact recovery (residue-level Jaccard within 4 Å)
    * PLIF Tanimoto similarity vs. crystal IFP (optional, requires ProLIF + obabel)

Output plots (in --out-dir):
    01_oracle_rmsd_cdf.png           — Oracle RMSD CDF, all tools
    02_oracle_rmsd_boxplot.png       — Oracle RMSD boxplot, all tools
    03_oracle_vs_top1_success.png    — Oracle vs Top-1 success bars
    04_rank_success_curve.png        — RMSD ≤ 2 Å success at each rank k
    04b_plif_recovery_by_rank.png    — Interaction-profile similarity to the crystal at
                                       each rank k: ProLIF IFP Tanimoto + native-contact
                                       recovery. AutoDock/DiffDock on native rank,
                                       EquiBind on gnina-affinity rank (benchmark only)
    05_cumulative_oracle_curve.png   — Best-of-top-k success as k grows
    06_top1_vs_oracle_scatter.png    — Top-1 RMSD vs oracle RMSD per pair
    07_oracle_rank_histogram.png     — Which rank does the oracle pose have?
    08_cross_tool_oracle_pairwise.png — Cross-tool oracle RMSD scatter matrix
    09a_accuracy_vs_validity_top1.png — PoseBusters-paper Fig.1-style bars for the
                                       TOP-1 (rank-1) pose (ranking tools): %RMSD ≤ 2 Å
                                       vs %(RMSD ≤ 2 Å & PB-valid)
    09b_accuracy_vs_validity_oracle.png — per variant, the RANK-1 pose (EquiBind: first
                                       generated pose) beside the ORACLE (best-of-N) pose,
                                       with each ML tool's RAW variant next to its best
                                       (DiffDock raw vs DiffDock*, EquiBind raw vs
                                       EquiBind*); light = RMSD ≤ 2 Å, dark = also PB-valid
    09c_optimization_raw_vs_best.png — grouped bars contrasting each tool's RAW
                                       representative pose with its BEST post-hoc-
                                       optimised variant (Vina: rank-1, no opt step;
                                       DiffDock: raw vs smina rank-1; EquiBind: first
                                       generated pose vs gnina-ranked) — shows how
                                       post-hoc optimisation recovers PB-validity
    09d_rank1_vs_topn.png            — grouped bars: best-of-top-d pose at several
                                       depths (top-1, top-15, top-30) for BOTH the raw
                                       and the smina-optimised run of each ML tool
                                       (ranking headroom + optimisation gain in one
                                       view; top-30 ≈ full oracle). Vina (confidence
                                       rank), DiffDock raw & smina-opt (confidence rank),
                                       EquiBind raw (generation order) & gnina-opt
                                       (gnina-affinity rank). Shows more poses raise
                                       accuracy but only
                                       optimisation recovers PB-validity
    09d_alt_scatter.png              — same data, non-bar: accuracy×validity scatter
                                       with the y=x "fully valid" diagonal; each variant
                                       a top-1→…→top-30 trajectory (gap below diagonal =
                                       accurate-but-invalid)
    09d_alt_slopegraph.png           — same data: accuracy | validity slopegraph across
                                       depths (slope = ranking headroom; flat deep
                                       segment = ceiling reached)
    09d_alt_dumbbell.png             — same data: per (variant, depth) validity●——○accuracy
                                       dumbbell; the segment is the accurate-but-invalid share
    09f_pbvalid_yield_boxplot.png    — "Search Effort, Post-Hoc Optimization and PoseBusters Validity":
                                       box/whisker of the PoseBusters-valid YIELD of the
                                       PRODUCED poses per variant (raw beside its smina /
                                       gnina optimized variants, grouped by tool). Each
                                       box = distribution across complexes of the % of a
                                       variant's poses that pass PoseBusters; smina/gnina
                                       refine the SAME poses (fixed denominator), so the
                                       raw→opt shift is recovered validity. Paired Wilcoxon
                                       signed-rank brackets (Holm across all raw→opt
                                       contrasts) + Friedman omnibus show the gains are
                                       significant. Tables: pbvalid_yield_by_variant.csv
                                       (per-variant pooled + per-complex distribution) and
                                       pbvalid_yield_per_complex.csv (the box source).
                                       09f_pbvalid_yield_report.txt holds the FULL stats
                                       (Friedman + Wilcoxon pairwise + raw→opt contrasts
                                       with effect sizes / CIs / Holm p) and the figure
                                       caption/notes (moved off the plot)
    09f_pbvalid_yield_boxplot_all_variants.png — same yield box/whisker but showing EVERY
                                       variant present in df_full (all EquiBind pocket ×
                                       refine [× clamp] variants + all DiffDock variants +
                                       AutoDock), independent of the best-variant / collapse
                                       flags; coloured per variant, no raw→opt brackets.
                                       Table: pbvalid_yield_by_variant_all.csv
    10_pb_valid_rmsd2_success_bars.png — Per-method success bars for the strict
                                       "≤ 2 Å & PB-valid" category alone (% of ALL
                                       complexes): Oracle (all tools) + Top-1
                                       (ranking tools)
    11_vs_posebusters_paper.png      — THIS study vs the PoseBusters paper
                                       (Benchmark set, primary pose): %RMSD ≤ 2 Å
                                       and %(RMSD ≤ 2 Å & PB-valid), our top-1
                                       (EquiBind*: best pose) beside the paper's
                                       published numbers (POSEBUSTERS_PAPER_BENCHMARK)
    12_rmsd2_vs_pbvalid_grouped.png  — The two criteria side by side per method
                                       (Oracle + Top-1), with the Δ pp drop = how
                                       much PB-validity costs each tool
    13_pbvalid_filter_influence.png  — Dumbbell per method (oracle): gap between
                                       RMSD ≤ 2 Å and RMSD ≤ 2 Å & PB-valid, with
                                       Δ pp and % of accurate poses retained
    14_twist_turn.png                — "Twisted & turned": per measure, box plots of
                                       how each tool re-orients/re-bends the ligand vs
                                       its crystal pose — translation (moved), rigid-
                                       body rotation (turned), best-fit RMSD / TFD /
                                       torsion deviations / strain (twisted); top-1 vs
                                       oracle pose, all tools
    15_oracle_rank_distribution.png  — How often the n-th ranked pose is the closest
                                       (oracle) pose: P(oracle at rank k) bars + the
                                       cumulative recovery within top-k (ranking tools;
                                       RMSD only — ignores PB-validity)
    15b_topk_recovery_validity.png   — Cumulative recovery within top-k with PB-validity
    15c_optimization_benefit_by_rank.png — Per-rank benefit of smina/gnina optimization
                                       vs the raw pose, DiffDock + EquiBind × (near-native,
                                       PB-valid, both); rank groups 1-10/11-20/21-30 shaded
    15d_optimization_benefit_by_group.png — Same as 15c but the +pp impact summarised as
                                       grouped bars over rank groups 1-10 / 11-20 / 21-30
                                       made explicit: per variant a solid (near-native,
                                       RMSD ≤ 2 Å) and dashed (near-native AND PB-valid)
                                       line — AutoDock, DiffDock raw & smina-opt, EquiBind
                                       raw (generation order) & gnina-opt (gnina affinity);
                                       the solid-dashed gap = accurate-but-invalid
    16a_pocket_targeting_by_rank.png — Pocket targeting (ranking tools): per rank
                                       k=1..N, the % of those rank-k poses in the
                                       validated (crystal) pocket plus the cumulative %
                                       of complexes with any of ranks 1..k in it (a pose
                                       counts as in-pocket when its centroid is within
                                       --pocket-cutoff Å of the crystal ligand centroid)
    16b_pocket_localization_summary.png — Per-tool summary (ranking tools): % of all
                                       top-N RANKED poses in the validated pocket beside
                                       the % of complexes whose ORACLE pose is in a
                                       DIFFERENT pocket than the top-N ranked set
    17_pb_test_waterfall.png         — PoseBusters paper-style failure waterfall, one
                                       panel per method: starting from every top-ranked
                                       pose (top-1 for ranking tools, best/oracle for
                                       unranked EquiBind), drop RMSD > 2 Å then each
                                       canonical PoseBusters test in turn — each red bar
                                       is the poses that test removes; green = passing all
    17b_pb_test_waterfall_topn.png   — Same cascade, but for ranking tools the
                                       representative pose is the best (lowest-RMSD) of the
                                       tool's top-N ranked poses instead of just rank-1
                                       (one pose per complex; EquiBind unchanged)
    17c_pb_waterfall_top1_vs_topn.png — Per ranking tool, the two cascade survival curves
                                       overlaid (rank-1 vs best-of-top-N): the shaded band
                                       is the complexes recovered by considering more
                                       ranked poses, with Δ callouts at RMSD ≤ 2 Å and at
                                       passing all tests. Includes a gnina-ranked EquiBind
                                       panel (its poses ordered by gnina affinity, since
                                       EquiBind has no native rank)
    17c2_pb_waterfall_top1_vs_top30.png — Same ranking-headroom comparison but against the
                                       best of the top-30 poses (a near-oracle ceiling, as
                                       the tools emit ~30 poses per complex)
    17f_pb_top30_recovery_by_test.png — The top-30 recovery in COUNTS, per cascade step: a
                                       bar at each step = extra complexes surviving with
                                       best-of-top-30 vs rank-1, colour-split into the
                                       RMSD ≤ 2 Å placement filter (blue) and the other
                                       physical-validity PoseBusters tests (amber), so the
                                       two categories are read off per test. Bars are drawn
                                       only where the top-30 pick makes a difference
    17c3_pb_waterfall_top1_vs_top30_raw.png / 17g_pb_top30_recovery_by_test_raw.png — The
                                       17c2 / 17f analyses repeated for the UN-REFINED
                                       variants (raw DiffDock, raw EquiBind geometry ranked
                                       by gnina score; AutoDock unchanged), for a direct
                                       raw-vs-refined comparison
    17h_pb_top{15,30}_recovery_raw_vs_refined.png — 17f + 17g in ONE figure: AutoDock once,
                                       then each tool's raw panel directly above its refined
                                       panel (DiffDock raw vs its collapsed refined variant —
                                       smina or gnina per --collapse-diffdock-variant;
                                       EquiBind raw/gnina), so the
                                       refinement benefit reads off panel-to-panel. Rendered at
                                       TWO pool depths — best-of-top-15 (= --top-n) and
                                       best-of-top-30 vs rank-1 — as two standalone figures.
                                       17f/17g/17h each annotate per panel a Wilson 95% CI on
                                       the pass-all rate + exact McNemar (Holm-adjusted) for
                                       the rank-1→top-N gain, and raw→refined McNemar on the
                                       refined panels; full payloads in the *_stats.json sidecar
                                       (keys top{15,30}_recovery_raw_vs_refined). Both 17h depths
                                       share ONE combined companion report
                                       17h_pb_recovery_raw_vs_refined_report.txt covering all
                                       THREE selections (top-1 · best-of-top-15 · best-of-top-30):
                                       the full paired-test details (a separate Holm family per
                                       figure) and a ';'-separated per-stage attrition table
                                       (poses failed at each cascade step + %, one selection block
                                       per depth)
    17d_pb_test_waterfall_raw_vs_best.png — Same cascade, top-1 pose, comparing the RAW pose
                                       with its best refined variant per method:
                                       DiffDock raw vs diffdock_smina (rank-1), EquiBind raw
                                       vs gnina (same prediction, ranked #1 by its gnina
                                       score), AutoDock once as a physics-docking reference;
                                       optimised panels annotate the net gain in poses
                                       passing every test vs their raw counterpart
    17e_pb_test_waterfall_raw_vs_best_validity.png — Same raw-vs-best panels, but WITHOUT
                                       the RMSD ≤ 2 Å gate: every ranked-1 pose faces the
                                       PoseBusters tests regardless of placement, exposing
                                       where poses fail on physical-validity grounds alone
                                       (the RMSD filter otherwise removes most poses first)
    17i_pb_rank_depth_added_pose_waterfall.png — The deepening question, per method (raw +
                                       best refined variant), WITHOUT the RMSD ≤ 2 Å gate:
                                       instead of one pose per complex, it pools the poses
                                       NEWLY considered when the rank cut widens from 1 to
                                       top-N (ranking positions 2..top-N; EquiBind gnina-
                                       ranked) and runs the PoseBusters cascade over that
                                       pool. Each panel's headline states the payoff — how
                                       many poses/complexes the deeper cut adds and how many
                                       complexes it recovers (rank-1 not PB-valid → an added
                                       pose PB-valid); red bars = added poses each test
                                       removes, green = added poses passing every test. CSVs
                                       pb_rank_depth_added_pose_waterfall.csv (cascade steps)
                                       + pb_rank_depth_added_pose_summary.csv (per-panel
                                       payoff + independent per-test fail counts).
    17i2_pb_rank_depth_added_pose_failbytest.png — Companion to 17i: of those rank-1→top-N
                                       ADDED poses, how many fail EACH PoseBusters test as an
                                       INDEPENDENT tally (a pose counted for every test it
                                       fails, so later tests aren't masked by earlier ones as
                                       in the cascade). Horizontal bars per method + count/%;
                                       CSV pb_rank_depth_added_pose_failbytest.csv.
    18_topn_within_thresholds_complex.png / _pose.png — Fine-grained RMSD sweep (split
                                       into two graphs): how many of the top-N ranked poses
                                       of AutoDock Vina / DiffDock (native rank) / EquiBind
                                       (gnina-affinity ranked) land within 1, 1.25, 1.5 … Å
                                       (--fine-rmsd-thresholds). _complex — % of complexes
                                       whose rank-1 pose (solid) / any of the top-N poses
                                       (best-of-top-N, dashed) is within each threshold;
                                       _pose — % of ALL pooled top-N ranked poses within
                                       each threshold. Denser than the 2 Å success line so
                                       the sub-2 Å accuracy of the top poses shows.
    18_topn_within_thresholds_pbvalid_complex.png / _pose.png — Same sweep, but a pose
                                       counts as a hit only if it is BOTH within the threshold
                                       AND PoseBusters-valid (combined near-native + physically
                                       valid, matching oracle_pb_valid_and_rmsd2_%). Denominators
                                       unchanged (all complexes / all top-N poses).
    18_topn_within_thresholds_kabsch_complex.png / _pose.png (+ _kabsch_pbvalid_*) — The
                                       same two sweeps, but distance is the best-fit (Kabsch)
                                       RMSD (bestfit_rmsd) instead of the as-placed RMSD: each
                                       pose is optimally superposed onto the crystal first, so
                                       the threshold measures conformer/shape fidelity only and
                                       ignores where the pose was docked. Companion CSVs
                                       topn_within_thresholds_kabsch[_pbvalid].csv. Both Kabsch
                                       figures use a HALVED grid (0–2.5 Å = half the 0–5 Å
                                       standard grid), since superposition collapses distances
                                       toward zero; the PB-valid one additionally emphasises a
                                       1 Å success line (form fidelity is a sub-Ångström question).
                                       NOTE (all fig-18 graphs): the legend now sits between the
                                       title and the plot, each tool's complex count n is stated
                                       once in the title, and the caption + full complex-panel
                                       statistics are written to a companion *_complex_report.txt
                                       / *_pose_report.txt (only the Wilson-CI whiskers stay on
                                       the plot).
    18_refinement_gain_vs_depth_pbvalid.png / _kabsch_pbvalid.png — Companion to the two PB-valid
                                       complex figures: how the raw→gnina-refined benefit evolves
                                       as the top-N pool deepens (N = 1, 2, 3 … top_n). Top panel
                                       = best-of-top-N success (% complexes within the figure's
                                       success line — 2 Å for RMSD, 1 Å for Kabsch — AND PB-valid)
                                       for the raw (dashed) vs gnina-refined (solid) variant of
                                       DiffDock / EquiBind, plus AutoDock Vina as a single native
                                       curve (no refinement step — its rise with N is the pure
                                       rank-depth gain, the reference the refined tools are read
                                       against); bottom panel = the refinement difference (gnina-
                                       refined − raw, pp) vs depth for the refined tools. CSVs
                                       refinement_gain_vs_depth[_kabsch]_pbvalid.csv + *_report.txt.
    18_topn_within_thresholds_pbvalid_depths.png / _kabsch_pbvalid_depths.png — Multi-depth
                                       version of the complex panel: the best-of-top-N curve
                                       at THREE pool depths (top-1 = rank-1, top-15, top-30)
                                       across the threshold sweep, for near-native (as-placed)
                                       RMSD and Kabsch RMSD respectively. Colour = tool, line
                                       style = depth. The value at the emphasised threshold
                                       (2 Å RMSD / 1 Å Kabsch) is called out on each curve; the
                                       full sweep of numbers lives in topn_within_thresholds_
                                       [kabsch_]pbvalid_depths.csv + *_report.txt. The report
                                       also carries the paired significance tests AT that
                                       threshold — Cochran's Q + pairwise McNemar (Holm)
                                       between tools per depth, and the nested ranking-headroom
                                       (recovered share + Wilson CI) across depths within each
                                       tool (sidecar keys within_by_depth[_kabsch]_pbvalid).
    19_within2_validity_dumbbell.png — Dumbbell per method: PB-valid share of near-native
                                       (RMSD ≤ 2 Å) poses over EVERY generated pose (circle) vs
                                       the oracle / min-RMSD pick only (diamond); each row labels
                                       its complex / total-pose / near-native counts. Shows
                                       whether the RMSD-greedy oracle pose is as physically valid
                                       as a typical near-native one. Rows are grouped tools →
                                       blind EquiBind → pocket-guided EquiBind, the guided
                                       fpocket/p2rank kept as separate rows but shaded + bracketed
                                       together (handed the pocket, so read apart). CSV
                                       (within2_validity_comparison.csv) has every variant.
                                       Complements fig 13 (oracle-only).
    20_form_fidelity.png             — Form fidelity of the successful poses (RMSD ≤ 2 Å
                                       AND PB-valid). The ≤ 2 Å line is scored on the
                                       IN-PLACE RMSD, which mixes placement and shape; this
                                       isolates the shape. (A) box of best-fit (Kabsch) RMSD
                                       per method (form error, placement removed); (B) % of
                                       those poses whose form is correct (best-fit ≤
                                       --form-ok-kabsch, default 1 Å); (C) 100 %-stacked split
                                       of each method's mean-square deviation into form vs
                                       placement — is the residual error conformation- or
                                       positioning-limited; (D) per-pose in-place RMSD vs
                                       best-fit RMSD, coloured by tool family (points on the
                                       x-axis = form perfect / placement-limited, points on
                                       the y = x bound = form-limited). Crystal-free sets
                                       (Orai) have no RMSD → figure/CSV are empty (no-op).
    20_form_fidelity__rank_depth_{pbvalid,within2}__*.png — Companion to fig 20: form fidelity
                                       across ranking depth (top-1 / top-5 / top-15 / top-30), as
                                       FOUR standalone graphs (__form_error box, __form_correct %,
                                       __mechanism composition, __form_vs_rank trend). Pooled
                                       per-pose, nested cohorts; ranking = native rank for
                                       AutoDock/DiffDock, gnina-affinity rank for EquiBind. Rendered
                                       under two selections for comparison: _pbvalid = all PB-valid
                                       poses (RMSD gate REMOVED), _within2 = PB-valid AND RMSD ≤ 2 Å.
                                       Summaries: form_fidelity_summary__rank_depth_{pbvalid,within2}.csv.
                                       Pose-level, so NOT comparable to the per-complex ≤ 2 Å summaries.
                                       Crystal-free sets (Orai) → no-op.
    20_form_fidelity__rank_depth_gate_impact__{form_error,form_correct}.png — Joint view of the two
                                       selections (lines only) so the IMPACT of adding the RMSD ≤ 2 Å
                                       criterion is the gap between paired curves (dashed = valid-only,
                                       solid = valid AND ≤ 2 Å) per method × depth, as two standalone
                                       graphs: __form_error (median form error) and __form_correct
                                       (% form-correct + ≤ 2 Å pose-retention per depth). Every point
                                       is value-labelled at each rank.
    20b_form_vs_success_count.png    — Method-level: NUMBER of ≤ 2 Å & PB-valid poses (x,
                                       log) vs their FORM error (y, median best-fit RMSD,
                                       IQR whiskers), one marker per method. Do tools that
                                       produce more successes also produce better-form ones?
    20d_form_vs_placement_by_family.png — Fig 20 panel D split into individual per-family
                                       graphs (AutoDock / DiffDock / EquiBind + combined),
                                       each coloured by mechanism region (r = form²/in-place²:
                                       placement-limited / mixed / form-limited) with the
                                       y = x bound and the r = 1/3, 2/3 rays.
    20e_form_vs_placement_clustering.png — Diagnostic: GMM (BIC-selected k) / k-means /
                                       HDBSCAN on the panel-D cloud, showing it has NO natural
                                       clusters (GMM tiles a gradient, HDBSCAN → mostly noise),
                                       which is why 20d segments by mechanism ratio instead.
                                       Skipped if scikit-learn is unavailable.
    20c_oracle_selection_comparison.png — Head-to-head of the two oracle selections:
                                       (A) complexes counted as a success per method under the
                                       nearest rule vs the validity-constrained rule (rescued
                                       count in red); (B) their median form error under each.
    20f_mechanism_examples_3d.png    — Didactic 3-D illustration (NOT from the data) of the three
                                       mechanism labels on a synthetic model ligand: top row each
                                       docked pose vs the crystal AS DOCKED (in-place RMSD), bottom
                                       row after best-fit superposition (form/Kabsch RMSD + r).
                                       Placement-limited collapses onto the crystal, form-limited
                                       stays split, mixed collapses partway. Also writes
                                       mechanism_examples.pdb (+ .pml) for a PyMOL screenshot:
                                       chain A = crystal, chain B = docked pose, resi 1/2/3 =
                                       placement/mixed/form, laid out side by side.
    20g_mechanism_examples_real_3d.png — Same idea but with REAL docked poses selected from
                                       per_pose_metrics.csv (one clean PB-valid exemplar per
                                       region), each vs its actual crystal ligand. Also writes
                                       mechanism_real_<region>.pdb (chain A = crystal/original,
                                       chain B = docked pose as docked, chain C = docked pose
                                       best-fit superposed on the crystal — shape-only),
                                       mechanism_real_examples.pml (grid view + 'as_docked' /
                                       'aligned' scenes to toggle B vs C),
                                       mechanism_real_<region>_with_receptor.pdb (receptor +
                                       crystal ligand resn CRY + docked pose resn DOK, in the
                                       binding pocket) with mechanism_real_with_receptor.pml, and
                                       mechanism_real_examples.csv (which pose was chosen per
                                       region). Needs --benchmark-dir (crystal SDFs) + the pose
                                       SDFs on disk; skips if unavailable.

    21_rank_concordance.png          — How faithfully does each tool's OWN ranking reproduce
                                       the oracle ordering of its poses (best pose → rank 1,
                                       2nd-best → rank 2, …)? Three panels: Crystal RMSD and
                                       Kabsch RMSD as grouped bars (% of complexes where the
                                       rank-k pose is the k-th best), PB-validity as achieved
                                       %-valid-by-rank vs the tool's overall-valid random baseline. AutoDock/
                                       DiffDock ranked by native confidence, EquiBind by gnina
                                       affinity. Writes rank_position_concordance.csv +
                                       rank_validity_by_rank.csv.
    21a/21b/21d_rank_vs_oracle_{crystal,kabsch,combined}.png — Individual line charts (more
                                       informative than 21's bars): median oracle rank (± IQR band)
                                       of the tool's rank-k pose vs k, with perfect-ranking diagonal
                                       and random baseline (≈ N/2). 21c_validity_by_rank.png = the
                                       PB-validity panel standalone. No stats box on any graph.
                                       Writes rank_vs_oracle.csv.
    22_rank_quality_whiskers.png     — Per-complex Kendall τ-b between each tool's ranking and
                                       the oracle order, as a grouped whisker plot: one group per
                                       criterion (Crystal RMSD, Kabsch RMSD, PB-validity, and a
                                       Combined valid + close + form criterion), one box per tool.
                                       τ = +1 perfect, 0 chance, −1 reversed. Writes rank_quality_tau.csv,
                                       rank_quality_tau_summary.csv, rank_quality_stats.csv.
    rank_quality_report.txt          — Full plain-text ranking-quality report: τ-b summary, the
                                       cross-tool paired tests (Friedman/Kendall-W + Wilcoxon-Holm
                                       pairwise — moved OFF the figures), rank-vs-oracle tables,
                                       identity-match, validity-by-rank, and the PB-valid complex mix.

    Oracle selection (fig 20 / 20b / 20d / 20e are each written TWICE):
      <name>.png              — RMSD-greedy oracle: the single nearest pose, kept only if it
                                is itself ≤ 2 Å AND PB-valid (matches oracle_pb_valid_and_rmsd2_%).
      <name>__valid_ceiling.png — validity-constrained oracle: the NEAREST pose that is ≤ 2 Å
                                AND PB-valid — the more generous sampling ceiling, which keeps a
                                complex whenever any valid near-native pose exists. Each figure's
                                title states which selection it uses.

PB-valid definition: a pose must pass EVERY canonical PoseBusters test
    (PB_CRITICAL_CHECKS = run_posebusters.CANONICAL_TEST_COLUMNS — the full 20-test
    intra- + intermolecular suite incl. internal energy and cofactor/water clashes).
    Identical to posebusters_validity_report.py / run_pandamap.py and to the
    PoseBusters paper, so every "valid" number in the study is directly comparable.

    pb_valid/                        — PB-valid variant of every RMSD figure (01–08),
                                       the twist/turn figure (14 + twist_turn_summary.csv)
                                       and the oracle/top1/rank CSVs, recomputed using
                                       only poses that pass all canonical PoseBusters
                                       tests. pb_valid_attrition.csv records, per
                                       method, how many pairs retain ≥1 valid pose
                                       (the denominator of these variant figures).

Output CSVs:
    per_pose_metrics.csv      — All metrics for every scored pose
    oracle_summary.csv        — Oracle stats per method (all 3 tools)
    top1_summary.csv          — Top-1 (rank-1) stats per method (ranking tools)
    per_rank_metrics.csv      — Per-rank success curves (ranking tools only)
    per_rank_ifp_recovery.csv — Per-rank interaction-fingerprint recovery vs crystal:
                                mean/median ProLIF IFP Tanimoto + native-contact
                                recovery at each rank k. rank_kind = native confidence
                                (AutoDock/DiffDock) or gnina affinity (EquiBind)
    oracle_rmsd_per_pair.csv  — Oracle RMSD pivot table (pairs × methods)
    comparison_vs_posebusters_paper.csv — This-study vs paper per method × metric
                                (this_study_% / posebusters_paper_% / delta_%)
    pbvalid_filter_influence.csv — Per method (oracle): rmsd_le_2A_%,
                                rmsd_le_2A_and_pb_valid_%, drop_pp, valid_retention_%
    twist_turn_summary.csv    — Per selection (top-1 / oracle) × method: median &
                                mean of centroid_dist, rot_angle_deg, bestfit_rmsd,
                                tfd, max/mean_torsion_dev_deg, n_torsions_flipped,
                                strain_energy (the "twisted & turned" decomposition)
    form_fidelity_summary.csv — Per method, over the near-native (≤ 2 Å) PB-valid
                                reps: n_success / success_%, best-fit (Kabsch) RMSD
                                median/mean/p90, tfd_median, form_correct_% (best-fit
                                ≤ --form-ok-kabsch), rms_inplace/rms_form/rms_placement
                                and form_share_of_error_% (>50 % ⇒ form-limited).
                                Also written as form_fidelity_summary__valid_ceiling.csv
                                for the validity-constrained oracle selection.
    oracle_selection_comparison.csv — Per method: current_nearest_rule (successes under
                                the RMSD-greedy oracle), relaxed_any_valid_le2A (under the
                                validity-constrained oracle) and rescuable (extra complexes
                                the nearest rule drops: closest pose ≤ 2 Å but invalid while
                                another ≤ 2 Å valid pose exists)
    oracle_rank_distribution.csv — Per ranking tool × rank k: pct_at_rank (% of
                                complexes whose oracle pose sits at rank k) + cum_pct
    pocket_localization.csv   — Per ranking tool: pct_topN_in_validated_pocket,
                                pct_oracle_diff_pocket_vs_topN (+ top1/any-topN/oracle
                                in-pocket context numbers, n_complexes, pocket_cutoff_A)
    pb_test_waterfall.csv     — Per method × cascade step: poses_removed / remaining /
                                remaining_pct for the PoseBusters failure waterfall
    topn_within_thresholds.csv — Per ranking tool × RMSD threshold (fine grid): how many
                                of the top-N ranked poses land within it — top1_within_%,
                                best_topN_within_% (any of top-N), pose_within_% (of all
                                pooled top-N poses), each with its raw count + n_pairs /
                                n_poses_topN denominators
    pose_validity_cascade.csv — Per docking variant (ALL of them): the pose-level cascade
                                produced → PoseBusters-valid → RMSD ≤ 2 Å → Kabsch RMSD < 1 Å
                                → (RMSD ≤ 2 Å & PB-valid & Kabsch < 1 Å), each stage as BOTH a
                                pose count and the number of distinct complexes still reached
                                (n_*_complexes / *_poses), plus the %-columns (pb_valid_%,
                                rmsd2_%, kabsch1_%, triple_of_rmsd2_%, triple_of_all_%). Kabsch
                                = best-fit / superposed RMSD (placement-free internal form;
                                threshold = --form-ok-kabsch, default 1 Å). The pure ≤ 2 Å &
                                PB-valid counts (paper headline, before the form filter) are
                                retained as the rmsd2_pbvalid_* columns. Pooled over every
                                produced pose — a sampling-volume view, NOT a per-complex
                                success rate. Companion human-readable table with grouped
                                headers + markdown: pose_validity_cascade_report.txt.

Per-pose cache (so re-runs that only change --top-n stay cheap)
    The heavy per-pose scoring is written to per_pose_metrics.csv alongside a
    per_pose_metrics.manifest.json fingerprint of the inputs it depends on (PB CSV
    signature, --ids-file, --split-equibind, --limit-pairs, --reference-convention —
    NOT --top-n, which only drives the cheap ranking aggregation). On the next run, if
    that fingerprint still matches, the cached metrics are reused and the per-pose
    scoring is skipped; pass --force to recompute regardless. Changing --top-n alone
    therefore reuses the cache. The manifest also records the reference convention
    at top level: --reuse-cache REFUSES a cache built under the other convention, and
    a plain run whose fingerprint differs ONLY in the convention refuses too (instead
    of silently recomputing) unless --force is given.

----------------------------------------------------------------------
Quick usage
----------------------------------------------------------------------
    python Scripts/Analysis/posebusters_pose_comparison.py

    python Scripts/Analysis/posebusters_pose_comparison.py --limit-pairs 20

    # Re-run a different ranking depth — reuses the cached per-pose metrics:
    python Scripts/Analysis/posebusters_pose_comparison.py --top-n 10

    # Force a full recompute of the per-pose metrics:
    python Scripts/Analysis/posebusters_pose_comparison.py --force

    # Score every pose against its nearest deposited ligand copy (PoseBusters
    # convention) instead of the single reference instance (the default):
    python Scripts/Analysis/posebusters_pose_comparison.py \\
        --reference-convention nearest --force --out-dir <new dir>

    python Scripts/Analysis/posebusters_pose_comparison.py \\
        --pb-csv  posebusters_results/benchmark/dock/posebusters_filtered_results.csv \\
        --benchmark-dir "Data/PoseBuster Benchmark Set" \\
        --out-dir posebusters_results/benchmark/dock/pose_comparison_report \\
        --top-n 5 --workers 8 --pocket-cutoff 6.0
"""

from __future__ import annotations
import argparse
import hashlib
import json
import math
import multiprocessing as mp
import os
import re
import subprocess
import sys
import tempfile
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path

import matplotlib
# Headless-safe backend: when NOT running under IPython/Jupyter and no backend
# was forced via MPLBACKEND, fall back to the non-interactive "Agg" backend.
# This script only saves figures (never shows them), so Agg is correct — and it
# avoids a Qt-backend segfault when the script is run as a CLI on a headless
# host. The notebook's inline backend (cell uses %run, in-process) is detected
# and left untouched so other cells still render normally.
if os.environ.get("MPLBACKEND") is None:
    try:
        from IPython import get_ipython
        _IN_IPYTHON = get_ipython() is not None
    except Exception:
        _IN_IPYTHON = False
    if not _IN_IPYTHON:
        matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
# Shared (A),(B),(C)… panel labeller for multi-panel figures.
from pocket_comparison_report import _label_panels  # noqa: E402
# Shared, unit-tested statistical helpers (Friedman/Kendall-W, paired-proportions
# with Cochran's Q + exact McNemar, Wilson CIs, Holm — see STATISTICAL_VALIDATION_PLAN.md).
import stats_utils as su  # noqa: E402
# Shared single-point method exclusion (--exclude-methods / --exclude-preset).
import method_filter as mf  # noqa: E402


def _json_default(o):
    """json.dump default: make numpy scalars/arrays serialisable in the stats sidecar."""
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, rdMolTransforms
from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate

try:
    from rdkit.Chem import TorsionFingerprints
    _HAS_TFD = True
except Exception:
    _HAS_TFD = False

RDLogger.DisableLog("rdApp.*")

try:
    import prolif as plf  # type: ignore
    _HAS_PROLIF = True
except Exception:
    _HAS_PROLIF = False

# Import posebusters' check_rmsd HERE — eagerly, in the parent, AFTER prolif.
# Order matters in this env: importing posebusters *before* prolif/MDAnalysis,
# or importing it lazily inside a forked worker (which already inherited
# MDAnalysis), both segfault. Importing it once in the parent after prolif is
# the only safe order; forked workers then inherit it and never re-import.
try:
    from posebusters.modules.rmsd import check_rmsd as _pb_check_rmsd
    _HAS_PB_RMSD = True
except Exception:
    _pb_check_rmsd = None
    _HAS_PB_RMSD = False

# Single source of truth for "PB-valid": reuse the pipeline's canonical PoseBusters
# test set so this script, posebusters_validity_report.py and run_pandamap.py agree
# EXACTLY on what counts as a valid pose — and so it matches the PoseBusters paper's
# validity definition (Buttenschoen, Morris & Deane, Chem. Sci. 2024, d3sc04185a:
# the full intra- + intermolecular suite incl. internal energy and cofactor/water
# clashes). run_posebusters lives one package over; add it to the path and import.
_PB_DIR = Path(__file__).resolve().parents[2] / "Scripts" / "Docking" / "Posebusters"
if str(_PB_DIR) not in sys.path:
    sys.path.insert(0, str(_PB_DIR))
try:
    from run_posebusters import CANONICAL_TEST_COLUMNS as _CANONICAL_PB_CHECKS  # noqa: E402
except Exception:
    _CANONICAL_PB_CHECKS = None

PLIF_INTERACTIONS = [
    "Hydrophobic", "HBDonor", "HBAcceptor", "PiStacking",
    "Anionic", "Cationic", "CationPi", "PiCation",
    "XBDonor", "XBAcceptor",
]

# ───────────────────────────────────────────────────────────────────
# Constants
# ───────────────────────────────────────────────────────────────────

_VDW = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98,
    "NA": 2.27, "MG": 1.73, "K": 2.75, "CA": 2.31, "MN": 2.05,
    "FE": 2.00, "CO": 2.00, "NI": 1.97, "CU": 1.96, "ZN": 1.39,
}
_DEFAULT_VDW = 1.70
CLASH_SCALE = 0.75
CONTACT_CUTOFF = 4.0
RMSD_THRESHOLDS = (1.0, 2.0, 5.0)
# Fine-grained RMSD grid (Å) for the "how close is the top-ranked pose" sweep (fig 18).
# The STANDARD (as-placed) RMSD runs 0 → 5 Å in 0.25 Å steps so the whole placement
# range is resolved; the Kabsch/best-fit view uses HALF that range (0 → 2.5 Å — see
# FINE_KABSCH_THRESHOLDS) because superposition removes placement error and collapses
# the distances toward zero, so a form-fidelity question is sub-Ångström. Overridable
# via --fine-rmsd-thresholds; the Kabsch grid always tracks it at half.
FINE_RMSD_THRESHOLDS = tuple(round(0.25 * i, 2) for i in range(21))    # 0, 0.25, … 5 Å
FINE_KABSCH_THRESHOLDS = tuple(round(t / 2.0, 4) for t in FINE_RMSD_THRESHOLDS)  # 0 … 2.5 Å
# Canonical docking-success line (Å): the in-place RMSD to the reference copy of
# the crystal ligand (see REFERENCE_CONVENTION) below which a pose is
# "near-native". Matches the 2 Å the paper and the rest of this report use.
NEAR_NATIVE_RMSD_A = 2.0

# ── Reference-ligand convention ───────────────────────────────────────────────
# Which deposited copy of the crystal ligand every crystal-referenced metric is
# measured against (module docstring, "Reference-ligand convention"):
#   instance — record 0 of <ID>_ligands.sdf (== <ID>_ligand.sdf), always. Default;
#              reproduces the frozen per-pose table exactly.
#   nearest  — per pose, the copy with the smallest symmetry-corrected in-place
#              RMSD (PoseBusters check_rmsd / paper convention); that one copy j*
#              drives every metric of the pose.
# ``REFERENCE_CONVENTION`` is set once from --reference-convention in main() and
# read by the text/label helpers below; process_pair receives the value through
# its work tuple, so forked workers never depend on this global.
REFERENCE_CONVENTIONS = ("instance", "nearest")
REFERENCE_CONVENTION = "instance"


def _set_reference_convention(conv: str) -> str:
    global REFERENCE_CONVENTION
    conv = str(conv or "instance")
    if conv not in REFERENCE_CONVENTIONS:
        raise ValueError(f"unknown reference convention {conv!r}; "
                         f"expected one of {REFERENCE_CONVENTIONS}")
    REFERENCE_CONVENTION = conv
    return conv


def _ref_noun() -> str:
    """Short noun for the reference ligand in labels: 'crystal' (instance) or
    'nearest deposited copy' (nearest). Under instance every label is unchanged."""
    return "crystal" if REFERENCE_CONVENTION == "instance" else "nearest deposited copy"


def _ref_ligand_phrase() -> str:
    """Longer noun phrase for prose: 'crystal ligand' or
    'nearest deposited copy of the crystal ligand'."""
    return ("crystal ligand" if REFERENCE_CONVENTION == "instance"
            else "nearest deposited copy of the crystal ligand")


def _ref_copy_label() -> str:
    """Label of the reference ligand in the mechanism-exemplar PDBs / legends."""
    return ("crystal (original)" if REFERENCE_CONVENTION == "instance"
            else "nearest deposited copy")


def _ref_convention_note() -> str:
    """One-line statement of the active convention for sidecar text files."""
    if REFERENCE_CONVENTION == "instance":
        return ("Reference convention: instance — every crystal-referenced metric is "
                "measured against the single reference instance (<ID>_ligand.sdf, "
                "record 0 of <ID>_ligands.sdf).")
    return ("Reference convention: nearest — every crystal-referenced metric is "
            "measured against the deposited copy of <ID>_ligands.sdf nearest to the "
            "pose (argmin symmetry-corrected in-place RMSD; PoseBusters check_rmsd "
            "convention), one copy per pose.")

# ── CSV numeric precision ────────────────────────────────────────────────────
# Percentage / percentage-point columns are written to the machine-readable CSVs
# with this many decimal places. 1–2 dp collapses near-ties that matter when
# ranking variants — e.g. AutoDock raw vs +gnina both print "1.7 % all" yet are
# 1.669800 % vs 1.680858 %. The fixed-width and markdown console tables format
# their own (usually 1 dp) precision independently, so this only widens the CSVs.
PCT_DECIMALS = 6


def _is_pct_name(name) -> bool:
    """True for a percentage / percentage-point column, recognised by the naming
    convention used throughout this script: the name ends in ``%`` or ``_pp`` (or
    ``_pp_...``), or contains ``pct`` / ``percent``."""
    n = str(name).lower()
    return (n.endswith("%") or n.endswith("_pp") or "_pp_" in n
            or "pct" in n or "percent" in n)


def _round_csv(df, pct=PCT_DECIMALS, other=2):
    """Return a copy of ``df`` rounded for CSV output: percentage / percentage-point
    columns (see :func:`_is_pct_name`) to ``pct`` decimals so near-ties stay
    distinguishable, and every other floating-point column to ``other`` (pass
    ``other=None`` to leave non-percentage columns at full precision). Integer and
    non-numeric columns pass through untouched."""
    if df is None or not hasattr(df, "columns"):
        return df
    out = df.copy()
    for c in out.columns:
        if (not pd.api.types.is_numeric_dtype(out[c])
                or pd.api.types.is_integer_dtype(out[c])
                or pd.api.types.is_bool_dtype(out[c])):
            continue
        dp = pct if _is_pct_name(c) else other
        if dp is not None:
            out[c] = out[c].round(dp)
    return out


def _write_csv(df, path, pct=PCT_DECIMALS, **kwargs):
    """``df.to_csv`` drop-in that writes percentage / percentage-point columns (see
    :func:`_is_pct_name`) with a FIXED ``pct`` decimal places, so the text always
    carries at least that many digits — plain float rounding drops trailing zeros
    (1.669800 → "1.6698"), which can dip below the promised precision on "round"
    values. Only percentage columns are reformatted (as strings, NaN → empty); every
    other column — counts, RMSDs, p-values — is written exactly as it stands, so tiny
    p-values keep their full precision. The frame is copied only when it actually has
    a percentage column, so large non-percentage tables are written without overhead."""
    cols = getattr(df, "columns", [])
    pct_cols = [c for c in cols
                if _is_pct_name(c) and pd.api.types.is_numeric_dtype(df[c])
                and not pd.api.types.is_bool_dtype(df[c])]
    if pct_cols:
        df = df.copy()
        fmt = "%." + str(int(pct)) + "f"
        for c in pct_cols:
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else fmt % v)
    df.to_csv(path, **kwargs)
# A near-native, PB-valid pose has the "correct form" when its best-fit (Kabsch)
# RMSD to the reference copy of the crystal ligand (the instance, or the pose's
# nearest copy j* under --reference-convention nearest) — heavy-atom RMSD AFTER
# optimal superposition, so translation and rotation are removed and only the
# internal conformation is compared — is within this many Å. Overridable via
# --form-ok-kabsch.
FORM_OK_KABSCH_A = 1.0
CENTROID_THRESHOLD = 4.0
# A pose counts as being in the experimentally validated pocket if its centroid is
# within this distance (Å) of the reference-copy centroid (``centroid_dist``, which
# under --reference-convention nearest is measured at the pose's nearest copy j*,
# so an alternate deposited site counts as validated for the poses nearest to it).
# Looser than CENTROID_THRESHOLD (which marks a near-correct PLACEMENT): a pocket
# spans ~10 Å, so a centroid within ~6 Å of the crystal sits in the same site, while
# anything farther is treated as a different (decoy) pocket. Overridable via
# --pocket-cutoff.
POCKET_CENTROID_CUTOFF = 6.0

# AutoDock-specific optimizer identities.  DiffDock and EquiBind deliberately
# retain their existing raw/smina/gnina semantics; the CNN-refinement alias is an
# AutoDock-only protocol and must remain separate from legacy ``gnina``.
_AUTODOCK_OPTIMIZERS = frozenset({"smina", "gnina", "gnina_refinement"})
_AUTODOCK_OPTIMIZER_SUFFIXES = ("_gnina_refinement", "_smina", "_gnina")

# Methods that produce an explicit pose ranking ordered by energy/confidence.
# AutoDock optimizer variants use their own ``optimized_rank``; the raw method
# uses Vina's ``autodock_rank``. Keeping every protocol here makes each aggregation
# group and rank the variants independently. Other tools (e.g. EquiBind) remain
# oracle-only unless a dedicated ranking adapter handles them.
RANKING_TOOLS = frozenset({
    "autodock", "autodock_smina", "autodock_gnina",
    "autodock_gnina_refinement",
    # Vinardo-scored AutoDock is the same Vina engine, different scoring function.
    "autodock_vinardo", "autodock_vinardo_smina", "autodock_vinardo_gnina",
    "autodock_vinardo_gnina_refinement",
    # MGLTools-ligand arm: same Vina engine + ranks (autodock_rank / optimized_rank).
    "autodock_mgltools", "autodock_mgltools_smina", "autodock_mgltools_gnina",
    "autodock_mgltools_gnina_refinement",
    # Exhaustiveness sweep on that same arm (the unsuffixed key above is its 32 point).
    "autodock_mgltools_exh18",
    "autodock_mgltools_exh64", "autodock_mgltools_exh64_gnina",
    "autodock_mgltools_exh92",
    "autodock_mgltools_exh128", "autodock_mgltools_exh128_gnina",
    # Uni-Dock (tiled) and Uni-Dock2 emit affinity-ordered poses (rank = pose order).
    "unidock", "unidock2",
    "diffdock",
})

DEFAULT_TOP_N = 5   # ranked poses per pair to consider in Part B

# PB-valid = a pose passes EVERY canonical PoseBusters test (the full intra- and
# intermolecular suite). Imported from run_posebusters so it stays identical to the
# rest of the pipeline (validity report, PandaMap) and to the paper; the literal
# fallback below is only used if that import fails. Only columns present in a given
# CSV are applied downstream, so 'mol'-mode CSVs (no protein/cofactor columns) still
# work. NOTE: this is the canonical 20-test set — it adds internal_energy and the
# cofactor/water clash & overlap checks that the older 12-test list omitted, so
# "PB-valid" here is now stricter and matches the paper's coral-bar definition.
PB_CRITICAL_CHECKS = list(_CANONICAL_PB_CHECKS) if _CANONICAL_PB_CHECKS else [
    "mol_pred_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
]

# Published reference numbers from the PoseBusters paper, for cross-checking ours.
# Source: Buttenschoen, Morris & Deane, "PoseBusters: AI-based docking methods fail
# to generate physically valid poses or generalise to novel sequences", Chem. Sci.
# 2024, 15, 3130-3139, DOI 10.1039/d3sc04185a — **PoseBusters Benchmark set (308
# complexes)**, single top-ranked predicted pose per target. From the results text:
#   "Gold (55%) and AutoDock Vina (58%) perform the best ... with ... and without
#    ... considering PB-validity"; "the best performing DL method, DiffDock (12%),
#    does not compete with the two standard docking methods"; "EquiBind, Uni-Mol
#    and TankBind generate almost no physically valid poses that pass all tests."
#   rmsd2       — % targets with RMSD ≤ 2 Å.  None = not stated unambiguously in the
#                 text; read it off Fig. 1 / the ESI table and fill in to plot it.
#   rmsd2_valid — % with RMSD ≤ 2 Å AND PB-valid (the paper's headline metric).
# This is a DIFFERENT experiment (their receptor prep / box / EquiBind variant), so
# treat the numbers as a literature sanity-check, not ground truth. Edit freely.
POSEBUSTERS_PAPER_BENCHMARK: dict[str, dict] = {
    "AutoDock Vina": {"rmsd2": 58.0, "rmsd2_valid": 58.0, "note": "text"},
    "DiffDock":      {"rmsd2": None, "rmsd2_valid": 12.0, "note": "valid 12% (text); RMSD≤2Å: see Fig.1"},
    "EquiBind":      {"rmsd2": None, "rmsd2_valid": 2.0,  "note": "paper: 'almost no valid poses' (<2%)"},
}


def _paper_method_name(method_key: str) -> str | None:
    """Map an internal method key to the matching PoseBusters-paper method name."""
    if method_key == "autodock":
        return "AutoDock Vina"
    if method_key == "diffdock":
        return "DiffDock"
    if method_key.startswith("equibind"):
        return "EquiBind"
    return None


# ── Method labelling / colours (EquiBind variant-aware) ─────────────
# EquiBind tags every pose "equibind_guided" in the raw CSV; the real flavour
# is encoded in the pose filename and unpacked into per-variant method labels
# (see _classify_equibind / _apply_equibind_split). TOOL_LABEL / TOOL_COLORS are
# therefore smart lookups (not plain dicts) so every existing ``.get(m, …)`` call
# site renders any variant — equibind_guided_smina, equibind_unguided, … — with a
# readable label and a stable colour, with no per-plot changes.

_BASE_COLORS = {
    "autodock": "#1f77b4",           # AutoDock Vina (blue)
    "autodock_vinardo": "#17becf",   # AutoDock Vinardo (cyan — AutoDock-family cousin)
    "unidock": "#9467bd",            # Uni-Dock tiled (purple)
    "unidock2": "#8c564b",           # Uni-Dock2 (brown)
    "diffdock": "#ff7f0e",           # DiffDock (orange)
}
_AD_VARIANT_COLORS = {
    "autodock_smina": "#6baed6", "autodock_gnina": "#08519c",
    "autodock_gnina_refinement": "#08306b",
    "autodock_vinardo_smina": "#9edae5", "autodock_vinardo_gnina": "#0e7c86",
    "autodock_vinardo_gnina_refinement": "#005f69",
    # MGLTools-ligand arm — OLIVE, ramped light→dark by search effort. Olive because
    # every other family hue in this file is taken and 09f draws them in one axes:
    # blue = AutoDock Vina, cyan = Vinardo, orange = DiffDock, purple = Uni-Dock,
    # brown = Uni-Dock2, and GREEN is EquiBind's family (_EQ_PALETTE below, which
    # _MethodColorMap cycles for any equibind* key). The previous greens here collided
    # with it badly — "#00441b" was byte-identical to equibind_unguided_raw and
    # "#78c679" was 5/255 from equibind_unguided_smina, both drawn in the same figure.
    # Within the ramp, lightness encodes exhaustiveness (house rule: tiers ramp, they
    # do not change hue); exh64's old orange is what made it read as its own engine.
    "autodock_mgltools_exh18": "#efeeb8",
    "autodock_mgltools": "#dcda6e",          # the exhaustiveness-32 point of the sweep
    "autodock_mgltools_exh64": "#c6c72e",
    "autodock_mgltools_exh92": "#a5a61f",
    "autodock_mgltools_exh128": "#8f9019",
    # Post-dock optimiser variants — darker tail of the same hue, so they read as
    # derived from their exhaustiveness point rather than as further search effort.
    "autodock_mgltools_smina": "#6b6c13",
    "autodock_mgltools_gnina": "#4e4f0e",
    "autodock_mgltools_gnina_refinement": "#33340a",
    "autodock_mgltools_exh64_gnina": "#5c5d10",
    "autodock_mgltools_exh128_gnina": "#787917",
}
# Green family for EquiBind variants (cycled if more than this many appear).
_EQ_PALETTE = ["#2ca02c", "#74c476", "#1b7837", "#a6dba0",
               "#006d2c", "#5aae61", "#00441b", "#c7e9c0"]
# Orange family for DiffDock smina/gnina optimizer variants (raw uses the base
# orange in _BASE_COLORS).
_DD_VARIANT_COLORS = {"diffdock_smina": "#d95f02", "diffdock_gnina": "#fdae6b"}

# Per-run display-label overrides (method key -> label). Populated in main() when
# --best-equibind-only is active, where the single retained EquiBind variant is
# shown as "EquiBind*". Honoured by _pretty_method, so every TOOL_LABEL.get(...)
# call site in the plots picks it up with no further changes.
_LABEL_OVERRIDES: dict[str, str] = {}

# Appended to every figure title. main() sets this while regenerating the
# PB-valid-only variant of each plot (poses passing all canonical PoseBusters
# tests), then resets it to "". Read by _vt() at each ax.set_title/suptitle call
# site so the variant figures are self-labelling with no per-plot signature change.
_PLOT_TITLE_SUFFIX = ""


def _vt(title: str) -> str:
    """Variant-aware title: append the active ``_PLOT_TITLE_SUFFIX`` (if any)."""
    return f"{title}{_PLOT_TITLE_SUFFIX}"


def _n_annot(summ: "pd.DataFrame | None", method: str, *, poses: bool = False) -> str:
    """Per-method count annotation for axis/tick labels.

    Returns ``'n=<pairs>'`` — the number of receptor-ligand pairs represented
    (the denominator) — and, when ``poses`` is set, also ``'· ≤<max> poses'``,
    the per-method pose cap the oracle selects from (it varies by method:
    Vina/DiffDock generate more poses per complex than EquiBind, and some
    complexes have fewer). Returns ``''`` if the method/columns are absent.
    """
    if summ is None or method not in getattr(summ, "index", []):
        return ""
    parts: list[str] = []
    if "n_pairs" in summ.columns and not pd.isna(summ.loc[method, "n_pairs"]):
        parts.append(f"n={int(summ.loc[method, 'n_pairs'])}")
    if (poses and "poses_per_pair_max" in summ.columns
            and not pd.isna(summ.loc[method, "poses_per_pair_max"])):
        parts.append(f"≤{int(summ.loc[method, 'poses_per_pair_max'])} poses")
    return " · ".join(parts)


def _eq_tokens(method: str) -> tuple[str | None, str | None, str | None]:
    """Split an ``equibind_*`` label into (pocket, refine, clamp) tokens."""
    pocket = refine = clamp = None
    for t in method.split("_")[1:]:          # drop the leading "equibind"
        if t in ("unguided", "fpocket", "p2rank", "guided"):
            pocket = t
        elif t in ("raw", "smina", "gnina"):
            refine = t
        elif t in ("clampON", "clampOFF"):
            clamp = t
    return pocket, refine, clamp


def _pretty_method(m: str) -> str:
    """Human-readable label for a method / EquiBind-variant key."""
    if m in _LABEL_OVERRIDES:
        return _LABEL_OVERRIDES[m]
    if m == "autodock":
        return "AutoDock Vina"
    if m == "autodock_smina":
        return "AutoDock Vina (smina-opt)"
    if m == "autodock_gnina":
        return "AutoDock Vina (gnina-opt)"
    if m == "autodock_gnina_refinement":
        return "AutoDock Vina (GNINA CNN-refinement)"
    if m == "autodock_vinardo":
        return "AutoDock Vinardo"
    if m == "autodock_vinardo_smina":
        return "AutoDock Vinardo (smina-opt)"
    if m == "autodock_vinardo_gnina":
        return "AutoDock Vinardo (gnina-opt)"
    if m == "autodock_vinardo_gnina_refinement":
        return "AutoDock Vinardo (GNINA CNN-refinement)"
    # The unsuffixed MGLTools key is the exhaustiveness-32 arm; name the setting so it
    # cannot be read as an exhaustiveness-agnostic baseline.
    if m == "autodock_mgltools":
        return "AutoDock Vina MGL-lig exh32"
    if m == "autodock_mgltools_smina":
        return "AutoDock Vina MGL-lig exh32 (smina-opt)"
    if m == "autodock_mgltools_gnina":
        return "AutoDock Vina MGL-lig exh32 (gnina-opt)"
    if m == "autodock_mgltools_gnina_refinement":
        return "AutoDock Vina MGL-lig exh32 (GNINA CNN-refinement)"
    if m == "autodock_mgltools_exh18":
        return "AutoDock Vina MGL-lig exh18"
    if m == "autodock_mgltools_exh64":
        return "AutoDock Vina MGL-lig exh64"
    if m == "autodock_mgltools_exh64_gnina":
        return "AutoDock Vina MGL-lig exh64 (gnina-opt)"
    if m == "autodock_mgltools_exh128_gnina":
        return "AutoDock Vina MGL-lig exh128 (gnina-opt)"
    if m == "autodock_mgltools_exh92":
        return "AutoDock Vina MGL-lig exh92"
    if m == "autodock_mgltools_exh128":
        return "AutoDock Vina MGL-lig exh128"
    if m == "unidock":
        return "Uni-Dock"
    if m == "unidock2":
        return "Uni-Dock2"
    if m == "diffdock":
        return "DiffDock"
    if m == "diffdock_smina":
        return "DiffDock (smina-opt)"
    if m == "diffdock_gnina":
        return "DiffDock (gnina-opt)"
    if not m.startswith("equibind"):
        return m
    pocket, refine, clamp = _eq_tokens(m)
    parts: list[str] = []
    if pocket:
        parts.append(pocket)
    if refine:
        parts.append({"smina": "smina-opt", "gnina": "gnina-opt"}.get(refine, "raw"))
    if clamp:
        parts.append("clamp on" if clamp == "clampON" else "clamp off")
    return f"EquiBind ({', '.join(parts)})" if parts else "EquiBind"


class _MethodLabelMap:
    """dict-like: ``.get(method, default)`` / ``[method]`` -> pretty label."""

    def get(self, key, default=None):
        return _pretty_method(str(key))

    def __getitem__(self, key):
        return _pretty_method(str(key))


class _MethodColorMap:
    """dict-like: stable colour per method. AutoDock/DiffDock are fixed; EquiBind
    variants draw from a green palette (common variants pre-seeded for stability).
    """

    _PRESEED = ["equibind_unguided", "equibind_fpocket", "equibind_p2rank",
                "equibind_unguided_gnina", "equibind_fpocket_gnina", "equibind_p2rank_gnina",
                "equibind_unguided_raw", "equibind_fpocket_raw", "equibind_p2rank_raw",
                "equibind_unguided_smina", "equibind_fpocket_smina", "equibind_p2rank_smina"]

    def __init__(self):
        self._cache = {m: _EQ_PALETTE[i % len(_EQ_PALETTE)]
                       for i, m in enumerate(self._PRESEED)}
        self._next = len(self._PRESEED)

    def get(self, key, default=None):
        key = str(key)
        if key in _BASE_COLORS:
            return _BASE_COLORS[key]
        if key in _AD_VARIANT_COLORS:
            return _AD_VARIANT_COLORS[key]
        if key in _DD_VARIANT_COLORS:
            return _DD_VARIANT_COLORS[key]
        if key.startswith("equibind"):
            if key not in self._cache:
                self._cache[key] = _EQ_PALETTE[self._next % len(_EQ_PALETTE)]
                self._next += 1
            return self._cache[key]
        return default

    def __getitem__(self, key):
        return self.get(key)


TOOL_LABEL = _MethodLabelMap()
TOOL_COLORS = _MethodColorMap()

# ───────────────────────────────────────────────────────────────────
# Helpers
# ───────────────────────────────────────────────────────────────────


def _to_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(bool)
    return (
        s.astype(str).str.strip().str.lower()
        .map({"true": True, "1": True, "1.0": True,
              "false": False, "0": False, "0.0": False, "nan": False, "": False})
        .fillna(False).astype(bool)
    )


def parse_rank(method: str, pose_file: str) -> int:
    """Extract pose rank from output filename (1 = top-ranked by the tool).

    Returns 999 for tools that provide no canonical ranking so those poses
    are never mistakenly selected as a valid top-k ranked pose.
    """
    name = Path(pose_file).name
    # AutoDock (Vina + Vinardo) and Uni-Dock / Uni-Dock2 all expand to per-pose
    # SDFs named "<stem>_model<N>.sdf" (N = the tool's affinity rank, 1 = best);
    # see run_posebusters._expand_{autodock,unidock,unidock2}_poses.
    if method.startswith("autodock") or method.startswith("unidock"):
        m = re.search(r"_model(\d+)\.sdf$", name)
        return int(m.group(1)) if m else 999
    # Every DiffDock flavour (base + smina/gnina-optimized variants, e.g.
    # diffdock_smina / diffdock_gnina) keeps the "rankNN" token in its filename,
    # so match on the prefix — otherwise the optimizer variants fall through to
    # 999 and, once --best-diffdock-only relabels the winner to "diffdock", the
    # Part-B ranking metrics see zero rank ≤ top_n poses.
    if method.startswith("diffdock"):
        m = re.search(r"rank(\d+)", name)
        return int(m.group(1)) if m else 999
    return 999  # unranked tool


def _positive_rank(value) -> int | None:
    """Return a positive integral rank, else ``None`` (NaN-safe)."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(number) or number < 1 or not number.is_integer():
        return None
    return int(number)


def _finite_number(value) -> float | None:
    """Return a finite float for CSV provenance, otherwise ``None``."""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _pose_rank(row, method: str, pose_file: str) -> int:
    """Resolve the ranking axis for a pose without crossing variant boundaries.

    Raw AutoDock poses use the immutable Vina ``autodock_rank``. Optimized
    AutoDock variants *only* use ``optimized_rank``; falling back to the source
    rank would silently undo gnina/smina re-ranking. Legacy raw rows may still
    use the converted ``_modelN.sdf`` filename fallback.
    """
    optimized_autodock = method.startswith("autodock") and method.endswith(
        _AUTODOCK_OPTIMIZER_SUFFIXES)
    if method.startswith("autodock") and not optimized_autodock:
        return _positive_rank(row.get("autodock_rank")) or parse_rank(method, pose_file)
    if optimized_autodock:
        return _positive_rank(row.get("optimized_rank")) or 999
    if method == "unidock":
        return _positive_rank(row.get("unidock_rank")) or parse_rank(method, pose_file)
    if method == "unidock2":
        return _positive_rank(row.get("unidock2_rank")) or parse_rank(method, pose_file)
    return parse_rank(method, pose_file)


def load_protein_heavy_atoms(pdb_path: Path) -> tuple[np.ndarray, list[str], list[tuple[str, int]]]:
    """Lightweight PDB parser → (coords [N×3], elements, residue_id list).

    Skips waters and hydrogens.  Residue id is (chain, resseq).
    """
    coords, elems, resids = [], [], []
    with open(pdb_path) as fh:
        for ln in fh:
            if not (ln.startswith("ATOM  ") or ln.startswith("HETATM")):
                continue
            resname = ln[17:20].strip()
            if resname in {"HOH", "WAT", "DOD"}:
                continue
            elem = (ln[76:78].strip() or ln[12:16].strip()[0]).upper()
            if elem == "H":
                continue
            try:
                x = float(ln[30:38]); y = float(ln[38:46]); z = float(ln[46:54])
            except ValueError:
                continue
            chain = ln[21].strip() or "A"
            try:
                resseq = int(ln[22:26])
            except ValueError:
                continue
            coords.append((x, y, z))
            elems.append(elem)
            resids.append((chain, resseq))
    return np.asarray(coords, dtype=np.float32), elems, resids


def load_first_mol(sdf_path: Path) -> Chem.Mol | None:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    for m in suppl:
        if m is not None:
            try:
                Chem.SanitizeMol(m)
            except Exception:
                try:
                    m.UpdatePropertyCache(strict=False)
                    Chem.GetSymmSSSR(m)
                except Exception:
                    return None
            return m
    return None


def load_all_mols(sdf_path: Path) -> list[Chem.Mol]:
    """Every record of an SDF (e.g. all deposited copies in <ID>_ligands.sdf), with
    the same per-record sanitize fallback as :func:`load_first_mol`. Records that
    fail both sanitization and the ring-perception fallback are skipped, so the
    returned index is NOT guaranteed to equal the file record index in that case."""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mols: list[Chem.Mol] = []
    for k, m in enumerate(suppl):
        if m is None:
            print(f"  WARNING: {sdf_path.name} record {k} could not be read; skipped "
                  "(copy indices below this record shift by one).")
            continue
        try:
            Chem.SanitizeMol(m)
        except Exception:
            try:
                m.UpdatePropertyCache(strict=False)
                Chem.GetSymmSSSR(m)
            except Exception:
                print(f"  WARNING: {sdf_path.name} record {k} failed sanitization and "
                      "ring perception; skipped (copy indices below this record shift by one).")
                continue
        mols.append(m)
    return mols


def reassign_template(pose: Chem.Mol, template: Chem.Mol) -> Chem.Mol | None:
    """Reassign bond orders from crystal template (handles PDBQT round-trips)."""
    try:
        return AssignBondOrdersFromTemplate(template, pose)
    except Exception:
        return None


def symmetry_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float:
    """Heavy-atom RMSD over all symmetry-equivalent matchings, NO superposition.

    The protein frame is shared across all docking methods so atom positions
    are directly comparable.  Enumerates substructure matches of ref in pose
    (captures topological symmetry) and returns the minimum coordinate RMSD.
    """
    pose = Chem.RemoveHs(pose)
    ref = Chem.RemoveHs(ref)
    # Bulk C-level coordinate accessor: ~70x faster than a per-atom
    # ``list(conf.GetAtomPosition(i))`` comprehension (the per-pose hotspot).
    pose_xyz = pose.GetConformer().GetPositions()
    ref_coords = ref.GetConformer().GetPositions()

    matches = pose.GetSubstructMatches(ref, uniquify=False, useChirality=False)
    if not matches:
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return float("nan")
        matches = [tuple(range(ref.GetNumAtoms()))]

    best = math.inf
    for mp_ in matches:
        pose_coords = pose_xyz[list(mp_)]
        diff = pose_coords - ref_coords
        rmsd = float(np.sqrt((diff * diff).sum() / len(ref_coords)))
        if rmsd < best:
            best = rmsd
    return best if best != math.inf else float("nan")


def _safe_symmetry_rmsd(pose: Chem.Mol, ref: Chem.Mol) -> float:
    """symmetry_rmsd that returns NaN instead of raising (per-copy scoring)."""
    try:
        return symmetry_rmsd(pose, ref)
    except Exception:
        return float("nan")


def centroid_distance(pose: Chem.Mol, ref: Chem.Mol) -> float:
    pose = Chem.RemoveHs(pose); ref = Chem.RemoveHs(ref)
    a = pose.GetConformer().GetPositions().mean(axis=0)
    b = ref.GetConformer().GetPositions().mean(axis=0)
    return float(np.linalg.norm(a - b))


def posebusters_rmsd(pose: Chem.Mol, ref: Chem.Mol,
                     threshold: float = 2.0) -> tuple[float, float, bool]:
    """PoseBusters' canonical ``check_rmsd`` of *pose* vs crystal *ref*.

    Returns (rmsd, kabsch_rmsd, within_threshold):
      * rmsd         — symmetry-corrected heavy-atom RMSD, NO superposition
                       (in-frame docking accuracy; comparable to symmetry_rmsd)
      * kabsch_rmsd  — same after optimal superposition (conformer similarity only)
      * within       — rmsd <= threshold (default 2.0 Å)
    NaN/False if PoseBusters is unavailable or the comparison fails (e.g. atom
    mismatch from a bad PDBQT round-trip).
    """
    if not _HAS_PB_RMSD:
        return float("nan"), float("nan"), False
    try:
        res = _pb_check_rmsd(Chem.RemoveHs(pose), Chem.RemoveHs(ref),
                             rmsd_threshold=threshold, heavy_only=True)["results"]
        return (float(res["rmsd"]), float(res.get("kabsch_rmsd", float("nan"))),
                bool(res["rmsd_within_threshold"]))
    except Exception:
        return float("nan"), float("nan"), False


def uff_strain(pose: Chem.Mol) -> float:
    """Strain = E(pose) − E(UFF-minimized pose).  Returns NaN on failure."""
    try:
        m = Chem.AddHs(pose, addCoords=True)
        ff = AllChem.UFFGetMoleculeForceField(m)
        if ff is None:
            return float("nan")
        e_now = ff.CalcEnergy()
        m_min = Chem.Mol(m)
        ff_min = AllChem.UFFGetMoleculeForceField(m_min)
        ff_min.Minimize(maxIts=200)
        e_min = ff_min.CalcEnergy()
        return float(e_now - e_min)
    except Exception:
        return float("nan")


def rigid_body_fit(pose: Chem.Mol, ref: Chem.Mol) -> tuple[float, float]:
    """Decompose pose-vs-crystal placement into a "turn" and a "twist".

    Returns ``(rotation_angle_deg, bestfit_rmsd)``:
      * rotation_angle_deg — rigid-body rotation of the fit ("turned"), taken
        from the matching with the lowest AS-PLACED RMSD (same criterion as
        ``symmetry_rmsd``); selecting by post-fit RMSD would let a symmetry-
        flipped relabelling of a near-perfect pose look like a ~180° turn.
      * bestfit_rmsd — lowest POST-superposition heavy-atom RMSD over all
        symmetry matchings ("twisted"): the irreducible conformational
        difference once translation+rotation are removed. This is the FALLBACK
        used only when PoseBusters' ``pb_kabsch_rmsd`` is unavailable; the caller
        prefers the PoseBusters value (see process_pair).
    Both NaN on atom mismatch.
    """
    pose = Chem.RemoveHs(pose); ref = Chem.RemoveHs(ref)
    # Bulk coordinate accessors (see symmetry_rmsd): avoid the slow per-atom
    # ``list(GetAtomPosition(i))`` comprehension.
    pose_xyz = pose.GetConformer().GetPositions()
    ref_coords = ref.GetConformer().GetPositions()
    matches = pose.GetSubstructMatches(ref, uniquify=False, useChirality=False)
    if not matches:
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return float("nan"), float("nan")
        matches = [tuple(range(ref.GetNumAtoms()))]

    qc = ref_coords - ref_coords.mean(axis=0)
    best_placed, best_angle = math.inf, float("nan")
    best_fit = math.inf
    for mp_ in matches:
        P = pose_xyz[list(mp_)]
        pcen = P - P.mean(axis=0)
        U, _, Vt = np.linalg.svd(pcen.T @ qc)       # Kabsch fit of pcen onto qc
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
        fit_rmsd = float(np.sqrt(((pcen @ R.T - qc) ** 2).sum() / len(qc)))
        if fit_rmsd < best_fit:
            best_fit = fit_rmsd
        placed = float(np.sqrt(((P - ref_coords) ** 2).sum() / len(ref_coords)))
        if placed < best_placed:
            best_placed = placed
            cos = (np.trace(R) - 1.0) / 2.0
            best_angle = float(math.degrees(math.acos(max(-1.0, min(1.0, cos)))))
    return best_angle, (best_fit if best_fit != math.inf else float("nan"))


def torsion_metrics(pose: Chem.Mol,
                    ref: Chem.Mol) -> tuple[float, float, float, int, int]:
    """How much the pose is internally *twisted* vs the crystal conformation.

    Maps the pose coordinates onto the crystal topology (one molecule, two
    conformers) so rotatable-bond dihedrals and the Torsion Fingerprint
    Deviation are directly comparable. Returns
    ``(tfd, max_torsion_dev_deg, mean_torsion_dev_deg, n_torsions_flipped,
    n_rot_bonds)``: ``tfd`` ∈ [0, 1] (RDKit, NaN if unavailable); torsion devs
    are |Δdihedral| folded to [0, 180]; a bond is 'flipped' when it moved > 90°.
    """
    nan = float("nan")
    try:
        pose = Chem.RemoveHs(pose); ref = Chem.RemoveHs(ref)
        match = pose.GetSubstructMatch(ref)          # match[i_ref] = pose atom idx
        if not match or len(match) != ref.GetNumAtoms():
            return nan, nan, nan, 0, 0
        work = Chem.Mol(ref)                          # crystal topology + conf 0
        pose_conf = pose.GetConformer()
        conf = Chem.Conformer(work.GetNumAtoms())
        for i in range(work.GetNumAtoms()):
            conf.SetAtomPosition(i, pose_conf.GetAtomPosition(match[i]))
        cid_pose = work.AddConformer(conf, assignId=True)
        cid_ref = work.GetConformer(0).GetId()

        patt = Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]")
        devs: list[float] = []
        for a, b in work.GetSubstructMatches(patt):
            na = [x.GetIdx() for x in work.GetAtomWithIdx(a).GetNeighbors()
                  if x.GetIdx() != b]
            nb = [x.GetIdx() for x in work.GetAtomWithIdx(b).GetNeighbors()
                  if x.GetIdx() != a]
            if not na or not nb:
                continue
            i, l = na[0], nb[0]
            dref = rdMolTransforms.GetDihedralDeg(work.GetConformer(cid_ref), i, a, b, l)
            dpose = rdMolTransforms.GetDihedralDeg(work.GetConformer(cid_pose), i, a, b, l)
            devs.append(abs((dpose - dref + 180.0) % 360.0 - 180.0))
        n_rot = len(devs)
        max_dev = max(devs) if devs else 0.0
        mean_dev = (sum(devs) / n_rot) if devs else 0.0
        n_flip = sum(1 for d in devs if d > 90.0)

        tfd = nan
        if _HAS_TFD:
            try:
                tfd = float(TorsionFingerprints.GetTFDBetweenConformers(
                    work, [cid_ref], [cid_pose])[0])
            except Exception:
                tfd = nan
        return tfd, max_dev, mean_dev, n_flip, n_rot
    except Exception:
        return nan, nan, nan, 0, 0


def _load_prolif_protein(pdb_path: Path) -> "plf.Molecule | None":
    if not _HAS_PROLIF:
        return None
    try:
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tf:
            tmp = tf.name
        subprocess.run(["obabel", str(pdb_path), "-O", tmp, "-h"],
                       check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL, timeout=120)
        rdmol = Chem.MolFromPDBFile(tmp, removeHs=False, sanitize=False,
                                    proximityBonding=True)
        try:
            os.unlink(tmp)
        except OSError:
            pass
        if rdmol is None:
            return None
        return plf.Molecule(rdmol)
    except Exception:
        return None


def _ligand_for_plif(mol: Chem.Mol) -> "plf.Molecule | None":
    if not _HAS_PROLIF or mol is None:
        return None
    try:
        m = Chem.AddHs(mol, addCoords=True)
        return plf.Molecule.from_rdkit(m)
    except Exception:
        return None


def compute_plif_recovery(crystal_mol,
                          pose_mols: list[Chem.Mol],
                          prot_mol: "plf.Molecule | None"):
    """Tanimoto similarity of pose IFP vs crystal IFP (single ProLIF run).

    ``crystal_mol`` is either ONE reference mol (legacy call: returns a flat list
    of length n_poses) or a LIST of deposited copies (returns an n_poses × n_copies
    nested list ``out[i][j]``; the caller selects column j* per pose). All copies
    are passed as the leading ligands of the same ProLIF run as the poses, so the
    pose fingerprints are computed once. A copy that ``_ligand_for_plif`` cannot
    convert — or a ``None`` placeholder, used by the instance convention to skip
    the copies it never scores — gets a NaN column instead of NaN-ing the whole
    pair; the run happens only if at least one copy converts.
    """
    single = not isinstance(crystal_mol, (list, tuple))
    ref_mols = [crystal_mol] if single else list(crystal_mol)
    n = len(pose_mols)
    n_ref = len(ref_mols)
    nan = float("nan")

    def _nan_out():
        if single:
            return [nan] * n
        return [[nan] * n_ref for _ in range(n)]

    if not _HAS_PROLIF or prot_mol is None:
        return _nan_out()

    ref_ligs = [_ligand_for_plif(m) if m is not None else None for m in ref_mols]
    ref_valid = [j for j, x in enumerate(ref_ligs) if x is not None]
    if not ref_valid:
        return _nan_out()

    pose_ligs = [_ligand_for_plif(m) for m in pose_mols]
    valid_idx = [i for i, x in enumerate(pose_ligs) if x is not None]
    valid_ligs = [ref_ligs[j] for j in ref_valid] + [pose_ligs[i] for i in valid_idx]

    try:
        fp = plf.Fingerprint(interactions=PLIF_INTERACTIONS)
        fp.run_from_iterable(valid_ligs, prot_mol, n_jobs=1, progress=False)
        bvs = fp.to_bitvectors()
    except Exception:
        return _nan_out()

    if not bvs:
        return _nan_out()
    n_lead = len(ref_valid)
    ref_bvs = {j: bvs[k] for k, j in enumerate(ref_valid)}
    out = _nan_out()
    for k, i in enumerate(valid_idx, start=n_lead):
        for j, ref_bv in ref_bvs.items():
            try:
                val = float(DataStructs.TanimotoSimilarity(ref_bv, bvs[k]))
            except Exception:
                val = nan
            if single:
                out[i] = val
            else:
                out[i][j] = val
    return out


def clash_and_contacts(pose: Chem.Mol,
                       prot_xyz: np.ndarray,
                       prot_elem: list[str],
                       prot_resid: list[tuple[str, int]]) -> tuple[int, set]:
    """Return (n_clashes, set_of_contact_residue_ids)."""
    pose = Chem.RemoveHs(pose)
    lig_xyz = pose.GetConformer().GetPositions()   # bulk accessor (see symmetry_rmsd)
    lig_elem = [a.GetSymbol().upper() for a in pose.GetAtoms()]

    if len(prot_xyz) == 0:
        return 0, set()
    lo = lig_xyz.min(0) - 6.0
    hi = lig_xyz.max(0) + 6.0
    mask = ((prot_xyz >= lo).all(1) & (prot_xyz <= hi).all(1))
    p_xyz = prot_xyz[mask]
    p_elem = [e for e, m in zip(prot_elem, mask) if m]
    p_res = [r for r, m in zip(prot_resid, mask) if m]
    if len(p_xyz) == 0:
        return 0, set()

    diff = lig_xyz[:, None, :] - p_xyz[None, :, :]
    d = np.sqrt((diff * diff).sum(-1))

    n_clash = 0
    for i, le in enumerate(lig_elem):
        rl = _VDW.get(le, _DEFAULT_VDW)
        for j, pe in enumerate(p_elem):
            rp = _VDW.get(pe, _DEFAULT_VDW)
            if d[i, j] < CLASH_SCALE * (rl + rp):
                n_clash += 1

    contact_mask = (d < CONTACT_CUTOFF).any(0)
    contacts = {p_res[j] for j in np.where(contact_mask)[0]}
    return n_clash, contacts


# ───────────────────────────────────────────────────────────────────
# Per-pair worker
# ───────────────────────────────────────────────────────────────────

@dataclass
class PoseRecord:
    method: str
    protein: str
    ligand: str
    pose_file: str
    pose_name: str
    rank: int
    rmsd: float
    pb_rmsd: float
    pb_kabsch_rmsd: float
    pb_rmsd_within_2A: bool
    centroid_dist: float
    strain_energy: float
    n_clashes: int
    contact_recovery: float
    plif_recovery: float
    pb_valid: bool
    # "Twisted & turned" decomposition of pose vs crystal: how the same ligand
    # was re-oriented (rot_angle_deg) and internally re-bent (bestfit_rmsd =
    # superposed conformer RMSD — PoseBusters' pb_kabsch_rmsd when available,
    # else our own Kabsch fit — plus tfd and torsion deviations) relative to its
    # experimental conformation. Pairs with centroid_dist (translation) above.
    rot_angle_deg: float
    bestfit_rmsd: float
    tfd: float
    max_torsion_dev_deg: float
    mean_torsion_dev_deg: float
    n_torsions_flipped: int
    n_rot_bonds: int
    # EquiBind variant provenance (None for autodock/diffdock), so per-pose
    # metric differences can be sliced by pocket source / clamp / re-search.
    pocket_source: str | None = None
    clamp_variant: str | None = None
    refine_variant: str | None = None
    smina_affinity: float | None = None
    gnina_affinity: float | None = None
    # AutoDock variant/ranking provenance. Both ranks remain available for a
    # gnina/smina pose, while ``rank`` above is the axis actually used by this
    # report (Vina rank for raw; optimized rank for optimized variants).
    optimizer: str | None = None
    autodock_rank: int | None = None
    optimized_rank: int | None = None
    rank_metric: str | None = None
    autodock_affinity: float | None = None
    minimized_affinity: float | None = None
    cnn_score: float | None = None
    cnn_affinity: float | None = None
    # Uni-Dock families expose affinity-ordered ranks directly from the PB CSV.
    unidock_rank: int | None = None
    unidock_affinity: float | None = None
    unidock2_rank: int | None = None
    unidock2_affinity: float | None = None
    # Reference-ligand convention provenance (module docstring, "Reference-ligand
    # convention"). ``n_copies`` deposited copies were found in <ID>_ligands.sdf;
    # ``ref_copy_index`` is the one whose coordinates match <ID>_ligand.sdf;
    # ``nearest_copy_index`` is the copy j* every primary metric above was
    # measured against (== ref_copy_index under ``instance``). The
    # ``*_ref_instance`` twins always hold the reference-instance value, so they
    # equal the primary columns under ``instance`` and give the sensitivity arm
    # under ``nearest``.
    reference_convention: str = "instance"
    n_copies: int = 1
    ref_copy_index: int = 0
    nearest_copy_index: int = 0
    nearest_copy_is_ref: bool = True
    rmsd_ref_instance: float = float("nan")
    centroid_dist_ref_instance: float = float("nan")
    pb_rmsd_ref_instance: float = float("nan")
    bestfit_rmsd_ref_instance: float = float("nan")


def _reference_copy_index(copies_h: list[Chem.Mol], crystal_h: Chem.Mol,
                          tol: float = 1e-3) -> int | None:
    """Index of the copy whose heavy-atom coordinates coincide with ``crystal_h``
    (symmetry-aware in-place RMSD ≤ ``tol`` Å); ``None`` if no copy matches."""
    for j, c in enumerate(copies_h):
        try:
            if c.GetNumAtoms() != crystal_h.GetNumAtoms():
                continue
            r = symmetry_rmsd(c, crystal_h)
        except Exception:
            continue
        if not math.isnan(r) and r <= tol:
            return j
    return None


def _nearest_copy_index(rmsds: list[float], ref_copy_index: int) -> int:
    """argmin over the per-copy RMSDs with NaN as +inf. Ties: the reference
    instance if it attains the minimum, else the lowest record index (np.argmin);
    all-NaN falls back to the reference instance."""
    arr = np.nan_to_num(np.asarray(rmsds, dtype=float), nan=np.inf)
    if not np.isfinite(arr).any():
        return int(ref_copy_index)
    j = int(np.argmin(arr))
    if arr[ref_copy_index] == arr[j]:
        return int(ref_copy_index)
    return j


def process_pair(args) -> list[dict]:
    # Work tuple: (pair_key, records, benchmark_dir, root[, reference_convention]).
    # The 4-tuple form is accepted for backwards compatibility and means "instance".
    if len(args) >= 5:
        pair_key, group_records, benchmark_dir, root, reference_convention = args[:5]
    else:
        pair_key, group_records, benchmark_dir, root = args
        reference_convention = "instance"
    reference_convention = str(reference_convention or "instance")
    if reference_convention not in REFERENCE_CONVENTIONS:
        raise ValueError(f"unknown reference convention {reference_convention!r}")
    nearest = reference_convention == "nearest"
    protein, ligand = pair_key
    out: list[dict] = []

    pdb_id = protein
    cdir = Path(benchmark_dir) / pdb_id
    crystal_sdf = cdir / f"{pdb_id}_ligand.sdf"
    copies_sdf = cdir / f"{pdb_id}_ligands.sdf"
    protein_pdb = cdir / f"{pdb_id}_protein.pdb"
    if not crystal_sdf.exists() or not protein_pdb.exists():
        return out

    crystal = load_first_mol(crystal_sdf)
    if crystal is None:
        return out
    crystal_h = Chem.RemoveHs(crystal)

    # Every deposited copy of the ligand (<ID>_ligands.sdf; the single file if
    # absent), bond orders re-assigned from the reference so the substructure
    # matching in symmetry_rmsd / rigid_body_fit sees one topology. The copy that
    # coincides with <ID>_ligand.sdf is the reference instance; it is REPLACED by
    # ``crystal_h`` itself so the instance convention (and the *_ref_instance
    # twins) are evaluated on exactly the mol the frozen table used, not on a
    # re-read of the same coordinates. If no copy matches (not the case on this
    # dataset) the reference is prepended as record 0.
    copies_raw = load_all_mols(copies_sdf) if copies_sdf.exists() else []
    if not copies_raw:
        copies_raw = [crystal]
    copies_h = [Chem.RemoveHs(reassign_template(m, crystal_h) or m) for m in copies_raw]
    ref_copy_index = _reference_copy_index(copies_h, crystal_h)
    # Written copy indices are FILE record indices of <ID>_ligands.sdf; when the
    # reference had to be prepended (not the case on this dataset) it is written
    # as -1, which _reference_copy_source maps back to <ID>_ligand.sdf.
    copy_index_offset = 0
    if ref_copy_index is None:
        print(f"  [{pdb_id}] WARNING: no record of {copies_sdf.name} matches "
              f"{crystal_sdf.name}; prepending the reference instance as copy 0.")
        copies_h = [crystal_h] + copies_h
        ref_copy_index = 0
        copy_index_offset = 1
    else:
        copies_h[ref_copy_index] = crystal_h
    n_copies = len(copies_h)

    prot_xyz, prot_elem, prot_resid = load_protein_heavy_atoms(protein_pdb)
    # Native contacts per copy (contact recovery is evaluated at the pose's copy
    # j*). The receptor is the PoseBusters-shipped <ID>_protein.pdb with every
    # chain, so alternate copies are always covered.
    native_contacts: list[set] = []
    n_native: list[int] = []
    for j, c in enumerate(copies_h):
        if not nearest and j != ref_copy_index:
            native_contacts.append(set()); n_native.append(1)
            continue
        _, nc = clash_and_contacts(c, prot_xyz, prot_elem, prot_resid)
        native_contacts.append(nc); n_native.append(len(nc) or 1)

    prot_plif = _load_prolif_protein(protein_pdb)

    loaded: list[tuple[dict, Chem.Mol, Path]] = []
    for rec in group_records:
        pose_path = Path(rec["pose_file"])
        if not pose_path.is_absolute():
            pose_path = root / pose_path
        if not pose_path.exists():
            continue
        pose = load_first_mol(pose_path)
        if pose is None:
            continue
        # gnina-refined EquiBind poses carry their minimised affinity as a
        # <gnina_affinity> SDF tag (smina poses use <smina_affinity>, which the
        # PoseBusters CSV already exposes). The CSV does NOT carry the gnina tag,
        # so read it straight off the pose mol here — this is the score EquiBind
        # is ranked on downstream. Done before reassign_template so the tag is
        # read from the untouched supplier mol.
        if rec.get("gnina_affinity") in (None, "", "nan") and pose.HasProp("gnina_affinity"):
            try:
                rec = {**rec, "gnina_affinity": float(pose.GetProp("gnina_affinity"))}
            except (ValueError, TypeError):
                pass
        pose = reassign_template(pose, crystal_h) or pose
        loaded.append((rec, pose, pose_path))

    # One ProLIF run: every copy this convention scores leads the pose list and
    # yields one column; column j* is selected per pose below. Under ``instance``
    # only the reference copy is passed (None placeholders keep the column
    # layout), which is exactly the legacy single-reference run.
    plif_refs = (copies_h if nearest
                 else [c if j == ref_copy_index else None for j, c in enumerate(copies_h)])
    plif_mat = compute_plif_recovery(plif_refs,
                                     [p for _, p, _ in loaded],
                                     prot_plif)

    def _bestfit_at(pose, ref, pb_kabsch):
        # "Twisted & turned" vs the reference conformation. For the best-fit
        # "twist" RMSD prefer PoseBusters' symmetry-corrected superposed RMSD
        # (pb_kabsch_rmsd); fall back to our own Kabsch fit only when PoseBusters
        # is unavailable. The rotation angle ("turn") is always ours — PoseBusters
        # does not expose it.
        try:
            rot_angle, self_bestfit = rigid_body_fit(pose, ref)
        except Exception:
            rot_angle, self_bestfit = float("nan"), float("nan")
        return rot_angle, (pb_kabsch if not math.isnan(pb_kabsch) else self_bestfit)

    for (rec, pose, pose_path), plif_row in zip(loaded, plif_mat):
        # Reference copy for this pose: j* = argmin in-place RMSD over the copies
        # (D1/D2 of the nearest-copy plan) or the reference instance (default).
        try:
            if nearest:
                # One copy failing (an unsanitisable alternate record) must not
                # NaN the whole pose: score each copy in its own try and let
                # _nearest_copy_index treat NaN as +inf.
                r_all = [_safe_symmetry_rmsd(pose, c) for c in copies_h]
                j_star = _nearest_copy_index(r_all, ref_copy_index)
                rmsd = r_all[j_star]
                rmsd_ref = r_all[ref_copy_index]
            else:
                j_star = ref_copy_index
                rmsd = symmetry_rmsd(pose, crystal_h)
                rmsd_ref = rmsd
            ref_j = copies_h[j_star]
            cdist = centroid_distance(pose, ref_j)
            cdist_ref = cdist if j_star == ref_copy_index else centroid_distance(pose, crystal_h)
        except Exception:
            j_star = ref_copy_index
            ref_j = crystal_h
            rmsd, cdist = float("nan"), float("nan")
            rmsd_ref, cdist_ref = float("nan"), float("nan")

        # PoseBusters' canonical RMSD vs the reference copy j* of the crystal ligand.
        pb_rmsd, pb_kabsch, pb_within = posebusters_rmsd(pose, ref_j)

        strain = uff_strain(pose)
        rot_angle, bestfit = _bestfit_at(pose, ref_j, pb_kabsch)
        if j_star == ref_copy_index:
            pb_rmsd_ref, bestfit_ref = pb_rmsd, bestfit
        else:
            pb_rmsd_ref, pb_kabsch_ref, _ = posebusters_rmsd(pose, crystal_h)
            _, bestfit_ref = _bestfit_at(pose, crystal_h, pb_kabsch_ref)
        tfd, max_td, mean_td, n_flip, n_rot = torsion_metrics(pose, ref_j)
        try:
            n_clash, contacts = clash_and_contacts(pose, prot_xyz, prot_elem, prot_resid)
            recov = len(native_contacts[j_star] & contacts) / n_native[j_star]
        except Exception:
            n_clash, recov = 0, float("nan")
        try:
            plif_val = plif_row[j_star]
        except Exception:
            plif_val = float("nan")

        out.append(asdict(PoseRecord(
            method=rec["docking_method"],
            protein=protein,
            ligand=ligand,
            pose_file=str(pose_path),
            pose_name=rec.get("pose_name", pose_path.stem),
            rank=_pose_rank(rec, rec["docking_method"], str(pose_path)),
            rmsd=rmsd,
            pb_rmsd=pb_rmsd,
            pb_kabsch_rmsd=pb_kabsch,
            pb_rmsd_within_2A=pb_within,
            centroid_dist=cdist,
            strain_energy=strain,
            n_clashes=int(n_clash),
            contact_recovery=float(recov),
            plif_recovery=float(plif_val),
            pb_valid=bool(rec["pb_valid"]),
            rot_angle_deg=float(rot_angle),
            bestfit_rmsd=float(bestfit),
            tfd=float(tfd),
            max_torsion_dev_deg=float(max_td),
            mean_torsion_dev_deg=float(mean_td),
            n_torsions_flipped=int(n_flip),
            n_rot_bonds=int(n_rot),
            pocket_source=rec.get("pocket_source"),
            clamp_variant=rec.get("clamp_variant"),
            refine_variant=rec.get("refine_variant"),
            smina_affinity=rec.get("smina_affinity"),
            gnina_affinity=rec.get("gnina_affinity"),
            optimizer=_col_value(rec, "optimizer"),
            autodock_rank=_positive_rank(rec.get("autodock_rank")),
            optimized_rank=_positive_rank(rec.get("optimized_rank")),
            rank_metric=_col_value(rec, "rank_metric"),
            autodock_affinity=_finite_number(rec.get("autodock_affinity")),
            minimized_affinity=_finite_number(rec.get("minimized_affinity")),
            cnn_score=_finite_number(rec.get("cnn_score")),
            cnn_affinity=_finite_number(rec.get("cnn_affinity")),
            unidock_rank=_positive_rank(rec.get("unidock_rank")),
            unidock_affinity=_finite_number(rec.get("unidock_affinity")),
            unidock2_rank=_positive_rank(rec.get("unidock2_rank")),
            unidock2_affinity=_finite_number(rec.get("unidock2_affinity")),
            reference_convention=reference_convention,
            n_copies=int(n_copies),
            ref_copy_index=int(ref_copy_index - copy_index_offset),
            nearest_copy_index=int(j_star - copy_index_offset),
            nearest_copy_is_ref=bool(j_star == ref_copy_index),
            rmsd_ref_instance=float(rmsd_ref),
            centroid_dist_ref_instance=float(cdist_ref),
            pb_rmsd_ref_instance=float(pb_rmsd_ref),
            bestfit_rmsd_ref_instance=float(bestfit_ref),
        )))
    return out


# ───────────────────────────────────────────────────────────────────
# Selection helpers
# ───────────────────────────────────────────────────────────────────


def _oracle_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (method, protein, ligand): the pose with minimum RMSD."""
    valid = df.dropna(subset=["rmsd"])
    if valid.empty:
        return valid
    idx = valid.groupby(["method", "protein", "ligand"])["rmsd"].idxmin()
    return valid.loc[idx.dropna()].reset_index(drop=True)


def _top1_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    """One row per (method, protein, ligand): the tool's rank-1 pose.

    Only applies to RANKING_TOOLS; other methods are excluded.
    """
    ranked = df[df["method"].isin(RANKING_TOOLS)]
    return (ranked.sort_values("rank")
                  .groupby(["method", "protein", "ligand"]).head(1)
                  .reset_index(drop=True))


def _best_topn_per_pair(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """One row per (method, protein, ligand): the lowest-RMSD pose among the
    tool's top-``top_n`` ranked poses (ranks 1..top_n).

    This is the "best of the top-N" a user would actually consider — it shares
    the top-1 denominator (one representative pose per complex), so a top-1 vs
    top-N comparison is apples-to-apples. Only applies to RANKING_TOOLS.
    """
    ranked = df[df["method"].isin(RANKING_TOOLS)]
    ranked = ranked[pd.to_numeric(ranked["rank"], errors="coerce") <= top_n]
    valid = ranked.dropna(subset=["rmsd"])
    if valid.empty:
        return valid.reset_index(drop=True)
    idx = valid.groupby(["method", "protein", "ligand"])["rmsd"].idxmin()
    return valid.loc[idx.dropna()].reset_index(drop=True)


# ───────────────────────────────────────────────────────────────────
# Aggregation
# ───────────────────────────────────────────────────────────────────


def aggregate_oracle(df: pd.DataFrame) -> pd.DataFrame:
    """Oracle (best-RMSD pose) stats per method — all tools included."""
    rows = []
    for method, sub in df.groupby("method"):
        oracle = _oracle_per_pair(sub)
        n_pairs = sub.groupby(["protein", "ligand"]).ngroups
        poses_pp = sub.groupby(["protein", "ligand"]).size()
        row = {"method": method, "n_pairs": n_pairs, "n_poses": len(sub),
               "poses_per_pair_max": int(poses_pp.max()) if len(poses_pp) else 0,
               "median_rmsd_oracle": float(oracle["rmsd"].median(skipna=True))}
        for thr in RMSD_THRESHOLDS:
            row[f"oracle_rmsd_le_{thr}A_%"] = (
                100 * float((oracle["rmsd"].dropna() <= thr).mean())
                if len(oracle) else float("nan"))
        row["oracle_centroid_le_4A_%"] = (
            100 * float((oracle["centroid_dist"].dropna() <= CENTROID_THRESHOLD).mean())
            if len(oracle) else float("nan"))
        row["oracle_pb_valid_and_rmsd2_%"] = (
            100 * float(((oracle["rmsd"] <= 2.0) & oracle["pb_valid"]).mean())
            if len(oracle) else float("nan"))

        # PoseBusters canonical RMSD oracle (best-pb_rmsd pose per pair).
        pbv = sub.dropna(subset=["pb_rmsd"])
        pb_oracle = (pbv.loc[pbv.groupby(["protein", "ligand"])["pb_rmsd"].idxmin()]
                     if len(pbv) else pbv)
        row["median_pb_rmsd_oracle"] = (
            float(pb_oracle["pb_rmsd"].median(skipna=True)) if len(pb_oracle) else float("nan"))
        for thr in RMSD_THRESHOLDS:
            row[f"oracle_pb_rmsd_le_{thr}A_%"] = (
                100 * float((pb_oracle["pb_rmsd"].dropna() <= thr).mean())
                if len(pb_oracle) else float("nan"))
        rows.append(row)
    return _round_csv(pd.DataFrame(rows).set_index("method"))


def aggregate_top1(df: pd.DataFrame) -> pd.DataFrame:
    """Top-1 (rank-1 pose) stats per method — ranking tools only."""
    rows = []
    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        top1 = _top1_per_pair(sub)
        n_pairs = sub.groupby(["protein", "ligand"]).ngroups
        row = {"method": method, "n_pairs": n_pairs,
               "median_rmsd_top1": float(top1["rmsd"].median(skipna=True))}
        for thr in RMSD_THRESHOLDS:
            row[f"top1_rmsd_le_{thr}A_%"] = (
                100 * float((top1["rmsd"].dropna() <= thr).mean())
                if len(top1) else float("nan"))
        row["top1_centroid_le_4A_%"] = (
            100 * float((top1["centroid_dist"].dropna() <= CENTROID_THRESHOLD).mean())
            if len(top1) else float("nan"))
        row["top1_pb_valid_and_rmsd2_%"] = (
            100 * float(((top1["rmsd"] <= 2.0) & top1["pb_valid"]).mean())
            if len(top1) else float("nan"))

        # PoseBusters canonical RMSD of the rank-1 pose.
        row["median_pb_rmsd_top1"] = float(top1["pb_rmsd"].median(skipna=True)) if len(top1) else float("nan")
        for thr in RMSD_THRESHOLDS:
            row[f"top1_pb_rmsd_le_{thr}A_%"] = (
                100 * float((top1["pb_rmsd"].dropna() <= thr).mean())
                if len(top1) else float("nan"))
        rows.append(row)
    return _round_csv(pd.DataFrame(rows).set_index("method"))


def aggregate_by_rank(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Per-rank success stats for ranking tools (rank k = 1 .. top_n).

    For each k:
      success_pct_rmsd_2.0A     — % of all pairs where the pose at rank k ≤ 2 Å RMSD
      mean_rmsd / median_rmsd   — average RMSD of the rank-k pose
      oracle_top_k_pct_rmsd_2.0A — % of pairs where any of top-k has RMSD ≤ 2 Å
    """
    rows = []
    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        sub_ranked = sub[sub["rank"] <= top_n]
        n_all_pairs = sub.groupby(["protein", "ligand"]).ngroups

        for k in range(1, top_n + 1):
            at_k = (sub_ranked[sub_ranked["rank"] == k]
                    .groupby(["protein", "ligand"]).first()
                    .reset_index())

            topk_valid = sub_ranked[sub_ranked["rank"] <= k].dropna(subset=["rmsd"])
            if len(topk_valid):
                oracle_idx = topk_valid.groupby(["protein", "ligand"])["rmsd"].idxmin()
                oracle_k = topk_valid.loc[oracle_idx.dropna()]
            else:
                oracle_k = topk_valid

            row = {
                "method": method,
                "rank": k,
                "n_pairs_total": n_all_pairs,
                "n_pairs_at_rank": len(at_k),
                "mean_rmsd": float(at_k["rmsd"].mean(skipna=True)),
                "median_rmsd": float(at_k["rmsd"].median(skipna=True)),
            }
            for thr in RMSD_THRESHOLDS:
                row[f"success_pct_rmsd_{thr}A"] = (
                    100 * float((at_k["rmsd"].dropna() <= thr).sum()) / n_all_pairs
                    if n_all_pairs > 0 else float("nan"))
                row[f"oracle_top_k_pct_rmsd_{thr}A"] = (
                    100 * float((oracle_k["rmsd"].dropna() <= thr).sum()) / n_all_pairs
                    if n_all_pairs > 0 else float("nan"))
            rows.append(row)
    return _round_csv(pd.DataFrame(rows))


def aggregate_ifp_by_rank(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Per-rank interaction-fingerprint recovery vs the crystal.

    Companion to :func:`aggregate_by_rank`, but scoring how well the pose at each
    rank reproduces the crystal ligand's *interaction profile* rather than its
    geometry. Each pose is placed on its tool's ranking axis via
    :func:`_effective_rank`: native confidence rank for the RANKING_TOOLS (AutoDock
    Vina, DiffDock) and gnina-affinity rank for EquiBind, which has no native ranking
    (only its gnina-optimised variants carry a gnina affinity, so only those are
    included; raw/smina EquiBind would fall back to arbitrary generation order and
    are skipped). For each tool and rank k = 1 .. top_n it takes the pose at rank k
    for every receptor-ligand pair and reports:
      mean/median_plif_recovery    — ProLIF interaction-fingerprint Tanimoto vs crystal
      mean/median_contact_recovery — fraction of the crystal's contact residues recovered
    ``rank_kind`` records which ranking produced the axis. A declining curve means the
    tool's higher-ranked poses match the crystal interactions better than its
    lower-ranked ones. ``plif_recovery`` needs the crystal IFP (ProLIF), so it is NaN
    on datasets without a crystal (e.g. Orai) and for the poses ProLIF cannot
    fingerprint; those rows are dropped from the plif statistics only
    (``contact_recovery`` is always defined). ``n_plif_at_rank`` reports how many
    poses back each plif cell.
    """
    have_plif = "plif_recovery" in df.columns
    have_contact = "contact_recovery" in df.columns
    has_gnina = "gnina_affinity" in df.columns
    rows = []
    for method, sub in df.groupby("method"):
        is_ranking = method in RANKING_TOOLS
        gnina_frac = (pd.to_numeric(sub["gnina_affinity"], errors="coerce").notna().mean()
                      if has_gnina else 0.0)
        # EquiBind is ranked by gnina affinity (its only ranking signal); include a
        # variant only if that signal is actually present (the gnina-optimised runs).
        use_gnina = (not is_ranking and str(method).startswith("equibind")
                     and gnina_frac >= 0.5)
        if not (is_ranking or use_gnina):
            continue
        rank_kind = "native confidence" if is_ranking else "gnina affinity"
        sub = sub.assign(_rk=pd.to_numeric(_effective_rank(sub), errors="coerce"))
        sub_ranked = sub[sub["_rk"] <= top_n]
        n_all_pairs = sub.groupby(["protein", "ligand"]).ngroups

        for k in range(1, top_n + 1):
            at_k = (sub_ranked[sub_ranked["_rk"] == k]
                    .groupby(["protein", "ligand"]).first()
                    .reset_index())
            row = {
                "method": method,
                "rank": k,
                "rank_kind": rank_kind,
                "n_pairs_total": n_all_pairs,
                "n_pairs_at_rank": len(at_k),
            }
            if have_plif:
                plif = at_k["plif_recovery"].dropna()
                row["n_plif_at_rank"] = int(len(plif))
                row["mean_plif_recovery"] = float(plif.mean()) if len(plif) else float("nan")
                row["median_plif_recovery"] = float(plif.median()) if len(plif) else float("nan")
            if have_contact:
                con = at_k["contact_recovery"].dropna()
                row["mean_contact_recovery"] = float(con.mean()) if len(con) else float("nan")
                row["median_contact_recovery"] = float(con.median()) if len(con) else float("nan")
            rows.append(row)
    return pd.DataFrame(rows).round(4)


def aggregate_topn_within_thresholds(
        df: pd.DataFrame, top_n: int,
        thresholds: tuple[float, ...] = FINE_RMSD_THRESHOLDS,
        eq_df: "pd.DataFrame | None" = None,
        pb_valid_only: bool = False,
        rmsd_col: str = "rmsd") -> pd.DataFrame:
    """Fine-grained RMSD sweep of the top-ranked poses.

    Answers "how many of the top-N ranked poses of AutoDock Vina / DiffDock /
    EquiBind are within 1, 1.25, 1.5 … Å of the crystal ligand" from three
    complementary angles. For each tool × threshold t (long-form, one row each):

      top1_within_%       — % of the tool's complexes whose RANK-1 pose is within
                            t Å (the everyday "did the top pick land close?")
      best_topN_within_%  — % of complexes where ANY of the first ``top_n`` ranked
                            poses is within t (best-of-top-N / oracle within the
                            ranked set — the ceiling a perfect rescorer could hit)
      pose_within_%       — % of ALL pooled rank ≤ N poses within t (the literal
                            "how many of the n top-ranked poses are within t")

    AutoDock Vina / DiffDock use their native rank; ``eq_df`` optionally supplies the
    gnina-optimised EquiBind poses (which have no native ranking) — they are ranked
    by gnina AFFINITY (most-negative = rank 1, ``_gnina_affinity_rank``) and appended
    as method ``equibind``. Denominators: the two complex-level columns use the tool's
    full pair count (``n_pairs``); the pose-level column uses the pooled top-N pose
    count (``n_poses_topN``). The distance metric is ``rmsd_col``: by default the
    symmetry-corrected heavy-atom ``rmsd`` (no superposition); pass
    ``rmsd_col="bestfit_rmsd"`` for the best-fit (Kabsch) RMSD — the same pose after
    optimal superposition, which isolates conformer/shape error from placement.
    DiffDock's duplicated top pose is de-duplicated first.
    """
    thresholds = tuple(thresholds)
    rows = []

    def _emit(method_label: str, sub: pd.DataFrame, rank_series) -> None:
        sub = sub.copy()
        sub["_rk"] = pd.to_numeric(rank_series, errors="coerce")
        # PB-valid gate: a pose counts toward a threshold only if it is ALSO
        # PoseBusters-valid (the combined near-native AND physically-valid criterion,
        # matching oracle_pb_valid_and_rmsd2_%). Complex/pose denominators are unchanged
        # (all of the tool's complexes / all its top-N poses) — the gate only restricts
        # which poses may count as a "hit".
        if pb_valid_only:
            sub["pb_valid"] = _to_bool(sub["pb_valid"])   # shared coercion — keeps
            #   _topn_within_frames' per-complex mirror byte-identical to this curve
        n_pairs = sub.groupby(["protein", "ligand"]).ngroups
        ranked = sub[sub["_rk"] <= top_n]
        top1 = (ranked[ranked["_rk"] == 1]
                .groupby(["protein", "ligand"]).first().reset_index())
        top1_hit = top1[top1["pb_valid"]] if pb_valid_only else top1
        rk_hit = ranked[ranked["pb_valid"]] if pb_valid_only else ranked
        top1_rmsd = top1_hit[rmsd_col].dropna()
        rk_valid = rk_hit.dropna(subset=[rmsd_col])
        best_topn = (rk_valid.groupby(["protein", "ligand"])[rmsd_col].min()
                     if len(rk_valid) else pd.Series(dtype=float))
        pose_rmsd = rk_hit[rmsd_col].dropna()
        n_poses = len(ranked[rmsd_col].dropna())   # denominator: ALL top-N poses
        for t in thresholds:
            top1_n = int((top1_rmsd <= t).sum())
            best_n = int((best_topn <= t).sum())
            pose_n = int((pose_rmsd <= t).sum())
            rows.append({
                "method": method_label,
                "rmsd_threshold_A": float(t),
                "n_pairs": n_pairs,
                "n_poses_topN": n_poses,
                "top1_within_n": top1_n,
                "top1_within_%": (100.0 * top1_n / n_pairs if n_pairs else float("nan")),
                "best_topN_within_n": best_n,
                "best_topN_within_%": (100.0 * best_n / n_pairs if n_pairs else float("nan")),
                "pose_within_n": pose_n,
                "pose_within_%": (100.0 * pose_n / n_poses if n_poses else float("nan")),
            })

    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        sub = sub.copy()
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        sub = _dedup_ranked_poses(sub)          # collapse DiffDock's duplicate top pose
        _emit(str(method), sub, sub["rank"])

    # EquiBind: no native ranking, so rank by gnina affinity (gnina-opt variant).
    if eq_df is not None and not eq_df.empty:
        _emit("equibind", eq_df, _gnina_affinity_rank(eq_df))

    return _round_csv(pd.DataFrame(rows))


# Pre-specified RMSD thresholds at which fig 18's complex-level panel gets paired
# tests — a SMALL fixed family (2 Å canonical docking success + 1 Å high-accuracy),
# NOT one test per cumulative sweep point. The 21 grid thresholds are nested (a hit
# at 1 Å is a hit at every larger t), so testing them all would multiply-count the
# same complexes; we pre-register these two and Holm-pool the pairwise McNemar family
# across them.
_TOPN_WITHIN_TEST_THRESHOLDS = (1.0, 2.0)


def _topn_within_frames(
        df: pd.DataFrame, top_n: int,
        eq_df: "pd.DataFrame | None" = None,
        pb_valid_only: bool = False,
        rmsd_col: str = "rmsd") -> dict:
    """Per method → per-complex frame feeding fig 18's paired tests.

    Mirrors :func:`aggregate_topn_within_thresholds` EXACTLY (same de-dup of
    DiffDock's duplicate top pose, same gnina-affinity ranking for EquiBind, same
    optional PB-valid gate) but returns, instead of the aggregate counts, a frame
    indexed by every ``(protein, ligand)`` complex of the method with two float
    columns:
      ``top1_rmsd`` — the rank-1 pose's heavy-atom RMSD, or NaN when the rank-1 pose
                      is absent or (under the gate) PoseBusters-invalid;
      ``best_rmsd`` — the min RMSD among the method's top-N ranked (gated) poses, or
                      NaN when none qualifies.
    A hit at threshold ``t`` is ``col <= t`` (NaN ≤ t is False, so non-qualifying
    complexes count as misses at every t — identical denominators to the plotted
    curves). Returns ``{method: frame}`` for ``autodock`` / ``diffdock`` / ``equibind``.
    """
    frames: dict = {}

    def _emit(label: str, sub: pd.DataFrame, rank_series) -> None:
        sub = sub.copy()
        sub["_rk"] = pd.to_numeric(rank_series, errors="coerce")
        if pb_valid_only:
            sub["pb_valid"] = _to_bool(sub["pb_valid"])
        idx = sub.groupby(["protein", "ligand"]).size().index   # every complex
        ranked = sub[sub["_rk"] <= top_n]
        top1 = (ranked[ranked["_rk"] == 1]
                .groupby(["protein", "ligand"]).first())
        top1_hit = top1[top1["pb_valid"]] if pb_valid_only else top1
        rk_hit = ranked[ranked["pb_valid"]] if pb_valid_only else ranked
        frame = pd.DataFrame(index=idx)
        frame["top1_rmsd"] = (pd.to_numeric(top1_hit[rmsd_col], errors="coerce")
                              .reindex(idx) if not top1_hit.empty
                              else pd.Series(index=idx, dtype=float))
        best = (rk_hit.dropna(subset=[rmsd_col])
                .groupby(["protein", "ligand"])[rmsd_col].min())
        frame["best_rmsd"] = best.reindex(idx) if len(best) else float("nan")
        frames[label] = frame

    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        sub = sub.copy()
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        sub = _dedup_ranked_poses(sub)          # collapse DiffDock's duplicate top pose
        _emit(str(method), sub, sub["rank"])
    if eq_df is not None and not eq_df.empty:
        _emit("equibind", eq_df, _gnina_affinity_rank(eq_df))
    return frames


def _stats_topn_within(
        df: pd.DataFrame, top_n: int,
        eq_df: "pd.DataFrame | None" = None,
        pb_valid_only: bool = False,
        test_thresholds: tuple[float, ...] = _TOPN_WITHIN_TEST_THRESHOLDS,
        rmsd_col: str = "rmsd",
        figure_base: str = "18_topn_within_thresholds",
        metric_word: str = "within t Å") -> dict | None:
    """Paired significance tests for fig 18's complex-level panel at PRE-SPECIFIED
    RMSD thresholds (default 1 Å + 2 Å).

    Unit = the ``(protein, ligand)`` complex (one rank-1 pose per complex, so no
    pose-level pseudo-replication), paired across tools on the same complexes.
    Three parts, all on per-complex booleans:

      * ``rank1_across_tools`` / ``best_topN_across_tools`` — at each tested t, the
        rank-1 (solid curve) and best-of-top-N (dashed curve) HIT outcome across
        tools: Cochran's Q omnibus (complete cases) + pairwise exact McNemar. The
        pairwise p-values of each family are pooled ACROSS the tested thresholds into
        ONE Holm family (nested thresholds → one small pre-registered family, not 21
        tests) and returned as ``rank1_pairwise_holm`` / ``best_topN_pairwise_holm``.
      * ``ranking_headroom`` — within a tool, best-of-top-N vs rank-1. The rank-1 hit
        set is a strict SUBSET of best-of-top-N (rank-1 ∈ top-N), so the discordance
        is one-directional and a McNemar p is degenerate (0.5**b); the gain is
        reported as an EFFECT SIZE (recovered share b/N with a Wilson 95 % CI), OUT of
        the p-value family — mirroring :func:`topk_recovery_stats`.

    Per-tool rates carry 95 % Wilson CIs on each tool's OWN denominator (matching the
    plotted percentages) to drive the figure's error bars. Returns the sidecar dict,
    or None when < 2 tools carry the per-complex frames.
    """
    frames = _topn_within_frames(df, top_n, eq_df=eq_df, pb_valid_only=pb_valid_only,
                                 rmsd_col=rmsd_col)
    methods = [m for m in ("autodock", "diffdock", "equibind") if m in frames]
    if len(methods) < 2:
        return None
    ts = [float(t) for t in test_thresholds]

    def _hits(frame: pd.DataFrame, col: str, t: float) -> dict:
        """complex-id → hit-at-t flag (paired keys shared across tools)."""
        rec = pd.to_numeric(frame[col], errors="coerce") <= t
        return {f"{p}||{l}": bool(v) for (p, l), v in zip(frame.index, rec)}

    def _across(col: str):
        """Cochran's Q per threshold (complete cases) + pairwise McNemar (raw p);
        Holm pooled across thresholds. Returns (per_threshold, pooled_pairwise)."""
        per_t, fam = [], []
        for t in ts:
            maps = {m: _hits(frames[m], col, t) for m in methods}
            common = sorted(set.intersection(*[set(mp) for mp in maps.values()]))
            n = len(common)
            rec = {"rmsd_threshold_A": t, "n_complete": n}
            if n >= _MIN_UNITS_STATS:
                M = np.array([[maps[m][c] for m in methods] for c in common], float)
                Q, pQ, dQ = su.cochran_q(M)
                rec["omnibus"] = {"Q": Q, "p": pQ, "df": dQ}
                for i in range(len(methods)):
                    for j in range(i + 1, len(methods)):
                        a, b = methods[i], methods[j]
                        va = np.fromiter((maps[a][c] for c in common), bool, n)
                        vb = np.fromiter((maps[b][c] for c in common), bool, n)
                        n10, n01, p = su.mcnemar_exact(va, vb)
                        # orient 'a' as the higher-rate condition
                        if int(vb.sum()) > int(va.sum()):
                            a, b, n10, n01 = b, a, n01, n10
                        fam.append({"rmsd_threshold_A": t, "a": a, "b": b,
                                    "a_wins": n10, "b_wins": n01, "n": n,
                                    "p_raw": float(p)})
            else:
                rec["omnibus"] = {"status": "n too small — exploratory", "n": n}
            per_t.append(rec)
        if fam:
            for pr, pa in zip(fam, su.holm([r["p_raw"] for r in fam])):
                pr["p_holm"] = float(pa)
                pr["star"] = su.p_stars(pa)
        return per_t, fam

    r1_per_t, r1_fam = _across("top1_rmsd")
    best_per_t, best_fam = _across("best_rmsd")

    headroom = []
    for m in methods:
        fr = frames[m]
        N = int(len(fr))
        t1_rmsd = pd.to_numeric(fr["top1_rmsd"], errors="coerce")
        bt_rmsd = pd.to_numeric(fr["best_rmsd"], errors="coerce")
        for t in ts:
            t1 = t1_rmsd <= t
            bt = bt_rmsd <= t
            k1, kb = int(t1.sum()), int(bt.sum())
            b = int((bt & ~t1).sum())          # recovered by a lower rank, missed by rank-1
            lo1, hi1 = su.wilson_ci(k1, N)
            lob, hib = su.wilson_ci(kb, N)
            loh, hih = su.wilson_ci(b, N)
            headroom.append({
                "method": m, "rmsd_threshold_A": t, "n": N,
                "top1_k": k1, "top1_rate": (k1 / N if N else float("nan")),
                "top1_ci": [lo1, hi1],
                "best_topN_k": kb, "best_topN_rate": (kb / N if N else float("nan")),
                "best_topN_ci": [lob, hib],
                "headroom_k": b, "headroom_share": (b / N if N else float("nan")),
                "headroom_ci": [loh, hih], "nested": True,
            })

    return {
        "figure": (f"{figure_base}_pbvalid_complex.png" if pb_valid_only
                   else f"{figure_base}_complex.png"),
        "metric": (f"rank-1 / best-of-top-N pose {metric_word} AND PoseBusters-valid"
                   if pb_valid_only else f"rank-1 / best-of-top-N pose {metric_word}"),
        "unit": ("(protein, ligand) complex — rank-1 / best-of-top-N hit boolean; "
                 "paired across tools on shared complexes"),
        "test": ("per-threshold Cochran's Q + pairwise exact McNemar (Holm pooled "
                 "across the pre-specified thresholds) on the across-tool rank-1 and "
                 "best-of-top-N hit; within-tool ranking headroom (best-of-top-N vs "
                 "rank-1) reported as a Wilson-CI effect size (nested → no p); 95% "
                 "Wilson CIs on every rate"),
        "test_thresholds_A": ts,
        "top_n": int(top_n),
        "methods": list(methods),
        "rank1_across_tools": r1_per_t,
        "rank1_pairwise_holm": r1_fam,
        "best_topN_across_tools": best_per_t,
        "best_topN_pairwise_holm": best_fam,
        "ranking_headroom": headroom,
    }


# Metric the variant-collapse selectors rank on: the combined docking-success criterion
# PB-Valid AND RMSD ≤ 2 Å (a pose must be BOTH near-native AND physically valid). Both this
# and the RMSD-only rate are produced by aggregate_oracle; the fallback covers a summary that
# somehow lacks the combined column.
BEST_VARIANT_METRIC = "oracle_pb_valid_and_rmsd2_%"
BEST_VARIANT_METRIC_FALLBACK = "oracle_rmsd_le_2.0A_%"


def _variant_rank_col(oracle_sum: pd.DataFrame) -> str | None:
    """The oracle-summary column to rank a tool's variants on: the combined PB-Valid AND
    RMSD ≤ 2 Å success rate, falling back to the RMSD-only rate, or None when neither is
    present."""
    if BEST_VARIANT_METRIC in oracle_sum.columns:
        return BEST_VARIANT_METRIC
    if BEST_VARIANT_METRIC_FALLBACK in oracle_sum.columns:
        return BEST_VARIANT_METRIC_FALLBACK
    return None


def _select_best_equibind(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind methods plus only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest PB-Valid AND RMSD ≤ 2 Å success rate
    (``oracle_pb_valid_and_rmsd2_%``, falling back to ``oracle_rmsd_le_2.0A_%``; see
    _variant_rank_col) from the oracle summary. Returns the filtered
    frame and the chosen variant's method key (``None`` if no EquiBind variant is
    present, in which case ``df`` is returned unchanged). AutoDock/DiffDock rows
    are always retained.
    """
    methods = df["method"].astype(str)
    eq_mask = methods.str.startswith("equibind")
    if not eq_mask.any():
        return df, None

    oracle_sum = aggregate_oracle(df)
    col = _variant_rank_col(oracle_sum)
    eq_keys = [m for m in oracle_sum.index if str(m).startswith("equibind")]
    scores = (oracle_sum.loc[eq_keys, col].astype(float).dropna()
              if col else pd.Series(dtype=float))
    best_variant = (str(scores.idxmax()) if not scores.empty
                    else sorted(methods[eq_mask].unique())[0])

    keep = (~eq_mask) | (methods == best_variant)
    return df[keep].reset_index(drop=True), best_variant


def _select_best_diffdock(df: pd.DataFrame, forced: str | None = None) -> tuple[pd.DataFrame, str | None]:
    """Keep non-DiffDock methods plus only the single best DiffDock variant.

    "Best" = the DiffDock optimizer variant (diffdock / diffdock_smina /
    diffdock_gnina) with the highest PB-Valid AND RMSD ≤ 2 Å success rate
    (``oracle_pb_valid_and_rmsd2_%``, falling back to ``oracle_rmsd_le_2.0A_%``; see
    _variant_rank_col). The winner is relabelled to the canonical
    ``diffdock`` method so it stays a ranking tool; the original variant key is
    returned (for the 'DiffDock*' label / print). AutoDock/EquiBind are retained.
    ``None`` when no DiffDock variant is present (``df`` unchanged)."""
    methods = df["method"].astype(str)
    dd_mask = methods.str.startswith("diffdock")
    if not dd_mask.any():
        return df, None

    present = sorted(methods[dd_mask].unique())
    if forced and forced in present:
        best_variant = forced
    else:
        if forced:
            print(f"  [collapse-diffdock-variant] {forced!r} not present among "
                  f"{present} — falling back to oracle ranking.")
        oracle_sum = aggregate_oracle(df)
        col = _variant_rank_col(oracle_sum)
        dd_keys = [m for m in oracle_sum.index if str(m).startswith("diffdock")]
        scores = (oracle_sum.loc[dd_keys, col].astype(float).dropna()
                  if col else pd.Series(dtype=float))
        best_variant = (str(scores.idxmax()) if not scores.empty else present[0])

    keep = (~dd_mask) | (methods == best_variant)
    out = df[keep].reset_index(drop=True)
    # Canonicalise the winner to 'diffdock' so it ranks/colours like the tool.
    winner = out["method"] == best_variant
    out.loc[winner, "method"] = "diffdock"
    # Re-derive rank from each pose filename now that the winner is a ranking
    # tool. A cached per_pose_metrics.csv scored the optimizer variants
    # (diffdock_smina / diffdock_gnina) under their own method key, which older
    # parse_rank did not recognise, so their cached rank is a flat 999. Recompute
    # it here (idempotent for the base 'diffdock', which already parsed) so the
    # Part-B ranking metrics see the real rank ≤ top_n poses without a --force
    # rescore of the whole per-pose pass.
    if "pose_file" in out.columns:
        out.loc[winner, "rank"] = [
            parse_rank("diffdock", str(p)) for p in out.loc[winner, "pose_file"]
        ]
    return out, best_variant


# Independent AutoDock arms. They ARE drawn in the headline, each as its own
# engine, but must never be pinned into the canonical "autodock" slot: doing so
# relabels them to "autodock" and silently redefines the reported AutoDock bar.
_SEPARATE_AUTODOCK_ARMS = ("autodock_mgltools",)

# The one MGLTools/ADFRsuite arm that IS allowed into the canonical "autodock"
# slot, because it is the reported pipeline rather than a control. Kept as an
# explicit allowlist rather than dropping the prefix above, so exh18/32/64/92 and
# the non-dominant gnina arms stay barred from silently redefining the AutoDock bar.
# The reported AutoDock pipeline, as a raw / optimised pair at ONE search effort so
# every "did optimisation help?" contrast is measured within a fixed exhaustiveness.
# Change these two names to move the dominant arm; the spec tables below and
# _oracle_variant_specs all read them rather than repeating the keys.
_DOMINANT_AUTODOCK_RAW = "autodock_mgltools_exh128"
_DOMINANT_AUTODOCK_OPT = "autodock_mgltools_exh128_gnina"
_PINNABLE_AUTODOCK_ARMS = (_DOMINANT_AUTODOCK_OPT,)


def _is_separate_autodock_arm(name: str) -> bool:
    """True for an AutoDock arm that renders as its own engine, not the Vina slot."""
    if str(name) in _PINNABLE_AUTODOCK_ARMS:
        return False
    return str(name).startswith(_SEPARATE_AUTODOCK_ARMS)


def _select_autodock_arm(df: pd.DataFrame,
                         forced: str | None = None) -> tuple[pd.DataFrame, str | None]:
    """Keep non-AutoDock methods plus a single pinned AutoDock Vina variant.

    Counterpart to :func:`_select_best_diffdock` for the AutoDock family, added
    because the presentation keep-set below hard-codes the raw ``autodock`` slot.
    Without this, every collapsed headline figure is drawn on raw Vina even when the
    chapter reports the gnina-rescored arm, and there was no way to say otherwise on
    the command line.

    ``forced`` names one of the Vina-scored variants (``autodock``,
    ``autodock_gnina``, ``autodock_gnina_refinement``). The winner is relabelled to
    the canonical ``autodock`` key so it keeps its ranking-tool behaviour and passes
    the keep-set; the original key is returned for the label and the provenance line.
    Vinardo variants are a separate scoring family and are left untouched, exactly as
    the keep-set already treats them. Returns ``(df, None)`` unchanged when no
    AutoDock variant is present.

    NOTE the relabel is rank-safe: for ``autodock_gnina`` the cached ``rank`` column
    is identical to ``optimized_rank``, i.e. it is already the gnina re-ranking, so
    no rank re-derivation is needed here (unlike the DiffDock optimizer variants,
    whose cached rank is a flat 999)."""
    methods = df["method"].astype(str)
    # Vina-scored AutoDock only. autodock_vinardo* is a different scoring function
    # and the keep-set carries it as its own engine column.
    ad_mask = (methods.str.startswith("autodock")
               & ~methods.str.startswith("autodock_vinardo")
               & ~methods.map(_is_separate_autodock_arm))
    if not ad_mask.any():
        return df, None

    present = sorted(methods[ad_mask].unique())
    if forced and _is_separate_autodock_arm(forced):
        # Pinning a sidecar arm here would relabel it to the canonical "autodock"
        # key and make it THE headline AutoDock bar. That arm is a converter /
        # search-effort control, not the reported pipeline — refuse explicitly
        # rather than silently redefining every collapsed figure.
        print(f"  [collapse-autodock-variant] {forced!r} is an independent AutoDock arm "
              "and cannot be pinned as the headline AutoDock slot — leaving it unchanged.")
        return df, None
    if not forced:
        return df, None
    if forced not in present:
        print(f"  [collapse-autodock-variant] {forced!r} not present among {present} "
              "— leaving the AutoDock slot unchanged.")
        return df, None

    keep = (~ad_mask) | (methods == forced)
    out = df[keep].reset_index(drop=True)
    out.loc[out["method"] == forced, "method"] = "autodock"
    return out, forced


def _select_presentation_tools(
        df: pd.DataFrame, equibind_variant: str | None) -> pd.DataFrame:
    """Keep the canonical presentation frame used by collapsed presentation plots.

    ``--collapse-plots-only`` writes its all-variant tables before calling the
    family selectors. Those selectors intentionally retain methods outside their
    own family, which used to be sufficient when raw AutoDock was the only
    non-DiffDock/non-EquiBind method. Once AutoDock optimizer/scoring variants and
    UniDock were added, those methods leaked into ``oracle_summary.csv`` and the
    headline plots. The collapsed presentation keeps the full-engine comparison:
    raw AutoDock (Vina) and raw AutoDock Vinardo (distinct scoring functions on the
    same engine), the canonicalized DiffDock winner, the selected EquiBind winner,
    and the single-variant engines Uni-Dock / Uni-Dock2. The post-hoc optimizer
    variants (gnina, refinement, smina) stay collapsed away — they are compared in
    their own per-family figures, not the cross-engine headline.
    """
    keep = {"autodock", "autodock_vinardo", "diffdock", "unidock", "unidock2",
            # Independent AutoDock arms: their own headline bars. They are still
            # barred from being pinned INTO the canonical "autodock" slot (see
            # _select_autodock_arm), so they can never redefine the AutoDock bar.
            # The MGLTools arm is represented here by its exhaustiveness-32 baseline
            # (the unsuffixed key) plus the pre-existing exh64 entry. The 18/92 arms are
            # deliberately NOT added: this frame feeds the cross-tool omnibus and the
            # Holm-corrected pairwise family, where four sampling points of ONE arm are
            # pseudoreplicates, not four engines — they would inflate the comparison
            # count and depress Kendall's W. The full sweep is carried by
            # oracle_summary_all_variants.csv and the 09f all-variants figure instead.
            # The exhaustiveness controls are NOT carried here any more. With the
            # dominant arm pinned into the "autodock" slot they would put three
            # points of ONE ladder on the cross-tool axis, which is the
            # pseudoreplication this comment already warned about. The sweep is
            # carried by oracle_summary_all_variants.csv and the 09f all-variants
            # figure instead.
            }
    if equibind_variant:
        keep.add(str(equibind_variant))
    return df[df["method"].astype(str).isin(keep)].reset_index(drop=True)


# ───────────────────────────────────────────────────────────────────
# Plots — Part A: oracle comparison (all tools)
# ───────────────────────────────────────────────────────────────────


# ───────────────────────────────────────────────────────────────────
# Statistical tests for the headline oracle / success / paper figures.
# Unit of analysis = one (protein, ligand) complex (poses are aggregated to the
# oracle / top-1 representative first); the design is paired across tools, so we
# use the shared paired bundles (Friedman/Kendall-W, Cochran-Q + exact McNemar,
# Wilson CIs) from stats_utils. Every caller wraps these in try/except so a stats
# failure degrades to the current test-free figure. Results are also written to
# posebusters_pose_comparison_stats.json in the out-dir. (See the plan doc.)
# ───────────────────────────────────────────────────────────────────
_MIN_UNITS_STATS = 5   # below this the paired tests are labelled 'exploratory'


def _stats_oracle_rmsd_paired(df: pd.DataFrame) -> dict | None:
    """Per-complex ORACLE (best) RMSD across tools → Friedman + Kendall's W with
    Wilcoxon signed-rank pairwise (Holm). Complete-case (protein, ligand) rows.

    Drives figures 01_oracle_rmsd_cdf and 02_oracle_rmsd_boxplot.
    """
    oracle = _oracle_per_pair(df)
    if oracle.empty:
        return None
    piv = oracle.pivot_table(index=["protein", "ligand"], columns="method",
                             values="rmsd")
    methods = list(piv.columns)
    complete = piv[methods].dropna()
    payload = {
        "figures": ["01_oracle_rmsd_cdf.png", "02_oracle_rmsd_boxplot.png"],
        "measure": "per-complex oracle (best) heavy-atom RMSD vs crystal (Å)",
        "unit": "(protein, ligand) complex — poses aggregated to the oracle first; complete cases across all tools",
        "test": "Friedman + Kendall's W omnibus; Wilcoxon signed-rank (Holm) pairwise",
        "methods": [str(m) for m in methods],
        "n_complete": int(len(complete)),
    }
    if len(complete) < _MIN_UNITS_STATS or len(methods) < 2:
        payload["status"] = "n too small — exploratory"
        return payload
    res = su.paired_continuous({str(m): complete[m].to_numpy(float) for m in methods})
    payload["status"] = "ok"
    payload["omnibus"] = res["omnibus"]
    payload["pairwise"] = res["pairwise"]
    payload["medians"] = res["medians"]
    return payload


def _annotate_oracle_rmsd_stats(ax, payload, x, y, ha, va) -> None:
    """Draw the oracle-RMSD Friedman/Kendall-W + pairwise-star box on an axis."""
    if not payload:
        return
    box = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.88)
    if payload.get("status") != "ok":
        ax.text(x, y, f"cross-tool test: {payload.get('status', 'n/a')} "
                f"(n={payload.get('n_complete', 0)})", transform=ax.transAxes,
                ha=ha, va=va, fontsize=7.3, color="0.35", bbox=box)
        return
    om = payload["omnibus"]
    lines = [f"Friedman χ²={om['chi2']:.1f}, {su.fmt_p(om['p'])} "
             f"(Kendall W={om['kendall_w']:.2f}, n={om['n']})"]
    for pr in payload["pairwise"]:
        a = TOOL_LABEL.get(pr["a"], pr["a"]); b = TOOL_LABEL.get(pr["b"], pr["b"])
        lines.append(f"{a} vs {b}: {su.p_stars(pr.get('p_holm'))} "
                     f"({su.fmt_p(pr.get('p_holm'))}, r={pr['rank_biserial']:+.2f})")
    ax.text(x, y, "\n".join(lines), transform=ax.transAxes, ha=ha, va=va,
            fontsize=7.3, color="0.20", bbox=box, linespacing=1.3)


def _stats_success_paired(df: pd.DataFrame, *, valid: bool) -> dict | None:
    """Per-complex success across tools (paired proportions) + within-tool
    oracle-vs-top-1 (exact McNemar).

    success = RMSD ≤ 2 Å  (valid=False, fig 03)
            = RMSD ≤ 2 Å AND PB-valid  (valid=True, fig 10)
    Unit = (protein, ligand) complex. Per-method Wilson CIs are on each method's
    own denominator (matching the plotted percentages); the cross-tool omnibus
    uses complete cases.
    """
    def _succ(sub):
        s = sub["rmsd"] <= 2.0
        if valid:
            s = s & sub["pb_valid"].astype(bool)
        return s.astype(float)

    oracle = _oracle_per_pair(df).copy()
    if oracle.empty:
        return None
    oracle["succ"] = _succ(oracle)
    orc = oracle.pivot_table(index=["protein", "ligand"], columns="method",
                             values="succ")
    top1 = _top1_per_pair(df).copy()
    if not top1.empty:
        top1["succ"] = _succ(top1)
        t1 = top1.pivot_table(index=["protein", "ligand"], columns="method",
                              values="succ")
    else:
        t1 = pd.DataFrame()

    payload = {
        "figure": ("10_pb_valid_rmsd2_success_bars.png" if valid
                   else "03_oracle_vs_top1_success.png"),
        "metric": ("RMSD ≤ 2 Å AND PB-valid" if valid else "RMSD ≤ 2 Å"),
        "unit": "(protein, ligand) complex — per-complex success boolean",
        "test": ("paired proportions across tools (Cochran's Q + exact McNemar, "
                 "Holm) on the oracle success; within-tool oracle-vs-top-1 exact "
                 "McNemar; 95% Wilson CIs"),
        "methods": [str(m) for m in orc.columns],
    }
    per = {}
    for m in orc.columns:
        col = orc[m].dropna(); k = int(col.sum()); n = int(col.size)
        lo, hi = su.wilson_ci(k, n)
        d = {"oracle_k": k, "oracle_n": n,
             "oracle_rate": (k / n if n else float("nan")),
             "oracle_lo": lo, "oracle_hi": hi}
        if not t1.empty and m in t1.columns:
            c2 = t1[m].dropna(); k2 = int(c2.sum()); n2 = int(c2.size)
            lo2, hi2 = su.wilson_ci(k2, n2)
            d.update(top1_k=k2, top1_n=n2,
                     top1_rate=(k2 / n2 if n2 else float("nan")),
                     top1_lo=lo2, top1_hi=hi2)
        per[str(m)] = d
    payload["per_method"] = per

    complete = orc[list(orc.columns)].dropna()
    if len(complete) >= _MIN_UNITS_STATS and orc.shape[1] >= 2:
        payload["cross_tool_oracle"] = su.paired_proportions(
            {str(m): complete[m].to_numpy(float) for m in orc.columns})
    else:
        payload["cross_tool_oracle"] = {"status": "n too small — exploratory",
                                        "n": int(len(complete))}

    wt = []
    for m in orc.columns:
        if t1.empty or m not in t1.columns:
            continue
        pair = pd.concat([orc[m].rename("o"), t1[m].rename("t")], axis=1).dropna()
        if len(pair) < _MIN_UNITS_STATS:
            wt.append({"method": str(m), "n": int(len(pair)),
                       "status": "n too small — exploratory"})
            continue
        n10, n01, p = su.mcnemar_exact(pair["o"].to_numpy(), pair["t"].to_numpy())
        wt.append({"method": str(m), "n": int(len(pair)), "oracle_wins": n10,
                   "top1_wins": n01, "p": float(p), "star": su.p_stars(p)})
    payload["within_tool_oracle_vs_top1"] = wt
    return payload


def _annotate_success_dumbbell(fig, ax, methods, oracle_vals, top1_vals, payload,
                               footnote_pairwise: bool = True) -> None:
    """Overlay Wilson-CI error bars, per-row within-tool McNemar stars, and a
    cross-tool omnibus footnote on the oracle-vs-top-1 dumbbell charts (03 / 10).

    ``footnote_pairwise`` inlines the significant cross-tool pairs into the
    footnote (fine for 3 tools). Set it False when many rows would make that list
    overflow (fig 09e's 5 variants) — the omnibus stays on the figure and the full
    pairwise matrix is recoverable from the stats sidecar.
    """
    if not payload:
        return
    yv = np.arange(len(methods))[::-1]
    per = payload.get("per_method", {})
    wt = {d["method"]: d for d in payload.get("within_tool_oracle_vs_top1", [])}
    for yi, m, ov, tv in zip(yv, methods, oracle_vals, top1_vals):
        d = per.get(str(m), {})
        if not pd.isna(ov) and "oracle_lo" in d:
            ax.errorbar(ov, yi,
                        xerr=[[max(0.0, ov - d["oracle_lo"] * 100)],
                              [max(0.0, d["oracle_hi"] * 100 - ov)]],
                        fmt="none", ecolor="0.30", elinewidth=1.1, capsize=2.5, zorder=2)
        if not pd.isna(tv) and "top1_lo" in d:
            ax.errorbar(tv, yi,
                        xerr=[[max(0.0, tv - d["top1_lo"] * 100)],
                              [max(0.0, d["top1_hi"] * 100 - tv)]],
                        fmt="none", ecolor="0.55", elinewidth=1.1, capsize=2.5, zorder=2)
        w = wt.get(str(m))
        if w and w.get("star") and not pd.isna(ov) and not pd.isna(tv):
            ax.text((ov + tv) / 2, yi + 0.24, f"{w['star']} ({su.fmt_p(w['p'])})",
                    ha="center", va="bottom", fontsize=7.3, color="#333",
                    fontweight="bold")
    ct = payload.get("cross_tool_oracle", {})
    if str(ct.get("status", "")).startswith("n too small"):
        note = f"Across-tool oracle test: n too small — exploratory (n={ct.get('n', 0)})"
    elif "omnibus" in ct:
        om = ct["omnibus"]
        if footnote_pairwise:
            sig = [f"{TOOL_LABEL.get(pr['a'], pr['a'])}>{TOOL_LABEL.get(pr['b'], pr['b'])} {pr['star']}"
                   for pr in ct.get("pairwise", []) if pr.get("star") not in ("", "ns")]
            tail = ("; " + ", ".join(sig)) if sig else ""
            unit = "tool"
        else:
            n_sig = sum(1 for pr in ct.get("pairwise", []) if pr.get("star") not in ("", "ns"))
            tail = (f"; {n_sig} of {len(ct.get('pairwise', []))} pairs significant "
                    "(Holm) — see stats sidecar") if ct.get("pairwise") else ""
            unit = "variant"
        note = f"Across-{unit} oracle success: Cochran's Q={om['Q']:.1f}, {su.fmt_p(om['p'])}{tail}"
    else:
        note = ""
    # Pairwise cross-tool detail (incl. non-significant) is in the stats sidecar.
    note = (note + "  ·  within-tool ●vs○: exact McNemar; bars = 95% Wilson CI").strip()
    fig.text(0.5, 0.012, note, ha="center", va="bottom", fontsize=7.0, color="0.30")


def _stats_vs_paper(df: pd.DataFrame, oracle_sum: pd.DataFrame,
                    top1_sum: pd.DataFrame) -> dict | None:
    """One-sample binomial tests of our primary-pose success rate vs the paper's
    FIXED published rate (fig 11). Selection mirrors plot_vs_posebusters_paper:
    top-1 for ranking tools, best pose for EquiBind. The complex sets differ, so
    this is descriptive — noted in the payload and the figure footnote.
    """
    from scipy.stats import binomtest
    methods = list(oracle_sum.index)
    chosen = [m for m in ("autodock", "diffdock") if m in methods]
    eq = [m for m in methods if m.startswith("equibind")]
    if eq:
        vcol = "oracle_pb_valid_and_rmsd2_%"
        best_eq = max(eq, key=lambda m: (float(oracle_sum.loc[m, vcol])
                                         if vcol in oracle_sum.columns else -1.0))
        chosen.append(best_eq)
    if not chosen:
        return None
    orc = _oracle_per_pair(df)
    t1 = _top1_per_pair(df)
    payload = {
        "figure": "11_vs_posebusters_paper.png",
        "test": "one-sample exact binomial vs the paper's fixed published rate",
        "note": ("the paper's value is a fixed published proportion (PoseBusters "
                 "Benchmark) with a different complex set from this run, so the "
                 "test is descriptive, not a like-for-like comparison"),
        "methods": {},
    }
    for m in chosen:
        ranking = m in RANKING_TOOLS
        sel_df = t1[t1["method"] == m] if ranking else orc[orc["method"] == m]
        n = int(len(sel_df))
        entry = {"selection": ("top-1" if ranking else "best pose*"), "n": n, "metrics": {}}
        ref = POSEBUSTERS_PAPER_BENCHMARK.get(_paper_method_name(m), {})
        for metric, mask, paperkey in (
            ("rmsd_le_2A", sel_df["rmsd"] <= 2.0, "rmsd2"),
            ("rmsd_le_2A_and_pb_valid",
             (sel_df["rmsd"] <= 2.0) & sel_df["pb_valid"].astype(bool), "rmsd2_valid"),
        ):
            k = int(mask.sum())
            lo, hi = su.wilson_ci(k, n)
            paper = ref.get(paperkey)
            md = {"k": k, "n": n,
                  "our_rate_%": (round(100 * k / n, PCT_DECIMALS) if n else None),
                  "wilson_lo_%": round(100 * lo, PCT_DECIMALS), "wilson_hi_%": round(100 * hi, PCT_DECIMALS),
                  "paper_%": paper}
            if paper is not None and n > 0 and 0.0 <= paper <= 100.0:
                bt = binomtest(k, n, paper / 100.0)
                md["binom_p"] = float(bt.pvalue); md["star"] = su.p_stars(bt.pvalue)
            else:
                md["binom_p"] = None; md["star"] = ""
            entry["metrics"][metric] = md
        payload["methods"][str(m)] = entry
    return payload


def plot_oracle_rmsd_cdf(df: pd.DataFrame, out: Path, stats: dict | None = None) -> None:
    """Oracle RMSD CDF — best pose per pair, all tools."""
    oracle = _oracle_per_pair(df)
    fig, ax = plt.subplots(figsize=(9, 5.8))
    for method, sub in oracle.groupby("method"):
        vals = sub["rmsd"].dropna().sort_values().values
        if len(vals) == 0:
            continue
        ys = np.arange(1, len(vals) + 1) / len(vals) * 100
        ax.plot(vals, ys,
                label=f"{TOOL_LABEL.get(method, method)} (n={len(vals)})",
                color=TOOL_COLORS.get(method), lw=2)
    ax.axvline(2.0, color="grey", ls="--", alpha=0.6, label="2 Å")
    ax.set_xlim(0, 15)
    ax.set_xlabel("Oracle (best) heavy-atom RMSD vs crystal (Å)")
    ax.set_ylabel("Cumulative % of pairs")
    ax.set_title(_vt("Oracle RMSD CDF — best pose per pair per method"))
    ax.legend(); ax.grid(alpha=0.3)
    # Footnote: define "cumulative % of pairs" and how each method's curve is built.
    footnote = "\n".join([
        "Cumulative % of pairs — per method, the oracle (lowest-RMSD) pose is taken for every receptor–ligand pair; at a given RMSD x the",
        "curve gives the percentage of that method's pairs whose oracle RMSD ≤ x (the method's oracle RMSDs sorted, plotted as i / n × 100).",
        "Each method uses its own pair count n (shown in the legend) as the denominator, so curves with different n stay comparable as percentages.",
    ])
    fig.text(0.5, 0.015, footnote, ha="center", va="bottom", fontsize=7.5,
             color="0.30", linespacing=1.35)
    # Paired cross-tool test on the per-complex oracle RMSD (Friedman + Kendall W,
    # Wilcoxon pairwise). Anchored lower-right, the empty corner of a RMSD CDF.
    try:
        _annotate_oracle_rmsd_stats(ax, stats, 0.985, 0.04, "right", "bottom")
    except Exception as exc:  # never let a stats failure break the figure
        print(f"  WARNING: oracle-RMSD CDF stats annotation skipped ({exc})")
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(out, dpi=160); plt.close(fig)


def plot_oracle_rmsd_box(df: pd.DataFrame, out: Path, stats: dict | None = None) -> None:
    """Oracle RMSD boxplot — best pose per pair, all tools.

    Outliers (points beyond the whiskers, i.e. > Q3 + 1.5×IQR) are hidden and
    the y-axis is capped at the tallest whisker so a few grossly mis-docked
    complexes don't compress every box into an unreadable sliver.
    """
    oracle = _oracle_per_pair(df)
    methods = sorted(oracle["method"].unique())
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    data = [oracle[oracle["method"] == m]["rmsd"].dropna().values for m in methods]
    bp = ax.boxplot(data,
                    tick_labels=[f"{TOOL_LABEL.get(m, m)}\n(n={len(d)})"
                                 for m, d in zip(methods, data)],
                    showmeans=True, patch_artist=True, showfliers=False)
    for patch, m in zip(bp["boxes"], methods):
        patch.set_facecolor(TOOL_COLORS.get(m, "grey"))
        patch.set_alpha(0.6)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    ax.axhline(2.0, color="grey", ls="--", alpha=0.6, label="2 Å")
    ax.set_ylabel("Oracle RMSD (Å)")
    ax.set_title(_vt("Distribution of oracle (best) RMSD per method"))
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    # Cap the y-axis at the tallest whisker (Q3 + 1.5×IQR) across methods so the
    # hidden outliers — and any outlier-driven mean marker — can't blow up the
    # scale and flatten the boxes.
    whisker_tops = []
    for d in data:
        if len(d) == 0:
            continue
        q1, q3 = np.percentile(d, [25, 75])
        inside = d[d <= q3 + 1.5 * (q3 - q1)]
        whisker_tops.append(float(inside.max()) if len(inside) else float(q3))
    if whisker_tops:
        ax.set_ylim(0, max(whisker_tops) * 1.08)
    fig.text(0.5, 0.012,
             "Outliers beyond the whiskers (> Q3 + 1.5×IQR) are hidden and the y-axis "
             "capped at the tallest whisker so the boxes stay legible; n is each "
             "method's full pair count.",
             ha="center", va="bottom", fontsize=7.5, color="0.30")
    # Same paired cross-tool test as the CDF (Friedman + Kendall W, Wilcoxon
    # pairwise). Anchored upper-right, above the low-RMSD boxes.
    try:
        _annotate_oracle_rmsd_stats(ax, stats, 0.985, 0.97, "right", "top")
    except Exception as exc:
        print(f"  WARNING: oracle-RMSD boxplot stats annotation skipped ({exc})")
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(out, dpi=160); plt.close(fig)


def _dumbbell(ax, labels, left, right, colors, *,
              left_name: str, right_name: str, value_fmt: str = "{:.0f}%",
              gap_fmt: str | None = None, gap_color: str = "#b22222",
              missing_note: str | None = None, xmax: float | None = None,
              label_pos: str = "beside"):
    """Horizontal dumbbell ("connected-dot") chart — a cleaner replacement for
    paired bars when the message is the *gap* between two per-row values.

    One row per ``labels`` entry (first label drawn at the top). A grey connector
    joins the ``left`` value (filled marker — the ceiling/larger quantity) and the
    ``right`` value (open marker) for that row, both in the row's ``colors`` entry.
    Rows whose ``right`` value is NaN draw only the left marker plus
    ``missing_note`` (e.g. EquiBind has no ranking). With ``gap_fmt`` the
    |left − right| gap is annotated on each connector. Returns legend handles so
    the caller controls legend placement.

    ``label_pos`` controls where the value labels go: ``"beside"`` (default) puts
    them just outside each marker horizontally; ``"below"`` hangs each value under
    its own marker, growing outward so even near-coincident markers stay legible
    (use when the two markers can sit close together, e.g. fig 09e).
    """
    from matplotlib.lines import Line2D
    y = np.arange(len(labels))[::-1]            # first label on top
    finite = [v for v in list(left) + list(right) if not pd.isna(v)]
    span = xmax if xmax else (max(finite) if finite else 1.0)
    dx = 0.012 * span

    for yi, lv, rv in zip(y, left, right):
        if not (pd.isna(lv) or pd.isna(rv)):
            ax.plot([min(lv, rv), max(lv, rv)], [yi, yi], color="#b8b8b8",
                    lw=2.6, zorder=1, solid_capstyle="round")
    for yi, lv, c in zip(y, left, colors):
        if not pd.isna(lv):
            ax.scatter(lv, yi, s=120, color=c, edgecolor="black",
                       linewidth=0.8, zorder=3)
    for yi, rv, c in zip(y, right, colors):
        if not pd.isna(rv):
            ax.scatter(rv, yi, s=120, facecolor="white", edgecolor=c,
                       linewidth=2.2, zorder=3)

    below = label_pos == "below"
    for yi, lv, rv in zip(y, left, right):
        if not (pd.isna(lv) or pd.isna(rv)):
            lo, hi = (lv, rv) if lv <= rv else (rv, lv)
            if below:
                # Each value hangs below its own marker, anchored on the far side
                # so the two texts grow apart — legible even when lo≈hi.
                ax.text(lo, yi - 0.20, value_fmt.format(lo), ha="right", va="top",
                        fontsize=8, fontweight="bold")
                ax.text(hi, yi - 0.20, value_fmt.format(hi), ha="left", va="top",
                        fontsize=8, fontweight="bold")
            else:
                ax.text(lo - dx, yi, value_fmt.format(lo), ha="right", va="center",
                        fontsize=8, fontweight="bold")
                ax.text(hi + dx, yi, value_fmt.format(hi), ha="left", va="center",
                        fontsize=8, fontweight="bold")
            if gap_fmt:
                ax.text((lo + hi) / 2, yi + 0.24, gap_fmt.format(abs(lv - rv)),
                        ha="center", va="bottom", fontsize=7.5,
                        color=gap_color, fontweight="bold")
        elif not pd.isna(lv):
            if below:
                ax.text(lv, yi - 0.20, value_fmt.format(lv), ha="left", va="top",
                        fontsize=8, fontweight="bold")
            else:
                ax.text(lv + dx, yi, value_fmt.format(lv), ha="left", va="center",
                        fontsize=8, fontweight="bold")
            if missing_note:
                ax.text(lv + dx, yi + 0.24, missing_note, ha="left", va="bottom",
                        fontsize=7, color="grey", style="italic")

    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, span * 1.12)
    ax.set_ylim(-0.6, len(labels) - 0.4)
    ax.grid(axis="x", alpha=0.3)
    ax.set_axisbelow(True)
    return [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555",
               markeredgecolor="black", markersize=11, label=left_name),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="#555555", markeredgewidth=2, markersize=11,
               label=right_name),
    ]


def plot_oracle_vs_top1_success(oracle_sum: pd.DataFrame,
                                top1_sum: pd.DataFrame,
                                out: Path,
                                stats: dict | None = None) -> None:
    """Oracle vs top-1 RMSD ≤ 2 Å success rate, as a dumbbell chart.

    Ranking tools (Vina, DiffDock) get a filled (oracle) and an open (top-1)
    marker joined by a connector whose length is the ranking-quality loss.
    Non-ranking tools (EquiBind) show only the oracle marker + "no ranking".
    """
    all_methods = list(oracle_sum.index)
    oracle_col = "oracle_rmsd_le_2.0A_%"
    top1_col = "top1_rmsd_le_2.0A_%"

    oracle_vals = [
        float(oracle_sum.loc[m, oracle_col]) if oracle_col in oracle_sum.columns else float("nan")
        for m in all_methods
    ]
    top1_vals = [
        float(top1_sum.loc[m, top1_col])
        if m in top1_sum.index and top1_col in top1_sum.columns
        else float("nan")
        for m in all_methods
    ]

    labels = [f"{TOOL_LABEL.get(m, m)}\n{_n_annot(oracle_sum, m, poses=True)}"
              for m in all_methods]
    colors = [TOOL_COLORS.get(m, "grey") for m in all_methods]
    fig, ax = plt.subplots(figsize=(9, 0.8 * len(all_methods) + 2.4))
    handles = _dumbbell(ax, labels, oracle_vals, top1_vals, colors,
                        left_name="Oracle (best of all docked poses)",
                        right_name="Top-1 (tool's rank-1 pose)",
                        value_fmt="{:.1f}%", missing_note="no ranking")
    ax.set_xlabel("% of receptor-ligand pairs with RMSD ≤ 2 Å")
    ax.set_title(_vt("Oracle vs Top-1 success — ranking quality gap\n"
                     "(connector length = ranking loss; ● oracle, ○ top-1)"))
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    # Paired proportion tests: cross-tool oracle omnibus + within-tool
    # oracle-vs-top-1 McNemar, with 95% Wilson CIs on each marker.
    if stats is not None:
        try:
            _annotate_success_dumbbell(fig, ax, all_methods, oracle_vals,
                                       top1_vals, stats)
            fig.tight_layout(rect=(0, 0.055, 1, 1))
        except Exception as exc:
            print(f"  WARNING: fig 03 stats annotation skipped ({exc})")
            fig.tight_layout()
    else:
        fig.tight_layout()
    fig.savefig(out, dpi=160); plt.close(fig)


def plot_pb_valid_success_bars(oracle_sum: pd.DataFrame,
                               top1_sum: pd.DataFrame,
                               out: Path,
                               stats: dict | None = None) -> None:
    """Per-method success dumbbell for the strict '≤ 2 Å & PB-valid' category.

    One first-class chart for the combined success rate that fig 09 only shows
    as an overlaid subset. For every method we draw, as a connected-dot pair:
      * Oracle — best-of-N pose per complex (all tools, incl. EquiBind)
      * Top-1  — the tool's rank-1 pose (ranking tools only; others get "no ranking")
    A pair counts as a success only if it has a pose that is BOTH RMSD ≤ 2 Å AND
    passes every canonical PoseBusters check. The denominator is ALL receptor-
    ligand complexes (an accurate-but-invalid pose counts as a failure) — the
    same numbers as ``{oracle,top1}_pb_valid_and_rmsd2_%`` in the summary CSVs.
    """
    all_methods = list(oracle_sum.index)
    oracle_col = "oracle_pb_valid_and_rmsd2_%"
    top1_col = "top1_pb_valid_and_rmsd2_%"

    oracle_vals = [
        float(oracle_sum.loc[m, oracle_col]) if oracle_col in oracle_sum.columns else float("nan")
        for m in all_methods
    ]
    top1_vals = [
        float(top1_sum.loc[m, top1_col])
        if m in top1_sum.index and top1_col in top1_sum.columns
        else float("nan")
        for m in all_methods
    ]

    labels = [f"{TOOL_LABEL.get(m, m)}\n{_n_annot(oracle_sum, m, poses=True)}"
              for m in all_methods]
    colors = [TOOL_COLORS.get(m, "grey") for m in all_methods]
    fig, ax = plt.subplots(figsize=(9, 0.8 * len(all_methods) + 2.4))
    handles = _dumbbell(ax, labels, oracle_vals, top1_vals, colors,
                        left_name="Oracle (best of all docked poses)",
                        right_name="Top-1 (tool's rank-1 pose)",
                        value_fmt="{:.1f}%", missing_note="no ranking")
    ax.set_xlabel("% of complexes with a pose ≤ 2 Å & PB-valid")
    ax.set_title(_vt("Accurate AND PoseBuster valid poses\n"
                     "(RMSD ≤ 2 Å and passes every PoseBusters check — "
                     "% of all complexes; ● oracle, ○ top-1)"))
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    # Same paired proportion tests as fig 03 on the strict ≤ 2 Å & PB-valid success.
    if stats is not None:
        try:
            _annotate_success_dumbbell(fig, ax, all_methods, oracle_vals,
                                       top1_vals, stats)
            fig.tight_layout(rect=(0, 0.055, 1, 1))
        except Exception as exc:
            print(f"  WARNING: fig 10 stats annotation skipped ({exc})")
            fig.tight_layout()
    else:
        fig.tight_layout()
    fig.savefig(out, dpi=160); plt.close(fig)


def plot_vs_posebusters_paper(oracle_sum: pd.DataFrame,
                              top1_sum: pd.DataFrame,
                              out: Path,
                              csv_out: Path | None = None,
                              stats: dict | None = None) -> None:
    """Side-by-side comparison: THIS study vs the published PoseBusters paper.

    For each shared method we plot our PRIMARY-prediction success next to the
    paper's reported value (PoseBusters Benchmark set), for the two headline
    metrics: %RMSD ≤ 2 Å and %(RMSD ≤ 2 Å & PB-valid). "Primary prediction" = the
    tool's top-1 (rank-1) pose for ranking tools (AutoDock Vina, DiffDock) — the
    same quantity the paper reports — and the best pose for EquiBind, which has no
    native ranking (marked '*', so it is an upper bound vs the paper's single
    output). PB-valid uses the full canonical PoseBusters suite (PB_CRITICAL_CHECKS),
    matching the paper. Reference numbers: POSEBUSTERS_PAPER_BENCHMARK (editable).
    Also writes a tidy CSV (this-study / paper / delta per method × metric).
    """
    methods = list(oracle_sum.index)
    chosen = [m for m in ("autodock", "diffdock") if m in methods]
    eq = [m for m in methods if m.startswith("equibind")]
    if eq:  # one EquiBind row — its best variant by oracle valid-success rate
        vcol = "oracle_pb_valid_and_rmsd2_%"
        best_eq = max(eq, key=lambda m: (float(oracle_sum.loc[m, vcol])
                                         if vcol in oracle_sum.columns else -1.0))
        chosen.append(best_eq)
    if not chosen:
        return

    rows = []
    for m in chosen:
        ranking = m in RANKING_TOOLS
        if ranking and m in top1_sum.index and "top1_rmsd_le_2.0A_%" in top1_sum.columns:
            our_rmsd2 = float(top1_sum.loc[m, "top1_rmsd_le_2.0A_%"])
            our_valid = float(top1_sum.loc[m, "top1_pb_valid_and_rmsd2_%"])
            sel = "top-1"
        else:
            our_rmsd2 = float(oracle_sum.loc[m, "oracle_rmsd_le_2.0A_%"])
            our_valid = float(oracle_sum.loc[m, "oracle_pb_valid_and_rmsd2_%"])
            sel = "best pose*"
        ref = POSEBUSTERS_PAPER_BENCHMARK.get(_paper_method_name(m), {})
        rows.append(dict(method_key=m, paper_method=_paper_method_name(m), selection=sel,
                         our_rmsd2=our_rmsd2, our_valid=our_valid,
                         paper_rmsd2=ref.get("rmsd2"), paper_valid=ref.get("rmsd2_valid"),
                         paper_note=ref.get("note", "")))

    labels = [f"{TOOL_LABEL.get(r['method_key'], r['method_key'])}\n"
              f"({r['selection']}, "
              f"{_n_annot(top1_sum if r['selection'] == 'top-1' else oracle_sum, r['method_key'])})"
              for r in rows]
    x = np.arange(len(rows)); w = 0.38
    panels = [("RMSD ≤ 2 Å (accuracy)", "our_rmsd2", "paper_rmsd2"),
              ("RMSD ≤ 2 Å & PB-valid (headline)", "our_valid", "paper_valid")]
    # metric-key per panel, to look up the one-sample binomial-test payload.
    _panel_metric = {"our_rmsd2": "rmsd_le_2A", "our_valid": "rmsd_le_2A_and_pb_valid"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    for ax, (title, ourk, paperk) in zip(axes, panels):
        our_vals = [r[ourk] for r in rows]
        paper_vals = [r[paperk] for r in rows]
        colors = [TOOL_COLORS.get(r["method_key"], "grey") for r in rows]
        ax.bar(x - w / 2, our_vals, w, color=colors, edgecolor="black",
               label="This study")
        ax.bar(x + w / 2, [v if v is not None else 0.0 for v in paper_vals], w,
               color="#9e9e9e", edgecolor="black", hatch="//", alpha=0.7,
               label="PoseBusters paper")
        # 95% Wilson CI whiskers + one-sample binomial star (our rate vs the
        # paper's fixed published rate) on each "this study" bar.
        if stats is not None:
            try:
                mkey = _panel_metric.get(ourk)
                for xi, r, v in zip(x, rows, our_vals):
                    md = (stats.get("methods", {}).get(str(r["method_key"]), {})
                          .get("metrics", {}).get(mkey))
                    if not md or pd.isna(v):
                        continue
                    ax.errorbar(xi - w / 2, v,
                                yerr=[[max(0.0, v - md["wilson_lo_%"])],
                                      [max(0.0, md["wilson_hi_%"] - v)]],
                                fmt="none", ecolor="black", elinewidth=1.1,
                                capsize=3, zorder=4)
                    if md.get("star"):
                        ax.text(xi - w / 2, min(v + 7, 103),
                                f"{md['star']}\n{su.fmt_p(md['binom_p'])}",
                                ha="center", va="bottom", fontsize=6.8,
                                color="#8b0000", fontweight="bold")
            except Exception as exc:
                print(f"  WARNING: fig 11 stats annotation skipped ({exc})")
        for xi, v in zip(x, our_vals):
            if not pd.isna(v):
                ax.text(xi - w / 2, v + 1, f"{v:.0f}%", ha="center",
                        fontsize=8, fontweight="bold")
        for xi, v in zip(x, paper_vals):
            if v is None:
                ax.text(xi + w / 2, 2, "n/a", ha="center", fontsize=7,
                        color="grey", style="italic", rotation=90)
            else:
                ax.text(xi + w / 2, v + 1, f"{v:.0f}%", ha="center",
                        fontsize=8, color="#444")
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_title(title, fontweight="bold")
        ax.set_ylim(0, 105); ax.grid(axis="y", alpha=0.3)
    axes[0].set_ylabel("% of receptor-ligand complexes")
    axes[0].legend(loc="upper right", fontsize=9)
    _label_panels(axes)
    fig.suptitle(_vt("This study vs PoseBusters paper — Benchmark set, primary pose\n"
                     "(ranking tools: our top-1 pose; EquiBind*: best pose, "
                     "no native ranking — upper bound)"),
                 fontsize=12, fontweight="bold")
    if stats is not None:
        try:
            fig.text(0.5, 0.012,
                     "Stars: one-sample exact binomial test of our rate vs the "
                     "paper's FIXED published rate (whiskers = 95% Wilson CI on "
                     "our rate). The complex sets differ (paper = PoseBusters "
                     "Benchmark), so this is descriptive, not like-for-like.",
                     ha="center", va="bottom", fontsize=7.2, color="0.30",
                     linespacing=1.3)
            fig.tight_layout(rect=(0, 0.05, 1, 0.91))
        except Exception as exc:
            print(f"  WARNING: fig 11 footnote skipped ({exc})")
            fig.tight_layout(rect=(0, 0, 1, 0.91))
    else:
        fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(out, dpi=160); plt.close(fig)

    if csv_out is not None:
        recs = []
        for r in rows:
            for metric, ours, paper in (
                ("rmsd_le_2A", r["our_rmsd2"], r["paper_rmsd2"]),
                ("rmsd_le_2A_and_pb_valid", r["our_valid"], r["paper_valid"]),
            ):
                ok = ours is not None and not pd.isna(ours)
                recs.append({
                    "method": r["paper_method"], "our_method_key": r["method_key"],
                    "metric": metric, "our_selection": r["selection"],
                    "this_study_%": round(ours, PCT_DECIMALS) if ok else None,
                    "posebusters_paper_%": paper,
                    "delta_%": (round(ours - paper, PCT_DECIMALS) if (ok and paper is not None) else None),
                    "paper_note": r["paper_note"],
                })
        _write_csv(pd.DataFrame(recs), csv_out, index=False)
        print(f"  wrote paper-comparison table → {csv_out.name}")


def plot_rmsd2_vs_pbvalid_grouped(oracle_sum: pd.DataFrame,
                                  top1_sum: pd.DataFrame,
                                  out: Path) -> None:
    """Per-method comparison of the two success criteria, as dumbbells.

    The companion to fig 09 (which OVERLAYS them): here "RMSD ≤ 2 Å" (filled) and
    "RMSD ≤ 2 Å & PB-valid" (open) are the two ends of a connector whose length is
    the drop (Δ pp = the accuracy lost to the physical-validity requirement),
    annotated on each connector. Two panels: Oracle (all tools) and Top-1
    (ranking tools).
    """
    panels = [("Oracle (best of all docked poses)", oracle_sum,
               "oracle_rmsd_le_2.0A_%", "oracle_pb_valid_and_rmsd2_%"),
              ("Top-1 (tool's rank-1 pose)", top1_sum,
               "top1_rmsd_le_2.0A_%", "top1_pb_valid_and_rmsd2_%")]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True)
    handles = None
    for ax, (title, summ, acol, vcol) in zip(axes, panels):
        if summ is None or summ.empty or acol not in summ.columns:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, color="grey", style="italic")
            ax.set_title(title, fontweight="bold"); ax.set_yticks([])
            continue
        methods = list(summ.index)
        a = [float(summ.loc[m, acol]) for m in methods]
        v = [float(summ.loc[m, vcol]) if vcol in summ.columns else float("nan")
             for m in methods]
        colors = [TOOL_COLORS.get(m, "grey") for m in methods]
        labels = [f"{TOOL_LABEL.get(m, m)}\n{_n_annot(summ, m, poses=acol.startswith('oracle'))}"
                  for m in methods]
        handles = _dumbbell(ax, labels, a, v, colors,
                            left_name="RMSD ≤ 2 Å",
                            right_name="RMSD ≤ 2 Å & PB-valid",
                            value_fmt="{:.0f}%", gap_fmt="−{:.0f} pp", xmax=100)
        ax.set_title(title, fontweight="bold")
        ax.set_xlabel("% of receptor-ligand complexes")
    if handles is not None:
        axes[0].legend(handles=handles, loc="lower right", fontsize=9)
    _label_panels(axes)
    fig.suptitle(_vt("Effect of the PB-valid filter — RMSD ≤ 2 Å vs RMSD ≤ 2 Å & PB-valid\n"
                     "(Δ pp on each connector = accurate poses lost to PoseBuster validity)"),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.90)); fig.savefig(out, dpi=160); plt.close(fig)


def plot_pbvalid_filter_influence(oracle_sum: pd.DataFrame,
                                  out: Path,
                                  csv_out: Path | None = None) -> None:
    """How much does requiring PB-validity cut each method? (oracle, all tools)

    A dumbbell per method: the circle is % RMSD ≤ 2 Å (accuracy), the diamond is
    % RMSD ≤ 2 Å & PB-valid; the bar between them is the drop caused by the
    physical-validity filter. Each row is labelled with the percentage-point drop
    and the retention rate (valid / accurate) — the fraction of accurate poses
    that survive PB-validity. Classical docking barely moves; DL methods lose most
    of their accurate-but-physically-invalid poses. Also writes the numbers as CSV.
    """
    from matplotlib.lines import Line2D
    acol, vcol = "oracle_rmsd_le_2.0A_%", "oracle_pb_valid_and_rmsd2_%"
    if oracle_sum is None or oracle_sum.empty or acol not in oracle_sum.columns:
        return
    methods = sorted(oracle_sum.index, key=lambda m: float(oracle_sum.loc[m, acol]),
                     reverse=True)
    a = [float(oracle_sum.loc[m, acol]) for m in methods]
    v = [float(oracle_sum.loc[m, vcol]) if vcol in oracle_sum.columns else float("nan")
         for m in methods]
    y = list(range(len(methods)))[::-1]   # top row = best accuracy

    fig, ax = plt.subplots(figsize=(9.5, max(3.5, 0.7 * len(methods) + 2)))
    for yi, m, av, vv in zip(y, methods, a, v):
        col = TOOL_COLORS.get(m, "grey")
        if not pd.isna(vv):
            ax.plot([vv, av], [yi, yi], color=col, lw=5, alpha=0.45, zorder=1)
            ax.scatter([vv], [yi], color=col, s=95, marker="D",
                       edgecolor="black", zorder=2)
        ax.scatter([av], [yi], color=col, s=95, edgecolor="black", zorder=2)
        ax.text(av + 1.5, yi, f"{av:.0f}%", va="center", fontsize=8)
        if not pd.isna(vv):
            ax.text(vv - 1.5, yi, f"{vv:.0f}%", va="center", ha="right",
                    fontsize=8, fontweight="bold")
            ret = (vv / av * 100) if av > 0 else float("nan")
            lbl = f"−{av - vv:.0f} pp" + ("" if math.isnan(ret) else f" · {ret:.0f}% kept")
            ax.text((av + vv) / 2, yi + 0.20, lbl, va="bottom", ha="center",
                    fontsize=7, color="#b22222", fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{TOOL_LABEL.get(m, m)}  ({_n_annot(oracle_sum, m, poses=True)})"
                        for m in methods])
    ax.set_xlabel("% of receptor-ligand complexes (oracle / best pose)")
    ax.set_xlim(0, (max(a) if a else 100) + 14)
    # Headroom so each row's "Δ pp · % kept" label (drawn at yi + 0.20) stays
    # inside the axes and never collides with the title.
    ax.set_ylim(-0.6, (len(methods) - 1) + 0.75)
    ax.set_title(_vt("Influence of the PB-valid filter (oracle, all tools)\n"
                     "circle = RMSD ≤ 2 Å · diamond = also PB-valid · bar = accuracy lost"),
                 fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=10, label="RMSD ≤ 2 Å"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=10, label="RMSD ≤ 2 Å & PB-valid"),
    ], loc="lower right", fontsize=9)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)

    if csv_out is not None:
        recs = []
        for m, av, vv in zip(methods, a, v):
            ok = not pd.isna(vv)
            recs.append({
                "method": TOOL_LABEL.get(m, m), "method_key": m,
                "rmsd_le_2A_%": round(av, PCT_DECIMALS),
                "rmsd_le_2A_and_pb_valid_%": round(vv, PCT_DECIMALS) if ok else None,
                "drop_pp": round(av - vv, PCT_DECIMALS) if ok else None,
                "valid_retention_%": round(vv / av * 100, PCT_DECIMALS) if (ok and av > 0) else None,
            })
        _write_csv(pd.DataFrame(recs), csv_out, index=False)
        print(f"  wrote PB-valid influence table → {csv_out.name}")


def _group_guided_equibind(method: str) -> str:
    """Merge EquiBind's two guided pocket front-ends (fpocket / p2rank) into one 'guided'
    bucket, preserving the refine/clamp axes; 'unguided' and non-EquiBind keys pass through.
    Lets the pocket-guided EquiBind variants — which are handed the binding site, so warrant
    reading apart from the blind methods — be considered as one group."""
    m = str(method)
    if not m.startswith("equibind"):
        return m
    pocket, refine, clamp = _eq_tokens(m)
    if pocket in ("fpocket", "p2rank"):
        pocket = "guided"
    return "_".join(["equibind", *[t for t in (pocket, refine, clamp) if t]])


def within2_validity_comparison(df: pd.DataFrame, thr: float = 2.0,
                                group_guided_equibind: bool = False) -> pd.DataFrame:
    """Per-method PB-validity among near-native (RMSD ≤ *thr* Å) poses, computed two ways.

    Contrasts a tool's near-native validity rate over (A) EVERY generated pose within *thr* Å
    of the crystal, versus (B) only the ORACLE pick — the single lowest-RMSD pose per complex,
    which the oracle selects on RMSD alone and so may be physically invalid. (A) answers "when
    this tool lands a pose near-native, how often is it also PB-valid?"; (B) restricts that to
    the best-RMSD pose, i.e. how often the oracle's near-native pick is valid — identical to
    ``oracle_pb_valid_and_rmsd2_% / oracle_rmsd_le_2.0A_%`` from the oracle summary. A gap
    (A − B) means the closest pose is systematically more/less valid than a typical near-native
    one.

    With *group_guided_equibind*, EquiBind's fpocket and p2rank poses are pooled into a single
    'guided' bucket per refine level (they are handed the pocket, so warrant reading apart from
    the blind methods); the oracle then picks the best pose across both pocket front-ends. One
    row per method, carrying n_complexes, total_poses and the near-native / oracle counts. Empty
    frame if the needed columns are absent (e.g. a crystal-free set, which has no RMSD)."""
    need = {"method", "protein", "ligand", "rmsd", "pb_valid"}
    if not need <= set(df.columns):
        return pd.DataFrame()
    d = df.dropna(subset=["rmsd"]).copy()
    if group_guided_equibind:
        d["method"] = d["method"].map(_group_guided_equibind)
    # Idempotent bool coercion (safe whether pb_valid is already bool or a "True"/"False" string).
    d["pb_valid"] = d["pb_valid"].map(lambda x: str(x).strip().lower() in ("true", "1"))
    rows = []
    for m, sub in d.groupby("method"):
        w2 = sub[sub["rmsd"] <= thr]                                        # all near-native poses
        orc = sub.loc[sub.groupby(["protein", "ligand"])["rmsd"].idxmin()]  # oracle pick per pair
        ow2 = orc[orc["rmsd"] <= thr]
        a_n, a_v = len(w2), int(w2["pb_valid"].sum())
        b_n, b_v = len(ow2), int(ow2["pb_valid"].sum())
        rows.append({
            "method": str(m),
            "n_complexes": int(sub.groupby(["protein", "ligand"]).ngroups),
            "total_poses": int(len(sub)),
            "near_native_poses": a_n, "near_native_valid": a_v,
            "near_native_pb_valid_%": round(100 * a_v / a_n, PCT_DECIMALS) if a_n else float("nan"),
            "oracle_near_native": b_n, "oracle_near_native_valid": b_v,
            "oracle_near_native_pb_valid_%": round(100 * b_v / b_n, PCT_DECIMALS) if b_n else float("nan"),
        })
    out = pd.DataFrame(rows).set_index("method")
    out["delta_pp"] = (out["near_native_pb_valid_%"]
                       - out["oracle_near_native_pb_valid_%"]).round(PCT_DECIMALS)
    return out


def _stats_within2_gap(df: pd.DataFrame, thr: float = 2.0) -> dict | None:
    """Within-method paired test for the near-native validity dumbbell (fig 19):
    is the oracle (min-RMSD) near-native pose systematically more/less PB-valid
    than a typical near-native pose of the SAME complex — i.e. is the diamond's
    offset from the circle real, or sampling noise?

    Per method, over (protein, ligand) complexes that have ≥1 near-native pose
    (RMSD ≤ *thr* Å; equivalently, whose oracle pick is itself near-native):
      pool_rate    = mean pb_valid over that complex's near-native poses (the circle)
      oracle_valid = pb_valid of the min-RMSD pose, 0/1                  (the diamond)
    Wilcoxon signed-rank on (oracle_valid − pool_rate) tests the gap; a complex
    with a single near-native pose contributes a zero difference and is dropped by
    the signed-rank (it carries no information about the gap). rank_biserial > 0 ⇒
    the closest pose is MORE valid than a typical near-native pose; < 0 ⇒ LESS
    (the "leftward diamond"). The pooled circle rate carries a **cluster (complex)
    bootstrap** CI — near-native poses are pooled and pseudoreplicated — while the
    diamond is one pose per complex, so a plain Wilson CI is exact. Mirrors the
    aggregation in :func:`within2_validity_comparison` so the numbers line up with
    the plotted percentages.
    """
    need = {"method", "protein", "ligand", "rmsd", "pb_valid"}
    if not need <= set(df.columns):
        return None
    d = df.dropna(subset=["rmsd"]).copy()
    d["pb_valid"] = d["pb_valid"].map(lambda x: str(x).strip().lower() in ("true", "1"))
    payload = {
        "figure": "19_within2_validity_dumbbell.png",
        "unit": "(protein, ligand) complex; near-native poses pooled with a cluster (complex) bootstrap",
        "test": (f"within-method paired Wilcoxon signed-rank of oracle-pose vs near-native-"
                 f"pool PB-validity (RMSD ≤ {thr:g} Å); cluster (complex) bootstrap 95% CI on "
                 "the pool rate, Wilson 95% CI on the oracle-pose rate"),
        "thr_A": float(thr),
        "per_method": {},
    }
    for m, sub in d.groupby("method"):
        w2 = sub[sub["rmsd"] <= thr]                                         # near-native poses
        if w2.empty:
            continue
        pair_id = (w2["protein"].astype(str) + "|" + w2["ligand"].astype(str)).to_numpy()
        circle_rate, c_lo, c_hi = su.cluster_bootstrap_ci(
            w2["pb_valid"].to_numpy(float), pair_id, statistic=np.mean)
        orc = sub.loc[sub.groupby(["protein", "ligand"])["rmsd"].idxmin()]   # oracle pick / complex
        ow2 = orc[orc["rmsd"] <= thr]
        b_n, b_v = int(len(ow2)), int(ow2["pb_valid"].sum())
        d_lo, d_hi = su.wilson_ci(b_v, b_n)
        pool = w2.groupby(["protein", "ligand"])["pb_valid"].mean()
        ora = ow2.set_index(["protein", "ligand"])["pb_valid"].astype(float)
        joined = pd.concat([pool.rename("pool"), ora.rename("oracle")], axis=1).dropna()
        rec = {"circle_rate": float(circle_rate), "circle_ci": [float(c_lo), float(c_hi)],
               "near_native_poses": int(len(w2)),
               "diamond_rate": (b_v / b_n if b_n else float("nan")),
               "diamond_ci": [float(d_lo), float(d_hi)], "oracle_near_native_n": b_n,
               "n_complexes_paired": int(len(joined))}
        if len(joined) >= _MIN_UNITS_STATS:
            rb, p, npair = su.wilcoxon_rankbiserial(
                joined["oracle"].to_numpy(float), joined["pool"].to_numpy(float))
            rec.update(gap_test="wilcoxon_signed_rank", rank_biserial=rb, p=p,
                       n_nonzero_pairs=int(npair), star=su.p_stars(p),
                       median_gap_pp=float(100 * (joined["oracle"] - joined["pool"]).median()))
        else:
            rec["status"] = "n too small — exploratory"
        payload["per_method"][str(m)] = rec
    return payload if payload["per_method"] else None


def plot_within2_validity_dumbbell(comp: pd.DataFrame, out: Path, thr: float = 2.0,
                                   stats: dict | None = None) -> None:
    """Dumbbell per method: PB-validity among near-native (RMSD ≤ *thr* Å) poses, computed over
    ALL generated poses (circle) vs the ORACLE pick only (diamond).

    The circle is the fraction of every within-*thr* Å pose that is PB-valid; the diamond is the
    same fraction restricted to the single lowest-RMSD pose per complex (what the oracle picks,
    on RMSD alone). The bar between them is the gap: when they coincide the oracle's near-native
    pick is as valid as a typical near-native pose; a leftward diamond means the closest pose is
    LESS valid than the pool. Complements fig 13 (which stays entirely at the oracle level).
    Marks are coloured by tool; the two measures are told apart by shape (+ legend), so identity
    never rests on colour alone. Rows are grouped pocket-free tools → blind (unguided) EquiBind →
    pocket-GUIDED EquiBind (fpocket / p2rank), the last kept as separate rows but shaded and
    bracketed together, since being handed the binding site warrants reading them apart. Skips
    methods with no near-native pose."""
    from matplotlib.lines import Line2D
    acol, bcol = "near_native_pb_valid_%", "oracle_near_native_pb_valid_%"
    if comp is None or comp.empty or acol not in comp.columns:
        return

    def _grp(m):                                        # 0 tools · 1 blind EquiBind · 2 guided
        m = str(m)
        if not m.startswith("equibind"):
            return 0
        return 1 if _eq_tokens(m)[0] == "unguided" else 2
    methods = sorted((m for m in comp.index if not pd.isna(comp.loc[m, acol])),
                     key=lambda m: (_grp(m), -float(comp.loc[m, acol])))
    if not methods:
        return
    y = list(range(len(methods)))[::-1]                 # first method = top row
    grps = [_grp(m) for m in methods]
    fig, ax = plt.subplots(figsize=(11.5, max(3.5, 0.62 * len(methods) + 2)))
    has_stats = bool((stats or {}).get("per_method"))
    star_x = 109.0                                       # right-hand significance column

    # Shade + bracket the pocket-guided EquiBind block; faint dividers between the three groups.
    guided_y = [yi for yi, g in zip(y, grps) if g == 2]
    if guided_y:
        lo, hi = min(guided_y) - 0.48, max(guided_y) + 0.48
        ax.axhspan(lo, hi, color="#5aae61", alpha=0.09, zorder=0)
        xb = 2.5
        ax.plot([xb, xb], [lo + 0.18, hi - 0.18], color="#2f6b34", lw=1.8, zorder=2)
        for yy in (lo + 0.18, hi - 0.18):
            ax.plot([xb, xb + 2.2], [yy, yy], color="#2f6b34", lw=1.8, zorder=2)
        ax.text(xb + 4.5, (lo + hi) / 2, "guided EquiBind\n(fpocket / p2rank)\n— pocket given",
                va="center", ha="left", fontsize=8.5, color="#2f6b34", fontweight="bold")
    for i in range(1, len(methods)):
        if grps[i] != grps[i - 1]:
            ax.axhline((y[i] + y[i - 1]) / 2, color="#cccccc", lw=0.8, zorder=0)

    for yi, m in zip(y, methods):
        col = TOOL_COLORS.get(m, "grey")
        av = float(comp.loc[m, acol])
        bv = float(comp.loc[m, bcol]) if not pd.isna(comp.loc[m, bcol]) else float("nan")
        if not pd.isna(bv):
            ax.plot([min(av, bv), max(av, bv)], [yi, yi], color=col, lw=5, alpha=0.45, zorder=1)
            ax.scatter([bv], [yi], color=col, s=110, marker="D", edgecolor="black", zorder=3)
            ax.text(bv, yi - 0.24, f"{bv:.0f}%", va="top", ha="center",
                    fontsize=8, fontweight="bold")
        ax.scatter([av], [yi], color=col, s=110, edgecolor="black", zorder=3)
        ax.text(av, yi + 0.24, f"{av:.0f}%", va="bottom", ha="center", fontsize=8)
        # 95% CI whiskers (circle: cluster/complex bootstrap · diamond: Wilson) + the
        # oracle-vs-pool gap significance star (paired Wilcoxon), from _stats_within2_gap.
        _st = (stats or {}).get("per_method", {}).get(str(m))
        if _st:
            cc = _st.get("circle_ci") or []
            if len(cc) == 2 and cc[0] == cc[0]:
                ax.errorbar(av, yi + 0.13, xerr=[[max(0.0, av - cc[0] * 100)],
                            [max(0.0, cc[1] * 100 - av)]], fmt="none", ecolor=col,
                            elinewidth=1.1, capsize=2.4, alpha=0.85, zorder=2)
            dc = _st.get("diamond_ci") or []
            if not pd.isna(bv) and len(dc) == 2 and dc[0] == dc[0]:
                ax.errorbar(bv, yi - 0.13, xerr=[[max(0.0, bv - dc[0] * 100)],
                            [max(0.0, dc[1] * 100 - bv)]], fmt="none", ecolor=col,
                            elinewidth=1.1, capsize=2.4, alpha=0.85, zorder=2)
            star = _st.get("star")
            if star:
                sig = star not in ("ns", "")
                ax.text(star_x, yi, star, ha="center", va="center",
                        fontsize=10 if sig else 8, fontweight="bold" if sig else "normal",
                        color="0.15" if sig else "0.55", zorder=4)
            elif _st.get("status"):
                ax.text(star_x, yi, "—", ha="center", va="center", fontsize=8,
                        color="0.6", zorder=4)
    if has_stats:
        # header sits just above the top spine (axes-fraction y), clear of the frame
        ax.text(star_x, 1.012, "oracle\nvs pool", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7.2, color="0.3",
                fontweight="bold", clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{TOOL_LABEL.get(m, m)}\n{int(comp.loc[m, 'n_complexes'])} complexes · "
         f"{int(comp.loc[m, 'total_poses']):,} poses · "
         f"{int(comp.loc[m, 'near_native_poses'])} ≤ {thr:g} Å"
         for m in methods], fontsize=8)
    ax.set_xlabel(f"PB-valid share of near-native poses (RMSD ≤ {thr:g} Å) (%)")
    ax.set_xlim(0, 116 if has_stats else 104)
    ax.set_ylim(-0.7, (len(methods) - 1) + 0.7)
    ax.set_title(_vt(f"PB-validity of near-native poses: all generated vs. oracle pick\n"
                     f"circle = all poses ≤ {thr:g} Å RMSD  ·  diamond = oracle (min-RMSD) pick"),
                 fontsize=12, fontweight="bold")
    ax.grid(axis="x", alpha=0.3)
    ax.legend(handles=[
        Line2D([0], [0], marker="o", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=10, label=f"all generated poses ≤ {thr:g} Å"),
        Line2D([0], [0], marker="D", color="w", markerfacecolor="grey",
               markeredgecolor="black", markersize=10, label="oracle pick (min-RMSD)"),
    ], loc="lower right", fontsize=9)
    if has_stats:
        fig.text(0.01, 0.008, "whiskers = 95% CI (circle: cluster/complex bootstrap · "
                 "diamond: Wilson).  ★ column = paired Wilcoxon signed-rank, oracle "
                 "(min-RMSD) pose vs near-native-pool PB-validity, per complex "
                 "(*** p<.001 · ** p<.01 · * p<.05 · ns; — = n too small).",
                 fontsize=6.6, color="0.4")
    fig.tight_layout(rect=(0, 0.03, 1, 1) if has_stats else None)
    fig.savefig(out, dpi=160); plt.close(fig)


def _draw_accuracy_validity_panel(ax, title: str, summ: pd.DataFrame,
                                  rmsd_col: str, valid_col: str, *,
                                  poses: bool, label_map: "dict | None" = None) -> None:
    """Draw one PoseBusters-paper Fig. 1-style accuracy-vs-validity panel on *ax*.

    Per method, one bar split into two overlaid parts:
      * the full (light) bar    — % of complexes with a pose RMSD ≤ 2 Å
      * the solid (dark) bar     — the subset that is ALSO PB-valid
                                  (RMSD ≤ 2 Å **and** passes every PoseBusters check)
    The PB-valid set is a subset of the RMSD-≤-2-Å set, so the dark bar is always
    contained in the light one; the visible light "cap" above it is the fraction of
    accurate-but-physically-invalid predictions — the paper's central point. Values
    come straight from the aggregation tables (``rmsd_col`` / ``valid_col``); set
    ``poses`` to also show each method's oracle pose cap in the tick labels.
    """
    if summ is None or summ.empty or rmsd_col not in summ.columns:
        ax.text(0.5, 0.5, "no data", ha="center", va="center",
                transform=ax.transAxes, color="grey", style="italic")
        ax.set_title(title, fontweight="bold")
        ax.set_xticks([])
        return

    methods = list(summ.index)
    rmsd_vals = [float(summ.loc[m, rmsd_col]) for m in methods]
    valid_vals = [float(summ.loc[m, valid_col])
                  if valid_col in summ.columns else float("nan")
                  for m in methods]
    x = np.arange(len(methods))
    colors = [TOOL_COLORS.get(m, "#888888") for m in methods]

    # Light "RMSD ≤ 2 Å" bar (accuracy) …
    ax.bar(x, rmsd_vals, width=0.62, color=colors, alpha=0.32,
           edgecolor=colors, linewidth=1.4, zorder=2)
    # … with the solid "RMSD ≤ 2 Å & PB-valid" subset overlaid on top.
    ax.bar(x, valid_vals, width=0.62, color=colors, edgecolor="white",
           linewidth=0.8, zorder=3)

    for xi, rv, vv in zip(x, rmsd_vals, valid_vals):
        if not math.isnan(rv):
            ax.text(xi, rv + 1.5, f"{rv:.0f}%", ha="center", va="bottom",
                    fontsize=10, fontweight="bold", zorder=4)
        if not math.isnan(vv) and vv > 4:
            ax.text(xi, vv / 2, f"{vv:.0f}%", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold", zorder=4)

    def _lbl(m):
        return (label_map.get(m) if label_map and m in label_map
                else TOOL_LABEL.get(m, m))
    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{_lbl(m)}\n{_n_annot(summ, m, poses=poses)}" for m in methods],
        rotation=25, ha="right", rotation_mode="anchor")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, 105)
    ax.grid(axis="y", alpha=0.3)


def _accuracy_validity_legend(fig) -> None:
    """The shared method-agnostic legend (grey swatches) for the two bar layers."""
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#888888", alpha=0.32, edgecolor="#888888",
              label="RMSD ≤ 2 Å (accurate)"),
        Patch(facecolor="#555555", edgecolor="white",
              label="RMSD ≤ 2 Å & PB-valid (accurate + PoseBuster valid)"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.02), fontsize=10)


def plot_accuracy_validity_top1(top1_sum: pd.DataFrame, out: Path) -> None:
    """Fig 09a — accuracy-vs-validity bars for the TOP-1 (rank-1) pose.

    Ranking tools only (AutoDock Vina, DiffDock): the light bar is the % of
    complexes whose rank-1 pose is RMSD ≤ 2 Å, the dark subset the fraction that
    is also PB-valid. Single panel (split out of the former two-panel fig 09).
    """
    fig, ax = plt.subplots(figsize=(7.0, 5.8))
    _draw_accuracy_validity_panel(
        ax, "Top-1 (tool's rank-1 pose)", top1_sum,
        "top1_rmsd_le_2.0A_%", "top1_pb_valid_and_rmsd2_%", poses=False)
    ax.set_ylabel("% of receptor-ligand complexes")
    _accuracy_validity_legend(fig)
    fig.suptitle("Docking accuracy vs. PoseBuster validity — top-1 pose\n"
                 "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.09)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _oracle_variant_specs(df_full: pd.DataFrame, forced_dd: str | None = None):
    """[(method_key, label, kind)] for the 09b variants — each ML tool's RAW variant
    beside its best-by-PB-valid&≤2Å variant. ``kind`` ∈ {'native','generation'} picks
    how the rank-1 representative is taken (native rank vs first generated pose).
    ``forced_dd`` pins the DiffDock 'best' variant (e.g. 'diffdock_gnina') instead of the
    oracle-ranked one, matching --collapse-diffdock-variant."""
    oa = aggregate_oracle(df_full)
    col = _variant_rank_col(oa)
    present = set(df_full["method"].astype(str))
    if col is None or oa.empty:
        return []

    def _best(prefix):
        ks = [m for m in oa.index if str(m).startswith(prefix)]
        sc = oa.loc[ks, col].astype(float).dropna()
        return str(sc.idxmax()) if not sc.empty else None

    def _opt(m):
        for t in ("smina", "gnina"):
            if t in str(m):
                return f"{t}-opt"
        return "raw"

    def _raw_of(best):
        if not best:
            return None
        if str(best).startswith("diffdock"):
            return "diffdock"
        pocket, _refine, _clamp = _eq_tokens(str(best))
        return f"equibind_{pocket}_raw" if pocket else "equibind_unguided_raw"

    dd_best, eb_best = _best("diffdock"), _best("equibind")
    if forced_dd and forced_dd in present:
        dd_best = forced_dd
    eb_raw = _raw_of(eb_best)
    specs = []

    def _add(key, label, kind):
        if key and key in present and key not in [s[0] for s in specs]:
            specs.append((key, label, kind))

    # AutoDock gets the same raw-beside-best treatment as the ML tools. The literal
    # "autodock" key used to be the only AutoDock entry here, so once the Meeko arms
    # were excluded the family vanished from 09b/09e entirely.
    _add(_DOMINANT_AUTODOCK_RAW, "AutoDock Vina", "native")
    _add(_DOMINANT_AUTODOCK_OPT, "AutoDock* (gnina-opt)", "native")
    _add("diffdock", "DiffDock (raw)", "native")
    if dd_best != "diffdock":
        _add(dd_best, f"DiffDock* ({_opt(dd_best)})", "native")
    _add(eb_raw, "EquiBind (raw)", "generation")
    if eb_best != eb_raw:
        _add(eb_best, f"EquiBind* ({_opt(eb_best)})", "generation")
    return specs


def plot_accuracy_validity_oracle(oracle_sum: pd.DataFrame, out: Path,
                                  df_full: "pd.DataFrame | None" = None,
                                  thr: float = 2.0,
                                  forced_dd: str | None = None) -> None:
    """Fig 09b — accuracy-vs-validity bars per variant: the RANK-1 pose (EquiBind has
    no ranking, so its FIRST generated pose) beside the ORACLE (best-of-N) pose, with
    each ML tool's raw variant next to its best variant. Light bar = RMSD ≤ 2 Å, dark
    subset = also PB-valid. Falls back to a single oracle bar per (collapsed) tool when
    ``df_full`` is absent or has no crystal RMSD.
    """
    specs = _oracle_variant_specs(df_full, forced_dd=forced_dd) if df_full is not None else []

    rows, color_key = [], {}
    for key, label, kind in specs:
        sub = df_full[df_full["method"].astype(str) == key].copy()
        if kind == "native":
            sub["_rk"] = pd.to_numeric(sub["rank"], errors="coerce")
            r1 = (sub[sub["_rk"] == 1].sort_values("pose_name")
                     .groupby(["protein", "ligand"]).first().reset_index())
            r1_sel = "rank-1"
        else:
            sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
            r1 = (sub.sort_values("_gi", kind="mergesort")
                     .groupby(["protein", "ligand"]).first().reset_index())
            r1_sel = "first pose"
        orc = _oracle_per_pair(sub)
        color_key[label] = key
        for sel, reps in [(r1_sel, r1), ("oracle", orc)]:
            d = _av_from_reps(reps, thr)
            rows.append({"variant": label, "selection": sel,
                         "acc": d.get("near_pct", float("nan")),
                         "valid": d.get("valid_pct", float("nan")), "n": d.get("n", 0)})
    comp = pd.DataFrame(rows)

    # ── fallback: no full frame / no crystal RMSD → the original single-bar panel ──
    if comp.empty or comp["acc"].isna().all():
        fig, ax = plt.subplots(figsize=(8.0, 5.8))
        _draw_accuracy_validity_panel(
            ax, "Oracle (best of all docked poses)", oracle_sum,
            "oracle_rmsd_le_2.0A_%", "oracle_pb_valid_and_rmsd2_%", poses=True)
        ax.set_ylabel("% of receptor-ligand complexes")
        _accuracy_validity_legend(fig)
        fig.suptitle("Docking accuracy vs. PoseBuster validity — oracle (best-of-N) pose\n"
                     "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.09)
        fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
        return

    # ── grouped bars, three spacing levels ──
    #   intra          — the two bars WITHIN a variant (rank-1 vs oracle), tightest
    #   inter_variant  — between a tool's raw and opt variants (kept close, same tool)
    #   inter_tool     — between different tools (widest)
    variants = [s[1] for s in specs]
    bar_w, intra, inter_variant, inter_tool = 0.56, 0.60, 0.72, 1.20
    positions, spans, tool_spans, plotted, x, prev_tool = [], {}, {}, [], 0.0, None
    for label in variants:
        tool = _fam_key(str(color_key[label]))
        g = comp[comp["variant"] == label]
        g = pd.concat([g[g["selection"] != "oracle"], g[g["selection"] == "oracle"]])  # rank-1 then oracle
        if prev_tool is not None:
            x += inter_tool if tool != prev_tool else inter_variant
        xs = []
        for k, (_, r) in enumerate(g.iterrows()):
            if k > 0:
                x += intra
            positions.append(x); xs.append(x); plotted.append(r)
        spans[label] = xs
        tool_spans.setdefault(tool, []).extend(xs)
        prev_tool = tool

    # Colour every bar by its TOOL's raw/base variant, not the per-variant tint, so
    # a tool's optimized bars (DiffDock*, EquiBind*) share the exact dark PB-valid
    # colour of its raw bars — the family is read off position + the tool bracket,
    # not hue. specs list the raw variant of each family first, so its colour is the
    # family base. (Local to 09b; 09c/09d already colour by family.)
    fam_base_color: dict[str, str] = {}
    for key, _label, _kind in specs:
        fam = _fam_key(str(key))
        fam_base_color.setdefault(fam, TOOL_COLORS.get(key, "#888888"))

    fig, ax = plt.subplots(figsize=(max(9.0, 1.45 * (max(positions) + 1.4)), 6.3))
    for xi, r in zip(positions, plotted):
        c = fam_base_color.get(_fam_key(str(color_key[r["variant"]])), "#888888")
        acc, val = float(r["acc"]), float(r["valid"])
        ax.bar(xi, acc, width=bar_w, color=c, alpha=0.32, edgecolor=c, linewidth=1.4, zorder=2)
        ax.bar(xi, val, width=bar_w, color=c, edgecolor="white", linewidth=0.8, zorder=3)
        if not math.isnan(acc):
            ax.text(xi, acc + 1.5, f"{acc:.1f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold", zorder=4)
        if not math.isnan(val) and val > 4:
            ax.text(xi, val / 2, f"{val:.1f}%", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold", zorder=4)
    ax.set_xticks(positions)
    ax.set_xticklabels([f"{r['selection']}\nn={int(r['n'])}" for r in plotted], fontsize=8)
    ax.tick_params(axis="x", length=0, pad=6)
    trans = ax.get_xaxis_transform()
    for label, xs in spans.items():
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - bar_w / 2, hi + bar_w / 2], [-0.16, -0.16], transform=trans,
                color="0.45", lw=1.0, clip_on=False, zorder=1)
        ax.text(xc, -0.18, label, transform=trans, ha="center", va="top",
                fontsize=9.5, fontweight="bold", clip_on=False)
    ax.set_ylabel("% of receptor-ligand complexes")
    ax.set_ylim(0, 105)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.grid(axis="y", alpha=0.3)
    _accuracy_validity_legend(fig)
    fig.suptitle("Docking accuracy vs. PoseBuster validity — rank-1 vs. oracle pose\n"
                 "raw vs. best variant · PoseBusters Benchmark  "
                 "(EquiBind rank-1 = first generated pose)",
                 fontsize=13, fontweight="bold", y=1.10)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── 09e: strict PB-valid success dumbbell over the 09b variants (fig-10 style) ──
#
# Takes 09b's per-variant data (each ML tool's raw vs best-optimised variant, plus
# AutoDock) and renders the strict "RMSD ≤ 2 Å & PB-valid" success as a fig-10
# dumbbell: ● oracle (best of all docked poses) vs ○ the pose you'd actually PICK
# by ranking. Unlike 09b — which uses EquiBind's FIRST generated pose as a stand-in
# rank-1 — EquiBind here is ranked by GNINA AFFINITY, the only real ranking it has.
# gnina scores exist only for the gnina-optimised variant (raw/smina EquiBind carry
# no gnina score), so raw EquiBind stays oracle-only, marked "no ranking". Same
# paired-proportion tests as figs 03/10: per-variant Wilson CIs, within-variant
# oracle-vs-rank-1 exact McNemar, and a cross-variant Cochran's Q + pairwise
# McNemar (Holm), rendered by the shared _annotate_success_dumbbell.


def _native_rank1_per_pair(sub: pd.DataFrame) -> pd.DataFrame:
    """One row per (protein, ligand): the tool's confidence rank-1 pose."""
    sub = sub.copy()
    sub["_rk"] = pd.to_numeric(sub["rank"], errors="coerce")
    return (sub[sub["_rk"] == 1].sort_values("pose_name")
               .groupby(["protein", "ligand"]).first().reset_index())


def _gnina_rank1_per_pair(sub: pd.DataFrame) -> pd.DataFrame:
    """One row per (protein, ligand): the best gnina-affinity pose (most negative;
    generation index breaks ties, NaN affinities sort last). This is the ranking
    gnina supplies to EquiBind, whose native output has no pose ranking."""
    sub = sub.copy()
    sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
    sub["_aff"] = pd.to_numeric(sub["gnina_affinity"], errors="coerce")
    return (sub.sort_values(["_aff", "_gi"], na_position="last", kind="mergesort")
               .groupby(["protein", "ligand"]).first().reset_index())


def _stats_pb_valid_variants(df_full: "pd.DataFrame | None", thr: float = 2.0,
                             forced_dd: str | None = None,
                             valid: bool = True) -> dict | None:
    """Per-09b-variant paired-proportion stats + plotted values for the fig-09e
    dumbbell. The payload is shaped exactly like ``_stats_success_paired`` so
    ``_annotate_success_dumbbell`` renders it directly. ``valid`` toggles the metric
    between '≤thr Å & PB-valid' (fig-10 style; default) and plain '≤thr Å'.
    """
    if df_full is None:
        return None
    specs = _oracle_variant_specs(df_full, forced_dd=forced_dd)
    if not specs:
        return None

    def _succ(reps: pd.DataFrame) -> pd.Series:
        s = reps["rmsd"] <= thr
        if valid:
            s = s & _to_bool(reps["pb_valid"])
        return s.astype(float)

    def _by_pair(reps: pd.DataFrame) -> pd.Series:
        r = reps.copy()
        r["_s"] = _succ(r)
        return r.set_index(["protein", "ligand"])["_s"]

    methods, per, orc_cols, wt = [], {}, {}, []
    for key, label, _kind in specs:
        sub = df_full[df_full["method"].astype(str) == key].dropna(subset=["rmsd"]).copy()
        if sub.empty:
            continue
        fam = _fam_key(str(key))
        has_gnina = ("gnina_affinity" in sub.columns
                     and pd.to_numeric(sub["gnina_affinity"], errors="coerce").notna().any())
        if fam in RANKING_TOOLS:
            r1, rank_note = _native_rank1_per_pair(sub), "confidence rank-1"
        elif has_gnina:
            r1, rank_note = _gnina_rank1_per_pair(sub), "gnina-affinity rank-1"
        else:
            r1, rank_note = None, "no ranking (raw)"

        o = _by_pair(_oracle_per_pair(sub))
        orc_cols[label] = o
        k, n = int(o.sum()), int(o.size)
        lo, hi = su.wilson_ci(k, n)
        d = {"method_key": key, "rank_note": rank_note, "n": n,
             "oracle_k": k, "oracle_n": n,
             "oracle_rate": (k / n if n else float("nan")),
             "oracle_lo": lo, "oracle_hi": hi}
        if r1 is not None and not r1.empty:
            t = _by_pair(r1)
            k2, n2 = int(t.sum()), int(t.size)
            lo2, hi2 = su.wilson_ci(k2, n2)
            d.update(top1_k=k2, top1_n=n2,
                     top1_rate=(k2 / n2 if n2 else float("nan")),
                     top1_lo=lo2, top1_hi=hi2)
            pair = pd.concat([o.rename("o"), t.rename("t")], axis=1).dropna()
            if len(pair) >= _MIN_UNITS_STATS:
                n10, n01, p = su.mcnemar_exact(pair["o"].to_numpy(), pair["t"].to_numpy())
                wt.append({"method": label, "n": int(len(pair)), "oracle_wins": n10,
                           "top1_wins": n01, "p": float(p), "star": su.p_stars(p)})
            else:
                wt.append({"method": label, "n": int(len(pair)),
                           "status": "n too small — exploratory"})
        methods.append(label)
        per[label] = d

    if not methods:
        return None

    payload = {
        "figure": "09e_pb_valid_variants_success_dumbbell.png",
        "metric": ("RMSD ≤ 2 Å AND PB-valid" if valid else "RMSD ≤ 2 Å"),
        "unit": "(protein, ligand) complex — per-complex success boolean",
        "ranking": ("Vina/DiffDock: native confidence rank-1; EquiBind: gnina-affinity "
                    "rank-1 (gnina-optimised variant only; raw/smina EquiBind carry no "
                    "gnina score → oracle-only)"),
        "test": ("paired proportions across variants (Cochran's Q + exact McNemar, "
                 "Holm) on the oracle success; within-variant oracle-vs-rank-1 exact "
                 "McNemar; 95% Wilson CIs"),
        "methods": methods,
        "per_method": per,
        "within_tool_oracle_vs_top1": wt,
    }
    orc = pd.DataFrame(orc_cols)
    complete = orc.dropna()
    if len(complete) >= _MIN_UNITS_STATS and orc.shape[1] >= 2:
        payload["cross_tool_oracle"] = su.paired_proportions(
            {m: complete[m].to_numpy(float) for m in methods})
    else:
        payload["cross_tool_oracle"] = {"status": "n too small — exploratory",
                                        "n": int(len(complete))}
    return payload


def plot_pb_valid_variants_dumbbell(df_full, out: Path, thr: float = 2.0,
                                    forced_dd: str | None = None,
                                    valid: bool = True,
                                    stats: dict | None = None) -> None:
    """Fig 09e — the strict '≤2 Å & PB-valid' success from 09b's per-variant data,
    drawn as a fig-10 dumbbell (● oracle vs ○ selected rank-1) with EquiBind ranked
    by gnina affinity. Paired-proportion tests (Wilson CIs, within-variant McNemar,
    cross-variant Cochran's Q) are annotated when available.
    """
    payload = stats if stats is not None else _stats_pb_valid_variants(
        df_full, thr, forced_dd, valid)
    if not payload or not payload.get("methods"):
        return
    labels = payload["methods"]
    per = payload["per_method"]
    oracle_vals = [per[l]["oracle_rate"] * 100 for l in labels]
    top1_vals = [(per[l]["top1_rate"] * 100 if per[l].get("top1_rate") is not None
                  else float("nan")) for l in labels]

    # Colour each row by its tool's raw/base hue (matches the 09b convention), so
    # a tool's raw and optimised rows read as the same colour, distinguished by label.
    fam_base: dict[str, str] = {}
    for l in labels:
        fam = _fam_key(str(per[l]["method_key"]))
        fam_base.setdefault(fam, TOOL_COLORS.get(per[l]["method_key"], "#888888"))
    colors = [fam_base[_fam_key(str(per[l]["method_key"]))] for l in labels]
    ylabels = [f"{l}\n{per[l]['rank_note']} · n={per[l]['n']}" for l in labels]

    metric_txt = "≤ 2 Å & PB-valid" if valid else "≤ 2 Å"
    fig, ax = plt.subplots(figsize=(10.6, 0.9 * len(labels) + 2.8))
    handles = _dumbbell(ax, ylabels, oracle_vals, top1_vals, colors,
                        left_name="Oracle (best of all docked poses)",
                        right_name="Selected rank-1 (confidence / gnina affinity)",
                        value_fmt="{:.1f}%", missing_note="no ranking",
                        label_pos="below")
    ax.set_xlabel(f"% of receptor-ligand complexes with a pose {metric_txt}")
    ax.set_title(_vt(f"Accurate & PoseBuster-valid poses per variant  ({metric_txt})\n"
                     "raw vs. gnina-optimised · EquiBind ranked by gnina affinity"),
                 fontsize=11.5)
    ax.legend(handles=handles, loc="lower right", fontsize=9)
    if payload.get("per_method"):
        try:
            _annotate_success_dumbbell(fig, ax, labels, oracle_vals, top1_vals,
                                       payload, footnote_pairwise=False)
            fig.tight_layout(rect=(0, 0.055, 1, 1))
        except Exception as exc:
            print(f"  WARNING: fig 09e stats annotation skipped ({exc})")
            fig.tight_layout()
    else:
        fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ── Raw vs best-variant — the effect of post-hoc optimization (fig 09c) ──────

_EQ_POSE_IDX_RE = re.compile(r"_(\d{2,3})__ref", re.IGNORECASE)


def _equibind_pose_index(name) -> int:
    """Generation index parsed from an EquiBind pose_name (…_003__refRAW.sdf → 3).

    EquiBind emits its poses in generation order as ``<pocket>_<NNN>__ref<TAG>.sdf``;
    the lowest index is the FIRST generated pose. Returns a large sentinel when the
    index can't be parsed so such poses sort last.
    """
    m = _EQ_POSE_IDX_RE.search(str(name))
    return int(m.group(1)) if m else 10 ** 6


def aggregate_optimization_raw_vs_best(df: pd.DataFrame, thr: float = 2.0) -> pd.DataFrame:
    """Per tool, contrast the RAW representative against post-hoc optimizer arms.

    ONE representative pose is chosen per complex, then aggregated over complexes:

      * AutoDock raw   — rank-1 pose under Vina's original ranking.
      * AutoDock rescore/refinement — rank-1 pose under each optimizer arm's
                                      own re-ranking.
      * DiffDock raw   — rank-1 pose of the un-optimised run (``diffdock``).
      * DiffDock best  — rank-1 pose of the smina-optimised run (``diffdock_smina``);
                         same DiffDock confidence ranking, refined geometry.
      * EquiBind raw   — the FIRST generated pose (lowest generation index) of the
                         blind/unguided run (``equibind_unguided_raw``): EquiBind has
                         no native ranking, so the first pose is what you'd take.
      * EquiBind best  — the gnina-optimised blind run (``equibind_unguided_gnina``),
                         RANKED by gnina affinity (best / most-negative = the pick);
                         gnina supplies both refined geometry AND a ranking.

    Requires the FULL (un-collapsed) per-pose frame — every DiffDock/EquiBind
    variant must still be present under its own method key. Missing variants are
    silently skipped. Returns one row per bar: accuracy (% RMSD ≤ ``thr``), validity
    (% RMSD ≤ ``thr`` AND PB-valid), n_complexes and the median RMSD of the picks.
    """
    def _rank1(sub: pd.DataFrame) -> pd.DataFrame:
        sub = sub.copy()
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        r1 = sub[sub["rank"] == 1].sort_values("pose_name")
        return r1.groupby(["protein", "ligand"]).first().reset_index()

    def _first_pose(sub: pd.DataFrame) -> pd.DataFrame:
        sub = sub.copy()
        sub["_pidx"] = sub["pose_name"].map(_equibind_pose_index)
        return (sub.sort_values("_pidx", kind="mergesort")
                   .groupby(["protein", "ligand"]).first().reset_index())

    def _gnina_ranked(sub: pd.DataFrame) -> pd.DataFrame:
        sub = sub.copy()
        sub["_aff"] = pd.to_numeric(sub["gnina_affinity"], errors="coerce")
        # Best affinity = most negative; NaN sorts last so an all-NaN complex
        # falls back to its first row.
        return (sub.sort_values("_aff", kind="mergesort", na_position="last")
                   .groupby(["protein", "ligand"]).first().reset_index())

    # (tool, role_label, is_best, method_key, picker)
    specs = [
        # AutoDock rows follow the dominant arm (ADFRsuite ligands, exhaustiveness
        # 128). Its ladder point has no smina or CNN-refine sibling, so the family
        # contributes raw + gnina only; the Meeko arms that used to sit here are
        # excluded from this report (see method_filter.json).
        ("autodock", "raw\n(Vina rank-1)",       False, _DOMINANT_AUTODOCK_RAW, _rank1),
        ("autodock", "CNN rescore\n(opt rank-1)", True, _DOMINANT_AUTODOCK_OPT, _rank1),
        ("diffdock", "raw\n(rank-1)",             False, "diffdock",                _rank1),
        ("diffdock", "gnina-opt\n(rank-1)",       True,  "diffdock_gnina",          _rank1),
        ("equibind", "raw\n(first pose)",         False, "equibind_unguided_raw",   _first_pose),
        ("equibind", "gnina-opt\n(gnina-ranked)", True,  "equibind_unguided_gnina", _gnina_ranked),
    ]
    acc_col = f"rmsd_le_{thr:g}A_%"
    val_col = f"pb_valid_and_rmsd{thr:g}_%"
    rows = []
    for tool, role_label, is_best, mkey, picker in specs:
        sub = df[df["method"].astype(str) == mkey]
        if sub.empty:
            continue
        reps = picker(sub).dropna(subset=["rmsd"])
        n = reps.groupby(["protein", "ligand"]).ngroups if len(reps) else 0
        if n == 0:
            continue
        near = reps["rmsd"] <= thr
        valid = near & _to_bool(reps["pb_valid"])
        rows.append({
            "tool": tool, "role": role_label.replace("\n", " "),
            "role_label": role_label, "is_best": is_best, "method_key": mkey,
            "n_complexes": n,
            acc_col: round(100 * float(near.mean()), 2),
            val_col: round(100 * float(valid.mean()), 2),
            "median_rmsd": round(float(reps["rmsd"].median()), 3),
        })
    return pd.DataFrame(rows)


def plot_optimization_raw_vs_best(agg: pd.DataFrame, out: Path,
                                  thr: float = 2.0) -> None:
    """Fig 09c — grouped accuracy-vs-validity bars contrasting each tool's RAW
    representative pose with its available post-hoc optimizer variants.

    Same light(accuracy) + dark(PB-valid subset) encoding as fig 09a/09b; bars are
    grouped by tool (colour = tool identity), raw on the left, followed by the
    optimizer arms, with a Δ callout of each validity gain. AutoDock uses each
    arm's explicit optimized ranking for its representative. See
    ``aggregate_optimization_raw_vs_best`` for the exact per-tool pose selection.
    """
    if agg is None or agg.empty:
        return
    acc_col = f"rmsd_le_{thr:g}A_%"
    val_col = f"pb_valid_and_rmsd{thr:g}_%"
    tool_order = ["autodock", "diffdock", "equibind"]
    color_key = {"autodock": "autodock", "diffdock": "diffdock",
                 "equibind": "equibind_unguided"}
    group_name = {"autodock": "AutoDock Vina", "diffdock": "DiffDock",
                  "equibind": "EquiBind (blind / unguided)"}

    agg = agg[agg["tool"].isin(tool_order)].copy()
    if agg.empty:
        return
    agg["_torder"] = agg["tool"].map({t: i for i, t in enumerate(tool_order)})
    agg = agg.sort_values(["_torder", "is_best"], kind="mergesort").reset_index(drop=True)

    # Lay out bar x-positions: within-group spacing + a wider gap between tools.
    intra, inter, bar_w = 0.80, 0.70, 0.60
    positions, spans, prev = [], {}, None
    x = 0.0
    for tool in agg["tool"]:
        if prev is not None and tool != prev:
            x += inter
        positions.append(x)
        spans.setdefault(tool, []).append(x)
        x += intra
        prev = tool
    agg["_x"] = positions

    fig, ax = plt.subplots(figsize=(11.0, 6.4))
    for _, r in agg.iterrows():
        c = TOOL_COLORS.get(color_key[r["tool"]], "#888888")
        xi, acc, val = r["_x"], float(r[acc_col]), float(r[val_col])
        ax.bar(xi, acc, width=bar_w, color=c, alpha=0.32, edgecolor=c,
               linewidth=1.4, zorder=2)
        ax.bar(xi, val, width=bar_w, color=c, edgecolor="white",
               linewidth=0.8, zorder=3)
        ax.text(xi, acc + 1.5, f"{acc:.0f}%", ha="center", va="bottom",
                fontsize=10, fontweight="bold", zorder=4)
        if val > 4:
            ax.text(xi, val / 2, f"{val:.0f}%", ha="center", va="center",
                    fontsize=9, color="white", fontweight="bold", zorder=4)

    # Raw → optimizer validity-gain callout above every optimized arm.
    for tool in spans:
        sub = agg[agg["tool"] == tool]
        raw, optimized = sub[~sub["is_best"]], sub[sub["is_best"]]
        if len(raw):
            raw_valid = float(raw[val_col].iloc[0])
            for _, best in optimized.iterrows():
                dv = float(best[val_col]) - raw_valid
                xb, yb = float(best["_x"]), float(best[acc_col])
                sign = "+" if dv >= 0 else ""
                ax.annotate(f"{sign}{dv:.0f} pp\nvalid vs raw", xy=(xb, yb),
                            xytext=(0, 15), textcoords="offset points",
                            ha="center", va="bottom", fontsize=8.5,
                            color="#1f7a1f" if dv >= 0 else "#a11d1d",
                            fontweight="bold")

    # Per-bar tick labels (role + n) and a second tier of centred tool names.
    ax.set_xticks(agg["_x"].tolist())
    ax.set_xticklabels([f"{r['role_label']}\nn={int(r['n_complexes'])}"
                        for _, r in agg.iterrows()], fontsize=8.5)
    ax.tick_params(axis="x", length=0, pad=6)
    trans = ax.get_xaxis_transform()
    for tool, xs in spans.items():
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - bar_w / 2, hi + bar_w / 2], [-0.155, -0.155],
                transform=trans, color="0.45", lw=1.0, clip_on=False, zorder=1)
        ax.text(xc, -0.175, group_name.get(tool, str(tool)),
                transform=trans, ha="center",
                va="top", fontsize=10.5, fontweight="bold", clip_on=False)

    ax.set_ylabel("% of receptor-ligand complexes")
    ax.set_ylim(0, 105)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.grid(axis="y", alpha=0.3)
    _accuracy_validity_legend(fig)
    fig.suptitle("Effect of post-hoc optimization — raw vs. optimizer variants\n"
                 "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.09)
    fig.text(0.5, -0.09,
             "Representative pose per complex — Vina raw, GNINA CNN-rescore, and "
             "GNINA CNN-refinement are shown independently using their rank-1 pose. "
             "DiffDock: rank-1 of the raw vs. smina-optimized run (same ranking). "
             "EquiBind (blind/unguided run, no native ranking): first generated pose "
             "vs. gnina-optimized, ranked by gnina affinity.",
             ha="center", va="top", fontsize=8, color="0.35", wrap=True)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── 09f: PoseBusters-valid YIELD of the produced poses — raw vs optimized ──
# Unlike 09c/09e (one representative pose per complex), this looks at EVERY produced
# pose: per complex, the share of a variant's poses that pass PoseBusters, then the
# distribution of that per-complex yield across complexes as a box/whisker. smina and
# gnina refine the SAME poses raw produced (post-hoc), so the denominator is fixed and
# the numerator — how many survive PoseBusters — is exactly what optimization moves.

# (tool, role, is_opt, method_key) — raw first, then the optimizer variants, so each
# tool group reads raw → smina-opt → gnina-opt left to right.
_PBVALID_YIELD_SPECS = [
    # AutoDock = the dominant arm (ADFRsuite ligands, exhaustiveness 128). No smina
    # or CNN-refine sibling exists at that ladder point, so the family is raw+gnina.
    ("autodock", "raw",       False, _DOMINANT_AUTODOCK_RAW),
    ("autodock", "gnina-opt", True,  _DOMINANT_AUTODOCK_OPT),
    ("diffdock", "raw",       False, "diffdock"),
    ("diffdock", "smina-opt", True,  "diffdock_smina"),
    ("diffdock", "gnina-opt", True,  "diffdock_gnina"),
    ("equibind", "raw",       False, "equibind_unguided_raw"),
    ("equibind", "smina-opt", True,  "equibind_unguided_smina"),
    ("equibind", "gnina-opt", True,  "equibind_unguided_gnina"),
]

# The raw variant each optimized variant of a tool is measured against (the paired
# baseline for the "did optimization help" contrasts).
_PBVALID_YIELD_RAW = {
    # AutoDock raw = the dominant arm's OWN raw counterpart (same ladder point),
    # so the raw->gnina gain is measured within one search effort.
    "autodock": _DOMINANT_AUTODOCK_RAW, "diffdock": "diffdock",
    "equibind": "equibind_unguided_raw",
}

# Figure caption / reading notes for 09f. These used to sit under the plot; per request
# they now live only in the companion report (09f_pbvalid_yield_report.txt) so the
# figure stays clean. Kept ASCII-friendly for plain-text viewers.
_PBVALID_YIELD_FIGURE_NOTES = (
    "Box/whisker of the PoseBusters-valid YIELD of the poses each variant PRODUCED, "
    "raw beside its post-hoc-optimized variants, grouped by tool.\n"
    "\n"
    "  - Each box = the distribution ACROSS COMPLEXES of that variant's per-complex "
    "yield: the % of the poses it produced for a complex that pass ALL PoseBusters "
    "tests. The heavy line is the median (also printed above each box); the open "
    "diamond is the mean.\n"
    "  - Raw / reference boxes are hatched (before optimization); post-hoc-optimized "
    "(smina / gnina) boxes are solid.\n"
    "  - smina and gnina REFINE THE SAME POSES the raw run produced (post-hoc), so raw "
    "and its optimized variants share the denominator (the pose count) and differ only "
    "in how many survive PoseBusters -- the raw->opt upward shift is exactly the "
    "validity that optimization recovers.\n"
    "  - AutoDock Vina has no ML-geometry optimization step, so it appears once as the "
    "physics reference (raw == final).\n"
    "  - Significance brackets: paired Wilcoxon signed-rank on the per-complex yield "
    "(raw vs each optimized variant), Holm-corrected across ALL raw->opt contrasts; the "
    "bracket label is the Hodges-Lehmann median gain in percentage points (pp) plus "
    "its significance stars.\n"
    "  - Tick labels: 'pool' = validity pooled over ALL produced poses of that variant "
    "(sum valid / sum produced, across every complex); n = number of complexes.\n"
    "  - Outliers beyond the whiskers (> Q3 + 1.5*IQR) are hidden so the boxes stay "
    "legible; the box, whiskers and printed median/mean use the full per-complex data."
)


def _drop_bare_diffdock_dupe(sub: pd.DataFrame, method_key: str) -> pd.DataFrame:
    """Drop raw DiffDock's duplicate top pose when pooling ALL produced poses.

    Raw DiffDock writes its rank-1 pose twice — a bare ``rankN.sdf`` beside the
    byte-identical ``rankN_confidence-*.sdf`` — so any denominator that pools every
    produced pose of the raw ``diffdock`` key (the pose-production cascade, the 09f
    produced-pose yield) would otherwise count the top pose twice and bias both the
    totals and the pooled percentages. The optimised variants (smina/gnina) are
    confidence-named only, so they never carry the bare copy; this is a no-op for
    every other method. Mirrors the guard in :func:`_positions_native`.
    """
    if method_key != "diffdock":
        return sub
    bare = sub["pose_name"].astype(str).str.contains(
        r"rank\d+\.sdf$", regex=True, case=False, na=False)
    return sub[~bare]


def pbvalid_yield_per_complex(df_full: pd.DataFrame,
                              specs=_PBVALID_YIELD_SPECS) -> pd.DataFrame:
    """Per (variant, complex): what fraction of the poses a variant PRODUCED are
    PoseBusters-valid.

    One row per (method_key, protein, ligand): ``n_poses`` produced, ``n_valid`` of
    them, and ``pct_valid`` = 100·n_valid/n_poses. ALL produced (scored) poses count
    — the denominator is what the variant generated for that complex, so a grossly
    distorted pose still counts against the yield (it just isn't valid). The optimizer
    variants (smina/gnina) refine the SAME poses raw produced, so raw and its optimized
    variants share the denominator and differ only in how many survive PoseBusters.

    Requires the FULL (un-collapsed) per-pose frame so every DiffDock/EquiBind variant
    is still present under its own method key; missing variants are skipped.
    """
    cols = ["tool", "role", "is_opt", "method_key", "protein", "ligand",
            "n_poses", "n_valid", "pct_valid"]
    frames = []
    for tool, role, is_opt, mkey in specs:
        sub = df_full[df_full["method"].astype(str) == mkey]
        sub = _drop_bare_diffdock_dupe(sub, mkey)
        if sub.empty:
            continue
        g = sub.assign(_v=_to_bool(sub["pb_valid"]).astype(int)).groupby(
            ["protein", "ligand"], sort=False)
        agg = g["_v"].agg(n_poses="size", n_valid="sum").reset_index()
        agg["pct_valid"] = 100.0 * agg["n_valid"] / agg["n_poses"]
        agg["tool"], agg["role"], agg["is_opt"], agg["method_key"] = \
            tool, role, is_opt, mkey
        frames.append(agg[cols])
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)


def aggregate_pbvalid_yield_by_variant(df_full: pd.DataFrame,
                                       specs=_PBVALID_YIELD_SPECS) -> pd.DataFrame:
    """Per-variant summary of the produced-pose PB-valid yield (for the CSV / table).

    One row per variant: n_complexes, total poses produced & valid, the POOLED
    validity rate (over all poses), and the per-complex yield distribution
    (mean / median / IQR). ``median_gain_over_raw_pp`` is the median per-complex yield
    minus the raw variant of the same tool (0 for raw / the physics reference) — the
    headline "how much optimization added" in percentage points.
    """
    per = pbvalid_yield_per_complex(df_full, specs)
    if per.empty:
        return per
    order = {mkey: i for i, (_, _, _, mkey) in enumerate(specs)}
    rows = []
    for tool, role, is_opt, mkey in specs:
        sub = per[per["method_key"] == mkey]
        if sub.empty:
            continue
        pct = sub["pct_valid"].to_numpy(float)
        tot_poses = int(sub["n_poses"].sum())
        n_cplx = int(len(sub))
        n_cplx_valid = int((sub["n_valid"] >= 1).sum())   # complexes with >=1 valid pose
        rows.append({
            "tool": tool, "role": role, "is_opt": is_opt, "method_key": mkey,
            "n_complexes": n_cplx,
            # Coverage: for how many complexes does the variant produce AT LEAST ONE
            # PB-valid pose. Distinct from n_complexes (complexes with ANY pose) — a
            # variant can be present for all 303 yet yield zero valid poses on some.
            "n_complexes_with_valid": n_cplx_valid,
            "complexes_with_valid_%": round(100.0 * n_cplx_valid / max(1, n_cplx), PCT_DECIMALS),
            "n_poses": tot_poses, "n_valid": int(sub["n_valid"].sum()),
            "pooled_pb_valid_%": round(100.0 * sub["n_valid"].sum()
                                       / max(1, tot_poses), PCT_DECIMALS),
            "per_complex_mean_%": round(float(np.mean(pct)), PCT_DECIMALS),
            "per_complex_median_%": round(float(np.median(pct)), PCT_DECIMALS),
            "per_complex_q1_%": round(float(np.percentile(pct, 25)), PCT_DECIMALS),
            "per_complex_q3_%": round(float(np.percentile(pct, 75)), PCT_DECIMALS),
        })
    out = pd.DataFrame(rows)
    med = dict(zip(out["method_key"], out["per_complex_median_%"]))
    def _gain(r):
        raw_key = _PBVALID_YIELD_RAW.get(r["tool"])
        if not r["is_opt"]:
            return 0.0
        if raw_key is None or raw_key not in med:
            # No baseline resolves for this tool — report "not computed", never 0.0.
            # A hard 0.0 here is indistinguishable from a genuine no-change result and
            # published exactly that for autodock_vinardo_gnina_refinement, whose real
            # gain is -6.7 pp. The MGLTools sweep has no single baseline by design: the
            # exh-32 optimiser box is measured against exh 32 and the exh-64 one against
            # exh 64, so a per-tool map cannot express it.
            return float("nan")
        return round(r["per_complex_median_%"] - med[raw_key], PCT_DECIMALS)
    out["median_gain_over_raw_pp"] = out.apply(_gain, axis=1)
    out["_order"] = out["method_key"].map(order)
    return out.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def _stats_pbvalid_yield(per: pd.DataFrame, specs=_PBVALID_YIELD_SPECS) -> dict | None:
    """Paired per-complex test that optimization raises the produced-pose PB-valid
    yield, run separately for each ML tool (DiffDock, EquiBind).

    For a tool, the per-complex yields of its raw / smina-opt / gnina-opt variants
    are a paired (same-complex) triple: Friedman + Kendall's W give the omnibus that
    the three differ, and one-sided-in-interpretation Wilcoxon signed-rank
    (rank-biserial effect size) gives the raw→smina and raw→gnina contrasts that are
    the actual "did optimization help?" question, with a Hodges–Lehmann median gain
    and bootstrap CI. Restricted to complexes scored under all three variants.

    All raw→opt contrasts (across both tools) form ONE Holm family so the
    significance stars survive multiplicity correction. Returns None on any shortfall.
    """
    if per is None or per.empty:
        return None
    tool_variants = {}
    for tool, role, is_opt, mkey in specs:
        if tool in _PBVALID_YIELD_RAW:
            tool_variants.setdefault(tool, []).append((role, mkey))
    out: dict = {}
    contrasts: list[dict] = []                 # pooled Holm family across tools
    for tool, variants in tool_variants.items():
        frames = []
        for role, mkey in variants:
            sub = per[per["method_key"] == mkey]
            if not sub.empty:
                frames.append(
                    sub.set_index(["protein", "ligand"])["pct_valid"].rename(role))
        labels = [f.name for f in frames]
        if len(frames) < 2 or "raw" not in labels:
            continue
        wide = pd.concat(frames, axis=1, join="inner").dropna()
        if len(wide) < 3:
            continue
        pc = su.paired_continuous(
            {l: wide[l].to_numpy(float) for l in labels}, labels=labels)
        raw = wide["raw"].to_numpy(float)
        tool_recs = []
        for role in labels:
            if role == "raw":
                continue
            opt = wide[role].to_numpy(float)
            rb, p, npair = su.wilcoxon_rankbiserial(opt, raw)
            gain, lo, hi = su.median_diff_ci(opt, raw, paired=True)
            rec = {"tool": tool, "variant": role, "n": int(len(wide)),
                   "median_raw_%": round(float(np.median(raw)), PCT_DECIMALS),
                   "median_opt_%": round(float(np.median(opt)), PCT_DECIMALS),
                   "median_gain_pp": round(float(gain), PCT_DECIMALS),
                   "gain_ci": [round(float(lo), PCT_DECIMALS), round(float(hi), PCT_DECIMALS)],
                   "rank_biserial": round(float(rb), 3) if rb == rb else None,
                   "p_raw": float(p) if p == p else None}
            tool_recs.append(rec)
            contrasts.append(rec)
        out[tool] = {"omnibus": pc["omnibus"], "medians": pc["medians"],
                     "pairwise": pc["pairwise"], "raw_vs_opt": tool_recs,
                     "n_complexes": int(len(wide))}
    if not out:
        return None
    valid_p = [c for c in contrasts if c["p_raw"] is not None]
    if valid_p:
        for c, pa in zip(valid_p, su.holm([c["p_raw"] for c in valid_p])):
            c["p_holm"] = float(pa)
            c["star"] = su.p_stars(pa)
    return out


def _pbvalid_yield_contrast(stats: dict | None, tool: str, variant: str) -> dict | None:
    """Look up one raw→opt contrast record from the 09f stats payload."""
    if not stats or tool not in stats:
        return None
    for rec in stats[tool].get("raw_vs_opt", []):
        if rec["variant"] == variant:
            return rec
    return None


def _write_pbvalid_yield_report(summary: "pd.DataFrame | None", stats: dict | None,
                                out_path: Path) -> None:
    """Write the FULL 09f statistics + the moved-off-figure notes to a text file.

    Companion to 09f_pbvalid_yield_boxplot.png / pbvalid_yield_by_variant.csv: the
    figure caption/reading notes (moved off the plot), the per-variant produced-pose
    PB-valid yield table, and the complete paired per-complex statistics — Friedman
    omnibus (+ Kendall's W), the full Wilcoxon pairwise family, and the raw→opt
    contrasts with rank-biserial effect sizes, Hodges–Lehmann median gains, 95%
    bootstrap CIs and Holm-adjusted p-values.
    """
    def _p(p) -> str:
        if p is None or p != p:
            return "n/a"
        if p <= 0:
            return "<1e-300"
        return f"{p:.3e}" if p < 1e-3 else f"{p:.4f}"

    def _f(v, nd=2) -> str:
        return "n/a" if v is None or v != v else f"{float(v):.{nd}f}"

    def _sf(v, nd=3) -> str:                       # None/NaN-safe SIGNED float
        return "n/a" if v is None or v != v else f"{float(v):+.{nd}f}"

    W = 84
    L: list[str] = []
    L.append("=" * W)
    L.append("Search Effort, Post-Hoc Optimization and PoseBusters Validity")
    L.append("Figure 09f — PoseBusters-valid yield of the produced poses (raw vs. optimized)")
    L.append("=" * W)
    L.append("")
    L.append("Source : per_pose_metrics.csv (df_full — all raw + smina/gnina variants)")
    L.append("Tables : pbvalid_yield_by_variant.csv (per-variant),")
    L.append("         pbvalid_yield_per_complex.csv (box source)")
    L.append("Figure : 09f_pbvalid_yield_boxplot.png")
    L.append("")
    L.append("-" * W)
    L.append("FIGURE NOTES  (moved off the plot)")
    L.append("-" * W)
    L.append(_PBVALID_YIELD_FIGURE_NOTES)
    L.append("")

    if summary is not None and not summary.empty:
        L.append("-" * W)
        L.append("PER-VARIANT SUMMARY  (produced-pose PoseBusters-valid yield)")
        L.append("-" * W)
        # cplx>=1v / cover% answer "for how many complexes does the variant produce AT
        # LEAST ONE PB-valid pose" — distinct from n_cplx (complexes with ANY pose).
        has_cover = "n_complexes_with_valid" in summary.columns
        hdr = (f"{'tool':9s} {'role':10s} {'method_key':24s} {'n_cplx':>6s} "
               + (f"{'cplx>=1v':>8s} {'cover%':>6s} " if has_cover else "")
               + f"{'n_poses':>7s} {'n_valid':>7s} {'pool%':>6s} {'mean%':>6s} "
               f"{'med%':>6s} {'Q1%':>6s} {'Q3%':>6s} {'gain_vs_raw_pp':>14s}")
        L.append(hdr)
        L.append("-" * len(hdr))
        for _, r in summary.iterrows():
            cov = (f"{int(r['n_complexes_with_valid']):>8d} "
                   f"{r['complexes_with_valid_%']:>6.1f} " if has_cover else "")
            L.append(f"{str(r['tool']):9s} {str(r['role']):10s} "
                     f"{str(r['method_key']):24s} {int(r['n_complexes']):>6d} " + cov
                     + f"{int(r['n_poses']):>7d} {int(r['n_valid']):>7d} "
                     f"{r['pooled_pb_valid_%']:>6.1f} {r['per_complex_mean_%']:>6.1f} "
                     f"{r['per_complex_median_%']:>6.1f} {r['per_complex_q1_%']:>6.1f} "
                     f"{r['per_complex_q3_%']:>6.1f} {r['median_gain_over_raw_pp']:>+14.1f}")
        L.append("")
        if has_cover:
            L.append("  cplx>=1v = complexes with >=1 PB-valid pose (of n_cplx); cover% = its rate")
        L.append("  pool%   = validity pooled over ALL produced poses (sum valid / sum produced)")
        L.append("  mean/med/Q1/Q3% = distribution across complexes of the per-complex yield")
        L.append("  gain_vs_raw_pp  = per-complex MEDIAN yield minus the tool's raw variant (pp)")
        L.append("  NOTE: n_cplx counts complexes with ANY pose; a variant can be present for all")
        L.append("        of them yet still yield ZERO valid poses on some (see cplx>=1v / cover%).")
        L.append("")

    L.append("-" * W)
    L.append("STATISTICAL RESULTS")
    L.append("-" * W)
    L.append("Unit of analysis : one value per complex = % of a variant's produced poses that")
    L.append("                   are PoseBusters-valid (correlated poses aggregated per complex).")
    L.append("Design           : paired — the SAME complexes across a tool's raw / smina-opt /")
    L.append("                   gnina-opt variants (listwise-complete complexes only).")
    L.append("Omnibus          : Friedman rank test (+ Kendall's W concordance), per tool.")
    L.append("Contrasts        : Wilcoxon signed-rank (matched-pairs rank-biserial effect size);")
    L.append("                   raw→opt gain = Hodges–Lehmann median of the paired differences")
    L.append("                   with a 95% percentile-bootstrap CI.")
    L.append("Multiplicity     : Holm–Bonferroni across ALL raw→opt contrasts (both tools pooled).")
    L.append("Significance      : *** p<.001   ** p<.01   * p<.05   ns = not significant.")
    L.append("")

    if not stats:
        L.append("No paired statistics available — optimized variants absent, or fewer than 3")
        L.append("complexes scored under all of a tool's raw/smina/gnina variants. Re-run the")
        L.append("report with all variants present (a --best-*-only / --collapse-plots-only")
        L.append("flag, or --diffdock-variant all --split-equibind).")
    else:
        pretty = {"diffdock": "DiffDock", "equibind": "EquiBind (blind / unguided)"}
        for tool in ("diffdock", "equibind"):
            blk = stats.get(tool)
            if not blk:
                continue
            om = blk.get("omnibus", {})
            n = blk.get("n_complexes", om.get("n", 0))
            L.append(f"### {pretty.get(tool, tool)}   (n = {n} complexes)")
            L.append("")
            meds = blk.get("medians", {})
            if meds:
                L.append("  Per-complex median yield:  "
                         + ",  ".join(f"{k} = {_f(v, 1)}%" for k, v in meds.items()))
            L.append("  Omnibus (raw vs smina-opt vs gnina-opt):")
            L.append(f"    Friedman χ²({om.get('df', '?')}) = {_f(om.get('chi2'), 2)}, "
                     f"p = {_p(om.get('p'))}, Kendall's W = {_f(om.get('kendall_w'), 3)}")
            L.append("")
            pw = blk.get("pairwise", [])
            if pw:
                L.append("  Pairwise Wilcoxon signed-rank (Holm within this tool's 3-pair family;")
                L.append("  rank-biserial > 0 ⇒ the FIRST-named variant has the higher yield):")
                for pr in pw:
                    npair = pr.get("n")
                    n_str = f", n={int(npair)}" if npair is not None and npair == npair else ""
                    L.append(f"    {str(pr['a']):>9s} vs {str(pr['b']):<9s} : "
                             f"r_rb = {_sf(pr.get('rank_biserial'))}, "
                             f"p_raw = {_p(pr.get('p_raw'))}, "
                             f"p_holm = {_p(pr.get('p_holm'))}  {pr.get('star', '')}{n_str}")
                L.append("")
            rvo = blk.get("raw_vs_opt", [])
            if rvo:
                L.append("  Raw → optimized (headline; Holm across all raw→opt contrasts):")
                for rec in rvo:
                    ci = rec.get("gain_ci", [None, None]) or [None, None]
                    L.append(f"    raw → {str(rec['variant']):<9s} : "
                             f"median {_f(rec.get('median_raw_%'), 1)}% → "
                             f"{_f(rec.get('median_opt_%'), 1)}%  (n={rec.get('n', '?')})")
                    L.append(f"        Hodges–Lehmann gain = {_sf(rec.get('median_gain_pp'), 1)} pp "
                             f"(95% CI [{_f(ci[0], 1)}, {_f(ci[1], 1)}])")
                    L.append(f"        rank-biserial = {_sf(rec.get('rank_biserial'))}, "
                             f"p_raw = {_p(rec.get('p_raw'))}, "
                             f"p_holm = {_p(rec.get('p_holm'))}  {rec.get('star', '')}")
                L.append("")
        if summary is not None and not summary.empty:
            ad = summary[summary["method_key"] == "autodock"]
            if not ad.empty:
                a = ad.iloc[0]
                L.append("### AutoDock Vina — physics reference (not tested)")
                L.append("")
                L.append("  Single variant, no ML-geometry optimization step (raw ≡ final):")
                L.append(f"    pooled valid = {a['pooled_pb_valid_%']:.1f}%, per-complex median "
                         f"= {a['per_complex_median_%']:.1f}% (n = {int(a['n_complexes'])}).")
                if "n_complexes_with_valid" in a.index:
                    L.append(f"    produces >=1 PB-valid pose for "
                             f"{int(a['n_complexes_with_valid'])} / {int(a['n_complexes'])} "
                             f"complexes ({a['complexes_with_valid_%']:.1f}%) — the rest yield "
                             f"zero valid poses despite having poses.")
                L.append("")

    L.append("=" * W)
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _all_variant_yield_specs(df_full: pd.DataFrame):
    """09f specs for EVERY variant present in ``df_full`` — all EquiBind pocket ×
    refine [× clamp] variants + all DiffDock variants + AutoDock — independent of the
    best-variant / collapse flags.

    Ordered AutoDock, DiffDock (raw → smina-opt → gnina-opt), then EquiBind grouped by
    pocket (unguided, fpocket, p2rank, …) and within a pocket raw → smina-opt → gnina-opt
    (clamp last). ``role`` carries a compact multi-line per-variant label for the ticks.
    """
    present = set(df_full["method"].astype(str).unique())
    specs = []
    for mkey, role in (("autodock", "raw"), ("autodock_smina", "smina-opt"),
                       ("autodock_gnina", "gnina-opt"),
                       ("autodock_gnina_refinement", "gnina CNN-refine")):
        if mkey in present:
            specs.append(("autodock", role, mkey != "autodock", mkey))
    # MGLTools-ligand arm — ONE group carrying the whole exhaustiveness sweep, with an
    # independent box per setting. Two things this ordering fixes:
    #   * the unsuffixed key IS the exhaustiveness-32 point (that tree docks at 32), so it
    #     is labelled "exh 32", not "raw". Labelling it "raw" conflated the ligand-prep arm
    #     with the sweep's baseline and left 18/92 with no place in the figure.
    #   * giving each setting its own top-level group would draw four underlines and read
    #     as four separate engines on the tool axis. It is one engine sampled at four
    #     search efforts, so it gets one underline and four boxes.
    # Post-dock optimiser variants stay attached to the exhaustiveness point they derive
    # from; only exh 32 has any, since the 18/64/92/128 arms are raw-only by construction.
    # Effort axis first and monotone (18 → 32 → 64 → 92 → 128), THEN the post-dock optimiser
    # boxes. Interleaving an optimiser box between two exhaustiveness points breaks the
    # only axis the group's underline claims, and does so where the optimiser's effect
    # runs the opposite way to its neighbours (exh32 99.8 % → +gnina 99.0 %).
    for mkey, role in (("autodock_mgltools_exh18", "exh 18"),
                       ("autodock_mgltools", "exh 32"),
                       ("autodock_mgltools_exh64", "exh 64"),
                       ("autodock_mgltools_exh92", "exh 92"),
                       ("autodock_mgltools_exh128", "exh 128"),
                       ("autodock_mgltools_smina", "exh 32\nsmina-opt"),
                       ("autodock_mgltools_gnina", "exh 32\ngnina-opt"),
                       ("autodock_mgltools_gnina_refinement", "exh 32\ngnina CNN-refine"),
                       ("autodock_mgltools_exh64_gnina", "exh 64\ngnina-opt"),
                       ("autodock_mgltools_exh128_gnina", "exh 128\ngnina-opt")):
        if mkey in present:
            # is_opt marks POST-DOCK optimisation, so a higher-exhaustiveness sibling is
            # not "opt" — it is still a raw Vina search, just a longer one.
            specs.append(("autodock_mgltools", role,
                          mkey.endswith(("_smina", "_gnina", "_gnina_refinement")), mkey))
    for mkey, role in (("autodock_vinardo", "raw"),
                       ("autodock_vinardo_smina", "smina-opt"),
                       ("autodock_vinardo_gnina", "gnina-opt"),
                       ("autodock_vinardo_gnina_refinement", "gnina CNN-refine")):
        if mkey in present:
            specs.append(("autodock_vinardo", role,
                          mkey != "autodock_vinardo", mkey))
    for mkey, role in (("diffdock", "raw"), ("diffdock_smina", "smina-opt"),
                       ("diffdock_gnina", "gnina-opt")):
        if mkey in present:
            specs.append(("diffdock", role, mkey != "diffdock", mkey))
    # Single-variant engines (no raw→smina→gnina family): one spec each.
    for mkey in ("unidock", "unidock2"):
        if mkey in present:
            specs.append((mkey, "raw", False, mkey))
    eq = [m for m in present if m.startswith("equibind")]
    pocket_order = {"unguided": 0, "fpocket": 1, "p2rank": 2, "guided": 3}
    refine_order = {None: 0, "raw": 0, "smina": 1, "gnina": 2}
    def _order(m):
        p, r, c = _eq_tokens(m)
        return (pocket_order.get(p, 9), refine_order.get(r, 9), c or "")
    for m in sorted(eq, key=_order):
        p, r, c = _eq_tokens(m)
        bits = [p or "equibind",
                {"smina": "smina-opt", "gnina": "gnina-opt"}.get(r, "raw")]
        if c:
            bits.append("clamp on" if c == "clampON" else "clamp off")
        specs.append(("equibind", "\n".join(bits), r in ("smina", "gnina"), m))
    return specs


def _cascade_label(mkey: str) -> str:
    """Table label for the pose-validity cascade, in the paper's variant style
    ('DiffDock (raw)', 'DiffDock + smina', 'EquiBind unguided (raw)', …)."""
    if mkey == "autodock":
        return "AutoDock Vina (raw)"
    if mkey == "autodock_smina":
        return "AutoDock Vina + smina"
    if mkey == "autodock_gnina":
        return "AutoDock Vina + gnina"
    if mkey == "autodock_gnina_refinement":
        return "AutoDock Vina + GNINA CNN-refinement"
    # exh32 is the unsuffixed arm; state the setting so the sweep's baseline is explicit.
    if mkey == "autodock_mgltools":
        return "AutoDock Vina MGLTools-lig exh32 (raw)"
    if mkey == "autodock_mgltools_smina":
        return "AutoDock Vina MGLTools-lig exh32 + smina"
    if mkey == "autodock_mgltools_gnina":
        return "AutoDock Vina MGLTools-lig exh32 + gnina"
    if mkey == "autodock_mgltools_gnina_refinement":
        return "AutoDock Vina MGLTools-lig exh32 + GNINA CNN-refinement"
    if mkey == "autodock_mgltools_exh18":
        return "AutoDock Vina MGLTools-lig exh18 (raw)"
    if mkey == "autodock_mgltools_exh64":
        return "AutoDock Vina MGLTools-lig exh64 (raw)"
    if mkey == "autodock_mgltools_exh64_gnina":
        return "AutoDock Vina MGLTools-lig exh64 + gnina"
    if mkey == "autodock_mgltools_exh128_gnina":
        return "AutoDock Vina MGLTools-lig exh128 + gnina"
    if mkey == "autodock_mgltools_exh92":
        return "AutoDock Vina MGLTools-lig exh92 (raw)"
    if mkey == "autodock_mgltools_exh128":
        return "AutoDock Vina MGLTools-lig exh128 (raw)"
    if mkey == "autodock_vinardo":
        return "AutoDock Vinardo (raw)"
    if mkey == "autodock_vinardo_smina":
        return "AutoDock Vinardo + smina"
    if mkey == "autodock_vinardo_gnina":
        return "AutoDock Vinardo + gnina"
    if mkey == "autodock_vinardo_gnina_refinement":
        return "AutoDock Vinardo + GNINA CNN-refinement"
    if mkey == "diffdock":
        return "DiffDock (raw)"
    if mkey == "diffdock_smina":
        return "DiffDock + smina"
    if mkey == "diffdock_gnina":
        return "DiffDock + gnina"
    if mkey == "unidock":
        return "Uni-Dock (tiled)"
    if mkey == "unidock2":
        return "Uni-Dock2"
    if mkey.startswith("equibind"):
        pocket, refine, clamp = _eq_tokens(mkey)
        pk = {"unguided": "unguided", "fpocket": "fpocket",
              "p2rank": "P2Rank", "guided": "guided"}.get(pocket or "", pocket or "")
        head = f"EquiBind {pk}".rstrip()
        suffix = {"smina": " + smina", "gnina": " + gnina"}.get(refine or "", " (raw)")
        lab = f"{head}{suffix}"
        if clamp:
            lab += " [clamp on]" if clamp == "clampON" else " [clamp off]"
        return lab
    return mkey


def aggregate_pose_validity_cascade(df_full: pd.DataFrame, specs=None,
                                    kabsch_thr: float = FORM_OK_KABSCH_A) -> pd.DataFrame:
    """Pose-level validity/accuracy cascade per docking variant.

    For every variant present in ``df_full`` this pools ALL of the poses it produced
    (the produced-pose denominator, the same one the 09f yield uses) and walks the
    cascade::

        produced  →  PoseBusters-valid  →  RMSD ≤ 2 Å  →  Kabsch RMSD < thr Å
                  →  RMSD ≤ 2 Å AND PB-valid AND Kabsch RMSD < thr Å

    Two placement-free "form" filters live in the middle/end: ``Kabsch RMSD < thr``
    (default 1 Å, ``kabsch_thr``) is the best-fit / superposed heavy-atom RMSD — the
    ligand's INTERNAL conformation regardless of where it was placed (PoseBusters'
    ``pb_kabsch_rmsd``, our own Kabsch fit as fallback, stored in ``bestfit_rmsd``).
    The final column is the strict apex: a pose that is near-native (≤ 2 Å), physically
    valid (PB), AND has the correct form (Kabsch < thr).

    Each stage is reported BOTH as a pose count (how many poses survive the filter)
    and as a complex count (how many distinct (protein, ligand) complexes still have
    at least one surviving pose). The pose count is the pooled sampling volume; the
    complex count is how many of the targets are actually covered at that stage —
    they answer different questions (a variant can pour thousands of valid poses into
    a handful of easy complexes, or spread a few valid poses across many).

    Percentage columns:
      pb_valid_%          — PB-valid poses / produced poses.
      rmsd2_%             — RMSD ≤ 2 Å poses / produced poses.
      kabsch1_%           — Kabsch RMSD < thr poses / produced poses.
      triple_of_rmsd2_%   — (≤ 2 Å & PB-valid & Kabsch < thr) / (≤ 2 Å) poses.
      triple_of_all_%     — (≤ 2 Å & PB-valid & Kabsch < thr) / produced poses.

    The pure ``rmsd2_pbvalid_*`` columns (≤ 2 Å & PB-valid, the previous headline that
    matches the PoseBusters paper) are RETAINED in the CSV for reference — the fixed-
    width / markdown tables show the stricter triple as the final category instead.

    Uses the FULL (un-collapsed) per-pose frame so every DiffDock/EquiBind variant is
    present under its own method key. Crystal RMSD and Kabsch RMSD that are NaN never
    satisfy their thresholds. Validity is the canonical PB-valid flag (``_to_bool``).
    """
    if specs is None:
        specs = _all_variant_yield_specs(df_full)
    cols = ["tool", "role", "is_opt", "method_key", "label",
            "n_complexes", "n_poses",
            "pb_valid_complexes", "pb_valid_poses", "pb_valid_%",
            "rmsd2_complexes", "rmsd2_poses", "rmsd2_%",
            "kabsch1_complexes", "kabsch1_poses", "kabsch1_%",
            # Pure ≤2Å & PB-valid (paper headline) — kept in the CSV for reference.
            "rmsd2_pbvalid_complexes", "rmsd2_pbvalid_poses",
            "rmsd2_pbvalid_of_rmsd2_%", "rmsd2_pbvalid_of_all_%",
            # Strict apex: ≤2Å & PB-valid & Kabsch < thr.
            "triple_complexes", "triple_poses",
            "triple_of_rmsd2_%", "triple_of_all_%"]
    rows = []
    for tool, role, is_opt, mkey in specs:
        sub = df_full[df_full["method"].astype(str) == mkey]
        sub = _drop_bare_diffdock_dupe(sub, mkey)
        if sub.empty:
            continue
        cplx = sub["protein"].astype(str) + "\x00" + sub["ligand"].astype(str)
        pbv = _to_bool(sub["pb_valid"])
        le2 = pd.to_numeric(sub["rmsd"], errors="coerce") <= 2.0
        # Best-fit (Kabsch) RMSD: bestfit_rmsd already coalesces PoseBusters'
        # pb_kabsch_rmsd with our own fit; fall back defensively for older caches.
        if "bestfit_rmsd" in sub.columns:
            kcol = sub["bestfit_rmsd"]
        elif "pb_kabsch_rmsd" in sub.columns:
            kcol = sub["pb_kabsch_rmsd"]
        else:
            kcol = pd.Series(np.nan, index=sub.index)
        kab = pd.to_numeric(kcol, errors="coerce") < kabsch_thr
        both = le2 & pbv
        triple = both & kab
        n_poses = int(len(sub))
        pb_p, rm_p = int(pbv.sum()), int(le2.sum())
        k1_p, bo_p, tr_p = int(kab.sum()), int(both.sum()), int(triple.sum())
        _nc = lambda mask: int(cplx[mask].nunique())
        rows.append({
            "tool": tool, "role": role, "is_opt": is_opt, "method_key": mkey,
            "label": _cascade_label(mkey),
            "n_complexes": int(cplx.nunique()), "n_poses": n_poses,
            # Percentages kept at full precision — rounded once, at display / CSV-write
            # time, so a single rate can't drift a digit via double rounding.
            "pb_valid_complexes": _nc(pbv), "pb_valid_poses": pb_p,
            "pb_valid_%": 100.0 * pb_p / max(1, n_poses),
            "rmsd2_complexes": _nc(le2), "rmsd2_poses": rm_p,
            "rmsd2_%": 100.0 * rm_p / max(1, n_poses),
            "kabsch1_complexes": _nc(kab), "kabsch1_poses": k1_p,
            "kabsch1_%": 100.0 * k1_p / max(1, n_poses),
            "rmsd2_pbvalid_complexes": _nc(both), "rmsd2_pbvalid_poses": bo_p,
            "rmsd2_pbvalid_of_rmsd2_%": (100.0 * bo_p / rm_p
                                         if rm_p else float("nan")),
            "rmsd2_pbvalid_of_all_%": 100.0 * bo_p / max(1, n_poses),
            "triple_complexes": _nc(triple), "triple_poses": tr_p,
            "triple_of_rmsd2_%": (100.0 * tr_p / rm_p if rm_p else float("nan")),
            "triple_of_all_%": 100.0 * tr_p / max(1, n_poses),
        })
    return pd.DataFrame(rows, columns=cols)


def _cascade_groups(kabsch_thr: float = FORM_OK_KABSCH_A):
    """Fixed-width column layout for the cascade text/markdown table: (group, [(header,
    width, key, kind)]) where kind 'i' = thousands-separated int, 'p' = 1-dp percent
    ('—' for NaN). The Kabsch threshold is folded into the two form-based headers."""
    kv = f"{kabsch_thr:g}"
    return [
        ("Produced", [("Cplx", 5, "n_complexes", "i"),
                      ("Poses", 7, "n_poses", "i")]),
        ("PoseBusters-valid", [("Cplx", 5, "pb_valid_complexes", "i"),
                               ("Poses", 7, "pb_valid_poses", "i"),
                               ("%", 6, "pb_valid_%", "p")]),
        ("RMSD ≤ 2 Å", [("Cplx", 5, "rmsd2_complexes", "i"),
                        ("Poses", 7, "rmsd2_poses", "i"),
                        ("%", 6, "rmsd2_%", "p")]),
        (f"Kabsch RMSD < {kv} Å", [("Cplx", 5, "kabsch1_complexes", "i"),
                                   ("Poses", 7, "kabsch1_poses", "i"),
                                   ("%", 6, "kabsch1_%", "p")]),
        (f"≤ 2 Å & PB-valid & Kabsch < {kv} Å",
         [("Cplx", 5, "triple_complexes", "i"),
          ("Poses", 7, "triple_poses", "i"),
          ("% of ≤2Å", 9, "triple_of_rmsd2_%", "p"),
          ("% all", 6, "triple_of_all_%", "p")]),
    ]


def _render_pose_validity_cascade_table(cascade: pd.DataFrame,
                                        kabsch_thr: float = FORM_OK_KABSCH_A) -> str:
    """Fixed-width, grouped-header rendering of the cascade DataFrame (a pose count and
    a complex count under each stage)."""
    # Wide enough for the longest registered _cascade_label ("AutoDock Vina MGLTools-lig
    # exh32 + GNINA CNN-refinement", 54 chars). The field pads but never truncates, so a
    # label longer than this pushes every numeric column right of the header rule.
    LBL_W, GAP, GGAP = 55, "  ", "   "
    groups = _cascade_groups(kabsch_thr)

    def _fmt(v, kind: str) -> str:
        if kind == "i":
            return f"{int(v):,}"
        if v is None or (isinstance(v, float) and v != v):
            return "—"
        return f"{float(v):.1f}"

    grp_hdr = " " * LBL_W
    sub_hdr = f"{'Docking variant':<{LBL_W}}"
    for gname, fields in groups:
        block_w = sum(w for _, w, _, _ in fields) + len(GAP) * (len(fields) - 1)
        grp_hdr += GGAP + gname.center(block_w)
        sub_hdr += GGAP + GAP.join(f"{h:>{w}}" for h, w, _, _ in fields)
    rule = "-" * len(sub_hdr)

    lines = [grp_hdr, sub_hdr, rule]
    for _, r in cascade.iterrows():
        cells = f"{str(r['label']):<{LBL_W}}"
        for _, fields in groups:
            cells += GGAP + GAP.join(
                f"{_fmt(r[key], kind):>{w}}" for _, w, key, kind in fields)
        lines.append(cells)
    return "\n".join(lines)


def _write_pose_validity_cascade_report(cascade: pd.DataFrame, out_path: Path,
                                        kabsch_thr: float = FORM_OK_KABSCH_A) -> str:
    """Write the pose-level validity/accuracy cascade as a fixed-width + markdown
    table with reading notes; return the fixed-width table (also printed by main)."""
    if cascade is None or cascade.empty:
        out_path.write_text("No variants present — pose validity cascade is empty.\n")
        return ""

    kv = f"{kabsch_thr:g}"
    table = _render_pose_validity_cascade_table(cascade, kabsch_thr)
    W = max(84, len(table.splitlines()[1]))

    L: list[str] = []
    L.append("=" * W)
    L.append("Pose-level validity & accuracy cascade  (per docking variant)")
    L.append("=" * W)
    L.append("")
    L.append("Source : per_pose_metrics.csv (df_full — ALL raw + smina/gnina variants)")
    L.append("Table  : pose_validity_cascade.csv")
    L.append("")
    L.append("Every variant's ENTIRE produced-pose pool is walked through the cascade")
    L.append(f"  produced -> PoseBusters-valid -> RMSD <= 2 A -> Kabsch RMSD < {kv} A")
    L.append(f"           -> RMSD <= 2 A AND PB-valid AND Kabsch RMSD < {kv} A,")
    L.append("reported at each stage as BOTH a pose count and the number of distinct")
    L.append("(protein, ligand) complexes that still have >= 1 surviving pose.")
    L.append("")
    L.append(table)
    L.append("")
    L.append("-" * W)
    L.append("COLUMN KEY")
    L.append("-" * W)
    L.append("Cplx / Poses      = complexes with >=1 qualifying pose  /  qualifying poses.")
    L.append("Produced          = every pose the variant generated (the denominator).")
    L.append("PoseBusters-valid = poses passing ALL canonical PoseBusters checks;")
    L.append("                    %   = PB-valid poses / produced poses.")
    L.append(f"RMSD <= 2 A       = poses within 2 A {_ref_noun()} RMSD (near-native placement);")
    L.append("                    %   = <=2 A poses / produced poses.")
    L.append(f"Kabsch RMSD < {kv} A   = poses whose best-fit / SUPERPOSED heavy-atom RMSD is")
    L.append("                    < the threshold — the ligand's internal CONFORMATION is")
    L.append("                    right regardless of where it was placed (PoseBusters'")
    L.append("                    pb_kabsch_rmsd, our own Kabsch fit as fallback);")
    L.append(f"                    %   = Kabsch < {kv} A poses / produced poses.")
    L.append(f"<= 2 A & PB-valid & Kabsch < {kv} A = the strict apex: poses that are")
    L.append("                    near-native AND PB-valid (the PoseBusters-paper headline)")
    L.append("                    AND have the correct form (Kabsch < threshold);")
    L.append("                    % of <=2A = of the <=2 A poses, how many clear all three;")
    L.append("                    % all     = triple-passing poses / produced poses.")
    L.append("")
    L.append("NOTES")
    L.append("  * Pose counts are POOLED across complexes (pseudoreplicated) — they are a")
    L.append("    sampling-volume view, not a per-complex success rate. The complex counts")
    L.append("    say how many of the targets each stage actually reaches.")
    L.append(f"  * RMSD is {_ref_noun()} (as-placed) RMSD; Kabsch is the placement-free best-fit")
    L.append("    RMSD. A pose with no computable value never satisfies that threshold.")
    if REFERENCE_CONVENTION != "instance":
        L.append("  * " + _ref_convention_note())
    L.append("  * The pure '<= 2 A & PB-valid' counts (paper headline, before the Kabsch")
    L.append("    form filter) are retained in pose_validity_cascade.csv as the")
    L.append("    rmsd2_pbvalid_* columns for reference.")
    L.append("  * A degenerate 'EquiBind guided' variant (1 complex) may appear last; it is a")
    L.append("    partial run, not a full benchmark variant.")
    L.append("")

    # Paste-ready markdown (one flat header row — markdown can't span group headers).
    def _mi(v):
        return f"{int(v):,}"

    def _mp(v):
        return "—" if v is None or (isinstance(v, float) and v != v) else f"{float(v):.1f}"

    md = [
        "MARKDOWN (paste-ready)",
        "",
        ("| Docking variant | Produced cplx | Produced poses | PB-valid cplx | "
         "PB-valid poses | PB-valid % | ≤2Å cplx | ≤2Å poses | ≤2Å % | "
         f"Kabsch<{kv}Å cplx | Kabsch<{kv}Å poses | Kabsch<{kv}Å % | "
         f"≤2Å&PB&Kabsch<{kv}Å cplx | ≤2Å&PB&Kabsch<{kv}Å poses | "
         "% of ≤2Å | % of all |"),
        ("|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|"),
    ]
    for _, r in cascade.iterrows():
        md.append(
            f"| {r['label']} | {_mi(r['n_complexes'])} | {_mi(r['n_poses'])} | "
            f"{_mi(r['pb_valid_complexes'])} | {_mi(r['pb_valid_poses'])} | {_mp(r['pb_valid_%'])} | "
            f"{_mi(r['rmsd2_complexes'])} | {_mi(r['rmsd2_poses'])} | {_mp(r['rmsd2_%'])} | "
            f"{_mi(r['kabsch1_complexes'])} | {_mi(r['kabsch1_poses'])} | {_mp(r['kabsch1_%'])} | "
            f"{_mi(r['triple_complexes'])} | {_mi(r['triple_poses'])} | "
            f"{_mp(r['triple_of_rmsd2_%'])} | {_mp(r['triple_of_all_%'])} |")
    L.extend(md)
    L.append("")

    out_path.write_text("\n".join(L))
    return table


def plot_pbvalid_yield_boxplot(per: pd.DataFrame, out: Path,
                               summary: "pd.DataFrame | None" = None,
                               stats: dict | None = None,
                               specs=_PBVALID_YIELD_SPECS,
                               *, draw_brackets: bool = True,
                               color_by_variant: bool = False,
                               figsize: tuple = (11.4, 6.6),
                               group_labels: "dict | None" = None,
                               legend_loc: str = "lower center",
                               legend_anchor: tuple = (0.5, 1.012),
                               legend_ncol: int = 3) -> None:
    """Fig 09f — box/whisker of the produced-pose PoseBusters-valid YIELD per variant,
    raw beside its optimized variants, grouped by tool.

    Each box is the distribution across complexes of that variant's per-complex % of
    produced poses that pass PoseBusters. Optimization (smina/gnina) refines the SAME
    poses raw produced, so the upward shift raw → opt is the validity it recovers.
    Raw boxes are hatched (before); optimized/reference boxes are solid (after). Paired
    Wilcoxon signed-rank brackets (Holm-corrected across all raw→opt contrasts) mark
    the significance and the median percentage-point gain. Outliers beyond the whiskers
    are hidden so the boxes stay legible.

    ``draw_brackets`` off, ``color_by_variant`` on and a wider ``figsize`` render the
    all-variants overview (every EquiBind pocket×refine variant), where per-tool raw→opt
    brackets no longer make sense; boxes are then coloured per variant instead of per tool.
    """
    if per is None or per.empty:
        return
    # Per-variant counts for the tick labels: # complexes with >=1 PB-valid pose, and the
    # # of PB-valid poses. The total n_complexes is the same for every variant (it lives in
    # the title).
    if (summary is not None and not summary.empty
            and "n_complexes_with_valid" in summary.columns):
        cover = dict(zip(summary["method_key"], summary["n_complexes_with_valid"]))
        n_valid = dict(zip(summary["method_key"], summary["n_valid"]))
    else:
        agg = {m: g for m, g in per.groupby("method_key")}
        cover = {m: int((g["n_valid"] >= 1).sum()) for m, g in agg.items()}
        n_valid = {m: int(g["n_valid"].sum()) for m, g in agg.items()}

    tool_color_key = {"autodock": "autodock", "autodock_vinardo": "autodock_vinardo",
                      "diffdock": "diffdock", "unidock": "unidock",
                      "unidock2": "unidock2", "equibind": "equibind_unguided"}
    group_name = group_labels or {"autodock": "AutoDock Vina", "diffdock": "DiffDock",
                                  "unidock": "Uni-Dock", "unidock2": "Uni-Dock2",
                                  "equibind": "EquiBind (blind / unguided)"}

    # Assemble the boxes present, in spec order, laying out x with a wider gap between
    # tools (same scheme as 09c).
    intra, inter, box_w = 0.90, 0.85, 0.62
    boxes, positions, spans = [], [], {}
    x, prev = 0.0, None
    xmap: dict[str, float] = {}
    for tool, role, is_opt, mkey in specs:
        data = per.loc[per["method_key"] == mkey, "pct_valid"].to_numpy(float)
        data = data[~np.isnan(data)]
        if data.size == 0:
            continue
        if prev is not None and tool != prev:
            x += inter
        q1, q3 = np.percentile(data, [25, 75])
        inside = data[data <= q3 + 1.5 * (q3 - q1)]
        whi = float(inside.max()) if inside.size else float(q3)
        color = (TOOL_COLORS.get(mkey, "#888888") if color_by_variant
                 else TOOL_COLORS.get(tool_color_key.get(tool, tool), "#888888"))
        boxes.append({"tool": tool, "role": role, "is_opt": is_opt, "mkey": mkey,
                      "data": data, "x": x, "n": int(data.size),
                      "median": float(np.median(data)), "whisker_hi": whi,
                      "color": color})
        positions.append(x); xmap[mkey] = x
        spans.setdefault(tool, []).append(x)
        x += intra; prev = tool
    if not boxes:
        return

    fig, ax = plt.subplots(figsize=figsize)
    bp = ax.boxplot([b["data"] for b in boxes], positions=positions, widths=box_w,
                    patch_artist=True, showmeans=True, showfliers=False,
                    medianprops=dict(color="#222222", lw=1.6),
                    meanprops=dict(marker="D", markerfacecolor="white",
                                   markeredgecolor="#222222", markersize=5))
    for patch, b in zip(bp["boxes"], boxes):
        patch.set_facecolor(b["color"])
        if b["is_opt"]:
            patch.set_alpha(0.72)
        else:                                   # raw / reference read as "before"
            patch.set_alpha(0.32)
            patch.set_hatch("///")
        patch.set_edgecolor(b["color"]); patch.set_linewidth(1.4)

    # MEDIAN per-complex yield printed at the LEFT of the median line, inside the box. The
    # text starts at the box's left edge so its white outline covers the median line's left
    # tip (otherwise the line pokes out as a stray dash that reads like a minus sign).
    import matplotlib.patheffects as _pe
    for b in boxes:
        ax.text(b["x"] - box_w / 2, b["median"], f"{b['median']:.1f}%",
                ha="left", va="center", fontsize=9.0, fontweight="bold",
                color="#111111", zorder=6,
                path_effects=[_pe.withStroke(linewidth=3.2, foreground="white")])

    # Significance brackets: raw → each optimized variant, per tool, stacked. Skipped in
    # the all-variants overview (multiple raws per tool → per-tool brackets don't apply).
    bracket_ceiling = 0.0
    if draw_brackets:
        group_top = {t: max(b["whisker_hi"] for b in boxes if b["tool"] == t)
                     for t in spans}
        for tool in spans:
            raw_key = _PBVALID_YIELD_RAW.get(tool)
            if raw_key not in xmap:
                continue
            opt_boxes = [b for b in boxes if b["tool"] == tool and b["is_opt"]]
            base = max(group_top[tool] + 3.0, 24.0)     # start just above the tallest box
            # Vertical spacing between stacked brackets must clear a full two-line label
            # (star + gain on line 1, p-value on line 2). At the 10 pt bracket font a
            # two-line label spans ~11 data units, so stack levels 13 units apart to keep
            # each label from colliding with the bracket line above it.
            for lvl, b in enumerate(sorted(opt_boxes, key=lambda z: z["x"])):
                rec = _pbvalid_yield_contrast(stats, tool, b["role"])
                y = base + lvl * 13.0
                x0, x1 = xmap[raw_key], b["x"]
                ax.plot([x0, x0, x1, x1], [y - 1.8, y, y, y - 1.8],
                        color="0.35", lw=1.1, clip_on=False, zorder=4)
                if rec is not None:
                    star = rec.get("star", su.p_stars(rec.get("p_holm", rec["p_raw"])))
                    pv = rec.get("p_holm", rec["p_raw"])
                    txt = (f"{star}  +{rec['median_gain_pp']:.1f} pp\n{su.fmt_p(pv)}"
                           if rec["median_gain_pp"] >= 0
                           else f"{star}  {rec['median_gain_pp']:.1f} pp\n{su.fmt_p(pv)}")
                else:
                    dv = b["median"] - next(z["median"] for z in boxes
                                            if z["mkey"] == raw_key)
                    txt = f"+{dv:.1f} pp"
                ax.text((x0 + x1) / 2, y + 0.4, txt, ha="center", va="bottom",
                        fontsize=10, color="#1f7a1f", fontweight="bold",
                        clip_on=False, zorder=5)
                bracket_ceiling = max(bracket_ceiling, y + 12.5)

    # Per-box tick labels + a second tier of tool names. The total n_complexes is the same
    # for every variant (it lives in the title). Line 2 = counts: complexes with >=1
    # PB-valid pose / number of PB-valid poses (e.g. DiffDock gnina-opt = 298/7276). Font
    # shrinks a touch when many variants are shown.
    _tick_fs = 9.6 if len(boxes) > 8 else 10.6
    ax.set_xticks(positions)
    ax.set_xticklabels(
        [f"{b['role']}\n{int(cover.get(b['mkey'], 0))}/{int(n_valid.get(b['mkey'], 0))}"
         for b in boxes], fontsize=_tick_fs)
    ax.tick_params(axis="x", length=0, pad=6)
    ax.tick_params(axis="y", labelsize=11.5)
    trans = ax.get_xaxis_transform()
    for tool, xs in spans.items():
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - box_w / 2, hi + box_w / 2], [-0.205, -0.205],
                transform=trans, color="0.45", lw=1.0, clip_on=False, zorder=1)
        # .get, not [] — a group key added upstream must degrade to a readable label,
        # never abort the whole figure run partway through (this is how the missing
        # MGLTools entry took out 09f and every figure after it).
        ax.text(xc, -0.225, group_name.get(tool, str(tool)),
                transform=trans, ha="center",
                va="top", fontsize=13.0, fontweight="bold", clip_on=False)

    ax.axhline(0, color="0.7", lw=0.8)
    ax.set_ylabel("PoseBusters-valid poses produced\n(% of a complex's poses, per complex)",
                  fontsize=13.0)
    # The y-axis is a percentage, so it ENDS AT 100 %: ticks stop at 100, the y-axis line
    # (left spine) is bounded to 0–100, and the top/right frame is dropped. Headroom above
    # 100 is kept to the minimum needed for the significance brackets (7-variant figure);
    # with no brackets the axis is capped just above 100 so the title sits close to the plot.
    ax.set_ylim(-3, (bracket_ceiling + 2.0) if (draw_brackets and bracket_ceiling) else 103.0)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_bounds(0, 100)
    ax.spines["bottom"].set_bounds(min(positions) - 0.75, max(positions) + 0.75)
    ax.set_xlim(min(positions) - 0.75, max(positions) + 0.75)
    ax.grid(axis="y", alpha=0.3)

    from matplotlib.patches import Patch
    ax.legend(handles=[
        Patch(facecolor="0.5", alpha=0.32, hatch="///", edgecolor="0.5",
              label="raw / reference (before optimization)"),
        Patch(facecolor="0.5", alpha=0.72, edgecolor="0.5",
              label="post-hoc optimized (smina / gnina)"),
        plt.Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
                   markeredgecolor="#222222", markersize=8, label="mean")],
        loc=legend_loc, bbox_to_anchor=legend_anchor, ncol=legend_ncol,
        fontsize=11, framealpha=0.92)

    # (The per-tool Friedman omnibus stats box was removed per request; the full
    # statistics live in the companion report 09f_pbvalid_yield_report.txt.)

    # Title carries the shared complex count (same n for every variant). When the legend
    # is placed above the plot (anchor y >= 1.0 → "between title and graph"), the top of
    # the figure is reserved for title + legend so they don't overlap the axes.
    _legend_above = legend_anchor[1] >= 1.0
    title_n = max((b["n"] for b in boxes), default=0)
    if _legend_above:
        # Legend rides just above the axes; reserve a thin top band for it so it clears
        # the plot, then drop the title to sit right on top of the legend. tight_layout
        # runs FIRST so the legend/axes settle, then the title is placed just above the
        # (already-drawn) legend so there is no empty gap between title and graph.
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        fig.canvas.draw()
        _leg = ax.get_legend()
        _leg_top = (_leg.get_window_extent(fig.canvas.get_renderer())
                    .transformed(fig.transFigure.inverted()).y1) if _leg else 0.90
        _title_y = min(0.995, _leg_top + 0.055)   # small gap above the legend
        fig.suptitle(_vt(f"Search Effort, Post-Hoc Optimization and PoseBusters Validity   "
                         f"(n = {title_n} complexes)"),
                     fontsize=16.5, fontweight="bold", y=_title_y, va="top")
    else:
        fig.suptitle(_vt(f"Search Effort, Post-Hoc Optimization and PoseBusters Validity   "
                         f"(n = {title_n} complexes)"),
                     fontsize=16.5, fontweight="bold", y=0.995)
        fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── Top-1 vs top-N poses, raw vs refined — ranking headroom + optimization (09d) ──


def _av_from_reps(reps: pd.DataFrame, thr: float) -> dict:
    """Accuracy/validity counts for ONE representative pose per complex.

    ``reps`` must already hold at most one row per (protein, ligand). Returns
    n_complexes, the count and % whose pose is RMSD ≤ ``thr`` (accurate), the count
    and % that are ALSO PB-valid (accurate + valid), and the median RMSD.
    """
    reps = reps.dropna(subset=["rmsd"])
    n = reps.groupby(["protein", "ligand"]).ngroups if len(reps) else 0
    if n == 0:
        return {"n": 0}
    near = reps["rmsd"] <= thr
    valid = near & _to_bool(reps["pb_valid"])
    return {"n": n,
            "n_near": int(near.sum()), "near_pct": round(100 * float(near.mean()), PCT_DECIMALS),
            "n_valid": int(valid.sum()), "valid_pct": round(100 * float(valid.mean()), PCT_DECIMALS),
            "median": round(float(reps["rmsd"].median()), 3)}


def _best_of_top_d(sub: pd.DataFrame, rank_series: pd.Series,
                   d: int) -> pd.DataFrame:
    """The lowest-RMSD pose among rank ≤ ``d`` — one row per complex (``d=1`` → the
    rank-1 pose). Accuracy is monotone non-decreasing in ``d`` (a min over a growing
    candidate set); ``rank_series`` is the variant's per-row 1-based ranking."""
    sub = sub.copy()
    sub["_rk"] = pd.to_numeric(rank_series, errors="coerce")
    cand = sub[sub["_rk"] <= d].dropna(subset=["rmsd"])
    if not len(cand):
        return cand
    return (cand.loc[cand.groupby(["protein", "ligand"])["rmsd"].idxmin()]
                .reset_index(drop=True))


def _gnina_affinity_rank(sub: pd.DataFrame) -> pd.Series:
    """Per-complex 1-based rank by gnina affinity (most-negative = rank 1). NaN
    affinities sort last, broken by generation index, so every pose gets a rank —
    this is the ranking EquiBind gains from gnina optimization. Indexed like ``sub``."""
    sub = sub.copy()
    sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
    sub["_aff"] = pd.to_numeric(sub["gnina_affinity"], errors="coerce")
    order = sub.sort_values(["_aff", "_gi"], na_position="last", kind="mergesort")
    order["_sr"] = order.groupby(["protein", "ligand"]).cumcount() + 1
    return order["_sr"].reindex(sub.index)


# Fig 09d bars: (tool colour key, variant label, method key, ranking label, rank kind).
# For DiffDock/EquiBind both the raw and the optimised run appear so raw vs. refined
# can be read off directly; needs the FULL (un-collapsed) frame.
_RANK1_TOPN_SPECS = [
    # AutoDock = the dominant arm (ADFRsuite ligands, exhaustiveness 128) and its own
    # raw counterpart, so raw vs refined still reads off directly.
    ("autodock", "AutoDock Vina",        _DOMINANT_AUTODOCK_RAW, "confidence rank",    "native"),
    ("autodock", "AutoDock (gnina-opt)", _DOMINANT_AUTODOCK_OPT,
     "optimized rank", "native"),
    ("diffdock", "DiffDock (raw)",       "diffdock",                "confidence rank",     "native"),
    ("diffdock", "DiffDock (gnina-opt)", "diffdock_gnina",          "confidence rank",     "native"),
    ("equibind", "EquiBind (raw)",       "equibind_unguided_raw",   "generation order",    "generation"),
    ("equibind", "EquiBind (gnina-opt)", "equibind_unguided_gnina", "gnina-affinity rank", "gnina"),
]


def aggregate_rank1_vs_topn(df: pd.DataFrame, depths,
                            thr: float = 2.0) -> pd.DataFrame:
    """Per tool variant, one bar per depth in ``depths`` (e.g. [1, 15, 30]) — the
    best (lowest-RMSD) pose among the variant's top-``d`` poses — on the same
    accuracy/validity axes as fig 09c, for BOTH the raw and the optimised run so
    raw-vs-refined reads off directly. Depth 1 = the rank-1 pose; deeper depths show
    the ranking headroom (top-30 ≈ the full-oracle ceiling for ~30-pose methods).

    "top-``d``" is the lowest-RMSD pose among the variant's first ``d`` poses under
    its ranking:
      * AutoDock Vina — its confidence ranking (single variant; reference).
      * DiffDock (raw) / (smina-opt) — DiffDock's confidence ranking on the raw vs.
        smina-refined geometry (``diffdock`` / ``diffdock_smina``).
      * EquiBind (raw) — GENERATION order (EquiBind has no native ranking, so top-d
        = best of the first d generated poses).
      * EquiBind (gnina-opt) — ranked by gnina AFFINITY (the ranking gnina supplies).

    Needs the FULL (un-collapsed) frame so every variant is present under its own
    method key; missing variants are silently skipped. Returns one row per bar.
    """
    acc_col = f"rmsd_le_{thr:g}A_%"
    val_col = f"pb_valid_and_rmsd{thr:g}_%"
    present = set(df["method"].astype(str))
    depths = sorted({int(x) for x in depths})
    rows = []

    def _emit(tool, variant, method, ranking, selection, depth, reps):
        d = _av_from_reps(reps, thr)
        if d.get("n", 0) == 0:
            return
        rows.append({
            "tool": tool, "variant": variant, "method_key": method,
            "ranking": ranking, "selection": selection, "depth": depth,
            "n_complexes": d["n"], "n_rmsd_le_2A": d["n_near"], acc_col: d["near_pct"],
            "n_pb_valid_and_rmsd2": d["n_valid"], val_col: d["valid_pct"],
            "median_rmsd": d["median"],
        })

    for tool, variant, mkey, ranking, kind in _RANK1_TOPN_SPECS:
        if mkey not in present:
            continue
        sub = df[df["method"].astype(str) == mkey].copy()
        if kind == "native":
            rk = sub["rank"]
        elif kind == "generation":
            rk = sub["pose_name"].map(_equibind_pose_index)
        else:  # gnina
            rk = _gnina_affinity_rank(sub)
        for d in depths:
            selection = "top-1" if d == 1 else f"top-{d}"
            _emit(tool, variant, mkey, ranking, selection, d,
                  _best_of_top_d(sub, rk, d))

    return pd.DataFrame(rows)


def plot_rank1_vs_topn(agg: pd.DataFrame, out: Path,
                       thr: float = 2.0) -> None:
    """Fig 09d — grouped accuracy-vs-validity bars: per tool variant, one bar per
    depth (top-1, top-15, top-30 …), with the RAW and the smina-optimised run of
    DiffDock/EquiBind side by side so ranking headroom (deeper top-d) and
    optimisation gain (raw → refined) read off together. Same light(accuracy) +
    dark(PB-valid subset) encoding as fig 09c; colour = tool. The Δ callout is the
    (always-positive) accuracy gain from top-1 to the deepest depth; validity is the
    dark bar — note it barely moves for the RAW runs (more poses find near-native
    geometry, but only optimisation makes it valid). See ``aggregate_rank1_vs_topn``.
    """
    if agg is None or agg.empty:
        return
    acc_col = f"rmsd_le_{thr:g}A_%"
    val_col = f"pb_valid_and_rmsd{thr:g}_%"
    variant_order = [s[1] for s in _RANK1_TOPN_SPECS]
    color_key = {"autodock": "autodock", "diffdock": "diffdock",
                 "equibind": "equibind_unguided"}

    agg = agg[agg["variant"].isin(variant_order)].copy()
    if agg.empty:
        return
    agg["_vo"] = agg["variant"].map({v: i for i, v in enumerate(variant_order)})
    agg = agg.sort_values(["_vo", "depth"], kind="mergesort").reset_index(drop=True)
    n_depths = agg["depth"].nunique()

    # Lay out: one bar per depth in each variant group, a small gap between variants
    # of the same tool, a larger gap between tools (raw/refined read as one tool).
    bar_w, intra, gap_same, gap_tool = 0.52, 0.62, 0.55, 1.05
    positions, spans, prev_v, prev_tool, x = [], {}, None, None, 0.0
    for _, r in agg.iterrows():
        v, tool = r["variant"], r["tool"]
        if prev_v is not None and v != prev_v:
            x += gap_tool if tool != prev_tool else gap_same
        positions.append(x)
        spans.setdefault(v, {"xs": [], "tool": tool})["xs"].append(x)
        x += intra
        prev_v, prev_tool = v, tool
    agg["_x"] = positions

    fig, ax = plt.subplots(figsize=(max(14.5, 1.05 * len(agg) + 3), 6.8))
    for _, r in agg.iterrows():
        c = TOOL_COLORS.get(color_key[r["tool"]], "#888888")
        xi, acc, val = r["_x"], float(r[acc_col]), float(r[val_col])
        ax.bar(xi, acc, width=bar_w, color=c, alpha=0.32, edgecolor=c,
               linewidth=1.4, zorder=2)
        ax.bar(xi, val, width=bar_w, color=c, edgecolor="white",
               linewidth=0.8, zorder=3)
        ax.text(xi, acc + 1.5, f"{acc:.0f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold", zorder=4)
        if val > 4:
            ax.text(xi, val / 2, f"{val:.0f}", ha="center", va="center",
                    fontsize=8, color="white", fontweight="bold", zorder=4)

    # Accuracy gain top-1 → deepest depth (provably ≥ 0), above each group's last bar.
    for v, info in spans.items():
        sub = agg[agg["variant"] == v].sort_values("depth")
        if len(sub) < 2:
            continue
        base, last = sub.iloc[0], sub.iloc[-1]
        da = float(last[acc_col]) - float(base[acc_col])
        xb, yb = float(last["_x"]), float(last[acc_col])
        ax.annotate(f"▲ +{da:.0f} pp ≤ 2 Å\n(top-1→top-{int(last['depth'])})",
                    xy=(xb, yb), xytext=(0, 12), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5, color="#1f7a1f",
                    fontweight="bold")

    ax.set_xticks(agg["_x"].tolist())
    ax.set_xticklabels([f"{r['selection']}\nn={int(r['n_complexes'])}"
                        for _, r in agg.iterrows()], fontsize=7.5)
    ax.tick_params(axis="x", length=0, pad=6)
    trans = ax.get_xaxis_transform()
    for v, info in spans.items():
        xs = info["xs"]
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - bar_w / 2, hi + bar_w / 2], [-0.14, -0.14],
                transform=trans, color="0.45", lw=1.0, clip_on=False, zorder=1)
        ax.text(xc, -0.155, v, transform=trans, ha="center", va="top",
                fontsize=9.5, fontweight="bold", clip_on=False)

    ax.set_ylabel("% of receptor-ligand complexes")
    ax.set_ylim(0, 105)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.grid(axis="y", alpha=0.3)
    _accuracy_validity_legend(fig)
    depth_lbl = ", ".join("top-" + str(d) for d in sorted(agg["depth"].unique()))
    fig.suptitle(f"Ranking headroom & optimization — {depth_lbl}, raw vs. refined\n"
                 "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.10)
    fig.text(0.5, -0.11,
             "Representative pose per complex — top-d = the lowest-RMSD pose among the "
             "variant's top-d poses (top-1 = the rank-1 pose; ~30 poses exist, so top-30 "
             "≈ the full oracle). Ranking used: Vina & DiffDock — native confidence rank "
             "(smina-opt keeps DiffDock's ranking on smina-refined geometry); EquiBind "
             "raw — generation order (no ranking); EquiBind gnina-opt — gnina affinity. "
             "'raw' = un-optimised geometry, 'gnina-opt' = gnina-refined for EquiBind, "
             "'smina-opt' = smina-refined for DiffDock.",
             ha="center", va="top", fontsize=8, color="0.35", wrap=True)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── 09d alternative (non-bar) views — same rank1_vs_topn aggregation ─────────
# Each variant is its own series here (a trajectory / slope / dumbbell row), so
# unlike the 09d bars (colour = tool, raw/refined by position) raw and refined get
# distinct shades of the tool colour to stay distinguishable.
_RANK1_TOPN_COLORS = {
    "AutoDock Vina":        "#1f77b4",
    # Without this entry the refined AutoDock series falls through to the "#888888"
    # default and reads as an un-keyed grey line rather than as the blue tool family.
    "AutoDock (gnina-opt)": "#08519c",
    "DiffDock (raw)":       "#ff7f0e",
    "DiffDock (gnina-opt)": "#d95f02",
    "EquiBind (raw)":       "#2ca02c",
    "EquiBind (gnina-opt)": "#1b7837",
}


def _r1tn_common(agg: pd.DataFrame, thr: float):
    """(variants in canonical order, sorted depths, acc_col, val_col, depth-label fn)."""
    acc_col = f"rmsd_le_{thr:g}A_%"
    val_col = f"pb_valid_and_rmsd{thr:g}_%"
    present = set(agg["variant"])
    variants = [s[1] for s in _RANK1_TOPN_SPECS if s[1] in present]
    depths = sorted({int(d) for d in agg["depth"].unique()})
    lbl = lambda d: "top-1" if d == 1 else f"top-{d}"
    return variants, depths, acc_col, val_col, lbl


def plot_rank1_vs_topn_scatter(agg: pd.DataFrame, out: Path, thr: float = 2.0) -> None:
    """Fig 09d-alt (scatter) — accuracy (x) vs. PB-validity (y) with the y=x "fully
    valid" diagonal; each variant is a top-1 → … → deepest-depth trajectory. The
    distance below the diagonal is the accurate-but-invalid gap; the arrow is the
    ranking headroom (→ more accurate, ↑ more valid)."""
    if agg is None or agg.empty:
        return
    from matplotlib.lines import Line2D
    variants, depths, acc_col, val_col, lbl = _r1tn_common(agg, thr)
    if not variants:
        return
    fig, ax = plt.subplots(figsize=(9.0, 8.4))
    lim = 95
    ax.plot([0, lim], [0, lim], ls="--", color="0.55", lw=1.2, zorder=1)
    ax.fill_between([0, lim], [0, lim], 0, color="#d1495b", alpha=0.05, zorder=0)
    ax.text(lim - 2, lim - 2, "fully valid (validity = accuracy)", rotation=45,
            rotation_mode="anchor", ha="right", va="bottom", fontsize=8.5,
            color="0.45", style="italic")
    ax.text(70, 26, "accurate but\nPB-invalid", ha="center", va="center",
            fontsize=10, color="#a83244", style="italic", alpha=0.8)
    for v in variants:
        c = _RANK1_TOPN_COLORS.get(v, "#888888")
        s = agg[agg["variant"] == v].sort_values("depth")
        xs = s[acc_col].astype(float).tolist()
        ys = s[val_col].astype(float).tolist()
        ax.plot(xs, ys, color=c, lw=1.8, alpha=0.85, zorder=3)
        if len(xs) >= 2:
            ax.annotate("", xy=(xs[-1], ys[-1]), xytext=(xs[-2], ys[-2]),
                        arrowprops=dict(arrowstyle="-|>", color=c, lw=1.8,
                                        shrinkA=3, shrinkB=3), zorder=3)
        for i, (a, b) in enumerate(zip(xs, ys)):
            first, last = i == 0, i == len(xs) - 1
            ax.scatter([a], [b], s=90 if last else 55,
                       facecolors="white" if first else c, edgecolors=c,
                       linewidths=2 if first else 1, zorder=4)
        dx, dy = (5, -9) if v.startswith("EquiBind") else (5, 6)
        ax.annotate(v, xy=(xs[-1], ys[-1]), xytext=(dx, dy),
                    textcoords="offset points", fontsize=9, fontweight="bold",
                    color=c, zorder=5)
    ax.set_xlim(-2, lim); ax.set_ylim(-2, lim); ax.set_aspect("equal")
    ax.set_xlabel("Accuracy — % of complexes with a pose RMSD ≤ 2 Å")
    ax.set_ylabel("Validity — % with a pose RMSD ≤ 2 Å AND PoseBusters-valid")
    ax.grid(alpha=0.25); ax.set_axisbelow(True)
    handles = [Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.35", mew=2,
                      ms=9, label=f"{lbl(depths[0])} (open)"),
               Line2D([0], [0], marker="o", ls="", color="0.35", ms=10,
                      label=f"{lbl(depths[-1])} (filled, arrow)"),
               Line2D([0], [0], color="0.35", lw=1.8,
                      label=" → ".join(lbl(d) for d in depths))]
    ax.legend(handles=handles, loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_title("Accuracy vs. PoseBuster validity — "
                 f"{' → '.join(lbl(d) for d in depths)} trajectory\n"
                 "closer to the diagonal = more of the accurate poses are valid",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_rank1_vs_topn_slopegraph(agg: pd.DataFrame, out: Path, thr: float = 2.0) -> None:
    """Fig 09d-alt (slopegraph) — two panels (accuracy | validity), each variant a
    line across the depths (raw dashed / refined solid, colour = variant). The slope
    is the ranking headroom; a flat deep segment means the ceiling is reached."""
    if agg is None or agg.empty:
        return
    from matplotlib.lines import Line2D
    variants, depths, acc_col, val_col, lbl = _r1tn_common(agg, thr)
    if not variants:
        return
    xpos = {d: i for i, d in enumerate(depths)}
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.6), sharey=True)
    for ax, (col, name) in zip(axes, [(acc_col, "Accuracy (RMSD ≤ 2 Å)"),
                                      (val_col, "Validity (RMSD ≤ 2 Å & PB-valid)")]):
        for v in variants:
            c = _RANK1_TOPN_COLORS.get(v, "#888888")
            s = agg[agg["variant"] == v].sort_values("depth")
            xs = [xpos[int(d)] for d in s["depth"]]
            ys = s[col].astype(float).tolist()
            ax.plot(xs, ys, color=c, lw=2.4, ls="--" if "raw" in v else "-",
                    marker="o", ms=6.5, zorder=3)
            ax.annotate(f"{ys[0]:.0f}%", (xs[0], ys[0]), xytext=(-6, 0),
                        textcoords="offset points", ha="right", va="center",
                        fontsize=8, color=c, fontweight="bold")
            ax.annotate(f"{ys[-1]:.0f}%  {v}", (xs[-1], ys[-1]), xytext=(6, 0),
                        textcoords="offset points", ha="left", va="center",
                        fontsize=8, color=c, fontweight="bold")
        ax.set_xlim(-0.6, len(depths) - 1 + 1.9); ax.set_ylim(-3, 95)
        ax.set_xticks(list(xpos.values()))
        ax.set_xticklabels([lbl(d) for d in depths], fontsize=10)
        ax.set_title(name, fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.25); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel("% of receptor-ligand complexes")
    fig.legend(handles=[Line2D([0], [0], color="0.35", lw=2.4, ls="-", label="refined (smina-opt)"),
                        Line2D([0], [0], color="0.35", lw=2.4, ls="--", label="raw")],
               loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.02), fontsize=9)
    fig.suptitle(f"Ranking headroom — {' → '.join(lbl(d) for d in depths)}  "
                 "(slope = gain from more poses; flat deep segment = ceiling reached)",
                 fontsize=13, fontweight="bold", y=1.06)
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_rank1_vs_topn_dumbbell(agg: pd.DataFrame, out: Path, thr: float = 2.0) -> None:
    """Fig 09d-alt (dumbbell) — one row per (variant, depth): validity ● —— ○ accuracy.
    The connecting segment IS the accurate-but-invalid share (labelled −N invalid)."""
    if agg is None or agg.empty:
        return
    from matplotlib.lines import Line2D
    variants, depths, acc_col, val_col, lbl = _r1tn_common(agg, thr)
    if not variants:
        return
    rows = []
    for v in variants:
        for d in depths:
            r = agg[(agg["variant"] == v) & (agg["depth"] == d)]
            if len(r):
                rows.append((v, lbl(d), float(r[acc_col].iloc[0]),
                             float(r[val_col].iloc[0])))
    rows = rows[::-1]
    ys = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(10.5, 0.42 * len(rows) + 2.2))
    for y, (v, sel, a, val) in zip(ys, rows):
        c = _RANK1_TOPN_COLORS.get(v, "#888888")
        ax.plot([val, a], [y, y], color=c, lw=2.5, alpha=0.5, zorder=2)
        ax.scatter([a], [y], s=90, facecolors="white", edgecolors=c, linewidths=2, zorder=3)
        ax.scatter([val], [y], s=90, color=c, edgecolors="white", linewidths=1, zorder=3)
        ax.annotate(f"{a:.0f}", (a, y), xytext=(7, 0), textcoords="offset points",
                    va="center", fontsize=8, color=c)
        gap = a - val
        if gap > 3:
            ax.annotate(f"−{gap:.0f} invalid", ((a + val) / 2, y), xytext=(0, 6),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=7, color="#a83244")
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{v} · {sel}" for v, sel, _, _ in rows], fontsize=8.5)
    ax.set_xlim(-2, 95); ax.set_xlabel("% of receptor-ligand complexes")
    ax.grid(axis="x", alpha=0.25); ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.legend(handles=[Line2D([0], [0], marker="o", ls="", color="0.35", ms=9,
                              label="validity (RMSD ≤ 2 Å & PB-valid)"),
                       Line2D([0], [0], marker="o", ls="", mfc="white", mec="0.35",
                              mew=2, ms=9, label="accuracy (RMSD ≤ 2 Å)")],
              loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_title("Accuracy vs. validity gap — the segment is the accurate-but-invalid share",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Plots — Part B: ranking quality (Vina + DiffDock only)
# ───────────────────────────────────────────────────────────────────


def plot_rank_success_curve(rank_df: pd.DataFrame, top_n: int, out: Path) -> None:
    """Success rate at each rank k (k = 1 .. top_n) for each ranking tool.

    A flat/declining curve means later ranks add no value; a step-like curve
    means rank-1 already captures most successes.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    col = "success_pct_rmsd_2.0A"
    for method, sub in rank_df.groupby("method"):
        sub = sub.sort_values("rank")
        n_pairs = int(sub["n_pairs_total"].iloc[0]) if "n_pairs_total" in sub else len(sub)
        ax.plot(sub["rank"], sub[col],
                marker="o", lw=2,
                label=f"{TOOL_LABEL.get(method, method)} (n={n_pairs})",
                color=TOOL_COLORS.get(method))
    ax.set_xticks(range(1, top_n + 1))
    ax.set_xlabel("Rank k — the tool's k-th ranked pose")
    ax.set_ylabel("% of pairs with RMSD ≤ 2 Å at exactly rank k")
    ax.set_title(_vt("Ranking quality — success rate at each pose rank\n"
                     "(rank 1 = tool's top-scored pose)"))
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_plif_recovery_by_rank(ifp_rank_df: pd.DataFrame, top_n: int, out: Path) -> None:
    """Interaction-profile similarity to the crystal at each pose rank k.

    Companion to :func:`plot_rank_success_curve` (which plots geometric success):
    one line per ranking tool, up to two panels sharing the rank axis —
      (A) ProLIF interaction-fingerprint Tanimoto similarity vs the crystal ligand
      (B) native-contact recovery (fraction of the crystal's contact residues the
          pose reproduces).
    A declining curve means the tool's top-ranked poses reproduce the crystal
    interaction profile better than its lower-ranked ones. Needs the crystal IFP
    (benchmark datasets only); no-ops when neither recovery column is present/finite.
    """
    if ifp_rank_df is None or ifp_rank_df.empty:
        return
    panels = [
        ("mean_plif_recovery", "o", "-",
         "Mean ProLIF interaction-fingerprint\nTanimoto similarity to crystal",
         "Interaction-fingerprint (ProLIF) recovery vs crystal by pose rank"),
        ("mean_contact_recovery", "s", "--",
         "Mean fraction of crystal contact\nresidues recovered",
         "Native-contact recovery vs crystal by pose rank"),
    ]
    panels = [p for p in panels
              if p[0] in ifp_rank_df.columns and ifp_rank_df[p[0]].notna().any()]
    if not panels:
        return
    fig, axes = plt.subplots(1, len(panels), figsize=(7.2 * len(panels), 5),
                             squeeze=False)
    axes = list(axes[0])
    for ax, (col, marker, ls, ylabel, title) in zip(axes, panels):
        for method, sub in ifp_rank_df.groupby("method"):
            sub = sub.sort_values("rank")
            if sub[col].notna().sum() == 0:
                continue
            n_pairs = int(sub["n_pairs_total"].iloc[0]) if "n_pairs_total" in sub else len(sub)
            ax.plot(sub["rank"], sub[col], marker=marker, ls=ls, lw=2,
                    label=f"{TOOL_LABEL.get(method, method)} (n={n_pairs})",
                    color=TOOL_COLORS.get(method))
        ax.set_xticks(range(1, top_n + 1))
        ax.set_xlabel("Rank k — the tool's k-th ranked pose")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0, 1)
        ax.set_title(_vt(title))
        ax.legend(); ax.grid(alpha=0.3)
    if len(axes) > 1:
        _label_panels(axes)
    fig.text(0.5, 0.005,
             "Rank basis: AutoDock Vina / DiffDock — native confidence rank;  "
             "EquiBind — gnina-affinity rank (most-negative = rank 1)",
             ha="center", va="bottom", fontsize=8, color="0.35")
    fig.tight_layout(rect=(0, 0.05, 1, 1)); fig.savefig(out, dpi=160); plt.close(fig)


def plot_cumulative_oracle_curve(rank_df: pd.DataFrame, top_n: int, out: Path) -> None:
    """Best-of-top-k success rate as k grows from 1 to top_n.

    Shows how much performance you recover by rescoring more poses, and how
    much headroom remains above the top-1 result.
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    col = "oracle_top_k_pct_rmsd_2.0A"
    for method, sub in rank_df.groupby("method"):
        sub = sub.sort_values("rank")
        n_pairs = int(sub["n_pairs_total"].iloc[0]) if "n_pairs_total" in sub else len(sub)
        ax.plot(sub["rank"], sub[col],
                marker="s", ls="--", lw=2,
                label=f"{TOOL_LABEL.get(method, method)} (n={n_pairs})",
                color=TOOL_COLORS.get(method))
    ax.set_xticks(range(1, top_n + 1))
    ax.set_xlabel("Top-k (best of first k ranked poses)")
    ax.set_ylabel("% of pairs where any of top-k has RMSD ≤ 2 Å")
    ax.set_title(_vt("Cumulative oracle success — best of top-k poses\n"
                     "(gap to 100 % = poses the method cannot produce within top-k)"))
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def _write_topn_within_report(out_path: Path, *, figure_name: str, csv_name: str,
                              title: str, panel_desc: str, notes_lines: list[str],
                              stats: dict | None, within_df: "pd.DataFrame | None",
                              top_n: int, success_line: float) -> None:
    """Write fig 18's moved-off-figure caption/notes AND the FULL complex-panel paired
    statistics to a companion text file (items 4 + 5 of the graph rework).

    Mirrors :func:`_write_pbvalid_yield_report`'s layout. The plot keeps only the visual
    Wilson-CI whiskers; every number behind them — per-threshold Cochran's Q, the full
    pairwise exact-McNemar (Holm-pooled) family for both the rank-1 and best-of-top-N
    curves, and the per-tool rates with 95 % Wilson CIs + ranking-headroom effect size —
    is tabulated here so the figure stays clean.
    """
    def _p(p) -> str:
        if p is None or p != p:
            return "n/a"
        if p <= 0:
            return "<1e-300"
        return f"{p:.3e}" if p < 1e-3 else f"{p:.4f}"

    def _pct(v) -> str:
        return "n/a" if v is None or v != v else f"{100 * float(v):5.1f}"

    W = 92
    L: list[str] = []
    L.append("=" * W)
    L.append(title)
    L.append("=" * W)
    L.append("")
    L.append(f"Figure : {figure_name}")
    L.append(f"Table  : {csv_name}")
    L.append(f"Source : per_pose_metrics.csv  (top-{top_n} ranked poses per tool)")
    L.append(f"Panel  : {panel_desc}")
    L.append("")
    L.append("-" * W)
    L.append("FIGURE NOTES  (moved off the plot)")
    L.append("-" * W)
    L.extend(notes_lines)
    L.append("")

    # Per-tool percentage table straight from the plotted curve, at the whole-Ångström
    # thresholds (a compact reproduction of what the figure's value labels showed).
    if within_df is not None and not within_df.empty:
        present = [m for m in ("autodock", "diffdock", "equibind")
                   if m in set(within_df["method"].unique())]
        ts_all = sorted(float(t) for t in within_df["rmsd_threshold_A"].unique())
        mark_ts = [t for t in ts_all if float(t).is_integer() and t >= 1.0] or ts_all[-4:]
        L.append("-" * W)
        L.append("TOP-RANKED ACCURACY  (% of each tool's complexes; rank-1 / best-of-top-N)")
        L.append("-" * W)
        hdr = f"{'threshold_A':>11s} " + " ".join(
            f"{TOOL_LABEL.get(m, m)[:13]:>15s}" for m in present)
        L.append(hdr + "    [rank-1 above · best-of-top-N below]")
        L.append("-" * len(hdr))
        for t in mark_ts:
            def _cell(m, col):
                row = within_df[(within_df["method"] == m)
                                & (within_df["rmsd_threshold_A"] == t)]
                return f"{float(row[col].iloc[0]):5.1f}" if not row.empty else "  n/a"
            L.append(f"{t:>11.3g} " + " ".join(f"{_cell(m, 'top1_within_%'):>15s}"
                                                for m in present))
            L.append(f"{'':>11s} " + " ".join(f"{_cell(m, 'best_topN_within_%'):>15s}"
                                                for m in present))
        L.append("")

    L.append("-" * W)
    L.append("STATISTICAL RESULTS  (complex-level paired tests)")
    L.append("-" * W)
    if not stats:
        L.append("No paired statistics available (fewer than two tools carried the")
        L.append("per-complex frames, or the stats step failed). See the figure's curves.")
        L.append("=" * W)
        out_path.write_text("\n".join(L) + "\n", encoding="utf-8")
        return

    L.append(f"Metric  : {stats.get('metric', '')}")
    L.append(f"Unit    : {stats.get('unit', '')}")
    L.append(f"Test    : {stats.get('test', '')}")
    ts = [float(t) for t in stats.get("test_thresholds_A", ())]
    L.append(f"Tested thresholds : {', '.join(f'{t:g} Å' for t in ts)}  "
             f"(emphasised success line = {success_line:g} Å)")
    L.append("Significance : *** p<.001   ** p<.01   * p<.05   ns = not significant.")
    L.append("")

    def _across_block(label, per_t_key, fam_key):
        per_t = {round(float(r["rmsd_threshold_A"]), 4): r
                 for r in stats.get(per_t_key, [])}
        fam = stats.get(fam_key, [])
        L.append(f"### {label} — across-tool HIT (paired on shared complexes)")
        for t in ts:
            rec = per_t.get(round(t, 4), {})
            om = rec.get("omnibus", {})
            n = rec.get("n_complete", "?")
            if "Q" in om:
                L.append(f"  @ {t:g} Å  (n={n}):  Cochran's Q({om.get('df', '?')}) "
                         f"= {om['Q']:.2f}, p = {_p(om.get('p'))}")
            else:
                L.append(f"  @ {t:g} Å  (n={n}):  {om.get('status', 'omnibus n/a')}")
            for pr in fam:
                if abs(float(pr["rmsd_threshold_A"]) - t) > 1e-6:
                    continue
                a = TOOL_LABEL.get(pr["a"], pr["a"]); b = TOOL_LABEL.get(pr["b"], pr["b"])
                L.append(f"      {a:>16s} vs {b:<16s}: "
                         f"wins {int(pr['a_wins'])}/{int(pr['b_wins'])}, "
                         f"p_raw = {_p(pr.get('p_raw'))}, p_holm = {_p(pr.get('p_holm'))}"
                         f"  {pr.get('star', '')}")
        L.append("")

    _across_block("Rank-1 (solid curve)", "rank1_across_tools", "rank1_pairwise_holm")
    _across_block("Best-of-top-N (dashed curve)", "best_topN_across_tools",
                  "best_topN_pairwise_holm")

    hr = stats.get("ranking_headroom", [])
    if hr:
        L.append("### Per-tool rates + 95% Wilson CIs, and ranking headroom "
                 "(best-of-top-N − rank-1)")
        L.append("  headroom is a NESTED effect size (rank-1 hits ⊂ best-of-top-N) — "
                 "reported as a share with a Wilson CI, not a p-value.")
        hdr = (f"  {'tool':>16s} {'thr_A':>6s} {'n':>4s} {'rank1%':>7s} "
               f"{'rank1 95%CI':>15s} {'bestN%':>7s} {'bestN 95%CI':>15s} "
               f"{'head_pp':>8s} {'head 95%CI':>15s}")
        L.append(hdr)
        L.append("  " + "-" * (len(hdr) - 2))
        for r in hr:
            t1c = r.get("top1_ci", [None, None]); btc = r.get("best_topN_ci", [None, None])
            hdc = r.get("headroom_ci", [None, None])
            L.append(f"  {TOOL_LABEL.get(r['method'], r['method'])[:16]:>16s} "
                     f"{float(r['rmsd_threshold_A']):>6g} {int(r['n']):>4d} "
                     f"{_pct(r['top1_rate']):>7s} "
                     f"[{_pct(t1c[0])},{_pct(t1c[1])}]".rjust(15) + " "
                     f"{_pct(r['best_topN_rate']):>7s} "
                     f"[{_pct(btc[0])},{_pct(btc[1])}]".rjust(15) + " "
                     f"{100 * float(r['headroom_share']):>8.1f} "
                     f"[{_pct(hdc[0])},{_pct(hdc[1])}]".rjust(15))
        L.append("")

    L.append("=" * W)
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _annotate_topn_within_stats(fig, ax, methods, payload, thr,
                                text_box: bool = True) -> None:
    """Overlay 95 % Wilson CI whiskers at the pre-specified thresholds and (only when
    ``text_box``) a Cochran's Q / pairwise-McNemar + ranking-headroom summary box on
    fig 18's complex panel (A). With ``text_box=False`` just the whiskers are drawn and
    the full summary is written to the companion report text file instead. Headline
    numbers are at the largest tested threshold; the full matrix lives in the sidecar."""
    if not payload:
        return
    ts = [float(t) for t in payload.get("test_thresholds_A", ())]
    if not ts:
        return

    def _rec_at(per_t, t):
        for r in per_t:
            if abs(float(r["rmsd_threshold_A"]) - t) < 1e-6:
                return r
        return None

    # per-tool, per-threshold rates + Wilson CIs (on each tool's own denominator).
    hr = {(r["method"], round(float(r["rmsd_threshold_A"]), 4)): r
          for r in payload.get("ranking_headroom", [])}
    for method in methods:
        c = TOOL_COLORS.get(method, "grey")
        for t in ts:
            r = hr.get((method, round(t, 4)))
            if not r:
                continue
            # rank-1 (solid) whisker left of the tick, best-of-top-N (dashed) right,
            # so a tool's two intervals never overprint.
            for col, ci, dx, a in (("top1_rate", "top1_ci", -0.035, 0.9),
                                    ("best_topN_rate", "best_topN_ci", 0.035, 0.55)):
                y = 100.0 * float(r[col]); lo, hi = (100.0 * float(v) for v in r[ci])
                ax.errorbar(t + dx, y, yerr=[[max(0.0, y - lo)], [max(0.0, hi - y)]],
                            fmt="none", ecolor=c, elinewidth=1.3, capsize=2.5,
                            alpha=a, zorder=5)

    if not text_box:
        return
    t_hi = max(ts)
    rec = _rec_at(payload.get("rank1_across_tools", []), t_hi) or {}
    n = rec.get("n_complete", "?")
    lines = [f"Paired tests · n={n} complexes · @ {t_hi:g} Å (1 Å in sidecar):"]
    om = rec.get("omnibus", {})
    if "Q" in om:
        lines.append(f"  rank-1 across tools: Cochran Q={om['Q']:.1f}, {su.fmt_p(om['p'])}")
    for pr in payload.get("rank1_pairwise_holm", []):
        if abs(float(pr["rmsd_threshold_A"]) - t_hi) > 1e-6:
            continue
        a = TOOL_LABEL.get(pr["a"], pr["a"]); b = TOOL_LABEL.get(pr["b"], pr["b"])
        lines.append(f"    {a} vs {b}: {su.p_stars(pr.get('p_holm'))} "
                     f"({su.fmt_p(pr.get('p_holm'))})")
    hs = [r for r in payload.get("ranking_headroom", [])
          if abs(float(r["rmsd_threshold_A"]) - t_hi) < 1e-6]
    if hs:
        parts = ", ".join(f"{TOOL_LABEL.get(r['method'], r['method'])} "
                          f"+{100 * float(r['headroom_share']):.1f}" for r in hs)
        lines.append(f"  ranking headroom (best−rank-1, pp): {parts}")
    lines.append("bars = 95% Wilson CI (rank-1 left · best-of-top-N right)")
    ax.text(0.02, 0.975, "\n".join(lines), transform=ax.transAxes, ha="left", va="top",
            fontsize=6.9, color="0.20", linespacing=1.35,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9), zorder=7)


def plot_topn_within_thresholds(within_df: pd.DataFrame, top_n: int,
                                thresholds: tuple[float, ...], out: Path,
                                pb_valid_only: bool = False,
                                stats: dict | None = None,
                                x_label: str = "RMSD threshold vs crystal ligand (Å)",
                                metric_note: str = ("RMSD = symmetry-corrected "
                                                    "heavy-atom, no superposition."),
                                success_line: float = 2.0,
                                csv_name: str | None = None) -> None:
    """How many of the top-ranked poses land within a fine RMSD grid.

    AutoDock Vina / DiffDock (native rank) + EquiBind (gnina-affinity ranked), two
    panels sharing the x-axis (distance threshold t ∈ {1, 1.25, 1.5, …} Å):
      (A) complex level — % of complexes whose RANK-1 pose is within t (solid) and
          % with ANY of the first ``top_n`` ranked poses within t (best-of-top-N,
          dashed). The vertical gap between a tool's two curves is the accuracy its
          ranking leaves on the table below that threshold.
      (B) pose level — % of ALL the tool's pooled top-N ranked poses within t.
    An emphasised dotted vertical line marks ``success_line`` Å (2 Å canonical docking
    success; 1 Å for the Kabsch/form-fidelity view).

    Layout (per the graph rework): the LEGEND sits between the title and the plot; each
    curve's complex count ``n`` is stated ONCE in the title, not per legend entry; and
    the full caption/notes + the complete paired statistics are written to a companion
    ``*_report.txt`` next to each PNG rather than printed on the figure — only the visual
    95 % Wilson-CI whiskers remain on panel A.
    """
    if within_df is None or within_df.empty:
        return
    present = set(within_df["method"].unique())
    methods = [m for m in ("autodock", "diffdock", "equibind") if m in present]
    if not methods:
        return
    thr = [float(t) for t in thresholds]

    def _label(method):
        lab = TOOL_LABEL.get(method, method)
        return lab + "*" if method == "equibind" else lab

    def _fmt(ax):
        ax.axvline(success_line, color="grey", ls=":", alpha=0.7, lw=1.2)
        ax.set_xticks(thr)
        ax.set_xticklabels([f"{t:g}" for t in thr], rotation=45, ha="right")
        ax.set_xlabel(x_label)
        ax.set_ylim(0, 100); ax.grid(alpha=0.3)

    rank_note = ("Rank: AutoDock/DiffDock native; EquiBind has none, so it is ranked by "
                 f"gnina affinity. Emphasised dotted line = {success_line:g} Å. " + metric_note)
    valid_tag  = " · PoseBusters-valid poses only" if pb_valid_only else ""
    pose_word  = "PB-valid pose" if pb_valid_only else "pose"
    title_word = "valid-pose" if pb_valid_only else "pose"
    valid_note = ("A pose counts only if it is within t Å AND PoseBusters-valid."
                  if pb_valid_only else "A pose counts if it is within t Å (no validity gate).")

    # complex count (item 4b: moved from the legend labels into the title line)
    n_pairs_by_m = {m: int(within_df[within_df["method"] == m]["n_pairs"].iloc[0])
                    for m in methods}
    uniq_n = sorted(set(n_pairs_by_m.values()))
    n_note = (f"n = {uniq_n[0]} complexes" if len(uniq_n) == 1
              else "n complexes: " + ", ".join(f"{TOOL_LABEL.get(m, m)} {n_pairs_by_m[m]}"
                                                for m in methods))

    def _finish_top_legend(fig, handles, labels, ncol, subtitle):
        """Legend BETWEEN the title and the axes (item 3): a horizontal legend anchored
        just under the (figure) title, above the plot."""
        fig.subplots_adjust(top=0.80, bottom=0.13, left=0.11, right=0.965)
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.895),
                   ncol=ncol, fontsize=8, framealpha=0.9, columnspacing=1.3,
                   handlelength=2.2, borderaxespad=0.2)
        fig.suptitle(_vt(subtitle), fontsize=12, fontweight="bold", y=0.995)

    # ── Graph 1 (complex level): rank-1 (solid) + best-of-top-N (dashed) ──
    figA, axA = plt.subplots(figsize=(8.6, 7.0))
    # Handles interleaved per tool (rank-1 then best-of-top-N): matplotlib fills the
    # horizontal legend column-major, so with ncol = #tools each column becomes one tool
    # (rank-1 on top, best-of-top-N below) instead of the two curve types scrambling.
    legA_h, legA_l = [], []
    for method in methods:
        sub = within_df[within_df["method"] == method].sort_values("rmsd_threshold_A")
        if sub.empty:
            continue
        color, label = TOOL_COLORS.get(method, "grey"), _label(method)
        (ln1,) = axA.plot(sub["rmsd_threshold_A"], sub["top1_within_%"], marker="o",
                          lw=2, color=color)
        (ln2,) = axA.plot(sub["rmsd_threshold_A"], sub["best_topN_within_%"], marker="s",
                          ls="--", lw=2, color=color, alpha=0.75)
        legA_h += [ln1, ln2]
        legA_l += [f"{label} — rank-1", f"{label} — best of top-{top_n}"]
    # Value labels at every whole-Ångström threshold: each curve's % where it crosses
    # that distance, in white boxes de-collided vertically so overlapping intercepts stay
    # legible. Interior marks label to the RIGHT of a reference line; the rightmost mark
    # labels to the LEFT to stay on-axis.
    mark_ts = [t for t in thr if float(t).is_integer() and t >= 1.0]

    def _mark_labels(t, side):
        dx, ha = (0.07, "left") if side == "right" else (-0.07, "right")
        marks = []
        for method in methods:
            row = within_df[(within_df["method"] == method)
                            & (within_df["rmsd_threshold_A"] == t)]
            if row.empty:
                continue
            c = TOOL_COLORS.get(method, "grey")
            for col in ("top1_within_%", "best_topN_within_%"):
                marks.append((float(row[col].iloc[0]), c))
        marks.sort(key=lambda m: m[0])
        _gap, _prev, placed = 5.0, None, []
        for _y, _c in marks:
            _ly = _y if _prev is None or _y - _prev >= _gap else _prev + _gap
            _prev = _ly
            placed.append([_ly, _y, _c])
        # Keep the de-collided stack inside the axes: when curves saturate near 100 %
        # (e.g. the Kabsch/non-PB-valid views) the upward push can drive the top label
        # past the axes ceiling and into the top legend band — slide the whole cluster
        # down by the overflow so every callout stays on-axis.
        _ceiling = 98.0
        _overflow = (max(p[0] for p in placed) - _ceiling) if placed else 0.0
        if _overflow > 0:
            for _p in placed:
                _p[0] -= _overflow
        for _ly, _y, _c in placed:
            axA.annotate(f"{_y:.1f}", xy=(t, _y), xytext=(t + dx, _ly),
                         va="center", ha=ha, fontsize=8, fontweight="bold", color=_c,
                         bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=_c,
                                   lw=0.8, alpha=0.9), zorder=6)

    for t in mark_ts:
        if abs(t - success_line) > 1e-9:   # success line already drawn emphasised by _fmt
            axA.axvline(t, color="0.75", ls=":", alpha=0.6, lw=1.0, zorder=1)
        _mark_labels(t, "left" if t == mark_ts[-1] else "right")
    _fmt(axA)
    axA.set_xlim(0, thr[-1] + (thr[-1] - thr[0]) * 0.04)   # x-axis starts at 0
    axA.set_xticks([0] + thr)
    axA.set_xticklabels(["0"] + [f"{t:g}" for t in thr], rotation=45, ha="right")
    axA.set_ylabel(f"% of complexes with a {pose_word} within the threshold")
    # Paired-test Wilson-CI whiskers only (the summary text goes to the report file).
    if stats:
        try:
            _annotate_topn_within_stats(figA, axA, methods, stats, thr, text_box=False)
        except Exception as _exc:
            print(f"  WARNING: fig 18 complex-panel whiskers failed ({_exc})")
    _finish_top_legend(
        figA, legA_h, legA_l, len(methods),
        f"Top-ranked {title_word} accuracy vs distance threshold  ·  {n_note}\n"
        f"(rank-1 solid · best-of-top-N dashed{valid_tag})")
    figA.savefig(out.with_stem(out.stem + "_complex"), dpi=160, bbox_inches="tight")
    plt.close(figA)

    # Companion report for the complex panel (items 4 + 5): full caption + statistics.
    notes_A = [
        "Denominator = each tool's complexes (n pairs).",
        "best-of-top-N = the closest of the tool's first N ranked poses per complex.",
        valid_note,
        rank_note,
    ]
    try:
        _write_topn_within_report(
            out.parent / (out.stem + "_complex_report.txt"),
            figure_name=out.with_stem(out.stem + "_complex").name,
            csv_name=csv_name or "(see topn_within_thresholds*.csv)",
            title=f"Fig 18 — top-ranked {title_word} accuracy within distance thresholds"
                  f" (complex level){valid_tag}",
            panel_desc="complex level — % of complexes with rank-1 (solid) and "
                       "best-of-top-N (dashed) pose within t",
            notes_lines=notes_A, stats=stats, within_df=within_df, top_n=top_n,
            success_line=success_line)
    except Exception as _exc:
        print(f"  WARNING: fig 18 complex report write failed ({_exc})")

    # ── Graph 2 (pose level): % of the pooled top-N poses within t ──
    figB, axB = plt.subplots(figsize=(8.6, 7.0))
    pose_counts = {m: int(within_df[within_df["method"] == m]["n_poses_topN"].iloc[0])
                   for m in methods}
    b_h, b_l = [], []
    for method in methods:
        sub = within_df[within_df["method"] == method].sort_values("rmsd_threshold_A")
        if sub.empty:
            continue
        color, label = TOOL_COLORS.get(method, "grey"), _label(method)
        (ln,) = axB.plot(sub["rmsd_threshold_A"], sub["pose_within_%"], marker="o", lw=2,
                         color=color)
        b_h.append(ln); b_l.append(label)
    _fmt(axB)
    axB.set_ylabel("% of the top-N ranked poses within the threshold"
                   + (" & PB-valid" if pb_valid_only else ""))
    pose_note = "; ".join(f"{TOOL_LABEL.get(m, m)} ≤{pose_counts[m]}" for m in methods)
    _finish_top_legend(
        figB, b_h, b_l, len(methods),
        f"How many of the top-{top_n} ranked poses are within t Å"
        + (" AND PB-valid" if pb_valid_only else "")
        + f"\n(pooled poses per tool: {pose_note})")
    figB.savefig(out.with_stem(out.stem + "_pose"), dpi=160, bbox_inches="tight")
    plt.close(figB)

    notes_B = [
        "Denominator = each tool's pooled rank ≤ N poses.",
        valid_note,
        rank_note,
    ]
    try:
        _write_topn_within_report(
            out.parent / (out.stem + "_pose_report.txt"),
            figure_name=out.with_stem(out.stem + "_pose").name,
            csv_name=csv_name or "(see topn_within_thresholds*.csv)",
            title=f"Fig 18 — pooled top-{top_n} pose accuracy within distance thresholds"
                  f" (pose level){valid_tag}",
            panel_desc="pose level — % of ALL pooled rank ≤ N poses within t",
            notes_lines=notes_B, stats=None, within_df=within_df, top_n=top_n,
            success_line=success_line)
    except Exception as _exc:
        print(f"  WARNING: fig 18 pose report write failed ({_exc})")


# Item 6 of the fig-18 rework: how the raw→gnina-refined benefit evolves as the top-N
# pool deepens. Each tool's raw variant is paired with its gnina-optimised variant on
# the SAME metric fig 18 uses (best-of-top-N pose within the success-line threshold,
# PB-valid-gated), swept over pool depth N = 1, 2, 3 … top_n. At the deepest N the pool
# is the whole set, so the gap there isolates the geometry/validity gain of refinement;
# at shallow N it also folds in the ranking gnina hands EquiBind (which has no native rank).
_REFINE_GAIN_PAIRS = [
    # (display, raw method key, refined method key, raw rank kind, refined rank kind)
    # AutoDock has NO ML-refinement step (raw key == refined key): it appears as a single
    # "native" best-of-top-N curve, and because rank depth is its only lever, the rise of
    # that curve with N is the pure rank-depth gain — the reference the refined tools are
    # read against.
    ("AutoDock Vina", "autodock",             "autodock",                "native",     "native"),
    ("DiffDock",      "diffdock",             "diffdock_gnina",          "native",     "native"),
    ("EquiBind",      "equibind_unguided_raw", "equibind_unguided_gnina", "generation", "gnina"),
]


def _rank_series_for(sub: pd.DataFrame, kind: str):
    """(possibly reordered frame, per-complex 1-based rank) under a ranking basis:
    ``native`` = DiffDock confidence rank (duplicate top pose de-duped);
    ``generation`` = EquiBind generation-order position;
    ``gnina`` = gnina-affinity rank (most-negative = 1)."""
    sub = sub.copy()
    if kind == "native":
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        sub = _dedup_ranked_poses(sub)
        return sub, sub["rank"]
    if kind == "generation":
        sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
        sub = sub.sort_values(["protein", "ligand", "_gi"], kind="mergesort")
        sub["_sr"] = sub.groupby(["protein", "ligand"]).cumcount() + 1
        return sub, sub["_sr"]
    return sub, _gnina_affinity_rank(sub)      # gnina


def aggregate_refinement_gain_vs_depth(
        df_full: pd.DataFrame, thr: float, depths,
        rmsd_col: str = "rmsd", pb_valid_only: bool = True) -> pd.DataFrame:
    """Per (tool, pool depth N): best-of-top-N success for the raw vs the gnina-refined
    variant and their difference (pp).

    A complex is a HIT at depth N if the MINIMUM ``rmsd_col`` among its first N ranked
    poses — PB-valid-gated when ``pb_valid_only`` — is ≤ ``thr``. The denominator is the
    complexes both variants share (paired), matching fig 18's "all of the tool's
    complexes" (a complex with zero qualifying poses counts as a miss). AutoDock Vina has
    no refinement step (raw key == refined key), so its raw and refined curves coincide
    (``gain_pp`` = 0, ``has_refinement`` = False) and it serves as the native rank-depth
    reference. Needs the FULL frame. Returns empty if no pair is available."""
    present = set(df_full["method"].astype(str))
    rows = []
    for disp, raw_key, ref_key, raw_kind, ref_kind in _REFINE_GAIN_PAIRS:
        if raw_key not in present or ref_key not in present:
            continue
        raw_sub, raw_rank = _rank_series_for(
            df_full[df_full["method"].astype(str) == raw_key], raw_kind)
        ref_sub, ref_rank = _rank_series_for(
            df_full[df_full["method"].astype(str) == ref_key], ref_kind)

        def _prep(sub, rank_series):
            sub = sub.copy()
            sub["_rk"] = pd.to_numeric(rank_series, errors="coerce")
            if pb_valid_only:
                sub = sub[_to_bool(sub["pb_valid"])]
            return sub.dropna(subset=[rmsd_col, "_rk"])

        raw_g, ref_g = _prep(raw_sub, raw_rank), _prep(ref_sub, ref_rank)
        # paired denominator = complexes each variant has ANY pose for (pre-gate).
        raw_ids = set(map(tuple, df_full[df_full["method"].astype(str) == raw_key]
                          [["protein", "ligand"]].drop_duplicates().to_numpy()))
        ref_ids = set(map(tuple, df_full[df_full["method"].astype(str) == ref_key]
                          [["protein", "ligand"]].drop_duplicates().to_numpy()))
        ids = raw_ids & ref_ids
        n_pairs = len(ids)
        if not n_pairs:
            continue

        def _hits(g, N):
            gg = g[g["_rk"] <= N]
            if gg.empty:
                return 0
            best = gg.groupby(["protein", "ligand"])[rmsd_col].min()
            best = best[[idx in ids for idx in best.index]]
            return int((best <= thr).sum())

        for N in depths:
            rh, fh = _hits(raw_g, N), _hits(ref_g, N)
            rows.append({
                "tool": disp, "depth": int(N), "n_pairs": n_pairs,
                "has_refinement": bool(raw_key != ref_key),
                "raw_success_%": round(100.0 * rh / n_pairs, PCT_DECIMALS),
                "refined_success_%": round(100.0 * fh / n_pairs, PCT_DECIMALS),
                "gain_pp": round(100.0 * (fh - rh) / n_pairs, PCT_DECIMALS),
                "raw_hit": rh, "refined_hit": fh,
            })
    return pd.DataFrame(rows)


def plot_refinement_gain_vs_depth(rec: pd.DataFrame, out: Path, top_n: int, thr: float,
                                  metric_label: str, pb_valid_only: bool = True,
                                  report_path: "Path | None" = None) -> None:
    """Companion to fig 18: the raw→gnina-refined benefit as the top-N pool deepens.

    Top panel — each tool's best-of-top-N success (% complexes within ``thr`` Å, on the
    parent figure's metric) for the raw (dashed) and gnina-refined (solid) variant, the
    band between them shaded. AutoDock Vina has no refinement step, so it is drawn as a
    single native curve (no dashed/band): its rise with N is the pure rank-depth gain,
    the reference the refined tools are read against. Bottom panel — the refinement
    difference (gnina-refined − raw, pp) for the refined tools (AutoDock has none). Legend
    between title and plot; the caption + per-depth table go to ``report_path``.
    """
    if rec is None or rec.empty:
        return
    tools = [t for t in ("AutoDock Vina", "DiffDock", "EquiBind") if t in set(rec["tool"])]
    if not tools:
        return
    colors = {"AutoDock Vina": TOOL_COLORS.get("autodock", "tab:blue"),
              "DiffDock": TOOL_COLORS.get("diffdock", "tab:orange"),
              "EquiBind": TOOL_COLORS.get("equibind", "tab:green")}
    depths = sorted(int(d) for d in rec["depth"].unique())
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(8.8, 8.6), sharex=True)
    handles, labels = [], []
    for tool in tools:
        s = rec[rec["tool"] == tool].sort_values("depth")
        c = colors.get(tool, "grey")
        has_ref = bool(s["has_refinement"].iloc[0]) if "has_refinement" in s.columns else True
        # Refined line (== the native line when the tool has no refinement step).
        (lref,) = axA.plot(s["depth"], s["refined_success_%"], marker="o", lw=2, color=c)
        if has_ref:
            (lraw,) = axA.plot(s["depth"], s["raw_success_%"], marker="s", ls="--", lw=2,
                               color=c, alpha=0.7)
            axA.fill_between(s["depth"], s["raw_success_%"], s["refined_success_%"],
                             color=c, alpha=0.12, zorder=0)
            axB.plot(s["depth"], s["gain_pp"], marker="o", lw=2, color=c)   # refinement Δ
            handles += [lref, lraw]
            labels += [f"{tool} — gnina-refined", f"{tool} — raw"]
        else:
            # AutoDock: no refinement — a single native curve whose rise with N is the
            # pure rank-depth gain (the reference the refined tools are read against).
            handles += [lref]
            labels += [f"{tool} — native (no refinement)"]
    axA.set_ylabel("% complexes with a best-of-top-N pose\n"
                   f"within {thr:g} Å ({metric_label})"
                   + (" & PB-valid" if pb_valid_only else ""))
    axA.set_ylim(0, max(5.0, float(rec[["raw_success_%", "refined_success_%"]]
                                   .to_numpy().max()) * 1.12))
    axA.grid(alpha=0.3)
    axB.axhline(0, color="0.6", lw=1.0)
    axB.set_ylabel("Refinement benefit\n(gnina-refined − raw, pp)")
    axB.set_xlabel("Top-N pool depth  (number of ranked poses considered, N = 1 … "
                   f"{top_n})")
    axB.grid(alpha=0.3)
    axB.set_xticks(depths)
    axB.set_xlim(depths[0] - 0.4, depths[-1] + 0.4)
    n_note = (f"n = {int(rec['n_pairs'].iloc[0])} complexes"
              if rec["n_pairs"].nunique() == 1
              else "n per tool: " + ", ".join(
                  f"{t} {int(rec[rec['tool'] == t]['n_pairs'].iloc[0])}" for t in tools))
    fig.subplots_adjust(top=0.85, bottom=0.09, left=0.12, right=0.965, hspace=0.13)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935),
               ncol=len(handles), fontsize=8, framealpha=0.9, columnspacing=1.2,
               handlelength=2.0, borderaxespad=0.2)
    fig.suptitle(_vt("Best-of-top-N success vs pool depth: refinement gain + rank-depth "
                     f"gain  ·  {n_note}\n(gnina-refined solid · raw dashed · AutoDock "
                     f"Vina native reference; success = within {thr:g} Å"
                     + (" & PB-valid" if pb_valid_only else "") + ")"),
                 fontsize=11.5, fontweight="bold", y=0.997)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)

    if report_path is not None:
        try:
            _write_refinement_gain_report(report_path, rec, out.name, thr, metric_label,
                                          pb_valid_only, top_n)
        except Exception as _exc:
            print(f"  WARNING: refinement-gain report write failed ({_exc})")


def _write_refinement_gain_report(out_path: Path, rec: pd.DataFrame, figure_name: str,
                                  thr: float, metric_label: str, pb_valid_only: bool,
                                  top_n: int) -> None:
    """Companion text for the raw→refined gain-vs-depth figure: caption + per-depth table."""
    W = 84
    L = ["=" * W,
         "Fig 18 companion — pose-refinement benefit as the top-N pool deepens",
         "=" * W, "",
         f"Figure : {figure_name}",
         f"Metric : best-of-top-N pose within {thr:g} Å ({metric_label})"
         + (" AND PoseBusters-valid" if pb_valid_only else ""),
         f"Depth  : N = 1 … {top_n} ranked poses (continuous sweep)",
         "",
         "-" * W, "FIGURE NOTES", "-" * W,
         "Each tool's raw variant is paired with its gnina-refined variant on the SAME",
         "complexes. A complex is a HIT at depth N when the closest of its first N ranked",
         "poses is within the threshold (PB-valid-gated when stated). Ranking basis:",
         "  DiffDock — native confidence rank for BOTH raw and gnina (refinement keeps it);",
         "  EquiBind — raw uses generation order; gnina uses gnina-affinity rank (EquiBind",
         "             has no native rank, so gnina supplies one — hence at shallow N the",
         "             benefit folds in ranking, while at the deepest N — the whole pool —",
         "             the gap is the pure geometry/validity gain of refinement.",
         "  AutoDock Vina — NO refinement step (raw == native): shown as a single native",
         "             curve. It has no refinement benefit; its rise in success as N grows",
         "             is the pure RANK-DEPTH gain (best-of-top-N − rank-1), the reference",
         "             the refined tools are read against.",
         "gain_pp = refined_success_% − raw_success_% (percentage points; 0 for AutoDock).",
         ""]
    for tool in [t for t in ("AutoDock Vina", "DiffDock", "EquiBind")
                 if t in set(rec["tool"])]:
        s = rec[rec["tool"] == tool].sort_values("depth")
        n = int(s["n_pairs"].iloc[0])
        has_ref = bool(s["has_refinement"].iloc[0]) if "has_refinement" in s.columns else True
        L.append("-" * W)
        if has_ref:
            L.append(f"{tool}   (n = {n} complexes)")
            L.append("-" * W)
            L.append(f"  {'depth N':>7s} {'raw %':>8s} {'refined %':>10s} {'gain (pp)':>10s}")
            for _, r in s.iterrows():
                L.append(f"  {int(r['depth']):>7d} {r['raw_success_%']:>8.1f} "
                         f"{r['refined_success_%']:>10.1f} {r['gain_pp']:>+10.1f}")
        else:
            top1 = float(s.iloc[0]["refined_success_%"])
            L.append(f"{tool}   (n = {n} complexes; native, no refinement)")
            L.append("-" * W)
            L.append(f"  {'depth N':>7s} {'success %':>10s} {'rank-depth gain vs N=1 (pp)':>28s}")
            for _, r in s.iterrows():
                succ = float(r["refined_success_%"])
                L.append(f"  {int(r['depth']):>7d} {succ:>10.1f} {succ - top1:>+28.1f}")
        L.append("")
    L.append("=" * W)
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")


# Multi-depth version of fig 18's complex panel: instead of rank-1 + best-of-top-N, show
# the best-of-top-N curve at SEVERAL pool depths (top-1 = rank-1, top-15, top-30) so the
# accuracy gained purely by deepening the ranked pool reads directly. Colour = tool,
# line style = pool depth. Per request the value-number callouts are omitted (the numbers
# live in the companion CSV + report text instead).
_MULTIDEPTH_POOL = (1, 15, 30)
_MULTIDEPTH_STYLE = {1: ("-", "o"), 15: ("--", "s"), 30: (":", "^")}   # depth → (ls, marker)


def aggregate_within_thresholds_by_depth(
        df: pd.DataFrame, depths, thresholds, eq_df: "pd.DataFrame | None" = None,
        pb_valid_only: bool = False, rmsd_col: str = "rmsd") -> pd.DataFrame:
    """Best-of-top-N within-threshold success at SEVERAL pool depths at once.

    For each depth d in ``depths`` this reuses :func:`aggregate_topn_within_thresholds`
    (top_n = d) and keeps its ``best_topN_within_%`` — the % of the tool's complexes whose
    CLOSEST of its first d ranked poses (PB-valid-gated when set) is within t. Depth 1 is
    the rank-1 curve. Returns long-form (method, depth, rmsd_threshold_A, within_%, n_pairs)."""
    frames = []
    for d in depths:
        a = aggregate_topn_within_thresholds(
            df, int(d), thresholds, eq_df=eq_df,
            pb_valid_only=pb_valid_only, rmsd_col=rmsd_col)
        if a.empty:
            continue
        a = a[["method", "rmsd_threshold_A", "n_pairs", "best_topN_within_%"]].copy()
        a = a.rename(columns={"best_topN_within_%": "within_%"})
        a["depth"] = int(d)
        frames.append(a)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def plot_within_thresholds_by_depth(
        within_df: pd.DataFrame, depths, thresholds, out: Path,
        pb_valid_only: bool = False,
        x_label: "str | None" = None,
        metric_note: str = "RMSD = symmetry-corrected heavy-atom, no superposition.",
        success_line: float = 2.0, csv_name: "str | None" = None,
        report_path: "Path | None" = None, stats: "dict | None" = None,
        data_table: "dict | None" = None) -> None:
    """Fig 18 multi-depth complex view: best-of-top-N accuracy vs distance threshold at
    pool depths top-1 / top-15 / top-30. Colour = tool, line style = depth. The value at
    the emphasised threshold (``success_line``: 2 Å RMSD / 1 Å Kabsch) is called out on
    each curve; the rest of the sweep's numbers go to the companion CSV + report text.
    Legend between title and plot: a tool colour key + a pool-depth line-style key.

    When *data_table* is given (``{"depths": [...], "metrics": [{"label", "by_tool"}]}``)
    a compact table of the at-threshold values for BOTH metrics (RMSD ≤ 2 Å and Kabsch
    < 1 Å) is drawn beneath the plot, so either figure carries both graphs' numbers."""
    if within_df is None or within_df.empty:
        return
    if x_label is None:
        x_label = f"RMSD threshold vs {_ref_ligand_phrase()} (Å)"
    from matplotlib.lines import Line2D
    present = set(within_df["method"].unique())
    methods = [m for m in ("autodock", "diffdock", "equibind") if m in present]
    if not methods:
        return
    thr = [float(t) for t in thresholds]
    depths = [int(d) for d in depths]
    _has_tab = bool(data_table and data_table.get("metrics"))

    def _style(d):
        return _MULTIDEPTH_STYLE.get(int(d), ("-", "o"))

    fig, ax = plt.subplots(figsize=(8.8, 9.7 if _has_tab else 7.0))
    for method in methods:
        c = TOOL_COLORS.get(method, "grey")
        for d in depths:
            s = within_df[(within_df["method"] == method)
                          & (within_df["depth"] == d)].sort_values("rmsd_threshold_A")
            if s.empty:
                continue
            ls, mk = _style(d)
            ax.plot(s["rmsd_threshold_A"], s["within_%"], ls=ls, marker=mk, ms=4,
                    lw=2, color=c, alpha=0.9)
    ax.axvline(success_line, color="grey", ls=":", alpha=0.7, lw=1.2)

    # Value callouts AT the emphasised threshold (per request): mark where each curve
    # crosses the success line and print its % there, colour = tool. Up to 9 labels are
    # de-collided vertically so they stay legible; a thin leader ties a nudged label
    # back to its true point.
    _pts = []
    for method in methods:
        c = TOOL_COLORS.get(method, "grey")
        for d in depths:
            s = within_df[(within_df["method"] == method) & (within_df["depth"] == d)
                          & (np.isclose(within_df["rmsd_threshold_A"], success_line,
                                        atol=1e-6))]
            if s.empty:
                continue
            y = float(s["within_%"].iloc[0])
            ax.plot([success_line], [y], marker=_style(d)[1], ms=8, color=c,
                    markeredgecolor="black", markeredgewidth=0.6, zorder=6)
            _pts.append((y, c, y))
    _pts.sort(key=lambda p: p[0])
    _last, _placed = -1e9, []
    for _y, _c, _val in _pts:                       # push overlapping labels apart
        _yy = max(_y, _last + 3.6)
        _placed.append((_yy, _c, _val))
        _last = _yy
    _xt = success_line + max(thr) * 0.012
    for _yy, _c, _val in _placed:
        ax.annotate(f"{_val:.1f}", xy=(success_line, _val), xytext=(_xt, _yy),
                    fontsize=7.5, fontweight="bold", color=_c, va="center", ha="left",
                    arrowprops=(dict(arrowstyle="-", color=_c, lw=0.5, alpha=0.6)
                                if abs(_yy - _val) > 1.0 else None), zorder=7)

    ax.set_xlim(0, thr[-1] + (thr[-1] - thr[0]) * 0.04)
    ax.set_xticks([0] + thr)
    ax.set_xticklabels(["0"] + [f"{t:g}" for t in thr], rotation=45, ha="right")
    ax.set_ylim(0, 100); ax.grid(alpha=0.3)
    ax.set_xlabel(x_label)
    pose_word = "PB-valid pose" if pb_valid_only else "pose"
    ax.set_ylabel(f"% of complexes with a {pose_word} within the threshold")

    n_by_m = {m: int(within_df[within_df["method"] == m]["n_pairs"].iloc[0])
              for m in methods}
    uniq = sorted(set(n_by_m.values()))
    n_note = (f"n = {uniq[0]} complexes" if len(uniq) == 1
              else "n: " + ", ".join(f"{TOOL_LABEL.get(m, m)} {n_by_m[m]}" for m in methods))

    def _tool_label(m):
        # One uniform rule for BOTH graphs: the two ML tools are each shown as a single
        # representative variant (DiffDock = the selected/collapsed variant; EquiBind =
        # the gnina-ranked unguided one) and carry a trailing '*'; AutoDock has no
        # variant, so none. Guarded with endswith so an already-starred label (e.g.
        # 'DiffDock*' from the collapse override) is never double-marked. This keeps the
        # RMSD and Kabsch depth figures' legends identical in every run mode.
        lab = TOOL_LABEL.get(m, m)
        if m in ("diffdock", "equibind") and not lab.endswith("*"):
            lab += "*"
        return lab

    tool_h = [Line2D([0], [0], color=TOOL_COLORS.get(m, "grey"), lw=2.6, marker="o", ms=4)
              for m in methods]
    depth_h = [Line2D([0], [0], color="0.35", lw=2.0, ls=_style(d)[0], marker=_style(d)[1],
                      ms=4) for d in depths]
    depth_l = [f"best of top-{d}" + (" (rank-1)" if d == 1 else "") for d in depths]
    # Two horizontal keys stacked between the title and the plot: tools (colour) above,
    # pool depths (line style) below. Reserve extra room at the bottom for the data table.
    fig.subplots_adjust(top=0.80, bottom=(0.335 if _has_tab else 0.14),
                        left=0.11, right=0.965)
    fig.legend(tool_h, [_tool_label(m) for m in methods], loc="upper center",
               bbox_to_anchor=(0.5, 0.915), ncol=len(methods), fontsize=8,
               framealpha=0.9, columnspacing=1.3, handlelength=2.0, borderaxespad=0.2)
    fig.legend(depth_h, depth_l, loc="upper center", bbox_to_anchor=(0.5, 0.865),
               ncol=len(depths), fontsize=8, framealpha=0.9, columnspacing=1.3,
               handlelength=2.6, borderaxespad=0.2)
    fig.suptitle(_vt("Best-of-top-N pose accuracy vs distance threshold  ·  "
                     f"{n_note}"),
                 fontsize=12, fontweight="bold", y=0.99)

    # ── Data table beneath the plot: at-threshold values for BOTH metrics ──
    # Same table on the RMSD and Kabsch figures, so either is self-contained. Rows are
    # grouped by tool (variant label colour-coded), one row per metric; columns = depths.
    if _has_tab:
        tdepths = [int(d) for d in data_table["depths"]]
        col_labels = ["Docking variant", "Metric"] + [f"top-{d}" for d in tdepths]
        cell_text, row_tools = [], []
        for tool in methods:                       # same tool order as the legend
            for met in data_table["metrics"]:
                vals = met.get("by_tool", {}).get(tool)
                if vals is None:
                    continue
                cell_text.append([_tool_label(tool), met["label"]]
                                 + [("—" if v != v else f"{v:.1f}") for v in vals])
                row_tools.append(tool)
        if cell_text:
            tab_ax = fig.add_axes([0.08, 0.045, 0.90, 0.215])
            tab_ax.axis("off")
            tbl = tab_ax.table(cellText=cell_text, colLabels=col_labels,
                               loc="center", cellLoc="center")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(8.5)
            tbl.scale(1, 1.45)
            ncol = len(col_labels)
            for (r_i, c_i), cell in tbl.get_celld().items():
                cell.set_edgecolor("0.8")
                if r_i == 0:                        # header row
                    cell.set_text_props(fontweight="bold")
                    cell.set_facecolor("#eef0f2")
                else:
                    tool = row_tools[r_i - 1]
                    if c_i == 0:                    # variant column: tool colour + bold
                        cell.get_text().set_color(TOOL_COLORS.get(tool, "black"))
                        cell.set_text_props(fontweight="bold")
                    if c_i <= 1:
                        cell.get_text().set_ha("left")
                        cell.PAD = 0.04
                    # zebra-shade each tool's block
                    if methods.index(tool) % 2 == 1:
                        cell.set_facecolor("#f7f8f9")
            fig.text(0.5, 0.028,
                     "% of complexes with a PoseBusters-valid pose within the threshold "
                     "(RMSD ≤ 2 Å = placement · Kabsch < 1 Å = internal form)",
                     ha="center", va="center", fontsize=7.5, color="0.35")

    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)

    if csv_name is not None:
        _write_csv(within_df, out.parent / csv_name, index=False)
    if report_path is not None:
        try:
            _write_within_by_depth_report(report_path, within_df, out.name, depths,
                                          methods, success_line, metric_note,
                                          pb_valid_only, stats=stats)
        except Exception as _exc:
            print(f"  WARNING: fig 18 multi-depth report write failed ({_exc})")


def _write_within_by_depth_report(out_path: Path, within_df: pd.DataFrame, figure_name: str,
                                  depths, methods, success_line: float, metric_note: str,
                                  pb_valid_only: bool, stats: "dict | None" = None) -> None:
    """Companion text for the multi-depth fig: the per-depth per-tool accuracy TABLE the
    plot omits at the whole-Ångström marks, plus (when *stats* is supplied) the paired
    significance tests at the emphasised threshold."""
    W = 90
    ts = sorted(float(t) for t in within_df["rmsd_threshold_A"].unique())
    mark_ts = [t for t in ts if float(t).is_integer() and t >= 1.0] or ts[-4:]
    L = ["=" * W,
         "Fig 18 multi-depth — best-of-top-N accuracy vs distance threshold",
         "=" * W, "",
         f"Figure : {figure_name}",
         f"Metric : best-of-top-N pose within t Å"
         + (" AND PoseBusters-valid" if pb_valid_only else ""),
         f"Pool depths : {', '.join(f'top-{int(d)}' for d in depths)}"
         f"   (emphasised success line = {success_line:g} Å)",
         "",
         "-" * W, "FIGURE NOTES", "-" * W,
         "Each curve = % of a tool's complexes whose CLOSEST of its first d ranked poses is",
         "within t Å (PB-valid-gated when stated). top-1 = the rank-1 pose. The on-plot value",
         "labels are omitted by request; the numbers at the whole-Ångström marks are below.",
         "Rank: AutoDock/DiffDock native; EquiBind ranked by gnina affinity. " + metric_note,
         ""]
    for d in depths:
        L.append("-" * W)
        L.append(f"best of top-{int(d)}  (% of each tool's complexes within t Å)")
        L.append("-" * W)
        hdr = f"  {'threshold_A':>11s} " + " ".join(
            f"{TOOL_LABEL.get(m, m)[:15]:>16s}" for m in methods)
        L.append(hdr)
        for t in mark_ts:
            def _cell(m):
                r = within_df[(within_df["method"] == m) & (within_df["depth"] == int(d))
                              & (within_df["rmsd_threshold_A"] == t)]
                return f"{float(r['within_%'].iloc[0]):.1f}" if not r.empty else "n/a"
            L.append(f"  {t:>11.3g} " + " ".join(f"{_cell(m):>16s}" for m in methods))
        L.append("")

    if stats:
        def _p(p) -> str:
            if p is None or p != p:
                return "n/a"
            return "<1e-300" if p <= 0 else (f"{p:.3e}" if p < 1e-3 else f"{p:.4f}")
        L.append("=" * W)
        L.append(f"STATISTICS  (at the emphasised threshold t = {stats['threshold_A']:g} Å"
                 + (", PoseBusters-valid poses only)" if pb_valid_only else ")"))
        L.append("=" * W)
        L.append("Unit = one (protein, ligand) complex; binary = the tool has a"
                 + (" PB-valid" if pb_valid_only else "") + " pose within t Å")
        L.append("among its best-of-top-d ranked poses. Paired across tools (same complexes);")
        L.append("NESTED across depths (top-1 ⊂ top-15 ⊂ top-30).")
        L.append("")
        L.append("(A) Between tools, at each pool depth — Cochran's Q omnibus (>=3 related")
        L.append("    binary conditions) + pairwise exact McNemar, Holm-corrected:")
        for b in stats.get("between_tools_by_depth", []):
            if "omnibus" not in b:
                L.append(f"    top-{b['depth']}: n={b.get('n', 0)} — {b.get('status', '')}")
                continue
            om = b["omnibus"]
            rates = b.get("rates", {})
            rt = ", ".join(f"{TOOL_LABEL.get(m, m)} {100 * rates[m][0]:.1f}%" for m in rates)
            L.append(f"    top-{b['depth']} (n={b['n']}): Q={om['Q']:.2f}, df={om['df']}, "
                     f"p={_p(om['p'])}   [{rt}]")
            for pr in b.get("pairwise", []):
                L.append(f"        {TOOL_LABEL.get(pr['a'], pr['a'])} > "
                         f"{TOOL_LABEL.get(pr['b'], pr['b'])}: "
                         f"discordant +{pr['a_wins']}/−{pr['b_wins']}, "
                         f"p_holm={_p(pr['p_holm'])} {pr['star']}")
        L.append("")
        L.append("(B) Ranking headroom within a tool — extra complexes a DEEPER pool recovers")
        L.append("    (nested → reported as recovered share with Wilson 95% CI), plus Cochran's")
        L.append("    Q across the three depths as the omnibus that depth matters at all:")
        for h in stats.get("depth_headroom", []):
            if "omnibus_across_depths" not in h:
                L.append(f"    {TOOL_LABEL.get(h['method'], h['method'])}: "
                         f"n={h.get('n', 0)} — {h.get('status', '')}")
                continue
            om = h["omnibus_across_depths"]
            rbd = h["rate_by_depth"]
            rt = " → ".join(f"top-{d} {100 * rbd[d]:.1f}%" for d in sorted(rbd))
            L.append(f"    {TOOL_LABEL.get(h['method'], h['method'])} (n={h['n']}): {rt};  "
                     f"Q={om['Q']:.2f}, df={om['df']}, p={_p(om['p'])}")
            for s in h["steps"]:
                lo, hi = s["recovered_ci"]
                L.append(f"        top-{s['from_depth']}→top-{s['to_depth']}: "
                         f"+{s['recovered_k']} complexes "
                         f"({100 * s['recovered_share']:.1f}%, 95% CI "
                         f"{100 * lo:.1f}–{100 * hi:.1f}%)")
        L.append("")
        L.append("*** p<.001   ** p<.01   * p<.05   ns = not significant.")
        L.append("Nested depth steps are one-directional (a deeper pool can only ADD hits), so")
        L.append("they carry an effect size (recovered share + CI), not a McNemar p-value.")
        L.append("")

    L.append("=" * W)
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _stats_within_by_depth(df: pd.DataFrame, depths, threshold: float,
                           eq_df: "pd.DataFrame | None" = None,
                           pb_valid_only: bool = False, rmsd_col: str = "rmsd",
                           figure_name: str = "", metric_label: str = "RMSD") -> "dict | None":
    """Paired significance tests for the multi-depth fig 18 view, AT the emphasised
    threshold (2 Å for as-placed RMSD, 1 Å for Kabsch RMSD).

    The datum is one binary per (protein, ligand) complex: does the tool have a
    (PB-valid) pose whose best-of-top-d RMSD is <= *threshold*. Two questions the graph
    poses, each with the test its data structure warrants:

      * between_tools_by_depth — at each pool depth, do the three tools differ? The same
        complexes are scored by every tool, so the data are PAIRED binary across >=3
        conditions → Cochran's Q omnibus + pairwise exact McNemar (Holm-corrected),
        via :func:`stats_utils.paired_proportions`.
      * depth_headroom — within a tool, does deepening the ranked pool (top-1 → top-15 →
        top-30) add recovery? The hit sets are NESTED (rank-1 ⊂ top-15 ⊂ top-30), so the
        discordance is one-directional and a McNemar p is degenerate; the gain is the
        recovered share of complexes (Wilson 95 % CI), with Cochran's Q across the three
        depths as the omnibus that depth matters at all.

    Returns the sidecar dict, or None when < 2 tools carry per-complex frames.
    """
    depths = sorted(int(d) for d in depths)
    tools = ("autodock", "diffdock", "equibind")
    hits: dict = {}
    present: set = set()
    for d in depths:
        frames = _topn_within_frames(df, d, eq_df=eq_df,
                                     pb_valid_only=pb_valid_only, rmsd_col=rmsd_col)
        hd = {}
        for m in tools:
            if m in frames:
                rec = pd.to_numeric(frames[m]["best_rmsd"], errors="coerce") <= threshold
                hd[m] = {f"{p}||{l}": bool(v)
                         for (p, l), v in zip(frames[m].index, rec)}
                present.add(m)
        hits[d] = hd
    methods = [m for m in tools if m in present]
    if len(methods) < 2:
        return None

    between = []
    for d in depths:
        hd = hits.get(d, {})
        mset = [m for m in methods if m in hd]
        if len(mset) < 2:
            continue
        common = sorted(set.intersection(*[set(hd[m]) for m in mset]))
        if len(common) < _MIN_UNITS_STATS:
            between.append({"depth": d, "n": len(common),
                            "status": "n too small — exploratory"})
            continue
        data = {m: [hd[m][c] for c in common] for m in mset}
        between.append({"depth": d, "n": len(common), **su.paired_proportions(data)})

    headroom = []
    for m in methods:
        dd = [d for d in depths if m in hits.get(d, {})]
        if len(dd) < 2:
            continue
        cids = sorted(set.intersection(*[set(hits[d][m]) for d in dd]))
        N = len(cids)
        if N < _MIN_UNITS_STATS:
            headroom.append({"method": m, "n": N, "status": "n too small — exploratory"})
            continue
        H = {d: np.fromiter((hits[d][m][c] for c in cids), bool, N) for d in dd}
        Q, pQ, dQ = su.cochran_q(np.array([H[d].astype(float) for d in dd]).T)
        pairs = list(zip(dd[:-1], dd[1:]))
        if len(dd) > 2:
            pairs.append((dd[0], dd[-1]))          # end-to-end headroom too
        steps = []
        for lo, hi in pairs:
            rec = int((H[hi] & ~H[lo]).sum())      # nested: only 0→1 transitions exist
            lo_ci, hi_ci = su.wilson_ci(rec, N)
            steps.append({"from_depth": lo, "to_depth": hi, "recovered_k": rec,
                          "recovered_share": rec / N,
                          "recovered_ci": [lo_ci, hi_ci], "nested": True})
        headroom.append({
            "method": m, "n": N,
            "k_by_depth": {int(d): int(H[d].sum()) for d in dd},
            "rate_by_depth": {int(d): float(H[d].mean()) for d in dd},
            "omnibus_across_depths": {"Q": Q, "p": pQ, "df": dQ},
            "steps": steps,
        })

    return {
        "figure": figure_name,
        "threshold_A": float(threshold),
        "metric": metric_label,
        "pb_valid_only": bool(pb_valid_only),
        "unit": "(protein, ligand) complex; paired across tools, nested across depths",
        "depths": depths,
        "between_tools_by_depth": between,
        "depth_headroom": headroom,
    }


def plot_top1_vs_oracle_scatter(df: pd.DataFrame, out: Path) -> None:
    """Scatter: top-1 RMSD (y-axis) vs oracle RMSD (x-axis) per pair.

    Points above the diagonal are pairs where rank-1 is worse than the
    method's best available pose — i.e. ranking failure.  The further above
    the diagonal, the more the ranking misleads.
    """
    ranking_methods = sorted(m for m in df["method"].unique() if m in RANKING_TOOLS)
    if not ranking_methods:
        return

    oracle = (_oracle_per_pair(df)
              .set_index(["method", "protein", "ligand"])["rmsd"]
              .rename("oracle_rmsd"))
    top1 = (_top1_per_pair(df)
             .set_index(["method", "protein", "ligand"])["rmsd"]
             .rename("top1_rmsd"))
    merged = pd.concat([oracle, top1], axis=1).dropna().reset_index()

    n = len(ranking_methods)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 5.2))
    if n == 1:
        axes = [axes]

    lim = 15
    for ax, method in zip(axes, ranking_methods):
        sub = merged[merged["method"] == method]
        ax.scatter(sub["oracle_rmsd"], sub["top1_rmsd"],
                   s=18, alpha=0.5,
                   color=TOOL_COLORS.get(method, "grey"),
                   edgecolors="none")
        ax.plot([0, lim], [0, lim], color="red", ls="--", alpha=0.7, lw=1.5,
                label="perfect ranking (rank 1 = oracle)")
        ax.axhline(2.0, color="grey", ls=":", alpha=0.5, lw=1)
        ax.axvline(2.0, color="grey", ls=":", alpha=0.5, lw=1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("Oracle RMSD (Å) — best achievable pose")
        ax.set_ylabel("Top-1 RMSD (Å) — rank-1 pose")
        ax.set_title(f"{TOOL_LABEL.get(method, method)} (n={len(sub)})")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)

        n_above = int((sub["top1_rmsd"] > sub["oracle_rmsd"] + 0.5).sum())
        pct_above = 100 * n_above / len(sub) if len(sub) else 0
        ax.text(0.04, 0.96,
                f"ranking worse than oracle\nby > 0.5 Å: {pct_above:.0f}% of pairs",
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    if n > 1:
        _label_panels(axes)
    fig.suptitle(_vt("Top-1 vs oracle RMSD — where does ranking fail?"), fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_oracle_rank_histogram(df: pd.DataFrame, top_n: int, out: Path) -> None:
    """Histogram: which rank does the oracle (best) pose receive per pair?

    A distribution concentrated at rank 1 indicates reliable scoring.
    Poses whose oracle rank exceeds top_n are collected in an overflow bin.
    """
    ranking_methods = sorted(m for m in df["method"].unique() if m in RANKING_TOOLS)
    if not ranking_methods:
        return

    oracle = _oracle_per_pair(df)
    oracle_ranked = oracle[oracle["method"].isin(RANKING_TOOLS)].copy()
    # Clamp ranks > top_n into an overflow bin labelled "> N"
    oracle_ranked["rank_bin"] = oracle_ranked["rank"].clip(upper=top_n + 1)

    n = len(ranking_methods)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.2), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, method in zip(axes, ranking_methods):
        sub = oracle_ranked[oracle_ranked["method"] == method]["rank_bin"].dropna()
        bins = np.arange(0.5, top_n + 2.5, 1)
        counts, _, patches = ax.hist(sub, bins=bins,
                                     color=TOOL_COLORS.get(method, "grey"),
                                     edgecolor="black", alpha=0.8)
        # Highlight rank-1 bar in a stronger colour
        patches[0].set_edgecolor("red"); patches[0].set_linewidth(2)
        pct_at_1 = 100 * (sub == 1).sum() / len(sub) if len(sub) else 0

        tick_labels = [str(k) for k in range(1, top_n + 1)] + [f"> {top_n}"]
        ax.set_xticks(range(1, top_n + 2))
        ax.set_xticklabels(tick_labels)
        ax.set_xlabel("Rank assigned to the oracle (best) pose")
        ax.set_ylabel("Number of receptor-ligand pairs")
        ax.set_title(f"{TOOL_LABEL.get(method, method)} (n={len(sub)})\n"
                     f"{pct_at_1:.0f}% of best poses ranked 1st")
        ax.grid(axis="y", alpha=0.3)

    if n > 1:
        _label_panels(axes)
    fig.suptitle(_vt("Oracle pose rank distribution — how reliably does rank 1 = best pose?"),
                 fontsize=12)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Plots — cross-tool comparison (oracle RMSD)
# ───────────────────────────────────────────────────────────────────


def plot_cross_tool_pairwise(df: pd.DataFrame, out: Path) -> None:
    """Per-pair oracle RMSD comparison across every tool pair.

    Diagonal: RMSD histogram for each tool.
    Upper triangle: delta histogram (row − col), showing which method is better.
    Lower triangle: scatter of oracle RMSDs between two tools.
    """
    oracle = _oracle_per_pair(df)
    pivot = oracle.pivot_table(index=["protein", "ligand"],
                               columns="method", values="rmsd")
    methods = list(pivot.columns)
    if len(methods) < 2:
        return
    n = len(methods)
    fig, axes = plt.subplots(n, n, figsize=(3.5 * n, 3.5 * n))
    if n == 1:
        axes = np.array([[axes]])

    for i, mi in enumerate(methods):
        for j, mj in enumerate(methods):
            ax = axes[i, j]
            if i == j:
                vals = pivot[mi].dropna()
                ax.hist(np.clip(vals, 0, 15), bins=30,
                        color=TOOL_COLORS.get(mi, "#888"),
                        edgecolor="black", alpha=0.8)
                ax.axvline(2.0, color="red", ls="--", alpha=0.6)
                ax.set_title(f"{TOOL_LABEL.get(mi, mi)} (n={len(vals)})", fontsize=10)
                ax.set_xlim(0, 15)
                ax.set_xlabel("Oracle RMSD (Å)")
            elif i < j:
                pair = pivot[[mi, mj]].dropna()
                if len(pair) == 0:
                    ax.set_visible(False); continue
                d = (pair[mi] - pair[mj]).values
                ax.hist(np.clip(d, -10, 10), bins=30,
                        color="#4c72b0", edgecolor="black")
                ax.axvline(0.0, color="red", ls="--", alpha=0.6)
                wins_i = int((pair[mi] < pair[mj]).sum())
                ax.set_title(
                    f"{TOOL_LABEL.get(mi, mi)[:8]} − {TOOL_LABEL.get(mj, mj)[:8]}\n"
                    f"row better in {wins_i}/{len(pair)}",
                    fontsize=8)
            else:
                pair = pivot[[mj, mi]].dropna()
                if len(pair) == 0:
                    ax.set_visible(False); continue
                ax.scatter(pair[mj], pair[mi], s=10, alpha=0.5,
                           color=TOOL_COLORS.get(mi, "grey"))
                lim = 15
                ax.plot([0, lim], [0, lim], color="red", ls="--", alpha=0.5)
                ax.set_xlim(0, lim); ax.set_ylim(0, lim)
                ax.set_xlabel(TOOL_LABEL.get(mj, mj), fontsize=9)
                ax.set_ylabel(TOOL_LABEL.get(mi, mi), fontsize=9)
            ax.grid(alpha=0.3)

    _label_panels(axes)
    fig.suptitle(_vt("Cross-tool comparison — per-pair oracle (best) RMSD"), fontsize=13)
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Plots — "twisted & turned" + oracle-rank distribution
# ───────────────────────────────────────────────────────────────────

# (column, axis label, which deformation it captures). The pose is compared to
# the crystal ligand: how far it MOVED (translation), how it TURNED (rigid-body
# rotation) and how it was TWISTED (internal conformation / torsions / strain).
TWIST_TURN_MEASURES: list[tuple[str, str, str]] = [
    ("centroid_dist",        "Translation (Å)",         "moved"),
    ("rot_angle_deg",        "Rotation (°)",            "turned"),
    ("bestfit_rmsd",         "Best-fit RMSD (Å)",       "twisted"),
    ("tfd",                  "TFD (0–1)",               "twisted"),
    ("max_torsion_dev_deg",  "Max torsion Δ (°)",       "twisted"),
    ("mean_torsion_dev_deg", "Mean torsion Δ (°)",      "twisted"),
    ("n_torsions_flipped",   "Torsions flipped (>90°)", "twisted"),
    ("strain_energy",        "UFF strain (kcal/mol)",   "twisted"),
]


def _twist_top1_per_pair(df: pd.DataFrame) -> pd.DataFrame:
    """Rank-1 pose per (method, protein, ligand) for the twist/turn figure.

    AutoDock/DiffDock use their native confidence rank. EquiBind has no native rank,
    so — matching how EquiBind is ranked elsewhere in the report (e.g. fig 18) — its
    poses are ranked by gnina AFFINITY (``_gnina_affinity_rank``, most-negative = rank 1)
    and the gnina-rank-1 pose is taken. Only EquiBind variants that actually carry gnina
    affinities get a row, so the box reflects a genuine gnina ranking rather than a
    generation-order fallback; a variant without them (raw / smina-opt) yields no top-1
    box, exactly as before.
    """
    parts = []
    ranked = df[df["method"].isin(RANKING_TOOLS)]
    if not ranked.empty:
        parts.append(ranked.sort_values("rank")
                           .groupby(["method", "protein", "ligand"]).head(1))
    eq = df[df["method"].astype(str).str.startswith("equibind")]
    for method, g in eq.groupby("method"):
        if pd.to_numeric(g["gnina_affinity"], errors="coerce").notna().any():
            g = g.copy()
            g["_sr"] = _gnina_affinity_rank(g)
            parts.append(g[g["_sr"] == 1].drop(columns="_sr"))
    if not parts:
        return df.iloc[0:0].reset_index(drop=True)
    return pd.concat(parts, ignore_index=True)


def _twist_best_topn_per_pair(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """Lowest-RMSD pose among each tool's top-``top_n`` RANKED poses, per complex.

    The 'oracle within the ranked shortlist' — the best a perfect rescorer could pick
    from the first ``top_n`` poses a user would actually look at, rather than the best of
    ALL docked poses. Ranking convention matches :func:`_twist_top1_per_pair`: native
    confidence rank for AutoDock/DiffDock, gnina-affinity rank for EquiBind (variants
    without gnina affinities are dropped, since they have no ranking to shortlist on).
    """
    parts = []
    ranked = df[df["method"].isin(RANKING_TOOLS)].copy()
    ranked["_rk"] = pd.to_numeric(ranked["rank"], errors="coerce")
    parts.append(ranked[ranked["_rk"] <= top_n])
    eq = df[df["method"].astype(str).str.startswith("equibind")]
    for method, g in eq.groupby("method"):
        if pd.to_numeric(g["gnina_affinity"], errors="coerce").notna().any():
            g = g.copy()
            g["_rk"] = _gnina_affinity_rank(g)
            parts.append(g[g["_rk"] <= top_n])
    pooled = pd.concat(parts, ignore_index=True) if parts else df.iloc[0:0]
    valid = pooled.dropna(subset=["rmsd"])
    if valid.empty:
        return valid.reset_index(drop=True)
    idx = valid.groupby(["method", "protein", "ligand"])["rmsd"].idxmin()
    return valid.loc[idx.dropna()].drop(columns="_rk").reset_index(drop=True)


def aggregate_twist_turn(df: pd.DataFrame, top_n: int = DEFAULT_TOP_N) -> pd.DataFrame:
    """Per-method median & mean of every twist/turn measure, for top-1 & the
    best-of-top-``top_n`` (ranked-shortlist oracle) selection.

    Long-form table (one row per selection × method). The companion plot draws
    the full distributions; this CSV gives the central tendencies for the text.
    """
    sels = {f"best_top{top_n}": _twist_best_topn_per_pair(df, top_n),
            "top1": _twist_top1_per_pair(df)}
    rows = []
    for sel_name, sub in sels.items():
        if sub.empty:
            continue
        for method, g in sub.groupby("method"):
            row = {"selection": sel_name, "method": method, "n_pairs": len(g)}
            for col, _label, _axis in TWIST_TURN_MEASURES:
                if col in g.columns:
                    row[f"{col}_median"] = round(float(g[col].median(skipna=True)), 3)
                    row[f"{col}_mean"] = round(float(g[col].mean(skipna=True)), 3)
            rows.append(row)
    return pd.DataFrame(rows)


def plot_twist_turn(df: pd.DataFrame, out: Path, top_n: int = DEFAULT_TOP_N,
                    pb_valid_only: bool = False) -> None:
    """Distributions of how each tool re-orients (turns) and re-bends (twists)
    the ligand relative to its crystal pose — top-1 vs best-of-top-N, per measure.

    One panel per measure; within a panel, each tool gets a solid 'best of top-N
    ranked' box (the lowest-RMSD pose among its first ``top_n`` ranked poses) and a
    hatched 'top-1' (rank-1) box. Ranking is the tool's native confidence rank for
    AutoDock/DiffDock and the gnina-affinity rank for EquiBind (which has no native
    ranking); EquiBind variants without gnina affinities show neither box. Outliers
    hidden so the bulk of each distribution stays readable.

    ``pb_valid_only`` only changes the footnote: pass True when ``df`` has already
    been filtered to PoseBusters-valid poses (the ``pb_valid/`` variant) so the note
    records that invalid poses were excluded and that "top-1" is therefore each
    tool's highest-ranked *valid* pose, not necessarily its true rank-1.
    """
    from matplotlib.patches import Patch
    oracle = _twist_best_topn_per_pair(df, top_n)
    top1 = _twist_top1_per_pair(df)
    methods = list(dict.fromkeys(oracle["method"]))
    if not methods:
        return

    ncol = 4
    nrow = int(np.ceil(len(TWIST_TURN_MEASURES) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4.0 * nrow))
    axes = np.atleast_1d(axes).ravel()

    for ax, (col, label, axis_kind) in zip(axes, TWIST_TURN_MEASURES):
        data, positions, facecolors, hatches = [], [], [], []
        tickpos, ticklabels = [], []
        for mi, method in enumerate(methods):
            o = oracle.loc[oracle["method"] == method, col].dropna() if col in oracle else pd.Series(dtype=float)
            # top1 holds a rank-1 row only for methods that HAVE a ranking (native for
            # AutoDock/DiffDock, gnina-affinity for EquiBind) — so an empty slice here
            # just means "no rank-1 box for this tool".
            t = (top1.loc[top1["method"] == method, col].dropna()
                 if col in top1 else pd.Series(dtype=float))
            c = TOOL_COLORS.get(method, "grey")
            if len(o):
                data.append(o.values); positions.append(mi - 0.18)
                facecolors.append(c); hatches.append(None)
            if len(t):
                data.append(t.values); positions.append(mi + 0.18)
                facecolors.append(c); hatches.append("//")
            tickpos.append(mi); ticklabels.append(TOOL_LABEL.get(method, method))
        if not data:
            ax.set_visible(False); continue
        bp = ax.boxplot(data, positions=positions, widths=0.32,
                        patch_artist=True, showfliers=False)
        for patch, c, h in zip(bp["boxes"], facecolors, hatches):
            patch.set_facecolor(c)
            patch.set_alpha(0.85 if h is None else 0.40)
            if h:
                patch.set_hatch(h)
        for med in bp["medians"]:
            med.set_color("black")
        ax.set_xticks(tickpos)
        ax.set_xticklabels(ticklabels, rotation=25, ha="right",
                           rotation_mode="anchor", fontsize=8)
        ax.set_title(f"{label}\n({axis_kind})", fontsize=10)
        ax.grid(axis="y", alpha=0.3)

    for ax in axes[len(TWIST_TURN_MEASURES):]:
        ax.set_visible(False)

    _label_panels(axes)
    fig.legend(handles=[
        Patch(facecolor="#888888", alpha=0.85,
              label=f"Best of top-{top_n} ranked poses"),
        Patch(facecolor="#888888", alpha=0.40, hatch="//", label="Top-1 (tool's rank-1 pose)"),
    ], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.015), fontsize=10)
    fig.suptitle(_vt("How the ligand is twisted & turned vs its crystal pose\n"
                     "(moved = translation · turned = rigid-body rotation · "
                     "twisted = internal conformation)"),
                 fontsize=13, fontweight="bold", y=1.07)
    footnote = (f"Both boxes are drawn from each tool's first {top_n} ranked poses. "
                "Ranking = native confidence rank (AutoDock, DiffDock); EquiBind has no native "
                "ranking, so it is ranked by gnina affinity.")
    if pb_valid_only:
        footnote += ("\nOnly PoseBusters-valid poses are considered here, so the top-1 box is "
                     "each tool's highest-ranked valid pose — not necessarily its true rank-1 pose.")
    fig.text(0.5, 0.005, footnote, ha="center", va="bottom", fontsize=8, color="0.30")
    fig.tight_layout(rect=(0, 0.02, 1, 0.97))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Plots — form fidelity of the near-native, physically-valid poses
# ───────────────────────────────────────────────────────────────────
# The strict success criterion (RMSD ≤ 2 Å AND PoseBusters-valid) is scored on
# the IN-PLACE RMSD, which conflates two errors: where the pose sits in the
# pocket (translation + rotation) and whether its internal conformation — its
# "form" — matches the crystal ligand. Conditioning on the poses we already call
# a success, this asks the leftover question: how good is the FORM, and is the
# residual deviation placement-limited (shape right, position off) or form-limited
# (the conformation itself is wrong)?
#
# Form error = best-fit (Kabsch) RMSD (``bestfit_rmsd`` — PoseBusters'
# pb_kabsch_rmsd, else our own Kabsch fit): heavy-atom RMSD after optimal
# superposition, translation + rotation removed. Placement error follows in
# quadrature from the exact per-pose identity in-place² = placement² + form²,
# so placement = √(max(0, in-place² − form²)).


def _near_native_valid_reps(df: pd.DataFrame,
                            rmsd_thr: float = NEAR_NATIVE_RMSD_A) -> pd.DataFrame:
    """Per (method, protein, ligand): the oracle (min in-place RMSD) pose, kept
    only when it is BOTH within ``rmsd_thr`` Å AND PoseBusters-valid.

    This is exactly the population behind ``oracle_pb_valid_and_rmsd2_%`` — the
    strict-success representatives — so the form analysis is conditioned on the
    same poses the report already counts as successful. Empty (→ callers no-op)
    for crystal-free sets where ``rmsd`` is undefined.
    """
    oracle = _oracle_per_pair(df)
    if oracle.empty:
        return oracle
    keep = (oracle["rmsd"] <= rmsd_thr) & oracle["pb_valid"].astype(bool)
    return oracle[keep].reset_index(drop=True)


def _best_valid_near_native_reps(df: pd.DataFrame,
                                 rmsd_thr: float = NEAR_NATIVE_RMSD_A) -> pd.DataFrame:
    """Per (method, protein, ligand): the NEAREST pose that is BOTH ≤ ``rmsd_thr`` Å
    AND PoseBusters-valid (validity-constrained oracle).

    Unlike :func:`_near_native_valid_reps` — which takes the single closest pose and
    only *then* checks validity, dropping the complex when that closest pose fails
    PoseBusters — this keeps a complex whenever ANY ≤ ``rmsd_thr`` Å valid pose
    exists, representing it by the closest such pose. It is the more generous
    sampling-ceiling reading ("can the tool produce a pose that is simultaneously
    near-native and valid?") and never counts fewer successes than the nearest rule.
    """
    ok = df.dropna(subset=["rmsd"])
    ok = ok[(ok["rmsd"] <= rmsd_thr) & ok["pb_valid"].astype(bool)]
    if ok.empty:
        return ok.reset_index(drop=True)
    idx = ok.groupby(["method", "protein", "ligand"])["rmsd"].idxmin()
    return ok.loc[idx.dropna()].reset_index(drop=True)


def _top1_valid_reps(df: pd.DataFrame,
                     rmsd_thr: float = NEAR_NATIVE_RMSD_A) -> pd.DataFrame:
    """Per (method, protein, ligand): the tool's RANK-1 pose, kept only if that
    single pose is itself ≤ ``rmsd_thr`` Å AND PoseBusters-valid.

    The realistic top-1 reading — the counterpart to the two oracle *ceilings*. The
    oracle selectors ask "can the tool put a near-native valid pose *somewhere* in
    its sample?"; this asks "does the pose you would actually take — its #1 ranked —
    succeed?". Effective rank is native confidence for RANKING_TOOLS and gnina-
    affinity for the blind tools (:func:`_effective_rank`), so every family is read
    at the pose a user is handed first. Because depth is 1 this yields at most one
    row per complex, exactly like the oracle selectors, so the panels stay
    apples-to-apples (same axes, same per-variant multiplicity — only the *pick*
    differs). Empty for crystal-free sets (``rmsd`` undefined → callers no-op).
    """
    return _valid_topd_poses(df, 1, rmsd_gate=rmsd_thr).reset_index(drop=True)


# The oracle selections + the realistic top-1 that the form-fidelity figures are
# rendered under. The note is stamped into each figure title so the graphs are
# self-identifying; plot_oracle_selection_comparison puts the two oracles head-to-head.
_SEL_NOTE_NEAREST = ("selection: RMSD-greedy oracle — the single nearest pose, "
                     "kept only if it is itself ≤ 2 Å AND PB-valid")
_SEL_NOTE_VALID = ("selection: validity-constrained oracle — the NEAREST pose that "
                   "is ≤ 2 Å AND PB-valid (more generous sampling ceiling)")
_SEL_NOTE_TOP1 = ("selection: realistic top-1 — each tool's RANK-1 pose, kept only "
                  "if that pose is itself ≤ 2 Å AND PB-valid (what you actually get first)")


def _form_components(reps: pd.DataFrame) -> pd.DataFrame:
    """Add ``form`` (best-fit/Kabsch RMSD), ``inplace`` (as-placed RMSD) and
    ``placement`` = √(max(0, in-place² − form²)) — the conformation-vs-positioning
    split of each pose's deviation from the crystal ligand."""
    reps = reps.copy()
    reps["form"] = pd.to_numeric(reps["bestfit_rmsd"], errors="coerce")
    inplace = reps["pb_rmsd"].where(reps["pb_rmsd"].notna(), reps["rmsd"])
    reps["inplace"] = pd.to_numeric(inplace, errors="coerce")
    resid = reps["inplace"] ** 2 - reps["form"] ** 2
    reps["placement"] = np.sqrt(resid.clip(lower=0))
    return reps


def aggregate_form_fidelity(df: pd.DataFrame,
                            form_ok: float = FORM_OK_KABSCH_A,
                            rmsd_thr: float = NEAR_NATIVE_RMSD_A,
                            selector=_near_native_valid_reps) -> pd.DataFrame:
    """Per method: among the near-native (≤ ``rmsd_thr`` Å) PB-valid poses, how
    well the internal conformation ("form") reproduces the crystal ligand.

    ``selector`` picks the per-complex representative — :func:`_near_native_valid_reps`
    (nearest pose, RMSD-greedy, the default) or :func:`_best_valid_near_native_reps`
    (nearest ≤ 2 Å valid pose, validity-constrained).

    Columns: ``n_pairs`` (all complexes the method covers), ``n_success`` /
    ``success_%`` (the strict-success reps; with the default selector this matches
    ``oracle_pb_valid_and_rmsd2_%``) and, over those reps, the best-fit (Kabsch) RMSD
    distribution, ``form_correct_%`` (best-fit ≤ ``form_ok``), TFD, and the quadrature
    split of the deviation into form vs placement (``form_share_of_error_%`` =
    mean form² / mean in-place² — the share of the residual mean-square deviation
    that is conformational rather than positional; > 50 % ⇒ form-limited).
    """
    reps = _form_components(selector(df, rmsd_thr))
    n_tot = df.drop_duplicates(["method", "protein", "ligand"]).groupby("method").size()
    rows = []
    for method in dict.fromkeys(df["method"]):
        n_pairs = int(n_tot.get(method, 0))
        g = reps[reps["method"] == method]
        n = len(g)
        row = {"method": method, "n_pairs": n_pairs, "n_success": n,
               "success_%": round(100 * n / n_pairs, PCT_DECIMALS) if n_pairs else float("nan")}
        if n:
            form = g["form"].dropna()
            ms_in = float((g["inplace"] ** 2).mean(skipna=True))
            ms_form = float((g["form"] ** 2).mean(skipna=True))
            row.update({
                "form_bestfit_median": round(float(form.median()), 3),
                "form_bestfit_mean": round(float(form.mean()), 3),
                "form_bestfit_p90": round(float(form.quantile(0.9)), 3),
                "tfd_median": round(float(g["tfd"].median(skipna=True)), 3),
                "form_correct_n": int((form <= form_ok).sum()),
                "form_correct_%": round(100 * float((form <= form_ok).mean()), PCT_DECIMALS),
                "rms_inplace": round(float(np.sqrt(ms_in)), 3),
                "rms_form": round(float(np.sqrt(ms_form)), 3),
                "rms_placement": round(float(np.sqrt((g["placement"] ** 2).mean(skipna=True))), 3),
                "form_share_of_error_%": round(100 * ms_form / ms_in, PCT_DECIMALS) if ms_in else float("nan"),
            })
        rows.append(row)
    return pd.DataFrame(rows)


def _fam_key(method: str) -> str:
    """Collapse a method to its tool family for the per-pose scatter colouring.

    Uni-Dock is matched EXPLICITLY and before the EquiBind fallback. This function used
    to end in an unconditional ``return "equibind"``, so every engine that was neither
    autodock* nor diffdock* was silently classified as EquiBind. Because
    :func:`_select_presentation_tools` deliberately keeps unidock and unidock2 in the
    collapsed presentation frame, Uni-Dock2 poses were being pooled into the EquiBind
    scatter of every family-pooled figure and CSV built from the whole-protein cache —
    more than doubling its rank-1 cohort and letting one complex contribute several
    "rank-1" points. Family-pooled plots iterate over the three tool families, so giving
    Uni-Dock its own key removes it from them rather than mislabelling it."""
    if method.startswith("autodock"):
        return "autodock"
    if method.startswith("diffdock"):
        return "diffdock"
    if method.startswith(("unidock", "uni_dock")):
        return "unidock"
    return "equibind"


def plot_form_fidelity(df: pd.DataFrame, out: Path,
                       form_ok: float = FORM_OK_KABSCH_A,
                       rmsd_thr: float = NEAR_NATIVE_RMSD_A,
                       selector=_near_native_valid_reps,
                       sel_note: str = _SEL_NOTE_NEAREST) -> None:
    """Form fidelity of the near-native, PB-valid poses (the strict successes).

    (A) box of best-fit (Kabsch) RMSD per method — how far the form is;
    (B) % of those poses whose form is correct (best-fit ≤ form_ok);
    (C) 100 %-stacked split of each method's mean-square deviation into form vs
        placement (is the residual error conformation- or positioning-limited);
    (D) per-pose in-place RMSD vs form (best-fit) RMSD, coloured by tool family —
        points near the x-axis are placement-limited (form perfect), points near
        the y = x bound are form-limited.

    ``selector`` chooses the per-complex representative and ``sel_note`` is stamped
    into the title so the figure is self-identifying (see _SEL_NOTE_*).
    """
    from matplotlib.patches import Patch
    reps = _form_components(selector(df, rmsd_thr))
    reps = reps.dropna(subset=["form", "inplace"])
    if reps.empty:
        return
    methods = [m for m in dict.fromkeys(reps["method"]) if (reps["method"] == m).any()]
    agg = aggregate_form_fidelity(df, form_ok, rmsd_thr, selector).set_index("method")

    # Split into two figures: (A+B) error & %-correct, (C+D) mechanism & per-pose cloud.
    out_ab = out.with_stem(out.stem + "_AB")
    out_cd = out.with_stem(out.stem + "_CD")
    suptitle = _vt(f"Form fidelity of the successful poses "
                   f"(RMSD ≤ {rmsd_thr:g} Å AND PoseBusters-valid)\n"
                   "how closely the docked ligand's internal conformation matches the "
                   f"crystal, once placement is removed\n{sel_note}")

    # ── Figure 1 — (A) form-error distribution + (B) % form-correct ─────────
    figAB, (axA, axB) = plt.subplots(1, 2, figsize=(13, 6))

    dataA, colsA, labelsA = [], [], []
    for m in methods:
        v = reps.loc[reps["method"] == m, "form"].dropna()
        if not len(v):
            continue
        dataA.append(v.values)
        colsA.append(TOOL_COLORS.get(m, "grey"))
        labelsA.append(f"{TOOL_LABEL.get(m, m)}\n(n={len(v)})")
    posA = range(len(dataA))
    bp = axA.boxplot(dataA, positions=list(posA), widths=0.6,
                     patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], colsA):
        patch.set_facecolor(c); patch.set_alpha(0.85)
    for med in bp["medians"]:
        med.set_color("black")
    axA.axhline(form_ok, ls="--", color="crimson", lw=1.2)
    axA.text(0.99, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
             va="bottom", ha="right", fontsize=8, transform=axA.get_yaxis_transform())
    axA.set_xticks(list(posA))
    axA.set_xticklabels(labelsA, rotation=25, ha="right", rotation_mode="anchor", fontsize=8)
    axA.set_ylabel(f"Form error — best-fit (Kabsch) RMSD to {_ref_noun()} (Å)")
    axA.set_title("How far is the form of each near-native, valid pose?", fontsize=11)
    axA.grid(axis="y", alpha=0.3)

    fc = [float(agg.loc[m, "form_correct_%"]) if m in agg.index else float("nan")
          for m in methods]
    colsB = [TOOL_COLORS.get(m, "grey") for m in methods]
    xb = range(len(methods))
    axB.bar(list(xb), fc, color=colsB, alpha=0.9)
    for x, m, v in zip(xb, methods, fc):
        if not np.isnan(v):
            n_ok = int(agg.loc[m, "form_correct_n"]); n_s = int(agg.loc[m, "n_success"])
            axB.text(x, v + 1, f"{v:.0f}%\n{n_ok}/{n_s}", ha="center", va="bottom", fontsize=8)
    axB.set_xticks(list(xb))
    axB.set_xticklabels([TOOL_LABEL.get(m, m) for m in methods], rotation=25,
                        ha="right", rotation_mode="anchor", fontsize=8)
    axB.set_ylim(0, 108)
    axB.set_ylabel(f"Poses with correct form (%)  [best-fit RMSD ≤ {form_ok:g} Å]")
    axB.set_title("How many near-native, valid poses get the form right?", fontsize=11)
    axB.grid(axis="y", alpha=0.3)

    _label_panels([axA, axB])
    figAB.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.02)
    figAB.tight_layout(rect=(0, 0, 1, 0.95))
    figAB.savefig(out_ab, dpi=160, bbox_inches="tight"); plt.close(figAB)

    # ── Figure 2 — (C) form-vs-placement share + (D) per-pose cloud ─────────
    figCD, (axC, axD) = plt.subplots(1, 2, figsize=(14, 6.5))

    share = [float(agg.loc[m, "form_share_of_error_%"]) if m in agg.index else float("nan")
             for m in methods]
    xc = range(len(methods))
    form_seg = [s if not np.isnan(s) else 0.0 for s in share]
    place_seg = [100 - s if not np.isnan(s) else 0.0 for s in share]
    axC.bar(list(xc), form_seg, color="#8856a7", label="Form (internal conformation)")
    axC.bar(list(xc), place_seg, bottom=form_seg, color="#c7c7c7",
            label="Placement (translation + rotation)")
    for x, m in zip(xc, methods):
        if m in agg.index and not np.isnan(agg.loc[m, "rms_inplace"]):
            axC.text(x, 101, f"{agg.loc[m,'rms_inplace']:.2f} Å", ha="center",
                     va="bottom", fontsize=7, color="#333333")
    axC.axhline(50, ls=":", color="black", lw=0.8)
    axC.set_xticks(list(xc))
    axC.set_xticklabels([TOOL_LABEL.get(m, m) for m in methods], rotation=25,
                        ha="right", rotation_mode="anchor", fontsize=8)
    axC.set_ylim(0, 108)
    axC.set_ylabel("Share of near-native deviation, mean-square (%)")
    axC.set_title("Is the residual error placement- or form-limited?\n"
                  "(bar labels = total root-mean-square deviation)", fontsize=11)
    axC.legend(loc="lower center", fontsize=8, ncol=2, frameon=False,
               bbox_to_anchor=(0.5, -0.02))

    fam_color = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind": "#2ca02c"}
    fam_label = {"autodock": "AutoDock", "diffdock": "DiffDock", "equibind": "EquiBind"}
    fams_present = []
    for fam in ("autodock", "diffdock", "equibind"):
        sub = reps[reps["method"].map(_fam_key) == fam]
        if sub.empty:
            continue
        fams_present.append(fam)
        axD.scatter(sub["inplace"], sub["form"], s=16, alpha=0.45,
                    color=fam_color[fam], edgecolors="none")
    lim = float(np.nanmax([reps["inplace"].max(), reps["form"].max(), rmsd_thr])) * 1.05
    axD.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
    axD.text(lim, lim, " y = x (all error is form)", color="grey", fontsize=8,
             va="top", ha="right")
    axD.axhline(form_ok, ls="--", color="crimson", lw=1)
    axD.axvline(rmsd_thr, ls=":", color="black", lw=0.8)
    axD.set_xlim(0, lim); axD.set_ylim(0, lim)
    axD.set_xlabel(f"In-place RMSD to {_ref_noun()} (Å)")
    axD.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)")
    axD.set_title("Placement- vs form-limited, pose by pose", fontsize=11)
    axD.grid(alpha=0.3)
    axD.legend(handles=[Patch(facecolor=fam_color[f], label=fam_label[f])
                        for f in fams_present], fontsize=8, frameon=False,
               loc="upper left")

    axC.set_title("(C)", loc="left", fontweight="bold", fontsize=13)
    axD.set_title("(D)", loc="left", fontweight="bold", fontsize=13)
    figCD.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.02)
    figCD.tight_layout(rect=(0, 0, 1, 0.95))
    figCD.savefig(out_cd, dpi=160, bbox_inches="tight"); plt.close(figCD)


# Mechanism regions for the form-vs-placement cloud: the share of a pose's
# mean-square deviation that is conformational, r = form² / in-place². The cloud
# is a continuous, bounded (form ≤ in-place) triangle with no natural clusters
# (GMM/BIC just tiles it; HDBSCAN calls most of it noise — see
# plot_form_placement_clustering), so the meaningful partition is this physically
# grounded segmentation, not unsupervised blob-clustering.
FORM_SHARE_BINS = (1 / 3, 2 / 3)
_MECH_LABELS = ("placement-limited", "mixed", "form-limited")
_MECH_COLORS = {"placement-limited": "#2c7fb8", "mixed": "#bdbdbd",
                "form-limited": "#8856a7"}
# Pose colours for the 3-D overlay illustrations (20f/20g): "mixed" grey is too
# low-contrast against the dark crystal, so use orange there (matches the PyMOL
# .pml colouring). The 2-D scatters keep the grey _MECH_COLORS convention.
_MECH_3D_COLORS = {"placement-limited": "#2c7fb8", "mixed": "#e6820e",
                   "form-limited": "#8856a7"}


def _mechanism_region(reps: pd.DataFrame) -> pd.Series:
    """Label each pose placement-/mixed/form-limited by r = form²/in-place² —
    the fraction of its deviation that is conformation rather than positioning."""
    r = (reps["form"] ** 2 / reps["inplace"] ** 2).clip(lower=0, upper=1)
    lo, hi = FORM_SHARE_BINS
    return pd.cut(r, [-0.01, lo, hi, 1.01], labels=list(_MECH_LABELS))


def plot_form_vs_success_count(df: pd.DataFrame, out: Path,
                               form_ok: float = FORM_OK_KABSCH_A,
                               rmsd_thr: float = NEAR_NATIVE_RMSD_A,
                               selector=_near_native_valid_reps,
                               sel_note: str = _SEL_NOTE_NEAREST) -> None:
    """Per method: NUMBER of successful (≤ ``rmsd_thr`` Å AND PB-valid) poses (x,
    log scale) vs the FORM ERROR of those poses (y = median best-fit/Kabsch RMSD,
    IQR whiskers). One marker per method; asks whether tools that produce more
    successful poses also produce better-form ones. Note the axes live at
    different levels — x is a per-method count, y summarises a per-pose
    distribution — so this is a method-level view, not the per-pose panel D."""
    reps = _form_components(selector(df, rmsd_thr)).dropna(subset=["form"])
    if reps.empty:
        return
    fig, ax = plt.subplots(figsize=(9.5, 7))
    xs, tops = [], []
    for method in dict.fromkeys(reps["method"]):
        g = reps.loc[reps["method"] == method, "form"]
        if not len(g):
            continue
        x, y = len(g), float(g.median())
        lo, hi = float(g.quantile(0.25)), float(g.quantile(0.75))
        xs.append(x); tops.append(hi)
        c = TOOL_COLORS.get(method, "grey")
        ax.errorbar(x, y, yerr=[[max(0, y - lo)], [max(0, hi - y)]], fmt="o", ms=9,
                    color=c, ecolor=c, elinewidth=1.2, capsize=3, alpha=0.9)
        ax.annotate(TOOL_LABEL.get(method, method), (x, y), textcoords="offset points",
                    xytext=(7, 4), fontsize=7.5)
    ax.axhline(form_ok, ls="--", color="crimson", lw=1)
    ax.text(0.99, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
            va="bottom", ha="right", fontsize=8, transform=ax.get_yaxis_transform())
    ax.set_xlabel(f"Number of RMSD ≤ {rmsd_thr:g} Å AND PoseBusters-valid poses "
                  "(per method)")
    ax.set_ylabel("Form error of those poses — median best-fit (Kabsch) RMSD (Å)\n"
                  "(whiskers = inter-quartile range)")
    ax.set_title(sel_note, fontsize=10, fontweight="bold")
    # y and x both start at 0 so the form-error magnitudes and the per-method counts
    # read against a true baseline (linear x with the actual counts as ticks — no
    # log / 10ⁿ scientific notation, and the origin is explicit).
    ax.set_ylim(0, max(tops + [form_ok]) * 1.15 if tops else None)
    if xs:
        ax.set_xlim(0, max(xs) * 1.12)
        ax.set_xticks(sorted(set(xs)))
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ── Form fidelity by ranking depth ───────────────────────────────────────────
# A companion reading of the form-fidelity figures that asks how the internal
# conformation of a method's poses holds up as the ranking net widens: the top-1,
# top-5, top-15 and top-30 ranked poses. Rendered under TWO selections so they can
# be read side by side:
#   pbvalid — all PB-VALID poses within top-d (the RMSD ≤ 2 Å gate is DROPPED);
#   within2 — only the PB-valid poses that are ALSO RMSD ≤ 2 Å (near-native).
# Ranking is the native confidence rank for AutoDock/DiffDock and the gnina-affinity
# rank for EquiBind (no native ranking) — the same rule every other ranked figure
# uses. Cohorts are pooled per-pose and nested (top-30 ⊇ top-15 ⊇ top-5 ⊇ top-1).
# Unlike plot_form_fidelity (one representative per complex), this is a pose-level
# view, so the numbers are not comparable to form_fidelity_summary.csv.
FORM_DEPTH_COHORTS = (1, 5, 15, 30)


def _depth_palette(depths) -> dict:
    """Single-hue Blues shades (light → dark = shallow → deep) keyed by depth, so
    the grouped bars/boxes read as an accumulating ranking net for any bucket set."""
    xs = np.linspace(0.32, 0.9, len(depths)) if len(depths) > 1 else [0.7]
    return {d: plt.cm.Blues(float(x)) for d, x in zip(depths, xs)}


def _effective_rank(sub: pd.DataFrame) -> pd.Series:
    """1-based per-complex rank for the form-depth cohorts: native confidence rank
    for RANKING_TOOLS, gnina-affinity rank (:func:`_gnina_affinity_rank`) for
    EquiBind, which has no native ranking. Indexed like ``sub`` (single method)."""
    method = str(sub["method"].iloc[0]) if len(sub) else ""
    if method in RANKING_TOOLS:
        return pd.to_numeric(sub["rank"], errors="coerce")
    return _gnina_affinity_rank(sub)


def _valid_topd_poses(df: pd.DataFrame, depth: int,
                      rmsd_gate: float | None = None) -> pd.DataFrame:
    """PoseBusters-valid poses within each method's top-``depth`` ranked poses,
    pooled per-pose (up to ``depth`` rows per complex); the effective rank is kept in
    ``eff_rank`` for the per-rank panel.

    ``rmsd_gate`` — when None (default) the RMSD ≤ 2 Å gate is NOT applied (validity
    only); when a float, poses must ALSO have in-place RMSD ≤ ``rmsd_gate`` (near-native).
    """
    parts = []
    for _, sub in df.groupby("method", sort=False):
        rk = _effective_rank(sub)
        keep_mask = (rk <= depth) & sub["pb_valid"].astype(bool)
        if rmsd_gate is not None:
            keep_mask &= pd.to_numeric(sub["rmsd"], errors="coerce") <= rmsd_gate
        keep = sub[keep_mask].copy()
        keep["eff_rank"] = rk.reindex(keep.index)
        parts.append(keep)
    if not parts:
        return df.iloc[0:0]
    return pd.concat(parts)


def aggregate_form_fidelity_by_depth(df: pd.DataFrame,
                                     form_ok: float = FORM_OK_KABSCH_A,
                                     depths=FORM_DEPTH_COHORTS,
                                     rmsd_gate: float | None = None) -> pd.DataFrame:
    """Per (method, ranking depth): among the method's PB-valid poses within the
    top-``d`` ranked, how well the internal conformation ("form") reproduces the
    crystal ligand. Pose-level cohorts, nested across depths. ``rmsd_gate`` gates on
    in-place RMSD (None = validity only; 2.0 = near-native AND valid) — see
    :func:`_valid_topd_poses`.

    ``n_valid_poses`` is the pooled pose count and ``n_complexes`` the distinct
    complexes they come from. ``form_correct_%`` = best-fit ≤ ``form_ok``. With the
    RMSD gate off the cohorts include PB-valid but far-mislocated poses (a blind tool
    can eject a ligand far from the pocket — internally valid, no clashes), so the
    placement summary uses the OUTLIER-ROBUST ``inplace_median`` and a per-pose
    error-mechanism composition (``pct_placement_limited`` / ``pct_mixed`` /
    ``pct_form_limited``, from r = form²/in-place²) rather than a mean-square share.
    """
    n_pairs = df.drop_duplicates(["method", "protein", "ligand"]).groupby("method").size()
    methods = list(dict.fromkeys(df["method"]))
    rows = []
    for d in sorted({int(x) for x in depths}):
        coh_all = _form_components(
            _valid_topd_poses(df, d, rmsd_gate)).dropna(subset=["form"])
        for method in methods:
            g = coh_all[coh_all["method"] == method]
            n = len(g)
            row = {"method": method, "depth": d,
                   "n_pairs": int(n_pairs.get(method, 0)),
                   "n_valid_poses": n,
                   "n_complexes": int(g.drop_duplicates(["protein", "ligand"]).shape[0])}
            if n:
                form = g["form"].dropna()
                mech = _mechanism_region(g)
                mech_pct = mech.value_counts(normalize=True).reindex(
                    list(_MECH_LABELS)).fillna(0.0) * 100
                row.update({
                    "form_bestfit_median": round(float(form.median()), 3),
                    "form_bestfit_mean": round(float(form.mean()), 3),
                    "form_bestfit_p90": round(float(form.quantile(0.9)), 3),
                    "tfd_median": round(float(g["tfd"].median(skipna=True)), 3),
                    "form_correct_n": int((form <= form_ok).sum()),
                    "form_correct_%": round(100 * float((form <= form_ok).mean()), PCT_DECIMALS),
                    "inplace_median": round(float(g["inplace"].median(skipna=True)), 3),
                    "pct_placement_limited": round(float(mech_pct["placement-limited"]), PCT_DECIMALS),
                    "pct_mixed": round(float(mech_pct["mixed"]), PCT_DECIMALS),
                    "pct_form_limited": round(float(mech_pct["form-limited"]), PCT_DECIMALS),
                })
            rows.append(row)
    return pd.DataFrame(rows)


def _save_titled(fig, path: Path, suptitle: str, fontsize: float = 10.5) -> None:
    """Save ``fig`` with a bold, multi-line ``suptitle`` sitting snug above the axes.

    The figure must be created with ``layout='constrained'`` — constrained_layout
    reserves EXACTLY the heading's height, so there is no empty band between the title
    and the plot (the earlier ``suptitle(y=1.0)`` + fixed ``tight_layout`` rect reserved
    a fixed 10 %, which showed up as a large gap). Keep title lines short — they are
    centred over the axes, not the wider figure — so the heading never overruns the plot."""
    if suptitle:
        fig.suptitle(suptitle, fontsize=fontsize, fontweight="bold")
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_form_fidelity_by_depth(df: pd.DataFrame, out: Path,
                                form_ok: float = FORM_OK_KABSCH_A,
                                depths=FORM_DEPTH_COHORTS,
                                rmsd_gate: float | None = None,
                                suptitle: str = "") -> None:
    """Form fidelity of each method's poses as the ranking net widens to the top-1 /
    top-5 / top-15 / top-30 ranked poses, as FOUR standalone single-panel figures:

      __form_error.png   — grouped box of form error (best-fit/Kabsch RMSD), method × depth
      __form_correct.png — % of those poses whose form is correct (best-fit ≤ ``form_ok``)
      __mechanism.png    — 100 %-stacked error-mechanism composition per depth
                           (placement-/mixed/form-limited by r = form²/in-place²)
      __form_vs_rank.png — form error vs rank (median + IQR band) per method

    ``rmsd_gate`` selects the cohort: None = all PB-valid poses (RMSD gate removed);
    2.0 = only the PB-valid poses that are also RMSD ≤ 2 Å (near-native). ``suptitle``
    is the (already variant-stamped) heading drawn on every panel so each file is
    self-identifying. No-op for crystal-free sets (no form error defined).
    """
    from matplotlib.patches import Patch
    depths = sorted({int(x) for x in depths})
    agg = aggregate_form_fidelity_by_depth(df, form_ok, depths, rmsd_gate)
    if agg.empty or "n_valid_poses" not in agg or not agg["n_valid_poses"].any():
        return
    methods = [m for m in dict.fromkeys(df["method"]) if (agg["method"] == m).any()]
    cohorts = {d: _form_components(_valid_topd_poses(df, d, rmsd_gate)).dropna(subset=["form"])
               for d in depths}
    if all(c.empty for c in cohorts.values()):
        return

    def _cell(m, d, col):
        r = agg[(agg["method"] == m) & (agg["depth"] == d)]
        if r.empty or col not in r or pd.isna(r[col].iloc[0]):
            return float("nan")
        return float(r[col].iloc[0])

    n_m = len(methods)
    grp_w = 0.82
    bar_w = grp_w / len(depths)
    offs = {d: (-grp_w / 2) + bar_w * (j + 0.5) for j, d in enumerate(depths)}
    xg = np.arange(n_m)
    palette = _depth_palette(depths)
    depth_legend = [Patch(facecolor=palette[d], edgecolor="none", label=f"top-{d}")
                    for d in depths]
    method_labels = [TOOL_LABEL.get(m, m) for m in methods]
    noun = "poses" if rmsd_gate is None else "near-native poses"

    def _save(fig, tag: str):
        _save_titled(fig, out.with_stem(out.stem + tag), suptitle)

    # ── Graph 1 — form-error boxes, grouped by depth ──────────────────────────
    figE, axE = plt.subplots(figsize=(9, 6.5), layout="constrained")
    for d in depths:
        coh = cohorts[d]
        data, pos = [], []
        for i, m in enumerate(methods):
            v = coh.loc[coh["method"] == m, "form"].dropna()
            if not len(v):
                continue
            data.append(v.values); pos.append(i + offs[d])
        if not data:
            continue
        bp = axE.boxplot(data, positions=pos, widths=bar_w * 0.9,
                         patch_artist=True, showfliers=False)
        for patch in bp["boxes"]:
            patch.set_facecolor(palette[d]); patch.set_alpha(0.9)
        for med in bp["medians"]:
            med.set_color("black")
    axE.axhline(form_ok, ls="--", color="crimson", lw=1.2)
    axE.text(0.99, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
             va="bottom", ha="right", fontsize=8, transform=axE.get_yaxis_transform())
    axE.set_xticks(xg); axE.set_xticklabels(method_labels, fontsize=9)
    axE.set_xlim(-0.5, n_m - 0.5); axE.set_ylim(bottom=0)
    axE.set_ylabel(f"Form error — best-fit (Kabsch) RMSD to {_ref_noun()} (Å)")
    axE.set_title(f"How far is the form of the {noun}, by ranking depth?", fontsize=11)
    axE.grid(axis="y", alpha=0.3)
    axE.legend(handles=depth_legend, title="ranking depth", fontsize=8,
               title_fontsize=8, frameon=False, loc="upper left")
    _save(figE, "__form_error")

    # ── Graph 2 — % form-correct, grouped by depth ────────────────────────────
    figC, axPct = plt.subplots(figsize=(9, 6.5), layout="constrained")
    for d in depths:
        xs = [i + offs[d] for i in range(n_m)]
        vals = [_cell(m, d, "form_correct_%") for m in methods]
        axPct.bar(xs, [v if not np.isnan(v) else 0 for v in vals], width=bar_w * 0.9,
                  color=palette[d])
        for x, m, v in zip(xs, methods, vals):
            if not np.isnan(v):
                n_ok = int(_cell(m, d, "form_correct_n"))
                n_tot = int(_cell(m, d, "n_valid_poses"))
                axPct.text(x, v + 1, f"{v:.0f}%\n{n_ok}/{n_tot}", ha="center",
                           va="bottom", fontsize=6.2)
    axPct.set_xticks(xg); axPct.set_xticklabels(method_labels, fontsize=9)
    axPct.set_xlim(-0.5, n_m - 0.5); axPct.set_ylim(0, 108)
    axPct.set_ylabel(f"{noun.capitalize()} with correct form (%)  "
                     f"[best-fit RMSD ≤ {form_ok:g} Å]")
    axPct.set_title(f"How many {noun} get the form right, by ranking depth?\n"
                    "(labels = correct / total poses)", fontsize=11)
    axPct.grid(axis="y", alpha=0.3)
    axPct.legend(handles=depth_legend, title="ranking depth", fontsize=8,
                 title_fontsize=8, frameon=False, loc="upper right")
    _save(figC, "__form_correct")

    # ── Graph 3 — 100%-stacked error-mechanism composition per (method, depth) ─
    figM, axM = plt.subplots(figsize=(9, 6.5), layout="constrained")
    for d in depths:
        xs = [i + offs[d] for i in range(n_m)]
        seg = {"placement-limited": [_cell(m, d, "pct_placement_limited") for m in methods],
               "mixed": [_cell(m, d, "pct_mixed") for m in methods],
               "form-limited": [_cell(m, d, "pct_form_limited") for m in methods]}
        bottom = np.zeros(n_m)
        for lab in _MECH_LABELS:
            vals = np.array([v if not np.isnan(v) else 0.0 for v in seg[lab]])
            axM.bar(xs, vals, width=bar_w * 0.9, bottom=bottom,
                    color=_MECH_COLORS[lab], edgecolor="white", linewidth=0.4)
            bottom += vals
        for x in xs:
            axM.text(x, 101, f"top-{d}", ha="center", va="bottom", fontsize=6,
                     rotation=90, color="#555555")
    axM.set_xticks(xg); axM.set_xticklabels(method_labels, fontsize=9)
    axM.set_xlim(-0.5, n_m - 0.5); axM.set_ylim(0, 113)
    axM.set_ylabel(f"Composition of the {noun} by error mechanism (%)")
    tail = ("\n(the far-off valid poses are placement-limited)" if rmsd_gate is None
            else "")
    axM.set_title("Is the error placement- or form-limited, by ranking depth?\n"
                  f"(r = form²/in-place²){tail}", fontsize=10.5)
    axM.legend(handles=[Patch(facecolor=_MECH_COLORS[l], label=l) for l in _MECH_LABELS],
               fontsize=8, ncol=3, frameon=False, loc="lower center",
               bbox_to_anchor=(0.5, -0.185))
    _save(figM, "__mechanism")

    # ── Graph 4 — form error vs rank (median + IQR band) per method ────────────
    # Each curve is labelled, where it crosses a ranking-depth line (top-1/5/15/30),
    # with the % form-correct (best-fit ≤ form_ok) of that cumulative top-d cohort.
    # Ranks with < min_rank_n poses are dropped: the near-native cohort thins to a
    # handful of poses at deep ranks, where a median over 1–2 poses is noise. PB-valid
    # cohorts stay large at every rank, so this only cleans the within2 curve.
    min_rank_n = 5
    figR, axR = plt.subplots(figsize=(10, 6), layout="constrained")
    deep = cohorts[max(depths)]
    dy_of = {m: (i - (len(methods) - 1) / 2) * 12 for i, m in enumerate(methods)}

    def _fc_at(m, d):
        r = agg[(agg["method"] == m) & (agg["depth"] == d)]
        if r.empty or "form_correct_%" not in r or pd.isna(r["form_correct_%"].iloc[0]):
            return float("nan")
        return float(r["form_correct_%"].iloc[0])

    for m in methods:
        sub = deep[deep["method"] == m]
        med, lo, hi, xr = [], [], [], []
        for k in range(1, max(depths) + 1):
            v = sub.loc[sub["eff_rank"] == k, "form"].dropna()
            if len(v) < min_rank_n:
                continue
            xr.append(k); med.append(float(v.median()))
            lo.append(float(v.quantile(0.25))); hi.append(float(v.quantile(0.75)))
        if not xr:
            continue
        c = TOOL_COLORS.get(m, "grey")
        axR.fill_between(xr, lo, hi, color=c, alpha=0.15, linewidth=0)
        axR.plot(xr, med, marker="o", ms=3.5, lw=1.8, color=c, label=TOOL_LABEL.get(m, m))
        mmed = dict(zip(xr, med))
        for d in depths:
            pct = _fc_at(m, d)
            if d not in mmed or np.isnan(pct):
                continue
            dy = dy_of[m]
            axR.annotate(f"{pct:.1f}%", (d, mmed[d]), textcoords="offset points",
                         xytext=(0, dy), ha="center",
                         va="bottom" if dy >= 0 else "top", fontsize=6.6,
                         fontweight="bold", color=c,
                         bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                   ec="none", alpha=0.6))
    axR.axhline(form_ok, ls="--", color="crimson", lw=1)
    axR.text(0.99, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
             va="bottom", ha="right", fontsize=8, transform=axR.get_yaxis_transform())
    for d in depths:
        axR.axvline(d, ls=":", color="#999999", lw=0.8)
        axR.text(d, 0.965, f"top-{d}", ha="center", va="top", fontsize=7,
                 color="#555555", transform=axR.get_xaxis_transform(),
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    axR.set_xlim(0.5, max(depths) + 0.5); axR.set_ylim(bottom=0)
    axR.set_xlabel("Rank (native for AutoDock/DiffDock, gnina-affinity for EquiBind)")
    axR.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)\n(median; band = IQR)")
    axR.set_title(f"Does the form of the {noun} degrade as the ranking net widens?",
                  fontsize=11)
    axR.grid(alpha=0.3)
    axR.legend(fontsize=8, frameon=False, loc="lower left")
    if rmsd_gate is not None:
        axR.text(0.99, 0.02, f"ranks with < {min_rank_n} near-native poses omitted",
                 transform=axR.transAxes, ha="right", va="bottom", fontsize=7,
                 color="#666666", style="italic")
    _save(figR, "__form_vs_rank")


def plot_form_fidelity_depth_gate_impact(df: pd.DataFrame, out: Path,
                                         form_ok: float = FORM_OK_KABSCH_A,
                                         depths=FORM_DEPTH_COHORTS,
                                         suptitle: str = "") -> None:
    """Joint comparison of the two rank-depth selections so the IMPACT of adding the
    RMSD ≤ 2 Å criterion reads directly as the gap between paired curves, per method
    and depth:
      dashed / open marker — PB-valid only (RMSD gate removed);
      solid  / filled marker — PB-valid AND RMSD ≤ 2 Å (near-native).

    Two standalone graphs (lines only, every point value-labelled at each depth):
      __form_error.png   — median form error vs depth (the gate pulls it DOWN);
      __form_correct.png — % form-correct vs depth (the gate pushes it UP), with the
                           gate's pose-retention (near-native / all-valid) per depth.
    No-op for crystal-free sets (no RMSD defined).
    """
    from matplotlib.lines import Line2D
    depths = sorted({int(x) for x in depths})
    agg_v = aggregate_form_fidelity_by_depth(df, form_ok, depths, rmsd_gate=None)
    agg_n = aggregate_form_fidelity_by_depth(df, form_ok, depths, rmsd_gate=NEAR_NATIVE_RMSD_A)
    if agg_v.empty or "n_valid_poses" not in agg_v or not agg_v["n_valid_poses"].any():
        return
    if not ("n_valid_poses" in agg_n and agg_n["n_valid_poses"].any()):
        return
    methods = [m for m in dict.fromkeys(df["method"]) if (agg_v["method"] == m).any()]
    x = np.arange(len(depths))

    def _get(agg, m, d, col):
        r = agg[(agg["method"] == m) & (agg["depth"] == d)]
        if r.empty or col not in r or pd.isna(r[col].iloc[0]):
            return float("nan")
        return float(r[col].iloc[0])

    def _ser(agg, m, col):
        return [_get(agg, m, d, col) for d in depths]

    # Per-method horizontal nudge (points) so coincident value labels at a depth
    # (methods whose curves cross) don't overprint — the small x-shift + colour keeps
    # each label read against its own line.
    dx_of = {m: (i - (len(methods) - 1) / 2) * 11 for i, m in enumerate(methods)}

    def _labels(ax, xs, ys, color, above, fmt, dx=0.0):
        dy, va = (6, "bottom") if above else (-6, "top")
        for xi, yi in zip(xs, ys):
            if np.isnan(yi):
                continue
            ax.annotate(fmt.format(yi), (xi, yi), textcoords="offset points",
                        xytext=(dx, dy), ha="center", va=va, fontsize=6.3,
                        color=color, fontweight="bold")

    method_handles = [Line2D([], [], color=TOOL_COLORS.get(m, "grey"), lw=3,
                             label=TOOL_LABEL.get(m, m)) for m in methods]
    style_handles = [
        Line2D([], [], color="0.35", ls="-", marker="o", ms=6,
               label="PB-valid AND RMSD ≤ 2 Å (near-native)"),
        Line2D([], [], color="0.35", ls="--", marker="o", mfc="white", ms=6,
               label="PB-valid only (RMSD gate removed)")]

    def _frame(ax):
        ax.set_xticks(x); ax.set_xticklabels([f"top-{d}" for d in depths], fontsize=9)
        ax.set_xlim(-0.35, len(depths) - 0.65)
        ax.set_xlabel("Ranking depth")
        ax.grid(axis="y", alpha=0.3)
        leg = ax.legend(handles=method_handles, title="method", fontsize=8,
                        title_fontsize=8, frameon=False, loc="lower left")
        ax.add_artist(leg)
        ax.legend(handles=style_handles, title="selection", fontsize=8,
                  title_fontsize=8, frameon=False, loc="lower right")

    def _save(fig, tag):
        _save_titled(fig, out.with_stem(out.stem + tag), suptitle)

    # ── Graph 1 — impact on form error (dashed = valid-only above, solid = ≤2Å below) ──
    figE, axE = plt.subplots(figsize=(10, 6.5), layout="constrained")
    for m in methods:
        c = TOOL_COLORS.get(m, "grey")
        yv, yn = _ser(agg_v, m, "form_bestfit_median"), _ser(agg_n, m, "form_bestfit_median")
        axE.plot(x, yv, ls="--", marker="o", mfc="white", ms=6, lw=1.6, color=c)
        axE.plot(x, yn, ls="-", marker="o", ms=6, lw=2.2, color=c)
        _labels(axE, x, yv, c, above=True, fmt="{:.2f}", dx=dx_of[m])
        _labels(axE, x, yn, c, above=False, fmt="{:.2f}", dx=dx_of[m])
    axE.axhline(form_ok, ls=":", color="crimson", lw=1)
    axE.text(0.01, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
             va="bottom", ha="left", fontsize=8, transform=axE.get_yaxis_transform())
    top = np.nanmax([_get(agg_v, m, d, "form_bestfit_median")
                     for m in methods for d in depths] + [form_ok])
    axE.set_ylim(0, top * 1.2)
    axE.set_ylabel(f"Form error — median best-fit (Kabsch) RMSD to {_ref_noun()} (Å)")
    axE.set_title("Median form error at each rank (value = Å)", fontsize=10.5)
    _frame(axE)
    _save(figE, "__form_error")

    # ── Graph 2 — impact on % form-correct (solid above, dashed below) + retention ──
    figC, axC = plt.subplots(figsize=(10, 6.5), layout="constrained")
    for m in methods:
        c = TOOL_COLORS.get(m, "grey")
        pv, pn = _ser(agg_v, m, "form_correct_%"), _ser(agg_n, m, "form_correct_%")
        axC.plot(x, pv, ls="--", marker="o", mfc="white", ms=6, lw=1.6, color=c)
        axC.plot(x, pn, ls="-", marker="o", ms=6, lw=2.2, color=c)
        _labels(axC, x, pn, c, above=True, fmt="{:.0f}%", dx=dx_of[m])
        _labels(axC, x, pv, c, above=False, fmt="{:.0f}%", dx=dx_of[m])
    for j, d in enumerate(depths):
        n_v = sum(int(_get(agg_v, m, d, "n_valid_poses") or 0) for m in methods)
        n_n = sum(int(_get(agg_n, m, d, "n_valid_poses") or 0) for m in methods)
        if n_v:
            axC.text(j, 112, f"≤2Å kept\n{n_n}/{n_v}\n({100*n_n/n_v:.0f}%)", ha="center",
                     va="bottom", fontsize=6.5, color="#444444")
    axC.set_ylim(0, 128)
    axC.set_ylabel(f"Poses with correct form (%)  [best-fit RMSD ≤ {form_ok:g} Å]")
    axC.set_title("% form-correct at each rank ('kept' = poses passing the ≤ 2 Å gate)",
                  fontsize=10.5)
    _frame(axC)
    _save(figC, "__form_correct")


def _stats_form_fidelity_gate_vs_rank(
        df: pd.DataFrame,
        form_ok: float = FORM_OK_KABSCH_A,
        depths=(1, 5, 10, 15),
        rmsd_gate: float = NEAR_NATIVE_RMSD_A,
        min_poses_corr: int = 3) -> "dict | None":
    """Paired significance tests for fig 20 ``…_gate_vs_rank`` (the three-panel
    form-fidelity-vs-rank companion). Crystal-only; returns ``None`` for crystal-free
    sets or when no method carries a form value.

    The unit of analysis is ALWAYS one ``(protein, ligand)`` complex — poses are
    aggregated to a per-complex representative or a per-complex summary BEFORE any
    test, so the figure's pooled-pose curves and top-d labels are never fed raw (no
    pose-level pseudo-replication). Four tests, each tied to a visual claim:

      ``gate_effect`` — (Panel C) does the RMSD ≤ 2 Å gate buy better internal
        geometry? The near-native set is a strict SUBSET of the PB-valid set, so a
        direct subset-vs-superset contrast is degenerate. Instead each tool's PB-valid
        poses are PARTITIONED into near-native (≤ ``rmsd_gate``) vs far (> gate); per
        complex holding both, the median best-fit (Kabsch) form RMSD of each group is
        compared with a paired Wilcoxon signed-rank + Hodges–Lehmann median-difference
        CI (far − near; > 0 ⇒ the gate removes worse-form poses).
      ``form_across_tools_{pbvalid,near_native}`` — (Panels A / B) do tools differ in
        the form error of the pose the ranker hands you first? Representative = the
        SHALLOWEST-ranked PB-valid (A) / near-native (B) pose per complex; Friedman +
        Kendall's W + pairwise Wilcoxon/Holm across the listwise-complete complexes.
        The comparison is conditional on each tool HAVING such a pose, so n shrinks to
        the shared set (reported as ``n_complete``).
      ``rank_trend`` — (all panels) does form error degrade with rank? Per complex,
        Kendall's τ of form vs effective rank over its PB-valid poses (≥
        ``min_poses_corr``); one-sample Wilcoxon of the per-complex τ vs 0 (repeated
        ranks within a complex are correlated, so the trend is tested on one slope per
        complex, never pooled). EquiBind (gnina-affinity rank, not a native confidence
        rank) is the pre-specified negative control.
      ``form_placement_coupling`` — (mechanism) is best-fit form even coupled to
        in-place RMSD? Per complex Spearman(form, in-place) over its PB-valid poses;
        one-sample Wilcoxon of the ρ vs 0. ρ ≈ 0 ⇒ the gate (an in-place criterion)
        can barely move form, so Panel C is small by construction.

    Each per-tool family (gate_effect / rank_trend / coupling) is Holm-corrected across
    the tools; the cross-tool panels carry their own Holm pairwise (via
    :func:`stats_utils.paired_continuous`). Full payload → the JSON sidecar; the figure
    shows only the headline lines.
    """
    max_d = max(int(x) for x in depths)
    coh = _form_components(_valid_topd_poses(df, max_d, None)).dropna(subset=["form"])
    if coh.empty:
        return None
    coh = coh.copy()
    coh["inplace_rmsd"] = pd.to_numeric(coh["inplace"], errors="coerce")
    coh["_near"] = pd.to_numeric(coh["rmsd"], errors="coerce") <= rmsd_gate
    methods = [m for m in dict.fromkeys(df["method"]) if (coh["method"] == m).any()]
    if not methods:
        return None

    def _holm_over(recs):
        """Add Holm-adjusted p + star to the records in ``recs`` that carry p_raw."""
        idx = [i for i, r in enumerate(recs) if "p_raw" in r]
        if not idx:
            return
        for i, pa in zip(idx, su.holm([recs[i]["p_raw"] for i in idx])):
            recs[i]["p_holm"] = float(pa)
            recs[i]["star"] = su.p_stars(pa)

    # ── gate_effect (Panel C): near vs far PB-valid poses, per-complex median form ──
    gate = []
    for m in methods:
        g = coh[coh["method"] == m]
        near = g[g["_near"]].groupby(["protein", "ligand"])["form"].median()
        far = g[~g["_near"]].groupby(["protein", "ligand"])["form"].median()
        common = near.index.intersection(far.index)
        rec = {"method": m, "n_complexes_both": int(len(common))}
        if len(common) >= _MIN_UNITS_STATS:
            nn = near.reindex(common).to_numpy(float)   # near-native (≤ gate)
            fn = far.reindex(common).to_numpy(float)    # far (> gate)
            rb, p, npair = su.wilcoxon_rankbiserial(fn, nn)
            est, lo, hi = su.median_diff_ci(fn, nn, paired=True)
            rec.update({
                "median_form_near": round(float(np.median(nn)), 3),
                "median_form_far": round(float(np.median(fn)), 3),
                "hl_median_diff_far_minus_near": round(float(est), 3),
                "hl_ci": [round(float(lo), 3), round(float(hi), 3)],
                "rank_biserial": round(float(rb), 3), "p_raw": float(p),
                "n": int(npair)})
        else:
            rec["status"] = "n too small — exploratory"
        gate.append(rec)
    _holm_over(gate)

    # ── form_across_tools (Panels A/B): shallowest-ranked rep per complex, paired ──
    def _across_tools(cohort):
        reps = {}
        for m in methods:
            sub = cohort[cohort["method"] == m]
            if sub.empty:
                continue
            reps[m] = (sub.sort_values("eff_rank")
                       .groupby(["protein", "ligand"]).first()["form"])
        labs = [m for m in methods if m in reps]
        if len(labs) < 2:
            return None
        common = sorted(set.intersection(*[set(reps[m].index) for m in labs]),
                        key=lambda t: (str(t[0]), str(t[1])))
        n = len(common)
        out = {"unit": "shallowest-ranked pose per complex; paired across tools",
               "methods": labs, "n_complete": n}
        if n >= _MIN_UNITS_STATS:
            data = {m: reps[m].reindex(common).to_numpy(float) for m in labs}
            res = su.paired_continuous(data, labels=labs)
            out.update({"omnibus": res["omnibus"], "pairwise": res["pairwise"],
                        "medians": {k: round(v, 3) for k, v in res["medians"].items()}})
        else:
            out["status"] = "n too small — exploratory"
            out["medians"] = {m: round(float(reps[m].reindex(common).median()), 3)
                              for m in labs}
        return out

    # ── rank_trend + form_placement_coupling: one slope/correlation per complex ──
    def _per_complex_corr(kind):
        recs = []
        for m in methods:
            g = coh[coh["method"] == m]
            coefs = []
            for _, cx in g.groupby(["protein", "ligand"]):
                if len(cx) < min_poses_corr:
                    continue
                if kind == "rank":
                    t, _, _ = su.kendall(cx["eff_rank"].to_numpy(float),
                                         cx["form"].to_numpy(float))
                else:
                    t, _, _ = su.spearman(cx["form"].to_numpy(float),
                                          cx["inplace_rmsd"].to_numpy(float))
                if t == t:                       # drop NaN (constant series)
                    coefs.append(float(t))
            rec = {"method": m, "n_complexes": len(coefs)}
            if len(coefs) >= _MIN_UNITS_STATS:
                arr = np.asarray(coefs, float)
                rb, p, npair = su.wilcoxon_rankbiserial(arr, np.zeros_like(arr))
                rec.update({"median_coef": round(float(np.median(arr)), 3),
                            "rank_biserial": round(float(rb), 3),
                            "p_raw": float(p), "n": int(npair)})
            else:
                rec["status"] = "n too small — exploratory"
            recs.append(rec)
        _holm_over(recs)
        return recs

    return {
        "figure": "20_form_fidelity__rank_depth_gate_vs_rank.png",
        "metric": f"best-fit (Kabsch) form RMSD of each tool's PB-valid poses "
                  f"within top-{max_d}",
        "unit": "one (protein, ligand) complex (poses aggregated per complex first)",
        "rmsd_gate_A": float(rmsd_gate), "form_ok_A": float(form_ok),
        "top_depth": max_d, "min_poses_corr": int(min_poses_corr),
        "methods": list(methods),
        "gate_effect": gate,
        "form_across_tools_pbvalid": _across_tools(coh),
        "form_across_tools_near_native": _across_tools(coh[coh["_near"]]),
        "rank_trend": _per_complex_corr("rank"),
        "form_placement_coupling": _per_complex_corr("coupling"),
    }


def _annotate_form_gate_vs_rank_stats(axA, axB, axC, payload, methods) -> None:
    """Overlay the headline paired tests on the three panels of fig 20's
    ``…_gate_vs_rank`` figure (full pairwise / per-tool detail is in the JSON
    sidecar). Panel A: cross-tool form-error omnibus + the per-complex rank trend;
    Panel B: cross-tool omnibus on the near-native reps; Panel C: the gate effect per
    tool + the form⊥placement coupling that frames it."""
    if not payload:
        return
    tl = lambda m: TOOL_LABEL.get(m, m)

    def _omni(block):
        if not block or "omnibus" not in block:
            return None
        om = block["omnibus"]
        w = om.get("kendall_w")
        return (f"Friedman {su.fmt_p(om.get('p'))}"
                + (f", W={w:.2f}" if isinstance(w, (int, float)) and w == w else "")
                + f" · n={block.get('n_complete', '?')}")

    def _trend_str(recs, sym):
        out = []
        by = {r["method"]: r for r in recs}
        for m in methods:
            r = by.get(m)
            if r and "median_coef" in r:
                out.append(f"{tl(m)} {sym}={r['median_coef']:+.2f}{r.get('star', '')}")
        return "; ".join(out)

    # ── Panel A — cross-tool form error (PB-valid) + rank trend ──
    A = payload.get("form_across_tools_pbvalid")
    linesA = ["Cross-tool form error — shallowest PB-valid pose / complex"]
    if _omni(A):
        linesA.append("  " + _omni(A))
    if A and A.get("pairwise"):
        linesA.append("  " + "; ".join(
            f"{tl(p['a'])}–{tl(p['b'])} {p.get('star', '')}" for p in A["pairwise"]))
    trend = _trend_str(payload.get("rank_trend", []), "τ")
    if trend:
        linesA.append(f"Rank trend (Kendall τ, form vs rank / complex vs 0): {trend}")
    axA.text(0.985, 0.03, "\n".join(linesA), transform=axA.transAxes, ha="right",
             va="bottom", fontsize=6.8, color="0.20", linespacing=1.3,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9), zorder=7)

    # ── Panel B — cross-tool form error on the near-native reps ──
    B = payload.get("form_across_tools_near_native")
    linesB = ["Cross-tool form error — shallowest near-native pose / complex"]
    if _omni(B):
        linesB.append("  " + _omni(B))
    if B and B.get("pairwise"):
        linesB.append("  " + "; ".join(
            f"{tl(p['a'])}–{tl(p['b'])} {p.get('star', '')}" for p in B["pairwise"]))
    axB.text(0.015, 0.975, "\n".join(linesB), transform=axB.transAxes, ha="left",
             va="top", fontsize=6.8, color="0.20", linespacing=1.3,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9), zorder=7)

    # ── Panel C — gate effect per tool (the headline) + form⊥placement framing ──
    linesC = ["Gate effect — median form: far (>2 Å) − near (≤2 Å) PB-valid, per complex"]
    by_g = {r["method"]: r for r in payload.get("gate_effect", [])}
    for m in methods:
        r = by_g.get(m)
        if r and "hl_median_diff_far_minus_near" in r:
            ci = r.get("hl_ci", [float("nan"), float("nan")])
            linesC.append(f"  {tl(m)} Δ={r['hl_median_diff_far_minus_near']:+.2f} Å "
                          f"[{ci[0]:+.2f}, {ci[1]:+.2f}] {r.get('star', '')} "
                          f"(n={r.get('n', '?')})")
    cpl = _trend_str(payload.get("form_placement_coupling", []), "ρ")
    if cpl:
        linesC.append(f"Form⊥placement: Spearman(form, in-place) / complex: {cpl}")
    axC.text(0.015, 0.975, "\n".join(linesC), transform=axC.transAxes, ha="left",
             va="top", fontsize=6.8, color="0.20", linespacing=1.3,
             bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9), zorder=7)


def plot_form_fidelity_gate_vs_rank(df: pd.DataFrame, out: Path,
                                    form_ok: float = FORM_OK_KABSCH_A,
                                    depths=(1, 5, 10, 15),
                                    min_n: int = 5,
                                    suptitle: str = "",
                                    stats: "dict | None" = None) -> None:
    """Per-rank companion to :func:`plot_form_fidelity_depth_gate_impact` — three
    panels stacked on a shared rank x-axis (1 … max(depths); default top-15):

      (A) form error (median best-fit/Kabsch RMSD + IQR band) vs rank for the
          PB-VALID poses — the RMSD ≤ 2 Å gate removed;
      (B) the same for the NEAR-NATIVE subset — PB-valid AND RMSD ≤ 2 Å;
      (C) the gate's effect at each rank = median(A) − median(B) per method — how
          many Å the ≤ 2 Å criterion pulls the median form error down.

    At each ranking-depth line (top-d) on each method's curve, the panels label:
      A / B — % form-correct (best-fit ≤ ``form_ok``) and ``n`` = distinct complexes
              of that cumulative top-d cohort;
      C     — the gate's lift in form-correct (near-native − PB-valid, percentage
              points) and ``n`` = near-native complexes.

    ``stats`` (from :func:`_stats_form_fidelity_gate_vs_rank`) adds the per-complex
    paired tests: cross-tool form-error omnibus on panels A/B, the per-complex rank
    trend, and the gate-effect (near vs far) + form⊥placement coupling on panel C.
    The figure shows only the headline lines; the full payload is in the JSON sidecar.
    No-op for crystal-free sets (no RMSD defined)."""
    from matplotlib.lines import Line2D
    depths = sorted({int(x) for x in depths})
    max_d = max(depths)
    coh_v = _form_components(_valid_topd_poses(df, max_d, None)).dropna(subset=["form"])
    coh_n = _form_components(
        _valid_topd_poses(df, max_d, NEAR_NATIVE_RMSD_A)).dropna(subset=["form"])
    if coh_v.empty:
        return
    methods = [m for m in dict.fromkeys(df["method"]) if (coh_v["method"] == m).any()]
    agg_v = aggregate_form_fidelity_by_depth(df, form_ok, depths, rmsd_gate=None)
    agg_n = aggregate_form_fidelity_by_depth(df, form_ok, depths,
                                             rmsd_gate=NEAR_NATIVE_RMSD_A)

    def _col(agg, m, d, col):
        r = agg[(agg["method"] == m) & (agg["depth"] == d)]
        if r.empty or col not in r or pd.isna(r[col].iloc[0]):
            return float("nan")
        return float(r[col].iloc[0])

    omitted = False

    def _per_rank(coh, m):
        """{rank: median/q25/q75} of form error for method m; ranks with < ``min_n``
        poses are dropped (the near-native cohort thins to a handful of poses at deep
        ranks, where a median over 1–2 poses is noise, not signal)."""
        nonlocal omitted
        sub = coh[coh["method"] == m]
        med, lo, hi = {}, {}, {}
        for k in range(1, max_d + 1):
            v = sub.loc[sub["eff_rank"] == k, "form"].dropna()
            if len(v) >= min_n:
                med[k] = float(v.median())
                lo[k] = float(v.quantile(0.25)); hi[k] = float(v.quantile(0.75))
            elif len(v):
                omitted = True
        return med, lo, hi

    def _edge(d):
        if d == depths[0]:
            return "left", 6
        if d == depths[-1]:
            return "right", -6
        return "center", 0

    def _stack_labels(ax, items, d, spacing=27):
        """``items`` = list of (y, colour, text). Anchor each label to its own curve
        at rank ``d`` but spread them vertically in y-order (lowest curve → label
        below, highest → above) so coincident/crossing curves never overprint."""
        items = sorted(items, key=lambda t: t[0])
        n = len(items)
        ha, hx = _edge(d)
        for rank, (y, color, text) in enumerate(items):
            dy = (rank - (n - 1) / 2) * spacing
            ax.annotate(text, (d, y), textcoords="offset points", xytext=(hx, dy),
                        ha=ha, va="center", fontsize=6.4, fontweight="bold", color=color,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                                  alpha=0.7))

    med_of = {"v": {}, "n": {}}  # per-method {rank: median} for the gap panel

    fig, (axA, axB, axC) = plt.subplots(3, 1, figsize=(11, 14), sharex=True,
                                        layout="constrained")

    def _draw(ax, coh, agg, key, label_line=True):
        top = form_ok
        for m in methods:
            c = TOOL_COLORS.get(m, "grey")
            med, lo, hi = _per_rank(coh, m)
            med_of[key][m] = med
            xr = sorted(med)
            if not xr:
                continue
            ax.fill_between(xr, [lo[k] for k in xr], [hi[k] for k in xr],
                            color=c, alpha=0.13, linewidth=0)
            ax.plot(xr, [med[k] for k in xr], marker="o", ms=3.2, lw=1.8, color=c)
            top = max(top, max(hi.values()))
        for d in depths:
            items = []
            for m in methods:
                med = med_of[key].get(m, {})
                pct = _col(agg, m, d, "form_correct_%")
                if d not in med or np.isnan(pct):
                    continue
                n = _col(agg, m, d, "n_complexes")
                items.append((med[d], TOOL_COLORS.get(m, "grey"),
                              f"{pct:.1f}%\nn={int(n)}"))
            _stack_labels(ax, items, d)
        ax.axhline(form_ok, ls="--", color="crimson", lw=1)
        if label_line:  # label the shared reference line once (panel A) to avoid clashes
            ax.text(0.012, form_ok, f"correct-form line ({form_ok:g} Å)", color="crimson",
                    va="bottom", ha="left", fontsize=7.5,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.7),
                    transform=ax.get_yaxis_transform())
        for d in depths:
            ax.axvline(d, ls=":", color="#999999", lw=0.8)
            ax.text(d, 0.98, f"top-{d}", ha="center", va="top", fontsize=7,
                    color="#555555", transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        ax.set_ylim(0, top * 1.16)
        ax.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)\n(median; band = IQR)")
        ax.grid(alpha=0.3)
        ax.tick_params(labelbottom=False)

    _draw(axA, coh_v, agg_v, "v")
    axA.set_title("(A) PB-valid poses — RMSD ≤ 2 Å gate removed   "
                  "(labels = % form-correct · n complexes at top-d)",
                  fontsize=10.5, loc="left")
    _draw(axB, coh_n, agg_n, "n", label_line=False)
    axB.set_title("(B) Near-native poses — PB-valid AND RMSD ≤ 2 Å   "
                  "(labels = % form-correct · n complexes at top-d)",
                  fontsize=10.5, loc="left")
    if omitted:
        axB.text(0.995, 0.03, f"ranks with < {min_n} near-native poses omitted",
                 transform=axB.transAxes, ha="right", va="bottom", fontsize=7,
                 color="#666666", style="italic")

    # ── Panel C — the gate's per-rank effect: median(PB-valid) − median(near-native) ──
    hi_c = 0.0
    gap_of = {}
    for m in methods:
        c = TOOL_COLORS.get(m, "grey")
        mv, mn = med_of["v"].get(m, {}), med_of["n"].get(m, {})
        xr = sorted(set(mv) & set(mn))
        if not xr:
            continue
        gap_of[m] = {k: mv[k] - mn[k] for k in xr}
        hi_c = max(hi_c, max(gap_of[m].values()))
        axC.plot(xr, [gap_of[m][k] for k in xr], marker="o", ms=3.2, lw=1.8, color=c)
    for d in depths:
        items = []
        for m in methods:
            if d not in gap_of.get(m, {}):
                continue
            pp = _col(agg_n, m, d, "form_correct_%") - _col(agg_v, m, d, "form_correct_%")
            n = _col(agg_n, m, d, "n_complexes")
            items.append((gap_of[m][d], TOOL_COLORS.get(m, "grey"),
                          f"{pp:+.1f}%\nn={int(n)}"))
        _stack_labels(axC, items, d)
    axC.axhline(0, color="#444444", lw=0.9)
    for d in depths:
        axC.axvline(d, ls=":", color="#999999", lw=0.8)
        axC.text(d, 0.98, f"top-{d}", ha="center", va="top", fontsize=7,
                 color="#555555", transform=axC.get_xaxis_transform(),
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    axC.set_ylim(0, hi_c * 1.25 + 0.05)
    axC.set_ylabel("Gate effect on form error (Å)\n"
                   "median(PB-valid) − median(near-native)")
    axC.set_title("(C) Gate effect per rank — curve = Å form-error removed;  "
                  "labels = form-correct lift (%) · n complexes", fontsize=10.5, loc="left")
    axC.grid(alpha=0.3)
    axC.set_xlim(0.5, max_d + 0.5)
    axC.set_xlabel("Rank (native for AutoDock/DiffDock, gnina-affinity for EquiBind)")

    handles = [Line2D([], [], color=TOOL_COLORS.get(m, "grey"), marker="o", lw=1.8,
                      label=TOOL_LABEL.get(m, m)) for m in methods]
    fig.legend(handles=handles, loc="outside lower center", ncol=len(methods),
               frameon=False, fontsize=9)
    if stats:
        try:
            _annotate_form_gate_vs_rank_stats(axA, axB, axC, stats, methods)
        except Exception as _exc:
            print(f"  WARNING: fig 20 gate-vs-rank stats annotation failed ({_exc})")
    _save_titled(fig, out, suptitle)


def _draw_form_placement_axes(ax, sub: pd.DataFrame, form_ok: float,
                              rmsd_thr: float, lim: float,
                              axis_fontsize: float | None = None,
                              tick_fontsize: float | None = None) -> None:
    """Scatter of in-place RMSD (x) vs form/best-fit RMSD (y) for ``sub``, points
    coloured by mechanism region, with the y = x bound and the r = 1/3, 2/3
    mechanism rays (form = √r · in-place). ``axis_fontsize``/``tick_fontsize``
    (default None → matplotlib default) enlarge the axis titles and tick labels."""
    reg = _mechanism_region(sub)
    for label in _MECH_LABELS:
        s = sub[reg == label]
        if len(s):
            ax.scatter(s["inplace"], s["form"], s=16, alpha=0.5,
                       color=_MECH_COLORS[label], edgecolors="none", label=label)
    ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
    for t in FORM_SHARE_BINS:
        ax.plot([0, lim], [0, math.sqrt(t) * lim], ls=":", color="black", lw=0.7)
    ax.axhline(form_ok, ls="--", color="crimson", lw=0.9)
    ax.axvline(rmsd_thr, ls=":", color="black", lw=0.6)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel(f"In-place RMSD to {_ref_noun()} (Å)", fontsize=axis_fontsize)
    ax.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)", fontsize=axis_fontsize)
    if tick_fontsize is not None:
        ax.tick_params(axis="both", labelsize=tick_fontsize)
    ax.grid(alpha=0.3)


def plot_form_vs_placement_by_family(df: pd.DataFrame, out: Path,
                                     form_ok: float = FORM_OK_KABSCH_A,
                                     rmsd_thr: float = NEAR_NATIVE_RMSD_A,
                                     selector=_near_native_valid_reps,
                                     sel_note: str = _SEL_NOTE_NEAREST) -> None:
    """Panel D of fig 20 split into individual per-family graphs, each coloured by
    mechanism region (placement-limited / mixed / form-limited). Points below the
    lower ray have the shape right and are merely mis-positioned; points near the
    y = x bound carry nearly all their error in the conformation itself."""
    reps = _form_components(selector(df, rmsd_thr)).dropna(subset=["form", "inplace"])
    if reps.empty:
        return
    lim = float(max(reps["inplace"].max(), reps["form"].max(), rmsd_thr)) * 1.05
    fams = [("autodock", "AutoDock"), ("diffdock", "DiffDock"), ("equibind", "EquiBind")]
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    ax_list = axes.ravel()
    for (fam, label), ax in zip(fams, ax_list[:3]):
        sub = reps[reps["method"].map(_fam_key) == fam]
        if sub.empty:
            ax.set_visible(False); continue
        _draw_form_placement_axes(ax, sub, form_ok, rmsd_thr, lim)
        ax.set_title(f"{label}  (n={len(sub)} successful poses)", fontsize=11)
        ax.legend(fontsize=8, frameon=False, loc="upper left", title="mechanism")
    axD = ax_list[3]
    _draw_form_placement_axes(axD, reps, form_ok, rmsd_thr, lim)
    axD.set_title(f"All tools combined  (n={len(reps)})", fontsize=11)
    axD.legend(fontsize=8, frameon=False, loc="upper left", title="mechanism")
    axD.text(lim * 0.98, lim * 0.02,
             "below lower ray: shape right, position off\n"
             "above upper ray: shape itself wrong",
             fontsize=7.5, ha="right", va="bottom", color="#333333")
    _label_panels(ax_list)
    fig.suptitle(_vt("Form vs placement of the successful poses, split by tool "
                     "family\ncoloured by mechanism (r = form² / in-place²): "
                     f"placement-limited · mixed · form-limited\n{sel_note}"),
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


FAR_FROM_RECEPTOR_CENTROID_A = 8.0  # a docked pose whose heavy-atom centroid sits
# farther than this from the crystal-ligand site is treated as ejected off the
# receptor (a different/decoy pocket). 8 Å is the repo's "same binding-site" scale
# (POCKET_RADIUS in pose_cluster_crystal_pocket_report.py / equibind_unguided_
# centroid_clusters.py) and sits in the density antimode between the in-pocket peak
# (0–2 Å) and the ejected blind-docking mode (>10 Å) — see centroid_dist (:833-837).


def plot_form_vs_placement_depth_filmstrip(
        df: pd.DataFrame, out: Path,
        form_ok: float = FORM_OK_KABSCH_A,
        rmsd_gate: float | None = NEAR_NATIVE_RMSD_A,
        depths=(1, 3, 5),
        centroid_max: float | None = None,
        axis_max: float | None = None,
        caption_out: Path | None = None,
        legend_loc: str = "panel",
        compact_titles: bool = False,
        title_mechanism_line: bool = True,
        title_fontsize: float = 9.0,
        axis_fontsize: float | None = None,
        tick_fontsize: float | None = None,
        legend_fontsize: float = 11.0,
        legend_title_fontsize: float = 12.0,
        suptitle_fontsize: float = 13.0,
        shared_axis_labels: bool = False) -> None:
    """Fig 20d read across ranking depth: rows = tool family, columns = the
    cumulative top-``d`` cohorts (top-1 ⊆ top-3 ⊆ top-5). Every cell is the same
    mechanism-coloured form-vs-placement scatter on SHARED axes, so the clouds are
    directly comparable across the row. A black ✕ marks each cohort's centroid
    (median in-place, median form) and an orange arrow runs from the top-1 centroid,
    making the DRIFT as the ranking net widens legible within each row: rightward =
    deeper ranks add more mis-positioned (placement-limited) valid poses; upward =
    the typical form degrades toward the y = x bound. The arrow's *shortness* is
    itself the result when widening the net barely moves the cloud. Each cell title
    carries ``n`` = PB-valid poses in the cohort and ``rec.`` = the distinct receptors
    (unique protein–ligand complexes) those poses span, so pose count and receptor
    breadth are read apart (a cohort can pile many poses onto few receptors).

    ``rmsd_gate`` gates the cohorts on in-place RMSD. Default 2 Å keeps the near-
    native 0–2 Å frame of the oracle 20d figures (the 'within2' reading). None =
    all PB-valid poses (RMSD gate removed), which by themselves sprawl to tens of Å.

    ``centroid_max`` (used with ``rmsd_gate=None``) removes poses that sit outside the
    crystal-site neighbourhood (they are NOT off the receptor — nearly all stay in van
    der Waals contact with the protein; they are simply away from the crystal site):
    poses whose docked centroid is > ``centroid_max`` Å from the crystal-ligand site
    (or has no centroid) are dropped, matching the code's out-of-pocket convention
    (``centroid_dist``). This trims the ejected blind-docking mode so the axes stay
    readable while the informative near-to-moderate spread the ≤ 2 Å gate hides is
    kept. No-op for crystal-free sets (no in-place RMSD / centroid defined).

    ``axis_max`` clips both axes to a fixed square (e.g. 8 Å) after the frame is
    otherwise computed, so the dense near-native core stays legible; poses that
    fall outside the frame are counted in the footnote, not silently dropped.

    ``caption_out`` (default None) — when a path is given, the descriptive notes are
    written to that text file (with the figure title as a heading) INSTEAD of being
    drawn as an on-figure footnote, and the freed bottom band is reclaimed by the
    panels. None keeps the footnote on the figure (unchanged default behaviour).

    ``legend_loc`` — "panel" (default) draws the mechanism legend inside the top-left
    panel; "top" draws a single figure-level horizontal legend below the title.
    ``compact_titles`` — False (default) keeps the full per-panel stat titles
    (n / rec. / place-mix-form %); True reduces each panel title to just
    "Family · top-d" and leaves the numbers to the sidecar tables.
    ``title_mechanism_line`` — True (default) keeps the second panel-title line
    ("place/mix/form = …%"); False drops it (keeping line 1: family · top-d ·
    n · rec.), for callers that move the mechanism breakdown into a sidecar table.
    Ignored when ``compact_titles`` is set.
    ``title_fontsize`` / ``axis_fontsize`` / ``tick_fontsize`` size the panel headings,
    axis titles and tick labels (None → matplotlib default for the latter two).
    ``legend_fontsize`` / ``legend_title_fontsize`` size the figure-level mechanism
    legend and its title (only used with ``legend_loc="top"``; the in-panel legend
    keeps its own compact size). ``suptitle_fontsize`` sizes the overall figure title.
    Defaults reproduce the historic hardcoded sizes so existing callers are unchanged.
    ``shared_axis_labels`` — False (default) keeps a per-panel x/y axis label on the
    outer panels; True clears them all and draws ONE figure-level x and y label
    (``fig.supxlabel`` / ``supylabel``) at ``axis_fontsize``. Recommended when
    ``axis_fontsize`` is large enough that the long per-panel y-label would exceed a
    single panel's height and overflow into the neighbouring row.
    """
    depths = sorted({int(d) for d in depths})
    if not depths:
        return
    fams = [("autodock", "AutoDock"), ("diffdock", "DiffDock"), ("equibind", "EquiBind")]
    raw = {d: _form_components(_valid_topd_poses(df, d, rmsd_gate)
                               ).dropna(subset=["form", "inplace"]) for d in depths}
    n_dropped = 0
    if centroid_max is not None:
        deepest_before = len(raw[depths[-1]])
        cohorts = {d: (c[c["centroid_dist"].le(centroid_max)]
                       if "centroid_dist" in c.columns else c) for d, c in raw.items()}
        n_dropped = deepest_before - len(cohorts[depths[-1]])
    else:
        cohorts = raw
    if all(c.empty for c in cohorts.values()):
        return
    vlim = rmsd_gate if rmsd_gate is not None else NEAR_NATIVE_RMSD_A
    raw_lim = max((max(c["inplace"].max(), c["form"].max())
                   for c in cohorts.values() if not c.empty), default=vlim)
    # Gated: snug 0–2 Å frame. Gate removed: round up to an integer-Å square frame
    # (the retained near-site poses still spread to ~15 Å in in-place RMSD).
    lim = (float(math.ceil(max(raw_lim, vlim))) if rmsd_gate is None
           else float(max(raw_lim, vlim)) * 1.05)
    # ``axis_max`` clips the frame to a fixed square so the dense near-native core is
    # readable; poses beyond it are counted for the footnote rather than dropped.
    n_clipped = 0
    if axis_max is not None:
        lim = float(axis_max)
        deep = cohorts[depths[-1]]
        n_clipped = int((deep["inplace"].gt(lim) | deep["form"].gt(lim)).sum())

    def _centroid(sub):
        return ((float(sub["inplace"].median()), float(sub["form"].median()))
                if len(sub) else None)

    d1 = depths[0]
    nrows, ncols = len(fams), len(depths)
    fig, axes = plt.subplots(nrows, ncols, figsize=(4.7 * ncols, 4.4 * nrows),
                             squeeze=False)
    for ri, (fam, flabel) in enumerate(fams):
        base = cohorts[d1][cohorts[d1]["method"].map(_fam_key) == fam]
        c1, n1 = _centroid(base), len(base)
        for ci, d in enumerate(depths):
            ax = axes[ri][ci]
            sub = cohorts[d][cohorts[d]["method"].map(_fam_key) == fam]
            if sub.empty:
                ax.set_visible(False); continue
            _draw_form_placement_axes(ax, sub, form_ok, vlim, lim,
                                      axis_fontsize=axis_fontsize,
                                      tick_fontsize=tick_fontsize)
            if ci > 0:
                ax.set_ylabel("")
            if ri < nrows - 1:
                ax.set_xlabel("")
            cc = _centroid(sub)
            if cc:
                if c1 and d != d1 and (abs(cc[0] - c1[0]) + abs(cc[1] - c1[1])) > 1e-3:
                    ax.annotate("", xy=cc, xytext=c1,
                                arrowprops=dict(arrowstyle="-|>", color="#e6820e",
                                                lw=2.0, alpha=0.95))
                ax.scatter(*cc, marker="X", s=90, color="black", zorder=6,
                           edgecolors="white", linewidths=1.0)
            if compact_titles:
                # Numbers (n / rec. / place-mix-form %) move to the sidecar tables;
                # the panel heading keeps only the tool family and ranking depth.
                ax.set_title(f"{flabel} · top-{d}", fontsize=title_fontsize,
                             fontweight="bold")
            else:
                n_rec = int(sub["protein"].nunique())  # distinct receptors (unique complexes)
                dn = len(sub) - n1
                extra = f"  (+{dn} poses vs top-{d1})" if d != d1 and dn else ""
                line1 = f"{flabel} · top-{d}   n={len(sub)} · {n_rec} rec.{extra}"
                if title_mechanism_line:
                    # Mechanism composition of this cohort (share of poses that are
                    # placement-limited / mixed / form-limited) — reading across a row
                    # shows the share shift as the ranking net widens. Callers that move
                    # this to a sidecar table pass title_mechanism_line=False.
                    _mc = _mechanism_region(sub).value_counts(normalize=True).reindex(
                        list(_MECH_LABELS)).fillna(0.0) * 100
                    ax.set_title(line1 + "\nplace/mix/form = "
                                 f"{_mc['placement-limited']:.0f}/"
                                 f"{_mc['mixed']:.0f}/{_mc['form-limited']:.0f}%",
                                 fontsize=title_fontsize)
                else:
                    ax.set_title(line1, fontsize=title_fontsize)
            if legend_loc == "panel" and ri == 0 and ci == 0:
                ax.legend(fontsize=7.5, frameon=False, loc="upper left", title="mechanism")
    if shared_axis_labels:
        # A single figure-level x/y label instead of one per outer panel: the long
        # y-label is taller than one panel at large fonts, so per-panel labels would
        # overflow into the neighbouring row. One shared label spans the whole grid.
        _lab_fs = axis_fontsize if axis_fontsize is not None else 12.0
        for ax in axes.ravel():
            if ax.get_visible():
                ax.set_xlabel(""); ax.set_ylabel("")
        fig.supxlabel(f"In-place RMSD to {_ref_noun()} (Å)", fontsize=_lab_fs)
        fig.supylabel("Form error — best-fit (Kabsch) RMSD (Å)", fontsize=_lab_fs)
    if rmsd_gate is not None:
        sel_line = ("cumulative PB-valid poses within each method's top-d "
                    f"(≤ {vlim:g} Å, near-native)")
    else:
        sel_line = "cumulative PB-valid poses within each method's top-d (RMSD gate removed)"

    # Short title; the descriptive detail that used to be stacked in it now sits in a
    # footnote so the heading rides just above the panels instead of floating off.
    title = _vt("Form vs placement across ranking depth — top-1 net widened to "
                f"top-{'/'.join(f'{d}' for d in depths[1:])}")
    fig.suptitle(title, fontsize=suptitle_fontsize, fontweight="bold", y=0.995)
    # ``legend_loc="top"`` lifts the mechanism key out of the top-left panel into a
    # single horizontal legend below the title; ``top`` then reserves that band so the
    # panels drop clear of it.
    top = 0.97
    if legend_loc == "top":
        from matplotlib.lines import Line2D
        handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=9,
                          markerfacecolor=_MECH_COLORS[l], markeredgecolor="none",
                          label=l) for l in _MECH_LABELS]
        fig.legend(handles=handles, title="mechanism", loc="upper center",
                   bbox_to_anchor=(0.5, 0.965), ncol=len(handles), frameon=False,
                   fontsize=legend_fontsize, title_fontsize=legend_title_fontsize,
                   handletextpad=0.4, columnspacing=1.6)
        top = 0.93
    notes = [f"Rows = tool family · columns = {sel_line}."]
    if rmsd_gate is None and centroid_max is not None:
        # NOT "off-receptor": the excluded poses were verified to sit in van der Waals
        # contact with the protein (autodock_gnina: 5,915 excluded PB-valid poses, median
        # closest protein contact 2.84 Å, max 3.53 Å). They are away from the CRYSTAL SITE,
        # which is a different claim — the caption must not imply they left the receptor.
        notes.append(f"Poses outside the crystal-site neighbourhood removed — docked centroid "
                     f"> {centroid_max:g} Å from the crystal-ligand site, not off the receptor "
                     f"({n_dropped} dropped at top-{depths[-1]}).")
    if n_clipped:
        notes.append(f"Axes clipped at {lim:g} Å for legibility — {n_clipped} pose(s) beyond "
                     f"the frame at top-{depths[-1]} not shown.")
    notes.append("Titles: n = PB-valid poses in the cohort; rec. = distinct receptors "
                 "(unique protein–ligand complexes) those poses span.")
    notes.append("Coloured by mechanism (r = form² / in-place²); black ✕ = cohort centroid "
                 "(median in-place, median form); orange arrow = drift from top-1.")
    if caption_out is not None:
        # Move the descriptive notes off the figure into a sidecar text file and let
        # the panels reclaim the band the footnote used to occupy.
        Path(caption_out).write_text(title + "\n\n" + "\n".join(notes) + "\n",
                                     encoding="utf-8")
        rect = (0, 0.0, 1, top)
    else:
        fig.text(0.5, 0.008, "\n".join(notes), ha="center", va="bottom",
                 fontsize=8, color="0.30", linespacing=1.35)
        rect = (0, 0.06, 1, top)
    fig.tight_layout(rect=rect)
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ── Honest statistics companion to the _pbvalid depth filmstrip ──────────────
# The filmstrip is a qualitative read of three confounds that its scatter cannot
# quantify on its own, so this ships the numbers alongside it:
#   1. COVERAGE — each per-tool cell is a DIFFERENT set of complexes (a tool only
#      appears where it produced a valid near-site pose), so raw cross-tool
#      medians are not apples-to-apples. We report valid-complex coverage and run
#      the cross-tool comparison PAIRED on the complexes all three tools cover.
#   2. r-COUPLING — r = form²/in-place² is inflated by low placement error, so a
#      high median r ("form-limited") can just mean placement is solved; we report
#      ABSOLUTE best-fit (form) RMSD as the primary form metric so r is never used
#      to rank tools on conformation.
#   3. POSE DIVERSITY — a flat depth profile can be genuine ranking robustness OR
#      near-duplicate poses (low sample diversity); the within-complex spread of a
#      tool's top-d in-place RMSD tells the two apart.
# Cohorts are reconstructed identically to the _pbvalid filmstrip (top-d PB-valid,
# RMSD gate removed, centroid-trimmed to the crystal-site neighbourhood). Crystal-free
# sets → {} (no-op).
_FAM_STATS = [("autodock", "AutoDock", "#1f77b4"),
              ("diffdock", "DiffDock", "#ff7f0e"),
              ("equibind", "EquiBind", "#2ca02c")]


def _holm(pvals: list[float]) -> list[float]:
    """Holm–Bonferroni step-down adjusted p-values (input order preserved)."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [0.0] * m
    run = 0.0
    for rank, idx in enumerate(order):
        run = max(run, (m - rank) * pvals[idx])
        adj[idx] = min(1.0, run)
    return adj


def _cliffs_delta(a, b) -> float:
    """Cliff's δ effect size: P(a>b) − P(a<b) ∈ [−1, 1]."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    if a.size == 0 or b.size == 0:
        return float("nan")
    gt = sum(np.sum(a > y) for y in b)
    lt = sum(np.sum(a < y) for y in b)
    return float((gt - lt) / (a.size * b.size))


def _sig_star(p: float) -> str:
    if p != p:            # NaN
        return ""
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"


def _boot_centroid_drift(sub1: pd.DataFrame, sub5: pd.DataFrame,
                         b_iter: int = 2000, seed: int = 0) -> dict:
    """Cluster (receptor) bootstrap 95 % CI on the top-1 → top-d centroid drift
    (Δmedian in-place, Δmedian form). Resamples RECEPTORS with replacement — the
    unit of independence — not poses, so the CI respects the pseudoreplication of
    many poses per complex."""
    def bank(sub):
        return {r: (g["inplace"].to_numpy(float), g["form"].to_numpy(float))
                for r, g in sub.groupby("protein")}
    b1, b5 = bank(sub1), bank(sub5)
    recs = np.array(sorted(set(b1) | set(b5)))
    if recs.size == 0:
        nan = float("nan")
        return dict(dx=nan, dy=nan, dx_ci=[nan, nan], dy_ci=[nan, nan])
    rng = np.random.default_rng(seed)
    dxs, dys = [], []
    for _ in range(b_iter):
        pick = rng.choice(recs, size=recs.size, replace=True)
        p1 = [r for r in pick if r in b1]
        p5 = [r for r in pick if r in b5]
        if not p1 or not p5:
            continue
        ip1 = np.concatenate([b1[r][0] for r in p1]); fm1 = np.concatenate([b1[r][1] for r in p1])
        ip5 = np.concatenate([b5[r][0] for r in p5]); fm5 = np.concatenate([b5[r][1] for r in p5])
        dxs.append(np.median(ip5) - np.median(ip1))
        dys.append(np.median(fm5) - np.median(fm1))
    dx = float(sub5["inplace"].median() - sub1["inplace"].median())
    dy = float(sub5["form"].median() - sub1["form"].median())
    return dict(dx=dx, dy=dy,
                dx_ci=[float(np.percentile(dxs, 2.5)), float(np.percentile(dxs, 97.5))],
                dy_ci=[float(np.percentile(dys, 2.5)), float(np.percentile(dys, 97.5))])


def _filmstrip_cohorts(df: pd.DataFrame, depths, centroid_max: float) -> dict:
    """The exact cohorts behind the _pbvalid filmstrip: PB-valid poses within each
    method's top-d (RMSD gate removed), outside-the-crystal-site (centroid > ``centroid_max``)
    trimmed, with ``fam``/``r`` columns added. ``eff_rank`` is carried for banding."""
    out = {}
    for d in sorted({int(x) for x in depths}):
        c = _form_components(_valid_topd_poses(df, d, None)).dropna(subset=["form", "inplace"])
        if "centroid_dist" in c.columns:
            c = c[c["centroid_dist"].le(centroid_max)]
        c = c.copy()
        c["fam"] = c["method"].map(_fam_key)
        c["r"] = (c["form"] ** 2 / c["inplace"] ** 2).clip(lower=0, upper=1)
        out[d] = c
    return out


def aggregate_filmstrip_statistics(df: pd.DataFrame, depths=(1, 3, 5),
                                   centroid_max: float = FAR_FROM_RECEPTOR_CENTROID_A,
                                   form_ok: float = FORM_OK_KABSCH_A) -> dict:
    """Statistics behind the _pbvalid depth filmstrip. Returns a dict of tidy
    DataFrames (``{}`` for crystal-free sets, or if SciPy is unavailable):

      per_tool_depth   — one row per (tool, depth): cumulative in-place & form
                         median/IQR, mechanism composition (%), median r, and
                         valid-complex COVERAGE of the benchmark universe.
      within_tool_trend— one row per (tool, metric∈{in-place, form}): the disjoint
                         rank-band medians (r1 / r2-3 / r4-5), a pose-level Spearman
                         rank-vs-value trend, the within-complex PAIRED Friedman
                         (the clustering-robust depth test), and the top-1→top-d
                         centroid drift with a receptor-bootstrap 95 % CI.
      crosstool_paired — one row per (depth, metric∈{in-place, form, r}): the tools
                         compared on the COMMON complexes only (Friedman + Wilcoxon
                         signed-rank post-hoc, Holm-adjusted, Cliff's δ).
      pose_diversity   — one row per tool at the deepest depth: within-complex spread
                         of in-place RMSD (mode-collapse vs genuine sample diversity).
    """
    try:
        from scipy import stats as ss
    except Exception as e:  # pragma: no cover
        print(f"  [filmstrip-stats] SciPy unavailable ({e}); skipping companion stats.")
        return {}
    depths = sorted({int(x) for x in depths})
    cohorts = _filmstrip_cohorts(df, depths, centroid_max)
    if all(c.empty for c in cohorts.values()):
        return {}  # crystal-free set — no in-place RMSD defined
    dmax = depths[-1]
    universe = int(df["protein"].nunique())
    covsets = {(f, d): set(cohorts[d].loc[cohorts[d]["fam"] == f, "protein"])
               for f, _, _ in _FAM_STATS for d in depths}

    # ── per (tool, depth): cumulative composition, medians, coverage ──
    rows = []
    for f, lab, _ in _FAM_STATS:
        for d in depths:
            sub = cohorts[d][cohorts[d]["fam"] == f]
            if sub.empty:
                continue
            mech = _mechanism_region(sub).astype(str)
            ncplx = int(sub["protein"].nunique())
            rows.append({
                "tool": lab, "depth": d, "n_poses": int(len(sub)),
                "n_valid_complexes": ncplx, "universe_complexes": universe,
                "coverage_pct": round(100 * ncplx / universe, PCT_DECIMALS) if universe else float("nan"),
                "inplace_median": round(float(sub["inplace"].median()), 3),
                "inplace_q1": round(float(sub["inplace"].quantile(0.25)), 3),
                "inplace_q3": round(float(sub["inplace"].quantile(0.75)), 3),
                "form_median": round(float(sub["form"].median()), 3),
                "form_q1": round(float(sub["form"].quantile(0.25)), 3),
                "form_q3": round(float(sub["form"].quantile(0.75)), 3),
                "median_r": round(float(sub["r"].median()), 3),
                "pct_placement_limited": round(100 * float((mech == "placement-limited").mean()), PCT_DECIMALS),
                "pct_mixed": round(100 * float((mech == "mixed").mean()), PCT_DECIMALS),
                "pct_form_limited": round(100 * float((mech == "form-limited").mean()), PCT_DECIMALS),
            })
    per_tool_depth = pd.DataFrame(rows)

    # ── within-tool depth trend (disjoint bands + paired Friedman + drift) ──
    top = cohorts[dmax]; base = cohorts[depths[0]]
    band_names = ["r1", "r2-3", "r4-5"]
    rows2, rows_div = [], []
    for f, lab, _ in _FAM_STATS:
        s5 = top[top["fam"] == f].copy()
        if s5.empty:
            continue
        er = pd.to_numeric(s5["eff_rank"], errors="coerce")
        s5["band"] = np.where(er <= 1, "r1", np.where(er <= 3, "r2-3", "r4-5"))
        s1 = base[base["fam"] == f]
        drift = _boot_centroid_drift(s1, s5)
        for metric in ("inplace", "form"):
            bmed = {b: (float(s5.loc[s5["band"] == b, metric].median())
                        if (s5["band"] == b).any() else float("nan")) for b in band_names}
            sr = ss.spearmanr(er, pd.to_numeric(s5[metric], errors="coerce"), nan_policy="omit")
            piv = s5.pivot_table(index="protein", columns="band", values=metric, aggfunc="median")
            have = [b for b in band_names if b in piv.columns]
            comp = piv[have].dropna() if len(have) == 3 else pd.DataFrame()
            fr_chi = fr_p = float("nan"); n_comp = 0
            if len(comp) >= 5:
                fr = ss.friedmanchisquare(*[comp[b].values for b in have])
                fr_chi, fr_p, n_comp = float(fr.statistic), float(fr.pvalue), int(len(comp))
            dv, ci = ((drift["dx"], drift["dx_ci"]) if metric == "inplace"
                      else (drift["dy"], drift["dy_ci"]))
            rows2.append({
                "tool": lab, "metric": metric,
                "band_r1_median": round(bmed["r1"], 3),
                "band_r2_3_median": round(bmed["r2-3"], 3),
                "band_r4_5_median": round(bmed["r4-5"], 3),
                "spearman_rank_rho": round(float(sr.statistic), 3),
                "spearman_rank_p": float(sr.pvalue),
                "friedman_chi2_paired": round(fr_chi, 2) if fr_chi == fr_chi else fr_chi,
                "friedman_p_paired": fr_p, "friedman_n_complexes": n_comp,
                f"drift_top{depths[0]}_to_top{dmax}": round(dv, 3),
                "drift_ci_lo": round(ci[0], 3), "drift_ci_hi": round(ci[1], 3),
            })
        # within-complex pose spread at the deepest depth (mode-collapse check)
        npose = s5.groupby("protein")["inplace"].size()
        multi = npose[npose >= 2].index
        sm = s5[s5["protein"].isin(multi)]
        spread = sm.groupby("protein")["inplace"].agg(lambda x: x.max() - x.min())
        sd = sm.groupby("protein")["inplace"].std()
        rows_div.append({
            "tool": lab, "depth": dmax,
            "n_complexes_ge2_valid_poses": int(len(multi)),
            "median_poses_per_complex": round(float(npose.median()), 1) if len(npose) else float("nan"),
            "median_within_complex_inplace_spread": round(float(spread.median()), 2) if len(spread) else float("nan"),
            "median_within_complex_inplace_sd": round(float(sd.median()), 2) if len(sd.dropna()) else float("nan"),
        })
    within_tool_trend = pd.DataFrame(rows2)
    pose_diversity = pd.DataFrame(rows_div)

    # ── cross-tool comparison on the COMMON complex set (paired) ──
    pair_idx = [(0, 1), (0, 2), (1, 2)]
    labels = [lab for _, lab, _ in _FAM_STATS]
    rows3 = []
    for d in depths:
        coh = cohorts[d]
        inter = set.intersection(*[covsets[(f, d)] for f, _, _ in _FAM_STATS])
        for metric in ("inplace", "form", "r"):
            wide = {}
            for f, lab, _ in _FAM_STATS:
                s = coh[(coh["fam"] == f) & (coh["protein"].isin(inter))]
                wide[lab] = s.groupby("protein")[metric].median()
            W = pd.DataFrame(wide)
            W = (W.loc[sorted(inter)].dropna() if inter else W.dropna())
            rec = {"depth": d, "metric": metric, "n_common": int(len(W))}
            for lab in labels:
                rec[f"median_{lab}"] = round(float(W[lab].median()), 3) if len(W) else float("nan")
            if len(W) >= 5:
                fr = ss.friedmanchisquare(*[W[lab].values for lab in labels])
                rec["friedman_chi2"] = round(float(fr.statistic), 2)
                rec["friedman_p"] = float(fr.pvalue)
                praw = [ss.wilcoxon(W[labels[i]].values, W[labels[j]].values).pvalue
                        for i, j in pair_idx]
                for (i, j), ph in zip(pair_idx, _holm(praw)):
                    rec[f"p_{labels[i]}_vs_{labels[j]}"] = ph
                    rec[f"cliffs_{labels[i]}_vs_{labels[j]}"] = round(
                        _cliffs_delta(W[labels[i]].values, W[labels[j]].values), 3)
            rows3.append(rec)
    crosstool_paired = pd.DataFrame(rows3)

    return {"per_tool_depth": per_tool_depth, "within_tool_trend": within_tool_trend,
            "crosstool_paired": crosstool_paired, "pose_diversity": pose_diversity}


def plot_filmstrip_statistics(stats: dict, out: Path,
                              depths=(1, 3, 5),
                              form_ok: float = FORM_OK_KABSCH_A) -> None:
    """Four-panel honest-statistics companion to the _pbvalid depth filmstrip.

    (A) valid-complex COVERAGE per tool (why the per-cell medians are not on a
        common set); (B) the depth drift of each tool's centroid in the
        form-vs-placement plane (placement degrades, form far less); (C) the
        cross-tool comparison PAIRED on the common complexes (tie at top-1, then
        DiffDock separates by top-d); (D) within-complex pose spread — whether a
        flat depth profile is ranking robustness or near-duplicate poses."""
    if not stats:
        return
    ptd = stats["per_tool_depth"]; xt = stats["crosstool_paired"]; pv = stats["pose_diversity"]
    depths = sorted({int(x) for x in depths})
    d0, dmax = depths[0], depths[-1]
    labels = [lab for _, lab, _ in _FAM_STATS]
    colors = {lab: col for _, lab, col in _FAM_STATS}
    present = [lab for lab in labels if (ptd["tool"] == lab).any()]
    x = np.arange(len(present))

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 11))
    (axA, axB), (axC, axD) = axes

    # (A) coverage — grouped bars, top-d0 vs top-dmax
    cov_depths = sorted({d0, dmax})
    w = 0.8 / len(cov_depths)
    for k, d in enumerate(cov_depths):
        vals = [float(ptd[(ptd["tool"] == lab) & (ptd["depth"] == d)]["coverage_pct"].iloc[0])
                if ((ptd["tool"] == lab) & (ptd["depth"] == d)).any() else 0.0 for lab in present]
        ncx = [int(ptd[(ptd["tool"] == lab) & (ptd["depth"] == d)]["n_valid_complexes"].iloc[0])
               if ((ptd["tool"] == lab) & (ptd["depth"] == d)).any() else 0 for lab in present]
        bars = axA.bar(x + (k - (len(cov_depths) - 1) / 2) * w, vals, w,
                       label=f"top-{d}", color=[colors[l] for l in present],
                       alpha=0.55 + 0.45 * k, edgecolor="white")
        for b, nc in zip(bars, ncx):
            axA.text(b.get_x() + b.get_width() / 2, b.get_height() + 1, str(nc),
                     ha="center", va="bottom", fontsize=8)
    uni = int(ptd["universe_complexes"].iloc[0])
    axA.set_ylim(0, 105); axA.set_xticks(x); axA.set_xticklabels(present)
    axA.set_ylabel("Valid-complex coverage (% of benchmark)")
    axA.set_title(f"(A) Coverage — each tool is scored on a different subset\n"
                  f"bars labelled with n valid complexes (universe = {uni})", fontsize=10.5)
    axA.legend(title="depth", fontsize=8, frameon=False, loc="upper right")
    axA.grid(axis="y", alpha=0.3)

    # (B) depth drift of the centroid in the form-vs-placement plane
    lim = 0.0
    for lab in present:
        g = ptd[ptd["tool"] == lab]
        lim = max(lim, g["inplace_median"].max(), g["form_median"].max())
    lim = float(max(lim * 1.15, 2.0))
    axB.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
    for t in FORM_SHARE_BINS:
        axB.plot([0, lim], [0, math.sqrt(t) * lim], ls=":", color="black", lw=0.7)
    axB.axhline(form_ok, ls="--", color="crimson", lw=0.9)
    for lab in present:
        g = ptd[ptd["tool"] == lab].sort_values("depth")
        xs, ys = g["inplace_median"].to_numpy(), g["form_median"].to_numpy()
        axB.plot(xs, ys, "-", color=colors[lab], lw=1.4, alpha=0.8)
        axB.scatter(xs[0], ys[0], s=44, color=colors[lab], zorder=5)      # top-d0 (dot)
        axB.scatter(xs[1:-1], ys[1:-1], s=22, color=colors[lab], zorder=5)
        axB.scatter(xs[-1], ys[-1], s=130, marker="X", color=colors[lab],
                    edgecolors="white", linewidths=1.0, zorder=6, label=lab)
    axB.set_xlim(0, lim); axB.set_ylim(0, lim)
    axB.set_xlabel(f"In-place RMSD to {_ref_noun()} — median (Å)")
    axB.set_ylabel("Form error — best-fit (Kabsch) RMSD, median (Å)")
    axB.set_title(f"(B) Depth drift of each tool's centroid (top-{d0} → top-{dmax})\n"
                  "rightward = placement degrades · upward = form degrades", fontsize=10.5)
    axB.legend(fontsize=8, frameon=False, loc="upper left",
               title=f"● = top-{d0}  ✕ = top-{dmax}")
    axB.grid(alpha=0.3)

    # (C) cross-tool, PAIRED on the common complexes — median in-place
    ci = xt[xt["metric"] == "inplace"]
    w = 0.8 / len(cov_depths)
    ymax = 0.0
    for k, d in enumerate(cov_depths):
        row = ci[ci["depth"] == d]
        if row.empty:
            continue
        row = row.iloc[0]
        vals = [float(row.get(f"median_{lab}", float("nan"))) for lab in present]
        ymax = max([ymax] + [v for v in vals if v == v])
        bars = axC.bar(x + (k - (len(cov_depths) - 1) / 2) * w, vals, w,
                       color=[colors[l] for l in present], alpha=0.55 + 0.45 * k,
                       edgecolor="white", label=f"top-{d} (n={int(row['n_common'])})")
        pad = row.get(f"p_{present[0]}_vs_{present[1]}", float("nan")) if len(present) >= 2 else float("nan")
        if pad == pad and len(present) >= 2:
            xa = x[0] + (k - (len(cov_depths) - 1) / 2) * w
            xb = x[1] + (k - (len(cov_depths) - 1) / 2) * w
            yb = max(vals[0], vals[1]) + ymax * 0.06
            axC.plot([xa, xa, xb, xb], [yb - ymax * 0.02, yb, yb, yb - ymax * 0.02],
                     lw=0.9, color="0.3")
            axC.text((xa + xb) / 2, yb, _sig_star(pad), ha="center", va="bottom", fontsize=9)
    axC.set_xticks(x); axC.set_xticklabels(present)
    axC.set_ylim(0, ymax * 1.3 if ymax else 1)
    axC.set_ylabel("Median in-place RMSD on shared complexes (Å)")
    _bracket = (f"bracket = {present[0]} vs {present[1]} (Wilcoxon, Holm)"
                if len(present) >= 2 else "Wilcoxon signed-rank, Holm-adjusted")
    axC.set_title(f"(C) Cross-tool PAIRED on common complexes\n{_bracket}", fontsize=10.5)
    axC.legend(fontsize=8, frameon=False, loc="upper left")
    axC.grid(axis="y", alpha=0.3)

    # (D) within-complex pose spread — mode-collapse check
    sp = [float(pv[pv["tool"] == lab]["median_within_complex_inplace_spread"].iloc[0])
          if (pv["tool"] == lab).any() else 0.0 for lab in present]
    sd = [float(pv[pv["tool"] == lab]["median_within_complex_inplace_sd"].iloc[0])
          if (pv["tool"] == lab).any() else 0.0 for lab in present]
    bars = axD.bar(x, sp, 0.6, color=[colors[l] for l in present], edgecolor="white")
    for b, s in zip(bars, sd):
        axD.text(b.get_x() + b.get_width() / 2, b.get_height(), f" SD {s:.2f}",
                 ha="center", va="bottom", fontsize=8)
    axD.set_xticks(x); axD.set_xticklabels(present)
    axD.set_ylabel(f"Within-complex spread of in-place RMSD\nacross top-{dmax} valid poses — median (Å)")
    axD.set_title(f"(D) Pose diversity within a complex (top-{dmax})\n"
                  "small spread = near-duplicate poses (flat ≠ better ranking)", fontsize=10.5)
    axD.grid(axis="y", alpha=0.3)

    _label_panels(np.array([axA, axB, axC, axD]))
    fig.suptitle(_vt("Honest statistics for the form-vs-placement depth filmstrip "
                     "(PB-valid, gate removed)\ncoverage · paired cross-tool test · "
                     "depth drift · pose diversity"),
                 fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_filmstrip_mechanism_share(stats: dict, out: Path, depths=(1, 3, 5)) -> None:
    """How the mechanism COMPOSITION of each tool's PB-valid cohort shifts as the
    ranking net widens (top-1 → top-3 → top-5). One panel per tool; within it a
    100 %-stacked bar per depth splits the cohort into placement-limited / mixed /
    form-limited (r = form²/in-place² bins), with the pose count n above each bar.
    Reads directly off ``per_tool_depth`` (the same %s stamped on the filmstrip
    panels). The deeper ranks add mis-positioned-but-valid poses, so the
    placement-limited share grows — steeply for AutoDock/EquiBind, barely for a
    tool whose ranking is depth-stable. ``{}`` / missing table → no-op."""
    if not stats or "per_tool_depth" not in stats:
        return
    ptd = stats["per_tool_depth"]
    if ptd.empty:
        return
    depths = sorted({int(x) for x in depths})
    labels = [lab for _, lab, _ in _FAM_STATS]
    present = [lab for lab in labels if (ptd["tool"] == lab).any()]
    if not present:
        return
    seg_cols = ["pct_placement_limited", "pct_mixed", "pct_form_limited"]
    fig, axes = plt.subplots(1, len(present), figsize=(4.4 * len(present), 5.0),
                             squeeze=False)
    axes = axes[0]
    for ax, lab in zip(axes, present):
        g = ptd[ptd["tool"] == lab].set_index("depth")
        ds = [d for d in depths if d in g.index]
        x = np.arange(len(ds))
        bottoms = np.zeros(len(ds))
        for seg, mech in zip(seg_cols, _MECH_LABELS):
            vals = np.array([float(g.loc[d, seg]) for d in ds])
            ax.bar(x, vals, 0.62, bottom=bottoms, color=_MECH_COLORS[mech],
                   edgecolor="white", label=mech)
            for xi, (v, b) in enumerate(zip(vals, bottoms)):
                if v >= 6:                       # skip labels on slivers
                    ax.text(xi, b + v / 2, f"{v:.0f}", ha="center", va="center",
                            fontsize=8.5, color="white" if mech != "mixed" else "0.25")
            bottoms += vals
        for xi, d in enumerate(ds):
            ax.text(xi, 101, f"n={int(g.loc[d, 'n_poses'])}", ha="center", va="bottom",
                    fontsize=8, color="0.35")
        ax.set_xticks(x); ax.set_xticklabels([f"top-{d}" for d in ds])
        ax.set_ylim(0, 108); ax.set_xlim(-0.6, len(ds) - 0.4)
        ax.set_title(lab, fontsize=11)
        if ax is axes[0]:
            ax.set_ylabel("Share of PB-valid poses (%)")
        ax.grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=8, frameon=False, loc="lower left", title="mechanism")
    fig.suptitle(_vt("Mechanism composition shifts as the ranking net widens "
                     "(top-1 → top-3 → top-5)\nshare of PB-valid poses that are "
                     "placement-limited · mixed · form-limited (r = form² / in-place²)"),
                 fontsize=12.5, fontweight="bold", y=1.02)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_form_placement_clustering(df: pd.DataFrame, out: Path,
                                   rmsd_thr: float = NEAR_NATIVE_RMSD_A,
                                   selector=_near_native_valid_reps,
                                   sel_note: str = _SEL_NOTE_NEAREST) -> None:
    """Why unsupervised blob-clustering is the wrong tool for the panel-D cloud:
    GMM (BIC-selected k) tiles a continuous gradient, k-means splits it by
    magnitude, HDBSCAN flags most points as noise. Rendered to justify the
    recommended mechanism-based segmentation. Skipped if scikit-learn is absent."""
    try:
        from sklearn.preprocessing import StandardScaler
        from sklearn.mixture import GaussianMixture
        from sklearn.cluster import KMeans, HDBSCAN
        from sklearn.metrics import silhouette_score
    except Exception as e:  # pragma: no cover
        print(f"  [clustering] scikit-learn unavailable ({e}); skipping 20e.")
        return
    reps = _form_components(selector(df, rmsd_thr)).dropna(subset=["form", "inplace"])
    if len(reps) < 20:
        return
    X = reps[["inplace", "form"]].to_numpy()
    Xs = StandardScaler().fit_transform(X)
    lim = float(max(X[:, 0].max(), X[:, 1].max(), rmsd_thr)) * 1.05
    bics = {k: GaussianMixture(k, covariance_type="full", random_state=0,
                               n_init=3).fit(Xs).bic(Xs) for k in range(1, 7)}
    kbest = min(bics, key=bics.get)
    gmm_lab = GaussianMixture(kbest, covariance_type="full", random_state=0,
                              n_init=3).fit(Xs).predict(Xs)
    km2 = KMeans(2, n_init=10, random_state=0).fit_predict(Xs)
    sil = silhouette_score(Xs, km2)
    hdb = HDBSCAN(min_cluster_size=30).fit(Xs).labels_
    noise = 100 * float((hdb == -1).mean())

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    cmap = plt.get_cmap("tab10")

    def _scat(ax, lab, title):
        for i, u in enumerate(sorted(set(lab))):
            s = reps[lab == u]
            col = "#cccccc" if u == -1 else cmap(i % 10)
            name = "noise" if u == -1 else f"cluster {u}"
            ax.scatter(s["inplace"], s["form"], s=14, alpha=0.5, color=col,
                       edgecolors="none", label=name)
        ax.plot([0, lim], [0, lim], ls="--", color="grey", lw=1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.grid(alpha=0.3)
        ax.set_xlabel("In-place RMSD (Å)")
        ax.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)")
        ax.set_title(title, fontsize=10.5)
        ax.legend(fontsize=7, frameon=False, loc="upper left")

    _scat(axes[0], gmm_lab, f"GMM · BIC picks k={kbest}\ncomponents tile a "
                            "continuous gradient (not real clusters)")
    _scat(axes[1], km2, f"k-means k=2 · silhouette {sil:.2f}\nsplits the cloud "
                        "by magnitude, no natural gap")
    _scat(axes[2], hdb, f"HDBSCAN · {noise:.0f}% flagged as noise\nno "
                        "density-separated clusters exist")
    _label_panels(axes)
    fig.suptitle(_vt("Clustering the form-vs-placement cloud: there are no natural "
                     "clusters — it is a bounded gradient\n(so segment by "
                     "mechanism ratio r = form²/in-place², not by blob-clustering)"
                     f"\n{sel_note}"),
                 fontsize=12.5, fontweight="bold", y=1.05)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def aggregate_oracle_selection_comparison(df: pd.DataFrame,
                                          rmsd_thr: float = NEAR_NATIVE_RMSD_A) -> pd.DataFrame:
    """Per method, how many complexes count as a success under each oracle
    selection, and how many the strict (nearest) rule needlessly drops.

    Columns:
      current_nearest_rule    — n complexes whose CLOSEST pose is ≤ rmsd_thr Å AND
                                PB-valid (the default figures' population).
      relaxed_any_valid_le2A  — n complexes with ANY ≤ rmsd_thr Å PB-valid pose
                                (the validity-constrained ceiling).
      rescuable               — relaxed − current: complexes whose closest pose is
                                ≤ rmsd_thr Å but INVALID while another ≤ rmsd_thr Å
                                valid pose exists (the only way the two rules differ,
                                since a >2 Å nearest pose means every pose is > 2 Å).
    """
    df = df.dropna(subset=["rmsd"])
    near = _near_native_valid_reps(df, rmsd_thr)
    valid = _best_valid_near_native_reps(df, rmsd_thr)
    n_near = near.groupby("method").size() if len(near) else pd.Series(dtype=int)
    n_valid = valid.groupby("method").size() if len(valid) else pd.Series(dtype=int)
    rows = []
    for method in dict.fromkeys(df["method"]):
        a, b = int(n_near.get(method, 0)), int(n_valid.get(method, 0))
        rows.append({"method": method, "current_nearest_rule": a,
                     "relaxed_any_valid_le2A": b, "rescuable": b - a})
    out = pd.DataFrame(rows)
    return (out[out["relaxed_any_valid_le2A"] > 0]
            .sort_values("current_nearest_rule", ascending=False)
            .reset_index(drop=True))


def plot_oracle_selection_comparison(df: pd.DataFrame, out: Path,
                                     form_ok: float = FORM_OK_KABSCH_A,
                                     rmsd_thr: float = NEAR_NATIVE_RMSD_A) -> None:
    """Head-to-head of the two oracle selections that feed the form-fidelity
    figures. (A) complexes counted as a success per method under each rule
    (nearest vs validity-constrained), rescued count annotated in red; (B) the
    median form error of those reps under each rule — do the rescued, slightly-
    farther valid poses change the typical shape quality?"""
    df = df.dropna(subset=["rmsd"])
    near = _form_components(_near_native_valid_reps(df, rmsd_thr))
    valid = _form_components(_best_valid_near_native_reps(df, rmsd_thr))
    if near.empty and valid.empty:
        return
    methods = [m for m in dict.fromkeys(df["method"])
               if (near["method"] == m).any() or (valid["method"] == m).any()]
    x = np.arange(len(methods)); w = 0.38
    n_near = [int((near["method"] == m).sum()) for m in methods]
    n_val = [int((valid["method"] == m).sum()) for m in methods]

    def _med(frame, m):
        v = frame.loc[frame["method"] == m, "form"].dropna()
        return float(v.median()) if len(v) else float("nan")
    f_near = [_med(near, m) for m in methods]
    f_val = [_med(valid, m) for m in methods]
    c_near, c_val = "#4477aa", "#66ccee"

    # Split into two standalone graphs: success counts, and median form error.
    out_counts = out.with_stem(out.stem + "_success_counts")
    out_form = out.with_stem(out.stem + "_form_error")
    suptitle = _vt("Oracle selection: RMSD-greedy (nearest pose) vs "
                   "validity-constrained (nearest ≤ 2 Å valid pose)")
    figw = max(8, 1.05 * len(methods))

    # ── Graph 1 — complexes qualifying under each selection ─────────────────
    figA, axA = plt.subplots(figsize=(figw, 6))
    axA.bar(x - w / 2, n_near, w, color=c_near, label="nearest rule (RMSD-greedy oracle)")
    axA.bar(x + w / 2, n_val, w, color=c_val,
            label="validity-constrained (nearest ≤ 2 Å valid pose)")
    for xi, a, b in zip(x, n_near, n_val):
        if b - a > 0:
            axA.annotate(f"+{b - a}", (xi + w / 2, b), textcoords="offset points",
                         xytext=(0, 2), ha="center", fontsize=8.5,
                         color="#cc3311", fontweight="bold")
    axA.set_xticks(x)
    axA.set_xticklabels([TOOL_LABEL.get(m, m) for m in methods], rotation=25,
                        ha="right", rotation_mode="anchor", fontsize=8)
    axA.set_ylabel("Complexes counted as a success")
    axA.set_title("How many complexes qualify under each oracle selection  "
                  "(red = extra complexes rescued)", fontsize=11)
    axA.legend(fontsize=8, frameon=False); axA.grid(axis="y", alpha=0.3)
    figA.suptitle(suptitle, fontsize=12, fontweight="bold")
    figA.tight_layout(rect=(0, 0, 1, 0.96))
    figA.savefig(out_counts, dpi=160, bbox_inches="tight"); plt.close(figA)

    # ── Graph 2 — median form error of the reps under each selection ────────
    figB, axB = plt.subplots(figsize=(figw, 6))
    axB.bar(x - w / 2, f_near, w, color=c_near, label="nearest rule")
    axB.bar(x + w / 2, f_val, w, color=c_val, label="validity-constrained")
    axB.axhline(form_ok, ls="--", color="crimson", lw=1)
    axB.text(0.99, form_ok, f" correct-form line ({form_ok:g} Å)", color="crimson",
             va="bottom", ha="right", fontsize=8, transform=axB.get_yaxis_transform())
    axB.set_xticks(x)
    axB.set_xticklabels([TOOL_LABEL.get(m, m) for m in methods], rotation=25,
                        ha="right", rotation_mode="anchor", fontsize=8)
    axB.set_ylabel("Median form error — best-fit (Kabsch) RMSD (Å)")
    axB.set_title("Median form error of the reps under each selection  "
                  "(does rescuing poses change shape quality?)", fontsize=11)
    axB.legend(fontsize=8, frameon=False); axB.grid(axis="y", alpha=0.3)
    figB.suptitle(suptitle, fontsize=12, fontweight="bold")
    figB.tight_layout(rect=(0, 0, 1, 0.96))
    figB.savefig(out_form, dpi=160, bbox_inches="tight"); plt.close(figB)


# ───────────────────────────────────────────────────────────────────
# Didactic 3-D illustration of the mechanism labels (+ PyMOL PDB)
# ───────────────────────────────────────────────────────────────────
# A synthetic 8-atom model ligand, its crystal reference and three hand-built
# docked poses that land cleanly in each region, purely to make the placement-vs-
# form distinction concrete (NOT from the docking data):
#   placement-limited — identical conformation, rigidly translated + rotated
#   form-limited       — a genuinely different conformation aligned onto the crystal
#                        (position ~correct, only the shape is wrong)
#   mixed              — a milder shape change plus a comparable positional offset
_MECH_REF_XYZ = np.array([
    [-1.70, 0.10, -1.00], [-0.60, -0.60, -0.10], [0.50, 0.40, 0.40], [1.60, -0.20, 1.00],
    [1.20, 1.30, -0.40], [0.10, 2.10, 0.60], [-1.20, 1.60, -0.20], [-2.00, 2.70, 0.80],
], float)
_MECH_BONDS = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
_MECH_ATOM_NAMES = ["N1", "C2", "C3", "C4", "C5", "C6", "C7", "O8"]
_MECH_ATOM_ELEMS = ["N", "C", "C", "C", "C", "C", "C", "O"]
# Order the illustration reads left→right / how the PDB numbers residues.
_MECH_ORDER = ["placement-limited", "mixed", "form-limited"]


def _rot_about(axis, theta: float) -> np.ndarray:
    """Rotation matrix (right-handed) of angle ``theta`` rad about ``axis``."""
    a = np.asarray(axis, float); a = a / np.linalg.norm(a)
    w = math.cos(theta / 2); x, y, z = -a * math.sin(theta / 2)
    return np.array([[w*w+x*x-y*y-z*z, 2*(x*y+w*z), 2*(x*z-w*y)],
                     [2*(x*y-w*z), w*w+y*y-x*x-z*z, 2*(y*z+w*x)],
                     [2*(x*z+w*y), 2*(y*z-w*x), w*w+z*z-x*x-y*y]])


def _kabsch_onto(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rigidly best-fit ``P`` onto ``Q`` (P aligned to Q's centroid + orientation)."""
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Pc @ (Vt.T @ np.diag([1.0, 1.0, d]) @ U.T).T + Q.mean(0)


def _pt_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    return float(np.sqrt(((P - Q) ** 2).sum() / len(P)))


def _build_mechanism_examples():
    """Return ``(ref, examples)`` where ``examples`` maps each region to a dict with
    the pose ``xyz`` and its ``inplace`` / ``form`` / ``placement`` RMSD and ``r``
    (= form²/in-place²). Geometry is hand-tuned so each lands in its own region."""
    ref = _MECH_REF_XYZ

    def _tors(c, idx, ang_deg, axis, piv):
        c = c.copy(); R = _rot_about(axis, math.radians(ang_deg))
        for i in idx:
            c[i] = piv + R @ (c[i] - piv)
        return c

    def _deform(head_deg, tail_deg):
        # bend the two ends in OPPOSITE senses so no rigid motion can undo it
        c = _tors(ref, [0, 1], head_deg, ref[2] - ref[1], ref[1])
        return _tors(c, [6, 7], tail_deg, ref[5] - ref[4], ref[5])

    poses = {
        "placement-limited": ref @ _rot_about([0.3, 0.4, 1.0], math.radians(13)).T
        + np.array([0.50, 0.35, 0.30]),
        "form-limited": _kabsch_onto(_deform(85, -92), ref) + np.array([0.10, 0.0, 0.0]),
    }
    base = _kabsch_onto(_deform(50, -56), ref)
    f = _pt_rmsd(base, ref)
    poses["mixed"] = base + np.array([f * 0.95, 0.0, 0.0])

    ex = {}
    for region, xyz in poses.items():
        ipv = _pt_rmsd(xyz, ref)
        fmv = _pt_rmsd(_kabsch_onto(xyz, ref), ref)
        ex[region] = {"xyz": xyz, "inplace": ipv, "form": fmv,
                      "placement": math.sqrt(max(0.0, ipv ** 2 - fmv ** 2)),
                      "r": (fmv / ipv) ** 2 if ipv else 0.0}
    return ref, ex


def plot_mechanism_examples_3d(out: Path, pdb_out=None) -> None:
    """Didactic 3-D illustration of the three mechanism labels on a model ligand.

    Top row: each docked pose vs the crystal AS DOCKED (in-place RMSD mixes position
    and shape). Bottom row: the SAME pair after best-fit (Kabsch) superposition,
    which removes position and leaves only the form error. Placement-limited
    collapses onto the crystal (its error was all position); form-limited stays
    split (its error is the conformation); mixed collapses partway. Optionally also
    writes a PyMOL-loadable PDB (+ .pml) via ``pdb_out``."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)
    from matplotlib.lines import Line2D
    ref, ex = _build_mechanism_examples()
    title_col = {"placement-limited": "#2c7fb8", "mixed": "#666666", "form-limited": "#8856a7"}

    def _draw(ax, pose, region, superimposed):
        P = _kabsch_onto(pose, ref) if superimposed else pose
        col = _MECH_3D_COLORS[region]
        for i, j in _MECH_BONDS:
            ax.plot(*zip(ref[i], ref[j]), color="#333333", lw=2.0, alpha=0.9)
            ax.plot(*zip(P[i], P[j]), color=col, lw=2.6, alpha=0.95)
        ax.scatter(*ref.T, color="#333333", s=55, depthshade=True)
        ax.scatter(*P.T, color=col, s=62, depthshade=True)
        ax.set_axis_off()
        allpts = np.vstack([ref, P]); c = allpts.mean(0)
        rng = (allpts.max(0) - allpts.min(0)).max() / 2 + 0.3
        ax.set_xlim(c[0]-rng, c[0]+rng); ax.set_ylim(c[1]-rng, c[1]+rng); ax.set_zlim(c[2]-rng, c[2]+rng)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        ax.view_init(elev=18, azim=-72)

    fig = plt.figure(figsize=(15, 9.5))
    for k, region in enumerate(_MECH_ORDER):
        e = ex[region]
        axT = fig.add_subplot(2, 3, k + 1, projection="3d")
        _draw(axT, e["xyz"], region, superimposed=False)
        axT.set_title(f"{region.upper()}\nas docked  ·  in-place RMSD = {e['inplace']:.2f} Å",
                      fontsize=11, color=title_col[region], fontweight="bold")
        axB = fig.add_subplot(2, 3, k + 4, projection="3d")
        _draw(axB, e["xyz"], region, superimposed=True)
        axB.set_title("after best-fit superposition\n"
                      f"form (Kabsch) RMSD = {e['form']:.2f} Å  ·  r = {e['r']:.2f}",
                      fontsize=10)
        if e["form"] < 0.05:
            axB.text2D(0.5, 0.04, "→ docked pose coincides with the crystal\n"
                       "(all error was position)", transform=axB.transAxes, ha="center",
                       va="bottom", fontsize=8.5, color=title_col[region], style="italic")

    handles = [Line2D([0], [0], color="#333333", lw=3, label="crystal reference")]
    handles += [Line2D([0], [0], color=_MECH_3D_COLORS[r], lw=3, label=f"docked pose ({r})")
                for r in _MECH_ORDER]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(_vt("What the mechanism labels mean — one crystal ligand vs three "
                     "docked poses\ntop: as docked (position + shape mixed)  ·  bottom: "
                     "after removing position (Kabsch), only shape error remains"),
                 fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)

    if pdb_out is not None:
        _write_mechanism_pdb(ref, ex, Path(pdb_out))


def _write_mechanism_pdb(ref: np.ndarray, ex: dict, pdb_out: Path) -> None:
    """Write the illustration to a PyMOL-loadable PDB: three example groups laid out
    side by side (offset in x), each with the crystal reference (chain A) and the
    docked pose (chain B); resi 1/2/3 = placement/mixed/form. CONECT bonds included.
    Also writes a sibling ``.pml`` with display settings for a clean screenshot."""
    lines, conect = [], []
    serial = 0
    span = float(ref[:, 0].max() - ref[:, 0].min()) + 12.0

    def _emit(coords, chain, resi, dx):
        nonlocal serial
        base = serial
        for i, (nm, el) in enumerate(zip(_MECH_ATOM_NAMES, _MECH_ATOM_ELEMS)):
            serial += 1
            x, y, z = coords[i] + np.array([dx, 0.0, 0.0])
            lines.append(f"ATOM  {serial:>5d} {nm:<4s} LIG {chain}{resi:>4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
        for a, b in _MECH_BONDS:
            conect.append(f"CONECT{base + a + 1:>5d}{base + b + 1:>5d}")

    for k, region in enumerate(_MECH_ORDER):
        dx = k * span
        _emit(ref, "A", k + 1, dx)                 # crystal reference
        _emit(ex[region]["xyz"], "B", k + 1, dx)   # docked pose
    header = [
        "REMARK   Illustrative mechanism examples — synthetic model ligand (not docking data)",
        "REMARK   chain A = crystal reference, chain B = docked pose",
        "REMARK   resi 1 = placement-limited, resi 2 = mixed, resi 3 = form-limited",
    ]
    pdb_out.write_text("\n".join(header + lines + conect) + "\nEND\n")

    pml = pdb_out.with_suffix(".pml")
    pml.write_text(
        f"# Screenshot helper for {pdb_out.name}\n"
        "# crystal (chain A) and pose (chain B) OVERLAP, so bond ONLY via CONECT\n"
        "# records (connect_mode 1) — otherwise PyMOL draws spurious A–B bonds.\n"
        "set connect_mode, 1\n"
        f"load {pdb_out.name}, mech\n"
        "hide everything\nshow sticks, mech\nshow spheres, mech\n"
        "set sphere_scale, 0.22\nset stick_radius, 0.13\n"
        "color grey40, chain A\n"
        "color marine,  (chain B and resi 1)\n"   # placement-limited
        "color orange,  (chain B and resi 2)\n"   # mixed
        "color purple,  (chain B and resi 3)\n"   # form-limited
        "bg_color white\nset ray_opaque_background, 0\nset orthoscopic, 1\n"
        "set label_size, 18\n"
        "label chain A and resi 1 and name C4, 'placement-limited'\n"
        "label chain A and resi 2 and name C4, 'mixed'\n"
        "label chain A and resi 3 and name C4, 'form-limited'\n"
        "orient\nray 1600, 900\n"
    )


# ───────────────────────────────────────────────────────────────────
# REAL mechanism examples — actual docked poses from the benchmark
# ───────────────────────────────────────────────────────────────────
# Unlike the synthetic illustration above, these are genuine (method, complex,
# pose) triples pulled from per_pose_metrics.csv that cleanly exemplify each
# region, rendered against their real crystal ligand + written to PyMOL PDBs.
_MECH_REAL_BANDS = {
    "placement-limited": (0.00, 0.12),
    "mixed": (0.44, 0.56),
    "form-limited": (0.85, 1.00),
}


def _select_real_mechanism_examples(metrics: pd.DataFrame, benchmark_dir, root,
                                    ip_range=(1.3, 3.2)) -> dict:
    """Pick one real exemplar pose per region from ``metrics``. Prefers PB-valid
    poses whose r = form²/in-place² sits near the band centre and whose in-place
    RMSD is large enough (``ip_range``) that the deviation is visible, and whose
    pose + crystal SDFs both exist on disk. Returns {region: row-dict}."""
    d = metrics.dropna(subset=["rmsd", "bestfit_rmsd"]).copy()
    d["inplace"] = d["pb_rmsd"].where(d["pb_rmsd"].notna(), d["rmsd"])
    d["form"] = pd.to_numeric(d["bestfit_rmsd"], errors="coerce")
    d = d[(d["form"] <= d["inplace"] + 1e-6) & d["inplace"].between(*ip_range)]
    d["r"] = (d["form"] ** 2 / d["inplace"] ** 2).clip(0, 1)
    bdir = Path(benchmark_dir)

    def _files_ok(row) -> bool:
        pf = Path(str(row["pose_file"]))
        if not pf.is_absolute():
            pf = Path(root) / pf
        cry, _idx = _reference_copy_source(bdir, row)
        return pf.exists() and cry.exists()

    picks = {}
    for region, (lo, hi) in _MECH_REAL_BANDS.items():
        c = d[(d["r"] >= lo) & (d["r"] < hi)]
        if c.empty:
            continue
        c = c[c.apply(_files_ok, axis=1)]
        if c.empty:
            continue
        mid = (lo + hi) / 2
        c = c.assign(_valid=c["pb_valid"].astype(bool), _rc=(c["r"] - mid).abs())
        # PB-valid first, then r nearest the band centre, then the larger offset
        c = c.sort_values(["_valid", "_rc", "inplace"], ascending=[False, True, False])
        picks[region] = c.iloc[0].to_dict()
    return picks


def _reference_copy_source(bdir: Path, row) -> tuple[Path, int]:
    """(sdf path, record index) of the reference ligand copy a per-pose row was
    scored against. Rows carrying ``nearest_copy_index`` that differs from
    ``ref_copy_index`` (nearest-copy convention, alternate copy chosen) point at
    that record of <ID>_ligands.sdf; every other row — including every row of an
    ``instance`` table — keeps <ID>_ligand.sdf, record 0, exactly as before."""
    prot = str(row["protein"])
    single = Path(bdir) / prot / f"{prot}_ligand.sdf"
    try:
        j = row["nearest_copy_index"] if "nearest_copy_index" in row else None
        j = int(j) if j is not None and not pd.isna(j) else None
        ref_j = row["ref_copy_index"] if "ref_copy_index" in row else 0
        ref_j = int(ref_j) if ref_j is not None and not pd.isna(ref_j) else 0
    except (TypeError, ValueError, KeyError):
        j, ref_j = None, 0
    if j is None or j == ref_j:
        return single, 0
    if int(j) < 0:
        return single, 0
    copies = Path(bdir) / prot / f"{prot}_ligands.sdf"
    if copies.exists():
        return copies, j
    return single, 0


def _load_reference_copy(sdf_path: Path, index: int) -> Chem.Mol | None:
    """Record ``index`` of ``sdf_path`` (all deposited copies) or, for index 0, the
    first record — the same loader the scoring used."""
    if int(index) <= 0:
        return load_first_mol(sdf_path)
    mols = load_all_mols(sdf_path)
    if int(index) < len(mols):
        return mols[int(index)]
    return None


def _example_geometry(exemplar: dict, benchmark_dir, root) -> dict | None:
    """Load the reference copy of the crystal ligand (the one the exemplar row was
    scored against, see :func:`_reference_copy_source`) + docked heavy-atom mols for
    one exemplar and best-fit the pose onto it. Returns dict of coords / bonds /
    elements (as docked and superimposed), or None on any load/align failure."""
    from rdkit.Chem import rdMolAlign
    prot = exemplar["protein"]
    pf = Path(str(exemplar["pose_file"]))
    if not pf.is_absolute():
        pf = Path(root) / pf
    cpath, cidx = _reference_copy_source(Path(benchmark_dir), exemplar)
    pose, crystal = load_first_mol(pf), _load_reference_copy(cpath, cidx)
    if pose is None or crystal is None:
        return None
    try:
        ph, ch = Chem.RemoveHs(pose), Chem.RemoveHs(crystal)
        c_xyz, p_xyz = ch.GetConformer().GetPositions(), ph.GetConformer().GetPositions()
        rms, xform, _ = rdMolAlign.GetBestAlignmentTransform(ph, ch)
        X = np.asarray(xform); R, t = X[:3, :3], X[:3, 3]
        p_sup = p_xyz @ R.T + t
    except Exception:
        return None
    return {
        "protein": prot, "method": exemplar["method"], "rank": int(exemplar["rank"]),
        "inplace": float(exemplar["inplace"]), "form": float(exemplar["form"]),
        "r": float(exemplar["r"]), "valid": bool(exemplar["pb_valid"]),
        "c_xyz": c_xyz, "p_xyz": p_xyz, "p_sup": p_sup,
        "c_bonds": [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in ch.GetBonds()],
        "p_bonds": [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in ph.GetBonds()],
        "c_el": [a.GetSymbol() for a in ch.GetAtoms()],
        "p_el": [a.GetSymbol() for a in ph.GetAtoms()],
    }


def _emit_mol_pdb(lines, conect, serial0, coords, elems, bonds, chain, resi,
                  resname="LIG") -> int:
    """Append ATOM + CONECT records for one molecule; return the new serial base."""
    for i, (xyz, el) in enumerate(zip(coords, elems)):
        s = serial0 + i + 1
        nm = f"{el}{i + 1}"[:4]
        x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
        lines.append(f"ATOM  {s:>5d} {nm:<4s} {resname:>3s} {chain}{resi:>4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
    for a, b in bonds:
        conect.append(f"CONECT{serial0 + a + 1:>5d}{serial0 + b + 1:>5d}")
    return serial0 + len(coords)


def _write_real_example_pdbs(geoms: dict, out_dir: Path) -> list:
    """Write one PDB per region — chain A = reference ligand copy ("crystal
    (original)" under the instance convention, "nearest deposited copy" under
    nearest), chain B = docked pose (as docked, real receptor frame), chain C =
    docked pose after best-fit superposition onto that reference (shape-only) —
    plus a master .pml with two toggleable scenes (as-docked / aligned) and a
    manifest CSV. Returns paths."""
    written = []
    for region, g in geoms.items():
        lines, conect, s = [], [], 0
        s = _emit_mol_pdb(lines, conect, s, g["c_xyz"], g["c_el"], g["c_bonds"], "A", 1)
        s = _emit_mol_pdb(lines, conect, s, g["p_xyz"], g["p_el"], g["p_bonds"], "B", 1)
        s = _emit_mol_pdb(lines, conect, s, g["p_sup"], g["p_el"], g["p_bonds"], "C", 1)
        header = [
            f"REMARK   {region.upper()} — real docked example",
            f"REMARK   {g['method']} on {g['protein']} (rank {g['rank']}, "
            f"PB-valid={g['valid']})",
            f"REMARK   in-place RMSD={g['inplace']:.2f} A  form(Kabsch)={g['form']:.2f} A  "
            f"r={g['r']:.2f}",
            f"REMARK   chain A = {_ref_copy_label()} ligand",
            "REMARK   chain B = docked pose, as docked (position + shape error)",
            "REMARK   chain C = docked pose after best-fit superposition (shape error only)",
        ]
        fn = out_dir / f"mechanism_real_{region.replace('-', '_')}.pdb"
        fn.write_text("\n".join(header + lines + conect) + "\nEND\n")
        written.append(fn)

    # master .pml — load each as its own object (grid), with two scenes to toggle
    objs = [(r, f"mechanism_real_{r.replace('-', '_')}") for r in geoms]
    col = {"placement-limited": "marine", "mixed": "orange", "form-limited": "purple"}
    pml = ["# Real docked mechanism examples",
           f"# chain A = {_ref_copy_label()} · chain B = docked as-docked · "
           "chain C = docked best-fit aligned",
           "# crystal and pose OVERLAP → bond only via CONECT (connect_mode 1).",
           "set connect_mode, 1"]
    for _region, obj in objs:
        pml.append(f"load {obj}.pdb, {obj}")
    pml += ["set grid_mode, 1", "set sphere_scale, 0.20", "set stick_radius, 0.14",
            "color grey50, chain A"]
    for region, obj in objs:
        pml.append(f"color {col.get(region, 'orange')}, ({obj} and (chain B or chain C))")
    pml += [
        "bg_color white", "set ray_opaque_background, 0",
        "# ---- toggle the two views with the scene buttons (or type the scene name) ----",
        "# scene as_docked : crystal (A) vs docked AS DOCKED (B) -> position + shape error",
        "# scene aligned   : crystal (A) vs docked BEST-FIT ALIGNED (C) -> ONLY the shape error",
        "hide everything",
        "show sticks, chain A or chain B", "show spheres, chain A or chain B",
        "orient", "scene as_docked, store",
        "hide everything",
        "show sticks, chain A or chain C", "show spheres, chain A or chain C",
        "scene aligned, store",
        "scene as_docked", "ray 1800, 700",
    ]
    pml_path = out_dir / "mechanism_real_examples.pml"
    pml_path.write_text("\n".join(pml) + "\n")
    written.append(pml_path)

    manifest = pd.DataFrame([{
        "region": r, "method": g["method"], "protein": g["protein"], "rank": g["rank"],
        "pb_valid": g["valid"], "inplace_rmsd_A": round(g["inplace"], 3),
        "form_kabsch_rmsd_A": round(g["form"], 3), "r_form_share": round(g["r"], 3),
    } for r, g in geoms.items()])
    man_path = out_dir / "mechanism_real_examples.csv"
    _write_csv(manifest, man_path, index=False)
    written.append(man_path)
    return written


def _emit_het_pdb(lines, conect, serial0, coords, elems, bonds, chain, resi,
                  resname) -> int:
    """Append HETATM + CONECT records for one ligand; return the new serial base."""
    for i, (xyz, el) in enumerate(zip(coords, elems)):
        s = serial0 + i + 1
        nm = f"{el}{i + 1}"[:4]
        x, y, z = float(xyz[0]), float(xyz[1]), float(xyz[2])
        lines.append(f"HETATM{s:>5d} {nm:<4s} {resname:>3s} {chain}{resi:>4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2s}")
    for a, b in bonds:
        conect.append(f"CONECT{serial0 + a + 1:>5d}{serial0 + b + 1:>5d}")
    return serial0 + len(coords)


def _write_receptor_context_pdbs(geoms: dict, benchmark_dir, out_dir: Path) -> list:
    """Per region, write receptor + reference ligand copy ("crystal (original)" or,
    under the nearest convention, the "nearest deposited copy") + docked pose in one
    PDB: the protein records verbatim, then the two ligands as HETATM (resn CRY /
    DOK) with serials continuing past the protein and chain IDs that don't collide
    with the protein's. Plus a master .pml (cartoon + sticks). Returns file paths."""
    written = []
    bdir = Path(benchmark_dir)
    for region, g in geoms.items():
        prot = g["protein"]
        ppath = bdir / str(prot) / f"{prot}_protein.pdb"
        if not ppath.exists():
            print(f"  [receptor-context] {ppath} missing; skipping {region}.")
            continue
        struct, rec_conect, max_serial, used_chains = [], [], 0, set()
        for ln in ppath.read_text().splitlines():
            rec = ln[:6].strip()
            if rec in ("END", "MASTER", "ENDMDL"):
                continue
            if rec == "CONECT":
                rec_conect.append(ln); continue
            if rec in ("ATOM", "HETATM"):
                try:
                    max_serial = max(max_serial, int(ln[6:11]))
                except ValueError:
                    pass
                ch = ln[21:22].strip()
                if ch:
                    used_chains.add(ch)
            struct.append(ln)
        free = [c for c in "XYZWVUTSR" if c not in used_chains]
        cry_ch = free[0] if free else "X"
        dok_ch = free[1] if len(free) > 1 else "Y"
        lig, lig_conect, s = [], [], max_serial
        s = _emit_het_pdb(lig, lig_conect, s, g["c_xyz"], g["c_el"], g["c_bonds"],
                          cry_ch, 900, "CRY")
        s = _emit_het_pdb(lig, lig_conect, s, g["p_xyz"], g["p_el"], g["p_bonds"],
                          dok_ch, 901, "DOK")
        header = [
            f"REMARK   {region.upper()} in receptor context — {g['method']} on {prot} "
            f"(rank {g['rank']}, PB-valid={g['valid']})",
            f"REMARK   in-place RMSD={g['inplace']:.2f} A  form(Kabsch)={g['form']:.2f} A  "
            f"r={g['r']:.2f}",
            f"REMARK   protein = original PDB chains {''.join(sorted(used_chains))} | "
            f"resn CRY = {_ref_copy_label()} ligand (chain {cry_ch}) | "
            f"resn DOK = docked pose (chain {dok_ch})",
        ]
        body = header + struct + lig + rec_conect + lig_conect + ["END"]
        fn = out_dir / f"mechanism_real_{region.replace('-', '_')}_with_receptor.pdb"
        fn.write_text("\n".join(body) + "\n")
        written.append(fn)

    objs = [(r, f"mechanism_real_{r.replace('-', '_')}_with_receptor") for r in geoms
            if (out_dir / f"mechanism_real_{r.replace('-', '_')}_with_receptor.pdb").exists()]
    if not objs:
        return written
    col = {"placement-limited": "marine", "mixed": "orange", "form-limited": "purple"}
    pml = ["# Receptor + crystal (original, resn CRY) + docked pose (resn DOK) per example"
           if REFERENCE_CONVENTION == "instance" else
           f"# Receptor + {_ref_copy_label()} (resn CRY) + docked pose (resn DOK) per example",
           "# grid_mode tiles the three; ligands overlap so drop spurious ligand bonds."]
    for _r, obj in objs:
        pml.append(f"load {obj}.pdb, {obj}")
    pml += ["set grid_mode, 1", "hide everything",
            "show cartoon, polymer", "color grey80, polymer",
            "set cartoon_transparency, 0.45",
            "show sticks, resn CRY+DOK", "show spheres, resn CRY+DOK",
            "set sphere_scale, 0.22, resn CRY+DOK",
            "unbond resn CRY, resn DOK", "unbond polymer, (resn CRY or resn DOK)",
            "color green, resn CRY"]
    for region, obj in objs:
        pml.append(f"color {col.get(region, 'orange')}, ({obj} and resn DOK)")
    pml += ["bg_color white", "set ray_opaque_background, 0",
            f"# green = {_ref_copy_label()}, coloured = docked pose; grey = receptor",
            "# to focus one example: set grid_mode, 0 ; disable all ; enable <object>",
            "orient resn CRY+DOK", "zoom resn CRY+DOK, 6", "ray 1800, 700"]
    pml_path = out_dir / "mechanism_real_with_receptor.pml"
    pml_path.write_text("\n".join(pml) + "\n")
    written.append(pml_path)
    return written


def plot_mechanism_examples_real(metrics: pd.DataFrame, benchmark_dir, out: Path,
                                 root=None) -> None:
    """Real docked poses that exemplify each mechanism label. Selects one clean
    exemplar per region from ``metrics``, renders a 3-D figure (top = as docked,
    bottom = after best-fit superposition) against the real crystal ligand, and
    writes per-region PyMOL PDBs (crystal + docked pose) + a manifest via
    :func:`_write_real_example_pdbs`. No-op if no exemplars/files are available."""
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    from matplotlib.lines import Line2D
    root = str(root or Path.cwd())
    picks = _select_real_mechanism_examples(metrics, benchmark_dir, root)
    geoms = {}
    for region in _MECH_ORDER:
        if region in picks:
            g = _example_geometry(picks[region], benchmark_dir, root)
            if g is not None:
                geoms[region] = g
    if not geoms:
        print("  [mechanism-real] no usable exemplars found; skipping 20g.")
        return
    order = [r for r in _MECH_ORDER if r in geoms]
    title_col = {"placement-limited": "#2c7fb8", "mixed": "#666666", "form-limited": "#8856a7"}

    def _draw(ax, g, region, superimposed):
        col = _MECH_3D_COLORS[region]
        P = g["p_sup"] if superimposed else g["p_xyz"]
        C = g["c_xyz"]
        for i, j in g["c_bonds"]:
            ax.plot(*zip(C[i], C[j]), color="#333333", lw=1.7, alpha=0.85)
        for i, j in g["p_bonds"]:
            ax.plot(*zip(P[i], P[j]), color=col, lw=2.2, alpha=0.9)
        ax.scatter(*C.T, color="#333333", s=14, depthshade=True)
        ax.scatter(*P.T, color=col, s=16, depthshade=True)
        ax.set_axis_off()
        allp = np.vstack([C, P]); ctr = allp.mean(0)
        rng = (allp.max(0) - allp.min(0)).max() / 2 + 0.8
        ax.set_xlim(ctr[0]-rng, ctr[0]+rng); ax.set_ylim(ctr[1]-rng, ctr[1]+rng); ax.set_zlim(ctr[2]-rng, ctr[2]+rng)
        try:
            ax.set_box_aspect((1, 1, 1))
        except Exception:
            pass
        ax.view_init(elev=16, azim=-70)

    n = len(order)
    fig = plt.figure(figsize=(5.2 * n, 9.5))
    for k, region in enumerate(order):
        g = geoms[region]
        axT = fig.add_subplot(2, n, k + 1, projection="3d")
        _draw(axT, g, region, superimposed=False)
        axT.set_title(f"{region.upper()}\n{TOOL_LABEL.get(g['method'], g['method'])} · "
                      f"{g['protein']} (rank {g['rank']})\nas docked · in-place RMSD "
                      f"= {g['inplace']:.2f} Å", fontsize=10, color=title_col[region],
                      fontweight="bold")
        axB = fig.add_subplot(2, n, k + 1 + n, projection="3d")
        _draw(axB, g, region, superimposed=True)
        axB.set_title("after best-fit superposition\n"
                      f"form (Kabsch) RMSD = {g['form']:.2f} Å  ·  r = {g['r']:.2f}",
                      fontsize=9.5)

    handles = [Line2D([0], [0], color="#333333", lw=3, label=f"{_ref_copy_label()} ligand")]
    handles += [Line2D([0], [0], color=_MECH_3D_COLORS[r], lw=3, label=f"docked pose ({r})")
                for r in order]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, -0.01))
    fig.suptitle(_vt("Real docked poses illustrating each mechanism label\n"
                     "top: as docked (position + shape)  ·  bottom: after best-fit "
                     "superposition (only shape error remains)"),
                 fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)

    _write_real_example_pdbs(geoms, out.parent)
    _write_receptor_context_pdbs(geoms, benchmark_dir, out.parent)
    print("  wrote real mechanism examples → 20g_mechanism_examples_real_3d.png + "
          f"{len(geoms)} ligand PDB(s) + {len(geoms)} receptor-context PDB(s) + manifest "
          f"[{', '.join(f'{r}:{geoms[r]['method']}/{geoms[r]['protein']}' for r in order)}]")


def aggregate_oracle_rank_distribution(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    """For each ranking tool: % of complexes whose oracle (best-RMSD) pose sits
    at exactly rank k, plus the running cumulative %. Ranks beyond ``top_n`` are
    pooled into a ``> top_n`` overflow bin so the columns always sum to 100 %.
    """
    oracle = _oracle_per_pair(df)
    oracle = oracle[oracle["method"].isin(RANKING_TOOLS)].copy()
    if oracle.empty:
        return pd.DataFrame()
    oracle["rank_bin"] = oracle["rank"].clip(upper=top_n + 1)
    rows = []
    for method, g in oracle.groupby("method"):
        n = len(g)
        counts = g["rank_bin"].value_counts()
        cum = 0.0
        for k in range(1, top_n + 2):
            pct = 100.0 * int(counts.get(k, 0)) / n if n else float("nan")
            cum += pct
            rows.append({
                "method": method,
                "rank": (str(k) if k <= top_n else f">{top_n}"),
                "n_pairs": n,
                "pct_at_rank": round(pct, PCT_DECIMALS),
                "cum_pct": round(cum, PCT_DECIMALS),
            })
    return pd.DataFrame(rows)


def plot_oracle_rank_distribution(dist: pd.DataFrame, top_n: int, out: Path) -> None:
    """How often the n-th ranked pose is the closest (oracle) one.

    Left: P(oracle pose is at exactly rank k) as grouped bars per tool.
    Right: cumulative — % of complexes whose oracle pose is recovered within
    the top-k (how much rescoring deeper buys you).
    """
    if dist is None or dist.empty:
        return
    methods = list(dict.fromkeys(dist["method"]))

    def _kpos(r) -> int:
        return top_n + 1 if str(r).startswith(">") else int(r)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    w = 0.8 / max(1, len(methods))
    for mi, method in enumerate(methods):
        g = dist[dist["method"] == method].copy()
        g["k"] = g["rank"].map(_kpos)
        g = g.sort_values("k")
        n = int(g["n_pairs"].iloc[0])
        offset = (mi - (len(methods) - 1) / 2) * w
        axes[0].bar(g["k"] + offset, g["pct_at_rank"], w,
                    label=f"{TOOL_LABEL.get(method, method)} (n={n})",
                    color=TOOL_COLORS.get(method), edgecolor="black", alpha=0.85)
        axes[1].plot(g["k"], g["cum_pct"], marker="o", lw=2,
                     label=f"{TOOL_LABEL.get(method, method)} (n={n})",
                     color=TOOL_COLORS.get(method))

    xt = list(range(1, top_n + 2))
    xl = [str(k) for k in range(1, top_n + 1)] + [f">{top_n}"]
    for ax in axes:
        ax.set_xticks(xt); ax.set_xticklabels(xl); ax.grid(axis="y", alpha=0.3)
    axes[0].set_xlabel("Rank k of the tool's pose")
    axes[0].set_ylabel("% of complexes where the oracle pose is at rank k")
    axes[0].set_title("Where does the best pose rank?")
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("Top-k")
    axes[1].set_ylabel("% of complexes whose oracle pose is within top-k")
    axes[1].set_title("Cumulative recovery of the oracle pose")
    axes[1].set_ylim(0, 105); axes[1].legend(fontsize=8)
    _label_panels(axes)
    fig.suptitle(_vt("How often is the n-th ranked pose the closest (oracle) pose?"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160); plt.close(fig)


# ── 15b: cumulative recovery within top-k — near-native vs. PB-valid ─────────
# Extends fig 15's recovery panel: raw DiffDock alongside DiffDock* (smina), and a
# validity-aware (dashed) line beside each variant's RMSD-only (solid) line. Fig 15
# ignores PB-validity (it's purely rank-of-min-RMSD); this makes validity explicit.
_TOPK_RECOVERY_SPECS = [
    # AutoDock = the dominant arm (ADFRsuite ligands, exhaustiveness 128) and its own
    # gnina-rescored counterpart, matching 09b/09d/09f and the 21/22 ranking-quality
    # table. The method keys here are the BENCHMARK ones; _topk_recovery_specs()
    # resolves them against the frame so an Orai frame still reads its literal
    # "autodock" key.
    ("AutoDock Vina",        _DOMINANT_AUTODOCK_RAW,    "native"),
    ("AutoDock (gnina-opt)", _DOMINANT_AUTODOCK_OPT,    "native"),
    ("DiffDock (raw)",       "diffdock",                "native"),
    # ADDED 2026-09-02. smina is the DiffDock variant the Results chapter carries
    # forward, so the table that reports DiffDock's optimisation should show it
    # beside the gnina one rather than only the arm that is not selected.
    ("DiffDock + smina",     "diffdock_smina",          "native"),
    ("DiffDock (gnina-opt)", "diffdock_gnina",          "native"),
    ("EquiBind (raw)",       "equibind_unguided_raw",   "generation"),
    ("EquiBind (gnina-opt)", "equibind_unguided_gnina", "gnina"),
]

# The key the AutoDock RAW slot falls back to when the dominant arm is not in the
# frame. This is the Meeko-ligand run, which is the reported AutoDock pipeline on the
# Orai datasets (their panels cannot be re-docked) and is excluded on the Benchmark.
_TOPK_AUTODOCK_FALLBACK = "autodock"


def _topk_recovery_specs(present) -> list:
    """:data:`_TOPK_RECOVERY_SPECS` with the AutoDock method key resolved against the
    methods actually in the frame, so ONE module-level table serves both datasets.

    The Benchmark whole-protein report drops every Meeko-ligand arm (method_filter
    preset ``meeko``), so a table hard-coding the literal ``autodock`` key silently
    loses the whole AutoDock family — no rows in topk_recovery_validity.csv, and the
    three pre-specified AutoDock pairs quietly gone from the Holm family. The Orai
    reports are the mirror image: ``autodock`` IS their reported arm and the dominant
    exh128 keys are absent. Resolving on presence rather than hard-coding either key
    keeps both right without a dataset flag.

    Only the RAW slot falls back. There is deliberately no gnina-opt fallback: an Orai
    frame can carry ``autodock_gnina``, and adopting it here would push a curve into
    artefacts that never had one.
    """
    present = {str(m) for m in present}
    if _DOMINANT_AUTODOCK_RAW in present or _TOPK_AUTODOCK_FALLBACK not in present:
        return list(_TOPK_RECOVERY_SPECS)
    return [(variant, _TOPK_AUTODOCK_FALLBACK if mkey == _DOMINANT_AUTODOCK_RAW else mkey,
             kind) for variant, mkey, kind in _TOPK_RECOVERY_SPECS]


# Pre-specified depths for the fig-15b significance tests. The 30 cumulative curve
# points are nested (recovered by k ⟹ recovered by k+1) and so heavily auto-
# correlated; testing all of them would be pseudo-replicated multiplicity. We fix a
# SMALL set of depths in advance — top-1 (the deployment pick), a short shortlist,
# and two deeper reads — and Holm-correct across the resulting family.
_TOPK_TEST_DEPTHS = (1, 5, 10, 15)

# Pre-specified between-method comparisons tested at each depth in _TOPK_TEST_DEPTHS.
# Each is a (baseline, comparison) pair of variant labels from _TOPK_RECOVERY_SPECS;
# the full family of (pair × depth) exact-McNemar p-values is Holm-corrected together.
_TOPK_TEST_PAIRS = [
    # ADDED 2026-09-02, the two within-tool contrasts this family was missing.
    #
    # AutoDock was absent while DiffDock and EquiBind were present, so the table
    # showed a raw-to-optimised effect for the two deep-learning tools and only the
    # optimised arm for the physics one. The contrast is not the same operation for
    # all three and the table note now says so: AutoDock's gnina step minimises AND
    # re-orders on CNN affinity, DiffDock keeps its confidence order, and EquiBind
    # has none to keep. Reporting all three is a comparison of what each deployed
    # pipeline gains from the same tool, which is what the chapter is about.
    ("AutoDock Vina",        "AutoDock (gnina-opt)"),      # rescoring effect (AutoDock)
    ("DiffDock (raw)",       "DiffDock + smina"),          # refinement effect, smina
    ("DiffDock (raw)",       "DiffDock (gnina-opt)"),      # refinement effect (DiffDock)
    ("EquiBind (raw)",       "EquiBind (gnina-opt)"),      # refinement effect (EquiBind)
    # -- the four above are WITHIN-tool; the rest are BETWEEN-tool. --
    # CORRECTED 2026-09-03. The cross-tool baseline is the gnina-rescored AutoDock
    # arm, not raw Vina. Appendix C declares the cross-tool family is "computed on
    # the gnina-rescored AutoDock Vina arm that the recovery tables below display",
    # and Table 21 prints that arm's rates (37.3/62.0/65.0/66.3). This list carried
    # raw Vina (31.4/50.8/58.7/61.1) in every commit of its history, so the script
    # tested a contrast the thesis never reported and never reproduced Table 21.
    # Restoring the declared arm reproduces the printed note exactly.
    ("AutoDock (gnina-opt)", "DiffDock (raw)"),            # physics vs raw DL
    ("AutoDock (gnina-opt)", "DiffDock (gnina-opt)"),      # physics vs refined DL
    ("DiffDock (gnina-opt)", "EquiBind (gnina-opt)"),      # DL vs DL (matched refinement)
    ("AutoDock (gnina-opt)", "EquiBind (gnina-opt)"),      # physics vs EquiBind
]

# Two Holm families, not one, changed 2026-09-02.
#
# The within-tool contrasts ask "does optimising this tool's own poses help?" and
# the between-tool ones ask "which tool wins?". Those are different questions, and
# pooling them means adding a tool's optimisation row makes every cross-tool verdict
# more conservative for no methodological reason. Each family is now corrected on
# its own and the table caption states which family a p-value belongs to.
_TOPK_WITHIN_TOOL_PAIRS = {
    ("AutoDock Vina",  "AutoDock (gnina-opt)"),
    ("DiffDock (raw)", "DiffDock + smina"),
    ("DiffDock (raw)", "DiffDock (gnina-opt)"),
    ("EquiBind (raw)", "EquiBind (gnina-opt)"),
}


def _topk_recovery_ranks(df: pd.DataFrame, thr: float = 2.0):
    """Per variant, the per-complex min rank at which a near-native (RMSD ≤ ``thr``) and
    a validity-aware (near-native AND PB-valid) pose first appears under the variant's
    ranking. Shared by :func:`aggregate_topk_recovery` (curve %) and
    :func:`topk_recovery_stats` (paired tests) so both read from identical ranking logic.

    Returns an ordered list of ``(variant, method_key, frame)`` where ``frame`` is
    indexed by ``(protein, ligand)`` over ALL the variant's complexes with float columns
    ``near_rank`` / ``valid_rank`` (NaN when the variant never yields such a pose — NaN
    ≤ k is False, so those complexes count as un-recovered at every depth). Needs the
    FULL frame so raw ``diffdock`` and its optimised variants are all present.
    """
    present = set(df["method"].astype(str))
    out = []
    for variant, mkey, kind in _topk_recovery_specs(present):
        if mkey not in present:
            continue
        sub = df[df["method"].astype(str) == mkey].copy()
        rk = (sub["rank"] if kind == "native"
              else sub["pose_name"].map(_equibind_pose_index) if kind == "generation"
              else _gnina_affinity_rank(sub))
        sub["_rk"] = pd.to_numeric(rk, errors="coerce")
        sub = sub.dropna(subset=["_rk", "rmsd"])
        idx = sub.groupby(["protein", "ligand"]).size().index    # every complex present
        if not len(idx):
            continue
        near = sub["rmsd"] <= thr
        valid = near & _to_bool(sub["pb_valid"])
        frame = pd.DataFrame(index=idx)
        frame["near_rank"] = (sub[near].groupby(["protein", "ligand"])["_rk"].min()
                              .reindex(idx))
        frame["valid_rank"] = (sub[valid].groupby(["protein", "ligand"])["_rk"].min()
                               .reindex(idx))
        out.append((variant, mkey, frame))
    return out


def aggregate_topk_recovery(df: pd.DataFrame, ks, thr: float = 2.0) -> pd.DataFrame:
    """Per variant × depth k: cumulative recovery within the variant's top-k poses.

    ``near_recovery_%`` = % of complexes with ≥ 1 top-k pose RMSD ≤ ``thr``;
    ``valid_recovery_%`` = % with ≥ 1 top-k pose RMSD ≤ ``thr`` AND PB-valid (the
    validity-aware line). Both are monotone non-decreasing in k and valid ≤ near, so
    the gap between them is the accurate-but-invalid share at depth k. ``near_k`` /
    ``valid_k`` carry the integer counts behind the percentages so callers (e.g. the
    Wilson-CI bands in :func:`plot_topk_recovery_validity`) needn't re-derive them from
    rounded percentages. Needs the FULL frame so raw ``diffdock`` and ``diffdock_smina``
    are both present.
    """
    ks = sorted({int(k) for k in ks})
    rows = []
    for variant, mkey, frame in _topk_recovery_ranks(df, thr):
        N = len(frame)
        near_rank, valid_rank = frame["near_rank"], frame["valid_rank"]
        for k in ks:
            nk = int((near_rank <= k).sum())
            vk = int((valid_rank <= k).sum())
            rows.append({
                "variant": variant, "method_key": mkey, "k": k, "n_complexes": N,
                "near_k": nk, "valid_k": vk,
                "near_recovery_%": round(100 * nk / N, PCT_DECIMALS),
                "valid_recovery_%": round(100 * vk / N, PCT_DECIMALS),
            })
    return pd.DataFrame(rows)


def topk_recovery_stats(df: pd.DataFrame, thr: float = 2.0,
                        test_ks=_TOPK_TEST_DEPTHS) -> dict | None:
    """Paired significance tests for fig 15b at PRE-SPECIFIED depths (default
    k ∈ {1, 5, 10, 15}) — a fixed, small number of tests rather than one per cumulative
    point. Unit = the receptor-ligand complex (recovered = ≥ 1 top-k pose meets the
    criterion), so there is no pose-level pseudo-replication, and every comparison is
    PAIRED across methods on the same complexes.

    Two families, both on the per-complex outcome:
      * ``between_method`` — exact McNemar on the NEAR-NATIVE outcome for each pre-
        specified (baseline, comparison) pair × depth in :data:`_TOPK_TEST_PAIRS`. ALL
        these p-values form ONE Holm family (``p_holm`` / ``star`` back-filled).
      * ``near_vs_valid_gap`` — within a variant, the near-native vs validity-aware gap.
        Because the validity-aware set is a strict SUBSET of near-native, discordance is
        one-directional and a McNemar p-value is degenerate (≈ 0.5**b); we therefore
        report the gap as an EFFECT SIZE — the accurate-but-invalid share (b/N) with a
        Wilson 95 % CI — and keep it OUT of the p-value family.

    Returns the sidecar dict, or None when < 2 variants are present.
    """
    test_ks = sorted({int(k) for k in test_ks})
    ranks = {variant: frame for variant, _mk, frame in _topk_recovery_ranks(df, thr)}
    if len(ranks) < 2:
        return None

    def _cmap(frame: pd.DataFrame, col: str, k: int) -> dict:
        """complex-id → recovered-by-depth-k flag (paired keys for _mcnemar_on_maps)."""
        rec = (frame[col] <= k)
        return {f"{p}||{l}": bool(v) for (p, l), v in zip(frame.index, rec)}

    between = []
    fams = {"within_tool": [], "between_tool": []}
    for a, b in _TOPK_TEST_PAIRS:
        if a not in ranks or b not in ranks:
            continue
        for k in test_ks:
            rec = _mcnemar_on_maps(_cmap(ranks[a], "near_rank", k),
                                   _cmap(ranks[b], "near_rank", k))
            if rec is None:
                continue
            which = ("within_tool" if (a, b) in _TOPK_WITHIN_TOOL_PAIRS
                     else "between_tool")
            rec = {"baseline": a, "comparison": b, "k": k, "family": which, **rec}
            between.append(rec)
            fams[which].append(rec)
    # Holm within each family separately; `family` is written out so a reader can
    # see which set of hypotheses a given p-value was corrected against.
    for name, fam in fams.items():
        if not fam:
            continue
        for rec, pa in zip(fam, su.holm([r["mcnemar_p"] for r in fam])):
            rec["p_holm"] = float(pa)
            rec["star"] = su.p_stars(pa)
            rec["family_size"] = len(fam)

    gap = []
    for variant, frame in ranks.items():
        N = len(frame)
        for k in test_ks:
            near = (frame["near_rank"] <= k)
            valid = (frame["valid_rank"] <= k)
            nk, vk = int(near.sum()), int(valid.sum())
            b = int((near & ~valid).sum())          # accurate but INVALID at depth k
            lo, hi = su.wilson_ci(b, N)
            gap.append({
                "variant": variant, "k": k, "n": N, "near_k": nk, "valid_k": vk,
                "near_rate": nk / N, "valid_rate": vk / N,
                "gap_pp": round(100 * (nk - vk) / N, PCT_DECIMALS),
                "accurate_invalid_k": b, "accurate_invalid_share": b / N,
                "accurate_invalid_ci": [lo, hi], "nested": True,
            })

    return {"test_depths": test_ks,
            "unit": "per-complex (≥1 top-k pose meets criterion); paired across methods",
            "between_method": between, "near_vs_valid_gap": gap}


def plot_topk_recovery_validity(rec: pd.DataFrame, out: Path, thr: float = 2.0,
                                stats: dict | None = None) -> None:
    """Fig 15b — cumulative recovery within top-k: one solid (near-native, RMSD ≤ 2 Å)
    and one dashed (validity-aware: near-native AND PB-valid) line per variant, colour
    = variant. Raw DiffDock and DiffDock* (smina) appear side by side; the gap between
    a variant's two lines is its accurate-but-invalid share at that depth. Unlike fig
    15 (rank of the min-RMSD pose, RMSD only), this brings PB-validity into the view.
    Every curve is annotated with its recovery % (one decimal) where it crosses each
    x-axis tick — k = 1, 5, 10, 15, 20, 25, 30 (dotted reference lines): near-native
    values left of each line, validity-aware right of it.

    Each curve carries a shaded Wilson 95 % CI band (from the ``near_k`` / ``valid_k``
    counts) — a proportion CI, so it stays valid at the extremes (EquiBind ≈ 1 %,
    AutoDock plateau ≈ 89 %) where a Gaussian/Wald interval would run past 0/100 %. The
    four pre-specified test depths k ∈ {1, 5, 10, 15} (see :func:`topk_recovery_stats`)
    get a darker reference line; ``stats`` (when supplied) drives only the footnote —
    the full paired-McNemar table lives in the CSV / JSON sidecar."""
    if rec is None or rec.empty:
        return
    from matplotlib.lines import Line2D
    variants = [v for v in (s[0] for s in _TOPK_RECOVERY_SPECS)
                if v in set(rec["variant"])]
    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    for variant in variants:
        g = rec[rec["variant"] == variant].sort_values("k")
        c = _RANK1_TOPN_COLORS.get(variant, "#888888")
        # Wilson 95% CI band per line (needs the integer counts; skip gracefully if a
        # cached CSV predates the near_k/valid_k columns). Near-native band a touch
        # stronger than the validity-aware one so overlapping variants stay legible.
        try:
            N = int(g["n_complexes"].iloc[0])
            for col, alpha in (("near_k", 0.13), ("valid_k", 0.07)):
                if col not in g.columns:
                    continue
                cis = [su.wilson_ci(int(kk), N) for kk in g[col]]
                ax.fill_between(g["k"], [100 * lo for lo, _ in cis],
                                [100 * hi for _, hi in cis],
                                color=c, alpha=alpha, lw=0, zorder=2)
        except Exception:
            pass
        ax.plot(g["k"], g["near_recovery_%"], color=c, lw=2.2, ls="-", zorder=3)
        ax.plot(g["k"], g["valid_recovery_%"], color=c, lw=2.2, ls="--", zorder=3)
        kmax = int(g["k"].max())
        ax.annotate(variant, xy=(kmax, float(g.loc[g["k"] == kmax, "near_recovery_%"].iloc[0])),
                    xytext=(6, 0), textcoords="offset points", fontsize=9, color=c,
                    fontweight="bold", va="center")
    # Read-off at EVERY x-axis tick (k = 1, 5, 10, 15, 20, 25, 30): annotate each curve's
    # recovery % (one decimal) where it crosses the tick, each a dotted reference line with
    # markers at both crossings. Near-native (solid, bold) is parked LEFT of each line and
    # validity-aware (dashed) RIGHT of it, so the two never collide. Within a column a small
    # vertical de-collision nudges labels that would otherwise overlap (e.g. DiffDock raw vs
    # smina near-native at k = 1, ~1.6 % apart). At the rightmost tick both columns go LEFT
    # so they clear the variant-name labels sitting to the right.
    kmax = int(rec["k"].max())
    MARK_KS = [k for k in (1, 5, 10, 15, 20, 25, 30) if k <= kmax]
    FLOOR, MIN_GAP = 1.8, 3.0    # keep labels off the x-axis; min vertical gap (% points)

    def _place_column(kx, items, side):
        """items: list of (y_true, color, bold); draw de-collided labels on one side."""
        dx, ha = (-5, "right") if side == "left" else (5, "left")
        prev = None
        for y_true, color, bold in sorted(items, key=lambda t: t[0]):
            y = max(y_true, FLOOR)
            if prev is not None and y - prev < MIN_GAP:
                y = prev + MIN_GAP           # nudge up off the label below it
            prev = y
            ax.annotate(f"{y_true:.1f}%", xy=(kx, y), xytext=(dx, 0),
                        textcoords="offset points", ha=ha, va="center", fontsize=7.5,
                        color=color, fontweight="bold" if bold else "normal", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.8))

    for mk in MARK_KS:
        is_test = mk in _TOPK_TEST_DEPTHS       # a pre-specified McNemar test depth
        ax.axvline(mk, color="0.45" if is_test else "0.72", ls=(0, (2, 3)),
                   lw=1.2 if is_test else 0.9, zorder=1)
        if is_test:
            ax.annotate("test k", xy=(mk, 100), xytext=(0, 1.5),
                        textcoords="offset points", ha="center", va="bottom",
                        fontsize=6, color="0.45", clip_on=False)
        near_col, valid_col = [], []
        for variant in variants:
            row = rec[(rec["variant"] == variant) & (rec["k"].astype(int) == mk)]
            if row.empty:
                continue
            c = _RANK1_TOPN_COLORS.get(variant, "#888888")
            near_y = float(row["near_recovery_%"].iloc[0])
            valid_y = float(row["valid_recovery_%"].iloc[0])
            ax.plot([mk, mk], [near_y, valid_y], marker="o", ms=3.5, ls="none",
                    color=c, zorder=5, mec="white", mew=0.5)
            near_col.append((near_y, c, True))    # solid = near-native (RMSD ≤ 2 Å)
            valid_col.append((valid_y, c, False))  # dashed = near-native AND PB-valid
        if mk == MARK_KS[-1]:                 # rightmost tick: both columns LEFT of the line
            _place_column(mk, near_col + valid_col, "left")
        else:
            _place_column(mk, near_col, "left")
            _place_column(mk, valid_col, "right")

    ax.set_xlabel("Top-k — within the tool's first k ranked poses")
    ax.set_ylabel("% of receptor-ligand complexes recovered")
    ax.set_ylim(0, 100)
    ax.set_xlim(0, kmax + 6)                  # start at 0 to give k = 1's LEFT labels room
    ax.set_xticks([k for k in (1, 5, 10, 15, 20, 25, 30) if k <= kmax])
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)   # dotted vlines are the x refs
    style = [Line2D([0], [0], color="0.35", lw=2.2, ls="-",
                    label="near-native (RMSD ≤ 2 Å)"),
             Line2D([0], [0], color="0.35", lw=2.2, ls="--",
                    label="near-native AND PB-valid (validity-aware)")]
    # Park the legend in the empty band between the AutoDock and DiffDock curves so
    # it clears every line and the right-hand variant labels.
    ax.legend(handles=style, loc="upper right", bbox_to_anchor=(0.99, 0.80),
              fontsize=9, framealpha=0.9, title="line = recovery of")
    ax.set_title("Cumulative recovery within top-k — near-native vs. PB-valid",
                 fontsize=12, fontweight="bold")
    depths = ", ".join(str(k) for k in (stats or {}).get("test_depths", _TOPK_TEST_DEPTHS))
    note = ("Shaded bands: Wilson 95% CI on each rate. Paired exact-McNemar tests "
            f"(Holm-corrected) at pre-specified depths k ∈ {{{depths}}}; full table → "
            "topk_recovery_stats.csv / …_stats.json.")
    fig.text(0.006, 0.006, note, fontsize=7, color="0.4", ha="left", va="bottom")
    fig.tight_layout(rect=(0, 0.025, 1, 1))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# 21/22 — Ranking quality: how faithfully does each tool's OWN pose ranking
# reproduce the oracle ordering of the poses it produced? AutoDock and DiffDock
# are ranked by their native confidence order; EquiBind has no native score, so
# its poses are ranked by gnina affinity (most-negative = rank 1), matching the
# convention of figs 09d/15b. Four quality criteria define the "oracle" order a
# pose SHOULD have (lower = better pose):
#   crystal  — as-placed heavy-atom RMSD to the crystal ligand
#   kabsch   — best-fit RMSD after optimal superposition (internal form only)
#   pbvalid  — PoseBusters validity (valid ranked above invalid)
#   combined — the deployment pick: PB-valid, near the crystal (as-placed RMSD)
#              AND well-formed (Kabsch RMSD). Valid poses first, then ranked by
#              placement+form (as-placed RMSD + Kabsch RMSD)
# Fig 21 shows position-by-position concordance (best pose → rank 1, 2nd-best →
# rank 2, …); fig 22 is the per-complex Kendall tau-b distribution as a whisker
# plot. Both need the FULL (un-collapsed) frame so the gnina EquiBind variant is
# present.
_RANK_QUALITY_SPECS = [
    # (colour key, display label, method key, ranking kind, ranking-axis label)
    # AutoDock = the dominant arm plus its raw counterpart. Keeping at least two
    # AutoDock entries here also keeps the cross-tool Friedman omnibus feasible,
    # which needs three or more tools.
    ("autodock", "AutoDock Vina",           _DOMINANT_AUTODOCK_RAW, "native", "confidence rank"),
    ("autodock_gnina", "AutoDock (gnina-opt)", _DOMINANT_AUTODOCK_OPT,
     "native", "optimized rank"),
    ("diffdock", "DiffDock",                "diffdock",                "native", "confidence rank"),
    ("equibind", "EquiBind (gnina-ranked)", "equibind_unguided_gnina", "gnina",  "gnina-affinity rank"),
]
# (key, display label) for the four ranking-quality criteria, in plot order.
_RANK_QUALITY_CRITERIA = [
    ("crystal",  "Crystal RMSD"),
    ("kabsch",   "Kabsch RMSD"),
    ("pbvalid",  "PB-validity"),
    ("combined", "Combined (valid + close + form)"),
]


def _rank_quality_badness(sub: pd.DataFrame, criterion: str) -> pd.Series:
    """Per-pose 'badness' (LOWER = better pose) under a ranking-quality criterion.

    A tool ranks well when its confidence order (rank 1 = highest confidence)
    tracks ascending badness. NaN badness marks a pose the criterion can't score
    (dropped downstream). ``combined`` is the deployment pick: PB-valid, near the
    crystal (as-placed RMSD) AND well-formed (Kabsch/best-fit RMSD). Its geometric
    badness sums the as-placed and form RMSDs, and a global penalty ≥ the frame's
    max sum pushes every PB-invalid pose below all valid poses within EVERY
    complex (valid first, then close-and-well-formed)."""
    rmsd = pd.to_numeric(sub["rmsd"], errors="coerce")
    if criterion == "crystal":
        return rmsd
    if criterion == "kabsch":
        return pd.to_numeric(sub["bestfit_rmsd"], errors="coerce")
    valid = _to_bool(sub["pb_valid"])
    if criterion == "pbvalid":
        return pd.Series(np.where(valid.to_numpy(), 0.0, 1.0), index=sub.index)
    if criterion == "combined":
        # geometric quality = placement (as-placed RMSD) + form (Kabsch RMSD)
        geo = rmsd + pd.to_numeric(sub["bestfit_rmsd"], errors="coerce")
        arr = geo.to_numpy()
        pen = (float(np.nanmax(arr)) + 1.0) if np.isfinite(arr).any() else 1.0
        return geo + np.where(valid.to_numpy(), 0.0, pen)
    raise ValueError(f"unknown ranking-quality criterion: {criterion}")


def aggregate_rank_quality_tau(df_full: pd.DataFrame,
                               specs=_RANK_QUALITY_SPECS) -> pd.DataFrame:
    """Per (tool, criterion, complex): Kendall tau-b between the tool's ranking
    (rank 1 = highest confidence) and the criterion badness (low = better pose).

    tau = +1 → the tool ranks its poses in perfect quality order; 0 → no better
    than random; −1 → reversed. A complex needs ≥ 3 poses and non-constant rank
    AND badness, so PB-validity tau is defined only where the pose set contains
    BOTH valid and invalid poses (a near-uniformly-valid tool contributes few
    complexes — reflected in the per-box n). Long form; one row per contributing
    complex. Needs the full frame for the gnina EquiBind variant."""
    present = set(df_full["method"].astype(str))
    rows = []
    for ckey, clabel, mkey, kind, _axis in specs:
        if mkey not in present:
            continue
        frame, rank = _rank_series_for(
            df_full[df_full["method"].astype(str) == mkey], kind)
        frame = frame.assign(_rk=pd.to_numeric(rank, errors="coerce"))
        for crit, _crit_label in _RANK_QUALITY_CRITERIA:
            work = frame.assign(
                _bad=pd.to_numeric(_rank_quality_badness(frame, crit), errors="coerce"))
            for (prot, lig), g in work.groupby(["protein", "ligand"]):
                gg = g.dropna(subset=["_rk", "_bad"])
                if len(gg) < 3 or gg["_rk"].nunique() < 2 or gg["_bad"].nunique() < 2:
                    continue
                tau, _p, _n = su.kendall(gg["_rk"].to_numpy(float),
                                         gg["_bad"].to_numpy(float))
                if not np.isfinite(tau):
                    continue
                rows.append({"tool_key": ckey, "tool": clabel, "criterion": crit,
                             "protein": prot, "ligand": lig,
                             "tau": round(float(tau), 4), "n_poses": int(len(gg))})
    return pd.DataFrame(rows)


def _summarize_rank_quality_tau(tau_df: pd.DataFrame) -> pd.DataFrame:
    """Compact per (tool, criterion) summary of the tau distribution: n complexes
    with a defined tau, median, IQR, and the % that are positive (ranking beats
    chance for that complex)."""
    rows = []
    for ckey, clabel, _mk, _kind, _ax in _RANK_QUALITY_SPECS:
        for crit, crit_label in _RANK_QUALITY_CRITERIA:
            s = tau_df[(tau_df["tool"] == clabel) & (tau_df["criterion"] == crit)]["tau"]
            s = s.dropna()
            if s.empty:
                continue
            rows.append({
                "tool": clabel, "criterion": crit, "criterion_label": crit_label,
                "n_complexes": int(len(s)),
                "median_tau": round(float(s.median()), 4),
                "q25_tau": round(float(s.quantile(0.25)), 4),
                "q75_tau": round(float(s.quantile(0.75)), 4),
                "pct_positive": round(100.0 * float((s > 0).mean()), PCT_DECIMALS),
            })
    return pd.DataFrame(rows)


def _stats_rank_quality(tau_df: pd.DataFrame) -> dict:
    """Per criterion: paired cross-tool test on the per-complex tau (the tools
    that scored the same complex are the paired units) — Friedman + Kendall's W
    omnibus and Wilcoxon signed-rank (Holm) pairwise, via the shared helper.
    Returns {criterion: paired_continuous(...)-dict (+ n_complete)} and skips a
    criterion with < 3 jointly-scored complexes or < 2 tools."""
    out = {}
    tools = [s[1] for s in _RANK_QUALITY_SPECS]
    for crit, _clabel in _RANK_QUALITY_CRITERIA:
        sub = tau_df[tau_df["criterion"] == crit]
        if sub.empty:
            continue
        wide = sub.pivot_table(index=["protein", "ligand"], columns="tool",
                               values="tau", aggfunc="first")
        cols = [t for t in tools if t in wide.columns]
        wide = wide[cols].dropna()
        if len(wide) < 3 or len(cols) < 2:
            out[crit] = {"insufficient": True, "n_complete": int(len(wide))}
            continue
        res = su.paired_continuous({t: wide[t].to_numpy(float) for t in cols}, labels=cols)
        res["n_complete"] = int(len(wide))
        out[crit] = res
    return out


def aggregate_rank_position_concordance(df_full: pd.DataFrame, specs, top_n: int):
    """Two frames feeding fig 21.

    diag  — for the two continuous criteria (crystal, kabsch) and each rank
      position k (1..top_n): the % of complexes whose pose at the tool's rank k
      is ALSO the k-th best pose by that metric (best → rank 1, 2nd-best →
      rank 2, …). Denominator = complexes with ≥ k metric-defined poses.
    valid — for PB-validity: the % of complexes whose tool-rank-k pose is
      PB-valid, plus the oracle ceiling (% of complexes with ≥ k valid poses =
      the best any ranking could do by front-loading valid poses). Identity
      matching is undefined for a binary label, so validity uses this rate view.
    """
    present = set(df_full["method"].astype(str))
    diag_rows, val_rows = [], []
    for ckey, clabel, mkey, kind, _axis in specs:
        if mkey not in present:
            continue
        frame, rank = _rank_series_for(
            df_full[df_full["method"].astype(str) == mkey], kind)
        frame = frame.assign(_rk=pd.to_numeric(rank, errors="coerce"))
        # continuous-metric diagonal (identity) agreement
        for crit, col in (("crystal", "rmsd"), ("kabsch", "bestfit_rmsd")):
            per_cx = list(frame.dropna(subset=["_rk", col]).groupby(["protein", "ligand"]))
            for k in range(1, top_n + 1):
                match = tot = 0
                for _, gg in per_cx:
                    if len(gg) < k:
                        continue
                    at_k = gg[gg["_rk"] == k]
                    if at_k.empty:
                        continue
                    order = gg.sort_values(col, kind="mergesort")
                    tot += 1
                    if order.iloc[k - 1]["pose_name"] == at_k.iloc[0]["pose_name"]:
                        match += 1
                if tot:
                    diag_rows.append({
                        "tool_key": ckey, "tool": clabel, "criterion": crit,
                        "rank": k, "pct_match": round(100.0 * match / tot, PCT_DECIMALS),
                        "n_pairs": tot})
        # PB-validity: achieved validity-by-rank vs the tool's OVERALL valid rate
        # (the random-ordering reference — a tool that ranks validity well sits
        # above it at rank 1 and declines below it at deep ranks). A front-loading
        # "ceiling" is NOT a valid per-position upper bound: front-loading exhausts
        # the valid poses early, so its deep-rank validity is a floor, not a ceiling.
        vframe = frame.assign(_v=_to_bool(frame["pb_valid"]))
        overall_valid = 100.0 * float(vframe["_v"].mean()) if len(vframe) else float("nan")
        for k in range(1, top_n + 1):
            at_k = vframe[vframe["_rk"] == k]
            npair = at_k.groupby(["protein", "ligand"]).ngroups
            if not npair:
                continue
            val_rows.append({
                "tool_key": ckey, "tool": clabel, "rank": k,
                "pct_valid": round(100.0 * float(at_k["_v"].mean()), PCT_DECIMALS),
                "pct_random": round(overall_valid, PCT_DECIMALS),
                "n_pairs": npair})
    return pd.DataFrame(diag_rows), pd.DataFrame(val_rows)


def plot_rank_position_concordance(diag_df: pd.DataFrame, valid_df: pd.DataFrame,
                                   top_n: int, out: Path) -> None:
    """Fig 21 — position-by-position ranking concordance, three panels.

    A/B (Crystal RMSD, Kabsch RMSD): grouped bars per tool of the % of complexes
    where the tool's rank-k pose is the k-th best pose by that metric (the
    literal "best → rank 1, 2nd-best → rank 2 …" question). C (PB-validity):
    achieved % PB-valid at each rank (solid) vs the front-loading ceiling
    (dashed) per tool."""
    if diag_df.empty and valid_df.empty:
        return
    ranks = list(range(1, top_n + 1))
    tools = [(ckey, clabel) for ckey, clabel, *_ in _RANK_QUALITY_SPECS]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

    def _bars(ax, crit, title):
        present = [(ck, cl) for ck, cl in tools
                   if not diag_df[(diag_df["tool"] == cl)
                                  & (diag_df["criterion"] == crit)].empty]
        w = 0.8 / max(1, len(present))
        for mi, (ck, cl) in enumerate(present):
            g = diag_df[(diag_df["tool"] == cl) & (diag_df["criterion"] == crit)]
            g = g.set_index("rank").reindex(ranks)
            offset = (mi - (len(present) - 1) / 2) * w
            ax.bar([r + offset for r in ranks], g["pct_match"], w,
                   label=cl, color=TOOL_COLORS.get(ck), edgecolor="black", alpha=0.85)
        ax.set_xticks(ranks)
        ax.set_xlabel("Rank position k (tool's own ranking)")
        ax.set_ylabel("% of complexes where rank-k pose\nis the k-th best pose")
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=8, title="ranked by native / gnina")

    _bars(axes[0], "crystal", "Crystal RMSD (as-placed)")
    _bars(axes[1], "kabsch", "Kabsch RMSD (best-fit, form only)")

    axC = axes[2]
    for ck, cl in tools:
        g = valid_df[valid_df["tool"] == cl].set_index("rank").reindex(ranks)
        if g["pct_valid"].isna().all():
            continue
        col = TOOL_COLORS.get(ck)
        axC.plot(ranks, g["pct_valid"], marker="o", lw=2, color=col, label=cl)
        axC.plot(ranks, g["pct_random"], ls=":", lw=1.2, color=col, alpha=0.7)
    axC.set_xticks(ranks)
    axC.set_ylim(0, 105)
    axC.set_xlabel("Rank position k (tool's own ranking)")
    axC.set_ylabel("% of rank-k poses that are PoseBusters-valid")
    axC.set_title("PoseBusters validity by rank")
    axC.grid(axis="y", alpha=0.3)
    axC.legend(fontsize=8)
    axC.text(0.98, 0.02, "dotted = random baseline (tool's overall valid rate)",
             transform=axC.transAxes, ha="right", va="bottom", fontsize=7.5,
             style="italic", color="0.35")

    _label_panels(axes)
    fig.suptitle(_vt("How faithfully does each tool rank its poses? "
                     "(best pose → rank 1, 2nd-best → rank 2, …)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_rank_quality_whiskers(tau_df: pd.DataFrame, out: Path) -> None:
    """Fig 22 — per-complex Kendall tau-b (tool ranking vs oracle order) as a
    grouped whisker plot: one group per criterion (Crystal RMSD, Kabsch RMSD,
    PB-validity, Combined), one box per tool. Higher = the tool ranks its poses
    closer to true quality order; the dashed line at 0 is chance. The cross-tool
    paired tests live in rank_quality_report.txt, not on the figure."""
    if tau_df.empty:
        return
    tools = [(ck, cl) for ck, cl, *_ in _RANK_QUALITY_SPECS]
    crits = _RANK_QUALITY_CRITERIA
    fig, ax = plt.subplots(figsize=(13, 5.6))
    group_w, box_w = 1.0, 0.8 / max(1, len(tools))
    positions, box_data, box_colors, n_labels = [], [], [], []
    centers = []
    for gi, (crit, clabel) in enumerate(crits):
        base = gi * group_w
        centers.append(base)
        for ti, (ck, cl) in enumerate(tools):
            s = tau_df[(tau_df["tool"] == cl) & (tau_df["criterion"] == crit)]["tau"].dropna()
            pos = base + (ti - (len(tools) - 1) / 2) * box_w
            positions.append(pos)
            box_data.append(s.to_numpy(float) if len(s) else np.array([np.nan]))
            box_colors.append(TOOL_COLORS.get(ck))
            n_labels.append((pos, len(s)))
    bp = ax.boxplot(box_data, positions=positions, widths=box_w * 0.9,
                    patch_artist=True, showfliers=False, showmeans=True,
                    meanprops=dict(marker="D", markerfacecolor="white",
                                   markeredgecolor="black", markersize=4),
                    medianprops=dict(color="black", lw=1.4))
    for patch, col in zip(bp["boxes"], box_colors):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
        patch.set_edgecolor("black")
    ax.axhline(0.0, ls="--", lw=1, color="0.4", zorder=0)
    for pos, n in n_labels:
        ax.text(pos, -1.06, f"n={n}", ha="center", va="top", fontsize=6.5, color="0.35")
    ax.set_xticks(centers)
    ax.set_xticklabels([cl for _c, cl in crits])
    ax.set_ylim(-1.12, 1.12)
    ax.set_ylabel("Kendall τ-b: tool ranking vs oracle order\n(+1 = perfect, 0 = chance, −1 = reversed)")
    ax.set_xlabel("Ranking-quality criterion (oracle ordering of the tool's own poses)")
    ax.grid(axis="y", alpha=0.3)
    # tool legend
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=TOOL_COLORS.get(ck),
                             edgecolor="black", alpha=0.75) for ck, _cl in tools]
    ax.legend(handles, [cl for _ck, cl in tools], fontsize=8,
              title="ranked by native confidence / gnina affinity", loc="upper right")

    fig.suptitle(_vt("Ranking quality per tool — does each tool order its poses "
                     "by true quality?"), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── 21a-d: rank-vs-oracle individual line charts + full text report ──────────────
# A more informative companion to the fig-21 identity-match bars: for each rank
# position k, the ORACLE rank (1 = best pose by the criterion) of the pose the tool
# actually placed at rank k. A perfect ranking sits on the diagonal (oracle rank = k);
# a useless one collapses onto the flat random baseline (≈ N/2). One standalone figure
# per criterion; the statistics go to rank_quality_report.txt (not onto the graphs).
_RVO_CRITERIA = [
    ("crystal",  "rmsd",         "Crystal RMSD"),
    ("kabsch",   "bestfit_rmsd", "Kabsch RMSD"),
    ("combined", None,           "Combined (valid + close + form)"),
]
# Depth for the rank-vs-oracle figures (21a/b/d): the tools produce ~30 poses per
# complex, so the axes span the full 1..30 rather than the report's --top-n cap.
_RVO_MAX_RANK = 30
# Depth for rank_quality_report.txt's per-rank tables (section 3 RANK vs ORACLE,
# 4 IDENTITY MATCH, 5 PB-VALIDITY BY RANK). Decoupled from --top-n (which caps the
# 21/21c figures) so the text report spans the full pose set. Kept equal to
# _RVO_MAX_RANK so section 3 (rvo) and sections 4/5 (identity/validity) share a depth.
_RANK_REPORT_DEPTH = _RVO_MAX_RANK


def aggregate_rank_vs_oracle(df_full: pd.DataFrame, specs, top_n: int) -> pd.DataFrame:
    """Per (tool, criterion, tool-rank k): the oracle rank (1 = best pose by that
    criterion) of the pose the tool placed at rank k, summarised over complexes
    (median + IQR + mean) with the random-ranking baseline (median (N+1)/2).
    Continuous criteria only (crystal, kabsch, combined) — PB-validity has no oracle
    ordering. Needs the full frame for the gnina EquiBind variant."""
    present = set(df_full["method"].astype(str))
    rows = []
    for ckey, clabel, mkey, kind, _axis in specs:
        if mkey not in present:
            continue
        frame, rank = _rank_series_for(
            df_full[df_full["method"].astype(str) == mkey], kind)
        frame = frame.assign(_rk=pd.to_numeric(rank, errors="coerce"))
        for crit, _col, crit_label in _RVO_CRITERIA:
            work = frame.assign(
                _bad=pd.to_numeric(_rank_quality_badness(frame, crit), errors="coerce"))
            work = work.dropna(subset=["_rk", "_bad"])
            by_k = {k: [] for k in range(1, top_n + 1)}
            half_n = []
            for _, g in work.groupby(["protein", "ligand"]):
                g = g.sort_values("_bad", kind="mergesort")
                g = g.assign(_o=np.arange(1, len(g) + 1))
                half_n.append((len(g) + 1) / 2.0)
                for k in range(1, top_n + 1):
                    r = g[g["_rk"] == k]
                    if not r.empty:
                        by_k[k].append(int(r["_o"].iloc[0]))
            base = float(np.median(half_n)) if half_n else float("nan")
            for k in range(1, top_n + 1):
                vals = by_k[k]
                if not vals:
                    continue
                a = np.asarray(vals, float)
                rows.append({
                    "tool_key": ckey, "tool": clabel, "criterion": crit,
                    "criterion_label": crit_label, "rank": k,
                    "median_oracle_rank": round(float(np.median(a)), 2),
                    "q25": round(float(np.percentile(a, 25)), 2),
                    "q75": round(float(np.percentile(a, 75)), 2),
                    "mean_oracle_rank": round(float(a.mean()), 2),
                    "n_pairs": len(vals),
                    "random_baseline": round(base, 2)})
    return pd.DataFrame(rows)


def plot_rank_vs_oracle_criterion(rvo: pd.DataFrame, crit: str, out: Path) -> None:
    """One standalone square rank-vs-rank figure: median oracle rank of the tool's
    rank-k pose vs rank k, per tool. Both axes are ranks — identical integer ticks
    starting at 1 — so the perfect-ranking diagonal (oracle rank = k) reads at 45°
    and the random level sits near the top edge (oracle rank ≈ N/2 ≥ top-N). A curve
    hugging the diagonal ranks well; one riding the ceiling is no better than random.
    Per-complex IQR is in rank_quality_report.txt. Points whose median oracle rank
    exceeds top-N clip at the ceiling (they are effectively random)."""
    sub = rvo[rvo["criterion"] == crit]
    if sub.empty:
        return
    label = sub["criterion_label"].iloc[0]
    ranks = sorted(int(r) for r in sub["rank"].unique())
    kmax = max(ranks)
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    ax.plot([1, kmax], [1, kmax], ls="--", lw=1.4, color="0.45", zorder=1,
            label="perfect ranking (oracle rank = k)")
    for ckey, clabel, *_rest in _RANK_QUALITY_SPECS:
        g = sub[sub["tool"] == clabel].set_index("rank").reindex(ranks)
        if g["median_oracle_rank"].isna().all():
            continue
        ax.plot(ranks, g["median_oracle_rank"], marker="o", lw=2,
                color=TOOL_COLORS.get(ckey), label=clabel, zorder=3)
    ax.set_xlim(1, kmax)
    ax.set_ylim(1, kmax)
    # both axes share the same rank ticks; stride out when the depth is large
    if kmax <= 16:
        ticks = list(range(1, kmax + 1))
    else:
        ticks = [1] + list(range(5, kmax + 1, 5))
        if kmax not in ticks:
            ticks.append(kmax)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Rank position k (tool's own ranking)")
    ax.set_ylabel("Oracle rank of the tool's rank-k pose\n(1 = best pose by this metric)")
    ax.grid(alpha=0.3)
    # legend ABOVE the plot (title sits above it)
    ax.legend(fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=2, borderaxespad=0.3, columnspacing=1.2,
              title="ranked by native confidence / gnina affinity")
    ax.set_title(_vt(f"Rank vs oracle — {label}"), pad=60)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_validity_by_rank_individual(valid_df: pd.DataFrame, top_n: int,
                                     out: Path) -> None:
    """Standalone PB-validity-by-rank line chart (achieved solid, ceiling dashed)."""
    if valid_df.empty:
        return
    ranks = list(range(1, top_n + 1))
    fig, ax = plt.subplots(figsize=(7.8, 5.4))
    for ckey, clabel, *_rest in _RANK_QUALITY_SPECS:
        g = valid_df[valid_df["tool"] == clabel].set_index("rank").reindex(ranks)
        if g["pct_valid"].isna().all():
            continue
        col = TOOL_COLORS.get(ckey)
        ax.plot(ranks, g["pct_valid"], marker="o", lw=2, color=col, label=clabel)
        ax.plot(ranks, g["pct_random"], ls=":", lw=1.2, color=col, alpha=0.7)
    ax.set_xticks(ranks)
    ax.set_xlim(1, top_n)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Rank position k (tool's own ranking)")
    ax.set_ylabel("% of rank-k poses that are PoseBusters-valid")
    ax.set_title(_vt("Rank vs oracle — PoseBusters validity by rank"))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper right",
              title="ranked by native confidence / gnina affinity")
    ax.text(0.98, 0.02, "dotted = random baseline (tool's overall valid rate);\n"
            "good ranking starts above it and declines toward it",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
            style="italic", color="0.35")
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def _validity_mix(df_full: pd.DataFrame, specs) -> pd.DataFrame:
    """Per tool: how many complexes are all-valid / all-invalid / mixed — explains
    why the PB-validity tau (fig 22) is defined only on the mixed complexes."""
    present = set(df_full["method"].astype(str))
    rows = []
    for _ck, clabel, mkey, _kind, _ax in specs:
        if mkey not in present:
            continue
        m = df_full[df_full["method"].astype(str) == mkey].copy()
        m["_v"] = _to_bool(m["pb_valid"])
        g = m.groupby(["protein", "ligand"])["_v"]
        nv, nt = g.sum(), g.size()
        rows.append({"tool": clabel, "pairs": int(g.ngroups),
                     "all_valid": int((nv == nt).sum()),
                     "all_invalid": int((nv == 0).sum()),
                     "mixed": int(((nv > 0) & (nv < nt)).sum())})
    return pd.DataFrame(rows)


def _write_rank_quality_report(path: Path, *, tau_df, summary_df, stats, rvo_df,
                               diag_df, valid_df, mix_df, top_n) -> None:
    """Full plain-text ranking-quality report — everything that used to be crammed
    onto the figures (and more), so the graphs stay clean."""
    L = []
    bar = "=" * 78
    L += [bar, "RANKING-QUALITY REPORT", bar,
          "How faithfully does each tool's OWN pose ranking reproduce the oracle",
          "ordering of the poses it produced? (benchmark set)", ""]
    L += ["Tools & ranking basis:"]
    for _ck, clabel, mkey, kind, axl in _RANK_QUALITY_SPECS:
        L.append(f"  - {clabel:24s}: {axl}  [{mkey}, kind={kind}]")
    L += ["",
          "Criteria (oracle order of a tool's poses; lower badness = better pose):",
          f"  - Crystal RMSD : as-placed heavy-atom RMSD to the {_ref_ligand_phrase()}",
          "  - Kabsch RMSD  : best-fit RMSD after optimal superposition (form only)",
          "  - PB-validity  : PoseBusters-valid ranked above invalid (binary)",
          "  - Combined     : deployment pick — valid first, then by placement+form",
          f"                   (as-placed {_ref_noun()} RMSD + Kabsch/best-fit RMSD)",
          ""]
    if REFERENCE_CONVENTION != "instance":
        L += ["  " + _ref_convention_note(), ""]

    L += [bar, "1. KENDALL tau-b SUMMARY  (per-complex tau; +1 perfect, 0 chance, -1 reversed)", bar]
    L.append(f"{'tool':24s} {'criterion':10s} {'n':>4s} {'median':>7s} "
             f"{'q25':>7s} {'q75':>7s} {'%>0':>5s}")
    for _r, row in summary_df.iterrows():
        L.append(f"{row['tool']:24s} {row['criterion']:10s} {int(row['n_complexes']):4d} "
                 f"{row['median_tau']:7.3f} {row['q25_tau']:7.3f} {row['q75_tau']:7.3f} "
                 f"{row['pct_positive']:5.0f}")
    L.append("")

    L += [bar, "2. CROSS-TOOL PAIRED TESTS  (per criterion; units = shared complexes)", bar]
    for crit, clabel in _RANK_QUALITY_CRITERIA:
        st = (stats or {}).get(crit) or {}
        L.append(f"[{clabel}]")
        if not st or st.get("insufficient"):
            L.append(f"  insufficient jointly-scored complexes (n={st.get('n_complete', 0)})")
            L.append("")
            continue
        om = st["omnibus"]
        L.append(f"  Friedman chi2={om['chi2']:.2f}, p={su.fmt_p(om['p'])}"
                 f"{su.p_stars(om['p'])}, Kendall's W={om['kendall_w']:.3f}, n={om['n']}")
        L.append(f"  medians: " + ", ".join(f"{k}={v:.3f}" for k, v in st["medians"].items()))
        for pr in st["pairwise"]:
            L.append(f"    {pr['a']:24s} vs {pr['b']:24s}  "
                     f"rank-biserial={pr['rank_biserial']:+.3f}  "
                     f"p_holm={su.fmt_p(pr['p_holm'])}{pr['star']}  (n={pr['n']})")
        L.append("")

    L += [bar, "3. RANK vs ORACLE  (oracle rank of the tool's rank-k pose; median [IQR];",
          f"   perfect = k, random ~ N/2; shown to rank {top_n}, figures 21a/b/d + "
          "rank_vs_oracle.csv go to rank 30)", bar]
    for crit, _col, crit_label in _RVO_CRITERIA:
        sub = rvo_df[rvo_df["criterion"] == crit]
        if sub.empty:
            continue
        base = float(sub["random_baseline"].median())
        L.append(f"[{crit_label}]   (random baseline ~ {base:.0f})")
        L.append("  " + "tool".ljust(24) + "".join(f"k{k:<9d}" for k in range(1, top_n + 1)))
        for _ck, clabel, *_ in _RANK_QUALITY_SPECS:
            g = sub[sub["tool"] == clabel].set_index("rank")
            cells = []
            for k in range(1, top_n + 1):
                if k in g.index:
                    cells.append(f"{int(g.loc[k, 'median_oracle_rank'])}"
                                 f"[{int(g.loc[k, 'q25'])}-{int(g.loc[k, 'q75'])}]".ljust(10))
                else:
                    cells.append("-".ljust(10))
            L.append("  " + clabel.ljust(24) + "".join(cells))
        L.append("")

    L += [bar, "4. IDENTITY MATCH  (% of complexes where rank-k pose IS the k-th best pose)", bar]
    for crit in ("crystal", "kabsch"):
        sub = diag_df[diag_df["criterion"] == crit]
        if sub.empty:
            continue
        L.append(f"[{crit}]")
        L.append("  " + "tool".ljust(24) + "".join(f"k{k:<6d}" for k in range(1, top_n + 1)))
        for _ck, clabel, *_ in _RANK_QUALITY_SPECS:
            g = sub[sub["tool"] == clabel].set_index("rank")
            cells = [(f"{g.loc[k, 'pct_match']:.0f}%".ljust(7) if k in g.index else "-".ljust(7))
                     for k in range(1, top_n + 1)]
            L.append("  " + clabel.ljust(24) + "".join(cells))
        L.append("")

    L += [bar, "5. PB-VALIDITY BY RANK  (% valid at rank k / tool overall valid rate =",
          "   random-order reference; good ranking is above it early, below it deep)", bar]
    L.append("  " + "tool".ljust(24) + "".join(f"k{k:<11d}" for k in range(1, top_n + 1)))
    for _ck, clabel, *_ in _RANK_QUALITY_SPECS:
        g = valid_df[valid_df["tool"] == clabel].set_index("rank")
        cells = [(f"{g.loc[k, 'pct_valid']:.0f}/{g.loc[k, 'pct_random']:.0f}".ljust(12)
                  if k in g.index else "-".ljust(12)) for k in range(1, top_n + 1)]
        L.append("  " + clabel.ljust(24) + "".join(cells))
    L.append("")

    L += [bar, "6. PB-VALIDITY COMPLEX MIX  (why the fig-22 PB-validity box n is small)", bar,
          "   tau needs BOTH valid and invalid poses in a complex; all-valid/all-invalid",
          "   complexes are constant -> tau undefined -> dropped (n = 'mixed').", ""]
    L.append(f"  {'tool':24s} {'pairs':>6s} {'all-valid':>10s} {'all-invalid':>12s} {'mixed(=n)':>10s}")
    for _r, row in mix_df.iterrows():
        L.append(f"  {row['tool']:24s} {int(row['pairs']):6d} {int(row['all_valid']):10d} "
                 f"{int(row['all_invalid']):12d} {int(row['mixed']):10d}")
    L.append("")

    Path(path).write_text("\n".join(L))


# ── Optimization benefit by rank — does refinement help early ranks more? ────────
# For DiffDock, smina/gnina keep the confidence rank on refined geometry, so rank r is
# the SAME underlying pose across raw/smina/gnina and the raw→optimized gap at each rank
# is optimization's per-rank benefit. EquiBind has no confidence rank for its raw poses,
# so poses are paired across raw/smina/gnina by GENERATION ORDER (the same generated
# sample in raw vs refined form); smina/gnina do add an affinity ranking, but generation
# order is the only axis shared with raw. Three criteria are tracked separately because
# they answer the "early vs late" question differently (RMSD barely moves, PB-validity
# jumps, the combined docking-success gain is largest at DiffDock's early ranks).
_OPT_BENEFIT_TOOLS = [
    # (tool label, x-axis label, rank kind, {optimizer: method_key})
    ("DiffDock", "DiffDock confidence rank", "native",
     {"raw": "diffdock", "smina": "diffdock_smina", "gnina": "diffdock_gnina"}),
    ("EquiBind", "EquiBind pose (generation order)", "generation",
     {"raw": "equibind_unguided_raw", "smina": "equibind_unguided_smina",
      "gnina": "equibind_unguided_gnina"}),
]
_OPT_BENEFIT_CRITERIA = [
    ("near_%", "Near-native (RMSD ≤ 2 Å)"),
    ("pbv_%",  "PoseBusters-valid"),
    ("both_%", "Docking success (≤ 2 Å AND PB-valid)"),
]
_OPT_RANK_GROUPS = [(1, 10), (11, 20), (21, 30)]
_OPT_C_RAW, _OPT_C_SMI, _OPT_C_GNI = "#f0a35e", "#d95f02", "#7f3b08"


def aggregate_optimization_benefit(df: pd.DataFrame, max_rank: int = 30,
                                   thr: float = 2.0) -> pd.DataFrame:
    """Per (tool, optimizer raw/smina/gnina, rank r): the fraction of poses at that rank
    that are near-native (RMSD ≤ ``thr``), PB-valid, and both, plus the pose count n.

    Rank r pairs the SAME underlying pose across the three optimizer variants: DiffDock
    uses its confidence rank (preserved through refinement); EquiBind — which has no raw
    confidence rank — uses generation-order position (the same generated sample in raw vs
    smina/gnina-refined form). The raw→optimized gap at each rank is optimization's
    per-rank benefit. Needs the FULL frame (all variants present)."""
    present = set(df["method"].astype(str))
    rows = []
    for tool, _xlab, kind, trio in _OPT_BENEFIT_TOOLS:
        for opt, mkey in trio.items():
            if mkey not in present:
                continue
            sub = df[df["method"].astype(str) == mkey].copy()
            if kind == "native":
                sub["_rk"] = pd.to_numeric(sub["rank"], errors="coerce")
                sub = sub.sort_values("_rk").drop_duplicates(
                    ["protein", "ligand", "_rk"])       # DiffDock writes rank1 twice
            else:                                        # EquiBind: generation-order pos
                sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
                sub = sub.sort_values(["protein", "ligand", "_gi"], kind="mergesort")
                sub["_rk"] = sub.groupby(["protein", "ligand"]).cumcount() + 1
            sub = sub.dropna(subset=["_rk"])
            sub = sub[sub["_rk"] <= max_rank]
            if sub.empty:
                continue
            near = pd.to_numeric(sub["rmsd"], errors="coerce") <= thr
            pbv = _to_bool(sub["pb_valid"])
            sub = sub.assign(_near=near, _pbv=pbv, _both=near & pbv)
            g = sub.groupby("_rk")
            agg = (g[["_near", "_pbv", "_both"]].mean() * 100).rename(
                columns={"_near": "near_%", "_pbv": "pbv_%", "_both": "both_%"})
            agg["n"] = g.size()
            agg.insert(0, "optimizer", opt)
            agg.insert(0, "tool", tool)
            agg.index.name = "rank"
            rows.append(agg.reset_index())
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _opt_pooled_rate(rec_v: pd.DataFrame, col: str, lo: int, hi: int) -> float:
    """Pose-count-weighted (pooled) rate over a rank group [lo, hi] from a per-rank slice
    of one (tool, optimizer): the fraction over ALL poses in the group."""
    d = rec_v[(rec_v["rank"] >= lo) & (rec_v["rank"] <= hi)]
    w = d["n"].to_numpy(dtype=float)
    v = d[col].to_numpy(dtype=float)
    tot = np.nansum(w)
    return float(np.nansum(v * w) / tot) if tot > 0 else float("nan")


def _opt_tools_in(rec: pd.DataFrame):
    seen = set(rec["tool"])
    return [t for (t, *_r) in _OPT_BENEFIT_TOOLS if t in seen]


# Map each 15c/15d criterion column to the endpoint tested in
# optimization_benefit_stats.py; the derived "both_%" is never tested.
_OPT_CRIT2EP = {"near_%": "near_native", "pbv_%": "pb_valid", "both_%": None}


def _annotate_15c_panel(ax, ranks, r, s, g, a_df, tool, col) -> None:
    """Overlay question-(a) paired-test results on a 15c panel: significance stars at
    the pre-specified tested ranks (raw→gnina exact McNemar, Holm) + a corner box with
    the primary rank-1 contrast. ``a_df`` is optimization_benefit_stats' pairwise table.
    The derived both_% column is labelled 'not tested'."""
    ep = _OPT_CRIT2EP.get(col)
    if ep is None:
        ax.text(0.97, 0.96, "derived (near ∧ PB-valid)\n— not tested",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.3,
                color="0.45", style="italic")
        return
    sub = a_df[(a_df["tool"] == tool) & (a_df["endpoint"] == ep)]
    if sub.empty:
        return
    rk_int = ranks.astype(int)
    arrs = [a for a in (s, g) if a is not None]
    top = np.nanmax(np.vstack(arrs), axis=0) if arrs else r
    span = ax.get_ylim()[1] - ax.get_ylim()[0]
    for _, rr in sub.iterrows():
        rk = int(rr["rank"])
        hit = np.where(rk_int == rk)[0]
        if not hit.size:
            continue
        star = str(rr["gnina_star"])
        ax.annotate(star, xy=(rk, float(top[hit[0]]) + 0.02 * span), ha="center",
                    va="bottom", fontsize=7.5 if star != "ns" else 6.0,
                    color=_OPT_C_GNI if star != "ns" else "0.55",
                    fontweight="bold", zorder=5)
    prim = sub[sub["primary"]]
    if not prim.empty:
        rr = prim.iloc[0]
        n_sig = int((sub["gnina_mcnemar_holm"] < 0.05).sum())
        line = (f"raw→gnina @r1: {rr['gnina_minus_raw_pp']:+.1f} pp  {rr['gnina_star']}\n"
                f"McNemar p_holm={su.fmt_p(rr['gnina_mcnemar_holm'])} · "
                f"{n_sig}/{len(sub)} tested ranks sig")
        yloc, vloc = (0.5, "center") if ep == "pb_valid" else (0.96, "top")
        ax.text(0.97, yloc, line, transform=ax.transAxes, ha="right", va=vloc,
                fontsize=6.4, color="0.12",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))


def _annotate_15d_panel(ax, b_df, tool, col) -> None:
    """Overlay question-(b) results on a 15d panel: the optimizer×rank interaction
    (gnina, cluster-bootstrap, BH-corrected) with a plain-language trend, plus the
    raw-only rank slope (DiffDock = confidence gradient; EquiBind = negative control).
    ``b_df`` is optimization_benefit_stats' interaction table."""
    ep = _OPT_CRIT2EP.get(col)
    if ep is None:
        ax.text(0.03, 0.97, "derived — not tested", transform=ax.transAxes,
                ha="left", va="top", fontsize=6.3, color="0.45", style="italic")
        return
    gn = b_df[(b_df["tool"] == tool) & (b_df["endpoint"] == ep)
              & (b_df["contrast"] == "gnina_vs_raw")]
    if gn.empty:
        return
    rr = gn.iloc[0]
    b3, star = rr["interaction_logodds_per_rank"], str(rr["star"])
    if star != "ns":
        trend = "benefit grows with rank" if b3 > 0 else "benefit shrinks with rank"
    else:
        trend = "benefit ~flat across ranks"
    ctl = "  (neg. control)" if tool == "EquiBind" else ""
    line = (f"opt×rank (gnina): {b3:+.3f}/rank\n"
            f"p_BH={su.fmt_p(rr['interaction_p_bh'])} {star} — {trend}\n"
            f"raw slope {rr['raw_rank_slope_logodds']:+.3f}/rank "
            f"{rr['raw_slope_star']}{ctl}")
    # top-left: these bars rise from 0 (leftmost group shortest for PB-valid), so the
    # upper-left corner clears them where a bottom anchor would sit inside the bars.
    ax.text(0.03, 0.97, line, transform=ax.transAxes, ha="left", va="top",
            fontsize=6.3, color="0.12",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.7", alpha=0.9))


def _load_opt_benefit_stats(out_dir: Path) -> dict | None:
    """Load the 15c/15d test results (written by optimization_benefit_stats.py) so a
    plain report run RE-ANNOTATES the figures instead of silently dropping the stats.
    Returns {"a": pairwise_df, "b": interaction_df} or None if not computed yet.

    The curves themselves are always recomputed from the current data; only the
    annotation layer is reused, so if the underlying poses changed since the stats were
    last computed, re-run optimization_benefit_stats.py to refresh the numbers."""
    ap = out_dir / "optimization_benefit_stats_pairwise.csv"
    bp = out_dir / "optimization_benefit_stats_interaction.csv"
    if not (ap.exists() and bp.exists()):
        return None
    try:
        a = pd.read_csv(ap)
        b = pd.read_csv(bp)
        if "primary" in a.columns:
            a["primary"] = a["primary"].astype(bool)
        return {"a": a, "b": b}
    except Exception as exc:                             # pragma: no cover
        print(f"  (15c/15d stats present but unreadable — drawing un-annotated: {exc})")
        return None


def plot_optimization_benefit_by_rank(rec: pd.DataFrame, out: Path,
                                      stats: dict | None = None) -> None:
    """Fig 15c — per-rank benefit of pose optimization (smina/gnina) over the raw pose:
    rows = DiffDock / EquiBind, cols = the three criteria. Each panel draws the raw /
    smina / gnina per-rank rates and shades the raw→smina gap as the benefit; the rank
    groups 1–10, 11–20, 21–30 are lightly shaded so the trend reads against them (the +pp
    per group is the companion fig 15d)."""
    if rec is None or rec.empty:
        return
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    tools = _opt_tools_in(rec)
    if not tools:
        return
    grp_colors = ["#2ca02c", "#1f77b4", "#e08214"]
    tool_xlab = {t: xl for (t, xl, *_r) in _OPT_BENEFIT_TOOLS}
    fig, axes = plt.subplots(len(tools), 3, figsize=(15.5, 4.7 * len(tools)),
                             squeeze=False)
    for ri, tool in enumerate(tools):
        sub = rec[rec["tool"] == tool]
        raw = sub[sub["optimizer"] == "raw"].set_index("rank").sort_index()
        smi = sub[sub["optimizer"] == "smina"].set_index("rank").sort_index()
        gni = sub[sub["optimizer"] == "gnina"].set_index("rank").sort_index()
        ranks = raw.index.to_numpy(dtype=float)
        for ci, (col, title) in enumerate(_OPT_BENEFIT_CRITERIA):
            ax = axes[ri][ci]
            r = raw[col].reindex(ranks).to_numpy(dtype=float)
            s = smi[col].reindex(ranks).to_numpy(dtype=float) if not smi.empty else None
            g = gni[col].reindex(ranks).to_numpy(dtype=float) if not gni.empty else None
            for (lo, hi), gc in zip(_OPT_RANK_GROUPS, grp_colors):
                ax.axvspan(lo - 0.5, hi + 0.5, color=gc, alpha=0.05, lw=0, zorder=0)
            if s is not None:
                ax.fill_between(ranks, r, s, where=(s >= r), color=_OPT_C_SMI,
                                alpha=0.15, lw=0)
            ax.plot(ranks, r, color=_OPT_C_RAW, lw=2.0, marker="o", ms=3, zorder=3)
            if s is not None:
                ax.plot(ranks, s, color=_OPT_C_SMI, lw=2.2, marker="o", ms=3, zorder=3)
            if g is not None:
                ax.plot(ranks, g, color=_OPT_C_GNI, lw=1.6, ls="--", zorder=3)
            ax.set_title(f"{tool} — {title}", fontsize=10, fontweight="bold")
            ax.set_xlabel(tool_xlab[tool])
            ax.set_ylabel("% of receptor-ligand complexes")
            ax.set_xlim(1, float(ranks.max()) if len(ranks) else 30)
            tops = [np.nanmax(a) for a in (r, s, g) if a is not None and a.size
                    and np.isfinite(np.nanmax(a))]
            ax.set_ylim(0, (max(tops) if tops else 1) * 1.15)
            ax.grid(alpha=0.3); ax.set_axisbelow(True)
            for xb in (10.5, 20.5):                     # rank-group dividers
                ax.axvline(xb, color="0.7", ls=(0, (1, 2)), lw=0.8, zorder=1)
            if ri == 0:                                 # group labels along the bottom
                for (lo, hi), gc in zip(_OPT_RANK_GROUPS, grp_colors):
                    ax.annotate(f"{lo}–{hi}", xy=((lo + hi) / 2, 0.02),
                                xycoords=("data", "axes fraction"), ha="center",
                                va="bottom", fontsize=7.5, color=gc, fontweight="bold")
            if stats is not None and stats.get("a") is not None:
                _annotate_15c_panel(ax, ranks, r, s, g, stats["a"], tool, col)
    handles = [
        Line2D([0], [0], color=_OPT_C_RAW, lw=2.0, marker="o", ms=3, label="raw (un-optimized)"),
        Line2D([0], [0], color=_OPT_C_SMI, lw=2.2, marker="o", ms=3, label="smina-opt"),
        Line2D([0], [0], color=_OPT_C_GNI, lw=1.6, ls="--", label="gnina-opt"),
        Patch(facecolor=_OPT_C_SMI, alpha=0.15, label="optimization benefit (raw → smina)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, frameon=False,
               bbox_to_anchor=(0.5, -0.01))
    fig.suptitle("Per-rank benefit of pose optimization (smina / gnina) vs the raw pose  "
                 "· rank groups 1–10 / 11–20 / 21–30 shaded",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_optimization_benefit_by_group(rec: pd.DataFrame, out: Path,
                                       stats: dict | None = None) -> None:
    """Fig 15d — optimization impact (+pp = optimized − raw, pose-count-weighted within
    each group) summarised over three rank groups (1–10, 11–20, 21–30). Rows = DiffDock /
    EquiBind, cols = the three criteria; grouped bars for smina and gnina, value-labelled
    in percentage points. Companion to fig 15c's per-rank curves."""
    if rec is None or rec.empty:
        return
    tools = _opt_tools_in(rec)
    if not tools:
        return
    labels = [f"{lo}–{hi}" for lo, hi in _OPT_RANK_GROUPS]
    x = np.arange(len(_OPT_RANK_GROUPS)); w = 0.38
    fig, axes = plt.subplots(len(tools), 3, figsize=(14.5, 4.5 * len(tools)),
                             squeeze=False, sharex=True)
    for ri, tool in enumerate(tools):
        sub = rec[rec["tool"] == tool]
        raw = sub[sub["optimizer"] == "raw"]
        smi = sub[sub["optimizer"] == "smina"]
        gni = sub[sub["optimizer"] == "gnina"]
        for ci, (col, title) in enumerate(_OPT_BENEFIT_CRITERIA):
            ax = axes[ri][ci]
            base = [_opt_pooled_rate(raw, col, lo, hi) for lo, hi in _OPT_RANK_GROUPS]
            sben = [(_opt_pooled_rate(smi, col, lo, hi) - b) if not smi.empty else np.nan
                    for (lo, hi), b in zip(_OPT_RANK_GROUPS, base)]
            gben = [(_opt_pooled_rate(gni, col, lo, hi) - b) if not gni.empty else np.nan
                    for (lo, hi), b in zip(_OPT_RANK_GROUPS, base)]
            bars = [(ax.bar(x - w / 2, sben, w, color=_OPT_C_SMI, label="smina-opt"), sben),
                    (ax.bar(x + w / 2, gben, w, color=_OPT_C_GNI, label="gnina-opt"), gben)]
            for rects, vals in bars:
                for rect, v in zip(rects, vals):
                    if np.isfinite(v):
                        ax.annotate(f"{v:+.1f}",
                                    xy=(rect.get_x() + rect.get_width() / 2, v),
                                    xytext=(0, 2 if v >= 0 else -2),
                                    textcoords="offset points", ha="center",
                                    va="bottom" if v >= 0 else "top", fontsize=7.5)
            ax.axhline(0, color="0.4", lw=0.8)
            ax.set_title(f"{tool} — {title}", fontsize=10, fontweight="bold")
            ax.set_xticks(x); ax.set_xticklabels(labels)
            if ri == len(tools) - 1:
                ax.set_xlabel("rank group")
            ax.set_ylabel("optimization impact (+pp vs raw)")
            ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
            ax.margins(y=0.20)
            if stats is not None and stats.get("b") is not None:
                _annotate_15d_panel(ax, stats["b"], tool, col)
    axes[0][0].legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    fig.suptitle("Optimization impact by rank group (+percentage points vs the raw pose)",
                 fontsize=13, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# Pocket localization — do the ranked poses target the validated pocket?
# ───────────────────────────────────────────────────────────────────


def _dedup_ranked_poses(sub: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate files of the same (complex, rank) to a single row.

    DiffDock writes its top pose twice — ``rank1.sdf`` and
    ``rank1_confidence-*.sdf`` are the same coordinates — so without this the
    rank-1 pose is counted twice (n=606 instead of 303 for 303 complexes).
    Keeps the first row per (protein, ligand, rank); since the duplicates are
    identical the kept metrics are unchanged.
    """
    return sub.sort_values("rank").drop_duplicates(
        subset=["protein", "ligand", "rank"], keep="first")


def aggregate_pocket_localization(df: pd.DataFrame, top_n: int,
                                  pocket_cutoff: float) -> pd.DataFrame:
    """Per ranking tool: how well do the top-N RANKED poses target the crystal
    (experimentally validated) pocket, and how often does the ORACLE (best-RMSD)
    pose sit in a different pocket than that ranked set?

    "In the validated pocket" = the pose centroid is within ``pocket_cutoff`` Å of
    the crystal ligand centroid (``centroid_dist`` ≤ cutoff). The crystal ligand
    defines the validated pocket; a pose outside the cutoff is taken to be in a
    different (decoy) pocket. Only RANKING_TOOLS are considered (EquiBind has no
    pose ranking). Columns:

      pct_topN_in_validated_pocket   — POSE level: % of all rank ≤ N poses whose
                                       centroid is in the validated pocket (how
                                       reliably the ranking aims at the right site)
      pct_oracle_diff_pocket_vs_topN — COMPLEX level: % of complexes whose oracle
                                       pose is in the validated pocket while NONE
                                       of the top-N ranked poses are — i.e. a
                                       better pose existed in the right pocket but
                                       the ranking sent its top-N elsewhere
      (supporting) pct_top1_in_validated_pocket / pct_any_topN_in_validated_pocket
                 / pct_oracle_in_validated_pocket — complex-level context numbers.
    """
    nan = float("nan")
    rows = []
    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        sub = sub.copy()
        # ``centroid_dist`` is measured at the pose's reference copy j* (the
        # instance by default; the nearest deposited copy under
        # --reference-convention nearest, decision D2). An any-copy variant would
        # switch here to a per-copy minimum instead of the j* value.
        cd = pd.to_numeric(sub["centroid_dist"], errors="coerce")
        sub["on_pocket"] = (cd <= pocket_cutoff).fillna(False)
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        sub = _dedup_ranked_poses(sub)
        topN = sub[sub["rank"] <= top_n]
        n_topN_poses = len(topN)
        pct_topN = 100 * float(topN["on_pocket"].mean()) if n_topN_poses else nan

        per_complex = []
        for _key, g in sub.groupby(["protein", "ligand"]):
            gtop = g[g["rank"] <= top_n]
            gr = g.dropna(subset=["rmsd"])
            if gtop.empty or gr.empty:
                continue
            oracle_on = bool(gr.loc[gr["rmsd"].idxmin(), "on_pocket"])
            top1_on = bool(gtop.sort_values("rank").iloc[0]["on_pocket"])
            any_on = bool(gtop["on_pocket"].any())
            per_complex.append((top1_on, any_on, oracle_on,
                                oracle_on and not any_on))
        pc = pd.DataFrame(per_complex,
                          columns=["top1_on", "any_on", "oracle_on", "oracle_diff"])
        nc = len(pc)
        rows.append({
            "method": method,
            "n_complexes": nc,
            "n_topN_poses": n_topN_poses,
            "top_n": top_n,
            "pocket_cutoff_A": pocket_cutoff,
            "pct_topN_in_validated_pocket": round(pct_topN, PCT_DECIMALS),
            "pct_top1_in_validated_pocket":
                round(100 * float(pc["top1_on"].mean()), PCT_DECIMALS) if nc else nan,
            "pct_any_topN_in_validated_pocket":
                round(100 * float(pc["any_on"].mean()), PCT_DECIMALS) if nc else nan,
            "pct_oracle_in_validated_pocket":
                round(100 * float(pc["oracle_on"].mean()), PCT_DECIMALS) if nc else nan,
            "pct_oracle_diff_pocket_vs_topN":
                round(100 * float(pc["oracle_diff"].mean()), PCT_DECIMALS) if nc else nan,
        })
    return pd.DataFrame(rows).set_index("method") if rows else pd.DataFrame()


def aggregate_pocket_localization_by_rank(df: pd.DataFrame, top_n: int,
                                          pocket_cutoff: float) -> pd.DataFrame:
    """Per ranking tool and per rank position k (1 … top_n): how well that single
    ranked pose targets the crystal (validated) pocket. Tidy/long columns:

      method, rank (k), n_poses (complexes with a rank-k pose),
      pct_in_pocket      — % of the rank-k poses whose centroid is ≤ cutoff Å
                           from the crystal ligand (the k-th pose on its own),
      pct_any_in_top_k   — % of complexes with ANY of ranks 1..k in the pocket
                           (cumulative: the payoff of keeping the first k poses).

    "In the validated pocket" matches ``aggregate_pocket_localization``
    (centroid_dist ≤ ``pocket_cutoff``; NaN distance counts as out-of-pocket).
    """
    nan = float("nan")
    rows = []
    for method, sub in df[df["method"].isin(RANKING_TOOLS)].groupby("method"):
        sub = sub.copy()
        cd = pd.to_numeric(sub["centroid_dist"], errors="coerce")
        sub["on_pocket"] = (cd <= pocket_cutoff).fillna(False)
        sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
        sub = _dedup_ranked_poses(sub)
        for k in range(1, top_n + 1):
            at_k = sub[sub["rank"] == k]
            n_k = len(at_k)
            top_k = sub[sub["rank"] <= k]
            any_by_complex = (top_k.groupby(["protein", "ligand"])["on_pocket"].any()
                              if len(top_k) else pd.Series(dtype=bool))
            rows.append({
                "method": method,
                "rank": k,
                "n_poses": n_k,
                "pct_in_pocket": round(100 * float(at_k["on_pocket"].mean()), PCT_DECIMALS)
                                 if n_k else nan,
                "pct_any_in_top_k": round(100 * float(any_by_complex.mean()), PCT_DECIMALS)
                                    if len(any_by_complex) else nan,
            })
    return pd.DataFrame(rows)


def plot_pocket_targeting_by_rank(by_rank: pd.DataFrame, top_n: int,
                                  pocket_cutoff: float, out: Path) -> None:
    """Per-rank pocket targeting for the ranking tools (fig 16a).

    For each rank position k = 1..N: the % of those rank-k poses that land in the
    validated pocket (solid), and the cumulative % of complexes with ANY of ranks
    1..k in the pocket (dashed). Shows how pocket targeting decays down the ranked
    list and how much keeping more poses buys. In-pocket = pose centroid ≤
    *pocket_cutoff* Å from the crystal ligand centroid.
    """
    if by_rank is None or by_rank.empty:
        return
    methods = list(pd.unique(by_rank["method"]))
    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    for m in methods:
        sub = by_rank[by_rank["method"] == m].sort_values("rank")
        if sub.empty:
            continue
        c = TOOL_COLORS.get(m, "grey")
        ax.plot(sub["rank"], sub["pct_in_pocket"], marker="o", lw=2, color=c,
                markeredgecolor="black", markeredgewidth=0.5,
                label=f"{TOOL_LABEL.get(m, m)} — rank-k pose")
        ax.plot(sub["rank"], sub["pct_any_in_top_k"], marker="s", lw=1.6, ls="--",
                color=c, alpha=0.85, markeredgecolor="black", markeredgewidth=0.4,
                label=f"{TOOL_LABEL.get(m, m)} — any of top-1..k")
    ax.set_xticks(range(1, top_n + 1))
    ax.set_xlabel("Rank position k  (1 = top-ranked pose)")
    ax.set_ylabel("Poses in the validated pocket (%)")
    ax.set_ylim(0, 105)
    ax.grid(alpha=0.3); ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="lower left")
    ax.set_title(_vt("Per-rank pocket targeting (ranking tools)\n"
                     "solid = the k-th pose alone · dashed = any of the first k\n"
                     f"(in-pocket ≤ {pocket_cutoff:g} Å from crystal ligand centroid)"),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def plot_pocket_localization_summary(summ: pd.DataFrame, top_n: int,
                                     pocket_cutoff: float, out: Path) -> None:
    """Per-tool pocket-localization summary dumbbell (fig 16b).

    Per ranking tool: % of all top-N ranked poses in the validated pocket (filled)
    vs % of complexes whose oracle pose sits in a DIFFERENT pocket than the top-N
    ranked set (open — a ranking failure: a better pose existed in the right pocket
    but the ranking sent its top-N elsewhere). In-pocket = pose centroid ≤
    *pocket_cutoff* Å from the crystal ligand centroid.
    """
    if summ is None or summ.empty:
        return
    methods = list(summ.index)
    colors = [TOOL_COLORS.get(m, "grey") for m in methods]
    v1 = [float(summ.loc[m, "pct_topN_in_validated_pocket"]) for m in methods]
    v2 = [float(summ.loc[m, "pct_oracle_diff_pocket_vs_topN"]) for m in methods]
    labels = [f"{TOOL_LABEL.get(m, m)}\n(n={int(summ.loc[m, 'n_complexes'])} complexes · "
              f"{int(summ.loc[m, 'n_topN_poses'])} top-{top_n} poses)" for m in methods]
    fig, ax = plt.subplots(figsize=(10, max(3.2, 1.25 * len(methods) + 2)))
    handles = _dumbbell(ax, labels, v1, v2, colors,
                        left_name=f"All top-{top_n} ranked poses in the validated "
                                  f"pocket (% of poses)",
                        right_name=f"Oracle pose in a different pocket than the "
                                   f"top-{top_n} ranked (% of complexes)",
                        value_fmt="{:.1f}%")
    ax.set_xlabel("Percentage")
    ax.legend(handles=handles, fontsize=8, loc="lower right")
    ax.set_title(_vt("Pocket localization summary — does the ranking target the "
                     "validated pocket?\n"
                     f"(in-pocket = pose centroid ≤ {pocket_cutoff:g} Å from the crystal "
                     f"ligand centroid)"), fontsize=12, fontweight="bold")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


# ───────────────────────────────────────────────────────────────────
# PoseBusters test-failure waterfall (paper Fig.-style cascade)
# ───────────────────────────────────────────────────────────────────

# Readable labels for the canonical PoseBusters test columns, phrased as the
# FAILURE mode (what removes a pose at that step) to mirror the paper's waterfall.
PB_TEST_LABELS: dict[str, str] = {
    "mol_pred_loaded": "Input cannot be loaded",
    "sanitization": "Sanitisation fails",
    "inchi_convertible": "InChI not convertible",
    "all_atoms_connected": "Bonds not preserved",
    "no_radicals": "Radicals present",
    "bond_lengths": "Bond lengths out of bounds",
    "bond_angles": "Bond angles out of bounds",
    "internal_steric_clash": "Internal steric clash",
    "aromatic_ring_flatness": "Deformed aromatic rings",
    "non-aromatic_ring_non-flatness": "Flat non-aromatic ring",
    "double_bond_flatness": "Deformed double bonds",
    "internal_energy": "Energy too high",
    "minimum_distance_to_protein": "Protein-ligand distance too small",
    "minimum_distance_to_organic_cofactors": "Distance to organic cofactors too small",
    "minimum_distance_to_inorganic_cofactors": "Distance to inorganic cofactors too small",
    "minimum_distance_to_waters": "Distance to waters too small",
    "volume_overlap_with_protein": "Volume overlap with protein",
    "volume_overlap_with_organic_cofactors": "Volume overlap with organic cofactors",
    "volume_overlap_with_inorganic_cofactors": "Volume overlap with inorganic cofactors",
    "volume_overlap_with_waters": "Volume overlap with waters",
}


def _load_pb_test_table(pb_csv: Path, split_equibind: bool = True,
                        diffdock_variant: str = "diffdock",
                        select_diffdock: bool = True):
    """Read the per-pose PoseBusters PASS/FAIL columns from the raw CSV.

    ``_build_pose_index`` drops the individual test columns (it keeps only the
    AND-ed ``pb_valid``), so the waterfall reads them back here. Returns
    ``(table, test_cols)`` where *table* is keyed by
    (method, protein, ligand, pose_name) — the same keys as per_pose_metrics —
    with one bool column per canonical test (True = pass). ``None`` if the CSV
    lacks the required identity columns.
    """
    raw = pd.read_csv(pb_csv, low_memory=False)
    key_cols = {"docking_method", "protein", "ligand", "pose_name"}
    if not key_cols.issubset(raw.columns):
        return None
    test_cols = [c for c in PB_CRITICAL_CHECKS if c in raw.columns]
    if not test_cols:
        return None
    raw = raw.copy()
    raw["docking_method"] = raw["docking_method"].astype(str).str.lower()
    if split_equibind:
        raw = _apply_equibind_split(raw)
    # Match _build_pose_index's AutoDock variant labels exactly; otherwise raw
    # and gnina rows share a method during waterfall joins/summaries.
    raw = _apply_autodock_split(raw)
    # Split + select the DiffDock variant identically to _build_pose_index, so
    # this table's method keys match per_pose_metrics (otherwise the selected
    # variant's rows fail to join in the waterfall).
    raw = _apply_diffdock_split(raw)
    if select_diffdock:
        raw = _select_diffdock_variant(raw, diffdock_variant)
    table = pd.DataFrame({
        "method": raw["docking_method"],
        "protein": raw["protein"],
        "ligand": raw["ligand"],
        "pose_name": raw["pose_name"],
    })
    for c in test_cols:
        table[c] = _to_bool(raw[c])
    table = table.drop_duplicates(subset=["method", "protein", "ligand", "pose_name"])
    return table, test_cols


def _pb_cascade_from_rep(rep: pd.DataFrame, test_table: pd.DataFrame,
                         test_cols: list[str], selection: str,
                         gate_rmsd: bool = True) -> dict | None:
    """Build one panel's sequential PoseBusters filter cascade from a
    representative-pose frame (one row per complex; must carry
    method/protein/ligand/pose_name/rmsd).

    Starts from every representative pose and, when ``gate_rmsd`` is True (the
    default), first drops those with RMSD > 2 Å; then applies each canonical
    PoseBusters test IN ORDER, removing at each step the poses that fail it AND
    passed every previous filter. With ``gate_rmsd=False`` the RMSD step is
    skipped entirely, so every pose faces the physical-validity tests regardless
    of placement — the cascade then shows where poses fail on PoseBusters grounds
    alone. Returns ``{N, selection, steps}`` (each step:
    label/kind/removed/remaining/before/after, kind ∈ {start, drop, milestone,
    end}) or ``None`` when the frame is empty / nothing joins the test table.
    Single source of truth for the cascade shared by the top-1/top-N and the
    raw-vs-best waterfalls.
    """
    rep = rep[["method", "protein", "ligand", "pose_name", "rmsd"]]
    if rep.empty:
        return None
    # Join on pose identity (protein, ligand, pose_name) — NOT method. A
    # representative frame may carry a collapsed method label (e.g. a best DiffDock
    # variant relabelled ``diffdock_smina`` → ``diffdock`` by _select_best_diffdock)
    # that no longer matches the test table's per-variant key, which would silently
    # miss every join and, via the fillna(True) below, count those poses as passing
    # every PoseBusters test. pose_name is globally unique per (protein, ligand)
    # across methods, so dropping method from the key is safe and correct.
    tt = test_table.drop(columns=["method"], errors="ignore")
    merged = rep.merge(tt, how="left", on=["protein", "ligand", "pose_name"])
    merged = merged.drop_duplicates(subset=["protein", "ligand", "pose_name"])
    N = len(merged)
    if N == 0:
        return None

    steps = [{"label": "All predictions", "kind": "start",
              "removed": 0, "remaining": N, "before": N, "after": N}]
    surv = pd.Series(True, index=merged.index)
    placement = pd.Series(True, index=merged.index)   # RMSD ≤ 2 Å status per complex

    if gate_rmsd:
        rmsd_ok = (pd.to_numeric(merged["rmsd"], errors="coerce") <= 2.0).fillna(False)
        placement = rmsd_ok
        before = int(surv.sum()); surv = surv & rmsd_ok; after = int(surv.sum())
        steps.append({"label": "RMSD > 2 Å", "kind": "drop",
                      "removed": before - after, "remaining": after,
                      "before": before, "after": after})
        steps.append({"label": "RMSD ≤ 2 Å", "kind": "milestone",
                      "removed": 0, "remaining": after, "before": after, "after": after})

    for col in test_cols:
        passed = (merged[col].fillna(True).astype(bool) if col in merged
                  else pd.Series(True, index=merged.index))
        before = int(surv.sum()); surv = surv & passed; after = int(surv.sum())
        steps.append({"label": PB_TEST_LABELS.get(col, col), "kind": "drop",
                      "removed": before - after, "remaining": after,
                      "before": before, "after": after})

    passing = int(surv.sum())
    steps.append({"label": "Passing all tests", "kind": "end",
                  "removed": 0, "remaining": passing,
                  "before": passing, "after": passing})
    # Per-complex outcomes keyed by "protein\x1fligand" so downstream paired tests
    # (McNemar: rank-1 vs top-N, raw vs refined) align the SAME complexes across
    # cascades. ``surv`` now holds passing-all-tests status; ``placement`` the RMSD
    # gate. One representative pose per complex, so each key is unique.
    ck = merged["protein"].astype(str) + "\x1f" + merged["ligand"].astype(str)
    pass_all = dict(zip(ck, surv.reindex(merged.index).to_numpy(dtype=bool)))
    rmsd_ok_map = dict(zip(ck, placement.reindex(merged.index).to_numpy(dtype=bool)))
    return {"N": N, "selection": selection, "steps": steps,
            "pass_all": pass_all, "rmsd_ok": rmsd_ok_map}


def aggregate_pb_waterfall(df: pd.DataFrame, test_table: pd.DataFrame,
                           test_cols: list[str], top_n: int,
                           rank_selection: str = "top1") -> dict:
    """Build, per method, the sequential PoseBusters filter cascade.

    Mirrors the PoseBusters paper's waterfall: start from every representative
    pose (see ``rank_selection`` for ranking tools; best/oracle pose for unranked
    tools such as EquiBind), drop those with RMSD > 2 Å, then apply each canonical
    PoseBusters test IN ORDER, removing at each step the poses that fail it AND
    passed every previous filter. The remainder at the end pass everything.

    ``rank_selection`` controls the representative pose for RANKING_TOOLS:
    ``"top1"`` (default) uses the rank-1 pose; ``"topn"`` uses the lowest-RMSD
    pose among the top-``top_n`` ranks (best-of-top-N — one pose per complex, so
    it shares the top-1 denominator). Unranked tools always use the oracle pose.

    EquiBind variants are collapsed to the single best one (``_select_best_equibind``)
    so the figure stays to a few panels. Returns ``{method: {N, selection, steps}}``
    where each step is a dict (label, kind, removed, remaining, before, after);
    *kind* ∈ {start, drop, milestone, end} and before/after are running counts.
    """
    wdf, _ = _select_best_equibind(df)
    if rank_selection == "topn":
        ranked_rep = _best_topn_per_pair(wdf, top_n)
        ranked_sel = f"best of top-{top_n} ranked poses"
    else:
        ranked_rep = _top1_per_pair(wdf)
        ranked_sel = "top-1 ranked pose"
    oracle = _oracle_per_pair(wdf)
    present = set(wdf["method"].astype(str))
    methods = [m for m in ("autodock", "diffdock") if m in present]
    methods += [m for m in sorted(present) if m not in RANKING_TOOLS]

    cascades: dict = {}
    for method in methods:
        if method in RANKING_TOOLS:
            rep, selection = ranked_rep[ranked_rep["method"] == method], ranked_sel
        else:
            rep, selection = oracle[oracle["method"] == method], \
                "best (oracle) pose — no ranking"
        cascade = _pb_cascade_from_rep(rep, test_table, test_cols, selection)
        if cascade is None:
            continue
        cascades[method] = cascade
    return cascades


def _waterfall_to_csv(cascades: dict, path: Path) -> None:
    rows = []
    for method, c in cascades.items():
        N = c["N"]
        for order, s in enumerate(c["steps"]):
            rows.append({
                "method": method, "selection": c["selection"], "order": order,
                "step": s["label"], "kind": s["kind"],
                "poses_removed": s["removed"], "poses_remaining": s["remaining"],
                "remaining_pct": round(100 * s["remaining"] / N, PCT_DECIMALS) if N else float("nan"),
            })
    _write_csv(pd.DataFrame(rows), path, index=False)


def plot_pb_waterfall(cascades: dict, out: Path,
                      sel_desc: str = "top-ranked poses per method") -> None:
    """Small-multiples PoseBusters waterfall — one panel per method, stacked so
    the long test labels are shared. Each red bar is the poses removed by that
    test (annotated −k); the teal milestone bar is the RMSD ≤ 2 Å running total
    and the green bar the poses passing every test.

    ``sel_desc`` names the representative pose in the suptitle (e.g.
    "best of top-15 ranked poses per method" for the top-N variant).
    """
    methods = list(cascades.keys())
    if not methods:
        return
    labels = [s["label"] for s in cascades[methods[0]]["steps"]]
    nx = len(labels)
    nrows = len(methods)
    fig, axes = plt.subplots(nrows, 1, sharex=True, squeeze=False,
                             figsize=(max(11, 0.62 * nx), 3.0 * nrows + 2.0))
    axes = list(axes[:, 0])
    teal, red, green = "#3a9d8f", "#d1495b", "#2a9d3f"

    for ax, method in zip(axes, methods):
        c = cascades[method]; N = c["N"]; steps = c["steps"]
        x = np.arange(nx)
        for xi, s in zip(x, steps):
            b = 100 * s["before"] / N if N else 0.0
            a = 100 * s["after"] / N if N else 0.0
            if s["kind"] == "start":
                ax.bar(xi, 100, color=teal, alpha=0.30, edgecolor=teal,
                       hatch="..", zorder=2)
                ax.text(xi, 101, f"{N}", ha="center", va="bottom",
                        fontsize=8, fontweight="bold")
            elif s["kind"] == "milestone":
                ax.bar(xi, a, color=teal, alpha=0.45, edgecolor=teal,
                       hatch="//", zorder=2)
                ax.text(xi, a + 1, f"{s['remaining']}", ha="center",
                        va="bottom", fontsize=8)
            elif s["kind"] == "end":
                ax.bar(xi, a, color=green, alpha=0.85, edgecolor="black", zorder=3)
                ax.text(xi, a + 1, f"{s['remaining']}\n({a:.1f}%)", ha="center",
                        va="bottom", fontsize=8, fontweight="bold")
            else:  # drop
                if s["removed"] > 0:
                    ax.bar(xi, b - a, bottom=a, color=red, alpha=0.9,
                           edgecolor="black", zorder=3)
                    ax.text(xi, b + 0.5, f"−{s['removed']}", ha="center",
                            va="bottom", fontsize=7, color=red, fontweight="bold")
                else:
                    ax.text(xi, a + 0.5, "0", ha="center", va="bottom",
                            fontsize=7, color="0.5")
        run = [100 * s["after"] / N if N else 0.0 for s in steps]
        ax.step(x, run, where="mid", color="0.4", lw=0.8, ls="--",
                alpha=0.7, zorder=1)
        ax.set_ylim(0, 110)
        ax.set_ylabel("% of poses\nremaining", fontsize=9)
        final_pct = 100 * steps[-1]["remaining"] / N if N else float("nan")
        ax.set_title(f"{TOOL_LABEL.get(method, method)} — {c['selection']}  "
                     f"(n={N}; passing all tests = {final_pct:.1f}%)",
                     fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)

    axes[-1].set_xticks(np.arange(nx))
    axes[-1].set_xticklabels(labels, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    fig.suptitle(_vt(f"PoseBusters test-failure waterfall — {sel_desc}\n"
                     "(sequential filter: each red bar = poses removed by "
                     "that test; green = poses passing every test)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


def _waterfall_survival(cascade: dict) -> tuple[list[str], list[float]]:
    """(step labels, running % remaining) for one method's cascade, each step's
    survivors as a % of that cascade's own start count N."""
    N = cascade["N"]
    labels = [s["label"] for s in cascade["steps"]]
    surv = [100 * s["after"] / N if N else float("nan") for s in cascade["steps"]]
    return labels, surv


def _rank_word(method: str) -> str:
    """Ranking-basis word for a panel: EquiBind variants have no native rank and
    are ordered by gnina affinity; every other method uses its native rank."""
    return "gnina rank" if str(method).startswith("equibind") else "rank"


def _panel_label(method: str, labels: dict | None = None) -> str:
    """Display label for a ranking-headroom / recovery panel, honouring an explicit
    per-figure override so a raw-variant figure can read 'DiffDock (raw)' where the
    refined figure reads 'DiffDock*'."""
    if labels and method in labels:
        return labels[method]
    if method == "equibind_gnina":
        return "EquiBind*"
    if method == "equibind_raw":
        return "EquiBind* (raw)"
    if method == "diffdock_raw":
        return "DiffDock (raw)"
    return TOOL_LABEL.get(method, method)


def plot_pb_waterfall_top1_vs_topn(cascades_top1: dict, cascades_topn: dict,
                                   top_n: int, out: Path,
                                   csv_out: Path | None = None,
                                   extra: tuple = (),
                                   methods: tuple | None = None,
                                   labels: dict | None = None) -> None:
    """Compare the PB cascade survival for top-1 vs best-of-top-N poses.

    One panel per ranking tool (AutoDock, DiffDock — native rank). Each panel
    overlays two survival curves — the % of complexes whose representative pose
    still survives after each sequential filter — for the rank-1 pose (solid grey)
    and the best-of-top-N pose (dashed green). The shaded band between them is the
    recovery from letting the ranker offer more poses; Δ callouts mark the gain at
    "RMSD ≤ 2 Å" and at "Passing all tests".

    ``extra`` lists additional method keys to append after AutoDock/DiffDock — used
    to add a gnina-ranked EquiBind panel, which stands in for EquiBind's missing
    native ranking by ordering poses on gnina affinity. ``methods`` overrides the
    panel set/order entirely (for the raw-variant figure); ``labels`` overrides the
    per-method display label.
    """
    wanted = methods if methods is not None else ("autodock", "diffdock", *extra)
    methods = [m for m in wanted
               if m in cascades_top1 and m in cascades_topn]
    if not methods:
        return
    grey, green = "#555555", "#2a9d3f"
    labels, _ = _waterfall_survival(cascades_top1[methods[0]])
    nx = len(labels)
    nrows = len(methods)
    fig, axes = plt.subplots(nrows, 1, sharex=True, squeeze=False,
                             figsize=(max(11, 0.62 * nx), 3.2 * nrows + 1.8))
    axes = list(axes[:, 0])
    x = np.arange(nx)
    csv_rows = []

    def _idx(labs: list[str], name: str) -> int | None:
        return labs.index(name) if name in labs else None

    for ax, method in zip(axes, methods):
        c1, cn = cascades_top1[method], cascades_topn[method]
        is_eq = str(method).startswith("equibind")
        rank_word = _rank_word(method)
        method_label = _panel_label(method, labels)
        _, s1 = _waterfall_survival(c1)
        ln, sn = _waterfall_survival(cn)
        ax.fill_between(x, s1, sn, step="mid", color=green, alpha=0.12, zorder=1)
        ax.step(x, s1, where="mid", color=grey, lw=1.8,
                label=f"{rank_word}-1 pose", zorder=3)
        ax.step(x, sn, where="mid", color=green, lw=1.8, ls="--",
                label=f"best of top-{top_n}" + (" (gnina)" if is_eq else ""),
                zorder=3)
        ax.plot(x, s1, "o", color=grey, ms=3.5, zorder=4)
        ax.plot(x, sn, "D", color=green, ms=3.5, zorder=4)

        for name in ("RMSD ≤ 2 Å", "Passing all tests"):
            i = _idx(labels, name)
            if i is None:
                continue
            d = sn[i] - s1[i]
            ax.annotate(f"+{d:.1f} pp", xy=(i, (s1[i] + sn[i]) / 2),
                        xytext=(4, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=8,
                        color=green, fontweight="bold")
            csv_rows.append({"method": method, "step": name,
                             "top1_remaining_%": round(s1[i], PCT_DECIMALS),
                             f"top{top_n}_remaining_%": round(sn[i], PCT_DECIMALS),
                             "gain_pp": round(d, PCT_DECIMALS)})

        ax.set_ylim(0, 110)
        ax.set_ylabel("% of complexes\nsurviving", fontsize=9)
        end1, endn = s1[-1], sn[-1]
        ax.set_title(f"{method_label} — passing all tests: "
                     f"{end1:.1f}% ({rank_word}-1) → {endn:.1f}% (top-{top_n}), "
                     f"+{endn - end1:.1f} pp",
                     fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    eq_note = ("  ·  EquiBind ranked by gnina affinity (no native rank)"
               if any(str(m).startswith("equibind") for m in methods) else "")
    fig.suptitle(_vt(f"Ranking headroom — rank-1 vs best-of-top-{top_n} pose "
                     f"through the PoseBusters cascade{eq_note}\n(green band = "
                     "complexes recovered by considering more ranked poses)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        _write_csv(pd.DataFrame(csv_rows), csv_out, index=False)


def _step_category(step: dict) -> tuple[str, str]:
    """(category, colour) for a cascade step in the top-30 recovery figure:
    the RMSD ≤ 2 Å placement filter, the physical-validity PoseBusters tests, or
    the final passing-all bar."""
    blue, amber, green, grey = "#3b6fb0", "#e08a2e", "#2a9d3f", "#9aa0a6"
    if step["label"] in ("RMSD > 2 Å", "RMSD ≤ 2 Å"):
        return "placement (RMSD ≤ 2 Å)", blue
    if step["kind"] == "end":
        return "passing all tests", green
    if step["kind"] == "start":
        return "start", grey
    return "other PoseBusters tests", amber


def _mcnemar_on_maps(map_a: dict | None, map_b: dict | None) -> dict | None:
    """Exact paired McNemar on two per-complex pass/fail maps (keyed by complex
    id), restricted to complexes present in BOTH. ``a`` is the baseline condition,
    ``b`` the comparison. Returns each rate + its Wilson 95% CI and the exact
    McNemar p, or None when the maps share no complex."""
    if not map_a or not map_b:
        return None
    keys = sorted(set(map_a) & set(map_b))
    n = len(keys)
    if n == 0:
        return None
    a = np.fromiter((bool(map_a[k]) for k in keys), dtype=bool, count=n)
    b = np.fromiter((bool(map_b[k]) for k in keys), dtype=bool, count=n)
    ka, kb = int(a.sum()), int(b.sum())
    n10, n01, p = su.mcnemar_exact(a, b)
    return {"n": n, "a_k": ka, "b_k": kb, "a_rate": ka / n, "b_rate": kb / n,
            "a_ci": su.wilson_ci(ka, n), "b_ci": su.wilson_ci(kb, n),
            "a_wins": int(n10), "b_wins": int(n01), "mcnemar_p": float(p)}


def _stats_top30_recovery(cascades_top1: dict, cascades_topn: dict, methods,
                          raw_key_map: dict | None = None) -> dict | None:
    """Paired significance tests for the top-30 recovery figures (17f/17g/17h).

    Per method panel, on the per-complex 'passing all PoseBusters tests' outcome
    (one representative pose per complex, so no pose-level pseudoreplication):
      * ``rank_gain``  — rank-1 vs best-of-top-N (the panel's coloured recovery
        bars): exact McNemar + Wilson 95% CIs on each pass-all rate. Two-sided —
        the top-N pick minimises RMSD, not PB-validity, so a complex can pass at
        rank-1 yet fail the top-N pick.
      * ``refine_gain`` — raw vs refined at BOTH rank-1 and top-N depths, but only
        for a refined panel whose raw counterpart is present in the same figure.
    Every p-value in the figure forms ONE family and is Holm-adjusted together.
    Returns ``{method: {...}}`` (each record carrying ``p_holm``/``star``) or None
    when no panel carries the per-complex maps.
    """
    raw_key_map = raw_key_map or {"diffdock": "diffdock_raw",
                                  "equibind_gnina": "equibind_raw"}
    out: dict = {}
    fam: list = []                 # records to back-fill Holm-adjusted p-values
    for m in methods:
        c1, cn = cascades_top1.get(m), cascades_topn.get(m)
        if not c1 or not cn:
            continue
        rank = _mcnemar_on_maps(c1.get("pass_all"), cn.get("pass_all"))
        if rank is None:
            continue
        entry: dict = {"rank_gain": rank}
        fam.append(rank)
        rk = raw_key_map.get(str(m))
        if rk and rk in cascades_top1 and rk in cascades_topn:
            refine: dict = {}
            for depth, casc in (("rank1", cascades_top1), ("topn", cascades_topn)):
                rec = _mcnemar_on_maps(casc[rk].get("pass_all"),
                                       casc[m].get("pass_all"))
                if rec is not None:
                    refine[depth] = rec
                    fam.append(rec)
            if refine:
                entry["refine_gain"] = refine
        out[str(m)] = entry
    if not out:
        return None
    for rec, pa in zip(fam, su.holm([r["mcnemar_p"] for r in fam])):
        rec["p_holm"] = float(pa)
        rec["star"] = su.p_stars(pa)
    return out


def plot_pb_top30_recovery_by_test(cascades_top1: dict, cascades_topn: dict,
                                   top_n: int, out: Path,
                                   csv_out: Path | None = None,
                                   extra: tuple = (),
                                   methods: tuple | None = None,
                                   labels: dict | None = None,
                                   stats: dict | None = None,
                                   only_differences: bool = True) -> None:
    """Per method, how many MORE complexes survive each cascade step when the
    representative pose is the best of the top-``top_n`` instead of rank-1.

    A count-space companion to 17c/17c2: at each step of the sequential PoseBusters
    cascade a bar spans the rank-1 survivor count (grey line) up to the
    best-of-top-N survivor count (green line) — the height is the extra complexes
    recovered at that step, labelled ``+k``. Bars are coloured by category so the
    two kinds of test are distinguished: the RMSD ≤ 2 Å placement filter (blue) vs
    each physical-validity PoseBusters test (amber), with the final passing-all bar
    in green. A per-panel summary states how many extra complexes the top-N pool
    places within 2 Å, how many of those extras are then lost to the other tests,
    and the net that pass everything.

    ``only_differences`` (default) draws a bar ONLY where the top-N pool changes the
    running surplus — the RMSD ≤ 2 Å placement bar, the final passing-all bar, and
    the individual PoseBusters tests that actually remove a different number of
    poses from the top-N pick than from rank-1 — so the many no-change carry-over
    steps stay bar-free. ``methods``/``labels`` override the panel set and labels
    (for the raw-variant figure).
    """
    wanted = methods if methods is not None else ("autodock", "diffdock", *extra)
    methods = [m for m in wanted
               if m in cascades_top1 and m in cascades_topn]
    if not methods:
        return
    grey, green = "#555555", "#2a9d3f"
    labels_x = [s["label"] for s in cascades_top1[methods[0]]["steps"]]
    nx = len(labels_x)
    nrows = len(methods)
    fig, axes = plt.subplots(nrows, 1, sharex=True, squeeze=False,
                             figsize=(max(11, 0.62 * nx), 3.3 * nrows + 2.0))
    axes = list(axes[:, 0])
    x = np.arange(nx)
    csv_rows = []

    # Alternating vertical bands at every x-tick, painted on each (shared-x) panel
    # so a cascade-step label on the bottom axis can be traced straight up through
    # all panels above it. Two faint tints alternate column-to-column; the low
    # zorder keeps them behind the grid, the gap bars and the survivor lines.
    band_colors = ("#eceff4", "#d9e1ec")
    for ax in axes:
        for xi in x:
            ax.axvspan(xi - 0.5, xi + 0.5, color=band_colors[int(xi) % 2],
                       alpha=0.7, lw=0, zorder=0)

    def _lbl_idx(name: str) -> int | None:
        return labels_x.index(name) if name in labels_x else None

    for ax, method in zip(axes, methods):
        c1, cn = cascades_top1[method], cascades_topn[method]
        steps1, stepsn = c1["steps"], cn["steps"]
        s1 = [s["after"] for s in steps1]
        sn = [s["after"] for s in stepsn]
        gaps = [sn[k] - s1[k] for k in range(nx)]
        is_eq = str(method).startswith("equibind")
        method_label = _panel_label(method, labels)
        y_top = max(sn) if sn else 1

        def _show_bar(k: int) -> bool:
            lbl, kind = steps1[k]["label"], steps1[k]["kind"]
            if kind == "start" or lbl == "RMSD > 2 Å":
                return False              # start / redundant with the milestone
            if lbl == "RMSD ≤ 2 Å" or kind == "end":
                return True               # placement bar / final passing-all bar
            if not only_differences:
                return True
            return gaps[k] != gaps[k - 1]  # a test that shifts the running surplus

        # Gap bars (rank-1 → top-N survivors), coloured by test category. Every
        # step is recorded to the CSV; only the difference-making steps get a bar.
        for xi in x:
            if xi == 0:
                continue
            gap = gaps[xi]
            cat, col = _step_category(steps1[xi])
            csv_rows.append({
                "method": method, "step": labels_x[xi], "category": cat,
                "rank1_survivors": s1[xi], f"top{top_n}_survivors": sn[xi],
                "extra_survivors": gap, "bar_shown": bool(_show_bar(xi))})
            if not _show_bar(xi):
                continue
            ax.bar(xi, gap, bottom=s1[xi], color=col, alpha=0.85,
                   edgecolor="white", linewidth=0.4, width=0.9, zorder=2)
            ax.text(xi, sn[xi] + y_top * 0.012, f"+{gap}", ha="center",
                    va="bottom", fontsize=7.5, fontweight="bold",
                    color="0.15", zorder=5)
        # Absolute survivor levels for context.
        ax.step(x, s1, where="mid", color=grey, lw=1.6, zorder=3,
                label=(_rank_word(method) + "-1") + " survivors")
        ax.step(x, sn, where="mid", color=green, lw=1.6, ls="--", zorder=3,
                label=f"best of top-{top_n} survivors")

        i_rmsd = _lbl_idx("RMSD ≤ 2 Å")
        i_end = _lbl_idx("Passing all tests")
        place_gain = sn[i_rmsd] - s1[i_rmsd] if i_rmsd is not None else 0
        final_gain = sn[i_end] - s1[i_end] if i_end is not None else 0
        lost = place_gain - final_gain
        p1 = 100 * s1[-1] / (c1["N"] or 1)     # passing-all %, rank-1
        p30 = 100 * sn[-1] / (cn["N"] or 1)    # passing-all %, top-N
        ax.set_ylim(0, y_top * 1.18)
        ax.set_ylabel("complexes\nsurviving", fontsize=9)
        # Title: passing-all rate and the rank-1 → top-N percentage-point gain.
        ax.set_title(
            f"{method_label} — passing all tests: {p1:.1f}% (rank-1) → {p30:.1f}% "
            f"(top-{top_n}), +{p30 - p1:.1f} pp   ·   n={c1['N']}",
            fontsize=9.5, fontweight="bold")
        # In-panel box: the count decomposition, and — on a refined panel whose raw
        # counterpart is also in the figure — the raw → refined percentage-point
        # lift at both rank-1 and top-N.
        note = (f"top-{top_n} recovers:  +{place_gain} within 2 Å  ·  −{lost} to "
                f"other tests  ·  +{final_gain} pass all")
        st = (stats or {}).get(str(method))
        raw_key = {"diffdock": "diffdock_raw",
                   "equibind_gnina": "equibind_raw"}.get(str(method))
        if raw_key and raw_key in cascades_top1 and raw_key in cascades_topn:
            rp1 = 100 * cascades_top1[raw_key]["steps"][-1]["after"] / (
                cascades_top1[raw_key]["N"] or 1)
            rp30 = 100 * cascades_topn[raw_key]["steps"][-1]["after"] / (
                cascades_topn[raw_key]["N"] or 1)
            # Star each raw → refined pp lift with its own paired McNemar result.
            ref = (st or {}).get("refine_gain", {})
            r1s = f" {ref['rank1']['star']}" if ref.get("rank1") else ""
            r30s = f" {ref['topn']['star']}" if ref.get("topn") else ""
            note += (f"\nraw → refined:  +{p1 - rp1:.1f} pp at rank-1{r1s}  ·  "
                     f"+{p30 - rp30:.1f} pp at top-{top_n}{r30s}")
        # Paired significance of the rank-1 → top-N passing-all gain (the panel's
        # coloured bars): Wilson 95% CI on each rate + exact McNemar, Holm-adjusted
        # across the figure. Absent when stats weren't computed (degrades silently).
        rg = st.get("rank_gain") if st else None
        if rg:
            lo1, hi1 = rg["a_ci"]; lon, hin = rg["b_ci"]
            pshow = rg.get("p_holm", rg["mcnemar_p"])
            note += (f"\npass-all % [95% CI]:  {100 * rg['a_rate']:.1f} "
                     f"[{100 * lo1:.0f}–{100 * hi1:.0f}] → {100 * rg['b_rate']:.1f} "
                     f"[{100 * lon:.0f}–{100 * hin:.0f}]   ·   McNemar "
                     f"{rg['star']} ({su.fmt_p(pshow)})")
        ax.text(0.5, 0.93, note, transform=ax.transAxes, ha="center", va="top",
                fontsize=8, color="0.12",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.55",
                          alpha=0.92), zorder=6)
        ax.grid(axis="y", alpha=0.25)

    # One combined legend (survivor lines + bar-colour key) in the empty top-right
    # of the first panel — avoids clobbering and keeps every panel uncluttered.
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], color="#555555", lw=1.6, label="rank-1 survivors"),
        Line2D([0], [0], color="#2a9d3f", lw=1.6, ls="--",
               label=f"best of top-{top_n} survivors"),
        Patch(facecolor="#3b6fb0", label="bar: RMSD ≤ 2 Å (placement)"),
        Patch(facecolor="#e08a2e", label="bar: other PoseBusters tests"),
        Patch(facecolor="#2a9d3f", label="bar: passing all tests"),
    ]
    # Placed low: the bars sit high (rank-1 → top-N survivor band), so the area
    # below the rank-1 line is empty in every panel.
    axes[0].legend(handles=handles, fontsize=7.5, loc="lower right",
                   framealpha=0.92, ncol=1)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels_x, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    # Title broken over four lines (each subtitle clause on its own line) so it
    # stays readable regardless of figure width, and pulled down close to the top
    # panel. NB: passing the title band to tight_layout via rect= does NOT work —
    # tight_layout leaves ~6% slack below the rect top instead of filling it, so the
    # title floats far above the graph. Instead run tight_layout() normally (for the
    # inter-panel + x-label spacing), then MEASURE the rendered title bottom and push
    # the subplot region up to just below it with subplots_adjust(top=…).
    title_y = 0.992
    suptitle = fig.suptitle(_vt(
        f"Extra complexes recovered by best-of-top-{top_n} vs rank-1,\n"
        "at each PoseBusters cascade step\n"
        f"(bar = extra survivors where the top-{top_n} pick makes a difference;\n"
        "blue = RMSD ≤ 2 Å placement, amber = other physical-validity tests, "
        "green = passing all tests)"),
        fontsize=13, fontweight="bold", y=title_y, va="top")
    fig.tight_layout()
    try:
        fig.canvas.draw()
        bb = suptitle.get_window_extent(renderer=fig.canvas.get_renderer())
        title_bottom = float(fig.transFigure.inverted().transform((0, bb.y0))[1])
        fig.subplots_adjust(top=min(0.97, title_bottom - 0.018))
    except Exception:
        fig.subplots_adjust(top=title_y - 4 * (0.253 / fig.get_figheight()) - 0.02)
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        _write_csv(pd.DataFrame(csv_rows), csv_out, index=False)


def _write_pb_recovery_report(out_path: Path, *, figure_names, csv_names, methods,
                              labels: dict | None, cascades_top1: dict,
                              depth_blocks) -> None:
    """ONE combined companion text report for the raw-vs-refined PoseBusters recovery
    figures across MULTIPLE pool depths — showing top-1 (rank-1 baseline),
    best-of-top-15 and best-of-top-30 together.

    ``cascades_top1`` is the rank-1 (top-1) cascade dict (depth-independent, so shown
    once). ``depth_blocks`` is an ordered list of ``{"depth": int,
    "cascades_topn": dict, "stats": dict | None}`` — one per best-of-top-N figure.

    Writes two blocks:
      1. STATISTICAL TESTS — a sub-section per figure/depth (each = exactly that
         figure's own Holm family, matching its on-panel annotations): exact two-sided
         McNemar (rank-1 → top-N pass-all gain, and raw → refined at rank-1 and top-N
         on the refined panels), each with both passing-all rates, Wilson 95 % CIs,
         the discordant-pair counts (b, c), the raw and Holm-adjusted p-values, and
         the significance star.
      2. PER-STAGE ATTRITION TABLE (';'-separated) — one row per
         method × selection × cascade stage, with one selection block per depth
         (top-1 shown once, then each best-of-top-N): the representative poses
         (complexes) removed AT that stage and the percentages (of the panel start N
         and of the poses entering the stage), plus the running survivor count.
    """
    def pm(m: str) -> str:
        if labels and m in labels:
            return labels[m]
        return _panel_label(m)

    def _mcnemar_lines(a_name: str, b_name: str, rec: dict, indent: str) -> list[str]:
        lo_a, hi_a = rec["a_ci"]; lo_b, hi_b = rec["b_ci"]
        pp = 100 * (rec["b_rate"] - rec["a_rate"])
        ls = [
            f"{indent}{a_name}: {rec['a_k']}/{rec['n']} = {100 * rec['a_rate']:.1f}% "
            f"[95% CI {100 * lo_a:.1f}–{100 * hi_a:.1f}]   →   "
            f"{b_name}: {rec['b_k']}/{rec['n']} = {100 * rec['b_rate']:.1f}% "
            f"[95% CI {100 * lo_b:.1f}–{100 * hi_b:.1f}]   (Δ {pp:+.1f} pp)",
            f"{indent}  discordant pairs: b={rec['a_wins']} ({a_name} only), "
            f"c={rec['b_wins']} ({b_name} only)  ·  "
            f"exact McNemar p={su.fmt_p(rec['mcnemar_p'])}"
            + (f", Holm-adjusted p={su.fmt_p(rec['p_holm'])}"
               if "p_holm" in rec else "")
            + (f"  {rec['star']}" if rec.get("star") else ""),
        ]
        return ls

    depths = [b["depth"] for b in depth_blocks]
    # Selections shown, in order: top-1 (once) then each best-of-top-N depth.
    selections = [("top-1", cascades_top1)] + [
        (f"best-of-top-{b['depth']}", b["cascades_topn"]) for b in depth_blocks]

    L: list[str] = []
    W = L.append
    bar = "=" * 78
    W(bar)
    for i, fn in enumerate(figure_names):
        W(f"{'Figures' if i == 0 else '       '} : {fn}")
    W("Data    : combined statistics + per-stage attrition for the raw-vs-refined")
    W("          PoseBusters recovery figures, at THREE pool depths —")
    W("          top-1 (rank-1 baseline) · "
      + " · ".join(f"best-of-top-{d}" for d in depths))
    for i, cn in enumerate(csv_names):
        W(f"{'Per-step CSVs' if i == 0 else '            '} : {cn}")
    W(bar)
    W("")
    W("WHAT THESE FIGURES SHOW")
    W("  For each method, ONE representative pose per complex is pushed through the")
    W("  sequential PoseBusters filter cascade: first drop poses with RMSD > 2 Å from")
    W(f"  the {_ref_ligand_phrase()}, then apply each canonical PoseBusters physical-validity")
    W("  test in order, removing at every step the poses that fail it having passed")
    W("  all previous ones. Three pose selections are compared per method:")
    W("    top-1 (rank-1)    = the tool's top-ranked pose")
    W("                        (EquiBind has no native rank → ranked by gnina affinity).")
    for d in depths:
        W(f"    best-of-top-{d:<6}= the lowest-RMSD pose among the tool's first {d} ranked")
        W("                        poses (one pose per complex → same denominator as top-1).")
    W("  'Failed at a stage' = the representative pose was removed by that stage's test.")
    W("  All counts are COMPLEXES (one representative pose each) — no pose-level")
    W("  pseudoreplication.")
    W("")
    W("PANELS  (n = complexes with a representative pose)")
    for m in methods:
        c = cascades_top1.get(m)
        if not c:
            for _, casc in selections:
                if m in casc:
                    c = casc[m]
                    break
        n = c["N"] if c else "?"
        W(f"    {pm(m):<26}  n = {n}")
    W("")
    W(bar)
    W("1) STATISTICAL TESTS")
    W(bar)
    W("  Outcome     : passing ALL canonical PoseBusters tests (binary, per complex).")
    W("  Design      : paired — the SAME complexes are compared across the two")
    W("                conditions; one representative pose per complex.")
    W("  rank gain   : rank-1  vs  best-of-top-N  passing-all rate. Tests whether")
    W("                deepening the ranked pool recovers more fully-valid poses.")
    W("  refine gain : raw variant  vs  refined variant passing-all rate, at BOTH the")
    W("                rank-1 and the best-of-top-N depth (refined panels only). Tests")
    W("                whether smina/gnina refinement raises pass-all.")
    W("  Test        : exact (binomial) two-sided McNemar on the discordant pairs")
    W("                (b, c); rates carry Wilson score 95% confidence intervals.")
    W("  Multiplicity: each best-of-top-N depth is a SEPARATE figure; every McNemar")
    W("                p-value WITHIN a figure forms one Holm–Bonferroni family (p_holm),")
    W("                matching that figure's on-panel stars. The two depths are NOT")
    W("                pooled into a single family. Stars: *** p<.001, ** p<.01,")
    W("                * p<.05, ns = not significant (on the Holm p).")
    W("")
    for blk in depth_blocks:
        d, stats = blk["depth"], blk["stats"]
        W(f"  ── rank-1 vs best-of-top-{d}  (own Holm family) "
          + "─" * max(0, 30 - len(str(d))))
        if not stats:
            W("    (no statistics available — the paired maps were empty)")
            W("")
            continue
        for m in methods:
            st = stats.get(str(m))
            if not st:
                continue
            W(f"    [{pm(m)}]")
            rg = st.get("rank_gain")
            if rg:
                W(f"      rank gain — top-1 vs best-of-top-{d}:")
                L.extend(_mcnemar_lines("top-1", f"top-{d}", rg, "        "))
            ref = st.get("refine_gain") or {}
            for dkey, dlab in (("rank1", "top-1"),
                               ("topn", f"best-of-top-{d}")):
                rec = ref.get(dkey)
                if rec:
                    W(f"      refine gain — raw vs refined @ {dlab}:")
                    L.extend(_mcnemar_lines("raw", "refined", rec, "        "))
            W("")
    W(bar)
    W("2) PER-STAGE ATTRITION TABLE  (';'-separated)")
    W(bar)
    W("  One row per method × selection × cascade stage; selections are top-1 (shown")
    W("  once) then each best-of-top-N. 'failed' = representative poses (complexes)")
    W("  removed at that stage; 'entering' = poses reaching it.")
    W("  pct_failed_of_start   = failed / N          (N = panel start count)")
    W("  pct_failed_of_entering= failed / entering   (stage-conditional failure rate)")
    W("  pct_remaining_of_start= remaining / N       (running survival)")
    W("  Stage kinds: start (all predictions) · drop (a filter) · milestone (RMSD ≤ 2 Å")
    W("  running total) · end (passing all tests).")
    W("")
    header = ["method", "selection", "stage_order", "stage", "kind",
              "entering", "failed", "remaining",
              "pct_failed_of_start", "pct_failed_of_entering",
              "pct_remaining_of_start"]
    W(";".join(header))
    for sel_name, casc in selections:
        for m in methods:
            c = casc.get(m)
            if not c:
                continue
            N = c["N"] or 1
            for order, s in enumerate(c["steps"]):
                entering = int(s["before"])
                failed = int(s["removed"])
                remaining = int(s["remaining"])
                pfs = 100 * failed / N
                pfe = (100 * failed / entering) if entering else float("nan")
                prs = 100 * remaining / N
                row = [
                    pm(m).replace(";", ","), sel_name, str(order),
                    str(s["label"]).replace(";", ","), s["kind"],
                    str(entering), str(failed), str(remaining),
                    f"{pfs:.2f}", (f"{pfe:.2f}" if entering else "NA"),
                    f"{prs:.2f}"]
                W(";".join(row))
    W("")
    out_path.write_text("\n".join(L) + "\n", encoding="utf-8")


def _dd_refine_label(variant: "str | None") -> str:
    """The 'DiffDock*' refined-panel label for figure 17h, reflecting which
    optimizer variant the report actually collapsed DiffDock to — the ``best_dd``
    returned by ``_select_best_diffdock`` (set by --collapse-diffdock-variant, or
    the oracle-best pick). The collapse relabels that winner to the ``diffdock``
    method key, so the panel data follow it; this keeps the label honest instead
    of hard-coding one refinement. Falls back to a generic 'refined' when the
    variant can't be resolved."""
    v = str(variant or "")
    if v.endswith("_smina"):
        return "DiffDock* (smina-refined)"
    if v.endswith("_gnina"):
        return "DiffDock* (gnina-refined)"
    if v == "diffdock":
        return "DiffDock* (raw)"
    return "DiffDock* (refined)"


def _emit_pb_recovery_raw_vs_refined(depth: int, *, df: pd.DataFrame,
                                     df_full: pd.DataFrame,
                                     test_table: pd.DataFrame, test_cols: list[str],
                                     cascades: dict, out_dir: Path,
                                     stats_sidecar: dict,
                                     dd_variant: "str | None" = None
                                     ) -> "dict | None":
    """Build + write the combined raw-vs-refined PoseBusters recovery figure (17h)
    at rank-1 vs best-of-top-``depth``, plus its per-step CSV, and fold the
    McNemar/Wilson payload into ``stats_sidecar`` under
    ``top{depth}_recovery_raw_vs_refined``.

    The companion text report is NOT written here — instead this returns the built
    cascades/stats so the caller can emit ONE combined report across all depths
    (top-1 · best-of-top-15 · best-of-top-30). Returns ``None`` if the figure could
    not be built (missing variants), else a dict with keys ``depth``, ``c1_all``,
    ``cn_all``, ``methods``, ``labels``, ``stats``, ``png``, ``csv``.

    ``cascades`` is the rank-1 (top-1) refined cascade dict (depth-independent, so
    reused across depths); the best-of-top-``depth`` refined cascades and both raw
    cascades are (re)built here so the figure is self-contained per depth.
    """
    cn_ref = aggregate_pb_waterfall(df, test_table, test_cols, depth,
                                    rank_selection="topn")
    if not cn_ref:
        return None
    c1 = dict(cascades)
    cn = dict(cn_ref)
    _inject_equibind_gnina_headroom(c1, cn, df_full, test_table, test_cols, depth)
    c1_raw, cn_raw = _build_raw_variant_cascades(df_full, test_table, test_cols, depth)
    c1_all = {**c1_raw, **c1}
    cn_all = {**cn_raw, **cn}
    combined_methods = tuple(
        m for m in ("autodock", "diffdock_raw", "diffdock",
                    "equibind_raw", "equibind_gnina")
        if m in c1_all and m in cn_all)
    if len(combined_methods) < 2:
        return None
    combined_labels = {
        "diffdock": _dd_refine_label(dd_variant),
        "diffdock_raw": "DiffDock (raw)",
        "equibind_gnina": "EquiBind* (gnina-refined)",
        "equibind_raw": "EquiBind* (raw)",
    }
    try:
        stats = _stats_top30_recovery(c1_all, cn_all, combined_methods)
    except Exception as exc:
        stats = None
        print(f"  WARNING: 17h top-{depth} recovery stats failed ({exc})")
    if stats:
        stats_sidecar[f"top{depth}_recovery_raw_vs_refined"] = stats
    png = out_dir / f"17h_pb_top{depth}_recovery_raw_vs_refined.png"
    csv = out_dir / f"pb_top{depth}_recovery_raw_vs_refined.csv"
    plot_pb_top30_recovery_by_test(
        c1_all, cn_all, depth, png, csv_out=csv,
        methods=combined_methods, labels=combined_labels, stats=stats)
    print(f"  wrote combined raw+refined top-{depth} recovery → {png.name} "
          f"(+ {csv.name})")
    return {"depth": depth, "c1_all": c1_all, "cn_all": cn_all,
            "methods": combined_methods, "labels": combined_labels,
            "stats": stats, "png": png.name, "csv": csv.name}


# ───────────────────────────────────────────────────────────────────
# Raw-vs-best-variant PoseBusters waterfall (17d)
# ───────────────────────────────────────────────────────────────────


def _pose_stem_no_refine(pose_name: str) -> str:
    """EquiBind pose id shared across refinement variants: strip the trailing
    ``__refRAW`` / ``__refSMINA`` / ``__refGNINA`` (and any ``.sdf``) so the raw
    and smina-refined copies of one prediction map to the same key."""
    return re.sub(r"__ref(?:raw|smina|gnina)(?:\.sdf)?$", "",
                  str(pose_name), flags=re.IGNORECASE)


def _rank1_rep(df: pd.DataFrame, method: str) -> pd.DataFrame:
    """The rank-1 pose per complex for a natively-ranked method key (AutoDock /
    DiffDock / its optimiser variants). One row per (protein, ligand)."""
    sub = df[df["method"].astype(str) == method].copy()
    sub = sub[pd.to_numeric(sub["rank"], errors="coerce").notna()]
    if sub.empty:
        return sub
    sub["rank"] = pd.to_numeric(sub["rank"], errors="coerce")
    return (sub.sort_values("rank")
               .groupby(["protein", "ligand"]).head(1).reset_index(drop=True))


def _best_equibind_gnina_pocket(df: pd.DataFrame) -> str | None:
    """The EquiBind pocket whose ``equibind_{pocket}_gnina`` variant has the
    highest oracle PB-Valid AND RMSD ≤ 2 Å rate (the same combined criterion the
    rest of the report ranks variants on). ``None`` if no gnina variant exists."""
    methods = df["method"].astype(str)
    pockets = sorted({m[len("equibind_"):-len("_gnina")] for m in methods.unique()
                      if m.startswith("equibind_") and m.endswith("_gnina")})
    if not pockets:
        return None
    best, best_score = None, -1.0
    for p in pockets:
        sub = df[methods == f"equibind_{p}_gnina"].dropna(subset=["rmsd"])
        if sub.empty:
            continue
        orc = sub.loc[sub.groupby(["protein", "ligand"])["rmsd"].idxmin()]
        score = float(((orc["rmsd"] <= 2.0)
                       & orc["pb_valid"].astype(bool)).mean())
        if score > best_score:
            best, best_score = p, score
    return best or pockets[0]


def _equibind_gnina_ranked_top1(df: pd.DataFrame, pocket: str
                                ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per complex, the EquiBind prediction ranked #1 by its gnina refinement
    score (lowest ``gnina_affinity`` among ``equibind_{pocket}_gnina`` poses),
    returned as ``(raw_rep, gnina_rep)`` — the SAME prediction in its raw and its
    gnina-refined geometry (one row per complex, restricted to complexes matched
    in both). EquiBind has no native ranking, so this is how a user would pick a
    single pose: keep the gnina-best one and compare its raw vs refined form.
    Empty frames when the pocket's raw/gnina variants are missing."""
    empty = df.iloc[0:0][["method", "protein", "ligand", "pose_name", "rmsd"]]
    gnina = df[df["method"].astype(str) == f"equibind_{pocket}_gnina"].copy()
    raw = df[df["method"].astype(str) == f"equibind_{pocket}_raw"].copy()
    if gnina.empty or raw.empty or "gnina_affinity" not in gnina:
        return empty, empty
    gnina = gnina.dropna(subset=["gnina_affinity"])
    if gnina.empty:
        return empty, empty
    gnina["_stem"] = gnina["pose_name"].map(_pose_stem_no_refine)
    raw["_stem"] = raw["pose_name"].map(_pose_stem_no_refine)
    idx = gnina.groupby(["protein", "ligand"])["gnina_affinity"].idxmin()
    gnina_top = gnina.loc[idx].drop_duplicates(["protein", "ligand"])
    key = ["protein", "ligand", "_stem"]
    raw_u = raw.drop_duplicates(key)
    merged = gnina_top.merge(raw_u[key + ["method", "pose_name", "rmsd"]],
                             on=key, how="inner", suffixes=("_gnina", "_raw"))
    cols = ["method", "protein", "ligand", "pose_name", "rmsd"]
    gnina_rep = merged.rename(columns={"method_gnina": "method",
                                       "pose_name_gnina": "pose_name",
                                       "rmsd_gnina": "rmsd"})[cols]
    raw_rep = merged.rename(columns={"method_raw": "method",
                                     "pose_name_raw": "pose_name",
                                     "rmsd_raw": "rmsd"})[cols]
    return raw_rep.reset_index(drop=True), gnina_rep.reset_index(drop=True)


def _equibind_gnina_rank_headroom(df: pd.DataFrame, pocket: str, top_n: int
                                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank EquiBind ``{pocket}_gnina`` poses by their gnina affinity (best =
    lowest) within each complex and return ``(top1_rep, topn_rep)``:
    ``top1_rep`` = the gnina rank-1 pose per complex; ``topn_rep`` = the
    lowest-RMSD pose among that complex's top-``top_n`` gnina-ranked poses.

    This is the gnina-score analogue of ``_top1_per_pair`` / ``_best_topn_per_pair``
    (which need a native ``rank``), letting the ranking-headroom figure (17c) carry
    an EquiBind panel even though EquiBind emits no native pose ranking. Both frames
    carry method/protein/ligand/pose_name/rmsd; empty when the variant is absent."""
    cols = ["method", "protein", "ligand", "pose_name", "rmsd"]
    sub = df[df["method"].astype(str) == f"equibind_{pocket}_gnina"].copy()
    sub = sub.dropna(subset=["gnina_affinity"])
    if sub.empty:
        return sub[cols], sub[cols]
    sub["_srank"] = (sub.groupby(["protein", "ligand"])["gnina_affinity"]
                        .rank(method="first", ascending=True))
    top1 = sub[sub["_srank"] == 1][cols].reset_index(drop=True)
    pool = sub[sub["_srank"] <= top_n].dropna(subset=["rmsd"])
    if pool.empty:
        return top1, top1
    topn = (pool.loc[pool.groupby(["protein", "ligand"])["rmsd"].idxmin()][cols]
                .reset_index(drop=True))
    return top1, topn


def _inject_equibind_gnina_headroom(cascades_top1: dict, cascades_topn: dict,
                                    df_full: pd.DataFrame, test_table: pd.DataFrame,
                                    test_cols: list[str], top_n: int) -> bool:
    """Add a gnina-ranked EquiBind panel (key ``equibind_gnina``) to the top-1 and
    best-of-top-N cascade dicts used by the ranking-headroom figure (17c).

    ``df_full`` must retain the ``equibind_{pocket}_gnina`` variant (the full
    per-pose frame). Mutates both dicts in place; returns True when a panel was
    added. Callers should pass shallow copies if the base cascades are reused
    elsewhere (17/17b)."""
    pocket = _best_equibind_gnina_pocket(df_full)
    if not pocket:
        return False
    top1_rep, topn_rep = _equibind_gnina_rank_headroom(df_full, pocket, top_n)
    c_top1 = _pb_cascade_from_rep(top1_rep, test_table, test_cols,
                                  "gnina rank-1 pose")
    c_topn = _pb_cascade_from_rep(topn_rep, test_table, test_cols,
                                  f"best of top-{top_n} gnina-ranked")
    if c_top1 is None or c_topn is None:
        return False
    cascades_top1["equibind_gnina"] = c_top1
    cascades_topn["equibind_gnina"] = c_topn
    return True


def _best_topn_rep(df: pd.DataFrame, method: str, top_n: int) -> pd.DataFrame:
    """The lowest-RMSD pose among a natively-ranked method's top-``top_n`` ranks,
    one row per complex (single-method analogue of _best_topn_per_pair)."""
    cols = ["method", "protein", "ligand", "pose_name", "rmsd"]
    s = df[df["method"].astype(str) == method].copy()
    s["rank"] = pd.to_numeric(s["rank"], errors="coerce")
    s = s[s["rank"] <= top_n].dropna(subset=["rmsd"])
    if s.empty:
        return s.reindex(columns=cols)
    return (s.loc[s.groupby(["protein", "ligand"])["rmsd"].idxmin()][cols]
             .reset_index(drop=True))


def _equibind_raw_gnina_ranked_headroom(df: pd.DataFrame, pocket: str, top_n: int
                                        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Like _equibind_gnina_rank_headroom, but returns the RAW EquiBind geometry
    while still ranking by the gnina refinement score (the raw poses carry no score
    of their own). For each complex the ``equibind_{pocket}_raw`` poses are ordered
    by their gnina counterpart's affinity (matched on pose stem); top1 = gnina
    rank-1 raw pose, topn = lowest-RMSD raw pose among the top-``top_n`` gnina ranks.
    This is the un-refined mirror of the gnina figures (raw geometry, gnina order)."""
    cols = ["method", "protein", "ligand", "pose_name", "rmsd"]
    gnina = df[df["method"].astype(str) == f"equibind_{pocket}_gnina"].dropna(
        subset=["gnina_affinity"]).copy()
    raw = df[df["method"].astype(str) == f"equibind_{pocket}_raw"].copy()
    if gnina.empty or raw.empty:
        return raw.reindex(columns=cols), raw.reindex(columns=cols)
    gnina["_stem"] = gnina["pose_name"].map(_pose_stem_no_refine)
    raw["_stem"] = raw["pose_name"].map(_pose_stem_no_refine)
    gnina["_sr"] = (gnina.groupby(["protein", "ligand"])["gnina_affinity"]
                    .rank(method="first", ascending=True))
    order = gnina[["protein", "ligand", "_stem", "_sr"]].drop_duplicates(
        ["protein", "ligand", "_stem"])
    rawm = raw.merge(order, on=["protein", "ligand", "_stem"], how="inner")
    top1 = rawm[rawm["_sr"] == 1][cols].reset_index(drop=True)
    pool = rawm[rawm["_sr"] <= top_n].dropna(subset=["rmsd"])
    if pool.empty:
        return top1, top1
    topn = (pool.loc[pool.groupby(["protein", "ligand"])["rmsd"].idxmin()][cols]
             .reset_index(drop=True))
    return top1, topn


def _build_raw_variant_cascades(df: pd.DataFrame, test_table: pd.DataFrame,
                                test_cols: list[str], top_n: int
                                ) -> tuple[dict, dict]:
    """Top-1 and best-of-top-N cascades for the UN-REFINED variants — the raw
    counterpart of the gnina figures (17c2/17f). Keys: ``autodock`` (physics,
    unchanged), ``diffdock_raw`` (raw DiffDock, native rank), ``equibind_raw`` (raw
    EquiBind geometry, ranked by gnina score). ``df`` must be the full per-variant
    frame. Returns ({}, {}) if the required variants are absent."""
    c1: dict = {}
    cn: dict = {}
    for key, method in (("autodock", "autodock"), ("diffdock_raw", "diffdock")):
        a = _pb_cascade_from_rep(_rank1_rep(df, method), test_table, test_cols,
                                 "rank-1 pose")
        b = _pb_cascade_from_rep(_best_topn_rep(df, method, top_n), test_table,
                                 test_cols, f"best of top-{top_n}")
        if a is not None and b is not None:
            c1[key], cn[key] = a, b
    pocket = _best_equibind_gnina_pocket(df)
    if pocket:
        t1, tn = _equibind_raw_gnina_ranked_headroom(df, pocket, top_n)
        a = _pb_cascade_from_rep(t1, test_table, test_cols,
                                 "raw pose · gnina rank-1")
        b = _pb_cascade_from_rep(tn, test_table, test_cols,
                                 f"best of top-{top_n} raw (gnina-ranked)")
        if a is not None and b is not None:
            c1["equibind_raw"], cn["equibind_raw"] = a, b
    return c1, cn


def _resolve_dd_best_variant(df: pd.DataFrame, forced: str | None) -> str:
    """The DiffDock optimiser variant to use as the 'best refined' DiffDock panel in
    the raw-vs-best waterfalls (17d/17e) and the rank-depth added-pose figures (17i),
    so they FOLLOW ``--collapse-diffdock-variant`` like the rest of the report instead
    of hard-coding one refiner.

    A 'best refined' panel needs a *refined* variant, so only the smina/gnina
    optimisers are eligible (raw ``diffdock`` is the counterpart it's compared
    against). Honours an explicit smina/gnina ``forced`` choice; otherwise (flag
    absent, or forced to raw) defaults to ``diffdock_gnina`` — the historical 17d
    pick and the report's usual oracle-best. Falls back to whatever refiner is
    present if the preferred one is missing."""
    present = set(df["method"].astype(str))
    refined = [v for v in ("diffdock_smina", "diffdock_gnina") if v in present]
    if not refined:
        return "diffdock_gnina"
    if forced in refined:
        return forced
    return "diffdock_gnina" if "diffdock_gnina" in refined else refined[0]


def _dd_refiner_word(dd_best: str) -> str:
    """Short label ('smina'/'gnina') for a DiffDock optimiser method key."""
    return "smina" if "smina" in str(dd_best) else "gnina"


def aggregate_pb_waterfall_raw_vs_best(df: pd.DataFrame, test_table: pd.DataFrame,
                                       test_cols: list[str],
                                       gate_rmsd: bool = True,
                                       dd_best: str = "diffdock_gnina") -> dict:
    """PoseBusters cascade, top-1 pose, for the RAW pose vs its best refined
    variant — as an ordered dict of panels ``{key: {..cascade.., title, role, pair,
    vs, opt}}``. ``dd_best`` selects the DiffDock refiner (follows
    ``--collapse-diffdock-variant``); EquiBind is always gnina-refined.

    ``df`` must retain every optimiser/refinement variant under its own method key
    (the full per-pose frame, before any best-variant collapse). Selection:

    * AutoDock — rank-1 pose. Physics docking has no smina/gnina re-optimisation,
      so it appears once as a reference (``role="reference"``).
    * DiffDock — rank-1 raw pose (``diffdock``) vs the rank-1 ``dd_best`` pose; the
      optimiser inherits the DiffDock rank, so the two panels are the same pose
      before/after refinement.
    * EquiBind — has no native ranking, so the representative is the prediction
      ranked #1 by its gnina refinement score, shown in its raw and its
      gnina-refined geometry (``_equibind_gnina_ranked_top1``).

    ``gate_rmsd`` (default True) drops RMSD > 2 Å before the PoseBusters tests;
    with ``gate_rmsd=False`` every ranked-1 pose faces the tests regardless of
    placement accuracy, so the cascade shows where poses fail on physical-validity
    grounds alone. Best panels carry ``vs`` = the key of their raw counterpart (so
    the plotter can annotate the net change in poses passing every test) and ``opt``
    = the refiner word for the panel label.
    """
    out: dict = {}
    dd_word = _dd_refiner_word(dd_best)

    ref = _pb_cascade_from_rep(_rank1_rep(df, "autodock"), test_table, test_cols,
                               "top-1 pose", gate_rmsd=gate_rmsd)
    if ref is not None:
        out["autodock"] = {**ref, "role": "reference", "pair": "autodock",
                           "title": "AutoDock Vina — top-1 pose "
                                    "(physics docking · no smina/gnina re-optimisation)"}

    dd_raw = _pb_cascade_from_rep(_rank1_rep(df, "diffdock"), test_table,
                                  test_cols, "raw pose", gate_rmsd=gate_rmsd)
    dd_best_c = _pb_cascade_from_rep(_rank1_rep(df, dd_best), test_table,
                                     test_cols, f"{dd_word}-refined pose",
                                     gate_rmsd=gate_rmsd)
    if dd_raw is not None:
        out["diffdock_raw"] = {**dd_raw, "role": "raw", "pair": "diffdock",
                               "title": "DiffDock — raw top-1 pose"}
    if dd_best_c is not None:
        e = {**dd_best_c, "role": "best", "pair": "diffdock", "opt": dd_word,
             "title": f"DiffDock + {dd_word} — best variant, top-1 pose"}
        if dd_raw is not None:
            e["vs"] = "diffdock_raw"
        out["diffdock_best"] = e

    pocket = _best_equibind_gnina_pocket(df)
    if pocket:
        eq_raw_rep, eq_gnina_rep = _equibind_gnina_ranked_top1(df, pocket)
        eq_raw = _pb_cascade_from_rep(eq_raw_rep, test_table, test_cols, "raw pose",
                                      gate_rmsd=gate_rmsd)
        eq_best = _pb_cascade_from_rep(eq_gnina_rep, test_table, test_cols,
                                       "gnina-refined pose", gate_rmsd=gate_rmsd)
        if eq_raw is not None:
            out["equibind_raw"] = {
                **eq_raw, "role": "raw", "pair": "equibind",
                "title": f"EquiBind ({pocket}) — raw pose · ranked #1 by gnina score"}
        if eq_best is not None:
            e = {**eq_best, "role": "best", "pair": "equibind", "opt": "gnina",
                 "title": f"EquiBind ({pocket}) + gnina — refined pose · "
                          "ranked #1 by gnina score"}
            if eq_raw is not None:
                e["vs"] = "equibind_raw"
            out["equibind_best"] = e
    return out


def plot_pb_waterfall_raw_vs_best(panels: dict, out: Path,
                                  csv_out: Path | None = None,
                                  gate_rmsd: bool = True) -> None:
    """Small-multiples PoseBusters waterfall comparing each method's RAW top-1
    pose with its best refined variant (DiffDock: smina; EquiBind: gnina).

    Same visual language as ``plot_pb_waterfall`` (red bar = poses a test removes;
    teal milestone = RMSD ≤ 2 Å survivors; green = poses passing every test), but
    panels are ordered raw-then-optimised per method and share a coloured left
    edge (grey = raw, green = optimised, neutral = the AutoDock reference). Each
    optimised panel is annotated with its net gain versus the paired raw panel.

    ``gate_rmsd`` must match the flag the cascades were built with: when False the
    panels have no RMSD step (every ranked-1 pose is tested for physical validity
    regardless of placement), and the labels/suptitle switch to PB-validity wording.
    """
    keys = list(panels.keys())
    if not keys:
        return
    pass_label = "passing all tests" if gate_rmsd else "PB-valid (any RMSD)"
    delta_label = "passing all tests" if gate_rmsd else "PB-valid"
    labels = [s["label"] for s in panels[keys[0]]["steps"]]
    nx = len(labels)
    nrows = len(keys)
    fig, axes = plt.subplots(nrows, 1, sharex=True, squeeze=False,
                             figsize=(max(11, 0.62 * nx), 2.7 * nrows + 2.0))
    axes = list(axes[:, 0])
    teal, red, green = "#3a9d8f", "#d1495b", "#2a9d3f"
    edge_col = {"reference": "#9aa0a6", "raw": "#c77b86", "best": green}
    csv_rows = []

    def _final_pct(c: dict) -> float:
        return 100 * c["steps"][-1]["remaining"] / c["N"] if c["N"] else float("nan")

    for ax, key in zip(axes, keys):
        c = panels[key]; N = c["N"]; steps = c["steps"]; role = c.get("role", "raw")
        x = np.arange(nx)
        for xi, s in zip(x, steps):
            b = 100 * s["before"] / N if N else 0.0
            a = 100 * s["after"] / N if N else 0.0
            if s["kind"] == "start":
                ax.bar(xi, 100, color=teal, alpha=0.30, edgecolor=teal,
                       hatch="..", zorder=2)
                ax.text(xi, 101, f"{N}", ha="center", va="bottom",
                        fontsize=8, fontweight="bold")
            elif s["kind"] == "milestone":
                ax.bar(xi, a, color=teal, alpha=0.45, edgecolor=teal,
                       hatch="//", zorder=2)
                ax.text(xi, a + 1, f"{s['remaining']}", ha="center",
                        va="bottom", fontsize=8)
            elif s["kind"] == "end":
                ax.bar(xi, a, color=green, alpha=0.85, edgecolor="black", zorder=3)
                ax.text(xi, a + 1, f"{s['remaining']}\n({a:.1f}%)", ha="center",
                        va="bottom", fontsize=8, fontweight="bold")
            else:  # drop
                if s["removed"] > 0:
                    ax.bar(xi, b - a, bottom=a, color=red, alpha=0.9,
                           edgecolor="black", zorder=3)
                    ax.text(xi, b + 0.5, f"−{s['removed']}", ha="center",
                            va="bottom", fontsize=7, color=red, fontweight="bold")
                else:
                    ax.text(xi, a + 0.5, "0", ha="center", va="bottom",
                            fontsize=7, color="0.5")
        run = [100 * s["after"] / N if N else 0.0 for s in steps]
        ax.step(x, run, where="mid", color="0.4", lw=0.8, ls="--",
                alpha=0.7, zorder=1)
        ax.set_ylim(0, 110)
        ax.set_ylabel("% of poses\nremaining", fontsize=9)
        final_pct = _final_pct(c)
        ax.set_title(f"{c.get('title', key)}   "
                     f"(n={N}; {pass_label} = {final_pct:.1f}%)",
                     fontsize=10, fontweight="bold")
        # Coloured left edge signalling raw / optimised / reference.
        ax.spines["left"].set_color(edge_col.get(role, "0.3"))
        ax.spines["left"].set_linewidth(3.2)
        ax.grid(axis="y", alpha=0.25)

        vs = c.get("vs")
        if vs in panels:
            raw_pct = _final_pct(panels[vs])
            d = final_pct - raw_pct
            accent = green if d >= 0 else red
            # Refiner word for the panel (follows --collapse-diffdock-variant for
            # DiffDock; EquiBind is always gnina). Falls back for older callers.
            opt_word = c.get("opt", "gnina" if c.get("pair") == "equibind" else "smina")
            ax.text(0.995, 0.975,
                    f"{opt_word} vs raw:  {'+' if d >= 0 else '−'}{abs(d):.1f} pp "
                    f"{delta_label}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=9.5, fontweight="bold", color=accent,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec=accent, alpha=0.92))
            csv_rows.append({"pair": c.get("pair"),
                             "raw_passing_all_%": round(raw_pct, PCT_DECIMALS),
                             "best_passing_all_%": round(final_pct, PCT_DECIMALS),
                             "gain_pp": round(d, PCT_DECIMALS)})

    axes[-1].set_xticks(np.arange(nx))
    axes[-1].set_xticklabels(labels, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    dd_word = panels.get("diffdock_best", {}).get("opt", "smina")
    if gate_rmsd:
        suptitle = ("PoseBusters test-failure waterfall — raw pose vs best "
                    f"refined variant, top-1 per method (DiffDock: {dd_word}; EquiBind: "
                    "gnina)\n(sequential filter: each red bar = poses removed by "
                    "that test; green = poses passing every test; left edge: grey = "
                    "raw, green = optimised)")
    else:
        suptitle = ("PoseBusters physical-validity waterfall — raw pose vs best "
                    f"refined variant, top-1 per method (DiffDock: {dd_word}; EquiBind: "
                    "gnina)\n(RMSD ≤ 2 Å NOT applied — every ranked-1 pose is tested; "
                    "each red bar = poses removed by that test; green = poses passing "
                    "every PoseBusters test, at any RMSD)")
    fig.suptitle(_vt(suptitle), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        _write_csv(pd.DataFrame(csv_rows), csv_out, index=False)


# ───────────────────────────────────────────────────────────────────
# Rank-depth 1→N ADDED-pose PoseBusters analysis (17i) — NO RMSD gate
# ───────────────────────────────────────────────────────────────────
# "How much does deepening the rank cut buy you, and at what physical-validity
# cost?" Instead of one representative pose per complex (17d/17e), this pools the
# poses a user NEWLY considers when the rank cut widens from 1 to N — ranking
# positions 2..N — and asks, WITHOUT the RMSD ≤ 2 Å placement gate: how many poses
# and complexes does the deeper cut add, and which PoseBusters tests do those added
# poses fail. Native rankers use their rank order; EquiBind (no native rank) is
# ordered by its gnina refinement score.


def _positions_native(df: pd.DataFrame, method: str, top_n: int
                      ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(rank1_rep, added_pool)`` for a natively-ranked method key, by rank ORDER:
    position 1 = the best-ranked pose per complex; the ADDED pool = positions
    2..``top_n`` — every pose gained by deepening the rank cut from 1 to ``top_n``.
    Rows carry method/protein/ligand/pose_name/rmsd/pb_valid.

    Raw DiffDock exports its top pose twice (a bare ``rankN.sdf`` copy beside the
    confidence-named file), so for the raw ``diffdock`` key the bare duplicate is
    dropped first — otherwise it would masquerade as a spurious extra added pose and
    the optimiser variants (confidence-named only) wouldn't align 1:1 by position."""
    cols = ["method", "protein", "ligand", "pose_name", "rmsd", "pb_valid"]
    s = df[df["method"].astype(str) == method].copy()
    if method == "diffdock":
        bare = s["pose_name"].astype(str).str.contains(
            r"rank\d+\.sdf$", regex=True, case=False, na=False)
        s = s[~bare]
    s["_r"] = pd.to_numeric(s["rank"], errors="coerce")
    s = s.dropna(subset=["_r"])
    if s.empty:
        e = s.reindex(columns=cols)
        return e, e
    s["_pos"] = s.groupby(["protein", "ligand"])["_r"].rank(method="first")
    r1 = s[s["_pos"] == 1][cols].reset_index(drop=True)
    add = s[(s["_pos"] >= 2) & (s["_pos"] <= top_n)][cols].reset_index(drop=True)
    return r1, add


def _positions_equibind(df: pd.DataFrame, pocket: str, top_n: int, geometry: str
                        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(rank1_rep, added_pool)`` for EquiBind, ranked by gnina affinity (lowest =
    best) — EquiBind has no native rank. ``geometry="gnina"`` returns the refined
    poses; ``geometry="raw"`` returns the SAME predictions' raw geometry (matched on
    pose stem, kept in the gnina rank order). ADDED pool = gnina positions 2..``top_n``.
    Same columns as ``_positions_native``; empty frames when the variant is absent."""
    cols = ["method", "protein", "ligand", "pose_name", "rmsd", "pb_valid"]
    gnina = df[df["method"].astype(str) == f"equibind_{pocket}_gnina"].copy()
    gnina = gnina.dropna(subset=["gnina_affinity"])
    if gnina.empty:
        e = gnina.reindex(columns=cols)
        return e, e
    gnina["_pos"] = gnina.groupby(["protein", "ligand"])["gnina_affinity"].rank(
        method="first", ascending=True)
    if geometry == "gnina":
        s = gnina
    else:                                   # raw geometry at the gnina-ranked slots
        gnina["_stem"] = gnina["pose_name"].map(_pose_stem_no_refine)
        raw = df[df["method"].astype(str) == f"equibind_{pocket}_raw"].copy()
        if raw.empty:
            e = raw.reindex(columns=cols)
            return e, e
        raw["_stem"] = raw["pose_name"].map(_pose_stem_no_refine)
        order = gnina[["protein", "ligand", "_stem", "_pos"]].drop_duplicates(
            ["protein", "ligand", "_stem"])
        s = raw.merge(order, on=["protein", "ligand", "_stem"], how="inner")
    r1 = s[s["_pos"] == 1][cols].reset_index(drop=True)
    add = s[(s["_pos"] >= 2) & (s["_pos"] <= top_n)][cols].reset_index(drop=True)
    return r1, add


def _added_pool_stats(rank1_rep: pd.DataFrame, added_pool: pd.DataFrame,
                      test_table: pd.DataFrame, test_cols: list[str]) -> dict:
    """Complex- and pose-level counts for one panel's rank-1 pose vs its
    positions-2..N added pool (NO RMSD gate; passing-all == PB-valid):

      * ``n_complex_rank1`` / ``n_complex_with_added`` — complexes with a rank-1
        pose / with ≥1 added pose;
      * ``rank1_pbvalid`` — complexes whose rank-1 pose is already PB-valid;
      * ``complexes_recovered`` — complexes whose rank-1 pose is NOT PB-valid but
        that gain a PB-valid pose somewhere in positions 2..N (the deepening payoff);
      * ``added_poses`` / ``added_pbvalid`` / ``added_fail_any`` — added-pose totals;
      * ``fail_by_test`` — the INDEPENDENT per-test failure count over the added pool
        (a pose is counted for EVERY test it fails, so — unlike the sequential
        cascade — later tests aren't masked by earlier ones)."""
    def _ck(d: pd.DataFrame) -> pd.Series:
        return d["protein"].astype(str) + "\x1f" + d["ligand"].astype(str)

    r1 = rank1_rep.drop_duplicates(["protein", "ligand"]).copy()
    r1["_c"] = _ck(r1)
    r1valid = dict(zip(r1["_c"], r1["pb_valid"].astype(bool)))
    add = added_pool.copy()
    n_added = len(add)
    stats = {
        "n_complex_rank1": int(len(r1)),
        "rank1_pbvalid": int(r1["pb_valid"].astype(bool).sum()) if len(r1) else 0,
        "added_poses": int(n_added),
        "added_pbvalid": 0, "added_fail_any": 0,
        "n_complex_with_added": 0, "complexes_recovered": 0,
        "fail_by_test": {},
    }
    if not n_added:
        return stats
    add["_c"] = _ck(add)
    addvalid_any = add.groupby("_c")["pb_valid"].any()
    stats["n_complex_with_added"] = int(add["_c"].nunique())
    stats["added_pbvalid"] = int(add["pb_valid"].astype(bool).sum())
    stats["added_fail_any"] = int(n_added - stats["added_pbvalid"])
    stats["complexes_recovered"] = int(sum(
        1 for c, av in addvalid_any.items() if av and not r1valid.get(c, False)))
    # Independent per-test failure tally: join the per-test PASS/FAIL columns on
    # pose identity (protein, ligand, pose_name) — NOT method — exactly as the
    # cascade does, then count explicit failures per test (missing join = pass).
    tt = test_table.drop(columns=["method"], errors="ignore")
    m = add.merge(tt, how="left", on=["protein", "ligand", "pose_name"])
    m = m.drop_duplicates(subset=["protein", "ligand", "pose_name"])
    for col in test_cols:
        failed = ((~m[col].fillna(True).astype(bool)).sum() if col in m else 0)
        stats["fail_by_test"][col] = int(failed)
    return stats


def aggregate_pb_rank_depth_added_pool(df: pd.DataFrame, test_table: pd.DataFrame,
                                       test_cols: list[str], top_n: int,
                                       dd_best: str = "diffdock_gnina") -> dict:
    """Ranking-depth 1→``top_n`` analysis WITHOUT the RMSD gate. Per method (raw +
    best refined variant), pools the NEWLY ADDED poses — ranking positions
    2..``top_n``, i.e. every pose gained by deepening the rank cut from 1 to
    ``top_n`` — and runs the sequential PoseBusters validity cascade over that pool
    (green = added poses PB-valid at any RMSD). Returns an ordered panels dict shaped
    like ``aggregate_pb_waterfall_raw_vs_best`` (each panel carries the cascade +
    role/pair/vs/title/short); each also carries a ``depth`` block with the
    added-pose / complexes-recovered counts and the independent per-test fail tally.

    The DiffDock refiner follows ``dd_best`` (set by ``--collapse-diffdock-variant``);
    EquiBind (no native rank) is ranked and refined by its gnina score. ``df`` must
    retain every optimiser/refinement variant under its own method key."""
    out: dict = {}
    dd_word = _dd_refiner_word(dd_best)
    plan = [
        ("autodock", ("native", "autodock"), "reference", "autodock",
         "AutoDock Vina",
         "AutoDock Vina — added poses, ranks 2–{n} (physics docking · no re-optimisation)"),
        ("diffdock_raw", ("native", "diffdock"), "raw", "diffdock",
         "DiffDock (raw)", "DiffDock — raw added poses, ranks 2–{n}"),
        ("diffdock_best", ("native", dd_best), "best", "diffdock",
         f"DiffDock + {dd_word}",
         f"DiffDock + {dd_word} — refined added poses, ranks 2–{{n}}"),
    ]
    pocket = _best_equibind_gnina_pocket(df)
    if pocket:
        plan += [
            ("equibind_raw", ("equibind", (pocket, "raw")), "raw", "equibind",
             f"EquiBind ({pocket}, raw)",
             f"EquiBind ({pocket}) — raw added poses · gnina ranks 2–{{n}}"),
            ("equibind_best", ("equibind", (pocket, "gnina")), "best", "equibind",
             f"EquiBind ({pocket}) + gnina",
             f"EquiBind ({pocket}) + gnina — refined added poses · gnina ranks 2–{{n}}"),
        ]
    for key, sel, role, pair, short, title in plan:
        kind, spec = sel
        if kind == "native":
            r1, add = _positions_native(df, spec, top_n)
        else:
            pk, geom = spec
            r1, add = _positions_equibind(df, pk, top_n, geom)
        if add.empty:
            continue
        cascade = _pb_cascade_from_rep(add, test_table, test_cols,
                                       f"added poses (ranks 2–{top_n})",
                                       gate_rmsd=False)
        if cascade is None:
            continue
        entry = {**cascade, "role": role, "pair": pair, "short": short,
                 "title": title.format(n=top_n),
                 "depth": _added_pool_stats(r1, add, test_table, test_cols)}
        if role == "best":
            entry["opt"] = "gnina" if pair == "equibind" else dd_word
            if f"{pair}_raw" in out:
                entry["vs"] = f"{pair}_raw"
        out[key] = entry
    return out


def plot_pb_rank_depth_added_waterfall(panels: dict, top_n: int, out: Path,
                                       csv_out: Path | None = None,
                                       summary_csv: Path | None = None) -> None:
    """Small-multiples PoseBusters waterfall over the poses ADDED by deepening the
    rank cut from 1 to ``top_n`` (positions 2..``top_n``), one panel per method
    (raw + best refined variant). No RMSD gate — every added pose faces the physical-
    validity tests. Each red bar = added poses that test removes; green = added poses
    passing every test (PB-valid). A per-panel headline states the deepening payoff:
    how many poses/complexes the deeper cut adds and how many complexes it recovers
    (rank-1 not PB-valid → an added pose PB-valid). Best panels also annotate the
    raw→refined change in added-pool PB-valid yield (refiner per --collapse-diffdock-variant)."""
    keys = list(panels.keys())
    if not keys:
        return
    labels = [s["label"] for s in panels[keys[0]]["steps"]]
    nx = len(labels)
    nrows = len(keys)
    fig, axes = plt.subplots(nrows, 1, sharex=True, squeeze=False,
                             figsize=(max(11, 0.62 * nx), 2.9 * nrows + 2.2))
    axes = list(axes[:, 0])
    teal, red, green = "#3a9d8f", "#d1495b", "#2a9d3f"
    edge_col = {"reference": "#9aa0a6", "raw": "#c77b86", "best": green}
    summ_rows = []

    def _final_pct(c: dict) -> float:
        return 100 * c["steps"][-1]["remaining"] / c["N"] if c["N"] else float("nan")

    for ax, key in zip(axes, keys):
        c = panels[key]; N = c["N"]; steps = c["steps"]; role = c.get("role", "raw")
        dep = c.get("depth", {})
        x = np.arange(nx)
        for xi, s in zip(x, steps):
            b = 100 * s["before"] / N if N else 0.0
            a = 100 * s["after"] / N if N else 0.0
            if s["kind"] == "start":
                ax.bar(xi, 100, color=teal, alpha=0.30, edgecolor=teal,
                       hatch="..", zorder=2)
                ax.text(xi, 101, f"{N}", ha="center", va="bottom",
                        fontsize=8, fontweight="bold")
            elif s["kind"] == "end":
                ax.bar(xi, a, color=green, alpha=0.85, edgecolor="black", zorder=3)
                ax.text(xi, a + 1, f"{s['remaining']}\n({a:.1f}%)", ha="center",
                        va="bottom", fontsize=8, fontweight="bold")
            else:  # drop
                if s["removed"] > 0:
                    ax.bar(xi, b - a, bottom=a, color=red, alpha=0.9,
                           edgecolor="black", zorder=3)
                    ax.text(xi, b + 0.5, f"−{s['removed']}", ha="center",
                            va="bottom", fontsize=7, color=red, fontweight="bold")
                else:
                    ax.text(xi, a + 0.5, "0", ha="center", va="bottom",
                            fontsize=7, color="0.5")
        run = [100 * s["after"] / N if N else 0.0 for s in steps]
        ax.step(x, run, where="mid", color="0.4", lw=0.8, ls="--",
                alpha=0.7, zorder=1)
        ax.set_ylim(0, 122)
        ax.set_ylabel("% of added\nposes remaining", fontsize=9)
        final_pct = _final_pct(c)
        ax.set_title(f"{c.get('title', key)}   (added poses n={N}; "
                     f"PB-valid = {final_pct:.1f}%)",
                     fontsize=10, fontweight="bold")
        ax.spines["left"].set_color(edge_col.get(role, "0.3"))
        ax.spines["left"].set_linewidth(3.2)
        ax.grid(axis="y", alpha=0.25)
        # The deepening payoff — the headline the figure is really about.
        head = (f"depth 1→{top_n}:  +{dep.get('added_poses', 0)} added poses across "
                f"{dep.get('n_complex_with_added', 0)} complexes   ·   "
                f"+{dep.get('complexes_recovered', 0)} complexes recovered "
                f"(rank-1 not PB-valid → an added pose PB-valid)")
        ax.text(0.5, 0.99, head, transform=ax.transAxes, ha="center", va="top",
                fontsize=8, color="0.12",
                bbox=dict(boxstyle="round,pad=0.3", fc="#eef3fa", ec="0.6",
                          alpha=0.95), zorder=6)
        vs = c.get("vs")
        if vs in panels:
            d = final_pct - _final_pct(panels[vs])
            accent = green if d >= 0 else red
            opt_word = c.get("opt", "gnina")
            ax.text(0.995, 0.03,
                    f"{opt_word} vs raw:  {'+' if d >= 0 else '−'}{abs(d):.1f} pp PB-valid",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=9, fontweight="bold", color=accent,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=accent,
                              alpha=0.92))
        row = {"panel": key, "role": role, "pair": c.get("pair"),
               "added_poses": dep.get("added_poses", 0),
               "complexes_with_added": dep.get("n_complex_with_added", 0),
               "complexes_recovered": dep.get("complexes_recovered", 0),
               "rank1_pbvalid_complexes": dep.get("rank1_pbvalid", 0),
               "added_pbvalid": dep.get("added_pbvalid", 0),
               "added_fail_any": dep.get("added_fail_any", 0),
               "added_pbvalid_pct": round(final_pct, PCT_DECIMALS)}
        for col, cnt in dep.get("fail_by_test", {}).items():
            row[f"failind_{col}"] = cnt
        summ_rows.append(row)

    axes[-1].set_xticks(np.arange(nx))
    axes[-1].set_xticklabels(labels, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    fig.suptitle(_vt(
        f"PoseBusters validity of the poses ADDED by deepening the rank cut 1 → {top_n}"
        f"\n(added poses = ranking positions 2–{top_n}; RMSD ≤ 2 Å NOT applied — pure "
        "physical validity; each red bar = added poses removed by that test; "
        "green = added poses passing every PoseBusters test)"),
        fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None:
        _waterfall_to_csv(panels, csv_out)
    if summary_csv is not None and summ_rows:
        _write_csv(pd.DataFrame(summ_rows), summary_csv, index=False)


def plot_pb_rank_depth_fail_by_test(panels: dict, top_n: int, out: Path,
                                    csv_out: Path | None = None) -> None:
    """Companion to the added-pose waterfall: of the poses ADDED by deepening the
    rank cut 1→``top_n``, how many fail EACH PoseBusters test — the INDEPENDENT tally
    (a pose is counted for every test it fails), so later tests aren't masked by
    earlier ones as they are in the sequential cascade. One panel per method (shared
    test rows), horizontal bars annotated with the count and its % of that panel's
    added poses; tests that never fail anywhere are dropped."""
    keys = list(panels.keys())
    if not keys:
        return
    totals: dict = {}
    for c in panels.values():
        for t, n in c.get("depth", {}).get("fail_by_test", {}).items():
            totals[t] = totals.get(t, 0) + n
    tests = [t for t in sorted(totals, key=lambda k: -totals[k]) if totals[t] > 0]
    if not tests:
        return
    labels_y = [PB_TEST_LABELS.get(t, t) for t in tests]
    y = np.arange(len(tests))[::-1]         # most-failed test at the top
    nrows = len(keys)
    bar_col = {"reference": "#6b7280", "raw": "#c77b86", "best": "#2a9d3f"}
    fig, axes = plt.subplots(1, nrows, sharey=True, squeeze=False,
                             figsize=(3.5 * nrows + 1.2, 0.46 * len(tests) + 2.6))
    axes = list(axes[0, :])
    csv_rows = []
    xmax = max((c.get("depth", {}).get("fail_by_test", {}).get(t, 0)
                for c in panels.values() for t in tests), default=1) or 1
    for ax, key in zip(axes, keys):
        c = panels[key]; dep = c.get("depth", {})
        fbt = dep.get("fail_by_test", {})
        n_added = dep.get("added_poses", 0)
        denom = n_added or 1
        vals = [fbt.get(t, 0) for t in tests]
        ax.barh(y, vals, color=bar_col.get(c.get("role"), "#888"),
                alpha=0.85, edgecolor="white", height=0.72)
        for yi, v in zip(y, vals):
            if v > 0:
                ax.text(v + xmax * 0.01, yi, f"{v} ({100 * v / denom:.0f}%)",
                        va="center", ha="left", fontsize=7)
        ax.set_xlim(0, xmax * 1.28)
        ax.set_title(f"{c.get('short', key)}\n(added n={n_added})",
                     fontsize=9.5, fontweight="bold")
        ax.set_xlabel("added poses failing test", fontsize=8)
        ax.tick_params(axis="x", labelsize=7)
        ax.grid(axis="x", alpha=0.25)
        for t in tests:
            csv_rows.append({"panel": key, "role": c.get("role"),
                             "test": t, "test_label": PB_TEST_LABELS.get(t, t),
                             "added_poses_failing": fbt.get(t, 0),
                             "added_poses": n_added})
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(labels_y, fontsize=8)
    fig.suptitle(_vt(
        f"Which PoseBusters tests the rank-1→{top_n} ADDED poses fail\n"
        f"(added poses = ranking positions 2–{top_n}; independent per-test tally — a "
        "pose is counted for every test it fails; RMSD not applied)"),
        fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        _write_csv(pd.DataFrame(csv_rows), csv_out, index=False)


# ───────────────────────────────────────────────────────────────────
# Index builder
# ───────────────────────────────────────────────────────────────────


# Provenance columns the PoseBusters CSV carries for EquiBind poses (written
# from the pose SDF tags by run_posebusters). Preferred over filename parsing
# and carried through to per_pose_metrics.csv for drill-down.
EQ_PROVENANCE_COLS = ("pocket_source", "clamp_variant", "refine_variant",
                      "smina_affinity", "gnina_affinity", "pocket_id")


def _col_value(row, key: str) -> str | None:
    """Read a provenance column from a row, treating empty/nan/unknown as absent."""
    v = row.get(key) if hasattr(row, "get") else None
    if v is None:
        return None
    s = str(v).strip()
    return s if s and s.lower() not in ("nan", "none", "unknown", "<na>") else None


def _classify_equibind(row) -> tuple[str, str | None, str | None]:
    """Map one EquiBind pose row to (pocket, refine, clamp) variant tokens.

    Prefers the provenance columns the PoseBusters CSV now carries
    (``pocket_source`` / ``refine_variant`` / ``clamp_variant``, written from the
    pose SDF tags); falls back to parsing the pose filename for older CSVs.

    pocket : "unguided" (blind) | "fpocket" | "p2rank" (| "guided" legacy)
    refine : "smina" (__refSMINA) | "gnina" (__refGNINA) | "raw" (__refRAW) | None (unsuffixed)
    clamp  : "clampON" (__clampON) | "clampOFF" (__clampOFF) | None
    """
    name = Path(str(row.get("pose_name", ""))).name.lower()   # strip "lig__prot/" prefix

    pocket = _col_value(row, "pocket_source")
    if pocket not in ("unguided", "fpocket", "p2rank"):
        if name.startswith("unguided"):
            pocket = "unguided"
        elif name.startswith("p2rank"):
            pocket = "p2rank"
        elif name.startswith("fpocket"):
            pocket = "fpocket"
        else:
            pocket = "guided"          # legacy / unrecognised → coarse bucket

    refine = _col_value(row, "refine_variant")
    if refine not in ("smina", "raw", "gnina"):
        refine = ("smina" if "__refsmina" in name
                  else "gnina" if "__refgnina" in name
                  else "raw" if "__refraw" in name else None)

    clamp = _col_value(row, "clamp_variant")
    if clamp not in ("clampON", "clampOFF"):
        clamp = ("clampOFF" if "__clampoff" in name
                 else "clampON" if "__clampon" in name else None)

    return pocket, refine, clamp


def _apply_equibind_split(df: pd.DataFrame) -> pd.DataFrame:
    """Refine ``docking_method`` for EquiBind rows into per-variant labels
    (e.g. ``equibind_fpocket_smina``) from the CSV provenance columns / filename."""
    is_eq = df["docking_method"].str.startswith("equibind")
    if not is_eq.any():
        return df
    eq_idx = df.index[is_eq]
    labels = []
    for i in eq_idx:
        pocket, refine, clamp = _classify_equibind(df.loc[i])
        lab = f"equibind_{pocket}"
        if refine:
            lab += f"_{refine}"
        if clamp:
            lab += f"_{clamp}"
        labels.append(lab)
    df.loc[eq_idx, "docking_method"] = labels
    return df


def _classify_autodock(row) -> str:
    """Return the independent optimizer identity for one AutoDock pose.

    Explicit PoseBusters provenance is authoritative. Method/path fallbacks are
    only for older result CSVs and deliberately recognize optimized directories
    so legacy optimized poses cannot be pooled with raw Vina.
    """
    opt = (_col_value(row, "optimizer") or "").lower()
    if opt in _AUTODOCK_OPTIMIZERS:
        return opt
    if opt in ("original", "raw", "native", "none"):
        return "original"
    method = str(row.get("docking_method", "")).strip().lower()
    if method.endswith("_gnina_refinement"):
        return "gnina_refinement"
    if method.endswith("_smina"):
        return "smina"
    if method.endswith("_gnina"):
        return "gnina"
    name = " ".join(str(row.get(k, "")) for k in ("pose_name", "pose_file")).lower()
    if "optimized_gnina_refinement" in name:
        return "gnina_refinement"
    if "optimized_smina" in name:
        return "smina"
    if "optimized_gnina" in name:
        return "gnina"
    return "original"


def _autodock_scoring_base(method: str) -> str:
    """The AutoDock scoring base of a method key — 'autodock' (Vina) or
    'autodock_vinardo' — with any existing optimizer suffix stripped, so the
    split stays idempotent and never merges Vinardo into Vina."""
    mm = str(method).strip().lower()
    for suf in _AUTODOCK_OPTIMIZER_SUFFIXES:
        if mm.endswith(suf):
            mm = mm[: -len(suf)]
            break
    # Return the base VERBATIM rather than folding onto autodock/autodock_vinardo.
    # Collapsing here silently merged the MGLTools-ligand and exhaustiveness-64
    # trees into the plain "autodock" rows, so the cascade reported one
    # "AutoDock Vina (raw)" row pooling three different preparations. Existing
    # keys are unaffected: "autodock"/"autodock_gnina" still reduce to "autodock"
    # and Vinardo is still never merged into Vina.
    return mm


def _apply_autodock_split(df: pd.DataFrame) -> pd.DataFrame:
    """Separate native and optimized AutoDock poses into stable method keys,
    preserving the scoring base (Vina vs Vinardo)."""
    is_ad = df["docking_method"].str.startswith("autodock")
    if not is_ad.any():
        return df
    ad_idx = df.index[is_ad]
    opts = [_classify_autodock(df.loc[i]) for i in ad_idx]
    bases = [_autodock_scoring_base(df.loc[i, "docking_method"]) for i in ad_idx]
    df.loc[ad_idx, "docking_method"] = [
        base if opt == "original" else f"{base}_{opt}"
        for base, opt in zip(bases, opts)
    ]
    # A stable explicit axis also survives into the pose index/cache even when
    # the input was a legacy CSV without an optimizer column.
    if "optimizer" not in df.columns:
        df["optimizer"] = pd.NA
    df.loc[ad_idx, "optimizer"] = opts
    return df


def _classify_diffdock(row) -> str | None:
    """Optimizer backend for one DiffDock pose: 'smina' | 'gnina' | None (original).

    Prefers the ``optimizer`` provenance column the PoseBusters CSV now carries
    (written by run_posebusters from each pose's ``optimized_<tool>/`` subfolder);
    falls back to the ``optimized_<tool>`` path component in ``pose_name``.
    """
    opt = _col_value(row, "optimizer")
    if opt in ("smina", "gnina"):
        return opt
    if opt == "original":
        return None
    name = str(row.get("pose_name", "")).lower()
    if "optimized_smina" in name:
        return "smina"
    if "optimized_gnina" in name:
        return "gnina"
    return None


def _apply_diffdock_split(df: pd.DataFrame) -> pd.DataFrame:
    """Relabel smina/gnina-optimised DiffDock poses' ``docking_method`` into
    ``diffdock_smina`` / ``diffdock_gnina`` (raw poses keep ``diffdock``).

    The optimised variants are deliberately NOT named ``diffdock`` so they fall
    outside RANKING_TOOLS and are treated oracle-only (like EquiBind): an
    optimised copy inherits its source pose's DiffDock rank, so re-ranking it
    would double-count rank-1. Keeping them separate also stops their RMSDs being
    pooled into raw DiffDock's oracle success rate."""
    is_dd = df["docking_method"].str.startswith("diffdock")
    if not is_dd.any():
        return df
    dd_idx = df.index[is_dd]
    labels = []
    for i in dd_idx:
        o = _classify_diffdock(df.loc[i])
        labels.append(f"diffdock_{o}" if o else "diffdock")
    df.loc[dd_idx, "docking_method"] = labels
    return df


def _select_diffdock_variant(df: pd.DataFrame, variant: str = "diffdock") -> pd.DataFrame:
    """Keep ONLY the chosen DiffDock optimizer variant and treat it as the
    canonical ``diffdock`` tool; drop the other variants.

    Run AFTER _apply_diffdock_split (methods are diffdock / diffdock_smina /
    diffdock_gnina). With the default ``"diffdock"`` this drops the optimised
    copies and keeps raw DiffDock. Selecting ``"diffdock_smina"`` /
    ``"diffdock_gnina"`` keeps that variant and relabels it to ``diffdock`` so it
    is ranked and compared like the docking tool (one DiffDock entry per run,
    mirroring the EquiBind variant selection)."""
    if variant == "all":               # keep every variant as a separate method
        return df
    is_dd = df["docking_method"].str.startswith("diffdock")
    if not is_dd.any():
        return df
    sel = df["docking_method"] == variant
    if not sel.any():
        print(f"  [warn] --diffdock-variant {variant!r} not present in this CSV; "
              f"available: {sorted(df.loc[is_dd, 'docking_method'].unique())}")
    df = df[~is_dd | sel].copy()
    df.loc[df["docking_method"].str.startswith("diffdock"), "docking_method"] = "diffdock"
    return df


def _load_allowed_ids(path: Path) -> set[str]:
    """Load '<PDBID>_<LIG>' complex ids (one per line; '#' comments ignored)."""
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def _build_pose_index(pb_csv: Path, split_equibind: bool = True,
                      diffdock_variant: str = "diffdock",
                      select_diffdock: bool = True) -> pd.DataFrame:
    df = pd.read_csv(pb_csv, low_memory=False)
    needed = {"docking_method", "protein", "ligand", "pose_file", "pose_name"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"PB CSV missing columns: {missing}")
    bool_df = pd.DataFrame({c: _to_bool(df[c])
                            for c in PB_CRITICAL_CHECKS if c in df.columns})
    df["pb_valid"] = bool_df.all(axis=1) if not bool_df.empty else False
    df["docking_method"] = df["docking_method"].astype(str).str.lower()
    if split_equibind:
        df = _apply_equibind_split(df)
    # Split AutoDock before any per-pair grouping: raw uses its original Vina
    # rank, optimized variants use their own optimized_rank downstream.
    df = _apply_autodock_split(df)
    # Split DiffDock optimizer variants, then (unless we keep all three for
    # --best-diffdock-only) keep only the selected one (raw 'diffdock' by default)
    # as the canonical DiffDock entry.
    df = _apply_diffdock_split(df)
    if select_diffdock:
        df = _select_diffdock_variant(df, diffdock_variant)
    keep = ["docking_method", "protein", "ligand", "pose_file", "pose_name", "pb_valid"]
    keep += [c for c in EQ_PROVENANCE_COLS if c in df.columns]
    # Preserve both AutoDock rank axes and their scoring provenance in the heavy
    # per-pose cache. This is what makes gnina/smina top-N use optimized_rank while
    # retaining the source Vina rank for paired before/after analysis.
    auto_cols = (
        "optimizer", "autodock_rank", "optimized_rank", "rank_metric",
        "autodock_affinity", "minimized_affinity", "cnn_score", "cnn_affinity",
    )
    for col in auto_cols:
        if col not in df.columns:
            df[col] = pd.NA
        if col not in keep:
            keep.append(col)
    for col in ("unidock_rank", "unidock_affinity", "unidock2_rank",
                "unidock2_affinity", "variant", "scoring"):
        if col in df.columns and col not in keep:
            keep.append(col)
    return df[keep]


# ───────────────────────────────────────────────────────────────────
# Per-pose cache (so re-runs that only change --top-n skip the heavy scoring)
# ───────────────────────────────────────────────────────────────────

# Columns a cached per_pose_metrics.csv must contain to be reusable by the current
# code (guards against reusing a CSV written by an older, narrower schema).
_CACHE_REQUIRED_COLS = {
    "method", "protein", "ligand", "pose_name", "rank", "rmsd",
    "centroid_dist", "pb_valid", "optimizer", "autodock_rank",
    "optimized_rank",
    # schema 6: reference-ligand convention provenance
    "reference_convention", "n_copies", "nearest_copy_index", "rmsd_ref_instance",
}


def _file_fingerprint(path: Path) -> dict:
    """Content fingerprint of a file: size + sha256 of its bytes.

    Keyed on *content*, not mtime, so the cache is NOT invalidated when the PB
    CSV is regenerated byte-identically upstream (a re-export bumps mtime but not
    content). Hashing ~100 MB takes well under a second — negligible next to the
    pair scoring it guards. Returns ``sha256=None`` if the file is unreadable.
    """
    try:
        st = path.stat()
    except OSError:
        return {"path": str(path), "size": None, "sha256": None}
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return {"path": str(path), "size": st.st_size, "sha256": None}
    return {"path": str(path), "size": st.st_size, "sha256": h.hexdigest()}


def _per_pose_signature(args) -> dict:
    """Fingerprint of the inputs the per-pose metrics depend on.

    The heavy per-pose scoring is a pure function of the PB CSV *content*, the id
    filter, the EquiBind split and --limit-pairs — but NOT of --top-n (that only
    drives the cheap downstream ranking aggregation). So a cache keyed on this
    signature is reusable across different --top-n values.
    """
    return {
        "schema": 6,  # 5: AutoDock variants + ranking provenance; 6: reference convention
        "pb_csv": _file_fingerprint(Path(args.pb_csv)),
        "ids_file": str(args.ids_file) if args.ids_file else None,
        "split_equibind": bool(args.split_equibind),
        "diffdock_variant": str(args.diffdock_variant),
        "best_diffdock_only": bool(args.best_diffdock_only),
        "limit_pairs": int(args.limit_pairs),
        "reference_convention": _args_reference_convention(args),
    }


def _args_reference_convention(args) -> str:
    return str(getattr(args, "reference_convention", None) or "instance")


def _manifest_reference_convention(man: dict) -> "str | None":
    """Convention a manifest was built under: the top-level key, else the one
    inside the signature, else ``None`` (legacy schema ≤ 5 manifest)."""
    if not isinstance(man, dict):
        return None
    conv = man.get("reference_convention")
    if conv is None:
        conv = (man.get("signature") or {}).get("reference_convention")
    if conv is None and isinstance(man.get("signature"), dict):
        # Every cache written before schema 6 was scored against the single
        # deposited instance; treat the missing key as that convention so the
        # refusals below apply to legacy caches too.
        return "instance"
    return str(conv) if conv is not None else None


def _read_cached_per_pose(per_pose_csv: Path):
    """Read + validate the cached per-pose CSV; ``None`` if missing/invalid."""
    if not per_pose_csv.exists():
        return None
    try:
        df = pd.read_csv(per_pose_csv, low_memory=False)
    except Exception:
        return None
    missing = sorted(c for c in _CACHE_REQUIRED_COLS if c not in df.columns)
    if missing:
        print("  cached per-pose metrics lack required columns — cannot reuse. "
              f"Missing: {', '.join(missing)} (a per_pose_metrics.csv written before "
              "schema 6 must be rebuilt with --force).")
        return None
    return df


def _load_cached_per_pose(args, sig: dict):
    """Return the cached per-pose DataFrame, or ``None`` to trigger a recompute.

    With ``--reuse-cache`` the existing per_pose_metrics.csv is loaded WITHOUT
    checking the input signature (use when only the plotting/aggregation code
    changed — the caller asserts the cached pairs are the ones they want). This
    bypasses ``--force`` too. The ONE thing --reuse-cache does check is the
    reference convention: a cache whose manifest records the other convention is
    refused (``args.cache_refused`` is set and ``None`` returned), because every
    crystal-referenced column would silently mean something else. Otherwise the
    cache is reused only when the manifest signature matches the current inputs;
    a signature that differs ONLY in the reference convention REFUSES rather than
    recomputing unless ``--force`` is given (schema / input mismatches recompute
    as before).
    """
    per_pose_csv = args.out_dir / "per_pose_metrics.csv"
    manifest = args.out_dir / "per_pose_metrics.manifest.json"
    want = _args_reference_convention(args)

    if getattr(args, "reuse_cache", False):
        if manifest.exists():
            try:
                man = json.loads(manifest.read_text())
            except Exception:
                man = {}
            have = _manifest_reference_convention(man)
            if have is not None and have != want:
                print(f"ERROR: --reuse-cache: the cache in {args.out_dir} was built under "
                      f"--reference-convention {have}, but this run requested {want}. "
                      "Refusing to reuse it (every crystal-referenced column would "
                      "change meaning). Point --out-dir at a cache built under "
                      f"{want}, or rebuild with --force.")
                args.cache_refused = f"reference_convention {have} != {want}"
                return None
        df = _read_cached_per_pose(per_pose_csv)
        if df is not None:
            print(f"--reuse-cache: loaded {per_pose_csv} ({len(df):,} rows) "
                  "WITHOUT checking the input signature — skipping pair scoring.")
        return df

    if args.force or not per_pose_csv.exists() or not manifest.exists():
        return None
    try:
        man = json.loads(manifest.read_text())
    except Exception:
        return None
    if man.get("signature") != sig:
        old_sig = man.get("signature") or {}
        if isinstance(old_sig, dict):
            diff = sorted(k for k in set(old_sig) | set(sig)
                          if old_sig.get(k) != sig.get(k))
        else:
            diff = ["signature"]
        have = _manifest_reference_convention(man)
        old_schema = old_sig.get("schema") if isinstance(old_sig, dict) else None
        # A cache whose CONVENTION or SCHEMA differs is never rebuilt implicitly:
        # both would silently replace a table every downstream sidecar depends
        # on with one that means something else (or that was built by another
        # generation of this script). Only --force may do that; every other
        # signature change (PB CSV content, id filter, variant flags) keeps the
        # long-standing recompute behaviour.
        if have != want or "reference_convention" in diff or old_schema != sig.get("schema"):
            others = [k for k in diff if k not in ("reference_convention", "schema")]
            why = []
            if have != want or "reference_convention" in diff:
                why.append(f"built under --reference-convention {have}, run requested {want}")
            if old_schema != sig.get("schema"):
                why.append(f"cache schema {old_schema} predates the current schema {sig.get('schema')}")
            if others:
                why.append("other inputs also differ: " + ", ".join(others))
            print(f"ERROR: per-pose cache in {args.out_dir}: " + "; ".join(why) +
                  ". Refusing to rebuild it implicitly; pass --force to rebuild in place "
                  "(back the directory up first) or point --out-dir at a new directory.")
            args.cache_refused = "; ".join(why)
            return None
        print("  cached per-pose metrics are stale (inputs changed: "
              f"{', '.join(diff)}) — recomputing.")
        return None
    df = _read_cached_per_pose(per_pose_csv)
    if df is None:
        return None
    print(f"Reusing cached per-pose metrics: {per_pose_csv} ({len(df):,} rows).")
    print(f"  inputs unchanged; previous top-n={man.get('top_n')}, "
          f"now top-n={args.top_n} (per-pose metrics are top-n independent). "
          f"Pass --force to recompute.")
    return df


def _write_per_pose_cache(args, df: pd.DataFrame, sig: dict) -> Path:
    per_pose_csv = args.out_dir / "per_pose_metrics.csv"
    _write_csv(df, per_pose_csv, index=False)
    (args.out_dir / "per_pose_metrics.manifest.json").write_text(
        json.dumps({"signature": sig, "top_n": int(args.top_n),
                    "n_rows": int(len(df)),
                    "reference_convention": _args_reference_convention(args)},
                   indent=2))
    return per_pose_csv


# ───────────────────────────────────────────────────────────────────
# Main
# ───────────────────────────────────────────────────────────────────


def _parse_threshold_list(s: str) -> tuple[float, ...]:
    """Parse a comma/space-separated list of RMSD thresholds into a sorted tuple."""
    vals = tuple(sorted({float(x) for x in str(s).replace(",", " ").split() if x}))
    if not vals:
        raise argparse.ArgumentTypeError("no thresholds parsed from %r" % s)
    return vals


def _write_rmsd_vs_pbvalid_table(oracle_sum: pd.DataFrame, top1_sum: pd.DataFrame,
                                 out: Path) -> pd.DataFrame:
    """Per-variant comparison of RMSD-only vs RMSD+PB-validity success.

    Puts the plain docking-accuracy criterion (RMSD-to-crystal ≤ 2 Å) next to the
    stricter one that ALSO requires the pose to pass PoseBusters (RMSD ≤ 2 Å AND
    PB-valid), for both the oracle (best) pose and the top-1 ranked pose, plus the
    ``cost`` of demanding validity (the percentage-point drop). Sourced from the
    columns aggregate_oracle/aggregate_top1 already compute, so the numbers match the
    12_rmsd2_vs_pbvalid_grouped / 13_pbvalid_filter_influence figures exactly."""
    cols: dict[str, pd.Series] = {}
    if "oracle_rmsd_le_2.0A_%" in oracle_sum:
        cols["oracle_rmsd_le_2A_%"] = oracle_sum["oracle_rmsd_le_2.0A_%"]
    if "oracle_pb_valid_and_rmsd2_%" in oracle_sum:
        cols["oracle_rmsd2_and_pbvalid_%"] = oracle_sum["oracle_pb_valid_and_rmsd2_%"]
    if not top1_sum.empty:
        if "top1_rmsd_le_2.0A_%" in top1_sum:
            cols["top1_rmsd_le_2A_%"] = top1_sum["top1_rmsd_le_2.0A_%"]
        if "top1_pb_valid_and_rmsd2_%" in top1_sum:
            cols["top1_rmsd2_and_pbvalid_%"] = top1_sum["top1_pb_valid_and_rmsd2_%"]
    tbl = pd.DataFrame(cols)
    if {"oracle_rmsd_le_2A_%", "oracle_rmsd2_and_pbvalid_%"} <= set(tbl.columns):
        tbl["oracle_pbvalid_cost_pp"] = (tbl["oracle_rmsd_le_2A_%"]
                                         - tbl["oracle_rmsd2_and_pbvalid_%"])
    if {"top1_rmsd_le_2A_%", "top1_rmsd2_and_pbvalid_%"} <= set(tbl.columns):
        tbl["top1_pbvalid_cost_pp"] = (tbl["top1_rmsd_le_2A_%"]
                                       - tbl["top1_rmsd2_and_pbvalid_%"])
    tbl.index.name = "method"
    _write_csv(_round_csv(tbl), out)
    return tbl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pb-csv", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "posebusters_filtered_results.csv"))
    ap.add_argument("--benchmark-dir", type=Path,
                    default=Path("Data/PoseBuster Benchmark Set"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "pose_comparison_report"))
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 1))
    ap.add_argument("--limit-pairs", type=int, default=0,
                    help="Process only the first N pairs (debug).")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Restrict analysis to the '<PDBID>_<LIG>' complex ids "
                         "listed in this file (one per line; '#' comments ok). "
                         "Pairs whose id is absent are dropped.")
    ap.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                    help="Ranked poses per pair to include in ranking analysis "
                         "(default %(default)s).")
    ap.add_argument("--fine-rmsd-thresholds", type=_parse_threshold_list,
                    default=FINE_RMSD_THRESHOLDS,
                    help="Comma/space-separated RMSD thresholds (Å) for the "
                         "top-ranked-pose within-distance sweep (figure 18 / "
                         "topn_within_thresholds.csv). Default: "
                         + ",".join(f"{t:g}" for t in FINE_RMSD_THRESHOLDS) + ".")
    ap.add_argument("--force", action="store_true",
                    help="Force a full recompute of the per-pose metrics even when "
                         "a cached per_pose_metrics.csv matching the inputs exists. "
                         "By default those cached results are reused (so re-running "
                         "with a different --top-n skips the heavy pose scoring, "
                         "since the per-pose metrics are top-n independent).")
    ap.add_argument("--reuse-cache", action="store_true",
                    help="Load the existing per_pose_metrics.csv and SKIP the heavy "
                         "pair scoring regardless of whether the input signature "
                         "matches (and regardless of --force). Use when only the "
                         "plotting/aggregation changed and you want to re-render the "
                         "figures from the cached scores. You are responsible for the "
                         "cache matching your intended inputs (e.g. the same "
                         "--ids-file). Errors out if no cache is present.")
    ap.add_argument("--form-ok-kabsch", type=float, default=FORM_OK_KABSCH_A,
                    help="Best-fit (Kabsch) RMSD (Å) below which a near-native, "
                         "PB-valid pose counts as having the CORRECT form / internal "
                         "conformation (figure 20 / form_fidelity_summary.csv). "
                         "Default %(default)s.")
    ap.add_argument("--pocket-cutoff", type=float, default=POCKET_CENTROID_CUTOFF,
                    help="Centroid distance (Å) to the reference copy of the crystal "
                         "ligand (per --reference-convention: the single instance, "
                         "or the pose's nearest deposited copy) within which a pose "
                         "counts as being in the experimentally validated pocket, for "
                         "the pocket-localization analysis (default %(default)s).")
    ap.add_argument("--reference-convention", choices=list(REFERENCE_CONVENTIONS),
                    default="instance",
                    help="Which deposited copy of the crystal ligand every "
                         "crystal-referenced per-pose metric (rmsd, pb_rmsd, "
                         "pb_kabsch_rmsd, pb_rmsd_within_2A, bestfit_rmsd, "
                         "centroid_dist, torsions, contact_recovery, plif_recovery) "
                         "is measured against. 'instance': record 0 of "
                         "<ID>_ligands.sdf (== <ID>_ligand.sdf), always — the "
                         "default, reproduces the frozen tables exactly. 'nearest': "
                         "per pose, the copy of <ID>_ligands.sdf with the smallest "
                         "symmetry-corrected in-place RMSD (PoseBusters check_rmsd "
                         "convention); that one copy drives every metric of the "
                         "pose. The convention is part of the per-pose cache "
                         "signature and manifest; a cache built under the other "
                         "convention is refused (see --reuse-cache / --force). "
                         "Default %(default)s.")
    ap.add_argument("--split-equibind", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Split EquiBind into fpocket/p2rank/unguided x "
                         "smina/gnina/raw x clamp variants from the CSV provenance "
                         "columns (filename fallback) (default: on).")
    ap.add_argument("--best-equibind-only", action="store_true",
                    help="Keep only the single best-performing EquiBind variant "
                         "(highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%) in "
                         "all summaries and plots, relabelled 'EquiBind*'. AutoDock/DiffDock are "
                         "unaffected. Most useful with --split-equibind (default).")
    ap.add_argument("--diffdock-variant", default="diffdock",
                    choices=("diffdock", "diffdock_smina", "diffdock_gnina", "all"),
                    help="Which DiffDock variant represents 'diffdock' in the "
                         "comparison: 'diffdock' (raw docking output, default) | "
                         "'diffdock_smina' | 'diffdock_gnina'. The chosen variant "
                         "is analysed as the canonical DiffDock tool; the others "
                         "are dropped (one DiffDock entry per run). 'all' keeps all "
                         "three as separate methods — use it to write an oracle_summary "
                         "that --best-diffdock-only (in the other reports) can rank.")
    ap.add_argument("--best-diffdock-only", action="store_true",
                    help="Keep only the single best-performing DiffDock optimizer "
                         "variant (raw/smina/gnina, highest PB-Valid AND RMSD ≤ 2 Å = "
                         "oracle_pb_valid_and_rmsd2_%%) in all summaries/plots, relabelled "
                         "'DiffDock*'. Scores all three then picks the best, so it OVERRIDES "
                         "--diffdock-variant. AutoDock/EquiBind are unaffected.")
    ap.add_argument("--collapse-diffdock-variant",
                    choices=("diffdock", "diffdock_smina", "diffdock_gnina"), default=None,
                    help="When --best-diffdock-only / --collapse-plots-only collapse DiffDock "
                         "to a single 'DiffDock*' variant, FORCE it to be this one instead of "
                         "the oracle-ranked best (does not affect the per-pose cache).")
    ap.add_argument("--collapse-autodock-variant",
                    help="Pin the AutoDock slot of the collapsed presentation frame to "
                         "one Vina-scored variant (autodock | autodock_gnina | "
                         "autodock_gnina_refinement) instead of always raw Vina. "
                         "Counterpart to --collapse-diffdock-variant. Takes effect under "
                         "--collapse-plots-only or --best-variants-only. Unset keeps the "
                         "previous raw-Vina behaviour, so existing runs are unchanged.")
    ap.add_argument("--best-variants-only", action="store_true",
                    help="Convenience umbrella: enable BOTH --best-equibind-only and "
                         "--best-diffdock-only, so ordinary summaries and headline "
                         "cross-tool plots show raw AutoDock, DiffDock*, and EquiBind*. "
                         "Explicit variant diagnostics and per_pose_metrics.csv still "
                         "keep every variant for drill-down.")
    ap.add_argument("--collapse-plots-only", action="store_true",
                    help="Keep EVERY variant in explicit all-variant tables and "
                         "diagnostics, while ordinary summaries and headline cross-tool "
                         "graphs use raw AutoDock plus one DiffDock and EquiBind winner. "
                         "Writes the full "
                         "all-variants tables (oracle_summary_all_variants.csv, "
                         "top1_summary_all_variants.csv, and rmsd_vs_pbvalid_all_variants.csv "
                         "— the RMSD vs RMSD+PB-validity comparison per variant). "
                         "oracle_summary.csv is collapsed; explicitly named all-variant, "
                         "raw-vs-refined, and optimization figures remain diagnostic. "
                         "Like --best-variants-only but the full tables are also written.")
    mf.add_method_filter_args(ap)
    args = ap.parse_args()
    _set_reference_convention(args.reference_convention)
    args.cache_refused = None
    if args.reference_convention != "instance":
        print(f"reference-convention: {args.reference_convention} — every "
              "crystal-referenced metric is measured against the pose's nearest "
              "deposited ligand copy (PoseBusters convention).")

    # --best-variants-only is a convenience umbrella for the two per-family "best"
    # filters. Apply it here, before anything reads those booleans: the per-pose
    # cache signature and the DiffDock variant selection both key off
    # best_diffdock_only, so it must be set prior to _per_pose_signature /
    # _build_pose_index below.
    if args.best_variants_only:
        args.best_equibind_only = True
        args.best_diffdock_only = True
        print("best-variants-only: enabling --best-equibind-only + --best-diffdock-only "
              "(graphs/summaries show AutoDock, DiffDock*, EquiBind* only).")
    # --collapse-plots-only scores every variant (like best-variants-only, so the full
    # tables can be written) then collapses to the best per tool only for the plots. It
    # sets both best-*-only so all variants stay in df; the all-variants tables are dumped
    # just before the collapse (see below).
    if args.collapse_plots_only:
        args.best_equibind_only = True
        args.best_diffdock_only = True
        print("collapse-plots-only: all variants kept in *_all_variants.csv and "
              "variant diagnostics; ordinary summaries + headline plots collapse to "
              "AutoDock, DiffDock*, EquiBind*.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    root = Path.cwd()

    # ── Per-pose metrics: reuse the cached results if the inputs are unchanged,
    #    otherwise run the (heavy) per-pose scoring and cache it. --top-n is NOT
    #    part of the cache key, so changing it alone reuses the cache. ──
    sig = _per_pose_signature(args)
    df = _load_cached_per_pose(args, sig)

    if df is None and getattr(args, "cache_refused", None):
        # The loader already printed the ERROR (reference-convention mismatch).
        # Exit non-zero so a harness / notebook stage cannot mistake the refusal
        # for a completed run.
        raise SystemExit(2)

    if df is None and args.reuse_cache:
        print(f"ERROR: --reuse-cache was set but no usable per-pose cache exists "
              f"in {args.out_dir}.\n       Run once without --reuse-cache to build "
              f"per_pose_metrics.csv, then re-run with it. A cache written before "
              "schema 6 (no reference_convention columns) must be rebuilt with --force "
              "into a NEW --out-dir or after backing the directory up.")
        raise SystemExit(2)

    if df is None:
        print(f"Loading pose index from: {args.pb_csv}")
        idx = _build_pose_index(args.pb_csv, split_equibind=args.split_equibind,
                                diffdock_variant=args.diffdock_variant,
                                select_diffdock=not args.best_diffdock_only)
        # Restrict to the official benchmark-set ids (the complex id is the
        # 'protein' column, == '<PDBID>_<LIG>' for the benchmark staging).
        if args.ids_file:
            allowed = _load_allowed_ids(args.ids_file)
            n0 = idx.groupby(["protein", "ligand"]).ngroups
            idx = idx[idx["protein"].astype(str).isin(allowed)].copy()
            n1 = idx.groupby(["protein", "ligand"]).ngroups
            print(f"  restricted to {len(allowed)} ids from {args.ids_file.name}: "
                  f"{n0} → {n1} pairs")

        grouped: dict[tuple[str, str], list[dict]] = {}
        for rec in idx.to_dict("records"):
            grouped.setdefault((rec["protein"], rec["ligand"]), []).append(rec)
        pairs = list(grouped.items())
        if args.limit_pairs:
            pairs = pairs[:args.limit_pairs]
        print(f"  pairs to process: {len(pairs):,} (workers={args.workers})")

        work = [((p, l), recs, str(args.benchmark_dir.resolve()), root,
                 args.reference_convention)
                for (p, l), recs in pairs]

        out_records: list[dict] = []
        if args.workers > 1:
            # chunksize=1: each pair is a big, roughly equal-cost task (~hundreds
            # of poses), so per-task IPC is negligible. Handing pairs out one at a
            # time keeps every worker busy right through the tail — with chunksize>1
            # the last chunks starve most workers (e.g. 40 pairs / chunksize 4 = 10
            # chunks leaves 21 of 31 workers idle), which is what capped CPU well
            # below 100%. One-at-a-time dispatch saturates all workers to the end.
            with mp.Pool(args.workers) as pool:
                for i, batch in enumerate(
                        pool.imap_unordered(process_pair, work, chunksize=1), 1):
                    out_records.extend(batch)
                    if i % 25 == 0 or i == len(work):
                        print(f"  [{i}/{len(work)}] pairs scored "
                              f"({len(out_records):,} poses)")
        else:
            for i, w in enumerate(work, 1):
                out_records.extend(process_pair(w))
                if i % 25 == 0 or i == len(work):
                    print(f"  [{i}/{len(work)}] pairs scored")

        if not out_records:
            print("No pose records produced — check input paths.")
            return

        df = pd.DataFrame(out_records)
        per_pose_csv = _write_per_pose_cache(args, df, sig)
        print(f"\nWrote per-pose metrics  → {per_pose_csv}  ({len(df):,} rows)")

    # ── Single-point method exclusion ────────────────────────────────
    # Both branches above converge here with a complete frame, and nothing has
    # consumed it yet, so this is the only place an exclusion has to be applied
    # to reach every downstream table, figure and statistic — including df_full,
    # the *_all_variants tables and the 09f cascade, each of which enumerates
    # variants independently and would otherwise need its own edit.
    # Deliberately AFTER _write_per_pose_cache: the cache keeps every variant, so
    # the exclusion is a presentation-time choice and a later --reuse-cache run is
    # not silently poisoned by it. (The cache signature does not record these
    # flags, exactly as it does not record the collapse flags.)
    df = mf.apply_method_filter(df, "method", args, label="pose-comparison",
                                out_dir=args.out_dir)

    # Summary of what we're analysing (works for cached & freshly scored df alike).
    methods_found = sorted(df["method"].astype(str).unique())
    print(f"  poses: {len(df):,}  "
          f"pairs: {df.groupby(['protein', 'ligand']).ngroups:,}  "
          f"methods: {methods_found}")
    ranking_found = [m for m in methods_found if m in RANKING_TOOLS]
    print(f"  ranking tools (Part B): {ranking_found or '(none found)'}")
    print(f"  top-n for ranking analysis: {args.top_n}")

    # ── PB-validity among near-native poses: all generated vs the oracle pick ──
    # Computed on the FULL (pre-collapse) df so the raw/smina/gnina contrast is preserved,
    # and labelled by true variant identity (the DiffDock*/EquiBind* overrides aren't set yet).
    _w2cmp = within2_validity_comparison(df)
    _w2_stats = None
    if not _w2cmp.empty:
        _write_csv(_w2cmp, args.out_dir / "within2_validity_comparison.csv")
        # Paired within-method test of the circle↔diamond gap (oracle vs near-native
        # pool PB-validity); degrades to the test-free figure on any failure.
        try:
            _w2_stats = _stats_within2_gap(df)
        except Exception as exc:
            print(f"  WARNING: fig 19 within2 gap stats failed ({exc})")
        # Figure keeps fpocket/p2rank as separate rows but brackets them as one guided block —
        # they are handed the pocket, so warrant reading apart from the blind methods.
        plot_within2_validity_dumbbell(
            _w2cmp, args.out_dir / "19_within2_validity_dumbbell.png", stats=_w2_stats)
        print(f"  wrote near-native validity comparison → within2_validity_comparison.csv "
              f"+ 19_within2_validity_dumbbell.png ({len(_w2cmp)} variants)")

    # ── --collapse-plots-only: dump the FULL (all-variants) summary tables BEFORE the
    #    best-per-tool collapse below, so the plots use only the best variant while the
    #    tables retain every EquiBind + DiffDock variant. Includes the RMSD vs
    #    RMSD+PB-validity comparison per variant (rmsd_vs_pbvalid_all_variants.csv).
    if args.collapse_plots_only:
        _oracle_all = aggregate_oracle(df)
        _top1_all = aggregate_top1(df)
        _write_csv(_oracle_all, args.out_dir / "oracle_summary_all_variants.csv")
        _write_csv(_top1_all, args.out_dir / "top1_summary_all_variants.csv")
        _write_rmsd_vs_pbvalid_table(
            _oracle_all, _top1_all, args.out_dir / "rmsd_vs_pbvalid_all_variants.csv")
        print(f"collapse-plots-only: wrote all-variants tables for {len(_oracle_all)} "
              "variants (oracle_summary_all_variants.csv, top1_summary_all_variants.csv, "
              "rmsd_vs_pbvalid_all_variants.csv) — headline plots use the canonical "
              "three-tool presentation.")

    # Keep the FULL per-variant frame (before any best-variant collapse) so the
    # raw-vs-best optimization figure (09c) can still see every DiffDock/EquiBind
    # variant under its own method key.
    df_full = df.copy()

    # ── Pose-level validity & accuracy cascade (every variant) ───────
    # Pools each variant's whole produced-pose set and walks produced → PB-valid →
    # RMSD ≤ 2 Å → Kabsch RMSD < 1 Å → (≤ 2 Å & PB-valid & Kabsch < 1 Å), reporting
    # BOTH a pose count and the number of complexes reached at each stage. Independent
    # of the best-variant/collapse flags (runs on df_full), so it always shows all
    # 13(+) variants. The Kabsch (best-fit) threshold follows --form-ok-kabsch.
    try:
        _cascade = aggregate_pose_validity_cascade(
            df_full, kabsch_thr=args.form_ok_kabsch)
        if not _cascade.empty:
            _write_csv(_round_csv(_cascade),
                args.out_dir / "pose_validity_cascade.csv", index=False)
            _cascade_tbl = _write_pose_validity_cascade_report(
                _cascade, args.out_dir / "pose_validity_cascade_report.txt",
                kabsch_thr=args.form_ok_kabsch)
            print("\nPose-level validity & accuracy cascade (all variants):")
            if _cascade_tbl:
                print(_cascade_tbl)
            print("  wrote pose_validity_cascade.csv + pose_validity_cascade_report.txt")
    except Exception as exc:
        print(f"  WARNING: pose validity cascade failed ({exc})")

    # ── Optional: restrict the report to the single best EquiBind variant ──
    # Done after the full per-pose dump (which keeps every variant) so only the
    # summaries and plots below are filtered.
    best_eq = None
    if args.best_equibind_only:
        df, best_eq = _select_best_equibind(df)
        if best_eq:
            _LABEL_OVERRIDES[best_eq] = "EquiBind*"
            print(f"best-equibind-only: '{best_eq}' is the top EquiBind variant "
                  "by PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'EquiBind*').")
        else:
            print("best-equibind-only: no EquiBind variants present — nothing filtered.")

    # ── Optional: restrict the report to the single best DiffDock variant ──
    # Scores all three optimizer variants (kept because select_diffdock was off),
    # collapses to the best, relabelled 'DiffDock*'.
    best_dd = None
    if args.best_diffdock_only:
        df, best_dd = _select_best_diffdock(df, forced=args.collapse_diffdock_variant)
        if best_dd:
            _LABEL_OVERRIDES["diffdock"] = "DiffDock*"
            if best_dd == args.collapse_diffdock_variant:
                print(f"best-diffdock-only: '{best_dd}' was pinned by "
                      "--collapse-diffdock-variant — keeping only it as 'diffdock' "
                      "(shown as 'DiffDock*').")
            else:
                print(f"best-diffdock-only: '{best_dd}' is the top DiffDock variant by "
                      "PB-Valid AND RMSD ≤ 2 Å — keeping only it as 'diffdock' "
                      "(shown as 'DiffDock*').")
        else:
            print("best-diffdock-only: no DiffDock variants present — nothing filtered.")

    # ── Aggregation ──────────────────────────────────────────────
    # The family selectors above preserve every unrelated method by design.
    # Under the two three-tool umbrella modes, finish the presentation collapse
    # explicitly so newly added AutoDock/Vinardo optimizer variants and UniDock do
    # not leak into ordinary summaries or headline plots. df_full and, for
    # --collapse-plots-only, the *_all_variants.csv tables remain complete.
    if args.collapse_plots_only or args.best_variants_only:
        # Pin the AutoDock slot BEFORE the keep-set runs, since the keep-set only
        # recognises the canonical 'autodock' key and would otherwise drop the
        # pinned optimizer variant outright.
        df, best_ad = _select_autodock_arm(df, forced=args.collapse_autodock_variant)
        if best_ad and best_ad != "autodock":
            # Star the pinned arm, matching the DiffDock*/EquiBind* convention that
            # _select_best_diffdock and _select_best_equibind already apply, so a
            # collapsed legend names all three engines the same way. (Before
            # 2026-08-20 this spelled the arm out as "AutoDock Vina + gnina", which
            # left the AutoDock slot as the only unstarred one in the legend.) The
            # refinement arm keeps its explicit name, since it is never the reported
            # pipeline and must not be mistaken for the starred default.
            _ad_label = {
                "autodock_gnina": "AutoDock*",
                "autodock_mgltools_exh128_gnina": "AutoDock*",
                "autodock_gnina_refinement": "AutoDock Vina + gnina refine",
            }.get(best_ad, f"AutoDock ({best_ad})")
            _LABEL_OVERRIDES["autodock"] = _ad_label
            print(f"collapse-autodock-variant: '{best_ad}' was pinned — keeping only it "
                  f"as 'autodock' (shown as '{_ad_label}').")
        df = _select_presentation_tools(df, best_eq)
        mode = ("collapse-plots-only" if args.collapse_plots_only
                else "best-variants-only")
        print(f"{mode}: presentation methods: "
              f"{sorted(df['method'].astype(str).unique())}"
              f"{f' (AutoDock slot = {best_ad})' if best_ad else ''}")

    oracle_sum = aggregate_oracle(df)
    _write_csv(oracle_sum, args.out_dir / "oracle_summary.csv")
    print("\nOracle summary (all tools):")
    print(oracle_sum.to_string())

    top1_sum = aggregate_top1(df)
    _write_csv(top1_sum, args.out_dir / "top1_summary.csv")
    if not top1_sum.empty:
        print("\nTop-1 summary (ranking tools):")
        print(top1_sum.to_string())

    rank_df = aggregate_by_rank(df, args.top_n)
    _write_csv(rank_df, args.out_dir / "per_rank_metrics.csv", index=False)

    ifp_rank_df = aggregate_ifp_by_rank(df, args.top_n)
    _write_csv(ifp_rank_df, args.out_dir / "per_rank_ifp_recovery.csv", index=False)
    print("  wrote per-rank interaction-fingerprint recovery → per_rank_ifp_recovery.csv")

    oracle_pivot = _oracle_per_pair(df).pivot_table(
        index=["protein", "ligand"], columns="method", values="rmsd")
    _write_csv(oracle_pivot, args.out_dir / "oracle_rmsd_per_pair.csv")

    # ── Part A: oracle comparison (all tools) ────────────────────
    # Statistical tests for the headline oracle / success / paper figures. Each
    # is computed once (unit = per-complex) and both annotated on the figure and
    # written to the stats sidecar. Any failure degrades to the test-free figure.
    _stats_sidecar: dict = {}
    try:
        _oracle_rmsd_stats = _stats_oracle_rmsd_paired(df)
        if _oracle_rmsd_stats:
            _stats_sidecar["oracle_rmsd_across_tools"] = _oracle_rmsd_stats
    except Exception as exc:
        _oracle_rmsd_stats = None
        print(f"  WARNING: oracle-RMSD paired stats failed ({exc})")
    try:
        _succ03_stats = _stats_success_paired(df, valid=False)
        if _succ03_stats:
            _stats_sidecar["success_rmsd2_oracle_vs_top1"] = _succ03_stats
    except Exception as exc:
        _succ03_stats = None
        print(f"  WARNING: fig 03 success stats failed ({exc})")
    try:
        _succ10_stats = _stats_success_paired(df, valid=True)
        if _succ10_stats:
            _stats_sidecar["success_rmsd2_pbvalid_oracle_vs_top1"] = _succ10_stats
    except Exception as exc:
        _succ10_stats = None
        print(f"  WARNING: fig 10 success stats failed ({exc})")
    try:
        _paper11_stats = _stats_vs_paper(df, oracle_sum, top1_sum)
        if _paper11_stats:
            _stats_sidecar["vs_posebusters_paper"] = _paper11_stats
    except Exception as exc:
        _paper11_stats = None
        print(f"  WARNING: fig 11 paper stats failed ({exc})")
    # fig 19 within-2Å validity gap (computed above, on the pre-collapse df)
    if _w2_stats:
        _stats_sidecar["within2_validity_dumbbell"] = _w2_stats
    # fig 15b top-k recovery — paired McNemar at pre-specified depths k ∈ {1,5,10,15}
    # (near-native) + the nested near-vs-valid gap as an effect size. On df_full so all
    # DiffDock/EquiBind variants are present; drives the 15b footnote + its own CSV.
    try:
        _topk_stats = topk_recovery_stats(df_full)
        if _topk_stats:
            _stats_sidecar["topk_recovery"] = _topk_stats
    except Exception as exc:
        _topk_stats = None
        print(f"  WARNING: fig 15b top-k recovery stats failed ({exc})")
    # fig 09e — strict ≤2Å & PB-valid success per 09b variant, oracle vs selected
    # rank-1 (EquiBind ranked by gnina affinity). Paired proportions across variants
    # + within-variant McNemar; needs df_full so every raw/opt variant is present.
    try:
        _variants09e_stats = _stats_pb_valid_variants(
            df_full, thr=2.0, forced_dd=args.collapse_diffdock_variant, valid=True)
        if _variants09e_stats:
            _stats_sidecar["pb_valid_variants_oracle_vs_rank1"] = _variants09e_stats
    except Exception as exc:
        _variants09e_stats = None
        print(f"  WARNING: fig 09e variant success stats failed ({exc})")

    print("\nGenerating plots …")
    plot_oracle_rmsd_cdf(df, args.out_dir / "01_oracle_rmsd_cdf.png",
                         stats=_oracle_rmsd_stats)
    plot_oracle_rmsd_box(df, args.out_dir / "02_oracle_rmsd_boxplot.png",
                         stats=_oracle_rmsd_stats)
    plot_oracle_vs_top1_success(oracle_sum, top1_sum,
                                args.out_dir / "03_oracle_vs_top1_success.png",
                                stats=_succ03_stats)
    # Accuracy-vs-validity (former two-panel fig 09), now split into 09a (top-1)
    # and 09b (oracle), plus 09c contrasting each tool's raw vs best-optimised pose.
    plot_accuracy_validity_top1(top1_sum,
                                args.out_dir / "09a_accuracy_vs_validity_top1.png")
    plot_accuracy_validity_oracle(oracle_sum,
                                  args.out_dir / "09b_accuracy_vs_validity_oracle.png",
                                  df_full=df_full,
                                  forced_dd=args.collapse_diffdock_variant)
    rvb = aggregate_optimization_raw_vs_best(df_full)
    if not rvb.empty:
        # role_label carries newlines for the two-line bar labels — keep the clean
        # single-line `role` column in the CSV.
        _write_csv(rvb.drop(columns=["role_label"]),
            args.out_dir / "optimization_raw_vs_best.csv", index=False)
        plot_optimization_raw_vs_best(
            rvb, args.out_dir / "09c_optimization_raw_vs_best.png")
        print("  wrote raw-vs-best optimization comparison → "
              "optimization_raw_vs_best.csv")
        # Guard: 09c needs the raw + optimised variant of each ML tool. If a run
        # built df without the smina/gnina DiffDock variants (e.g. a plain default
        # run, no --diffdock-variant all), the smina-opt bar is silently absent.
        if not (rvb["method_key"] == "diffdock_gnina").any() \
                and (rvb["method_key"] == "diffdock").any():
            print("  WARNING: 09c is missing DiffDock's gnina-opt bar — df lacks "
                  "diffdock_gnina (re-run with --diffdock-variant all).")

    # 09f — box/whisker of the produced-pose PB-valid YIELD per variant (raw beside
    # its smina/gnina optimized variants): how much post-hoc optimization raises the
    # share of a variant's poses that pass PoseBusters. Per complex the % of produced
    # poses that are PB-valid; the box is the distribution across complexes. Uses
    # df_full so every raw/smina/gnina variant is present under its own method key.
    # Paired per-complex Wilcoxon (raw→opt, Holm across all contrasts) + Friedman
    # omnibus annotate the figure and feed the JSON sidecar.
    _yield_per = pbvalid_yield_per_complex(df_full)
    if not _yield_per.empty:
        _write_csv(_yield_per, args.out_dir / "pbvalid_yield_per_complex.csv", index=False)
        _yield_sum = aggregate_pbvalid_yield_by_variant(df_full)
        _write_csv(_yield_sum, args.out_dir / "pbvalid_yield_by_variant.csv", index=False)
        try:
            _yield_stats = _stats_pbvalid_yield(_yield_per)
            if _yield_stats:
                _stats_sidecar["pbvalid_yield_by_variant"] = _yield_stats
        except Exception as _exc:
            _yield_stats = None
            print(f"  WARNING: fig 09f PB-valid yield stats failed ({_exc})")
        plot_pbvalid_yield_boxplot(
            _yield_per, args.out_dir / "09f_pbvalid_yield_boxplot.png",
            summary=_yield_sum, stats=_yield_stats)
        # Full statistics + the (moved-off-figure) caption/notes → companion text file.
        # Guarded so a formatting hiccup can never abort the remaining figures/outputs.
        try:
            _write_pbvalid_yield_report(
                _yield_sum, _yield_stats,
                args.out_dir / "09f_pbvalid_yield_report.txt")
        except Exception as _exc:
            print(f"  WARNING: fig 09f report text write failed ({_exc})")
        print("  wrote produced-pose PB-valid yield (raw vs optimized) → "
              "pbvalid_yield_by_variant.csv (+ per-complex + 09f box/whisker + "
              "09f_pbvalid_yield_report.txt full stats/notes)")
        _opt_present = [s[3] for s in _PBVALID_YIELD_SPECS if s[2]]
        if not any((_yield_sum["method_key"] == k).any() for k in _opt_present):
            print("  WARNING: 09f has no optimized variants — df_full lacks the "
                  "smina/gnina variants (re-run with --diffdock-variant all "
                  "--split-equibind, or a --best-*-only / --collapse-plots-only flag).")

        # 09f (all variants) — the same produced-pose PB-valid yield box/whisker, but
        # showing EVERY variant present in df_full (all EquiBind pocket × refine [× clamp]
        # variants + all DiffDock variants + AutoDock), independent of the best-variant /
        # collapse flags. Coloured per variant, no raw→opt brackets (multiple raws per
        # tool). Companion table pbvalid_yield_by_variant_all.csv.
        _all_specs = _all_variant_yield_specs(df_full)
        _yield_per_all = pbvalid_yield_per_complex(df_full, _all_specs)
        if not _yield_per_all.empty:
            _yield_sum_all = aggregate_pbvalid_yield_by_variant(df_full, _all_specs)
            _write_csv(_yield_sum_all,
                args.out_dir / "pbvalid_yield_by_variant_all.csv", index=False)
            _n_variants = _yield_sum_all["method_key"].nunique()
            plot_pbvalid_yield_boxplot(
                _yield_per_all,
                args.out_dir / "09f_pbvalid_yield_boxplot_all_variants.png",
                summary=_yield_sum_all, stats=None, specs=_all_specs,
                draw_brackets=False, color_by_variant=True,
                figsize=(max(13.0, 1.15 * _n_variants + 3.0), 7.2),
                group_labels={"autodock": "AutoDock Vina",
                              "autodock_vinardo": "AutoDock Vinardo",
                              "autodock_mgltools":
                                  "AutoDock Vina MGL-lig (exhaustiveness sweep)",
                              "diffdock": "DiffDock",
                              "unidock": "Uni-Dock",
                              "unidock2": "Uni-Dock2",
                              "equibind": "EquiBind (all pocket × refine variants)"},
                # AutoDock's 100 % median label sits under a top-left legend here (tight
                # y-axis), so move the legend to the empty top-right corner (in-axes,
                # single column — the anchor y < 1.0 keeps it off the title band).
                legend_loc="upper right", legend_anchor=(0.992, 0.995),
                legend_ncol=1)
            print(f"  wrote 09f ALL-variants overview ({_n_variants} variants) → "
                  "09f_pbvalid_yield_boxplot_all_variants.png + "
                  "pbvalid_yield_by_variant_all.csv")

    # 09d — best-of-top-d at several depths (top-1, top-N, top-30), raw vs refined
    # (ranking headroom + optimization). Uses df_full so BOTH the raw and the
    # optimised run of DiffDock (smina) / EquiBind (gnina) are present under their
    # own method keys.
    r1tn_depths = sorted({1, args.top_n, 30})
    r1tn = aggregate_rank1_vs_topn(df_full, r1tn_depths)
    if not r1tn.empty:
        _write_csv(r1tn, args.out_dir / "rank1_vs_topn.csv", index=False)
        plot_rank1_vs_topn(r1tn, args.out_dir / "09d_rank1_vs_topn.png")
        # Non-bar alternative views of the same data (trajectory / slope / gap).
        plot_rank1_vs_topn_scatter(r1tn, args.out_dir / "09d_alt_scatter.png")
        plot_rank1_vs_topn_slopegraph(r1tn, args.out_dir / "09d_alt_slopegraph.png")
        plot_rank1_vs_topn_dumbbell(r1tn, args.out_dir / "09d_alt_dumbbell.png")
        print("  wrote top-d (raw vs refined) ranking comparison → rank1_vs_topn.csv "
              "(+ 09d bar / scatter / slopegraph / dumbbell views)")
        missing = [s[1] for s in _RANK1_TOPN_SPECS
                   if s[2] not in set(r1tn["method_key"])]
        if missing:
            print(f"  WARNING: 09d is missing variant(s) {missing} — df lacks those "
                  "method keys; the corresponding optimizer results are not present "
                  "in this report input.")

    # 09e — 09b's per-variant data as a fig-10 dumbbell (● oracle vs ○ selected
    # rank-1) for the strict ≤2Å & PB-valid metric, EquiBind ranked by gnina affinity.
    plot_pb_valid_variants_dumbbell(
        df_full, args.out_dir / "09e_pb_valid_variants_success_dumbbell.png",
        thr=2.0, forced_dd=args.collapse_diffdock_variant, valid=True,
        stats=_variants09e_stats)

    plot_pb_valid_success_bars(oracle_sum, top1_sum,
                               args.out_dir / "10_pb_valid_rmsd2_success_bars.png",
                               stats=_succ10_stats)
    plot_vs_posebusters_paper(oracle_sum, top1_sum,
                              args.out_dir / "11_vs_posebusters_paper.png",
                              args.out_dir / "comparison_vs_posebusters_paper.csv",
                              stats=_paper11_stats)
    # Persist the full numeric stats for the headline figures (recoverable numbers).
    try:
        with open(args.out_dir / "posebusters_pose_comparison_stats.json", "w") as _fh:
            json.dump(_stats_sidecar, _fh, indent=2, default=_json_default)
        print("  wrote statistical sidecar → posebusters_pose_comparison_stats.json")
    except Exception as exc:
        print(f"  WARNING: could not write stats sidecar ({exc})")
    plot_rmsd2_vs_pbvalid_grouped(oracle_sum, top1_sum,
                                  args.out_dir / "12_rmsd2_vs_pbvalid_grouped.png")
    plot_pbvalid_filter_influence(oracle_sum,
                                  args.out_dir / "13_pbvalid_filter_influence.png",
                                  args.out_dir / "pbvalid_filter_influence.csv")

    # ── "Twisted & turned": how the docked ligand deviates from its crystal
    #    conformation (translation / rotation / internal twisting), for the
    #    top-1 and oracle pose of each tool ───────────────────────────────
    _write_csv(aggregate_twist_turn(df, args.top_n),
        args.out_dir / "twist_turn_summary.csv", index=False)
    plot_twist_turn(df, args.out_dir / "14_twist_turn.png", top_n=args.top_n)
    print("  wrote twist/turn summary → twist_turn_summary.csv")

    # ── Form fidelity: of the poses we already call a success (≤ 2 Å AND
    #    PB-valid), how good is the internal conformation, and is the residual
    #    error placement- or form-limited? (crystal-free sets → empty, no-op) ──
    # Form fidelity, rendered under BOTH oracle selections (suffix "" = default):
    #   nearest       — RMSD-greedy oracle (matches oracle_pb_valid_and_rmsd2_%)
    #   valid-nearest — validity-constrained ceiling (nearest ≤ 2 Å PB-valid pose)
    # Each figure — 20 (4-panel overview), 20b (count vs form), 20d (panel D split
    # by family, mechanism-coloured), 20e (clustering diagnostic) — is stamped with
    # its selection in the title.
    form_selections = [
        ("", _near_native_valid_reps, _SEL_NOTE_NEAREST),
        ("__valid_ceiling", _best_valid_near_native_reps, _SEL_NOTE_VALID),
    ]
    for suffix, selector, note in form_selections:
        _write_csv(aggregate_form_fidelity(df, args.form_ok_kabsch, selector=selector),
            args.out_dir / f"form_fidelity_summary{suffix}.csv", index=False)
        plot_form_fidelity(df, args.out_dir / f"20_form_fidelity{suffix}.png",
                           args.form_ok_kabsch, selector=selector, sel_note=note)
        plot_form_vs_success_count(
            df, args.out_dir / f"20b_form_vs_success_count{suffix}.png",
            args.form_ok_kabsch, selector=selector, sel_note=note)
        plot_form_vs_placement_by_family(
            df, args.out_dir / f"20d_form_vs_placement_by_family{suffix}.png",
            args.form_ok_kabsch, selector=selector, sel_note=note)
        plot_form_placement_clustering(
            df, args.out_dir / f"20e_form_vs_placement_clustering{suffix}.png",
            selector=selector, sel_note=note)
    # Realistic top-1 counterpart to the two oracle ceilings (20d only): each tool's
    # RANK-1 pose, kept if it is itself a ≤ 2 Å PB-valid success — what the user is
    # handed first, not what the sampler could produce at any depth. Same axes /
    # mechanism colouring as the oracle 20d figures, so the three read side by side.
    plot_form_vs_placement_by_family(
        df, args.out_dir / "20d_form_vs_placement_by_family__top1.png",
        args.form_ok_kabsch, selector=_top1_valid_reps, sel_note=_SEL_NOTE_TOP1)
    # Filmstrip of 20d across ranking depth (rows = family, cols = top-1/3/5): how the
    # near-native valid cloud drifts as the ranking net widens, centroid-arrow annotated.
    plot_form_vs_placement_depth_filmstrip(
        df, args.out_dir / "20d_form_vs_placement_by_family__depth_filmstrip.png",
        args.form_ok_kabsch, rmsd_gate=NEAR_NATIVE_RMSD_A, depths=(1, 3, 5))
    # Same filmstrip with the ≤ 2 Å gate REMOVED (all PB-valid poses), away-from-site
    # outliers dropped (docked centroid > FAR_FROM_RECEPTOR_CENTROID_A from the native
    # site) so the axes stay readable while the near-to-moderate spread the gate hides
    # is exposed — the deeper ranks visibly add mis-positioned but still-valid poses.
    plot_form_vs_placement_depth_filmstrip(
        df, args.out_dir / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid.png",
        args.form_ok_kabsch, rmsd_gate=None, depths=(1, 3, 5),
        centroid_max=FAR_FROM_RECEPTOR_CENTROID_A, axis_max=8.0)
    # Honest-statistics companion to the _pbvalid filmstrip: the numbers the scatter
    # can only gesture at (differential coverage, r-coupling, paired cross-tool tests,
    # pose diversity). Tidy CSVs + a 4-panel figure; no-op for crystal-free sets.
    _fs_stats = aggregate_filmstrip_statistics(
        df, depths=(1, 3, 5), centroid_max=FAR_FROM_RECEPTOR_CENTROID_A,
        form_ok=args.form_ok_kabsch)
    if _fs_stats:
        for _key in ("per_tool_depth", "within_tool_trend", "crosstool_paired", "pose_diversity"):
            _write_csv(_fs_stats[_key],
                args.out_dir / f"filmstrip_stats__{_key}.csv", index=False)
        plot_filmstrip_statistics(
            _fs_stats,
            args.out_dir / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__stats.png",
            depths=(1, 3, 5), form_ok=args.form_ok_kabsch)
        # Mechanism-composition shift across depth (share of placement/mixed/form-
        # limited poses growing from top-1 → top-3 → top-5, per tool).
        plot_filmstrip_mechanism_share(
            _fs_stats,
            args.out_dir / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid__mechanism_share.png",
            depths=(1, 3, 5))
        print("  wrote filmstrip honest-statistics → filmstrip_stats__{per_tool_depth,"
              "within_tool_trend,crosstool_paired,pose_diversity}.csv + "
              "20d_...__depth_filmstrip_pbvalid__stats.png + __mechanism_share.png")
    # Companion reading across ranking depth (top-1 / top-5 / top-15 / top-30), each
    # panel a standalone figure, rendered under BOTH selections so they read side by
    # side: pbvalid (all PB-valid poses, RMSD gate removed — "how does form hold up as
    # the net widens?") and within2 (PB-valid AND RMSD ≤ 2 Å — the near-native subset).
    depth_modes = [
        ("_pbvalid", None,
         "Form fidelity of the PoseBusters-valid poses by ranking depth\n"
         "(RMSD ≤ 2 Å gate removed; top-1 / top-5 / top-15 / top-30)\n"
         "selection: all PB-valid poses within each method's top-d ranked"),
        ("_within2", NEAR_NATIVE_RMSD_A,
         "Form fidelity of the near-native, PoseBusters-valid poses\n"
         "by ranking depth (top-1 / top-5 / top-15 / top-30)\n"
         "selection: PB-valid AND RMSD ≤ 2 Å poses within top-d ranked"),
    ]
    for suffix, gate, sup in depth_modes:
        _write_csv(aggregate_form_fidelity_by_depth(df, args.form_ok_kabsch, rmsd_gate=gate),
            args.out_dir / f"form_fidelity_summary__rank_depth{suffix}.csv", index=False)
        plot_form_fidelity_by_depth(
            df, args.out_dir / f"20_form_fidelity__rank_depth{suffix}.png",
            args.form_ok_kabsch, rmsd_gate=gate, suptitle=_vt(sup))
    # Joint view: the two selections overlaid (lines only) so the IMPACT of the ≤ 2 Å
    # criterion is the gap between the paired curves (per method × depth).
    plot_form_fidelity_depth_gate_impact(
        df, args.out_dir / "20_form_fidelity__rank_depth_gate_impact.png",
        args.form_ok_kabsch, suptitle=_vt(
            "Impact of the RMSD ≤ 2 Å criterion on form fidelity, by ranking depth\n"
            "solid = PB-valid AND ≤ 2 Å (near-native)\n"
            "dashed = PB-valid only (RMSD gate removed)"))
    # Per-rank companion to the gate-impact view: form error vs rank for the PB-valid
    # (A) and near-native (B) cohorts stacked, plus the per-rank gate effect (C =
    # A − B), each curve labelled with the % form-correct at every top-d marker.
    # Paired per-complex tests (gate effect near-vs-far, cross-tool form omnibus,
    # per-complex rank trend, form⊥placement coupling) drive the panel annotations
    # + the JSON sidecar; degrades to the test-free figure if the stats raise.
    try:
        _gate_vr_stats = _stats_form_fidelity_gate_vs_rank(
            df, args.form_ok_kabsch, depths=(1, 5, 10, 15))
        if _gate_vr_stats:
            _stats_sidecar["form_fidelity_gate_vs_rank"] = _gate_vr_stats
    except Exception as _exc:
        _gate_vr_stats = None
        print(f"  WARNING: fig 20 gate-vs-rank stats failed ({_exc})")
    plot_form_fidelity_gate_vs_rank(
        df, args.out_dir / "20_form_fidelity__rank_depth_gate_vs_rank.png",
        args.form_ok_kabsch, depths=(1, 5, 10, 15), stats=_gate_vr_stats, suptitle=_vt(
            "Form fidelity vs rank (top-15) — impact of the RMSD ≤ 2 Å criterion\n"
            "(A) PB-valid   (B) near-native (PB-valid AND ≤ 2 Å)   (C) gate effect\n"
            "curves = median best-fit RMSD (band = IQR); labels = % form-correct · n complexes"))
    if _gate_vr_stats:
        try:
            with open(args.out_dir / "posebusters_pose_comparison_stats.json", "w") as _fh:
                json.dump(_stats_sidecar, _fh, indent=2, default=_json_default)
            print("  updated statistical sidecar with fig-20 gate-vs-rank tests")
        except Exception as _exc:
            print(f"  WARNING: could not update stats sidecar ({_exc})")
    # Head-to-head of the two selections: table (how many complexes each rule
    # keeps + how many the nearest rule needlessly drops) and figure (20c).
    _write_csv(aggregate_oracle_selection_comparison(df),
        args.out_dir / "oracle_selection_comparison.csv", index=False)
    plot_oracle_selection_comparison(
        df, args.out_dir / "20c_oracle_selection_comparison.png", args.form_ok_kabsch)
    # Didactic 3-D illustration of what the mechanism labels mean (data-independent)
    # + a PyMOL-loadable PDB (chain A = crystal, chain B = docked pose).
    plot_mechanism_examples_3d(args.out_dir / "20f_mechanism_examples_3d.png",
                               args.out_dir / "mechanism_examples.pdb")
    # Real docked poses that exemplify each mechanism label (from the data) + per-
    # region PDBs holding the crystal (original) + docked ligand for PyMOL.
    plot_mechanism_examples_real(df, args.benchmark_dir,
                                 args.out_dir / "20g_mechanism_examples_real_3d.png",
                                 root=root)
    print("  wrote form-fidelity summaries (nearest + valid_ceiling + rank_depth "
          "pbvalid/within2) + oracle_selection_comparison.csv (+ 20/20b/20c/20d/20e figures, "
          "both selections; + 20d __top1 realistic-rank-1 + __depth_filmstrip top-1/3/5 drift "
          "(gated + __depth_filmstrip_pbvalid gate-removed, centroid-trimmed to the "
          "crystal-site neighbourhood); "
          "+ 20_form_fidelity__rank_depth_{pbvalid,within2} standalone panels, "
          "top-1/5/15/30; + 20_form_fidelity__rank_depth_gate_impact joint ≤2Å-criterion view) "
          "+ 20f mechanism illustration & mechanism_examples.pdb")

    # ── Part B: ranking quality (Vina + DiffDock only) ───────────
    if not rank_df.empty:
        plot_rank_success_curve(rank_df, args.top_n,
                                args.out_dir / "04_rank_success_curve.png")
        plot_plif_recovery_by_rank(ifp_rank_df, args.top_n,
                                   args.out_dir / "04b_plif_recovery_by_rank.png")
        plot_cumulative_oracle_curve(rank_df, args.top_n,
                                     args.out_dir / "05_cumulative_oracle_curve.png")
        plot_top1_vs_oracle_scatter(df, args.out_dir / "06_top1_vs_oracle_scatter.png")
        plot_oracle_rank_histogram(df, args.top_n,
                                   args.out_dir / "07_oracle_rank_histogram.png")
        # How often is the n-th ranked pose the closest (oracle) one?
        odist = aggregate_oracle_rank_distribution(df, args.top_n)
        _write_csv(odist, args.out_dir / "oracle_rank_distribution.csv", index=False)
        plot_oracle_rank_distribution(odist, args.top_n,
                                      args.out_dir / "15_oracle_rank_distribution.png")
        print("  wrote oracle-rank distribution → oracle_rank_distribution.csv")

        # 15b — cumulative recovery within top-k, near-native (solid) vs. PB-valid
        # (dashed), raw DiffDock alongside DiffDock* (smina). Uses df_full so both
        # DiffDock variants are present; fig 15 itself ignores PB-validity.
        topk_rec = aggregate_topk_recovery(df_full, range(1, 31))
        if not topk_rec.empty:
            _write_csv(topk_rec, args.out_dir / "topk_recovery_validity.csv", index=False)
            plot_topk_recovery_validity(
                topk_rec, args.out_dir / "15b_topk_recovery_validity.png",
                stats=_topk_stats)
            print("  wrote top-k recovery (near-native vs PB-valid) → "
                  "topk_recovery_validity.csv")
            # Tidy CSV of the pre-specified paired McNemar tests (k ∈ {1,5,10,15}).
            if _topk_stats and _topk_stats.get("between_method"):
                _bm = pd.DataFrame(_topk_stats["between_method"])
                _bm["baseline_rate_%"] = (100 * _bm["a_rate"]).round(PCT_DECIMALS)
                _bm["comparison_rate_%"] = (100 * _bm["b_rate"]).round(PCT_DECIMALS)
                _bm = _bm.rename(columns={"a_wins": "baseline_only_wins",
                                          "b_wins": "comparison_only_wins",
                                          "n": "n_complexes"})
                # `family` and `family_size` are written out so a reader can see
                # which set of hypotheses each p-value was corrected against.
                _cols = ["baseline", "comparison", "k", "n_complexes",
                         "baseline_rate_%", "comparison_rate_%",
                         "baseline_only_wins", "comparison_only_wins",
                         "mcnemar_p", "p_holm", "star"]
                _cols += [c for c in ("family", "family_size") if c in _bm.columns]
                _write_csv(_bm[_cols],
                    args.out_dir / "topk_recovery_stats.csv", index=False)
                print("  wrote top-k paired McNemar tests → topk_recovery_stats.csv")

        # 21/22 — Ranking quality: how faithfully does each tool's OWN ranking
        # reproduce the oracle ordering of its poses? AutoDock/DiffDock use native
        # confidence; EquiBind (no native score) is ranked by gnina affinity. Fig 21
        # = per-position concordance (best → rank 1, 2nd-best → rank 2 …) under
        # crystal RMSD, Kabsch RMSD and PB-validity; fig 22 = per-complex Kendall
        # tau-b whisker plot over the four criteria (incl. combined valid+close).
        # Uses df_full so the gnina-ranked EquiBind variant is present.
        rq_tau = aggregate_rank_quality_tau(df_full)
        if not rq_tau.empty:
            _write_csv(rq_tau, args.out_dir / "rank_quality_tau.csv", index=False)
            _write_csv(_summarize_rank_quality_tau(rq_tau),
                args.out_dir / "rank_quality_tau_summary.csv", index=False)
            try:
                rq_stats = _stats_rank_quality(rq_tau)
            except Exception as _exc:
                rq_stats = None
                print(f"  WARNING: rank-quality cross-tool stats failed ({_exc})")
            if rq_stats:
                _rq_rows = [{"criterion": crit, **pr}
                            for crit, _cl in _RANK_QUALITY_CRITERIA
                            for pr in (rq_stats.get(crit) or {}).get("pairwise", [])]
                if _rq_rows:
                    _write_csv(pd.DataFrame(_rq_rows),
                        args.out_dir / "rank_quality_stats.csv", index=False)
            plot_rank_quality_whiskers(
                rq_tau, args.out_dir / "22_rank_quality_whiskers.png")
            diag_df, valid_df = aggregate_rank_position_concordance(
                df_full, _RANK_QUALITY_SPECS, args.top_n)
            _write_csv(diag_df, args.out_dir / "rank_position_concordance.csv", index=False)
            _write_csv(valid_df, args.out_dir / "rank_validity_by_rank.csv", index=False)
            plot_rank_position_concordance(
                diag_df, valid_df, args.top_n,
                args.out_dir / "21_rank_concordance.png")
            # 21a-d — individual rank-vs-oracle line charts (oracle rank of each
            # rank-k pose vs k; more informative than the identity-match bars) + a
            # full text report. The graphs carry NO embedded stats box — the paired
            # tests and every table live in rank_quality_report.txt.
            rvo_df = aggregate_rank_vs_oracle(df_full, _RANK_QUALITY_SPECS, _RVO_MAX_RANK)
            _write_csv(rvo_df, args.out_dir / "rank_vs_oracle.csv", index=False)
            plot_rank_vs_oracle_criterion(
                rvo_df, "crystal", args.out_dir / "21a_rank_vs_oracle_crystal.png")
            plot_rank_vs_oracle_criterion(
                rvo_df, "kabsch", args.out_dir / "21b_rank_vs_oracle_kabsch.png")
            plot_validity_by_rank_individual(
                valid_df, args.top_n, args.out_dir / "21c_validity_by_rank.png")
            plot_rank_vs_oracle_criterion(
                rvo_df, "combined", args.out_dir / "21d_rank_vs_oracle_combined.png")
            # The text report's per-rank tables (sections 3-5) span the full pose depth
            # (_RANK_REPORT_DEPTH), independent of --top-n which governs the 21/21c
            # figures + their backing CSVs. rvo already runs to _RVO_MAX_RANK; recompute
            # identity-match + validity-by-rank at the deeper depth for the report only.
            report_depth = max(int(args.top_n), _RANK_REPORT_DEPTH)
            if report_depth > args.top_n:
                diag_rep, valid_rep = aggregate_rank_position_concordance(
                    df_full, _RANK_QUALITY_SPECS, report_depth)
            else:
                diag_rep, valid_rep = diag_df, valid_df
            _write_rank_quality_report(
                args.out_dir / "rank_quality_report.txt",
                tau_df=rq_tau, summary_df=_summarize_rank_quality_tau(rq_tau),
                stats=rq_stats, rvo_df=rvo_df, diag_df=diag_rep, valid_df=valid_rep,
                mix_df=_validity_mix(df_full, _RANK_QUALITY_SPECS), top_n=report_depth)
            print("  wrote ranking-quality figures → 21_rank_concordance.png, "
                  "21a/21b/21c/21d_*.png, 22_rank_quality_whiskers.png "
                  "(+ rank_quality_tau*.csv, rank_vs_oracle.csv, "
                  "rank_position_concordance.csv, rank_validity_by_rank.csv, "
                  "rank_quality_report.txt)")

        # 15c/15d — per-rank benefit of pose optimization (smina/gnina) over the raw
        # pose, for DiffDock (confidence rank) and EquiBind (generation order): does
        # refinement help the early ranks more? 15c = per-rank curves, 15d = the +pp
        # impact summarised over rank groups 1-10 / 11-20 / 21-30. Uses df_full so raw +
        # smina + gnina variants of both tools are present.
        opt_benefit = aggregate_optimization_benefit(df_full, max_rank=30)
        if not opt_benefit.empty:
            _write_csv(opt_benefit, args.out_dir / "optimization_benefit_by_rank.csv",
                               index=False)
            # Re-annotate from the sidecar (optimization_benefit_stats.py) if present, so
            # a plain report run keeps the 15c/15d stats instead of overwriting them plain.
            _ob_stats = _load_opt_benefit_stats(args.out_dir)
            plot_optimization_benefit_by_rank(
                opt_benefit, args.out_dir / "15c_optimization_benefit_by_rank.png",
                stats=_ob_stats)
            plot_optimization_benefit_by_group(
                opt_benefit, args.out_dir / "15d_optimization_benefit_by_group.png",
                stats=_ob_stats)
            print("  wrote per-rank optimization benefit → "
                  "optimization_benefit_by_rank.csv"
                  + ("  (+15c/15d stats annotations)" if _ob_stats else ""))

        # How many of the top-N ranked poses land within 1, 1.25, 1.5 … Å?
        # EquiBind (gnina-ranked) is added from df_full since the collapsed df keeps
        # only the best EquiBind variant (which may not be the gnina one).
        eq_gnina = df_full[df_full["method"].astype(str) == "equibind_unguided_gnina"]
        _eq = eq_gnina if not eq_gnina.empty else None
        _depths = tuple(range(1, args.top_n + 1))    # item-6 pool sweep N = 1 … top_n
        within_df = aggregate_topn_within_thresholds(
            df, args.top_n, args.fine_rmsd_thresholds, eq_df=_eq)
        if not within_df.empty:
            _write_csv(within_df, args.out_dir / "topn_within_thresholds.csv", index=False)
            # Paired complex-level tests (Cochran's Q + pairwise McNemar, Holm; Wilson
            # CIs; ranking headroom) at the pre-specified 1 Å + 2 Å thresholds — drives
            # panel A's annotations and the JSON sidecar. Same (df, eq_df) as the curves.
            try:
                _topn18_stats = _stats_topn_within(
                    df, args.top_n, eq_df=_eq, pb_valid_only=False)
                if _topn18_stats:
                    _stats_sidecar["topn_within_thresholds"] = _topn18_stats
            except Exception as _exc:
                _topn18_stats = None
                print(f"  WARNING: fig 18 within-threshold stats failed ({_exc})")
            plot_topn_within_thresholds(within_df, args.top_n,
                                        args.fine_rmsd_thresholds,
                                        args.out_dir / "18_topn_within_thresholds.png",
                                        stats=_topn18_stats,
                                        csv_name="topn_within_thresholds.csv")
            thr_str = ", ".join(f"{t:g}" for t in args.fine_rmsd_thresholds)
            print(f"\nTop-ranked poses within RMSD thresholds ({thr_str} Å):")
            print(within_df.pivot_table(index="method", columns="rmsd_threshold_A",
                                        values="top1_within_%").to_string())
            print("  wrote within-threshold sweep → topn_within_thresholds.csv")

        # 18 PB-valid variant — same sweep but a pose only counts as a hit if it
        # is BOTH within t Å AND PoseBusters-valid (combined near-native + valid).
        within_pbv = aggregate_topn_within_thresholds(
            df, args.top_n, args.fine_rmsd_thresholds, eq_df=_eq, pb_valid_only=True)
        if not within_pbv.empty:
            _write_csv(within_pbv, args.out_dir / "topn_within_thresholds_pbvalid.csv",
                              index=False)
            try:
                _topn18v_stats = _stats_topn_within(
                    df, args.top_n, eq_df=_eq, pb_valid_only=True)
                if _topn18v_stats:
                    _stats_sidecar["topn_within_thresholds_pbvalid"] = _topn18v_stats
            except Exception as _exc:
                _topn18v_stats = None
                print(f"  WARNING: fig 18 PB-valid within-threshold stats failed ({_exc})")
            plot_topn_within_thresholds(
                within_pbv, args.top_n, args.fine_rmsd_thresholds,
                args.out_dir / "18_topn_within_thresholds_pbvalid.png",
                pb_valid_only=True, stats=_topn18v_stats,
                csv_name="topn_within_thresholds_pbvalid.csv")
            print("  wrote PB-valid within-threshold sweep → "
                  "topn_within_thresholds_pbvalid.csv")
            # Item 6 companion — raw→gnina-refined benefit as the top-N pool deepens,
            # on the as-placed RMSD at the 2 Å success line (same metric as this figure).
            try:
                _rg_rmsd = aggregate_refinement_gain_vs_depth(
                    df_full, thr=2.0, depths=_depths, rmsd_col="rmsd", pb_valid_only=True)
                if not _rg_rmsd.empty:
                    _write_csv(_rg_rmsd, args.out_dir / "refinement_gain_vs_depth_pbvalid.csv",
                                    index=False)
                    plot_refinement_gain_vs_depth(
                        _rg_rmsd, args.out_dir / "18_refinement_gain_vs_depth_pbvalid.png",
                        args.top_n, thr=2.0, metric_label="RMSD", pb_valid_only=True,
                        report_path=args.out_dir
                        / "18_refinement_gain_vs_depth_pbvalid_report.txt")
                    print("  wrote raw-vs-refined gain-vs-depth (RMSD) → "
                          "18_refinement_gain_vs_depth_pbvalid.png")
            except Exception as _exc:
                print(f"  WARNING: fig 18 RMSD refinement-gain companion failed ({_exc})")

        # 18 Kabsch variant — the SAME sweep, but distance is the best-fit (Kabsch)
        # RMSD (``bestfit_rmsd`` = PoseBusters' pb_kabsch_rmsd, else our own Kabsch fit)
        # instead of the as-placed RMSD. Kabsch superposes each pose onto the crystal
        # first, so the threshold measures conformer/shape fidelity ONLY (ignores where
        # the pose was docked). Plain + PB-valid, each → _complex and _pose panels.
        # Both Kabsch variants use the HALVED grid — half of every --fine-rmsd-thresholds
        # value (0 → 2.5 Å when the standard grid is 0 → 5 Å) — since superposition
        # collapses the distances toward zero, making form fidelity a sub-Ångström
        # question. The PB-valid variant additionally uses a 1 Å emphasised success line
        # (and 0.5/1 Å paired tests); the plain variant keeps the 2 Å line.
        _KABSCH_XLABEL = "Best-fit (Kabsch) RMSD threshold vs crystal ligand (Å)"
        _KABSCH_NOTE = ("Kabsch RMSD = heavy-atom RMSD after optimal superposition onto "
                        "the crystal ligand (conformer/shape error only; ignores placement).")
        _KABSCH_METRIC_WORD = "within t Å (best-fit/Kabsch RMSD)"
        _kab_half = tuple(round(t / 2.0, 4) for t in args.fine_rmsd_thresholds)
        for (_pbv, _csv_name, _png_name, _sidecar_key,
             _kab_thr, _kab_line, _kab_test) in (
                (False, "topn_within_thresholds_kabsch.csv",
                 "18_topn_within_thresholds_kabsch.png",
                 "topn_within_thresholds_kabsch",
                 _kab_half, 2.0, (1.0, 2.0)),
                (True, "topn_within_thresholds_kabsch_pbvalid.csv",
                 "18_topn_within_thresholds_kabsch_pbvalid.png",
                 "topn_within_thresholds_kabsch_pbvalid",
                 _kab_half, 1.0, (0.5, 1.0))):
            within_kab = aggregate_topn_within_thresholds(
                df, args.top_n, _kab_thr, eq_df=_eq,
                pb_valid_only=_pbv, rmsd_col="bestfit_rmsd")
            if within_kab.empty:
                continue
            _write_csv(within_kab, args.out_dir / _csv_name, index=False)
            try:
                _kab_stats = _stats_topn_within(
                    df, args.top_n, eq_df=_eq, pb_valid_only=_pbv,
                    test_thresholds=_kab_test, rmsd_col="bestfit_rmsd",
                    figure_base="18_topn_within_thresholds_kabsch",
                    metric_word=_KABSCH_METRIC_WORD)
                if _kab_stats:
                    _stats_sidecar[_sidecar_key] = _kab_stats
            except Exception as _exc:
                _kab_stats = None
                print(f"  WARNING: fig 18 Kabsch within-threshold stats failed ({_exc})")
            plot_topn_within_thresholds(
                within_kab, args.top_n, _kab_thr,
                args.out_dir / _png_name, pb_valid_only=_pbv, stats=_kab_stats,
                x_label=_KABSCH_XLABEL, metric_note=_KABSCH_NOTE,
                success_line=_kab_line, csv_name=_csv_name)
            print(f"  wrote Kabsch within-threshold sweep → {_csv_name}")
            # Item 6 companion for the Kabsch PB-valid figure — raw→refined benefit vs
            # pool depth on the best-fit (Kabsch) RMSD at the 1 Å form-fidelity line.
            if _pbv:
                try:
                    _rg_kab = aggregate_refinement_gain_vs_depth(
                        df_full, thr=1.0, depths=_depths, rmsd_col="bestfit_rmsd",
                        pb_valid_only=True)
                    if not _rg_kab.empty:
                        _write_csv(_rg_kab,
                            args.out_dir / "refinement_gain_vs_depth_kabsch_pbvalid.csv",
                            index=False)
                        plot_refinement_gain_vs_depth(
                            _rg_kab,
                            args.out_dir / "18_refinement_gain_vs_depth_kabsch_pbvalid.png",
                            args.top_n, thr=1.0, metric_label="Kabsch RMSD",
                            pb_valid_only=True,
                            report_path=args.out_dir
                            / "18_refinement_gain_vs_depth_kabsch_pbvalid_report.txt")
                        print("  wrote raw-vs-refined gain-vs-depth (Kabsch) → "
                              "18_refinement_gain_vs_depth_kabsch_pbvalid.png")
                except Exception as _exc:
                    print(f"  WARNING: fig 18 Kabsch refinement-gain companion failed "
                          f"({_exc})")

        # Multi-depth version (per request): best-of-top-1 / top-15 / top-30 curves, with
        # the value at the emphasised threshold called out on each curve. Colour = tool,
        # line style = pool depth. The full sweep + paired stats go to the companion CSV +
        # report text (NO on-figure data table, per request).
        _md_specs = (
            ("RMSD ≤ 2 Å", args.fine_rmsd_thresholds, "rmsd", 2.0,
             "RMSD threshold vs crystal ligand (Å)",
             "RMSD = symmetry-corrected heavy-atom, no superposition.",
             "18_topn_within_thresholds_pbvalid_depths.png",
             "topn_within_thresholds_pbvalid_depths.csv"),
            ("Kabsch < 1 Å", _kab_half, "bestfit_rmsd", 1.0, _KABSCH_XLABEL, _KABSCH_NOTE,
             "18_topn_within_thresholds_kabsch_pbvalid_depths.png",
             "topn_within_thresholds_kabsch_pbvalid_depths.csv"),
        )
        for (_mlab, _md_thr, _md_col, _md_line, _md_xlab, _md_note,
             _md_png, _md_csv) in _md_specs:
            try:
                _md = aggregate_within_thresholds_by_depth(
                    df, _MULTIDEPTH_POOL, _md_thr, eq_df=_eq,
                    pb_valid_only=True, rmsd_col=_md_col)
            except Exception as _exc:
                print(f"  WARNING: fig 18 multi-depth ({_md_png}) aggregate failed ({_exc})")
                continue
            if _md.empty:
                continue
            _md_stats = None
            try:
                _md_stats = _stats_within_by_depth(
                    df, _MULTIDEPTH_POOL, _md_line, eq_df=_eq,
                    pb_valid_only=True, rmsd_col=_md_col, figure_name=_md_png,
                    metric_label=("Kabsch RMSD" if _md_col == "bestfit_rmsd" else "RMSD"))
                if _md_stats:
                    _skey = ("within_by_depth_"
                             + ("kabsch_" if _md_col == "bestfit_rmsd" else "") + "pbvalid")
                    _stats_sidecar[_skey] = _md_stats
            except Exception as _e2:
                print(f"  WARNING: fig 18 multi-depth stats ({_md_png}) failed ({_e2})")
            try:
                _md_out = args.out_dir / _md_png
                plot_within_thresholds_by_depth(
                    _md, _MULTIDEPTH_POOL, _md_thr, _md_out, pb_valid_only=True,
                    x_label=_md_xlab, metric_note=_md_note, success_line=_md_line,
                    csv_name=_md_csv, stats=_md_stats,
                    report_path=_md_out.parent / (_md_out.stem + "_report.txt"))
                print("  wrote multi-depth (top-1/15/30) within-threshold → " + _md_png)
            except Exception as _exc:
                print(f"  WARNING: fig 18 multi-depth ({_md_png}) failed ({_exc})")

        # Re-persist the sidecar: the headline dump ran before fig 18's within-threshold
        # tests existed, and the top-30 re-persist below is gated on a deeper block that
        # need not run — fold the fig-18 payloads in now (idempotent superset write).
        if any(k.startswith(("topn_within_thresholds", "within_by_depth"))
               for k in _stats_sidecar):
            try:
                with open(args.out_dir / "posebusters_pose_comparison_stats.json",
                          "w") as _fh:
                    json.dump(_stats_sidecar, _fh, indent=2, default=_json_default)
                print("  updated statistical sidecar with fig-18 within-threshold tests")
            except Exception as exc:
                print(f"  WARNING: could not update stats sidecar ({exc})")
    else:
        print("  Skipping ranking plots (no ranking tool data found).")

    # ── Cross-tool ───────────────────────────────────────────────
    plot_cross_tool_pairwise(df, args.out_dir / "08_cross_tool_oracle_pairwise.png")

    # ── Pocket localization (ranking tools): do the top-N ranked poses target
    #    the experimentally validated (crystal) pocket, and is the oracle pose in
    #    a different pocket than that ranked set? ──────────────────
    pocket_sum = aggregate_pocket_localization(df, args.top_n, args.pocket_cutoff)
    pocket_by_rank = aggregate_pocket_localization_by_rank(
        df, args.top_n, args.pocket_cutoff)
    if not pocket_sum.empty:
        _write_csv(pocket_sum, args.out_dir / "pocket_localization.csv")
        if not pocket_by_rank.empty:
            _write_csv(pocket_by_rank, args.out_dir / "pocket_localization_by_rank.csv",
                                  index=False)
        plot_pocket_targeting_by_rank(pocket_by_rank, args.top_n, args.pocket_cutoff,
                                      args.out_dir / "16a_pocket_targeting_by_rank.png")
        plot_pocket_localization_summary(pocket_sum, args.top_n, args.pocket_cutoff,
                                         args.out_dir / "16b_pocket_localization_summary.png")
        print("\nPocket localization (ranking tools, "
              f"in-pocket ≤ {args.pocket_cutoff:g} Å):")
        print(pocket_sum[["pct_topN_in_validated_pocket",
                          "pct_oracle_diff_pocket_vs_topN"]].to_string())
        print("  wrote pocket localization → pocket_localization.csv")
    else:
        print("  Skipping pocket localization (no ranking-tool data).")

    # ── PoseBusters test-failure waterfall (paper-style cascade), per method ──
    pb_tests = _load_pb_test_table(args.pb_csv, split_equibind=args.split_equibind,
                                   diffdock_variant=args.diffdock_variant,
                                   select_diffdock=not args.best_diffdock_only)
    if pb_tests is not None:
        test_table, test_cols = pb_tests
        cascades = aggregate_pb_waterfall(df, test_table, test_cols, args.top_n)
        if cascades:
            _waterfall_to_csv(cascades, args.out_dir / "pb_test_waterfall.csv")
            plot_pb_waterfall(cascades, args.out_dir / "17_pb_test_waterfall.png")
            print("  wrote PoseBusters test waterfall → pb_test_waterfall.csv")

            # Same cascade, but the representative pose is the best of the tool's
            # top-N ranked poses (ranking tools only) — 17b — plus a top-1 vs
            # top-N recovery comparison — 17c.
            cascades_topn = aggregate_pb_waterfall(
                df, test_table, test_cols, args.top_n, rank_selection="topn")
            if cascades_topn:
                _waterfall_to_csv(cascades_topn,
                                  args.out_dir / "pb_test_waterfall_topn.csv")
                plot_pb_waterfall(
                    cascades_topn, args.out_dir / "17b_pb_test_waterfall_topn.png",
                    sel_desc=f"best of top-{args.top_n} ranked poses per method")
                print("  wrote top-N PoseBusters test waterfall → "
                      "pb_test_waterfall_topn.csv")
                # 17c — ranking headroom (rank-1 vs best-of-top-N), now with a
                # gnina-ranked EquiBind panel (EquiBind has no native rank, so its
                # poses are ordered by gnina affinity). Built on shallow copies so
                # 17/17b's cascades stay untouched; EquiBind gnina data comes from
                # df_full (the collapsed df no longer carries the gnina variant).
                c1, cn = dict(cascades), dict(cascades_topn)
                has_eq = _inject_equibind_gnina_headroom(
                    c1, cn, df_full, test_table, test_cols, args.top_n)
                plot_pb_waterfall_top1_vs_topn(
                    c1, cn, args.top_n,
                    args.out_dir / "17c_pb_waterfall_top1_vs_topn.png",
                    csv_out=args.out_dir / "pb_waterfall_top1_vs_topn.csv",
                    extra=("equibind_gnina",) if has_eq else ())
                print("  wrote top-1 vs top-N recovery → "
                      "pb_waterfall_top1_vs_topn.csv")

                # 17c2 — same ranking-headroom comparison, but against the best of
                # the top-30 poses (a near-oracle ceiling, since the tools emit
                # ~30 poses per complex). Top-1 is rank-depth-independent, so it is
                # reused; only the best-of-top-N depth changes.
                TOP30 = 30
                cn30 = aggregate_pb_waterfall(
                    df, test_table, test_cols, TOP30, rank_selection="topn")
                if cn30:
                    c1_30, cn_30 = dict(cascades), dict(cn30)
                    has_eq30 = _inject_equibind_gnina_headroom(
                        c1_30, cn_30, df_full, test_table, test_cols, TOP30)
                    ex30 = ("equibind_gnina",) if has_eq30 else ()
                    plot_pb_waterfall_top1_vs_topn(
                        c1_30, cn_30, TOP30,
                        args.out_dir / "17c2_pb_waterfall_top1_vs_top30.png",
                        csv_out=args.out_dir / "pb_waterfall_top1_vs_top30.csv",
                        extra=ex30)
                    print("  wrote top-1 vs top-30 recovery → "
                          "17c2_pb_waterfall_top1_vs_top30.png")

                    # 17f — the same top-30 recovery, but expressed per cascade
                    # step in counts: how many MORE complexes survive each test
                    # with best-of-top-30 vs rank-1, colour-split into the RMSD ≤ 2 Å
                    # placement filter and the other physical-validity tests.
                    try:
                        stats17f = _stats_top30_recovery(
                            c1_30, cn_30, ("autodock", "diffdock", *ex30))
                    except Exception as exc:
                        stats17f = None
                        print(f"  WARNING: 17f recovery stats failed ({exc})")
                    if stats17f:
                        _stats_sidecar["top30_recovery_refined"] = stats17f
                    plot_pb_top30_recovery_by_test(
                        c1_30, cn_30, TOP30,
                        args.out_dir / "17f_pb_top30_recovery_by_test.png",
                        csv_out=args.out_dir / "pb_top30_recovery_by_test.csv",
                        extra=ex30, stats=stats17f)
                    print("  wrote top-30 per-test recovery → "
                          "17f_pb_top30_recovery_by_test.png")

                    # 17c3 / 17g — the SAME top-30 analysis for the UN-REFINED
                    # variants: AutoDock (physics, unchanged), raw DiffDock, and raw
                    # EquiBind geometry (ranked by gnina score, since raw poses have
                    # no score). Lets raw be compared directly against the refined
                    # 17c2/17f above.
                    c1_raw, cn_raw = _build_raw_variant_cascades(
                        df_full, test_table, test_cols, TOP30)
                    raw_methods = tuple(
                        m for m in ("autodock", "diffdock_raw", "equibind_raw")
                        if m in c1_raw and m in cn_raw)
                    if raw_methods:
                        plot_pb_waterfall_top1_vs_topn(
                            c1_raw, cn_raw, TOP30,
                            args.out_dir / "17c3_pb_waterfall_top1_vs_top30_raw.png",
                            csv_out=args.out_dir / "pb_waterfall_top1_vs_top30_raw.csv",
                            methods=raw_methods)
                        try:
                            stats17g = _stats_top30_recovery(
                                c1_raw, cn_raw, raw_methods)
                        except Exception as exc:
                            stats17g = None
                            print(f"  WARNING: 17g recovery stats failed ({exc})")
                        if stats17g:
                            _stats_sidecar["top30_recovery_raw"] = stats17g
                        plot_pb_top30_recovery_by_test(
                            c1_raw, cn_raw, TOP30,
                            args.out_dir / "17g_pb_top30_recovery_by_test_raw.png",
                            csv_out=args.out_dir / "pb_top30_recovery_by_test_raw.csv",
                            methods=raw_methods, stats=stats17g)
                        print("  wrote un-refined (raw) top-30 recovery → "
                              "17c3_pb_waterfall_top1_vs_top30_raw.png, "
                              "17g_pb_top30_recovery_by_test_raw.png")

                        # 17h — raw AND refined in ONE figure: AutoDock (physics,
                        # once) then each tool's raw panel directly above its
                        # refined panel (DiffDock: the collapsed best_dd variant —
                        # smina/gnina per --collapse-diffdock-variant; EquiBind:
                        # gnina), so the
                        # refinement benefit reads off panel-to-panel. Rendered at TWO
                        # pool depths — best-of-top-{args.top_n} and best-of-top-30 vs
                        # rank-1 — each as its own PNG + per-step CSV. The rank-1 (top-1)
                        # cascade is depth-independent, so ALL THREE depths (top-1 ·
                        # top-{args.top_n} · top-30) are then written to ONE combined
                        # companion text report (full paired McNemar/Wilson stats — a
                        # separate Holm family per figure — plus a ';'-separated per-stage
                        # attrition table with one selection block per depth).
                        _recovery_blocks = []
                        for _depth in sorted({args.top_n, TOP30}):
                            _res = _emit_pb_recovery_raw_vs_refined(
                                _depth, df=df, df_full=df_full,
                                test_table=test_table, test_cols=test_cols,
                                cascades=cascades, out_dir=args.out_dir,
                                stats_sidecar=_stats_sidecar, dd_variant=best_dd)
                            if _res:
                                _recovery_blocks.append(_res)
                        if _recovery_blocks:
                            # top-1 cascades are identical across depths; take the first.
                            _rep = args.out_dir / \
                                "17h_pb_recovery_raw_vs_refined_report.txt"
                            try:
                                _write_pb_recovery_report(
                                    _rep,
                                    figure_names=[b["png"] for b in _recovery_blocks],
                                    csv_names=[b["csv"] for b in _recovery_blocks],
                                    methods=_recovery_blocks[0]["methods"],
                                    labels=_recovery_blocks[0]["labels"],
                                    cascades_top1=_recovery_blocks[0]["c1_all"],
                                    depth_blocks=[
                                        {"depth": b["depth"],
                                         "cascades_topn": b["cn_all"],
                                         "stats": b["stats"]}
                                        for b in _recovery_blocks])
                                print("  wrote combined top-1/top-"
                                      f"{args.top_n}/top-{TOP30} recovery report → "
                                      f"{_rep.name}")
                            except Exception as exc:
                                print(f"  WARNING: combined 17h report failed ({exc})")

            # Re-persist the sidecar: the headline dump above ran before the
            # top-30 recovery tests existed, so fold their McNemar/Wilson payloads
            # in now (idempotent superset write; degrades silently on failure).
            if any(k.startswith("top30_recovery") for k in _stats_sidecar):
                try:
                    with open(args.out_dir / "posebusters_pose_comparison_stats.json",
                              "w") as _fh:
                        json.dump(_stats_sidecar, _fh, indent=2, default=_json_default)
                    print("  updated statistical sidecar with top-30 recovery tests")
                except Exception as exc:
                    print(f"  WARNING: could not update stats sidecar ({exc})")

            # Raw pose vs best refined variant, top-1 per method — 17d. Uses
            # df_full (every optimiser/refinement variant under its own method key,
            # before the best-variant collapse) so raw and refined poses can be
            # cascaded side by side. The DiffDock refiner follows
            # --collapse-diffdock-variant (so 17d/17e/17i agree with the collapsed
            # report); EquiBind, having no native rank, is ranked by its gnina score.
            _dd_best_key = _resolve_dd_best_variant(
                df_full, args.collapse_diffdock_variant)
            rvb = aggregate_pb_waterfall_raw_vs_best(
                df_full, test_table, test_cols, dd_best=_dd_best_key)
            if rvb:
                _waterfall_to_csv(rvb, args.out_dir / "pb_test_waterfall_raw_vs_best.csv")
                plot_pb_waterfall_raw_vs_best(
                    rvb, args.out_dir / "17d_pb_test_waterfall_raw_vs_best.png",
                    csv_out=args.out_dir / "pb_waterfall_raw_vs_best_summary.csv")
                print("  wrote raw-vs-best PoseBusters waterfall → "
                      "17d_pb_test_waterfall_raw_vs_best.png, "
                      "pb_test_waterfall_raw_vs_best.csv")

                # Same comparison, but WITHOUT the RMSD ≤ 2 Å gate — every ranked-1
                # pose faces the PoseBusters tests regardless of placement, so the
                # cascade shows where poses fail on physical-validity grounds
                # alone (the RMSD filter otherwise removes most poses first) — 17e.
                rvb_pb = aggregate_pb_waterfall_raw_vs_best(
                    df_full, test_table, test_cols, gate_rmsd=False,
                    dd_best=_dd_best_key)
                if rvb_pb:
                    _waterfall_to_csv(
                        rvb_pb,
                        args.out_dir / "pb_test_waterfall_raw_vs_best_validity.csv")
                    plot_pb_waterfall_raw_vs_best(
                        rvb_pb,
                        args.out_dir / "17e_pb_test_waterfall_raw_vs_best_validity.png",
                        csv_out=args.out_dir / "pb_waterfall_raw_vs_best_validity_summary.csv",
                        gate_rmsd=False)
                    print("  wrote raw-vs-best PB-validity waterfall (no RMSD gate) "
                          "→ 17e_pb_test_waterfall_raw_vs_best_validity.png")

                # Rank-depth 1→top_n ADDED-pose validity — 17i. The poses newly
                # considered when the rank cut widens from 1 to top_n (ranking
                # positions 2..top_n), WITHOUT the RMSD ≤ 2 Å gate: how many poses /
                # complexes the deeper cut adds, and which PoseBusters tests those
                # added poses fail (waterfall + independent per-test tally).
                added = aggregate_pb_rank_depth_added_pool(
                    df_full, test_table, test_cols, args.top_n,
                    dd_best=_dd_best_key)
                if added:
                    plot_pb_rank_depth_added_waterfall(
                        added, args.top_n,
                        args.out_dir / "17i_pb_rank_depth_added_pose_waterfall.png",
                        csv_out=args.out_dir / "pb_rank_depth_added_pose_waterfall.csv",
                        summary_csv=args.out_dir
                        / "pb_rank_depth_added_pose_summary.csv")
                    plot_pb_rank_depth_fail_by_test(
                        added, args.top_n,
                        args.out_dir / "17i2_pb_rank_depth_added_pose_failbytest.png",
                        csv_out=args.out_dir
                        / "pb_rank_depth_added_pose_failbytest.csv")
                    print(f"  wrote rank-depth 1→{args.top_n} added-pose validity → "
                          "17i_pb_rank_depth_added_pose_waterfall.png (+ 17i2 "
                          "fail-by-test; pb_rank_depth_added_pose_"
                          "{waterfall,summary,failbytest}.csv)")
            else:
                print("  Skipping raw-vs-best waterfall (need diffdock raw + smina "
                      "and equibind raw + gnina variants; run with "
                      "--diffdock-variant all).")
        else:
            print("  Skipping PB waterfall (no representative poses to cascade).")
    else:
        print("  Skipping PB waterfall (raw CSV lacks per-test columns).")

    # ── PB-valid variants ────────────────────────────────────────
    # Re-emit every RMSD figure (01–08), the twist/turn figure (14) and their summary
    # CSVs using ONLY poses that pass all canonical PoseBusters tests (pb_valid). Each
    # graph then describes
    # the physically-valid poses alone; the denominator becomes "pairs for which
    # the method produced ≥1 valid pose" (attrition reported below + per-method in
    # pb_valid_attrition.csv). The all-complexes "valid & RMSD≤2 Å" headline stays
    # in 09a/09b_accuracy_vs_validity_*.png and the oracle/top1 summary CSVs.
    global _PLOT_TITLE_SUFFIX
    df_valid = df[df["pb_valid"]].copy()
    n_valid = len(df_valid)
    print(f"\nPB-valid poses: {n_valid:,}/{len(df):,} "
          f"({100 * n_valid / len(df):.1f}%) pass all canonical PoseBusters tests")
    if df_valid.empty:
        print("  No PB-valid poses — skipping PB-valid plot variants.")
    else:
        valid_dir = args.out_dir / "pb_valid"
        valid_dir.mkdir(parents=True, exist_ok=True)

        # Per-method validity attrition — makes the changed denominator explicit.
        attr_rows = []
        for method, g in df.groupby("method"):
            pt = g.groupby(["protein", "ligand"]).ngroups
            gv = g[g["pb_valid"]]
            pv = gv.groupby(["protein", "ligand"]).ngroups if len(gv) else 0
            attr_rows.append({
                "method": method, "pairs_total": pt, "pairs_with_valid_pose": pv,
                "pairs_with_valid_pose_%": round(100 * pv / pt, PCT_DECIMALS) if pt else 0.0,
                "poses_total": len(g), "poses_valid": int(g["pb_valid"].sum()),
                "poses_valid_%": round(100 * g["pb_valid"].mean(), PCT_DECIMALS) if len(g) else 0.0,
            })
        attr = pd.DataFrame(attr_rows).set_index("method")
        _write_csv(attr, valid_dir / "pb_valid_attrition.csv")
        print("\nPB-valid attrition by method "
              "(denominator of the pb_valid/ figures = pairs_with_valid_pose):")
        print(attr.to_string())

        # Headings are intentionally NOT suffixed with a "PB-valid only" tag —
        # the pb_valid/ figures are distinguished by their output folder and by
        # their (lower) per-method n= annotations, not by the title.
        _PLOT_TITLE_SUFFIX = ""
        try:
            oracle_sum_v = aggregate_oracle(df_valid)
            _write_csv(oracle_sum_v, valid_dir / "oracle_summary.csv")
            has_ranking_valid = bool(df_valid["method"].isin(RANKING_TOOLS).any())
            top1_sum_v = aggregate_top1(df_valid) if has_ranking_valid else pd.DataFrame()
            if not top1_sum_v.empty:
                _write_csv(top1_sum_v, valid_dir / "top1_summary.csv")
            rank_df_v = aggregate_by_rank(df_valid, args.top_n)
            ifp_rank_df_v = aggregate_ifp_by_rank(df_valid, args.top_n)
            if not rank_df_v.empty:
                _write_csv(rank_df_v, valid_dir / "per_rank_metrics.csv", index=False)
                _write_csv(ifp_rank_df_v, valid_dir / "per_rank_ifp_recovery.csv", index=False)

            print("\nGenerating PB-valid plot variants …")
            plot_oracle_rmsd_cdf(df_valid, valid_dir / "01_oracle_rmsd_cdf.png")
            plot_oracle_rmsd_box(df_valid, valid_dir / "02_oracle_rmsd_boxplot.png")
            plot_oracle_vs_top1_success(oracle_sum_v, top1_sum_v,
                                        valid_dir / "03_oracle_vs_top1_success.png")
            if not rank_df_v.empty:
                plot_rank_success_curve(rank_df_v, args.top_n,
                                        valid_dir / "04_rank_success_curve.png")
                plot_plif_recovery_by_rank(ifp_rank_df_v, args.top_n,
                                           valid_dir / "04b_plif_recovery_by_rank.png")
                plot_cumulative_oracle_curve(rank_df_v, args.top_n,
                                             valid_dir / "05_cumulative_oracle_curve.png")
                plot_top1_vs_oracle_scatter(df_valid,
                                            valid_dir / "06_top1_vs_oracle_scatter.png")
                plot_oracle_rank_histogram(df_valid, args.top_n,
                                           valid_dir / "07_oracle_rank_histogram.png")
            else:
                print("  Skipping PB-valid ranking plots (no valid ranking-tool poses).")
            plot_cross_tool_pairwise(df_valid,
                                     valid_dir / "08_cross_tool_oracle_pairwise.png")
            # Twist/turn decomposition restricted to PB-valid poses. The "twisted"
            # panels (strain energy, torsion Δ, TFD) are exactly the measures that
            # blow up for physically-invalid conformations, so a valid-only view is
            # worth having beside the unfiltered top-level 14_twist_turn.png. Ranking
            # is unchanged; because df_valid drops invalid poses first, "top-1" here
            # is each tool's highest-ranked *valid* pose, not necessarily its rank-1.
            _write_csv(aggregate_twist_turn(df_valid, args.top_n),
                valid_dir / "twist_turn_summary.csv", index=False)
            plot_twist_turn(df_valid, valid_dir / "14_twist_turn.png", top_n=args.top_n,
                            pb_valid_only=True)
        finally:
            _PLOT_TITLE_SUFFIX = ""
        print(f"\nPB-valid figures + CSVs → {valid_dir.resolve()}")

    print(f"\nAll figures + CSVs → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
