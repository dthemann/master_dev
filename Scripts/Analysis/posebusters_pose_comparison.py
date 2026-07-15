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
    variants) down to the single best-performing one each (highest PB-Valid AND
    RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%) so every summary/plot compares just AutoDock, one
    DiffDock curve and one EquiBind curve. The retained variants are relabelled
    "EquiBind*" / "DiffDock*" in every legend/axis. --best-variants-only is a
    convenience umbrella that turns on both. All off by default; AutoDock is never
    touched and per_pose_metrics.csv always carries every variant for drill-down.

Metrics computed per pose:
    * Symmetry-corrected heavy-atom RMSD vs. crystal ligand (no superposition)
    * Centroid distance to crystal ligand (translation — how far it "moved")
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
    17h_pb_top30_recovery_raw_vs_refined.png — 17f + 17g in ONE figure: AutoDock once, then
                                       each tool's raw panel directly above its refined
                                       panel (DiffDock raw/smina*, EquiBind raw/gnina), so the
                                       refinement benefit reads off panel-to-panel
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

Per-pose cache (so re-runs that only change --top-n stay cheap)
    The heavy per-pose scoring is written to per_pose_metrics.csv alongside a
    per_pose_metrics.manifest.json fingerprint of the inputs it depends on (PB CSV
    signature, --ids-file, --split-equibind, --limit-pairs — NOT --top-n, which only
    drives the cheap ranking aggregation). On the next run, if that fingerprint still
    matches, the cached metrics are reused and the per-pose scoring is skipped; pass
    --force to recompute regardless. Changing --top-n alone therefore reuses the cache.

----------------------------------------------------------------------
Quick usage
----------------------------------------------------------------------
    python Scripts/Analysis/posebusters_pose_comparison.py

    python Scripts/Analysis/posebusters_pose_comparison.py --limit-pairs 20

    # Re-run a different ranking depth — reuses the cached per-pose metrics:
    python Scripts/Analysis/posebusters_pose_comparison.py --top-n 10

    # Force a full recompute of the per-pose metrics:
    python Scripts/Analysis/posebusters_pose_comparison.py --force

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
# Fine-grained RMSD grid (Å) for the "how close is the top-ranked pose" sweep —
# a denser set than RMSD_THRESHOLDS, stepped by 0.25 Å around the canonical 2 Å
# docking-success line, so the sub-2 Å accuracy of the ranking tools' top poses
# is resolved (how many of the top-N ranked poses land within 1, 1.25, 1.5 … Å).
# Overridable via --fine-rmsd-thresholds.
FINE_RMSD_THRESHOLDS = (1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0)
# Canonical docking-success line (Å): the in-place RMSD-to-crystal below which a
# pose is "near-native". Matches the 2 Å the paper and the rest of this report use.
NEAR_NATIVE_RMSD_A = 2.0
# A near-native, PB-valid pose has the "correct form" when its best-fit (Kabsch)
# RMSD to the crystal ligand — heavy-atom RMSD AFTER optimal superposition, so
# translation and rotation are removed and only the internal conformation is
# compared — is within this many Å. Overridable via --form-ok-kabsch.
FORM_OK_KABSCH_A = 1.0
CENTROID_THRESHOLD = 4.0
# A pose counts as being in the experimentally validated pocket if its centroid is
# within this distance (Å) of the crystal ligand centroid. Looser than
# CENTROID_THRESHOLD (which marks a near-correct PLACEMENT): a pocket spans ~10 Å,
# so a centroid within ~6 Å of the crystal sits in the same site, while anything
# farther is treated as a different (decoy) pocket. Overridable via --pocket-cutoff.
POCKET_CENTROID_CUTOFF = 6.0

# Tools that produce pose rankings ordered by energy/confidence.
# All other tools (e.g. EquiBind) are oracle-only.
RANKING_TOOLS = frozenset({"autodock", "diffdock"})

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

_BASE_COLORS = {"autodock": "#1f77b4", "diffdock": "#ff7f0e"}
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
    if method == "autodock":
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


def compute_plif_recovery(crystal_mol: Chem.Mol,
                          pose_mols: list[Chem.Mol],
                          prot_mol: "plf.Molecule | None") -> list[float]:
    """Tanimoto similarity of pose IFP vs crystal IFP (single ProLIF run)."""
    n = len(pose_mols)
    if not _HAS_PROLIF or prot_mol is None:
        return [float("nan")] * n

    crystal_lig = _ligand_for_plif(crystal_mol)
    if crystal_lig is None:
        return [float("nan")] * n

    pose_ligs = [_ligand_for_plif(m) for m in pose_mols]
    valid_idx = [i for i, x in enumerate(pose_ligs) if x is not None]
    valid_ligs = [crystal_lig] + [pose_ligs[i] for i in valid_idx]

    try:
        fp = plf.Fingerprint(interactions=PLIF_INTERACTIONS)
        fp.run_from_iterable(valid_ligs, prot_mol, n_jobs=1, progress=False)
        bvs = fp.to_bitvectors()
    except Exception:
        return [float("nan")] * n

    if not bvs:
        return [float("nan")] * n
    crystal_bv = bvs[0]
    out = [float("nan")] * n
    for k, i in enumerate(valid_idx, start=1):
        try:
            out[i] = float(DataStructs.TanimotoSimilarity(crystal_bv, bvs[k]))
        except Exception:
            out[i] = float("nan")
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


def process_pair(args) -> list[dict]:
    pair_key, group_records, benchmark_dir, root = args
    protein, ligand = pair_key
    out: list[dict] = []

    pdb_id = protein
    cdir = Path(benchmark_dir) / pdb_id
    crystal_sdf = cdir / f"{pdb_id}_ligand.sdf"
    protein_pdb = cdir / f"{pdb_id}_protein.pdb"
    if not crystal_sdf.exists() or not protein_pdb.exists():
        return out

    crystal = load_first_mol(crystal_sdf)
    if crystal is None:
        return out
    crystal_h = Chem.RemoveHs(crystal)

    prot_xyz, prot_elem, prot_resid = load_protein_heavy_atoms(protein_pdb)
    _, native_contacts = clash_and_contacts(crystal_h, prot_xyz, prot_elem, prot_resid)
    n_native = len(native_contacts) or 1

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

    plif_vals = compute_plif_recovery(crystal_h,
                                      [p for _, p, _ in loaded],
                                      prot_plif)

    for (rec, pose, pose_path), plif_val in zip(loaded, plif_vals):
        try:
            rmsd = symmetry_rmsd(pose, crystal_h)
            cdist = centroid_distance(pose, crystal_h)
        except Exception:
            rmsd, cdist = float("nan"), float("nan")

        # PoseBusters' canonical RMSD vs the experimental crystal ligand.
        pb_rmsd, pb_kabsch, pb_within = posebusters_rmsd(pose, crystal_h)

        strain = uff_strain(pose)
        # "Twisted & turned" vs the crystal conformation. For the best-fit
        # "twist" RMSD prefer PoseBusters' symmetry-corrected superposed RMSD
        # (pb_kabsch_rmsd); fall back to our own Kabsch fit only when PoseBusters
        # is unavailable. The rotation angle ("turn") is always ours — PoseBusters
        # does not expose it.
        try:
            rot_angle, self_bestfit = rigid_body_fit(pose, crystal_h)
        except Exception:
            rot_angle, self_bestfit = float("nan"), float("nan")
        bestfit = pb_kabsch if not math.isnan(pb_kabsch) else self_bestfit
        tfd, max_td, mean_td, n_flip, n_rot = torsion_metrics(pose, crystal_h)
        try:
            n_clash, contacts = clash_and_contacts(pose, prot_xyz, prot_elem, prot_resid)
            recov = len(native_contacts & contacts) / n_native
        except Exception:
            n_clash, recov = 0, float("nan")

        out.append(asdict(PoseRecord(
            method=rec["docking_method"],
            protein=protein,
            ligand=ligand,
            pose_file=str(pose_path),
            pose_name=rec.get("pose_name", pose_path.stem),
            rank=parse_rank(rec["docking_method"], str(pose_path)),
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
    return pd.DataFrame(rows).set_index("method").round(2)


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
    return pd.DataFrame(rows).set_index("method").round(2)


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
    return pd.DataFrame(rows).round(2)


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
        pb_valid_only: bool = False) -> pd.DataFrame:
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
    count (``n_poses_topN``). RMSD is the symmetry-corrected heavy-atom ``rmsd`` (no
    superposition). DiffDock's duplicated top pose is de-duplicated first.
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
            sub["pb_valid"] = sub["pb_valid"].map(
                lambda x: str(x).strip().lower() in ("true", "1"))
        n_pairs = sub.groupby(["protein", "ligand"]).ngroups
        ranked = sub[sub["_rk"] <= top_n]
        top1 = (ranked[ranked["_rk"] == 1]
                .groupby(["protein", "ligand"]).first().reset_index())
        top1_hit = top1[top1["pb_valid"]] if pb_valid_only else top1
        rk_hit = ranked[ranked["pb_valid"]] if pb_valid_only else ranked
        top1_rmsd = top1_hit["rmsd"].dropna()
        rk_valid = rk_hit.dropna(subset=["rmsd"])
        best_topn = (rk_valid.groupby(["protein", "ligand"])["rmsd"].min()
                     if len(rk_valid) else pd.Series(dtype=float))
        pose_rmsd = rk_hit["rmsd"].dropna()
        n_poses = len(ranked["rmsd"].dropna())   # denominator: ALL top-N poses
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

    return pd.DataFrame(rows).round(2)


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


def _select_best_diffdock(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
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

    oracle_sum = aggregate_oracle(df)
    col = _variant_rank_col(oracle_sum)
    dd_keys = [m for m in oracle_sum.index if str(m).startswith("diffdock")]
    scores = (oracle_sum.loc[dd_keys, col].astype(float).dropna()
              if col else pd.Series(dtype=float))
    best_variant = (str(scores.idxmax()) if not scores.empty
                    else sorted(methods[dd_mask].unique())[0])

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


# ───────────────────────────────────────────────────────────────────
# Plots — Part A: oracle comparison (all tools)
# ───────────────────────────────────────────────────────────────────


def plot_oracle_rmsd_cdf(df: pd.DataFrame, out: Path) -> None:
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
    fig.tight_layout(rect=(0, 0.12, 1, 1))
    fig.savefig(out, dpi=160); plt.close(fig)


def plot_oracle_rmsd_box(df: pd.DataFrame, out: Path) -> None:
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
    fig.tight_layout(rect=(0, 0.07, 1, 1)); fig.savefig(out, dpi=160); plt.close(fig)


def _dumbbell(ax, labels, left, right, colors, *,
              left_name: str, right_name: str, value_fmt: str = "{:.0f}%",
              gap_fmt: str | None = None, gap_color: str = "#b22222",
              missing_note: str | None = None, xmax: float | None = None):
    """Horizontal dumbbell ("connected-dot") chart — a cleaner replacement for
    paired bars when the message is the *gap* between two per-row values.

    One row per ``labels`` entry (first label drawn at the top). A grey connector
    joins the ``left`` value (filled marker — the ceiling/larger quantity) and the
    ``right`` value (open marker) for that row, both in the row's ``colors`` entry.
    Rows whose ``right`` value is NaN draw only the left marker plus
    ``missing_note`` (e.g. EquiBind has no ranking). With ``gap_fmt`` the
    |left − right| gap is annotated on each connector. Returns legend handles so
    the caller controls legend placement.
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

    for yi, lv, rv in zip(y, left, right):
        if not (pd.isna(lv) or pd.isna(rv)):
            lo, hi = (lv, rv) if lv <= rv else (rv, lv)
            ax.text(lo - dx, yi, value_fmt.format(lo), ha="right", va="center",
                    fontsize=8, fontweight="bold")
            ax.text(hi + dx, yi, value_fmt.format(hi), ha="left", va="center",
                    fontsize=8, fontweight="bold")
            if gap_fmt:
                ax.text((lo + hi) / 2, yi + 0.24, gap_fmt.format(abs(lv - rv)),
                        ha="center", va="bottom", fontsize=7.5,
                        color=gap_color, fontweight="bold")
        elif not pd.isna(lv):
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
                                out: Path) -> None:
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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_pb_valid_success_bars(oracle_sum: pd.DataFrame,
                               top1_sum: pd.DataFrame,
                               out: Path) -> None:
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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_vs_posebusters_paper(oracle_sum: pd.DataFrame,
                              top1_sum: pd.DataFrame,
                              out: Path,
                              csv_out: Path | None = None) -> None:
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
                    "this_study_%": round(ours, 2) if ok else None,
                    "posebusters_paper_%": paper,
                    "delta_%": (round(ours - paper, 2) if (ok and paper is not None) else None),
                    "paper_note": r["paper_note"],
                })
        pd.DataFrame(recs).to_csv(csv_out, index=False)
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
                "rmsd_le_2A_%": round(av, 2),
                "rmsd_le_2A_and_pb_valid_%": round(vv, 2) if ok else None,
                "drop_pp": round(av - vv, 2) if ok else None,
                "valid_retention_%": round(vv / av * 100, 1) if (ok and av > 0) else None,
            })
        pd.DataFrame(recs).to_csv(csv_out, index=False)
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
            "near_native_pb_valid_%": round(100 * a_v / a_n, 2) if a_n else float("nan"),
            "oracle_near_native": b_n, "oracle_near_native_valid": b_v,
            "oracle_near_native_pb_valid_%": round(100 * b_v / b_n, 2) if b_n else float("nan"),
        })
    out = pd.DataFrame(rows).set_index("method")
    out["delta_pp"] = (out["near_native_pb_valid_%"]
                       - out["oracle_near_native_pb_valid_%"]).round(2)
    return out


def plot_within2_validity_dumbbell(comp: pd.DataFrame, out: Path, thr: float = 2.0) -> None:
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
    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{TOOL_LABEL.get(m, m)}\n{int(comp.loc[m, 'n_complexes'])} complexes · "
         f"{int(comp.loc[m, 'total_poses']):,} poses · "
         f"{int(comp.loc[m, 'near_native_poses'])} ≤ {thr:g} Å"
         for m in methods], fontsize=8)
    ax.set_xlabel(f"PB-valid share of near-native poses (RMSD ≤ {thr:g} Å) (%)")
    ax.set_xlim(0, 104)
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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


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


def _oracle_variant_specs(df_full: pd.DataFrame):
    """[(method_key, label, kind)] for the 09b variants — each ML tool's RAW variant
    beside its best-by-PB-valid&≤2Å variant. ``kind`` ∈ {'native','generation'} picks
    how the rank-1 representative is taken (native rank vs first generated pose)."""
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
    eb_raw = _raw_of(eb_best)
    specs = []

    def _add(key, label, kind):
        if key and key in present and key not in [s[0] for s in specs]:
            specs.append((key, label, kind))

    _add("autodock", "AutoDock Vina", "native")
    _add("diffdock", "DiffDock (raw)", "native")
    if dd_best != "diffdock":
        _add(dd_best, f"DiffDock* ({_opt(dd_best)})", "native")
    _add(eb_raw, "EquiBind (raw)", "generation")
    if eb_best != eb_raw:
        _add(eb_best, f"EquiBind* ({_opt(eb_best)})", "generation")
    return specs


def plot_accuracy_validity_oracle(oracle_sum: pd.DataFrame, out: Path,
                                  df_full: "pd.DataFrame | None" = None,
                                  thr: float = 2.0) -> None:
    """Fig 09b — accuracy-vs-validity bars per variant: the RANK-1 pose (EquiBind has
    no ranking, so its FIRST generated pose) beside the ORACLE (best-of-N) pose, with
    each ML tool's raw variant next to its best variant. Light bar = RMSD ≤ 2 Å, dark
    subset = also PB-valid. Falls back to a single oracle bar per (collapsed) tool when
    ``df_full`` is absent or has no crystal RMSD.
    """
    specs = _oracle_variant_specs(df_full) if df_full is not None else []

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

    fig, ax = plt.subplots(figsize=(max(9.0, 1.45 * (max(positions) + 1.4)), 6.3))
    for xi, r in zip(positions, plotted):
        c = TOOL_COLORS.get(color_key[r["variant"]], "#888888")
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
    """Per tool, contrast the RAW (un-optimised) representative pose against the
    BEST post-hoc-optimised variant, on the same accuracy/validity axes as fig 09.

    ONE representative pose is chosen per complex, then aggregated over complexes:

      * AutoDock Vina  — rank-1 pose (single variant; no ML-geometry optimisation,
                         so raw ≡ final — shown once as the physics reference).
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
        ("autodock", "final\n(rank-1)",           False, "autodock",                _rank1),
        ("diffdock", "raw\n(rank-1)",             False, "diffdock",                _rank1),
        ("diffdock", "smina-opt\n(rank-1)",       True,  "diffdock_smina",          _rank1),
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
    representative pose with its BEST post-hoc-optimised variant.

    Same light(accuracy) + dark(PB-valid subset) encoding as fig 09a/09b; bars are
    grouped by tool (colour = tool identity), raw on the left, optimised on the
    right, with a Δ callout of the validity gain. AutoDock Vina has no ML-geometry
    optimisation step, so it appears once as the physics reference. See
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

    # Raw → best validity-gain callout above each optimised bar.
    for tool in spans:
        sub = agg[agg["tool"] == tool]
        raw, best = sub[~sub["is_best"]], sub[sub["is_best"]]
        if len(raw) and len(best):
            dv = float(best[val_col].iloc[0]) - float(raw[val_col].iloc[0])
            xb, yb = float(best["_x"].iloc[0]), float(best[acc_col].iloc[0])
            ax.annotate(f"▲ +{dv:.0f} pp\nvalid vs raw", xy=(xb, yb),
                        xytext=(0, 15), textcoords="offset points",
                        ha="center", va="bottom", fontsize=8.5,
                        color="#1f7a1f", fontweight="bold")

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
        ax.text(xc, -0.175, group_name[tool], transform=trans, ha="center",
                va="top", fontsize=10.5, fontweight="bold", clip_on=False)

    ax.set_ylabel("% of receptor-ligand complexes")
    ax.set_ylim(0, 105)
    ax.set_xlim(min(positions) - 0.7, max(positions) + 0.7)
    ax.grid(axis="y", alpha=0.3)
    _accuracy_validity_legend(fig)
    fig.suptitle("Effect of post-hoc optimization — raw vs. best variant\n"
                 "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.09)
    fig.text(0.5, -0.09,
             "Representative pose per complex — Vina: rank-1 (no optimization step). "
             "DiffDock: rank-1 of the raw vs. smina-optimized run (same ranking). "
             "EquiBind (blind/unguided run, no native ranking): first generated pose "
             "vs. gnina-optimized, ranked by gnina affinity.",
             ha="center", va="top", fontsize=8, color="0.35", wrap=True)
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
            "n_near": int(near.sum()), "near_pct": round(100 * float(near.mean()), 2),
            "n_valid": int(valid.sum()), "valid_pct": round(100 * float(valid.mean()), 2),
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
    ("autodock", "AutoDock Vina",        "autodock",                "confidence rank",     "native"),
    ("diffdock", "DiffDock (raw)",       "diffdock",                "confidence rank",     "native"),
    ("diffdock", "DiffDock (smina-opt)", "diffdock_smina",          "confidence rank",     "native"),
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
    "DiffDock (raw)":       "#ff7f0e",
    "DiffDock (smina-opt)": "#d95f02",
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


def plot_topn_within_thresholds(within_df: pd.DataFrame, top_n: int,
                                thresholds: tuple[float, ...], out: Path,
                                pb_valid_only: bool = False) -> None:
    """How many of the top-ranked poses land within a fine RMSD grid.

    AutoDock Vina / DiffDock (native rank) + EquiBind (gnina-affinity ranked), two
    panels sharing the x-axis (RMSD threshold t ∈ {1, 1.25, 1.5, …} Å):
      (A) complex level — % of complexes whose RANK-1 pose is within t (solid) and
          % with ANY of the first ``top_n`` ranked poses within t (best-of-top-N,
          dashed). The vertical gap between a tool's two curves is the accuracy its
          ranking leaves on the table below that threshold.
      (B) pose level — % of ALL the tool's pooled top-N ranked poses within t, the
          literal "how many of the n top-ranked poses are within t Å".
    Dotted vertical line marks the 2 Å canonical docking-success threshold.
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
        return lab + " (gnina-ranked)" if method == "equibind" else lab

    def _fmt(ax):
        ax.axvline(2.0, color="grey", ls=":", alpha=0.7, lw=1.2)
        ax.set_xticks(thr)
        ax.set_xticklabels([f"{t:g}" for t in thr], rotation=45, ha="right")
        ax.set_xlabel("RMSD threshold vs crystal ligand (Å)")
        ax.set_ylim(0, 100); ax.grid(alpha=0.3)

    rank_note = ("Rank: AutoDock/DiffDock native; EquiBind has none, so it is ranked by "
                 "gnina affinity. Dotted line = 2 Å. RMSD = symmetry-corrected heavy-atom, "
                 "no superposition.")
    span = f"{thr[0]:g}–{thr[-1]:g} Å"
    # PB-valid variant: a pose only counts as a hit if it is ALSO PoseBusters-valid.
    valid_tag  = " — PoseBusters-valid poses only" if pb_valid_only else ""
    pose_word  = "PB-valid pose" if pb_valid_only else "pose"
    title_word = "valid-pose" if pb_valid_only else "pose"
    valid_note = (" A pose counts only if it is within t Å AND PoseBusters-valid."
                  if pb_valid_only else "")

    # ── Graph 1 (complex level): rank-1 (solid) + best-of-top-N (dashed) ──
    figA, axA = plt.subplots(figsize=(8.2, 6.0))
    for method in methods:
        sub = within_df[within_df["method"] == method].sort_values("rmsd_threshold_A")
        if sub.empty:
            continue
        color, label = TOOL_COLORS.get(method, "grey"), _label(method)
        n_pairs = int(sub["n_pairs"].iloc[0])
        axA.plot(sub["rmsd_threshold_A"], sub["top1_within_%"], marker="o", lw=2,
                 color=color, label=f"{label} — top-1 (n={n_pairs})")
        axA.plot(sub["rmsd_threshold_A"], sub["best_topN_within_%"], marker="s", ls="--",
                 lw=2, color=color, alpha=0.75, label=f"{label} — best of top-{top_n}")
    # Value labels at every whole-Ångström threshold (1, 2, 3, 4, 5 Å): each curve's
    # % where it crosses that distance, in white boxes de-collided vertically so
    # overlapping intercepts stay legible. Interior marks label to the RIGHT of a dotted
    # reference line; the rightmost mark labels to the LEFT to stay on-axis. The 2 Å
    # canonical success line (drawn emphasised by _fmt) is one of these marks.
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
        _gap, _prev = 5.0, None
        for _y, _c in marks:
            _ly = _y if _prev is None or _y - _prev >= _gap else _prev + _gap
            _prev = _ly
            axA.annotate(f"{_y:.1f}", xy=(t, _y), xytext=(t + dx, _ly),
                         va="center", ha=ha, fontsize=8, fontweight="bold", color=_c,
                         bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=_c,
                                   lw=0.8, alpha=0.9), zorder=6)

    for t in mark_ts:
        if t != 2.0:                       # 2 Å reference line already drawn by _fmt
            axA.axvline(t, color="0.75", ls=":", alpha=0.6, lw=1.0, zorder=1)
        _mark_labels(t, "left" if t == mark_ts[-1] else "right")
    _fmt(axA)
    axA.set_xlim(0, thr[-1] + (thr[-1] - thr[0]) * 0.04)   # x-axis starts at 0
    axA.set_xticks([0] + thr)
    axA.set_xticklabels(["0"] + [f"{t:g}" for t in thr], rotation=45, ha="right")
    axA.set_ylabel(f"% of complexes with a {pose_word} within the threshold")
    axA.set_title(f"Top-ranked {title_word} accuracy vs distance threshold\n"
                  "(rank-1 solid, best-of-top-N dashed)", fontweight="bold")
    axA.legend(fontsize=8)
    if not pb_valid_only:   # PB-valid complex figure is shown title-less
        figA.suptitle(_vt(f"Top-ranked {title_word} accuracy across {span}{valid_tag} \n"
                          f"(AutoDock Vina, DiffDock, EquiBind; top-{top_n})"),
                      fontsize=12, fontweight="bold")
    figA.text(0.5, 0.005, "Denominator = each tool's complexes (n pairs); best-of-top-N = "
              "the closest of the tool's first N ranked poses per complex." + valid_note
              + " " + rank_note,
              ha="center", va="bottom", fontsize=7.5, color="0.30", wrap=True)
    figA.tight_layout(rect=(0, 0.05, 1, 0.96 if not pb_valid_only else 1.0))
    figA.savefig(out.with_stem(out.stem + "_complex"), dpi=160, bbox_inches="tight")
    plt.close(figA)

    # ── Graph 2 (pose level): % of the pooled top-N poses within t ──
    figB, axB = plt.subplots(figsize=(8.2, 6.0))
    for method in methods:
        sub = within_df[within_df["method"] == method].sort_values("rmsd_threshold_A")
        if sub.empty:
            continue
        color, label = TOOL_COLORS.get(method, "grey"), _label(method)
        n_poses = int(sub["n_poses_topN"].iloc[0])
        axB.plot(sub["rmsd_threshold_A"], sub["pose_within_%"], marker="o", lw=2,
                 color=color, label=f"{label} (≤{n_poses} poses)")
    _fmt(axB)
    axB.set_ylabel(f"% of the top-N ranked poses within the threshold"
                   + (" & PB-valid" if pb_valid_only else ""))
    axB.set_title(f"How many of the top-{top_n} ranked poses are within t Å"
                  + (" AND PB-valid" if pb_valid_only else ""),
                  fontweight="bold")
    axB.legend(fontsize=8)
    figB.suptitle(_vt(f"Top-{top_n} pooled-pose accuracy across {span}{valid_tag} "
                      f"(AutoDock Vina, DiffDock, EquiBind)"),
                  fontsize=12, fontweight="bold")
    figB.text(0.5, 0.005, "Denominator = each tool's pooled rank ≤ N poses."
              + valid_note + " " + rank_note,
              ha="center", va="bottom", fontsize=7.5, color="0.30", wrap=True)
    figB.tight_layout(rect=(0, 0.05, 1, 0.96))
    figB.savefig(out.with_stem(out.stem + "_pose"), dpi=160, bbox_inches="tight")
    plt.close(figB)


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
               "success_%": round(100 * n / n_pairs, 2) if n_pairs else float("nan")}
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
                "form_correct_%": round(100 * float((form <= form_ok).mean()), 2),
                "rms_inplace": round(float(np.sqrt(ms_in)), 3),
                "rms_form": round(float(np.sqrt(ms_form)), 3),
                "rms_placement": round(float(np.sqrt((g["placement"] ** 2).mean(skipna=True))), 3),
                "form_share_of_error_%": round(100 * ms_form / ms_in, 1) if ms_in else float("nan"),
            })
        rows.append(row)
    return pd.DataFrame(rows)


def _fam_key(method: str) -> str:
    """Collapse a method to its tool family for the per-pose scatter colouring."""
    if method.startswith("autodock"):
        return "autodock"
    if method.startswith("diffdock"):
        return "diffdock"
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
    axA.set_ylabel("Form error — best-fit (Kabsch) RMSD to crystal (Å)")
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
    axD.set_xlabel("In-place RMSD to crystal (Å)")
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
                    "form_correct_%": round(100 * float((form <= form_ok).mean()), 2),
                    "inplace_median": round(float(g["inplace"].median(skipna=True)), 3),
                    "pct_placement_limited": round(float(mech_pct["placement-limited"]), 1),
                    "pct_mixed": round(float(mech_pct["mixed"]), 1),
                    "pct_form_limited": round(float(mech_pct["form-limited"]), 1),
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
    axE.set_ylabel("Form error — best-fit (Kabsch) RMSD to crystal (Å)")
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
    axE.set_ylabel("Form error — median best-fit (Kabsch) RMSD to crystal (Å)")
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


def plot_form_fidelity_gate_vs_rank(df: pd.DataFrame, out: Path,
                                    form_ok: float = FORM_OK_KABSCH_A,
                                    depths=(1, 5, 10, 15),
                                    min_n: int = 5,
                                    suptitle: str = "") -> None:
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
    _save_titled(fig, out, suptitle)


def _draw_form_placement_axes(ax, sub: pd.DataFrame, form_ok: float,
                              rmsd_thr: float, lim: float) -> None:
    """Scatter of in-place RMSD (x) vs form/best-fit RMSD (y) for ``sub``, points
    coloured by mechanism region, with the y = x bound and the r = 1/3, 2/3
    mechanism rays (form = √r · in-place)."""
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
    ax.set_xlabel("In-place RMSD to crystal (Å)")
    ax.set_ylabel("Form error — best-fit (Kabsch) RMSD (Å)")
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
        axis_max: float | None = None) -> None:
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

    ``centroid_max`` (used with ``rmsd_gate=None``) removes off-receptor outliers:
    poses whose docked centroid is > ``centroid_max`` Å from the crystal-ligand site
    (or has no centroid) are dropped, matching the code's out-of-pocket convention
    (``centroid_dist``). This trims the ejected blind-docking mode so the axes stay
    readable while the informative near-to-moderate spread the ≤ 2 Å gate hides is
    kept. No-op for crystal-free sets (no in-place RMSD / centroid defined).

    ``axis_max`` clips both axes to a fixed square (e.g. 8 Å) after the frame is
    otherwise computed, so the dense near-native core stays legible; poses that
    fall outside the frame are counted in the footnote, not silently dropped.
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
    # (the retained on-receptor poses still spread to ~15 Å in in-place RMSD).
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
            _draw_form_placement_axes(ax, sub, form_ok, vlim, lim)
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
            n_rec = int(sub["protein"].nunique())  # distinct receptors (unique complexes)
            dn = len(sub) - n1
            extra = f"  (+{dn} poses vs top-{d1})" if d != d1 and dn else ""
            ax.set_title(f"{flabel} · top-{d}   n={len(sub)} · {n_rec} rec.{extra}",
                         fontsize=9.5)
            if ri == 0 and ci == 0:
                ax.legend(fontsize=7.5, frameon=False, loc="upper left", title="mechanism")
    if rmsd_gate is not None:
        sel_line = ("cumulative PB-valid poses within each method's top-d "
                    f"(≤ {vlim:g} Å, near-native)")
    else:
        sel_line = "cumulative PB-valid poses within each method's top-d (RMSD gate removed)"

    # Short title; the descriptive detail that used to be stacked in it now sits in a
    # footnote so the heading rides just above the panels instead of floating off.
    fig.suptitle(_vt(
        "Form vs placement across ranking depth — top-1 net widened to "
        f"top-{'/'.join(f'{d}' for d in depths[1:])}"),
        fontsize=13, fontweight="bold", y=0.995)
    notes = [f"Rows = tool family · columns = {sel_line}."]
    if rmsd_gate is None and centroid_max is not None:
        notes.append(f"Off-receptor poses removed — docked centroid > {centroid_max:g} Å from "
                     f"the crystal-ligand site ({n_dropped} dropped at top-{depths[-1]}).")
    if n_clipped:
        notes.append(f"Axes clipped at {lim:g} Å for legibility — {n_clipped} pose(s) beyond "
                     f"the frame at top-{depths[-1]} not shown.")
    notes.append("Titles: n = PB-valid poses in the cohort; rec. = distinct receptors "
                 "(unique protein–ligand complexes) those poses span.")
    notes.append("Coloured by mechanism (r = form² / in-place²); black ✕ = cohort centroid "
                 "(median in-place, median form); orange arrow = drift from top-1.")
    fig.text(0.5, 0.008, "\n".join(notes), ha="center", va="bottom",
             fontsize=8, color="0.30", linespacing=1.35)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
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
        cry = bdir / str(row["protein"]) / f"{row['protein']}_ligand.sdf"
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


def _example_geometry(exemplar: dict, benchmark_dir, root) -> dict | None:
    """Load the crystal + docked heavy-atom mols for one exemplar and best-fit the
    pose onto the crystal. Returns dict of coords / bonds / elements (as docked and
    superimposed), or None on any load/align failure."""
    from rdkit.Chem import rdMolAlign
    prot = exemplar["protein"]
    pf = Path(str(exemplar["pose_file"]))
    if not pf.is_absolute():
        pf = Path(root) / pf
    cpath = Path(benchmark_dir) / str(prot) / f"{prot}_ligand.sdf"
    pose, crystal = load_first_mol(pf), load_first_mol(cpath)
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
    """Write one PDB per region — chain A = crystal (original) ligand, chain B =
    docked pose (as docked, real receptor frame), chain C = docked pose after
    best-fit superposition onto the crystal (shape-only) — plus a master .pml with
    two toggleable scenes (as-docked / aligned) and a manifest CSV. Returns paths."""
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
            "REMARK   chain A = crystal (original) ligand",
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
           "# chain A = crystal (original) · chain B = docked as-docked · "
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
    manifest.to_csv(man_path, index=False)
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
    """Per region, write receptor + crystal (original) ligand + docked pose in one
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
            f"resn CRY = crystal (original) ligand (chain {cry_ch}) | "
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
    pml = ["# Receptor + crystal (original, resn CRY) + docked pose (resn DOK) per example",
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
            "# green = crystal (original), coloured = docked pose; grey = receptor",
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

    handles = [Line2D([0], [0], color="#333333", lw=3, label="crystal (original) ligand")]
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
                "pct_at_rank": round(pct, 2),
                "cum_pct": round(cum, 2),
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
    ("AutoDock Vina",        "autodock",                "native"),
    ("DiffDock (raw)",       "diffdock",                "native"),
    ("DiffDock (smina-opt)", "diffdock_smina",          "native"),
    ("EquiBind (raw)",       "equibind_unguided_raw",   "generation"),
    ("EquiBind (gnina-opt)", "equibind_unguided_gnina", "gnina"),
]


def aggregate_topk_recovery(df: pd.DataFrame, ks, thr: float = 2.0) -> pd.DataFrame:
    """Per variant × depth k: cumulative recovery within the variant's top-k poses.

    ``near_recovery_%`` = % of complexes with ≥ 1 top-k pose RMSD ≤ ``thr``;
    ``valid_recovery_%`` = % with ≥ 1 top-k pose RMSD ≤ ``thr`` AND PB-valid (the
    validity-aware line). Both are monotone non-decreasing in k and valid ≤ near, so
    the gap between them is the accurate-but-invalid share at depth k. Needs the FULL
    frame so raw ``diffdock`` and ``diffdock_smina`` are both present.
    """
    ks = sorted({int(k) for k in ks})
    present = set(df["method"].astype(str))
    rows = []
    for variant, mkey, kind in _TOPK_RECOVERY_SPECS:
        if mkey not in present:
            continue
        sub = df[df["method"].astype(str) == mkey].copy()
        rk = (sub["rank"] if kind == "native"
              else sub["pose_name"].map(_equibind_pose_index) if kind == "generation"
              else _gnina_affinity_rank(sub))
        sub["_rk"] = pd.to_numeric(rk, errors="coerce")
        sub = sub.dropna(subset=["_rk", "rmsd"])
        N = sub.groupby(["protein", "ligand"]).ngroups
        if not N:
            continue
        near = sub["rmsd"] <= thr
        valid = near & _to_bool(sub["pb_valid"])
        near_rank = sub[near].groupby(["protein", "ligand"])["_rk"].min()
        valid_rank = sub[valid].groupby(["protein", "ligand"])["_rk"].min()
        for k in ks:
            rows.append({
                "variant": variant, "method_key": mkey, "k": k, "n_complexes": N,
                "near_recovery_%": round(100 * int((near_rank <= k).sum()) / N, 2),
                "valid_recovery_%": round(100 * int((valid_rank <= k).sum()) / N, 2),
            })
    return pd.DataFrame(rows)


def plot_topk_recovery_validity(rec: pd.DataFrame, out: Path, thr: float = 2.0) -> None:
    """Fig 15b — cumulative recovery within top-k: one solid (near-native, RMSD ≤ 2 Å)
    and one dashed (validity-aware: near-native AND PB-valid) line per variant, colour
    = variant. Raw DiffDock and DiffDock* (smina) appear side by side; the gap between
    a variant's two lines is its accurate-but-invalid share at that depth. Unlike fig
    15 (rank of the min-RMSD pose, RMSD only), this brings PB-validity into the view.
    Every curve is annotated with its recovery % (one decimal) where it crosses each
    x-axis tick — k = 1, 5, 10, 15, 20, 25, 30 (dotted reference lines): near-native
    values left of each line, validity-aware right of it."""
    if rec is None or rec.empty:
        return
    from matplotlib.lines import Line2D
    variants = [v for v in (s[0] for s in _TOPK_RECOVERY_SPECS)
                if v in set(rec["variant"])]
    fig, ax = plt.subplots(figsize=(9.8, 6.2))
    for variant in variants:
        g = rec[rec["variant"] == variant].sort_values("k")
        c = _RANK1_TOPN_COLORS.get(variant, "#888888")
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
        ax.axvline(mk, color="0.6", ls=(0, (2, 3)), lw=1.0, zorder=1)
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
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


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


def plot_optimization_benefit_by_rank(rec: pd.DataFrame, out: Path) -> None:
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


def plot_optimization_benefit_by_group(rec: pd.DataFrame, out: Path) -> None:
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
            "pct_topN_in_validated_pocket": round(pct_topN, 2),
            "pct_top1_in_validated_pocket":
                round(100 * float(pc["top1_on"].mean()), 2) if nc else nan,
            "pct_any_topN_in_validated_pocket":
                round(100 * float(pc["any_on"].mean()), 2) if nc else nan,
            "pct_oracle_in_validated_pocket":
                round(100 * float(pc["oracle_on"].mean()), 2) if nc else nan,
            "pct_oracle_diff_pocket_vs_topN":
                round(100 * float(pc["oracle_diff"].mean()), 2) if nc else nan,
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
                "pct_in_pocket": round(100 * float(at_k["on_pocket"].mean()), 2)
                                 if n_k else nan,
                "pct_any_in_top_k": round(100 * float(any_by_complex.mean()), 2)
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

    if gate_rmsd:
        rmsd_ok = (pd.to_numeric(merged["rmsd"], errors="coerce") <= 2.0).fillna(False)
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
    return {"N": N, "selection": selection, "steps": steps}


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
                "remaining_pct": round(100 * s["remaining"] / N, 2) if N else float("nan"),
            })
    pd.DataFrame(rows).to_csv(path, index=False)


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
        return "EquiBind (gnina-ranked)"
    if method == "equibind_raw":
        return "EquiBind (raw · gnina-ranked)"
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
                             "top1_remaining_%": round(s1[i], 2),
                             f"top{top_n}_remaining_%": round(sn[i], 2),
                             "gain_pp": round(d, 2)})

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
        pd.DataFrame(csv_rows).to_csv(csv_out, index=False)


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


def plot_pb_top30_recovery_by_test(cascades_top1: dict, cascades_topn: dict,
                                   top_n: int, out: Path,
                                   csv_out: Path | None = None,
                                   extra: tuple = (),
                                   methods: tuple | None = None,
                                   labels: dict | None = None,
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
        raw_key = {"diffdock": "diffdock_raw",
                   "equibind_gnina": "equibind_raw"}.get(str(method))
        if raw_key and raw_key in cascades_top1 and raw_key in cascades_topn:
            rp1 = 100 * cascades_top1[raw_key]["steps"][-1]["after"] / (
                cascades_top1[raw_key]["N"] or 1)
            rp30 = 100 * cascades_topn[raw_key]["steps"][-1]["after"] / (
                cascades_topn[raw_key]["N"] or 1)
            note += (f"\nraw → refined:  +{p1 - rp1:.1f} pp at rank-1  ·  "
                     f"+{p30 - rp30:.1f} pp at top-{top_n}")
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
    fig.suptitle(_vt(f"Extra complexes recovered by best-of-top-{top_n} vs rank-1, "
                     "at each PoseBusters cascade step\n(bar = extra survivors "
                     "where the top-30 pick makes a difference; blue = RMSD ≤ 2 Å "
                     "placement, amber = other physical-validity tests, green = "
                     "passing all tests)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        pd.DataFrame(csv_rows).to_csv(csv_out, index=False)


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


def aggregate_pb_waterfall_raw_vs_best(df: pd.DataFrame, test_table: pd.DataFrame,
                                       test_cols: list[str],
                                       gate_rmsd: bool = True) -> dict:
    """PoseBusters cascade, top-1 pose, for the RAW pose vs its best refined
    variant — as an ordered dict of panels ``{key: {..cascade.., title, role, pair,
    vs}}``. DiffDock is refined/ranked with smina; EquiBind with gnina.

    ``df`` must retain every optimiser/refinement variant under its own method key
    (the full per-pose frame, before any best-variant collapse). Selection:

    * AutoDock — rank-1 pose. Physics docking has no smina/gnina re-optimisation,
      so it appears once as a reference (``role="reference"``).
    * DiffDock — rank-1 raw pose (``diffdock``) vs rank-1 smina pose
      (``diffdock_smina``); the optimiser inherits the DiffDock rank, so the two
      panels are the same pose before/after refinement.
    * EquiBind — has no native ranking, so the representative is the prediction
      ranked #1 by its gnina refinement score, shown in its raw and its
      gnina-refined geometry (``_equibind_gnina_ranked_top1``).

    ``gate_rmsd`` (default True) drops RMSD > 2 Å before the PoseBusters tests;
    with ``gate_rmsd=False`` every ranked-1 pose faces the tests regardless of
    placement accuracy, so the cascade shows where poses fail on physical-validity
    grounds alone. Best panels carry ``vs`` = the key of their raw counterpart so
    the plotter can annotate the net change in poses passing every test.
    """
    out: dict = {}

    ref = _pb_cascade_from_rep(_rank1_rep(df, "autodock"), test_table, test_cols,
                               "top-1 pose", gate_rmsd=gate_rmsd)
    if ref is not None:
        out["autodock"] = {**ref, "role": "reference", "pair": "autodock",
                           "title": "AutoDock Vina — top-1 pose "
                                    "(physics docking · no smina/gnina re-optimisation)"}

    dd_raw = _pb_cascade_from_rep(_rank1_rep(df, "diffdock"), test_table,
                                  test_cols, "raw pose", gate_rmsd=gate_rmsd)
    dd_best = _pb_cascade_from_rep(_rank1_rep(df, "diffdock_smina"), test_table,
                                   test_cols, "smina-refined pose", gate_rmsd=gate_rmsd)
    if dd_raw is not None:
        out["diffdock_raw"] = {**dd_raw, "role": "raw", "pair": "diffdock",
                               "title": "DiffDock — raw top-1 pose"}
    if dd_best is not None:
        e = {**dd_best, "role": "best", "pair": "diffdock",
             "title": "DiffDock + smina — best variant, top-1 pose"}
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
            e = {**eq_best, "role": "best", "pair": "equibind",
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
            # EquiBind is refined/ranked with gnina, DiffDock with smina.
            opt_word = "gnina" if c.get("pair") == "equibind" else "smina"
            ax.text(0.995, 0.975,
                    f"{opt_word} vs raw:  {'+' if d >= 0 else '−'}{abs(d):.1f} pp "
                    f"{delta_label}",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=9.5, fontweight="bold", color=accent,
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec=accent, alpha=0.92))
            csv_rows.append({"pair": c.get("pair"),
                             "raw_passing_all_%": round(raw_pct, 1),
                             "best_passing_all_%": round(final_pct, 1),
                             "gain_pp": round(d, 1)})

    axes[-1].set_xticks(np.arange(nx))
    axes[-1].set_xticklabels(labels, rotation=45, ha="right",
                             rotation_mode="anchor", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    if gate_rmsd:
        suptitle = ("PoseBusters test-failure waterfall — raw pose vs best "
                    "refined variant, top-1 per method (DiffDock: smina; EquiBind: "
                    "gnina)\n(sequential filter: each red bar = poses removed by "
                    "that test; green = poses passing every test; left edge: grey = "
                    "raw, green = optimised)")
    else:
        suptitle = ("PoseBusters physical-validity waterfall — raw pose vs best "
                    "refined variant, top-1 per method (DiffDock: smina; EquiBind: "
                    "gnina)\n(RMSD ≤ 2 Å NOT applied — every ranked-1 pose is tested; "
                    "each red bar = poses removed by that test; green = poses passing "
                    "every PoseBusters test, at any RMSD)")
    fig.suptitle(_vt(suptitle), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    if csv_out is not None and csv_rows:
        pd.DataFrame(csv_rows).to_csv(csv_out, index=False)


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
    # Split DiffDock optimizer variants, then (unless we keep all three for
    # --best-diffdock-only) keep only the selected one (raw 'diffdock' by default)
    # as the canonical DiffDock entry.
    df = _apply_diffdock_split(df)
    if select_diffdock:
        df = _select_diffdock_variant(df, diffdock_variant)
    keep = ["docking_method", "protein", "ligand", "pose_file", "pose_name", "pb_valid"]
    keep += [c for c in EQ_PROVENANCE_COLS if c in df.columns]
    if "optimizer" in df.columns:          # DiffDock optimizer provenance, for drill-down
        keep.append("optimizer")
    return df[keep]


# ───────────────────────────────────────────────────────────────────
# Per-pose cache (so re-runs that only change --top-n skip the heavy scoring)
# ───────────────────────────────────────────────────────────────────

# Columns a cached per_pose_metrics.csv must contain to be reusable by the current
# code (guards against reusing a CSV written by an older, narrower schema).
_CACHE_REQUIRED_COLS = {"method", "protein", "ligand", "pose_name", "rank",
                        "rmsd", "centroid_dist", "pb_valid"}


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
        "schema": 4,  # bumped: per-pose rows now carry gnina_affinity (EquiBind ranking)
        "pb_csv": _file_fingerprint(Path(args.pb_csv)),
        "ids_file": str(args.ids_file) if args.ids_file else None,
        "split_equibind": bool(args.split_equibind),
        "diffdock_variant": str(args.diffdock_variant),
        "best_diffdock_only": bool(args.best_diffdock_only),
        "limit_pairs": int(args.limit_pairs),
    }


def _read_cached_per_pose(per_pose_csv: Path):
    """Read + validate the cached per-pose CSV; ``None`` if missing/invalid."""
    if not per_pose_csv.exists():
        return None
    try:
        df = pd.read_csv(per_pose_csv, low_memory=False)
    except Exception:
        return None
    if not _CACHE_REQUIRED_COLS.issubset(df.columns):
        print("  cached per-pose metrics lack required columns — cannot reuse.")
        return None
    return df


def _load_cached_per_pose(args, sig: dict):
    """Return the cached per-pose DataFrame, or ``None`` to trigger a recompute.

    With ``--reuse-cache`` the existing per_pose_metrics.csv is loaded WITHOUT
    checking the input signature (use when only the plotting/aggregation code
    changed — the caller asserts the cached pairs are the ones they want). This
    bypasses ``--force`` too. Otherwise the cache is reused only when the manifest
    signature matches the current inputs.
    """
    per_pose_csv = args.out_dir / "per_pose_metrics.csv"
    manifest = args.out_dir / "per_pose_metrics.manifest.json"

    if getattr(args, "reuse_cache", False):
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
        print("  cached per-pose metrics are stale (inputs changed) — recomputing.")
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
    df.to_csv(per_pose_csv, index=False)
    (args.out_dir / "per_pose_metrics.manifest.json").write_text(
        json.dumps({"signature": sig, "top_n": int(args.top_n),
                    "n_rows": int(len(df))}, indent=2))
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
    tbl.round(2).to_csv(out)
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
                    help="Centroid distance (Å) to the crystal ligand within which "
                         "a pose counts as being in the experimentally validated "
                         "pocket, for the pocket-localization analysis "
                         "(default %(default)s).")
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
    ap.add_argument("--best-variants-only", action="store_true",
                    help="Convenience umbrella: enable BOTH --best-equibind-only and "
                         "--best-diffdock-only, so every summary and plot shows only the "
                         "single best EquiBind and best DiffDock variant (i.e. AutoDock, "
                         "DiffDock*, EquiBind*) instead of all 18 EquiBind + 3 DiffDock "
                         "variants. per_pose_metrics.csv still keeps every variant for "
                         "drill-down.")
    ap.add_argument("--collapse-plots-only", action="store_true",
                    help="Keep EVERY variant in the summary TABLES but show only the "
                         "single best variant per tool in the GRAPHS. Writes the full "
                         "all-variants tables (oracle_summary_all_variants.csv, "
                         "top1_summary_all_variants.csv, and rmsd_vs_pbvalid_all_variants.csv "
                         "— the RMSD vs RMSD+PB-validity comparison per variant) before "
                         "collapsing to AutoDock, DiffDock*, EquiBind* for every plot and "
                         "oracle_summary.csv. Like --best-variants-only but tables stay full.")
    args = ap.parse_args()

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
        print("collapse-plots-only: all variants kept in the *_all_variants.csv tables; "
              "plots + oracle_summary.csv collapse to AutoDock, DiffDock*, EquiBind*.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    root = Path.cwd()

    # ── Per-pose metrics: reuse the cached results if the inputs are unchanged,
    #    otherwise run the (heavy) per-pose scoring and cache it. --top-n is NOT
    #    part of the cache key, so changing it alone reuses the cache. ──
    sig = _per_pose_signature(args)
    df = _load_cached_per_pose(args, sig)

    if df is None and args.reuse_cache:
        print(f"ERROR: --reuse-cache was set but no usable per-pose cache exists "
              f"in {args.out_dir}.\n       Run once without --reuse-cache to build "
              f"per_pose_metrics.csv, then re-run with it.")
        return

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

        work = [((p, l), recs, str(args.benchmark_dir.resolve()), root)
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
    if not _w2cmp.empty:
        _w2cmp.to_csv(args.out_dir / "within2_validity_comparison.csv")
        # Figure keeps fpocket/p2rank as separate rows but brackets them as one guided block —
        # they are handed the pocket, so warrant reading apart from the blind methods.
        plot_within2_validity_dumbbell(
            _w2cmp, args.out_dir / "19_within2_validity_dumbbell.png")
        print(f"  wrote near-native validity comparison → within2_validity_comparison.csv "
              f"+ 19_within2_validity_dumbbell.png ({len(_w2cmp)} variants)")

    # ── --collapse-plots-only: dump the FULL (all-variants) summary tables BEFORE the
    #    best-per-tool collapse below, so the plots use only the best variant while the
    #    tables retain every EquiBind + DiffDock variant. Includes the RMSD vs
    #    RMSD+PB-validity comparison per variant (rmsd_vs_pbvalid_all_variants.csv).
    if args.collapse_plots_only:
        _oracle_all = aggregate_oracle(df)
        _top1_all = aggregate_top1(df)
        _oracle_all.to_csv(args.out_dir / "oracle_summary_all_variants.csv")
        _top1_all.to_csv(args.out_dir / "top1_summary_all_variants.csv")
        _write_rmsd_vs_pbvalid_table(
            _oracle_all, _top1_all, args.out_dir / "rmsd_vs_pbvalid_all_variants.csv")
        print(f"collapse-plots-only: wrote all-variants tables for {len(_oracle_all)} "
              "variants (oracle_summary_all_variants.csv, top1_summary_all_variants.csv, "
              "rmsd_vs_pbvalid_all_variants.csv) — plots below use the best per tool.")

    # Keep the FULL per-variant frame (before any best-variant collapse) so the
    # raw-vs-best optimization figure (09c) can still see every DiffDock/EquiBind
    # variant under its own method key.
    df_full = df.copy()

    # ── Optional: restrict the report to the single best EquiBind variant ──
    # Done after the full per-pose dump (which keeps every variant) so only the
    # summaries and plots below are filtered.
    if args.best_equibind_only:
        df, best_variant = _select_best_equibind(df)
        if best_variant:
            _LABEL_OVERRIDES[best_variant] = "EquiBind*"
            print(f"best-equibind-only: '{best_variant}' is the top EquiBind variant "
                  "by PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'EquiBind*').")
        else:
            print("best-equibind-only: no EquiBind variants present — nothing filtered.")

    # ── Optional: restrict the report to the single best DiffDock variant ──
    # Scores all three optimizer variants (kept because select_diffdock was off),
    # collapses to the best, relabelled 'DiffDock*'.
    if args.best_diffdock_only:
        df, best_dd = _select_best_diffdock(df)
        if best_dd:
            _LABEL_OVERRIDES["diffdock"] = "DiffDock*"
            print(f"best-diffdock-only: '{best_dd}' is the top DiffDock variant by "
                  "PB-Valid AND RMSD ≤ 2 Å — keeping only it as 'diffdock' (shown as 'DiffDock*').")
        else:
            print("best-diffdock-only: no DiffDock variants present — nothing filtered.")

    # ── Aggregation ──────────────────────────────────────────────
    oracle_sum = aggregate_oracle(df)
    oracle_sum.to_csv(args.out_dir / "oracle_summary.csv")
    print("\nOracle summary (all tools):")
    print(oracle_sum.to_string())

    top1_sum = aggregate_top1(df)
    top1_sum.to_csv(args.out_dir / "top1_summary.csv")
    if not top1_sum.empty:
        print("\nTop-1 summary (ranking tools):")
        print(top1_sum.to_string())

    rank_df = aggregate_by_rank(df, args.top_n)
    rank_df.to_csv(args.out_dir / "per_rank_metrics.csv", index=False)

    ifp_rank_df = aggregate_ifp_by_rank(df, args.top_n)
    ifp_rank_df.to_csv(args.out_dir / "per_rank_ifp_recovery.csv", index=False)
    print("  wrote per-rank interaction-fingerprint recovery → per_rank_ifp_recovery.csv")

    oracle_pivot = _oracle_per_pair(df).pivot_table(
        index=["protein", "ligand"], columns="method", values="rmsd")
    oracle_pivot.to_csv(args.out_dir / "oracle_rmsd_per_pair.csv")

    # ── Part A: oracle comparison (all tools) ────────────────────
    print("\nGenerating plots …")
    plot_oracle_rmsd_cdf(df, args.out_dir / "01_oracle_rmsd_cdf.png")
    plot_oracle_rmsd_box(df, args.out_dir / "02_oracle_rmsd_boxplot.png")
    plot_oracle_vs_top1_success(oracle_sum, top1_sum,
                                args.out_dir / "03_oracle_vs_top1_success.png")
    # Accuracy-vs-validity (former two-panel fig 09), now split into 09a (top-1)
    # and 09b (oracle), plus 09c contrasting each tool's raw vs best-optimised pose.
    plot_accuracy_validity_top1(top1_sum,
                                args.out_dir / "09a_accuracy_vs_validity_top1.png")
    plot_accuracy_validity_oracle(oracle_sum,
                                  args.out_dir / "09b_accuracy_vs_validity_oracle.png",
                                  df_full=df_full)
    rvb = aggregate_optimization_raw_vs_best(df_full)
    if not rvb.empty:
        # role_label carries newlines for the two-line bar labels — keep the clean
        # single-line `role` column in the CSV.
        rvb.drop(columns=["role_label"]).to_csv(
            args.out_dir / "optimization_raw_vs_best.csv", index=False)
        plot_optimization_raw_vs_best(
            rvb, args.out_dir / "09c_optimization_raw_vs_best.png")
        print("  wrote raw-vs-best optimization comparison → "
              "optimization_raw_vs_best.csv")
        # Guard: 09c needs the raw + optimised variant of each ML tool. If a run
        # built df without the smina/gnina DiffDock variants (e.g. a plain default
        # run, no --diffdock-variant all), the smina-opt bar is silently absent.
        if not (rvb["method_key"] == "diffdock_smina").any() \
                and (rvb["method_key"] == "diffdock").any():
            print("  WARNING: 09c is missing DiffDock's smina-opt bar — df lacks "
                  "diffdock_smina (re-run with --diffdock-variant all).")

    # 09d — best-of-top-d at several depths (top-1, top-N, top-30), raw vs refined
    # (ranking headroom + optimization). Uses df_full so BOTH the raw and the
    # optimised run of DiffDock (smina) / EquiBind (gnina) are present under their
    # own method keys.
    r1tn_depths = sorted({1, args.top_n, 30})
    r1tn = aggregate_rank1_vs_topn(df_full, r1tn_depths)
    if not r1tn.empty:
        r1tn.to_csv(args.out_dir / "rank1_vs_topn.csv", index=False)
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
                  "method keys (re-run with --diffdock-variant all --split-equibind).")

    plot_pb_valid_success_bars(oracle_sum, top1_sum,
                               args.out_dir / "10_pb_valid_rmsd2_success_bars.png")
    plot_vs_posebusters_paper(oracle_sum, top1_sum,
                              args.out_dir / "11_vs_posebusters_paper.png",
                              args.out_dir / "comparison_vs_posebusters_paper.csv")
    plot_rmsd2_vs_pbvalid_grouped(oracle_sum, top1_sum,
                                  args.out_dir / "12_rmsd2_vs_pbvalid_grouped.png")
    plot_pbvalid_filter_influence(oracle_sum,
                                  args.out_dir / "13_pbvalid_filter_influence.png",
                                  args.out_dir / "pbvalid_filter_influence.csv")

    # ── "Twisted & turned": how the docked ligand deviates from its crystal
    #    conformation (translation / rotation / internal twisting), for the
    #    top-1 and oracle pose of each tool ───────────────────────────────
    aggregate_twist_turn(df, args.top_n).to_csv(
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
        aggregate_form_fidelity(df, args.form_ok_kabsch, selector=selector).to_csv(
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
    # Same filmstrip with the ≤ 2 Å gate REMOVED (all PB-valid poses), off-receptor
    # outliers dropped (docked centroid > FAR_FROM_RECEPTOR_CENTROID_A from the native
    # site) so the axes stay readable while the near-to-moderate spread the gate hides
    # is exposed — the deeper ranks visibly add mis-positioned but still-valid poses.
    plot_form_vs_placement_depth_filmstrip(
        df, args.out_dir / "20d_form_vs_placement_by_family__depth_filmstrip_pbvalid.png",
        args.form_ok_kabsch, rmsd_gate=None, depths=(1, 3, 5),
        centroid_max=FAR_FROM_RECEPTOR_CENTROID_A, axis_max=8.0)
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
        aggregate_form_fidelity_by_depth(df, args.form_ok_kabsch, rmsd_gate=gate).to_csv(
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
    plot_form_fidelity_gate_vs_rank(
        df, args.out_dir / "20_form_fidelity__rank_depth_gate_vs_rank.png",
        args.form_ok_kabsch, depths=(1, 5, 10, 15), suptitle=_vt(
            "Form fidelity vs rank (top-15) — impact of the RMSD ≤ 2 Å criterion\n"
            "(A) PB-valid   (B) near-native (PB-valid AND ≤ 2 Å)   (C) gate effect\n"
            "curves = median best-fit RMSD (band = IQR); labels = % form-correct · n complexes"))
    # Head-to-head of the two selections: table (how many complexes each rule
    # keeps + how many the nearest rule needlessly drops) and figure (20c).
    aggregate_oracle_selection_comparison(df).to_csv(
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
          "(gated + __depth_filmstrip_pbvalid gate-removed, off-receptor centroid-trimmed); "
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
        odist.to_csv(args.out_dir / "oracle_rank_distribution.csv", index=False)
        plot_oracle_rank_distribution(odist, args.top_n,
                                      args.out_dir / "15_oracle_rank_distribution.png")
        print("  wrote oracle-rank distribution → oracle_rank_distribution.csv")

        # 15b — cumulative recovery within top-k, near-native (solid) vs. PB-valid
        # (dashed), raw DiffDock alongside DiffDock* (smina). Uses df_full so both
        # DiffDock variants are present; fig 15 itself ignores PB-validity.
        topk_rec = aggregate_topk_recovery(df_full, range(1, 31))
        if not topk_rec.empty:
            topk_rec.to_csv(args.out_dir / "topk_recovery_validity.csv", index=False)
            plot_topk_recovery_validity(
                topk_rec, args.out_dir / "15b_topk_recovery_validity.png")
            print("  wrote top-k recovery (near-native vs PB-valid) → "
                  "topk_recovery_validity.csv")

        # 15c/15d — per-rank benefit of pose optimization (smina/gnina) over the raw
        # pose, for DiffDock (confidence rank) and EquiBind (generation order): does
        # refinement help the early ranks more? 15c = per-rank curves, 15d = the +pp
        # impact summarised over rank groups 1-10 / 11-20 / 21-30. Uses df_full so raw +
        # smina + gnina variants of both tools are present.
        opt_benefit = aggregate_optimization_benefit(df_full, max_rank=30)
        if not opt_benefit.empty:
            opt_benefit.to_csv(args.out_dir / "optimization_benefit_by_rank.csv",
                               index=False)
            plot_optimization_benefit_by_rank(
                opt_benefit, args.out_dir / "15c_optimization_benefit_by_rank.png")
            plot_optimization_benefit_by_group(
                opt_benefit, args.out_dir / "15d_optimization_benefit_by_group.png")
            print("  wrote per-rank optimization benefit → "
                  "optimization_benefit_by_rank.csv")

        # How many of the top-N ranked poses land within 1, 1.25, 1.5 … Å?
        # EquiBind (gnina-ranked) is added from df_full since the collapsed df keeps
        # only the best EquiBind variant (which may not be the gnina one).
        eq_gnina = df_full[df_full["method"].astype(str) == "equibind_unguided_gnina"]
        within_df = aggregate_topn_within_thresholds(
            df, args.top_n, args.fine_rmsd_thresholds,
            eq_df=eq_gnina if not eq_gnina.empty else None)
        if not within_df.empty:
            within_df.to_csv(args.out_dir / "topn_within_thresholds.csv", index=False)
            plot_topn_within_thresholds(within_df, args.top_n,
                                        args.fine_rmsd_thresholds,
                                        args.out_dir / "18_topn_within_thresholds.png")
            thr_str = ", ".join(f"{t:g}" for t in args.fine_rmsd_thresholds)
            print(f"\nTop-ranked poses within RMSD thresholds ({thr_str} Å):")
            print(within_df.pivot_table(index="method", columns="rmsd_threshold_A",
                                        values="top1_within_%").to_string())
            print("  wrote within-threshold sweep → topn_within_thresholds.csv")

        # 18 PB-valid variant — same sweep but a pose only counts as a hit if it
        # is BOTH within t Å AND PoseBusters-valid (combined near-native + valid).
        within_pbv = aggregate_topn_within_thresholds(
            df, args.top_n, args.fine_rmsd_thresholds,
            eq_df=eq_gnina if not eq_gnina.empty else None, pb_valid_only=True)
        if not within_pbv.empty:
            within_pbv.to_csv(args.out_dir / "topn_within_thresholds_pbvalid.csv",
                              index=False)
            plot_topn_within_thresholds(
                within_pbv, args.top_n, args.fine_rmsd_thresholds,
                args.out_dir / "18_topn_within_thresholds_pbvalid.png",
                pb_valid_only=True)
            print("  wrote PB-valid within-threshold sweep → "
                  "topn_within_thresholds_pbvalid.csv")
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
        pocket_sum.to_csv(args.out_dir / "pocket_localization.csv")
        if not pocket_by_rank.empty:
            pocket_by_rank.to_csv(args.out_dir / "pocket_localization_by_rank.csv",
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
                    plot_pb_top30_recovery_by_test(
                        c1_30, cn_30, TOP30,
                        args.out_dir / "17f_pb_top30_recovery_by_test.png",
                        csv_out=args.out_dir / "pb_top30_recovery_by_test.csv",
                        extra=ex30)
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
                        plot_pb_top30_recovery_by_test(
                            c1_raw, cn_raw, TOP30,
                            args.out_dir / "17g_pb_top30_recovery_by_test_raw.png",
                            csv_out=args.out_dir / "pb_top30_recovery_by_test_raw.csv",
                            methods=raw_methods)
                        print("  wrote un-refined (raw) top-30 recovery → "
                              "17c3_pb_waterfall_top1_vs_top30_raw.png, "
                              "17g_pb_top30_recovery_by_test_raw.png")

                        # 17h — raw AND refined in ONE figure: AutoDock (physics,
                        # once) then each tool's raw panel directly above its
                        # refined panel (DiffDock: smina; EquiBind: gnina), so the
                        # refinement benefit reads off panel-to-panel. Reuses the
                        # (already verified) cascades.
                        c1_all = {**c1_raw, **c1_30}
                        cn_all = {**cn_raw, **cn_30}
                        combined_methods = tuple(
                            m for m in ("autodock", "diffdock_raw", "diffdock",
                                        "equibind_raw", "equibind_gnina")
                            if m in c1_all and m in cn_all)
                        combined_labels = {
                            "diffdock": "DiffDock* (smina-refined)",
                            "diffdock_raw": "DiffDock (raw)",
                            "equibind_gnina": "EquiBind (gnina-refined)",
                            "equibind_raw": "EquiBind (raw)",
                        }
                        if len(combined_methods) >= 2:
                            plot_pb_top30_recovery_by_test(
                                c1_all, cn_all, TOP30,
                                args.out_dir / "17h_pb_top30_recovery_raw_vs_refined.png",
                                csv_out=args.out_dir / "pb_top30_recovery_raw_vs_refined.csv",
                                methods=combined_methods, labels=combined_labels)
                            print("  wrote combined raw+refined top-30 recovery → "
                                  "17h_pb_top30_recovery_raw_vs_refined.png")

            # Raw pose vs best refined variant, top-1 per method — 17d. Uses
            # df_full (every optimiser/refinement variant under its own method key,
            # before the best-variant collapse) so raw and refined poses can be
            # cascaded side by side. DiffDock is refined/ranked with smina; EquiBind,
            # having no native rank, is ranked by its gnina refinement score.
            rvb = aggregate_pb_waterfall_raw_vs_best(df_full, test_table, test_cols)
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
                    df_full, test_table, test_cols, gate_rmsd=False)
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
                "pairs_with_valid_pose_%": round(100 * pv / pt, 1) if pt else 0.0,
                "poses_total": len(g), "poses_valid": int(g["pb_valid"].sum()),
                "poses_valid_%": round(100 * g["pb_valid"].mean(), 1) if len(g) else 0.0,
            })
        attr = pd.DataFrame(attr_rows).set_index("method")
        attr.to_csv(valid_dir / "pb_valid_attrition.csv")
        print("\nPB-valid attrition by method "
              "(denominator of the pb_valid/ figures = pairs_with_valid_pose):")
        print(attr.to_string())

        # Headings are intentionally NOT suffixed with a "PB-valid only" tag —
        # the pb_valid/ figures are distinguished by their output folder and by
        # their (lower) per-method n= annotations, not by the title.
        _PLOT_TITLE_SUFFIX = ""
        try:
            oracle_sum_v = aggregate_oracle(df_valid)
            oracle_sum_v.to_csv(valid_dir / "oracle_summary.csv")
            has_ranking_valid = bool(df_valid["method"].isin(RANKING_TOOLS).any())
            top1_sum_v = aggregate_top1(df_valid) if has_ranking_valid else pd.DataFrame()
            if not top1_sum_v.empty:
                top1_sum_v.to_csv(valid_dir / "top1_summary.csv")
            rank_df_v = aggregate_by_rank(df_valid, args.top_n)
            ifp_rank_df_v = aggregate_ifp_by_rank(df_valid, args.top_n)
            if not rank_df_v.empty:
                rank_df_v.to_csv(valid_dir / "per_rank_metrics.csv", index=False)
                ifp_rank_df_v.to_csv(valid_dir / "per_rank_ifp_recovery.csv", index=False)

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
            aggregate_twist_turn(df_valid, args.top_n).to_csv(
                valid_dir / "twist_turn_summary.csv", index=False)
            plot_twist_turn(df_valid, valid_dir / "14_twist_turn.png", top_n=args.top_n,
                            pb_valid_only=True)
        finally:
            _PLOT_TITLE_SUFFIX = ""
        print(f"\nPB-valid figures + CSVs → {valid_dir.resolve()}")

    print(f"\nAll figures + CSVs → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
