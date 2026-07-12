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
      refine   __refSMINA -> "smina" (smina-optimized);  __refRAW -> "raw" control
      clamp    __clampON / __clampOFF                    (centroid-clamp study)
    e.g. fpocket_raw_p001_pose01__refSMINA.sdf -> "equibind_fpocket_smina". Axes
    not present collapse to "equibind_fpocket" / "equibind_p2rank" /
    "equibind_unguided". All EquiBind variants stay oracle-only (Part B is still
    AutoDock/DiffDock, which alone rank poses). The per-pose provenance
    (pocket_source, clamp_variant, refine_variant, smina_affinity) is also carried
    into per_pose_metrics.csv. NOTE: refine_mode='on' rewrites the pose in place
    with no suffix/tag, so the smina split is only visible for refine_mode='both'.

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
    05_cumulative_oracle_curve.png   — Best-of-top-k success as k grows
    06_top1_vs_oracle_scatter.png    — Top-1 RMSD vs oracle RMSD per pair
    07_oracle_rank_histogram.png     — Which rank does the oracle pose have?
    08_cross_tool_oracle_pairwise.png — Cross-tool oracle RMSD scatter matrix
    09a_accuracy_vs_validity_top1.png — PoseBusters-paper Fig.1-style bars for the
                                       TOP-1 (rank-1) pose (ranking tools): %RMSD ≤ 2 Å
                                       vs %(RMSD ≤ 2 Å & PB-valid)
    09b_accuracy_vs_validity_oracle.png — same bars for the ORACLE (best-of-N) pose
                                       (all tools incl. EquiBind)
    09c_optimization_raw_vs_best.png — grouped bars contrasting each tool's RAW
                                       representative pose with its BEST post-hoc-
                                       optimised variant (Vina: rank-1, no opt step;
                                       DiffDock: raw vs smina rank-1; EquiBind: first
                                       generated pose vs smina-ranked) — shows how
                                       smina optimisation recovers PB-validity
    09d_rank1_vs_topn.png            — grouped bars: best-of-top-d pose at several
                                       depths (top-1, top-15, top-30) for BOTH the raw
                                       and the smina-optimised run of each ML tool
                                       (ranking headroom + optimisation gain in one
                                       view; top-30 ≈ full oracle). Vina (confidence
                                       rank), DiffDock raw & smina-opt (confidence rank),
                                       EquiBind raw (generation order) & smina-opt
                                       (smina-affinity rank). Shows more poses raise
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
                                       made explicit: per variant a solid (near-native,
                                       RMSD ≤ 2 Å) and dashed (near-native AND PB-valid)
                                       line — AutoDock, DiffDock raw & smina-opt, EquiBind
                                       raw (generation order) & smina-opt (smina affinity);
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
                                       passing all tests
    18_topn_within_thresholds.png    — Fine-grained RMSD sweep: how many of the top-N
                                       ranked poses of AutoDock Vina / DiffDock (native
                                       rank) / EquiBind (smina-affinity ranked) land
                                       within 1, 1.25, 1.5 … Å of the crystal
                                       (--fine-rmsd-thresholds). Two panels — (A) % of
                                       complexes whose rank-1 pose (solid) / any of the
                                       top-N poses (best-of-top-N, dashed) is within each
                                       threshold; (B) % of ALL pooled top-N ranked poses
                                       within each threshold. Denser than the 2 Å success
                                       line so the sub-2 Å accuracy of the top poses shows.
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

    pb_valid/                        — PB-valid variant of every RMSD figure (01–08)
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
                "equibind_unguided_smina", "equibind_fpocket_smina", "equibind_p2rank_smina",
                "equibind_unguided_raw", "equibind_fpocket_raw", "equibind_p2rank_raw"]

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


def aggregate_topn_within_thresholds(
        df: pd.DataFrame, top_n: int,
        thresholds: tuple[float, ...] = FINE_RMSD_THRESHOLDS,
        eq_df: "pd.DataFrame | None" = None) -> pd.DataFrame:
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
    smina-optimised EquiBind poses (which have no native ranking) — they are ranked
    by smina AFFINITY (most-negative = rank 1, ``_smina_affinity_rank``) and appended
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
        n_pairs = sub.groupby(["protein", "ligand"]).ngroups
        ranked = sub[sub["_rk"] <= top_n]
        top1 = (ranked[ranked["_rk"] == 1]
                .groupby(["protein", "ligand"]).first().reset_index())
        top1_rmsd = top1["rmsd"].dropna()
        rk_valid = ranked.dropna(subset=["rmsd"])
        best_topn = (rk_valid.groupby(["protein", "ligand"])["rmsd"].min()
                     if len(rk_valid) else pd.Series(dtype=float))
        pose_rmsd = ranked["rmsd"].dropna()
        n_poses = len(pose_rmsd)
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

    # EquiBind: no native ranking, so rank by smina affinity (smina-opt variant).
    if eq_df is not None and not eq_df.empty:
        _emit("equibind", eq_df, _smina_affinity_rank(eq_df))

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
                                  poses: bool) -> None:
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

    ax.set_xticks(x)
    ax.set_xticklabels(
        [f"{TOOL_LABEL.get(m, m)}\n{_n_annot(summ, m, poses=poses)}" for m in methods],
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


def plot_accuracy_validity_oracle(oracle_sum: pd.DataFrame, out: Path) -> None:
    """Fig 09b — accuracy-vs-validity bars for the ORACLE (best-of-N) pose.

    All tools including EquiBind: the light bar is the % of complexes whose best
    docked pose is RMSD ≤ 2 Å, the dark subset the fraction that is also PB-valid.
    Single panel (split out of the former two-panel fig 09).
    """
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    _draw_accuracy_validity_panel(
        ax, "Oracle (best of all docked poses)", oracle_sum,
        "oracle_rmsd_le_2.0A_%", "oracle_pb_valid_and_rmsd2_%", poses=True)
    ax.set_ylabel("% of receptor-ligand complexes")
    _accuracy_validity_legend(fig)
    fig.suptitle("Docking accuracy vs. PoseBuster validity — oracle (best-of-N) pose\n"
                 "PoseBusters Benchmark", fontsize=13, fontweight="bold", y=1.09)
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
      * EquiBind best  — the smina-optimised blind run (``equibind_unguided_smina``),
                         RANKED by smina affinity (best / most-negative = the pick);
                         smina supplies both refined geometry AND a ranking.

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

    def _smina_ranked(sub: pd.DataFrame) -> pd.DataFrame:
        sub = sub.copy()
        sub["_aff"] = pd.to_numeric(sub["smina_affinity"], errors="coerce")
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
        ("equibind", "smina-opt\n(smina-ranked)", True,  "equibind_unguided_smina", _smina_ranked),
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
             "vs. smina-optimized, ranked by smina affinity.",
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


def _smina_affinity_rank(sub: pd.DataFrame) -> pd.Series:
    """Per-complex 1-based rank by smina affinity (most-negative = rank 1). NaN
    affinities sort last, broken by generation index, so every pose gets a rank —
    this is the ranking EquiBind gains from smina optimization. Indexed like ``sub``."""
    sub = sub.copy()
    sub["_gi"] = sub["pose_name"].map(_equibind_pose_index)
    sub["_aff"] = pd.to_numeric(sub["smina_affinity"], errors="coerce")
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
    ("equibind", "EquiBind (smina-opt)", "equibind_unguided_smina", "smina-affinity rank", "smina"),
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
      * EquiBind (smina-opt) — ranked by smina AFFINITY (the ranking smina supplies).

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
        else:  # smina
            rk = _smina_affinity_rank(sub)
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
             "raw — generation order (no ranking); EquiBind smina-opt — smina affinity. "
             "'raw' = un-optimised geometry, 'smina-opt' = smina-refined.",
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
    "EquiBind (smina-opt)": "#1b7837",
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
                                thresholds: tuple[float, ...], out: Path) -> None:
    """How many of the top-ranked poses land within a fine RMSD grid.

    AutoDock Vina / DiffDock (native rank) + EquiBind (smina-affinity ranked), two
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

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.5, 5.8))
    for method in methods:
        sub = within_df[within_df["method"] == method].sort_values("rmsd_threshold_A")
        if sub.empty:
            continue
        color = TOOL_COLORS.get(method, "grey")
        label = TOOL_LABEL.get(method, method)
        if method == "equibind":
            label += " (smina-ranked)"
        n_pairs = int(sub["n_pairs"].iloc[0])
        n_poses = int(sub["n_poses_topN"].iloc[0])
        axA.plot(sub["rmsd_threshold_A"], sub["top1_within_%"],
                 marker="o", lw=2, color=color,
                 label=f"{label} — top-1 (n={n_pairs})")
        axA.plot(sub["rmsd_threshold_A"], sub["best_topN_within_%"],
                 marker="s", ls="--", lw=2, color=color, alpha=0.75,
                 label=f"{label} — best of top-{top_n}")
        axB.plot(sub["rmsd_threshold_A"], sub["pose_within_%"],
                 marker="o", lw=2, color=color,
                 label=f"{label} (≤{n_poses} poses)")

    for ax in (axA, axB):
        ax.axvline(2.0, color="grey", ls=":", alpha=0.7, lw=1.2)
        ax.set_xticks(thr)
        ax.set_xticklabels([f"{t:g}" for t in thr], rotation=45, ha="right")
        ax.set_xlabel("RMSD threshold vs crystal ligand (Å)")
        ax.set_ylim(0, 100)
        ax.grid(alpha=0.3)
    axA.set_ylabel("% of complexes with a pose within the threshold")
    axA.set_title("Top-ranked pose accuracy vs distance threshold\n"
                  "(rank-1 solid, best-of-top-N dashed)")
    axA.legend(fontsize=8)
    axB.set_ylabel("% of the top-N ranked poses within the threshold")
    axB.set_title(f"How many of the top-{top_n} ranked poses are within t Å")
    axB.legend(fontsize=8)

    _label_panels([axA, axB])
    fig.suptitle(_vt(f"Top-ranked pose accuracy across {thr[0]:g}–{thr[-1]:g} Å "
                     f"(AutoDock Vina, DiffDock, EquiBind; top-{top_n})"),
                 fontsize=13, fontweight="bold")
    footnote = "\n".join([
        "Denominators — (A) top-1 & best-of-top-N are % of each tool's complexes (n pairs); (B) is % of that tool's pooled rank ≤ N poses.",
        "Best-of-top-N = the closest of the tool's first N ranked poses per complex (the oracle within the ranked set). Rank: AutoDock/DiffDock",
        "native; EquiBind has none, so it is ranked by smina affinity. Dotted line = 2 Å. RMSD = symmetry-corrected heavy-atom, no superposition.",
    ])
    fig.text(0.5, 0.015, footnote, ha="center", va="bottom", fontsize=7.5,
             color="0.30", linespacing=1.35)
    fig.tight_layout(rect=(0, 0.11, 1, 0.95))
    fig.savefig(out, dpi=160); plt.close(fig)


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


def aggregate_twist_turn(df: pd.DataFrame) -> pd.DataFrame:
    """Per-method median & mean of every twist/turn measure, for top-1 & oracle.

    Long-form table (one row per selection × method). The companion plot draws
    the full distributions; this CSV gives the central tendencies for the text.
    """
    sels = {"oracle": _oracle_per_pair(df), "top1": _top1_per_pair(df)}
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


def plot_twist_turn(df: pd.DataFrame, out: Path) -> None:
    """Distributions of how each tool re-orients (turns) and re-bends (twists)
    the ligand relative to its crystal pose — top-1 vs oracle, per measure.

    One panel per measure; within a panel, each tool gets a solid 'oracle' box
    and (ranking tools only) a hatched 'top-1' box. Outliers hidden so the
    bulk of each distribution stays readable.
    """
    from matplotlib.patches import Patch
    oracle = _oracle_per_pair(df)
    top1 = _top1_per_pair(df)
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
            t = (top1.loc[top1["method"] == method, col].dropna()
                 if (method in RANKING_TOOLS and col in top1) else pd.Series(dtype=float))
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
        Patch(facecolor="#888888", alpha=0.85, label="Oracle (best of all docked poses)"),
        Patch(facecolor="#888888", alpha=0.40, hatch="//", label="Top-1 (tool's rank-1 pose)"),
    ], loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.015), fontsize=10)
    fig.suptitle(_vt("How the ligand is twisted & turned vs its crystal pose\n"
                     "(moved = translation · turned = rigid-body rotation · "
                     "twisted = internal conformation)"),
                 fontsize=13, fontweight="bold", y=1.07)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
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


# The two oracle selections the form-fidelity figures are rendered under (both by
# default). The note is stamped into each figure title so the graphs are
# self-identifying; plot_oracle_selection_comparison puts them head-to-head.
_SEL_NOTE_NEAREST = ("selection: RMSD-greedy oracle — the single nearest pose, "
                     "kept only if it is itself ≤ 2 Å AND PB-valid")
_SEL_NOTE_VALID = ("selection: validity-constrained oracle — the NEAREST pose that "
                   "is ≤ 2 Å AND PB-valid (more generous sampling ceiling)")


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
    ("EquiBind (smina-opt)", "equibind_unguided_smina", "smina"),
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
              else _smina_affinity_rank(sub))
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
    15 (rank of the min-RMSD pose, RMSD only), this brings PB-validity into the view."""
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
    ax.set_xlabel("Top-k — within the tool's first k ranked poses")
    ax.set_ylabel("% of receptor-ligand complexes recovered")
    ax.set_ylim(0, 100)
    ax.set_xlim(1, int(rec["k"].max()) + 6)
    ax.set_xticks([k for k in (1, 5, 10, 15, 20, 25, 30) if k <= int(rec["k"].max())])
    ax.grid(alpha=0.3); ax.set_axisbelow(True)
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
        rep = rep[["method", "protein", "ligand", "pose_name", "rmsd"]]
        if rep.empty:
            continue
        merged = rep.merge(test_table, how="left",
                           on=["method", "protein", "ligand", "pose_name"])
        merged = merged.drop_duplicates(
            subset=["method", "protein", "ligand", "pose_name"])
        N = len(merged)
        if N == 0:
            continue

        steps = [{"label": "All predictions", "kind": "start",
                  "removed": 0, "remaining": N, "before": N, "after": N}]
        surv = pd.Series(True, index=merged.index)

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
        cascades[method] = {"N": N, "selection": selection, "steps": steps}
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
                ax.text(xi, a + 1, f"{s['remaining']}\n({a:.0f}%)", ha="center",
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
                     f"(n={N}; passing all tests = {final_pct:.0f}%)",
                     fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)

    axes[-1].set_xticks(np.arange(nx))
    axes[-1].set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
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


def plot_pb_waterfall_top1_vs_topn(cascades_top1: dict, cascades_topn: dict,
                                   top_n: int, out: Path,
                                   csv_out: Path | None = None) -> None:
    """Compare the PB cascade survival for top-1 vs best-of-top-N poses.

    One panel per RANKING tool (EquiBind has no ranking, so it has no top-N pick
    to compare). Each panel overlays two survival curves — the % of complexes
    whose representative pose still survives after each sequential filter — for
    the rank-1 pose (solid grey) and the best-of-top-N pose (dashed green). The
    shaded band between them is the recovery from letting the ranker offer more
    poses; Δ callouts mark the gain at "RMSD ≤ 2 Å" and at "Passing all tests".
    """
    methods = [m for m in ("autodock", "diffdock")
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
        _, s1 = _waterfall_survival(c1)
        ln, sn = _waterfall_survival(cn)
        ax.fill_between(x, s1, sn, step="mid", color=green, alpha=0.12, zorder=1)
        ax.step(x, s1, where="mid", color=grey, lw=1.8, label="rank-1 pose",
                zorder=3)
        ax.step(x, sn, where="mid", color=green, lw=1.8, ls="--",
                label=f"best of top-{top_n}", zorder=3)
        ax.plot(x, s1, "o", color=grey, ms=3.5, zorder=4)
        ax.plot(x, sn, "D", color=green, ms=3.5, zorder=4)

        for name in ("RMSD ≤ 2 Å", "Passing all tests"):
            i = _idx(labels, name)
            if i is None:
                continue
            d = sn[i] - s1[i]
            ax.annotate(f"+{d:.0f} pp", xy=(i, (s1[i] + sn[i]) / 2),
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
        ax.set_title(f"{TOOL_LABEL.get(method, method)} — passing all tests: "
                     f"{end1:.0f}% (rank-1) → {endn:.0f}% (top-{top_n}), "
                     f"+{endn - end1:.0f} pp",
                     fontsize=10, fontweight="bold")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8, loc="upper right", framealpha=0.9)

    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels(labels, rotation=90, ha="center", fontsize=8)
    if nrows > 1:
        _label_panels(axes)
    fig.suptitle(_vt(f"Ranking headroom — rank-1 vs best-of-top-{top_n} pose "
                     "through the PoseBusters cascade\n(green band = complexes "
                     "recovered by considering more ranked poses)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
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
                      "smina_affinity", "pocket_id")


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
        "schema": 3,
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
                    help="Split EquiBind into fpocket/p2rank/unguided x smina/raw "
                         "x clamp variants from the CSV provenance columns "
                         "(filename fallback) (default: on).")
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
                                  args.out_dir / "09b_accuracy_vs_validity_oracle.png")
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
    # smina-optimised run of DiffDock/EquiBind are present under their own method keys.
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
    aggregate_twist_turn(df).to_csv(args.out_dir / "twist_turn_summary.csv", index=False)
    plot_twist_turn(df, args.out_dir / "14_twist_turn.png")
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
    # Head-to-head of the two selections: table (how many complexes each rule
    # keeps + how many the nearest rule needlessly drops) and figure (20c).
    aggregate_oracle_selection_comparison(df).to_csv(
        args.out_dir / "oracle_selection_comparison.csv", index=False)
    plot_oracle_selection_comparison(
        df, args.out_dir / "20c_oracle_selection_comparison.png", args.form_ok_kabsch)
    print("  wrote form-fidelity summaries (nearest + valid_ceiling) + "
          "oracle_selection_comparison.csv (+ 20/20b/20c/20d/20e figures, both selections)")

    # ── Part B: ranking quality (Vina + DiffDock only) ───────────
    if not rank_df.empty:
        plot_rank_success_curve(rank_df, args.top_n,
                                args.out_dir / "04_rank_success_curve.png")
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

        # How many of the top-N ranked poses land within 1, 1.25, 1.5 … Å?
        # EquiBind (smina-ranked) is added from df_full since the collapsed df keeps
        # only the best EquiBind variant (which may not be the smina one).
        eq_smina = df_full[df_full["method"].astype(str) == "equibind_unguided_smina"]
        within_df = aggregate_topn_within_thresholds(
            df, args.top_n, args.fine_rmsd_thresholds,
            eq_df=eq_smina if not eq_smina.empty else None)
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
                plot_pb_waterfall_top1_vs_topn(
                    cascades, cascades_topn, args.top_n,
                    args.out_dir / "17c_pb_waterfall_top1_vs_topn.png",
                    csv_out=args.out_dir / "pb_waterfall_top1_vs_topn.csv")
                print("  wrote top-1 vs top-N recovery → "
                      "pb_waterfall_top1_vs_topn.csv")
        else:
            print("  Skipping PB waterfall (no representative poses to cascade).")
    else:
        print("  Skipping PB waterfall (raw CSV lacks per-test columns).")

    # ── PB-valid variants ────────────────────────────────────────
    # Re-emit every RMSD figure (01–08) and its summary CSVs using ONLY poses that
    # pass all canonical PoseBusters tests (pb_valid). Each graph then describes
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
            if not rank_df_v.empty:
                rank_df_v.to_csv(valid_dir / "per_rank_metrics.csv", index=False)

            print("\nGenerating PB-valid plot variants …")
            plot_oracle_rmsd_cdf(df_valid, valid_dir / "01_oracle_rmsd_cdf.png")
            plot_oracle_rmsd_box(df_valid, valid_dir / "02_oracle_rmsd_boxplot.png")
            plot_oracle_vs_top1_success(oracle_sum_v, top1_sum_v,
                                        valid_dir / "03_oracle_vs_top1_success.png")
            if not rank_df_v.empty:
                plot_rank_success_curve(rank_df_v, args.top_n,
                                        valid_dir / "04_rank_success_curve.png")
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
        finally:
            _PLOT_TITLE_SUFFIX = ""
        print(f"\nPB-valid figures + CSVs → {valid_dir.resolve()}")

    print(f"\nAll figures + CSVs → {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
