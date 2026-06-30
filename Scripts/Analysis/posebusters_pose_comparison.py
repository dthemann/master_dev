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

Best-EquiBind filter (--best-equibind-only)
    Collapse the many EquiBind variants down to the single best-performing one
    (highest oracle_rmsd_le_2.0A_%) so every summary/plot compares AutoDock,
    DiffDock and just one EquiBind curve. The retained variant is relabelled
    "EquiBind*" in every legend/axis. Off by default; AutoDock/DiffDock are never
    touched and per_pose_metrics.csv still carries all variants.

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
    09_accuracy_vs_validity_bars.png — PoseBusters-paper Fig.1-style bars:
                                       %RMSD ≤ 2 Å vs %(RMSD ≤ 2 Å & PB-valid),
                                       Top-1 and Oracle panels
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
                                       cumulative recovery within top-k (ranking tools)
    16_pocket_localization.png       — Pocket targeting (ranking tools), two panels:
                                       (A) per rank k=1..N, the % of those rank-k poses
                                       in the validated (crystal) pocket plus the
                                       cumulative % of complexes with any of ranks 1..k
                                       in it; (B) per-tool summary — % of all top-N
                                       RANKED poses in the validated pocket beside the %
                                       of complexes whose ORACLE pose is in a DIFFERENT
                                       pocket than the top-N ranked set (a pose counts
                                       as in-pocket when its centroid is within
                                       --pocket-cutoff Å of the crystal ligand centroid)
    17_pb_test_waterfall.png         — PoseBusters paper-style failure waterfall, one
                                       panel per method: starting from every top-ranked
                                       pose (top-1 for ranking tools, best/oracle for
                                       unranked EquiBind), drop RMSD > 2 Å then each
                                       canonical PoseBusters test in turn — each red bar
                                       is the poses that test removes; green = passing all

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
    oracle_rank_distribution.csv — Per ranking tool × rank k: pct_at_rank (% of
                                complexes whose oracle pose sits at rank k) + cum_pct
    pocket_localization.csv   — Per ranking tool: pct_topN_in_validated_pocket,
                                pct_oracle_diff_pocket_vs_topN (+ top1/any-topN/oracle
                                in-pocket context numbers, n_complexes, pocket_cutoff_A)
    pb_test_waterfall.csv     — Per method × cascade step: poses_removed / remaining /
                                remaining_pct for the PoseBusters failure waterfall

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
        elif t in ("raw", "smina"):
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
        parts.append("smina-opt" if refine == "smina" else "raw")
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
    if method == "diffdock":
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
    pose_conf = pose.GetConformer()
    ref_conf = ref.GetConformer()
    ref_coords = np.array([list(ref_conf.GetAtomPosition(i))
                           for i in range(ref.GetNumAtoms())])

    matches = pose.GetSubstructMatches(ref, uniquify=False, useChirality=False)
    if not matches:
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return float("nan")
        matches = [tuple(range(ref.GetNumAtoms()))]

    best = math.inf
    for mp_ in matches:
        pose_coords = np.array([list(pose_conf.GetAtomPosition(j)) for j in mp_])
        diff = pose_coords - ref_coords
        rmsd = float(np.sqrt((diff * diff).sum() / len(ref_coords)))
        if rmsd < best:
            best = rmsd
    return best if best != math.inf else float("nan")


def centroid_distance(pose: Chem.Mol, ref: Chem.Mol) -> float:
    pose = Chem.RemoveHs(pose); ref = Chem.RemoveHs(ref)
    pc = pose.GetConformer(); rc = ref.GetConformer()
    a = np.mean([list(pc.GetAtomPosition(i)) for i in range(pose.GetNumAtoms())], axis=0)
    b = np.mean([list(rc.GetAtomPosition(i)) for i in range(ref.GetNumAtoms())], axis=0)
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
    pc = pose.GetConformer(); rc = ref.GetConformer()
    ref_coords = np.array([list(rc.GetAtomPosition(i))
                           for i in range(ref.GetNumAtoms())])
    matches = pose.GetSubstructMatches(ref, uniquify=False, useChirality=False)
    if not matches:
        if pose.GetNumAtoms() != ref.GetNumAtoms():
            return float("nan"), float("nan")
        matches = [tuple(range(ref.GetNumAtoms()))]

    qc = ref_coords - ref_coords.mean(axis=0)
    best_placed, best_angle = math.inf, float("nan")
    best_fit = math.inf
    for mp_ in matches:
        P = np.array([list(pc.GetAtomPosition(j)) for j in mp_])
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
    pc = pose.GetConformer()
    lig_xyz = np.array([list(pc.GetAtomPosition(i)) for i in range(pose.GetNumAtoms())])
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


def _select_best_equibind(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind methods plus only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest oracle RMSD ≤ 2 Å success rate
    (``oracle_rmsd_le_2.0A_%`` from the oracle summary). Returns the filtered
    frame and the chosen variant's method key (``None`` if no EquiBind variant is
    present, in which case ``df`` is returned unchanged). AutoDock/DiffDock rows
    are always retained.
    """
    methods = df["method"].astype(str)
    eq_mask = methods.str.startswith("equibind")
    if not eq_mask.any():
        return df, None

    oracle_sum = aggregate_oracle(df)
    col = "oracle_rmsd_le_2.0A_%"
    eq_keys = [m for m in oracle_sum.index if str(m).startswith("equibind")]
    scores = (oracle_sum.loc[eq_keys, col].astype(float).dropna()
              if col in oracle_sum.columns else pd.Series(dtype=float))
    best_variant = (str(scores.idxmax()) if not scores.empty
                    else sorted(methods[eq_mask].unique())[0])

    keep = (~eq_mask) | (methods == best_variant)
    return df[keep].reset_index(drop=True), best_variant


def _select_best_diffdock(df: pd.DataFrame) -> tuple[pd.DataFrame, str | None]:
    """Keep non-DiffDock methods plus only the single best DiffDock variant.

    "Best" = the DiffDock optimizer variant (diffdock / diffdock_smina /
    diffdock_gnina) with the highest oracle RMSD ≤ 2 Å success rate
    (``oracle_rmsd_le_2.0A_%``). The winner is relabelled to the canonical
    ``diffdock`` method so it stays a ranking tool; the original variant key is
    returned (for the 'DiffDock*' label / print). AutoDock/EquiBind are retained.
    ``None`` when no DiffDock variant is present (``df`` unchanged)."""
    methods = df["method"].astype(str)
    dd_mask = methods.str.startswith("diffdock")
    if not dd_mask.any():
        return df, None

    oracle_sum = aggregate_oracle(df)
    col = "oracle_rmsd_le_2.0A_%"
    dd_keys = [m for m in oracle_sum.index if str(m).startswith("diffdock")]
    scores = (oracle_sum.loc[dd_keys, col].astype(float).dropna()
              if col in oracle_sum.columns else pd.Series(dtype=float))
    best_variant = (str(scores.idxmax()) if not scores.empty
                    else sorted(methods[dd_mask].unique())[0])

    keep = (~dd_mask) | (methods == best_variant)
    out = df[keep].reset_index(drop=True)
    # Canonicalise the winner to 'diffdock' so it ranks/colours like the tool.
    out.loc[out["method"] == best_variant, "method"] = "diffdock"
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


def plot_accuracy_validity_bars(oracle_sum: pd.DataFrame,
                                top1_sum: pd.DataFrame,
                                out: Path) -> None:
    """PoseBusters-paper Fig. 1-style accuracy-vs-validity bars.

    For every method we draw, per panel, one bar split into two overlaid parts:
      * the full (light) bar    — % of complexes with a pose RMSD ≤ 2 Å
      * the solid (dark) bar    — the subset that is ALSO PB-valid
                                  (RMSD ≤ 2 Å **and** passes every PoseBusters check)
    Because the PB-valid set is a subset of the RMSD-≤-2-Å set, the dark bar is
    always contained in the light one; the visible light "cap" above it is the
    fraction of accurate-but-physically-invalid predictions — the paper's central
    point. Two panels:
      * Top-1   — the tool's rank-1 pose (ranking tools only: Vina, DiffDock)
      * Oracle  — the best-of-N pose per complex (all tools, incl. EquiBind)

    The numbers come straight from the aggregation tables
    (``{top1,oracle}_rmsd_le_2.0A_%`` and ``{top1,oracle}_pb_valid_and_rmsd2_%``).
    """
    panels = [
        ("Top-1 (tool's rank-1 pose)", top1_sum,
         "top1_rmsd_le_2.0A_%", "top1_pb_valid_and_rmsd2_%"),
        ("Oracle (best of all docked poses)", oracle_sum,
         "oracle_rmsd_le_2.0A_%", "oracle_pb_valid_and_rmsd2_%"),
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    for ax, (title, summ, rmsd_col, valid_col) in zip(axes, panels):
        if summ is None or summ.empty or rmsd_col not in summ.columns:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, color="grey", style="italic")
            ax.set_title(title, fontweight="bold")
            ax.set_xticks([])
            continue

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
            [f"{TOOL_LABEL.get(m, m)}\n{_n_annot(summ, m, poses=rmsd_col.startswith('oracle'))}"
             for m in methods],
            rotation=25, ha="right", rotation_mode="anchor")
        ax.set_title(title, fontweight="bold")
        ax.set_ylim(0, 105)
        ax.grid(axis="y", alpha=0.3)

    axes[0].set_ylabel("% of receptor-ligand complexes")

    # Method-agnostic legend (grey swatches): the two bar layers.
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor="#888888", alpha=0.32, edgecolor="#888888",
              label="RMSD ≤ 2 Å (accurate)"),
        Patch(facecolor="#555555", edgecolor="white",
              label="RMSD ≤ 2 Å & PB-valid (accurate + PoseBuster valid)"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 1.02), fontsize=10)
    _label_panels(axes)
    fig.suptitle("Docking accuracy vs. PoseBuster validity — PoseBusters Benchmark",
                 fontsize=14, fontweight="bold", y=1.10)
    fig.tight_layout()
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


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


def plot_pocket_localization(summ: pd.DataFrame, by_rank: pd.DataFrame, top_n: int,
                             pocket_cutoff: float, out: Path) -> None:
    """Pocket localization for the ranking tools, in two panels.

    (A) Per-rank targeting — for each rank position k = 1..N, the % of those rank-k
        poses that land in the validated pocket (solid), and the cumulative % of
        complexes with ANY of ranks 1..k in the pocket (dashed). Shows how pocket
        targeting decays down the ranked list and how much keeping more poses buys.
    (B) Per-tool summary dumbbell — % of all top-N poses in the validated pocket
        (filled) vs % of complexes whose oracle pose sits in a DIFFERENT pocket
        than the top-N ranked set (open: a ranking failure).

    "In the validated pocket" = pose centroid ≤ ``pocket_cutoff`` Å from the
    crystal ligand centroid.
    """
    if summ is None or summ.empty:
        return
    methods = list(summ.index)
    colors = [TOOL_COLORS.get(m, "grey") for m in methods]

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.8),
                             gridspec_kw={"width_ratios": [1.45, 1.0]})

    # ── Panel A: per-rank pocket targeting (the requested per-top-pose view) ──
    axA = axes[0]
    if by_rank is not None and not by_rank.empty:
        for m in methods:
            sub = by_rank[by_rank["method"] == m].sort_values("rank")
            if sub.empty:
                continue
            c = TOOL_COLORS.get(m, "grey")
            axA.plot(sub["rank"], sub["pct_in_pocket"], marker="o", lw=2, color=c,
                     markeredgecolor="black", markeredgewidth=0.5,
                     label=f"{TOOL_LABEL.get(m, m)} — rank-k pose")
            axA.plot(sub["rank"], sub["pct_any_in_top_k"], marker="s", lw=1.6,
                     ls="--", color=c, alpha=0.85, markeredgecolor="black",
                     markeredgewidth=0.4,
                     label=f"{TOOL_LABEL.get(m, m)} — any of top-1..k")
    axA.set_xticks(range(1, top_n + 1))
    axA.set_xlabel("Rank position k  (1 = top-ranked pose)")
    axA.set_ylabel("% in the validated pocket")
    axA.set_ylim(0, 105)
    axA.set_title("Per-rank pocket targeting\n(solid = the k-th pose alone · "
                  "dashed = any of the first k)", fontsize=10, fontweight="bold")
    axA.grid(alpha=0.3); axA.set_axisbelow(True)
    axA.legend(fontsize=7.5, loc="lower left")

    # ── Panel B: per-tool summary dumbbell (the original complex-level view) ──
    v1 = [float(summ.loc[m, "pct_topN_in_validated_pocket"]) for m in methods]
    v2 = [float(summ.loc[m, "pct_oracle_diff_pocket_vs_topN"]) for m in methods]
    labels = [f"{TOOL_LABEL.get(m, m)}\n(n={int(summ.loc[m, 'n_complexes'])} complexes · "
              f"{int(summ.loc[m, 'n_topN_poses'])} top-{top_n} poses)" for m in methods]
    handles = _dumbbell(axes[1], labels, v1, v2, colors,
                        left_name=f"All top-{top_n} ranked poses in the validated "
                                  f"pocket (% of poses)",
                        right_name=f"Oracle pose in a different pocket than the "
                                   f"top-{top_n} ranked (% of complexes)",
                        value_fmt="{:.1f}%")
    axes[1].set_xlabel("Percentage")
    axes[1].set_title("Per-tool summary", fontsize=10, fontweight="bold")
    axes[1].legend(handles=handles, fontsize=7.5, loc="lower right")

    _label_panels(axes)
    fig.suptitle(_vt("Pocket localization — does the ranking target the validated "
                     "pocket?\n"
                     f"(in-pocket = pose centroid ≤ {pocket_cutoff:g} Å from the "
                     f"crystal ligand centroid)"),
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


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
                           test_cols: list[str], top_n: int) -> dict:
    """Build, per method, the sequential PoseBusters filter cascade.

    Mirrors the PoseBusters paper's waterfall: start from every representative
    pose (top-1 ranked pose for ranking tools; best/oracle pose for unranked
    tools such as EquiBind), drop those with RMSD > 2 Å, then apply each canonical
    PoseBusters test IN ORDER, removing at each step the poses that fail it AND
    passed every previous filter. The remainder at the end pass everything.

    EquiBind variants are collapsed to the single best one (``_select_best_equibind``)
    so the figure stays to a few panels. Returns ``{method: {N, selection, steps}}``
    where each step is a dict (label, kind, removed, remaining, before, after);
    *kind* ∈ {start, drop, milestone, end} and before/after are running counts.
    """
    wdf, _ = _select_best_equibind(df)
    top1 = _top1_per_pair(wdf)
    oracle = _oracle_per_pair(wdf)
    present = set(wdf["method"].astype(str))
    methods = [m for m in ("autodock", "diffdock") if m in present]
    methods += [m for m in sorted(present) if m not in RANKING_TOOLS]

    cascades: dict = {}
    for method in methods:
        if method in RANKING_TOOLS:
            rep, selection = top1[top1["method"] == method], "top-1 ranked pose"
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


def plot_pb_waterfall(cascades: dict, out: Path) -> None:
    """Small-multiples PoseBusters waterfall — one panel per method, stacked so
    the long test labels are shared. Each red bar is the poses removed by that
    test (annotated −k); the teal milestone bar is the RMSD ≤ 2 Å running total
    and the green bar the poses passing every test.
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
    fig.suptitle(_vt("PoseBusters test-failure waterfall — top-ranked poses per "
                     "method\n(sequential filter: each red bar = poses removed by "
                     "that test; green = poses passing every test)"),
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)


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
    refine : "smina" (__refSMINA) | "raw" (__refRAW) | None (unsuffixed)
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
    if refine not in ("smina", "raw"):
        refine = ("smina" if "__refsmina" in name
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
                         "(highest oracle_rmsd_le_2.0A_%%) in all summaries and "
                         "plots, relabelled 'EquiBind*'. AutoDock/DiffDock are "
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
                         "variant (raw/smina/gnina, highest oracle_rmsd_le_2.0A_%%) "
                         "in all summaries/plots, relabelled 'DiffDock*'. Scores all "
                         "three then picks the best, so it OVERRIDES --diffdock-variant. "
                         "AutoDock/EquiBind are unaffected.")
    args = ap.parse_args()

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
            with mp.Pool(args.workers) as pool:
                for i, batch in enumerate(
                        pool.imap_unordered(process_pair, work, chunksize=4), 1):
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

    # ── Optional: restrict the report to the single best EquiBind variant ──
    # Done after the full per-pose dump (which keeps every variant) so only the
    # summaries and plots below are filtered.
    if args.best_equibind_only:
        df, best_variant = _select_best_equibind(df)
        if best_variant:
            _LABEL_OVERRIDES[best_variant] = "EquiBind*"
            print(f"best-equibind-only: '{best_variant}' is the top EquiBind variant "
                  f"by oracle_rmsd_le_2.0A_% — keeping only it (shown as 'EquiBind*').")
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
                  f"oracle_rmsd_le_2.0A_% — keeping only it as 'diffdock' (shown as 'DiffDock*').")
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
    plot_accuracy_validity_bars(oracle_sum, top1_sum,
                                args.out_dir / "09_accuracy_vs_validity_bars.png")
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
        plot_pocket_localization(pocket_sum, pocket_by_rank, args.top_n,
                                 args.pocket_cutoff,
                                 args.out_dir / "16_pocket_localization.png")
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
    # in 09_accuracy_vs_validity_bars.png and the oracle/top1 summary CSVs.
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
