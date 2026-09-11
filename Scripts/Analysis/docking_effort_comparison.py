#!/usr/bin/env python3
"""Compare the computational effort of docking the PoseBusters benchmark set
across AutoDock Vina, DiffDock, and EquiBind (best configuration).

For every benchmark receptor-ligand pair this reads the per-complex timing each
tool already wrote to disk, joins it with the per-pose PoseBusters verdicts, and
reports two complementary views of cost, per method:

  1. WALL-CLOCK (time-to-result): how long you actually wait. Comparable across
     tools as elapsed time, but NOT as identical hardware work (AutoDock saturates
     many CPU threads; DiffDock/EquiBind use a GPU).
  2. RESOURCE-SECONDS (hardware work): the actual compute consumed, split by the
     hardware that did it —
       * CPU-core-seconds  = serial-equivalent CPU work (parallelism removed)
       * GPU-seconds        = time a GPU was occupied
     CPU-core-seconds and GPU-seconds are DIFFERENT currencies — a GPU-second is
     not a CPU-core-second — so they are reported side by side, never summed.

Both views are given overall and per docked pose, for *all generated* poses and
for only the poses that *survive the PoseBusters tests* (``pb_valid``).

A third figure (docking_effort_by_quality_tier.png) turns the "per-pose cost" into a
box-and-whiskers distribution across four NESTED quality tiers — all poses, PB-valid,
PB-valid & near-native (rmsd ≤ 2 Å), and PB-valid & near-native & correct form
(bestfit_rmsd ≤ 1 Å) — the same per-pose endpoints the headline accuracy report scores
on, overridable with --rmsd-column / --kabsch-column. For each complex the docking wall-clock is amortised over the
poses that reach a tier (wall_s ÷ poses-in-tier), so a box reads "wall-clock seconds of
docking effort spent per pose of that quality"; the cost climbs as the bar tightens and
the yield drops. Paired cross-tool tests per tier land in effort_by_quality_stats.json.

``--batch-comparison`` adds a CROSS-dataset figure (docking_effort_batch_vs_perrun.png)
that isolates the wall-clock batch docking saves. It builds BOTH campaigns and, per tool,
contrasts two ways to dock the same pairs holding marginal compute fixed: a fresh run per
(receptor, ligand) that reloads the tool each time (per-run) vs one invocation per receptor
that loads it once and streams its ligands (batch). The saving = L·(N_dockings − N_receptors)
with L the tool's fixed startup — free of the receptor-size / sample-count differences
between the campaigns because only the reload count varies. It is large only where a receptor
is reused (the Orai frames, ~300 ligands each) AND the tool has a real startup to amortise:
the benchmark gives every ligand its own crystal (reuse 1×) so its saving is 0; on Orai,
EquiBind's tiny marginal makes batching cut its cost ~57%, DiffDock saves hours but only ~5%
(the huge receptor dominates), and AutoDock saves nothing (Vina reloads the receptor per
ligand anyway). L is measured for EquiBind (pipeline phase timing) and --load-* tunable.

Timing sources (per complex):
  * AutoDock : Dockings/Benchmark/<id>/<prep>/docking_summary.json
               -> overall.total_time_seconds              (CPU wall; Vina, N threads)
               + optimization_log.csv for --autodock-refine
                 smina/gnina/gnina_refinement
               CPU-core-seconds = wall × N threads (--autodock-cpu, default 32);
               optimizer work is added on its recorded CPU/GPU device. Under
               --autodock-gnina-accounting device-occupancy a GPU-run gnina stage is
               charged the measured union of its provenance-sidecar intervals (wall_s
               and gpu_s) plus --autodock-gnina-host-cores per invocation-second of CPU.
  * Uni-Dock : Dockings/unidock_results_full_protein_vina_scoring/<id>/
               docking_summary.json -> elapsed_s, accepted only with the matching
               v2 .unidock_done + run_manifest.json commit. Tiled and Uni-Dock2
               remain separate methods. Tiled elapsed is reported as GPU occupancy;
               Uni-Dock2 exposes end-to-end wall time but no CPU/GPU phase split,
               so its resource-second fields are left unknown (NaN), never zero.
  * DiffDock : Dockings/Benchmark_DiffDock/<id>/docking_log.csv (docking wall, GPU)
               + optimization_log.csv (post-hoc rescoring, --diffdock-refine)
               GPU-seconds = docking wall + gnina rescoring (serial, CUDA CNN, GPU
               because the benchmark ran gnina_use_gpu true); with --diffdock-refine
               smina — or --no-diffdock-gnina-gpu — the rescoring is CPU instead.
  * EquiBind : Dockings/Benchmark_Equibind/<id>/pipeline_summary.json
               -> global_timing.pipeline_wall_time_s       (wall)
               -> all_results[].{prep,dock,post,refine}_time_s  (per-pose work)
               GPU-seconds = Σ dock_time_s (inference) + Σ refine_time_s only if gnina
               ran on the GPU (--eq-gnina-gpu). The benchmark ran EquiBind's gnina with
               --no_gpu (gnina_use_gpu false), so gnina refine is CPU-core-seconds;
               CPU-core-seconds = Σ (prep+post)_time_s (+ refine unless it was GPU gnina).

EquiBind "best configuration" (default **unguided + gnina**): the benchmark run
produced the full clamp×refine×pocket fan-out; of those, unguided + gnina gives the
best PoseBusters validity (~52% vs ~40% for smina). We therefore (a) count only that
variant's poses, (b) sum *only those poses'* per-pose resource-seconds (GPU/CPU), and
(c) estimate its wall-clock as those poses' OWN compute spread across the pipeline's
worker pool (config.n_parallel_workers) — a pure best-config-only estimate that does
NOT depend on the rest of the fan-out. The real full-run (all-variant) wall is kept in
the CSV (total_wall_full_h) for transparency. The refiner is selectable with
--eq-refine (smina is CPU; gnina is CPU too unless --eq-gnina-gpu, matching the
benchmark's gnina_use_gpu false); DiffDock's rescoring is selected in parallel with
--diffdock-refine / --diffdock-gnina-gpu (its gnina ran on the GPU).

All methods are compared on the COMMON set of complexes that have usable per-complex
timing for *every* method: complexes whose timing was lost to a skip/empty re-run
(elapsed_time_s / total_time_seconds == 0) are dropped per method, then intersected,
so total effort, pose counts and validity are like-for-like rather than spanning
different subsets.

Two docking campaigns are supported via ``--dataset`` (defaults + readers switch with it):
  * ``benchmark`` (default) — PoseBusters ligands in their NATIVE crystals. One timing file
    per complex; per_pose_metrics.csv carries crystal RMSD, so all four quality tiers apply.
  * ``orai_benchmark`` — the same ligands docked into Orai receptor frames. Each tool writes
    ONE aggregated timing log (AutoDock batch docking_log_*.csv; DiffDock docking_log.csv +
    optimization_log.csv; EquiBind a single pipeline_summary.json with per-pose all_results).
    pb_valid is derived from the PoseBusters test columns (reusing
    posebusters_validity_report.load_and_score); there is NO crystal, so only the ``all`` and
    ``PB-valid`` quality tiers exist. Defaults to the smina refiner for both DiffDock and
    EquiBind (matching the run configs). EquiBind's globally-decoupled pipeline records fresh
    timing only for the frames re-docked in that run (others cached at ~0 s), so the common
    timed set is restricted to those frames — surfaced in the printed drop diagnostics.

Run in the analysis env (conda env ``vina`` — pandas + matplotlib):

    # native-crystal benchmark
    python Scripts/Analysis/docking_effort_comparison.py \
        --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"

    # Orai × benchmark ligands (aggregated logs, crystal-free)
    python Scripts/Analysis/docking_effort_comparison.py --dataset orai_benchmark
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)
_PROJECT_ROOT = _HERE.parents[1]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
try:                                                # shared, unit-tested stats helpers
    import stats_utils as su                        # noqa: E402
except Exception:                                   # pragma: no cover - stats degrade off
    su = None
try:                                                # reuse the shared panel-labeller
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                   # pragma: no cover - fallback
    def _label_panels(axes, fontsize=12):
        flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
        for i, ax in enumerate(flat):
            if ax is not None and getattr(ax, "get_visible", lambda: True)():
                ax.text(-0.08, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                        fontsize=fontsize, fontweight="bold", va="bottom", ha="right")

# Method label -> (pretty name, bar colour, device used). The label + device are
# rewritten per run in main() to reflect the selected refine variant. Bar colours
# match the per-tool palette of 09f_pbvalid_yield_boxplot.png (posebusters_pose_comparison
# TOOL_COLORS): AutoDock #1f77b4, DiffDock #ff7f0e, EquiBind #2ca02c.
METHODS = {
    "autodock":  ("AutoDock Vina", "#1f77b4", "CPU"),
    "unidock":   ("Uni-Dock (tiled)", "#9467bd", "GPU"),
    # The current Uni-Dock2 summary records end-to-end wall time but not its
    # receptor-prep CPU / docking-GPU phase split. Keep resource work unknown
    # rather than inventing a CPU/GPU allocation.
    "unidock2":  ("Uni-Dock2", "#8c564b", "CPU+GPU (unpartitioned)"),
    "diffdock":  ("DiffDock (gnina)", "#ff7f0e", "GPU"),
    "equibind":  ("EquiBind (unguided+gnina)", "#2ca02c", "GPU+CPU"),
}

# AutoDock-specific post-processing identities. ``gnina`` retains its legacy
# empirical-minimization + CNN-rescore meaning; ``gnina_refinement`` is timed and
# selected independently. DiffDock/EquiBind continue to use only smina/gnina.
_AUTODOCK_REFINERS = frozenset({"smina", "gnina", "gnina_refinement"})
_AUTODOCK_GNINA_REFINERS = frozenset({"gnina", "gnina_refinement"})
# How the GPU-side AutoDock gnina stage is billed (see _autodock_effort). 'process-sum' is
# the historical accounting: gpu_s is the UNDIVIDED sum of the per-invocation elapsed times
# in optimization_log.csv. On the exh128 arm those invocations ran optimize_workers=16 deep
# on ONE device, so the sum is ~14x the time the device was actually busy (10.27 h against a
# measured 0.71 h). 'device-occupancy' charges the measured per-complex union of the
# optimizer intervals read from the provenance sidecars, and bills the stage's host CPU
# at DEFAULT_GNINA_HOST_CORES instead of dropping it.
GNINA_ACCOUNTINGS = ("process-sum", "device-occupancy")
DEFAULT_GNINA_ACCOUNTING = "process-sum"      # byte-preserves every committed effort tree
DEFAULT_GNINA_HOST_CORES = 2.1                 # recorded as measured in the exh128 gnina run config (optimize_cpu comment)
_GNINA_ACCOUNTING_MODE = DEFAULT_GNINA_ACCOUNTING   # set from args in main; read by labels
_ORAI_PROCESS_SUM_WARNED = False
_ORAI_PROCESS_SUM_HIT = False       # set when an Orai row was charged process-sum under a device-occupancy request


def _set_gnina_accounting_mode(args) -> None:
    global _GNINA_ACCOUNTING_MODE
    _GNINA_ACCOUNTING_MODE = getattr(args, "autodock_gnina_accounting", DEFAULT_GNINA_ACCOUNTING)


def _gnina_accounting_label() -> str:
    """Provenance label for the stats sidecars: the mode APPLIED, not merely requested."""
    if _GNINA_ACCOUNTING_MODE == "device-occupancy" and _ORAI_PROCESS_SUM_HIT:
        return "device-occupancy (benchmark path); Orai AutoDock gnina term: process-sum"
    return _GNINA_ACCOUNTING_MODE


def _warn_orai_process_sum_once() -> None:
    global _ORAI_PROCESS_SUM_WARNED, _ORAI_PROCESS_SUM_HIT
    _ORAI_PROCESS_SUM_HIT = True
    if not _ORAI_PROCESS_SUM_WARNED:
        print("NOTE: --autodock-gnina-accounting device-occupancy applies to the benchmark path only. "
              "Orai timing comes from one aggregated optimiser log with no per-invocation sidecars, "
              "so its AutoDock gnina term stays process-sum.", file=sys.stderr)
        _ORAI_PROCESS_SUM_WARNED = True
# EquiBind per-pose timing fields, classified by the hardware that does the work.
# Inference is GPU; conformer prep / UFF / IO are CPU. The REFINE step's hardware
# depends on the tool AND the pipeline's gnina_use_gpu flag: smina is always CPU;
# gnina is a CUDA-capable CNN but is run with --no_gpu (CPU) unless gnina_use_gpu is
# set. The benchmark ran EquiBind's gnina on CPU (gnina_use_gpu: false) and DiffDock's
# gnina on GPU (gnina_use_gpu: true) — hence the per-pipeline --*-gnina-gpu flags.
EQ_INFER_FIELD = "dock_time_s"                     # EquiBind model inference (GPU)
EQ_HOST_FIELDS = ("prep_time_s", "post_time_s")    # RDKit conformers / UFF / IO (CPU)

# resource-bar colours (by hardware, shared across methods)
CPU_COLOR, GPU_COLOR = "#8C8C8C", "#E1A100"

# ── Shared house style (matches orai_tool_agreement_compare.png) ──────────────
# Colour ALWAYS encodes the tool (AutoDock blue / DiffDock orange / EquiBind green,
# the METHODS palette). A secondary 2-category split within a tool (generated vs valid,
# CPU vs GPU work, per-run vs batch) is drawn as HATCH over the tool colour, never a
# second hue; ordinal quality tiers use a lightness ramp of the tool hue. Legends
# describe that split with NEUTRAL-GREY swatches and sit centred between the (sub)title
# and the axes — the tool identity is read from the x-axis + the consistent colour.
HATCH_SECONDARY = "////"          # hatch marking the "second" category of a within-tool split
LEGEND_GREY = "#9E9E9E"           # neutral swatch colour for the hatch/shade legends


def _tint(color, frac):
    """Blend a colour toward white by ``frac`` in [0,1] (0 = original, 1 = white)."""
    import matplotlib.colors as mcolors
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * frac, g + (1 - g) * frac, b + (1 - b) * frac)


def _legend_between(fig, handles, labels, y=0.9, ncol=None, fontsize=9, title=None):
    """Horizontal legend centred between the suptitle and the axes (reference layout)."""
    return fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y),
                      ncol=ncol or len(labels), frameon=True, fontsize=fontsize,
                      title=title, columnspacing=1.4)

# ── Batch-vs-per-run model (--batch-comparison) ───────────────────────────────
# Compares two ways to dock a set of receptor×ligand pairs, holding the marginal
# docking compute fixed and varying only how often the tool is (re)started:
#   * per-run : every pair is a fresh invocation that reloads the tool → pays the
#               fixed startup L each time  → total = M + L·N_dockings
#   * batch   : one invocation per RECEPTOR streams all its ligands, loading the tool
#               once → total = M + L·N_receptors
# The saving = L·(N_dockings − N_receptors) is the ONLY term that differs, so it is
# free of the receptor-size / sample-count confounds that dominate the raw per-pose
# cost across the two campaigns. It is large only when a receptor is reused (the Orai
# frames, ~300 ligands each) AND the tool has a real startup to amortise; the benchmark
# gives every ligand its own crystal (reuse 1×), so its saving is 0 by construction.
#
# Fixed per-run load L (seconds), per tool. All are --load-* overridable.
#   autodock : 0.0  — Vina has no neural model; even a "batch" loops Vina once per ligand,
#              so its receptor pdbqt load is NOT amortised. The one batchable cost (the
#              Meeko pdbqt prep) is upstream and absent from the docking timing → ~0 saving.
#   diffdock : 15.0 — ESM embeddings + score/confidence checkpoints, loaded once per
#              subprocess. ESTIMATE, bounded above by the smallest observed per-run docking
#              wall (~18.5 s at 30 samples), which is load + minimal inference.
#   equibind : MEASURED at runtime from the benchmark per-complex pipeline phase timing
#              (phase1 prep + [phase2 dock − pose inference]); ~4.7 s. Falls back to
#              DEFAULT_LOAD_EQUIBIND when the phase timing can't be read.
DEFAULT_LOAD_AUTODOCK = 0.0
DEFAULT_LOAD_DIFFDOCK = 15.0
DEFAULT_LOAD_EQUIBIND = 4.7
# per-run vs batch strategy colours (slow-you-wait red vs amortised green), shared panels.
PERRUN_COLOR, BATCH_COLOR = "#C44E52", "#55A868"
# Generic tool labels for the cross-dataset batch figure (refiner-agnostic — the two
# campaigns use different DiffDock refiners, and batching structure is refiner-independent).
BASE_NAMES = {
    "autodock": "AutoDock Vina", "unidock": "Uni-Dock (tiled)",
    "unidock2": "Uni-Dock2", "diffdock": "DiffDock", "equibind": "EquiBind",
}

# ── Quality tiers for the "effort to obtain a pose of quality X" view ─────────
# Nested per-pose subsets, each computed from the PoseBusters metrics already in
# per_pose_metrics.csv. The effort-by-quality figure amortises a complex's docking
# wall-clock over the poses that reach a tier (wall_s / poses-in-tier), so a box reads
# "seconds of docking effort spent per pose of that quality" — it climbs as the bar
# tightens and the yield falls.
#   all    every generated pose
#   valid  survives all PoseBusters tests   (pb_valid)
#   near2  valid AND placed near-native     (--rmsd-column   <= 2 Å)
#   form1  near2 AND correct internal form  (--kabsch-column <= 1 Å, the report's
#          --form-ok-kabsch default). Best-fit RMSD is a lower bound on placement RMSD,
#          so this only adds a conformer-fidelity constraint on top of near-nativeness.
NEAR2_RMSD_A = 2.0
FORM1_KABSCH_A = 1.0
# Which per-pose columns those thresholds are read from. These MUST match the headline
# accuracy report (posebusters_validity_report.py: --rmsd-column default "rmsd",
# --form-ok-kabsch scored on "bestfit_rmsd") or the cost view answers a different
# question than the accuracy view and a cost sub-cohort can report MORE qualifying poses
# than the full accuracy cohort. PoseBusters' own pb_rmsd / pb_kabsch_rmsd stay
# selectable via --rmsd-column / --kabsch-column but are strictly more permissive.
DEFAULT_RMSD_COLUMN = "rmsd"
DEFAULT_KABSCH_COLUMN = "bestfit_rmsd"
# (tier key, per-complex count column, x-tick / legend label, colour). Looser → stricter.
# The full four-tier ladder needs a crystal reference (pb_rmsd / pb_kabsch_rmsd). The
# Orai datasets have none, so main() swaps in ORAI_QUALITY_TIERS (all + PB-valid only).
QUALITY_TIERS = [
    ("all",   "generated", "all poses",                              "#C6DBEF"),
    ("valid", "pb_valid",  "PB-valid",                               "#6BAED6"),
    ("near2", "near2",     "PB-valid & RMSD ≤ 2 Å",                  "#2171B5"),
    ("form1", "form1",     "PB-valid & RMSD ≤ 2 Å & Kabsch ≤ 1 Å",   "#08306B"),
]
# Crystal-free datasets (Orai): no RMSD → only the placement-agnostic tiers survive.
ORAI_QUALITY_TIERS = [
    ("all",   "generated", "all poses", "#C6DBEF"),
    ("valid", "pb_valid",  "PB-valid",  "#2171B5"),
]

# Short dataset name spliced into figure/console titles; rewritten in main() per --dataset.
DATASET_LABEL = "PoseBuster benchmark"

# Per-dataset defaults for the paths + refiner picks that differ between the native-crystal
# benchmark and the Orai runs. CLI args for these default to None and are filled from here in
# main() (so an explicit flag still wins). ``title`` becomes DATASET_LABEL in figures/console.
#   benchmark      : one timing file per complex, per_pose_metrics w/ crystal RMSD;
#                    DiffDock=smina (oracle collapse pick), EquiBind=gnina (best PB validity).
#   orai_benchmark : aggregated per-tool logs, PoseBusters CSV (no crystal), smina refiner.
DATASET_DEFAULTS = {
    "benchmark": {
        "title": "PoseBuster benchmark",
        "per_pose_csv": "posebusters_results/benchmark/dock/"
                        "pose_comparison_report/per_pose_metrics.csv",
        "autodock_dir": "Dockings/Benchmark",
        "unidock_dir": "Dockings/unidock_results_full_protein_vina_scoring",
        "unidock2_dir": "Dockings/unidock2_results_full_protein_vina_scoring",
        "diffdock_dir": "Dockings/Benchmark_DiffDock",
        "equibind_dir": "Dockings/Benchmark_Equibind",
        "out_dir": "posebusters_results/benchmark/docking_effort",
        "ids_file": "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt",
        "autodock_refine": "raw",
        "diffdock_refine": "smina",             # Benchmark DiffDock canonical variant: smina (matches oracle collapse pick)
        "eq_refine": "gnina",
    },
    "orai_benchmark": {
        "title": "Orai × PoseBuster benchmark ligands",
        "per_pose_csv": "posebusters_results/orai_benchmark/dock/"
                        "posebusters_filtered_results.csv",
        "autodock_dir": "Dockings/Orai_Benchmark",
        "unidock_dir": None,
        "unidock2_dir": None,
        "diffdock_dir": "Dockings/Orai_Benchmark_DiffDock",
        "equibind_dir": "Dockings/Orai_Benchmark_Equibind",
        "out_dir": "posebusters_results/orai_benchmark/docking_effort",
        "ids_file": None,                       # no id filter — all frame×ligand combos
        "autodock_refine": "raw",
        "diffdock_refine": "smina",             # Orai DiffDock ran optimization: smina
        "eq_refine": "smina",                   # Orai EquiBind ran refine_tool: smina
    },
}


# Every arg whose default lives in DATASET_DEFAULTS rather than in the parser. Kept in one
# place because _apply_dataset_defaults fills it and _clone_dataset_args clears it: two
# hand-maintained copies of this tuple would drift, and a key present in the clear list but
# missing from the fill list leaves a clone with a None path.
PER_DATASET_ARGS = ("per_pose_csv", "autodock_dir", "unidock_dir", "unidock2_dir",
                    "diffdock_dir", "equibind_dir",
                    "out_dir", "ids_file", "autodock_refine", "diffdock_refine", "eq_refine")

# Where _record_explicit_dataset_args stashes its snapshot on the args namespace.
_EXPLICIT_KEYS_ATTR = "_explicit_dataset_args"
_EXPLICIT_DS_ATTR = "_explicit_dataset"


def _apply_dataset_defaults(args) -> None:
    """Fill any dataset-specific arg the user left unset (None) from DATASET_DEFAULTS."""
    d = DATASET_DEFAULTS[args.dataset]
    for k in PER_DATASET_ARGS:
        if getattr(args, k, None) is None:
            setattr(args, k, d[k])


def _record_explicit_dataset_args(args) -> None:
    """Snapshot which per-dataset flags the user actually typed, and for which campaign.

    MUST run before _apply_dataset_defaults. These args all parse to None so the defaults
    can be filled per --dataset; once filled, an explicitly passed path is indistinguishable
    from a dataset default, and _clone_dataset_args can no longer tell which of its resets
    would be discarding a user instruction. An empty string counts as explicit — the
    ``--unidock-dir ""`` idiom that suppresses an engine is falsy, not unset."""
    setattr(args, _EXPLICIT_DS_ATTR, args.dataset)
    setattr(args, _EXPLICIT_KEYS_ATTR,
            frozenset(k for k in PER_DATASET_ARGS + ("autodock_method",)
                      if getattr(args, k, None) is not None))

# ── Orai complex-id normalisation ─────────────────────────────────────────────
# The Orai per-pose CSV keys a complex by (frame, ligand); the aggregated timing logs
# key it by tool-specific names — ligand carries a "_ligand_start_conf" / "_start_conf"
# suffix and the frame an EquiBind "_cleaned" suffix. Normalise both to the canonical
# "<PDBID_LIG>__<frame>" so pose counts and timing join on the same key.
def _norm_lig(s: str) -> str:
    s = re.sub(r"_(?:ligand_)?start_conf$", "", str(s))
    return re.sub(r"_ligand$", "", s)


def _norm_frame(s: str) -> str:
    return re.sub(r"_cleaned$", "", str(s))


def _orai_cid(lig: str, frame: str) -> str:
    return f"{_norm_lig(lig)}__{_norm_frame(frame)}"


# ════════════════════════════════════════════════════════════════════════
# Per-complex timing readers (one per tool)
# ════════════════════════════════════════════════════════════════════════

def _autodock_time(cid: str, root: Path, prep: str) -> Optional[float]:
    p = root / cid / prep / "docking_summary.json"
    if not p.exists():
        return None
    try:
        t = float(json.loads(p.read_text())["overall"]["total_time_seconds"])
    except Exception:
        return None
    # A skip/empty re-run overwrites the original timed summary with
    # total_combinations=0 / total_time_seconds=0; treat that as "no timing"
    # (mirrors the DiffDock reader) rather than a real 0-second dock.
    return t if t > 0 else None


def _autodock_opt_time(cid: str, root: Path, prep: str, tool: str,
                       success_only: bool = False) -> Optional[float]:
    """Serial optimizer wall time for one AutoDock complex and committed tool.

    Failed attempts are included because they consumed real compute; collectors
    separately exclude their absent poses. A zero-only log usually means a cache
    reuse overwrote the original timing, so it is treated as missing rather than
    advertised as a free optimization.
    """
    path = root / cid / prep / "optimization_log.csv"
    if not path.is_file():
        return None
    try:
        # dtype=str: every column used here is coerced explicitly below, and type
        # inference is unsafe on these logs — a sha256 that happens to read as
        # "<digits>e<huge exponent>" sends pandas' C float parser into a
        # multi-hour scaling loop (observed on 7F51_BA7).
        rows = pd.read_csv(path, low_memory=False, dtype=str)
    except Exception:
        return None
    if not {"tool", "elapsed_time_s", "status"}.issubset(rows.columns):
        return None
    selected = rows[rows["tool"].astype(str).str.lower() == str(tool).lower()]
    statuses = selected["status"].fillna("").astype(str).str.strip().str.lower()
    # Failed rows consumed real compute, but cannot establish the existence of
    # the requested optimized variant. Once at least one committed-success row
    # exists, retain all its failed-attempt time in the measured total below.
    if not statuses.eq("success").any():
        return None
    if success_only:                      # rows that wrote a provenance sidecar
        selected = selected[statuses.eq("success").to_numpy()]
    elapsed = pd.to_numeric(selected["elapsed_time_s"], errors="coerce")
    elapsed = elapsed[np.isfinite(elapsed) & (elapsed >= 0)]
    total = float(elapsed.sum()) if len(elapsed) else 0.0
    return total if total > 0 else None


# Two top-level scalars out of a sidecar that also carries a full input-identity block.
# A byte scan is ~2x faster than json.loads over the exh128 arm's 9,126 files; a file that
# does not yield BOTH groups falls back to json.loads before being called unreadable, so a
# writer-side formatting change degrades to slow, never to wrong.
_PROV_SCAN = re.compile(rb'"created_at"\s*:\s*"([^"]+)"'
                        rb'|"optimizer_elapsed_time_s"\s*:\s*(-?[0-9][0-9.eE+-]*)')
_PROV_EPOCH = datetime(2000, 1, 1)   # naive local stamps differenced against a naive origin
_GNINA_OCCUPANCY_CACHE: Dict[tuple, Dict[str, tuple]] = {}


def _interval_union_s(intervals) -> float:
    """Total length of the union of half-open [start, end) intervals, in seconds."""
    if not intervals:
        return 0.0
    items = sorted(intervals)
    total, cur_s, cur_e = 0.0, items[0][0], items[0][1]
    for s, e in items[1:]:
        if s > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
        elif e > cur_e:
            cur_e = e
    return total + (cur_e - cur_s)


def _autodock_gnina_occupancy(root: Path, prep: str, tool: str) -> Dict[str, Tuple[float, float, int]]:
    """MEASURED per-complex device occupancy of one AutoDock optimizer tool.

    Returns ``{cid: (occupancy_s, process_s, n)}``: the union of that complex's own
    optimizer intervals, their undivided sum, and the invocation count.

    Read from the provenance sidecars, NOT from optimization_log.csv. The CSV's
    ``timestamp`` is a batch-WRITE stamp applied to every row of a complex after the
    whole batch returns (measured spread across a complex's 30 rows: 0.57 ms), so a
    union taken on it collapses to max(elapsed) and reports 0.39 h against the sidecars'
    0.71 h on the exh128 arm. The sidecar's ``created_at`` is stamped per invocation
    BEFORE execution (run_autodock.py, _build_optimizer_provenance) and
    ``optimizer_elapsed_time_s`` wraps only the subprocess.

    Complexes did not overlap each other on that arm: the sum of the 303 per-complex
    unions equals the union of all their intervals to 1e-6, so this is an exact additive
    decomposition of global occupancy and stays a true occupancy figure after the
    common-set restriction.
    """
    key = (str(Path(root).resolve()), str(prep), str(tool))
    hit = _GNINA_OCCUPANCY_CACHE.get(key)
    if hit is not None:
        return hit
    per: Dict[str, list] = {}
    for f in Path(root).glob(f"*/{prep}/docking/optimized_{tool}/*.provenance.json"):
        cid = f.relative_to(root).parts[0]
        created = elapsed = None
        try:
            raw = f.read_bytes()
            for m in _PROV_SCAN.finditer(raw):
                if m.group(1) is not None and created is None:
                    created = m.group(1).decode()
                elif m.group(2) is not None and elapsed is None:
                    elapsed = float(m.group(2))
                if created is not None and elapsed is not None:
                    break
            if created is None or elapsed is None:
                d = json.loads(raw)
                created, elapsed = d["created_at"], float(d["optimizer_elapsed_time_s"])
            t0 = (datetime.fromisoformat(created) - _PROV_EPOCH).total_seconds()
        except Exception:
            continue                       # unreadable sidecar: leave that invocation out
        if elapsed < 0:
            continue
        per.setdefault(cid, []).append((t0, t0 + elapsed))
    out = {cid: (_interval_union_s(iv), sum(e - s for s, e in iv), len(iv))
           for cid, iv in per.items()}
    _GNINA_OCCUPANCY_CACHE[key] = out
    return out


def _autodock_effort(
    cid: str,
    root: Path,
    prep: str,
    refine: str,
    docking_cpu: int,
    optimizer_cpu: int,
    gnina_gpu: bool = False,
    optimizer_workers: Optional[int] = None,
    gnina_accounting: str = DEFAULT_GNINA_ACCOUNTING,
    gnina_host_cores: float = DEFAULT_GNINA_HOST_CORES,
) -> Optional[tuple]:
    """Return ``(wall_s, gpu_s, cpu_core_s)`` for an AutoDock variant.

    ``gnina_accounting`` decides what gpu_s IS for a GPU-run gnina stage: 'process-sum'
    (default, historical) the undivided Σ of per-invocation elapsed times; 'device-occupancy'
    the measured per-complex union of those invocations' intervals, with the stage's host
    CPU billed at ``gnina_host_cores`` per invocation-second instead of dropped.
    """
    dock = _autodock_time(cid, root, prep)
    if dock is None:
        return None
    dock_cpu_s = dock * max(int(docking_cpu), 1)
    if refine not in _AUTODOCK_REFINERS:
        return dock, 0.0, dock_cpu_s
    opt = _autodock_opt_time(cid, root, prep, refine)
    if opt is None:
        return None
    # Σ per-pose optimizer elapsed IS the optimizer's wall clock only while the poses of a
    # complex were minimised one at a time. Arms run with optimize_workers > 1 overlap those
    # intervals, so the sum becomes compute time and the divisor recovers an approximate wall
    # clock. Under 'process-sum' the resource-seconds are left un-divided (sixteen concurrent
    # gnina calls bill sixteen calls' worth of GPU); under 'device-occupancy' gpu_s is the
    # MEASURED union of the intervals instead, so no divisor is involved at all.
    opt_wall = opt / max(int(optimizer_workers), 1) if optimizer_workers else opt
    if refine in _AUTODOCK_GNINA_REFINERS and gnina_gpu:
        if gnina_accounting == "device-occupancy":
            occ = _autodock_gnina_occupancy(root, prep, refine).get(cid)
            if occ is None:                # no sidecars: never fabricate an occupancy
                print(f"NOTE: {cid} {refine}: no provenance sidecars; complex dropped from the "
                      "timed set under --autodock-gnina-accounting device-occupancy", file=sys.stderr)
                return None
            occupancy_s, process_s, n = occ
            # Same invocations, two records. Only committed-success rows write a sidecar, so
            # the sidecar sum is compared with the success-only CSV sum (rounded to 2 dp); a
            # gap wider than rounding means the two sources describe different sets (stale
            # optimized_* dir, partial re-run). Failed attempts still bill CPU below.
            opt_ok = _autodock_opt_time(cid, root, prep, refine, success_only=True) or 0.0
            if abs(process_s - opt_ok) > 0.01 * max(n, 1) + 1e-6:
                raise SystemExit(
                    f"ERROR: {cid} {refine}: provenance sidecars sum to {process_s:.3f} s over "
                    f"{n} invocations but optimization_log.csv success rows sum to {opt_ok:.3f} s. "
                    "Under --autodock-gnina-accounting device-occupancy the two must describe "
                    "the same invocations.")
            # wall_s and gpu_s are the MEASURED union (the /optimizer_workers idealisation is
            # not needed once the stage's occupancy is in hand). The host cores each concurrent
            # gnina held are billed as CPU instead of dropped; opt (the full CSV sum) multiplies
            # the cores because it is the record that includes failed attempts' real CPU.
            return (dock + occupancy_s, occupancy_s,
                    dock_cpu_s + opt * float(gnina_host_cores))
        if gnina_accounting != "process-sum":
            raise ValueError(f"unknown gnina accounting {gnina_accounting!r}; "
                             f"expected one of {GNINA_ACCOUNTINGS}")
        return dock + opt_wall, opt, dock_cpu_s
    return dock + opt_wall, 0.0, dock_cpu_s + opt * max(int(optimizer_cpu), 1)


def _json_dict(path: Path) -> Optional[dict]:
    try:
        value = json.loads(Path(path).read_text())
    except Exception:
        return None
    return value if isinstance(value, dict) else None


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _current_hashed_inputs(provenance: dict) -> bool:
    """Require every fingerprinted docking input to still be the same file."""
    try:
        inputs = provenance["inputs"]
        if not isinstance(inputs, dict):
            return False
        for name in ("receptor", "ligand", "box"):
            item = inputs[name]
            raw = item["path"]
            path = Path(raw).expanduser()
            if (not isinstance(raw, str) or not path.is_absolute()
                    or str(path.resolve()) != raw or not path.is_file()
                    or isinstance(item.get("size"), bool)
                    or int(item["size"]) != path.stat().st_size
                    or item["sha256"] != _sha256_file(path)):
                return False
    except (KeyError, TypeError, ValueError, OSError):
        return False
    return True


_UNIDOCK2_SENTINEL_SCHEMA = 3
_UNIDOCK2_PROVENANCE_KEYS = (
    "schema_version", "engine", "tool", "inputs", "center", "size",
    "effective_config", "effective_config_sha256", "driver_execution",
)


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _current_unidock2_tool_identity(tool: object) -> bool:
    """Rehash the executable artifacts and Python source tree recorded by v3."""
    expected_keys = {
        "path", "version", "launcher_sha256", "tleap_path", "ambpdb_path",
        "artifacts", "python_source_tree",
    }
    if not isinstance(tool, dict) or set(tool) != expected_keys:
        return False
    if not isinstance(tool.get("version"), str) or not tool["version"].strip():
        return False
    try:
        artifacts = tool["artifacts"]
        if not isinstance(artifacts, list) or not artifacts:
            return False
        artifact_paths: set[str] = set()
        for item in artifacts:
            if not isinstance(item, dict) or set(item) != {"path", "size", "sha256"}:
                return False
            raw = item["path"]
            path = Path(raw).expanduser()
            if (not isinstance(raw, str) or not path.is_absolute()
                    or str(path.resolve()) != raw or not path.is_file()
                    or isinstance(item["size"], bool) or not isinstance(item["size"], int)
                    or item["size"] != path.stat().st_size
                    or not _is_sha256(item["sha256"])
                    or item["sha256"] != _sha256_file(path)
                    or raw in artifact_paths):
                return False
            artifact_paths.add(raw)
        for key in ("path", "tleap_path", "ambpdb_path"):
            raw = tool[key]
            path = Path(raw).expanduser()
            if (not isinstance(raw, str) or not path.is_absolute()
                    or str(path.resolve()) != raw or raw not in artifact_paths):
                return False
        if (not _is_sha256(tool["launcher_sha256"])
                or tool["launcher_sha256"] != _sha256_file(tool["path"])):
            return False

        source_tree = tool["python_source_tree"]
        if (not isinstance(source_tree, dict)
                or set(source_tree) != {"path", "file_count", "size", "sha256"}):
            return False
        source_root = Path(source_tree["path"]).expanduser()
        if (not source_root.is_absolute()
                or str(source_root.resolve()) != source_tree["path"]
                or not source_root.is_dir()):
            return False
        from Scripts.Docking import run_unidock2 as unidock2
        return unidock2._json_equal(
            source_tree, unidock2._python_tree_provenance(source_root))
    except (ImportError, KeyError, OSError, RuntimeError, TypeError, ValueError):
        return False


def _validated_tiled_timing_commit(cdir: Path, summary: dict) -> bool:
    """Validate the normalized three-record tiled Uni-Dock commit."""
    marker = _json_dict(cdir / ".unidock_done")
    manifest = _json_dict(cdir / "run_manifest.json")
    if not marker or not manifest:
        return False
    if (marker.get("schema_version") != 2
            or marker.get("engine") != "unidock-tiled"
            or marker.get("status") != "complete"
            or manifest.get("schema_version") != 2
            or manifest.get("status") != "complete"
            or str(summary.get("status") or "").strip().lower()
            not in {"success", "done"}):
        return False
    for record in (summary, marker, manifest):
        if record.get("commit_status") != "complete":
            return False
    keys = (
        "fingerprint", "generation_id", "output_file", "output_sha256",
        "num_poses", "n_subboxes",
    )
    if any(marker.get(key) in (None, "") for key in keys):
        return False
    if any(not (marker.get(key) == manifest.get(key) == summary.get(key))
           for key in keys):
        return False
    try:
        output_name = marker["output_file"]
        if (not isinstance(output_name, str) or Path(output_name).name != output_name
                or output_name != f"{cdir.name}_unidock_out.pdbqt"):
            return False
        output = cdir / output_name
        if marker["output_sha256"] != _sha256_file(output):
            return False
        from Scripts.Docking import run_unidock_tiled as tiled
        poses = tiled.parse_out_pdbqt(output)
        affinities = [pose.affinity for pose in poses]
        if (len(poses) != int(marker["num_poses"]) or len(poses) < 1
                or int(marker["n_subboxes"]) < 1
                or affinities != sorted(affinities)):
            return False
        provenance = manifest["provenance"]
        canonical = json.dumps(
            provenance, sort_keys=True, separators=(",", ":"), allow_nan=False)
        if hashlib.sha256(canonical.encode()).hexdigest() != marker["fingerprint"]:
            return False
        if (provenance.get("schema_version") != 2
                or not isinstance(provenance.get("engine"), dict)
                or len(provenance.get("subboxes") or []) != int(marker["n_subboxes"])
                or not _current_hashed_inputs(provenance)):
            return False
    except (ImportError, KeyError, OSError, TypeError, ValueError):
        return False
    return True


def _validated_unidock2_timing_commit(cdir: Path, cid: str, summary: dict) -> bool:
    """Validate Uni-Dock2's authoritative schema-3 sentinel and timed summary."""
    marker_path = cdir / f"{cid}_unidock2_completion.json"
    marker = _json_dict(marker_path)
    if not marker:
        return False
    output = cdir / f"{cid}_unidock2_out.sdf"
    prepared = cdir / f"{cid}_unidock2_receptor_prepared.pdb"
    if (marker.get("schema_version") != _UNIDOCK2_SENTINEL_SCHEMA
            or marker.get("engine") != "unidock2"
            or marker.get("status") != "success"
            or marker.get("output_file") != output.name
            or marker.get("prepared_receptor_file") != prepared.name):
        return False
    try:
        provenance = {key: marker[key] for key in _UNIDOCK2_PROVENANCE_KEYS}
        effective = marker["effective_config"]
        effective_json = json.dumps(
            effective, sort_keys=True, separators=(",", ":"), allow_nan=False)
        if hashlib.sha256(effective_json.encode()).hexdigest() != marker[
                "effective_config_sha256"]:
            return False
        canonical = json.dumps(
            provenance, sort_keys=True, separators=(",", ":"), allow_nan=False)
        if hashlib.sha256(canonical.encode()).hexdigest() != marker["fingerprint"]:
            return False
        if not _current_unidock2_tool_identity(marker["tool"]):
            return False
        lock_path = marker["driver_execution"]["gpu_lock_file"]
        if (set(marker["driver_execution"]) != {"gpu_lock_file"}
                or not isinstance(lock_path, str) or not Path(lock_path).is_absolute()
                or str(Path(lock_path).expanduser().resolve()) != lock_path
                or not _current_hashed_inputs(provenance)):
            return False

        advanced = effective["Advanced"]
        settings = effective["Settings"]
        if settings.get("search_mode") != "free" or settings.get("task") != "screen":
            return False
        from Scripts.Docking import run_unidock2 as unidock2
        validation_cfg = {
            "num_pose": int(advanced["num_pose"]),
            "energy_range": float(advanced["energy_range"]),
        }
        if not unidock2._completion_matches(
                output, marker_path, marker["fingerprint"], provenance, validation_cfg):
            return False
        validation = unidock2._validate_docking_output(output, validation_cfg)
        num_poses = int(marker["num_poses"])
        best = float(marker["best_affinity"])
        if (num_poses != validation["num_poses"]
                or not np.isfinite(best)
                or abs(best - validation["best_affinity"]) > 0.002
                or marker["output_sha256"] != _sha256_file(output)):
            return False
        if (not _is_sha256(marker.get("receptor_prmtop_sha256"))
                or not _is_sha256(marker.get("receptor_inpcrd_sha256"))
                or int(marker["input_protein_heavy_atoms"]) < 1
                or not np.isfinite(float(marker["completed_unix_s"]))
                or float(marker["completed_unix_s"]) <= 0):
            return False

        # The sentinel is commit authority; these fields only bind its generation
        # to the separate summary from which elapsed_s is being attributed.
        if (str(summary.get("status") or "").strip().lower() != "success"
                or summary.get("engine") != "unidock2"
                or summary.get("fingerprint") != marker["fingerprint"]
                or summary.get("generation_id") != marker["generation_id"]
                or int(summary.get("num_poses")) != num_poses
                or abs(float(summary.get("best_affinity")) - best) > 0.002
                or Path(str(summary.get("output_sdf"))).expanduser().resolve()
                != output.resolve()
                or Path(str(summary.get("completion_manifest"))).expanduser().resolve()
                != marker_path.resolve()
                or Path(str(summary.get("prepared_receptor_pdb"))).expanduser().resolve()
                != prepared.resolve()
                or summary.get("published_output_current") is not True
                or summary.get("inputs_changed_during_run") is not False):
            return False
    except (ImportError, KeyError, OSError, TypeError, ValueError):
        return False
    return True


def _unidock_effort(cid: str, root: Path, engine: str) -> Optional[tuple]:
    """Committed end-to-end effort for tiled Uni-Dock or Uni-Dock2.

    The timing is accepted only when the production ``docking_summary.json``
    agrees with the engine's completion manifest. This prevents a failed/stale
    summary from being paired with the currently collected poses.

    Returns ``(wall_s, gpu_s, cpu_core_s)``. Tiled Uni-Dock runs its sequential
    searches on the GPU, so its elapsed time is reported as GPU occupancy (the
    summary does not split its small planning/merge CPU overhead). Uni-Dock2's
    end-to-end timing includes substantial receptor preparation plus GPU docking
    but exposes no phase timing; both resource currencies are therefore NaN while
    the measured wall time remains usable.
    """
    cdir = Path(root) / cid
    summary = _json_dict(cdir / "docking_summary.json")
    if not summary:
        return None
    try:
        wall = float(summary.get("elapsed_s"))
    except (TypeError, ValueError):
        return None
    if not np.isfinite(wall) or wall <= 0:
        return None

    if engine == "unidock":
        if not _validated_tiled_timing_commit(cdir, summary):
            return None
        return wall, wall, 0.0

    if engine == "unidock2":
        if not _validated_unidock2_timing_commit(cdir, cid, summary):
            return None
        return wall, float("nan"), float("nan")

    raise ValueError(f"unknown Uni-Dock engine: {engine}")


def _diffdock_time(cid: str, root: Path) -> Optional[float]:
    """Real per-complex DiffDock docking wall-clock from docking_log.csv.

    Use docking_log.csv's elapsed_time_s (recorded on the FIRST real dock of each
    complex, e.g. the 2026-05-30 production run). Multiple rows (appended skip-runs
    write 0.0) -> take the max, i.e. the one real dock. Returns None when no positive
    time is recorded (a skip-run whose original production timing was overwritten).

    docking_summary.json is NOT consulted: its overall.total_time_seconds is
    co-zeroed on every skip-run (skipped:1 -> 0.0), so it can never recover a
    complex the log has already lost.
    """
    log = root / cid / "docking_log.csv"
    if log.exists():
        try:
            import csv
            with open(log) as f:
                ts = [float(r.get("elapsed_time_s") or 0) for r in csv.DictReader(f)]
            mx = max(ts) if ts else 0.0
            if mx > 0:
                return mx
        except Exception:
            pass
    return None


def _diffdock_opt_time(cid: str, root: Path, tool: str) -> Optional[float]:
    """Total post-hoc {tool} rescoring/minimization time from optimization_log.csv.

    DiffDock's ``optimized_<tool>`` variants are produced by minimizing each DiffDock
    pose with smina (CPU) or gnina (CUDA CNN, GPU or CPU depending on gnina_use_gpu).
    The rescoring is serial per pose (verified: timestamp span ≈ Σ elapsed), so the
    summed per-pose elapsed_time_s is the added wall-clock. Returns None if the log
    has no rows for ``tool``.
    """
    p = root / cid / "optimization_log.csv"
    if not p.exists():
        return None
    try:
        import csv
        with open(p) as f:
            ts = [float(r.get("elapsed_time_s") or 0.0)
                  for r in csv.DictReader(f) if r.get("tool") == tool]
    except Exception:
        return None
    return sum(ts) if ts else None


def _diffdock_effort(cid: str, root: Path, refine: str,
                     gnina_gpu: bool = True) -> Optional[tuple]:
    """Return (wall_s, gpu_s, cpu_core_s) for the requested DiffDock variant.

    ``raw``   = DiffDock docking only (GPU).
    ``smina`` = docking (GPU) + smina minimization (CPU).
    ``gnina`` = docking (GPU) + gnina CNN minimization (GPU if gnina_gpu else CPU).
    None if the docking wall — or, for a rescored variant, its optimization time —
    is missing.
    """
    dock = _diffdock_time(cid, root)          # docking wall (GPU)
    if dock is None:
        return None
    if refine in ("smina", "gnina"):
        opt = _diffdock_opt_time(cid, root, refine)
        if opt is None:
            return None
        if refine == "gnina" and gnina_gpu:
            return dock + opt, dock + opt, 0.0        # docking + gnina both GPU
        return dock + opt, dock, opt                  # docking GPU; smina / --no_gpu gnina on CPU
    return dock, dock, 0.0                            # raw docking only


def _equibind_time(cid: str, root: Path, eq_mode: str, eq_refine: str,
                   eq_clamp: str, gnina_gpu: bool = False) -> Optional[tuple]:
    """Return (wall_full_s, wall_best_s, best_gpu_s, best_cpu_core_s) for one complex.

    best_gpu_s / best_cpu_core_s are the summed per-pose resource-seconds of the
    best-config (unguided + ``eq_refine`` + clampOFF) poses ONLY, split by hardware:
    inference is GPU, prep/post are CPU, refine follows the tool AND its device
    (smina -> CPU; gnina -> GPU only if gnina_gpu, else CPU via --no_gpu). wall_best
    estimates the wall to run *only* that best config — those poses' compute spread
    across config.n_parallel_workers — so it depends solely on the best-config poses,
    not the rest of the fan-out. wall_full is the real full-run (all-variant) wall,
    kept only for transparency.
    """
    p = root / cid / "pipeline_summary.json"
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
    except Exception:
        return None
    wall_full = float(d.get("global_timing", {}).get("pipeline_wall_time_s", np.nan))
    if not np.isfinite(wall_full):
        return None
    rows = d.get("all_results") or []
    n_workers = float(d.get("config", {}).get("n_parallel_workers") or 32) or 32.0

    refine_on_gpu = (eq_refine == "gnina") and gnina_gpu       # gnina only, and only if run on GPU

    def gpu_s(r):
        g = float(r.get(EQ_INFER_FIELD) or 0.0)                 # inference (GPU)
        if refine_on_gpu:
            g += float(r.get("refine_time_s") or 0.0)          # gnina CNN (GPU)
        return g

    def cpu_s(r):
        c = sum(float(r.get(k) or 0.0) for k in EQ_HOST_FIELDS)  # prep/post (CPU)
        if not refine_on_gpu:
            c += float(r.get("refine_time_s") or 0.0)           # smina (CPU)
        return c

    best = [r for r in rows if r.get("mode") == eq_mode
            and r.get("refine_variant") == eq_refine
            and r.get("clamp_variant") == eq_clamp]
    best_gpu = sum(gpu_s(r) for r in best)
    best_cpu = sum(cpu_s(r) for r in best)
    best_c = best_gpu + best_cpu
    # Wall to run ONLY the best config: its OWN compute-seconds spread across the
    # pipeline's worker pool (config.n_parallel_workers). Depends solely on the
    # unguided+<refine> poses — the rest of the fan-out (other modes / raw / smina)
    # no longer enters. Idealised (assumes perfect n-worker overlap; the full run's
    # measured overlap was a touch lower, so this is a slight lower bound).
    wall_best = best_c / n_workers if n_workers > 0 else np.nan
    return wall_full, wall_best, best_gpu, best_cpu


# ════════════════════════════════════════════════════════════════════════
# Orai × Benchmark timing readers (aggregated, one log per tool)
# ════════════════════════════════════════════════════════════════════════
# Unlike the benchmark (one timing file per complex), the Orai run writes ONE
# aggregated log per tool at the toolchain root, keyed by combo_name / ligand+frame.
# Each reader preloads its whole log into a {canonical-cid -> effort} dict once.

def _orai_autodock_timing(root: Path) -> Dict[str, float]:
    """{cid -> Vina wall-clock seconds} from the batch docking_log_*.csv.

    A combo's real timing is the max positive elapsed_time_s across its rows (skip/empty
    re-runs append 0.0). Combos whose only rows are 0 (original timing overwritten) are
    dropped, mirroring the benchmark AutoDock reader.
    """
    import csv
    out: Dict[str, float] = {}
    for p in sorted(root.glob("docking_log_*.csv")):
        try:
            with open(p) as f:
                for r in csv.DictReader(f):
                    t = float(r.get("elapsed_time_s") or 0.0)
                    if t <= 0:
                        continue
                    k = _orai_cid(r.get("ligand_name", ""), r.get("protein_name", ""))
                    out[k] = max(out.get(k, 0.0), t)
        except Exception:                                    # pragma: no cover
            continue
    return out


def _orai_autodock_opt_timing(root: Path, tool: str) -> Dict[str, float]:
    """Optimizer seconds by Orai ``cid`` from all AutoDock optimization logs."""
    out: Dict[str, float] = {}
    successful: set[str] = set()
    for path in sorted(root.glob("**/optimization_log.csv")):
        try:
            # dtype=str for the same reason as _autodock_opt_time: sha256 columns
            # can look like huge-exponent floats and stall the C parser.
            rows = pd.read_csv(path, low_memory=False, dtype=str)
        except Exception:
            continue
        if not {"tool", "elapsed_time_s", "status"}.issubset(rows.columns):
            continue
        rows = rows[rows["tool"].astype(str).str.lower() == str(tool).lower()]
        for record in rows.to_dict("records"):
            elapsed = pd.to_numeric(pd.Series([record.get("elapsed_time_s")]),
                                    errors="coerce").iloc[0]
            if not np.isfinite(elapsed) or elapsed < 0:
                continue
            ligand = record.get("ligand_name", "")
            protein = record.get("protein_name", "")
            if not ligand or not protein:
                combo = str(record.get("combo_name") or "")
                if "__" not in combo:
                    continue
                ligand, protein = combo.split("__", 1)
            cid = _orai_cid(str(ligand), str(protein))
            out[cid] = out.get(cid, 0.0) + float(elapsed)
            if str(record.get("status") or "").strip().lower() == "success":
                successful.add(cid)
    return {
        cid: seconds for cid, seconds in out.items()
        if cid in successful and seconds > 0
    }


def _orai_diffdock_timing(root: Path, refine: str,
                          gnina_gpu: bool = True) -> Dict[str, tuple]:
    """{cid -> (wall_s, gpu_s, cpu_core_s)} for the requested DiffDock variant.

    Reads the aggregated docking_log.csv (docking wall, GPU) and, for a rescored
    variant, optimization_log.csv (summed per-pose smina/gnina time). Same currency
    split as the benchmark _diffdock_effort: raw / GPU-gnina add no CPU work; smina and
    CPU (--no_gpu) gnina are CPU-core-seconds. A combo missing its rescoring time drops.
    """
    import csv
    dock: Dict[str, float] = {}
    try:
        with open(root / "docking_log.csv") as f:
            for r in csv.DictReader(f):
                t = float(r.get("elapsed_time_s") or 0.0)
                if t <= 0:
                    continue
                k = _orai_cid(r.get("ligand_name", ""), r.get("protein_name", ""))
                dock[k] = max(dock.get(k, 0.0), t)
    except Exception:                                        # pragma: no cover
        return {}
    if refine not in ("smina", "gnina"):
        return {k: (v, v, 0.0) for k, v in dock.items()}     # raw docking only (GPU)
    opt: Dict[str, float] = {}
    try:
        with open(root / "optimization_log.csv") as f:
            for r in csv.DictReader(f):
                if r.get("tool") != refine:
                    continue
                k = _orai_cid(r.get("ligand_name", ""), r.get("protein_name", ""))
                opt[k] = opt.get(k, 0.0) + float(r.get("elapsed_time_s") or 0.0)
    except Exception:                                        # pragma: no cover
        return {}
    out: Dict[str, tuple] = {}
    for k, dv in dock.items():
        if k not in opt:                                     # no rescoring time → drop
            continue
        ov = opt[k]
        if refine == "gnina" and gnina_gpu:
            out[k] = (dv + ov, dv + ov, 0.0)                 # docking + gnina both GPU
        else:
            out[k] = (dv + ov, dv, ov)                       # docking GPU; smina / CPU gnina on CPU
    return out


def _orai_equibind_timing(root: Path, eq_mode: str, eq_refine: str,
                          eq_clamp: Optional[str], gnina_gpu: bool = False) -> Dict[str, tuple]:
    """{cid -> (wall_s, wall_s, gpu_s, cpu_core_s)} from the single pipeline_summary.json.

    The Orai EquiBind run is a single variant (no clamp/refine fan-out), so it writes one
    aggregated summary whose ``all_results`` hold per-pose resource-seconds. Those are
    SERIAL work-seconds (Σ post/refine ÷ n_parallel_workers ≈ the measured phase wall), so
    the per-complex wall is the summed compute ÷ the worker pool — the same best-config
    estimate the benchmark reader uses. wall_full == wall (the one variant IS the whole run).
    Split by hardware: inference (dock_time_s) is GPU; prep/post are CPU; refine follows the
    tool AND device (smina → CPU; gnina → GPU only if gnina_gpu). Combos with zero recorded
    compute (poses cached from a prior run) are dropped, mirroring the benchmark's skip-drop.
    """
    p = root / "pipeline_summary.json"
    if not p.exists():
        return {}
    try:
        d = json.loads(p.read_text())
    except Exception:                                        # pragma: no cover
        return {}
    n_workers = float(d.get("config", {}).get("n_parallel_workers") or 32) or 32.0
    refine_on_gpu = (eq_refine == "gnina") and gnina_gpu
    agg: Dict[str, list] = {}
    for r in d.get("all_results") or []:
        if (r.get("mode") != eq_mode or r.get("refine_variant") != eq_refine
                or r.get("clamp_variant") != eq_clamp):
            continue
        gpu = float(r.get(EQ_INFER_FIELD) or 0.0)            # inference (GPU)
        cpu = sum(float(r.get(k) or 0.0) for k in EQ_HOST_FIELDS)  # prep/post (CPU)
        rt = float(r.get("refine_time_s") or 0.0)
        if refine_on_gpu:
            gpu += rt
        else:
            cpu += rt
        k = _orai_cid(r.get("ligand_name", ""), r.get("protein_name", ""))
        a = agg.setdefault(k, [0.0, 0.0])
        a[0] += gpu
        a[1] += cpu
    out: Dict[str, tuple] = {}
    for k, (g, c) in agg.items():
        tot = g + c
        if tot <= 0:                                         # cached / untimed combo → drop
            continue
        wall = tot / n_workers if n_workers > 0 else np.nan
        out[k] = (wall, wall, g, c)
    return out


# ════════════════════════════════════════════════════════════════════════
# Aggregation
# ════════════════════════════════════════════════════════════════════════

def _per_pose_methods(args) -> Dict[str, str]:
    """Map each method key to its ``method`` value in per_pose_metrics.csv, per the
    selected variant (EquiBind mode/refine, DiffDock refine)."""
    ad_refine = getattr(args, "autodock_refine", "raw") or "raw"
    ad = "autodock" if ad_refine == "raw" else f"autodock_{ad_refine}"
    # The derivation above only ever spells the four Meeko labels, so the ADFRsuite
    # ladder (autodock_mgltools_exh128_gnina and its rungs) is unreachable through
    # --autodock-refine alone: --autodock-prep swaps the directory segment, not the
    # per-pose method key. An explicit label overrides the derivation.
    ad = getattr(args, "autodock_method", None) or ad
    dd = "diffdock" if args.diffdock_refine == "raw" else f"diffdock_{args.diffdock_refine}"
    return {
        "autodock": ad,
        "unidock": "unidock",
        "unidock2": "unidock2",
        "diffdock": dd,
        "equibind": f"equibind_{args.eq_mode}_{args.eq_refine}",   # clamp-off drops the clamp token
    }


def _endpoint_cols(args) -> tuple:
    """(rmsd_column, kabsch_column) for the near2/form1 tiers — canonical by default."""
    return (getattr(args, "rmsd_column", None) or DEFAULT_RMSD_COLUMN,
            getattr(args, "kabsch_column", None) or DEFAULT_KABSCH_COLUMN)


def _per_pose_optimizers(args) -> Dict[str, str]:
    """Additional provenance filter for tools whose method label stays stable."""
    refine = getattr(args, "autodock_refine", "raw") or "raw"
    return {"autodock": "original" if refine == "raw" else refine}


def _filter_optimizer_variant(sub: pd.DataFrame, key: str,
                              optimizer_by_method: Optional[Dict[str, str]]) -> pd.DataFrame:
    wanted = (optimizer_by_method or {}).get(key)
    if not wanted:
        return sub
    aliases = {"raw": "original", "native": "original", "none": "original"}
    want = aliases.get(str(wanted).strip().lower(), str(wanted).strip().lower())
    if "optimizer" in sub.columns:
        actual = sub["optimizer"].fillna("original").astype(str).str.strip().str.lower()
    elif "pose_file" in sub.columns:
        paths = sub["pose_file"].fillna("").astype(str).str.lower().str.replace("\\", "/", regex=False)
        actual = pd.Series("original", index=sub.index, dtype=object)
        actual.loc[paths.str.contains("/optimized_smina/")] = "smina"
        actual.loc[paths.str.contains("/optimized_gnina_refinement/")] = "gnina_refinement"
        actual.loc[paths.str.contains("/optimized_gnina/")] = "gnina"
    else:
        # Without either explicit provenance or an optimizer path component the
        # requested variant cannot be isolated; fail closed instead of pooling.
        return sub.iloc[0:0]
    actual = actual.map(lambda value: aliases.get(value, value))
    return sub[actual == want]


def load_pose_counts(csv: Path, ids: Optional[set],
                     per_pose_method: Dict[str, str],
                     optimizer_by_method: Optional[Dict[str, str]] = None,
                     rmsd_column: str = DEFAULT_RMSD_COLUMN,
                     kabsch_column: str = DEFAULT_KABSCH_COLUMN) -> Dict[str, pd.DataFrame]:
    """Per method: a DataFrame indexed by complex id with generated + pb_valid counts.

    ``rmsd_column`` / ``kabsch_column`` pick the near-nativeness and internal-form
    endpoints. They default to the CANONICAL columns the headline accuracy tables score
    on (posebusters_validity_report.py --rmsd-column / --form-ok-kabsch), so the cost
    view and the accuracy view answer the same question. PoseBusters' own ``pb_rmsd`` /
    ``pb_kabsch_rmsd`` remain selectable but are systematically more permissive
    (pb_rmsd <= rmsd pose-by-pose), which inflates the qualifying-pose counts.
    """
    df = pd.read_csv(csv, low_memory=False)
    if "pb_valid" in df.columns:
        df["pb_valid"] = df["pb_valid"].astype(str).str.lower().isin(("true", "1", "1.0"))
    else:
        df["pb_valid"] = False
    missing = [c for c in (rmsd_column, kabsch_column) if c not in df.columns]
    if missing:
        raise SystemExit(
            f"ERROR: per-pose CSV {csv} has no column(s) {missing}. Available RMSD-like "
            f"columns: {[c for c in df.columns if 'rmsd' in c.lower()]}. Pick one with "
            "--rmsd-column / --kabsch-column (silently treating them as NaN would zero "
            "out every near-native pose).")
    # Per-pose quality-tier membership (nested subsets) for the effort-by-quality view.
    # NaN RMSDs (e.g. a target with no crystal) compare False, so such poses simply drop
    # out of near2/form1 rather than poisoning the count.
    def _numcol(name):
        return (pd.to_numeric(df[name], errors="coerce") if name in df.columns
                else pd.Series(np.nan, index=df.index))
    df["near2"] = df["pb_valid"] & (_numcol(rmsd_column) <= NEAR2_RMSD_A)
    df["form1"] = df["near2"] & (_numcol(kabsch_column) <= FORM1_KABSCH_A)
    df["cid"] = df["protein"].astype(str)
    if ids is not None:
        df = df[df["cid"].isin(ids)].copy()
    out = {}
    for key, per_pose_name in per_pose_method.items():
        sub = df[df["method"] == per_pose_name]
        sub = _filter_optimizer_variant(sub, key, optimizer_by_method)
        if per_pose_name == "diffdock" and "pose_file" in sub.columns:
            # Base DiffDock writes its top pose under TWO filenames: rankN.sdf and
            # rankN_confidence-X.sdf (byte-identical). Drop the no-confidence copies
            # so each physical pose is counted once. (The smina/gnina rescored variants
            # carry a tool suffix, so they have no such twin.)
            base = sub["pose_file"].astype(str).map(os.path.basename)
            sub = sub[~base.str.match(r"rank\d+\.sdf$")]
        g = sub.groupby("cid").agg(generated=("pb_valid", "size"),
                                   pb_valid=("pb_valid", "sum"),
                                   near2=("near2", "sum"),
                                   form1=("form1", "sum"))
        out[key] = g
    return out


def load_pose_counts_orai(csv_path: Path, ids: Optional[set],
                          per_pose_method: Dict[str, str],
                          optimizer_by_method: Optional[Dict[str, str]] = None) -> Dict[str, pd.DataFrame]:
    """Orai per-method pose counts (generated + pb_valid) indexed by canonical cid.

    The Orai PoseBusters CSV carries no ``pb_valid`` / ``method`` columns and no crystal
    RMSD. Reuse posebusters_validity_report.load_and_score — the single source of truth for
    deriving pb_valid from the canonical PoseBusters test set AND splitting docking_method
    into per-variant labels (equibind_unguided_smina, diffdock_smina, …) — then key each
    complex by (frame, ligand). near2/form1 stay 0: without a crystal those tiers don't exist
    (main() swaps in ORAI_QUALITY_TIERS so they are never plotted).
    """
    try:
        import sys as _sys
        _here = str(Path(__file__).resolve().parent)
        if _here not in _sys.path:
            _sys.path.insert(0, _here)
        from posebusters_validity_report import load_and_score
    except Exception as e:                                   # pragma: no cover
        raise SystemExit(
            "ERROR: could not import posebusters_validity_report.load_and_score "
            f"({e!r}). The Orai path needs it to derive pb_valid + method labels; run in "
            "the 'vina' conda env where run_posebusters (PoseBusters) is importable.")
    df = load_and_score(Path(csv_path))
    df["cid"] = [_orai_cid(l, p) for l, p in zip(df["ligand"].astype(str),
                                                 df["protein"].astype(str))]
    if ids is not None:
        df = df[df["cid"].isin(ids)].copy()
    out = {}
    for key, per_pose_name in per_pose_method.items():
        sub = df[df["docking_method"] == per_pose_name]
        sub = _filter_optimizer_variant(sub, key, optimizer_by_method)
        if per_pose_name == "diffdock" and "pose_file" in sub.columns:
            # Raw DiffDock writes its top pose twice (rankN.sdf + rankN_confidence-X.sdf);
            # drop the no-confidence twin so each pose counts once (rescored variants are
            # suffix-tagged and have no twin). Mirrors load_pose_counts.
            base = sub["pose_file"].astype(str).map(os.path.basename)
            sub = sub[~base.str.match(r"rank\d+\.sdf$")]
        g = sub.groupby("cid").agg(generated=("pb_valid", "size"),
                                   pb_valid=("pb_valid", "sum"))
        g["near2"] = 0
        g["form1"] = 0
        out[key] = g
    return out


def _collect_rows_benchmark(args, counts) -> tuple:
    """Per-method per-complex timing+count rows for the benchmark (one file per complex)."""
    ad_root = Path(args.autodock_dir)
    ud_root = Path(args.unidock_dir) if args.unidock_dir else None
    ud2_root = Path(args.unidock2_dir) if args.unidock2_dir else None
    dd_root = Path(args.diffdock_dir)
    eb_root = Path(args.equibind_dir)
    rows: List[dict] = []
    per_method_timed = {k: 0 for k in METHODS}
    eq_processed = eq_best = 0
    for key in METHODS:
        for cid, crow in counts[key].iterrows():
            wall = wall_full = cpu_core_s = gpu_s = None
            if key == "autodock":
                effort = _autodock_effort(
                    cid, ad_root, args.autodock_prep,
                    getattr(args, "autodock_refine", "raw"), args.autodock_cpu,
                    getattr(args, "autodock_opt_cpu", 1),
                    getattr(args, "autodock_gnina_gpu", False),
                    getattr(args, "autodock_optimizer_workers", None),
                    getattr(args, "autodock_gnina_accounting", DEFAULT_GNINA_ACCOUNTING),
                    getattr(args, "autodock_gnina_host_cores", DEFAULT_GNINA_HOST_CORES),
                )
                if effort is not None:
                    wall, gpu_s, cpu_core_s = effort
                    wall_full = wall
            elif key in {"unidock", "unidock2"}:
                root = ud_root if key == "unidock" else ud2_root
                effort = _unidock_effort(cid, root, key) if root is not None else None
                if effort is not None:
                    wall, gpu_s, cpu_core_s = effort
                    wall_full = wall
            elif key == "diffdock":
                dd = _diffdock_effort(cid, dd_root, args.diffdock_refine, args.diffdock_gnina_gpu)
                if dd is not None:
                    wall, gpu_s, cpu_core_s = dd
                    wall_full = wall
            elif key == "equibind":
                eq_processed += 1
                t = _equibind_time(cid, eb_root, args.eq_mode, args.eq_refine, args.eq_clamp,
                                   args.eq_gnina_gpu)
                if t:
                    wall_full, wall, gpu_s, cpu_core_s = t[0], t[1], t[2], t[3]
                    if (gpu_s + cpu_core_s) > 0:
                        eq_best += 1
            if wall is None or not np.isfinite(wall) or wall <= 0:
                continue
            per_method_timed[key] += 1
            rows.append({
                "method": key, "cid": cid,
                "wall_s": float(wall), "wall_full_s": float(wall_full),
                "cpu_core_s": float(cpu_core_s), "gpu_s": float(gpu_s),
                "generated": int(crow["generated"]), "pb_valid": int(crow["pb_valid"]),
                "near2": int(crow["near2"]), "form1": int(crow["form1"]),
            })
    return rows, per_method_timed, eq_processed, eq_best


def _collect_rows_orai(args, counts) -> tuple:
    """Per-method per-complex timing+count rows for Orai (aggregated logs, preloaded once)."""
    ad = _orai_autodock_timing(Path(args.autodock_dir))
    ud_root = Path(args.unidock_dir) if args.unidock_dir else None
    ud2_root = Path(args.unidock2_dir) if args.unidock2_dir else None
    ad_refine = getattr(args, "autodock_refine", "raw")
    ad_opt = (_orai_autodock_opt_timing(Path(args.autodock_dir), ad_refine)
              if ad_refine in _AUTODOCK_REFINERS else {})
    dd = _orai_diffdock_timing(Path(args.diffdock_dir), args.diffdock_refine,
                               args.diffdock_gnina_gpu)
    eb = _orai_equibind_timing(Path(args.equibind_dir), args.eq_mode, args.eq_refine,
                               args.eq_clamp, args.eq_gnina_gpu)
    rows: List[dict] = []
    per_method_timed = {k: 0 for k in METHODS}
    eq_processed = eq_best = 0
    for key in METHODS:
        for cid, crow in counts[key].iterrows():
            wall = wall_full = cpu_core_s = gpu_s = None
            if key == "autodock":
                w = ad.get(cid)
                if w is not None:
                    opt = ad_opt.get(cid, 0.0) if ad_refine in _AUTODOCK_REFINERS else 0.0
                    if ad_refine in _AUTODOCK_REFINERS and opt <= 0:
                        continue
                    wall = wall_full = w + opt
                    cpu_core_s = w * args.autodock_cpu
                    if (ad_refine in _AUTODOCK_GNINA_REFINERS
                            and getattr(args, "autodock_gnina_gpu", False)):
                        # Orai has ONE aggregated optimiser log and no per-invocation sidecars,
                        # so device-occupancy cannot be measured here. Stay process-sum, say so.
                        if getattr(args, "autodock_gnina_accounting",
                                   DEFAULT_GNINA_ACCOUNTING) == "device-occupancy":
                            _warn_orai_process_sum_once()
                        gpu_s = opt
                    else:
                        gpu_s = 0.0
                        cpu_core_s += opt * max(int(getattr(args, "autodock_opt_cpu", 1)), 1)
            elif key in {"unidock", "unidock2"}:
                root = ud_root if key == "unidock" else ud2_root
                effort = _unidock_effort(cid, root, key) if root is not None else None
                if effort is not None:
                    wall, gpu_s, cpu_core_s = effort
                    wall_full = wall
            elif key == "diffdock":
                e = dd.get(cid)
                if e is not None:
                    wall, gpu_s, cpu_core_s = e
                    wall_full = wall
            elif key == "equibind":
                eq_processed += 1
                e = eb.get(cid)
                if e is not None:
                    wall_full, wall, gpu_s, cpu_core_s = e
                    if (gpu_s + cpu_core_s) > 0:
                        eq_best += 1
            if wall is None or not np.isfinite(wall) or wall <= 0:
                continue
            per_method_timed[key] += 1
            rows.append({
                "method": key, "cid": cid,
                "wall_s": float(wall), "wall_full_s": float(wall_full),
                "cpu_core_s": float(cpu_core_s), "gpu_s": float(gpu_s),
                "generated": int(crow["generated"]), "pb_valid": int(crow["pb_valid"]),
                "near2": int(crow["near2"]), "form1": int(crow["form1"]),
            })
    return rows, per_method_timed, eq_processed, eq_best


def build_table(args) -> tuple:
    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    # Benchmark: per-complex timing files + per_pose_metrics (crystal RMSD present).
    # Orai: aggregated per-tool timing logs + PoseBusters CSV (no crystal → pb_valid derived).
    if args.dataset == "orai_benchmark":
        counts = load_pose_counts_orai(
            Path(args.per_pose_csv), ids, _per_pose_methods(args), _per_pose_optimizers(args))
        rows, per_method_timed, eq_processed, eq_best = _collect_rows_orai(args, counts)
    else:
        counts = load_pose_counts(
            Path(args.per_pose_csv), ids, _per_pose_methods(args), _per_pose_optimizers(args),
            rmsd_column=_endpoint_cols(args)[0], kabsch_column=_endpoint_cols(args)[1])
        rows, per_method_timed, eq_processed, eq_best = _collect_rows_benchmark(args, counts)

    if eq_processed and not eq_best:
        raise SystemExit(
            "ERROR: the EquiBind best-config filter matched no poses for any complex "
            f"(--eq-mode={args.eq_mode} --eq-refine={args.eq_refine} --eq-clamp={args.eq_clamp}).\n"
            "       Check these against the mode / refine_variant / clamp_variant values in "
            "pipeline_summary.json — clamp_variant is None for clamp-off runs, so leave "
            "--eq-clamp unset (do NOT pass clampOFF).")
    pc = pd.DataFrame(rows)
    # Inside build_table so BOTH call sites are covered — main() and the
    # cross-dataset path, which re-enters with a cloned namespace. Note the
    # method column here holds engine SLOT keys (METHODS), not per-pose variant
    # keys, so only whole engines can be dropped at this level.
    pc = mf.apply_method_filter(pc, "method", args, label="docking-effort")

    # Compare methods on the COMMON set of complexes that have usable timing for
    # ALL methods, so total effort / pose counts / validity are like-for-like
    # (methods otherwise cover different subsets: skip-runs drop out per method).
    common = set()
    dropped_to_common = {}
    if not pc.empty:
        present = [k for k in METHODS if (pc["method"] == k).any()]
        common = set.intersection(*(set(pc.loc[pc["method"] == k, "cid"]) for k in present))
        dropped_to_common = {k: per_method_timed[k] - len(common) for k in present}
        pc = pc[pc["cid"].isin(common)].copy()
    meta = {"common_n": len(common), "per_method_timed": per_method_timed,
            "dropped_to_common": dropped_to_common}

    summary_rows = []
    for key, (pretty, _c, device) in METHODS.items():
        m = pc[pc["method"] == key]
        if m.empty:
            continue
        tot_wall = float(m["wall_s"].sum())
        # A missing CPU/GPU phase split is unknown, not zero work. In particular
        # Uni-Dock2 currently records only end-to-end wall time.
        tot_cpu = (float(m["cpu_core_s"].sum())
                   if m["cpu_core_s"].notna().all() else float("nan"))
        tot_gpu = (float(m["gpu_s"].sum())
                   if m["gpu_s"].notna().all() else float("nan"))
        n_gen = int(m["generated"].sum())
        n_val = int(m["pb_valid"].sum())
        # nested crystal tiers (0 for crystal-free Orai): near2 = PB-valid & RMSD ≤ 2 Å,
        # form1 = near2 & Kabsch ≤ 1 Å. Used for the per-quality resource panels.
        n_near2 = int(m["near2"].sum())
        n_form1 = int(m["form1"].sum())
        summary_rows.append({
            "method": key, "name": pretty, "device": device,
            "n_complexes": int(len(m)),
            "total_wall_h": round(tot_wall / 3600.0, 4),
            "total_wall_full_h": round(float(m["wall_full_s"].sum()) / 3600.0, 4),
            "total_cpu_core_h": round(tot_cpu / 3600.0, 4),
            "total_gpu_h": round(tot_gpu / 3600.0, 4),
            "poses_generated": n_gen,
            "poses_pb_valid": n_val,
            "poses_near2": n_near2,
            "poses_form1": n_form1,
            "pb_valid_rate_pct": round(100.0 * n_val / n_gen, 2) if n_gen else np.nan,
            # wall-clock per pose
            "wall_s_per_generated": round(tot_wall / n_gen, 4) if n_gen else np.nan,
            "wall_s_per_pb_valid": round(tot_wall / n_val, 4) if n_val else np.nan,
            # resource-seconds per pose (per generated / PB-valid / near-native / correct-form)
            "cpu_core_s_per_generated": round(tot_cpu / n_gen, 4) if n_gen else np.nan,
            "gpu_s_per_generated": round(tot_gpu / n_gen, 4) if n_gen else np.nan,
            "cpu_core_s_per_pb_valid": round(tot_cpu / n_val, 4) if n_val else np.nan,
            "gpu_s_per_pb_valid": round(tot_gpu / n_val, 4) if n_val else np.nan,
            "cpu_core_s_per_near2": round(tot_cpu / n_near2, 4) if n_near2 else np.nan,
            "gpu_s_per_near2": round(tot_gpu / n_near2, 4) if n_near2 else np.nan,
            "cpu_core_s_per_form1": round(tot_cpu / n_form1, 4) if n_form1 else np.nan,
            "gpu_s_per_form1": round(tot_gpu / n_form1, 4) if n_form1 else np.nan,
        })
    summ = pd.DataFrame(summary_rows)
    return pc, summ, meta


# ════════════════════════════════════════════════════════════════════════
# Statistics (paired across methods on the common timed set)
# ════════════════════════════════════════════════════════════════════════

_MIN_PAIRED = 3          # need at least this many paired complexes for a real test

# ── cost basis ──────────────────────────────────────────────────────────────
# Two ways to charge a pipeline for the time it spent, selected by --basis.
#
#   elapsed  wall_s as each reader built it. This is NOT one quantity across the
#            three arms: AutoDock's is vina elapsed plus its gnina sum (divided by
#            --autodock-optimizer-workers when that flag is given), DiffDock's is a
#            measured elapsed time, and EquiBind's is best-config compute divided by
#            the pipeline's own n_parallel_workers. Historical default; every
#            committed effort directory reproduces byte-for-byte under it.
#
#   charged  cpu_core_s / cpu_threads + gpu_s, applied identically to every arm.
#            Device occupancy: a processor stage that saturates all cpu_threads is
#            charged at its elapsed-equivalent, a GPU stage at its gpu_s (a true
#            device-occupancy time only under --autodock-gnina-accounting
#            device-occupancy), and consecutive stages sum. This is a time, not a resource-second,
#            so it does not violate the "CPU-core-s and GPU-s are never summed" rule
#            that governs the resource panels.
#
# The charged numerator reads only cpu_core_s and gpu_s, and neither is divided by
# --autodock-optimizer-workers (that flag rescales only wall_s), so the charged basis is
# arithmetically invariant to it under either gnina accounting. That invariance is not, on
# its own, an argument for the basis: under --autodock-gnina-accounting process-sum (the
# default) AutoDock's gpu_s is the undivided sum of ~9,000 gnina process elapsed times that
# ran sixteen-deep on ONE device, so it counts the same device-second up to sixteen times
# and is invariant to the worker count only because it ignores it. Under device-occupancy
# gpu_s is the measured union of each complex's own optimizer intervals, so concurrency is
# removed by measurement rather than by a divisor, and the stage's host CPU is billed at
# --autodock-gnina-host-cores instead of dropped.
COST_BASES = ("elapsed", "charged")
DEFAULT_COST_BASIS = "elapsed"
DEFAULT_COST_CPU_THREADS = 32


def _cost_numerator(m: pd.DataFrame, basis: str = DEFAULT_COST_BASIS,
                    cpu_threads: int = DEFAULT_COST_CPU_THREADS) -> np.ndarray:
    """Per-complex cost numerator in seconds for the rows of one method.

    ``elapsed`` returns wall_s untouched, so the default path is bit-identical to the
    pre-basis code. ``charged`` returns cpu_core_s / cpu_threads + gpu_s.
    """
    if basis == "charged":
        thr = max(int(cpu_threads), 1)
        cpu = pd.to_numeric(m["cpu_core_s"], errors="coerce").to_numpy(float)
        gpu = pd.to_numeric(m["gpu_s"], errors="coerce").to_numpy(float)
        return np.nan_to_num(cpu, nan=0.0) / thr + np.nan_to_num(gpu, nan=0.0)
    if basis != "elapsed":
        raise ValueError(f"unknown cost basis {basis!r}; expected one of {COST_BASES}")
    return m["wall_s"].to_numpy(float)


def _basis_label(basis: str) -> str:
    """Axis-label fragment naming the currency, for figures and sidecar units."""
    return ("Charged device-occupancy seconds" if basis == "charged"
            else "Wall-clock seconds")


def compute_effort_stats(pc: pd.DataFrame, summ: pd.DataFrame, meta: dict,
                         seed: int = 0, basis: str = DEFAULT_COST_BASIS,
                         cpu_threads: int = DEFAULT_COST_CPU_THREADS) -> dict:
    """Paired stats over the COMMON timed set (one value per complex per method).

    All methods dock the SAME complexes, so ``pc`` (already restricted to the
    common set in ``build_table``) is a fully paired design once pivoted to
    cid × method. Everything here degrades gracefully: if stats_utils is missing,
    too few paired complexes exist, or a test raises, the corresponding block is
    skipped with a note and the figures fall back to their test-free form.

    Returns a JSON-serialisable dict (also written to effort_stats.json):
      unit, common_n, wall_distribution{omnibus, pairwise, medians, median_ratios},
      validity{omnibus, pairwise, rates_paired, rates_pooled}, notes[].
    """
    _basis_unit = ("per-complex charged device-occupancy seconds "
                   f"(cpu_core_s / {int(cpu_threads)} + gpu_s)" if basis == "charged"
                   else "per-complex wall-clock")
    out = {
        "unit": f"{_basis_unit} / per-complex 'produced >=1 pb_valid pose' "
                f"boolean; paired across methods on the common timed set ({DATASET_LABEL}). "
                "Pooled per-pose validity rates are shown descriptively with Wilson CIs.",
        "cost_basis": basis,
        "cpu_threads": int(cpu_threads),
        "gnina_accounting": _gnina_accounting_label(),
        "common_n": int(meta.get("common_n", 0)),
        "min_paired": _MIN_PAIRED,
        "notes": [],
    }
    if su is None:
        out["notes"].append("stats_utils unavailable — all tests skipped (figures unannotated).")
        return out
    if pc is None or pc.empty:
        out["notes"].append("no per-complex effort rows — tests skipped.")
        return out

    order = [k for k in METHODS if (pc["method"] == k).any()]
    pretty = {k: METHODS[k][0] for k in order}

    # ── wall-clock distribution: Friedman + Wilcoxon/Holm + median-ratio CIs ──
    try:
        _cost = pc.assign(_cost_s=_cost_numerator(pc, basis, cpu_threads))
        wall = _cost.pivot_table(index="cid", columns="method", values="_cost_s")
        wall = wall[[k for k in order if k in wall.columns]].dropna()
        n_pair = int(len(wall))
        wd = {"n_paired": n_pair}
        if n_pair >= _MIN_PAIRED and len(order) >= 2:
            data = {pretty[k]: wall[k].to_numpy(float) for k in order}
            res = su.paired_continuous(data)
            wd["omnibus"] = res["omnibus"]
            wd["pairwise"] = res["pairwise"]
            wd["medians_s"] = res["medians"]
            # per-complex median ratio (paired) for each pair, higher / lower median
            ratios = []
            for k1 in order:
                for k2 in order:
                    if k1 >= k2:
                        continue
                    a, b = wall[k1].to_numpy(float), wall[k2].to_numpy(float)
                    m = (a > 0) & (b > 0)
                    if m.sum() < _MIN_PAIRED:
                        continue
                    # orient hi/lo by median so the ratio reads >= 1
                    if np.median(a[m]) >= np.median(b[m]):
                        hi, lo, hn, ln = a[m], b[m], pretty[k1], pretty[k2]
                    else:
                        hi, lo, hn, ln = b[m], a[m], pretty[k2], pretty[k1]
                    est, clo, chi = su.bootstrap_ci(hi / lo, np.median, seed=seed)
                    ratios.append({"slower": hn, "faster": ln,
                                   "median_ratio": est, "ci": [clo, chi], "n": int(m.sum())})
            wd["median_ratios"] = ratios
        else:
            wd["skipped"] = f"only {n_pair} paired complexes (<{_MIN_PAIRED}) — descriptive only."
            out["notes"].append("wall distribution: " + wd["skipped"])
        out["wall_distribution"] = wd
    except Exception as e:                                       # pragma: no cover
        out["notes"].append(f"wall distribution stats failed: {e!r}")

    # ── validity: per-complex 'produced a pb_valid pose' -> paired proportions ──
    try:
        val = pc.pivot_table(index="cid", columns="method", values="pb_valid")
        val = val[[k for k in order if k in val.columns]].dropna()
        n_pair = int(len(val))
        vd = {"n_paired": n_pair}
        if n_pair >= _MIN_PAIRED and len(order) >= 2:
            has_valid = {pretty[k]: (val[k].to_numpy(float) > 0).astype(int) for k in order}
            res = su.paired_proportions(has_valid)
            vd["omnibus"] = res["omnibus"]
            vd["pairwise"] = res["pairwise"]
            vd["rate_has_valid"] = {lab: {"rate": r, "ci": [lo, hi]}
                                    for lab, (r, lo, hi) in res["rates"].items()}
        else:
            vd["skipped"] = f"only {n_pair} paired complexes (<{_MIN_PAIRED}) — descriptive only."
            out["notes"].append("validity: " + vd["skipped"])
        # pooled per-pose validity rate (what the bars actually show) + Wilson CI
        pooled = {}
        for _, r in summ.iterrows():
            k = r["method"]
            n_gen, n_val = int(r["poses_generated"]), int(r["poses_pb_valid"])
            lo, hi = su.wilson_ci(n_val, n_gen)
            pooled[pretty.get(k, k)] = {"poses_generated": n_gen, "poses_pb_valid": n_val,
                                        "rate": (n_val / n_gen) if n_gen else float("nan"),
                                        "ci": [lo, hi],
                                        "note": "pooled over correlated poses — CI is "
                                                "descriptive, not the unit of inference"}
        vd["rate_pooled_per_pose"] = pooled
        out["validity"] = vd
    except Exception as e:                                       # pragma: no cover
        out["notes"].append(f"validity stats failed: {e!r}")

    return out


def _jsonable(o):
    """Recursively coerce numpy scalars / NaNs so a stats dict is json.dump-able."""
    if isinstance(o, dict):
        return {k: _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, float) and o != o:
        return None
    return o


def compute_effort_by_quality_stats(pc: pd.DataFrame, seed: int = 0,
                                    endpoint: Optional[tuple] = None,
                                    basis: str = DEFAULT_COST_BASIS,
                                    cpu_threads: int = DEFAULT_COST_CPU_THREADS) -> dict:
    """Paired cross-tool stats on the amortised docking cost per quality tier.

    For each nested tier the per-complex cost is wall_s(complex) / (poses reaching the
    tier); complexes yielding 0 such poses drop out of that tier. Per tier we run a
    paired Friedman + Wilcoxon/Holm across the tools that share a qualifying complex
    (plus bootstrapped paired median-ratio CIs), and report per-tool descriptive
    median/IQR over each tool's OWN qualifying complexes (what the boxes show). Degrades
    to descriptive-only when stats_utils is missing or too few complexes are paired.
    """
    _num = (f"charged device-occupancy seconds (cpu_core_s / {int(cpu_threads)} + gpu_s)"
            if basis == "charged" else "docking wall-clock seconds")
    _expr = (f"(cpu_core_s / {int(cpu_threads)} + gpu_s) / poses-in-tier"
             if basis == "charged" else "wall_s / poses-in-tier")
    out = {"unit": f"per-complex {_num} amortised over the poses "
                   f"reaching each quality tier ({_expr}); paired across "
                   f"tools on the common timed set ({DATASET_LABEL}). Descriptive per-tool "
                   "median/IQR are over each tool's own qualifying complexes; the paired "
                   "tests use only complexes every tool populates at that tier.",
           "cost_basis": basis, "cpu_threads": int(cpu_threads),
           "gnina_accounting": _gnina_accounting_label(),
           "min_paired": _MIN_PAIRED, "tiers": {}, "notes": []}
    # Self-document which per-pose endpoint the near2 / form1 tiers were scored on, so a
    # reader can tell at a glance whether this cost view matches the headline accuracy view.
    _rc, _kc = endpoint or (DEFAULT_RMSD_COLUMN, DEFAULT_KABSCH_COLUMN)
    out["endpoint"] = {
        "near2_rmsd_column": _rc, "near2_rmsd_threshold_A": NEAR2_RMSD_A,
        "form1_kabsch_column": _kc, "form1_kabsch_threshold_A": FORM1_KABSCH_A,
        "matches_headline_accuracy_report": bool(_rc == DEFAULT_RMSD_COLUMN
                                                 and _kc == DEFAULT_KABSCH_COLUMN),
    }
    if pc is None or pc.empty:
        out["notes"].append("no per-complex effort rows — tests skipped.")
        return out
    order = [k for k in METHODS if (pc["method"] == k).any()]
    pretty = {k: METHODS[k][0] for k in order}

    def _tool_cost(k, tcol):
        """Per-complex cost Series (cost numerator / poses-in-tier), count>0 only.

        Numerator is wall_s under the elapsed basis and cpu_core_s/threads + gpu_s
        under the charged one; see _cost_numerator.
        """
        m = pc[pc["method"] == k]
        cnt = m[tcol].to_numpy(float)
        keep = cnt > 0
        return pd.Series(_cost_numerator(m, basis, cpu_threads)[keep] / cnt[keep],
                         index=m["cid"].to_numpy()[keep])

    for tkey, tcol, _lab, _c in QUALITY_TIERS:
        cost = {k: _tool_cost(k, tcol) for k in order}
        td = {"per_method": {}}
        for k in order:
            v = cost[k].to_numpy(float)
            td["per_method"][pretty[k]] = {
                "n_complexes": int(v.size),
                "median_s": float(np.median(v)) if v.size else None,
                "q1_s": float(np.percentile(v, 25)) if v.size else None,
                "q3_s": float(np.percentile(v, 75)) if v.size else None,
            }
        wide = pd.DataFrame(cost).dropna()           # complexes every tool populates
        n = int(len(wide))
        td["n_paired"] = n
        if su is not None and n >= _MIN_PAIRED and len(order) >= 2:
            # Wrap per tier so one tier's failure (e.g. Friedman needs >=3 groups when
            # only two tools have timing) degrades THAT tier to descriptive-only rather
            # than aborting the whole function and losing every tier's sidecar + brackets.
            try:
                data = {pretty[k]: wide[k].to_numpy(float) for k in order}
                res = su.paired_continuous(data)
                td["omnibus"] = res["omnibus"]
                td["pairwise"] = res["pairwise"]
                td["medians_paired_s"] = res["medians"]
                ratios = []
                for a in range(len(order)):
                    for b in range(a + 1, len(order)):
                        x1, x2 = wide[order[a]].to_numpy(float), wide[order[b]].to_numpy(float)
                        if np.median(x1) >= np.median(x2):
                            hi, lo, hn, ln = x1, x2, pretty[order[a]], pretty[order[b]]
                        else:
                            hi, lo, hn, ln = x2, x1, pretty[order[b]], pretty[order[a]]
                        est, clo, chi = su.bootstrap_ci(hi / lo, np.median, seed=seed)
                        ratios.append({"costlier": hn, "cheaper": ln,
                                       "median_ratio": est, "ci": [clo, chi], "n": n})
                td["median_ratios"] = ratios
            except Exception as e:                       # e.g. <3 groups for Friedman
                td["skipped"] = f"paired test failed ({e!r}) — descriptive only."
                out["notes"].append(f"tier {tkey}: " + td["skipped"])
        else:
            reason = ("stats_utils unavailable" if su is None
                      else f"only {n} paired complexes (<{_MIN_PAIRED})")
            td["skipped"] = reason + " — descriptive only."
            out["notes"].append(f"tier {tkey}: " + td["skipped"])
        out["tiers"][tkey] = td
    return out


def _pretty_pos(order):
    """{pretty method name -> x position} for placing star brackets."""
    return {METHODS[k][0]: i for i, k in enumerate(order)}


def _sig_brackets(ax, pairwise, pos, log=False, show_ns=True):
    """Stack star brackets over paired comparisons without overplotting.

    pairwise : list of {a, b, star} (pretty names); pos maps name->x index.
    Reserves headroom at the top of the axis, then draws each bracket at a fixed
    axis-fraction height so log/linear scales are handled uniformly.
    """
    pairs = []
    for pr in pairwise:
        star = pr.get("star", "")
        if not star or (star == "ns" and not show_ns):
            continue
        if pr["a"] not in pos or pr["b"] not in pos:
            continue
        i, j = sorted((pos[pr["a"]], pos[pr["b"]]))
        pairs.append((i, j, star))
    if not pairs:
        return
    pairs.sort(key=lambda t: (t[1] - t[0], t[0]))               # nested first
    # headroom so the highest bracket sits above the data
    frac = max(0.55, 1.0 - 0.09 * (len(pairs) + 1))
    lo, hi = ax.get_ylim()
    if log and lo > 0 and hi > 0:
        llo, lhi = np.log10(lo), np.log10(hi)
        ax.set_ylim(lo, 10 ** (llo + (lhi - llo) / frac))
    else:
        ax.set_ylim(lo, lo + (hi - lo) / frac)
    trans = ax.get_xaxis_transform()
    for lvl, (i, j, star) in enumerate(pairs):
        y = 0.80 + 0.075 * lvl
        ax.plot([i, i, j, j], [y, y + 0.012, y + 0.012, y], transform=trans,
                color="black", lw=0.8, clip_on=False)
        ax.text((i + j) / 2.0, y + 0.016, star, transform=trans, ha="center",
                va="bottom", fontsize=8, clip_on=False)


def _omnibus_note(omni, kind):
    """Compact omnibus subtitle line for a title, reusing su.fmt_p/p_stars."""
    if not omni or su is None:
        return ""
    p = omni.get("p")
    if p is None or p != p:
        return ""
    if kind == "friedman":
        return (f"Friedman chi2={omni.get('chi2', float('nan')):.1f}, {su.fmt_p(p)} "
                f"{su.p_stars(p)}, Kendall W={omni.get('kendall_w', float('nan')):.2f} "
                f"(n={omni.get('n', '?')} paired)")
    return (f"Cochran Q={omni.get('Q', float('nan')):.1f}, {su.fmt_p(p)} "
            f"{su.p_stars(p)} (n paired)")


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _bars(ax, x, ys, colors, fmt, devices=None, ymax_mult=1.25):
    ax.bar(x, ys, color=colors)
    ymax = max([v for v in ys if v == v] + [0]) or 1
    for i, (xi, y) in enumerate(zip(x, ys)):
        if y == y:
            lab = fmt.format(y) + (f"\n[{devices[i]}]" if devices else "")
            ax.text(xi, y + 0.01 * ymax, lab, ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * ymax_mult)


def _save_panels(panels, panel_dir, plt, figsize=(6.2, 4.8)):
    """Render each (name, draw_fn) as its own single-panel PNG (no (a)/(b) label).

    Panels whose draw_fn accepts ``standalone`` get it True here so they render their
    own between-title-and-axes legend (in the combined grids that legend is figure-level).
    """
    import inspect
    panel_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for nm, draw in panels:
        f, a = plt.subplots(figsize=figsize)
        if "standalone" in inspect.signature(draw).parameters:
            draw(a, standalone=True)
        else:
            draw(a)
        f.tight_layout()
        pp = panel_dir / f"{nm}.png"
        f.savefig(pp, dpi=150, bbox_inches="tight"); plt.close(f)
        written.append(pp)
    return written


def make_wall_figure(pc, summ, out_dir, stats=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    colors = [METHODS[k][1] for k in order]
    devices = [METHODS[k][2] for k in order]
    # neutral-grey swatches for the generated/survive split (colour still = tool)
    _counts_labels = ["generated", "survive PoseBusters"]
    _counts_handles = [Patch(facecolor=LEGEND_GREY, hatch=HATCH_SECONDARY, edgecolor="white",
                             label=_counts_labels[0]),
                       Patch(facecolor=LEGEND_GREY, label=_counts_labels[1])]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order))
    stats = stats or {}
    pos = _pretty_pos(order)
    _wd = stats.get("wall_distribution", {}) or {}
    _vd = stats.get("validity", {}) or {}

    def _style(ax):
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25)

    def p_overall(ax):
        _bars(ax, x, s["total_wall_h"].to_numpy(float), colors, "{:.2f} h", devices)
        ax.set_title("Overall computational effort (wall-clock to dock the benchmark)\n"
                     "EquiBind = best-config-only estimate (full fan-out wall in CSV)")
        ax.set_ylabel("Total wall-clock time (hours)"); _style(ax)

    def p_distribution(ax):
        data = [pc[pc["method"] == k]["wall_s"].to_numpy(float) for k in order]
        bp = ax.boxplot(data, positions=x, widths=0.6, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.75)
        ax.set_yscale("log")
        title = ("Per-complex docking wall-clock time (distribution)\n"
                 "EquiBind = best-config-only estimate (compute ÷ workers)")
        note = _omnibus_note(_wd.get("omnibus"), "friedman")
        if note:
            title += "\n" + note + "  [paired Friedman; Wilcoxon/Holm brackets]"
        ax.set_title(title)
        ax.set_ylabel("Wall-clock time per complex (seconds, log scale)"); _style(ax)
        try:
            if _wd.get("pairwise"):
                _sig_brackets(ax, _wd["pairwise"], pos, log=True)
        except Exception as _e:                                  # pragma: no cover
            print(f"  [warn] wall-distribution brackets skipped: {_e!r}")

    def p_counts(ax, standalone=False):
        w = 0.38
        gen = s["poses_generated"].to_numpy(float)
        val = s["poses_pb_valid"].to_numpy(float)
        # colour = tool; hatch marks "generated", solid marks "survive" (house style)
        ax.bar(x - w / 2, gen, w, color=colors, hatch=HATCH_SECONDARY, edgecolor="white")
        ax.bar(x + w / 2, val, w, color=colors)
        for i in range(len(order)):
            rate = s["pb_valid_rate_pct"].to_numpy(float)[i]
            if rate == rate:
                ax.text(x[i] + w / 2, val[i], f"{rate:.1f}%", ha="center",
                        va="bottom", fontsize=8)
        title = "Poses generated vs. poses that survive\nthe PoseBusters tests"
        note = _omnibus_note(_vd.get("omnibus"), "cochran")
        if note:
            title += ("\nper-complex 'produced a valid pose': " + note
                      + " [McNemar/Holm brackets]")
        ax.set_title(title, pad=26 if standalone else None)
        ax.set_ylabel("Number of poses"); _style(ax)
        if standalone:
            ax.legend(handles=_counts_handles, loc="lower center",
                      bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, frameon=True)
        try:
            if _vd.get("pairwise"):
                _sig_brackets(ax, _vd["pairwise"], pos, log=False)
        except Exception as _e:                                  # pragma: no cover
            print(f"  [warn] poses-count brackets skipped: {_e!r}")

    def p_per_generated(ax):
        _bars(ax, x, s["wall_s_per_generated"].to_numpy(float), colors, "{:.2f} s")
        ax.set_title("Timing per generated pose")
        ax.set_ylabel("Wall-clock time per generated pose (seconds)"); _style(ax)

    def p_per_valid(ax):
        _bars(ax, x, s["wall_s_per_pb_valid"].to_numpy(float), colors, "{:.2f} s")
        ax.set_title("Timing per pose that survives the PoseBusters tests")
        ax.set_ylabel("Wall-clock time per valid pose (seconds)"); _style(ax)

    def p_rate(ax):
        _bars(ax, x, s["pb_valid_rate_pct"].to_numpy(float), colors, "{:.1f}%", ymax_mult=1.3)
        title = "PoseBusters validity rate\n(fraction of generated poses that survive)"
        note = _omnibus_note(_vd.get("omnibus"), "cochran")
        if note:
            title += ("\npaired on per-complex 'produced a valid pose': " + note)
        ax.set_title(title)
        ax.set_ylabel("Valid poses (percent of generated)")
        ax.set_ylim(0, min(100, ax.get_ylim()[1])); _style(ax)
        try:
            if _vd.get("pairwise"):
                _sig_brackets(ax, _vd["pairwise"], pos, log=False)
        except Exception as _e:                                  # pragma: no cover
            print(f"  [warn] validity-rate brackets skipped: {_e!r}")

    panels = [
        ("wall_overall_total", p_overall),
        ("wall_per_complex_distribution", p_distribution),
        ("wall_poses_generated_vs_valid", p_counts),
        ("wall_per_generated_pose", p_per_generated),
        ("wall_per_valid_pose", p_per_valid),
        ("wall_validity_rate", p_rate),
    ]

    # combined 2×3 grid (kept — the notebook displays this file)
    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    for a, (_n, draw) in zip(ax.ravel(), panels):
        draw(a)
    _label_panels(ax)
    fig.suptitle("Docking WALL-CLOCK (time-to-result) — AutoDock Vina vs DiffDock vs "
                 f"EquiBind (best config), {DATASET_LABEL} (common timed set).\n"
                 "AutoDock = CPU (multi-thread), DiffDock/EquiBind = single GPU — comparable "
                 "as elapsed time, not as identical-hardware work.", fontsize=12, y=0.995)
    _legend_between(fig, _counts_handles, _counts_labels, y=0.945, ncol=2, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.915))
    p = out_dir / "docking_effort_comparison.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    return [p] + _save_panels(panels, out_dir / "panels", plt)


def make_resource_figure(summ, out_dir):
    """CPU-core-seconds vs GPU-seconds — the true hardware work (different currencies)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    tool_colors = [METHODS[k][1] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order)); w = 0.38

    def _res_handles(labels):                          # neutral swatches: solid=CPU, hatch=GPU
        return [Patch(facecolor=LEGEND_GREY, label=labels[0]),
                Patch(facecolor=LEGEND_GREY, hatch=HATCH_SECONDARY, edgecolor="white",
                      label=labels[1])]

    def grouped(ax, cpu, gpu, fmt, labels=("CPU-core-seconds", "GPU-seconds"), standalone=False):
        # colour = tool; CPU bar solid, GPU bar hatched (house style)
        ax.bar(x - w / 2, cpu, w, color=tool_colors)
        ax.bar(x + w / 2, gpu, w, color=tool_colors, hatch=HATCH_SECONDARY, edgecolor="white")
        finite_vals = [float(v) for v in list(cpu) + list(gpu) if np.isfinite(v)]
        ymax = max(finite_vals + [0.0]) or 1
        for off, vals in ((-w / 2, cpu), (w / 2, gpu)):
            for xi, v in zip(x + off, vals):
                if np.isfinite(v):
                    ax.text(xi, v + 0.01 * ymax, fmt.format(v), ha="center",
                            va="bottom", fontsize=8, rotation=0)
        ax.set_ylim(0, ymax * 1.25)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25)
        if standalone:
            ax.legend(handles=_res_handles(labels), loc="lower center",
                      bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=8, frameon=True)

    def p_overall(ax, standalone=False):
        grouped(ax, s["total_cpu_core_h"].to_numpy(float), s["total_gpu_h"].to_numpy(float),
                "{:.1f}", labels=("CPU-core-hours", "GPU-hours"), standalone=standalone)
        ax.set_title("Overall hardware effort\n(total resource consumed)",
                     pad=26 if standalone else None)
        ax.set_ylabel("Resource used (hours): CPU-core-hours | GPU-hours")

    def p_per_generated(ax, standalone=False):
        grouped(ax, s["cpu_core_s_per_generated"].to_numpy(float),
                s["gpu_s_per_generated"].to_numpy(float), "{:.2f}", standalone=standalone)
        ax.set_title("Hardware effort per generated pose", pad=26 if standalone else None)
        ax.set_ylabel("Resource per generated pose (seconds)")

    def p_per_valid(ax, standalone=False):
        grouped(ax, s["cpu_core_s_per_pb_valid"].to_numpy(float),
                s["gpu_s_per_pb_valid"].to_numpy(float), "{:.2f}", standalone=standalone)
        ax.set_title("Hardware effort per pose that survives\nthe PoseBusters tests",
                     pad=26 if standalone else None)
        ax.set_ylabel("Resource per valid pose (seconds)")

    def p_per_near2(ax, standalone=False):
        grouped(ax, s["cpu_core_s_per_near2"].to_numpy(float),
                s["gpu_s_per_near2"].to_numpy(float), "{:.2f}", standalone=standalone)
        ax.set_title("Hardware effort per pose that survives PoseBusters\n"
                     "AND is near-native (RMSD ≤ 2 Å)", pad=26 if standalone else None)
        ax.set_ylabel("Resource per near-native valid pose (seconds)")

    def p_per_form1(ax, standalone=False):
        grouped(ax, s["cpu_core_s_per_form1"].to_numpy(float),
                s["gpu_s_per_form1"].to_numpy(float), "{:.2f}", standalone=standalone)
        ax.set_title("Hardware effort per pose that survives PoseBusters, is\n"
                     "near-native (RMSD ≤ 2 Å) AND correct form (Kabsch ≤ 1 Å)",
                     pad=26 if standalone else None)
        ax.set_ylabel("Resource per near-native correct-form pose (seconds)")

    # Core 1×3 panels (the notebook displays the combined file). The crystal-only quality
    # tiers (near-native, correct-form) exist only when a crystal reference is present
    # (benchmark, not Orai) → added as standalone panels when any tool populates them.
    core_panels = [
        ("resource_overall_total", p_overall),
        ("resource_per_generated_pose", p_per_generated),
        ("resource_per_valid_pose", p_per_valid),
    ]
    extra_panels = []
    if (s.get("poses_near2", pd.Series(0.0, index=s.index)).fillna(0) > 0).any():
        extra_panels.append(("resource_per_near_native_valid_pose", p_per_near2))
    if (s.get("poses_form1", pd.Series(0.0, index=s.index)).fillna(0) > 0).any():
        extra_panels.append(("resource_per_near_native_correct_form_pose", p_per_form1))

    # combined 1×3 grid (kept — the notebook displays this file)
    fig, ax = plt.subplots(1, 3, figsize=(18, 6.4))
    for a, (_n, draw) in zip(ax, core_panels):
        draw(a)
    _label_panels(ax)
    _dd_nm = METHODS["diffdock"][0]
    _eq_nm = METHODS["equibind"][0]
    fig.suptitle("Docking HARDWARE EFFORT — CPU-core-seconds vs GPU-seconds (the actual work, "
                 f"parallelism removed) — {DATASET_LABEL}.\nA GPU-second is NOT a CPU-core-second — "
                 "the two bars are different currencies and are never summed. AutoDock is CPU-only; "
                 f"{_dd_nm} docks on the GPU (rescoring on GPU or CPU per refiner); {_eq_nm} uses the "
                 "GPU only for inference — conformer prep/UFF and its refine run on CPU, so it is "
                 "CPU-dominated.", fontsize=12, y=0.995)
    _legend_between(fig, _res_handles(("CPU-core work (solid)", "GPU work (hatched)")),
                    ["CPU-core work (solid)", "GPU work (hatched)"], y=0.86, ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.82))
    p = out_dir / "docking_resource_effort.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    return [p] + _save_panels(core_panels + extra_panels, out_dir / "panels", plt)


def make_time_per_pose_figure(summ, out_dir):
    """Wall-clock time to produce ONE pose vs ONE PoseBusters-valid pose, per toolchain,
    with each bar split into the share of that time delivered by CPU vs GPU work.

    Two grouped bars per method (per generated pose, per surviving pose); each bar is
    stacked into CPU-delivered and GPU-delivered wall-seconds, apportioned by the
    method's CPU-core-second : GPU-second ratio (so the stack total is the per-pose time).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    tool_colors = [METHODS[k][1] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order)); w = 0.38

    cpu_h = s["total_cpu_core_h"].to_numpy(float)
    gpu_h = s["total_gpu_h"].to_numpy(float)
    tot_h = cpu_h + gpu_h
    known_split = np.isfinite(cpu_h) & np.isfinite(gpu_h) & (tot_h > 0)
    gpu_frac = np.divide(gpu_h, tot_h, out=np.zeros_like(gpu_h), where=known_split)
    cpu_frac = np.divide(cpu_h, tot_h, out=np.zeros_like(cpu_h), where=known_split)
    unknown_frac = (~known_split).astype(float)

    metrics = [(s["wall_s_per_generated"].to_numpy(float), -w / 2),
               (s["wall_s_per_pb_valid"].to_numpy(float), +w / 2)]

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    allv = [v for wall, _ in metrics for v in wall if v == v]
    ymax = (max(allv) if allv else 1) or 1
    for wall, off in metrics:
        cpu_part = cpu_frac * wall
        gpu_part = gpu_frac * wall
        unknown_part = unknown_frac * wall
        # colour = tool; CPU share solid, GPU share hatched (house style)
        ax.bar(x + off, cpu_part, w, color=tool_colors)
        ax.bar(x + off, gpu_part, w, bottom=cpu_part, color=tool_colors,
               hatch=HATCH_SECONDARY, edgecolor="white")
        ax.bar(x + off, unknown_part, w, bottom=cpu_part + gpu_part,
               color=tool_colors, hatch="xx", edgecolor="white")
        for xi, wv in zip(x + off, wall):
            if wv == wv:
                ax.text(xi, wv + 0.012 * ymax, f"{wv:.2f} s",
                        ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * 1.24)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("Wall-clock time per pose (seconds)")
    ax.set_title("Wall-clock time per pose, split into CPU- vs GPU-delivered work\n"
                 "left bar = per generated pose; right bar = per pose surviving PoseBusters",
                 pad=30)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(handles=[Patch(facecolor=LEGEND_GREY, label="CPU-delivered"),
                       Patch(facecolor=LEGEND_GREY, hatch=HATCH_SECONDARY, edgecolor="white",
                             label="GPU-delivered"),
                       Patch(facecolor=LEGEND_GREY, hatch="xx", edgecolor="white",
                             label="CPU/GPU split unavailable")],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=9, frameon=True)
    fig.tight_layout()
    p = out_dir / "docking_time_per_pose.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p


def _fmt_secs(v):
    """Compact seconds label that keeps sub-second per-pose costs legible (0.22, not 0).

    A plain ``{:.0f}`` collapses AutoDock's ~0.22 s/pose median to "0" on a log-scaled
    cost axis where zero is impossible — so scale the precision to the magnitude.
    """
    v = float(v)
    if v >= 10:
        return f"{v:.0f}"
    if v >= 1:
        return f"{v:.1f}"
    return f"{v:.2f}"


def make_effort_by_quality_figure(pc, out_dir, stats=None,
                                  basis=DEFAULT_COST_BASIS,
                                  cpu_threads=DEFAULT_COST_CPU_THREADS):
    """Box-and-whiskers of docking effort per pose reaching each quality tier, per tool.

    The cost numerator follows --basis: measured wall_s (elapsed) or
    cpu_core_s / cpu_threads + gpu_s (charged device occupancy).

    For every complex the tool's docking wall-clock is amortised over the poses that
    reach a tier (wall_s / poses-in-tier); each box is the distribution of that
    per-complex cost across the complexes that yield >=1 such pose. As the quality bar
    tightens (all poses -> PB-valid -> +RMSD<=2 Å -> +Kabsch<=1 Å) fewer poses qualify,
    so the seconds spent per good pose climb.

    Emits a combined grouped figure (tools on x, one box per tier) plus one single-panel
    PNG per tier (tools side by side, with the paired cross-tool test drawn on it).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [k for k in METHODS if (pc["method"] == k).any()]
    if not order:
        return []
    names = [METHODS[k][0] for k in order]
    tool_colors = [METHODS[k][1] for k in order]
    stats = stats or {}
    tier_stats = stats.get("tiers", {}) or {}
    pos = _pretty_pos(order)
    x = np.arange(len(order), dtype=float)

    def _costs(k, tcol):
        """Per-complex cost array (cost numerator / poses-in-tier), count>0 only."""
        m = pc[pc["method"] == k]
        cnt = m[tcol].to_numpy(float)
        keep = cnt > 0
        return _cost_numerator(m, basis, cpu_threads)[keep] / cnt[keep]

    # ── combined: grouped by tool, one box per tier ──────────────────────────
    # colour = tool (hue); quality tier = lightness ramp of that hue (light=loose→dark=strict)
    n_tier = len(QUALITY_TIERS)
    tier_off = (np.arange(n_tier) - (n_tier - 1) / 2.0) * 0.205
    tint_fracs = np.linspace(0.62, 0.0, n_tier)      # tier 0 lightest → last tier full tool colour
    bw = 0.17
    fig, ax = plt.subplots(figsize=(12.5, 6.8))
    med_labels = []                                  # (x, median, n) drawn after scaling
    for j, (tkey, tcol, _lab, _c) in enumerate(QUALITY_TIERS):
        for i, k in enumerate(order):
            v = _costs(k, tcol)
            if v.size == 0:
                continue
            posx = x[i] + tier_off[j]
            boxc = _tint(tool_colors[i], tint_fracs[j])
            bp = ax.boxplot([v], positions=[posx], widths=bw, patch_artist=True,
                            showfliers=True, medianprops=dict(color="black", lw=1.1),
                            flierprops=dict(marker="o", ms=2.5, mfc=boxc, mec="none",
                                            alpha=0.30))
            bp["boxes"][0].set_facecolor(boxc); bp["boxes"][0].set_alpha(0.95)
            med_labels.append((posx, float(np.median(v)), int(v.size)))
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(names, fontsize=10)
    ax.set_xlim(-0.7, len(order) - 0.3)
    _combo_y = ("Charged device-occupancy seconds per pose reaching the tier (log scale)"
                if basis == "charged"
                else "Docking wall-clock seconds per pose reaching the tier (log scale)")
    ax.set_ylabel(_combo_y)
    ax.grid(axis="y", which="both", alpha=0.2)
    trans = ax.get_xaxis_transform()
    for posx, med, nn in med_labels:
        ax.text(posx, med, _fmt_secs(med), ha="center", va="bottom", fontsize=6,
                rotation=90, color="0.15")
        ax.text(posx, 0.008, f"{nn}", transform=trans, ha="center", va="bottom",
                fontsize=5.5, color="0.4", clip_on=False)
    ax.text(-0.68, 0.008, "n cplx:", transform=trans, ha="left", va="bottom",
            fontsize=5.5, color="0.4", clip_on=False)
    # neutral-grey ramp shows the tier ordinal (hue itself = tool, read from the x-axis)
    handles = [Patch(facecolor=_tint("#3a3a3a", tint_fracs[j]), label=lab)
               for j, (_tk, _tc, lab, _c) in enumerate(QUALITY_TIERS)]
    if basis == "charged":
        _combo_title = (
            "Docking effort per pose by quality tier — charged device-occupancy seconds "
            f"to obtain one pose of increasing quality, per tool ({DATASET_LABEL})\n"
            f"per complex: (CPU-core-s ÷ {int(cpu_threads)} + GPU-s) ÷ poses reaching the "
            "tier; box = spread over complexes yielding ≥1 such pose (n below). One basis "
            "for every arm, so the three columns are directly comparable.")
    else:
        _combo_title = (
            "Docking effort per pose by quality tier — wall-clock seconds to obtain "
            f"one pose of increasing quality, per tool ({DATASET_LABEL})\n"
            "per complex: docking wall-clock ÷ poses reaching the tier; box = spread "
            "over complexes yielding ≥1 such pose (n below). AutoDock = CPU "
            "(multi-thread), DiffDock/EquiBind = single GPU — comparable as elapsed "
            "time, not identical-hardware work.")
    ax.set_title(_combo_title, fontsize=10, pad=46)
    ax.legend(handles=handles, title="quality tier — shade light→dark  (hue = tool)",
              fontsize=8, title_fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0),
              ncol=len(QUALITY_TIERS), framealpha=0.9)
    fig.tight_layout()
    p = out_dir / "docking_effort_by_quality_tier.png"
    fig.savefig(p, dpi=140, bbox_inches="tight"); plt.close(fig)

    # ── one single-panel PNG per tier (tools adjacent → paired-test brackets) ──
    def _panel(tkey, tcol, lab):
        def draw(a):
            drawn = [(i, _costs(order[i], tcol)) for i in range(len(order))]
            drawn = [(i, d) for i, d in drawn if d.size]
            if not drawn:
                a.text(0.5, 0.5, "no qualifying poses", ha="center", va="center",
                       transform=a.transAxes)
                return
            bp = a.boxplot([d for _i, d in drawn], positions=[x[i] for i, _d in drawn],
                           widths=0.55, patch_artist=True, showfliers=True,
                           medianprops=dict(color="black", lw=1.2),
                           flierprops=dict(marker="o", ms=3, mec="none", alpha=0.4))
            for (i, _d), patch in zip(drawn, bp["boxes"]):
                patch.set_facecolor(tool_colors[i]); patch.set_alpha(0.8)
            a.set_yscale("log")
            a.set_xticks(x); a.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
            a.set_xlim(-0.6, len(order) - 0.4)
            a.set_ylabel(f"{_basis_label(basis)} per qualifying pose (log)")
            a.grid(axis="y", which="both", alpha=0.2)
            tr = a.get_xaxis_transform()
            for i, d in drawn:
                a.text(x[i], float(np.median(d)), f"{_fmt_secs(np.median(d))} s",
                       ha="center", va="bottom", fontsize=8)
                a.text(x[i], 0.01, f"n={d.size}", transform=tr, ha="center", va="bottom",
                       fontsize=7, color="0.4")
            ts = tier_stats.get(tkey, {}) or {}
            title = f"Docking effort per pose — {lab}"
            if basis == "charged":
                # Name the currency on the panel itself. Under the elapsed basis the
                # title is left exactly as it was, so committed outputs reproduce.
                _tail = (" + GPU-s, gnina at measured device occupancy"
                         if _GNINA_ACCOUNTING_MODE == "device-occupancy" and not _ORAI_PROCESS_SUM_HIT
                         else " + GPU-s, one basis for every arm")
                title += f"\ncharged device occupancy: CPU-core-s / {int(cpu_threads)}" + _tail
            note = _omnibus_note(ts.get("omnibus"), "friedman")
            if note:
                title += "\n" + note + "  [paired Friedman; Wilcoxon/Holm]"
            a.set_title(title, fontsize=9)
            try:
                if ts.get("pairwise"):
                    _sig_brackets(a, ts["pairwise"], pos, log=True)
            except Exception as _e:                                  # pragma: no cover
                print(f"  [warn] effort-by-quality brackets ({tkey}) skipped: {_e!r}")
        return draw

    panels = [(f"effort_by_quality_{tkey}", _panel(tkey, tcol, lab))
              for tkey, tcol, lab, _c in QUALITY_TIERS]
    written = _save_panels(panels, out_dir / "panels", plt, figsize=(6.4, 5.0))
    return [p] + written


# ════════════════════════════════════════════════════════════════════════
# Batch-vs-per-run comparison (cross-dataset: benchmark per-run vs Orai batch)
# ════════════════════════════════════════════════════════════════════════
# Unlike every figure above (single dataset), this models BOTH campaigns to isolate
# the time that batch docking saves by amortising each tool's fixed startup over the
# ligands sharing a receptor. See the "Batch-vs-per-run model" block near the top for
# the L·(N_dockings − N_receptors) derivation.

def _measure_equibind_load(benchmark_equibind_dir) -> Optional[float]:
    """Median measured EquiBind per-run fixed load (seconds) from the benchmark run.

    Each benchmark complex writes its own pipeline_summary.json, so its global_timing
    captures the once-per-run cost that batch mode amortises: phase1 prep (receptor /
    ligand setup) plus the model load buried in phase2 (= phase2_dock − the poses' own
    inference, Σ dock_time_s). Returns the median over complexes, or None if unreadable
    (caller then falls back to DEFAULT_LOAD_EQUIBIND).
    """
    d = Path(benchmark_equibind_dir)
    if not d.exists():
        return None
    loads = []
    for p in d.glob("*/pipeline_summary.json"):
        try:
            js = json.loads(p.read_text())
        except Exception:
            continue
        gt = js.get("global_timing", {})
        p2 = gt.get("phase2_dock_s")
        if p2 is None:
            continue
        # Σ inference over UNIQUE poses (dock_time_s is duplicated across refine variants).
        seen = {}
        for r in js.get("all_results") or []:
            seen[(r.get("mode"), r.get("pocket_unique_id"), r.get("pose_num"))] = \
                float(r.get("dock_time_s") or 0.0)
        model_load = max(float(p2) - sum(seen.values()), 0.0)
        loads.append(model_load + float(gt.get("phase1_prep_s") or 0.0))
    return float(np.median(loads)) if loads else None


def _receptor_of(cid: str) -> str:
    """Receptor identity of a complex id, for counting distinct receptor-loads.

    Orai cids are ``<ligand>__<frame>`` → the frame is the receptor (≈300 ligands share
    each of the 4 frames). Benchmark cids are ``<PDBID>_<CCD>`` → the PDB id is the
    receptor (each is a unique native crystal, so reuse is 1×).
    """
    cid = str(cid)
    return cid.split("__")[-1] if "__" in cid else cid.split("_")[0]


def _clone_dataset_args(base, dataset: str):
    """A fresh args Namespace for one dataset: keep the user's shared flags, but reset the
    per-dataset path/refiner fields to their DATASET_DEFAULTS so each campaign reads its own
    logs. The Orai EquiBind run used gnina (its smina variant has no timed poses), so force
    the EquiBind refiner to gnina there — matching the benchmark default and the actual run.

    The reset is scoped to the campaign the user was NOT addressing. A blanket reset silently
    discarded the --autodock-dir / --per-pose-csv / --autodock-method the caller typed for the
    benchmark, so the cross-dataset figures costed the retired crystal-boxed Meeko campaign
    (7,337 poses / 1.19 h) while the very same invocation's effort_summary.csv reported the
    whole-protein exh128+gnina arm (8,976 poses / 5.37 h) — two AutoDock campaigns, one label,
    eleven seconds apart. Flags typed for THIS dataset survive; the other campaign still gets
    its own defaults, so an Orai clone can never inherit a benchmark-only path or method."""
    a = argparse.Namespace(**vars(base))
    a.dataset = dataset
    # out_dir is deliberately never kept: it is an output location, not a data source, and
    # a clone only ever reads. Keeping it would aim a clone at the primary run's directory.
    keep = frozenset()
    if getattr(base, _EXPLICIT_DS_ATTR, None) == dataset:
        keep = frozenset(getattr(base, _EXPLICIT_KEYS_ATTR, ())) - {"out_dir"}
    for k in PER_DATASET_ARGS:
        if k not in keep:
            setattr(a, k, None)
    # Not in the loop above because DATASET_DEFAULTS carries no such key: it is cleared
    # outright rather than refilled. The Orai per-pose CSV only ever spells docking_method
    # 'autodock', so a benchmark-only label inherited from the base args would match nothing
    # there and AutoDock would vanish from the batch / individual-vs-batch figures without a
    # warning. That is exactly the case `keep` excludes it from.
    if "autodock_method" not in keep:
        a.autodock_method = None
    _apply_dataset_defaults(a)
    if dataset == "orai_benchmark":
        a.eq_refine = "gnina"
    return a


def _dataset_batch_rows(a) -> pd.DataFrame:
    """Per-complex timing+count rows (pre-common-set) for one dataset, for the batch model.

    Uses each tool's OWN timed set — NOT the cross-tool common set build_table restricts to —
    because per-run vs batch is compared WITHIN a tool, so one tool's missing complexes must
    not shrink another's campaign.
    """
    ids = None
    if a.ids_file and Path(a.ids_file).exists():
        ids = {ln.strip() for ln in Path(a.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}
    if a.dataset == "orai_benchmark":
        counts = load_pose_counts_orai(
            Path(a.per_pose_csv), ids, _per_pose_methods(a), _per_pose_optimizers(a))
        rows, *_ = _collect_rows_orai(a, counts)
    else:
        counts = load_pose_counts(
            Path(a.per_pose_csv), ids, _per_pose_methods(a), _per_pose_optimizers(a),
            rmsd_column=_endpoint_cols(a)[0], kabsch_column=_endpoint_cols(a)[1])
        rows, *_ = _collect_rows_benchmark(a, counts)
    return pd.DataFrame(rows)


def _batch_model_records(pc: pd.DataFrame, dataset: str, loads: Dict[str, float]) -> List[dict]:
    """Per-tool per-run vs batch effort for one dataset (see the module-top model block).

    marginal M is the measured per-docking wall stripped of any startup it already bundles:
    a per-run campaign's wall (benchmark AutoDock/DiffDock) includes L on every pair, so L is
    subtracted back out before re-adding it per strategy; EquiBind's wall is a load-free
    compute÷workers estimate and Orai's batch wall is already amortised, so neither subtracts.
    """
    per_run_dataset = (dataset != "orai_benchmark")     # benchmark docked one fresh run per pair
    out = []
    for key in METHODS:
        m = pc[pc["method"] == key]
        if m.empty:
            continue
        L = float(loads.get(key, 0.0))
        n_dock = int(len(m))
        n_recep = int(m["cid"].map(_receptor_of).nunique())
        n_valid = int(m["pb_valid"].sum())
        n_gen = int(m["generated"].sum())
        # measured wall bundles L only for a per-run campaign of a model-reloading tool.
        bundles_L = per_run_dataset and key in ("autodock", "diffdock")
        marg = m["wall_s"].to_numpy(float) - (L if bundles_L else 0.0)
        M = float(np.clip(marg, 0.0, None).sum())
        per_run_s = M + L * n_dock
        batch_s = M + L * n_recep
        saving_s = per_run_s - batch_s               # == L·(n_dock − n_recep), confound-free
        out.append({
            "dataset": dataset, "method": key, "name": BASE_NAMES.get(key, key),
            "load_s": L, "n_dockings": n_dock, "n_receptors": n_recep,
            "reuse_factor": (n_dock / n_recep) if n_recep else np.nan,
            "poses_generated": n_gen, "poses_pb_valid": n_valid,
            "per_run_wall_h": per_run_s / 3600.0, "batch_wall_h": batch_s / 3600.0,
            "saving_wall_h": saving_s / 3600.0,
            "saving_pct": (100.0 * saving_s / per_run_s) if per_run_s > 0 else 0.0,
            "per_run_s_per_valid": (per_run_s / n_valid) if n_valid else np.nan,
            "batch_s_per_valid": (batch_s / n_valid) if n_valid else np.nan,
        })
    return out


def make_batch_vs_perrun_figure(base_args, out_dir: Path):
    """Cross-dataset figure: how much wall-clock batch docking saves vs a fresh run per pair.

    Builds the benchmark (per-run, unique crystal per ligand) and Orai (batch, 4 shared frames)
    campaigns, models each tool's per-run vs batch total under a fixed marginal compute, and
    draws total-hours and per-PB-valid-pose panels per dataset. Returns [combined_png, *panels]
    and writes batch_vs_perrun.csv. Empty [] if neither campaign could be built.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # Fixed per-run load L per tool (EquiBind measured from the benchmark phase timing).
    eq_L = base_args.load_equibind
    eq_src = "user"
    if eq_L is None:
        bench_eb = DATASET_DEFAULTS["benchmark"]["equibind_dir"]
        eq_L = _measure_equibind_load(bench_eb)
        eq_src = "measured" if eq_L is not None else "default"
        if eq_L is None:
            eq_L = DEFAULT_LOAD_EQUIBIND
    loads = {"autodock": base_args.load_autodock,
             "diffdock": base_args.load_diffdock,
             "equibind": float(eq_L)}
    load_src = {"autodock": "estimate", "diffdock": "estimate", "equibind": eq_src}

    # ── build both campaigns (skip a dataset whose logs are absent) ───────────
    datasets = [("benchmark", "PoseBuster benchmark\n(per-run: a unique crystal per ligand)"),
                ("orai_benchmark", "Orai × benchmark\n(batch: 4 shared receptor frames)")]
    records, titles = {}, {}
    for ds, ds_title in datasets:
        try:
            a = _clone_dataset_args(base_args, ds)
            pc = _dataset_batch_rows(a)
            if pc.empty:
                print(f"  [warn] batch comparison: no timed rows for {ds} — skipped.")
                continue
            records[ds] = _batch_model_records(pc, ds, loads)
            titles[ds] = ds_title
        except Exception as e:                                   # pragma: no cover
            print(f"  [warn] batch comparison: could not build {ds} ({e!r}) — skipped.")
    if not records:
        print("  [warn] batch comparison: no datasets built — figure skipped.")
        return []

    # ── flat CSV of every number behind the figure ───────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    flat = [r for ds in records for r in records[ds]]
    df_out = pd.DataFrame(flat)
    for k, v in loads.items():
        df_out.loc[df_out["method"] == k, "load_source"] = load_src[k]
    df_out.to_csv(out_dir / "batch_vs_perrun.csv", index=False)

    # ── panel drawer: grouped per-run vs batch bars, one group per tool ───────
    _batch_labels = ["per-run (reload each pair)", "batch (reload once per receptor)"]
    _batch_handles = [Patch(facecolor=LEGEND_GREY, label=_batch_labels[0]),
                      Patch(facecolor=LEGEND_GREY, hatch=HATCH_SECONDARY, edgecolor="white",
                            label=_batch_labels[1])]

    def _panel(ds, ycol, ylabel, fmt, title, annotate):
        recs = records[ds]
        names = [r["name"] for r in recs]
        bar_colors = [METHODS[r["method"]][1] for r in recs]     # colour = tool
        x = np.arange(len(recs)); w = 0.38

        def draw(ax, standalone=False):
            pr = np.array([r["per_run_" + ycol] for r in recs], float)
            ba = np.array([r["batch_" + ycol] for r in recs], float)
            # colour = tool; per-run solid, batch hatched (house style)
            ax.bar(x - w / 2, pr, w, color=bar_colors)
            ax.bar(x + w / 2, ba, w, color=bar_colors, hatch=HATCH_SECONDARY, edgecolor="white")
            ymax = max([v for v in np.r_[pr, ba] if v == v] + [0]) or 1
            for xi, v in zip(x - w / 2, pr):
                if v == v:
                    ax.text(xi, v + 0.01 * ymax, fmt.format(v), ha="center", va="bottom", fontsize=7)
            for xi, v in zip(x + w / 2, ba):
                if v == v:
                    ax.text(xi, v + 0.01 * ymax, fmt.format(v), ha="center", va="bottom", fontsize=7)
            if annotate:                                          # saving callout per tool
                for i, r in enumerate(recs):
                    d_h = r["saving_wall_h"]
                    if d_h < 1e-6:                                # explain WHY nothing is saved
                        why = ("no reuse" if r["reuse_factor"] <= 1.001    # unique receptors
                               else "no model load")                       # reuse, but L≈0 (AutoDock)
                        txt = f"{why} → 0 saved"
                    elif ycol == "wall_h":
                        txt = f"−{d_h:.2f} h ({r['saving_pct']:.0f}%)"
                    else:
                        d = r["per_run_s_per_valid"] - r["batch_s_per_valid"]
                        txt = f"−{d:.2f} s ({r['saving_pct']:.0f}%)"
                    ax.text(x[i], max(pr[i], ba[i]) + 0.06 * ymax, txt, ha="center",
                            va="bottom", fontsize=7, color="0.25", fontweight="bold")
            ax.set_ylim(0, ymax * 1.32)
            ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
            ax.set_ylabel(ylabel); ax.set_title(title, fontsize=10, pad=26 if standalone else None)
            ax.grid(axis="y", alpha=0.25)
            if standalone:
                ax.legend(handles=_batch_handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
                          ncol=2, fontsize=7.5, frameon=True)
        return draw

    order_ds = [d for d, _ in datasets if d in records]
    panels = []
    for ds in order_ds:
        panels.append((f"batch_total_hours_{ds}",
                       _panel(ds, "wall_h", "Total docking wall-clock (hours)", "{:.2f}",
                              titles[ds] + "\ntotal wall-clock: per-run vs batch", True)))
    for ds in order_ds:
        panels.append((f"batch_cost_per_valid_{ds}",
                       _panel(ds, "s_per_valid", "Wall-clock seconds per PB-valid pose", "{:.2f}",
                              titles[ds] + "\ncost per PB-valid pose: per-run vs batch", True)))

    # combined grid: rows = metric (total h, per valid), cols = dataset
    ncol = len(order_ds)
    fig, axes = plt.subplots(2, ncol, figsize=(6.4 * ncol, 10.4), squeeze=False)
    grid = [(0, ds, "wall_h", "Total docking wall-clock (hours)", "{:.2f}", "total wall-clock: per-run vs batch")
            for ds in order_ds] + \
           [(1, ds, "s_per_valid", "Wall-clock seconds per PB-valid pose", "{:.2f}",
             "cost per PB-valid pose: per-run vs batch") for ds in order_ds]
    for row, ds, ycol, ylab, fmt, sub in grid:
        col = order_ds.index(ds)
        _panel(ds, ycol, ylab, fmt, titles[ds] + "\n" + sub, True)(axes[row][col])
    _label_panels(axes)
    _Ls = "  ".join(f"{BASE_NAMES[k].split(' ')[0]} L={loads[k]:.1f}s [{load_src[k]}]"
                    for k in METHODS if k in loads)
    fig.suptitle(
        "Batch docking vs a fresh run per pair — wall-clock saved by loading each tool once "
        "per receptor instead of once per (receptor, ligand).\n"
        "Marginal docking compute is held fixed; only the reload count changes, so the saving = "
        "L·(dockings − receptors) is free of the receptor-size / sample-count differences between "
        "the campaigns.\nFixed per-run load  " + _Ls
        + ".  AutoDock loops Vina per ligand (receptor reloaded anyway) → no batch saving.",
        fontsize=11, y=0.995)
    _legend_between(fig, _batch_handles, _batch_labels, y=0.9, ncol=2, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.87))
    p = out_dir / "docking_effort_batch_vs_perrun.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    # ── console summary ──────────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("BATCH vs PER-RUN pose-production cost (load amortised across ligands sharing a receptor)")
    print("=" * 100)
    print("Fixed per-run load L:  " + "  ".join(
        f"{BASE_NAMES[k]} = {loads[k]:.1f}s [{load_src[k]}]" for k in METHODS if k in loads))
    for ds in order_ds:
        print(f"\n{titles[ds].splitlines()[0]} — {titles[ds].splitlines()[1]}")
        print(f"  {'tool':<26}{'dockings':>9}{'recep':>7}{'reuse':>8}"
              f"{'per-run h':>11}{'batch h':>9}{'saved h':>9}{'saved%':>8}")
        for r in records[ds]:
            print(f"  {r['name']:<26}{r['n_dockings']:>9}{r['n_receptors']:>7}"
                  f"{r['reuse_factor']:>7.0f}x{r['per_run_wall_h']:>11.2f}{r['batch_wall_h']:>9.2f}"
                  f"{r['saving_wall_h']:>9.2f}{r['saving_pct']:>7.0f}%")
    print(f"\n  CSV written to: {out_dir / 'batch_vs_perrun.csv'}")

    return [p] + _save_panels(panels, out_dir / "panels", plt, figsize=(6.6, 5.2))


def make_dataset_comparison_figure(base_args, out_dir: Path):
    """Empirical INDIVIDUAL- vs BATCH-docking timing, per tool (a cross-campaign figure).

    The benchmark campaign docks each (receptor, ligand) on its OWN native crystal — one fresh
    run per complex. The Orai × benchmark campaign batch-docks the same ligands into shared
    receptor frames (~300 ligands per frame), so a receptor's setup is amortised across its
    ligands. This reads BOTH campaigns' actual per-complex timing (via ``build_table``, each on
    its own cross-tool common timed set) and contrasts wall-clock cost per pose — the headline
    panel is seconds per PB-VALID pose.

    Caveat surfaced in the figure: the campaigns also differ in the RECEPTOR docked (small
    crystals vs large Orai frames), so this empirical gap mixes batching with receptor size;
    ``--batch-comparison`` isolates the pure reload saving by holding marginal compute fixed.
    Both campaigns are put on the same footing here (DiffDock=smina, EquiBind=unguided+gnina).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    ds_meta = [("benchmark",      "Benchmark (docked individually)"),
               ("orai_benchmark", "Orai × benchmark (batch-docked)")]
    summ_by_ds = {}
    for ds, _lab in ds_meta:
        a = _clone_dataset_args(base_args, ds)
        try:
            _pc, summ, _meta = build_table(a)
        except SystemExit as e:                                  # EquiBind filter etc.
            print(f"  [warn] dataset-comparison: {ds} build failed ({e}).")
            continue
        if summ is not None and not summ.empty:
            summ_by_ds[ds] = summ.set_index("method")
    present = [ds for ds, _ in ds_meta if ds in summ_by_ds]
    if len(present) < 2:
        print("  [warn] dataset-comparison needs BOTH campaigns timed; skipping.")
        return None
    b_ds, o_ds = "benchmark", "orai_benchmark"
    order = [k for k in METHODS if all(k in summ_by_ds[ds].index for ds in present)]
    if not order:
        print("  [warn] dataset-comparison: no tool timed in both campaigns; skipping.")
        return None
    names = [METHODS[k][0] for k in order]
    tool_colors = [METHODS[k][1] for k in order]
    x = np.arange(len(order)); w = 0.38

    # neutral swatches: solid = individual (benchmark), hatched = batch (Orai) — colour still = tool
    cmp_labels = ["individual (benchmark crystal)", "batch (Orai frame)"]
    cmp_handles = [Patch(facecolor=LEGEND_GREY, label=cmp_labels[0]),
                   Patch(facecolor=LEGEND_GREY, hatch=HATCH_SECONDARY, edgecolor="white",
                         label=cmp_labels[1])]

    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for k in order:
        rb, ro = summ_by_ds[b_ds].loc[k], summ_by_ds[o_ds].loc[k]
        rows.append({
            "method": k, "name": METHODS[k][0],
            "benchmark_complexes": int(rb["n_complexes"]), "orai_complexes": int(ro["n_complexes"]),
            "benchmark_s_per_generated": rb["wall_s_per_generated"],
            "orai_s_per_generated": ro["wall_s_per_generated"],
            "benchmark_s_per_pb_valid": rb["wall_s_per_pb_valid"],
            "orai_s_per_pb_valid": ro["wall_s_per_pb_valid"],
            "benchmark_total_wall_h": rb["total_wall_h"], "orai_total_wall_h": ro["total_wall_h"],
        })
    pd.DataFrame(rows).to_csv(out_dir / "dataset_timing_comparison.csv", index=False)

    metrics = [
        ("wall_s_per_pb_valid",  "Wall-clock seconds per PB-valid pose",  "{:.2f}"),
        ("wall_s_per_generated", "Wall-clock seconds per generated pose", "{:.2f}"),
        ("total_wall_h",         "Total campaign wall-clock (hours)",     "{:.2f}"),
    ]

    def _panel(col, ylabel, fmt):
        def draw(ax, standalone=False):
            bv = np.array([summ_by_ds[b_ds].loc[k, col] for k in order], float)
            ov = np.array([summ_by_ds[o_ds].loc[k, col] for k in order], float)
            ax.bar(x - w / 2, bv, w, color=tool_colors)
            ax.bar(x + w / 2, ov, w, color=tool_colors, hatch=HATCH_SECONDARY, edgecolor="white")
            allv = [v for v in np.r_[bv, ov] if v == v]
            ymax = (max(allv) if allv else 1) or 1
            for off, vals in ((-w / 2, bv), (w / 2, ov)):
                for xi, v in zip(x + off, vals):
                    if v == v:
                        ax.text(xi, v + 0.01 * ymax, fmt.format(v), ha="center",
                                va="bottom", fontsize=8)
            ax.set_ylim(0, ymax * 1.28)
            ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
            ax.set_ylabel(ylabel); ax.grid(axis="y", alpha=0.25)
            if standalone:
                ax.legend(handles=cmp_handles, loc="lower center", bbox_to_anchor=(0.5, 1.0),
                          ncol=2, fontsize=8, frameon=True)
        return draw

    panels = [(f"dataset_cmp_{col}", _panel(col, ylabel, fmt)) for col, ylabel, fmt in metrics]

    fig, ax = plt.subplots(1, 3, figsize=(18, 6.2))
    for a, (_n, draw) in zip(ax, panels):
        draw(a)
    _label_panels(ax)
    fig.suptitle(
        "Docking cost per pose: each complex docked INDIVIDUALLY (benchmark, native crystals) "
        "vs ligands BATCH-docked into shared Orai frames.\n"
        "DiffDock = smina, EquiBind = unguided+gnina in both campaigns; each on its own common "
        "timed set. The receptor also differs (crystal vs Orai frame), so this gap mixes batching "
        "with receptor size — --batch-comparison isolates the pure reload saving.",
        fontsize=12, y=0.995)
    _legend_between(fig, cmp_handles, cmp_labels, y=0.87, ncol=2)
    fig.tight_layout(rect=(0, 0, 1, 0.83))
    p = out_dir / "docking_effort_individual_vs_batch.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    # ── console summary (answers "per-pose cost for PB-valid poses, in seconds") ──
    print("\n" + "=" * 100)
    print("INDIVIDUAL (benchmark crystal) vs BATCH (Orai frame) docking cost per pose — seconds")
    print("=" * 100)
    print(f"  {'tool':<26}{'per generated pose':>24}{'per PB-valid pose':>24}")
    print(f"  {'':<26}{'individual / batch':>24}{'individual / batch':>24}")
    for k in order:
        rb, ro = summ_by_ds[b_ds].loc[k], summ_by_ds[o_ds].loc[k]
        print(f"  {METHODS[k][0]:<26}"
              f"{rb['wall_s_per_generated']:>11.2f} /{ro['wall_s_per_generated']:>10.2f}"
              f"{rb['wall_s_per_pb_valid']:>13.2f} /{ro['wall_s_per_pb_valid']:>9.2f}")
    print(f"\n  CSV written to: {out_dir / 'dataset_timing_comparison.csv'}")

    return [p] + _save_panels(panels, out_dir / "panels", plt, figsize=(6.6, 5.2))


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dataset", default="benchmark",
                    choices=tuple(DATASET_DEFAULTS),
                    help="Which docking campaign to compare. 'benchmark' = PoseBusters ligands "
                         "in their native crystals (per-complex timing, crystal RMSD tiers). "
                         "'orai_benchmark' = the same ligands docked into Orai frames "
                         "(aggregated per-tool timing logs, no crystal → all/PB-valid tiers only). "
                         "Path / refiner args left unset default per dataset (see DATASET_DEFAULTS).")
    # Path + refiner args default to None so _apply_dataset_defaults can fill them per --dataset;
    # an explicit flag still overrides. (autodock-prep/-cpu, eq-mode/-clamp, gnina-gpu flags are
    # shared across datasets, so they keep concrete defaults.)
    ap.add_argument("--per-pose-csv", default=None)
    ap.add_argument("--basis", default=DEFAULT_COST_BASIS, choices=COST_BASES,
                    help="Currency for every per-pose cost view (the by-quality panels, "
                         "their stats sidecar, and the paired wall-clock block of "
                         "effort_stats.json). 'elapsed' divides wall_s as each reader "
                         "built it, which is NOT one quantity across arms: AutoDock's is "
                         "vina elapsed plus its gnina sum, DiffDock's is measured, and "
                         "EquiBind's is best-config compute over its own worker count. "
                         "'charged' divides cpu_core_s / --cpu-threads + gpu_s for every "
                         "arm alike. It is a true device-occupancy time only when the "
                         "AutoDock gnina term is charged at its measured occupancy "
                         "(--autodock-gnina-accounting device-occupancy); it does not break "
                         "the rule that CPU-core-s and GPU-s are never summed. Neither term "
                         "is divided by --autodock-optimizer-workers, so 'charged' is "
                         "invariant to that flag while 'elapsed' is not; whether that "
                         "invariance is meaningful depends on the gnina accounting "
                         "(process-sum ignores the stage's concurrency, device-occupancy "
                         "measures it). Default "
                         f"'{DEFAULT_COST_BASIS}', under which every committed effort "
                         "directory reproduces unchanged.")
    ap.add_argument("--cpu-threads", type=int, default=DEFAULT_COST_CPU_THREADS,
                    help="Divisor turning CPU-core-seconds into an elapsed-equivalent for "
                         "--basis charged: the number of hardware threads a saturating "
                         f"processor stage occupies. Default {DEFAULT_COST_CPU_THREADS}, "
                         "matching the workstation. Ignored under --basis elapsed.")
    ap.add_argument("--rmsd-column", default=DEFAULT_RMSD_COLUMN,
                    help="Per-pose column scored for the near-native tier (RMSD ≤ 2 Å). "
                         f"Default '{DEFAULT_RMSD_COLUMN}' — the SAME endpoint the headline "
                         "accuracy tables use (posebusters_validity_report.py --rmsd-column), "
                         "so cost and accuracy stay comparable. 'pb_rmsd' reproduces the "
                         "legacy, systematically more permissive PoseBusters endpoint.")
    ap.add_argument("--kabsch-column", default=DEFAULT_KABSCH_COLUMN,
                    help="Per-pose column scored for the internal-form tier (≤ 1 Å). "
                         f"Default '{DEFAULT_KABSCH_COLUMN}', matching the headline report's "
                         "--form-ok-kabsch endpoint; 'pb_kabsch_rmsd' is the legacy choice.")
    ap.add_argument("--autodock-dir", default=None)
    ap.add_argument("--autodock-prep", default="meeko", choices=("meeko", "mgl_tools"),
                    help="Benchmark AutoDock prep subfolder (ignored for orai_benchmark, whose "
                         "timing is a single batch log).")
    ap.add_argument("--autodock-cpu", type=int, default=32,
                    help="CPU threads Vina used per complex (for CPU-core-seconds). "
                         "Benchmark and Orai runs both used 32.")
    ap.add_argument("--autodock-refine", default=None,
                    choices=("raw", "smina", "gnina", "gnina_refinement"),
                    help="Which OPTIMIZER's elapsed time is added to the Vina search wall: "
                         "raw = none, otherwise the matching optimization_log.csv rows. It also "
                         "supplies the per-pose method label unless --autodock-method overrides "
                         "it. Default: raw.")
    ap.add_argument("--autodock-method", default=None,
                    help="Exact per-pose 'method' label of the AutoDock arm to cost, e.g. "
                         "autodock_mgltools_exh128_gnina. Unset, the label is derived from "
                         "--autodock-refine, which only ever reaches the four Meeko arms "
                         "(autodock, autodock_{smina,gnina,gnina_refinement}) — pairing an "
                         "ADFRsuite tree with a derived Meeko label would join this arm's "
                         "timings to another arm's pose counts. Campaign-scoped: the "
                         "cross-dataset figures keep it on the campaign it was typed for and "
                         "clear it on the other (see _clone_dataset_args).")
    ap.add_argument("--autodock-opt-cpu", type=int, default=1,
                    help="CPU threads used by each AutoDock smina/CPU-gnina optimization "
                         "(for optimizer CPU-core-seconds).")
    ap.add_argument("--autodock-optimizer-workers", type=int, default=None,
                    help="Concurrency divisor for the AutoDock optimizer term: the arm's "
                         "optimize_workers. Σ per-pose elapsed equals wall clock only for a "
                         "serial optimizer; at optimize_workers=16 (the exh128 gnina arm) that "
                         "sum is compute time and overstates the wall by roughly N×, so the "
                         "wall term is divided by N. An approximation, not a measurement — real "
                         "speed-up is sublinear. Unset = today's serial accounting, so every "
                         "committed effort directory reproduces unchanged. Resource-seconds are "
                         "never divided. Benchmark-only; Orai timing takes a different path.")
    ap.add_argument("--autodock-gnina-gpu", action=argparse.BooleanOptionalAction,
                    default=False,
                    help="AutoDock gnina/gnina_refinement used the GPU. By default elapsed "
                         "time is classified as CPU work; set this for gnina_use_gpu=true.")
    ap.add_argument("--autodock-gnina-accounting", default=DEFAULT_GNINA_ACCOUNTING,
                    choices=GNINA_ACCOUNTINGS,
                    help="How the GPU-side AutoDock gnina stage is billed. 'process-sum' "
                         "(default) charges gpu_s the UNDIVIDED sum of the per-invocation "
                         "elapsed times in optimization_log.csv; on the exh128 arm those ran "
                         "optimize_workers=16 deep on one device, so the sum is ~14x the time "
                         "the device was busy, and every committed effort directory reproduces "
                         "unchanged under it. 'device-occupancy' charges the MEASURED union of "
                         "each complex's own optimizer intervals, read from the "
                         "optimized_<tool> provenance sidecars (created_at + "
                         "optimizer_elapsed_time_s), and additionally bills the stage's host "
                         "CPU at --autodock-gnina-host-cores instead of dropping it, and sets the "
                         "stage's wall_s to that union as well. Applies "
                         "only to gnina/gnina_refinement with --autodock-gnina-gpu; the CPU "
                         "branch already bills both. Benchmark-only; the Orai path stays "
                         "process-sum and prints a note.")
    ap.add_argument("--autodock-gnina-host-cores", type=float,
                    default=DEFAULT_GNINA_HOST_CORES,
                    help="CPU cores one concurrent GPU-gnina invocation occupies, for its "
                         f"CPU-core-seconds. Default {DEFAULT_GNINA_HOST_CORES}, the value the "
                         "exh128 gnina run config records as measured (16 x 2.1 = 33.6 against "
                         "32 hardware threads: the pass saturated the processor). Distinct "
                         "from --autodock-opt-cpu, the REQUESTED --cpu of a CPU-branch "
                         "optimization. Ignored under --autodock-gnina-accounting process-sum.")
    ap.add_argument("--unidock-dir", default=None,
                    help="Tiled Uni-Dock result root. A complex is timed only when its "
                         "summary, v2 JSON sentinel, and run manifest agree. Dataset default: "
                         "Dockings/unidock_results_full_protein_vina_scoring (benchmark).")
    ap.add_argument("--unidock2-dir", default=None,
                    help="Uni-Dock2 result root. A complex is timed only when its summary "
                         "matches the committed completion manifest. Dataset default: "
                         "Dockings/unidock2_results_full_protein_vina_scoring (benchmark).")
    ap.add_argument("--diffdock-dir", default=None)
    ap.add_argument("--diffdock-refine", default=None, choices=("raw", "smina", "gnina"),
                    help="DiffDock variant to time: raw docking, or docking + smina/gnina "
                         "post-hoc rescoring (adds optimization_log.csv time). Default: smina "
                         "(both datasets — the canonical oracle-collapse variant).")
    ap.add_argument("--diffdock-gnina-gpu", action=argparse.BooleanOptionalAction, default=True,
                    help="DiffDock ran gnina's CNN on the GPU (benchmark & orai: gnina_use_gpu true). "
                         "Use --no-diffdock-gnina-gpu if it was run with --no_gpu (CPU).")
    ap.add_argument("--equibind-dir", default=None)
    ap.add_argument("--eq-mode", default="unguided")
    ap.add_argument("--eq-refine", default=None, choices=("raw", "smina", "gnina"),
                    help="EquiBind refine variant. Default: gnina (benchmark, best PB validity) / "
                         "smina (orai_benchmark, the only refiner that run timed). smina is CPU.")
    ap.add_argument("--eq-gnina-gpu", action=argparse.BooleanOptionalAction, default=False,
                    help="EquiBind ran gnina's CNN on the GPU. Benchmark used --no_gpu "
                         "(gnina_use_gpu false, GPU busy with inference), so gnina refine is "
                         "CPU work by default; pass --eq-gnina-gpu if it was run on the GPU.")
    ap.add_argument("--eq-clamp", default=None,
                    help="EquiBind clamp variant (clampON/clampOFF). Leave unset for "
                         "clamp-off runs, whose poses carry no clamp token "
                         "(clamp_variant=None in pipeline_summary.json).")
    ap.add_argument("--ids-file", default=None,
                    help="Complex-id allowlist (one per line). Default: the benchmark id file / "
                         "no filter for orai_benchmark.")
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--no-plot", action="store_true")
    # ── empirical individual- vs batch-docking comparison (cross-dataset) ─────
    ap.add_argument("--dataset-comparison", action="store_true",
                    help="Also emit the empirical individual-vs-batch figure "
                         "(docking_effort_individual_vs_batch.png): reads BOTH the benchmark "
                         "(each complex docked on its own crystal) and Orai (ligands batch-docked "
                         "into shared frames) campaigns and contrasts wall-clock cost per pose — "
                         "headline seconds per PB-valid pose. Builds both datasets regardless of "
                         "--dataset. (Empirical: the receptor also differs, so it mixes batching "
                         "with receptor size; --batch-comparison isolates the pure reload saving.)")
    ap.add_argument("--dataset-comparison-out-dir", default=None,
                    help="Output dir for the individual-vs-batch figure/CSV "
                         "(default: posebusters_results/docking_effort_individual_vs_batch).")
    # ── batch-vs-per-run comparison (cross-dataset) ──────────────────────────
    ap.add_argument("--batch-comparison", action="store_true",
                    help="Also emit the cross-dataset batch-vs-per-run figure: models BOTH the "
                         "benchmark (per-run, unique crystal per ligand) and Orai (batch, shared "
                         "frames) campaigns to isolate the wall-clock batching saves by loading "
                         "each tool once per receptor instead of once per pair. Builds both "
                         "datasets regardless of --dataset.")
    ap.add_argument("--batch-out-dir", default=None,
                    help="Output dir for the batch-comparison figure/CSV "
                         "(default: posebusters_results/docking_effort_batch_vs_perrun).")
    ap.add_argument("--load-autodock", type=float, default=DEFAULT_LOAD_AUTODOCK,
                    help="Fixed per-run load L (s) for AutoDock in the batch model. Default 0: "
                         "Vina has no NN model and reloads the receptor per ligand even in a batch, "
                         "so batching saves it nothing at the docking level.")
    ap.add_argument("--load-diffdock", type=float, default=DEFAULT_LOAD_DIFFDOCK,
                    help="Fixed per-run load L (s) for DiffDock (ESM embeddings + model "
                         "checkpoints). Default 15 (estimate, under the ~18.5 s smallest observed "
                         "per-run wall).")
    ap.add_argument("--load-equibind", type=float, default=None,
                    help="Fixed per-run load L (s) for EquiBind. Default: measured from the "
                         "benchmark pipeline phase timing (~4.7 s); pass a value to override.")
    mf.add_method_filter_args(ap)
    args = ap.parse_args(argv)
    _set_gnina_accounting_mode(args)
    _record_explicit_dataset_args(args)      # before the fill, while None still means "unset"
    _apply_dataset_defaults(args)

    # Dataset-wide title + quality-tier ladder (crystal-free Orai loses the RMSD/Kabsch tiers).
    global DATASET_LABEL, QUALITY_TIERS
    DATASET_LABEL = DATASET_DEFAULTS[args.dataset]["title"]
    if args.dataset == "orai_benchmark":
        QUALITY_TIERS = ORAI_QUALITY_TIERS

    # Rewrite method labels + device tags to reflect the selected variants and gnina device.
    _ad_variant_label = {
        "gnina_refinement": "GNINA CNN-refinement",
    }.get(args.autodock_refine, args.autodock_refine)
    _ad_name = ("AutoDock Vina" if args.autodock_refine == "raw"
                else f"AutoDock Vina + {_ad_variant_label}")
    _ad_dev = ("CPU+GPU" if args.autodock_refine in _AUTODOCK_GNINA_REFINERS
               and args.autodock_gnina_gpu
               else "CPU")
    METHODS["autodock"] = (_ad_name, METHODS["autodock"][1], _ad_dev)
    if args.diffdock_refine == "raw" or (args.diffdock_refine == "gnina" and args.diffdock_gnina_gpu):
        _dd_dev = "GPU"                                   # docking (+ GPU gnina) only
    else:
        _dd_dev = "GPU+CPU"                               # docking GPU + CPU rescoring (smina / --no_gpu gnina)
    _dd_name = "DiffDock" if args.diffdock_refine == "raw" else f"DiffDock ({args.diffdock_refine})"
    METHODS["diffdock"] = (_dd_name, METHODS["diffdock"][1], _dd_dev)
    METHODS["equibind"] = (f"EquiBind ({args.eq_mode}+{args.eq_refine})",
                           METHODS["equibind"][1], "GPU+CPU")  # inference GPU + prep/post (+CPU gnina) CPU

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pc, summ, meta = build_table(args)
    if summ.empty:
        print("No data found — check the docking dirs and per-pose CSV path.")
        return 1

    pc.to_csv(out_dir / "per_complex_effort.csv", index=False)
    summ.to_csv(out_dir / "effort_summary.csv", index=False)

    # ── paired stats over the common timed set (degrades gracefully) ─────────
    stats = {}
    try:
        stats = compute_effort_stats(pc, summ, meta, basis=args.basis,
                                     cpu_threads=args.cpu_threads)
        (out_dir / "effort_stats.json").write_text(
            json.dumps(_jsonable(stats), indent=2))
        print(f"  Stats sidecar:       {out_dir / 'effort_stats.json'}")
    except Exception as e:
        print(f"  [warn] effort stats failed ({e!r}); figures fall back to test-free.")
        stats = {}

    # ── effort-by-quality-tier paired stats (degrades gracefully) ────────────
    q_stats = {}
    try:
        q_stats = compute_effort_by_quality_stats(pc, endpoint=_endpoint_cols(args),
                                                  basis=args.basis,
                                                  cpu_threads=args.cpu_threads)
        (out_dir / "effort_by_quality_stats.json").write_text(
            json.dumps(_jsonable(q_stats), indent=2))
        print(f"  Stats sidecar:       {out_dir / 'effort_by_quality_stats.json'}")
    except Exception as e:
        print(f"  [warn] effort-by-quality stats failed ({e!r}); figure falls back to test-free.")
        q_stats = {}

    # ── console report ───────────────────────────────────────────────────
    print("\n" + "=" * 100)
    _eq_cfg = f"{args.eq_mode}+{args.eq_refine}" + (f"+{args.eq_clamp}" if args.eq_clamp else "")
    print(f"DOCKING COMPUTATIONAL EFFORT — {DATASET_LABEL} "
          f"(EquiBind best config: {_eq_cfg})")
    print("=" * 100)
    _drops = ", ".join(f"{k} {v}" for k, v in meta.get("dropped_to_common", {}).items() if v)
    print(f"Compared on the COMMON set of {meta['common_n']} complexes timed for all methods"
          + (f" (dropped to reach it — {_drops})." if _drops else "."))
    print("WALL-CLOCK (time-to-result):")
    print(f"  {'method':<26}{'dev':>8}{'cplx':>6}{'total h':>9}{'gen':>9}"
          f"{'valid':>8}{'valid%':>8}{'s/gen':>8}{'s/valid':>9}")
    for _, r in summ.iterrows():
        print(f"  {r['name']:<26}{r['device']:>8}{r['n_complexes']:>6}{r['total_wall_h']:>9.2f}"
              f"{r['poses_generated']:>9,}{r['poses_pb_valid']:>8,}{r['pb_valid_rate_pct']:>7.1f}%"
              f"{r['wall_s_per_generated']:>8.2f}{r['wall_s_per_pb_valid']:>9.2f}")
    print("\nHARDWARE EFFORT (resource-seconds; CPU-core-s and GPU-s are different currencies):")
    print(f"  {'method':<26}{'CPU core-h':>12}{'GPU-h':>9}"
          f"{'CPUs/valid':>12}{'GPUs/valid':>12}")
    for _, r in summ.iterrows():
        print(f"  {r['name']:<26}{r['total_cpu_core_h']:>12.2f}{r['total_gpu_h']:>9.2f}"
              f"{r['cpu_core_s_per_pb_valid']:>12.2f}{r['gpu_s_per_pb_valid']:>12.2f}")
    print(f"\n  AutoDock docking CPU-core-seconds = Vina wall × {args.autodock_cpu} threads. "
          "DiffDock docking is GPU (CPU prep not separately recorded).")
    if args.autodock_refine != "raw":
        _ad_opt_dev = ("GPU" if args.autodock_refine in _AUTODOCK_GNINA_REFINERS
                       and args.autodock_gnina_gpu
                       else f"CPU × {args.autodock_opt_cpu} threads")
        _ad_opt_workers = getattr(args, "autodock_optimizer_workers", None)
        _ad_opt_conc = ("serial" if not _ad_opt_workers
                        else f"÷{int(_ad_opt_workers)} for optimize_workers concurrency")
        print(f"  AutoDock timing = Vina wall + {args.autodock_refine} optimization "
              f"(optimization_log.csv, {_ad_opt_conc}; {_ad_opt_dev}).")
    if "unidock" in set(summ["method"]):
        print("  Tiled Uni-Dock timing = committed summary elapsed_s; reported as GPU "
              "occupancy because planning/merge phase timing is not split out.")
    if "unidock2" in set(summ["method"]):
        print("  Uni-Dock2 timing = committed end-to-end summary elapsed_s. Its CPU prep / "
              "GPU docking split is unavailable, so resource-second fields are NaN.")
    if args.diffdock_refine != "raw":
        _dd_dev_note = ("GPU (CUDA CNN)" if (args.diffdock_refine == "gnina" and args.diffdock_gnina_gpu)
                        else "CPU")
        print(f"  DiffDock timing = docking wall + {args.diffdock_refine} rescoring "
              f"(optimization_log.csv, serial; {_dd_dev_note}).")
    if args.eq_refine == "gnina":
        print(f"  EquiBind gnina refine ran on {'GPU' if args.eq_gnina_gpu else 'CPU (--no_gpu)'} "
              "(gnina_use_gpu " + ("true" if args.eq_gnina_gpu else "false") + ").")
    eqfull = summ.loc[summ['method'] == 'equibind', 'total_wall_full_h'].values
    eqfull_str = f"{eqfull[0]:.2f}" if len(eqfull) else "n/a"
    if args.dataset == "orai_benchmark":
        print(f"  EquiBind wall = per-complex compute (unguided+{args.eq_refine}) ÷ worker pool; the "
              "Orai run had no fan-out, so full-run wall = best-config wall. Its globally-decoupled "
              "pipeline recorded fresh timing only for the frames re-docked this run (others were "
              "cached at ~0 s) — the common timed set is restricted to those frames.")
    else:
        print(f"  EquiBind wall = best-config (unguided+{args.eq_refine}) compute ÷ worker pool "
              f"(pure best-config-only estimate; full fan-out run wall = {eqfull_str} h).")
    print(f"\n  CSVs written to: {out_dir}/")

    if not args.no_plot:
        wall_files = make_wall_figure(pc, summ, out_dir, stats)
        res_files = make_resource_figure(summ, out_dir)
        time_fig = make_time_per_pose_figure(summ, out_dir)
        quality_files = make_effort_by_quality_figure(pc, out_dir, q_stats,
                                                      basis=args.basis,
                                                      cpu_threads=args.cpu_threads)
        print(f"  Figure (wall-clock): {wall_files[0]}")
        print(f"  Figure (resources):  {res_files[0]}")
        print(f"  Figure (time/pose):  {time_fig}")
        if quality_files:
            print(f"  Figure (by quality): {quality_files[0]}")
        panels = wall_files[1:] + res_files[1:] + (quality_files[1:] if quality_files else [])
        print(f"  Individual panels ({len(panels)}) in {out_dir / 'panels'}/:")
        for pp in panels:
            print(f"      - {pp.name}")

    # ── empirical individual-vs-batch figure (builds both campaigns itself) ───
    if args.dataset_comparison:
        cmp_out = Path(args.dataset_comparison_out_dir or
                       "posebusters_results/docking_effort_individual_vs_batch")
        try:
            cmp_files = make_dataset_comparison_figure(args, cmp_out)
            if cmp_files:
                print(f"\n  Figure (individual vs batch): {cmp_files[0]}")
                print(f"  Individual panels ({len(cmp_files) - 1}) in {cmp_out / 'panels'}/")
        except Exception as e:
            print(f"  [warn] individual-vs-batch comparison failed ({e!r}).")

    # ── cross-dataset batch-vs-per-run figure (builds both campaigns itself) ──
    if args.batch_comparison:
        batch_out = Path(args.batch_out_dir or
                         "posebusters_results/docking_effort_batch_vs_perrun")
        try:
            batch_files = make_batch_vs_perrun_figure(args, batch_out)
            if batch_files:
                print(f"\n  Figure (batch vs per-run): {batch_files[0]}")
                print(f"  Individual panels ({len(batch_files) - 1}) in {batch_out / 'panels'}/")
        except Exception as e:
            print(f"  [warn] batch-vs-per-run comparison failed ({e!r}).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
