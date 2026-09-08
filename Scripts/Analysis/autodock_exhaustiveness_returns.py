#!/usr/bin/env python3
"""Marginal return on compute for the RAW AutoDock Vina exhaustiveness ladder.

Answers one question: as exhaustiveness rises, how much extra docking success does
each additional hour of search buy, and where does that return stop being worth
paying for? Measured separately on each of the three gates and on their conjunction,
because they saturate at very different points:

  (a) RMSD <= 2 A        near-native PLACEMENT      column `rmsd`          (<= 2.0)
  (b) PoseBusters-valid  physical plausibility      column `pb_valid`
  (c) Kabsch < 1 A       internal conformer FORM    column `bestfit_rmsd`  (< 1.0, STRICT)
  (d) triple             (a) AND (b) AND (c) on ONE pose

Gate expressions are copied verbatim from the published cascade
(`aggregate_pose_validity_cascade`, posebusters_pose_comparison.py:4601-4613) and the
script ASSERTS that it reproduces that table before computing anything else. A drift
there means every number below it is unmoored from the thesis, so it is a hard error.

TWO STRANDS, ONE LADDER. The search-effort ladder is RAW arms only, and every cost and
marginal-return number here is computed on it alone. Where a rung also has a
gnina-rescored counterpart, that is carried as a SECOND STRAND at the same
exhaustiveness — never as a further rung — because it shares the raw arm's poses,
receptor and ligands and can only reorder what the search already found. Folding it
into the ladder would make the optimiser's gain look like extra sampling, which is the
easiest way to misread this data. `--no-optimizer-strand` drops it entirely.

    Search hours and re-ranking hours are NEVER combined. The first is 32-thread Vina
    wall time on CPU; the second is summed single-process optimiser time on GPU. They
    are different work on different devices and are not commensurable, so section 9
    reports the re-ranking cost separately and in COMPUTE hours — which, unlike wall,
    does not depend on how many optimiser workers a given pass happened to use.

WHAT COSTS WHAT
    The cost of an arm is the measured Vina-subprocess wall-clock second, summed over
    the same complexes the gates use, at a stated per-arm thread allocation. The
    allocation is part of the cost DEFINITION, not a normaliser: core-seconds
    (wall x cpus_per_worker) are a constant x32 rescale for every arm except exh32,
    where the x20 is the sole reason the first marginal cost comes out negative. No
    marginal return in this report uses them.

    Machine state is reported as a two-way arm effect on log wall time, fitted over the
    balanced complex x arm panel. That is symmetric across arms. A leave-one-arm-out
    "model" is not: it holds one arm out of its own fit and then calls the residual an
    inflation, which flags whichever arm you chose to exclude. This script does not do
    that, does not exclude any run-time window, and never substitutes a modelled cost
    for a measured one.

Usage (run from the repository root):

    /home/manndo/anaconda3/envs/vina/bin/python \\
        Scripts/Analysis/autodock_exhaustiveness_returns.py \\
        --per-pose-csv posebusters_results/benchmark_full_protein_vina_scoring/dock/pose_comparison_report/per_pose_metrics.csv \\
        --out-dir posebusters_results/autodock_exhaustiveness_returns

KNOWN TRAPS, all of them load-bearing
  * Depth changes the ANSWER, not just the precision. At rank-1 the ladder stops
    climbing after exh64; from depth 2 onward the top rung leads at every depth. Extra
    sampling finds poses Vina's own scoring function then fails to promote, so a
    rank-1-only curve and an oracle-only curve support opposite recommendations. Both
    are always reported.
  * The batch effect is the same size as the treatment effect and is perfectly aliased
    with it. Each arm is ONE batch, run on a different day, with no replication and no
    background-load telemetry. The residual arm-effect spread is ~19 pp, and the fitted
    per-rung wall step it has to be read against is not constant: large at the bottom of
    the ladder, but only ~15 % and ~14 % on the top two rungs. Where the knee is claimed,
    the batch effect is LARGER than the step it has to resolve. This is the dominant
    limitation of the cost side and the report prints it whether or not anything looks
    wrong.
  * exh32 is the `autodock_mgltools` arm, not a separate "raw" arm, and it is the one
    rung that ran at a different thread allocation (20 threads against 32). On wall
    time it interpolates within ~3 % of its neighbours, so it stays in the ladder; on
    core-seconds it is the only arm that can flip a sign, so core-seconds are excluded.
  * Arm-level `wall_clock_s` in exhaustiveness_arm_status.json covers 308 complexes, is
    absent for exh32, and is non-monotone. Cost is always re-summed per complex from
    docking_log_mgl_tools.csv over the analysed cohort.
  * `bestfit_rmsd` is a conformer deviation measured AFTER optimal superposition. Most
    poses clearing gate (c) are nowhere near the crystal site. Gate (c) is a diagnostic
    and is never plotted as a standalone success curve.
  * PoseBusters validity is saturated across this ladder and its residual gap is a
    hydrogen-handling artefact, not pose repair. Flatness is reported as a bounded
    equivalence, never as "no effect" read off a non-significant test.
  * Degeneracy is a property of a gate x DEPTH cell, never of a gate. PoseBusters
    validity currently has zero discordant complexes at top-3, top-10, top-15 and the
    oracle, but one at rank-1 and one at top-5, so those two cells are SATURATED and
    carry a live contrast rather than being excluded. Which cells are degenerate also
    moves when a rung is added, so it is computed per cell at run time and never
    inferred from a sampled depth.
  * The reciprocal (hours per additional complex) divides by the YIELD delta, which is
    a handful of complexes and is negative at rank-1 on the top rung. It is emitted
    only when the bootstrap numerator stays positive; otherwise the cell refuses.
  * The knee is inferentially claimable in ONE pre-registered cell only. Taking the
    largest rung whose interval clears zero across a family of cells places a knee at
    the top rung far more often than its nominal rate.

Statistical detail lives in the `*_stats.txt` sidecars; figures carry descriptive
content only, per the project convention.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1]
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import stats_utils as su  # noqa: E402
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)

try:
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                        # pragma: no cover
    def _label_panels(axes, fontsize: int = 13) -> None:
        import string
        flat = list(axes.ravel()) if hasattr(axes, "ravel") else (
            list(axes) if isinstance(axes, (list, tuple)) else [axes])
        i = 0
        for ax in flat:
            if ax is None or not ax.get_visible():
                continue
            ax.set_title(f"({string.ascii_uppercase[i % 26]})", loc="left",
                         fontweight="bold", fontsize=fontsize)
            i += 1

# ── Analysis axes ────────────────────────────────────────────────────────────────
GATES = ["rmsd2", "pb_valid", "kabsch1", "triple"]
GATE_LABEL = {
    "rmsd2": "Root-mean-square deviation ≤ 2 Å",
    "pb_valid": "PoseBusters-valid",
    "kabsch1": "Best-fit (Kabsch) deviation < 1 Å",
    "triple": "≤ 2 Å & PoseBusters-valid & Kabsch < 1 Å",
}
GATE_SHORT = {"rmsd2": "RMSD≤2Å", "pb_valid": "PB-valid",
              "kabsch1": "Kabsch<1Å", "triple": "triple"}
# Gate (c) alone is a conformer metric, not a placement metric, so it is never one of
# the plotted yield curves. Gate (b) is saturated. Both remain in every table.
HEADLINE_GATES = ["rmsd2", "triple"]

DEPTHS = [1, 3, 5, 10, 15, 30]
SUPPORT_DEPTHS = [3, 5, 10, 15]

# Colour encodes the exhaustiveness RUNG in every figure here, and nothing else. Gates
# are carried by marker shape, line style and hatch, so the two axes never share a
# channel.
#
# Thesis house style: hue encodes the TOOL, and an ordinal variable WITHIN one tool is a
# lightness ramp of that tool's hue. Every arm here is AutoDock, so the ramp is built
# around AutoDock blue #1f77b4 — two tints above it and two shades below, using the same
# tint(c, f) = c + (1 - c) * f helper as docking_effort_comparison.py, mirrored as
# shade(c, f) = c * (1 - f), at f = 0.42 and 0.22. The middle rung is literally the tool
# colour. Legend swatches that stand for the GATE rather than for a rung use LEGEND_GREY,
# per the house rule that a secondary split keeps its identity out of the colour channel.
#
# Light -> dark with search effort. Measured L* runs 69.6 / 59.1 / 48.0 / 37.7 / 27.8, so
# the adjacent steps are 10.4 / 11.2 / 10.3 / 9.9 and the ladder stays ordered in print
# greyscale. Separation is carried by luminance alone, which no dichromacy can collapse.
# exh32 is the lightest rung ever drawn as a line or a bar fill (exh18 appears only as an
# edged marker) and holds 3.26:1 on white, clearing the 3:1 line-legibility bar.
ARM_COLOR = {18: "#7db0d4", 32: "#5095c4", 64: "#1f77b4",
             92: "#185d8c", 128: "#124568"}
_RAMP_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "autodock_exh", [ARM_COLOR[e] for e in (18, 32, 64, 92, 128)])
_RAMP_EXH_SPAN = (16.0, 128.0)
# Neutral ink for marker edges and bar outlines. The fills now span L* 70 to 28, so no
# single separator suits both ends: marker edges take white (best on the dark rungs) and
# bar outlines keep #333333 (best on the light rungs).
EDGE_LIGHT = "white"
EDGE_DARK = "#333333"
# Neutral swatch for legend entries that denote the GATE, which applies at every rung and
# must therefore not borrow any one rung's fill.
LEGEND_GREY = "#9E9E9E"
# Reserved for the "this arm is not comparable" flag. Carried by ring weight and black
# ink rather than by a second hue, because under house style a hue would read as another
# tool — #d95f02 in particular is DiffDock's variant colour elsewhere in this thesis.
FLAG_COLOR = "#000000"
_OPTIMIZER_SUFFIXES = ("_gnina_refinement", "_smina", "_gnina")

# The two strands the report carries. RAW is the search-effort ladder and is the ONLY
# one the cost / marginal-return machinery runs on: a GPU re-ranking pass is a
# different kind of work on a different device, and folding it into a search-effort
# cost axis is exactly the conflation this script exists to prevent. The rescored
# strand is instead compared against raw at equal exhaustiveness.
RAW_STRAND = "raw"
OPTIMIZER_STRAND = "gnina"
OPTIMIZER_SUFFIX = "_gnina"
STRAND_LABEL = {RAW_STRAND: "raw Vina", OPTIMIZER_STRAND: "+ gnina re-rank"}

RMSD_OK_A = 2.0
FORM_OK_KABSCH_A = 1.0              # posebusters_pose_comparison.py:724
BOOT_SEED = 20260828

# The knee is a maximum over a family of intervals, so it is only claimable where the
# cell was fixed before looking. Everything else is descriptive.
KNEE_CELL = ("triple", 0)           # 0 == oracle depth

# The published cascade rows this script must reproduce, or it has drifted from the
# thesis. Keyed FIRST by the reference convention the per-pose table was scored under
# (`reference_convention` column of per_pose_metrics.csv; a table that predates the
# column is `instance`), THEN by exhaustiveness; values are the 15 cascade columns in
# table order (CASCADE_PIN_ORDER). The `instance` rows are the published ones. Under
# `nearest` (every RMSD-derived gate evaluated against the closest deposited copy of the
# ligand, plan v2 D1/D2) there is no published table yet, so the entry is None: the
# self-check CSV is still written but the regression abort is skipped with a WARNING
# unless pins are supplied with --cascade-pins-from (plan step 3.1a).
CASCADE_PINS = {
    "instance": {
        18: (303, 8986, 303, 8974, 99.9, 126, 160, 1.8, 209, 2774, 30.9, 98, 112, 70.0, 1.2),
        92: (303, 8984, 303, 8954, 99.7, 194, 290, 3.2, 221, 2828, 31.5, 159, 212, 73.1, 2.4),
    },
    "nearest": None,
}
DEFAULT_REFERENCE_CONVENTION = "instance"
# The rungs a pin set covers, and the hub method keys they are read from when pins are
# loaded from a pose_validity_cascade.csv (`<method_prefix>_exh<E>`).
PINNED_EXHAUSTIVENESS = (18, 92)
# Cascade columns in pin order. The hub's CSV spells the percentage columns `*_%`.
CASCADE_PIN_ORDER = [
    "n_complexes", "n_poses", "pb_valid_complexes", "pb_valid_poses",
    "pb_valid_pct", "rmsd2_complexes", "rmsd2_poses", "rmsd2_pct",
    "kabsch1_complexes", "kabsch1_poses", "kabsch1_pct",
    "triple_complexes", "triple_poses", "triple_of_rmsd2_pct",
    "triple_of_all_pct"]

DEFAULT_PB_CONFIG = ("Scripts/Docking/Posebusters/"
                     "posebusters_benchmark_full_protein_config.yaml")
# TRAP: this default is the full-protein tree, NOT the canonical matched_equibind tree
# the REGENERATE.md stage command passes explicitly. Never invoke the script bare.
DEFAULT_PER_POSE = ("posebusters_results/benchmark_full_protein_vina_scoring/dock/"
                    "pose_comparison_report/per_pose_metrics.csv")
DEFAULT_OUT = "posebusters_results/autodock_exhaustiveness_returns"

CORE_SECOND_WARNING = (
    "Core-seconds are reported for reference only. They are wall x cpus_per_worker, "
    "i.e. a constant rescale for every arm sharing the modal thread allocation and a "
    "different one for any arm that does not. On that basis the first marginal cost of "
    "this ladder is NEGATIVE, which is multiplier arithmetic, not a measurement. No "
    "marginal return in this report uses core-seconds.")


# ── Helpers ──────────────────────────────────────────────────────────────────────
def _to_bool(s: pd.Series) -> pd.Series:
    """The cascade's validity coercion (posebusters_pose_comparison.py:1035-1046).

    Reimplemented rather than imported so this script does not pull a 14k-line module
    in for one helper. NaN maps to False, so a missing measurement FAILS its gate
    instead of being dropped — the pose still counts in the denominator.
    """
    if s.dtype == bool:
        return s
    if pd.api.types.is_numeric_dtype(s):
        return s.fillna(0).astype(bool)
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0"))


def _depth_label(d: int) -> str:
    return "oracle (all poses)" if d == 0 else f"top-{d}" if d > 1 else "rank-1"


def _color(e: int) -> str:
    """Colour for one exhaustiveness rung.

    Pinned for the rungs the campaign actually ran, so a rung keeps its colour when the
    ladder is run as a subset — sampling by position in the present ladder would recolour
    exh64 the moment `--min-exhaustiveness` dropped a rung below it. An unpinned rung is
    interpolated on log2(exhaustiveness) across the same span, which keeps it both in
    ladder order and inside the legible part of the ramp.
    """
    if e in ARM_COLOR:
        return ARM_COLOR[e]
    lo, hi = _RAMP_EXH_SPAN
    t = (np.log2(max(e, 1)) - np.log2(lo)) / (np.log2(hi) - np.log2(lo))
    return mcolors.to_hex(_RAMP_CMAP(float(np.clip(t, 0.0, 1.0))))


def _fw(rows: list[list[str]], headers: list[str], aligns: str | None = None) -> str:
    """Fixed-width table, the house rendering for a console/report block."""
    cols = len(headers)
    w = [max([len(str(headers[i]))] + [len(str(r[i])) for r in rows]) for i in range(cols)] \
        if rows else [len(str(h)) for h in headers]
    aligns = aligns or (">" * cols)
    def _row(cells):
        return "  ".join(
            f"{str(cells[i]):{'>' if aligns[i] == '>' else '<'}{w[i]}}"
            for i in range(cols))
    out = [_row(headers), "-" * len(_row(headers))]
    out += [_row(r) for r in rows]
    return "\n".join(out)


# ── Arm discovery ────────────────────────────────────────────────────────────────
def _exh_from_method(mkey: str, cfg_exh: int | None) -> int | None:
    """Exhaustiveness a method key stands for. The config value wins; the key's own
    ``_exhNN`` suffix is the fallback, and the unsuffixed MGLTools arm is 32."""
    if cfg_exh is not None:
        return int(cfg_exh)
    m = re.search(r"_exh(\d+)$", mkey)
    if m:
        return int(m.group(1))
    return 32 if mkey == "autodock_mgltools" else None


def _find_arm_config(tree: Path, root: Path) -> Path | None:
    """Locate the docking config whose ``output_dir`` IS this arm's tree.

    Matching on the declared output_dir rather than on a filename convention means a
    renamed config is still found, and a config pointing somewhere else is never
    silently attributed to this arm.

    A tree can legitimately have MORE than one config: the post-dock optimiser writes
    into its own namespace inside the same tree, so the gnina arm's config declares the
    same output_dir as the raw arm's. This function serves raw arms only, so it prefers
    the config declaring ``optimization: none`` and falls back to the first match.
    Picking by sort order instead would work today only by the accident that '.' sorts
    before '_', and would report the optimiser's settings for the raw arm the moment a
    config were renamed.
    """
    want = str(tree).rstrip("/")
    matches = []
    for cand in sorted((root / "Scripts" / "Docking").glob("*.yaml")):
        try:
            cfg = yaml.safe_load(cand.read_text()) or {}
        except Exception:
            continue
        out = cfg.get("output_dir")
        if out and str(root / str(out).rstrip("/")) == want:
            matches.append((cand, str(cfg.get("optimization", "none")).strip().lower()))
    if not matches:
        return None
    for cand, opt in matches:
        if opt in ("none", "", "null", "false"):
            return cand
    return matches[0][0]


def discover_arms(per_pose: pd.DataFrame, pb_config: Path, root: Path,
                  method_prefix: str) -> list[dict]:
    """Every raw exhaustiveness arm, read out of the PoseBusters arm registry.

    The ladder comes from the config's ``docking_directories`` block, so a rung added
    there (exh128 already is) joins this analysis the moment its poses reach
    per_pose_metrics.csv — no code change and no hardcoded list.
    """
    cfg = yaml.safe_load(pb_config.read_text()) or {}
    trees = cfg.get("docking_directories", {}) or {}
    scored = set(per_pose["method"].astype(str).unique())

    arms: list[dict] = []
    for mkey, rel in sorted(trees.items()):
        if not mkey.startswith(method_prefix) or mkey.endswith(_OPTIMIZER_SUFFIXES):
            continue
        tree = root / str(rel)
        cfg_path = _find_arm_config(tree, root)
        arm_cfg = yaml.safe_load(cfg_path.read_text()) if cfg_path else {}
        arm_cfg = arm_cfg or {}
        exh = _exh_from_method(mkey, arm_cfg.get("exhaustiveness"))
        if exh is None:
            continue
        # An arm that ran an optimiser inline still has clean Vina search timing (the
        # two are serialised), but its METHOD KEY must still be the raw one.
        arms.append({
            "method_key": mkey, "exhaustiveness": int(exh), "tree": tree,
            "config": cfg_path,
            "cpus_per_worker": arm_cfg.get("cpus_per_worker"),
            "num_modes": arm_cfg.get("num_modes"),
            "energy_range": arm_cfg.get("energy_range"),
            "seed": arm_cfg.get("seed"),
            "scoring_function": arm_cfg.get("scoring_function"),
            "optimization": str(arm_cfg.get("optimization", "none")).strip().lower(),
            "has_poses": mkey in scored,
            "has_timing": bool(list(tree.glob("*/*/docking_log_*.csv"))),
        })
    return sorted(arms, key=lambda a: a["exhaustiveness"])


def discover_optimizer_arms(raw_arms: list[dict], per_pose: pd.DataFrame,
                            root: Path,
                            cohort: set[str] | None = None) -> list[dict]:
    """The post-dock re-ranked counterpart of each raw rung, where one was scored.

    These are NOT separate docking trees. The optimiser writes into a namespace inside
    the raw arm's own tree, so a rung's rescored variant appears in per_pose_metrics.csv
    as ``<raw method key>_gnina`` while sharing the raw arm's poses, receptor and
    ligands. Registering it as its own PoseBusters directory would make that tree be
    scanned twice and double-count every pose.

    Cost is read from the tree's optimization_log.csv and reported as COMPUTE hours (the
    sum of per-pose optimiser seconds), restricted to the analysed cohort. Wall clock is
    unusable because the passes ran at 1 to 16 workers.

    Compute is NOT comparable across that boundary either, and the report says so. On
    this campaign the resourcing changed with the rung: measured per-pose optimiser time
    is roughly 1.5x higher at 16 workers than at 1, while doubling exhaustiveness at
    fixed resourcing moves it only a few percent — and under heavy GPU sharing the
    per-complex spread collapses, so the number starts measuring the scheduler rather
    than the pose. Each row is therefore the cost of THAT pass as it was run, never a
    point on a cost-versus-exhaustiveness curve. The resourcing is printed beside it so
    the confound is visible rather than implied.
    """
    scored = set(per_pose["method"].astype(str).unique())
    out: list[dict] = []
    for raw in raw_arms:
        mkey = f"{raw['method_key']}{OPTIMIZER_SUFFIX}"
        if mkey not in scored:
            continue
        rows = 0
        compute_s = 0.0
        tools: set[str] = set()
        # Restricted to the analysed cohort. The trees hold 308 complexes while every
        # other number in this report is on the paired subset, and the extras are not
        # cheap ones — mixing the bases would overstate the pass by a few percent.
        for log in raw["tree"].glob("*/*/optimization_log.csv"):
            cid = log.parents[1].name
            if cohort is not None and cid not in cohort:
                continue
            try:
                with open(log, newline="") as fh:
                    for r in csv.DictReader(fh):
                        rows += 1
                        tools.add(str(r.get("tool", "")))
                        try:
                            compute_s += float(r.get("elapsed_time_s") or 0.0)
                        except (TypeError, ValueError):
                            pass
            except Exception:
                continue
        # The gnina config is a sibling of the raw one and declares the worker count
        # the pass actually used; without it, assume the historical serial default.
        workers = 1
        opt_cpu = None
        gpu_slots = None
        opt_cfg = None
        if raw["config"] is not None:
            cand = raw["config"].with_name(raw["config"].stem + "_gnina.yaml")
            if not cand.exists() and raw["config"].stem.endswith("_mgltools"):
                cand = raw["config"]           # the exh32 arm optimises in-place
            if cand.exists():
                opt_cfg = cand
                try:
                    c = yaml.safe_load(cand.read_text()) or {}
                    workers = max(1, int(c.get("optimize_workers", 1) or 1))
                    opt_cpu = c.get("optimize_cpu")
                    gpu_slots = c.get("gpu_max_concurrent", 1)
                except Exception:
                    pass
        out.append({**raw, "method_key": mkey, "strand": OPTIMIZER_STRAND,
                    "raw_method_key": raw["method_key"], "optimizer_config": opt_cfg,
                    "optimizer_rows": rows, "optimizer_tools": sorted(t for t in tools if t),
                    "optimizer_compute_h": compute_s / 3600.0,
                    "optimize_workers": workers, "optimize_cpu": opt_cpu,
                    "gpu_max_concurrent": gpu_slots,
                    "has_poses": True})
    return sorted(out, key=lambda a: a["exhaustiveness"])


def load_timing(tree: Path) -> pd.DataFrame:
    """Per-complex Vina-subprocess wall time for one arm.

    `elapsed_time_s` brackets the Vina subprocess only (run_autodock.py:1556-1592):
    preparation, pose validation and any optimiser fall outside it. The complex id is
    the directory name, which is exactly the `protein` value in per_pose_metrics.csv —
    joining on protein + '__' + ligand instead matches nothing at all.
    """
    rows = []
    for log in sorted(tree.glob("*/*/docking_log_*.csv")):
        cid = log.parents[1].name
        try:
            df = pd.read_csv(log)
        except Exception:
            continue
        for _, r in df.iterrows():
            rows.append({
                "protein": cid,
                "elapsed_time_s": pd.to_numeric(r.get("elapsed_time_s"), errors="coerce"),
                "status": str(r.get("status", "")),
                "timestamp": r.get("timestamp"),
                "log_exhaustiveness": pd.to_numeric(r.get("exhaustiveness"), errors="coerce"),
                "log_num_modes": pd.to_numeric(r.get("num_modes"), errors="coerce"),
                "cpu_model": r.get("cpu_model"),
                "lig_heavy_atoms": pd.to_numeric(r.get("lig_heavy_atoms"), errors="coerce"),
                "lig_rotatable_bonds": pd.to_numeric(r.get("lig_rotatable_bonds"), errors="coerce"),
                "prot_num_atoms": pd.to_numeric(r.get("prot_num_atoms"), errors="coerce"),
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values("timestamp").drop_duplicates("protein", keep="last")
        out["run_order"] = out["timestamp"].rank(method="first").astype(int)
    return out


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_input_parity(arms: list[dict], complexes: list[str], full: bool) -> dict:
    """Assert the ladder is a single-variable comparison.

    Exhaustiveness is only the independent variable if every arm docked byte-identical
    ligands into byte-identical receptors inside byte-identical boxes. Ligand and box
    are always hashed; receptors only under ``full``, because they are ~1 MB each.
    """
    kinds = {"ligand": "_staging/ligands/pdbqt/*.pdbqt",
             "box": "_staging/receptors/pdbqt/*.box.txt"}
    if full:
        kinds["receptor"] = "_staging/receptors/pdbqt/*.pdbqt"
    report = {}
    for kind, pat in kinds.items():
        mismatched, partial, checked = [], [], 0
        for cid in complexes:
            digests, n_have = set(), 0
            for arm in arms:
                files = sorted((arm["tree"] / cid).glob(pat))
                if kind == "receptor":
                    files = [f for f in files if not f.name.endswith(".box.txt")]
                if files:
                    n_have += 1
                    digests.add("|".join(_sha(f) for f in files))
            if not digests:
                continue
            checked += 1
            if len(digests) > 1:
                mismatched.append(cid)
            elif n_have < len(arms):
                # One arm agreeing with itself is not evidence that the arms agree.
                # Certifying parity from a single staged copy would be the strongest
                # claim in this report resting on the weakest evidence.
                partial.append(cid)
        report[kind] = {"checked": checked, "n_arms": len(arms),
                        "n_mismatched": len(mismatched), "mismatched": mismatched[:20],
                        "n_partial": len(partial), "partial": partial[:20],
                        "identical": len(mismatched) == 0 and not partial}
    return report


# ── Gates ────────────────────────────────────────────────────────────────────────
def apply_gates(df: pd.DataFrame) -> pd.DataFrame:
    """Attach the four per-POSE gate booleans, using the published cascade's columns
    and comparison operators (posebusters_pose_comparison.py:4601-4613)."""
    out = df.copy()
    out["rmsd"] = pd.to_numeric(out["rmsd"], errors="coerce")
    kcol = "bestfit_rmsd" if "bestfit_rmsd" in out.columns else "pb_kabsch_rmsd"
    out["_kabsch"] = pd.to_numeric(out[kcol], errors="coerce")
    out["pb_valid"] = _to_bool(out["pb_valid"])
    out["rmsd2"] = out["rmsd"] <= RMSD_OK_A
    out["kabsch1"] = out["_kabsch"] < FORM_OK_KABSCH_A       # STRICT, as in the cascade
    # The triple gate lives on ONE pose. Gating the pose pool three times separately
    # and AND-ing the answers would count three different poses as a single success.
    out["triple"] = out["rmsd2"] & out["pb_valid"] & out["kabsch1"]
    return out


def audit_gate_inputs(df: pd.DataFrame) -> list[dict]:
    """Invariants that must hold for the gates to mean what the report says.

    These are warnings, not gates: each names a way the input could silently stop
    matching the published cascade's assumptions.
    """
    checks = []

    def add(name, ok, detail):
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    n_exploded = int((df["rmsd"] >= 1000).sum())
    add("no exploded-pose RMSD sentinels", n_exploded == 0,
        f"{n_exploded} poses with rmsd >= 1000 (max rmsd {df['rmsd'].max():.2f} A). "
        "The published cascade carries no exploded-pose guard; inside the triple gate "
        "such a guard is unreachable anyway, since rmsd <= 2 implies rmsd < 1000.")
    n_edge = int((df["_kabsch"] == FORM_OK_KABSCH_A).sum())
    add("strict-vs-nonstrict Kabsch is a no-op", n_edge == 0,
        f"{n_edge} poses sit exactly on the {FORM_OK_KABSCH_A} A boundary")
    if "pb_kabsch_rmsd" in df.columns:
        gap = float(np.nanmax(np.abs(pd.to_numeric(df["pb_kabsch_rmsd"], errors="coerce")
                                     - df["_kabsch"])))
        add("bestfit_rmsd == pb_kabsch_rmsd", gap == 0.0, f"max |difference| = {gap}")
    nan_counts = {c: int(df[c].isna().sum()) for c in ("rmsd", "_kabsch")}
    add("no NaN in the gate columns", all(v == 0 for v in nan_counts.values()),
        f"{nan_counts}; a NaN FAILS its gate rather than dropping the pose")
    # These two are RAW-strand invariants. A rescored arm ranks on the optimiser's CNN
    # affinity by design, so `rank` there is optimized_rank and Vina's affinity ordering
    # is deliberately broken — running these checks over both strands would report the
    # re-ranking working as a defect.
    if "autodock_rank" in df.columns:
        ar = pd.to_numeric(df["autodock_rank"], errors="coerce")
        frac = float((ar == df["rank"]).mean())
        add("rank is the raw Vina rank", frac == 1.0,
            f"rank == autodock_rank in {100 * frac:.2f}% of poses")
    if "optimized_rank" in df.columns:
        frac = float(pd.to_numeric(df["optimized_rank"], errors="coerce").isna().mean())
        add("no optimiser ranking present", frac == 1.0,
            f"optimized_rank is NaN in {100 * frac:.2f}% of poses")
    if "autodock_affinity" in df.columns:
        aff = pd.to_numeric(df["autodock_affinity"], errors="coerce")
        tmp = df.assign(_aff=aff).sort_values(["method", "protein", "rank"])
        bad = int(tmp.groupby(["method", "protein"])["_aff"]
                  .apply(lambda s: bool((s.diff().dropna() < -1e-9).any())).sum())
        add("rank follows Vina's affinity ordering", bad == 0,
            f"{bad} complex-arm pairs where affinity improves as rank worsens")
    dup = int(df.duplicated(subset=["method", "protein", "ligand", "rank"]).sum())
    add("no duplicate (arm, complex, rank) rows", dup == 0, f"{dup} duplicates")
    return checks


def cascade_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Reproduce the published pose-validity cascade row for every arm."""
    rows = []
    for exh, sub in df.groupby("exhaustiveness"):
        cplx = sub["protein"].astype(str) + "\x00" + sub["ligand"].astype(str)
        n = len(sub)
        rec = {"exhaustiveness": exh, "n_complexes": cplx.nunique(), "n_poses": n}
        for gate in ("pb_valid", "rmsd2", "kabsch1", "triple"):
            m = sub[gate].values
            rec[f"{gate}_complexes"] = int(cplx[m].nunique())
            rec[f"{gate}_poses"] = int(m.sum())
            rec[f"{gate}_pct"] = 100.0 * m.sum() / max(1, n)
        rm = int(sub["rmsd2"].sum())
        rec["triple_of_rmsd2_pct"] = (100.0 * rec["triple_poses"] / rm
                                      if rm else float("nan"))
        rec["triple_of_all_pct"] = 100.0 * rec["triple_poses"] / max(1, n)
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("exhaustiveness").reset_index(drop=True)


def detect_reference_convention(per_pose: pd.DataFrame) -> tuple[str, bool]:
    """(convention, declared): the reference convention the per-pose table was scored
    under, and whether the table itself declares it.

    A table written by the hub under --reference-convention carries one constant
    `reference_convention` column; a table that predates that column is by definition
    the single-instance convention. A mixed column is refused: the gates would then
    compare poses scored against different references inside one ladder.
    """
    if "reference_convention" not in per_pose.columns:
        return DEFAULT_REFERENCE_CONVENTION, False
    vals = sorted(set(per_pose["reference_convention"].dropna().astype(str).str.strip()))
    if len(vals) != 1:
        raise SystemExit(f"per_pose_metrics.csv carries a mixed reference_convention "
                         f"column {vals}; the ladder needs one convention throughout")
    conv = vals[0]
    if conv not in CASCADE_PINS:
        raise SystemExit(f"unknown reference_convention {conv!r}; this script knows "
                         f"{sorted(CASCADE_PINS)}")
    return conv, True


def load_cascade_pins(path: Path, method_prefix: str) -> dict[int, tuple]:
    """Read the pinned rungs out of a hub `pose_validity_cascade.csv`.

    Rows are matched on `method_key` == `<method_prefix>_exh<E>` for every rung in
    PINNED_EXHAUSTIVENESS, and the values are rounded exactly as assert_cascade rounds
    (percentages to one decimal, counts to int), so a pin set read from the canonical
    hub table reproduces CASCADE_PINS["instance"] to the digit.
    """
    if not path.is_file():
        raise SystemExit(f"--cascade-pins-from {path}: file not found. Under the nearest "
                         "convention the pins come from the hub's rebuilt "
                         "pose_validity_cascade.csv (plan 3.1a); omit the flag to run "
                         "unpinned with a WARNING.")
    tbl = pd.read_csv(path)
    if "method_key" not in tbl.columns:
        raise SystemExit(f"--cascade-pins-from {path}: no method_key column; expected a "
                         "hub pose_validity_cascade.csv")
    hub_cols = [c.replace("_pct", "_%") for c in CASCADE_PIN_ORDER]
    missing = [c for c in hub_cols if c not in tbl.columns]
    if missing:
        raise SystemExit(f"--cascade-pins-from {path}: missing columns {missing}")
    pins = {}
    for e in PINNED_EXHAUSTIVENESS:
        key = f"{method_prefix}_exh{e}"
        row = tbl[tbl["method_key"].astype(str) == key]
        if row.empty:
            raise SystemExit(f"--cascade-pins-from {path}: no row with method_key "
                             f"{key!r}; the pinned rungs are {PINNED_EXHAUSTIVENESS}")
        row = row.iloc[0]
        pins[e] = tuple(round(float(row[c]), 1) if c.endswith("_%") else int(row[c])
                        for c in hub_cols)
    return pins


def assert_cascade(casc: pd.DataFrame, pins_by_exh: dict[int, tuple]) -> list[str]:
    """Hard regression check against the pinned cascade table.

    Percentages are rounded once, at comparison time, exactly as the cascade rounds
    once at display time — comparing rounded intermediates can drift a digit.
    """
    problems = []
    order = CASCADE_PIN_ORDER
    for exh, pins in pins_by_exh.items():
        row = casc[casc.exhaustiveness == exh]
        if row.empty:
            continue
        row = row.iloc[0]
        for col, want in zip(order, pins):
            got = row[col]
            got = round(float(got), 1) if col.endswith("_pct") else int(got)
            if got != want:
                problems.append(f"exh{exh} {col}: got {got}, published {want}")
    return problems


def complex_matrix(df: pd.DataFrame, gate: str, depth: int, arms_exh: list[int],
                   complexes: list[str]) -> pd.DataFrame:
    """n_complexes x n_arms 0/1 frame: did this arm land at least one qualifying pose
    within the first `depth` Vina-ranked poses (depth 0 == every pose it produced)?"""
    sub = df if depth == 0 else df[df["rank"] <= depth]
    piv = (sub.groupby(["protein", "exhaustiveness"])[gate].any()
           .unstack("exhaustiveness"))
    return piv.reindex(index=complexes, columns=arms_exh).fillna(False).astype(bool)


# ── Cost ─────────────────────────────────────────────────────────────────────────
def two_way_arm_effects(cost: pd.DataFrame, arms_exh: list[int]) -> dict:
    """Machine-state diagnostic: residual arm effects on log wall time.

    Fits log t_ia = mu + alpha_i(complex) + beta_a(arm) by complex-demeaning, then
    removes the fitted exhaustiveness trend from beta. What is left is the part of an
    arm's cost that its exhaustiveness does not explain — i.e. the day it ran on.

    This is symmetric across arms by construction. A leave-one-arm-out fit is not: the
    held-out arm is the only one not constrained by its own residual, so whichever arm
    you exclude comes back looking inflated. That procedure is deliberately not used.
    """
    cols = [e for e in arms_exh if e in cost.columns]
    sub = cost[cols].dropna()
    if len(sub) < 10 or len(cols) < 3:
        return {"ok": False, "reason": "need >= 3 arms and >= 10 shared complexes"}
    L = np.log(sub.values)
    L = L - L.mean(axis=1, keepdims=True)          # remove the complex effect
    beta = L.mean(axis=0)                          # raw arm effects, log units
    x = np.log(np.asarray(cols, float))
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, beta, rcond=None)
    resid = beta - A @ coef                        # arm effect net of the exh trend
    return {"ok": True, "arms": cols,
            "arm_effect_pct": {int(e): float(100 * (np.exp(r) - 1))
                               for e, r in zip(cols, resid)},
            "spread_pp": float(100 * (np.exp(resid.max()) - np.exp(resid.min()))),
            "log_log_slope": float(coef[1]),
            "n_complexes": int(len(sub))}


def build_cost_table(arms: list[dict], timing: dict[int, pd.DataFrame],
                     complexes: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Per-complex and per-arm cost on the analysed cohort, in every currency."""
    cost = pd.DataFrame(index=complexes)
    for arm in arms:
        t = timing.get(arm["exhaustiveness"])
        if t is None or t.empty:
            continue
        cost[arm["exhaustiveness"]] = (t.set_index("protein")["elapsed_time_s"]
                                       .reindex(complexes))
    cost = cost.dropna(axis=1, how="all")

    cpus = {a["exhaustiveness"]: a["cpus_per_worker"] for a in arms}
    vals = [cpus.get(e) for e in cost.columns if cpus.get(e) is not None]
    modal = int(pd.Series(vals).mode().iloc[0]) if vals else None
    # Machine state is estimated on the arms that share a thread allocation. An arm
    # that ran on a different number of threads carries a systematic penalty, not a
    # machine-state effect, and folding it in would bend the exhaustiveness trend the
    # residuals are measured against.
    modal_arms = [a["exhaustiveness"] for a in arms
                  if modal is None or cpus.get(a["exhaustiveness"]) == modal]
    effects = two_way_arm_effects(cost, modal_arms)

    rows = []
    prev_e = None
    for arm in arms:
        e = arm["exhaustiveness"]
        if e not in cost.columns:
            rows.append({"exhaustiveness": e, "method_key": arm["method_key"],
                         "timing_available": False})
            prev_e = e
            continue
        col = cost[e].dropna()
        s = col.sort_values()
        k = int(0.1 * s.size)
        ratio = np.nan
        if prev_e is not None and prev_e in cost.columns:
            pair = cost[[prev_e, e]].dropna()
            ratio = float((pair[e] / pair[prev_e]).median())
        rows.append({
            "exhaustiveness": e, "method_key": arm["method_key"],
            "timing_available": True, "n_complexes": int(col.size),
            "cpus_per_worker": cpus.get(e),
            "total_wall_s": float(col.sum()),
            "total_wall_h": float(col.sum()) / 3600.0,
            "mean_wall_s": float(col.mean()),
            "median_wall_s": float(col.median()),
            "trimmed_mean_10pct_wall_s": float(s.iloc[k:s.size - k].mean()),
            "iqr_lo_s": float(col.quantile(0.25)),
            "iqr_hi_s": float(col.quantile(0.75)),
            "paired_median_ratio_to_prev": ratio,
            # Diagnostic only. Never an axis, never a denominator — see the warning.
            "total_core_h": (float(col.sum()) * cpus[e] / 3600.0
                             if cpus.get(e) else np.nan),
            "arm_effect_pct": effects.get("arm_effect_pct", {}).get(e, np.nan),
            "thread_allocation_modal": (modal is None or cpus.get(e) == modal),
            "cost_comparable_wall": True,
            "cost_comparable_core": (modal is not None and cpus.get(e) == modal),
        })
        prev_e = e
    return cost, pd.DataFrame(rows), {"modal_cpus": modal, "arm_effects": effects}


# ── Marginal return ──────────────────────────────────────────────────────────────
# A ratio needs a denominator that cannot change sign. These are the two thresholds
# that decide whether a cost delta is a cost at all. On this ladder they separate
# cleanly: the usable segments put 0.00-0.01 % of resamples at or below zero, the
# inverted one puts 76 %.
SIGN_CHANGE_MAX = 0.005          # refuse a ratio if the denominator flips this often
COST_FLOOR_H = 0.0               # a cost delta must be positive to be divided by


def _paired_index_draws(n: int, n_boot: int, seed: int) -> np.ndarray:
    """One shared set of complex-level resamples.

    Every quantity that must stay paired — the yield delta, the cost delta, and both
    segments of the diminishing-returns contrast — is evaluated on these same draws.
    """
    return np.random.default_rng(seed).integers(0, n, size=(n_boot, n))


def _nan_aware_ci(boots: np.ndarray, alpha: float) -> tuple[float, float, float]:
    """Percentile interval over the finite draws.

    A handful of degenerate draws (a resample whose cost delta collapses) must not
    poison the whole interval, which is what a plain np.percentile over an array
    containing NaN would do — it returns NaN and the cell silently reports nothing.
    Whether that handful is tolerable is decided by the sign-change guard, not here.
    """
    finite = boots[np.isfinite(boots)]
    if finite.size == 0:
        return np.nan, np.nan, np.nan
    lo, hi = np.percentile(finite, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(np.median(finite)), float(lo), float(hi)


def marginal_return(y_lo, y_hi, c_lo, c_hi, n_boot: int, seed: int,
                    alpha: float = 0.05) -> dict:
    """Net complexes recovered per extra wall-hour, with a paired complex bootstrap.

    Numerator and denominator are recomputed on the SAME resampled complexes. Holding
    the cost fixed would be wrong twice: it is measured per complex, and it correlates
    with difficulty, because the complexes extra sampling rescues are the slow ones.

    dY is a NET difference — at rank-1 the top segment has gains and losses that partly
    cancel, and reporting only the gains would overstate it.

    BOTH ends can fail. The denominator fails when a rung's measured cost is not above
    the rung below it, which happens on this ladder: the arms are single unreplicated
    batches run on different days, and the machine-state spread between days is as
    large as the cost step between rungs. A rate per hour computed against a
    non-positive or near-zero cost delta is not a measurement, so the cell refuses.
    """
    dy = np.asarray(y_hi, float) - np.asarray(y_lo, float)
    dc = (np.asarray(c_hi, float) - np.asarray(c_lo, float)) / 3600.0
    n = len(dy)
    draws = _paired_index_draws(n, n_boot, seed)
    num_b, den_b = dy[draws].sum(axis=1), dc[draws].sum(axis=1)
    # Negated comparison on purpose: `den_b <= floor` scores a NaN draw as False, i.e.
    # as GOOD, so a complex with no timing row would sail past both guards and publish
    # an interval computed from the handful of draws that happened to miss it.
    p_den_bad = float((~(den_b > COST_FLOOR_H)).mean())
    total_dc = float(dc.sum())

    if not total_dc > COST_FLOOR_H:
        est = lo = hi = np.nan
        status = (f"REFUSED: cost delta not positive ({total_dc:+.3f} wall-h) — the "
                  "rung above did not measurably cost more than the rung below")
    elif p_den_bad > SIGN_CHANGE_MAX:
        est = lo = hi = np.nan
        status = (f"REFUSED: cost delta changes sign in {100 * p_den_bad:.1f}% of "
                  "resamples")
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(den_b > COST_FLOOR_H, num_b / den_b, np.nan)
        est, lo, hi = _nan_aware_ci(ratio, alpha)
        status = "ok"

    gained = int((np.asarray(y_hi, bool) & ~np.asarray(y_lo, bool)).sum())
    lost = int((~np.asarray(y_hi, bool) & np.asarray(y_lo, bool)).sum())
    return {"delta_yield": float(dy.sum()), "gained": gained, "lost": lost,
            "delta_cost_h": total_dc,
            "rate_per_h": (float(dy.sum() / total_dc)
                           if total_dc > COST_FLOOR_H else np.nan),
            "rate_boot": est, "rate_lo": lo, "rate_hi": hi, "rate_status": status,
            "p_numerator_nonpositive": float((num_b <= 0).mean()),
            "p_denominator_nonpositive": p_den_bad}


def reciprocal_or_refusal(mr: dict, guard: float = 0.05) -> tuple[float, str]:
    """Wall-hours per additional complex recovered — or a typed refusal.

    Guarded at BOTH ends. The denominator here is the YIELD delta: a handful of
    complexes out of a few hundred, and negative at rank-1 on the top rung. A
    percentile interval for a ratio whose denominator can cross zero is not an
    interval, so the cell refuses rather than printing a confident-looking number.
    """
    if mr["rate_status"] != "ok":
        return float("nan"), mr["rate_status"]
    if mr["delta_yield"] <= 0:
        return float("nan"), "REFUSED: net yield delta is not positive"
    if mr["p_numerator_nonpositive"] > guard:
        return float("nan"), (f"REFUSED: bootstrap yield delta <= 0 in "
                              f"{100 * mr['p_numerator_nonpositive']:.1f}% of resamples")
    return mr["delta_cost_h"] / mr["delta_yield"], "ok"


def diminishing_contrast(seg_first: dict, seg_last: dict, n_boot: int, seed: int,
                         alpha: float = 0.05) -> dict:
    """Is the LAST segment's return per hour lower than the FIRST segment's?

    This is the diminishing-returns test. It targets the economic claim directly, needs
    no functional form, and works on a short ladder — where a fitted saturation curve
    would have one or two residual degrees of freedom and a ceiling estimate spanning
    most of the feasible range.

    Both segments are recomputed on the SAME resampled complexes, so complex-level
    pairing always holds. When the two chosen segments happen to be ADJACENT they also
    share a rung, which enters the two rates with opposite signs and makes the contrast
    conservative. When they are disjoint — which is the usual case once three or more
    segments have a usable cost — that term is absent and the interval is effectively
    the independent-segments interval. The report states which case it is in.
    """
    n = len(seg_first["dy"])
    draws = _paired_index_draws(n, n_boot, seed)
    da = seg_first["dc"][draws].sum(axis=1) / 3600.0
    db = seg_last["dc"][draws].sum(axis=1) / 3600.0
    ok = (da > COST_FLOOR_H) & (db > COST_FLOOR_H)
    with np.errstate(divide="ignore", invalid="ignore"):
        boots = np.where(ok,
                         seg_first["dy"][draws].sum(axis=1) / da
                         - seg_last["dy"][draws].sum(axis=1) / db,
                         np.nan)
    p_bad = float((~ok).mean())
    if p_bad > SIGN_CHANGE_MAX:
        return {"contrast": np.nan, "lo": np.nan, "hi": np.nan, "diminishing": False,
                "status": f"REFUSED: a segment's cost delta changes sign in "
                          f"{100 * p_bad:.1f}% of resamples"}
    est, lo, hi = _nan_aware_ci(boots, alpha)
    return {"contrast": est, "lo": lo, "hi": hi,
            "diminishing": bool(lo == lo and lo > 0), "status": "ok"}


# ── Saturation ───────────────────────────────────────────────────────────────────
def saturation_cell(mat: pd.DataFrame, arms_exh: list[int],
                    margins=(0.05, 0.025)) -> dict:
    """Classify one gate x depth cell, and bound it when it is flat.

    DEGENERATE: not one complex anywhere on the ladder changes its answer, so there is
    no test and the honest output is a bound. SATURATED: it moves, but by less than
    this design can resolve. INFORMATIVE: the omnibus rejects.

    The minimum detectable difference is computed from the LARGEST discordance on the
    ladder, not from each contrast's own. Using a contrast's own discordance is
    circular — discordance is small exactly where the rung adds nothing, so the
    diagnostic would certify the design as most sensitive precisely where it is blind.
    """
    n = len(mat)
    rates = {int(e): 100.0 * float(mat[e].mean()) for e in arms_exh}
    pairs = list(zip(arms_exh, arms_exh[1:]))
    disc = {}
    for lo_e, hi_e in pairs:
        a, b = mat[lo_e].values, mat[hi_e].values
        disc[(lo_e, hi_e)] = int(((a & ~b) | (~a & b)).sum())
    total_disc = sum(disc.values())

    q_stat = q_p = q_df = np.nan
    if total_disc > 0:
        q_stat, q_p, q_df = su.cochran_q(mat[arms_exh].values.astype(int))

    mde_pp = np.nan
    if total_disc > 0:
        lo_e, hi_e = max(disc, key=disc.get)
        try:
            pw = su.mcnemar_power(mat[lo_e].values, mat[hi_e].values)
            mde_pp = 100.0 * float(pw.get("mde", np.nan))
        except Exception:
            mde_pp = np.nan

    obs_range = max(rates.values()) - min(rates.values())
    if total_disc == 0:
        flag = "DEGENERATE"
    elif mde_pp == mde_pp and obs_range < mde_pp:
        flag = "SATURATED"
    elif q_p == q_p and q_p < 0.05:
        flag = "INFORMATIVE"
    else:
        flag = "SATURATED"

    equiv = {}
    if flag != "INFORMATIVE":
        a, b = mat[arms_exh[0]].values, mat[arms_exh[-1]].values
        for m in margins:
            t = su.tost_paired_proportions(a, b, margin=m)
            equiv[m] = {"diff_pp": 100 * t["diff"], "lo_pp": 100 * t["lo"],
                        "hi_pp": 100 * t["hi"], "equivalent": bool(t["equivalent"])}
    return {"n": n, "rates_pct": rates, "observable_range_pp": obs_range,
            "n_discordant_total": total_disc,
            "n_discordant_max_contrast": max(disc.values()) if disc else 0,
            "mde_design_pp": mde_pp, "cochran_q": q_stat, "cochran_p": q_p,
            "cochran_df": q_df, "flag": flag, "equivalence": equiv}


# ── Figures ──────────────────────────────────────────────────────────────────────
def fig_yield_vs_cost(yield_df, cost_tbl, arms_exh, out: Path, n_complexes: int,
                      headline_depths: list[int]) -> None:
    """Yield against measured search cost, one panel per headline depth."""
    cost_by_e = cost_tbl[cost_tbl.timing_available].set_index("exhaustiveness")
    fig, axes = plt.subplots(1, len(headline_depths),
                             figsize=(5.7 * len(headline_depths), 5.0), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, depth in zip(axes, headline_depths):
        for gate, marker, ls in (("rmsd2", "o", "-"), ("triple", "s", "--")):
            pts = []
            for e in arms_exh:
                row = yield_df[(yield_df.exhaustiveness == e) & (yield_df.gate == gate)
                               & (yield_df.depth == depth)]
                if row.empty or e not in cost_by_e.index:
                    continue
                r = row.iloc[0]
                pts.append((cost_by_e.loc[e, "total_wall_h"], 100 * r["rate"], e,
                            bool(cost_by_e.loc[e, "thread_allocation_modal"]),
                            100 * r["wilson_lo"], 100 * r["wilson_hi"]))
            if len(pts) > 1:
                ax.plot([p[0] for p in pts], [p[1] for p in pts], ls,
                        color="#555555", lw=1.1, zorder=1)
            # Rungs can sit almost on top of each other in x when the cost ladder is
            # not monotone, so a label that always goes above would overprint the
            # neighbouring marker. Nudge the lower of any close pair downwards.
            xs = [p[0] for p in pts]
            span = (max(xs) - min(xs)) or 1.0
            for i, (c, y, e, modal, lo, hi) in enumerate(pts):
                close = [j for j, p in enumerate(pts)
                         if j != i and abs(p[0] - c) < 0.06 * span]
                below = any(pts[j][1] > y for j in close)
                ax.errorbar(c, y, yerr=[[y - lo], [hi - y]], fmt="none",
                            ecolor="#999999", elinewidth=1, capsize=2, zorder=2)
                ax.plot(c, y, marker=marker, ms=11, zorder=3, color=_color(e),
                        markeredgecolor=EDGE_LIGHT if modal else FLAG_COLOR,
                        markeredgewidth=1.0 if modal else 2.4)
                ax.annotate(f"{e}", (c, y), textcoords="offset points",
                            xytext=(0, -18 if below else 11), ha="center",
                            fontsize=8.5, zorder=4)
        ax.set_xlabel("Measured AutoDock Vina search time\n"
                      "(wall-clock hours over the whole set)")
        ax.set_title(_depth_label(depth), fontsize=10, pad=16)
        ax.grid(alpha=0.25, lw=0.6)
    axes[0].set_ylabel("Complexes with at least one qualifying pose\n"
                       f"(% of {n_complexes})")
    # Legend swatches use a mid-ramp colour so they read as "any rung", not as a
    # particular one: in these figures colour never encodes the gate.
    mid = LEGEND_GREY
    handles = [
        plt.Line2D([], [], marker="o", ls="-", color="#555555", markerfacecolor=mid,
                   markeredgecolor=EDGE_LIGHT, ms=9, label=GATE_LABEL["rmsd2"]),
        plt.Line2D([], [], marker="s", ls="--", color="#555555", markerfacecolor=mid,
                   markeredgecolor=EDGE_LIGHT, ms=9, label=GATE_LABEL["triple"]),
        plt.Line2D([], [], marker="o", ls="none", markerfacecolor=mid,
                   markeredgecolor=FLAG_COLOR, markeredgewidth=2.4, ms=9,
                   label="ran at a different thread allocation"),
        plt.Line2D([], [], ls="none", marker=None,
                   label="marker fill = exhaustiveness (light → dark with search effort)"),
    ]
    _label_panels(axes)
    fig.legend(handles=handles, loc="lower center", ncol=1, frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, -0.11))
    fig.suptitle("Docking success against the compute it cost — "
                 "raw AutoDock Vina exhaustiveness ladder", fontsize=11.5)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_marginal_return(mr_df, out: Path, headline_depths: list[int]) -> None:
    """Return per extra hour for each step of the ladder, with bootstrap intervals."""
    sub = mr_df[mr_df.gate.isin(HEADLINE_GATES)].copy()
    if sub.empty:
        return
    fig, axes = plt.subplots(1, len(headline_depths),
                             figsize=(5.7 * len(headline_depths), 4.7), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, depth in zip(axes, headline_depths):
        d = sub[sub.depth == depth]
        segs = list(dict.fromkeys(d["segment"]))
        width = 0.38
        # Colour is the rung the step arrives at; the gate is carried by hatch alone.
        # Colouring the two gates instead would put gate and exhaustiveness in the same
        # channel, and the two hexes previously used for it were literally the exh64 and
        # exh92 arm colours.
        seg_color = [_color(int(d[d.segment == s]["exh_hi"].iloc[0])) for s in segs]
        for gi, gate in enumerate(HEADLINE_GATES):
            g = d[d.gate == gate].set_index("segment").reindex(segs)
            x = np.arange(len(segs)) + (gi - 0.5) * width
            vals = g["rate_per_h"].values.astype(float)
            lo = np.abs(vals - g["rate_lo"].values.astype(float))
            hi = np.abs(g["rate_hi"].values.astype(float) - vals)
            ax.bar(x, vals, width, color=seg_color, edgecolor=EDGE_DARK, lw=0.7,
                   hatch=None if gate == "rmsd2" else "////")
            ax.errorbar(x, vals, yerr=[lo, hi], fmt="none", ecolor=EDGE_DARK,
                        elinewidth=1.1, capsize=3)
        ax.axhline(0, color="#333333", lw=0.9)
        # A segment with no measurable cost delta has no bar. Say so on the panel:
        # an empty slot otherwise reads as a measured zero.
        for si, seg in enumerate(segs):
            if d[d.segment == seg]["rate_per_h"].notna().any():
                continue
            ax.annotate("cost delta\nnot measurable", (si, 0),
                        textcoords="offset points", xytext=(0, 14), ha="center",
                        fontsize=7.5, color="#7a7a7a", style="italic")
        ax.set_xticks(np.arange(len(segs)))
        ax.set_xticklabels([s.replace("->", "→") for s in segs], fontsize=9)
        ax.set_xlabel("Step up the exhaustiveness ladder")
        ax.set_title(_depth_label(depth), fontsize=10, pad=16)
        ax.grid(axis="y", alpha=0.25, lw=0.6)
    axes[0].set_ylabel("Additional complexes recovered\nper extra wall-clock hour")
    _label_panels(axes)
    handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=LEGEND_GREY, edgecolor=EDGE_DARK,
                      label=GATE_LABEL["rmsd2"]),
        plt.Rectangle((0, 0), 1, 1, facecolor=LEGEND_GREY, edgecolor=EDGE_DARK,
                      hatch="////", label=GATE_LABEL["triple"]),
        plt.Rectangle((0, 0), 1, 1, facecolor="none", edgecolor="none",
                      label="bar fill = the rung the step arrives at"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.09))
    fig.suptitle("What each step up the ladder buys per hour "
                 "(whiskers are a paired complex bootstrap)", fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_raw_vs_rescored(strand_df, out: Path, headline_depths: list[int],
                        n_complexes: int) -> None:
    """Raw against gnina-rescored at each rung — deliberately with NO cost axis.

    Search hours and re-ranking hours are measured on different devices in different
    units, so putting both strands on one cost axis would invite a comparison the data
    cannot support. The x axis is exhaustiveness, which both strands genuinely share.
    """
    if strand_df.empty:
        return
    fig, axes = plt.subplots(1, len(headline_depths),
                             figsize=(5.7 * len(headline_depths), 4.7), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, depth in zip(axes, headline_depths):
        d = strand_df[strand_df.depth == depth]
        for gate, marker, ls in (("rmsd2", "o", "-"), ("triple", "s", "--")):
            g = d[d.gate == gate].sort_values("exhaustiveness")
            if g.empty:
                continue
            x = np.arange(len(g))
            for col, lw, fill in (("k_raw", 1.6, "white"), ("k_rescored", 2.2, None)):
                y = 100 * g[col].values / n_complexes
                ax.plot(x, y, ls, marker=marker, ms=9, lw=lw, color="#555555",
                        markerfacecolor=[_color(int(e)) for e in g.exhaustiveness][0]
                        if fill is None else fill,
                        markeredgecolor="#333333", zorder=3)
                # Raw sits ABOVE rescored so a coincident pair still shows both: at
                # oracle depth the two strands are near-identical and equal zorder would
                # hide the raw marker entirely.
                for xi, yi, e in zip(x, y, g.exhaustiveness):
                    ax.plot(xi, yi, marker=marker, ms=9,
                            zorder=6 if fill == "white" else 4,
                            color=_color(int(e)) if fill is None else "white",
                            markeredgecolor="#333333", markeredgewidth=1.2)
            ax.set_xticks(np.arange(len(g)))
            ax.set_xticklabels([f"exh{int(e)}" for e in g.exhaustiveness])
        ax.set_xlabel("Exhaustiveness rung")
        ax.set_title(_depth_label(depth), fontsize=10, pad=16)
        ax.grid(alpha=0.25, lw=0.6)
    axes[0].set_ylabel("Complexes with at least one qualifying pose\n"
                       f"(% of {n_complexes})")
    _label_panels(axes)
    handles = [
        plt.Line2D([], [], marker="o", ls="none", markerfacecolor="white",
                   markeredgecolor="#333333", ms=9, label="raw Vina ranking"),
        plt.Line2D([], [], marker="o", ls="none", markerfacecolor=LEGEND_GREY,
                   markeredgecolor="#333333", ms=9, label="+ gnina re-rank"),
        plt.Line2D([], [], marker="o", ls="-", color="#555555", markerfacecolor="none",
                   ms=0, label=GATE_LABEL["rmsd2"]),
        plt.Line2D([], [], marker="s", ls="--", color="#555555", markerfacecolor="none",
                   ms=0, label=GATE_LABEL["triple"]),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.10))
    fig.suptitle("Re-ranking a fixed pose pool, at each search effort", fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def fig_depth_sweep(sweep, out: Path) -> None:
    """Every rung's yield gain against ranking depth, one panel per headline gate.

    Sweeping all rungs rather than a chosen one avoids picking the segment after
    seeing which looks most interesting, and it is the panel that shows the depth
    dependence is a property of the ladder, not of one step.
    """
    if sweep.empty:
        return
    fig, axes = plt.subplots(1, len(HEADLINE_GATES),
                             figsize=(5.9 * len(HEADLINE_GATES), 4.5), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, gate in zip(axes, HEADLINE_GATES):
        g = sweep[sweep.gate == gate]
        for seg, sg in g.groupby("segment", sort=False):
            sg = sg.sort_values("depth")
            e_hi = int(sg["exh_hi"].iloc[0])
            ax.plot(sg["depth"], sg["delta_complexes"], marker="o", ms=3.8,
                    color=_color(e_hi), lw=1.7, markeredgecolor=EDGE_LIGHT,
                    markeredgewidth=0.5, label=seg.replace("->", " → "))
            sig = sg[sg["p_raw"] < 0.05]
            ax.plot(sig["depth"], sig["delta_complexes"], marker="o", ms=8.5,
                    ls="none", markerfacecolor="none", markeredgecolor=EDGE_DARK,
                    markeredgewidth=1.1)
        ax.axhline(0, color="#333333", lw=0.9)
        ax.set_xlabel("Ranking depth\n(number of top-ranked poses inspected)")
        ax.set_title(GATE_LABEL[gate], fontsize=10, pad=16)
        ax.grid(alpha=0.25, lw=0.6)
    axes[0].set_ylabel("Complexes gained over the rung below\n(of 303)")
    _label_panels(axes)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels), frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, -0.07),
               title="ringed markers: unadjusted McNemar p < 0.05", title_fontsize=8)
    fig.suptitle("What each rung buys, as a function of how many poses are inspected",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.93))
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


# ── Main ─────────────────────────────────────────────────────────────────────────
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv", type=Path, default=Path(DEFAULT_PER_POSE),
                    help="per_pose_metrics.csv from posebusters_pose_comparison.py. "
                         "Supplies every gate column and the Vina rank.")
    ap.add_argument("--pb-config", type=Path, default=Path(DEFAULT_PB_CONFIG),
                    help="PoseBusters run config. Its docking_directories block IS the "
                         "arm registry, so a rung added there joins this analysis "
                         "automatically once its poses are scored.")
    ap.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT))
    ap.add_argument("--method-prefix", default="autodock_mgltools",
                    help="Restricts the ladder to one preparation arm. The sweep ran on "
                         "MGLTools/ADFRsuite ligands, so admitting the Meeko-ligand arm "
                         "would change two variables at once.")
    ap.add_argument("--min-exhaustiveness", type=int, default=0)
    ap.add_argument("--max-exhaustiveness", type=int, default=10 ** 6)
    def _n_boot(v: str) -> int:
        n = int(v)
        if n < 1000:
            raise argparse.ArgumentTypeError(
                f"--n-boot must be at least 1000; {n} resamples cannot carry a 95% "
                "percentile interval. At n=1 the interval collapses onto one draw and "
                "the knee guard then tests that single draw.")
        return n

    ap.add_argument("--n-boot", type=_n_boot, default=5000,
                    help="Paired complex-level bootstrap resamples (minimum 1000).")
    ap.add_argument("--seed", type=int, default=BOOT_SEED)
    ap.add_argument("--knee-cell", default=None,
                    help="GATE:DEPTH for the one cell in which a knee is claimed "
                         "inferentially (depth 0 = oracle). Requires "
                         "--knee-cell-justification, which is echoed into the report. "
                         "Default is the pre-registered triple gate at oracle depth; "
                         "picking the best-looking cell after the fact inflates the "
                         "error rate roughly six-fold.")
    ap.add_argument("--knee-cell-justification", default=None)
    ap.add_argument("--no-optimizer-strand", action="store_true",
                    help="Analyse the raw search ladder only. By default a rescored "
                         "strand is added for every rung whose <key>_gnina variant is "
                         "scored; it shares the raw arm's poses, so it is a second "
                         "strand at the same exhaustiveness, never a further rung.")
    ap.add_argument("--trend-supplement", action="store_true",
                    help="Add a Cochran-Armitage trend test with log2(exhaustiveness) "
                         "scores. Off by default: on a paired design it is valid but "
                         "extremely conservative, so a null result there is a floor, "
                         "not evidence of no trend.")
    ap.add_argument("--verify-input-parity", action="store_true",
                    help="sha256 the staged ligand and box per complex per arm, proving "
                         "exhaustiveness is the only variable. Add --full-parity to "
                         "hash the receptors too.")
    ap.add_argument("--full-parity", action="store_true")
    ap.add_argument("--cascade-pins-from", type=Path, default=None,
                    help="A hub pose_validity_cascade.csv whose <method-prefix>_exh18 "
                         "and _exh92 rows become the cascade pins for this run's "
                         "reference convention, replacing the built-in CASCADE_PINS "
                         "entry. Required for a hard self-check under the nearest "
                         "convention, which has no published rows yet (plan 3.1a).")
    ap.add_argument("--write-pins", type=Path, default=None,
                    help="Write the pins this run compared against (and where they came "
                         "from) to this JSON path.")
    mf.add_method_filter_args(ap)
    args = ap.parse_args(argv)

    root = _ROOT
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    per_pose_path = (args.per_pose_csv if args.per_pose_csv.is_absolute()
                     else root / args.per_pose_csv)
    pb_config = args.pb_config if args.pb_config.is_absolute() else root / args.pb_config

    knee_gate, knee_depth = KNEE_CELL
    if args.knee_cell:
        if not args.knee_cell_justification:
            raise SystemExit("--knee-cell requires --knee-cell-justification: moving the "
                             "confirmatory cell after seeing the data is the selection "
                             "effect this pre-registration exists to prevent.")
        g, _, d = args.knee_cell.partition(":")
        knee_gate, knee_depth = g.strip(), int(d)

    print(f"reading {per_pose_path}")
    head = pd.read_csv(per_pose_path, nrows=0)
    want = ["method", "protein", "ligand", "rank", "rmsd", "pb_valid", "bestfit_rmsd",
            "pb_kabsch_rmsd", "autodock_rank", "optimized_rank", "autodock_affinity", "optimizer",
            "reference_convention"]
    per_pose = pd.read_csv(per_pose_path,
                           usecols=[c for c in want if c in head.columns],
                           low_memory=False)

    # Which deposited ligand copy the RMSD-derived gate columns were scored against.
    # Only a table that declares the convention (or a run that supplies pins) says so in
    # its outputs: a table from before the column exists is instance by definition and
    # its outputs stay byte-identical to the published run.
    convention, conv_declared = detect_reference_convention(per_pose)
    conv_note = (f"{convention} (per_pose_metrics.csv column reference_convention)"
                 if conv_declared else
                 f"{convention} (table predates the reference_convention column)")
    if conv_declared:
        print(f"reference convention: {conv_note}")
    if args.cascade_pins_from is not None:
        pins_path = (args.cascade_pins_from if args.cascade_pins_from.is_absolute()
                     else root / args.cascade_pins_from)
        cascade_pins = load_cascade_pins(pins_path, args.method_prefix)
        pins_source = str(pins_path)
        print(f"cascade pins: {convention} convention, read from {pins_path}")
    else:
        cascade_pins = CASCADE_PINS[convention]
        pins_source = "built-in CASCADE_PINS[%r]" % convention if cascade_pins else None
        if cascade_pins is None:
            print(f"  WARNING: no built-in cascade pins for the {convention!r} "
                  f"convention and no --cascade-pins-from given; the cascade self-check "
                  f"CSV will be written but NOT asserted, so this run is not pinned "
                  f"to any published cascade")
    announce_conv = conv_declared or args.cascade_pins_from is not None
    if args.write_pins is not None:
        wp = args.write_pins if args.write_pins.is_absolute() else root / args.write_pins
        wp.parent.mkdir(parents=True, exist_ok=True)
        wp.write_text(json.dumps({
            "reference_convention": convention,
            "convention_declared_by_table": conv_declared,
            "per_pose_csv": str(per_pose_path),
            "pins_source": pins_source,
            "pin_order": CASCADE_PIN_ORDER,
            "pins": ({str(e): list(v) for e, v in cascade_pins.items()}
                     if cascade_pins else None),
        }, indent=2) + "\n")
        print(f"wrote pins to {wp}")

    # Before discover_arms: the ladder registry is built by enumerating method
    # keys off this frame, so filtering later would leave excluded rungs in the
    # registry and in its duplicate-rung guard.
    per_pose = mf.apply_method_filter(per_pose, "method", args,
                                      label="exhaustiveness-returns",
                                      out_dir=args.out_dir)
    arms_all = discover_arms(per_pose, pb_config, root, args.method_prefix)
    for a in arms_all:
        if a["optimization"] not in ("none", "", "null", "false"):
            print(f"  note: {a['method_key']} declares optimization="
                  f"{a['optimization']!r}; only its RAW poses and its Vina search timing "
                  f"are used (the optimiser runs after the timed interval)")
    arms = [a for a in arms_all if a["has_poses"] and a["has_timing"]
            and args.min_exhaustiveness <= a["exhaustiveness"] <= args.max_exhaustiveness]
    pending = [a for a in arms_all if a["has_timing"] and not a["has_poses"]]
    if len(arms) < 2:
        raise SystemExit(f"need at least two scored arms with timing; found "
                         f"{[a['method_key'] for a in arms]}")

    # The ladder is indexed BY exhaustiveness, so two arms at one setting would be
    # pooled into a single rung and silently mix preparation or scoring arms. A wider
    # --method-prefix is the way this happens.
    seen: dict[int, str] = {}
    for a in arms:
        e = a["exhaustiveness"]
        if e in seen:
            raise SystemExit(
                f"--method-prefix {args.method_prefix!r} admits two arms at "
                f"exhaustiveness {e} ({seen[e]} and {a['method_key']}). This ladder is "
                f"indexed by exhaustiveness, so a duplicate rung would pool different "
                f"arms into one row. Narrow the prefix — the default "
                f"'autodock_mgltools' gives one arm per setting.")
        seen[e] = a["method_key"]

    arms_exh = [a["exhaustiveness"] for a in arms]
    print(f"ladder: {', '.join('exh%d' % e for e in arms_exh)}")
    for a in pending:
        print(f"  pending: exh{a['exhaustiveness']} is docked and timed but not scored "
              f"yet — it joins once {a['method_key']} reaches per_pose_metrics.csv")

    # The rescored counterpart of each rung, where one exists. It shares the raw arm's
    # poses and inputs, so it is a SECOND STRAND at the same exhaustiveness rather than
    # a new rung — the ladder axis and the re-ranking axis are kept orthogonal.
    opt_arms = [] if args.no_optimizer_strand else discover_optimizer_arms(
        arms, per_pose, root)   # cost re-read below, once the cohort is known
    if opt_arms:
        print(f"rescored strand: {', '.join('exh%d' % a['exhaustiveness'] for a in opt_arms)}"
              f"  (+{OPTIMIZER_SUFFIX})")
        for a in opt_arms:
            if a["optimizer_rows"] == 0:
                print(f"  warning: exh{a['exhaustiveness']}{OPTIMIZER_SUFFIX} has scored "
                      f"poses but no optimization_log.csv — its cost cannot be reported")

    keys = {a["method_key"]: a["exhaustiveness"] for a in arms}
    strands = {a["method_key"]: RAW_STRAND for a in arms}
    for a in opt_arms:
        keys[a["method_key"]] = a["exhaustiveness"]
        strands[a["method_key"]] = OPTIMIZER_STRAND
    df = per_pose[per_pose["method"].isin(keys)].copy()
    df["exhaustiveness"] = df["method"].map(keys)
    df["strand"] = df["method"].map(strands)
    df = apply_gates(df)
    # The ladder's own cohort is defined by the RAW strand; a rescored arm inherits it.
    raw_df = df[df["strand"] == RAW_STRAND]

    sets = {e: set(g["protein"].astype(str)) for e, g in raw_df.groupby("exhaustiveness")}
    complexes = sorted(set.intersection(*sets.values()))
    dropped_cplx = sorted(set.union(*sets.values()) - set(complexes))
    df = df[df["protein"].astype(str).isin(complexes)]
    raw_df = df[df["strand"] == RAW_STRAND]
    n_cplx = len(complexes)
    print(f"paired cohort: {n_cplx} complexes"
          + (f"  (dropped {len(dropped_cplx)} not present in every arm)"
             if dropped_cplx else ""))

    # The cohort is only known now, so the optimiser cost is re-read against it —
    # discovery above ran before the paired set existed.
    if opt_arms:
        opt_arms = discover_optimizer_arms(arms, per_pose, root, cohort=set(complexes))

    audit = audit_gate_inputs(raw_df)
    for c in audit:
        c["strand"] = RAW_STRAND
    # The rescored strand has its own invariants: it must rank on the optimiser's score,
    # which is precisely what the raw checks above would flag as broken.
    opt_df = df[df["strand"] == OPTIMIZER_STRAND]
    if not opt_df.empty:
        if "optimized_rank" in opt_df.columns:
            orank = pd.to_numeric(opt_df["optimized_rank"], errors="coerce")
            frac = float((orank == opt_df["rank"]).mean())
            audit.append({"check": "rescored rank is the optimiser rank",
                          "ok": bool(frac == 1.0), "strand": OPTIMIZER_STRAND,
                          "detail": f"rank == optimized_rank in {100 * frac:.2f}% of poses"})
        if "optimizer" in opt_df.columns:
            tools = sorted(set(opt_df["optimizer"].dropna().astype(str)))
            audit.append({"check": "rescored poses carry an optimiser tag",
                          "ok": bool(tools and "original" not in tools),
                          "strand": OPTIMIZER_STRAND,
                          "detail": f"optimizer values: {tools}"})
        shared = (raw_df.groupby("exhaustiveness").size().reindex(
                      sorted(opt_df["exhaustiveness"].unique()))
                  == opt_df.groupby("exhaustiveness").size())
        audit.append({"check": "rescored arm has the same pose count as its raw rung",
                      "ok": bool(shared.all()), "strand": OPTIMIZER_STRAND,
                      "detail": f"per-rung equality: {shared.to_dict()}"})
    if announce_conv:
        # Which copy of the deposited ligand the gate columns were scored against. Not a
        # pass/fail check; it records the convention beside the invariants it governs.
        audit.append({"check": "reference convention of the gate columns", "ok": True,
                      "strand": RAW_STRAND,
                      "detail": f"{conv_note}; rmsd / bestfit_rmsd / pb_valid are "
                                f"evaluated against the {convention} deposited copy"})
    for c in audit:
        if not c["ok"]:
            print(f"  WARNING [{c['strand']}/{c['check']}]: {c['detail']}")
    pd.DataFrame(audit).to_csv(out_dir / "exh_gate_input_audit.csv", index=False)

    # ── cascade regression assertion ────────────────────────────────────────────
    casc = cascade_rows(raw_df)
    casc.to_csv(out_dir / "exh_cascade_selfcheck.csv", index=False)
    checked: list[int] = []
    skipped: list[int] = []
    if cascade_pins is None:
        print(f"  WARNING: cascade self-check written to exh_cascade_selfcheck.csv but "
              f"NOT asserted — no pins for the {convention!r} convention")
    else:
        problems = assert_cascade(casc, cascade_pins)
        if problems:
            raise SystemExit("CASCADE REGRESSION — this script no longer reproduces the "
                             "published pose_validity_cascade table, so every number below "
                             "it would be unmoored from the thesis:\n  "
                             + "\n  ".join(problems))
        # Name only the pins that were actually compared. assert_cascade skips a pin
        # whose arm is absent, so announcing the whole pin list would claim a
        # regression check that did not run.
        checked = sorted(e for e in cascade_pins if (casc.exhaustiveness == e).any())
        skipped = sorted(set(cascade_pins) - set(checked))
        if checked:
            print("cascade self-check: reproduces the published rows for "
                  + ", ".join("exh%d" % e for e in checked))
        if skipped:
            print("  WARNING: pinned arm(s) "
                  + ", ".join("exh%d" % e for e in skipped)
                  + " are absent from this ladder, so their pins were NOT evaluated — "
                    "this run is not pinned to the published cascade")

    # depths: oracle is derived, never assumed to be 30
    max_rank = int(raw_df["rank"].max())
    depths = [d for d in DEPTHS if d < max_rank] + [0]
    headline_depths = [1, 0]
    oracle_is_top30 = max_rank == 30

    # ── timing and cost ─────────────────────────────────────────────────────────
    timing = {}
    for a in arms:
        t = load_timing(a["tree"])
        if not t.empty:
            bad = t[t["status"].str.lower() != "success"]
            if len(bad):
                print(f"  warning: exh{a['exhaustiveness']} has {len(bad)} non-success "
                      f"docking rows; excluded from the cost")
                t = t[t["status"].str.lower() == "success"]
            logged = t["log_exhaustiveness"].dropna().unique()
            if len(logged) == 1 and int(logged[0]) != a["exhaustiveness"]:
                raise SystemExit(f"exh{a['exhaustiveness']} tree logs exhaustiveness="
                                 f"{int(logged[0])}: config and run disagree")
        timing[a["exhaustiveness"]] = t

    cost, cost_tbl, cost_meta = build_cost_table(arms, timing, complexes)
    cost_tbl.to_csv(out_dir / "exh_cost_per_arm.csv", index=False)
    long_cost = (cost.reset_index().rename(columns={"index": "protein"})
                 .melt(id_vars="protein", var_name="exhaustiveness",
                       value_name="elapsed_time_s"))
    long_cost.to_csv(out_dir / "exh_cost_per_complex.csv", index=False)
    cost_arms = [e for e in arms_exh if e in cost.columns]

    parity = None
    if args.verify_input_parity:
        print("hashing staged inputs to prove exhaustiveness is the only variable ...")
        parity = verify_input_parity(arms, complexes, full=args.full_parity)
        for kind, r in parity.items():
            state = ("identical" if r["identical"]
                     else f"MISMATCH on {r['n_mismatched']}" if r["n_mismatched"]
                     else f"UNPROVEN — {r['n_partial']} complexes staged by fewer than "
                          f"all {r['n_arms']} arms")
            print(f"  {kind:9s} {state} across {r['n_arms']} arms "
                  f"({r['checked']} complexes)")

    # ── yield ───────────────────────────────────────────────────────────────────
    mats: dict[tuple[str, int], pd.DataFrame] = {}
    yrows = []
    for gate in GATES:
        for depth in depths:
            m = complex_matrix(raw_df, gate, depth, arms_exh, complexes)
            mats[(gate, depth)] = m
            for e in arms_exh:
                k = int(m[e].sum())
                lo, hi = su.wilson_ci(k, n_cplx)
                yrows.append({"exhaustiveness": e, "gate": gate, "depth": depth,
                              "depth_label": _depth_label(depth), "k_complexes": k,
                              "n_complexes": n_cplx, "rate": k / n_cplx,
                              "wilson_lo": lo, "wilson_hi": hi,
                              "is_headline_depth": depth in headline_depths})
    yield_df = pd.DataFrame(yrows)
    yield_df.to_csv(out_dir / "exh_yield_by_depth.csv", index=False)

    prows = []
    for gate in GATES:
        for e, g in raw_df.groupby("exhaustiveness"):
            est, lo, hi = su.cluster_bootstrap_ci(
                g[gate].values.astype(float), g["protein"].values,
                n_boot=2000, seed=args.seed)
            prows.append({"exhaustiveness": e, "gate": gate, "n_poses": len(g),
                          "pose_share": est, "ci_lo": lo, "ci_hi": hi})
    pd.DataFrame(prows).to_csv(out_dir / "exh_pose_level_share.csv", index=False)

    # ── saturation ──────────────────────────────────────────────────────────────
    srows = []
    for gate in GATES:
        for depth in depths:
            c = saturation_cell(mats[(gate, depth)], arms_exh)
            rec = {"gate": gate, "depth": depth, "depth_label": _depth_label(depth),
                   "flag": c["flag"], "n": c["n"],
                   "rate_min_pct": min(c["rates_pct"].values()),
                   "rate_max_pct": max(c["rates_pct"].values()),
                   "observable_range_pp": c["observable_range_pp"],
                   "n_discordant_total": c["n_discordant_total"],
                   "n_discordant_max_contrast": c["n_discordant_max_contrast"],
                   "mde_design_pp": c["mde_design_pp"],
                   "cochran_q": c["cochran_q"], "cochran_p": c["cochran_p"]}
            for m, v in c["equivalence"].items():
                # No dots in column names: DataFrame.itertuples() renames any column
                # that is not a valid identifier to a positional _N, and a getattr
                # lookup on the intended name then silently returns the default.
                tag = f"tost_{100 * m:g}pp".replace(".", "p")
                rec[f"{tag}_diff_pp"] = v["diff_pp"]
                rec[f"{tag}_lo_pp"] = v["lo_pp"]
                rec[f"{tag}_hi_pp"] = v["hi_pp"]
                rec[f"{tag}_equivalent"] = v["equivalent"]
            srows.append(rec)
    sat_df = pd.DataFrame(srows)
    sat_df.to_csv(out_dir / "exh_saturation_diagnostics.csv", index=False)
    degenerate = {(r.gate, r.depth) for r in sat_df.itertuples() if r.flag == "DEGENERATE"}

    # ── adjacent-rung yield contrasts, with the multiplicity families ───────────
    crows = []
    for gate in GATES:
        for depth in depths:
            m = mats[(gate, depth)]
            for lo_e, hi_e in zip(arms_exh, arms_exh[1:]):
                a, b = m[lo_e].values, m[hi_e].values
                n10, n01, p = su.mcnemar_exact(a, b)
                ci = su.newcombe_paired_diff_ci(a, b)
                crows.append({
                    "gate": gate, "depth": depth, "depth_label": _depth_label(depth),
                    "segment": f"exh{lo_e}->exh{hi_e}", "exh_lo": lo_e, "exh_hi": hi_e,
                    "k_lo": int(a.sum()), "k_hi": int(b.sum()),
                    "delta_complexes": int(b.sum() - a.sum()),
                    "gained": int((b & ~a).sum()), "lost": int((a & ~b).sum()),
                    "delta_pp": 100 * ci["diff"], "ci_lo_pp": 100 * ci["lo"],
                    "ci_hi_pp": 100 * ci["hi"], "n_discordant": n10 + n01, "p_raw": p,
                    "degenerate_cell": (gate, depth) in degenerate})
    contrasts = pd.DataFrame(crows)
    # Family A, confirmatory: informative gates x headline depths, Holm.
    # Degenerate cells are excluded BEFORE the correction — mcnemar_exact returns
    # p = 1.0 (not NaN) with zero discordant pairs, and holm only auto-drops NaN, so
    # leaving them in would inflate the family size and deflate every real contrast.
    famA = contrasts[(contrasts.gate.isin(HEADLINE_GATES))
                     & (contrasts.depth.isin(headline_depths))
                     & (~contrasts.degenerate_cell)]
    famB = contrasts[~contrasts.index.isin(famA.index) & (~contrasts.degenerate_cell)]
    contrasts["family"] = np.where(contrasts.index.isin(famA.index), "A-confirmatory",
                                   np.where(contrasts.index.isin(famB.index),
                                            "B-exploratory", "C-degenerate"))
    contrasts["p_adjusted"] = np.nan
    if len(famA):
        contrasts.loc[famA.index, "p_adjusted"] = su.holm(famA["p_raw"].values)
    if len(famB):
        contrasts.loc[famB.index, "p_adjusted"] = su.bh_fdr(famB["p_raw"].values)
    contrasts.to_csv(out_dir / "exh_yield_contrasts.csv", index=False)

    # ── marginal return ─────────────────────────────────────────────────────────
    segments = list(zip(cost_arms, cost_arms[1:]))
    mrows = []
    for gate in GATES:
        for depth in depths:
            m = mats[(gate, depth)]
            for lo_e, hi_e in segments:
                mr = marginal_return(m[lo_e].values, m[hi_e].values,
                                     cost[lo_e].values, cost[hi_e].values,
                                     n_boot=args.n_boot, seed=args.seed)
                recip, note = reciprocal_or_refusal(mr)
                d_exh = hi_e - lo_e
                is_knee_cell = (gate == knee_gate and depth == knee_depth)
                mrows.append({
                    "gate": gate, "depth": depth, "depth_label": _depth_label(depth),
                    "segment": f"exh{lo_e}->exh{hi_e}", "exh_lo": lo_e, "exh_hi": hi_e,
                    "delta_exhaustiveness": d_exh, **mr,
                    "delta_cost_h_per_exh_unit": mr["delta_cost_h"] / d_exh,
                    "wall_h_per_extra_complex": recip, "reciprocal_status": note,
                    "family": "D-confirmatory" if is_knee_cell else "D-descriptive",
                    "inference_claimable": is_knee_cell,
                    "is_headline": gate in HEADLINE_GATES and depth in headline_depths})
    mr_df = pd.DataFrame(mrows)
    mr_df.to_csv(out_dir / "exh_marginal_return.csv", index=False)

    # ── diminishing returns ─────────────────────────────────────────────────────
    # The contrast needs two segments whose cost delta is actually a cost. A segment
    # whose measured cost went DOWN cannot anchor a return-per-hour comparison, so the
    # last usable segment is used and named explicitly rather than silently assumed to
    # be the top of the ladder.
    usable = [(lo, hi) for lo, hi in segments
              if (cost[hi] - cost[lo]).sum() / 3600.0 > COST_FLOOR_H]
    unusable = [s for s in segments if s not in usable]
    for lo_e, hi_e in unusable:
        print(f"  note: segment exh{lo_e}->exh{hi_e} has a measured cost delta of "
              f"{(cost[hi_e] - cost[lo_e]).sum() / 3600.0:+.3f} wall-h and cannot "
              f"anchor a return-per-hour statistic")
    drows = []
    if len(usable) >= 2:
        first, last = usable[0], usable[-1]
        for gate in GATES:
            for depth in depths:
                m = mats[(gate, depth)]
                packs = [{"dy": m[hi].values.astype(float) - m[lo].values.astype(float),
                          "dc": cost[hi].values - cost[lo].values}
                         for lo, hi in (first, last)]
                d = diminishing_contrast(packs[0], packs[1], n_boot=args.n_boot,
                                         seed=args.seed)
                drows.append({
                    "gate": gate, "depth": depth, "depth_label": _depth_label(depth),
                    "first_segment": f"exh{first[0]}->exh{first[1]}",
                    "last_segment": f"exh{last[0]}->exh{last[1]}",
                    "rate_first_per_h": packs[0]["dy"].sum() / (packs[0]["dc"].sum() / 3600),
                    "rate_last_per_h": packs[1]["dy"].sum() / (packs[1]["dc"].sum() / 3600),
                    **d,
                    "inference_claimable": (gate == knee_gate and depth == knee_depth),
                    "is_headline": gate in HEADLINE_GATES and depth in headline_depths})
    dim_df = pd.DataFrame(drows)
    dim_df.to_csv(out_dir / "exh_diminishing_returns.csv", index=False)

    # ── depth sweep, every rung ─────────────────────────────────────────────────
    # Swept for all segments rather than a chosen one: picking the segment after
    # seeing which looks most interesting is exactly the selection this report
    # otherwise takes care to avoid.
    swrows = []
    for depth in range(1, max_rank + 1):
        mats_d = {g: complex_matrix(raw_df, g, depth, arms_exh, complexes) for g in GATES}
        for gate in GATES:
            m = mats_d[gate]
            for lo_e, hi_e in zip(arms_exh, arms_exh[1:]):
                a, b = m[lo_e].values, m[hi_e].values
                n10, n01, p = su.mcnemar_exact(a, b)
                swrows.append({"gate": gate, "depth": depth,
                               "segment": f"exh{lo_e}->exh{hi_e}",
                               "exh_lo": lo_e, "exh_hi": hi_e,
                               "k_lo": int(a.sum()), "k_hi": int(b.sum()),
                               "delta_complexes": int(b.sum() - a.sum()),
                               "n_discordant": n10 + n01, "p_raw": p})
    sweep = pd.DataFrame(swrows)
    sweep.to_csv(out_dir / "exh_depth_sweep.csv", index=False)

    # ── raw versus rescored, at equal exhaustiveness ────────────────────────────
    # Paired on the same complexes AND the same underlying poses: re-ranking reorders a
    # fixed pool, so any rank-1 gain is the CNN promoting a pose the Vina score buried,
    # never a pose the search did not find. At oracle depth the two strands can only
    # differ through the optimiser's minimisation, which moves geometry slightly.
    srows = []
    opt_by_exh = {a["exhaustiveness"]: a for a in opt_arms}
    for e, oa in sorted(opt_by_exh.items()):
        for gate in GATES:
            for depth in depths:
                mr = complex_matrix(raw_df, gate, depth, [e], complexes)[e].values
                mo = complex_matrix(df[df["strand"] == OPTIMIZER_STRAND], gate, depth,
                                    [e], complexes)[e].values
                n10, n01, pv = su.mcnemar_exact(mr, mo)
                ci = su.newcombe_paired_diff_ci(mr, mo)
                srows.append({
                    "exhaustiveness": e, "gate": gate, "depth": depth,
                    "depth_label": _depth_label(depth),
                    "k_raw": int(mr.sum()), "k_rescored": int(mo.sum()),
                    "delta_complexes": int(mo.sum() - mr.sum()),
                    "gained": int((mo & ~mr).sum()), "lost": int((mr & ~mo).sum()),
                    "delta_pp": 100 * ci["diff"], "ci_lo_pp": 100 * ci["lo"],
                    "ci_hi_pp": 100 * ci["hi"], "n_discordant": n10 + n01, "p_raw": pv,
                    "optimizer_compute_h": oa["optimizer_compute_h"],
                    "optimize_workers": oa["optimize_workers"],
                    "optimize_cpu": oa.get("optimize_cpu"),
                    "gpu_max_concurrent": oa.get("gpu_max_concurrent"),
                    "is_headline": gate in HEADLINE_GATES and depth in headline_depths})
    strand_df = pd.DataFrame(srows)
    if strand_df.empty:
        # A previous run's strand outputs must not survive a run that has no strand,
        # or the directory asserts an analysis this report explicitly denies making.
        for stale in (out_dir / "exh_raw_vs_rescored.csv",
                      fig_dir / "exh_raw_vs_rescored.png",
                      fig_dir / "exh_raw_vs_rescored_stats.txt"):
            if stale.exists():
                stale.unlink()
                print(f"  removed stale {stale.name} (no rescored strand in this run)")
    if not strand_df.empty:
        # Confirmatory family: headline gates x headline depths x rungs, Holm-corrected.
        fam = strand_df[strand_df.is_headline & (strand_df.n_discordant > 0)]
        strand_df["p_adjusted"] = np.nan
        if len(fam):
            strand_df.loc[fam.index, "p_adjusted"] = su.holm(fam["p_raw"].values)
        strand_df.to_csv(out_dir / "exh_raw_vs_rescored.csv", index=False)

    trend = pd.DataFrame()
    if args.trend_supplement:
        trows = []
        for gate in GATES:
            for depth in depths:
                m = mats[(gate, depth)]
                ks = [int(m[e].sum()) for e in arms_exh]
                try:
                    z, p, _ = su.cochran_armitage(ks, [n_cplx] * len(arms_exh),
                                                  scores=np.log2(arms_exh))
                except Exception:
                    z, p = np.nan, np.nan
                trows.append({"gate": gate, "depth": depth, "z": z, "p_raw": p,
                              "label": "CONSERVATIVE — a null here is a floor, not "
                                       "evidence of no trend"})
        trend = pd.DataFrame(trows)
        trend.to_csv(out_dir / "exh_trend_supplement.csv", index=False)

    # ── figures ─────────────────────────────────────────────────────────────────
    fig_yield_vs_cost(yield_df, cost_tbl, arms_exh,
                      fig_dir / "exh_yield_vs_cost.png", n_cplx, headline_depths)
    fig_marginal_return(mr_df, fig_dir / "exh_marginal_return.png", headline_depths)
    fig_depth_sweep(sweep, fig_dir / "exh_depth_sweep.png")
    fig_raw_vs_rescored(strand_df, fig_dir / "exh_raw_vs_rescored.png",
                        headline_depths, n_cplx)

    # ── report ──────────────────────────────────────────────────────────────────
    L: list[str] = []
    A = L.append
    A("Marginal return on compute — raw AutoDock Vina exhaustiveness ladder")
    A("#" * 78)
    A(f"generated by Scripts/Analysis/{Path(__file__).name}")
    A("")
    A("SOURCES")
    A(f"  per-pose metrics   {per_pose_path}")
    if announce_conv:
        A(f"  reference conv.    {conv_note}")
        A(f"  cascade pins       {pins_source or 'NONE — self-check not asserted'}")
    A(f"  arm registry       {pb_config}")
    for a in arms:
        A(f"  exh{a['exhaustiveness']:<4d} timing     {a['tree']}/*/*/docking_log_*.csv")
        if a["config"]:
            A(f"  exh{a['exhaustiveness']:<4d} config     {a['config']}")
    A("")
    A("SCOPE")
    A(f"  {n_cplx} complexes covered by every arm; Vina ranks 1..{max_rank}.")
    if opt_arms:
        A(f"  Two strands: the RAW search ladder (all cost and marginal-return numbers) "
          f"and its")
        A(f"  gnina-rescored counterpart at {len(opt_arms)} rung(s), compared in section 9. "
          f"A rescored arm")
        A("  shares its raw rung's poses, so it is never treated as a further rung.")
    else:
        A("  RAW poses only. Post-dock optimiser variants are a different axis and are "
          "excluded.")
    if oracle_is_top30:
        A(f"  NOTE: num_modes = {max_rank} in every arm, so top-{max_rank} and the "
          f"oracle are the SAME slice.")
        A("  They are one measurement, not two, and are not independent evidence of a "
          "plateau.")
    for a in pending:
        A(f"  PENDING: exh{a['exhaustiveness']} is docked and timed but not yet scored. "
          f"It joins this")
        A(f"  ladder automatically once {a['method_key']} appears in per_pose_metrics.csv.")
    if dropped_cplx:
        A(f"  {len(dropped_cplx)} complexes are not present in every arm and were "
          f"dropped from both the")
        A(f"  yield and the cost: {', '.join(dropped_cplx[:8])}"
          + (" ..." if len(dropped_cplx) > 8 else ""))
    if parity:
        for kind, r in parity.items():
            if r["identical"]:
                msg = (f"identical across all {r['n_arms']} arms on {r['checked']} "
                       f"complexes — exhaustiveness is the only variable")
            elif r["n_mismatched"]:
                msg = (f"MISMATCH on {r['n_mismatched']} complexes — this ladder is NOT "
                       f"a single-variable comparison")
            else:
                msg = (f"UNPROVEN — {r['n_partial']} complexes are staged by fewer than "
                       f"all {r['n_arms']} arms, so agreement was not testable there")
            A(f"  input parity ({kind}): {msg}")
    A("")

    if cascade_pins is None:
        A(f"1. CASCADE SELF-CHECK — NOT ASSERTED: no cascade pins for the {convention!r} "
          f"convention")
        A("   (pass --cascade-pins-from <pose_validity_cascade.csv> of the matching "
          "hub run to pin it)")
    elif args.cascade_pins_from is not None:
        A(f"1. CASCADE SELF-CHECK — reproduces {pins_source}")
    else:
        A("1. CASCADE SELF-CHECK — reproduces pose_validity_cascade_report.txt")
    A("")
    rows = []
    for r in casc.itertuples():
        d = r._asdict()
        rows.append([f"exh{int(r.exhaustiveness)}", f"{int(r.n_complexes)}",
                     f"{int(r.n_poses):,}", f"{int(r.pb_valid_complexes)}",
                     f"{int(r.pb_valid_poses):,}", f"{d['pb_valid_pct']:.1f}",
                     f"{int(r.rmsd2_complexes)}", f"{int(r.rmsd2_poses):,}",
                     f"{d['rmsd2_pct']:.1f}", f"{int(r.kabsch1_complexes)}",
                     f"{int(r.kabsch1_poses):,}", f"{d['kabsch1_pct']:.1f}",
                     f"{int(r.triple_complexes)}", f"{int(r.triple_poses):,}",
                     f"{d['triple_of_rmsd2_pct']:.1f}", f"{d['triple_of_all_pct']:.1f}"])
    A("   " + _fw(rows, ["arm", "cplx", "poses", "cplx", "poses", "%", "cplx", "poses",
                         "%", "cplx", "poses", "%", "cplx", "poses", "%of≤2Å", "%all"],
                  "<" + ">" * 15).replace("\n", "\n   "))
    A("        produced   PoseBusters-valid      RMSD ≤ 2 Å        Kabsch < 1 Å"
      "         triple gate")
    A("")
    A("   Complex counts are 'at least one pose passing'; every percentage is "
      "pose-level.")
    A("")

    A("2. WHAT EACH ARM COST")
    A("")
    A(f"   Measured Vina-subprocess wall time, summed over the same {n_cplx} complexes")
    A("   the gates use.")
    A("")
    rows = []
    for r in cost_tbl.itertuples():
        if not getattr(r, "timing_available", False):
            rows.append([f"exh{int(r.exhaustiveness)}"] + ["—"] * 7)
            continue
        rows.append([
            f"exh{int(r.exhaustiveness)}",
            f"{int(r.cpus_per_worker)}" if r.cpus_per_worker else "—",
            f"{r.total_wall_h:.3f}", f"{r.median_wall_s:.1f}",
            f"{r.trimmed_mean_10pct_wall_s:.1f}",
            f"{r.paired_median_ratio_to_prev:.3f}"
            if r.paired_median_ratio_to_prev == r.paired_median_ratio_to_prev else "—",
            f"{r.arm_effect_pct:+.1f}%" if r.arm_effect_pct == r.arm_effect_pct else "—",
            f"{r.total_core_h:.1f}" if r.total_core_h == r.total_core_h else "—"])
    A("   " + _fw(rows, ["arm", "threads", "wall h", "median s", "trim-mean s",
                         "×prev", "machine", "core h"],
                  "<" + ">" * 7).replace("\n", "\n   "))
    A("")
    A("   'median s' and 'trim-mean s' are per complex. '×prev' is the paired "
      "per-complex")
    A("   median ratio to the rung below — a robust step size that a single slow "
      "complex")
    A("   cannot invert, unlike the total.")
    A("")
    eff = cost_meta["arm_effects"]
    if eff.get("ok"):
        A("   'machine' is the residual arm effect on log wall time after removing the "
          "per-complex")
        A("   effect and the fitted exhaustiveness trend — the part of an arm's cost "
          "its search")
        A(f"   effort does not explain. Spread across the ladder: "
          f"{eff['spread_pp']:.1f} percentage points.")
        A("")
        A("   THIS IS THE DOMINANT LIMITATION OF THE COST SIDE. Each arm is one batch, "
          "run on")
        A("   one day, with no replication and no background-load telemetry. The batch "
          "effect is")
        A("   the same size as the treatment effect and is perfectly aliased with it. "
          "No cost")
        A("   number here separates 'exhaustiveness 92 is expensive' from 'the machine "
          "was busy")
        A("   on the afternoon exhaustiveness 92 ran'.")
        A("")
    inverted = [(lo, hi) for lo, hi in segments
                if (cost[hi] - cost[lo]).sum() / 3600.0 <= COST_FLOOR_H]
    if inverted:
        A("   COST INVERSION — the ladder is not monotone in measured wall time:")
        for lo_e, hi_e in inverted:
            pair = cost[[lo_e, hi_e]].dropna()
            delta = pair[hi_e] - pair[lo_e]
            d = float(delta.sum()) / 3600.0
            med = float((pair[hi_e] / pair[lo_e]).median())
            slower = int((delta > 0).sum())
            try:
                from scipy.stats import binomtest, wilcoxon
                p_sign = float(binomtest(slower, len(pair)).pvalue)
                p_wil = float(wilcoxon(pair[hi_e], pair[lo_e]).pvalue)
                # fmt_p already carries its own '<' for tiny values, so the label must
                # not add a second one.
                def _p(v):
                    s = su.fmt_p(v)
                    return f"p{s}" if s.startswith("<") else f"p={s}"
                tests = f" (sign test {_p(p_sign)}, Wilcoxon {_p(p_wil)})"
            except Exception:
                tests = ""
            # How much of the total is carried by its largest contributors: if dropping
            # a handful of complexes flips the sign, the total is not describing the arm.
            k = 3
            flipped = float(delta.drop(delta.nsmallest(k).index).sum()) / 3600.0
            A(f"     exh{lo_e} → exh{hi_e}: total {d:+.3f} wall-h, yet {slower} of "
              f"{len(pair)} complexes are")
            A(f"     individually SLOWER at exh{hi_e}{tests} and the paired per-complex "
              f"median ratio")
            A(f"     is {med:.3f}. Dropping the {k} largest single contributors turns "
              f"the total to {flipped:+.3f} wall-h.")
            A(f"     The sum inverts; the location does not. A few complexes that ran "
              f"anomalously slowly")
            A(f"     in the exh{lo_e} batch carry it.")
        A("")
        A("     No return-per-hour is computed across an inverted segment. The cost "
          "delta is not a")
        A("     cost, and dividing by it would produce a confident number with no "
          "referent. This is")
        A("     the batch confound above, showing up as an arithmetic impossibility "
          "rather than as")
        A("     a subtle bias — the same algorithm cannot get cheaper by searching "
          "harder.")
        A("")
    for line in textwrap.wrap(CORE_SECOND_WARNING, 84):
        A("   " + line)
    A("")
    A("   Wall time also grows far more slowly than exhaustiveness"
      + (f" (log-log slope {eff['log_log_slope']:.2f} against 1.0 for proportional)"
         if eff.get("ok") else "") + ",")
    A("   because the Monte-Carlo runs are spread over the thread pool and a large "
      "part of each")
    A("   whole-protein complex's cost is exhaustiveness-independent setup. Wall time "
      "at a fixed")
    A("   allocation therefore FLATTERS high exhaustiveness relative to the compute it "
      "consumes.")
    A("   Exhaustiveness is never used as a proxy for cost here, in either direction.")
    A("")

    A("3. YIELD — complexes with at least one qualifying pose")
    A("")
    for depth in headline_depths:
        A(f"   ranking depth: {_depth_label(depth)}"
          + ("   (what a user actually consumes)" if depth == 1
             else "   (what the search could deliver with a perfect re-ranker)"))
        rows = []
        for e in arms_exh:
            r = [f"exh{e}"]
            for gate in GATES:
                y = yield_df[(yield_df.exhaustiveness == e) & (yield_df.gate == gate)
                             & (yield_df.depth == depth)].iloc[0]
                r.append(f"{int(y.k_complexes):3d}  {100 * y.rate:5.1f}%")
            rows.append(r)
        A("   " + _fw(rows, ["arm"] + [GATE_SHORT[g] for g in GATES],
                      "<" + ">" * len(GATES)).replace("\n", "\n   "))
        A("")

    A("4. SATURATION — which gates can answer the question at all")
    A("")
    rows = []
    for r in sat_df[sat_df.depth.isin(headline_depths)].itertuples():
        lo = getattr(r, "tost_5pp_lo_pp", np.nan)
        hi = getattr(r, "tost_5pp_hi_pp", np.nan)
        ok = getattr(r, "tost_5pp_equivalent", None)
        bound = (f"[{lo:+.2f}, {hi:+.2f}] pp"
                 + ("  equivalent at 5 pp" if ok else "  NOT equivalent at 5 pp")
                 if lo == lo else "—")
        rows.append([GATE_SHORT[r.gate], _depth_label(r.depth), r.flag,
                     f"{r.observable_range_pp:.2f}", str(int(r.n_discordant_total)),
                     f"{r.mde_design_pp:.2f}" if r.mde_design_pp == r.mde_design_pp else "—",
                     su.fmt_p(r.cochran_p) if r.cochran_p == r.cochran_p else "—",
                     bound])
    A("   " + _fw(rows, ["gate", "depth", "state", "range pp", "discordant",
                         "MDE pp", "Cochran Q p", "first-vs-last bound"],
                  "<<<>>>><").replace("\n", "\n   "))
    A("")
    A("   DEGENERATE means not one complex on the whole ladder changes its answer, so "
      "there is")
    A("   no test to run and the honest output is a bound. A non-significant McNemar "
      "would bound")
    A("   nothing by itself. The state is computed per gate AND depth: the same gate "
      "can be")
    A("   degenerate at one depth and live at another, so it is never inferred from a "
      "sampled")
    A("   depth. The minimum detectable difference is fixed from the LARGEST "
      "discordance on the")
    A("   ladder, not from each contrast's own — using a contrast's own discordance "
      "would certify")
    A("   the design as most sensitive exactly where it is blind.")
    A("")
    A("   PoseBusters validity is flat here because it is saturated by a hydrogen-"
      "handling")
    A("   artefact in the validation path, not because exhaustiveness repairs or "
      "breaks poses.")
    # Stated from the data rather than hardcoded: whether any form step survives
    # correction depends on which rungs are on the ladder.
    kb = contrasts[(contrasts.gate == "kabsch1") & (contrasts.p_adjusted < 0.05)]
    A("   The form gate (Kabsch) improves across the sweep as a whole but flattens near "
      "the top.")
    if kb.empty:
        A("   No single rung-to-rung form step is individually resolvable after "
          "correction.")
    else:
        which = ", ".join(f"{r.segment} at {_depth_label(r.depth)}"
                          for r in kb.sort_values("p_adjusted").head(4).itertuples())
        A(f"   {len(kb)} of {int((contrasts.gate == 'kabsch1').sum())} form steps "
          f"survive correction ({which}"
          + (", ..." if len(kb) > 4 else "") + ");")
        A("   the rest are not individually resolvable.")
    A("")

    A("5. MARGINAL RETURN PER HOUR")
    A("")
    # Refusals are frequent and verbose, so tables carry a short code and each table
    # prints only the legend entries it actually used.
    refusal_codes: dict[str, str] = {}
    used: set[str] = set()

    def _code(status: str) -> str:
        if status == "ok":
            return ""
        key = status.split("—")[0].strip()
        if key not in refusal_codes:
            refusal_codes[key] = f"R{len(refusal_codes) + 1}"
        used.add(key)
        return refusal_codes[key]

    def _legend() -> None:
        if not used:
            return
        A("")
        for text in sorted(used, key=lambda t: refusal_codes[t]):
            A(f"     {refusal_codes[text]} = {text}")
        used.clear()

    rows = []
    for r in mr_df[mr_df.is_headline].itertuples():
        rate = f"{r.rate_per_h:+.2f}" if r.rate_per_h == r.rate_per_h else _code(r.rate_status)
        ci = (f"[{r.rate_lo:+.2f}, {r.rate_hi:+.2f}]"
              if r.rate_lo == r.rate_lo else "—")
        recip = (f"{60 * r.wall_h_per_extra_complex:.0f} min"
                 if r.wall_h_per_extra_complex == r.wall_h_per_extra_complex
                 else _code(r.reciprocal_status))
        rows.append([GATE_SHORT[r.gate], _depth_label(r.depth), r.segment,
                     f"{r.delta_yield:+.0f}", f"{r.gained}/{r.lost}",
                     f"{r.delta_cost_h:+.3f}", rate, ci, recip])
    A("   " + _fw(rows, ["gate", "depth", "step", "Δcplx", "gain/loss", "Δh",
                         "cplx/h", "95% CI", "per extra cplx"],
                  "<<<>>>>>>").replace("\n", "\n   "))
    _legend()
    A("")
    A("   Δcplx is a NET difference; 'gain/loss' shows the two directions separately, "
      "because")
    A("   at rank-1 the top step both gains and loses complexes and they partly cancel.")
    A("   The interval is a paired complex-level bootstrap that recomputes the cost on "
      "the SAME")
    A("   resampled complexes — holding it fixed would understate the variance and "
      "break the")
    A("   pairing, since the complexes extra sampling rescues are the slow ones.")
    A("")
    A("   'per extra cplx' is the wall time bought per additional recovered complex. It "
      "refuses")
    A("   rather than printing a number when the yield delta it divides by is not "
      "reliably")
    A("   positive. Where it does print, it is a plug-in quotient (Δ wall-h / Δcplx) "
      "carrying NO")
    A("   interval: its denominator has mass near zero, so a percentile interval on it "
      "would not")
    A("   be an interval. The bootstrap decides whether the cell prints, not what it "
      "prints, so")
    A("   conservative claims should rest on complexes-per-hour rather than on its "
      "reciprocal.")
    A("")

    A("6. DO RETURNS DIMINISH?")
    A("")
    if dim_df.empty:
        A("   Not testable: fewer than three rungs with timing.")
    else:
        rows = []
        for r in dim_df[dim_df.is_headline].itertuples():
            rows.append([GATE_SHORT[r.gate], _depth_label(r.depth),
                         f"{r.rate_first_per_h:+.2f}", f"{r.rate_last_per_h:+.2f}",
                         f"{r.contrast:+.2f}" if r.contrast == r.contrast else "—",
                         (f"[{r.lo:+.2f}, {r.hi:+.2f}]" if r.lo == r.lo
                          else _code(r.status)),
                         "YES" if r.diminishing else "not resolved",
                         "confirmatory" if r.inference_claimable else "descriptive"])
        A("   " + _fw(rows, ["gate", "depth", "first cplx/h", "last cplx/h", "contrast",
                             "95% CI", "diminishing", "family"],
                      "<<>>>><<").replace("\n", "\n   "))
        _legend()
        A("")
        first_lbl = dim_df.iloc[0]["first_segment"].replace("->", " → ")
        last_lbl = dim_df.iloc[0]["last_segment"].replace("->", " → ")
        A(f"   The contrast is (return per hour on {first_lbl}) minus (return per hour "
          f"on {last_lbl}),")
        A("   both recomputed on the SAME resampled complexes. A CI above zero is the "
          "diminishing-")
        A("   returns claim.")
        shares = (dim_df.iloc[0]["first_segment"].split("->")[1]
                  == dim_df.iloc[0]["last_segment"].split("->")[0])
        if shares:
            A("   The two segments abut, so they share a rung that enters the two rates "
              "with opposite")
            A("   signs and makes the contrast conservative rather than liberal.")
        else:
            A("   The two segments are disjoint on this ladder, so the pairing is at "
              "the COMPLEX level")
            A("   only and the interval is effectively the independent-segments "
              "interval. No shared-rung")
            A("   term widens it.")
        A("")
        A("   No minimum detectable difference, sample-size argument or type-I "
          "simulation exists")
        A("   for this statistic. Its SIGN and interval are the claim; its MAGNITUDE "
          "moves with the")
        A("   cost basis and should never be quoted without naming that basis.")
        A("")

    if not sweep.empty:
        A("7. WHERE EACH RUNG PAYS — yield gain against ranking depth")
        A("")
        rows = []
        for gate in HEADLINE_GATES:
            for seg, g in sweep[sweep.gate == gate].groupby("segment", sort=False):
                g = g.sort_values("depth")
                pos = g[g.delta_complexes > 0]
                sig = g[(g.p_raw < 0.05) & (g.delta_complexes > 0)]
                rows.append([GATE_SHORT[gate], seg,
                             f"{int(g[g.depth == 1].delta_complexes.iloc[0]):+d}",
                             str(int(pos.depth.min())) if len(pos) else "never",
                             str(int(sig.depth.min())) if len(sig) else "never",
                             f"{int(g[g.depth == max_rank].delta_complexes.iloc[0]):+d}"])
        A("   " + _fw(rows, ["gate", "step", "Δ at rank-1", "first depth positive",
                             "first depth p<0.05", f"Δ at depth {max_rank}"],
                      "<<>>>>").replace("\n", "\n   "))
        A("")
        A("   This is the number a reader needs before choosing an exhaustiveness. A "
          "rung's verdict")
        A("   depends on how many poses get inspected: extra sampling finds poses that "
          "Vina's own")
        A("   scoring function then fails to promote to rank 1, so a step can be worth "
          "nothing at")
        A("   rank-1 and clearly worth paying for a few ranks down.")
        A("")
        A("   This section is a YIELD statement only. It carries no cost, so it stands "
          "whether or")
        A("   not a segment's cost delta was measurable.")
        A("")

    knee_rows = mr_df[(mr_df.gate == knee_gate) & (mr_df.depth == knee_depth)]
    A("8. KNEE")
    A("")
    if knee_rows.empty:
        A("   Not evaluated: the pre-registered cell is not in this run.")
    else:
        # A rung only counts as reached if its incoming step both had a usable cost and
        # a return interval clear of zero. A step whose cost could not be measured
        # cannot promote the knee.
        knee = None
        for r in knee_rows.sort_values("exh_hi").itertuples():
            # rate_status is "ok" only when the cost delta was measurable, so this is
            # exactly the cost half of the invariant and survives any rewording of the
            # refusal message. reciprocal_status would wrongly also catch YIELD refusals.
            if r.rate_status != "ok":
                continue
            if r.rate_lo == r.rate_lo and r.rate_lo > 0:
                knee = r.exh_hi
        A(f"   Pre-registered cell: {GATE_SHORT[knee_gate]} at "
          f"{_depth_label(knee_depth)}.")
        if args.knee_cell_justification:
            A(f"   Cell overridden. Justification: {args.knee_cell_justification}")
        A(f"   Knee = exhaustiveness {knee} — the largest rung whose incoming "
          f"marginal-return" if knee else
          "   No rung's incoming marginal-return interval excludes zero in this cell.")
        if knee:
            A("   interval excludes zero, in that cell.")
        A("")
        A("   Claimed in ONE cell only, fixed before looking. Taking the largest rung "
          "whose")
        A("   interval clears zero across many gate x depth cells is a maximum over a "
          "family of")
        A("   uncorrected intervals and places a knee at the top rung far more often "
          "than its")
        A("   nominal rate. Every other cell in this report is descriptive.")
        A("")
        A("   The knee is bounded above by the largest arm in the ladder. "
          "'Exhaustiveness beyond")
        A(f"   {arms_exh[-1]} does not help' is NOT supported by anything here; only "
          f"'no benefit was")
        A("   detectable between the top two rungs at this depth' is.")
        A("")
        A("   Budget crossing, for a reader with their own threshold:")
        rows = [[r.segment, f"{r.delta_yield:+.0f}", f"{r.delta_cost_h:+.3f}",
                 (f"{60 * r.wall_h_per_extra_complex:.0f} min"
                  if r.wall_h_per_extra_complex == r.wall_h_per_extra_complex
                  else _code(r.reciprocal_status))]
                for r in knee_rows.itertuples()]
        A("   " + _fw(rows, ["step", "Δcplx", "Δ wall h", "wall time per extra cplx"],
                      "<>>>").replace("\n", "\n   "))
        _legend()
        A("")

    if not strand_df.empty:
        A("9. RAW VERSUS gnina RE-RANKED, at equal exhaustiveness")
        A("")
        A("   A second STRAND, not a further rung. The rescored arm shares the raw arm's")
        A("   poses, receptor and ligands — re-ranking reorders a fixed pool, so a rank-1")
        A("   gain is the CNN promoting a pose the Vina score buried, never a pose the")
        A("   search did not find.")
        A("")
        for depth in headline_depths:
            A(f"   ranking depth: {_depth_label(depth)}")
            rows = []
            for gate in HEADLINE_GATES:
                for r in strand_df[(strand_df.gate == gate)
                                   & (strand_df.depth == depth)].itertuples():
                    rows.append([GATE_SHORT[gate], f"exh{int(r.exhaustiveness)}",
                                 f"{r.k_raw:3d} → {r.k_rescored:3d}",
                                 f"{r.delta_complexes:+d}",
                                 f"{r.gained}/{r.lost}", f"{r.delta_pp:+.2f}",
                                 f"[{r.ci_lo_pp:+.2f}, {r.ci_hi_pp:+.2f}]",
                                 su.fmt_p(r.p_raw),
                                 su.fmt_p(r.p_adjusted) if r.p_adjusted == r.p_adjusted
                                 else "—"])
            A("   " + _fw(rows, ["gate", "rung", "raw → rescored", "Δ", "gain/loss",
                                 "Δ pp", "95% CI", "p", "Holm"],
                          "<<>>>>>>>").replace("\n", "\n   "))
            A("")
        A("   COST OF THE RE-RANKING PASS (compute hours, summed per-pose optimiser time,")
        A("   restricted to the same complexes as every other number in this report)")
        rows = []
        for e, oa in sorted(opt_by_exh.items()):
            rows.append([f"exh{e}", f"{oa['optimizer_compute_h']:.2f}",
                         f"{oa['optimize_workers']}",
                         str(oa.get("optimize_cpu") or "—"),
                         str(oa.get("gpu_max_concurrent") or "—"),
                         f"{oa['optimizer_rows']:,}",
                         ", ".join(oa["optimizer_tools"]) or "—"])
        A("   " + _fw(rows, ["rung", "compute h", "workers", "cpu", "gpu slots",
                             "poses", "tool"], "<>>>>><").replace("\n", "\n   "))
        A("")
        resourcing = {(oa["optimize_workers"], oa.get("optimize_cpu"))
                      for oa in opt_by_exh.values()}
        A("   Compute hours, NOT wall clock: these passes ran at different worker counts,")
        A("   so their wall times differ by more than an order of magnitude for reasons")
        A("   that have nothing to do with the method.")
        if len(resourcing) > 1:
            A("")
            A("   BUT COMPUTE IS NOT COMPARABLE ACROSS THESE ROWS EITHER. The resourcing")
            A("   changed with the rung — see the workers / cpu / gpu-slot columns above —")
            A("   and on this campaign that change moves measured per-pose optimiser time")
            A("   by far more than the rung does. Under heavy GPU sharing the per-pose")
            A("   number also stops tracking the pose at all and starts tracking the")
            A("   scheduler, because every process waits its turn on a slot. READ EACH ROW")
            A("   AS THE COST OF THAT PASS AS IT WAS RUN. This column is NOT a")
            A("   cost-versus-exhaustiveness curve and no trend should be read across it.")
        A("")
        A("   It is also NOT comparable with the search hours in section 2 — that is")
        A("   32-thread Vina wall time on CPU, this is summed optimiser time on a GPU.")
        A("")
        A("   The re-ranking is NOT folded into the marginal-return ladder. Search effort")
        A("   and re-ranking are different kinds of work on different hardware; a single")
        A("   cost axis spanning both would compare quantities that are not commensurable.")
        A("")

    A("10. WHAT THIS ANALYSIS DOES NOT CLAIM")
    A("")
    refusals = [
        "No saturation ceiling, and no 'compute needed to reach 90% of maximum yield'. "
        "A saturating fit on this few rungs returns a ceiling whose interval spans most "
        "of the feasible range and an implied cost beyond the most expensive arm ever "
        "measured. That is extrapolation, not a result.",
        f"No claim that exhaustiveness above {arms_exh[-1]} is useless. The ladder "
        "simply stops there.",
        "No knee outside the one pre-registered cell, and no elbow from a curvature "
        "method. With at most two interior rungs, a curvature elbow returns the middle "
        "rung by construction.",
        "No 'no effect' conclusions from a non-significant test. Flat gates are "
        "reported as equivalence bounds instead.",
        "Gate (c) alone is never a docking success rate. Best-fit deviation is measured "
        "after optimal superposition, so it says the ligand's internal geometry is "
        "right, not that the ligand is in the right place.",
        "No marginal cost in core-hours. On this ladder that currency is non-monotone "
        "and its first step is negative.",
        "No modelled cost is ever substituted for a measured one, no run-time window is "
        "excluded, and no arm is called inflated on the strength of a fit it was held "
        "out of.",
        "No separation of search effort from batch. Each arm ran once, on its own day, "
        "and the machine-state spread is as large as the effect being measured.",
        "No single cost axis spanning both strands. Search hours are 32-thread Vina wall "
        "time on CPU; re-ranking hours are summed single-process optimiser time on GPU. "
        "They are not commensurable and are never added, divided or plotted together.",
        "No claim that re-ranking cost rises with exhaustiveness. The optimiser passes "
        "were run at different worker, CPU and GPU-slot settings, and that difference "
        "moves per-pose time by more than the rung difference it is confounded with.",
        "No claim that a rescored rung is a point on the exhaustiveness ladder. It shares "
        "its raw rung's poses, so it can only reorder them — at oracle depth the two "
        "strands differ solely through the optimiser's minimisation.",
    ]
    for i, claim in enumerate(refusals, 1):
        A(f"   ({i}) {claim}")
    A("")

    report = "\n".join(L)
    (out_dir / "exh_returns_report.txt").write_text(report + "\n")
    print()
    print(report)

    for fig_name, frames in (
            ("exh_yield_vs_cost",
             {"yield (headline depths)": yield_df[yield_df.is_headline_depth],
              "cost per arm": cost_tbl,
              "saturation": sat_df[sat_df.depth.isin(headline_depths)]}),
            ("exh_marginal_return",
             {"marginal return": mr_df[mr_df.is_headline],
              "diminishing returns": dim_df,
              "yield contrasts (family A, Holm)":
                  contrasts[contrasts.family == "A-confirmatory"]}),
            ("exh_depth_sweep", {"depth sweep": sweep}),
            ("exh_raw_vs_rescored", {"raw vs rescored": strand_df})):
        if not (fig_dir / f"{fig_name}.png").exists():
            continue
        blocks = [f"{fig_name}.png — statistics sidecar", "#" * 78,
                  f"generated by Scripts/Analysis/{Path(__file__).name}",
                  f"source: {per_pose_path}"]
        if announce_conv:
            blocks.append(f"reference convention: {conv_note}")
        blocks += ["",
                  "Depths are nested and monotone — a complex passing at depth 1 passes "
                  "at every deeper",
                  "depth — so the exploratory family's members are strongly positively "
                  "dependent. BH is",
                  "valid under positive dependence, but the effective number of "
                  "independent tests is far",
                  "below the nominal family size. Cochran's Q is reported uncorrected "
                  "as a gatekeeper;",
                  "no alpha is allocated across the hierarchy.", ""]
        for name, frame in frames.items():
            blocks.append(f"[{name}]  n_rows={len(frame)}")
            blocks.append(frame.to_string(index=False) if len(frame) else "  (empty)")
            blocks.append("")
        (fig_dir / f"{fig_name}_stats.txt").write_text("\n".join(blocks) + "\n")

    manifest = {
        "n_complexes": n_cplx, "max_rank": max_rank,
        "oracle_equals_top30": oracle_is_top30,
        "arms": [{k: (str(v) if isinstance(v, Path) else v) for k, v in a.items()}
                 for a in arms],
        "pending_arms": [a["method_key"] for a in pending],
        "optimizer_arms": [{"method_key": a["method_key"],
                            "exhaustiveness": a["exhaustiveness"],
                            "optimizer_compute_h": a["optimizer_compute_h"],
                            "optimize_workers": a["optimize_workers"],
                            "optimize_cpu": a.get("optimize_cpu"),
                            "gpu_max_concurrent": a.get("gpu_max_concurrent"),
                            "optimizer_rows": a["optimizer_rows"]} for a in opt_arms],
        "cost_arms": cost_arms, "headline_depths": headline_depths,
        "modal_cpus_per_worker": cost_meta["modal_cpus"],
        "arm_effects_pct": cost_meta["arm_effects"].get("arm_effect_pct"),
        "arm_effect_spread_pp": cost_meta["arm_effects"].get("spread_pp"),
        "knee_cell": [knee_gate, knee_depth],
        "knee_cell_justification": args.knee_cell_justification,
        "input_parity": parity, "n_boot": args.n_boot, "seed": args.seed,
    }
    if announce_conv:
        manifest["reference_convention"] = convention
        manifest["reference_convention_declared_by_table"] = conv_declared
        manifest["cascade_pins_source"] = pins_source
        manifest["cascade_pins_checked"] = checked
    (out_dir / "exh_run_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str) + "\n")

    print()
    for f in sorted(out_dir.glob("*.csv")) + sorted(out_dir.glob("*.txt")) \
            + sorted(out_dir.glob("*.json")) + sorted(fig_dir.glob("*")):
        print(f"wrote {f.relative_to(root) if f.is_relative_to(root) else f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
