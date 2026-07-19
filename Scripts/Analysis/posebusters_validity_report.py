"""Visualize PoseBusters validity per docking tool & receptor-ligand pair.

Reads posebusters_filtered_results.csv and produces:
    1. Bar chart: total / valid pose counts per docking method.
    2. Heatmap: valid-pose counts per receptor-ligand x method — 02 easiest (top-N) and
       02b hardest (bottom-N) complexes.
    3. Cleveland dot plot: top-N receptor-ligand pairs comparing valid poses per method.
    4. Boxplot: distribution of per-pair valid poses per method.
    5. Heatmap: per-check pass rate per method (column labels carry each variant's n).
    6/7. Bar chart: PB-validity per variant, one figure per tool (EquiBind / DiffDock);
         7b: the same as a base-method × optimizer validity matrix (heatmap).
    8. (--rmsd-sweep-max) Sweep: PB-valid poses vs RMSD-to-crystal cutoff, as two figures
       — 8a cumulative count and 8b yield (% of the method's poses).
    9. Waterfall: per-check pose attrition for the top variant of each tool
       (AutoDock / best DiffDock / best EquiBind, best = highest PB-valid fraction).
    + CSV summaries next to the plots: summary_per_tool.csv (per method:
      total/valid/valid_fraction, plus rmsd2_* over the ≤ 2 Å near-native subset for
      crystal sets) and check_failure_rates_per_variant.csv (per variant: n_poses,
      survival %, and per-check failure %).

A pose is considered "physically valid" when ALL of PoseBusters' canonical
pass/fail tests succeed — the same 20-test PB-Valid criterion used by
run_posebusters.py (imported from there so the two never drift). This includes
the cofactor / metal / water clash checks, not just protein-side checks.

EquiBind variant split
----------------------
EquiBind emits several flavours of pose per complex, all tagged "equibind_guided"
in the raw CSV. This report unpacks every variant axis so the flavours can be
compared directly (toggle with --no-split-equibind). It reads the provenance
columns the PoseBusters CSV now carries (``pocket_source`` / ``clamp_variant`` /
``refine_variant``, written from each pose's SDF tags) and falls back to parsing
the pose filename for older CSVs that lack them:

  pocket axis  (pocket_source)     unguided   -> blind EquiBind (no pocket)
                                   fpocket    -> guided by an fpocket pocket
                                   p2rank     -> guided by a p2rank pocket
  re-search    (refine_variant)    smina (__refSMINA) | gnina (__refGNINA) | raw (__refRAW) | None
  clamp axis   (clamp_variant)     clampON (__clampON) | clampOFF (__clampOFF) | None

So e.g. ``fpocket_raw_p001_pose01__clampOFF__refSMINA.sdf`` becomes the method
``equibind_fpocket_smina_clampOFF``. fpocket and p2rank are reported as distinct
methods (not merged into one "guided" bucket). Axes that aren't present (the
default single-variant runs) simply collapse to ``equibind_fpocket`` /
``equibind_p2rank`` / ``equibind_unguided``.
NOTE: refine_mode='on' rewrites the pose in place with no suffix/tag, so an
unsuffixed pose is reported under its base label (the smina split is only
visible for refine_mode='both' runs, which emit __refRAW/__refSMINA/__refGNINA).

Quick usage:
    # Run with defaults
    python Scripts/Analysis/posebusters_validity_report.py

    # Run with custom input/output locations
    python Scripts/Analysis/posebusters_validity_report.py \
        --csv posebusters_results/benchmark/dock/posebusters_filtered_results.csv \
        --out-dir posebusters_results/benchmark/dock/validity_report

    # Keep EquiBind as a single "equibind_guided" bucket (old behaviour)
    python Scripts/Analysis/posebusters_validity_report.py --no-split-equibind

    # Show only the single best EquiBind variant (highest PB-Valid AND RMSD ≤ 2 Å =
    # oracle_pb_valid_and_rmsd2_%, read from --oracle-summary), relabelled "EquiBind*". For
    # crystal-free sets (Orai) the default --oracle-summary borrows the benchmark ranking.
    python Scripts/Analysis/posebusters_validity_report.py --best-equibind-only

    # Keep only poses within 5 Å of the crystal (RMSD joined from per_pose_metrics.csv,
    # written by posebusters_pose_comparison.py — run that first). Every CSV/figure
    # then reflects the filtered set, and the figures are titled accordingly. Under
    # --max-rmsd the per-method summary reports validity over both populations side by
    # side: generated_poses / generated_valid_poses / generated_valid_fraction (ALL
    # poses produced) and within_rmsd_poses / valid_poses / valid_fraction (the kept
    # subset within the cutoff).
    python Scripts/Analysis/posebusters_validity_report.py --max-rmsd 5

    # Keep every variant in the CSV tables / printed summary but collapse each tool to
    # a single variant ('DiffDock*'/'EquiBind*') in the main comparison figures (01-05
    # + the RMSD sweep). Also emits a per-tool variant-comparison figure (06 EquiBind,
    # 07 DiffDock). Pair with a selector to choose the kept variant, or let it default
    # to each tool's best on PB-Valid AND RMSD ≤ 2 Å.
    python Scripts/Analysis/posebusters_validity_report.py --collapse-plots-only
    python Scripts/Analysis/posebusters_validity_report.py \
        --collapse-plots-only --diffdock-variant smina --best-equibind-only

    # Sweep figure (08): how the PB-valid pose count grows as the RMSD-to-crystal cutoff
    # is relaxed from 0 to 10 A (crystal sets only; joined from --per-pose-metrics).
    python Scripts/Analysis/posebusters_validity_report.py --rmsd-sweep-max 10

    # A per-variant PB check failure-rate table is written on every run
    # (check_failure_rates_per_variant.csv): rows = PB checks, columns = variants,
    # values = %% of that variant's poses that fail the check (worst checks first).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Shared, unit-tested statistical helpers (paired McNemar/Cochran's Q, Friedman +
# Kendall's W, Wilson CI, Holm/BH). Never reimplement a test — import them.
import os as _os, sys as _sys  # noqa: E402
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import stats_utils as su  # noqa: E402

# Single source of truth: reuse the pipeline's canonical PoseBusters test set and
# its bool-coercion so this report and run_posebusters.py agree exactly on what
# counts as a valid pose. (run_posebusters lives one package over; add it to the
# path and import the constant + coercer rather than duplicating a check list.)
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PB_DIR = _PROJECT_ROOT / "Scripts" / "Docking" / "Posebusters"
for _p in (str(_PROJECT_ROOT), str(_PB_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)
from run_posebusters import (  # noqa: E402
    CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool,
)

# Boolean columns; True == pass. The canonical 20-test PB-Valid set (the full
# criterion, incl. cofactor/metal/water clash checks). Only the columns actually
# present in a given CSV are used (a 'mol'-mode CSV omits the protein/cofactor
# columns), so the report degrades gracefully instead of erroring.
CRITICAL_CHECKS: list[str] = list(CANONICAL_TEST_COLUMNS)

# ───────────────────── method labelling / ordering ──────────────────────
# The grouping key is the ``docking_method`` column. For EquiBind it is
# refined into per-variant labels (see _classify_equibind / _apply_equibind_split).
# Everything below is data-driven so any subset of variants renders correctly.

# Token orderings used to sort EquiBind variants consistently.
# "guided" is kept for back-compat with older, unsplit CSVs.
_POCKET_ORDER = {"unguided": 0, "fpocket": 1, "p2rank": 2, "guided": 3}
_REFINE_ORDER = {None: 0, "raw": 1, "smina": 2, "gnina": 3}
_CLAMP_ORDER = {None: 0, "clampON": 1, "clampOFF": 2}

# Green family for EquiBind variants (cycled if more than this many appear).
_EQ_PALETTE = ["#2ca02c", "#74c476", "#1b7837", "#a6dba0",
               "#006d2c", "#5aae61", "#00441b", "#c7e9c0"]
# Orange family for DiffDock smina/gnina optimizer variants (raw uses the base
# orange in _BASE_COLORS; these are the lighter/darker shades for the variants).
_DD_PALETTE = ["#ffbb78", "#d95f02", "#fdae6b", "#a63603"]
_BASE_COLORS = {"autodock": "#1f77b4", "diffdock": "#ff7f0e"}

# Per-run display-label overrides (method key -> label). Populated in main() when
# --best-equibind-only is active, where the single retained EquiBind variant is
# shown as "EquiBind*". Honoured by _pretty_method, so every plot/axis label picks
# it up with no further changes.
_LABEL_OVERRIDES: dict[str, str] = {}

# Appended to every figure title when --max-rmsd is active, so a filtered report's
# figures can never be mistaken for the unfiltered one. Set once in main().
_RMSD_FILTER_NOTE: str = ""


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
    """Human-readable label for a method / EquiBind- or DiffDock-variant key."""
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


def _method_sort_key(m: str):
    if m == "autodock":
        return (0, 0, 0, 0)
    if m.startswith("diffdock"):
        # raw DiffDock first, then smina-/gnina-optimised variants.
        return (1, {"diffdock": 0, "diffdock_smina": 1, "diffdock_gnina": 2}.get(m, 3), 0, 0)
    if m.startswith("equibind"):
        p, r, c = _eq_tokens(m)
        return (2, _POCKET_ORDER.get(p, 9), _REFINE_ORDER.get(r, 0),
                _CLAMP_ORDER.get(c, 0))
    return (3, 0, 0, 0)


def _method_order(df: pd.DataFrame) -> list[str]:
    """All methods present, in a stable, readable order."""
    methods = [str(m) for m in pd.unique(df["docking_method"])]
    return sorted(methods, key=lambda m: (_method_sort_key(m), m))


def _method_colors(order: list[str]) -> dict[str, object]:
    colors: dict[str, object] = {}
    eq = [m for m in order if m.startswith("equibind")]
    # DiffDock optimizer variants share an orange family (base orange = raw).
    dd = [m for m in order if m.startswith("diffdock") and m not in _BASE_COLORS]
    leftover = plt.cm.tab10.colors
    j = 0
    for m in order:
        if m in _BASE_COLORS:
            colors[m] = _BASE_COLORS[m]
        elif m in eq:
            colors[m] = _EQ_PALETTE[eq.index(m) % len(_EQ_PALETTE)]
        elif m in dd:
            colors[m] = _DD_PALETTE[dd.index(m) % len(_DD_PALETTE)]
        else:
            colors[m] = leftover[j % len(leftover)]
            j += 1
    return colors


def _present_checks(df: pd.DataFrame) -> list[str]:
    """Canonical checks actually present in *df* (in canonical order)."""
    return [c for c in CRITICAL_CHECKS if c in df.columns]


def _bool_checks(df: pd.DataFrame, checks: list[str]) -> pd.DataFrame:
    """Bool DataFrame for *checks*, coerced identically to run_posebusters.

    Uses the pipeline's ``coerce_test_cols_to_bool`` so mixed bool/float/string
    PB outputs map the same way here as in the verdict — in particular an absent
    cofactor/water check (NaN) is treated as a pass (not-applicable), matching
    the pipeline rather than the old NaN->False behaviour.
    """
    sub = df[checks].copy()
    coerce_test_cols_to_bool(sub, checks)
    return sub


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
    """Refine ``docking_method`` for EquiBind rows into per-variant labels and
    add the parsed ``eq_pocket`` / ``eq_refine`` / ``eq_clamp`` axis columns."""
    df["eq_pocket"] = pd.NA
    df["eq_refine"] = pd.NA
    df["eq_clamp"] = pd.NA

    is_eq = df["docking_method"].str.startswith("equibind")
    if not is_eq.any():
        return df
    if "pose_name" not in df.columns:
        print("[warn] 'pose_name' column absent — cannot split EquiBind variants; "
              "leaving 'equibind_guided' as a single bucket.")
        return df

    eq_idx = df.index[is_eq]
    cls = [_classify_equibind(df.loc[i]) for i in eq_idx]
    labels = []
    for pocket, refine, clamp in cls:
        lab = f"equibind_{pocket}"
        if refine:
            lab += f"_{refine}"
        if clamp:
            lab += f"_{clamp}"
        labels.append(lab)

    df.loc[eq_idx, "docking_method"] = labels
    df.loc[eq_idx, "eq_pocket"] = [c[0] for c in cls]
    df.loc[eq_idx, "eq_refine"] = [c[1] for c in cls]
    df.loc[eq_idx, "eq_clamp"] = [c[2] for c in cls]
    return df


def _classify_diffdock(row) -> str | None:
    """Optimizer backend for one DiffDock pose: 'smina' | 'gnina' | None (original).

    Prefers the ``optimizer`` provenance column the PoseBusters CSV now carries
    (written by run_posebusters from each pose's ``optimized_<tool>/`` subfolder);
    falls back to the ``optimized_<tool>`` path component in ``pose_name`` for
    older CSVs.
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
    per-optimizer labels (``diffdock_smina`` / ``diffdock_gnina``) and add a
    ``dd_optimizer`` axis column; raw (un-optimised) DiffDock poses keep
    ``diffdock``. Mirrors _apply_equibind_split so optimised vs original poses
    are never pooled into a single DiffDock pass rate."""
    df["dd_optimizer"] = pd.NA
    is_dd = df["docking_method"].str.startswith("diffdock")
    if not is_dd.any():
        return df
    dd_idx = df.index[is_dd]
    opts = [_classify_diffdock(df.loc[i]) for i in dd_idx]
    df.loc[dd_idx, "docking_method"] = [f"diffdock_{o}" if o else "diffdock" for o in opts]
    df.loc[dd_idx, "dd_optimizer"] = [o or "original" for o in opts]
    return df


def _load_allowed_ids(path: Path) -> set[str]:
    """Load '<PDBID>_<LIG>' complex ids (one per line; '#' comments ignored)."""
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def load_and_score(csv_path: Path, split_equibind: bool = True) -> pd.DataFrame:
    df = pd.read_csv(csv_path, low_memory=False)
    checks = _present_checks(df)
    if not checks:
        raise ValueError(
            "No canonical PoseBusters check columns found in CSV: "
            f"{CRITICAL_CHECKS}"
        )
    missing = [c for c in CRITICAL_CHECKS if c not in df.columns]
    if missing:
        print(f"[warn] {len(missing)} canonical check(s) absent from CSV "
              f"(skipped — likely a 'mol'-mode run): {missing}")

    df["pb_valid"] = _bool_checks(df, checks).all(axis=1)
    df["pair"] = df["protein"].astype(str) + " / " + df["ligand"].astype(str)
    df["docking_method"] = df["docking_method"].astype(str).str.lower()
    df["method_raw"] = df["docking_method"]
    if split_equibind:
        df = _apply_equibind_split(df)
    # Always separate DiffDock optimizer variants — pooling optimised poses with
    # the raw DiffDock output would corrupt the per-method pass rate. No-op when
    # no optimised poses are present (every row stays ``diffdock``).
    df = _apply_diffdock_split(df)
    return df


# Default location of the oracle summary written by posebusters_pose_comparison.py
# (the only place the oracle_rmsd_le_2.0A_% metric exists). Datasets without their
# own crystal (e.g. Orai) borrow the benchmark ranking from here. Mirrors
# run_pandamap / pandamap_interaction_report DEFAULT_ORACLE_SUMMARY.
DEFAULT_ORACLE_SUMMARY = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/oracle_summary.csv")

# The PoseBusters filtered-results CSV carries no RMSD; the only place per-pose
# RMSD-to-crystal exists is per_pose_metrics.csv, written by
# posebusters_pose_comparison.py. --max-rmsd joins that column in from here (same
# dir as DEFAULT_ORACLE_SUMMARY). Crystal-free sets (Orai) have no RMSD, so
# --max-rmsd is not meaningful there and errors out rather than emptying the report.
DEFAULT_PER_POSE_METRICS = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")

# RMSD columns in per_pose_metrics.csv that --max-rmsd may filter on:
#   rmsd           — symmetry-corrected heavy-atom RMSD, NO superposition (docking
#                    accuracy; the column the oracle RMSD ≤ 2 Å metric uses). Default.
#   pb_rmsd        — PoseBusters' canonical check_rmsd (also no superposition).
#   pb_kabsch_rmsd — after optimal superposition (conformer similarity, not placement).
_RMSD_COLUMNS = ("rmsd", "pb_rmsd", "pb_kabsch_rmsd")


def apply_rmsd_filter(df: pd.DataFrame, metrics_csv: Path, max_rmsd: float,
                      rmsd_col: str, out_dir: Path) -> pd.DataFrame:
    """Keep only poses whose RMSD-to-crystal (``rmsd_col``) is ≤ ``max_rmsd`` Å.

    The PoseBusters filtered-results CSV has no RMSD, so the per-pose RMSD-to-crystal
    is joined in from *metrics_csv* (per_pose_metrics.csv, written by
    posebusters_pose_comparison.py) on the globally-unique (protein, ligand,
    pose_name) key — sidestepping any method-label differences between the two files.
    Poses with no RMSD value — absent from the metrics file or scored NaN — cannot be
    confirmed ≤ threshold and are dropped and counted, the same as poses above it.
    Writes a per-method rmsd_filter_summary.csv next to the plots.

    Exits with a clear message (rather than silently emptying the report) when the
    metrics file is missing, lacks ``rmsd_col``/the join keys, or matches no pose —
    e.g. a crystal-free set (Orai) where RMSD-to-crystal is undefined.
    """
    if not metrics_csv or not Path(metrics_csv).exists():
        sys.exit(f"[rmsd-filter] --max-rmsd {max_rmsd} requested but per-pose metrics "
                 f"not found at {metrics_csv}. Run posebusters_pose_comparison.py first "
                 f"(it writes per_pose_metrics.csv), or pass --per-pose-metrics.")

    metrics = pd.read_csv(metrics_csv, low_memory=False)
    if rmsd_col not in metrics.columns:
        sys.exit(f"[rmsd-filter] column '{rmsd_col}' not in {metrics_csv} "
                 f"(available: {[c for c in _RMSD_COLUMNS if c in metrics.columns]}).")
    key = ["protein", "ligand", "pose_name"]
    if any(k not in metrics.columns for k in key):
        sys.exit(f"[rmsd-filter] {metrics_csv} lacks join columns {key}.")

    m = metrics[key + [rmsd_col]].copy()
    for k in key:
        m[k] = m[k].astype(str)
    n_dup = int(m.duplicated(subset=key).sum())
    if n_dup:
        print(f"[rmsd-filter] warning: {n_dup} duplicate (protein,ligand,pose_name) "
              "rows in the metrics file — keeping the first of each.")
        m = m.drop_duplicates(subset=key)

    d = df.copy()
    for k in key:
        d[k] = d[k].astype(str)
    n_in = len(d)
    # left + validate='m:1' keeps every report row exactly once (never multiplies).
    d = d.merge(m.rename(columns={rmsd_col: "rmsd_to_crystal"}), on=key,
                how="left", validate="m:1")
    d["rmsd_to_crystal"] = pd.to_numeric(d["rmsd_to_crystal"], errors="coerce")

    d["_has_rmsd"] = d["rmsd_to_crystal"].notna()
    d["_kept"] = d["_has_rmsd"] & (d["rmsd_to_crystal"] <= max_rmsd)
    n_with = int(d["_has_rmsd"].sum())
    n_kept = int(d["_kept"].sum())
    n_no_rmsd = n_in - n_with
    n_over = n_with - n_kept

    if n_kept == 0:
        sys.exit(f"[rmsd-filter] no pose has {rmsd_col} ≤ {max_rmsd} Å joined from "
                 f"{metrics_csv} (matched {n_with}/{n_in}). Nothing to report — is this "
                 "the right metrics file, and a crystal-bearing set?")

    rep = (d.groupby("docking_method")
             .agg(total_in=("_kept", "size"),
                  with_rmsd=("_has_rmsd", "sum"),
                  kept=("_kept", "sum")))
    rep["dropped_no_rmsd"] = rep["total_in"] - rep["with_rmsd"]
    rep["dropped_over_thresh"] = rep["with_rmsd"] - rep["kept"]
    out_dir.mkdir(parents=True, exist_ok=True)
    rep.to_csv(out_dir / "rmsd_filter_summary.csv")

    print(f"[rmsd-filter] source: {metrics_csv} (column '{rmsd_col}', ≤ {max_rmsd:g} Å)")
    print(f"[rmsd-filter] poses in: {n_in:,} | with RMSD: {n_with:,} | kept: {n_kept:,} "
          f"| dropped >thr: {n_over:,} | dropped no-RMSD: {n_no_rmsd:,}")
    if n_no_rmsd:
        print(f"[rmsd-filter] note: {n_no_rmsd:,} pose(s) had no RMSD in the metrics "
              "file (unscored, or absent — e.g. DiffDock optimizer variants the "
              "comparison collapsed) and were dropped. See rmsd_filter_summary.csv.")

    return d[d["_kept"]].drop(columns=["_has_rmsd", "_kept"]).reset_index(drop=True)


def _resolve_variant_oracle(oracle_csv: Path) -> Path:
    """Prefer the ``*_all_variants.csv`` sibling of an oracle summary for per-variant ranking.

    posebusters_pose_comparison.py writes oracle_summary.csv with each tool COLLAPSED to a
    single row (``diffdock`` / ``equibind_unguided_gnina`` — its own best variant, relabelled
    to the bare tool name). This report, by contrast, splits DiffDock/EquiBind into per-variant
    labels (``diffdock_smina`` / ``diffdock_gnina`` / …). Ranking those split labels against the
    collapsed file matches only the bare ``diffdock`` row, so every other variant is dropped
    (``m in scores.index`` is False for it) and raw DiffDock "wins" by being the sole candidate.

    The comparison also writes an ``oracle_summary_all_variants.csv`` sibling that keeps one row
    per variant. When it exists we rank against that instead, so every variant is actually
    compared. Falls back to the given path (older runs / crystal-free sets that lack the sibling);
    a path already ending in ``_all_variants.csv`` is returned unchanged.
    """
    if not oracle_csv:
        return oracle_csv
    p = Path(oracle_csv)
    if p.name.endswith("_all_variants.csv"):
        return p
    sibling = p.with_name(f"{p.stem}_all_variants{p.suffix}")
    return sibling if sibling.exists() else p


# The metric the best-variant selectors rank on: the combined docking-success criterion
# PB-Valid AND RMSD ≤ 2 Å (a pose must be BOTH near-native AND physically valid), written
# per variant by posebusters_pose_comparison.py. This is what makes "best" mean a pose that
# is usable, not merely close to the crystal (raw DiffDock, say, is often near-native but
# clashing). Older summaries that predate the combined column fall back to the RMSD-only rate.
BEST_VARIANT_METRIC = "oracle_pb_valid_and_rmsd2_%"
BEST_VARIANT_METRIC_FALLBACK = "oracle_rmsd_le_2.0A_%"


def _variant_ranking_scores(oracle_csv: Path, tag: str) -> tuple[pd.Series | None, str | None]:
    """Per-variant ranking scores and the metric column used, read from *oracle_csv*.

    Ranks on ``oracle_pb_valid_and_rmsd2_%`` (PB-Valid AND RMSD ≤ 2 Å — the combined
    docking-success criterion), falling back to ``oracle_rmsd_le_2.0A_%`` when a summary
    predates the combined column. Returns (scores, metric) or (None, None) — with a ``[tag]``
    note — when the file carries neither column. The caller has already resolved *oracle_csv*
    to its ``*_all_variants.csv`` sibling (see _resolve_variant_oracle)."""
    osum = pd.read_csv(oracle_csv, index_col=0)
    col = (BEST_VARIANT_METRIC if BEST_VARIANT_METRIC in osum.columns
           else BEST_VARIANT_METRIC_FALLBACK if BEST_VARIANT_METRIC_FALLBACK in osum.columns
           else None)
    if col is None:
        print(f"  [{tag}] neither '{BEST_VARIANT_METRIC}' nor "
              f"'{BEST_VARIANT_METRIC_FALLBACK}' in {oracle_csv} — keeping all variants.")
        return None, None
    return pd.to_numeric(osum[col], errors="coerce"), col


def select_best_equibind(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind rows + only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest PB-Valid AND RMSD ≤ 2 Å rate
    (``oracle_pb_valid_and_rmsd2_%``, falling back to ``oracle_rmsd_le_2.0A_%`` for older
    summaries — see _variant_ranking_scores) in *oracle_csv*, the oracle_summary.csv written by
    posebusters_pose_comparison.py (its per-variant ``*_all_variants.csv`` sibling is preferred
    when present, see _resolve_variant_oracle, so the split variant labels are ranked rather than
    silently collapsed to the one ``equibind_*`` row the summary carries). For datasets without
    their own crystal (Orai),
    point --oracle-summary at the benchmark summary to borrow its ranking. The best
    variant is chosen among the variants actually present in *df*; AutoDock/DiffDock
    rows are always kept. Returns (filtered_df, best_variant_key); best is None and
    df unchanged when no EquiBind variant is present, the summary is missing/unreadable,
    or no present variant has a score.
    """
    methods = df["docking_method"].astype(str)
    eq_mask = methods.str.startswith("equibind")
    if not eq_mask.any():
        return df, None
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-equibind-only] oracle summary not found at {oracle_csv} — "
              "keeping all EquiBind variants.")
        return df, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-equibind-only")
    if scores is None:
        return df, None
    eq_present = sorted(methods[eq_mask].unique())
    cand = {m: float(scores[m]) for m in eq_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-equibind-only] none of the present EquiBind variants have a "
              f"score in {oracle_csv} — keeping all EquiBind variants.")
        return df, None

    best = max(cand, key=cand.get)
    keep = (~eq_mask) | (methods == best)
    return df[keep].reset_index(drop=True), best


def select_best_diffdock(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-DiffDock rows + only the single best DiffDock optimizer variant.

    "Best" = the DiffDock variant (diffdock / diffdock_smina / diffdock_gnina) with the
    highest PB-Valid AND RMSD ≤ 2 Å rate (``oracle_pb_valid_and_rmsd2_%``, falling back to
    ``oracle_rmsd_le_2.0A_%`` for older summaries — see _variant_ranking_scores) in *oracle_csv*
    (its per-variant ``*_all_variants.csv`` sibling is used when present — see
    _resolve_variant_oracle — so the split smina/gnina labels are actually ranked, not silently
    dropped, which would otherwise leave bare ``diffdock`` as the only candidate). Mirrors
    select_best_equibind; AutoDock/EquiBind rows are always kept. Returns
    (filtered_df, best_variant_key), or the inputs unchanged with ``None`` when no
    DiffDock variant is present, the summary is missing/unreadable, or no present
    variant has a score."""
    methods = df["docking_method"].astype(str)
    dd_mask = methods.str.startswith("diffdock")
    if not dd_mask.any():
        return df, None
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-diffdock-only] oracle summary not found at {oracle_csv} — "
              "keeping all DiffDock variants.")
        return df, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-diffdock-only")
    if scores is None:
        return df, None
    dd_present = sorted(methods[dd_mask].unique())
    cand = {m: float(scores[m]) for m in dd_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-diffdock-only] none of the present DiffDock variants have a "
              f"score in {oracle_csv} — keeping all DiffDock variants.")
        return df, None

    best = max(cand, key=cand.get)
    keep = (~dd_mask) | (methods == best)
    return df[keep].reset_index(drop=True), best


def select_diffdock_variant(df: pd.DataFrame, variant: str) -> tuple[pd.DataFrame, str | None]:
    """Keep non-DiffDock rows + only the DiffDock poses from optimizer *variant*.

    Explicit, crystal-free alternative to --best-diffdock-only: it names the variant
    directly (raw/original | smina | gnina) instead of ranking by PB-Valid AND RMSD ≤ 2 Å
    (which needs a crystal the Orai receptors lack). AutoDock/EquiBind rows are
    always kept. Returns (filtered_df, kept_label), or the inputs unchanged with
    ``None`` when no DiffDock pose carries that optimizer."""
    want = "original" if str(variant).lower() in ("raw", "original") else str(variant).lower()
    methods = df["docking_method"].astype(str)
    dd_mask = methods.str.startswith("diffdock")
    if not dd_mask.any():
        return df, None
    keep_dd = dd_mask & (df["dd_optimizer"].astype(str) == want)
    if not keep_dd.any():
        present = sorted(df.loc[dd_mask, "dd_optimizer"].dropna().astype(str).unique())
        print(f"  [diffdock-variant] no DiffDock poses with optimizer '{want}' "
              f"(present: {present or 'none'}) — keeping all DiffDock variants.")
        return df, None
    kept_label = methods[keep_dd].iloc[0]
    return df[(~dd_mask) | keep_dd].reset_index(drop=True), kept_label


def select_equibind_variant(df: pd.DataFrame, spec: str) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind rows + only the EquiBind poses whose variant matches *spec*.

    Explicit, crystal-free counterpart to --best-equibind-only (no oracle needed).
    *spec* tokens (split on '_' or '/') are matched against the parsed pocket /
    refine / clamp axes, e.g. 'unguided_gnina' or 'unguided/gnina/clampOFF'; a token
    subset constrains only the axes it names. AutoDock/DiffDock rows are always
    kept. Returns (filtered_df, kept_label) — kept_label is the retained variant
    when *spec* resolves to exactly one, else ``None`` (rows still filtered)."""
    tokens = [t for t in str(spec).strip().lower().replace("/", "_").split("_") if t]
    methods = df["docking_method"].astype(str)
    eq_mask = methods.str.startswith("equibind")
    if not tokens or not eq_mask.any():
        return df, None
    if not {"eq_pocket", "eq_refine", "eq_clamp"} <= set(df.columns):
        print("  [equibind-variant] EquiBind axes not parsed (needs --split-equibind) "
              "— keeping all EquiBind variants.")
        return df, None
    axis_tok = (df[["eq_pocket", "eq_refine", "eq_clamp"]]
                .astype(str).apply(lambda r: {v.lower() for v in r}, axis=1))
    match = eq_mask & axis_tok.apply(lambda s: all(t in s for t in tokens))
    if not match.any():
        present = sorted(methods[eq_mask].unique())
        print(f"  [equibind-variant] no EquiBind poses match '{spec}' "
              f"(present: {present}) — keeping all EquiBind variants.")
        return df, None
    df = df[(~eq_mask) | match].reset_index(drop=True)
    kept = sorted(df.loc[df["docking_method"].astype(str).str.startswith("equibind"),
                         "docking_method"].astype(str).unique())
    return df, (kept[0] if len(kept) == 1 else None)


def _analyzed_count_col(summary: pd.DataFrame) -> str:
    """Name of the column holding the pose count validity is scored against.

    ``within_rmsd_poses`` for an --max-rmsd run (the kept subset), else the plain
    ``total_poses``. Lets the plots stay agnostic to whether the filter was on.
    """
    return "within_rmsd_poses" if "within_rmsd_poses" in summary.columns else "total_poses"


def per_tool_summary(df: pd.DataFrame, order: list[str],
                     generated_counts: pd.Series | None = None,
                     generated_valid_counts: pd.Series | None = None,
                     rmsd2_counts: pd.Series | None = None,
                     rmsd2_valid_counts: pd.Series | None = None) -> pd.DataFrame:
    """Per-method pose/validity summary.

    ``generated_counts`` / ``generated_valid_counts`` (poses produced, and PB-valid
    among them, per method *before* any --max-rmsd filter, indexed by
    ``docking_method``) are threaded in from main() only for a filtered run. When
    present the summary reports the full generated population and the RMSD-limited
    subset side by side:

      generated_poses  generated_valid_poses  generated_valid_fraction  ← ALL poses
      within_rmsd_poses  valid_poses  valid_fraction                    ← ≤ cutoff subset

    so validity is available both over every pose produced and over just the
    near-crystal ones actually scored. Without the filter a single ``total_poses`` /
    ``valid_poses`` / ``valid_fraction`` triple is emitted, exactly as before (the two
    populations would be identical, so the split would be redundant).

    ``rmsd2_counts`` / ``rmsd2_valid_counts`` (poses within NEAR_NATIVE_RMSD_A Å of the
    crystal, and PB-valid among them, per ``docking_method``) add three more columns —
    ``rmsd2_poses`` / ``rmsd2_valid_poses`` / ``rmsd2_valid_fraction`` — reporting the
    same total/valid/fraction triple over just the near-native (≤ 2 Å) subset. Passed
    only for crystal-bearing sets; omitted (columns absent) for crystal-free ones.
    """
    grp = (
        df.groupby("docking_method")
        .agg(analyzed_poses=("pb_valid", "size"),
             valid_poses=("pb_valid", "sum"))
        .reindex(order)
        .dropna(how="all")
    )
    if generated_counts is not None:
        grp = grp.rename(columns={"analyzed_poses": "within_rmsd_poses"})
        # Every method in the (post-filter) summary is a subset of the pre-filter
        # population, so the reindex never introduces a NaN → safe int cast.
        gen = generated_counts.reindex(grp.index).astype("int64")
        grp.insert(0, "generated_poses", gen)
        if generated_valid_counts is not None:
            gval = generated_valid_counts.reindex(grp.index).astype("int64")
            grp.insert(1, "generated_valid_poses", gval)
            grp.insert(2, "generated_valid_fraction", gval / gen)
        analyzed_col = "within_rmsd_poses"
    else:
        grp = grp.rename(columns={"analyzed_poses": "total_poses"})
        analyzed_col = "total_poses"
    grp["valid_fraction"] = grp["valid_poses"] / grp[analyzed_col]

    if rmsd2_counts is not None:
        # Near-native subset: poses within NEAR_NATIVE_RMSD_A Å of the crystal, and
        # PB-valid among them. Methods with none in-cutoff reindex to 0 → fraction NaN.
        r2 = rmsd2_counts.reindex(grp.index).fillna(0).astype("int64")
        grp["rmsd2_poses"] = r2
        if rmsd2_valid_counts is not None:
            r2v = rmsd2_valid_counts.reindex(grp.index).fillna(0).astype("int64")
            grp["rmsd2_valid_poses"] = r2v
            grp["rmsd2_valid_fraction"] = r2v / r2.replace(0, np.nan)
    return grp


def per_pair_tool_matrix(df: pd.DataFrame, value: str, order: list[str]) -> pd.DataFrame:
    """Pivot table: rows=pair, cols=method. value in {'valid', 'total', 'fraction'}."""
    if value == "valid":
        mat = df.pivot_table(index="pair", columns="docking_method",
                             values="pb_valid", aggfunc="sum", fill_value=0)
    elif value == "total":
        mat = df.pivot_table(index="pair", columns="docking_method",
                             values="pb_valid", aggfunc="size", fill_value=0)
    elif value == "fraction":
        valid = df.pivot_table(index="pair", columns="docking_method",
                               values="pb_valid", aggfunc="sum", fill_value=0)
        total = df.pivot_table(index="pair", columns="docking_method",
                               values="pb_valid", aggfunc="size", fill_value=0)
        mat = (valid / total.replace(0, np.nan)).fillna(0)
    else:
        raise ValueError(value)

    cols = [c for c in order if c in mat.columns]
    return mat[cols]


def equibind_axis_breakdown(df: pd.DataFrame) -> pd.DataFrame:
    """Valid-fraction broken down by EquiBind (pocket × refine × clamp) axes."""
    if not {"eq_pocket", "eq_refine", "eq_clamp"}.issubset(df.columns):
        return pd.DataFrame()        # --no-split-equibind: axes weren't parsed
    eq = df[df["docking_method"].str.startswith("equibind")]
    if eq.empty:
        return pd.DataFrame()
    axis = (
        eq.groupby(["eq_pocket", "eq_refine", "eq_clamp"], dropna=False)
        .agg(total_poses=("pb_valid", "size"), valid_poses=("pb_valid", "sum"))
    )
    axis["valid_fraction"] = axis["valid_poses"] / axis["total_poses"]
    return axis


# ─────────────────────── statistics (paired, per-complex) ───────────────────
# HARD RULE: aggregate the tool's many correlated poses to ONE value per complex
# FIRST, then run the paired test across tools/variants. The unit of analysis is the
# receptor-ligand ``pair`` (Benchmark ≈ 303 complexes; Orai×JKU ≈ 12 (frame,ligand)
# pairs). All tests degrade gracefully (try/except at the call site) — a stats failure
# just falls back to the current, test-free figure.

# Below this many complete paired units the omnibus/pairwise tests are drawn but flagged
# 'exploratory (small n)' (Orai×JKU's real n ≈ 12); below _MIN_UNITS_FOR_TEST we skip the
# test entirely (nothing to say) but still draw the figure.
_EXPLORATORY_MAX_UNITS = 20
_MIN_UNITS_FOR_TEST = 5


def _usable_methods(df: pd.DataFrame, methods: list[str], min_frac: float = 0.5) -> list[str]:
    """Drop methods with structural missingness before a listwise-complete paired test.

    A legacy/optional variant present for only a handful of complexes (e.g. the
    ``equibind_guided`` back-compat bucket, 1/303) would otherwise collapse the
    listwise-complete n to a handful. Keep only methods covering ≥ *min_frac* of the
    best-covered method's complexes; order is preserved."""
    tot = df.pivot_table(index="pair", columns="docking_method",
                         values="pb_valid", aggfunc="size")
    cov = tot.notna().sum()
    if cov.empty:
        return list(methods)
    mx = float(cov.max())
    return [m for m in methods if m in cov.index and float(cov[m]) >= min_frac * mx]


def _pair_valid_count_matrix(df: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    """Per-complex PB-valid pose COUNT, one column per method (NaN where a tool
    produced no pose for that pair — informative missingness, kept out of the paired
    test rather than counted as a 0)."""
    total = df.pivot_table(index="pair", columns="docking_method",
                           values="pb_valid", aggfunc="size")
    valid = df.pivot_table(index="pair", columns="docking_method",
                           values="pb_valid", aggfunc="sum")
    valid = valid.where(total.notna())
    cols = [m for m in methods if m in valid.columns]
    return valid[cols]


def _pair_any_valid_matrix(df: pd.DataFrame, methods: list[str]) -> pd.DataFrame:
    """Per-complex boolean '≥1 PB-valid pose' (0/1), one column per method; NaN where
    the tool produced no pose for that pair."""
    cnt = _pair_valid_count_matrix(df, methods)
    return (cnt >= 1).astype(float).where(cnt.notna())


def _pair_check_anypass_matrix(df: pd.DataFrame, methods: list[str],
                               bool_checks: pd.DataFrame, check: str) -> pd.DataFrame:
    """Per-complex boolean 'the tool has ≥1 pose passing *check*' (0/1) per method;
    NaN where the tool produced no pose for that pair. Aggregates poses to the complex
    first (max over the tool's poses)."""
    work = pd.DataFrame({"pair": df["pair"].to_numpy(),
                         "docking_method": df["docking_method"].to_numpy(),
                         "_pass": bool_checks[check].to_numpy().astype(float)})
    anypass = work.groupby(["pair", "docking_method"])["_pass"].max().unstack("docking_method")
    cols = [m for m in methods if m in anypass.columns]
    return anypass[cols]


def _json_safe(o):
    """Recursively coerce numpy / NaN into JSON-serialisable Python."""
    if isinstance(o, dict):
        return {str(k): _json_safe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_json_safe(v) for v in o]
    if isinstance(o, np.ndarray):
        return _json_safe(o.tolist())
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        o = float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, float):
        return None if math.isnan(o) else o
    return o


def _exploratory(n: int) -> bool:
    return n < _EXPLORATORY_MAX_UNITS


def _tag(n: int) -> str:
    return "  — exploratory (small n)" if _exploratory(n) else ""


def _fmt_cont_omnibus(res: dict) -> str | None:
    """One-line Friedman + Kendall's W caption for a paired_continuous result."""
    om = res.get("omnibus", {})
    p, W, chi2, n, dfree = (om.get("p"), om.get("kendall_w"), om.get("chi2"),
                            om.get("n"), om.get("df"))
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return None
    return (f"Friedman χ²({dfree})={chi2:.1f}, p={su.fmt_p(p)} {su.p_stars(p)}; "
            f"Kendall W={W:.2f}; n={n}{_tag(int(n))}")


def _fmt_prop_omnibus(res: dict, label: str = "≥1 valid pose") -> str | None:
    """One-line Cochran's Q caption for a paired_proportions result."""
    om = res.get("omnibus", {})
    p, Q, dfq, n = om.get("p"), om.get("Q"), om.get("df"), res.get("n")
    if p is None or (isinstance(p, float) and math.isnan(p)):
        return None
    return (f"{label}: Cochran Q({dfq})={Q:.1f}, p={su.fmt_p(p)} "
            f"{su.p_stars(p)}; n={n}{_tag(int(n))}")


def _pairwise_caption(pairwise: list[dict], label_fn=_pretty_method,
                      max_show: int = 3) -> str | None:
    """Compact caption of the most-significant pairwise contrasts (Holm-adjusted)."""
    if not pairwise:
        return None
    items = sorted(pairwise, key=lambda d: (d.get("p_holm", 1.0) if
                   d.get("p_holm") == d.get("p_holm") else 1.0))
    parts = []
    for d in items[:max_show]:
        star = d.get("star", "")
        if star in ("", "n/a"):
            continue
        parts.append(f"{label_fn(d['a'])} vs {label_fn(d['b'])} {star}")
    return "pairwise (Holm): " + "; ".join(parts) if parts else None


def _draw_star_brackets(ax, method_to_x: dict, pairwise: list[dict],
                        y0: float, step: float, max_brackets: int = 6) -> float:
    """Draw significance brackets between bars for the significant pairwise contrasts.

    Returns the top y reached (so the caller can extend the y-limit). Only drawn for a
    small number of bars — the caller gates on bar count to avoid overplotting."""
    sig = [d for d in pairwise if d.get("star") in ("*", "**", "***")
           and d["a"] in method_to_x and d["b"] in method_to_x]
    sig.sort(key=lambda d: abs(method_to_x[d["a"]] - method_to_x[d["b"]]))
    sig = sig[:max_brackets]
    top = y0
    for lvl, d in enumerate(sig):
        xa, xb = method_to_x[d["a"]], method_to_x[d["b"]]
        y = y0 + lvl * step
        ax.plot([xa, xa, xb, xb], [y, y + step * 0.25, y + step * 0.25, y],
                color="black", lw=0.8, clip_on=False)
        ax.text((xa + xb) / 2, y + step * 0.25, d["star"], ha="center",
                va="bottom", fontsize=8)
        top = max(top, y + step * 0.5)
    return top


def _variant_paired_stats(df_full: pd.DataFrame, prefix: str,
                          n_units: int) -> dict | None:
    """Paired per-complex stats across the variants of one tool (*prefix*).

    Returns {any_valid: paired_proportions, valid_count: paired_continuous,
    variants:[...], n} or None when the tool has <2 variants / too few units."""
    variants = sorted(m for m in df_full["docking_method"].astype(str).unique()
                      if m.startswith(prefix))
    variants = _usable_methods(df_full, variants)
    if len(variants) < 2:
        return None
    cnt = _pair_valid_count_matrix(df_full, variants).dropna()
    if len(cnt) < _MIN_UNITS_FOR_TEST:
        return {"variants": variants, "n": int(len(cnt)),
                "note": "n too small — exploratory"}
    anyv = (cnt >= 1).astype(float)
    out: dict = {"variants": variants, "n": int(len(cnt))}
    try:
        out["any_valid"] = su.paired_proportions({v: anyv[v].to_numpy() for v in variants})
    except Exception as e:  # pragma: no cover - degrade gracefully
        print(f"  [stats] {prefix} variant any-valid failed: {e}")
    # Friedman needs >=3 conditions; with 2 variants the McNemar in any_valid is the test.
    if len(variants) >= 3:
        try:
            out["valid_count"] = su.paired_continuous(
                {v: cnt[v].to_numpy() for v in variants})
        except Exception as e:  # pragma: no cover
            print(f"  [stats] {prefix} variant valid-count failed: {e}")
    return out


def compute_validity_stats(df_plot: pd.DataFrame, df_full: pd.DataFrame,
                           order_plot: list[str], order_full: list[str],
                           dataset_label: str) -> dict:
    """All paired, per-complex validity statistics for the report.

    Aggregates the correlated poses to one value per receptor-ligand ``pair`` FIRST,
    then runs the paired tests across tools (df_plot / order_plot) and across each
    tool's variants (df_full / order_full). Returns a JSON-ready dict; individual
    tests are wrapped so one failure never sinks the rest."""
    n_units = int(df_plot["pair"].nunique())
    # Drop structurally-sparse methods (e.g. the 1-complex legacy 'equibind_guided'
    # bucket) so listwise completeness across the remaining tools doesn't collapse n.
    usable_plot = _usable_methods(df_plot, order_plot)
    stats: dict = {"dataset": dataset_label, "unit_of_analysis": "receptor-ligand pair",
                   "n_units": n_units,
                   "exploratory": _exploratory(n_units),
                   "min_units_for_test": _MIN_UNITS_FOR_TEST,
                   "methods_tested": usable_plot,
                   "methods_dropped_sparse": [m for m in order_plot if m not in usable_plot]}

    # ── Across tools: per-complex valid-pose count + ≥1-valid (figs 01/03/04) ──
    cnt = _pair_valid_count_matrix(df_plot, usable_plot).dropna()
    stats["n_complete_units"] = int(len(cnt))
    if len(cnt) >= _MIN_UNITS_FOR_TEST and cnt.shape[1] >= 2:
        if cnt.shape[1] >= 3:            # Friedman needs >=3 tools
            try:
                stats["valid_count"] = su.paired_continuous(
                    {m: cnt[m].to_numpy() for m in cnt.columns})
            except Exception as e:  # pragma: no cover
                print(f"  [stats] valid-count paired_continuous failed: {e}")
        try:
            anyv = (cnt >= 1).astype(float)
            stats["any_valid"] = su.paired_proportions(
                {m: anyv[m].to_numpy() for m in anyv.columns})
        except Exception as e:  # pragma: no cover
            print(f"  [stats] any-valid paired_proportions failed: {e}")
    else:
        stats["note_tools"] = "n too small — exploratory"

    # ── Per-check family: '≥1 pose passes check c' across tools, BH over checks (fig 05) ──
    checks = _present_checks(df_plot)
    if checks and cnt.shape[1] >= 2:
        bchecks = _bool_checks(df_plot, checks)
        per_check: dict = {}
        pvals, keys = [], []
        for chk in checks:
            try:
                m = _pair_check_anypass_matrix(df_plot, usable_plot, bchecks, chk).dropna()
                if len(m) < _MIN_UNITS_FOR_TEST or m.shape[1] < 2:
                    continue
                res = su.paired_proportions({c: m[c].to_numpy() for c in m.columns})
                per_check[chk] = res
                pvals.append(res["omnibus"]["p"])
                keys.append(chk)
            except Exception as e:  # pragma: no cover
                print(f"  [stats] per-check '{chk}' failed: {e}")
        if pvals:
            q = su.bh_fdr(pvals)
            sig = []
            for k, p, qq in zip(keys, pvals, q):
                per_check[k]["omnibus"]["p_bh"] = float(qq)
                per_check[k]["omnibus"]["star_bh"] = su.p_stars(qq)
                if qq == qq and qq < 0.05:
                    sig.append(k)
            stats["per_check_significant_bh"] = sig
        stats["per_check"] = per_check

    # ── Variants are paired (same poses minimized): raw vs smina vs gnina (06/07/07b) ──
    for prefix, key in (("equibind", "equibind_variants"),
                        ("diffdock", "diffdock_variants")):
        try:
            res = _variant_paired_stats(df_full, prefix, n_units)
            if res is not None:
                stats[key] = res
        except Exception as e:  # pragma: no cover
            print(f"  [stats] {prefix} variant stats failed: {e}")

    return stats


# ───────────────────────────── plots ─────────────────────────────

def plot_per_tool(summary: pd.DataFrame, out: Path, stats: dict | None = None) -> None:
    count_col = _analyzed_count_col(summary)
    has_gen = "generated_poses" in summary.columns
    analyzed_label = "Poses within cutoff" if count_col == "within_rmsd_poses" else "Total poses"
    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(summary)), 5))
    x = np.arange(len(summary))
    analyzed = summary[count_col].astype(float)
    valid = summary["valid_poses"].astype(float)

    if has_gen:
        # Three bars per method: all poses produced (generated) → the subset within
        # the RMSD cutoff (analyzed/scored) → the PB-valid subset of that.
        w = 0.28
        gen = summary["generated_poses"].astype(float)
        ax.bar(x - w, gen, w, label="Generated poses",
               color="#7f7f7f", edgecolor="black")
        ax.bar(x, analyzed, w, label=analyzed_label,
               color="#bdbdbd", edgecolor="black")
        ax.bar(x + w, valid, w, label="PB-Valid poses",
               color="#2ca02c", edgecolor="black")
        for i, (g, tot, val) in enumerate(zip(gen, analyzed, valid)):
            ax.text(i - w, g, f"{int(g)}", ha="center", va="bottom", fontsize=8)
            ax.text(i, tot, f"{int(tot)}", ha="center", va="bottom", fontsize=8)
            pct = (val / tot * 100) if tot else 0
            ax.text(i + w, val, f"{int(val)}\n({pct:.1f}%)",
                    ha="center", va="bottom", fontsize=8)
    else:
        w = 0.38
        ax.bar(x - w / 2, analyzed, w, label=analyzed_label,
               color="#bdbdbd", edgecolor="black")
        ax.bar(x + w / 2, valid, w, label="PB-Valid poses",
               color="#2ca02c", edgecolor="black")
        for i, (tot, val) in enumerate(zip(analyzed, valid)):
            pct = (val / tot * 100) if tot else 0
            ax.text(i + w / 2, val, f"{int(val)}\n({pct:.1f}%)",
                    ha="center", va="bottom", fontsize=9)
            ax.text(i - w / 2, tot, f"{int(tot)}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([_pretty_method(t) for t in summary.index],
                       rotation=20, ha="right")
    ax.set_ylabel("Number of poses")
    title = ("PoseBusters Benchmark — Pose validity per docking method\n"
             "(valid = passes all canonical PB checks)" + _RMSD_FILTER_NOTE)
    # Paired per-complex tests (aggregated to one value per receptor-ligand pair first):
    # valid-pose COUNT (Friedman + Kendall W) and ≥1-valid-pose (Cochran Q).
    if stats:
        try:
            caps = [c for c in (_fmt_cont_omnibus(stats.get("valid_count", {})),
                                _fmt_prop_omnibus(stats.get("any_valid", {}))) if c]
            if caps:
                # Each omnibus on its own line so the caption never overflows the width;
                # per-complex pairwise stars are drawn on fig 04 (matched y-axis).
                title += "\nper-complex paired: " + caps[0]
                title += "".join("\n" + c for c in caps[1:])
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_per_tool annotation failed: {e}")
    ax.set_title(title, fontsize=10)
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_heatmap(valid_mat: pd.DataFrame, out: Path, top_n: int = 60,
                 worst: bool = False, pose_scale: int = 31) -> None:
    """Heatmap of valid-pose counts per receptor-ligand pair × method.

    Pairs are ranked by their total valid poses across methods and the extreme *top_n*
    are shown: the best-docking pairs by default, or the *worst*-docking pairs (fewest
    valid poses — the hardest complexes) when ``worst=True``. The colour range is fixed
    to 0..*pose_scale* (default 31, the max poses produced per complex) so the top (02)
    and bottom (02b) heatmaps share one scale and are directly comparable."""
    mat = valid_mat.copy()
    mat["__sum"] = mat.sum(axis=1)
    mat = mat.sort_values("__sum", ascending=worst).head(top_n).drop(columns="__sum")
    mat.columns = [_pretty_method(c) for c in mat.columns]

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(mat.columns)), max(8, 0.22 * len(mat))))
    sns.heatmap(mat, annot=True, fmt=".0f", cmap="YlGnBu", vmin=0, vmax=pose_scale,
                cbar_kws={"label": f"Valid poses (shared 0–{pose_scale} scale)"},
                linewidths=0.3, linecolor="white", ax=ax)
    ax.set_xlabel("Docking method")
    ax.set_ylabel("Receptor / Ligand")
    kind = "worst" if worst else "top"
    ax.set_title(f"Valid poses per receptor-ligand pair ({kind} {len(mat)} — "
                 f"{'hardest' if worst else 'easiest'} complexes)" + _RMSD_FILTER_NOTE)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_grouped_bars(valid_mat: pd.DataFrame, out: Path, colors: dict,
                      top_n: int = 30, stats: dict | None = None) -> None:
    # Cleveland dot plot: one row per receptor-ligand pair, one coloured dot per
    # tool, joined by a faint connector. Replaces a dense top-N grouped-bar wall
    # and complements the per-pair validity heatmap (fig 02).
    mat = valid_mat.copy()
    mat["__sum"] = mat.sum(axis=1)
    mat = mat.sort_values("__sum", ascending=False).head(top_n).drop(columns="__sum")

    pairs = mat.index.tolist()
    tools = list(mat.columns)
    y = np.arange(len(pairs))[::-1]   # highest-total pair on top

    fig, ax = plt.subplots(figsize=(max(8, 0.9 * len(tools) + 5.0),
                                    max(5, 0.32 * len(pairs) + 1.5)))
    for yi, p in zip(y, pairs):
        vals = [mat.loc[p, t] for t in tools]
        ax.plot([min(vals), max(vals)], [yi, yi], color="#dddddd", lw=2,
                zorder=1, solid_capstyle="round")
    for t in tools:
        ax.scatter(mat[t].to_numpy(), y, s=45, color=colors.get(t, None),
                   edgecolor="black", linewidth=0.4, label=_pretty_method(t), zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(pairs, fontsize=8)
    ax.set_xlabel("Valid poses")
    title = (f"Valid poses per docking method — top {len(pairs)} receptor-ligand pairs"
             + _RMSD_FILTER_NOTE)
    # Whole-set paired test (all complexes, not just the top-N shown): does the tool a
    # complex gets ≥1 valid pose from differ? Cochran's Q + pairwise McNemar/Holm.
    if stats:
        try:
            cap = _fmt_prop_omnibus(stats.get("any_valid", {}))
            pw = _pairwise_caption((stats.get("any_valid") or {}).get("pairwise", []))
            if cap:
                title += "\n(all complexes) " + cap
            if pw:
                title += "\n" + pw
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_grouped_bars annotation failed: {e}")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_validity_distribution(df: pd.DataFrame, out: Path, order: list[str],
                               stats: dict | None = None) -> None:
    """Distribution of per-pair valid-pose counts by method (boxplot)."""
    valid = df.pivot_table(index="pair", columns="docking_method",
                           values="pb_valid", aggfunc="sum", fill_value=0)
    cols = [c for c in order if c in valid.columns]
    valid = valid[cols]

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(cols)), 5))
    # Set tick labels via the Axes (not boxplot's renamed 'labels'/'tick_labels'
    # kwarg) so this works across Matplotlib versions without deprecation noise.
    bp = ax.boxplot([valid[c].values for c in cols], showmeans=True)
    # Whiskers (and their caps) in a distinct colour so they stand out from the boxes.
    for el in bp["whiskers"] + bp["caps"]:
        el.set_color("#d62728")
        el.set_linewidth(1.5)
    # Total poses assessed per variant, shown under each label.
    counts = df.groupby("docking_method").size()
    ax.set_xticks(range(1, len(cols) + 1))
    ax.set_xticklabels([f"{_pretty_method(c)}\n(n = {int(counts.get(c, 0)):,})" for c in cols],
                       rotation=20, ha="right")
    ax.set_ylabel("Valid poses per receptor-ligand pair")
    title = ("Distribution of valid poses across receptor-ligand pairs\n"
             "(n = poses assessed per variant; whiskers/caps in red)"
             + _RMSD_FILTER_NOTE)
    # Paired per-complex valid-pose COUNT: Friedman + Kendall's W, Wilcoxon post-hoc (Holm).
    if stats:
        try:
            vc = stats.get("valid_count")
            cap = _fmt_cont_omnibus(vc or {})
            if cap:
                title += "\nper-complex paired: " + cap
            if vc and len(cols) <= 5:
                m2x = {m: i + 1 for i, m in enumerate(cols)}
                top = float(np.nanmax(valid.to_numpy())) * 1.05 + 0.5
                ytop = _draw_star_brackets(ax, m2x, vc.get("pairwise", []),
                                           top, max(0.5, top * 0.06))
                ax.set_ylim(top=max(ax.get_ylim()[1], ytop * 1.05))
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_validity_distribution annotation failed: {e}")
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_check_passrate(df: pd.DataFrame, out: Path, order: list[str],
                        stats: dict | None = None) -> None:
    """Per-check pass rate per method (which PB criterion is failing most)."""
    checks = _present_checks(df)
    rows = []
    for tool, sub in df.groupby("docking_method"):
        b = _bool_checks(sub, checks)
        for chk in checks:
            rows.append({"tool": tool, "check": chk,
                         "pass_rate": b[chk].mean()})
    pr = pd.DataFrame(rows).pivot(index="check", columns="tool", values="pass_rate")
    cols = [c for c in order if c in pr.columns]
    pr = pr[cols]
    # Column labels carry each variant's pose count (n = …) so the pass rates can be
    # read against how many poses they were computed over.
    counts = df.groupby("docking_method").size()
    pr.columns = [f"{_pretty_method(c)}\n(n = {int(counts.get(c, 0)):,})" for c in cols]

    fig, ax = plt.subplots(figsize=(max(8, 1.2 * len(pr.columns)), 6))
    sns.heatmap(pr * 100, annot=True, fmt=".1f", cmap="RdYlGn", vmin=0, vmax=100,
                cbar_kws={"label": "Pass rate (%)"}, ax=ax)
    title = ("Per-check pass rate (%) per docking method (n = poses per variant)"
             + _RMSD_FILTER_NOTE)
    # Per-complex paired test on EACH check ('≥1 pose passes it'): Cochran's Q, BH-corrected
    # across the checks. Star (★) the checks where tools genuinely differ — this is what
    # substantiates the 'EquiBind fails bond angles/lengths' claim.
    if stats:
        try:
            sig = set(stats.get("per_check_significant_bh", []))
            per_check = stats.get("per_check", {})
            if per_check:
                ylabels = []
                for chk in pr.index:                       # raw check names
                    star = "★ " if chk in sig else ""
                    ylabels.append(star + str(chk))
                ax.set_yticklabels(ylabels, rotation=0)
                title += ("\n★ tools differ (per-complex ≥1-pose-passes; "
                          "Cochran Q, BH q<0.05)")
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_check_passrate annotation failed: {e}")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("PoseBusters critical check")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def _variant_short_label(m: str, prefix: str) -> str:
    """Short, override-free label for one variant of the *prefix* tool.

    Deliberately ignores ``_LABEL_OVERRIDES`` (the 'DiffDock*'/'EquiBind*' stars) so
    per-variant comparison figures always show the true optimizer/pocket identity of
    every bar — including the variant chosen as 'best' elsewhere."""
    if prefix == "diffdock":
        return {"diffdock": "original", "diffdock_smina": "smina-opt",
                "diffdock_gnina": "gnina-opt"}.get(m, m)
    if prefix == "equibind":
        pocket, refine, clamp = _eq_tokens(m)
        parts: list[str] = []
        if pocket:
            parts.append(pocket)
        if refine:
            parts.append({"smina": "smina-opt", "gnina": "gnina-opt"}.get(refine, "raw"))
        if clamp:
            parts.append("clamp on" if clamp == "clampON" else "clamp off")
        return ", ".join(parts) if parts else "base"
    return m


def plot_tool_variant_comparison(summary: pd.DataFrame, out: Path, colors: dict, *,
                                 prefix: str, tool_label: str,
                                 subtitle: str = "", stats: dict | None = None) -> bool:
    """Bar chart comparing PB-validity (%) across the variants of one docking tool.

    Filters *summary* to methods whose key starts with *prefix* (``diffdock`` /
    ``equibind``) and draws one bar per variant, annotated with the pass rate and
    valid/total counts. Returns False (nothing written) when the tool has < 2 variants —
    there is nothing to compare. Labels bypass the 'best'-variant star so every
    optimizer/pocket flavour stays identifiable."""
    sub = summary[summary.index.str.startswith(prefix)]
    if len(sub) < 2:
        return False
    sub = sub.loc[sorted(sub.index, key=lambda m: (_method_sort_key(m), m))]

    count_col = _analyzed_count_col(summary)
    labels = [_variant_short_label(m, prefix) for m in sub.index]
    x = np.arange(len(sub))
    fig, ax = plt.subplots(figsize=(max(7, 1.6 * len(sub)), 5))
    bars = ax.bar(x, sub["valid_fraction"] * 100, 0.6,
                  color=[colors.get(m, "#777777") for m in sub.index],
                  edgecolor="black")
    for rect, (tot, val) in zip(bars, zip(sub[count_col], sub["valid_poses"])):
        pct = (val / tot * 100) if tot else 0
        ax.text(rect.get_x() + rect.get_width() / 2, pct,
                f"{pct:.1f}%\n{int(val)}/{int(tot)}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("PB-Valid poses (%)")
    ax.set_ylim(0, max(5, sub["valid_fraction"].max() * 100 * 1.25))
    title = f"{tool_label} — PB-validity per variant"
    if subtitle:
        title += "\n" + subtitle
    title += _RMSD_FILTER_NOTE
    # Variants are PAIRED (the same complexes' poses minimized): per-complex ≥1-valid
    # (Cochran's Q) + raw-vs-smina-vs-gnina pairwise McNemar/Holm brackets. Validates
    # gnina > smina.
    vstats = (stats or {}).get(f"{prefix}_variants")
    if vstats:
        try:
            av = vstats.get("any_valid")
            vc = vstats.get("valid_count")
            caps = [c for c in (_fmt_prop_omnibus(av or {}),
                                _fmt_cont_omnibus(vc or {}) and
                                ("valid-count " + _fmt_cont_omnibus(vc))) if c]
            if caps:
                title += "\nper-complex paired: " + caps[0]
                title += "".join("\n" + c for c in caps[1:])
            # Prefer the valid-COUNT pairwise for the brackets — the ≥1-valid level
            # saturates (both smina & gnina almost always get one), so gnina>smina only
            # shows in the count. Fall back to ≥1-valid pairwise when count is absent (k=2).
            pw_src = (vc or av or {}).get("pairwise", [])
            if pw_src and len(sub) <= 6:
                m2x = {m: i for i, m in enumerate(sub.index)}
                top = float(sub["valid_fraction"].max() * 100) * 1.10 + 2
                ytop = _draw_star_brackets(ax, m2x, pw_src,
                                           top, max(2.0, top * 0.06))
                ax.set_ylim(top=max(ax.get_ylim()[1], ytop * 1.05))
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_tool_variant_comparison annotation failed: {e}")
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def plot_equibind_variants(summary: pd.DataFrame, out: Path, colors: dict,
                           stats: dict | None = None) -> bool:
    """Valid fraction (%) per EquiBind variant. Returns False if < 2 variants.

    Thin wrapper over plot_tool_variant_comparison kept for its original call site."""
    return plot_tool_variant_comparison(
        summary, out, colors, prefix="equibind", tool_label="EquiBind",
        subtitle="(fpocket / p2rank / unguided, smina re-search, centroid clamp)",
        stats=stats)


def _variant_matrix_cell(m: str) -> tuple[str, str] | None:
    """Map a variant key to (row, column) for the base-method × optimizer validity
    matrix, or None for methods with no place in it. EquiBind and DiffDock share a
    raw/smina/gnina optimizer axis, so both slot into the same three columns. AutoDock has
    no optimizer axis, so — for completeness — it occupies the 'raw / original' column
    only (its smina/gnina cells stay blank)."""
    col_of = {"": "raw / original", "raw": "raw / original", "original": "raw / original",
              "smina": "smina", "gnina": "gnina"}
    if m == "autodock":
        return "AutoDock Vina", "raw / original"
    if m.startswith("diffdock"):
        opt = m[len("diffdock"):].lstrip("_")            # "", "smina", "gnina"
        return "DiffDock", col_of.get(opt, opt)
    if m.startswith("equibind"):
        pocket, refine, clamp = _eq_tokens(m)
        row = f"EquiBind · {pocket or '?'}"
        if clamp:
            row += f" · clamp {'on' if clamp == 'clampON' else 'off'}"
        return row, col_of.get(refine or "", "raw / original")
    return None


def plot_variant_validity_matrix(summary: pd.DataFrame, out: Path,
                                 stats: dict | None = None) -> bool:
    """Heatmap of PB-validity (%) per variant — the matrix view of the per-tool variant
    bar charts (figs 06/07). Rows are the base method (AutoDock, DiffDock, and EquiBind by
    pocket × clamp); columns are the shared optimizer axis (raw/original, smina, gnina);
    each cell is that variant's PB-valid fraction. AutoDock has no optimizer variants, so
    only its raw/original cell is filled. Returns False when fewer than two variant cells
    exist (nothing worth a matrix)."""
    cells: dict[str, dict[str, float]] = {}
    for m in summary.index:
        rc = _variant_matrix_cell(str(m))
        if rc is None:
            continue
        r, c = rc
        cells.setdefault(r, {})[c] = float(summary.loc[m, "valid_fraction"]) * 100.0
    if not cells:
        return False
    col_order = ["raw / original", "smina", "gnina"]
    mat = pd.DataFrame(cells).T
    mat = mat.reindex(columns=[c for c in col_order if c in mat.columns])
    lead = [r for r in ("AutoDock Vina", "DiffDock") if r in mat.index]
    row_order = lead + sorted(r for r in mat.index if r not in lead)
    mat = mat.reindex(row_order)
    if int(mat.notna().values.sum()) < 2:
        return False

    fig, ax = plt.subplots(figsize=(max(6.0, 1.7 * mat.shape[1] + 3.0),
                                    max(3.5, 0.62 * mat.shape[0] + 1.6)))
    sns.heatmap(mat, annot=True, fmt=".1f", cmap="YlGn", vmin=0, vmax=100,
                linewidths=0.4, linecolor="white", mask=mat.isna(),
                cbar_kws={"label": "PB-Valid poses (%)"}, ax=ax)
    ax.set_xlabel("Optimizer / re-search")
    ax.set_ylabel("Base method")
    title = ("PB-validity (%) per variant — base method × optimizer" + _RMSD_FILTER_NOTE)
    # Paired per-complex omnibus across each tool's optimizer variants (Cochran's Q on
    # ≥1-valid); the raw-vs-smina-vs-gnina pairwise stars live in the 06/07 bar charts.
    if stats:
        try:
            caps = []
            for prefix, name in (("equibind", "EquiBind"), ("diffdock", "DiffDock")):
                v = stats.get(f"{prefix}_variants", {})
                c = _fmt_prop_omnibus((v or {}).get("any_valid", {}))
                if c:
                    caps.append(f"{name} {c}")
            if caps:
                title += "\nper-complex paired: " + "\n".join(caps)
        except Exception as e:  # pragma: no cover
            print(f"  [stats] plot_variant_validity_matrix annotation failed: {e}")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def _join_rmsd_to_crystal(df: pd.DataFrame, metrics_csv: Path, rmsd_col: str,
                          label: str = "rmsd-sweep") -> pd.DataFrame | None:
    """Return *df* with a numeric ``rmsd_to_crystal`` column joined from *metrics_csv*
    (per_pose_metrics.csv) on (protein, ligand, pose_name), or None when the metrics
    file / column / join keys are unavailable. Unlike apply_rmsd_filter this never
    exits — callers (the RMSD sweep, the ≤2 Å summary columns) degrade gracefully (e.g.
    crystal-free Orai) instead of aborting the run. *label* only tags the skip message."""
    if not metrics_csv or not Path(metrics_csv).exists():
        print(f"  [{label}] per-pose metrics not found at {metrics_csv} — skipping.")
        return None
    metrics = pd.read_csv(metrics_csv, low_memory=False)
    key = ["protein", "ligand", "pose_name"]
    if rmsd_col not in metrics.columns or any(k not in metrics.columns for k in key):
        print(f"  [{label}] {metrics_csv} lacks '{rmsd_col}' or join keys {key} — skipping.")
        return None
    m = metrics[key + [rmsd_col]].copy()
    for k in key:
        m[k] = m[k].astype(str)
    m = m.drop_duplicates(subset=key)
    d = df.copy()
    for k in key:
        d[k] = d[k].astype(str)
    d = d.merge(m.rename(columns={rmsd_col: "rmsd_to_crystal"}), on=key,
                how="left", validate="m:1")
    d["rmsd_to_crystal"] = pd.to_numeric(d["rmsd_to_crystal"], errors="coerce")
    return d


# Near-native RMSD-to-crystal cutoff for the rmsd2_* summary columns (the canonical
# docking-success threshold; matches the oracle_rmsd_le_2.0A_% / pb_rmsd_within_2A metrics).
NEAR_NATIVE_RMSD_A = 2.0


def _rmsd_le2_counts(df: pd.DataFrame, metrics_csv: Path, rmsd_col: str,
                     thresh: float = NEAR_NATIVE_RMSD_A):
    """Per-method (poses, PB-valid poses) with RMSD-to-crystal ≤ *thresh* Å.

    Reuses the ``rmsd_to_crystal`` column already joined by --max-rmsd when present,
    otherwise joins it from *metrics_csv*. Returns (None, None) for a crystal-free set
    (no metrics / no RMSD) or when nothing lands within the cutoff, so the caller simply
    omits the rmsd2_* summary columns."""
    if "rmsd_to_crystal" in df.columns:
        d = df
    else:
        d = _join_rmsd_to_crystal(df, metrics_csv, rmsd_col, label="rmsd≤2")
        if d is None:
            return None, None
    r = pd.to_numeric(d["rmsd_to_crystal"], errors="coerce")
    near = d[r <= thresh]
    if near.empty:
        print(f"  [rmsd≤2] no pose within {thresh:g} Å of the crystal — "
              "omitting rmsd2_* summary columns.")
        return None, None
    return (near.groupby("docking_method").size(),
            near.groupby("docking_method")["pb_valid"].sum())


def plot_validity_vs_rmsd(df: pd.DataFrame, metrics_csv: Path, rmsd_col: str,
                          max_cutoff: float, steps: int, order: list[str],
                          colors: dict, out_count: Path, out_yield: Path) -> list[str]:
    """Two separate sweeps of how PB-validity accumulates as the RMSD-to-crystal cutoff
    is relaxed from 0 to *max_cutoff* Å, written as two files:

    out_count : cumulative COUNT of PB-valid poses with RMSD ≤ τ, per method.
    out_yield : the same as a %% of each method's total poses (a comparable 'yield' curve).

    RMSD-to-crystal is joined from *metrics_csv*; returns [] (no files) for a set with no
    usable RMSD (e.g. crystal-free Orai), else the basenames of the figures written."""
    d = _join_rmsd_to_crystal(df, metrics_csv, rmsd_col)
    if d is None:
        return []
    d = d[d["rmsd_to_crystal"].notna()]
    if d.empty:
        print("  [rmsd-sweep] no pose has an RMSD-to-crystal value — skipping figure.")
        return []

    taus = np.linspace(0.0, float(max_cutoff), max(2, int(steps)))
    present = [m for m in order if (d["docking_method"] == m).any()]
    curves = {}
    for method in present:
        sub = d[d["docking_method"] == method]
        rvalid = np.sort(sub.loc[sub["pb_valid"], "rmsd_to_crystal"].to_numpy())
        cum = np.searchsorted(rvalid, taus, side="right").astype(float)
        curves[method] = (cum, len(sub))

    note = _RMSD_FILTER_NOTE.replace("\n", " ")

    # Figure 1 — cumulative count.
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    for method in present:
        cum, _ = curves[method]
        ax.plot(taus, cum, label=_pretty_method(method), color=colors.get(method), lw=1.9)
    ax.set_xlabel(f"RMSD-to-crystal cutoff ({rmsd_col}, Å)")
    ax.set_ylabel("PB-valid poses within cutoff (cumulative count)")
    ax.set_title("Valid poses accumulated vs RMSD cutoff" + note)
    ax.set_xlim(0, float(max_cutoff))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_count, dpi=160)
    plt.close(fig)

    # Figure 2 — yield (% of the method's poses).
    fig, ax = plt.subplots(figsize=(8.5, 5.4))
    for method in present:
        cum, total = curves[method]
        ax.plot(taus, (cum / total * 100.0) if total else cum,
                label=_pretty_method(method), color=colors.get(method), lw=1.9)
    ax.set_xlabel(f"RMSD-to-crystal cutoff ({rmsd_col}, Å)")
    ax.set_ylabel("PB-valid poses within cutoff (% of the method's poses)")
    ax.set_title("Valid-pose yield vs RMSD cutoff" + note)
    ax.set_xlim(0, float(max_cutoff))
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(out_yield, dpi=160)
    plt.close(fig)

    return [Path(out_count).name, Path(out_yield).name]


def _prettify_check(c: str) -> str:
    """Shorten a PoseBusters check column name for a compact table header."""
    return (str(c).replace("minimum_distance_to_", "min_dist_")
             .replace("volume_overlap_with_", "vol_ovlp_")
             .replace("non-aromatic", "non-arom")
             .replace("_", " "))


def check_failure_table(df: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    """Per-variant PoseBusters failure-rate table (variants × checks).

    Rows are the methods/variants (in *order*); the first two columns are ``n_poses``
    (total poses that variant produced) and ``survival_%`` (poses passing ALL checks —
    the PB-valid rate), followed by one column per PB check giving the %% of that
    variant's poses that FAIL it (100 − pass rate), worst-failing check leftmost. This is
    the tabular answer to 'how many poses, how many survive, and which PB test fails by
    how much, per variant' — the counterpart to the pass-rate heatmap (fig 05)."""
    checks = _present_checks(df)
    methods = [m for m in order if (df["docking_method"] == m).any()]
    if not checks or not methods:
        return pd.DataFrame()
    data, n_poses, survival = {}, {}, {}
    for m in methods:
        sub = df[df["docking_method"] == m]
        b = _bool_checks(sub, checks)
        data[m] = {chk: (1.0 - float(b[chk].mean())) * 100.0 for chk in checks}
        n_poses[m] = len(sub)
        survival[m] = float(sub["pb_valid"].mean()) * 100.0
    # rows = variant (method order); check columns worst-failing first, with the
    # per-variant pose count and overall survival (%) pinned in front.
    tbl = pd.DataFrame(data).T.reindex(index=methods, columns=checks)
    tbl = tbl[tbl.max(axis=0).sort_values(ascending=False).index]
    tbl.insert(0, "survival_%", pd.Series(survival))
    tbl.insert(0, "n_poses", pd.Series(n_poses))
    tbl.index.name = "docking_method"
    return tbl


def _top_variant_by_validity(df: pd.DataFrame, prefix: str) -> str | None:
    """The variant of the *prefix* tool with the highest PB-valid fraction (its 'best'
    variant for a validity/failure view), or None when the tool is absent."""
    sub = df[df["docking_method"].astype(str).str.startswith(prefix)]
    if sub.empty:
        return None
    return str(sub.groupby("docking_method")["pb_valid"].mean().idxmax())


def _top_variant_for_waterfall(df: pd.DataFrame, prefix: str,
                               oracle_csv: Path) -> str | None:
    """Top variant of *prefix* for the failure waterfall (fig 09).

    Prefers the combined-metric best — PB-Valid AND RMSD ≤ 2 Å, ranked from the oracle
    summary via _variant_ranking_scores — so the waterfall highlights the SAME 'best'
    variant the collapsed comparison figures do (select_best_*). Falls back to the highest
    PB-valid fraction when no oracle score is available: a crystal-free set (Orai with no
    borrowed ranking) or a single-variant tool (AutoDock) with nothing to rank."""
    present = sorted({m for m in df["docking_method"].astype(str) if m.startswith(prefix)})
    if not present:
        return None
    if len(present) == 1:
        return present[0]
    if oracle_csv and Path(oracle_csv := _resolve_variant_oracle(oracle_csv)).exists():
        scores, _ = _variant_ranking_scores(oracle_csv, f"waterfall:{prefix}")
        if scores is not None:
            cand = {m: float(scores[m]) for m in present
                    if m in scores.index and pd.notna(scores[m])}
            if cand:
                return max(cand, key=cand.get)
    return _top_variant_by_validity(df, prefix)


def _waterfall_steps(bool_df: pd.DataFrame, checks: list[str]):
    """Pareto attrition steps for one variant's poses.

    Orders *checks* by descending independent fail rate, then walks them accumulating the
    *marginal* fraction of poses each removes among those still passing every earlier
    check. Returns (steps, valid_pct) where steps is a list of
    (check, running_before_%, running_after_%) for checks that actually drop poses, and
    valid_pct is the final PB-valid %. Order-dependent by construction (a pose failing
    several checks is charged to the first) — which is what makes the biggest culprit the
    tallest step."""
    n = len(bool_df)
    if n == 0:
        return [], 100.0
    fail_rate = {c: 1.0 - float(bool_df[c].mean()) for c in checks}
    ordered = [c for c in sorted(checks, key=lambda c: -fail_rate[c]) if fail_rate[c] > 0]
    remaining = pd.Series(True, index=bool_df.index)
    running = 100.0
    steps = []
    for c in ordered:
        new_remaining = remaining & bool_df[c]
        drop = (int(remaining.sum()) - int(new_remaining.sum())) / n * 100.0
        remaining = new_remaining
        if drop <= 0:
            continue
        steps.append((c, running, running - drop))
        running -= drop
    return steps, running


def plot_failure_waterfall(df: pd.DataFrame, methods: list[str], out: Path) -> bool:
    """Per-check pose-attrition waterfall for each variant in *methods* (one panel each).

    Each panel starts at 100% of that variant's poses and drops one red step per failing
    PB check (Pareto order — biggest culprit first), ending at the green PB-valid %. A
    step's height is the marginal % of poses that check removes, so the panel reads
    directly as 'which PB test the poses fail, and by how much'. Returns False when no
    checks or none of *methods* are present."""
    checks = _present_checks(df)
    methods = [m for m in methods if (df["docking_method"] == m).any()]
    if not checks or not methods:
        return False
    fig, axes = plt.subplots(1, len(methods), figsize=(6.3 * len(methods), 5.4),
                             squeeze=False)
    for idx, method in enumerate(methods):
        ax = axes[0][idx]
        sub = df[df["docking_method"] == method]
        b = _bool_checks(sub, checks)
        steps, valid_pct = _waterfall_steps(b, checks)

        labels = ["All poses"] + [_prettify_check(c) for c, _, _ in steps] + ["PB-Valid"]
        ax.bar(0, 100.0, 0.62, color="#bdbdbd", edgecolor="black")
        ax.text(0, 100.0, f"n={len(sub)}", ha="center", va="bottom", fontsize=8)
        prev_top = 100.0
        for i, (c, top, bot) in enumerate(steps, start=1):
            ax.bar(i, top - bot, 0.62, bottom=bot, color="#d62728", edgecolor="black")
            ax.text(i, top, f"-{top - bot:.1f}", ha="center", va="bottom", fontsize=7)
            ax.plot([i - 1 + 0.31, i - 0.31], [prev_top, top], color="#888",
                    ls="--", lw=0.7)
            prev_top = bot
        fx = len(steps) + 1
        ax.bar(fx, valid_pct, 0.62, color="#2ca02c", edgecolor="black")
        ax.plot([fx - 1 + 0.31, fx - 0.31], [prev_top, valid_pct], color="#888",
                ls="--", lw=0.7)
        ax.text(fx, valid_pct, f"{valid_pct:.1f}%", ha="center", va="bottom", fontsize=9)

        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_ylim(0, 112)
        if idx == 0:
            ax.set_ylabel("% of the variant's poses still passing")
        ax.set_title(_pretty_method(method), fontsize=10)
        ax.grid(axis="y", alpha=0.3)
        ax.set_axisbelow(True)

    fig.suptitle("PoseBusters failure waterfall — marginal % of poses lost per check "
                 "(Pareto order)" + _RMSD_FILTER_NOTE.replace("\n", " "), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


# ───────────────────────────── main ─────────────────────────────

def main() -> None:
    # How to run this file from a notebook cell:
    #   !python Scripts/Analysis/posebusters_validity_report.py
    #
    # Helpful optional flags:
    #   --csv               Path to PoseBusters filtered-results CSV.
    #   --out-dir           Folder where plots and summary CSVs are saved.
    #   --top-n-heatmap     Number of receptor-ligand pairs shown in heatmap.
    #   --top-n-bars        Number of receptor-ligand pairs shown in grouped bars.
    #   --no-split-equibind Treat EquiBind as one bucket (skip variant split).
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path,
                    default=Path("posebusters_results/benchmark/dock/"
                                 "posebusters_filtered_results.csv"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("posebusters_results/benchmark/dock/validity_report"))
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Restrict the report to the '<PDBID>_<LIG>' complex ids "
                         "listed in this file (one per line; '#' comments ok). "
                         "Pairs whose id is absent are dropped.")
    ap.add_argument("--top-n-heatmap", type=int, default=60)
    ap.add_argument("--top-n-bars", type=int, default=30)
    ap.add_argument("--split-equibind", action=argparse.BooleanOptionalAction,
                    default=True,
                    help="Split EquiBind into fpocket/p2rank/unguided x raw/smina/gnina "
                         "x clamp variants from the CSV provenance columns "
                         "(filename fallback) (default: on).")
    ap.add_argument("--best-equibind-only", action="store_true",
                    help="Keep only the single best-performing EquiBind variant "
                         "(highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%, "
                         "read from --oracle-summary) in all plots/CSVs, relabelled "
                         "'EquiBind*'. AutoDock/DiffDock are unaffected. For datasets without "
                         "a crystal (Orai), the default --oracle-summary borrows the benchmark "
                         "ranking.")
    ap.add_argument("--oracle-summary", type=Path, default=DEFAULT_ORACLE_SUMMARY,
                    help="oracle_summary.csv from posebusters_pose_comparison.py used "
                         "to pick the best EquiBind/DiffDock variant for "
                         "--best-equibind-only / --best-diffdock-only. Its per-variant "
                         "'*_all_variants.csv' sibling is used automatically when present "
                         "(required to rank the split smina/gnina variants; the collapsed "
                         "file lists only one row per tool) (default: %(default)s).")
    ap.add_argument("--best-diffdock-only", action="store_true",
                    help="Keep only the single best-performing DiffDock optimizer "
                         "variant (raw/smina/gnina, highest PB-Valid AND RMSD ≤ 2 Å = "
                         "oracle_pb_valid_and_rmsd2_%%, read from --oracle-summary) in all "
                         "plots/CSVs, relabelled 'DiffDock*'. AutoDock/EquiBind are unaffected.")
    ap.add_argument("--diffdock-variant", choices=["raw", "original", "smina", "gnina"],
                    default=None,
                    help="Keep only this DiffDock optimizer variant in all plots/CSVs, "
                         "relabelled 'DiffDock*'. Explicit, crystal-free alternative to "
                         "--best-diffdock-only (names the variant instead of ranking by "
                         "PB-Valid AND RMSD ≤ 2 Å). AutoDock/EquiBind are unaffected.")
    ap.add_argument("--equibind-variant", default=None, metavar="SPEC",
                    help="Keep only the EquiBind variant matching SPEC (tokens split on "
                         "'_' or '/', matched against the pocket/refine/clamp axes, e.g. "
                         "'unguided_gnina' or 'unguided/gnina/clampOFF'); relabelled "
                         "'EquiBind*' when it resolves to one variant. Explicit, "
                         "crystal-free counterpart to --best-equibind-only. "
                         "AutoDock/DiffDock are unaffected.")
    ap.add_argument("--max-rmsd", type=float, default=None, metavar="A",
                    help="Keep only poses whose RMSD-to-crystal is ≤ this many Å "
                         "(e.g. 5.0) before scoring/plotting. Off by default. RMSD is "
                         "joined in from --per-pose-metrics; poses with no RMSD value "
                         "are dropped (and counted). Not meaningful for crystal-free "
                         "sets (Orai).")
    ap.add_argument("--per-pose-metrics", type=Path, default=DEFAULT_PER_POSE_METRICS,
                    help="per_pose_metrics.csv from posebusters_pose_comparison.py — "
                         "the source of the per-pose RMSD used by --max-rmsd "
                         "(default: %(default)s).")
    ap.add_argument("--rmsd-column", choices=_RMSD_COLUMNS, default="rmsd",
                    help="Which per_pose_metrics RMSD column --max-rmsd / --rmsd-sweep-max "
                         "use (default: %(default)s — symmetry-corrected, no "
                         "superposition, matching the oracle RMSD ≤ 2 Å metric).")
    ap.add_argument("--collapse-plots-only", action="store_true",
                    help="Keep every variant in the CSV tables and the printed summary, "
                         "but collapse each tool to a single variant in the main "
                         "comparison figures (01-05 + the RMSD sweep), shown as "
                         "'DiffDock*'/'EquiBind*'. Pairs with --best-*-only / "
                         "--diffdock-variant / --equibind-variant to choose which variant "
                         "is kept; with none of those it keeps each tool's best variant on "
                         "PB-Valid AND RMSD ≤ 2 Å (oracle_pb_valid_and_rmsd2_%%). Also emits a "
                         "per-tool variant comparison figure for every tool that has >1 variant.")
    ap.add_argument("--rmsd-sweep-max", type=float, default=None, metavar="A",
                    help="Also draw a sweep figure (08) of how PB-validity accumulates as "
                         "the RMSD-to-crystal cutoff is relaxed from 0 to this many Å "
                         "(e.g. 10). Uses --per-pose-metrics / --rmsd-column; skipped "
                         "(with a note) for crystal-free sets (Orai).")
    ap.add_argument("--rmsd-sweep-steps", type=int, default=60,
                    help="Number of cutoff samples in the --rmsd-sweep-max curve "
                         "(default: %(default)s).")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = load_and_score(args.csv, split_equibind=args.split_equibind)

    # Restrict to the official benchmark-set ids (the complex id is the 'protein'
    # column, which equals '<PDBID>_<LIG>' for the benchmark staging).
    if args.ids_file:
        allowed = _load_allowed_ids(args.ids_file)
        n0 = df["protein"].nunique()
        df = df[df["protein"].astype(str).isin(allowed)].copy()
        print(f"Restricted to {len(allowed)} ids from {args.ids_file.name}: "
              f"{n0} → {df['protein'].nunique()} receptor-ligand pairs")

    # Optionally keep only poses within --max-rmsd Å of the crystal. Done before the
    # best-variant selection and all scoring so every downstream CSV/figure reflects
    # the filtered set (which is also flagged in each figure title via _RMSD_FILTER_NOTE).
    generated_counts: pd.Series | None = None
    generated_valid_counts: pd.Series | None = None
    if args.max_rmsd is not None:
        global _RMSD_FILTER_NOTE
        # Poses produced per method *before* the RMSD filter — count and PB-valid count
        # — so the summary can show validity over ALL generated poses next to the
        # within-cutoff subset. Captured now (after any ids-file restriction, before the
        # filter): the docking_method labels are already final (equibind/diffdock splits
        # applied in load_and_score), pb_valid is already scored on the full set, and the
        # later best-*-only selection only drops rows, so reindexing to the final method
        # order in per_tool_summary is safe.
        generated_counts = df.groupby("docking_method").size()
        generated_valid_counts = df.groupby("docking_method")["pb_valid"].sum()
        df = apply_rmsd_filter(df, args.per_pose_metrics, args.max_rmsd,
                               args.rmsd_column, args.out_dir)
        _RMSD_FILTER_NOTE = (f"\n(poses filtered to {args.rmsd_column} "
                             f"≤ {args.max_rmsd:g} Å vs crystal)")

    # ── Variant selection ────────────────────────────────────────────────────
    # Two frames: ``df_full`` keeps every variant (drives the CSV tables, the printed
    # summary, the per-tool variant-comparison figures and the failure-rate table);
    # ``df_plot`` may be collapsed to one variant per tool for the main comparison
    # figures. Without --collapse-plots-only the two are the same object, so the
    # selectors shrink tables and plots alike — exactly the legacy behaviour.
    df_full = df

    def _select(frame: pd.DataFrame) -> pd.DataFrame:
        """Apply the requested variant selectors to *frame* (populating the
        'DiffDock*'/'EquiBind*' label overrides). With --collapse-plots-only and no
        explicit selector, falls back to keeping each tool's best variant on
        PB-Valid AND RMSD ≤ 2 Å."""
        dd_selected = bool(args.best_diffdock_only or args.diffdock_variant)
        eq_selected = bool(args.best_equibind_only or args.equibind_variant)
        if args.best_equibind_only:
            frame, best = select_best_equibind(frame, args.oracle_summary)
            if best:
                _LABEL_OVERRIDES[best] = "EquiBind*"
                print(f"best-equibind-only: '{best}' is the top EquiBind variant by "
                      "PB-Valid AND RMSD ≤ 2 Å — collapsing to it (shown as 'EquiBind*').")
        if args.best_diffdock_only:
            frame, best = select_best_diffdock(frame, args.oracle_summary)
            if best:
                _LABEL_OVERRIDES[best] = "DiffDock*"
                print(f"best-diffdock-only: '{best}' is the top DiffDock variant by "
                      "PB-Valid AND RMSD ≤ 2 Å — collapsing to it (shown as 'DiffDock*').")
        if args.diffdock_variant:
            frame, kept = select_diffdock_variant(frame, args.diffdock_variant)
            if kept:
                _LABEL_OVERRIDES[kept] = "DiffDock*"
                print(f"diffdock-variant: keeping only '{kept}' (shown as 'DiffDock*').")
        if args.equibind_variant:
            frame, kept = select_equibind_variant(frame, args.equibind_variant)
            if kept:
                _LABEL_OVERRIDES[kept] = "EquiBind*"
                print(f"equibind-variant: keeping only '{kept}' (shown as 'EquiBind*').")
        # Under --collapse-plots-only, collapse each tool that wasn't given an explicit
        # selector above to its best variant on PB-Valid AND RMSD ≤ 2 Å, so both tools
        # collapse by default (a lone --diffdock-variant no longer leaves EquiBind expanded).
        # Needs an oracle; crystal-free sets (Orai) should pin the variant explicitly instead.
        if args.collapse_plots_only and not dd_selected:
            frame, best_dd = select_best_diffdock(frame, args.oracle_summary)
            if best_dd:
                _LABEL_OVERRIDES[best_dd] = "DiffDock*"
                print(f"collapse-plots-only: keeping DiffDock's best on PB-Valid & RMSD ≤ 2 Å "
                      f"'{best_dd}' in the plots (shown as 'DiffDock*').")
        if args.collapse_plots_only and not eq_selected:
            frame, best_eq = select_best_equibind(frame, args.oracle_summary)
            if best_eq:
                _LABEL_OVERRIDES[best_eq] = "EquiBind*"
                print(f"collapse-plots-only: keeping EquiBind's best on PB-Valid & RMSD ≤ 2 Å "
                      f"'{best_eq}' in the plots (shown as 'EquiBind*').")
        return frame

    df_plot = _select(df_full)
    if not args.collapse_plots_only:
        # Legacy: no table/plot split — the selection (if any) applies everywhere.
        df_full = df_plot

    # Orderings / colours: the full set drives tables + variant/failure figures; the
    # (possibly collapsed) plot set drives the main comparison figures. In the common
    # non-collapse case the two frames are the same object, so reuse the computation.
    order_full = _method_order(df_full)
    colors_full = _method_colors(order_full)
    if df_plot is df_full:
        order_plot, colors_plot = order_full, colors_full
    else:
        order_plot = _method_order(df_plot)
        colors_plot = _method_colors(order_plot)

    # Near-native (≤ 2 Å RMSD-to-crystal) subset counts per method — poses and PB-valid
    # among them — added to the summary as rmsd2_poses / rmsd2_valid_poses /
    # rmsd2_valid_fraction. Reuses --max-rmsd's join if present, else joins
    # per_pose_metrics; (None, None) for crystal-free sets (Orai) → columns omitted.
    rmsd2_counts, rmsd2_valid_counts = _rmsd_le2_counts(
        df_full, args.per_pose_metrics, args.rmsd_column)

    # ── Tables (always the full, all-variants set) ───────────────────────────
    summary_full = per_tool_summary(df_full, order_full,
                                    generated_counts, generated_valid_counts,
                                    rmsd2_counts, rmsd2_valid_counts)
    valid_mat = per_pair_tool_matrix(df_full, "valid", order_full)
    total_mat = per_pair_tool_matrix(df_full, "total", order_full)
    frac_mat = per_pair_tool_matrix(df_full, "fraction", order_full)

    summary_full.to_csv(args.out_dir / "summary_per_tool.csv")
    valid_mat.to_csv(args.out_dir / "valid_poses_per_pair_tool.csv")
    total_mat.to_csv(args.out_dir / "total_poses_per_pair_tool.csv")
    frac_mat.round(3).to_csv(args.out_dir / "valid_fraction_per_pair_tool.csv")

    axis = equibind_axis_breakdown(df_full)
    if not axis.empty:
        axis.round(4).to_csv(args.out_dir / "equibind_axis_breakdown.csv")

    # Per-variant PB check failure-rate table (what fails, and by how much, per variant).
    fail_tbl = check_failure_table(df_full, order_full)
    if not fail_tbl.empty:
        fail_tbl.round(2).to_csv(args.out_dir / "check_failure_rates_per_variant.csv")

    # ── Main comparison figures (df_plot: collapsed under --collapse-plots-only) ──
    if df_plot is df_full:
        summary_plot, valid_mat_plot = summary_full, valid_mat
    else:
        summary_plot = per_tool_summary(df_plot, order_plot,
                                        generated_counts, generated_valid_counts,
                                        rmsd2_counts, rmsd2_valid_counts)
        valid_mat_plot = per_pair_tool_matrix(df_plot, "valid", order_plot)

    # ── Paired, per-complex validity statistics (aggregate poses to one value per
    # receptor-ligand pair FIRST). Always-on and cheap; wrapped so any failure degrades
    # to the current, test-free figures. Numbers are also written to validity_stats.json.
    dataset_label = args.out_dir.parent.parent.name if len(args.out_dir.parts) >= 2 \
        else str(args.out_dir.name)
    stats: dict | None = None
    try:
        stats = compute_validity_stats(df_plot, df_full, order_plot, order_full,
                                       dataset_label)
        (args.out_dir / "validity_stats.json").write_text(
            json.dumps(_json_safe(stats), indent=2, default=str))
        note = (" (exploratory — small n)" if stats.get("exploratory") else "")
        print(f"Stats: paired per-complex tests over n={stats.get('n_complete_units')} "
              f"complete units{note}; wrote validity_stats.json")
    except Exception as e:  # pragma: no cover - never crash the pipeline
        print(f"[warn] statistics step failed ({e}); figures drawn without tests.")
        stats = None

    figures: list[str] = []
    plot_per_tool(summary_plot, args.out_dir / "01_per_tool_validity.png", stats=stats)
    plot_heatmap(valid_mat_plot, args.out_dir / "02_valid_heatmap_top.png",
                 top_n=args.top_n_heatmap)
    plot_heatmap(valid_mat_plot, args.out_dir / "02b_valid_heatmap_bottom.png",
                 top_n=args.top_n_heatmap, worst=True)
    plot_grouped_bars(valid_mat_plot, args.out_dir / "03_valid_grouped_bars_top.png",
                      colors_plot, top_n=args.top_n_bars, stats=stats)
    plot_validity_distribution(df_plot,
                               args.out_dir / "04_valid_per_pair_distribution.png",
                               order_plot, stats=stats)
    plot_check_passrate(df_plot, args.out_dir / "05_per_check_passrate.png", order_plot,
                        stats=stats)
    figures += ["01_per_tool_validity.png", "02_valid_heatmap_top.png",
                "02b_valid_heatmap_bottom.png", "03_valid_grouped_bars_top.png",
                "04_valid_per_pair_distribution.png", "05_per_check_passrate.png"]

    # ── Per-tool variant-comparison figures (from the full, all-variants set) ──
    # Drawn for every tool that actually has >1 variant. Under --collapse-plots-only (or
    # a plain default run) df_full keeps every variant, so these carry the detail the
    # collapsed main figures drop; with a global selector and no collapse only one
    # variant is left, so they self-skip.
    if plot_equibind_variants(summary_full,
                              args.out_dir / "06_equibind_variant_validity.png",
                              colors_full, stats=stats):
        figures.append("06_equibind_variant_validity.png")
    if plot_tool_variant_comparison(
            summary_full, args.out_dir / "07_diffdock_variant_validity.png",
            colors_full, prefix="diffdock", tool_label="DiffDock",
            subtitle="(original vs smina- / gnina-optimised)", stats=stats):
        figures.append("07_diffdock_variant_validity.png")
    # Matrix view of 06/07: PB-validity (%) per variant as base method × optimizer.
    if plot_variant_validity_matrix(
            summary_full, args.out_dir / "07b_variant_validity_matrix.png", stats=stats):
        figures.append("07b_variant_validity_matrix.png")

    # ── RMSD-relaxation sweep (opt-in via --rmsd-sweep-max; crystal sets only) ──
    # Written as two separate figures: 08a cumulative count, 08b yield (%).
    if args.rmsd_sweep_max is not None:
        figures += plot_validity_vs_rmsd(
            df_plot, args.per_pose_metrics, args.rmsd_column,
            args.rmsd_sweep_max, args.rmsd_sweep_steps, order_plot, colors_plot,
            args.out_dir / "08a_validity_vs_rmsd_count.png",
            args.out_dir / "08b_validity_vs_rmsd_yield.png")

    # ── Failure waterfall for the top variant of each tool ──
    # AutoDock, best DiffDock, best EquiBind side by side — each variant's poses cascade
    # from 100% down through the checks they fail to the PB-valid %. "Best" is the same
    # PB-Valid AND RMSD ≤ 2 Å pick the collapsed comparison figures use (falls back to
    # highest PB-valid fraction for crystal-free sets), so the report names one 'best' variant.
    # Honour an explicit --diffdock-variant pin so the waterfall highlights the
    # SAME DiffDock variant as the collapsed figures (otherwise it re-ranks via the
    # oracle and can disagree — e.g. showing smina while the plots show gnina).
    def _diffdock_wf_variant():
        if args.diffdock_variant:
            opt = ("original" if str(args.diffdock_variant).lower() in ("raw", "original")
                   else str(args.diffdock_variant).lower())
            want = "diffdock" if opt == "original" else f"diffdock_{opt}"
            if (df_full["docking_method"].astype(str) == want).any():
                return want
        return _top_variant_for_waterfall(df_full, "diffdock", args.oracle_summary)
    top_variants = [v for v in (
        _top_variant_for_waterfall(df_full, "autodock", args.oracle_summary),
        _diffdock_wf_variant(),
        _top_variant_for_waterfall(df_full, "equibind", args.oracle_summary),
    ) if v]
    if top_variants and plot_failure_waterfall(
            df_full, top_variants, args.out_dir / "09_failure_waterfall_top.png"):
        figures.append("09_failure_waterfall_top.png")
        print(f"failure waterfall: top variant per tool = {', '.join(top_variants)}")

    # ── Console summary ──────────────────────────────────────────────────────
    print(f"Total poses scored        : {len(df_full):,}")
    print(f"Total PB-valid poses      : {int(df_full['pb_valid'].sum()):,} "
          f"({df_full['pb_valid'].mean() * 100:.2f}%)")
    print(f"Receptor-ligand pairs     : {df_full['pair'].nunique()}")
    print(f"Methods / variants found  : {len(order_full)} "
          f"({', '.join(order_full)})")
    if args.collapse_plots_only:
        print(f"Plots collapsed to        : {', '.join(order_plot)}")
    print()
    print("── Per-method summary (all variants) ─────────────")
    print(summary_full.to_string(float_format=lambda x: f"{x:.3f}"))
    if generated_counts is not None:
        print(f"  (generated_* = count / PB-valid / valid-fraction over ALL poses "
              f"produced; within_rmsd_poses + valid_poses + valid_fraction = the same "
              f"over the subset within {args.rmsd_column} ≤ {args.max_rmsd:g} Å)")
    if rmsd2_counts is not None:
        print(f"  (rmsd2_* = count / PB-valid / valid-fraction over poses within "
              f"{args.rmsd_column} ≤ {NEAR_NATIVE_RMSD_A:g} Å of the crystal — the "
              f"near-native subset)")
    if not fail_tbl.empty:
        print()
        print("── PB check failure rate (%) per variant "
              "(n_poses, survival %, then worst checks) ──")
        meta = ["n_poses", "survival_%"]
        check_cols = [c for c in fail_tbl.columns if c not in meta]
        worst = [c for c in check_cols if fail_tbl[c].max() > 0.0][:10]
        shown = fail_tbl[meta + worst].rename(columns={c: _prettify_check(c) for c in worst})
        fmt = {c: (lambda v: f"{v:6.1f}") for c in shown.columns}
        fmt["n_poses"] = lambda v: f"{int(v):>7d}"
        print(shown.to_string(formatters=fmt))
        print("  survival_% = poses passing ALL checks; "
              "full table (all checks): check_failure_rates_per_variant.csv")
    if not axis.empty:
        print()
        print("── EquiBind axis breakdown (pocket × refine × clamp) ──")
        print(axis.to_string(float_format=lambda x: f"{x:.3f}"))
    print()
    print(f"Wrote CSVs and {len(figures)} figures to: {args.out_dir.resolve()}")
    print(f"  figures: {', '.join(figures)}")


if __name__ == "__main__":
    main()
