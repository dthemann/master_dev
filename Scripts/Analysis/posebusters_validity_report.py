"""Visualize PoseBusters validity per docking tool & receptor-ligand pair.

Reads posebusters_filtered_results.csv and produces:
    1. Bar chart: total / valid pose counts per docking method.
    2. Heatmap: valid-pose counts per receptor-ligand x method — 02 easiest (top-N) and
       02b hardest (bottom-N) complexes.
    3. Cleveland dot plot: top-N receptor-ligand pairs comparing valid poses per method.
    4. Boxplot: distribution of per-pair valid poses per method.
    5. Heatmap: per-check pass rate per method (column labels carry each variant's n).
    6/7. Bar chart: PB-validity per variant, one figure per tool
         (EquiBind / AutoDock / DiffDock);
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

    # Show only the single best EquiBind variant (highest docking-success gate =
    # oracle_pb_valid_and_rmsd2_%, near-native AND PB-valid, read from --oracle-summary),
    # relabelled "EquiBind*". Because that gate contains PB-validity, the collapsed across-tool
    # comparison is reported descriptively (see validity_stats.txt). For crystal-free sets
    # (Orai) the default --oracle-summary borrows the benchmark ranking.
    python Scripts/Analysis/posebusters_validity_report.py --best-equibind-only

    # Keep only poses within 5 Å of the crystal (RMSD joined from per_pose_metrics.csv,
    # written by posebusters_pose_comparison.py — run that first). Every CSV/figure
    # then reflects the filtered set, and the figures are titled accordingly. Under
    # --max-rmsd the per-method summary reports validity over both populations side by
    # side: generated_poses / generated_valid_poses / generated_valid_fraction (ALL
    # poses produced) and within_rmsd_poses / valid_poses / valid_fraction (the kept
    # subset within the cutoff).
    python Scripts/Analysis/posebusters_validity_report.py --max-rmsd 5

    # Keep every variant in the CSV tables / printed summary but collapse EACH engine to a
    # single best variant ('AutoDock*'/'DiffDock*'/'EquiBind*') in the main comparison
    # figures (01-05 + the RMSD sweep) — symmetric, none pinned to raw. Also emits per-tool
    # variant-comparison figures (06 EquiBind, 06a/06b AutoDock Vina/Vinardo, 07 DiffDock).
    # Pair with a selector to choose the kept variant, or let it default to each engine's
    # best on the docking-success gate (near-native AND PB-valid → collapsed across-tool
    # comparison is descriptive, not inferential).
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
# counts as a valid pose. (run_posebusters lives one package over; add it to the path
# and import the constant + coercer rather than duplicating a check list.)
#
# The import is DEFERRED to the first real scoring call (_load_pb_constants), not done
# at module scope: importing run_posebusters pulls in the whole PoseBusters runtime,
# which only exists in the docking ('vina') conda env. Doing it at import time would
# make even `python posebusters_validity_report.py --help` crash in any other shell.
# So argparse/--help stay dependency-free and only an actual run needs PoseBusters.
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PB_DIR = _PROJECT_ROOT / "Scripts" / "Docking" / "Posebusters"
for _p in (str(_PROJECT_ROOT), str(_PB_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# Populated lazily by _load_pb_constants(). CRITICAL_CHECKS is the canonical dock-mode
# PB-Valid set (the full criterion, incl. cofactor/metal/water clash checks);
# MOLECULE_TEST_COLUMNS is the molecule-only subset for 'mol'-mode CSVs.
CANONICAL_TEST_COLUMNS: tuple = ()
MOLECULE_TEST_COLUMNS: tuple = ()
CRITICAL_CHECKS: list[str] = []
coerce_test_cols_to_bool = None          # set by _load_pb_constants
expected_test_columns = None             # set by _load_pb_constants


def _load_pb_constants() -> None:
    """Import the canonical PoseBusters columns + coercers from run_posebusters.

    Deferred so `--help` works without the PoseBusters runtime; raises a clear,
    actionable error the first time real scoring needs it and it is unavailable."""
    global CANONICAL_TEST_COLUMNS, MOLECULE_TEST_COLUMNS, CRITICAL_CHECKS
    global coerce_test_cols_to_bool, expected_test_columns
    if CRITICAL_CHECKS:
        return
    try:
        from run_posebusters import (
            CANONICAL_TEST_COLUMNS as _canon,
            MOLECULE_TEST_COLUMNS as _mol,
            coerce_test_cols_to_bool as _coerce,
            expected_test_columns as _expected,
        )
    except Exception as e:  # PoseBusters runtime absent in this interpreter
        raise ImportError(
            "posebusters_validity_report needs the PoseBusters runtime "
            "(run_posebusters) to score poses; run it in the 'vina' conda env "
            f"(e.g. `conda run -n vina python …`). Original import error: {e}"
        ) from e
    CANONICAL_TEST_COLUMNS = _canon
    MOLECULE_TEST_COLUMNS = _mol
    coerce_test_cols_to_bool = _coerce
    expected_test_columns = _expected
    CRITICAL_CHECKS = list(_canon)

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
# Blue family for AutoDock smina/GNINA optimizer variants (raw uses the base
# blue in _BASE_COLORS). Optimizer variants are separate methods: combining
# them with the Vina output would mix different geometries and ranking axes.
_AD_PALETTE = ["#6baed6", "#08519c", "#08306b"]
# Orange family for DiffDock smina/gnina optimizer variants (raw uses the base
# orange in _BASE_COLORS; these are the lighter/darker shades for the variants).
_DD_PALETTE = ["#ffbb78", "#d95f02", "#fdae6b", "#a63603"]
_BASE_COLORS = {
    "autodock": "#1f77b4",           # AutoDock Vina (blue)
    "autodock_vinardo": "#17becf",   # AutoDock Vinardo (cyan)
    "unidock": "#9467bd",            # Uni-Dock tiled (purple)
    "unidock2": "#8c564b",           # Uni-Dock2 (brown)
    "diffdock": "#ff7f0e",           # DiffDock (orange)
}

# Per-run display-label overrides (method key -> label). Populated in main() when
# --best-equibind-only is active, where the single retained EquiBind variant is
# shown as "EquiBind*". Honoured by _pretty_method, so every plot/axis label picks
# it up with no further changes.
_LABEL_OVERRIDES: dict[str, str] = {}

# Appended to every figure title when --max-rmsd is active, so a filtered report's
# figures can never be mistaken for the unfiltered one. Set once in main().
_RMSD_FILTER_NOTE: str = ""

# AutoDock alone has the additional CNN-refinement protocol.  The legacy
# ``gnina`` ID remains empirical minimization + CNN rescore; keeping the longer
# suffix first prevents ``gnina_refinement`` paths/methods from being mistaken
# for the legacy GNINA variant.  DiffDock and EquiBind retain their original
# raw/smina/gnina axes.
_AUTODOCK_OPTIMIZERS = frozenset({"smina", "gnina", "gnina_refinement"})
_AUTODOCK_OPTIMIZER_SUFFIXES = ("_gnina_refinement", "_smina", "_gnina")


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


def _method_sort_key(m: str):
    # Bucket 0 groups the Vina-family engines (AutoDock Vina then Vinardo, then
    # Uni-Dock / Uni-Dock2); bucket 1 = DiffDock, bucket 2 = EquiBind.
    if m.startswith("autodock_vinardo"):
        return (0, 10 + {"autodock_vinardo": 0, "autodock_vinardo_smina": 1,
                         "autodock_vinardo_gnina": 2,
                         "autodock_vinardo_gnina_refinement": 3}.get(m, 4), 0, 0)
    if m.startswith("autodock"):
        # Raw Vina first, then its smina/gnina post-optimization variants.
        return (0, {"autodock": 0, "autodock_smina": 1,
                    "autodock_gnina": 2,
                    "autodock_gnina_refinement": 3}.get(m, 4), 0, 0)
    if m == "unidock":
        return (0, 20, 0, 0)
    if m == "unidock2":
        return (0, 21, 0, 0)
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
    ad = [m for m in order if m.startswith("autodock") and m not in _BASE_COLORS]
    # DiffDock optimizer variants share an orange family (base orange = raw).
    dd = [m for m in order if m.startswith("diffdock") and m not in _BASE_COLORS]
    leftover = plt.cm.tab10.colors
    j = 0
    for m in order:
        if m in _BASE_COLORS:
            colors[m] = _BASE_COLORS[m]
        elif m in ad:
            colors[m] = _AD_PALETTE[ad.index(m) % len(_AD_PALETTE)]
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


_WARNED_UNKNOWN_AUTODOCK: set[str] = set()


def _warn_unknown_autodock(token: str) -> None:
    """One-time warning that an unrecognized AutoDock optimizer/method token was seen."""
    if token and token not in _WARNED_UNKNOWN_AUTODOCK:
        _WARNED_UNKNOWN_AUTODOCK.add(token)
        print(f"[warn] unrecognized AutoDock optimizer/method token '{token}' — bucketed "
              "as 'unknown' (kept out of the raw Vina denominator, not silently pooled).")


def _classify_autodock(row) -> str:
    """AutoDock optimizer axis, including the distinct GNINA refinement mode.

    Current PoseBusters exports carry an explicit ``optimizer`` field.  The
    path/method fallbacks keep older exports readable without ever pooling an
    ``optimized_*`` pose into the raw Vina bucket. An explicit-but-unrecognized
    optimizer value, or an ``autodock_<suffix>`` method whose suffix matches no known
    optimizer, is bucketed as ``"unknown"`` (with a one-time warning) rather than
    silently mislabelled ``"original"`` (raw Vina) — a new/typo'd optimizer must not
    inflate the raw denominator.
    """
    opt = (_col_value(row, "optimizer") or "").lower()
    if opt in _AUTODOCK_OPTIMIZERS:
        return opt
    if opt in ("original", "raw", "native", "none"):
        return "original"
    if opt:                                   # explicit optimizer field, unknown value
        _warn_unknown_autodock(opt)
        return "unknown"
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
    # No optimizer signal from the field or the pose path. If the method key still
    # carries an unrecognized suffix beyond its Vina/Vinardo scoring base, that is an
    # unknown variant — flag it rather than folding it into raw Vina.
    base = _autodock_scoring_base(method)
    suffix = method[len(base):].lstrip("_") if method.startswith(base) else ""
    if suffix:
        _warn_unknown_autodock(suffix)
        return "unknown"
    return "original"


def _autodock_scoring_base(method: str) -> str:
    """AutoDock scoring base of a method key — 'autodock' (Vina) or
    'autodock_vinardo' — with any optimizer suffix stripped, so the split is
    idempotent and never merges Vinardo into Vina."""
    mm = str(method).strip().lower()
    for suf in _AUTODOCK_OPTIMIZER_SUFFIXES:
        if mm.endswith(suf):
            mm = mm[: -len(suf)]
            break
    return "autodock_vinardo" if mm.startswith("autodock_vinardo") else "autodock"


def _apply_autodock_split(df: pd.DataFrame) -> pd.DataFrame:
    """Give every AutoDock optimizer its own stable method and axis label,
    preserving the scoring base (Vina vs Vinardo)."""
    df["ad_optimizer"] = pd.NA
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
    df.loc[ad_idx, "ad_optimizer"] = opts
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


def _resolve_check_schema(df: pd.DataFrame) -> tuple[list[str], str]:
    """Resolve which COMPLETE PoseBusters schema *df* carries → (checks, mode).

    A pose is PB-valid only when it passes the FULL mode-appropriate check set — the
    same rule run_posebusters enforces via require_test_columns. Accepting any nonempty
    subset (the old behaviour) would silently score a truncated CSV missing 21 of 22
    checks as '100% PB-valid'. So require every column of one known mode:
      dock : the 22 canonical checks (protein/cofactor/metal/water clashes included)
      mol  : the 12 molecule-only checks (a legitimate 'mol'-mode PoseBusters run)
    and raise if the CSV matches neither completely."""
    dock_missing = [c for c in CANONICAL_TEST_COLUMNS if c not in df.columns]
    if not dock_missing:
        return list(CANONICAL_TEST_COLUMNS), "dock"
    mol_missing = [c for c in MOLECULE_TEST_COLUMNS if c not in df.columns]
    if not mol_missing:
        print("[info] CSV carries the molecule-only PB schema (no protein/cofactor "
              "clash checks) — scoring molecule-valid, not full dock PB-valid.")
        return list(MOLECULE_TEST_COLUMNS), "mol"
    present = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
    raise ValueError(
        "PoseBusters CSV matches neither the complete dock schema "
        f"({len(dock_missing)} of {len(CANONICAL_TEST_COLUMNS)} checks missing, e.g. "
        f"{dock_missing[:5]}) nor the complete molecule schema ({len(mol_missing)} of "
        f"{len(MOLECULE_TEST_COLUMNS)} missing). Present canonical checks: {present}. "
        "A partial schema must never be scored as PB-valid — re-run PoseBusters so "
        "every mode-specific check column is present, or pass a complete CSV."
    )


def load_and_score(csv_path: Path, split_equibind: bool = True) -> pd.DataFrame:
    _load_pb_constants()
    df = pd.read_csv(csv_path, low_memory=False)
    checks, _mode = _resolve_check_schema(df)
    df["pb_valid"] = _bool_checks(df, checks).all(axis=1)
    df["pair"] = df["protein"].astype(str) + " / " + df["ligand"].astype(str)
    df["docking_method"] = df["docking_method"].astype(str).str.lower()
    df["method_raw"] = df["docking_method"]
    if split_equibind:
        df = _apply_equibind_split(df)
    # AutoDock's post-optimization rows have a distinct geometry and, for gnina,
    # a distinct optimized ranking. Keep them out of the native-Vina denominator.
    df = _apply_autodock_split(df)
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


def _attach_rmsd_to_crystal(df: pd.DataFrame, metrics: pd.DataFrame,
                            rmsd_col: str) -> pd.DataFrame:
    """Return a copy of *df* with a numeric ``rmsd_to_crystal`` column joined from the
    already-loaded *metrics* frame (per_pose_metrics.csv).

    The join is method-AWARE whenever *metrics* carries a per-pose ``method`` column and
    *df* carries ``docking_method``: it keys on (protein, ligand, pose_name, method).
    This matters because the full-protein per_pose_metrics reuse the SAME ``…_poseN``
    pose_name across autodock / autodock_vinardo / unidock2 (with genuinely different
    RMSDs) — a method-blind (protein, ligand, pose_name) key collapses those rows via
    drop_duplicates so every Vina-family variant inherits the first method's RMSD. Both
    files label the variants identically (load_and_score has already relabelled
    ``docking_method`` to match ``method``), so the method-aware key lines up 1:1.

    Rows whose (protein, ligand, pose_name, method) has no metrics match fall back to the
    method-blind key (first row per pose), so pose sets that never reused a pose_name
    across methods — the pocket benchmark, and the DiffDock/EquiBind families here — join
    exactly as before. Unmatched poses end up NaN and are handled by the caller."""
    key = ["protein", "ligand", "pose_name"]
    d = df.copy()
    for k in key:
        d[k] = d[k].astype(str)
    d["rmsd_to_crystal"] = np.nan
    matched = pd.Series(False, index=d.index)

    if "method" in metrics.columns and "docking_method" in d.columns:
        mkey = key + ["docking_method"]
        m = metrics[key + ["method", rmsd_col]].copy()
        for k in key:
            m[k] = m[k].astype(str)
        m["docking_method"] = m["method"].astype(str).str.lower()
        m = m.drop_duplicates(subset=mkey)[mkey + [rmsd_col]]
        left = d[mkey].copy()
        left["docking_method"] = left["docking_method"].astype(str).str.lower()
        j = left.merge(m.rename(columns={rmsd_col: "_r"}), on=mkey,
                       how="left", validate="m:1")
        d["rmsd_to_crystal"] = pd.to_numeric(j["_r"].to_numpy(), errors="coerce")
        matched = d["rmsd_to_crystal"].notna()

    need = ~matched
    if need.any():
        m2 = metrics[key + [rmsd_col]].copy()
        for k in key:
            m2[k] = m2[k].astype(str)
        # Method-blind fallback for poses the method-aware key didn't match — but ONLY
        # for pose keys that are UNAMBIGUOUS in the metrics file (i.e. (protein, ligand,
        # pose_name) occurs under a single method). The full-protein metrics reuse the
        # same pose_name across autodock / autodock_vinardo / unidock2 with genuinely
        # different RMSDs; a method-blind "first row per pose" pick there is order-
        # dependent and would assign another method's RMSD. Those ambiguous keys are left
        # NaN (dropped and counted) rather than silently mis-joined.
        grp_size = m2.groupby(key, sort=False)[rmsd_col].transform("size")
        m2_unique = m2[grp_size == 1]
        if not m2_unique.empty:
            j2 = d.loc[need, key].merge(m2_unique.rename(columns={rmsd_col: "_r"}),
                                        on=key, how="left", validate="m:1")
            d.loc[need, "rmsd_to_crystal"] = pd.to_numeric(j2["_r"].to_numpy(),
                                                           errors="coerce")
        # A method present in the poses but absent from the metrics file can only ever
        # reach this fallback — the pocket-default per_pose_metrics lacks autodock_vinardo
        # / unidock2, so a full-protein report run without an explicit --per-pose-metrics
        # would otherwise inherit AutoDock's RMSD for those rows. Warn loudly and leave the
        # ambiguous ones unmatched instead.
        if "method" in metrics.columns and "docking_method" in d.columns:
            df_methods = set(d["docking_method"].astype(str).str.lower().unique())
            met_methods = set(metrics["method"].astype(str).str.lower().unique())
            missing_methods = sorted(df_methods - met_methods)
            if missing_methods:
                still_nan = int(d.loc[need, "rmsd_to_crystal"].isna().sum())
                print(f"  [rmsd-join] {len(missing_methods)} pose method(s) have NO rows "
                      f"in the metrics file (e.g. {missing_methods[:4]}); their RMSD can "
                      f"only come from the method-blind fallback, and {still_nan} pose(s) "
                      "with a reused pose_name are left unmatched rather than assigned "
                      "another method's RMSD. Pass the matching --per-pose-metrics.")

    d["rmsd_to_crystal"] = pd.to_numeric(d["rmsd_to_crystal"], errors="coerce")
    return d


def apply_rmsd_filter(df: pd.DataFrame, metrics_csv: Path, max_rmsd: float,
                      rmsd_col: str, out_dir: Path) -> pd.DataFrame:
    """Keep only poses whose RMSD-to-crystal (``rmsd_col``) is ≤ ``max_rmsd`` Å.

    The PoseBusters filtered-results CSV has no RMSD, so the per-pose RMSD-to-crystal
    is joined in from *metrics_csv* (per_pose_metrics.csv, written by
    posebusters_pose_comparison.py) via _attach_rmsd_to_crystal — a method-aware
    (protein, ligand, pose_name, method) join, because the full-protein metrics reuse
    the same pose_name across autodock / autodock_vinardo / unidock2.
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

    n_in = len(df)
    d = _attach_rmsd_to_crystal(df, metrics, rmsd_col)

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


# The metric the best-variant selectors rank on: the docking-success GATE
# ``oracle_pb_valid_and_rmsd2_%`` — a pose must be near-native (RMSD ≤ 2 Å) AND PB-valid.
# On a single-site crystal RMSD ≤ 2 Å already implies the correct binding site (centroid
# ≤ 4 Å), so this is the full site → placement → validity 3-step success gate for every
# crystal-bearing dataset here. "Best" therefore means the most USABLE variant, not merely
# the one closest to the crystal (raw DiffDock is often near-native but clashing). Falls
# back to the RMSD-only rate for summaries predating the combined column.
#
# IMPORTANT: because this gate CONTAINS PB-validity, a variant selected on it and then
# compared on PB-validity is post-selection. The collapsed across-tool comparison is
# therefore reported as DESCRIPTIVE, not inferential (compute_validity_stats /
# write_stats_sidecar flag it via ``post_selection``). The per-variant tests, run on the
# full unselected set, are unaffected and stay inferential.
BEST_VARIANT_METRIC = "oracle_pb_valid_and_rmsd2_%"
BEST_VARIANT_METRIC_FALLBACK = "oracle_rmsd_le_2.0A_%"


def _variant_ranking_scores(oracle_csv: Path, tag: str) -> tuple[pd.Series | None, str | None]:
    """Per-variant ranking scores and the metric column used, read from *oracle_csv*.

    Ranks on the docking-success gate ``oracle_pb_valid_and_rmsd2_%`` (near-native AND
    PB-valid), falling back to ``oracle_rmsd_le_2.0A_%`` (RMSD-only) for summaries predating
    the combined column. Returns (scores, metric) or (None, None) — with a ``[tag]`` note —
    when the file carries neither column. The caller has already resolved *oracle_csv* to its
    ``*_all_variants.csv`` sibling (see _resolve_variant_oracle)."""
    osum = pd.read_csv(oracle_csv, index_col=0)
    col = (BEST_VARIANT_METRIC if BEST_VARIANT_METRIC in osum.columns
           else BEST_VARIANT_METRIC_FALLBACK if BEST_VARIANT_METRIC_FALLBACK in osum.columns
           else None)
    if col is None:
        print(f"  [{tag}] neither '{BEST_VARIANT_METRIC}' nor "
              f"'{BEST_VARIANT_METRIC_FALLBACK}' in {oracle_csv} — keeping all variants.")
        return None, None
    if col == BEST_VARIANT_METRIC_FALLBACK:
        print(f"  [{tag}] success-gate metric '{BEST_VARIANT_METRIC}' absent; falling back to "
              f"'{BEST_VARIANT_METRIC_FALLBACK}' (RMSD ≤ 2 Å only) for variant selection.")
    return pd.to_numeric(osum[col], errors="coerce"), col


def _select_best_family(df: pd.DataFrame, oracle_csv: Path, family_fn,
                        tag: str, keep_msg: str) -> tuple[pd.DataFrame, str | None]:
    """Keep all rows outside a tool family + only the single best in-family variant.

    "Best" = the in-family variant with the highest docking-success-gate rate (near-native
    AND PB-valid; see BEST_VARIANT_METRIC). *family_fn* maps a method key to True when it
    belongs to the family; the per-variant ``*_all_variants.csv`` sibling of *oracle_csv* is
    used when present (see _resolve_variant_oracle) so split smina/gnina/pocket labels are
    actually ranked rather than collapsed to one row. Returns (filtered_df, best_key), or the
    inputs unchanged with None when the family is absent, the oracle is missing/unreadable, or
    no present variant has a score."""
    methods = df["docking_method"].astype(str)
    mask = methods.map(family_fn).astype(bool)
    if not mask.any():
        return df, None
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [{tag}] oracle summary not found at {oracle_csv} — {keep_msg}")
        return df, None
    scores, _ = _variant_ranking_scores(oracle_csv, tag)
    if scores is None:
        return df, None
    present = sorted(methods[mask].unique())
    cand = {m: float(scores[m]) for m in present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print(f"  [{tag}] none of the present variants have a score in {oracle_csv} — {keep_msg}")
        return df, None
    best = max(cand, key=cand.get)
    keep = (~mask) | (methods == best)
    return df[keep].reset_index(drop=True), best


def select_best_equibind(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind rows + only the single best EquiBind variant (docking-success gate:
    near-native AND PB-valid). For crystal-free sets (Orai) point --oracle-summary at the
    benchmark summary to borrow its ranking. AutoDock/DiffDock rows are always kept. Returns
    (filtered_df, best_variant_key) — best None / df unchanged when no EquiBind variant is
    present, the summary is missing, or no variant has a score."""
    return _select_best_family(df, oracle_csv, lambda m: m.startswith("equibind"),
                               "best-equibind-only", "keeping all EquiBind variants.")


def select_best_diffdock(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-DiffDock rows + only the single best DiffDock optimizer variant (docking-
    success gate: near-native AND PB-valid). AutoDock/EquiBind rows are always kept."""
    return _select_best_family(df, oracle_csv, lambda m: m.startswith("diffdock"),
                               "best-diffdock-only", "keeping all DiffDock variants.")


def select_best_autodock(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-AutoDock-Vina rows + only the single best AutoDock Vina variant (raw vs
    smina-/gnina-optimised / GNINA-refinement) by the docking-success gate. Symmetric with the
    DiffDock/EquiBind selectors so AutoDock is not pinned to raw in the collapsed headline;
    Vinardo, Uni-Dock and other engines are untouched."""
    return _select_best_family(
        df, oracle_csv,
        lambda m: m.startswith("autodock") and _autodock_scoring_base(m) == "autodock",
        "best-autodock-only", "keeping all AutoDock Vina variants.")


def select_best_autodock_vinardo(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-Vinardo rows + only the single best AutoDock Vinardo variant by the docking-
    success gate (symmetric counterpart to select_best_autodock for the Vinardo family)."""
    return _select_best_family(df, oracle_csv, lambda m: m.startswith("autodock_vinardo"),
                               "best-autodock-vinardo-only", "keeping all AutoDock Vinardo variants.")


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


def _engine_family(m: str) -> str:
    """Coarse docking-engine bucket for a method key — one 'engine' per headline bar."""
    m = str(m)
    if m.startswith("autodock_vinardo"):
        return "autodock_vinardo"
    if m.startswith("autodock"):
        return "autodock"
    if m.startswith("diffdock"):
        return "diffdock"
    if m.startswith("equibind"):
        return "equibind"
    return m                       # unidock, unidock2, or any single-variant engine


def _select_presentation_tools(df: pd.DataFrame) -> pd.DataFrame:
    """Reduce a collapsed plot frame to ONE variant per docking engine.

    The caller (--collapse-plots-only) has already collapsed each multi-variant family to
    its best variant on the docking-success gate (near-native AND PB-valid) and starred it
    in _LABEL_OVERRIDES. Here we keep exactly one representative per
    engine so the headline compares each tool at its selected best — symmetric across
    AutoDock (Vina and Vinardo), Uni-Dock, Uni-Dock2, DiffDock and EquiBind, with no engine
    pinned to raw and none silently dropped. For a family with no oracle-based selection
    (crystal-free set / no explicit pin) fall back to its raw base variant; single-variant
    engines pass through unchanged."""
    methods = df["docking_method"].astype(str)
    present = list(dict.fromkeys(methods))
    starred = set(_LABEL_OVERRIDES)
    fams: dict[str, list[str]] = {}
    for m in present:
        fams.setdefault(_engine_family(m), []).append(m)
    keep_labels: set[str] = set()
    for fam, variants in fams.items():
        fam_starred = [v for v in variants if v in starred]
        if fam_starred:
            keep_labels.update(fam_starred)          # oracle-selected best for this engine
        elif len(variants) == 1:
            keep_labels.add(variants[0])             # single-variant engine — pass through
        elif fam in variants:
            keep_labels.add(fam)                     # unselected multi-variant → raw base
        else:
            keep_labels.update(variants)             # no base label present → keep all
    return df[methods.isin(keep_labels)].reset_index(drop=True)


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
    agg = df.groupby("docking_method").agg(
        analyzed_poses=("pb_valid", "size"), valid_poses=("pb_valid", "sum"))

    if generated_counts is not None:
        # Under --max-rmsd, keep methods that generated poses but had ZERO within the
        # cutoff (generated_poses = N, within_rmsd_poses = 0) instead of letting them
        # vanish from the summary because the post-filter frame no longer contains them.
        idx = [m for m in order if (m in agg.index) or (m in generated_counts.index)]
        for m in generated_counts.index:
            if m not in idx:
                idx.append(m)
        grp = agg.reindex(idx)
        grp["analyzed_poses"] = grp["analyzed_poses"].fillna(0)
        grp["valid_poses"] = grp["valid_poses"].fillna(0)
        grp = grp.rename(columns={"analyzed_poses": "within_rmsd_poses"})
        gen = generated_counts.reindex(grp.index).fillna(0).astype("int64")
        grp.insert(0, "generated_poses", gen)
        if generated_valid_counts is not None:
            gval = generated_valid_counts.reindex(grp.index).fillna(0).astype("int64")
            grp.insert(1, "generated_valid_poses", gval)
            grp.insert(2, "generated_valid_fraction", gval / gen.replace(0, np.nan))
        analyzed_col = "within_rmsd_poses"
    else:
        grp = agg.reindex(order).dropna(how="all").rename(
            columns={"analyzed_poses": "total_poses"})
        analyzed_col = "total_poses"
    grp[analyzed_col] = grp[analyzed_col].astype("int64")
    grp["valid_poses"] = grp["valid_poses"].astype("int64")
    # A method with 0 poses in the analysed population has an undefined validity fraction
    # (NaN), not 0/0 → 0; keep it NaN so it reads as "no data", not "all invalid".
    grp["valid_fraction"] = grp["valid_poses"] / grp[analyzed_col].replace(0, np.nan)

    # Per-method complex coverage: how many receptor-ligand pairs the method actually
    # produced a pose for, out of the dataset total. Makes a valid_fraction computed over
    # a subset of complexes impossible to misread as full-set performance — a method that
    # skipped the hard complexes shows pairs_covered < pairs_total here.
    if "pair" in df.columns:
        cov = df.groupby("docking_method")["pair"].nunique().reindex(grp.index).fillna(0)
        grp["pairs_covered"] = cov.astype("int64")
        grp["pairs_total"] = int(df["pair"].nunique())

    if rmsd2_counts is not None:
        # Near-native subset: poses within NEAR_NATIVE_RMSD_A Å of the crystal, and
        # PB-valid among them. An empty (but non-None) rmsd2_counts means "RMSD available,
        # zero within 2 Å" → all 0 (distinct from the crystal-free case, which omits these
        # columns entirely). Methods with none in-cutoff reindex to 0 → fraction NaN.
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


def _variant_paired_stats(df_full: pd.DataFrame, variants: list[str],
                          n_units: int) -> dict | None:
    """Paired per-complex stats across a tool's *variants* (an explicit method list).

    Returns {any_valid: paired_proportions, valid_count: paired_continuous,
    variants:[...], n} or None when the tool has <2 usable variants / too few units."""
    variants = _usable_methods(df_full, sorted(variants))
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
                           dataset_label: str, post_selection: bool = False) -> dict:
    """All paired, per-complex validity statistics for the report.

    Aggregates the correlated poses to one value per receptor-ligand ``pair`` FIRST,
    then runs the paired tests across tools (df_plot / order_plot) and across each
    tool's variants (df_full / order_full). Returns a JSON-ready dict; individual
    tests are wrapped so one failure never sinks the rest.

    *post_selection* True means each tool's variant in ``df_plot`` was chosen on the
    docking-success gate (which contains PB-validity), so the across-tool contrasts are
    post-selection and are reported DESCRIPTIVELY, not as unbiased inference (the sidecar
    banners them). The per-variant tests (on ``df_full``, unselected) stay inferential."""
    n_units = int(df_plot["pair"].nunique())
    # Drop structurally-sparse methods (e.g. the 1-complex legacy 'equibind_guided'
    # bucket) so listwise completeness across the remaining tools doesn't collapse n.
    usable_plot = _usable_methods(df_plot, order_plot)
    stats: dict = {"dataset": dataset_label, "unit_of_analysis": "receptor-ligand pair",
                   "n_units": n_units,
                   "exploratory": _exploratory(n_units),
                   "min_units_for_test": _MIN_UNITS_FOR_TEST,
                   "across_tools_post_selection": bool(post_selection),
                   "methods_tested": usable_plot,
                   "methods_dropped_sparse": [m for m in order_plot if m not in usable_plot]}

    # ── Across tools: per-complex valid-pose count + ≥1-valid (figs 01/03/04) ──
    cnt_full = _pair_valid_count_matrix(df_plot, usable_plot)
    cnt = cnt_full.dropna()
    stats["n_complete_units"] = int(len(cnt))
    # Complexes listwise-deleted from the paired tests because a tested method produced no
    # pose there. Surfaced so a shrink from n_units → n_complete_units (a method skipping
    # hard complexes) is visible in the sidecar, not silent.
    dropped_units = [p for p in cnt_full.index if p not in cnt.index]
    stats["incomplete_units_dropped"] = int(len(dropped_units))
    stats["incomplete_units"] = [str(p) for p in dropped_units[:50]]
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

    # ── Variants are paired (same poses minimized): raw vs smina vs gnina, per engine
    # family (06/06a/06b/07/07b). AutoDock Vina and Vinardo are tested separately so the
    # Vina figure is never mixed with Vinardo. ──
    methods_full = [str(m) for m in pd.unique(df_full["docking_method"])]
    variant_groups = {
        "equibind_variants": [m for m in methods_full if m.startswith("equibind")],
        "diffdock_variants": [m for m in methods_full if m.startswith("diffdock")],
        "autodock_variants": [m for m in methods_full if _engine_family(m) == "autodock"],
        "autodock_vinardo_variants":
            [m for m in methods_full if m.startswith("autodock_vinardo")],
    }
    for key, variants in variant_groups.items():
        try:
            res = _variant_paired_stats(df_full, variants, n_units)
            if res is not None:
                stats[key] = res
        except Exception as e:  # pragma: no cover
            print(f"  [stats] {key} stats failed: {e}")

    return stats


def write_stats_sidecar(stats: dict | None, out_dir: Path) -> None:
    """Write validity_stats.txt — the inferential results that used to be annotated onto the
    figure panels. Per the project's 'stats off panels into sidecars' convention, the figures
    stay purely descriptive (pose counts / pooled % / per-pair distributions) and every test
    lives here, each labelled with the exact per-complex estimand and unit of analysis it
    tests — so a per-complex '≥1-valid' or 'valid-count' test is never mistaken for the
    pose-level rate a panel plots. Numbers are also in validity_stats.json."""
    if not stats:
        return
    lines: list[str] = []
    w = lines.append

    def w_pairwise(pairwise, indent="    "):
        """One significant Holm-adjusted contrast per line (most significant first), with a
        tally of the non-significant ones so the block stays readable for many methods."""
        sig = [d for d in pairwise if d.get("star") in ("*", "**", "***")]
        ns = len(pairwise) - len(sig)
        sig.sort(key=lambda d: (d.get("p_holm", 1.0)
                                if d.get("p_holm") == d.get("p_holm") else 1.0))
        for d in sig:
            w(f"{indent}{_pretty_method(d['a'])} vs {_pretty_method(d['b'])}  "
              f"{d['star']}  (Holm p={su.fmt_p(d.get('p_holm'))})")
        if ns:
            w(f"{indent}(+{ns} further pairwise contrast(s) n.s.)")

    w("PoseBusters validity — inferential statistics (companion to the figures)")
    w(f"dataset: {stats.get('dataset', '?')}")
    w(f"unit of analysis: {stats.get('unit_of_analysis', 'receptor-ligand pair')}")
    w(f"n units (total): {stats.get('n_units', '?')}   |   "
      f"n complete units used in paired tests: {stats.get('n_complete_units', '?')}")
    if stats.get("exploratory"):
        w("NOTE: exploratory — small n (below the pre-registered threshold).")
    dropped = stats.get("incomplete_units_dropped", 0)
    if dropped:
        w(f"listwise-deleted complexes (a tested method produced no pose there): {dropped}")
        names = stats.get("incomplete_units", [])
        if names:
            w("  " + ", ".join(names) + ("  …" if dropped > len(names) else ""))
    if stats.get("methods_dropped_sparse"):
        w("methods dropped as structurally sparse (<50% coverage): "
          + ", ".join(stats["methods_dropped_sparse"]))
    w("")

    w("== Across tools (per-complex; poses aggregated to one value per pair first) ==")
    if stats.get("across_tools_post_selection"):
        w("*** DESCRIPTIVE ONLY — each tool's variant was selected on the docking-success")
        w("    gate (which contains PB-validity), so these across-tool contrasts are")
        w("    post-selection: read them as effect sizes, NOT as unbiased inference. The")
        w("    per-variant tests below (run on the full, unselected set) are unaffected. ***")
    vc = stats.get("valid_count")
    if vc:
        w("valid-pose COUNT per complex — Friedman + Kendall's W:")
        w("  " + (_fmt_cont_omnibus(vc) or "n/a"))
        w_pairwise(vc.get("pairwise", []), indent="  ")
    av = stats.get("any_valid")
    if av:
        w("≥1 PB-valid pose per complex — Cochran's Q:")
        w("  " + (_fmt_prop_omnibus(av) or "n/a"))
        w_pairwise(av.get("pairwise", []), indent="  ")
    if not (vc or av):
        w(stats.get("note_tools", "  (no across-tool test available)"))

    if stats.get("per_check") is not None:
        w("")
        w("== Per-check family (per-complex '≥1 pose passes check', Cochran Q, BH over checks) ==")
        sig = stats.get("per_check_significant_bh", [])
        w("checks where tools differ (BH q<0.05): " + (", ".join(sig) if sig else "none"))

    for key, lab in (("equibind_variants", "EquiBind variants"),
                     ("diffdock_variants", "DiffDock variants"),
                     ("autodock_variants", "AutoDock Vina variants"),
                     ("autodock_vinardo_variants", "AutoDock Vinardo variants")):
        v = stats.get(key)
        if not v:
            continue
        w("")
        w(f"== {lab} (paired, per-complex) ==")
        w(f"variants: {', '.join(v.get('variants', []))}   (n={v.get('n', '?')})")
        if v.get("note"):
            w("  " + str(v["note"]))
        av2 = v.get("any_valid")
        if av2:
            c = _fmt_prop_omnibus(av2)
            if c:
                w("  ≥1-valid: " + c)
            w_pairwise(av2.get("pairwise", []))
        vc2 = v.get("valid_count")
        if vc2:
            c = _fmt_cont_omnibus(vc2)
            if c:
                w("  valid-count: " + c)
            w_pairwise(vc2.get("pairwise", []))

    w("")
    w("Estimand note: the figure panels are descriptive only. The tests above are paired,")
    w("per-complex (one value per receptor-ligand pair before testing). Do not read a bar's")
    w("height as the tested quantity — e.g. fig 05 plots pose-level per-check pass rates,")
    w("while its test is the per-complex '≥1 pose passes' contrast reported here.")
    (out_dir / "validity_stats.txt").write_text("\n".join(lines) + "\n")


# ───────────────────────────── plots ─────────────────────────────

def plot_per_tool(summary: pd.DataFrame, out: Path,
                  title_prefix: str = "PoseBusters Benchmark") -> None:
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
    # Descriptive title only — the paired per-complex tests live in validity_stats.txt
    # (write_stats_sidecar), not on the panel.
    title = (f"{title_prefix} — Pose validity per docking method\n"
             "(valid = passes all canonical PB checks)" + _RMSD_FILTER_NOTE)
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
                      top_n: int = 30) -> None:
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
    # Descriptive title only — the whole-set paired ≥1-valid test lives in validity_stats.txt.
    title = (f"Valid poses per docking method — top {len(pairs)} receptor-ligand pairs"
             + _RMSD_FILTER_NOTE)
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_validity_distribution(df: pd.DataFrame, out: Path, order: list[str]) -> None:
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
    # Descriptive title only — the paired per-complex valid-count test lives in
    # validity_stats.txt (no significance brackets are drawn on the panel).
    title = ("Distribution of valid poses across receptor-ligand pairs\n"
             "(n = poses assessed per variant; whiskers/caps in red)"
             + _RMSD_FILTER_NOTE)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_check_passrate(df: pd.DataFrame, out: Path, order: list[str]) -> None:
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
    # Descriptive panel only: the heatmap shows POSE-LEVEL per-check pass rates. The
    # per-complex '≥1 pose passes' test (a different estimand) that used to ★-mark rows is
    # in validity_stats.txt, so a significant per-complex contrast is never conflated with
    # the pose-level rate plotted here.
    title = ("Per-check pass rate (%) per docking method (n = poses per variant)"
             + _RMSD_FILTER_NOTE)
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
    if str(m).startswith("autodock"):
        base = _autodock_scoring_base(m)
        opt = m[len(base):].lstrip("_") if str(m).startswith(base) else ""
        return {"": "native", "smina": "smina-opt", "gnina": "gnina-opt",
                "gnina_refinement": "GNINA CNN-refine"}.get(opt, opt or "native")
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
                                 prefix: str, tool_label: str, subtitle: str = "",
                                 select=None) -> bool:
    """Bar chart comparing PB-validity (%) across the variants of one docking tool.

    Selects *summary* rows via *select* (a method-key predicate) when given, else by
    ``str.startswith(prefix)``. The predicate lets the caller separate the AutoDock Vina
    family from Vinardo (a bare ``startswith("autodock")`` would sweep Vinardo into a
    "Vina" figure). Draws one bar per variant, annotated with the pass rate and valid/total
    counts. Returns False (nothing written) when the tool has < 2 variants. Labels bypass
    the 'best'-variant star so every optimizer/pocket flavour stays identifiable. Inferential
    tests live in validity_stats.txt, not on this panel."""
    if select is not None:
        sub = summary[[bool(select(str(m))) for m in summary.index]]
    else:
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
    # Descriptive title only — the paired raw-vs-smina-vs-gnina per-complex tests (which
    # validate gnina > smina) live in validity_stats.txt, not as brackets on this panel.
    title = f"{tool_label} — PB-validity per variant"
    if subtitle:
        title += "\n" + subtitle
    title += _RMSD_FILTER_NOTE
    ax.set_title(title, fontsize=10)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def plot_equibind_variants(summary: pd.DataFrame, out: Path, colors: dict) -> bool:
    """Valid fraction (%) per EquiBind variant. Returns False if < 2 variants.

    Thin wrapper over plot_tool_variant_comparison kept for its original call site."""
    return plot_tool_variant_comparison(
        summary, out, colors, prefix="equibind", tool_label="EquiBind",
        subtitle="(fpocket / p2rank / unguided, smina re-search, centroid clamp)")


def _variant_matrix_cell(m: str) -> tuple[str, str] | None:
    """Map a variant key to (row, column) for the base-method × optimizer validity
    matrix, or None for methods with no place in it. EquiBind and DiffDock share a
    raw/smina/gnina optimizer axis. AutoDock additionally gets a dedicated
    GNINA-refinement column; DiffDock and EquiBind keep their existing axes."""
    col_of = {"": "raw / original", "raw": "raw / original", "original": "raw / original",
              "smina": "smina", "gnina": "gnina",
              "gnina_refinement": "gnina refinement"}
    if m.startswith("autodock"):
        base = "autodock_vinardo" if m.startswith("autodock_vinardo") else "autodock"
        opt = m[len(base):].lstrip("_")
        row = "AutoDock Vinardo" if base.endswith("vinardo") else "AutoDock Vina"
        return row, col_of.get(opt, opt)
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


def plot_variant_validity_matrix(summary: pd.DataFrame, out: Path) -> bool:
    """Heatmap of PB-validity (%) per variant — the matrix view of the per-tool variant
    bar charts (figs 06/07). Rows are the base method (AutoDock, DiffDock, and EquiBind by
    pocket × clamp); columns are the optimizer axis (raw/original, smina, gnina,
    plus AutoDock-only GNINA refinement);
    each cell is that variant's PB-valid fraction. Returns False when fewer than two
    variant cells exist (nothing worth a matrix)."""
    # Accumulate valid/total pose COUNTS per (row, col) rather than assigning the
    # fraction directly: two distinct variants can legitimately map to the same cell
    # (e.g. the unsuffixed ``equibind_fpocket`` and the explicit-raw ``equibind_fpocket_raw``
    # both mean "fpocket · raw / original"). A plain assignment let the second silently
    # overwrite the first; pooling their counts gives the correct combined rate instead.
    count_col = _analyzed_count_col(summary)
    cell_valid: dict[str, dict[str, float]] = {}
    cell_total: dict[str, dict[str, float]] = {}
    for m in summary.index:
        rc = _variant_matrix_cell(str(m))
        if rc is None:
            continue
        r, c = rc
        cell_valid.setdefault(r, {}).setdefault(c, 0.0)
        cell_total.setdefault(r, {}).setdefault(c, 0.0)
        cell_valid[r][c] += float(summary.loc[m, "valid_poses"])
        cell_total[r][c] += float(summary.loc[m, count_col])
    cells = {r: {c: (cell_valid[r][c] / cell_total[r][c] * 100.0
                     if cell_total[r][c] else np.nan)
                 for c in cols} for r, cols in cell_valid.items()}
    if not cells:
        return False
    col_order = ["raw / original", "smina", "gnina", "gnina refinement"]
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
    # Descriptive title only — the per-complex omnibus across each engine's optimizer
    # variants lives in validity_stats.txt.
    title = ("PB-validity (%) per variant — base method × optimizer" + _RMSD_FILTER_NOTE)
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return True


def _join_rmsd_to_crystal(df: pd.DataFrame, metrics_csv: Path, rmsd_col: str,
                          label: str = "rmsd-sweep") -> pd.DataFrame | None:
    """Return *df* with a numeric ``rmsd_to_crystal`` column joined from *metrics_csv*
    (per_pose_metrics.csv) via _attach_rmsd_to_crystal — a method-aware (protein, ligand,
    pose_name, method) join — or None when the metrics file / column / join keys are
    unavailable. Unlike apply_rmsd_filter this never exits — callers (the RMSD sweep, the
    ≤2 Å summary columns) degrade gracefully (e.g. crystal-free Orai) instead of aborting
    the run. *label* only tags the skip message."""
    if not metrics_csv or not Path(metrics_csv).exists():
        print(f"  [{label}] per-pose metrics not found at {metrics_csv} — skipping.")
        return None
    metrics = pd.read_csv(metrics_csv, low_memory=False)
    key = ["protein", "ligand", "pose_name"]
    if rmsd_col not in metrics.columns or any(k not in metrics.columns for k in key):
        print(f"  [{label}] {metrics_csv} lacks '{rmsd_col}' or join keys {key} — skipping.")
        return None
    return _attach_rmsd_to_crystal(df, metrics, rmsd_col)


# Near-native RMSD-to-crystal cutoff for the rmsd2_* summary columns (the canonical
# docking-success threshold; matches the oracle_rmsd_le_2.0A_% / pb_rmsd_within_2A metrics).
NEAR_NATIVE_RMSD_A = 2.0


def _rmsd_le2_counts(df: pd.DataFrame, metrics_csv: Path, rmsd_col: str,
                     thresh: float = NEAR_NATIVE_RMSD_A):
    """Per-method (poses, PB-valid poses) within *thresh* Å of the crystal.

    Joins RMSD-to-crystal fresh from *metrics_csv*; *df* should be the UNFILTERED frame so
    the near-native (≤ *thresh* Å) subset is measured over every generated pose. Reusing a
    pre-filtered --max-rmsd column would make the ≤ 2 Å count a no-op when the filter cutoff
    is itself < 2 Å — the columns would be labelled ≤ 2 Å while actually holding the ≤ cutoff
    subset.

    Returns (poses, valid_poses) per docking_method, distinguishing the two zero cases the
    caller renders differently:
      * no RMSD data at all (crystal-free set / metrics missing) → (None, None): columns omitted.
      * RMSD available but nothing within *thresh* → EMPTY Series (not None): the caller then
        still reports rmsd2_poses = 0, so a genuine 0 % near-native is not shown like N/A."""
    d = _join_rmsd_to_crystal(df, metrics_csv, rmsd_col, label="rmsd≤2")
    if d is None:
        return None, None
    r = pd.to_numeric(d["rmsd_to_crystal"], errors="coerce")
    near = d[r <= thresh]
    if near.empty:
        print(f"  [rmsd≤2] RMSD available but no pose within {thresh:g} Å of the crystal "
              "— reporting rmsd2_* columns as 0 (not omitting them).")
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
    d_all = _join_rmsd_to_crystal(df, metrics_csv, rmsd_col)
    if d_all is None:
        return []
    # Yield denominator = each method's FULL pose count (BEFORE dropping poses with no
    # joined RMSD), so the curve is a true "% of the method's poses" (as the axis says),
    # not "% of the RMSD-matched poses" — poses absent from the metrics file must count
    # against the denominator, not silently shrink it.
    method_totals = d_all.groupby("docking_method").size()
    d = d_all[d_all["rmsd_to_crystal"].notna()]
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
        curves[method] = (cum, int(method_totals.get(method, len(sub))))

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


# Conditionally-produced outputs: written only for some configs (a variant figure that
# needs >1 variant, an --max-rmsd summary, the RMSD sweep, the failure waterfall, the
# stats sidecar). Cleared at the start of every run so a reused out-dir can never mix a
# stale artifact from a different config into the current report. Always-overwritten
# outputs (01–05, summary_per_tool.csv, validity_stats.json …) are not listed.
_CONDITIONAL_OUTPUTS = (
    "06_equibind_variant_validity.png",
    "06a_autodock_variant_validity.png",
    "06b_autodock_vinardo_variant_validity.png",
    "07_diffdock_variant_validity.png",
    "07b_variant_validity_matrix.png",
    "08a_validity_vs_rmsd_count.png",
    "08b_validity_vs_rmsd_yield.png",
    "09_failure_waterfall_top.png",
    "equibind_axis_breakdown.csv",
    "check_failure_rates_per_variant.csv",
    "rmsd_filter_summary.csv",
    "validity_stats.txt",
)


def _clear_stale_outputs(out_dir: Path) -> None:
    """Unlink this report's conditionally-produced artifacts in *out_dir* (if present)."""
    for name in _CONDITIONAL_OUTPUTS:
        f = out_dir / name
        if f.exists():
            f.unlink()


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
                         "(highest docking-success gate = oracle_pb_valid_and_rmsd2_%%, "
                         "near-native AND PB-valid, read from --oracle-summary) in all "
                         "plots/CSVs, relabelled 'EquiBind*'. AutoDock/DiffDock are unaffected. "
                         "For datasets without a crystal (Orai), the default --oracle-summary "
                         "borrows the benchmark ranking.")
    ap.add_argument("--oracle-summary", type=Path, default=DEFAULT_ORACLE_SUMMARY,
                    help="oracle_summary.csv from posebusters_pose_comparison.py used "
                         "to pick the best EquiBind/DiffDock variant for "
                         "--best-equibind-only / --best-diffdock-only. Its per-variant "
                         "'*_all_variants.csv' sibling is used automatically when present "
                         "(required to rank the split smina/gnina variants; the collapsed "
                         "file lists only one row per tool) (default: %(default)s).")
    ap.add_argument("--best-diffdock-only", action="store_true",
                    help="Keep only the single best-performing DiffDock optimizer "
                         "variant (raw/smina/gnina, highest docking-success gate = "
                         "oracle_pb_valid_and_rmsd2_%%, near-native AND PB-valid, read from "
                         "--oracle-summary) in all plots/CSVs, relabelled 'DiffDock*'. "
                         "AutoDock/EquiBind are unaffected.")
    ap.add_argument("--diffdock-variant", choices=["raw", "original", "smina", "gnina"],
                    default=None,
                    help="Keep only this DiffDock optimizer variant in all plots/CSVs, "
                         "relabelled 'DiffDock*'. Explicit, crystal-free alternative to "
                         "--best-diffdock-only (names the variant instead of ranking by "
                         "the docking-success gate). AutoDock/EquiBind are unaffected.")
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
                    help="Keep every variant in the CSV tables and the printed summary, but "
                         "show ONE variant per engine in the main comparison figures (01-05 + "
                         "the RMSD sweep), each at its best variant (relabelled "
                         "'AutoDock*'/'DiffDock*'/'EquiBind*'), symmetric across engines with "
                         "none pinned to raw. Pairs with --best-*-only / --diffdock-variant / "
                         "--equibind-variant to choose which variant is kept; with none of "
                         "those it selects each engine's best by the docking-success gate "
                         "(oracle_pb_valid_and_rmsd2_%%, near-native AND PB-valid) — so the "
                         "collapsed across-tool comparison is DESCRIPTIVE, not inferential "
                         "(see validity_stats.txt). Also emits a per-tool variant comparison "
                         "figure for every tool with >1 variant.")
    ap.add_argument("--rmsd-sweep-max", type=float, default=None, metavar="A",
                    help="Also draw a sweep figure (08) of how PB-validity accumulates as "
                         "the RMSD-to-crystal cutoff is relaxed from 0 to this many Å "
                         "(e.g. 10). Uses --per-pose-metrics / --rmsd-column; skipped "
                         "(with a note) for crystal-free sets (Orai).")
    ap.add_argument("--rmsd-sweep-steps", type=int, default=60,
                    help="Number of cutoff samples in the --rmsd-sweep-max curve "
                         "(default: %(default)s).")
    args = ap.parse_args()

    # Conflicting per-tool selectors must fail loudly, not resolve by silent apply-order
    # (an explicit --*-variant would otherwise be discarded after --best-*-only ran).
    if args.best_diffdock_only and args.diffdock_variant:
        ap.error("--best-diffdock-only and --diffdock-variant are mutually exclusive; "
                 "choose one DiffDock selector.")
    if args.best_equibind_only and args.equibind_variant:
        ap.error("--best-equibind-only and --equibind-variant are mutually exclusive; "
                 "choose one EquiBind selector.")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    # Reusing an output dir must not leave stale, conditionally-produced artifacts from a
    # different config (e.g. a prior --max-rmsd run's rmsd_filter_summary.csv, or a variant
    # figure that self-skips this run) masquerading as current output.
    _clear_stale_outputs(args.out_dir)
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
    # The pre-filter frame is retained so the near-native (≤ 2 Å) rmsd2_* columns are
    # always measured over ALL generated poses, never a < 2 Å --max-rmsd subset.
    df_unfiltered = df
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
        'DiffDock*'/'EquiBind*'/'AutoDock*' label overrides). With --collapse-plots-only and
        no explicit selector, collapses EVERY multi-variant engine to its best variant on the
        docking-success gate (near-native AND PB-valid) — the most usable variant per engine,
        symmetric (no engine pinned to raw). Because the gate contains PB-validity, the
        collapsed across-tool comparison is reported descriptively, not inferentially."""
        dd_selected = bool(args.best_diffdock_only or args.diffdock_variant)
        eq_selected = bool(args.best_equibind_only or args.equibind_variant)
        if args.best_equibind_only:
            frame, best = select_best_equibind(frame, args.oracle_summary)
            if best:
                _LABEL_OVERRIDES[best] = "EquiBind*"
                print(f"best-equibind-only: '{best}' is the top EquiBind variant by the "
                      "docking-success gate — collapsing to it (shown as 'EquiBind*').")
        if args.best_diffdock_only:
            frame, best = select_best_diffdock(frame, args.oracle_summary)
            if best:
                _LABEL_OVERRIDES[best] = "DiffDock*"
                print(f"best-diffdock-only: '{best}' is the top DiffDock variant by the "
                      "docking-success gate — collapsing to it (shown as 'DiffDock*').")
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
        # Under --collapse-plots-only, collapse EVERY engine that wasn't given an explicit
        # selector to its best variant on the docking-success gate (near-native AND PB-valid).
        # AutoDock (Vina and Vinardo) get the same treatment as DiffDock/EquiBind, so no engine
        # is pinned to raw. Needs an oracle; crystal-free sets (Orai) should pin the ML variants
        # explicitly (AutoDock has no crystal ranking there and stays raw).
        if args.collapse_plots_only:
            for selector, star, on in (
                (select_best_diffdock, "DiffDock*", not dd_selected),
                (select_best_equibind, "EquiBind*", not eq_selected),
                (select_best_autodock, "AutoDock*", True),
                (select_best_autodock_vinardo, "AutoDock Vinardo*", True),
            ):
                if not on:
                    continue
                frame, best = selector(frame, args.oracle_summary)
                if best:
                    _LABEL_OVERRIDES[best] = star
                    print(f"collapse-plots-only: keeping {star.rstrip('*')}'s best on the "
                          f"docking-success gate '{best}' in the plots (shown as '{star}').")
        return frame

    df_plot = _select(df_full)
    if args.collapse_plots_only:
        df_plot = _select_presentation_tools(df_plot)
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
        df_unfiltered, args.per_pose_metrics, args.rmsd_column)

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
    # Friendly title prefix for fig 01 — derived from the dataset, no longer hardcoded
    # "PoseBusters Benchmark" for every input (F13).
    dataset_title = {
        "benchmark": "PoseBusters Benchmark",
        "benchmark_full_protein_vina_scoring": "PoseBusters Benchmark (full protein)",
        "orai_benchmark": "Orai (benchmark ligands)",
        "orai": "Orai",
    }.get(dataset_label, dataset_label.replace("_", " ").title())
    # Across-tool contrasts are post-selection (descriptive) whenever a tool's variant was
    # data-driven-selected on the docking-success gate: --collapse-plots-only auto-collapse or
    # --best-*-only. Explicit --*-variant pins are pre-declared, but a collapse run still
    # auto-selects the other engines, so flag the whole across-tool block conservatively.
    post_selection = bool(args.collapse_plots_only or args.best_equibind_only
                          or args.best_diffdock_only)
    stats: dict | None = None
    try:
        stats = compute_validity_stats(df_plot, df_full, order_plot, order_full,
                                       dataset_label, post_selection=post_selection)
        (args.out_dir / "validity_stats.json").write_text(
            json.dumps(_json_safe(stats), indent=2, default=str))
        # Inferential tests live in a companion sidecar, not on the figure panels.
        write_stats_sidecar(stats, args.out_dir)
        note = (" (exploratory — small n)" if stats.get("exploratory") else "")
        print(f"Stats: paired per-complex tests over n={stats.get('n_complete_units')} "
              f"complete units{note}; wrote validity_stats.json + validity_stats.txt")
    except Exception as e:  # pragma: no cover - never crash the pipeline
        print(f"[warn] statistics step failed ({e}); figures drawn without the sidecar.")
        stats = None

    figures: list[str] = []
    plot_per_tool(summary_plot, args.out_dir / "01_per_tool_validity.png",
                  title_prefix=dataset_title)
    plot_heatmap(valid_mat_plot, args.out_dir / "02_valid_heatmap_top.png",
                 top_n=args.top_n_heatmap)
    plot_heatmap(valid_mat_plot, args.out_dir / "02b_valid_heatmap_bottom.png",
                 top_n=args.top_n_heatmap, worst=True)
    plot_grouped_bars(valid_mat_plot, args.out_dir / "03_valid_grouped_bars_top.png",
                      colors_plot, top_n=args.top_n_bars)
    plot_validity_distribution(df_plot,
                               args.out_dir / "04_valid_per_pair_distribution.png",
                               order_plot)
    plot_check_passrate(df_plot, args.out_dir / "05_per_check_passrate.png", order_plot)
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
                              colors_full):
        figures.append("06_equibind_variant_validity.png")
    # AutoDock Vina and Vinardo are SEPARATE figures — a bare "autodock" prefix would sweep
    # the Vinardo variants into a figure titled "AutoDock Vina" (F7). Each selects its own
    # scoring-base family explicitly.
    if plot_tool_variant_comparison(
            summary_full, args.out_dir / "06a_autodock_variant_validity.png",
            colors_full, prefix="autodock", tool_label="AutoDock Vina",
            subtitle="(native Vina vs smina- / gnina-optimised and re-ranked)",
            select=lambda m: m.startswith("autodock") and _autodock_scoring_base(m) == "autodock"):
        figures.append("06a_autodock_variant_validity.png")
    if plot_tool_variant_comparison(
            summary_full, args.out_dir / "06b_autodock_vinardo_variant_validity.png",
            colors_full, prefix="autodock_vinardo", tool_label="AutoDock Vinardo",
            subtitle="(native Vinardo vs smina- / gnina-optimised and re-ranked)",
            select=lambda m: m.startswith("autodock_vinardo")):
        figures.append("06b_autodock_vinardo_variant_validity.png")
    if plot_tool_variant_comparison(
            summary_full, args.out_dir / "07_diffdock_variant_validity.png",
            colors_full, prefix="diffdock", tool_label="DiffDock",
            subtitle="(original vs smina- / gnina-optimised)"):
        figures.append("07_diffdock_variant_validity.png")
    # Matrix view of 06/06a/06b/07: PB-validity (%) per variant as base method × optimizer.
    if plot_variant_validity_matrix(
            summary_full, args.out_dir / "07b_variant_validity_matrix.png"):
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
    # docking-success-gate pick the collapsed comparison figures use (falls back to highest
    # PB-valid fraction for crystal-free sets), so the report names one 'best' variant.
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
