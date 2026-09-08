"""Interaction-difference report from PandaMap fingerprints.

Consumes the CSVs written by ``run_pandamap.py``
(``pandamap_pose_summary.csv``, ``pandamap_interactions.csv``,
``crystal_interactions.csv``) and produces analyses that show **how the
protein-ligand interactions differ between docking methods** and, at a
**chemical level**, which bond/interaction types form for which receptor-ligand
chemotypes.

A. Cross-method interaction differences
   01 per-method interaction-type profile (mean count per pose)         [bar]
   02 interaction-type x method heatmap + total-interactions boxplot
   03 residue hot-spot heatmap (top residues x method, contact freq.)
   04 native-interaction recovery vs the crystal ligand (precision/recall/F1) [bar+CDF]
   05 fingerprint similarity (Jaccard): method<->method and method<->crystal

B. Chemical-level analyses (what bonds for what receptor-ligand types)
   06 interaction type vs ligand chemotype (halogen bonds vs #halogens;
      pi/cation-pi vs #aromatic rings; ionic/salt-bridge vs formal charge;
      H-bonds vs HBD+HBA) — RDKit descriptors from the crystal ligand
   07 residue-class preference per interaction type (aromatic / +charged /
      -charged / polar / hydrophobic)
   08 ligand-atom element vs interaction type

C. Deeper crystal-vs-pose comparison (extends 04/05; needs crystal fingerprint)
   09 native recall & precision resolved per interaction type (fig 04 split by type)
   10 matched / missed / spurious contact decomposition per method (stacked counts)
   11 per-residue native-contact recovery + per-residue spurious (hallucinated) contacts
   12 strict typed recall vs loose residue-contact recall (right residue, wrong bond)
   13 native-F1 as a pose selector: best-in-set (oracle) vs the tool's top-ranked pose
   14 geometry vs chemistry: RMSD-to-crystal vs native-interaction F1 per pose
      (joins the posebusters_pose_comparison per_pose_metrics.csv RMSD on
       (method, protein, ligand, pose_name); Spearman rho per method)
   15 per-pair crystal-vs-pose typed-contact overlays (best/worst/divergent case studies)

Proposals not built here (sketch for later): per-pair radar plots; ligand
chemotype clustering (nucleotide/cofactor vs drug-like) then interaction
profiles per cluster; 3D PandaMap / PyMOL exports.

Best-EquiBind filter (--best-equibind-only)
   Collapse the many EquiBind variants down to the single best-performing one
   (highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%, read from the
   posebusters_pose_comparison oracle_summary.csv via --oracle-summary) so every chart/CSV
   compares AutoDock,
   DiffDock, crystal and just one EquiBind series. The retained variant is
   relabelled "EquiBind*" in every legend/axis. Off by default.

Config (--config, optional)
   Reuses the same pandamap_config.yaml as run_pandamap.py. When given, it
   supplies --in-dir (from the config's output_dir), --benchmark-dir,
   best_equibind_only and oracle_summary, so one config can drive both stages.
   Any CLI flag overrides the matching config value.

Run under the vina env (needs rdkit):
   /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
       --in-dir pandamap_results/benchmark
   # focused report (best EquiBind variant only):
   /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
       --in-dir pandamap_results/benchmark --best-equibind-only
   # config-driven (mirrors run_pandamap.py): picks up output_dir + best_equibind_only:
   /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/pandamap_interaction_report.py \
       --config Scripts/Analysis/pandamap_config.yaml
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
if os.environ.get("MPLBACKEND") is None:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.patches
import seaborn as sns

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
# Shared (A),(B),(C)… panel labeller for multi-panel figures.
from pocket_comparison_report import _label_panels  # noqa: E402
import method_filter as mf  # noqa: E402  (shared single-point method exclusion)

# Shared, unit-tested statistical helpers (paired Friedman/Wilcoxon, Cochran's Q /
# McNemar, G-test of independence, Spearman, Wilson/bootstrap CIs, Holm/BH). Never
# reimplement a test — import them. Degrade to the current test-free figures if the
# module can't be imported (kept optional so the pipeline never hard-fails on stats).
try:
    import stats_utils as su  # noqa: E402
    _HAVE_STATS = True
except Exception as _stats_exc:  # pragma: no cover
    su = None
    _HAVE_STATS = False
    print(f"  [stats] stats_utils unavailable ({_stats_exc}) — figures drawn without tests.")

# 16 PandaMap interaction types (stable column order).
INTERACTION_TYPES = [
    "hydrogen_bonds", "carbon_pi", "pi_pi_stacking", "donor_pi", "amide_pi",
    "hydrophobic", "ionic", "halogen_bonds", "cation_pi", "metal_coordination",
    "salt_bridge", "covalent", "alkyl_pi", "attractive_charge", "pi_cation",
    "repulsion",
]

# Sodium's element symbol is the literal string "NA", which pandas reads as a
# missing value by default. That silently deletes every sodium coordination
# contact from ``lig_atom_element`` before figure 08 ever sees it (fig 08 drops
# null elements), so Na is the one coordinating metal missing from the panel.
# Read the residue-level frames with the default token list off and only the
# empty field treated as missing, so element symbols stay literal while a field
# PandaMap could not fill still arrives as NaN.
_CSV_NA_VALUES = [""]


def read_pandamap_csv(path, **kwargs) -> pd.DataFrame:
    """``read_csv`` that keeps element symbols literal (see ``_CSV_NA_VALUES``)."""
    return pd.read_csv(path, keep_default_na=False, na_values=_CSV_NA_VALUES, **kwargs)


# Readable axis labels for the interaction types (avoid raw snake_case on figures).
INTERACTION_LABELS = {
    "hydrogen_bonds": "Hydrogen bonds", "carbon_pi": "Carbon–π",
    "pi_pi_stacking": "π–π stacking", "donor_pi": "Donor–π", "amide_pi": "Amide–π",
    "hydrophobic": "Hydrophobic", "ionic": "Ionic", "halogen_bonds": "Halogen bonds",
    "cation_pi": "Cation–π", "metal_coordination": "Metal coordination",
    "salt_bridge": "Salt bridge", "covalent": "Covalent", "alkyl_pi": "Alkyl–π",
    "attractive_charge": "Attractive charge", "pi_cation": "π–cation",
    "repulsion": "Repulsion",
}


def pretty_itype(t: str) -> str:
    return INTERACTION_LABELS.get(t, t.replace("_", " ").capitalize())

# Residue classes (mirrors panda_maps.ipynb cell 19).
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}
POS_CHARGED = {"ARG", "LYS", "HIS"}
NEG_CHARGED = {"ASP", "GLU"}
POLAR = {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "HIS"}
HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY"}


def residue_class(resname: str) -> str:
    r = str(resname).upper()
    if r in POS_CHARGED and r != "HIS":
        return "+charged"
    if r in NEG_CHARGED:
        return "-charged"
    if r in AROMATIC:
        return "aromatic"
    if r in POLAR:
        return "polar"
    if r in HYDROPHOBIC:
        return "hydrophobic"
    return "other"


# ── method labelling (mirrors the other report scripts, kept local) ────────
_POCKET_ORDER = {"unguided": 0, "fpocket": 1, "p2rank": 2, "guided": 3}

# Per-run display-label overrides (method key -> label). Populated in main() when
# --best-equibind-only is active, where the single retained EquiBind variant is
# shown as "EquiBind*". Honoured by pretty_method, so every legend / tick label
# across the charts picks it up with no further changes.
_LABEL_OVERRIDES: dict[str, str] = {}


def _eq_tokens(method: str):
    pocket = refine = clamp = None
    for t in method.split("_")[1:]:
        if t in ("unguided", "fpocket", "p2rank", "guided"):
            pocket = t
        elif t in ("raw", "smina", "gnina"):
            refine = t
        elif t in ("clampON", "clampOFF"):
            clamp = t
    return pocket, refine, clamp


def pretty_method(m: str) -> str:
    if m in _LABEL_OVERRIDES:
        return _LABEL_OVERRIDES[m]
    # ADFRsuite/MGLTools exhaustiveness ladder. Without this the raw method key reaches the
    # axis of every figure, because the exact-match chain below only knows the Meeko keys.
    if m.startswith("autodock_mgltools"):
        base, opt = m, ""
        for suf in ("_gnina", "_smina"):
            if base.endswith(suf):
                base, opt = base[: -len(suf)], f" ({suf[1:]}-opt)"
                break
        rest = base[len("autodock_mgltools"):].lstrip("_")
        exh = rest[3:] if rest.startswith("exh") else "32"
        return f"AutoDock Vina exh{exh}{opt}"
    if m == "autodock":
        return "AutoDock Vina"
    if m == "autodock_smina":
        return "AutoDock Vina (smina-opt)"
    if m == "autodock_gnina":
        return "AutoDock Vina (gnina-opt)"
    if m == "autodock_vinardo":
        return "AutoDock Vinardo"
    if m == "autodock_vinardo_smina":
        return "AutoDock Vinardo (smina-opt)"
    if m == "autodock_vinardo_gnina":
        return "AutoDock Vinardo (gnina-opt)"
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
    if m == "crystal":
        return "Crystal (native)"
    if not m.startswith("equibind"):
        return m
    pocket, refine, clamp = _eq_tokens(m)
    parts = [p for p in (pocket,
                         "smina-opt" if refine == "smina" else
                         "gnina-opt" if refine == "gnina" else
                         "raw" if refine else None,
                         "clamp on" if clamp == "clampON" else "clamp off" if clamp else None)
             if p]
    return f"EquiBind ({', '.join(parts)})" if parts else "EquiBind"


def method_sort_key(m: str):
    if m == "crystal":
        return (-1, 0, 0)
    # Bucket 0 = Vina-family engines: AutoDock Vina (+opt), Vinardo (+opt),
    # Uni-Dock, Uni-Dock2 — raw first, then smina/gnina optimizer variants.
    if m.startswith("autodock_vinardo"):
        return (0, 10 + {"autodock_vinardo": 0, "autodock_vinardo_smina": 1,
                         "autodock_vinardo_gnina": 2}.get(m, 3), 0)
    if m.startswith("autodock"):
        return (0, {"autodock": 0, "autodock_smina": 1, "autodock_gnina": 2}.get(m, 3), 0)
    if m == "unidock":
        return (0, 20, 0)
    if m == "unidock2":
        return (0, 21, 0)
    if m.startswith("diffdock"):
        return (1, {"diffdock": 0, "diffdock_smina": 1, "diffdock_gnina": 2}.get(m, 3), 0)
    if m.startswith("equibind"):
        p, r, c = _eq_tokens(m)
        return (2, _POCKET_ORDER.get(p, 9), {None: 0, "raw": 1, "smina": 2}.get(r, 0))
    return (3, 0, 0)


def ordered_methods(methods) -> list[str]:
    return sorted([m for m in methods if m != "crystal"], key=lambda m: (method_sort_key(m), m))


# ── docking-tool colour palette (one palette for the whole script) ───────────
# Every figure colours a docking tool the same way, matching the PoseBusters
# pose_comparison report: AutoDock Vina = blue, DiffDock = orange (its smina/gnina
# optimiser variants a darker/lighter orange), and every EquiBind variant a shade
# of green. Mirrors posebusters_pose_comparison.TOOL_COLORS; kept local (like
# pretty_method / _eq_tokens above) so this report needs no import of that module.
_BASE_COLORS = {
    "autodock": "#1f77b4",           # AutoDock Vina (blue)
    "autodock_vinardo": "#17becf",   # AutoDock Vinardo (cyan)
    "unidock": "#9467bd",            # Uni-Dock tiled (purple)
    "unidock2": "#8c564b",           # Uni-Dock2 (brown)
    "diffdock": "#ff7f0e",           # DiffDock (orange)
}
_DD_VARIANT_COLORS = {"diffdock_smina": "#d95f02", "diffdock_gnina": "#fdae6b"}
# Blue/cyan family for AutoDock Vina + Vinardo gnina/smina optimizer variants
# (raw uses the base blue/cyan in _BASE_COLORS).
_AD_VARIANT_COLORS = {
    "autodock_smina": "#6baed6", "autodock_gnina": "#08519c",
    "autodock_vinardo_smina": "#9edae5", "autodock_vinardo_gnina": "#0e7c86",
}
# Green family for EquiBind variants (cycled if more than this many appear).
_EQ_PALETTE = ["#2ca02c", "#74c476", "#1b7837", "#a6dba0",
               "#006d2c", "#5aae61", "#00441b", "#c7e9c0"]
# Common EquiBind variants pre-seeded so a variant keeps its shade regardless of
# how many appear / which figure is drawn first.
_EQ_PRESEED = ["equibind_unguided", "equibind_fpocket", "equibind_p2rank",
               "equibind_unguided_gnina", "equibind_fpocket_gnina", "equibind_p2rank_gnina",
               "equibind_unguided_raw", "equibind_fpocket_raw", "equibind_p2rank_raw",
               "equibind_unguided_smina", "equibind_fpocket_smina", "equibind_p2rank_smina"]


class _MethodColorMap:
    """Stable colour per docking-tool method key. AutoDock/DiffDock are fixed;
    EquiBind variants draw from a green palette (encounter order, cached), so a
    given variant reads the same green in every figure of a run."""

    def __init__(self):
        self._cache = {m: _EQ_PALETTE[i % len(_EQ_PALETTE)]
                       for i, m in enumerate(_EQ_PRESEED)}
        self._next = len(_EQ_PRESEED)

    def get(self, key, default="#7f7f7f"):
        key = str(key)
        if key in _BASE_COLORS:
            return _BASE_COLORS[key]
        if key in _AD_VARIANT_COLORS:
            return _AD_VARIANT_COLORS[key]
        if key in _DD_VARIANT_COLORS:
            return _DD_VARIANT_COLORS[key]
        if key == "crystal":
            return "#555555"
        if key.startswith("equibind"):
            if key not in self._cache:
                self._cache[key] = _EQ_PALETTE[self._next % len(_EQ_PALETTE)]
                self._next += 1
            return self._cache[key]
        return default

    def __getitem__(self, key):
        return self.get(key)


TOOL_COLORS = _MethodColorMap()


# Default location of the oracle summary written by posebusters_pose_comparison.py
# (the only place the per-variant success metrics exist), used to rank EquiBind/DiffDock
# variants for --best-*-only (by PB-Valid AND RMSD ≤ 2 Å) when neither --oracle-summary nor a
# config value is given. Mirrors run_pandamap.DEFAULT_ORACLE_SUMMARY.
DEFAULT_ORACLE_SUMMARY = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/oracle_summary.csv")

# Per-pose RMSD-to-crystal table written by posebusters_pose_comparison.py — the only
# place per-pose RMSD lives (PandaMap carries none). Used by fig 14 to cross geometry
# against interaction recovery. Sits beside the oracle summary; join key is
# (method, protein, ligand, pose_name).
DEFAULT_PER_POSE_METRICS = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")


def _resolve_variant_oracle(oracle_csv: Path) -> Path:
    """Prefer the ``*_all_variants.csv`` sibling of an oracle summary for per-variant ranking.

    posebusters_pose_comparison.py writes oracle_summary.csv with each tool COLLAPSED to a
    single row (``diffdock`` / ``equibind_unguided_gnina`` — its own best variant, relabelled
    to the bare tool name). This report splits DiffDock/EquiBind into per-variant labels
    (``diffdock_smina`` / ``diffdock_gnina`` / …). Ranking those split labels against the
    collapsed file matches only the bare ``diffdock`` row, so every other variant is dropped
    (``m in scores.index`` is False for it) and raw DiffDock "wins" by being the sole candidate.

    The comparison also writes an ``oracle_summary_all_variants.csv`` sibling that keeps one row
    per variant; when it exists we rank against that instead, so every variant is compared. Falls
    back to the given path (older runs / crystal-free sets that lack the sibling); a path already
    ending in ``_all_variants.csv`` is returned unchanged.
    """
    if not oracle_csv:
        return oracle_csv
    p = Path(oracle_csv)
    if p.name.endswith("_all_variants.csv"):
        return p
    sibling = p.with_name(f"{p.stem}_all_variants{p.suffix}")
    return sibling if sibling.exists() else p


# Metric the best-variant selectors rank on: the combined docking-success criterion
# PB-Valid AND RMSD ≤ 2 Å (a pose must be BOTH near-native AND physically valid), written
# per variant by posebusters_pose_comparison.py. Older summaries that predate the combined
# column fall back to the RMSD-only rate.
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


def select_best_equibind(summary: pd.DataFrame, inter: pd.DataFrame,
                         oracle_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    """Drop every EquiBind variant except the single best-performing one.

    "Best" = the EquiBind variant with the highest PB-Valid AND RMSD ≤ 2 Å rate
    (``oracle_pb_valid_and_rmsd2_%``, falling back to ``oracle_rmsd_le_2.0A_%`` for older
    summaries — see _variant_ranking_scores) in *oracle_csv* — the oracle_summary.csv written
    by posebusters_pose_comparison.py (the only place that metric exists; PandaMap carries no
    RMSD). The chosen variant is selected among the EquiBind variants actually present here.
    AutoDock/DiffDock/crystal rows are always kept.

    Returns the filtered (summary, inter) frames and the chosen variant's method
    key, or the inputs unchanged with ``None`` when filtering can't be applied
    (no EquiBind variants present, missing/unreadable oracle summary, or no
    overlap between the present variants and the ranking metric).
    """
    eq_present = sorted({m for m in summary["method"].astype(str).unique()
                         if m.startswith("equibind")})
    if not eq_present:
        return summary, inter, None
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-equibind-only] oracle summary not found at {oracle_csv} — "
              "keeping all EquiBind variants.")
        return summary, inter, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-equibind-only")
    if scores is None:
        return summary, inter, None
    cand = {m: float(scores[m]) for m in eq_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-equibind-only] none of the present EquiBind variants have a "
              f"score in {oracle_csv} — keeping all EquiBind variants.")
        return summary, inter, None

    best = max(cand, key=cand.get)

    def _keep(df: pd.DataFrame) -> pd.DataFrame:
        m = df["method"].astype(str)
        return df[(~m.str.startswith("equibind")) | (m == best)].reset_index(drop=True)

    return _keep(summary), _keep(inter), best


def select_equibind_variant(summary: pd.DataFrame, inter: pd.DataFrame,
                            spec: str) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    """Drop every EquiBind variant except the one(s) matching *spec*.

    Explicit, oracle-free counterpart to select_best_equibind (mirrors
    posebusters_validity_report.select_equibind_variant). *spec* tokens (split on
    '_' or '/') are matched against the label's pocket/refine/clamp tokens, e.g.
    'gnina' or 'unguided_gnina'; a token subset constrains only the axes it names.
    AutoDock/DiffDock/crystal rows are always kept. Returns the filtered
    (summary, inter) frames and the kept variant key when *spec* resolves to
    exactly one (else ``None``, with rows still filtered / unchanged when nothing
    matches)."""
    tokens = [t for t in str(spec).strip().lower().replace("/", "_").split("_") if t]
    eq_present = sorted({m for m in summary["method"].astype(str).unique()
                         if m.startswith("equibind")})
    if not tokens or not eq_present:
        return summary, inter, None
    matched = [m for m in eq_present
               if all(t in set(m.lower().split("_")[1:]) for t in tokens)]
    if not matched:
        print(f"  [equibind-variant] no EquiBind variant matches '{spec}' "
              f"(present: {eq_present}) — keeping all EquiBind variants.")
        return summary, inter, None

    def _keep(df: pd.DataFrame) -> pd.DataFrame:
        m = df["method"].astype(str)
        return df[(~m.str.startswith("equibind")) | (m.isin(matched))].reset_index(drop=True)

    return _keep(summary), _keep(inter), (matched[0] if len(matched) == 1 else None)


def select_best_diffdock(summary: pd.DataFrame, inter: pd.DataFrame,
                         oracle_csv: Path,
                         pin: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    """Drop every DiffDock optimizer variant except one.

    When *pin* is given ('raw' | 'smina' | 'gnina') the named variant is kept
    directly (``diffdock`` for 'raw', else ``diffdock_<pin>``), bypassing the
    oracle entirely — use this to force gnina regardless of the benchmark ranking.

    Otherwise "best" = the DiffDock variant (diffdock / diffdock_smina / diffdock_gnina) with
    the highest PB-Valid AND RMSD ≤ 2 Å rate (``oracle_pb_valid_and_rmsd2_%``, falling back to
    ``oracle_rmsd_le_2.0A_%`` for older summaries — see _variant_ranking_scores) in *oracle_csv*.
    Mirrors select_best_equibind; AutoDock/EquiBind/crystal rows are always kept. Returns
    the filtered (summary, inter) frames and the kept variant key, or the inputs
    unchanged with ``None`` when no DiffDock variant is present, a pinned variant is
    absent, the summary is missing/unreadable, or no present variant has a score."""
    dd_present = sorted({m for m in summary["method"].astype(str).unique()
                         if m.startswith("diffdock")})
    if not dd_present:
        return summary, inter, None
    if pin:
        target = "diffdock" if pin == "raw" else f"diffdock_{pin}"
        if target not in dd_present:
            print(f"  [diffdock-variant={pin}] '{target}' not among present DiffDock "
                  f"variants {dd_present} — keeping all DiffDock variants.")
            return summary, inter, None

        def _keep_pin(df: pd.DataFrame) -> pd.DataFrame:
            m = df["method"].astype(str)
            return df[(~m.str.startswith("diffdock")) | (m == target)].reset_index(drop=True)

        return _keep_pin(summary), _keep_pin(inter), target
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-diffdock-only] oracle summary not found at {oracle_csv} — "
              "keeping all DiffDock variants.")
        return summary, inter, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-diffdock-only")
    if scores is None:
        return summary, inter, None
    cand = {m: float(scores[m]) for m in dd_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-diffdock-only] none of the present DiffDock variants have a "
              f"score in {oracle_csv} — keeping all DiffDock variants.")
        return summary, inter, None

    best = max(cand, key=cand.get)

    def _keep(df: pd.DataFrame) -> pd.DataFrame:
        m = df["method"].astype(str)
        return df[(~m.str.startswith("diffdock")) | (m == best)].reset_index(drop=True)

    return _keep(summary), _keep(inter), best


# ── pose-selection cap (re-enforce run_pandamap's top-N on the loaded data) ──

def cap_top_n_per_combo(summary: pd.DataFrame, inter: pd.DataFrame,
                        n: "int | None") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep only each method's top-*n* poses per (protein, ligand), ranked by pose_rank.

    run_pandamap already selects the top-N poses per method × complex before mapping,
    but it writes ``pandamap_pose_summary.csv`` / ``pandamap_interactions.csv`` in
    resume/append mode. A results directory reused across runs therefore accumulates
    MORE than N poses per method × complex (and stale variants/pairs from older runs),
    which would silently inflate every figure past the stated top-N scope. This
    re-applies the cap on the loaded frames so the report honours that scope no matter
    how the CSVs were built. It only limits the per-combo COUNT (by pose_rank); it does
    not re-derive validity or pair scope — regenerate with ``--overwrite`` for that.

    No-ops when *n* is unset, the frames are empty, or nothing exceeds N (e.g. a freshly
    overwritten, already-capped set). Prints a note only when it actually drops rows.
    """
    key = ["method", "protein", "ligand", "pose_name"]
    if not n or summary.empty or "pose_rank" not in summary.columns:
        return summary, inter
    s = summary.copy()
    s["_rank"] = pd.to_numeric(s["pose_rank"], errors="coerce").fillna(999)
    keep_keys = (s.drop_duplicates(key)
                 .sort_values(["_rank", "pose_name"])
                 .groupby(["method", "protein", "ligand"], sort=False)
                 .head(n)[key].drop_duplicates())
    summary2 = summary.merge(keep_keys, on=key, how="inner").drop_duplicates(key).reset_index(drop=True)
    inter2 = (inter.merge(keep_keys, on=key, how="inner")
              if not inter.empty and set(key).issubset(inter.columns) else inter)
    dropped = len(summary) - len(summary2)
    if dropped:
        print(f"  [top-{n} cap] dropped {dropped} accumulated pose row(s) beyond the top {n} "
              f"per method × complex (summary CSV was built in resume/append mode).")
    return summary2, inter2


# ── ligand chemotype descriptors (RDKit, from the crystal ligand) ──────────

def ligand_descriptors(sdf_path: Path) -> dict | None:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, Crippen
    mol = Chem.MolFromMolFile(str(sdf_path), sanitize=True)
    if mol is None:
        mol = Chem.MolFromMolFile(str(sdf_path), sanitize=False)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return None
    halo = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ("F", "Cl", "Br", "I"))
    return {
        "mw": round(Descriptors.MolWt(mol), 1),
        "n_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "n_halogen": halo,
        "formal_charge": Chem.GetFormalCharge(mol),
        "logp": round(Crippen.MolLogP(mol), 2),
        "n_heavy": mol.GetNumHeavyAtoms(),
    }


def build_descriptor_table(pairs, benchmark_dir: Path | None, summary: pd.DataFrame) -> pd.DataFrame:
    """One descriptor row per (protein, ligand) from the crystal SDF (fallback: a pose)."""
    rows = []
    pose_lookup = (summary.dropna(subset=["pose_file"])
                   .groupby(["protein", "ligand"])["pose_file"].first().to_dict())
    for (protein, ligand) in pairs:
        sdf = None
        if benchmark_dir:
            cand = benchmark_dir / protein / f"{protein}_ligand.sdf"
            if cand.exists():
                sdf = cand
        if sdf is None:
            pf = pose_lookup.get((protein, ligand))
            if pf and Path(pf).exists() and str(pf).lower().endswith(".sdf"):
                sdf = Path(pf)
        if sdf is None:
            continue
        d = ligand_descriptors(sdf)
        if d:
            rows.append({"protein": protein, "ligand": ligand, **d})
    return pd.DataFrame(rows)


# ── fingerprints ───────────────────────────────────────────────────────────

def pose_fingerprint(df_pose: pd.DataFrame, mode: str = "typed") -> set:
    """Interaction fingerprint of one pose's residue-level rows.

    mode='typed'   -> {(interaction_type, resname, resnum, chain)}
    mode='contact' -> {(resname, resnum, chain)}  (residue contacted, any type)
    """
    if mode == "contact":
        return set(zip(df_pose["resname"], df_pose["resnum"], df_pose["chain"]))
    return set(zip(df_pose["interaction_type"], df_pose["resname"],
                   df_pose["resnum"], df_pose["chain"]))


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return np.nan
    u = len(a | b)
    return len(a & b) / u if u else np.nan


# ── nearest-copy crystal reference (opt-in; nearest-copy endpoint program) ──
# The canonical crystal_interactions.csv fingerprints only the reference copy (record 0
# of <ID>_ligands.sdf) and the per-pose metrics table carries no copy index, so every
# function below falls back to today's ``crystal.groupby(["protein", "ligand"])`` and the
# report is byte-identical. When run_pandamap ran with ``crystal_copies: all`` (crystal
# rows carry ``copy_index``) AND the per-pose metrics table was built under
# ``--reference-convention nearest`` (rows carry ``nearest_copy_index``), each pose is
# scored against the crystal fingerprint of the deposited copy it is nearest to (plan
# v2 D2), joined on (method, protein, ligand, pose_name). A top-k UNION fingerprint is
# scored against the union of the fingerprints of the copies its poses selected (plan
# v2 §1.8). A pose absent from the metrics table falls back to the reference copy and
# is counted.

class CrystalCopySelector:
    """Per-pose crystal reference under the nearest-copy convention."""

    KEY = ["method", "protein", "ligand", "pose_name"]

    def __init__(self, crystal: pd.DataFrame, nearest: dict, ref_index: dict):
        ci = crystal["copy_index"].fillna(0).astype(int)
        self.copies = {(p, l, k): g for (p, l, k), g in
                       crystal.assign(_k=ci).groupby(["protein", "ligand", "_k"])}
        self.pairs = {(p, l) for (p, l, _k) in self.copies}
        self.nearest = nearest            # (method, protein, ligand, pose_name) -> copy k
        self.ref_index = ref_index        # (protein, ligand) -> reference copy index (0)
        self.fallback: set = set()        # pose keys that had no nearest_copy_index
        self._fp: dict = {}

    @classmethod
    def build(cls, crystal: pd.DataFrame, metrics_csv) -> "CrystalCopySelector | None":
        """None unless BOTH ``copy_index`` (crystal) and ``nearest_copy_index`` (metrics)."""
        if crystal.empty or "copy_index" not in crystal.columns:
            return None
        if not metrics_csv or not Path(metrics_csv).exists():
            return None
        head = pd.read_csv(metrics_csv, nrows=0)
        if "nearest_copy_index" not in head.columns:
            return None
        cols = cls.KEY + ["nearest_copy_index"] + (["ref_copy_index"] if "ref_copy_index" in head.columns else [])
        mt = pd.read_csv(metrics_csv, usecols=cols, low_memory=False).dropna(subset=["nearest_copy_index"])
        mt = mt.drop_duplicates(cls.KEY)
        nearest = {tuple(k): int(v) for *k, v in
                   mt[cls.KEY + ["nearest_copy_index"]].itertuples(index=False, name=None)}
        ref_index = {}
        if "ref_copy_index" in mt.columns:
            r = mt.dropna(subset=["ref_copy_index"]).drop_duplicates(["protein", "ligand"])
            ref_index = {(p, l): int(v) for p, l, v in
                         r[["protein", "ligand", "ref_copy_index"]].itertuples(index=False, name=None)}
        return cls(crystal, nearest, ref_index)

    # ----- copy resolution -------------------------------------------------
    def ref_copy(self, p, l) -> int:
        return self.ref_index.get((p, l), 0)

    def copy_of(self, method, p, l, pose_name) -> int:
        k = self.nearest.get((method, p, l, pose_name))
        if k is None:
            self.fallback.add((method, p, l, pose_name))
            return self.ref_copy(p, l)
        return k

    def frame(self, p, l, k) -> pd.DataFrame:
        """Crystal rows of copy k (empty frame when that copy has no interaction rows)."""
        g = self.copies.get((p, l, k))
        if g is None:
            return pd.DataFrame(columns=["interaction_type", "resname", "resnum", "chain"])
        return g

    def fp(self, p, l, k, mode="typed") -> set:
        key = (p, l, k, mode)
        if key not in self._fp:
            self._fp[key] = pose_fingerprint(self.frame(p, l, k), mode)
        return self._fp[key]

    def has_pair(self, p, l) -> bool:
        return (p, l) in self.pairs

    # ----- per-pose / per-pose-set references -------------------------------
    def pose_fp(self, method, p, l, pose_name, mode="typed") -> set:
        return self.fp(p, l, self.copy_of(method, p, l, pose_name), mode)

    def copies_for(self, method, p, l, pose_names) -> list[int]:
        """Distinct copy indices (sorted) the given poses of one method × pair select."""
        return sorted({self.copy_of(method, p, l, n) for n in pose_names})

    def union_fp(self, method, p, l, pose_names, mode="typed") -> set:
        """Union of the fingerprints of every copy the pose set selects (plan v2 §1.8)."""
        out: set = set()
        for k in self.copies_for(method, p, l, pose_names):
            out |= self.fp(p, l, k, mode)
        return out

    def report(self) -> str:
        n_multi = sum(1 for (p, l) in self.pairs
                      if len({k for (pp, ll, k) in self.copies if (pp, ll) == (p, l)}) > 1)
        alt = sum(1 for (m, p, l, n), k in self.nearest.items() if k != self.ref_copy(p, l))
        return (f"nearest-copy crystal reference: {len(self.pairs)} pairs, {n_multi} with >1 copy "
                f"fingerprinted; {len(self.nearest)} poses carry a nearest_copy_index "
                f"({alt} point at an alternate copy); {len(self.fallback)} profiled poses without "
                f"an index fell back to the reference copy "
                f"(in {len({(p, l) for (_m, p, l, _n) in self.fallback})} pairs)")


# ── shared figure finalisation + per-figure text sidecar ────────────────────
# run_pandamap keeps only each method's top-N ranked poses per complex
# (poses_per_combo / select_top_n), so the whole interaction analysis is built
# from that top-N set. Rather than stamp that scope (or any second title line)
# onto the image, every figure carries a single companion "<stem>.txt" holding
# what used to sit on the figure: its title, any second/subtitle line, the
# pose-selection-scope footer (plus per-figure footnotes) and the statistical
# results. So each graph keeps a one-line title and no footer, with the detail
# one file away. The scope string is set once in main() from poses_per_combo.
_POSE_SCOPE_CAPTION: str = ""


def set_pose_scope_caption(n: "int | None") -> None:
    """Record the top-N pose-selection scope written to every figure's .txt footer."""
    global _POSE_SCOPE_CAPTION
    _POSE_SCOPE_CAPTION = (
        f"Built from each method's top {n} ranked pose(s) per complex "
        f"(fewer where a method has under {n} valid poses)." if n else "")


# Per-figure text buffer, keyed by the figure's output-path string. _title /
# _note_subtitle / _note_footer / _write_stats_txt append here; _flush_fig_txt
# (called from main) writes one "<stem>.txt" per figure. Buffering (rather than
# writing on the spot) keeps a single, consistently ordered file regardless of
# whether a plot finalises the figure before or after recording its stats, and
# keeps re-runs idempotent (the buffer is rebuilt fresh each process).
_FIG_TXT: "dict[str, dict]" = {}


def _fig_entry(out: Path) -> dict:
    return _FIG_TXT.setdefault(
        str(Path(out)), {"title": None, "subtitle": [], "footer": [], "stats": []})


def _note_subtitle(out: Path, line) -> None:
    """Queue a figure's moved-off subtitle line(s) for its .txt sidecar."""
    e = _fig_entry(out)
    for ln in ([line] if isinstance(line, str) else list(line)):
        ln = " ".join(str(ln).split())
        if ln and ln not in e["subtitle"]:
            e["subtitle"].append(ln)


def _note_footer(out: Path, foot) -> None:
    """Queue a figure's footer / footnote line(s) for its .txt sidecar."""
    e = _fig_entry(out)
    for ln in ([foot] if isinstance(foot, str) else list(foot)):
        ln = " ".join(str(ln).split())
        if ln and ln not in e["footer"]:
            e["footer"].append(ln)


def _title(out: Path, text: str) -> str:
    r"""Keep only a figure title's first line on the image; move the rest to the .txt.

    Wrap every figure-level ``set_title`` / ``suptitle`` in this: it records the
    full title and any ``\n``-separated continuation as the sidecar subtitle, and
    returns the single first line to hand to matplotlib. Guarantees no figure
    carries a two-line heading."""
    parts = str(text).split("\n")
    _fig_entry(out)["title"] = parts[0].strip()
    _note_subtitle(out, [p for p in parts[1:] if p.strip()])
    return parts[0].strip()


def _finalize_fig(fig, out: Path, dpi: int = 160, caption: bool = True) -> None:
    """Save + close. The pose-selection-scope footer is recorded to the figure's
    "<stem>.txt" sidecar (never stamped on the image), so every graph stays
    footer-free. ``caption`` is accepted for call-site compatibility and no longer
    changes the image; the scope is always recorded to the sidecar."""
    if _POSE_SCOPE_CAPTION:
        _note_footer(out, _POSE_SCOPE_CAPTION)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


def _flush_fig_txt(out_dir: Path) -> None:
    """Write one "<stem>.txt" per figure: title, subtitle, footer, then stats.

    Everything moved off the figures (second title lines, footers) plus the
    statistical-test results, one human-readable file beside each .png."""
    written = 0
    for path, e in _FIG_TXT.items():
        body: list[str] = []
        if e["title"]:
            body += [f"Figure: {e['title']}"]
        if e["subtitle"]:
            body += ["", "Subtitle:"] + [f"  {s}" for s in e["subtitle"]]
        if e["footer"]:
            body += ["", "Footer:"] + [f"  {f}" for f in e["footer"]]
        for st in e["stats"]:
            body += ["", st["header"], "=" * len(st["header"]), ""] + st["lines"]
        text = "\n".join(body).strip()
        if not text:
            continue
        try:
            Path(path).with_suffix(".txt").write_text(text + "\n")
            written += 1
        except Exception as ex:
            print(f"  [fig-txt] could not write sidecar for {Path(path).name}: {ex}")
    print(f"  [fig-txt] wrote {written} figure text sidecar(s).")


# ── statistical tests (stats_utils) ─────────────────────────────────────────
# Every test is wrapped by its caller in try/except so a stats failure degrades
# to the current test-free figure and never crashes the pipeline. Numeric results
# accumulate in _STATS and are written to interaction_stats.json at the end of the
# run (recoverable sidecar next to the figures). UNIT OF ANALYSIS: correlated poses
# are aggregated to ONE value per complex (protein, ligand) before any paired test;
# the pooled-pose tests (14, and the categorical association G-tests 03/07/08) say
# so in their annotation / JSON 'unit' field.
_STATS: dict = {}
_MIN_UNITS = 3            # fewest complexes for a paired test; below -> 'exploratory'
_EXPLORATORY_MAX = 8      # n below this is flagged small/exploratory in captions


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


def _record(section: str, payload) -> None:
    _STATS[section] = payload


def _write_stats_json(out_dir: Path) -> None:
    if not _STATS:
        return
    try:
        (out_dir / "interaction_stats.json").write_text(
            json.dumps(_json_safe(_STATS), indent=2))
        print(f"  [stats] wrote interaction_stats.json ({len(_STATS)} section(s)).")
    except Exception as e:
        print(f"  [stats] could not write interaction_stats.json: {e}")


def _fmt_p_txt(p) -> str:
    """'p=<value> <stars>' for the plain-text stats sidecars."""
    if p is None or (isinstance(p, float) and p != p):
        return "p=n/a"
    val = su.fmt_p(p) if _HAVE_STATS else f"{p:.3g}"
    star = su.p_stars(p) if _HAVE_STATS else ""
    return f"p={val}" + (f" {star}" if star else "")


def _write_stats_txt(fig_out: Path, header: str, lines) -> None:
    """Record a figure's statistical test results for its sidecar ``<stem>.txt``.

    Buffers into the per-figure text sidecar (flushed by ``_flush_fig_txt``) so the
    stats share one file with the figure's moved-off subtitle and footer. The
    figures themselves carry no test annotations (significance stars, p-values,
    test statistics, effect sizes, CI numbers) — this file holds them instead.
    ``lines`` may nest lists (flattened) and contain ``None`` entries (skipped).
    Nothing is recorded when there is no content. Never raises."""
    flat: list[str] = []

    def _add(x):
        if x is None:
            return
        if isinstance(x, (list, tuple)):
            for y in x:
                _add(y)
        else:
            flat.append(str(x))

    _add(lines)
    if not flat:
        return
    _fig_entry(fig_out)["stats"].append({"header": str(header), "lines": flat})


def _sig_star(s) -> str:
    """Keep only a real significance star ('*','**','***'); drop 'ns'/''/None."""
    return s if s in ("*", "**", "***") else ""


def _ntag(n: int) -> str:
    return "  — exploratory (small n)" if (n or 0) < _EXPLORATORY_MAX else ""


def _paired_count_matrix(summary: pd.DataFrame, order, col: str):
    """complex × tool matrix of the per-complex MEAN of *col* (listwise complete)."""
    if col not in summary.columns:
        return None, []
    wide = (summary.groupby(["protein", "ligand", "method"])[col].mean()
            .unstack("method"))
    methods = [m for m in order if m in wide.columns]
    wide = wide[methods].dropna()
    return wide, methods


def _paired_continuous_entry(wide, methods, unit: str) -> dict:
    """Run su.paired_continuous on a complex×tool matrix, packaged for the sidecar."""
    if wide is None or wide.empty:
        return {"unit": unit, "n_units": 0, "note": "no complete complexes"}
    n, k = wide.shape
    entry = {"unit": unit, "n_units": int(n),
             "methods": [pretty_method(m) for m in methods], "method_keys": list(methods)}
    # Friedman (via su.paired_continuous -> scipy.friedmanchisquare) requires >=3 tool
    # series; with exactly 2 it raises, so bail with a note for k<3 (a 2-tool set only
    # arises when a dataset lacks a tool or is pinned down to two variants).
    if n < _MIN_UNITS or k < 3 or not _HAVE_STATS:
        entry["note"] = f"n={n} too small — exploratory (no test run)" if n < _MIN_UNITS else \
            ("stats_utils unavailable" if not _HAVE_STATS else
             f"only {k} tool(s) — Friedman omnibus needs >=3")
        return entry
    res = su.paired_continuous(wide, labels=list(methods))
    entry["omnibus"] = res["omnibus"]
    entry["medians"] = {pretty_method(m): res["medians"][m] for m in methods}
    entry["pairwise"] = [{**pw, "a": pretty_method(pw["a"]), "b": pretty_method(pw["b"])}
                         for pw in res["pairwise"]]
    return entry


def _interaction_count_stats(summary: pd.DataFrame, order, present_types) -> dict:
    """Paired per-complex interaction-count comparison across tools (figs 01 / 02b).

    Per interaction type (and the total): per-complex MEAN count per tool ->
    su.paired_continuous (Friedman + Wilcoxon pairwise), with su.bh_fdr across the
    interaction-type family. Idempotent — safe to call from both plot functions.
    """
    per_type = {}
    for t in present_types:
        wide, methods = _paired_count_matrix(summary, order, t)
        per_type[t] = _paired_continuous_entry(wide, methods,
                                               "per-complex mean count per tool")
    ptypes = [t for t in present_types if "omnibus" in per_type[t]
              and per_type[t]["omnibus"].get("p") == per_type[t]["omnibus"].get("p")]
    if ptypes and _HAVE_STATS:
        qs = su.bh_fdr([per_type[t]["omnibus"]["p"] for t in ptypes])
        for t, q in zip(ptypes, qs):
            per_type[t]["omnibus"]["p_bh"] = float(q)
            per_type[t]["omnibus"]["star_bh"] = su.p_stars(q)
    twide, tmethods = _paired_count_matrix(summary, order, "total_interactions")
    total = _paired_continuous_entry(twide, tmethods, "per-complex mean total per tool")
    payload = {"family": "16 PandaMap interaction types (BH-FDR)",
               "per_type": per_type, "total": total}
    _record("01_02b_interaction_counts", payload)
    return payload


def _rank_trend_stats(df: pd.DataFrame, order, value_cols, section: str, unit: str) -> dict:
    """Per-complex rank-trend test for the by-rank figures (04c / 10b).

    For each tool and each value column, compute the per-complex Kendall τ between
    pose_rank and the value (one slope per complex — correlated poses are never
    pooled), then a one-sample Wilcoxon signed-rank of those τ's against 0 with a
    matched-pairs rank-biserial effect size. Interpretation: a **negative** median τ
    on a quality metric (F1, precision, recall, matched contacts) means the tool's
    own ranking orders pose quality (rank 1 is best); a **positive** τ on
    'missed'/'spurious' means the same. BH-FDR across the whole (tool × column)
    family. EquiBind's rank axis is gnina-affinity (a re-scoring proxy, not a native
    confidence), so read its trend as the negative-control-ish baseline.

    Records the full result under `section` and returns {method: {col: entry}}.
    """
    methods = [m for m in order if m in set(df["method"])]
    out: dict = {m: {} for m in methods}
    if not _HAVE_STATS:
        _record(section, {"unit": unit, "note": "stats_utils unavailable"})
        return out
    flat_p, flat_key = [], []
    for m in methods:
        sub = df[df["method"] == m]
        for col in value_cols:
            taus = []
            for _, g in sub.groupby(["protein", "ligand"], sort=False):
                if g["pose_rank"].nunique() < 3:       # need ≥3 ranked poses for a τ
                    continue
                tau, _p, _n = su.kendall(g["pose_rank"].to_numpy(float),
                                         g[col].to_numpy(float))
                if tau == tau:
                    taus.append(tau)
            entry = {"n_complexes": len(taus)}
            if len(taus) >= _MIN_UNITS:
                arr = np.asarray(taus, float)
                rb, p, npair = su.wilcoxon_rankbiserial(arr, np.zeros_like(arr))
                entry.update({"median_tau": float(np.median(arr)),
                              "rank_biserial": rb, "p_raw": p, "n_pairs": npair})
                if p == p:
                    flat_p.append(p); flat_key.append((m, col))
            else:
                entry["note"] = "n too small — exploratory"
            out[m][col] = entry
    if flat_p:
        for (m, col), q in zip(flat_key, su.bh_fdr(flat_p)):
            out[m][col]["p_bh"] = float(q)
            out[m][col]["star_bh"] = su.p_stars(q)
    _record(section, {
        "unit": unit,
        "family": "per-complex Kendall tau(rank, value) -> one-sample Wilcoxon vs 0 (BH-FDR)",
        "n_units_note": "one Kendall tau per complex; tool needs >=3 ranked poses in that complex",
        "per_method": {pretty_method(m): out[m] for m in methods}})
    return out


# ── 01: per-method interaction profile ─────────────────────────────────────

def plot_profile(summary: pd.DataFrame, order, out: Path) -> None:
    # Cleveland dot plot: one row per interaction type, one coloured dot per
    # method. Replaces a 16×N grouped-bar wall — far easier to compare methods
    # within a type, and complements the absolute-magnitude heatmap (fig 02).
    present = [t for t in INTERACTION_TYPES if t in summary.columns and summary[t].sum() > 0]
    means = summary.groupby("method")[present].mean().reindex(order).dropna(how="all")
    methods = list(means.index)
    mcolors = {m: TOOL_COLORS.get(m) for m in methods}
    y = np.arange(len(present))[::-1]
    fig, ax = plt.subplots(figsize=(9, max(4.0, 0.46 * len(present) + 1.4)))
    for yi, t in zip(y, present):
        vals = [means.loc[m, t] for m in methods]
        ax.plot([min(vals), max(vals)], [yi, yi], color="#dddddd", lw=2,
                zorder=1, solid_capstyle="round")
    for m in methods:
        ax.scatter(means.loc[m, present].to_numpy(), y, s=55, color=mcolors[m],
                   edgecolor="black", linewidth=0.4, label=pretty_method(m), zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(present)
    ax.set_xlabel("Mean count per pose")
    ax.set_ylabel("Interaction type")
    ax.set_title(_title(out, "Protein-ligand interaction profile by docking method (PandaMap)"))
    ax.legend(fontsize=7, ncol=2)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    # Stats (per-complex mean count per tool -> Friedman across tools per interaction
    # type, BH-FDR across the 16 types) are written to the sidecar .txt, not drawn.
    try:
        st = _interaction_count_stats(summary, order, present)
        pt = st["per_type"]
        n_u = st["total"].get("n_units")
        lines = ["Test: Friedman across tools of the per-complex mean count per pose, one "
                 "test per interaction type; BH-FDR across the interaction-type family.",
                 f"n = {n_u} complexes (listwise-complete).", ""]
        for t in present:
            omn = pt.get(t, {}).get("omnibus", {})
            p = omn.get("p")
            if p is not None and p == p:
                bh = (f"; BH q={su.fmt_p(omn['p_bh'])} {omn.get('star_bh', '')}".rstrip()
                      if omn.get("p_bh") is not None else "")
                lines.append(f"  {pretty_itype(t)}: chi2({omn.get('df')})={omn.get('chi2'):.2f}, "
                             f"{_fmt_p_txt(omn['p'])}, Kendall W={omn.get('kendall_w'):.2f}, "
                             f"n={omn.get('n')}{bh}")
            else:
                lines.append(f"  {pretty_itype(t)}: {pt.get(t, {}).get('note', 'n/a')}")
        _write_stats_txt(out, "01 — Protein-ligand interaction profile by docking method", lines)
    except Exception as e:
        print(f"  [stats] 01 interaction-profile tests skipped: {e}")
    fig.tight_layout(); _finalize_fig(fig, out)


def plot_type_heatmap_box(summary: pd.DataFrame, order, out_heat: Path, out_box: Path) -> None:
    present = [t for t in INTERACTION_TYPES if t in summary.columns and summary[t].sum() > 0]
    means = summary.groupby("method")[present].mean().reindex(order).dropna(how="all")
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(present)), max(5, 0.5 * len(means))))
    sns.heatmap(means, annot=True, fmt=".1f", cmap="viridis",
                yticklabels=[pretty_method(m) for m in means.index],
                cbar_kws={"label": "mean / pose"}, ax=ax)
    ax.set_title(_title(out_heat, "Mean interactions per pose — type × method"))
    fig.tight_layout(); _finalize_fig(fig, out_heat)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(means)), 5))
    data = [summary[summary["method"] == m]["total_interactions"].values for m in means.index]
    bp = ax.boxplot(data, showmeans=True, patch_artist=True)
    for patch, m in zip(bp["boxes"], means.index):   # colour each box by its tool
        patch.set_facecolor(TOOL_COLORS.get(m)); patch.set_alpha(0.85)
    for med in bp["medians"]:
        med.set_color("#333333")
    ax.set_xticks(range(1, len(means) + 1))
    ax.set_xticklabels([pretty_method(m) for m in means.index], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylabel("Total interactions per pose")
    ax.set_title(_title(out_box, "Distribution of total interactions per pose"))
    ax.grid(axis="y", alpha=0.3)
    # Stats (per-complex mean total per tool -> Friedman omnibus + Wilcoxon pairwise)
    # are written to the sidecar .txt, not drawn on the figure.
    try:
        st = _interaction_count_stats(summary, order, present)
        total = st.get("total") or {}
        omn = total.get("omnibus")
        lines = ["Test: Friedman across tools of the per-complex mean total interactions per "
                 "pose (paired), with Wilcoxon signed-rank pairwise (Holm-adjusted).", ""]
        if omn and omn.get("p") == omn.get("p"):
            lines.append(f"Omnibus: chi2({omn['df']})={omn['chi2']:.2f}, {_fmt_p_txt(omn['p'])}, "
                         f"Kendall W={omn['kendall_w']:.2f}, n={omn['n']} complexes.")
        else:
            lines.append(f"Omnibus: {total.get('note', 'n/a')}")
        if total.get("medians"):
            lines += ["", "Per-tool median (per-complex mean total):"]
            lines += [f"  {m}: {v:.2f}" for m, v in total["medians"].items()]
        if total.get("pairwise"):
            lines += ["", "Pairwise (Wilcoxon signed-rank, Holm):"]
            for pw in total["pairwise"]:
                rb = (f", rank-biserial={pw['rank_biserial']:+.2f}"
                      if pw.get("rank_biserial") is not None else "")
                lines.append(f"  {pw['a']} vs {pw['b']}: {_fmt_p_txt(pw.get('p_holm'))} (Holm){rb}")
        _write_stats_txt(out_box, "02b — Distribution of total interactions per pose", lines)
    except Exception as e:
        print(f"  [stats] 02b total-box tests skipped: {e}")
    fig.tight_layout(); _finalize_fig(fig, out_box)


# ── 03: residue hot-spots ───────────────────────────────────────────────────

def _residue_hotspot_stats(inter: pd.DataFrame, order, top_res) -> tuple[dict | None, dict, int]:
    """Association + paired tests for the residue hot-spot map (fig 03).

    * G-test of independence on the (residue × tool) pose-contact-count table
      (Cramér's V + adjusted residuals) — a categorical association test on pooled
      pose contacts.
    * Per-residue McNemar across tools on the PER-COMPLEX contact boolean (did the
      tool contact this residue in any of its poses for that complex?), via
      su.paired_proportions, with su.bh_fdr across the top residues.
    Returns (gtest_dict|None, {res: entry}, n_complexes).
    """
    pc = inter.copy()
    pc["res"] = pc["resname"].astype(str) + pc["resnum"].astype(str)
    methods = [m for m in order if m in set(pc["method"])]
    top_res = [r for r in top_res]

    # G-test on residue × tool pose-contact counts
    gstat = None
    contact = (pc[pc["res"].isin(top_res)]
               .groupby(["method", "res"])["pose_name"].nunique().rename("n").reset_index())
    tab = (contact.pivot(index="res", columns="method", values="n")
           .reindex(index=top_res, columns=methods).fillna(0.0))
    if _HAVE_STATS and tab.shape[0] >= 2 and tab.shape[1] >= 2 and tab.to_numpy().sum() > 0:
        try:
            g = su.gtest_independence(tab.to_numpy(float))
            gstat = {"G": g["G"], "p": g["p"], "df": g["df"], "cramers_v": g["cramers_v"],
                     "n_low_expected": g["n_low_expected"],
                     "rows": top_res, "cols": [pretty_method(m) for m in methods],
                     "residuals": g["residuals"]}
        except Exception as e:
            print(f"  [stats] 03 G-test failed: {e}")

    # Per-residue McNemar across tools on the per-complex contact boolean
    rb = pc.groupby(["method", "protein", "ligand"])["res"].apply(set)
    present: dict[str, set] = {m: set() for m in methods}
    for (m, p, l) in rb.index:
        if m in present:
            present[m].add((p, l))
    common = (sorted(set.intersection(*[present[m] for m in methods]))
              if methods and all(present.values()) else [])
    res_stats: dict[str, dict] = {}
    omni_p, keys = [], []
    for r in top_res:
        entry = {"n_units": len(common)}
        if _HAVE_STATS and len(common) >= _MIN_UNITS and len(methods) >= 2:
            cols = {m: np.array([1 if r in rb.get((m, p, l), set()) else 0
                                 for (p, l) in common]) for m in methods}
            try:
                pp = su.paired_proportions(cols, labels=methods)
                entry["omnibus"] = pp["omnibus"]
                entry["rates"] = {pretty_method(m): pp["rates"][m] for m in methods}
                entry["pairwise"] = [{**pw, "a": pretty_method(pw["a"]),
                                      "b": pretty_method(pw["b"])} for pw in pp["pairwise"]]
                if pp["omnibus"].get("p") == pp["omnibus"].get("p"):
                    omni_p.append(pp["omnibus"]["p"]); keys.append(r)
            except Exception as e:
                entry["error"] = str(e)
        else:
            entry["note"] = "n too small — exploratory"
        res_stats[r] = entry
    if omni_p and _HAVE_STATS:
        for r, q in zip(keys, su.bh_fdr(omni_p)):
            res_stats[r]["omnibus"]["p_bh"] = float(q)
            res_stats[r]["omnibus"]["star_bh"] = su.p_stars(q)
    _record("03_residue_hotspots", {
        "unit": "G-test: pooled pose-contact counts; McNemar: per-complex contact boolean",
        "family": "per-residue McNemar across the top residues (BH-FDR)",
        "n_complexes": len(common), "gtest": gstat, "per_residue": res_stats})
    return gstat, res_stats, len(common)


def plot_residue_hotspots(inter: pd.DataFrame, summary: pd.DataFrame, order,
                          out: Path, top_n: int = 25) -> pd.DataFrame:
    inter = inter.copy()
    inter["res"] = inter["resname"].astype(str) + inter["resnum"].astype(str)
    # poses per method (denominator)
    poses_per_method = summary.groupby("method")["pose_name"].nunique()
    # poses contacting each residue per method
    contact = (inter.groupby(["method", "res"])["pose_name"].nunique()
               .rename("n").reset_index())
    contact["frac"] = contact.apply(lambda r: r["n"] / poses_per_method.get(r["method"], np.nan), axis=1)
    mat = contact.pivot_table(index="res", columns="method", values="frac", fill_value=0.0)
    mat = mat.reindex(columns=[m for m in order if m in mat.columns])
    top = mat.sum(axis=1).sort_values(ascending=False).head(top_n).index
    mat = mat.loc[top]
    # Stats (per-residue McNemar across tools + a residue×tool G-test of independence)
    # are computed here and written to the sidecar .txt, not stamped on the figure.
    gstat = None
    res_stats: dict = {}
    try:
        gstat, res_stats, _n = _residue_hotspot_stats(inter, order, list(mat.index))
    except Exception as e:
        print(f"  [stats] 03 residue-hotspot tests skipped: {e}")
    fig, ax = plt.subplots(figsize=(max(8, 0.9 * mat.shape[1]), max(6, 0.32 * len(mat))))
    sns.heatmap(mat, annot=False, cmap="rocket_r", vmin=0, vmax=1,
                xticklabels=[pretty_method(m) for m in mat.columns],
                yticklabels=[str(r) for r in mat.index],
                cbar_kws={"label": "fraction of poses contacting residue"}, ax=ax)
    ax.set_title(_title(out, f"Residue hot-spots — top {len(mat)} contacted residues × method"),
                 fontsize=12)
    ax.set_ylabel("Residue"); ax.set_xlabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(); _finalize_fig(fig, out)
    try:
        lines = ["Residue hot-spots — fraction of a tool's poses contacting each top residue.", ""]
        if gstat and gstat.get("p") == gstat.get("p"):
            lines += ["Test 1 — G-test of independence on the residue × tool pose-contact "
                      "counts (pooled pose contacts):",
                      f"  G={gstat['G']:.2f}, df={gstat['df']}, {_fmt_p_txt(gstat['p'])}, "
                      f"Cramér's V={gstat['cramers_v']:.2f}, cells with expected<5: "
                      f"{gstat['n_low_expected']}.", ""]
        lines.append("Test 2 — per-residue Cochran's Q / McNemar across tools on the "
                     "per-complex contact boolean (BH-FDR across residues):")
        for r in mat.index:
            e = res_stats.get(r, {})
            omn = e.get("omnibus", {})
            p = omn.get("p")
            if p is not None and p == p:
                bh = (f", BH q={su.fmt_p(omn['p_bh'])} {omn.get('star_bh', '')}".rstrip()
                      if omn.get("p_bh") is not None else "")
                lines.append(f"  {r}: Q={omn.get('Q'):.2f}, df={omn.get('df')}, "
                             f"{_fmt_p_txt(omn['p'])}{bh}  (n={e.get('n_units')} complexes)")
            else:
                lines.append(f"  {r}: {e.get('note', 'n/a')}")
        _write_stats_txt(out, "03 — Residue hot-spots", lines)
    except Exception as e:
        print(f"  [stats] 03 residue-hotspot txt skipped: {e}")
    return mat


# ── 04/05: native recovery + fingerprint similarity ─────────────────────────

def native_recovery(inter: pd.DataFrame, crystal: pd.DataFrame, order,
                    out_dir: Path, copy_sel: "CrystalCopySelector | None" = None) -> pd.DataFrame:
    """Per-pose precision/recall/F1 of interactions vs the crystal fingerprint.

    ``copy_sel`` (nearest-copy convention) scores each pose against the crystal copy it
    is nearest to; ``None`` keeps the per-pair reference fingerprint.
    """
    cryst_fp = {}
    if copy_sel is None:
        for (p, l), g in crystal.groupby(["protein", "ligand"]):
            cryst_fp[(p, l)] = pose_fingerprint(g, "typed")

    recs = []
    for (method, p, l, pose), g in inter.groupby(["method", "protein", "ligand", "pose_name"]):
        if copy_sel is not None:
            ref = copy_sel.pose_fp(method, p, l, pose) if copy_sel.has_pair(p, l) else None
        else:
            ref = cryst_fp.get((p, l))
        if not ref:
            continue
        fp = pose_fingerprint(g, "typed")
        inter_n = len(fp & ref)
        prec = inter_n / len(fp) if fp else np.nan
        rec = inter_n / len(ref) if ref else np.nan
        f1 = (2 * prec * rec / (prec + rec)) if (prec and rec) else 0.0
        recs.append({"method": method, "protein": p, "ligand": l, "pose_name": pose,
                     "precision": prec, "recall": rec, "f1": f1,
                     "jaccard_vs_crystal": jaccard(fp, ref)})
    rdf = pd.DataFrame(recs)
    if rdf.empty:
        return rdf
    rdf.to_csv(out_dir / "native_recovery_per_pose.csv", index=False)

    summ = rdf.groupby("method")[["precision", "recall", "f1", "jaccard_vs_crystal"]].mean()
    summ = summ.reindex([m for m in order if m in summ.index])
    summ.round(3).to_csv(out_dir / "native_recovery_summary.csv")

    # Stats: cluster-bootstrap 95% CI (clusters = complex) on each mean bar, plus a
    # paired per-complex F1 comparison across tools (Friedman + Wilcoxon pairwise).
    ci: dict = {}
    f1_stats = None
    try:
        for m in summ.index:
            sub = rdf[rdf["method"] == m]
            clusters = (sub["protein"].astype(str) + "|" + sub["ligand"].astype(str)).to_numpy()
            for col in ["precision", "recall", "f1"]:
                vals = sub[col].to_numpy(float)
                if _HAVE_STATS and np.isfinite(vals).sum() >= _MIN_UNITS \
                        and len(np.unique(clusters)) >= _MIN_UNITS:
                    ci[(m, col)] = su.cluster_bootstrap_ci(vals, clusters)
        f1w = rdf.groupby(["protein", "ligand", "method"])["f1"].mean().unstack("method")
        f1_methods = [m for m in order if m in f1w.columns]
        f1w = f1w[f1_methods].dropna()
        f1_stats = _paired_continuous_entry(f1w, f1_methods,
                                            "per-complex mean native-F1 per tool")
        _record("04_native_recovery", {
            "paired_f1": f1_stats,
            "bootstrap_ci": {f"{m}|{c}": list(v) for (m, c), v in ci.items()}})
    except Exception as e:
        print(f"  [stats] 04 native-recovery tests skipped: {e}")

    # bar chart of mean precision/recall/F1 (cluster-bootstrap 95% CIs + the paired
    # F1 Friedman go to the sidecar .txt, not onto the figure)
    fig, ax = plt.subplots(figsize=(max(9, 1.2 * len(summ)), 5.5))
    x = np.arange(len(summ)); w = 0.25
    for i, col in enumerate(["precision", "recall", "f1"]):
        ax.bar(x + (i - 1) * w, summ[col].to_numpy(float), w, label=col)
    ax.set_xticks(x)
    ax.set_xticklabels([pretty_method(m) for m in summ.index], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylabel("score"); ax.set_ylim(0, 1)
    ax.set_title(_title(out_dir / "04_native_recovery.png",
                        "Native-interaction recovery vs crystal (mean over poses)\n"
                        "fingerprint = (interaction_type, residue)"), fontsize=12)
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "04_native_recovery.png")
    try:
        lines = ["Native-interaction recovery vs crystal — fingerprint = (interaction type, "
                 "residue).", "",
                 "Per-tool mean precision / recall / F1 with cluster-bootstrap 95% CI "
                 "(clusters = complex):"]
        for m in summ.index:
            parts = []
            for col in ["precision", "recall", "f1"]:
                c = ci.get((m, col))
                if c is not None and c[1] == c[1]:
                    parts.append(f"{col} {summ.loc[m, col]:.3f} [{c[1]:.3f}, {c[2]:.3f}]")
                else:
                    parts.append(f"{col} {summ.loc[m, col]:.3f}")
            lines.append(f"  {pretty_method(m)}: " + "; ".join(parts))
        omn = (f1_stats or {}).get("omnibus")
        if omn and omn.get("p") == omn.get("p"):
            lines += ["", "Paired Friedman across tools (per-complex mean native-F1): "
                          f"chi2({omn['df']})={omn['chi2']:.2f}, {_fmt_p_txt(omn['p'])}, "
                          f"Kendall W={omn['kendall_w']:.2f}, n={omn['n']} complexes."]
            for pw in (f1_stats or {}).get("pairwise", []):
                rb = (f", rank-biserial={pw['rank_biserial']:+.2f}"
                      if pw.get("rank_biserial") is not None else "")
                lines.append(f"  {pw['a']} vs {pw['b']}: {_fmt_p_txt(pw.get('p_holm'))} (Holm){rb}")
        _write_stats_txt(out_dir / "04_native_recovery.png",
                         "04 — Native-interaction recovery", lines)
    except Exception as e:
        print(f"  [stats] 04 native-recovery txt skipped: {e}")

    # Stats for 04b: the 'total-miss' rate = fraction of poses recovering ZERO native
    # typed contacts (F1≈0 — literally the CDF's y-intercept), per method with a Wilson
    # 95% CI. That rate is pose-level (matches the plotted curve) and so pseudoreplicated,
    # so also carry a per-complex miss-fraction paired Friedman across tools (the
    # inference-grade comparison) + a cluster-robust (complex-block) bootstrap CI of that
    # SAME pose-level rate (resamples whole complexes → a cluster-honest interval, not a
    # per-complex-fraction mean; its point estimate equals the pose-level rate).
    miss: dict = {}
    miss_paired: dict = {}
    try:
        for m in summ.index:
            v = rdf[rdf["method"] == m]["f1"].to_numpy(float)
            v = v[np.isfinite(v)]
            k = int((v <= 1e-9).sum()); n = int(v.size)
            lo, hi = su.wilson_ci(k, n) if (_HAVE_STATS and n) else (np.nan, np.nan)
            miss[m] = {"pose_miss_rate": (k / n if n else np.nan), "k": k,
                       "n_poses": n, "wilson_ci": [lo, hi]}
            sub = rdf[rdf["method"] == m]
            clusters = (sub["protein"].astype(str) + "|" + sub["ligand"].astype(str)).to_numpy()
            zvals = (sub["f1"].to_numpy(float) <= 1e-9).astype(float)
            if _HAVE_STATS and len(np.unique(clusters)) >= _MIN_UNITS:
                miss[m]["miss_rate_cluster_robust_ci"] = list(su.cluster_bootstrap_ci(zvals, clusters))
        mc = (rdf.assign(_miss=(rdf["f1"] <= 1e-9).astype(float))
              .groupby(["protein", "ligand", "method"])["_miss"].mean().unstack("method"))
        mm = [m for m in order if m in mc.columns]
        miss_paired = _paired_continuous_entry(mc[mm].dropna(), mm,
                                               "per-complex fraction of poses with F1≈0")
        _record("04b_native_recovery_cdf", {
            "unit": "pose-level total-miss rate (F1≈0) + Wilson CI; per-complex miss-fraction paired",
            "per_method_pose_level": {pretty_method(m): miss[m] for m in summ.index},
            "per_complex_paired": miss_paired,
            "notes": "F1≈0 = a pose recovering zero native typed contacts (the CDF y-intercept). "
                     "The pose-level Wilson CI treats poses as independent (pseudoreplication); "
                     "the per-complex paired Friedman is the inference-grade tool comparison."})
    except Exception as e:
        print(f"  [stats] 04b CDF total-miss tests skipped: {e}")

    # F1 CDF per method (pose-level F1≈0 rate + Wilson CI and the paired Friedman go
    # to the sidecar .txt, not into the legend / title)
    fig, ax = plt.subplots(figsize=(9, 5))
    for m in summ.index:
        vals = np.sort(rdf[rdf["method"] == m]["f1"].dropna().values)
        if not len(vals):
            continue
        ax.plot(vals, np.arange(1, len(vals) + 1) / len(vals) * 100, lw=2,
                color=TOOL_COLORS.get(m), label=pretty_method(m))
    ax.set_xlabel("F1 of native interaction recovery"); ax.set_ylabel("cumulative % of poses")
    ax.set_title(_title(out_dir / "04b_native_recovery_cdf.png",
                        "Native-interaction recovery — F1 CDF"), fontsize=12)
    ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "04b_native_recovery_cdf.png")
    try:
        lines = ["Total-miss rate = fraction of poses recovering ZERO native typed contacts "
                 "(F1≈0, the CDF's y-intercept).", "",
                 "Per-tool pose-level F1≈0 rate with Wilson 95% CI "
                 "(poses treated as independent — pseudoreplicated screening bound):"]
        for m in summ.index:
            e = miss.get(m)
            if e and e.get("n_poses"):
                lo, hi = e["wilson_ci"]
                extra = ""
                cr = e.get("miss_rate_cluster_robust_ci")
                if cr:
                    extra = f"; cluster-robust CI [{cr[1] * 100:.0f}–{cr[2] * 100:.0f}%]"
                lines.append(f"  {pretty_method(m)}: {e['pose_miss_rate'] * 100:.1f}% "
                             f"({e['k']}/{e['n_poses']}) [Wilson {lo * 100:.0f}–{hi * 100:.0f}%]{extra}")
        omn = (miss_paired or {}).get("omnibus")
        if omn and omn.get("p") == omn.get("p"):
            lines += ["", "Per-complex miss-fraction paired Friedman across tools "
                          "(inference-grade comparison): "
                          f"chi2({omn['df']})={omn['chi2']:.2f}, {_fmt_p_txt(omn['p'])}, "
                          f"Kendall W={omn['kendall_w']:.2f}, n={omn['n']} complexes."]
        _write_stats_txt(out_dir / "04b_native_recovery_cdf.png",
                         "04b — Native-interaction recovery F1 CDF", lines)
    except Exception as e:
        print(f"  [stats] 04b CDF txt skipped: {e}")
    return summ


def plot_native_recovery_by_rank(per_pose: pd.DataFrame, order, out: Path,
                                 top_n: int = 5) -> None:
    """Evolution of native-recovery precision / recall / F1 with pose rank, per method.

    Rank-resolved companion to fig 04: three panels (precision, recall, F1), x = pose
    rank k (1 = the tool's top-ranked pose), one line per toolchain. run_pandamap ranks
    AutoDock and DiffDock by native confidence rank and EquiBind's gnina-optimised poses
    by gnina affinity, so all three toolchains evolve across ranks. Any tool that still
    lacks a real rank (e.g. a run without the gnina ranking) is drawn as a flat dashed
    baseline at its pose-set mean instead. Poses are already PB-valid only (pb_valid_only
    in the fingerprint run). Writes native_recovery_by_rank.csv alongside the figure.
    """
    if per_pose is None or per_pose.empty or "pose_rank" not in per_pose.columns:
        return
    df = per_pose.copy()
    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="coerce")
    present = set(df["method"])
    methods = [m for m in order if m in present]
    if not methods:
        return
    # Colour each toolchain by its BASE tool so the palette matches
    # 09f_pbvalid_yield_boxplot (AutoDock blue, DiffDock orange, EquiBind green) —
    # the variant suffix (smina / gnina / pocket) is ignored here, mirroring 09f's
    # tool_color_key {autodock, diffdock, equibind->equibind_unguided}.
    def _base_tool_color(m: str) -> str:
        if m.startswith("diffdock"):
            return TOOL_COLORS.get("diffdock")
        if m.startswith("equibind"):
            return TOOL_COLORS.get("equibind_unguided")
        return TOOL_COLORS.get(m)
    colors = {m: _base_tool_color(m) for m in order}
    scores = [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1")]

    rows, ranked, unranked = [], {}, {}
    for m in methods:
        sub = df[df["method"] == m]
        real = sub[(sub["pose_rank"] >= 1) & (sub["pose_rank"] <= top_n)]
        if real["pose_rank"].nunique() >= 2:          # a genuine rank axis
            per = real.groupby("pose_rank")[["precision", "recall", "f1"]].mean()
            cnt = real.groupby("pose_rank").size()
            ranked[m] = per
            for k, r in per.iterrows():
                rows.append({"method": m, "rank": int(k), "n": int(cnt[k]),
                             "precision": round(r["precision"], 4),
                             "recall": round(r["recall"], 4), "f1": round(r["f1"], 4)})
        else:                                          # unranked (pose_rank 999)
            mean = sub[["precision", "recall", "f1"]].mean()
            unranked[m] = mean
            rows.append({"method": m, "rank": "unranked", "n": int(len(sub)),
                         "precision": round(mean["precision"], 4),
                         "recall": round(mean["recall"], 4), "f1": round(mean["f1"], 4)})
    pd.DataFrame(rows).to_csv(out.parent / "native_recovery_by_rank.csv", index=False)

    # Stats: per-complex Kendall τ(pose rank, score) → one-sample Wilcoxon vs 0, per
    # tool, BH-FDR across the tool×metric family. A negative τ means the tool's own
    # ranking puts better poses first. Only the genuinely-ranked tools are tested.
    trend: dict = {}
    try:
        trend_df = df[(df["pose_rank"] >= 1) & (df["pose_rank"] <= top_n)]
        trend = _rank_trend_stats(trend_df, [m for m in methods if m in ranked],
                                  ["precision", "recall", "f1"],
                                  "04c_native_recovery_by_rank",
                                  "per-complex Kendall tau(pose rank, score); one per complex")
    except Exception as e:
        print(f"  [stats] 04c rank-trend tests skipped: {e}")

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 12), sharex=True, sharey=True)
    for ax, (col, title) in zip(axes, scores):
        for m in methods:
            c = colors.get(m)
            if m in ranked:
                per = ranked[m]
                ax.plot(per.index, per[col], marker="o", lw=2, color=c,
                        label=pretty_method(m))
            else:
                ax.axhline(unranked[m][col], ls="--", lw=1.6, color=c, alpha=0.85,
                           label=f"{pretty_method(m)} (unranked)")
        ax.set_xticks(range(1, top_n + 1))
        ax.set_ylabel("Mean score vs crystal\n(typed interaction recovery)")
        ax.set_title(title)
        ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    axes[-1].set_xlabel("Pose rank k (1 = tool's top-ranked pose)")
    _label_panels(axes)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=11, ncol=len(labels),
               loc="upper center", bbox_to_anchor=(0.5, 0.945), frameon=False)
    fig.suptitle(_title(out,
                        "Native-interaction recovery vs crystal by pose rank\n"
                        "fingerprint = (interaction_type, residue), PB-valid poses only"),
                 fontsize=13, y=0.98)
    _note_footer(out,
                 "Rank basis: AutoDock / DiffDock native confidence rank, EquiBind "
                 "gnina-affinity rank (a tool lacking a rank is shown as a flat baseline)")
    fig.tight_layout(rect=(0, 0.02, 1, 0.925))
    _finalize_fig(fig, out, caption=False)
    try:
        tlines = ["Rank trend: per-complex Kendall tau(pose rank, score) -> one-sample "
                  "Wilcoxon signed-rank vs 0, per tool, BH-FDR across the tool x metric family.",
                  "tau<0 => the tool's own ranking puts better poses first "
                  "(one Kendall tau per complex; a tool needs >=3 ranked poses in that complex).",
                  ""]
        for m in [mm for mm in methods if mm in ranked]:
            tlines.append(f"{pretty_method(m)}:")
            for col, name in (("precision", "Precision"), ("recall", "Recall"), ("f1", "F1")):
                e = (trend.get(m) or {}).get(col, {})
                if "median_tau" in e:
                    bh = (f", BH q={su.fmt_p(e['p_bh'])} {e.get('star_bh', '')}".rstrip()
                          if e.get("p_bh") is not None else "")
                    tlines.append(f"  {name}: median tau={e['median_tau']:+.2f}, "
                                  f"rank-biserial={e.get('rank_biserial'):+.2f}, "
                                  f"{_fmt_p_txt(e.get('p_raw'))}{bh}  "
                                  f"(n={e.get('n_complexes')} complexes)")
                else:
                    tlines.append(f"  {name}: {e.get('note', 'n/a')}")
        _write_stats_txt(out, "04c — Native-interaction recovery by pose rank", tlines)
    except Exception as e:
        print(f"  [stats] 04c rank-trend txt skipped: {e}")


def plot_fingerprint_similarity(inter: pd.DataFrame, crystal: pd.DataFrame, order, out: Path,
                                depth: int = 1, show_caption: bool = True,
                                copy_sel: "CrystalCopySelector | None" = None) -> None:
    """Mean cross-method Jaccard of the top-`depth` pose(s) per pair, + vs crystal.

    ``depth=1`` uses each method's single best (lowest pose_rank) pose. ``depth>1`` uses
    the UNION of that method's top-`depth` poses' typed fingerprints per complex — the
    interaction repertoire the tool covers across its top ranked poses — so the matrix
    compares interaction coverage rather than a single pose.

    ``copy_sel`` (nearest-copy convention): the crystal fingerprint a method's top-`depth`
    union is compared with is the UNION of the fingerprints of the copies those poses are
    nearest to (plan v2 §1.8); crystal-vs-crystal stays the reference copy. ``None`` keeps
    the per-pair reference fingerprint.
    """
    if copy_sel is None:
        cryst_fp = {(p, l): pose_fingerprint(g, "typed") for (p, l), g in crystal.groupby(["protein", "ligand"])}
    else:
        cryst_fp = {(p, l): copy_sel.fp(p, l, copy_sel.ref_copy(p, l)) for (p, l) in copy_sel.pairs}
    # keep the top-`depth` poses per (method, pair) by pose_rank; the union of their
    # interaction rows (pose_fingerprint builds a set over all rows) is the depth-`depth`
    # coverage fingerprint. depth=1 → the single best pose (original behaviour).
    keep = (inter[["method", "protein", "ligand", "pose_name", "pose_rank"]]
            .drop_duplicates()
            .sort_values("pose_rank")
            .groupby(["method", "protein", "ligand"], sort=False)
            .head(depth))
    sel = inter.merge(keep[["method", "protein", "ligand", "pose_name"]],
                      on=["method", "protein", "ligand", "pose_name"], how="inner")
    fp_by = {}
    cryst_for = {}      # nearest-copy: (method, p, l) -> union fp of the copies its poses selected
    for (method, p, l), g in sel.groupby(["method", "protein", "ligand"], sort=False):
        fp_by[(method, p, l)] = pose_fingerprint(g, "typed")   # union over the top-depth poses
        if copy_sel is not None and copy_sel.has_pair(p, l):
            cryst_for[(method, p, l)] = copy_sel.union_fp(method, p, l, g["pose_name"].unique())
    methods = ordered_methods({m for (m, _, _) in fp_by})
    labels = methods + (["crystal"] if cryst_fp else [])
    mat = np.full((len(labels), len(labels)), np.nan)
    pairs = {(p, l) for (_, p, l) in fp_by}

    def _cryst(other, p, l):
        # the crystal fingerprint paired with method `other` at (p, l)
        if copy_sel is not None and other != "crystal":
            return cryst_for.get((other, p, l))
        return cryst_fp.get((p, l))

    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            vals = []
            for (p, l) in pairs:
                fa = _cryst(b, p, l) if a == "crystal" else fp_by.get((a, p, l))
                fb = _cryst(a, p, l) if b == "crystal" else fp_by.get((b, p, l))
                if fa is not None and fb is not None:
                    vals.append(jaccard(fa, fb))
            mat[i, j] = np.nanmean(vals) if vals else np.nan
    fig, ax = plt.subplots(figsize=(max(7, 0.7 * len(labels)), max(6, 0.6 * len(labels))))
    # symmetric matrix — show the lower triangle + diagonal only (drop redundant upper half)
    mask = np.triu(np.ones((len(labels), len(labels)), dtype=bool), k=1)
    sns.heatmap(mat, annot=True, fmt=".2f", cmap="mako", mask=mask,
                xticklabels=[pretty_method(m) for m in labels],
                yticklabels=[pretty_method(m) for m in labels],
                vmin=0, vmax=1, cbar_kws={"label": "mean Jaccard"}, ax=ax)
    pose_desc = "best pose per pair" if depth == 1 else f"union of top-{depth} poses per pair"
    ax.set_title(_title(out, f"Interaction-fingerprint similarity ({pose_desc})"))
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(); _finalize_fig(fig, out, caption=show_caption)


# ── 06/07/08: chemical-level analyses ───────────────────────────────────────

def plot_chem_chemotype(summary: pd.DataFrame, desc: pd.DataFrame, out_dir: Path) -> None:
    if desc.empty:
        print("  [skip] chemotype plots — no ligand descriptors")
        return
    df = summary.merge(desc, on=["protein", "ligand"], how="inner")
    if df.empty:
        print("  [skip] chemotype plots — no descriptor join")
        return
    panels = [
        ("n_halogen", "halogen_bonds", "# halogens", "halogen bonds"),
        ("n_aromatic_rings", "pi_pi_stacking", "# aromatic rings", "π–π stacking"),
        ("n_aromatic_rings", "cation_pi", "# aromatic rings", "cation–π"),
        ("formal_charge", "salt_bridge", "ligand formal charge", "salt bridges"),
        ("formal_charge", "ionic", "ligand formal charge", "ionic"),
        ("hba_hbd", "hydrogen_bonds", "HBD + HBA", "hydrogen bonds"),
    ]
    df["hba_hbd"] = df.get("hbd", 0) + df.get("hba", 0)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, (xcol, ycol, xlab, ylab) in zip(axes.flat, panels):
        if xcol not in df.columns or ycol not in df.columns:
            ax.set_visible(False); continue
        sub = df[[xcol, ycol]].dropna()
        if sub.empty:
            ax.set_visible(False); continue
        grp = sub.groupby(xcol)[ycol]
        ax.scatter(sub[xcol] + np.random.uniform(-0.12, 0.12, len(sub)), sub[ycol],
                   s=12, alpha=0.3, color="#4c72b0")
        m = grp.mean()
        ax.plot(m.index, m.values, "o-", color="#c44e52", lw=2, label="mean")
        ax.set_xlabel(xlab); ax.set_ylabel(f"{ylab} / pose")
        ax.set_title(f"{ylab} vs {xlab}"); ax.grid(alpha=0.3); ax.legend(fontsize=8)
    _label_panels(axes)
    fig.suptitle(_title(out_dir / "06_interaction_vs_chemotype.png",
                        "Chemical level — interaction type vs ligand chemotype"), fontsize=14)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "06_interaction_vs_chemotype.png")
    # correlation table (Pearson CSV kept for backward compatibility)
    cols = ["mw", "n_aromatic_rings", "hbd", "hba", "n_halogen", "formal_charge", "logp"]
    cols = [c for c in cols if c in df.columns]
    itypes = [t for t in INTERACTION_TYPES if t in df.columns and df[t].sum() > 0]
    corr = df[cols + itypes].corr().loc[cols, itypes]
    corr.round(3).to_csv(out_dir / "chem_descriptor_interaction_corr.csv")
    # Heatmap: switch to Spearman (interaction counts are non-normal) computed at the
    # per-COMPLEX unit (descriptor is per complex; pose counts averaged per complex), so
    # n = complexes rather than pooled poses. BH-FDR across the whole grid; significant
    # cells (q<0.05) are starred.
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(itypes)), max(4, 0.6 * len(cols))))
    try:
        if not _HAVE_STATS or not cols or not itypes:
            raise RuntimeError("stats_utils unavailable or empty grid")
        per_complex = df.groupby(["protein", "ligand"])[cols + itypes].mean()
        n_cx = len(per_complex)
        rho = pd.DataFrame(index=cols, columns=itypes, dtype=float)
        pmat = np.full((len(cols), len(itypes)), np.nan)
        for i, c in enumerate(cols):
            for j, t in enumerate(itypes):
                r, p, _ = su.spearman(per_complex[c].to_numpy(float),
                                      per_complex[t].to_numpy(float))
                rho.iloc[i, j] = r
                pmat[i, j] = p
        q = su.bh_fdr(pmat.ravel()).reshape(pmat.shape)
        annot = np.empty(rho.shape, dtype=object)
        for i in range(len(cols)):
            for j in range(len(itypes)):
                r = rho.iloc[i, j]
                annot[i, j] = "" if r != r else f"{r:.2f}"
        sns.heatmap(rho.astype(float), annot=annot, fmt="", cmap="coolwarm", center=0,
                    vmin=-1, vmax=1, annot_kws={"fontsize": 8}, ax=ax)
        ax.set_title(_title(out_dir / "06b_descriptor_corr.png",
                            "Spearman correlation (per complex): ligand descriptor × interaction type\n"
                            f"(n={n_cx} complexes)"), fontsize=10)
        _record("06b_descriptor_corr", {
            "unit": "per-complex mean interaction count vs per-complex descriptor",
            "family": "descriptor × interaction-type grid (BH-FDR)", "n_complexes": n_cx,
            "descriptors": cols, "interaction_types": itypes,
            "rho": rho.astype(float).to_dict(), "p": pmat, "q_bh": q})
        sig = []
        for i, c in enumerate(cols):
            for j, t in enumerate(itypes):
                qq = q[i, j]
                if qq == qq and qq < 0.05:
                    sig.append((abs(float(rho.iloc[i, j])), c, t,
                                float(rho.iloc[i, j]), pmat[i, j], qq))
        sig.sort(reverse=True)
        txt_lines = [f"Test: per-complex Spearman correlation, ligand descriptor × interaction "
                     f"type (n={n_cx} complexes).",
                     f"BH-FDR across the {len(cols)}×{len(itypes)} grid; significant cells "
                     "(q<0.05), strongest first:", ""]
        if sig:
            for _, c, t, r, p, qq in sig:
                txt_lines.append(f"  {c} × {pretty_itype(t)}: rho={r:+.2f}, {_fmt_p_txt(p)}, "
                                 f"BH q={su.fmt_p(qq)}")
        else:
            txt_lines.append("  (no cell significant at FDR<0.05)")
        _write_stats_txt(out_dir / "06b_descriptor_corr.png",
                         "06b — Ligand descriptor × interaction-type correlation", txt_lines)
    except Exception as e:
        print(f"  [stats] 06b Spearman grid skipped ({e}) — Pearson heatmap.")
        ax.clear()
        sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", center=0,
                    vmin=-1, vmax=1, ax=ax)
        ax.set_title(_title(out_dir / "06b_descriptor_corr.png",
                            "Correlation: ligand descriptor × interaction type"))
    fig.tight_layout(); _finalize_fig(fig, out_dir / "06b_descriptor_corr.png")


def plot_residue_class(inter: pd.DataFrame, out: Path) -> None:
    df = inter.copy()
    df["res_class"] = df["resname"].map(residue_class)
    present = [t for t in INTERACTION_TYPES if t in set(df["interaction_type"])]
    tab = (df.groupby(["interaction_type", "res_class"]).size()
           .rename("n").reset_index()
           .pivot(index="interaction_type", columns="res_class", values="n").fillna(0))
    tab = tab.reindex([t for t in present if t in tab.index])
    frac = tab.div(tab.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(frac))))
    frac.plot(kind="barh", stacked=True, ax=ax, colormap="Set2")
    ax.set_xlabel("fraction of interactions"); ax.set_ylabel("")
    title = "Residue-class preference per interaction type"
    # Stats (G-test of independence, interaction type × residue class — a categorical
    # association test on pose-pooled counts) go to the sidecar .txt, not the title.
    try:
        if _HAVE_STATS and tab.shape[0] >= 2 and tab.shape[1] >= 2 and tab.to_numpy().sum() > 0:
            g = su.gtest_independence(tab.to_numpy(float))
            _record("07_residue_class", {
                "unit": "pooled interaction rows (categorical association test)",
                "G": g["G"], "p": g["p"], "df": g["df"], "cramers_v": g["cramers_v"],
                "n_low_expected": g["n_low_expected"],
                "rows": list(tab.index), "cols": list(tab.columns), "residuals": g["residuals"]})
            _write_stats_txt(out, "07 — Residue-class preference per interaction type", [
                "Test: G-test of independence, interaction type × residue class "
                "(counts pooled over poses — a categorical association test).", "",
                f"G={g['G']:.2f}, df={g['df']}, {_fmt_p_txt(g['p'])}, "
                f"Cramér's V={g['cramers_v']:.2f}, cells with expected<5: {g['n_low_expected']}."])
    except Exception as e:
        print(f"  [stats] 07 residue-class G-test skipped: {e}")
    ax.set_title(_title(out, title), fontsize=10 if "\n" in title else 12)
    ax.legend(title="residue class", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout(); _finalize_fig(fig, out)
    tab.to_csv(out.with_suffix(".csv"))


def plot_ligand_element(inter: pd.DataFrame, out: Path) -> None:
    df = inter.dropna(subset=["lig_atom_element"]).copy()
    if df.empty:
        return
    df["elem"] = df["lig_atom_element"].astype(str).str.upper().replace(
        {"F": "halogen", "CL": "halogen", "BR": "halogen", "I": "halogen"})
    present = [t for t in INTERACTION_TYPES if t in set(df["interaction_type"])]
    tab = (df.groupby(["interaction_type", "elem"]).size().rename("n").reset_index()
           .pivot(index="interaction_type", columns="elem", values="n").fillna(0))
    tab = tab.reindex([t for t in present if t in tab.index])
    fig, ax = plt.subplots(figsize=(10, max(5, 0.5 * len(tab))))
    tab.div(tab.sum(axis=1), axis=0).plot(kind="barh", stacked=True, ax=ax, colormap="tab10")
    ax.set_xlabel("fraction of interactions"); ax.set_ylabel("")
    title = "Ligand-atom element per interaction type"
    # Stats (G-test of independence, interaction type × ligand-atom element — a
    # categorical association test on pose-pooled counts) go to the sidecar .txt.
    try:
        if _HAVE_STATS and tab.shape[0] >= 2 and tab.shape[1] >= 2 and tab.to_numpy().sum() > 0:
            g = su.gtest_independence(tab.to_numpy(float))
            _record("08_ligand_element", {
                "unit": "pooled interaction rows (categorical association test)",
                "G": g["G"], "p": g["p"], "df": g["df"], "cramers_v": g["cramers_v"],
                "n_low_expected": g["n_low_expected"],
                "rows": list(tab.index), "cols": list(tab.columns), "residuals": g["residuals"]})
            _write_stats_txt(out, "08 — Ligand-atom element per interaction type", [
                "Test: G-test of independence, interaction type × ligand-atom element "
                "(counts pooled over poses — a categorical association test).", "",
                f"G={g['G']:.2f}, df={g['df']}, {_fmt_p_txt(g['p'])}, "
                f"Cramér's V={g['cramers_v']:.2f}, cells with expected<5: {g['n_low_expected']}."])
    except Exception as e:
        print(f"  [stats] 08 ligand-element G-test skipped: {e}")
    ax.set_title(_title(out, title), fontsize=10 if "\n" in title else 12)
    ax.legend(title="ligand atom", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout(); _finalize_fig(fig, out)
    tab.to_csv(out.with_suffix(".csv"))


# ── 09-15: deeper crystal-vs-pose comparison ───────────────────────────────
# These extend the three native-reference figures (04/04b/05) with residue- and
# interaction-type-resolved recovery, a matched/missed/spurious decomposition, a
# geometry-vs-chemistry (RMSD-vs-F1) cross-check, an interaction-based pose
# selector, and per-pair fingerprint overlays. All consume the same residue-level
# ``inter`` / ``crystal`` frames; only 14 also needs the pose_comparison RMSD table.

def _method_colors(order) -> dict:
    """Per-tool colour map from the shared docking-tool palette (see TOOL_COLORS)."""
    return {m: TOOL_COLORS.get(m) for m in order}


def build_recovery_detail(inter: pd.DataFrame, crystal: pd.DataFrame,
                          copy_sel: "CrystalCopySelector | None" = None):
    """Per-pose typed & residue-contact recovery detail vs the crystal fingerprint.

    ``copy_sel`` (nearest-copy convention) scores each pose against the crystal copy it
    is nearest to; the extra ``crystal_copy`` column records that copy index.

    Returns ``(per_pose, per_type)``:
      * ``per_pose``  — one row per (method, protein, ligand, pose_name, pose_rank):
        typed matched/missed/spurious counts (tp/fn/fp), precision/recall/F1,
        Jaccard, and the residue-level (any-type) counterparts (tp_res,
        recall_contact) so the strict-vs-loose gap (fig 12) can be measured.
      * ``per_type`` — one row per (…pose…, interaction_type): typed tp/fn/fp
        restricted to that interaction type, plus the native prevalence, feeding
        the interaction-type-resolved recall/precision (fig 09).
    """
    cryst_typed, cryst_contact, cryst_by_type = {}, {}, {}

    def _by_type(t):
        by = {}
        for (it, rn, rnum, ch) in t:
            by.setdefault(it, set()).add((rn, rnum, ch))
        return by

    if copy_sel is None:
        for (p, l), g in crystal.groupby(["protein", "ligand"]):
            t = pose_fingerprint(g, "typed")
            cryst_typed[(p, l)] = t
            cryst_contact[(p, l)] = pose_fingerprint(g, "contact")
            cryst_by_type[(p, l)] = _by_type(t)

    has_rank = "pose_rank" in inter.columns
    keys = ["method", "protein", "ligand", "pose_name"] + (["pose_rank"] if has_rank else [])
    pose_rows, type_rows = [], []
    for k, g in inter.groupby(keys, sort=False):
        if has_rank:
            method, p, l, pose, rank = k
        else:
            (method, p, l, pose), rank = k, np.nan
        if copy_sel is not None:
            if not copy_sel.has_pair(p, l):
                continue
            kcopy = copy_sel.copy_of(method, p, l, pose)
            ref = copy_sel.fp(p, l, kcopy, "typed")
            ref_c = copy_sel.fp(p, l, kcopy, "contact")
            ref_by_type = cryst_by_type.setdefault((p, l, kcopy), _by_type(ref))
        else:
            kcopy = None
            ref = cryst_typed.get((p, l))
            if ref is None:
                continue
            ref_c = cryst_contact.get((p, l), set())
            ref_by_type = cryst_by_type.get((p, l), {})
        fp = pose_fingerprint(g, "typed")
        fp_c = pose_fingerprint(g, "contact")
        tp, fp_n, fn = len(fp & ref), len(fp - ref), len(ref - fp)
        prec = tp / len(fp) if fp else np.nan
        rec = tp / len(ref) if ref else np.nan
        f1 = (2 * prec * rec / (prec + rec)) if (prec and rec) else 0.0
        tp_c = len(fp_c & ref_c)
        pose_rows.append({
            "method": method, "protein": p, "ligand": l, "pose_name": pose,
            "pose_rank": rank, "n_native": len(ref), "n_pred": len(fp),
            "tp": tp, "fn": fn, "fp": fp_n, "precision": prec, "recall": rec,
            "f1": f1, "jaccard_vs_crystal": jaccard(fp, ref),
            "n_native_res": len(ref_c), "tp_res": tp_c,
            "recall_contact": (tp_c / len(ref_c) if ref_c else np.nan),
            **({"crystal_copy": kcopy} if copy_sel is not None else {}),
        })
        pose_by_type: dict[str, set] = {}
        for (it, rn, rnum, ch) in fp:
            pose_by_type.setdefault(it, set()).add((rn, rnum, ch))
        for it in set(pose_by_type) | set(ref_by_type):
            ps, rs = pose_by_type.get(it, set()), ref_by_type.get(it, set())
            type_rows.append({
                "method": method, "protein": p, "ligand": l, "pose_name": pose,
                "pose_rank": rank, "interaction_type": it, "tp": len(ps & rs),
                "fn": len(rs - ps), "fp": len(ps - rs), "n_native_type": len(rs),
            })
    return pd.DataFrame(pose_rows), pd.DataFrame(type_rows)


def plot_type_resolved_recovery(per_type: pd.DataFrame, order, out: Path,
                                depth: "int | None" = None) -> None:
    """09 — native recall & precision resolved per interaction type (fig 04 split by type).

    ``depth`` bounds which ranked poses contribute: ``depth=1`` uses only each method's
    top-ranked pose; ``depth=5`` pools its top-5; ``None`` pools all. n = the crystal's
    native contacts of that type across complexes (pose-count-independent).
    """
    if per_type.empty:
        print("  [skip] type-resolved recovery — no data"); return
    if depth is not None and "pose_rank" in per_type.columns:
        keep = (per_type[["method", "protein", "ligand", "pose_name", "pose_rank"]]
                .drop_duplicates().sort_values("pose_rank")
                .groupby(["method", "protein", "ligand"], sort=False).head(depth))
        per_type = per_type.merge(keep[["method", "protein", "ligand", "pose_name"]],
                                  on=["method", "protein", "ligand", "pose_name"], how="inner")
    agg = (per_type.groupby(["method", "interaction_type"])[["tp", "fn", "fp"]].sum()
           .reset_index())
    agg["recall"] = agg["tp"] / (agg["tp"] + agg["fn"]).replace(0, np.nan)
    agg["precision"] = agg["tp"] / (agg["tp"] + agg["fp"]).replace(0, np.nan)
    # Wilson 95% CI per cell (recall: k=tp, n=tp+fn; precision: k=tp, n=tp+fp).
    if _HAVE_STATS:
        rc = [su.wilson_ci(int(t), int(t + f)) for t, f in zip(agg["tp"], agg["fn"])]
        pc = [su.wilson_ci(int(t), int(t + p)) for t, p in zip(agg["tp"], agg["fp"])]
        agg["recall_lo"] = [c[0] for c in rc]; agg["recall_hi"] = [c[1] for c in rc]
        agg["precision_lo"] = [c[0] for c in pc]; agg["precision_hi"] = [c[1] for c in pc]
    # native prevalence per COMPLEX (dedup poses) so n is independent of the depth bucket
    native = (per_type.drop_duplicates(["protein", "ligand", "interaction_type"])
              .groupby("interaction_type")["n_native_type"].sum())
    types = [t for t in INTERACTION_TYPES if native.get(t, 0) > 0]
    methods = [m for m in order if m in set(agg["method"])]
    rec = agg.pivot(index="interaction_type", columns="method", values="recall").reindex(index=types, columns=methods)
    prc = agg.pivot(index="interaction_type", columns="method", values="precision").reindex(index=types, columns=methods)

    # Cell annotation: the point value only. The Wilson 95% CI per cell moves to the
    # sidecar .txt (written below) so the heatmap cells stay uncluttered.
    def _cell_annot(metric: str):
        vpiv = agg.pivot(index="interaction_type", columns="method", values=metric).reindex(index=types, columns=methods)
        out = np.empty(vpiv.shape, dtype=object)
        for r in range(vpiv.shape[0]):
            for c in range(vpiv.shape[1]):
                v = vpiv.iloc[r, c]
                out[r, c] = "" if pd.isna(v) else f"{v:.2f}"
        return out
    rec_ann, prc_ann = _cell_annot("recall"), _cell_annot("precision")
    ylabels = [f"{pretty_itype(t)}  (n={int(native.get(t, 0))})" for t in types]
    # low-n / saturated flags kept as sidecar metadata (no longer stamped on the figure)
    LOWN = 200
    lown = {t for t in types if int(native.get(t, 0)) < LOWN}
    satur = {t for t in types
             if (agg["interaction_type"] == t).any()
             and (agg.loc[agg["interaction_type"] == t, "recall"].fillna(0) >= 0.999).all()
             and (agg.loc[agg["interaction_type"] == t, "precision"].fillna(0) >= 0.999).all()}
    if _HAVE_STATS:
        _record(out.stem, {
            "unit": "pooled typed contacts per (interaction type, tool): recall=tp/(tp+fn), "
                    "precision=tp/(tp+fp), Wilson 95% CI per cell",
            "depth": depth, "low_n_flag_threshold_native_contacts": LOWN,
            "flagged_low_n": sorted(lown), "flagged_saturated": sorted(satur),
            "cells": [{"interaction_type": r.interaction_type, "method": pretty_method(r.method),
                       "tp": int(r.tp), "fn": int(r.fn), "fp": int(r.fp),
                       "recall": (None if r.recall != r.recall else float(r.recall)),
                       "recall_ci": [r.recall_lo, r.recall_hi],
                       "precision": (None if r.precision != r.precision else float(r.precision)),
                       "precision_ci": [r.precision_lo, r.precision_hi]}
                      for r in agg.itertuples()],
            "notes": "Pooled over poses → each contact treated as independent "
                     "(pseudoreplication); read the CIs as optimistic screening bounds. "
                     "n in the row label = the crystal's native contacts of that type "
                     "(independent-unit count). flagged_low_n = types with n<LOWN (CI genuinely "
                     "wide); flagged_saturated = recall=precision=1.00 for every tool "
                     "(trivially recovered)."})
    fig, axes = plt.subplots(
        1, 2, figsize=(max(11, 2.4 * len(methods) + 6), max(5.5, 0.6 * len(types) + 2)))
    for i, (ax, mat, ann, name) in enumerate(((axes[0], rec, rec_ann, "Recall"),
                                              (axes[1], prc, prc_ann, "Precision"))):
        sns.heatmap(mat, annot=ann, fmt="", annot_kws={"fontsize": 8}, cmap="viridis",
                    vmin=0, vmax=1, linewidths=0.5, linecolor="white",
                    xticklabels=[pretty_method(m) for m in mat.columns],
                    yticklabels=(ylabels if i == 0 else False),   # labels once, on the left
                    cbar_kws={"label": "fraction", "shrink": 0.8}, ax=ax)
        ax.set_title(name, loc="center", pad=10, fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("Interaction type   (n = the crystal's native contacts of that type)"
                      if i == 0 else "")
        plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")
        plt.setp(ax.get_yticklabels(), rotation=0)
    _label_panels(axes)
    pose_desc = ("top-1 pose" if depth == 1
                 else f"top-{depth} poses pooled" if depth else "all poses pooled")
    fig.suptitle(_title(out,
                        "Native-interaction recovery resolved by interaction type  —  " + pose_desc + "\n"
                        "(A) recall = native contacts reproduced;   (B) precision = predicted contacts that are native"),
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _finalize_fig(fig, out)
    agg.round(4).to_csv(out.with_suffix(".csv"), index=False)
    try:
        pose_desc_txt = ("top-1 pose" if depth == 1
                         else f"top-{depth} poses pooled" if depth else "all poses pooled")
        have_ci = "recall_lo" in agg.columns

        def _cell(r, val, lo, hi):
            if pd.isna(r[val]):
                return "n/a"
            if have_ci and not pd.isna(r[lo]):
                return f"{r[val]:.2f} [{r[lo]:.2f}, {r[hi]:.2f}]"
            return f"{r[val]:.2f}"

        lines = [f"Native-interaction recovery resolved by interaction type ({pose_desc_txt}).",
                 "Pooled typed contacts per (interaction type, tool): recall=tp/(tp+fn), "
                 "precision=tp/(tp+fp), with Wilson 95% CI per cell.",
                 "Poses pooled -> each contact treated as independent (pseudoreplication); "
                 "read the CIs as optimistic screening bounds.",
                 "n (per type) = the crystal's native contacts of that type across complexes.", ""]
        for t in types:
            lines.append(f"{pretty_itype(t)} (n={int(native.get(t, 0))}):")
            for m in methods:
                row = agg[(agg["interaction_type"] == t) & (agg["method"] == m)]
                if row.empty:
                    continue
                r = row.iloc[0]
                lines.append(f"  {pretty_method(m)}: "
                             f"recall {_cell(r, 'recall', 'recall_lo', 'recall_hi')}; "
                             f"precision {_cell(r, 'precision', 'precision_lo', 'precision_hi')} "
                             f"(tp={int(r['tp'])}, fn={int(r['fn'])}, fp={int(r['fp'])})")
        if lown:
            lines += ["", f"Low-n interaction types (native contacts < {LOWN}; CIs genuinely "
                          "wide): " + ", ".join(pretty_itype(t) for t in sorted(lown))]
        if satur:
            lines += ["", "Saturated types (recall=precision=1.00 for every tool — trivially "
                          "recovered): " + ", ".join(pretty_itype(t) for t in sorted(satur))]
        _write_stats_txt(out, f"09 — Type-resolved native recovery ({pose_desc_txt})", lines)
    except Exception as e:
        print(f"  [stats] 09 type-resolved txt skipped: {e}")


def plot_contact_decomposition(per_pose: pd.DataFrame, order, out: Path,
                               depth: int = 5) -> None:
    """10 — mean matched / missed / spurious contacts per pose, per method (stacked).

    ``depth`` bounds which ranked poses feed the average: ``depth=1`` scores only each
    method's top-ranked pose (the pick you'd actually use); ``depth=5`` averages over the
    method's top-5 poses (typical-pose behaviour). Poses are the lowest ``depth`` by
    pose_rank per (method, complex).
    """
    if per_pose.empty:
        print("  [skip] contact decomposition — no data"); return
    sel = (per_pose.sort_values("pose_rank")
           .groupby(["method", "protein", "ligand"], sort=False).head(depth))
    methods = [m for m in order if m in set(sel["method"])]
    means = sel.groupby("method")[["tp", "fn", "fp"]].mean().reindex(methods)
    x = np.arange(len(methods))
    fig, ax = plt.subplots(figsize=(max(7, 1.4 * len(methods)), 5.5))
    ax.bar(x, means["tp"], 0.6, label="matched (recovered native contact)", color="#55A868")
    ax.bar(x, means["fn"], 0.6, bottom=means["tp"], label="missed (native contact absent)", color="#C44E52")
    ax.bar(x, means["fp"], 0.6, bottom=means["tp"] + means["fn"],
           label="spurious (non-native contact invented)", color="#8172B3")
    ax.set_xticks(x); ax.set_xticklabels([pretty_method(m) for m in methods], rotation=20, ha="right")
    ax.set_ylabel("Mean number of typed contacts per pose")
    pose_desc = "top-1 pose" if depth == 1 else f"mean over top-{depth} poses"
    ax.set_title(_title(out,
                        f"Matched / missed / spurious interaction contacts vs the crystal — {pose_desc}\n"
                        "(typed contact = interaction type + residue)"))
    ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out)
    means.round(3).to_csv(out.with_suffix(".csv"))


def plot_contact_decomposition_by_rank(per_pose: pd.DataFrame, order, out: Path,
                                       top_n: int = 5) -> None:
    """10b — matched / missed / spurious contacts vs pose rank, as lines.

    Line version of fig 10: x = pose rank k (1..top_n), y = mean typed contacts per pose.
    Colour encodes the tool (AutoDock / DiffDock / EquiBind); line style encodes the contact
    class (matched solid / missed dashed / spurious dotted). Two legends (tool, contact class)
    keep the tools×3 lines legible.
    """
    if per_pose.empty or "pose_rank" not in per_pose.columns:
        print("  [skip] contact decomposition by rank — no data"); return
    df = per_pose.copy()
    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="coerce")
    df = df[df["pose_rank"].between(1, top_n)]
    methods = [m for m in order if m in set(df["method"])]
    if not methods:
        return
    # Stats: per-complex Kendall τ(pose rank, count) → one-sample Wilcoxon vs 0, per tool,
    # BH-FDR across the tool × {matched, missed, spurious} family. Matched (tp) declining
    # with rank (τ<0) — or missed/spurious rising (τ>0) — means the ranking orders quality.
    trend: dict = {}
    try:
        trend = _rank_trend_stats(df, order, ["tp", "fn", "fp"],
                                  "10b_contact_decomposition_by_rank",
                                  "per-complex Kendall tau(pose rank, contact count); one per complex")
    except Exception as e:
        print(f"  [stats] 10b rank-trend tests skipped: {e}")
    classes = [("tp", "matched (recovered native)", "-"),
               ("fn", "missed (native absent)", "--"),
               ("fp", "spurious (non-native invented)", ":")]
    colors = _method_colors(order)   # AutoDock blue, DiffDock orange, EquiBind green
    fig, ax = plt.subplots(figsize=(11.5, 6))
    for m in methods:
        color = colors.get(m)
        agg = df[df["method"] == m].groupby("pose_rank")[["tp", "fn", "fp"]].mean()
        for col, _, ls in classes:
            ax.plot(agg.index, agg[col], ls=ls, color=color, lw=2, marker="o", ms=4)
    ax.set_xticks(range(1, top_n + 1))
    ax.set_xlabel("Pose rank k (1 = tool's top-ranked pose)")
    ax.set_ylabel("Mean number of typed contacts per pose")
    ax.set_ylim(bottom=0); ax.grid(alpha=0.3); ax.set_axisbelow(True)
    # Title at the very top; the two legends sit between it and the plot as two
    # centered horizontal rows (Tool above Contact class) so neither overlaps the other.
    fig.suptitle(_title(out,
                        "Matched / missed / spurious interaction contacts vs the crystal, by pose rank\n"
                        "(typed contact = interaction type + residue)"), fontsize=11, y=0.98)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.74, bottom=0.10)
    from matplotlib.lines import Line2D
    tool_handles = [Line2D([0], [0], color=colors.get(m), lw=3) for m in methods]
    class_handles = [Line2D([0], [0], color="0.35", lw=2, ls=ls) for _, _, ls in classes]
    fig.legend(tool_handles, [pretty_method(m) for m in methods], title="Tool",
               loc="upper center", bbox_to_anchor=(0.5, 0.905), ncol=min(len(methods), 6),
               fontsize=9, frameon=True)
    fig.legend(class_handles, [lab for _, lab, _ in classes], title="Contact class",
               loc="upper center", bbox_to_anchor=(0.5, 0.83), ncol=len(classes),
               fontsize=9, frameon=True)
    _finalize_fig(fig, out, caption=False)
    (df.groupby(["method", "pose_rank"])[["tp", "fn", "fp"]].mean().round(3)
     .to_csv(out.with_suffix(".csv")))
    try:
        tlines = ["Rank trend: per-complex Kendall tau(pose rank, contact count) -> one-sample "
                  "Wilcoxon signed-rank vs 0, per tool, BH-FDR across the tool x class family.",
                  "matched (tp) tau<0 (falls with rank) or missed/spurious (fn/fp) tau>0 (rises) "
                  "=> the tool's ranking orders pose quality.", ""]
        classlabels = [("tp", "matched"), ("fn", "missed"), ("fp", "spurious")]
        for m in methods:
            tlines.append(f"{pretty_method(m)}:")
            for col, name in classlabels:
                e = (trend.get(m) or {}).get(col, {})
                if "median_tau" in e:
                    bh = (f", BH q={su.fmt_p(e['p_bh'])} {e.get('star_bh', '')}".rstrip()
                          if e.get("p_bh") is not None else "")
                    tlines.append(f"  {name}: median tau={e['median_tau']:+.2f}, "
                                  f"rank-biserial={e.get('rank_biserial'):+.2f}, "
                                  f"{_fmt_p_txt(e.get('p_raw'))}{bh}  "
                                  f"(n={e.get('n_complexes')} complexes)")
                else:
                    tlines.append(f"  {name}: {e.get('note', 'n/a')}")
        _write_stats_txt(out, "10b — Contact decomposition by pose rank", tlines)
    except Exception as e:
        print(f"  [stats] 10b rank-trend txt skipped: {e}")


def plot_residue_confusion(inter: pd.DataFrame, crystal: pd.DataFrame, order,
                           out: Path, top_n: int = 20,
                           copy_sel: "CrystalCopySelector | None" = None) -> None:
    """11 — per-residue recovery of native contacts and per-residue spurious contacts.

    Residues are pooled by resname+resnum across complexes (as in fig 03). Left panel:
    fraction of a method's poses that recover each top native hot-spot residue (only
    counting poses whose own crystal contacts that residue). Right panel: the top
    residues each method contacts that its crystal does NOT (hallucinated contacts),
    as a fraction of that method's poses.
    """
    def _res(df):
        return df["resname"].astype(str) + df["resnum"].astype(str)
    pc = inter.copy(); pc["res"] = _res(pc)
    pose_res = (pc.groupby(["method", "protein", "ligand", "pose_name"])["res"]
                .apply(set).reset_index())
    cryst = crystal.copy(); cryst["res"] = _res(cryst)
    if copy_sel is None:
        native_by_pair = cryst.groupby(["protein", "ligand"])["res"].apply(set).to_dict()
        native_freq = cryst.drop_duplicates(["protein", "ligand", "res"])["res"].value_counts()
        native_by_copy = None
    else:
        # nearest-copy: a pose's native residue set is that of the copy it is nearest to;
        # the hot-spot ranking pools, per complex, the residues native in ANY copy some
        # profiled pose selected (a residue counts once per complex).
        ck = cryst["copy_index"].fillna(0).astype(int)
        native_by_copy = cryst.assign(_k=ck).groupby(["protein", "ligand", "_k"])["res"].apply(set).to_dict()
        pose_res["_k"] = [copy_sel.copy_of(m, p, l, n) for m, p, l, n in
                          pose_res[["method", "protein", "ligand", "pose_name"]].itertuples(index=False, name=None)]
        used = pose_res[["protein", "ligand", "_k"]].drop_duplicates()
        used = used[[copy_sel.has_pair(p, l) for p, l in zip(used["protein"], used["ligand"])]]
        cu = cryst.assign(_k=ck).merge(used, on=["protein", "ligand", "_k"], how="inner")
        native_freq = cu.drop_duplicates(["protein", "ligand", "res"])["res"].value_counts()
    methods = [m for m in order if m in set(pose_res["method"])]
    poses_per_method = pose_res.groupby("method").size()

    recov_num: dict = {}; recov_den: dict = {}; spur_num: dict = {}
    for _, r in pose_res.iterrows():
        m = r["method"]
        if native_by_copy is not None:
            nat = native_by_copy.get((r["protein"], r["ligand"], r["_k"]), set())
        else:
            nat = native_by_pair.get((r["protein"], r["ligand"]), set())
        contacted = r["res"]
        for res in nat:                       # recovery denominator = poses whose pair has res native
            recov_den[(m, res)] = recov_den.get((m, res), 0) + 1
            if res in contacted:
                recov_num[(m, res)] = recov_num.get((m, res), 0) + 1
        for res in contacted - nat:           # spurious = contacted but not native for this pair
            spur_num[(m, res)] = spur_num.get((m, res), 0) + 1

    top_native = list(native_freq.head(top_n).index)
    rec_mat = pd.DataFrame(index=top_native, columns=methods, dtype=float)
    for res in top_native:
        for m in methods:
            den = recov_den.get((m, res), 0)
            rec_mat.loc[res, m] = (recov_num.get((m, res), 0) / den) if den else np.nan

    spur_tot = {}
    for (m, res), n in spur_num.items():
        spur_tot[res] = spur_tot.get(res, 0) + n
    top_spur = [res for res, _ in sorted(spur_tot.items(), key=lambda kv: kv[1], reverse=True)[:top_n]]
    spur_mat = pd.DataFrame(index=top_spur, columns=methods, dtype=float)
    for res in top_spur:
        for m in methods:
            spur_mat.loc[res, m] = spur_num.get((m, res), 0) / poses_per_method.get(m, np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(max(10, 1.7 * len(methods) + 6), max(6, 0.35 * top_n + 1.5)))
    sns.heatmap(rec_mat.astype(float), annot=True, fmt=".2f", cmap="Greens", vmin=0, vmax=1,
                xticklabels=[pretty_method(m) for m in methods],
                cbar_kws={"label": "fraction of poses recovering the contact"}, ax=axes[0])
    axes[0].set_title("Native recovery", loc="center", pad=16)
    axes[0].set_ylabel("Residue (name + number)"); axes[0].set_xlabel("")
    if top_spur:
        sns.heatmap(spur_mat.astype(float), annot=True, fmt=".2f", cmap="Reds", vmin=0,
                    xticklabels=[pretty_method(m) for m in methods],
                    cbar_kws={"label": "fraction of poses making the non-native contact"}, ax=axes[1])
        axes[1].set_title("Spurious contacts", loc="center", pad=16)
    else:
        axes[1].set_visible(False)
    axes[1].set_ylabel("Residue (name + number)"); axes[1].set_xlabel("")
    for ax in axes:
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    _label_panels(axes)
    fig.suptitle(_title(out,
                        f"Per-residue native-contact recovery (top {len(top_native)} hot-spot residues) "
                        f"and hallucinated non-native contacts (top {len(top_spur)})"), fontsize=13)
    fig.tight_layout(); _finalize_fig(fig, out)
    rec_mat.round(3).to_csv(out.with_name(out.stem + "_recovery.csv"))
    spur_mat.round(3).to_csv(out.with_name(out.stem + "_spurious.csv"))


def plot_typed_vs_loose(per_pose: pd.DataFrame, order, out: Path) -> None:
    """12 — strict typed recall vs loose residue-contact recall (right residue, wrong bond)."""
    if per_pose.empty:
        print("  [skip] typed-vs-loose gap — no data"); return
    methods = [m for m in order if m in set(per_pose["method"])]
    strict = per_pose.groupby("method")["recall"].mean().reindex(methods)
    loose = per_pose.groupby("method")["recall_contact"].mean().reindex(methods)
    # Stats: per-complex paired Wilcoxon of loose vs strict recall (within-pose paired,
    # aggregated to one value per complex), BH-FDR across tools. Star flags a tool whose
    # loose>strict gap is significant.
    stars: dict = {}
    try:
        agg = per_pose.groupby(["method", "protein", "ligand"])[["recall", "recall_contact"]].mean()
        per_m, ps, keys = {}, [], []
        for m in methods:
            try:
                sub = agg.xs(m, level="method")
            except KeyError:
                continue
            xr = sub["recall"].to_numpy(float); yv = sub["recall_contact"].to_numpy(float)
            mask = np.isfinite(xr) & np.isfinite(yv)
            if _HAVE_STATS and mask.sum() >= _MIN_UNITS:
                rb, p, npair = su.wilcoxon_rankbiserial(yv[mask], xr[mask])  # loose - strict
                per_m[m] = {"rank_biserial": rb, "p_raw": p, "n": int(mask.sum())}
                if p == p:
                    ps.append(p); keys.append(m)
        if ps and _HAVE_STATS:
            for m, q in zip(keys, su.bh_fdr(ps)):
                per_m[m]["p_bh"] = float(q); per_m[m]["star_bh"] = su.p_stars(q)
                stars[m] = su.p_stars(q)
        _record("12_typed_vs_loose", {
            "unit": "per-complex mean recall, loose vs strict (Wilcoxon signed-rank)",
            "family": "per-tool Wilcoxon (BH-FDR)",
            "per_method": {pretty_method(m): v for m, v in per_m.items()}})
        lines = ["Test: per-complex paired Wilcoxon signed-rank of loose vs strict native "
                 "recall (aggregated to one value per complex), BH-FDR across tools.",
                 "gap = loose − strict mean recall (poses hitting the native residue but "
                 "forming a different interaction type).", ""]
        for m in methods:
            gap = loose.get(m, np.nan) - strict.get(m, np.nan)
            v = per_m.get(m)
            if v:
                lines.append(f"  {pretty_method(m)}: gap={gap:+.2f}, "
                             f"rank-biserial={v['rank_biserial']:+.2f}, "
                             f"{_fmt_p_txt(v.get('p_bh', v.get('p_raw')))} (BH), n={v['n']}")
            else:
                lines.append(f"  {pretty_method(m)}: gap={gap:+.2f}  (no test — small n)")
        _write_stats_txt(out, "12 — Typed vs loose native recall", lines)
    except Exception as e:
        print(f"  [stats] 12 typed-vs-loose test skipped: {e}")
    x = np.arange(len(methods)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(methods)), 5.5))
    ax.bar(x - w / 2, strict.values, w, label="typed recall (interaction type + residue must match)", color="#4C72B0")
    ax.bar(x + w / 2, loose.values, w, label="loose recall (residue contacted, any interaction type)", color="#DD8452")
    for xi, m in enumerate(methods):
        s, ll = strict.values[xi], loose.values[xi]
        if np.isfinite(s) and np.isfinite(ll):
            ax.annotate(f"gap {ll - s:.2f}", (xi, ll),
                        textcoords="offset points", xytext=(0, 4), ha="center",
                        fontsize=7, color="#333333")
    ax.set_xticks(x); ax.set_xticklabels([pretty_method(m) for m in methods], rotation=20, ha="right")
    ax.set_ylabel("Mean recall of native contacts"); ax.set_ylim(0, 1)
    ax.set_title(_title(out,
                        "Right residue, wrong chemistry: typed vs loose native recall\n"
                        "(gap = poses that hit the native residue but form a different interaction type)"),
                 fontsize=11)
    ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out)


def plot_native_f1_oracle(per_pose: pd.DataFrame, order, out: Path) -> None:
    """13 — best-achievable native-F1 in a pose set vs the tool's own top-ranked pose."""
    if per_pose.empty:
        print("  [skip] native-F1 oracle — no data"); return
    rows = []
    for (m, p, l), g in per_pose.groupby(["method", "protein", "ligand"]):
        oracle = g["f1"].max()
        top1 = g.sort_values(["pose_rank", "pose_name"])["f1"].iloc[0]
        rows.append({"method": m, "oracle_f1": oracle, "top1_f1": top1})
    d = pd.DataFrame(rows)
    methods = [m for m in order if m in set(d["method"])]
    oracle = d.groupby("method")["oracle_f1"].mean().reindex(methods)
    top1 = d.groupby("method")["top1_f1"].mean().reindex(methods)
    # Stats: per-complex paired Wilcoxon of oracle vs top-1 native-F1 (one value per
    # complex), BH-FDR across tools. Star flags a tool whose own ranking leaves a
    # significant recovery gap.
    stars: dict = {}
    try:
        per_m, ps, keys = {}, [], []
        for m in methods:
            sub = d[d["method"] == m]
            o = sub["oracle_f1"].to_numpy(float); t1 = sub["top1_f1"].to_numpy(float)
            mask = np.isfinite(o) & np.isfinite(t1)
            if _HAVE_STATS and mask.sum() >= _MIN_UNITS:
                rb, p, npair = su.wilcoxon_rankbiserial(o[mask], t1[mask])
                per_m[m] = {"rank_biserial": rb, "p_raw": p, "n": int(mask.sum())}
                if p == p:
                    ps.append(p); keys.append(m)
        if ps and _HAVE_STATS:
            for m, q in zip(keys, su.bh_fdr(ps)):
                per_m[m]["p_bh"] = float(q); per_m[m]["star_bh"] = su.p_stars(q)
                stars[m] = su.p_stars(q)
        _record("13_native_f1_oracle", {
            "unit": "per-complex oracle vs top-1 native-F1 (Wilcoxon signed-rank)",
            "family": "per-tool Wilcoxon (BH-FDR)",
            "per_method": {pretty_method(m): v for m, v in per_m.items()}})
        lines = ["Test: per-complex paired Wilcoxon signed-rank of oracle (best pose in set) "
                 "vs the tool's own top-1 native-F1 (one value per complex), BH-FDR across tools.",
                 "gap = recovery the tool's own ranking leaves on the table.", ""]
        for m in methods:
            o = oracle.get(m, np.nan); t1 = top1.get(m, np.nan)
            v = per_m.get(m)
            base = f"  {pretty_method(m)}: oracle F1={o:.2f}, top-1 F1={t1:.2f}, gap={o - t1:+.2f}"
            if v:
                lines.append(base + f", rank-biserial={v['rank_biserial']:+.2f}, "
                             f"{_fmt_p_txt(v.get('p_bh', v.get('p_raw')))} (BH), n={v['n']}")
            else:
                lines.append(base + "  (no test — small n)")
        _write_stats_txt(out, "13 — Native-F1 oracle vs top-ranked pose", lines)
    except Exception as e:
        print(f"  [stats] 13 native-F1 oracle test skipped: {e}")
    x = np.arange(len(methods)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(9, 2.0 * len(methods)), 5.5))
    ax.bar(x - w / 2, oracle.values, w, label="best native-F1 pose in the set (interaction oracle)", color="#55A868")
    ax.bar(x + w / 2, top1.values, w, label="native-F1 of the tool's own top-ranked pose", color="#937860")
    ax.set_xticks(x); ax.set_xticklabels([pretty_method(m) for m in methods], rotation=20, ha="right")
    ax.set_ylabel("Mean native-interaction recovery F1"); ax.set_ylim(0, 1)
    ax.set_title(_title(out,
                        "Interaction-recovery F1: oracle vs the tool's top-ranked pose\n"
                        "gap = recovery the tool's own ranking leaves on the table"), fontsize=11)
    ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out)
    d.groupby("method")[["oracle_f1", "top1_f1"]].mean().round(3).to_csv(out.with_suffix(".csv"))


def plot_rmsd_vs_recovery(per_pose: pd.DataFrame, metrics_csv: Path, order, out: Path) -> None:
    """14 — geometry vs chemistry: RMSD-to-crystal vs native-interaction F1, per pose.

    Joins the per-pose native recovery to the pose_comparison RMSD table on
    (method, protein, ligand, pose_name) — the identity key that survives across the
    two stages. Coverage per method is logged (AutoDock has some complexes with no
    RMSD row; EquiBind/DiffDock join on the exact variant PandaMap fingerprinted).
    """
    if per_pose.empty or not metrics_csv or not Path(metrics_csv).exists():
        print(f"  [skip] RMSD-vs-recovery — metrics table not found ({metrics_csv})"); return
    mt = pd.read_csv(metrics_csv, low_memory=False)
    need = {"method", "protein", "ligand", "pose_name", "rmsd"}
    if need - set(mt.columns):
        print(f"  [skip] RMSD-vs-recovery — {need - set(mt.columns)} missing from {metrics_csv}"); return
    key = ["method", "protein", "ligand", "pose_name"]
    mt = mt[key + ["rmsd"]].dropna(subset=["rmsd"]).drop_duplicates(key)
    j = per_pose.merge(mt, on=key, how="inner")
    if j.empty:
        print("  [skip] RMSD-vs-recovery — no (method,protein,ligand,pose_name) overlap"); return
    methods = [m for m in order if m in set(j["method"])]
    colors = _method_colors(order)
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    # Stats: per-method Spearman rho (RMSD vs native-F1) WITH a p-value, BH-FDR across
    # methods. Unit = per pose (pooled) — this is inherently a pose-by-pose relationship;
    # flagged in the JSON 'unit' field. A star on the legend = FDR<0.05.
    per_m: dict = {}
    ps, keys = [], []
    for m in methods:
        sub = j[j["method"] == m]
        cov = len(sub); tot = int((per_pose["method"] == m).sum())
        if _HAVE_STATS and len(sub) > 2:
            rho, p, nn = su.spearman(sub["rmsd"].to_numpy(float), sub["f1"].to_numpy(float))
        else:
            rho, p, nn = np.nan, np.nan, len(sub)
        per_m[m] = {"rho": rho, "p_raw": p, "n_poses": nn, "coverage": f"{cov}/{tot}"}
        if p == p:
            ps.append(p); keys.append(m)
        ax.scatter(sub["rmsd"], sub["f1"], s=12, alpha=0.35, color=colors.get(m), edgecolor="none")
    if ps and _HAVE_STATS:
        for m, q in zip(keys, su.bh_fdr(ps)):
            per_m[m]["p_bh"] = float(q); per_m[m]["star_bh"] = su.p_stars(q)
    legend = [f"{pretty_method(m)}  (n={per_m[m]['coverage']} poses)" for m in methods]
    # Test 3 — do the tools differ in how tightly RMSD couples to recovery?
    # Fisher r-to-z on each method's Spearman rho, pairwise, Holm across the pairs.
    # Bonett-Wright SE for Spearman: var(z)=(1+rho^2/2)/(n-3). The methods share
    # complexes (different poses each), so this treats them as independent samples —
    # a screening approximation flagged in the sidecar; the coupling difference here
    # dwarfs that dependence.
    fisher: list = []
    fz = {m: (math.atanh(per_m[m]["rho"]), per_m[m]["n_poses"]) for m in methods
          if per_m[m]["rho"] == per_m[m]["rho"] and abs(per_m[m]["rho"]) < 1
          and (per_m[m]["n_poses"] or 0) > 3}
    fzp: list = []
    for a, b in itertools.combinations([m for m in methods if m in fz], 2):
        za, na = fz[a]; zb, nb = fz[b]
        se = math.sqrt((1 + per_m[a]["rho"] ** 2 / 2) / (na - 3)
                       + (1 + per_m[b]["rho"] ** 2 / 2) / (nb - 3))
        z = (za - zb) / se
        p = math.erfc(abs(z) / math.sqrt(2))          # two-sided
        fisher.append({"pair": [pretty_method(a), pretty_method(b)],
                       "delta_rho": per_m[a]["rho"] - per_m[b]["rho"],
                       "z": z, "p_raw": p})
        fzp.append(p)
    if fzp and _HAVE_STATS:
        for e, q in zip(fisher, su.holm(fzp)):
            e["p_holm"] = float(q); e["star_holm"] = su.p_stars(q)

    # Test 4 — 2 Å threshold framing: is native-F1 higher for sub-2 Å poses?
    # Mann-Whitney U on F1 (RMSD<2 vs >=2) + Cliff's delta per method, Holm across methods.
    thr: dict = {}
    thr_keys, thr_ps = [], []
    for m in methods:
        sub = j[j["method"] == m]
        lo = sub.loc[sub["rmsd"] < 2.0, "f1"].to_numpy(float)
        hi = sub.loc[sub["rmsd"] >= 2.0, "f1"].to_numpy(float)
        mw = su.mannwhitney_cliffs(lo, hi) if (_HAVE_STATS and lo.size and hi.size) else None
        if mw:
            thr[m] = {"n_lt2": mw["n_a"], "n_ge2": mw["n_b"],
                      "median_f1_lt2": mw["medians"][0], "median_f1_ge2": mw["medians"][1],
                      "U": mw["U"], "p_raw": mw["p"],
                      "cliffs_delta": mw["cliffs_delta"], "delta_ci": mw["delta_ci"]}
            thr_keys.append(m); thr_ps.append(mw["p"])
    if thr_ps and _HAVE_STATS:
        for m, q in zip(thr_keys, su.holm(thr_ps)):
            thr[m]["p_holm"] = float(q); thr[m]["star_holm"] = su.p_stars(q)

    _record("14_rmsd_vs_recovery", {
        "unit": "per pose (pooled): RMSD-to-crystal vs native-F1",
        "family": "per-method Spearman (BH-FDR) + between-method Fisher-z on rho (Holm) "
                  "+ 2A-threshold Mann-Whitney/Cliff's delta (Holm)",
        "per_method": {pretty_method(m): per_m[m] for m in methods},
        "between_method_fisher_z": fisher,
        "threshold_2A": {pretty_method(m): thr[m] for m in thr},
        "notes": "Fisher-z treats methods as independent samples (they share complexes "
                 "but different poses) — screening approximation. Split at RMSD=2 Å; "
                 "Cliff's delta = P(F1|<2Å > F1|>=2Å) - P(<)."})

    ax.axvline(2.0, color="#888888", ls="--", lw=1)
    ax.annotate("2 Å", (2.05, 0.985), color="#666666", fontsize=8, ha="left", va="top")
    ax.set_xlabel("RMSD to crystal ligand (Å)")
    ax.set_ylabel("Native-interaction recovery F1")
    ax.set_title(_title(out, "Geometry vs chemistry: does low RMSD predict recovering native contacts?"),
                 fontsize=12)
    # Cap the RMSD axis at 10 Å — far-off decoy poses beyond this only compress the
    # informative 0–10 Å range (per-method Spearman ρ above is still over all poses).
    ax.set_xlim(0, 10); ax.set_ylim(-0.02, 1.02)
    ax.legend(legend, fontsize=8, loc="upper right"); ax.grid(alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out, caption=False)
    j[key + ["rmsd", "f1", "recall", "precision"]].to_csv(out.with_suffix(".csv"), index=False)
    # Compact per-figure stats sidecar (tidy long form) — the three tests behind the
    # annotations, so the panel's numbers are reproducible without parsing the pooled
    # interaction_stats.json.
    srows: list = []
    for m in methods:
        e = per_m[m]
        srows.append({"test": "spearman_rmsd_vs_f1", "subject": pretty_method(m),
                      "n": e.get("n_poses"), "estimate_name": "spearman_rho",
                      "estimate": e.get("rho"), "p_raw": e.get("p_raw"),
                      "p_adj": e.get("p_bh"), "p_adj_method": "BH-FDR",
                      "stars": e.get("star_bh", "")})
    for e in fisher:
        srows.append({"test": "fisher_z_between_rho",
                      "subject": f"{e['pair'][0]} vs {e['pair'][1]}",
                      "estimate_name": "delta_rho", "estimate": e["delta_rho"], "z": e["z"],
                      "p_raw": e["p_raw"], "p_adj": e.get("p_holm"), "p_adj_method": "Holm",
                      "stars": e.get("star_holm", "")})
    for m in methods:
        if m in thr:
            t = thr[m]
            srows.append({"test": "mannwhitney_f1_by_2A", "subject": pretty_method(m),
                          "n": t["n_lt2"] + t["n_ge2"], "n_lt2": t["n_lt2"], "n_ge2": t["n_ge2"],
                          "median_f1_lt2": t["median_f1_lt2"], "median_f1_ge2": t["median_f1_ge2"],
                          "estimate_name": "cliffs_delta", "estimate": t["cliffs_delta"],
                          "cliffs_delta_lo": t["delta_ci"][0], "cliffs_delta_hi": t["delta_ci"][1],
                          "U": t["U"], "p_raw": t["p_raw"], "p_adj": t.get("p_holm"),
                          "p_adj_method": "Holm", "stars": t.get("star_holm", "")})
    pd.DataFrame(srows).to_csv(out.with_name(out.stem + "_stats.csv"), index=False)
    try:
        tlines = ["Unit: per pose (pooled) — RMSD-to-crystal vs native-interaction F1.",
                  "Fisher-z treats methods as independent samples (they share complexes but "
                  "different poses) — a screening approximation.", "",
                  "Per-method Spearman correlation (BH-FDR across methods):"]
        for m in methods:
            e = per_m[m]
            bh = (f", BH q={su.fmt_p(e['p_bh'])} {e.get('star_bh', '')}".rstrip()
                  if e.get("p_bh") is not None else "")
            tlines.append(f"  {pretty_method(m)}: rho={e['rho']:+.2f}, "
                          f"{_fmt_p_txt(e.get('p_raw'))}{bh}  (n={e['coverage']} poses)")
        if fisher:
            tlines += ["", "Between-method coupling (Fisher r-to-z on rho, Holm):"]
            for e in fisher:
                tlines.append(f"  {e['pair'][0]} vs {e['pair'][1]}: Δrho={e['delta_rho']:+.2f}, "
                              f"z={e['z']:+.2f}, {_fmt_p_txt(e.get('p_holm'))} (Holm)")
        if thr:
            tlines += ["", "2 Å threshold (Mann-Whitney U on F1, <2Å vs ≥2Å; Cliff's delta, Holm):"]
            for m in methods:
                if m not in thr:
                    continue
                t = thr[m]
                tlines.append(f"  {pretty_method(m)}: median F1 {t['median_f1_lt2']:.2f} "
                              f"(<2Å, n={t['n_lt2']}) vs {t['median_f1_ge2']:.2f} "
                              f"(≥2Å, n={t['n_ge2']}); Cliff's delta={t['cliffs_delta']:+.2f} "
                              f"[{t['delta_ci'][0]:+.2f}, {t['delta_ci'][1]:+.2f}], "
                              f"{_fmt_p_txt(t.get('p_holm'))} (Holm)")
        _write_stats_txt(out, "14 — Geometry vs chemistry (RMSD vs native-F1)", tlines)
    except Exception as e:
        print(f"  [stats] 14 txt skipped: {e}")
    print(f"  RMSD-vs-recovery: joined {len(j)} poses "
          f"({', '.join(f'{m}={int((j.method==m).sum())}' for m in methods)}).")


def plot_pair_overlays(inter: pd.DataFrame, crystal: pd.DataFrame, per_pose: pd.DataFrame,
                       order, out: Path, n_pairs: int = 4, max_rows: int = 26,
                       copy_sel: "CrystalCopySelector | None" = None) -> None:
    """15 — per-complex crystal-vs-pose typed-contact overlays (case studies).

    Picks a spread of complexes (best / worst / most-divergent mean native-F1 across
    methods) and, for each, shows a present/absent grid: rows = native + method
    contacts, columns = crystal and each method's best (top-ranked) pose. Cells are
    coloured native-only / recovered / spurious / absent.
    """
    if per_pose.empty or crystal.empty:
        print("  [skip] pair overlays — no data"); return
    methods = [m for m in order if m in set(per_pose["method"])]
    pair_f1 = per_pose.groupby(["protein", "ligand"])["f1"].mean().dropna()
    spread = per_pose.groupby(["protein", "ligand"])["f1"].agg(lambda s: s.max() - s.min())
    if pair_f1.empty:
        print("  [skip] pair overlays — no per-pair F1"); return
    picks, seen = [], set()
    for cand in ([pair_f1.idxmax(), pair_f1.idxmin(), spread.idxmax()] +
                 list(pair_f1.sort_values().index[len(pair_f1) // 2: len(pair_f1) // 2 + 2])):
        if cand not in seen:
            picks.append(cand); seen.add(cand)
        if len(picks) >= n_pairs:
            break

    if copy_sel is None:
        cryst_typed = {(p, l): pose_fingerprint(g, "typed") for (p, l), g in crystal.groupby(["protein", "ligand"])}
    else:
        cryst_typed = {}   # nearest-copy: filled per picked pair from the shown poses' copies
    best_pose = (inter.sort_values(["pose_rank", "pose_name"])
                 .drop_duplicates(["method", "protein", "ligand"]))  # top-ranked pose id per (method,pair)

    cmap = matplotlib.colors.ListedColormap(["#F2F2F2", "#C44E52", "#55A868", "#8172B3"])
    #                                          absent      native-only recovered   spurious
    fig, axes = plt.subplots(1, len(picks), figsize=(max(4.5 * len(picks), 9), max(6, 0.28 * max_rows + 2)))
    axes = np.atleast_1d(axes)
    for ax, (p, l) in zip(axes, picks):
        cols, fps, refs = ["crystal"], [None], [None]
        for m in methods:
            row = best_pose[(best_pose["method"] == m) & (best_pose["protein"] == p) & (best_pose["ligand"] == l)]
            if row.empty:
                continue
            pose = row.iloc[0]["pose_name"]
            g = inter[(inter["method"] == m) & (inter["protein"] == p) &
                      (inter["ligand"] == l) & (inter["pose_name"] == pose)]
            cols.append(m); fps.append(pose_fingerprint(g, "typed"))
            # nearest-copy: each shown pose is judged against the copy it is nearest to
            refs.append(copy_sel.pose_fp(m, p, l, pose) if copy_sel is not None else None)
        if copy_sel is None:
            ref = cryst_typed.get((p, l), set())
        else:
            # the crystal column shows the union of the copies the shown poses selected
            ref = set().union(*[r for r in refs if r is not None]) if any(r is not None for r in refs) else set()
        fps[0] = ref
        refs = [ref if r is None else r for r in refs]
        contacts = sorted(set().union(*fps), key=lambda c: (c[0] not in {t for (t, *_ ) in ref}, str(c)))
        contacts = contacts[:max_rows]
        M = np.zeros((len(contacts), len(cols)))
        for j, (cn, fp, ref_j) in enumerate(zip(cols, fps, refs)):
            for i, c in enumerate(contacts):
                if c not in fp:
                    M[i, j] = 0
                elif cn == "crystal":
                    M[i, j] = 1                       # native (shown in crystal column)
                elif c in ref_j:
                    M[i, j] = 2                       # recovered native
                else:
                    M[i, j] = 3                       # spurious
        ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=3)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([pretty_method(c) for c in cols], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(contacts)))
        ax.set_yticklabels([f"{it}:{rn}{rnum}" for (it, rn, rnum, ch) in contacts], fontsize=6)
        ax.set_title(str(p), fontsize=9)
        _note_subtitle(out, f"{p}: mean native-F1 {pair_f1.get((p, l), np.nan):.2f}")
    handles = [matplotlib.patches.Patch(color=c, label=lbl) for c, lbl in
               (("#C44E52", "native contact (crystal)"), ("#55A868", "recovered by pose"),
                ("#8172B3", "spurious (non-native)"), ("#F2F2F2", "absent"))]
    # lift the legend into its own band (bbox y=0.05); the bottom 10% reserved by the
    # tight_layout rect below keeps it clear of the panels.
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.05),
               ncol=4, fontsize=8, frameon=False)
    _label_panels(axes)
    fig.suptitle(_title(out, "Crystal-vs-pose interaction overlays (top-ranked pose per method)"),
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.10, 1, 1)); _finalize_fig(fig, out)


# ── main ────────────────────────────────────────────────────────────────────

def _load_allowed_ids(path: Path) -> set[str]:
    """Load '<PDBID>_<LIG>' complex ids (one per line; '#' comments ignored)."""
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", "-c", type=Path, default=None,
                    help="Optional pandamap_config.yaml (the file run_pandamap.py "
                         "uses). Supplies in-dir (from its output_dir), benchmark-dir, "
                         "best_equibind_only and oracle_summary. Any CLI flag below "
                         "overrides the matching config value.")
    ap.add_argument("--in-dir", type=Path, default=None,
                    help="run_pandamap output dir (has pandamap_*.csv). "
                         "Default: config output_dir.")
    ap.add_argument("--out-dir", type=Path, default=None, help="default <in-dir>/report")
    ap.add_argument("--benchmark-dir", type=Path, default=None,
                    help="for crystal-ligand chemotype descriptors "
                         "(default: config benchmark_dir, else 'Data/PoseBuster Benchmark Set').")
    ap.add_argument("--top-residues", type=int, default=25)
    ap.add_argument("--poses-per-combo", type=int, default=None,
                    help="Top-N ranked poses per method × complex that the fingerprints "
                         "were built from (run_pandamap's poses_per_combo). Shown as a "
                         "footer on every figure and used to cap the depth-pooled charts. "
                         "Default: config 'poses_per_combo', else inferred from the data.")
    ap.add_argument("--best-equibind-only", action="store_true", default=None,
                    help="Keep only the single best-performing EquiBind variant "
                         "(highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%, "
                         "read from --oracle-summary) in all charts/CSVs, relabelled "
                         "'EquiBind*'. AutoDock/DiffDock/crystal are unaffected. "
                         "Overrides config 'best_equibind_only'.")
    ap.add_argument("--equibind-variant", default=None, metavar="SPEC",
                    help="Pin EquiBind to the variant matching SPEC (tokens split on '_' or "
                         "'/', e.g. 'gnina' or 'unguided_gnina') instead of oracle-ranking it; "
                         "kept in all charts/CSVs under its explicit variant label. Takes precedence "
                         "over --best-equibind-only. Overrides config 'equibind_variant'.")
    ap.add_argument("--best-diffdock-only", action="store_true", default=None,
                    help="Keep only the single best-performing DiffDock optimizer "
                         "variant (highest PB-Valid AND RMSD ≤ 2 Å = "
                         "oracle_pb_valid_and_rmsd2_%%, read from --oracle-summary) in all "
                         "charts/CSVs, relabelled 'DiffDock*'. AutoDock/EquiBind/crystal are "
                         "unaffected. "
                         "Overrides config 'best_diffdock_only'.")
    ap.add_argument("--diffdock-variant", choices=("raw", "smina", "gnina"), default=None,
                    help="Pin DiffDock to this optimizer variant (e.g. 'gnina') instead of "
                         "oracle-ranking it; kept in all charts/CSVs under its explicit variant label. "
                         "Takes precedence over --best-diffdock-only. Overrides config "
                         "'diffdock_variant'.")
    ap.add_argument("--oracle-summary", type=Path, default=None,
                    help="oracle_summary.csv from posebusters_pose_comparison.py, used "
                         "to pick the best EquiBind variant for --best-equibind-only. "
                         "Overrides config 'oracle_summary' "
                         f"(default: {DEFAULT_ORACLE_SUMMARY}).")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Restrict the report to the '<PDBID>_<LIG>' complex ids "
                         "listed in this file (one per line; '#' comments ok). "
                         "Overrides config 'ids_file'.")
    ap.add_argument("--per-pose-metrics", type=Path, default=None,
                    help="per_pose_metrics.csv from posebusters_pose_comparison.py "
                         "(per-pose RMSD-to-crystal) for the geometry-vs-recovery figure. "
                         f"Default: sibling of --oracle-summary, else {DEFAULT_PER_POSE_METRICS}.")
    mf.add_method_filter_args(ap)
    args = ap.parse_args()

    # Optional config (same YAML as run_pandamap.py); CLI flags override its values.
    cfg = None
    if args.config:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        from run_pandamap import load_config
        cfg = load_config(args.config)

    in_dir = args.in_dir or (cfg.output_dir if cfg else None)
    if in_dir is None:
        ap.error("--in-dir is required (or pass --config with an 'output_dir').")
    benchmark_dir = (args.benchmark_dir or (cfg.benchmark_dir if cfg else None)
                     or Path("Data/PoseBuster Benchmark Set"))
    best_equibind_only = (args.best_equibind_only if args.best_equibind_only is not None
                          else (cfg.best_equibind_only if cfg else False))
    equibind_variant = (args.equibind_variant
                        or (getattr(cfg, "equibind_variant", None) if cfg else None))
    best_diffdock_only = (args.best_diffdock_only if args.best_diffdock_only is not None
                          else (getattr(cfg, "best_diffdock_only", False) if cfg else False))
    diffdock_variant = (args.diffdock_variant
                        or (getattr(cfg, "diffdock_variant", None) if cfg else None))
    oracle_summary = (args.oracle_summary or (cfg.oracle_summary if cfg else None)
                      or DEFAULT_ORACLE_SUMMARY)
    ids_file = args.ids_file or (cfg.ids_file if cfg else None)
    # per-pose RMSD table for fig 14: CLI, else the oracle summary's sibling, else default.
    if args.per_pose_metrics:
        per_pose_metrics = args.per_pose_metrics
    else:
        sib = Path(oracle_summary).with_name("per_pose_metrics.csv") if oracle_summary else None
        per_pose_metrics = sib if (sib and sib.exists()) else DEFAULT_PER_POSE_METRICS

    out_dir = args.out_dir or (in_dir / "report")
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = read_pandamap_csv(in_dir / "pandamap_pose_summary.csv", low_memory=False)
    inter = read_pandamap_csv(in_dir / "pandamap_interactions.csv", low_memory=False)
    crystal_path = in_dir / "crystal_interactions.csv"
    crystal = (read_pandamap_csv(crystal_path, low_memory=False)
               if crystal_path.exists() else pd.DataFrame())
    # Nearest-copy crystal reference (opt-in, see CrystalCopySelector). Needs BOTH a
    # per-copy crystal file (copy_index) and a nearest-convention metrics table
    # (nearest_copy_index). A per-copy crystal file read against an instance-convention
    # table is reduced to its reference copy, which is exactly the canonical file.
    copy_sel = None
    if not crystal.empty and "copy_index" in crystal.columns:
        copy_sel = CrystalCopySelector.build(crystal, per_pose_metrics)
        if copy_sel is None:
            _ck = crystal["copy_index"].fillna(0).astype(int)
            crystal = crystal[_ck == 0].drop(columns=["copy_index"]).reset_index(drop=True)
            print(f"crystal_interactions.csv carries copy_index but {per_pose_metrics} has no "
                  "nearest_copy_index — using the reference copy only (single-instance convention).")
        else:
            print(f"crystal_interactions.csv carries copy_index and {per_pose_metrics} carries "
                  "nearest_copy_index — every pose is scored against its nearest deposited copy; "
                  "top-k unions against the union of the selected copies (nearest-copy convention).")

    # Drop whole methods from every chart/CSV. The fingerprints stay on disk; this only
    # controls what is presented, so a tool can be generated once and then left out of a
    # particular report (e.g. an engine a supervisor asked to keep out of the write-up)
    # without discarding its data or re-running PandaMap.
    # Three frames share one exclusion. Only the summary reports it and writes the
    # sidecar; the other two would repeat the same block for no new information.
    # crystal_interactions.csv is optional and the interaction frames are empty for
    # a run with no fingerprints, so each is filtered only when it can be.
    summary = mf.apply_method_filter(summary, "method", args,
                                     label="pandamap", out_dir=out_dir)
    _pats = mf.resolve_patterns(args)

    def _also_filter(frame, name):
        """crystal_interactions.csv is optional and both frames are empty when a
        run produced no fingerprints, so filter each only when it can be."""
        if frame.empty or "method" not in frame.columns:
            return frame
        return mf.apply_patterns(frame, "method", _pats, label=f"pandamap/{name}")

    inter = _also_filter(inter, "inter")
    crystal = _also_filter(crystal, "crystal")

    # Restrict to the official benchmark-set ids (complex id == 'protein' column,
    # which equals '<PDBID>_<LIG>' for the benchmark staging).
    if ids_file:
        allowed = _load_allowed_ids(ids_file)
        def _restrict(d):
            return (d[d["protein"].astype(str).isin(allowed)].copy()
                    if not d.empty and "protein" in d.columns else d)
        def _npairs(d):
            return d.drop_duplicates(["protein", "ligand"]).shape[0] if not d.empty else 0
        p0 = _npairs(summary)
        summary, inter, crystal = _restrict(summary), _restrict(inter), _restrict(crystal)
        print(f"Restricted to {len(allowed)} ids from {Path(ids_file).name}: "
              f"{p0} → {_npairs(summary)} pairs")

    # Optionally restrict every chart/CSV to one EquiBind variant — a pinned
    # variant (--equibind-variant) takes precedence over the oracle-best pick.
    if equibind_variant:
        summary, inter, best_eq = select_equibind_variant(summary, inter, equibind_variant)
        if best_eq:
            print(f"equibind-variant={equibind_variant}: keeping only '{best_eq}' "
                  f"(shown as '{pretty_method(best_eq)}').")
    elif best_equibind_only:
        summary, inter, best_eq = select_best_equibind(summary, inter, oracle_summary)
        if best_eq:
            _LABEL_OVERRIDES[best_eq] = "EquiBind*"
            print(f"best-equibind-only: '{best_eq}' is the top EquiBind variant by "
                  "PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'EquiBind*').")

    # Optionally restrict every chart/CSV to one DiffDock variant — a pinned
    # optimizer (--diffdock-variant) takes precedence over the oracle-best pick.
    if diffdock_variant:
        summary, inter, best_dd = select_best_diffdock(summary, inter, oracle_summary,
                                                       pin=diffdock_variant)
        if best_dd:
            print(f"diffdock-variant={diffdock_variant}: keeping only '{best_dd}' "
                  f"(shown as '{pretty_method(best_dd)}').")
    elif best_diffdock_only:
        summary, inter, best_dd = select_best_diffdock(summary, inter, oracle_summary)
        if best_dd:
            _LABEL_OVERRIDES[best_dd] = "DiffDock*"
            print(f"best-diffdock-only: '{best_dd}' is the top DiffDock variant by "
                  "PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'DiffDock*').")

    # Pose-selection scope: run_pandamap keeps only each method's top-N ranked poses
    # per complex (poses_per_combo / select_top_n). Determine N (config, else CLI,
    # else the data's own max) so every figure can carry a footer stating the analysis
    # is built from those top-N poses, and the depth-pooled figures never advertise
    # more poses than the run kept.
    poses_per_combo = (args.poses_per_combo
                       or (getattr(cfg, "poses_per_combo", None) if cfg else None))
    if not poses_per_combo and not summary.empty:
        _pc = summary.groupby(["method", "protein", "ligand"])["pose_name"].nunique()
        poses_per_combo = int(_pc.max()) if len(_pc) else None
    # Re-enforce the cap on the loaded data: the summary CSV is written in resume/append
    # mode, so a reused results dir can accumulate MORE than N poses per method × complex.
    # No-ops on a freshly overwritten, already-capped set.
    summary, inter = cap_top_n_per_combo(summary, inter, poses_per_combo)
    set_pose_scope_caption(poses_per_combo)
    # Depth for the "top-5" pooled / by-rank figures, capped at the actual per-complex
    # top-N (min(5, N)): a top-3 run shows "top-3" instead of an impossible "top-5",
    # while a run that kept ≥5 poses is unchanged.
    depth_pool = min(5, poses_per_combo) if poses_per_combo else 5

    order = ordered_methods(summary["method"].unique())
    for _m in order:                 # seed the shared tool-colour cache in canonical
        TOOL_COLORS.get(_m)          # order so every EquiBind variant's green is stable
    print(f"Loaded {len(summary)} poses, {len(inter)} interaction rows, "
          f"{len(crystal)} crystal rows. Methods: {order}")
    if poses_per_combo:
        print(f"Pose-selection scope: top {poses_per_combo} ranked poses per method × complex"
              f"  (depth-pooled / by-rank figures use top {depth_pool}).")

    # A. cross-method
    plot_profile(summary, order, out_dir / "01_interaction_profile.png")
    plot_type_heatmap_box(summary, order,
                          out_dir / "02_type_heatmap.png", out_dir / "02b_total_boxplot.png")
    plot_residue_hotspots(inter, summary, order, out_dir / "03_residue_hotspots.png", args.top_residues)
    if not crystal.empty:
        if copy_sel is not None:
            # the selector was built on the unfiltered crystal frame; rebuild it on the
            # method-filtered / id-restricted one so its pair set matches the report's
            copy_sel = CrystalCopySelector(crystal, copy_sel.nearest, copy_sel.ref_index)
        rec = native_recovery(inter, crystal, order, out_dir, copy_sel=copy_sel)
        if not rec.empty:
            print("\nNative-interaction recovery (mean per method):")
            print(rec.round(3).to_string())
        # top-1 pose (base name, kept) and top-5 union bucket
        plot_fingerprint_similarity(inter, crystal, order,
                                    out_dir / "05_fingerprint_similarity.png", depth=1,
                                    copy_sel=copy_sel)
        plot_fingerprint_similarity(inter, crystal, order,
                                    out_dir / "05_fingerprint_similarity_top5.png",
                                    depth=depth_pool, show_caption=False, copy_sel=copy_sel)

        # Deeper crystal-vs-pose comparison (figs 09-15). One residue-level detail
        # pass feeds the type-resolved, decomposition, gap and oracle figures.
        per_pose, per_type = build_recovery_detail(inter, crystal, copy_sel=copy_sel)
        if not per_pose.empty:
            per_pose.round(4).to_csv(out_dir / "recovery_detail_per_pose.csv", index=False)
        plot_native_recovery_by_rank(per_pose, order,
                                     out_dir / "04c_native_recovery_by_rank.png",
                                     top_n=depth_pool)
        # top-1 pose and top-5 (base name, kept) buckets
        plot_type_resolved_recovery(per_type, order,
                                    out_dir / "09_type_resolved_recovery_top1.png", depth=1)
        plot_type_resolved_recovery(per_type, order,
                                    out_dir / "09_type_resolved_recovery.png", depth=depth_pool)
        # top-1 pose and top-5 (base name, kept) buckets
        plot_contact_decomposition(per_pose, order,
                                   out_dir / "10_contact_decomposition_top1.png", depth=1)
        plot_contact_decomposition(per_pose, order,
                                   out_dir / "10_contact_decomposition.png", depth=depth_pool)
        plot_contact_decomposition_by_rank(per_pose, order,
                                           out_dir / "10b_contact_decomposition_by_rank.png",
                                           top_n=depth_pool)
        plot_residue_confusion(inter, crystal, order,
                               out_dir / "11_residue_confusion.png", args.top_residues,
                               copy_sel=copy_sel)
        plot_typed_vs_loose(per_pose, order, out_dir / "12_typed_vs_loose_gap.png")
        plot_native_f1_oracle(per_pose, order, out_dir / "13_native_f1_oracle.png")
        plot_rmsd_vs_recovery(per_pose, per_pose_metrics, order,
                              out_dir / "14_rmsd_vs_recovery.png")
        plot_pair_overlays(inter, crystal, per_pose, order,
                           out_dir / "15_pair_overlays.png", copy_sel=copy_sel)
        if copy_sel is not None:
            msg = copy_sel.report()
            print("  " + msg)
            (out_dir / "crystal_copy_selection.txt").write_text(msg + "\n")
    else:
        print("  [skip] native recovery / similarity + figs 09-15 — no crystal_interactions.csv")

    # B. chemical level
    pairs = sorted({(p, l) for p, l in zip(summary["protein"], summary["ligand"])})
    bdir = benchmark_dir if benchmark_dir and benchmark_dir.exists() else None
    desc = build_descriptor_table(pairs, bdir, summary)
    if not desc.empty:
        desc.to_csv(out_dir / "ligand_descriptors.csv", index=False)
    plot_chem_chemotype(summary, desc, out_dir)
    plot_residue_class(inter, out_dir / "07_residue_class.png")
    plot_ligand_element(inter, out_dir / "08_ligand_element.png")

    # Recoverable numeric sidecar for every test annotated on the figures above.
    _write_stats_json(out_dir)
    # One "<stem>.txt" per figure: its one-line title, moved-off subtitle, footer
    # and the statistical results (everything no longer drawn on the image).
    _flush_fig_txt(out_dir)

    print(f"\nReport written to: {out_dir.resolve()}")
    for f in sorted(out_dir.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
