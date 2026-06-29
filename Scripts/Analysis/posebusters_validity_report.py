"""Visualize PoseBusters validity per docking tool & receptor-ligand pair.

Reads posebusters_filtered_results.csv and produces:
    1. Bar chart: total / valid pose counts per docking method.
    2. Heatmap: valid-pose counts per receptor-ligand x method.
    3. Stacked bar: top-N receptor-ligand pairs comparing valid poses per method.
    4. Boxplot: distribution of per-pair valid poses per method.
    5. Heatmap: per-check pass rate per method.
    6. (EquiBind only) Bar chart: valid fraction per EquiBind variant.
    + CSV summaries written next to the plots.

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
  re-search    (refine_variant)    smina (__refSMINA) | raw (__refRAW) | None
  clamp axis   (clamp_variant)     clampON (__clampON) | clampOFF (__clampOFF) | None

So e.g. ``fpocket_raw_p001_pose01__clampOFF__refSMINA.sdf`` becomes the method
``equibind_fpocket_smina_clampOFF``. fpocket and p2rank are reported as distinct
methods (not merged into one "guided" bucket). Axes that aren't present (the
default single-variant runs) simply collapse to ``equibind_fpocket`` /
``equibind_p2rank`` / ``equibind_unguided``.
NOTE: refine_mode='on' rewrites the pose in place with no suffix/tag, so an
unsuffixed pose is reported under its base label (the smina split is only
visible for refine_mode='both' runs, which emit __refRAW/__refSMINA).

Quick usage:
    # Run with defaults
    python Scripts/Analysis/posebusters_validity_report.py

    # Run with custom input/output locations
    python Scripts/Analysis/posebusters_validity_report.py \
        --csv posebusters_results/benchmark/dock/posebusters_filtered_results.csv \
        --out-dir posebusters_results/benchmark/dock/validity_report

    # Keep EquiBind as a single "equibind_guided" bucket (old behaviour)
    python Scripts/Analysis/posebusters_validity_report.py --no-split-equibind

    # Show only the single best EquiBind variant (highest oracle_rmsd_le_2.0A_%,
    # read from --oracle-summary), relabelled "EquiBind*". For crystal-free sets
    # (Orai) the default --oracle-summary borrows the benchmark ranking.
    python Scripts/Analysis/posebusters_validity_report.py --best-equibind-only
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
_REFINE_ORDER = {None: 0, "raw": 1, "smina": 2}
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
        parts.append("smina-opt" if refine == "smina" else "raw")
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


def select_best_equibind(df: pd.DataFrame, oracle_csv: Path) -> tuple[pd.DataFrame, str | None]:
    """Keep non-EquiBind rows + only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest oracle RMSD ≤ 2 Å success rate
    (``oracle_rmsd_le_2.0A_%``) in *oracle_csv* — the oracle_summary.csv written by
    posebusters_pose_comparison.py. For datasets without their own crystal (Orai),
    point --oracle-summary at the benchmark summary to borrow its ranking. The best
    variant is chosen among the variants actually present in *df*; AutoDock/DiffDock
    rows are always kept. Returns (filtered_df, best_variant_key); best is None and
    df unchanged when no EquiBind variant is present, the summary is missing/unreadable,
    or no present variant has an oracle score.
    """
    methods = df["docking_method"].astype(str)
    eq_mask = methods.str.startswith("equibind")
    if not eq_mask.any():
        return df, None
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-equibind-only] oracle summary not found at {oracle_csv} — "
              "keeping all EquiBind variants.")
        return df, None

    col = "oracle_rmsd_le_2.0A_%"
    osum = pd.read_csv(oracle_csv, index_col=0)
    if col not in osum.columns:
        print(f"  [best-equibind-only] '{col}' missing from {oracle_csv} — "
              "keeping all EquiBind variants.")
        return df, None

    scores = pd.to_numeric(osum[col], errors="coerce")
    eq_present = sorted(methods[eq_mask].unique())
    cand = {m: float(scores[m]) for m in eq_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-equibind-only] none of the present EquiBind variants have an "
              f"oracle score in {oracle_csv} — keeping all EquiBind variants.")
        return df, None

    best = max(cand, key=cand.get)
    keep = (~eq_mask) | (methods == best)
    return df[keep].reset_index(drop=True), best


def per_tool_summary(df: pd.DataFrame, order: list[str]) -> pd.DataFrame:
    grp = (
        df.groupby("docking_method")
        .agg(total_poses=("pb_valid", "size"),
             valid_poses=("pb_valid", "sum"))
        .reindex(order)
        .dropna(how="all")
    )
    grp["valid_fraction"] = grp["valid_poses"] / grp["total_poses"]
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


# ───────────────────────────── plots ─────────────────────────────

def plot_per_tool(summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(max(8, 1.4 * len(summary)), 5))
    x = np.arange(len(summary))
    w = 0.38
    ax.bar(x - w / 2, summary["total_poses"], w, label="Total poses",
           color="#bdbdbd", edgecolor="black")
    ax.bar(x + w / 2, summary["valid_poses"], w, label="PB-Valid poses",
           color="#2ca02c", edgecolor="black")

    for i, (tot, val) in enumerate(zip(summary["total_poses"], summary["valid_poses"])):
        pct = (val / tot * 100) if tot else 0
        ax.text(i + w / 2, val, f"{int(val)}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=9)
        ax.text(i - w / 2, tot, f"{int(tot)}", ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels([_pretty_method(t) for t in summary.index],
                       rotation=20, ha="right")
    ax.set_ylabel("Number of poses")
    ax.set_title("PoseBusters Benchmark — Pose validity per docking method\n"
                 "(valid = passes all canonical PB checks)")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_heatmap(valid_mat: pd.DataFrame, out: Path, top_n: int = 60) -> None:
    mat = valid_mat.copy()
    mat["__sum"] = mat.sum(axis=1)
    mat = mat.sort_values("__sum", ascending=False).head(top_n).drop(columns="__sum")
    mat.columns = [_pretty_method(c) for c in mat.columns]

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(mat.columns)), max(8, 0.22 * len(mat))))
    sns.heatmap(mat, annot=True, fmt=".0f", cmap="YlGnBu", cbar_kws={"label": "Valid poses"},
                linewidths=0.3, linecolor="white", ax=ax)
    ax.set_xlabel("Docking method")
    ax.set_ylabel("Receptor / Ligand")
    ax.set_title(f"Valid poses per receptor-ligand pair (top {len(mat)})")
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
    ax.set_title(f"Valid poses per docking method — top {len(pairs)} receptor-ligand pairs")
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
    ax.boxplot([valid[c].values for c in cols], showmeans=True)
    ax.set_xticks(range(1, len(cols) + 1))
    ax.set_xticklabels([_pretty_method(c) for c in cols], rotation=20, ha="right")
    ax.set_ylabel("Valid poses per receptor-ligand pair")
    ax.set_title("Distribution of valid poses across receptor-ligand pairs")
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
    pr = pr[[c for c in order if c in pr.columns]]
    pr.columns = [_pretty_method(c) for c in pr.columns]

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(pr.columns)), 6))
    sns.heatmap(pr * 100, annot=True, fmt=".1f", cmap="RdYlGn", vmin=0, vmax=100,
                cbar_kws={"label": "Pass rate (%)"}, ax=ax)
    ax.set_title("Per-check pass rate (%) per docking method")
    ax.set_xlabel("")
    ax.set_ylabel("PoseBusters critical check")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    plt.close(fig)


def plot_equibind_variants(summary: pd.DataFrame, out: Path, colors: dict) -> bool:
    """Valid fraction (%) per EquiBind variant. Returns False if no variants."""
    eq = summary[summary.index.str.startswith("equibind")]
    if eq.empty:
        return False

    labels = [_pretty_method(m).replace("EquiBind ", "").strip("()") or "EquiBind"
              for m in eq.index]
    x = np.arange(len(eq))
    fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(eq)), 5))
    bars = ax.bar(x, eq["valid_fraction"] * 100, 0.6,
                  color=[colors.get(m, "#2ca02c") for m in eq.index],
                  edgecolor="black")
    for rect, (tot, val) in zip(bars, zip(eq["total_poses"], eq["valid_poses"])):
        pct = (val / tot * 100) if tot else 0
        ax.text(rect.get_x() + rect.get_width() / 2, pct,
                f"{pct:.1f}%\n{int(val)}/{int(tot)}",
                ha="center", va="bottom", fontsize=9)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("PB-Valid poses (%)")
    ax.set_ylim(0, max(5, eq["valid_fraction"].max() * 100 * 1.25))
    ax.set_title("EquiBind — PB-validity per variant\n"
                 "(fpocket / p2rank / unguided, smina re-search, centroid clamp)")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
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
                    help="Split EquiBind into fpocket/p2rank/unguided x smina/raw "
                         "x clamp variants from the CSV provenance columns "
                         "(filename fallback) (default: on).")
    ap.add_argument("--best-equibind-only", action="store_true",
                    help="Keep only the single best-performing EquiBind variant "
                         "(highest oracle_rmsd_le_2.0A_%%, read from --oracle-summary) "
                         "in all plots/CSVs, relabelled 'EquiBind*'. AutoDock/DiffDock "
                         "are unaffected. For datasets without a crystal (Orai), the "
                         "default --oracle-summary borrows the benchmark ranking.")
    ap.add_argument("--oracle-summary", type=Path, default=DEFAULT_ORACLE_SUMMARY,
                    help="oracle_summary.csv from posebusters_pose_comparison.py used "
                         "to pick the best EquiBind variant for --best-equibind-only "
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

    # Optionally restrict every plot/CSV to the single best EquiBind variant.
    if args.best_equibind_only:
        df, best_eq = select_best_equibind(df, args.oracle_summary)
        if best_eq:
            _LABEL_OVERRIDES[best_eq] = "EquiBind*"
            print(f"best-equibind-only: '{best_eq}' is the top EquiBind variant by "
                  f"oracle_rmsd_le_2.0A_% — keeping only it (shown as 'EquiBind*').")

    order = _method_order(df)
    colors = _method_colors(order)

    summary = per_tool_summary(df, order)
    valid_mat = per_pair_tool_matrix(df, "valid", order)
    total_mat = per_pair_tool_matrix(df, "total", order)
    frac_mat = per_pair_tool_matrix(df, "fraction", order)

    summary.to_csv(args.out_dir / "summary_per_tool.csv")
    valid_mat.to_csv(args.out_dir / "valid_poses_per_pair_tool.csv")
    total_mat.to_csv(args.out_dir / "total_poses_per_pair_tool.csv")
    frac_mat.round(3).to_csv(args.out_dir / "valid_fraction_per_pair_tool.csv")

    axis = equibind_axis_breakdown(df)
    if not axis.empty:
        axis.round(4).to_csv(args.out_dir / "equibind_axis_breakdown.csv")

    plot_per_tool(summary, args.out_dir / "01_per_tool_validity.png")
    plot_heatmap(valid_mat, args.out_dir / "02_valid_heatmap_top.png",
                 top_n=args.top_n_heatmap)
    plot_grouped_bars(valid_mat, args.out_dir / "03_valid_grouped_bars_top.png",
                      colors, top_n=args.top_n_bars)
    plot_validity_distribution(df, args.out_dir / "04_valid_per_pair_distribution.png",
                               order)
    plot_check_passrate(df, args.out_dir / "05_per_check_passrate.png", order)
    has_eq_plot = plot_equibind_variants(
        summary, args.out_dir / "06_equibind_variant_validity.png", colors)

    n_figs = 6 if has_eq_plot else 5
    print(f"Total poses scored        : {len(df):,}")
    print(f"Total PB-valid poses      : {int(df['pb_valid'].sum()):,} "
          f"({df['pb_valid'].mean() * 100:.2f}%)")
    print(f"Receptor-ligand pairs     : {df['pair'].nunique()}")
    print(f"Methods / variants found  : {len(order)} "
          f"({', '.join(order)})")
    print()
    print("── Per-method summary ───────────────────────────")
    print(summary.to_string(float_format=lambda x: f"{x:.3f}"))
    if not axis.empty:
        print()
        print("── EquiBind axis breakdown (pocket × refine × clamp) ──")
        print(axis.to_string(float_format=lambda x: f"{x:.3f}"))
    print()
    print(f"Wrote CSVs and {n_figs} figures to: {args.out_dir.resolve()}")


if __name__ == "__main__":
    main()
