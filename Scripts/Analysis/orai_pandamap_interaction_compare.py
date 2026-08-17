#!/usr/bin/env python
"""Compare PandaMap protein–ligand interactions on Orai: Experimental (JKU) vs Benchmark ligands.

Companion to ``pandamap_interaction_report.py`` (which runs once per dataset and
writes a single-dataset ``report/``, kept untouched by this script) and a sibling of
``orai_cross_tool_agreement_compare.py`` / ``orai_pbvalid_tm_share_compare.py``. Here we
read the cached PandaMap fingerprints each dataset already produced
(``pandamap_pose_summary.csv`` + ``pandamap_interactions.csv``) and overlay the two Orai
ligand categories:

    Experimental ligands  →  "Exp. Ligands"       (Orai × JKU:  pandamap_results/orai_jku)
    Benchmark ligands     →  "Benchmark ligands"  (Orai × Bench: pandamap_results/orai_benchmark)

**Scope — PB-valid poses outside the transmembrane.** Both PandaMap runs are built from
``posebusters_filtered_results.no_tm.csv`` (kept = placed AND outside the TM1 pore, see
orai_transmembrane_exclusion.py) with ``pb_valid_only: true`` (config), so every fingerprinted
pose is already PoseBusters-valid AND non-transmembrane. This script additionally re-asserts
that scope defensively: it filters on the ``pb_valid`` flag where present and re-caps to the
top-N poses per (tool, complex), so a results dir that accumulated stale rows in resume/append
mode cannot silently widen the population.

**Same receptor both sides.** Both ligand sets dock into the SAME four Orai1 MD frames, so the
Orai residue axis (Glu106 selectivity filter, Val102 gate, Arg91 basic funnel …) is directly
comparable across the two sets — the residue hot-spot and fingerprint-overlap views exploit this.

**Toolchains.** One variant per tool, matched across the two datasets (both now use the same
selection): AutoDock Vina, DiffDock (smina-optimised), EquiBind (gnina-optimised). EquiBind is
kept but flagged: it has very few PB-valid non-TM poses on either side (it dumps ligands in the
pore), so its Experimental distribution rests on a handful of poses — read it as exploratory.

Four comparison views (each a figure; inferential tests go to the companion ``*_stats.txt``,
the panels carry descriptive content only — project convention):

  1. fig_total_interactions_compare.png   per-tool total contacts per pose, Exp vs Benchmark.
  2. fig_type_profile_compare.png         per-tool mean count per pose of each interaction type.
  3. fig_residue_hotspots_compare.png     per-tool Orai residue hot-spots (contact frequency),
                                          Exp vs Benchmark, on a shared residue axis.
  4. fig_fingerprint_overlap_compare.png  per-tool overlap of the two sets' characteristic Orai
                                          contact fingerprints (shared / exp-only / bench-only +
                                          Jaccard), residue-level and typed.

Plus CSVs (pose totals, type profile, residue hot-spots, fingerprint overlap) and a stats txt.

Run under the vina env (scipy + matplotlib; pandamap not needed here):
    /home/manndo/anaconda3/envs/vina/bin/python \
        Scripts/Analysis/orai_pandamap_interaction_compare.py
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as _pe
from matplotlib.patches import Patch

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

from pose_topn import top_n_keys                       # noqa: E402  top-N pose cap

# ── shared constants (reuse the single-dataset report's, fall back to local copies
#    so this script still runs if that heavy module can't be imported) ────────────
try:
    from pandamap_interaction_report import (
        INTERACTION_TYPES, INTERACTION_LABELS, pretty_itype,
    )
except Exception:                                   # pragma: no cover
    INTERACTION_TYPES = [
        "hydrogen_bonds", "carbon_pi", "pi_pi_stacking", "donor_pi", "amide_pi",
        "hydrophobic", "ionic", "halogen_bonds", "cation_pi", "metal_coordination",
        "salt_bridge", "covalent", "alkyl_pi", "attractive_charge", "pi_cation",
        "repulsion",
    ]
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

try:                                                # shared, unit-tested stats helpers
    import stats_utils as su
    _HAVE_SU = True
except Exception:                                   # pragma: no cover
    su = None
    _HAVE_SU = False


# ── canonical tools + palette ────────────────────────────────────────────────────
# Palette matches 09f_pbvalid_yield_dominant_orai_benchmark_vs_experimental.png (and the
# orai_pbvalid_tm_share_compare whisker): box/bar FILL encodes the docking TOOL — AutoDock
# blue, DiffDock orange, EquiBind green — and the two ligand SETS differ ONLY by fill STYLE
# (Benchmark solid, Experimental hatched). Benchmark is drawn first (left / solid).
CANONICAL_TOOLS = ["autodock", "diffdock", "equibind"]
TOOL_LABEL = {"autodock": "AutoDock Vina", "diffdock": "DiffDock", "equibind": "EquiBind"}
TOOL_COLORS = {"autodock": "#1f77b4", "diffdock": "#ff7f0e", "equibind": "#2ca02c"}

# Per-dataset fill STYLE (the fill COLOUR comes from the tool, never the dataset).
DATASET_STYLE = {
    "bench": {"label": "Orai × Benchmark ligands",    "short": "Benchmark",
              "alpha": 0.85, "hatch": ""},
    "exp":   {"label": "Orai × Experimental ligands", "short": "Experimental",
              "alpha": 0.45, "hatch": "///"},
}
CAT_ORDER = ["bench", "exp"]          # Benchmark first (left / solid), matching 09f


# ── method → canonical tool + variant selection ──────────────────────────────────

def _tool_of(method: str) -> Optional[str]:
    m = str(method)
    for t in CANONICAL_TOOLS:
        if m.startswith(t):
            return t
    return None


def _variant_tokens(method: str) -> List[str]:
    """Tokens after the base tool name, e.g. 'diffdock_smina' → ['smina'],
    'equibind_unguided_gnina' → ['unguided','gnina']."""
    return str(method).split("_")[1:]


def _pick_method(methods_present: List[str], tool: str, want: Optional[str]) -> Optional[str]:
    """Choose the single method key for *tool* preferring the variant token *want*.

    AutoDock has no variant. For DiffDock/EquiBind: keep the candidate whose variant
    tokens include *want* (e.g. 'smina' / 'gnina'); if none match but exactly one
    candidate exists, use it (the datasets each carry one variant per tool); otherwise
    fall back to the lexicographically-first candidate."""
    cands = sorted(m for m in methods_present if _tool_of(m) == tool)
    if not cands:
        return None
    if tool == "autodock" or not want:
        return cands[0]
    matched = [m for m in cands if want in _variant_tokens(m)]
    if matched:
        return matched[0]
    # Wanted variant absent: prefer a non-raw variant (one with an optimiser/refine
    # token) over the bare/raw base tool, and warn — never silently fall back to raw.
    non_raw = [m for m in cands if _variant_tokens(m) and "raw" not in _variant_tokens(m)]
    pick = (non_raw or cands)[0]
    print(f"  [variant] {tool}: wanted '{want}' not present among {cands}; using '{pick}'.")
    return pick


def _coerce_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0"))


def _cap_top_n(summary: pd.DataFrame, inter: pd.DataFrame, n: int
               ) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """Defensively keep only each method's top-*n* poses per (protein, ligand) by pose_rank.

    Mirrors pandamap_interaction_report.cap_top_n_per_combo: PandaMap writes its CSVs in
    resume/append mode, so a reused dir can hold more than N poses per (tool, complex).
    No-op when nothing exceeds N. Returns (summary, inter, n_dropped)."""
    key = ["method", "protein", "ligand", "pose_name"]
    if not n or summary.empty or "pose_rank" not in summary.columns:
        return summary, inter, 0
    s = summary.copy()
    s["_rank"] = pd.to_numeric(s["pose_rank"], errors="coerce").fillna(999)
    keep_keys = (s.drop_duplicates(key)
                 .sort_values(["_rank", "pose_name"])
                 .groupby(["method", "protein", "ligand"], sort=False)
                 .head(n)[key].drop_duplicates())
    summary2 = (summary.merge(keep_keys, on=key, how="inner")
                .drop_duplicates(key).reset_index(drop=True))
    inter2 = (inter.merge(keep_keys, on=key, how="inner")
              if not inter.empty and set(key).issubset(inter.columns) else inter)
    return summary2, inter2, len(summary) - len(summary2)


def _apply_topn_gnina(summary: pd.DataFrame, inter: pd.DataFrame, keys: set
                      ) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """Keep only poses whose (protein, ligand, pose-file-basename) is in *keys* — the
    gnina-energy top-N allowlist. Matches on pose identity (path-independent) because the
    PandaMap summary stores EquiBind's pose_rank as a 999 placeholder, so a rank threshold
    can't select its top-N; the allowlist carries the real per-tool rank."""
    def _k(df):
        return list(zip(df["protein"].astype(str), df["ligand"].astype(str),
                        df["pose_file"].astype(str).map(lambda p: Path(str(p)).name)))
    before = len(summary)
    summary = summary[[k in keys for k in _k(summary)]].reset_index(drop=True)
    if not inter.empty and {"protein", "ligand", "pose_file"}.issubset(inter.columns):
        inter = inter[[k in keys for k in _k(inter)]].reset_index(drop=True)
    return summary, inter, before - len(summary)


def load_dataset(key: str, root: Path, diffdock_variant: str, equibind_variant: str,
                 top_n: int, posebusters_csv: Optional[Path] = None,
                 top_n_poses: int = 0) -> dict:
    """Load one PandaMap dataset, select one variant per tool, re-assert the
    PB-valid + top-N scope, and tag each pose/interaction with its canonical tool.

    ``top_n_poses`` > 0 additionally restricts every tool to its top-N ranked poses per
    (protein, ligand) using the NATIVE rank (AutoDock=Vina, DiffDock=confidence, EquiBind=
    gnina energy), read from ``posebusters_csv``. This makes the interaction figure use the
    same top-N pose budget as the yield / transmembrane / pose-cluster figures."""
    root = Path(root)
    sp = root / "pandamap_pose_summary.csv"
    ip = root / "pandamap_interactions.csv"
    if not sp.exists() or not ip.exists():
        raise FileNotFoundError(f"missing PandaMap CSVs under {root} "
                                f"({sp.name} / {ip.name}) — run run_pandamap.py first")
    summary = pd.read_csv(sp, low_memory=False)
    inter = pd.read_csv(ip, low_memory=False)

    notes: List[str] = []
    # PB-valid filter (defensive; both configs already set pb_valid_only:true).
    if "pb_valid" in summary.columns:
        before = len(summary)
        pbv = _coerce_bool(summary["pb_valid"])
        summary = summary[pbv].reset_index(drop=True)
        if "pb_valid" in inter.columns:
            inter = inter[_coerce_bool(inter["pb_valid"])].reset_index(drop=True)
        if before != len(summary):
            notes.append(f"pb_valid filter: {before} → {len(summary)} poses")
    else:
        notes.append("no pb_valid column (older run); relying on the run's "
                     "pb_valid_only:true source scope")

    methods_present = sorted(summary["method"].astype(str).unique())
    chosen: Dict[str, str] = {}
    for tool, want in (("autodock", None), ("diffdock", diffdock_variant),
                       ("equibind", equibind_variant)):
        mk = _pick_method(methods_present, tool, want)
        if mk is not None:
            chosen[tool] = mk
    keep_methods = set(chosen.values())
    summary = summary[summary["method"].astype(str).isin(keep_methods)].reset_index(drop=True)
    inter = inter[inter["method"].astype(str).isin(keep_methods)].reset_index(drop=True)

    # canonical-tool column
    summary["tool"] = summary["method"].map(_tool_of)
    inter["tool"] = inter["method"].map(_tool_of)

    summary, inter, dropped = _cap_top_n(summary, inter, top_n)
    if dropped:
        notes.append(f"top-{top_n} cap dropped {dropped} accumulated pose row(s)")

    # Native top-N gnina/confidence/Vina rank cap (equalises the pose budget with the yield /
    # TM / pose-cluster figures). PandaMap stores EquiBind pose_rank=999, so we filter by a
    # pose-identity allowlist carrying the real rank instead of thresholding pose_rank.
    if top_n_poses and top_n_poses > 0:
        if posebusters_csv and Path(posebusters_csv).exists():
            keys = top_n_keys(posebusters_csv, top_n_poses)
            summary, inter, dropped2 = _apply_topn_gnina(summary, inter, keys)
            notes.append(f"top-{top_n_poses}-poses gnina cap: dropped {dropped2} pose row(s) "
                         f"(ranked on {Path(posebusters_csv).name})")
        else:
            notes.append(f"top-{top_n_poses}-poses cap REQUESTED but posebusters_csv missing "
                         f"({posebusters_csv}) — skipped")

    # normalise residue key (resname + resnum) for the hot-spot / fingerprint views
    inter = inter[inter["resname"].notna() & inter["resnum"].notna()].copy()
    inter["resnum_i"] = pd.to_numeric(inter["resnum"], errors="coerce")
    inter = inter[inter["resnum_i"].notna()].copy()
    inter["res"] = (inter["resname"].astype(str).str.upper()
                    + inter["resnum_i"].astype(int).astype(str))

    return {
        "key": key, "label": DATASET_STYLE[key]["label"], "root": root,
        "summary": summary, "inter": inter, "chosen": chosen,
        "methods_present": methods_present, "notes": notes,
        "mtime": datetime.fromtimestamp(sp.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S"),
    }


# ── small helpers ────────────────────────────────────────────────────────────────

def _pretty_res(res: str) -> str:
    """'GLU106' → 'Glu106' (3-letter code title-cased, number kept)."""
    i = 0
    while i < len(res) and not res[i].isdigit():
        i += 1
    return res[:i].capitalize() + res[i:]


def _poses_per_tool(summary: pd.DataFrame, tool: str) -> int:
    return int(summary[summary["tool"] == tool]["pose_name"].nunique())


def _per_complex_means(summary: pd.DataFrame, tool: str, col: str) -> np.ndarray:
    """Per-(protein,ligand) mean of *col* for one tool (the pseudoreplication-free unit)."""
    s = summary[summary["tool"] == tool]
    if s.empty or col not in s.columns:
        return np.array([])
    v = s.groupby(["protein", "ligand"])[col].mean().to_numpy()
    return v[~np.isnan(v)]


def _label_panels(axes, fontsize: int = 12) -> None:
    flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
    for i, ax in enumerate([a for a in flat if a is not None]):
        ax.text(-0.06, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                fontsize=fontsize, fontweight="bold", va="bottom", ha="right")


def _set_legend_handles() -> list:
    """Neutral-grey SET legend — solid = Benchmark, hatched = Experimental. The fill
    COLOUR on the bars encodes the tool, so the legend conveys STYLE only, matching the
    09f_pbvalid_yield_dominant / orai_pbvalid_tm_share_compare convention."""
    return [Patch(facecolor="0.5", alpha=DATASET_STYLE[k]["alpha"], edgecolor="0.3",
                  hatch=DATASET_STYLE[k]["hatch"] or None, label=DATASET_STYLE[k]["label"])
            for k in CAT_ORDER]


def _legend_fontsize(fig_w: float, ref_w: float = 10.5, ref_fs: float = 9.0) -> float:
    """Legend font size that keeps the legend as visually LARGE as the
    fig_total_interactions reference (fontsize 9 on its ~10.5-inch-wide figure).
    Font points are absolute, so a wider multi-panel figure shrinks a fixed-point
    legend once it is scaled to fit; scaling the size by figure width restores parity."""
    return round(max(ref_fs, ref_fs * fig_w / ref_w), 1)


def _group_x_axis(ax, spans: dict, box_w: float, y_line: float = -0.15,
                  y_label: float = -0.175, fontsize: int = 11) -> None:
    """Draw the 09f-style second x-axis level: an underline bracket under each tool's
    boxes plus the tool name below it (the per-box tick labels carry the set + n=)."""
    trans = ax.get_xaxis_transform()
    for tool, xs in spans.items():
        lo, hi, xc = min(xs), max(xs), sum(xs) / len(xs)
        ax.plot([lo - box_w / 2, hi + box_w / 2], [y_line, y_line], transform=trans,
                color="0.45", lw=1.0, clip_on=False, zorder=1)
        ax.text(xc, y_label, TOOL_LABEL[tool], transform=trans, ha="center", va="top",
                fontsize=fontsize, fontweight="bold", clip_on=False)


# ── View 1: total interactions per pose, Exp vs Benchmark, per tool ──────────────

def fig_total_interactions(datasets: Dict[str, dict], tools: List[str], out: Path):
    fig, ax = plt.subplots(figsize=(max(8.5, 3.0 * len(tools) + 1.5), 6.6))
    rng = np.random.default_rng(0)
    box_w = 0.34
    gap = box_w + 0.05
    boxes, spans = [], {}
    for ti, tool in enumerate(tools):
        xs = []
        for di, key in enumerate(CAT_ORDER):
            ds = datasets[key]; st = DATASET_STYLE[key]
            v = pd.to_numeric(ds["summary"][ds["summary"]["tool"] == tool]["total_interactions"],
                              errors="coerce").dropna().to_numpy()
            pos = ti + (di - 0.5) * gap
            boxes.append(dict(tool=tool, pos=pos, data=v, n=int(v.size), color=TOOL_COLORS[tool],
                              alpha=st["alpha"], hatch=st["hatch"], short=st["short"]))
            xs.append(pos)
        spans[tool] = xs
    bp = ax.boxplot([b["data"] if b["n"] else np.array([np.nan]) for b in boxes],
                    positions=[b["pos"] for b in boxes], widths=box_w, patch_artist=True,
                    showmeans=True, showfliers=False, manage_ticks=False,
                    medianprops=dict(color="#222222", lw=1.6),
                    meanprops=dict(marker="D", markerfacecolor="white",
                                   markeredgecolor="#222222", markersize=5))
    for patch, b in zip(bp["boxes"], boxes):
        patch.set_facecolor(b["color"]); patch.set_alpha(b["alpha"])
        patch.set_edgecolor(b["color"]); patch.set_linewidth(1.4)
        if b["hatch"]:
            patch.set_hatch(b["hatch"])
    for b in boxes:
        vv = b["data"]
        if vv.size:
            jit = b["pos"] + (rng.random(vv.size) - 0.5) * box_w * 0.6
            ax.scatter(jit, vv, s=9, color=b["color"], edgecolor="white", lw=0.25,
                       alpha=0.5 if vv.size > 40 else 0.9, zorder=3)
    # Pose-level median value printed on each box so the descriptive medians the text
    # quotes are legible on the figure (the inferential tests aggregate to the per-complex
    # mean; that unit is reported in the stats file and the prose, not on the plot).
    for b in boxes:
        if b["n"]:
            med = float(np.median(b["data"]))
            ax.text(b["pos"], med, f"{med:g}", ha="center", va="center",
                    fontsize=8, fontweight="bold", color="#111111", zorder=6,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                              edgecolor="none", alpha=0.82))
    # two-level x-axis: per-box tick = set + n complexes-of-poses; tool group label below
    ax.set_xticks([b["pos"] for b in boxes])
    ax.set_xticklabels([f"{b['short']}\nn={b['n']}" for b in boxes], fontsize=8.6)
    ax.tick_params(axis="x", length=0, pad=6)
    _group_x_axis(ax, spans, box_w, y_line=-0.15, y_label=-0.175, fontsize=11)
    ax.set_ylabel("Total protein–ligand interactions per pose (PandaMap)")
    ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    ax.margins(x=0.03)
    ax.grid(axis="y", alpha=0.3); ax.set_axisbelow(True)
    handles = _set_legend_handles() + [
        plt.Line2D([0], [0], marker="D", color="none", markerfacecolor="white",
                   markeredgecolor="#222222", markersize=6, label="mean")]
    fig.suptitle("Orai contact count per pose — Benchmark vs Experimental ligands\n"
                 "(PoseBusters-valid poses outside the transmembrane)", fontsize=12.5, y=0.995)
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), fontsize=9,
               frameon=True, bbox_to_anchor=(0.5, 0.905))
    fig.subplots_adjust(top=0.83, bottom=0.17)
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


# ── View 2: interaction-type profile, Exp vs Benchmark, per tool ─────────────────

def _present_types(datasets: Dict[str, dict]) -> List[str]:
    present = []
    for t in INTERACTION_TYPES:
        tot = 0.0
        for ds in datasets.values():
            if t in ds["summary"].columns:
                tot += pd.to_numeric(ds["summary"][t], errors="coerce").fillna(0).sum()
        if tot > 0:
            present.append(t)
    return present


def fig_type_profile(datasets: Dict[str, dict], tools: List[str], out: Path):
    types = _present_types(datasets)
    if not types:
        print("  [warn] fig_type_profile: no interaction types present — skipping")
        return
    # order types by overall mean-per-pose (pooled across tools/datasets), most common first
    order_key = {}
    for t in types:
        vals = []
        for ds in datasets.values():
            if t in ds["summary"].columns:
                vals.append(pd.to_numeric(ds["summary"][t], errors="coerce").mean())
        order_key[t] = np.nanmean(vals) if vals else 0.0
    types = sorted(types, key=lambda t: order_key[t], reverse=True)
    y = np.arange(len(types))[::-1]

    # vertical stack: one panel per tool (rows), shared residue/type axis kept per-panel so
    # each tool's own x-scale is preserved (tools differ by an order of magnitude in count).
    n_panels = len(tools)
    per_h = max(3.2, 0.32 * len(types) + 1.2)
    reserve = 0.6                                    # inches reserved at top for title + legend
    fig_h = per_h * n_panels + reserve
    fig, axes = plt.subplots(n_panels, 1, figsize=(8.2, fig_h))
    axes = np.atleast_1d(axes)
    bh = 0.38
    for ax, tool in zip(axes, tools):
        col = TOOL_COLORS[tool]; ns = {}
        for di, key in enumerate(CAT_ORDER):
            ds = datasets[key]; st = DATASET_STYLE[key]
            s = ds["summary"][ds["summary"]["tool"] == tool]
            n = int(s["pose_name"].nunique()); ns[key] = n
            means = [pd.to_numeric(s[t], errors="coerce").mean() if (n and t in s.columns) else np.nan
                     for t in types]
            offs = (di - 0.5) * (bh + 0.02)
            ax.barh(y + offs, np.nan_to_num(means), height=bh, color=col, alpha=st["alpha"],
                    hatch=st["hatch"] or None, edgecolor=col, lw=0.7)
        ax.set_title(f"{TOOL_LABEL[tool]}\nBenchmark n={ns['bench']} · Experimental n={ns['exp']}",
                     fontsize=10)
        ax.set_xlabel("Mean count per pose")
        ax.set_yticks(y)
        ax.set_yticklabels([pretty_itype(t) for t in types], fontsize=9)
        ax.set_ylabel("Interaction type")
        ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    _label_panels(axes)
    fig.suptitle("Interaction-type profile on Orai — Benchmark vs Experimental ligands, per tool\n"
                 "(mean count per pose; PoseBusters-valid poses outside the transmembrane)",
                 fontsize=12.5, y=1 - 0.28 / fig_h)
    fig.legend(handles=_set_legend_handles(), loc="upper center", ncol=2,
               fontsize=_legend_fontsize(8.2),
               frameon=True, bbox_to_anchor=(0.5, 1 - 0.85 / fig_h))
    fig.tight_layout(rect=(0, 0, 1, 1 - reserve / fig_h))
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


# ── residue hot-spot frequency table (shared by fig 3 + stats) ───────────────────

def _residue_freq(ds: dict, tool: str) -> pd.Series:
    """fraction of a (dataset, tool)'s poses that contact each Orai residue."""
    inter = ds["inter"]; s = inter[inter["tool"] == tool]
    n = _poses_per_tool(ds["summary"], tool)
    if n == 0 or s.empty:
        return pd.Series(dtype=float)
    contact = s.groupby("res")["pose_name"].nunique()
    return (contact / n).sort_values(ascending=False)


def _top_shared_residues(datasets: Dict[str, dict], tools: List[str], top_n: int) -> List[str]:
    """Top-N Orai residues by pooled contact frequency across all tools + both datasets."""
    score: Dict[str, float] = {}
    for ds in datasets.values():
        for tool in tools:
            fr = _residue_freq(ds, tool)
            for res, f in fr.items():
                score[res] = score.get(res, 0.0) + float(f)
    ordered = sorted(score, key=lambda r: score[r], reverse=True)
    return ordered[:top_n]


# ── View 3: Orai residue hot-spots, Exp vs Benchmark, per tool ───────────────────

def fig_residue_hotspots(datasets: Dict[str, dict], tools: List[str], out: Path,
                         top_n: int = 15):
    residues = _top_shared_residues(datasets, tools, top_n)
    if not residues:
        print("  [warn] fig_residue_hotspots: no contacted residues — skipping")
        return
    y = np.arange(len(residues))[::-1]
    freqs = {key: {tool: _residue_freq(datasets[key], tool) for tool in tools}
             for key in CAT_ORDER}
    # vertical stack: one panel per tool (rows), shared 0–100% residue-contact x-axis.
    n_panels = len(tools)
    per_h = max(3.2, 0.34 * len(residues) + 1.2)
    reserve = 0.6                                    # inches reserved at top for title + legend
    fig_h = per_h * n_panels + reserve
    fig, axes = plt.subplots(n_panels, 1, figsize=(8.2, fig_h))
    axes = np.atleast_1d(axes)
    bh = 0.38
    for ax, tool in zip(axes, tools):
        col = TOOL_COLORS[tool]; ns = {}
        for di, key in enumerate(CAT_ORDER):
            st = DATASET_STYLE[key]
            fr = freqs[key][tool]
            n = _poses_per_tool(datasets[key]["summary"], tool); ns[key] = n
            vals = [100 * float(fr.get(r, 0.0)) for r in residues]
            offs = (di - 0.5) * (bh + 0.02)
            ax.barh(y + offs, vals, height=bh, color=col, alpha=st["alpha"],
                    hatch=st["hatch"] or None, edgecolor=col, lw=0.7)
        ax.set_title(f"{TOOL_LABEL[tool]}\nBenchmark n={ns['bench']} · Experimental n={ns['exp']}",
                     fontsize=10)
        ax.set_xlabel("Poses contacting residue (%)")
        ax.set_xlim(0, 100)
        ax.set_yticks(y)
        ax.set_yticklabels([_pretty_res(r) for r in residues], fontsize=9)
        ax.set_ylabel("Orai1 residue")
        ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    _label_panels(axes)
    fig.suptitle("Orai residue hot-spots — Benchmark vs Experimental ligands, per tool\n"
                 "(fraction of a tool's poses contacting each residue; same receptor both sets)",
                 fontsize=12.5, y=1 - 0.28 / fig_h)
    fig.legend(handles=_set_legend_handles(), loc="upper center", ncol=2,
               fontsize=_legend_fontsize(8.2),
               frameon=True, bbox_to_anchor=(0.5, 1 - 0.85 / fig_h))
    fig.tight_layout(rect=(0, 0, 1, 1 - reserve / fig_h))
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


# ── fingerprint overlap (shared by fig 4 + stats) ────────────────────────────────

def _char_residue_set(ds: dict, tool: str, min_freq: float) -> set:
    """Residues contacted in at least *min_freq* fraction of a (dataset, tool)'s poses —
    a 'characteristic contact' set robust to the large pose-count imbalance (Bench ≫ Exp)."""
    fr = _residue_freq(ds, tool)
    return set(fr[fr >= min_freq].index)


def _typed_fingerprint(ds: dict, tool: str) -> set:
    """Set of (interaction_type, residue) contacts observed for a (dataset, tool)."""
    s = ds["inter"][ds["inter"]["tool"] == tool]
    if s.empty:
        return set()
    return set(zip(s["interaction_type"].astype(str), s["res"].astype(str)))


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return np.nan
    u = len(a | b)
    return len(a & b) / u if u else np.nan


def _overlap_rows(datasets: Dict[str, dict], tools: List[str], min_freq: float) -> List[dict]:
    rows = []
    for tool in tools:
        exp_res = _char_residue_set(datasets["exp"], tool, min_freq)
        bench_res = _char_residue_set(datasets["bench"], tool, min_freq)
        exp_typ = _typed_fingerprint(datasets["exp"], tool)
        bench_typ = _typed_fingerprint(datasets["bench"], tool)
        rows.append({
            "tool": tool,
            "n_exp_poses": _poses_per_tool(datasets["exp"]["summary"], tool),
            "n_bench_poses": _poses_per_tool(datasets["bench"]["summary"], tool),
            "shared": len(exp_res & bench_res),
            "exp_only": len(exp_res - bench_res),
            "bench_only": len(bench_res - exp_res),
            "jaccard_res": _jaccard(exp_res, bench_res),
            "jaccard_typed": _jaccard(exp_typ, bench_typ),
            "shared_res": sorted(exp_res & bench_res),
        })
    return rows


def fig_fingerprint_overlap(datasets: Dict[str, dict], tools: List[str], out: Path,
                            min_freq: float):
    rows = _overlap_rows(datasets, tools, min_freq)
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(10, 10.5),
                                   gridspec_kw={"height_ratios": [1.0, 1.4]})
    y = np.arange(len(tools))[::-1]

    # (A) per-tool residue-set decomposition. Each row is in its TOOL colour; the three
    # set-membership segments are told apart by fill STYLE — Benchmark-only solid,
    # Experimental-only hatched (the 09f set convention), Shared a brighter solid.
    seg_style = {"exp_only":   dict(alpha=0.55, hatch="///"),
                 "shared":     dict(alpha=0.92, hatch=""),
                 "bench_only": dict(alpha=0.50, hatch="")}
    for yi, r in zip(y, rows):
        col = TOOL_COLORS[r["tool"]]; left = 0.0
        for part in ("exp_only", "shared", "bench_only"):
            w = r[part]; ss = seg_style[part]
            axA.barh(yi, w, left=left, color=col, alpha=ss["alpha"], hatch=ss["hatch"] or None,
                     edgecolor="white", lw=0.9, zorder=2)
            if w:
                axA.text(left + w / 2, yi, str(w), ha="center", va="center", fontsize=8.5,
                         color="#111111", fontweight="bold",
                         path_effects=[_pe.withStroke(linewidth=2.4, foreground="white")])
            left += w
        j = r["jaccard_res"]
        axA.text(left + 0.4, yi, f"Jaccard {j:.2f}" if j == j else "Jaccard n/a",
                 ha="left", va="center", fontsize=8.5, color="#333333")
    axA.set_yticks(y); axA.set_yticklabels([TOOL_LABEL[t] for t in tools], fontsize=10)
    axA.set_xlabel(f"Number of characteristic Orai residues (contacted in ≥ {min_freq:.0%} of a tool's poses)")
    axA.set_title("Overlap of characteristic contact fingerprints", fontsize=11)
    axA.set_xlim(0, max(6, max((r["exp_only"] + r["shared"] + r["bench_only"]) for r in rows) * 1.28))
    axA.grid(axis="x", alpha=0.3); axA.set_axisbelow(True)

    # (B) Jaccard per tool: residue-level vs typed-contact, tool-coloured bars with a
    # two-level x-axis (residue / typed sub-ticks under each tool's group label).
    box_w = 0.34; gap = box_w + 0.06
    spans = {}
    tick_pos, tick_lab = [], []
    for ti, tool in enumerate(tools):
        r = next(rr for rr in rows if rr["tool"] == tool)
        xs = []
        for gi, (jkey, sub, alpha) in enumerate((("jaccard_res", "residue", 0.85),
                                                 ("jaccard_typed", "typed", 0.42))):
            pos = ti + (gi - 0.5) * gap
            v = 100 * (r[jkey] if r[jkey] == r[jkey] else 0.0)
            axB.bar(pos, v, box_w, color=TOOL_COLORS[tool], alpha=alpha,
                    edgecolor=TOOL_COLORS[tool], lw=1.0, zorder=2)
            axB.text(pos, v + 1.4, f"{v:.0f}%" if r[jkey] == r[jkey] else "n/a",
                     ha="center", va="bottom", fontsize=8)
            tick_pos.append(pos); tick_lab.append(sub); xs.append(pos)
        spans[tool] = xs
    axB.set_xticks(tick_pos); axB.set_xticklabels(tick_lab, fontsize=8.4)
    axB.tick_params(axis="x", length=0, pad=5)
    _group_x_axis(axB, spans, box_w, y_line=-0.11, y_label=-0.135, fontsize=10)
    axB.set_ylabel("Fingerprint overlap — Jaccard (%)"); axB.set_ylim(0, 108)
    axB.set_title("How similarly the two sets engage Orai\n(residue-level vs typed contacts)",
                  fontsize=10.5)
    axB.grid(axis="y", alpha=0.3); axB.set_axisbelow(True)

    _label_panels([axA, axB])
    # single figure-level legend for the set-membership decomposition (grey = style only)
    over_handles = [Patch(facecolor="0.5", alpha=seg_style["shared"]["alpha"], edgecolor="0.3",
                          label="Shared by both sets"),
                    Patch(facecolor="0.5", alpha=seg_style["bench_only"]["alpha"], edgecolor="0.3",
                          label="Benchmark only"),
                    Patch(facecolor="0.5", alpha=seg_style["exp_only"]["alpha"], edgecolor="0.3",
                          hatch="///", label="Experimental only")]
    fig.legend(handles=over_handles, loc="upper center", ncol=3,
               fontsize=_legend_fontsize(10), frameon=True, bbox_to_anchor=(0.5, 1 - 0.80 / 10.5))
    fig.suptitle("Orai contact-fingerprint overlap — Benchmark vs Experimental ligands, per tool\n"
                 f"(same receptor both sets; characteristic residues = contacted in ≥ {min_freq:.0%} "
                 "of a tool's poses)", fontsize=12.5, y=1 - 0.25 / 10.5)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.6 / 10.5))
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print(f"  wrote {out}")


# ── CSV exports ──────────────────────────────────────────────────────────────────

def write_csvs(datasets: Dict[str, dict], tools: List[str], types: List[str],
               residues: List[str], overlap_rows: List[dict], out_dir: Path):
    # pose totals (per pose)
    tot = []
    for key in CAT_ORDER:
        ds = datasets[key]
        s = ds["summary"][["tool", "protein", "ligand", "pose_name", "total_interactions"]].copy()
        s.insert(0, "dataset", key); s.insert(1, "category", ds["label"])
        tot.append(s)
    pd.concat(tot, ignore_index=True).to_csv(out_dir / "pandamap_pose_totals.csv", index=False)

    # type profile (per dataset × tool × type: mean count/pose)
    prof = []
    for key in CAT_ORDER:
        ds = datasets[key]
        for tool in tools:
            s = ds["summary"][ds["summary"]["tool"] == tool]
            n = int(s["pose_name"].nunique())
            for t in types:
                if t in s.columns:
                    prof.append(dict(dataset=key, category=ds["label"], tool=tool, n_poses=n,
                                     interaction_type=t,
                                     mean_per_pose=float(pd.to_numeric(s[t], errors="coerce").mean())
                                     if n else np.nan))
    pd.DataFrame(prof).to_csv(out_dir / "pandamap_type_profile.csv", index=False)

    # residue hot-spots (per dataset × tool × residue: contact frequency)
    hot = []
    for key in CAT_ORDER:
        ds = datasets[key]
        for tool in tools:
            fr = _residue_freq(ds, tool)
            n = _poses_per_tool(ds["summary"], tool)
            for res in residues:
                hot.append(dict(dataset=key, category=ds["label"], tool=tool, residue=res,
                                n_poses=n, n_contacting=int(round(fr.get(res, 0.0) * n)),
                                frac_poses=float(fr.get(res, 0.0))))
    pd.DataFrame(hot).to_csv(out_dir / "pandamap_residue_hotspots.csv", index=False)

    # fingerprint overlap (per tool)
    ov = [{k: v for k, v in r.items() if k != "shared_res"} for r in overlap_rows]
    for r, src in zip(ov, overlap_rows):
        r["shared_residues"] = ";".join(_pretty_res(x) for x in src["shared_res"])
    pd.DataFrame(ov).to_csv(out_dir / "pandamap_fingerprint_overlap.csv", index=False)
    print(f"  wrote 4 CSVs → {out_dir}")


# ── companion stats (inferential tests kept OFF the panels) ──────────────────────

def _fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and p != p):
        return "p=n/a"
    return f"p={su.fmt_p(p)} {su.p_stars(p)}".rstrip() if _HAVE_SU else f"p={p:.3g}"


def write_stats(datasets: Dict[str, dict], tools: List[str], types: List[str],
                residues: List[str], overlap_rows: List[dict], min_freq: float,
                out: Path):
    L: List[str] = []
    A = L.append
    A("PandaMap protein–ligand interactions on Orai — Experimental (JKU) vs Benchmark ligands")
    A("#" * 90)
    A("")
    A("Scope: PoseBusters-VALID poses OUTSIDE the transmembrane (both PandaMap runs are built")
    A("from posebusters_filtered_results.no_tm.csv + pb_valid_only:true). Same four Orai1 MD")
    A("frames both sets → the residue axis is directly comparable across ligand categories.")
    A("Toolchains (one variant per tool, matched across datasets): "
      + ", ".join(f"{TOOL_LABEL[t]}={datasets['bench']['chosen'].get(t, 'n/a')}" for t in tools) + ".")
    A("")
    for key in CAT_ORDER:
        ds = datasets[key]
        per_tool = ", ".join(f"{TOOL_LABEL[t]} {_poses_per_tool(ds['summary'], t)}p/"
                             f"{ds['summary'][ds['summary'].tool == t].groupby(['protein','ligand']).ngroups}c"
                             for t in tools)
        A(f"  {ds['label']:<20} {ds['root']}  (mtime {ds['mtime']})")
        A(f"       poses/complexes per tool: {per_tool}")
        for nt in ds["notes"]:
            A(f"       note: {nt}")
    A("")
    A("EquiBind caveat: it places most ligands in the TM pore, so its PB-valid non-TM pose count")
    A("is tiny on BOTH sides (esp. Experimental) — its distributions are EXPLORATORY, read n.")
    A("Unit of analysis: unpaired (different ligand sets). Pose-level medians are pseudoreplicated")
    A("(poses within a complex are correlated); the tests below aggregate to the per-complex mean")
    A("first, so the unit is the (frame × ligand) complex, not the pose.")
    A("")

    # 1. total interactions per pose
    A("=" * 90)
    A("1. Total interactions per pose — Benchmark vs Experimental, per tool")
    A("=" * 90)
    A("   Mann–Whitney U + Cliff's δ on the per-complex MEAN total (δ>0 ⇒ Benchmark higher).")
    praw, pack = [], []
    for tool in tools:
        b = _per_complex_means(datasets["bench"]["summary"], tool, "total_interactions")
        e = _per_complex_means(datasets["exp"]["summary"], tool, "total_interactions")
        if _HAVE_SU and b.size and e.size:
            res = su.mannwhitney_cliffs(b, e)
            pack.append((tool, b, e, res)); praw.append(res["p"])
        else:
            pack.append((tool, b, e, None))
    adj = su.holm(praw) if (_HAVE_SU and praw) else praw
    ai = 0
    for tool, b, e, res in pack:
        # pooled pose-level medians (descriptive)
        bp = pd.to_numeric(datasets["bench"]["summary"].query("tool==@tool")["total_interactions"],
                           errors="coerce").dropna()
        ep = pd.to_numeric(datasets["exp"]["summary"].query("tool==@tool")["total_interactions"],
                           errors="coerce").dropna()
        head = (f"   {TOOL_LABEL[tool]:<15} pose-median Bench {bp.median() if len(bp) else float('nan'):.1f} "
                f"(n={len(bp)}) vs Exp {ep.median() if len(ep) else float('nan'):.1f} (n={len(ep)})")
        if res is None:
            A(head + "   [no test — empty complex set on one side]"); continue
        pa = adj[ai]; ai += 1
        A(head + f"   | per-complex med Bench {np.median(b):.1f} (nc={b.size}) vs "
          f"Exp {np.median(e):.1f} (nc={e.size})  {_fmt_p(res['p'])} (Holm {_fmt_p(pa)}) "
          f"Cliff δ={res['cliffs_delta']:+.2f}")
    A("")

    # 2. interaction-type profile
    A("=" * 90)
    A("2. Interaction-type profile — Benchmark vs Experimental, per tool")
    A("=" * 90)
    A("   Per (tool, type): Mann–Whitney U + Cliff's δ on the per-complex mean count per pose,")
    A("   BH-FDR across the interaction-type family within each tool. Only types present shown.")
    for tool in tools:
        A(f"   ── {TOOL_LABEL[tool]} ──")
        rows, praw2 = [], []
        for t in types:
            b = _per_complex_means(datasets["bench"]["summary"], tool, t)
            e = _per_complex_means(datasets["exp"]["summary"], tool, t)
            if _HAVE_SU and b.size and e.size and (b.sum() + e.sum()) > 0:
                res = su.mannwhitney_cliffs(b, e)
                rows.append((t, b, e, res)); praw2.append(res["p"])
            else:
                rows.append((t, b, e, None))
        qs = su.bh_fdr(praw2) if (_HAVE_SU and praw2) else praw2
        qi = 0
        for t, b, e, res in rows:
            if res is None:
                continue
            q = qs[qi]; qi += 1
            A(f"       {pretty_itype(t):<18} mean/pose Bench {np.mean(b):.2f} vs Exp {np.mean(e):.2f}"
              f"   {_fmt_p(res['p'])} (BH {_fmt_p(q)}) Cliff δ={res['cliffs_delta']:+.2f}")
    A("")

    # 3. residue hot-spots
    A("=" * 90)
    A("3. Orai residue hot-spots — contact frequency, Benchmark vs Experimental, per tool")
    A("=" * 90)
    A("   Contact frequency = fraction of a (tool, dataset)'s poses contacting the residue.")
    A("   Two-proportion Fisher exact per residue (POSE-LEVEL — pseudoreplicated, descriptive only;")
    A("   Holm across the shown residues within each tool). Δpp = Benchmark% − Experimental%.")
    for tool in tools:
        A(f"   ── {TOOL_LABEL[tool]} ──")
        nb = _poses_per_tool(datasets["bench"]["summary"], tool)
        ne = _poses_per_tool(datasets["exp"]["summary"], tool)
        frb = _residue_freq(datasets["bench"], tool)
        fre = _residue_freq(datasets["exp"], tool)
        rows, praw3 = [], []
        for res in residues:
            kb = int(round(frb.get(res, 0.0) * nb)); ke = int(round(fre.get(res, 0.0) * ne))
            if _HAVE_SU and nb and ne:
                tp = su.two_proportion_test(kb, nb, ke, ne)
                rows.append((res, kb, ke, tp)); praw3.append(tp.get("p"))
            else:
                rows.append((res, kb, ke, None))
        # praw3 has exactly one entry per computed test, in row order — adjust it
        # directly so it stays aligned with the `j` counter below.
        adj3 = su.holm(praw3) if (_HAVE_SU and praw3) else praw3
        j = 0
        for res, kb, ke, tp in rows:
            bpc = 100 * kb / nb if nb else float("nan")
            epc = 100 * ke / ne if ne else float("nan")
            line = f"       {_pretty_res(res):<9} Bench {bpc:5.1f}% ({kb}/{nb}) vs Exp {epc:5.1f}% ({ke}/{ne})  Δ={bpc-epc:+5.1f}pp"
            if tp is not None:
                pa = adj3[j]; j += 1
                line += f"  {_fmt_p(tp.get('p'))} (Holm {_fmt_p(pa)})"
            A(line)
    A("")

    # 4. fingerprint overlap
    A("=" * 90)
    A("4. Contact-fingerprint overlap — Benchmark vs Experimental, per tool")
    A("=" * 90)
    A(f"   Characteristic Orai residue set = residues contacted in ≥ {min_freq:.0%} of a tool's poses.")
    A("   Jaccard(res) on those sets; Jaccard(typed) on (interaction_type, residue) contacts.")
    A(f"   CAVEAT — set-size asymmetry: the ≥{min_freq:.0%} threshold is an ABSOLUTE count of")
    A("   ~ceil(freq·n_poses) poses, i.e. only a handful for Experimental (36 poses/tool) vs")
    A("   hundreds for Benchmark (~3400/tool). The huge Benchmark set spreads its contacts across")
    A("   ~300 diverse ligands, so FEWER residues clear its threshold; the small, homogeneous")
    A("   Experimental set (3 ligands × 4 frames, docked repeatedly) concentrates contacts, so MORE")
    A("   residues clear its low absolute threshold. Hence 'exp-only' can EXCEED 'bench-only' — a")
    A("   property of set size/homogeneity, NOT of ligand promiscuity. Read the SHARED set")
    A("   (threshold-robust: frequent in BOTH) and the per-residue rows (section 3), not Jaccard alone.")
    for r in overlap_rows:
        tool = r["tool"]
        shared = ", ".join(_pretty_res(x) for x in r["shared_res"]) or "(none)"
        A(f"   {TOOL_LABEL[tool]:<15} shared {r['shared']} | exp-only {r['exp_only']} | "
          f"bench-only {r['bench_only']}   Jaccard(res)={r['jaccard_res']:.2f} "
          f"Jaccard(typed)={r['jaccard_typed']:.2f}   (poses: Exp {r['n_exp_poses']}, Bench {r['n_bench_poses']})")
        A(f"       shared residues: {shared}")
    A("")
    A("=" * 90)
    A("CAVEATS")
    A("=" * 90)
    A("• Same 4 Orai frames both sets, but DIFFERENT ligands → comparisons are unpaired.")
    A("• Experimental = 3 ligands × 4 frames; EquiBind has ~a handful of PB-valid non-TM poses")
    A("  (both sets) → EquiBind rows and any Experimental cell with tiny n are EXPLORATORY.")
    A("• Section 3 Fisher tests are pose-level and pseudoreplicated (poses within a complex are")
    A("  correlated); they OVERSTATE significance and are descriptive support for the figures only.")
    A("• Panels carry descriptive content only; the inferential tests live here (project convention).")
    out.write_text("\n".join(L) + "\n")
    print(f"  wrote {out}")


# ── driver ───────────────────────────────────────────────────────────────────────

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exp-dir", default="pandamap_results/orai_jku", type=Path,
                    help="PandaMap output dir for the Experimental (JKU) ligand run.")
    ap.add_argument("--bench-dir", default="pandamap_results/orai_benchmark", type=Path,
                    help="PandaMap output dir for the Benchmark ligand run.")
    ap.add_argument("--out-dir", default="pandamap_results/orai_interaction_compare", type=Path,
                    help="Where the *_compare figures / CSVs / stats are written.")
    ap.add_argument("--diffdock-variant", default="smina",
                    help="DiffDock optimiser variant to select in both datasets (default smina).")
    ap.add_argument("--equibind-variant", default="gnina",
                    help="EquiBind refine variant to select in both datasets (default gnina).")
    ap.add_argument("--top-n", type=int, default=3,
                    help="Defensive cap: top-N poses per (tool, complex) by pose_rank (default 3).")
    ap.add_argument("--top-n-poses", type=int, default=0,
                    help="If >0, also restrict every tool to its top-N NATIVE-rank poses per "
                         "(protein, ligand) — AutoDock=Vina, DiffDock=confidence, EquiBind=gnina "
                         "energy — to match the yield / TM / pose-cluster figures. Ranks are read "
                         "from --exp-posebusters-csv / --bench-posebusters-csv (EquiBind's "
                         "PandaMap pose_rank is a 999 placeholder, so a rank threshold can't do it).")
    ap.add_argument("--exp-posebusters-csv",
                    default="posebusters_results/orai_jku/dock/posebusters_filtered_results.csv",
                    type=Path, help="Full per-pose PoseBusters CSV for the Experimental run "
                                    "(source of the native per-tool ranks for --top-n-poses).")
    ap.add_argument("--bench-posebusters-csv",
                    default="posebusters_results/orai_benchmark/dock/posebusters_filtered_results.csv",
                    type=Path, help="Full per-pose PoseBusters CSV for the Benchmark run.")
    ap.add_argument("--residue-top-n", type=int, default=15,
                    help="Number of top Orai residues on the hot-spot figure (default 15).")
    ap.add_argument("--overlap-min-freq", type=float, default=0.10,
                    help="A residue is 'characteristic' for a (tool, dataset) if contacted in at "
                         "least this fraction of its poses (default 0.10). Feeds fig 4 + stats.")
    ap.add_argument("--tools", nargs="*", default=None,
                    help="Subset/reorder tools (choices: autodock diffdock equibind). "
                         "Default: all present in BOTH datasets.")
    args = ap.parse_args(argv)

    for p in (args.exp_dir, args.bench_dir):
        if not (Path(p) / "pandamap_pose_summary.csv").exists():
            print(f"ERROR: no pandamap_pose_summary.csv under {p} — run run_pandamap.py first "
                  f"(config: pandamap_orai_{'jku' if 'jku' in str(p) else 'benchmark'}_config.yaml).")
            return 1

    datasets = {
        "exp": load_dataset("exp", args.exp_dir, args.diffdock_variant,
                            args.equibind_variant, args.top_n,
                            args.exp_posebusters_csv, args.top_n_poses),
        "bench": load_dataset("bench", args.bench_dir, args.diffdock_variant,
                              args.equibind_variant, args.top_n,
                              args.bench_posebusters_csv, args.top_n_poses),
    }

    # tools present in BOTH datasets (keep canonical order), unless overridden
    if args.tools:
        tools = [t for t in args.tools if t in CANONICAL_TOOLS]
    else:
        tools = [t for t in CANONICAL_TOOLS
                 if t in set(datasets["exp"]["summary"]["tool"])
                 and t in set(datasets["bench"]["summary"]["tool"])]
    if not tools:
        print("ERROR: no tool is present in both datasets — nothing to compare.")
        return 1

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Comparing PandaMap interactions on Orai (tools: {', '.join(tools)}):")
    for key in CAT_ORDER:
        ds = datasets[key]
        counts = ", ".join(f"{TOOL_LABEL[t]}={_poses_per_tool(ds['summary'], t)}" for t in tools)
        print(f"  {ds['label']:<20} chosen={ds['chosen']}  poses: {counts}")

    types = _present_types(datasets)
    residues = _top_shared_residues(datasets, tools, args.residue_top_n)
    overlap_rows = _overlap_rows(datasets, tools, args.overlap_min_freq)

    fig_total_interactions(datasets, tools, out_dir / "fig_total_interactions_compare.png")
    fig_type_profile(datasets, tools, out_dir / "fig_type_profile_compare.png")
    fig_residue_hotspots(datasets, tools, out_dir / "fig_residue_hotspots_compare.png",
                         args.residue_top_n)
    fig_fingerprint_overlap(datasets, tools, out_dir / "fig_fingerprint_overlap_compare.png",
                            args.overlap_min_freq)

    write_csvs(datasets, tools, types, residues, overlap_rows, out_dir)
    if _HAVE_SU:
        write_stats(datasets, tools, types, residues, overlap_rows, args.overlap_min_freq,
                    out_dir / "orai_pandamap_interaction_compare_stats.txt")
    else:
        print("  [stats] stats_utils unavailable — companion stats skipped.")
    print("done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
