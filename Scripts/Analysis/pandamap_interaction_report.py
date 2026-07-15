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

# 16 PandaMap interaction types (stable column order).
INTERACTION_TYPES = [
    "hydrogen_bonds", "carbon_pi", "pi_pi_stacking", "donor_pi", "amide_pi",
    "hydrophobic", "ionic", "halogen_bonds", "cation_pi", "metal_coordination",
    "salt_bridge", "covalent", "alkyl_pi", "attractive_charge", "pi_cation",
    "repulsion",
]

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
        elif t in ("raw", "smina"):
            refine = t
        elif t in ("clampON", "clampOFF"):
            clamp = t
    return pocket, refine, clamp


def pretty_method(m: str) -> str:
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
    if m == "crystal":
        return "Crystal (native)"
    if not m.startswith("equibind"):
        return m
    pocket, refine, clamp = _eq_tokens(m)
    parts = [p for p in (pocket,
                         "smina-opt" if refine == "smina" else "raw" if refine else None,
                         "clamp on" if clamp == "clampON" else "clamp off" if clamp else None)
             if p]
    return f"EquiBind ({', '.join(parts)})" if parts else "EquiBind"


def method_sort_key(m: str):
    if m == "crystal":
        return (-1, 0, 0)
    if m == "autodock":
        return (0, 0, 0)
    if m.startswith("diffdock"):
        return (1, {"diffdock": 0, "diffdock_smina": 1, "diffdock_gnina": 2}.get(m, 3), 0)
    if m.startswith("equibind"):
        p, r, c = _eq_tokens(m)
        return (2, _POCKET_ORDER.get(p, 9), {None: 0, "raw": 1, "smina": 2}.get(r, 0))
    return (3, 0, 0)


def ordered_methods(methods) -> list[str]:
    return sorted([m for m in methods if m != "crystal"], key=lambda m: (method_sort_key(m), m))


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


# ── shared figure finalisation (pose-selection scope footer) ────────────────
# run_pandamap keeps only each method's top-N ranked poses per complex
# (poses_per_combo / select_top_n), so the whole interaction analysis is built
# from that top-N set. Every saved figure carries a uniform footer stating this,
# so each graph makes the pose-selection scope explicit. Set once in main() from
# the config's poses_per_combo (or inferred from the loaded data).
_POSE_SCOPE_CAPTION: str = ""


def set_pose_scope_caption(n: "int | None") -> None:
    """Record the top-N pose-selection scope stamped as a footer on every figure."""
    global _POSE_SCOPE_CAPTION
    _POSE_SCOPE_CAPTION = (
        f"Built from each method's top {n} ranked pose(s) per complex "
        f"(fewer where a method has under {n} valid poses)." if n else "")


def _finalize_fig(fig, out: Path, dpi: int = 160) -> None:
    """Stamp the shared pose-selection-scope footer (if set), then save + close."""
    if _POSE_SCOPE_CAPTION:
        fig.text(0.5, 0.006, _POSE_SCOPE_CAPTION, ha="center", va="bottom",
                 fontsize=7, color="#666666", style="italic")
    fig.savefig(out, dpi=dpi)
    plt.close(fig)


# ── 01: per-method interaction profile ─────────────────────────────────────

def plot_profile(summary: pd.DataFrame, order, out: Path) -> None:
    # Cleveland dot plot: one row per interaction type, one coloured dot per
    # method. Replaces a 16×N grouped-bar wall — far easier to compare methods
    # within a type, and complements the absolute-magnitude heatmap (fig 02).
    present = [t for t in INTERACTION_TYPES if t in summary.columns and summary[t].sum() > 0]
    means = summary.groupby("method")[present].mean().reindex(order).dropna(how="all")
    methods = list(means.index)
    cmap = plt.get_cmap("tab10")
    mcolors = {m: cmap(i % 10) for i, m in enumerate(methods)}
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
    ax.set_title("Protein-ligand interaction profile by docking method (PandaMap)")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(axis="x", alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out)


def plot_type_heatmap_box(summary: pd.DataFrame, order, out_heat: Path, out_box: Path) -> None:
    present = [t for t in INTERACTION_TYPES if t in summary.columns and summary[t].sum() > 0]
    means = summary.groupby("method")[present].mean().reindex(order).dropna(how="all")
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(present)), max(5, 0.5 * len(means))))
    sns.heatmap(means, annot=True, fmt=".1f", cmap="viridis",
                yticklabels=[pretty_method(m) for m in means.index],
                cbar_kws={"label": "mean / pose"}, ax=ax)
    ax.set_title("Mean interactions per pose — type × method")
    fig.tight_layout(); _finalize_fig(fig, out_heat)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(means)), 5))
    data = [summary[summary["method"] == m]["total_interactions"].values for m in means.index]
    ax.boxplot(data, showmeans=True)
    ax.set_xticks(range(1, len(means) + 1))
    ax.set_xticklabels([pretty_method(m) for m in means.index], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylabel("Total interactions per pose")
    ax.set_title("Distribution of total interactions per pose")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); _finalize_fig(fig, out_box)


# ── 03: residue hot-spots ───────────────────────────────────────────────────

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
    fig, ax = plt.subplots(figsize=(max(8, 0.9 * mat.shape[1]), max(6, 0.32 * len(mat))))
    sns.heatmap(mat, annot=False, cmap="rocket_r", vmin=0, vmax=1,
                xticklabels=[pretty_method(m) for m in mat.columns],
                cbar_kws={"label": "fraction of poses contacting residue"}, ax=ax)
    ax.set_title(f"Residue hot-spots — top {len(mat)} contacted residues × method")
    ax.set_ylabel("Residue"); ax.set_xlabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(); _finalize_fig(fig, out)
    return mat


# ── 04/05: native recovery + fingerprint similarity ─────────────────────────

def native_recovery(inter: pd.DataFrame, crystal: pd.DataFrame, order,
                    out_dir: Path) -> pd.DataFrame:
    """Per-pose precision/recall/F1 of interactions vs the crystal fingerprint."""
    cryst_fp = {}
    for (p, l), g in crystal.groupby(["protein", "ligand"]):
        cryst_fp[(p, l)] = pose_fingerprint(g, "typed")

    recs = []
    for (method, p, l, pose), g in inter.groupby(["method", "protein", "ligand", "pose_name"]):
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

    # bar chart of mean precision/recall/F1
    fig, ax = plt.subplots(figsize=(max(9, 1.2 * len(summ)), 5.5))
    x = np.arange(len(summ)); w = 0.25
    for i, col in enumerate(["precision", "recall", "f1"]):
        ax.bar(x + (i - 1) * w, summ[col].values, w, label=col)
    ax.set_xticks(x)
    ax.set_xticklabels([pretty_method(m) for m in summ.index], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylabel("score"); ax.set_ylim(0, 1)
    ax.set_title("Native-interaction recovery vs crystal (mean over poses)\n"
                 "fingerprint = (interaction_type, residue)")
    ax.legend(); ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "04_native_recovery.png")

    # F1 CDF per method
    fig, ax = plt.subplots(figsize=(8, 5))
    for m in summ.index:
        vals = np.sort(rdf[rdf["method"] == m]["f1"].dropna().values)
        if len(vals):
            ax.plot(vals, np.arange(1, len(vals) + 1) / len(vals) * 100, lw=2, label=pretty_method(m))
    ax.set_xlabel("F1 of native interaction recovery"); ax.set_ylabel("cumulative % of poses")
    ax.set_title("Native-interaction recovery — F1 CDF"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "04b_native_recovery_cdf.png")
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
    colors = _method_colors(order)
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

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=True, sharey=True)
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
        ax.set_xlabel("Pose rank k (1 = tool's top-ranked pose)")
        ax.set_title(title)
        ax.set_ylim(0, 1); ax.grid(alpha=0.3)
    axes[0].set_ylabel("Mean score vs crystal\n(typed interaction recovery)")
    axes[0].legend(fontsize=7, loc="lower left")
    _label_panels(axes)
    fig.suptitle("Native-interaction recovery vs crystal by pose rank  —  "
                 "fingerprint = (interaction_type, residue), PB-valid poses only\n"
                 "Rank basis: AutoDock / DiffDock native confidence rank, EquiBind "
                 "gnina-affinity rank (a tool lacking a rank is shown as a flat baseline)")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    _finalize_fig(fig, out)


def plot_fingerprint_similarity(inter: pd.DataFrame, crystal: pd.DataFrame, order, out: Path,
                                depth: int = 1) -> None:
    """Mean cross-method Jaccard of the top-`depth` pose(s) per pair, + vs crystal.

    ``depth=1`` uses each method's single best (lowest pose_rank) pose. ``depth>1`` uses
    the UNION of that method's top-`depth` poses' typed fingerprints per complex — the
    interaction repertoire the tool covers across its top ranked poses — so the matrix
    compares interaction coverage rather than a single pose.
    """
    cryst_fp = {(p, l): pose_fingerprint(g, "typed") for (p, l), g in crystal.groupby(["protein", "ligand"])}
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
    for (method, p, l), g in sel.groupby(["method", "protein", "ligand"], sort=False):
        fp_by[(method, p, l)] = pose_fingerprint(g, "typed")   # union over the top-depth poses
    methods = ordered_methods({m for (m, _, _) in fp_by})
    labels = methods + (["crystal"] if cryst_fp else [])
    mat = np.full((len(labels), len(labels)), np.nan)
    pairs = {(p, l) for (_, p, l) in fp_by}
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            vals = []
            for (p, l) in pairs:
                fa = cryst_fp.get((p, l)) if a == "crystal" else fp_by.get((a, p, l))
                fb = cryst_fp.get((p, l)) if b == "crystal" else fp_by.get((b, p, l))
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
    ax.set_title(f"Interaction-fingerprint similarity ({pose_desc})")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(); _finalize_fig(fig, out)


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
    fig.suptitle("Chemical level — interaction type vs ligand chemotype", fontsize=14)
    fig.tight_layout(); _finalize_fig(fig, out_dir / "06_interaction_vs_chemotype.png")
    # correlation table
    cols = ["mw", "n_aromatic_rings", "hbd", "hba", "n_halogen", "formal_charge", "logp"]
    cols = [c for c in cols if c in df.columns]
    itypes = [t for t in INTERACTION_TYPES if t in df.columns and df[t].sum() > 0]
    corr = df[cols + itypes].corr().loc[cols, itypes]
    corr.round(3).to_csv(out_dir / "chem_descriptor_interaction_corr.csv")
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(itypes)), max(4, 0.6 * len(cols))))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", center=0, vmin=-1, vmax=1, ax=ax)
    ax.set_title("Correlation: ligand descriptor × interaction type")
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
    ax.set_title("Residue-class preference per interaction type")
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
    ax.set_title("Ligand-atom element per interaction type")
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
    """Stable per-method colour map (tab10), shared across the new figures."""
    cmap = plt.get_cmap("tab10")
    return {m: cmap(i % 10) for i, m in enumerate(order)}


def build_recovery_detail(inter: pd.DataFrame, crystal: pd.DataFrame):
    """Per-pose typed & residue-contact recovery detail vs the crystal fingerprint.

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
    for (p, l), g in crystal.groupby(["protein", "ligand"]):
        t = pose_fingerprint(g, "typed")
        cryst_typed[(p, l)] = t
        cryst_contact[(p, l)] = pose_fingerprint(g, "contact")
        by = {}
        for (it, rn, rnum, ch) in t:
            by.setdefault(it, set()).add((rn, rnum, ch))
        cryst_by_type[(p, l)] = by

    has_rank = "pose_rank" in inter.columns
    keys = ["method", "protein", "ligand", "pose_name"] + (["pose_rank"] if has_rank else [])
    pose_rows, type_rows = [], []
    for k, g in inter.groupby(keys, sort=False):
        if has_rank:
            method, p, l, pose, rank = k
        else:
            (method, p, l, pose), rank = k, np.nan
        ref = cryst_typed.get((p, l))
        if ref is None:
            continue
        fp = pose_fingerprint(g, "typed")
        fp_c = pose_fingerprint(g, "contact")
        ref_c = cryst_contact.get((p, l), set())
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
        })
        pose_by_type: dict[str, set] = {}
        for (it, rn, rnum, ch) in fp:
            pose_by_type.setdefault(it, set()).add((rn, rnum, ch))
        ref_by_type = cryst_by_type.get((p, l), {})
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
    # native prevalence per COMPLEX (dedup poses) so n is independent of the depth bucket
    native = (per_type.drop_duplicates(["protein", "ligand", "interaction_type"])
              .groupby("interaction_type")["n_native_type"].sum())
    types = [t for t in INTERACTION_TYPES if native.get(t, 0) > 0]
    methods = [m for m in order if m in set(agg["method"])]
    rec = agg.pivot(index="interaction_type", columns="method", values="recall").reindex(index=types, columns=methods)
    prc = agg.pivot(index="interaction_type", columns="method", values="precision").reindex(index=types, columns=methods)
    ylabels = [f"{pretty_itype(t)}  (n={int(native.get(t, 0))})" for t in types]
    fig, axes = plt.subplots(
        1, 2, figsize=(max(11, 2.4 * len(methods) + 6), max(5.5, 0.6 * len(types) + 2)))
    for i, (ax, mat, name) in enumerate(((axes[0], rec, "Recall"), (axes[1], prc, "Precision"))):
        sns.heatmap(mat, annot=True, fmt=".2f", annot_kws={"fontsize": 9}, cmap="viridis",
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
    fig.suptitle("Native-interaction recovery resolved by interaction type  —  " + pose_desc + "\n"
                 "(A) recall = native contacts reproduced;   (B) precision = predicted contacts that are native",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    _finalize_fig(fig, out)
    agg.round(4).to_csv(out.with_suffix(".csv"), index=False)


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
    ax.set_title(f"Matched / missed / spurious interaction contacts vs the crystal — {pose_desc}\n"
                 "(typed contact = interaction type + residue)")
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
    ax.set_title("Matched / missed / spurious interaction contacts vs the crystal, by pose rank\n"
                 "(typed contact = interaction type + residue)")
    # reserve a right margin so the two legends sit fully inside the figure
    fig.subplots_adjust(left=0.08, right=0.68, top=0.88, bottom=0.11)
    from matplotlib.lines import Line2D
    tool_handles = [Line2D([0], [0], color=colors.get(m), lw=3) for m in methods]
    class_handles = [Line2D([0], [0], color="0.35", lw=2, ls=ls) for _, _, ls in classes]
    leg1 = ax.legend(tool_handles, [pretty_method(m) for m in methods], title="Tool",
                     loc="upper left", bbox_to_anchor=(1.03, 1.0), fontsize=9)
    ax.add_artist(leg1)
    ax.legend(class_handles, [lab for _, lab, _ in classes], title="Contact class",
              loc="upper left", bbox_to_anchor=(1.03, 0.5), fontsize=9)
    _finalize_fig(fig, out)
    (df.groupby(["method", "pose_rank"])[["tp", "fn", "fp"]].mean().round(3)
     .to_csv(out.with_suffix(".csv")))


def plot_residue_confusion(inter: pd.DataFrame, crystal: pd.DataFrame, order,
                           out: Path, top_n: int = 20) -> None:
    """11 — per-residue recovery of native contacts and per-residue spurious contacts.

    Residues are pooled by resname+resnum across complexes (as in fig 03). Left panel:
    fraction of a method's poses that recover each top native hot-spot residue (only
    counting poses whose own crystal contacts that residue). Right panel: the top
    residues each method contacts that its crystal does NOT (hallucinated contacts),
    as a fraction of that method's poses.
    """
    def _res(df):
        return df["resname"].astype(str) + df["resnum"].astype(str)
    cryst = crystal.copy(); cryst["res"] = _res(cryst)
    native_by_pair = cryst.groupby(["protein", "ligand"])["res"].apply(set).to_dict()
    native_freq = cryst.drop_duplicates(["protein", "ligand", "res"])["res"].value_counts()

    pc = inter.copy(); pc["res"] = _res(pc)
    pose_res = (pc.groupby(["method", "protein", "ligand", "pose_name"])["res"]
                .apply(set).reset_index())
    methods = [m for m in order if m in set(pose_res["method"])]
    poses_per_method = pose_res.groupby("method").size()

    recov_num: dict = {}; recov_den: dict = {}; spur_num: dict = {}
    for _, r in pose_res.iterrows():
        m = r["method"]; nat = native_by_pair.get((r["protein"], r["ligand"]), set())
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
    fig.suptitle(f"Per-residue native-contact recovery (top {len(top_native)} hot-spot residues) "
                 f"and hallucinated non-native contacts (top {len(top_spur)})", fontsize=13)
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
    x = np.arange(len(methods)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(methods)), 5.5))
    ax.bar(x - w / 2, strict.values, w, label="typed recall (interaction type + residue must match)", color="#4C72B0")
    ax.bar(x + w / 2, loose.values, w, label="loose recall (residue contacted, any interaction type)", color="#DD8452")
    for xi, (s, ll) in enumerate(zip(strict.values, loose.values)):
        if np.isfinite(s) and np.isfinite(ll):
            ax.annotate(f"gap {ll - s:.2f}", (xi, ll), textcoords="offset points",
                        xytext=(0, 4), ha="center", fontsize=7, color="#333333")
    ax.set_xticks(x); ax.set_xticklabels([pretty_method(m) for m in methods], rotation=20, ha="right")
    ax.set_ylabel("Mean recall of native contacts"); ax.set_ylim(0, 1)
    ax.set_title("Right residue, wrong chemistry: typed vs loose native recall\n"
                 "(gap = poses that hit the native residue but form a different interaction type)")
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
    x = np.arange(len(methods)); w = 0.38
    fig, ax = plt.subplots(figsize=(max(9, 2.0 * len(methods)), 5.5))
    ax.bar(x - w / 2, oracle.values, w, label="best native-F1 pose in the set (interaction oracle)", color="#55A868")
    ax.bar(x + w / 2, top1.values, w, label="native-F1 of the tool's own top-ranked pose", color="#937860")
    ax.set_xticks(x); ax.set_xticklabels([pretty_method(m) for m in methods], rotation=20, ha="right")
    ax.set_ylabel("Mean native-interaction recovery F1"); ax.set_ylim(0, 1)
    ax.set_title("Interaction-recovery F1: oracle vs the tool's top-ranked pose\n"
                 "gap = recovery the tool's own ranking leaves on the table", fontsize=11)
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
    legend = []
    for m in methods:
        sub = j[j["method"] == m]
        cov = len(sub); tot = int((per_pose["method"] == m).sum())
        rho = sub[["rmsd", "f1"]].corr(method="spearman").iloc[0, 1] if len(sub) > 2 else np.nan
        ax.scatter(sub["rmsd"], sub["f1"], s=12, alpha=0.35, color=colors.get(m), edgecolor="none")
        legend.append(f"{pretty_method(m)}  (n={cov}/{tot}, Spearman ρ={rho:.2f})")
    ax.axvline(2.0, color="#888888", ls="--", lw=1)
    ax.annotate("2 Å", (2.0, 0.02), color="#666666", fontsize=8, ha="left")
    ax.set_xlabel("RMSD to crystal ligand (Å)")
    ax.set_ylabel("Native-interaction recovery F1")
    ax.set_title("Geometry vs chemistry: does low RMSD predict recovering native contacts?")
    # Cap the RMSD axis at 10 Å — far-off decoy poses beyond this only compress the
    # informative 0–10 Å range (per-method Spearman ρ above is still over all poses).
    ax.set_xlim(0, 10); ax.set_ylim(-0.02, 1.02)
    ax.legend(legend, fontsize=8, loc="upper right"); ax.grid(alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout(); _finalize_fig(fig, out)
    j[key + ["rmsd", "f1", "recall", "precision"]].to_csv(out.with_suffix(".csv"), index=False)
    print(f"  RMSD-vs-recovery: joined {len(j)} poses "
          f"({', '.join(f'{m}={int((j.method==m).sum())}' for m in methods)}).")


def plot_pair_overlays(inter: pd.DataFrame, crystal: pd.DataFrame, per_pose: pd.DataFrame,
                       order, out: Path, n_pairs: int = 4, max_rows: int = 26) -> None:
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

    cryst_typed = {(p, l): pose_fingerprint(g, "typed") for (p, l), g in crystal.groupby(["protein", "ligand"])}
    best_pose = (inter.sort_values(["pose_rank", "pose_name"])
                 .drop_duplicates(["method", "protein", "ligand"]))  # top-ranked pose id per (method,pair)

    cmap = matplotlib.colors.ListedColormap(["#F2F2F2", "#C44E52", "#55A868", "#8172B3"])
    #                                          absent      native-only recovered   spurious
    fig, axes = plt.subplots(1, len(picks), figsize=(max(4.5 * len(picks), 9), max(6, 0.28 * max_rows + 2)))
    axes = np.atleast_1d(axes)
    for ax, (p, l) in zip(axes, picks):
        ref = cryst_typed.get((p, l), set())
        cols, fps = ["crystal"], [ref]
        for m in methods:
            row = best_pose[(best_pose["method"] == m) & (best_pose["protein"] == p) & (best_pose["ligand"] == l)]
            if row.empty:
                continue
            pose = row.iloc[0]["pose_name"]
            g = inter[(inter["method"] == m) & (inter["protein"] == p) &
                      (inter["ligand"] == l) & (inter["pose_name"] == pose)]
            cols.append(m); fps.append(pose_fingerprint(g, "typed"))
        contacts = sorted(set().union(*fps), key=lambda c: (c[0] not in {t for (t, *_ ) in ref}, str(c)))
        contacts = contacts[:max_rows]
        M = np.zeros((len(contacts), len(cols)))
        for j, (cn, fp) in enumerate(zip(cols, fps)):
            for i, c in enumerate(contacts):
                if c not in fp:
                    M[i, j] = 0
                elif cn == "crystal":
                    M[i, j] = 1                       # native (shown in crystal column)
                elif c in ref:
                    M[i, j] = 2                       # recovered native
                else:
                    M[i, j] = 3                       # spurious
        ax.imshow(M, aspect="auto", cmap=cmap, vmin=0, vmax=3)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels([pretty_method(c) for c in cols], rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(contacts)))
        ax.set_yticklabels([f"{it}:{rn}{rnum}" for (it, rn, rnum, ch) in contacts], fontsize=6)
        ax.set_title(f"{p}\nmean native-F1 {pair_f1.get((p, l), np.nan):.2f}", fontsize=9)
    handles = [matplotlib.patches.Patch(color=c, label=lbl) for c, lbl in
               (("#C44E52", "native contact (crystal)"), ("#55A868", "recovered by pose"),
                ("#8172B3", "spurious (non-native)"), ("#F2F2F2", "absent"))]
    # lift the legend into its own band (bbox y=0.05) and reserve the bottom 10% so it
    # clears the pose-scope footer that _finalize_fig stamps at the very bottom edge.
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.05),
               ncol=4, fontsize=8, frameon=False)
    _label_panels(axes)
    fig.suptitle("Crystal-vs-pose interaction overlays (top-ranked pose per method)", fontsize=13)
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
                         "kept in all charts/CSVs, relabelled 'EquiBind*'. Takes precedence "
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
                         "oracle-ranking it; kept in all charts/CSVs, relabelled 'DiffDock*'. "
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

    summary = pd.read_csv(in_dir / "pandamap_pose_summary.csv", low_memory=False)
    inter = pd.read_csv(in_dir / "pandamap_interactions.csv", low_memory=False)
    crystal_path = in_dir / "crystal_interactions.csv"
    crystal = pd.read_csv(crystal_path, low_memory=False) if crystal_path.exists() else pd.DataFrame()

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
            _LABEL_OVERRIDES[best_eq] = "EquiBind*"
            print(f"equibind-variant={equibind_variant}: keeping only '{best_eq}' "
                  "(shown as 'EquiBind*').")
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
            _LABEL_OVERRIDES[best_dd] = "DiffDock*"
            print(f"diffdock-variant={diffdock_variant}: keeping only '{best_dd}' "
                  "(shown as 'DiffDock*').")
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
        rec = native_recovery(inter, crystal, order, out_dir)
        if not rec.empty:
            print("\nNative-interaction recovery (mean per method):")
            print(rec.round(3).to_string())
        # top-1 pose (base name, kept) and top-5 union bucket
        plot_fingerprint_similarity(inter, crystal, order,
                                    out_dir / "05_fingerprint_similarity.png", depth=1)
        plot_fingerprint_similarity(inter, crystal, order,
                                    out_dir / "05_fingerprint_similarity_top5.png", depth=depth_pool)

        # Deeper crystal-vs-pose comparison (figs 09-15). One residue-level detail
        # pass feeds the type-resolved, decomposition, gap and oracle figures.
        per_pose, per_type = build_recovery_detail(inter, crystal)
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
                               out_dir / "11_residue_confusion.png", args.top_residues)
        plot_typed_vs_loose(per_pose, order, out_dir / "12_typed_vs_loose_gap.png")
        plot_native_f1_oracle(per_pose, order, out_dir / "13_native_f1_oracle.png")
        plot_rmsd_vs_recovery(per_pose, per_pose_metrics, order,
                              out_dir / "14_rmsd_vs_recovery.png")
        plot_pair_overlays(inter, crystal, per_pose, order,
                           out_dir / "15_pair_overlays.png")
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

    print(f"\nReport written to: {out_dir.resolve()}")
    for f in sorted(out_dir.glob("*.png")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
