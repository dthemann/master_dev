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


def select_best_diffdock(summary: pd.DataFrame, inter: pd.DataFrame,
                         oracle_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame, str | None]:
    """Drop every DiffDock optimizer variant except the single best-performing one.

    "Best" = the DiffDock variant (diffdock / diffdock_smina / diffdock_gnina) with the
    highest PB-Valid AND RMSD ≤ 2 Å rate (``oracle_pb_valid_and_rmsd2_%``, falling back to
    ``oracle_rmsd_le_2.0A_%`` for older summaries — see _variant_ranking_scores) in *oracle_csv*.
    Mirrors select_best_equibind; AutoDock/EquiBind/crystal rows are always kept. Returns
    the filtered (summary, inter) frames and the chosen variant key, or the inputs
    unchanged with ``None`` when no DiffDock variant is present, the summary is
    missing/unreadable, or no present variant has a score."""
    dd_present = sorted({m for m in summary["method"].astype(str).unique()
                         if m.startswith("diffdock")})
    if not dd_present:
        return summary, inter, None
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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


def plot_type_heatmap_box(summary: pd.DataFrame, order, out_heat: Path, out_box: Path) -> None:
    present = [t for t in INTERACTION_TYPES if t in summary.columns and summary[t].sum() > 0]
    means = summary.groupby("method")[present].mean().reindex(order).dropna(how="all")
    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(present)), max(5, 0.5 * len(means))))
    sns.heatmap(means, annot=True, fmt=".1f", cmap="viridis",
                yticklabels=[pretty_method(m) for m in means.index],
                cbar_kws={"label": "mean / pose"}, ax=ax)
    ax.set_title("Mean interactions per pose — type × method")
    fig.tight_layout(); fig.savefig(out_heat, dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(max(8, 1.1 * len(means)), 5))
    data = [summary[summary["method"] == m]["total_interactions"].values for m in means.index]
    ax.boxplot(data, showmeans=True)
    ax.set_xticks(range(1, len(means) + 1))
    ax.set_xticklabels([pretty_method(m) for m in means.index], rotation=45,
                       ha="right", rotation_mode="anchor", fontsize=8)
    ax.set_ylabel("Total interactions per pose")
    ax.set_title("Distribution of total interactions per pose")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout(); fig.savefig(out_box, dpi=160); plt.close(fig)


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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
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
    fig.tight_layout(); fig.savefig(out_dir / "04_native_recovery.png", dpi=160); plt.close(fig)

    # F1 CDF per method
    fig, ax = plt.subplots(figsize=(8, 5))
    for m in summ.index:
        vals = np.sort(rdf[rdf["method"] == m]["f1"].dropna().values)
        if len(vals):
            ax.plot(vals, np.arange(1, len(vals) + 1) / len(vals) * 100, lw=2, label=pretty_method(m))
    ax.set_xlabel("F1 of native interaction recovery"); ax.set_ylabel("cumulative % of poses")
    ax.set_title("Native-interaction recovery — F1 CDF"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(out_dir / "04b_native_recovery_cdf.png", dpi=160); plt.close(fig)
    return summ


def plot_fingerprint_similarity(inter: pd.DataFrame, crystal: pd.DataFrame, order, out: Path) -> None:
    """Mean cross-method Jaccard of the best (top-rank) pose per pair, + vs crystal."""
    cryst_fp = {(p, l): pose_fingerprint(g, "typed") for (p, l), g in crystal.groupby(["protein", "ligand"])}
    # best pose fingerprint per (method, pair) = lowest pose_rank
    best = inter.sort_values("pose_rank").groupby(["method", "protein", "ligand", "pose_name"], sort=False)
    fp_by = {}
    for (method, p, l, pose), g in best:
        fp_by.setdefault((method, p, l), pose_fingerprint(g, "typed"))  # first (top-rank) only
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
    sns.heatmap(mat, annot=True, fmt=".2f", cmap="mako",
                xticklabels=[pretty_method(m) for m in labels],
                yticklabels=[pretty_method(m) for m in labels],
                vmin=0, vmax=1, cbar_kws={"label": "mean Jaccard"}, ax=ax)
    ax.set_title("Interaction-fingerprint similarity (best pose per pair)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)


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
    fig.tight_layout(); fig.savefig(out_dir / "06_interaction_vs_chemotype.png", dpi=160); plt.close(fig)
    # correlation table
    cols = ["mw", "n_aromatic_rings", "hbd", "hba", "n_halogen", "formal_charge", "logp"]
    cols = [c for c in cols if c in df.columns]
    itypes = [t for t in INTERACTION_TYPES if t in df.columns and df[t].sum() > 0]
    corr = df[cols + itypes].corr().loc[cols, itypes]
    corr.round(3).to_csv(out_dir / "chem_descriptor_interaction_corr.csv")
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(itypes)), max(4, 0.6 * len(cols))))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", center=0, vmin=-1, vmax=1, ax=ax)
    ax.set_title("Correlation: ligand descriptor × interaction type")
    fig.tight_layout(); fig.savefig(out_dir / "06b_descriptor_corr.png", dpi=160); plt.close(fig)


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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
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
    fig.tight_layout(); fig.savefig(out, dpi=160); plt.close(fig)
    tab.to_csv(out.with_suffix(".csv"))


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
    ap.add_argument("--best-equibind-only", action="store_true", default=None,
                    help="Keep only the single best-performing EquiBind variant "
                         "(highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%, "
                         "read from --oracle-summary) in all charts/CSVs, relabelled "
                         "'EquiBind*'. AutoDock/DiffDock/crystal are unaffected. "
                         "Overrides config 'best_equibind_only'.")
    ap.add_argument("--best-diffdock-only", action="store_true", default=None,
                    help="Keep only the single best-performing DiffDock optimizer "
                         "variant (highest PB-Valid AND RMSD ≤ 2 Å = "
                         "oracle_pb_valid_and_rmsd2_%%, read from --oracle-summary) in all "
                         "charts/CSVs, relabelled 'DiffDock*'. AutoDock/EquiBind/crystal are "
                         "unaffected. "
                         "Overrides config 'best_diffdock_only'.")
    ap.add_argument("--oracle-summary", type=Path, default=None,
                    help="oracle_summary.csv from posebusters_pose_comparison.py, used "
                         "to pick the best EquiBind variant for --best-equibind-only. "
                         "Overrides config 'oracle_summary' "
                         f"(default: {DEFAULT_ORACLE_SUMMARY}).")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Restrict the report to the '<PDBID>_<LIG>' complex ids "
                         "listed in this file (one per line; '#' comments ok). "
                         "Overrides config 'ids_file'.")
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
    best_diffdock_only = (args.best_diffdock_only if args.best_diffdock_only is not None
                          else (getattr(cfg, "best_diffdock_only", False) if cfg else False))
    oracle_summary = (args.oracle_summary or (cfg.oracle_summary if cfg else None)
                      or DEFAULT_ORACLE_SUMMARY)
    ids_file = args.ids_file or (cfg.ids_file if cfg else None)

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

    # Optionally restrict every chart/CSV to the single best EquiBind variant.
    if best_equibind_only:
        summary, inter, best_eq = select_best_equibind(summary, inter, oracle_summary)
        if best_eq:
            _LABEL_OVERRIDES[best_eq] = "EquiBind*"
            print(f"best-equibind-only: '{best_eq}' is the top EquiBind variant by "
                  "PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'EquiBind*').")

    # Optionally restrict every chart/CSV to the single best DiffDock variant.
    if best_diffdock_only:
        summary, inter, best_dd = select_best_diffdock(summary, inter, oracle_summary)
        if best_dd:
            _LABEL_OVERRIDES[best_dd] = "DiffDock*"
            print(f"best-diffdock-only: '{best_dd}' is the top DiffDock variant by "
                  "PB-Valid AND RMSD ≤ 2 Å — keeping only it (shown as 'DiffDock*').")

    order = ordered_methods(summary["method"].unique())
    print(f"Loaded {len(summary)} poses, {len(inter)} interaction rows, "
          f"{len(crystal)} crystal rows. Methods: {order}")

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
        plot_fingerprint_similarity(inter, crystal, order, out_dir / "05_fingerprint_similarity.png")
    else:
        print("  [skip] native recovery / similarity — no crystal_interactions.csv")

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
