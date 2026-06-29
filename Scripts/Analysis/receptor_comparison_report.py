#!/usr/bin/env python3
"""Compare the Orai MD-snapshot receptors against the PoseBusters benchmark receptors.

The Orai cross-docking targets are four MD snapshots of one hexameric Ca²⁺ channel,
whereas the benchmark targets are 308 diverse, mostly-soluble cognate drug targets.
This report quantifies *how atypical* Orai is as a docking target on the only axes
that are meaningfully comparable across unrelated proteins:

  1. PROTEIN SIZE / COMPOSITION — heavy-atom count, number of chains (oligomeric
     state), residue count. (Fold / sequence are NOT comparable and are not shown.)
  2. POCKET PROPERTIES — fpocket + p2rank pockets, already computed for both: pocket
     count, top/maximum volume, top/maximum druggability, and the best pocket scores.
     This asks the key question: do Orai's predicted pockets resemble typical
     benchmark binding sites, or atypical channel cavities (large, low-druggability)?

Two data subtleties handled here:
  * The Orai PDBs store the six subunits in the **segid** column (MONA…MONF), not the
    chain column (all "M"); the benchmark notebook's `protein_features` keys residues
    on (chain, resnum) and so mis-counts Orai as 1 chain / 223 residues. This script
    counts subunits from chain-id OR segid OR residue-number resets → the correct
    6 chains / 1338 residues.
  * Benchmark PDBs carry no hydrogens, so total-atom counts are confounded; the size
    comparison uses **heavy atoms** only.

Benchmark descriptors are read from the cached RDKit feature table
(``ligand_protein_features.csv``); Orai descriptors are computed from the original
structures (default ``Data/Receptors/Original``). Pockets come from
``pocket_results/{fpocket_results,p2rank_results}`` for both.

Run in the analysis env (conda env ``vina`` — pandas + matplotlib):

    python Scripts/Analysis/receptor_comparison_report.py
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from pocket_comparison_report import parse_fpocket, parse_p2rank  # noqa: E402
try:
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                   # pragma: no cover
    def _label_panels(axes, fontsize=12):
        flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
        for i, ax in enumerate(flat):
            if ax is not None and getattr(ax, "get_visible", lambda: True)():
                ax.text(-0.08, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                        fontsize=fontsize, fontweight="bold", va="bottom", ha="right")

ORAI_COLORS = ["#C44E52", "#8C2D2F", "#E07B7D", "#9B1B30"]


# ════════════════════════════════════════════════════════════════════════
# Orai protein descriptors (robust to chain-id vs segid vs flattened numbering)
# ════════════════════════════════════════════════════════════════════════

def orai_protein_features(pdb: Path) -> Optional[dict]:
    lines = [ln for ln in Path(pdb).read_text(errors="ignore").splitlines()
             if ln.startswith(("ATOM", "HETATM"))]
    if not lines:
        return None
    chain_ids = {ln[21] for ln in lines if ln[21].strip()}
    seg_ids = {ln[72:76].strip() for ln in lines if ln[72:76].strip()}
    use_chain = len(chain_ids) > 1
    use_seg = (not use_chain) and len(seg_ids) > 1

    heavy = 0
    res = set()
    reset_sub, prev = 0, None
    subunits = set()
    for ln in lines:
        el = ln[76:78].strip()
        name = ln[12:16].strip()
        is_h = el == "H" or (el == "" and name[:1] == "H")
        if not is_h:
            heavy += 1
        if name != "CA":
            continue
        try:
            rn = int(ln[22:26])
        except ValueError:
            continue
        if use_chain:
            sub = ln[21]
        elif use_seg:
            sub = ln[72:76].strip()
        else:                                  # flattened — new subunit when resnum resets
            if prev is not None and rn < prev:
                reset_sub += 1
            sub = reset_sub
        prev = rn
        subunits.add(sub)
        res.add((sub, rn))
    return dict(prot_heavy_atoms=heavy, n_chains=len(subunits), n_residues=len(res),
                n_cofactor_types=0)


# ════════════════════════════════════════════════════════════════════════
# Pocket properties (fpocket + p2rank) -> per-receptor summary
# ════════════════════════════════════════════════════════════════════════

def pocket_summary(fp_out: Path, pr_csv: Path) -> dict:
    fp = parse_fpocket(fp_out) if fp_out.is_dir() else []
    pr = parse_p2rank(pr_csv) if pr_csv.exists() else []

    def g(plist, key, agg):
        vals = [float(p[key]) for p in plist if p.get(key) is not None]
        return float(agg(vals)) if vals else np.nan

    out = {
        "n_fpocket": len(fp), "n_p2rank": len(pr),
        "top1_volume": (fp[0].get("volume") if fp else np.nan),
        "max_volume": g(fp, "volume", max),
        "top1_druggability": (fp[0].get("druggability") if fp else np.nan),
        "max_druggability": g(fp, "druggability", max),
        "top1_fpocket_score": (fp[0].get("score") if fp else np.nan),
        "top_p2rank_score": (pr[0].get("score") if pr else np.nan),
    }
    return out


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _overlay(ax, bench, orai_vals, title, xlabel, log=False, bins=30):
    v = pd.to_numeric(pd.Series(bench), errors="coerce").dropna()
    v = v[v > 0] if log else v
    if v.empty:
        ax.set_visible(False); return
    edges = (np.logspace(np.log10(v.min()), np.log10(v.max()), bins) if log
             else np.linspace(v.min(), v.max(), bins))
    ax.hist(v, bins=edges, color="#bcbcbc", alpha=0.85, label="benchmark receptors")
    ov = [x for x in orai_vals if x == x]
    for i, x in enumerate(ov):
        ax.axvline(x, color=ORAI_COLORS[i % len(ORAI_COLORS)], lw=1.6,
                   label="Orai snapshots" if i == 0 else None)
    if log:
        ax.set_xscale("log")
    rep = np.median(ov) if ov else np.nan
    pct = float((v < rep).mean() * 100) if rep == rep else np.nan
    ax.set_title(f"{title}\nOrai ≈ {rep:.0f} ({pct:.0f}th percentile of benchmark)"
                 if rep == rep else title)
    ax.set_xlabel(xlabel); ax.set_ylabel("Number of benchmark receptors")
    ax.legend(fontsize=7)


def fig_descriptors(bench: pd.DataFrame, orai: pd.DataFrame, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    specs = [
        ("prot_heavy_atoms", "Heavy atoms", False),
        ("n_chains", "Number of chains (oligomeric state)", False),
        ("n_residues", "Number of residues", False),
        ("n_cofactor_types", "Number of cofactor types", False),
    ]
    for a, (col, lab, log) in zip(ax.ravel(), specs):
        if col in bench:
            _overlay(a, bench[col], orai[col].tolist() if col in orai else [], lab, lab, log=log)
    # (E) chemical-size scatter: residues vs heavy atoms
    ax[1, 1].scatter(bench["n_residues"], bench["prot_heavy_atoms"], s=14, color="#bcbcbc",
                     alpha=0.7, label="benchmark")
    ax[1, 1].scatter(orai["n_residues"], orai["prot_heavy_atoms"], s=90, marker="*",
                     color="#C44E52", edgecolors="k", label="Orai", zorder=5)
    ax[1, 1].set_xlabel("Number of residues"); ax[1, 1].set_ylabel("Heavy atoms")
    ax[1, 1].set_title("Receptor size: residues vs heavy atoms")
    ax[1, 1].legend(fontsize=7)
    # (F) oligomeric state vs size
    ax[1, 2].scatter(bench["n_chains"], bench["n_residues"], s=14, color="#bcbcbc",
                     alpha=0.7, label="benchmark")
    ax[1, 2].scatter(orai["n_chains"], orai["n_residues"], s=90, marker="*",
                     color="#C44E52", edgecolors="k", label="Orai", zorder=5)
    ax[1, 2].set_xlabel("Number of chains"); ax[1, 2].set_ylabel("Number of residues")
    ax[1, 2].set_title("Oligomeric state vs residue count")
    ax[1, 2].legend(fontsize=7)
    for a in ax.ravel():
        a.grid(alpha=0.2)
    _label_panels(ax)
    fig.suptitle("Orai receptors vs PoseBusters benchmark receptors — size & composition "
                 f"(n={len(bench)} benchmark, {len(orai)} Orai snapshots). "
                 "Heavy atoms used (benchmark PDBs carry no hydrogens).", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "receptor_descriptor_comparison.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def fig_pockets(bench: pd.DataFrame, orai: pd.DataFrame, out_dir: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(18, 10))
    specs = [
        ("n_fpocket", "Number of fpocket pockets", False),
        ("top1_volume", "Top fpocket pocket volume (Å³)", True),
        ("max_volume", "Largest pocket volume (Å³)", True),
        ("top1_druggability", "Top pocket druggability score", False),
        ("max_druggability", "Best druggability score (any pocket)", False),
        ("top_p2rank_score", "Top p2rank pocket score", False),
    ]
    for a, (col, lab, log) in zip(ax.ravel(), specs):
        if col in bench:
            _overlay(a, bench[col], orai[col].tolist() if col in orai else [], lab, lab, log=log)
    _label_panels(ax)
    fig.suptitle("Orai receptors vs PoseBusters benchmark receptors — predicted pocket "
                 f"properties (fpocket + p2rank; n={len(bench)} benchmark, {len(orai)} Orai). "
                 "Do Orai's pockets look like druggable benchmark sites?", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "receptor_pocket_comparison.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def _pct(bench: pd.Series, val: float) -> float:
    v = pd.to_numeric(bench, errors="coerce").dropna()
    return float((v < val).mean() * 100) if len(v) and val == val else np.nan


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--features-csv",
                    default="PoseBusters_Benchmark_Analysis/ligand_protein_features.csv")
    ap.add_argument("--orai-dir", default="Data/Receptors/Original")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--p2rank-dir", default="pocket_results/p2rank_results")
    ap.add_argument("--ids-file", default="Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt")
    ap.add_argument("--out-dir", default="posebusters_results/benchmark/receptor_comparison")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    fp_dir, pr_dir = Path(args.fpocket_dir), Path(args.p2rank_dir)

    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    # ── benchmark: descriptors (cache) + pockets ─────────────────────────
    bench = pd.read_csv(args.features_csv)
    if ids is not None:
        bench = bench[bench["entry"].isin(ids)].copy()
    print(f"Benchmark receptors: {len(bench)}  | gathering pockets ...")
    prows = []
    for e in bench["entry"]:
        prows.append({"entry": e, **pocket_summary(
            fp_dir / f"{e}_protein_out", pr_dir / f"{e}_protein.pdb_predictions.csv")})
    bench = bench.merge(pd.DataFrame(prows), on="entry", how="left")

    # ── Orai: descriptors (originals) + pockets (cleaned) ────────────────
    orows = []
    for pdb in sorted(Path(args.orai_dir).glob("Orai*.pdb")):
        frame = pdb.stem
        feat = orai_protein_features(pdb)
        if feat is None:
            continue
        pk = pocket_summary(fp_dir / f"{frame}_cleaned_out",
                            pr_dir / f"{frame}_cleaned.pdb_predictions.csv")
        orows.append({"frame": frame, **feat, **pk})
    orai = pd.DataFrame(orows)
    if orai.empty:
        print("No Orai receptors found.")
        return 1
    print(f"Orai snapshots: {len(orai)}")

    bench.to_csv(out_dir / "benchmark_receptor_properties.csv", index=False)
    orai.to_csv(out_dir / "orai_receptor_properties.csv", index=False)

    # ── console insights ─────────────────────────────────────────────────
    descr_cols = ["prot_heavy_atoms", "n_chains", "n_residues"]
    pock_cols = ["n_fpocket", "top1_volume", "max_volume", "top1_druggability",
                 "max_druggability", "top_p2rank_score"]
    print("\n" + "=" * 78)
    print("ORAI vs BENCHMARK RECEPTORS")
    print("=" * 78)
    print(f"  {'property':<24}{'Orai (median)':>16}{'bench median':>16}{'Orai pct':>11}")
    for c in descr_cols + pock_cols:
        if c not in bench:
            continue
        ov = pd.to_numeric(orai[c], errors="coerce").median()
        bm = pd.to_numeric(bench[c], errors="coerce").median()
        print(f"  {c:<24}{ov:>16.2f}{bm:>16.2f}{_pct(bench[c], ov):>10.0f}%")

    # headline interpretation
    od = pd.to_numeric(orai["max_druggability"], errors="coerce").median()
    print(f"\n  Orai best-pocket druggability {od:.2f} "
          f"(benchmark median {pd.to_numeric(bench['max_druggability'],errors='coerce').median():.2f}); "
          f"largest cavity {pd.to_numeric(orai['max_volume'],errors='coerce').median():.0f} Å³ "
          f"at the {_pct(bench['max_volume'], pd.to_numeric(orai['max_volume'],errors='coerce').median()):.0f}th pct.")
    print(f"\n  CSVs written to: {out_dir}/")

    summary = {
        "n_benchmark": int(len(bench)), "n_orai": int(len(orai)),
        "orai_descriptor_percentiles": {c: _pct(bench[c], pd.to_numeric(orai[c], errors="coerce").median())
                                        for c in descr_cols if c in bench},
        "orai_pocket_percentiles": {c: _pct(bench[c], pd.to_numeric(orai[c], errors="coerce").median())
                                    for c in pock_cols if c in bench},
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))

    if not args.no_plot:
        print(f"  Figure: {fig_descriptors(bench, orai, out_dir)}")
        print(f"  Figure: {fig_pockets(bench, orai, out_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
