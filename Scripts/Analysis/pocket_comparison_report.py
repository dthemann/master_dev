#!/usr/bin/env python3
"""Compare fpocket vs p2rank pocket predictions across the benchmark receptors.

For every benchmark protein (``Data/PoseBuster Benchmark Set/<id>/<id>_protein.pdb``)
both fpocket and p2rank have already been run into ``pocket_results/`` by the
"Benchmark Set || Equibind" notebook cell. This script consumes those results and
answers:

  1. Overlap   — for the top-N pockets of each program, how many agree in 3D
                 (optimal one-to-one center matching within a distance cutoff),
                 per receptor and overall.
  2. Distance  — how close/far the top pockets are: top1↔top1 distance, nearest-
                 neighbour distances, and a rank×rank mean-distance matrix.
  3. Clustering— pool both programs' top-N centers and cluster by distance, so a
                 site is "consensus" (both programs), "fpocket-only" or
                 "p2rank-only" regardless of rank.
  4. Rank vs   — matched pairs keep their (fpocket_rank, p2rank_rank): a rank
     distance    crossover matrix shows e.g. top-rank fpocket aligning with a
                 3rd-rank p2rank pocket.
  5. Attributes— which attribute each program ranks by (Spearman rank↔attribute):
                 fpocket Score/Druggability/Volume/AlphaSpheres,
                 p2rank score/probability/sas_points/surf_atoms.
  6. Crystal   — how far each program's pockets sit from the experimentally
     accuracy    validated ligand, using two proximity metrics: DCC (pocket
                 centre → crystal ligand centroid) and DCA (pocket centre →
                 nearest ligand heavy atom, the fpocket / p2rank benchmark
                 metric). Reports the top-1 and best-of-top-N distance per
                 receptor and overall, with within-threshold "finds the true
                 site" rates for fpocket vs p2rank. Also resolves this by rank
                 (median + IQR distance of each ranked pocket) and a rank ×
                 distance-threshold hit-rate matrix (how often the rank-n pocket
                 lands within 3/5/7 Å, per --crystal-thresholds). Adds a Top-N
                 success curve, a success-vs-threshold sweep (Wilson 95% CI),
                 and a headline success-rate table (crystal_success_rate_table.csv).

Outputs CSVs, a summary JSON, and PNG figures under --out-dir. Every graph is
written as its own standalone PNG (no multi-panel grids).

Run in the analysis env (conda env ``vina`` has scipy/sklearn/matplotlib):

    python Scripts/Analysis/pocket_comparison_report.py --top-n 3 --match-thr 5.0
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import string
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ════════════════════════════════════════════════════════════════════════
# Parsing — self-contained so the script runs without the equibind package
# ════════════════════════════════════════════════════════════════════════

# fpocket descriptors we keep (info.txt key -> short column name).
_FPOCKET_ATTRS = {
    "Score": "score",
    "Druggability Score": "druggability",
    "Number of Alpha Spheres": "alpha_spheres",
    "Total SASA": "total_sasa",
    "Volume": "volume",
    "Hydrophobicity score": "hydrophobicity",
    "Polarity score": "polarity",
    "Flexibility": "flexibility",
}


def _load_ids(path: Path) -> set:
    """Load '<PDBID>_<CCD>' ids (one per line; '#' comments / blanks ignored)."""
    return {ln.strip() for ln in path.read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def _pdb_centroid(pdb: Path) -> Optional[Tuple[float, float, float]]:
    coords = []
    try:
        with open(pdb) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        coords.append((float(line[30:38]), float(line[38:46]),
                                       float(line[46:54])))
                    except (ValueError, IndexError):
                        continue
    except OSError:
        return None
    if not coords:
        return None
    c = np.asarray(coords).mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])


def _sdf_centroid(sdf: Path) -> Optional[Tuple[float, float, float]]:
    """Centroid of the first molecule in a V2000 SDF (the crystal ligand).

    Parses the atom block directly (columns x/y/z in 10-char fields) so the
    script stays dependency-free — no RDKit needed just to get a centre.
    """
    try:
        lines = sdf.read_text().splitlines()
    except OSError:
        return None
    if len(lines) < 4:
        return None
    try:
        natoms = int(lines[3][:3])          # V2000 counts line
    except (ValueError, IndexError):
        return None
    coords = []
    for ln in lines[4:4 + natoms]:
        try:
            coords.append((float(ln[0:10]), float(ln[10:20]), float(ln[20:30])))
        except (ValueError, IndexError):
            continue
    if not coords:
        return None
    c = np.asarray(coords).mean(axis=0)
    return float(c[0]), float(c[1]), float(c[2])


def _sdf_heavy_atoms(sdf: Path) -> Optional[np.ndarray]:
    """(N, 3) array of the crystal ligand's heavy-atom coordinates (V2000 SDF).

    Like ``_sdf_centroid`` but keeps every atom (dropping hydrogens) so the
    pocket-centre → nearest-ligand-atom distance (DCA) can be computed. DCA is
    the metric the fpocket / p2rank benchmarks report; unlike the centroid
    distance (DCC) it does not penalise a large pocket whose centre drifts from
    the ligand centre-of-mass. Returns None if the atom block is unreadable.
    """
    try:
        lines = sdf.read_text().splitlines()
    except OSError:
        return None
    if len(lines) < 4:
        return None
    try:
        natoms = int(lines[3][:3])          # V2000 counts line
    except (ValueError, IndexError):
        return None
    coords = []
    for ln in lines[4:4 + natoms]:
        try:
            x, y, z = float(ln[0:10]), float(ln[10:20]), float(ln[20:30])
        except (ValueError, IndexError):
            continue
        elem = ln[31:34].strip() if len(ln) > 31 else ""   # V2000 atom symbol
        if elem.upper() == "H":                             # heavy atoms only
            continue
        coords.append((x, y, z))
    if not coords:
        return None
    return np.asarray(coords, dtype=float)


def _dca(center, atoms: np.ndarray) -> float:
    """Distance from a pocket centre to the nearest ligand heavy atom (DCA)."""
    return float(np.linalg.norm(atoms - np.asarray(center, dtype=float),
                                axis=1).min())


def parse_fpocket(out_dir: Path) -> List[dict]:
    """Parse one ``*_out`` fpocket dir into ranked pocket dicts."""
    info_files = list(out_dir.glob("*_info.txt"))
    pockets_dir = out_dir / "pockets"
    if not info_files or not pockets_dir.exists():
        return []
    # attributes per pocket number
    attrs: Dict[int, dict] = defaultdict(dict)
    cur = None
    with open(info_files[0]) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            m = re.match(r"Pocket\s+(\d+)\s*:", s)
            if m:
                cur = int(m.group(1))
                continue
            if cur is None or ":" not in s:
                continue
            key, _, val = s.partition(":")
            key = key.strip()
            short = _FPOCKET_ATTRS.get(key)
            if short:
                try:
                    attrs[cur][short] = float(val.strip())
                except ValueError:
                    pass
    pockets: List[dict] = []
    for num, a in attrs.items():
        atm = pockets_dir / f"pocket{num}_atm.pdb"
        center = _pdb_centroid(atm) if atm.exists() else None
        if center is None:
            continue
        pockets.append({"program": "fpocket", "pocket_id": num, "center": center,
                        **a})
    # fpocket numbers pockets by descending Score; rank by score to be safe.
    pockets.sort(key=lambda p: p.get("score", 0.0), reverse=True)
    for i, p in enumerate(pockets, 1):
        p["rank"] = i
    return pockets


def parse_p2rank(pred_csv: Path) -> List[dict]:
    pockets: List[dict] = []
    try:
        with open(pred_csv) as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                r = {k.strip(): (v.strip() if isinstance(v, str) else v)
                     for k, v in row.items()}
                try:
                    pockets.append({
                        "program": "p2rank",
                        "rank": int(r["rank"]),
                        "pocket_id": int(r["rank"]),
                        "center": (float(r["center_x"]), float(r["center_y"]),
                                   float(r["center_z"])),
                        "score": float(r["score"]),
                        "probability": float(r["probability"]),
                        "sas_points": float(r.get("sas_points") or 0),
                        "surf_atoms": float(r.get("surf_atoms") or 0),
                    })
                except (KeyError, TypeError, ValueError):
                    continue
    except OSError:
        return []
    pockets.sort(key=lambda p: p["rank"])
    return pockets


# ════════════════════════════════════════════════════════════════════════
# Geometry helpers
# ════════════════════════════════════════════════════════════════════════

def _dist_matrix(a: List[dict], b: List[dict]) -> np.ndarray:
    pa = np.asarray([p["center"] for p in a])
    pb = np.asarray([p["center"] for p in b])
    return np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2)


def match_pockets(fp: List[dict], pr: List[dict], thr: float):
    """Optimal one-to-one center matching (Hungarian). Returns list of
    (fpocket, p2rank, distance, matched_bool)."""
    from scipy.optimize import linear_sum_assignment
    if not fp or not pr:
        return []
    D = _dist_matrix(fp, pr)
    ri, ci = linear_sum_assignment(D)
    out = []
    for i, j in zip(ri, ci):
        d = float(D[i, j])
        out.append((fp[i], pr[j], d, d <= thr))
    return out


# ════════════════════════════════════════════════════════════════════════
# Per-receptor analysis
# ════════════════════════════════════════════════════════════════════════

def analyze_receptor(name: str, fp_all: List[dict], pr_all: List[dict],
                     top_n: int, match_thr: float, consensus_thr: float) -> dict:
    fp = fp_all[:top_n]
    pr = pr_all[:top_n]
    n_eff = min(len(fp), len(pr))

    matches = match_pockets(fp, pr, match_thr)
    n_matched = sum(1 for *_, ok in matches if ok)
    overlap = (n_matched / n_eff) if n_eff else np.nan

    # top1 / nearest-neighbour distances
    top1_dist = np.nan
    if fp and pr:
        top1_dist = float(np.linalg.norm(
            np.subtract(fp[0]["center"], pr[0]["center"])))
    nn_fp_to_pr = np.nan
    if fp and pr:
        D = _dist_matrix(fp, pr)
        nn_fp_to_pr = float(D.min(axis=1).mean())

    # distance-based clustering of the pooled top-N (consensus detection)
    pooled = fp + pr
    consensus = fponly = pronly = 0
    if pooled:
        progs = [p["program"] for p in pooled]
        if len(pooled) == 1:
            labels = np.zeros(1, dtype=int)
        else:
            from sklearn.cluster import AgglomerativeClustering
            labels = AgglomerativeClustering(
                n_clusters=None, distance_threshold=consensus_thr,
                linkage="single").fit_predict(
                    np.asarray([p["center"] for p in pooled]))
        for lab in set(labels.tolist()):
            members = {progs[k] for k in range(len(progs)) if labels[k] == lab}
            if "fpocket" in members and "p2rank" in members:
                consensus += 1
            elif "fpocket" in members:
                fponly += 1
            else:
                pronly += 1
    total_clusters = consensus + fponly + pronly
    consensus_frac = (consensus / total_clusters) if total_clusters else np.nan

    return {
        "protein": name,
        "n_fpocket": len(fp_all), "n_p2rank": len(pr_all),
        "top_n_used": n_eff,
        "n_matched": n_matched, "overlap": round(overlap, 4) if n_eff else np.nan,
        "top1_dist": round(top1_dist, 3) if fp and pr else np.nan,
        "mean_nn_fp_to_pr": round(nn_fp_to_pr, 3) if fp and pr else np.nan,
        "consensus_sites": consensus, "fpocket_only_sites": fponly,
        "p2rank_only_sites": pronly,
        "consensus_frac": round(consensus_frac, 4) if total_clusters else np.nan,
        "_matches": matches,  # kept for aggregation, dropped before CSV
    }


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _label_panels(axes, fontsize: int = 13) -> None:
    """Annotate each visible subplot with (A), (B), (C)... in reading order
    (top-left → bottom-right) as a bold left-aligned title above the panel.

    This report now writes every graph as its own standalone PNG and no longer
    calls this itself, but it is still imported by the multi-panel figures in
    pose_cluster_crystal_pocket_report.py. Accepts a Matplotlib Axes array
    (2D/1D), a list of Axes, or a single Axes; hidden axes are skipped.
    """
    if hasattr(axes, "ravel"):
        flat = list(axes.ravel())
    elif isinstance(axes, (list, tuple)):
        flat = list(axes)
    else:
        flat = [axes]
    i = 0
    for ax in flat:
        if ax is None:
            continue
        try:
            if not ax.get_visible():
                continue
        except Exception:
            pass
        ax.set_title(f"({string.ascii_uppercase[i % 26]})", loc="left",
                     fontweight="bold", fontsize=fontsize)
        i += 1


# Consistent programme colours across every figure: fpocket = blue,
# p2rank = orange, the consensus of the two = purple. Matches the seaborn-muted
# scheme used by the other analysis figures in this repo (blue/orange/purple is a
# colourblind-safe triad).
_C_FP = "#4C72B0"       # fpocket
_C_PR = "#DD8452"       # p2rank
_C_CONS = "#8172B3"     # consensus / both programs
_C_NEG = "#C44E52"      # negative-correlation red
_C_ECDF = "#55A868"     # single-series ECDF green


def _ecdf_xy(vals):
    """Return (sorted values, cumulative fraction) for a step ECDF."""
    srt = np.sort(np.asarray(vals, dtype=float))
    y = np.arange(1, len(srt) + 1) / len(srt)
    return srt, y


def _rank_dist_arrays(crystal_rows, top_n, metric="dist"):
    """(n_receptors × top_n) arrays of each ranked pocket's distance to the
    crystal ligand, for fpocket and p2rank. ``metric`` selects the stored
    distance family: "dist" = centroid distance (DCC), "dca" = centre → nearest
    ligand atom (DCA). NaN where a program had < rank pockets (or that metric was
    unavailable). Shared by the summary (hit rates) and the figures."""
    fp = np.array([[r.get(f"fp_rank{k+1}_{metric}", np.nan) for k in range(top_n)]
                   for r in crystal_rows], dtype=float)
    pr = np.array([[r.get(f"pr_rank{k+1}_{metric}", np.nan) for k in range(top_n)]
                   for r in crystal_rows], dtype=float)
    return fp, pr


def _rank_hit_rates(arr, thresholds):
    """For each rank (row of ``arr``), the fraction of receptors whose pocket at
    that rank is within each threshold. Returns (rate_matrix, count_matrix),
    both (top_n × len(thresholds)); NaN/0 for a rank with no data."""
    top_n = arr.shape[1] if arr.ndim == 2 else 0
    rate = np.full((top_n, len(thresholds)), np.nan)
    cnt = np.zeros((top_n, len(thresholds)), dtype=int)
    for ri in range(top_n):
        valid = arr[:, ri][~np.isnan(arr[:, ri])]
        if not len(valid):
            continue
        for ti, thr in enumerate(thresholds):
            rate[ri, ti] = float((valid <= thr).mean())
            cnt[ri, ti] = int((valid <= thr).sum())
    return rate, cnt


def _wilson_ci(k, n, z=1.96):
    """Wilson score 95% CI for a binomial success rate k/n → (lo, hi)."""
    if not n:
        return (float("nan"), float("nan"))
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def _topn_success(arr, thr):
    """Cumulative Top-N success from a (n_receptors × top_n) distance array: for
    each N, the fraction of receptors whose BEST pocket among ranks 1..N is
    within ``thr``. Returns (frac[top_n], hits[top_n], n[top_n]); receptors with
    no data at all ranks 1..N are dropped from that N's denominator."""
    top_n = arr.shape[1] if arr.ndim == 2 else 0
    frac = np.full(top_n, np.nan)
    hits = np.zeros(top_n, dtype=int)
    ns = np.zeros(top_n, dtype=int)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for N in range(1, top_n + 1):
            best = np.nanmin(arr[:, :N], axis=1)       # nearest of ranks 1..N
            valid = best[~np.isnan(best)]
            if len(valid):
                hits[N - 1] = int((valid <= thr).sum())
                ns[N - 1] = len(valid)
                frac[N - 1] = hits[N - 1] / len(valid)
    return frac, hits, ns


def _save_one(fig, ax, out_dir, name, saved, grid=True):
    import matplotlib.pyplot as plt
    if grid and ax is not None:
        ax.grid(alpha=0.3); ax.set_axisbelow(True)
    fig.tight_layout()
    p = out_dir / name
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); saved.append(p)


def _donut(ax, values, colors, labels, centre_big, centre_small, pctfmt):
    """Draw a labelled donut on ``ax``; returns the wedge handles."""
    wedges, _, autotexts = ax.pie(
        values, colors=colors, startangle=90, counterclock=False,
        autopct=pctfmt, pctdistance=0.78,
        wedgeprops=dict(width=0.42, edgecolor="white", linewidth=2))
    for t in autotexts:
        t.set_color("white"); t.set_fontweight("bold"); t.set_fontsize(10)
    ax.text(0, 0.12, centre_big, ha="center", va="center",
            fontsize=19, fontweight="bold")
    ax.text(0, -0.18, centre_small, ha="center", va="center",
            fontsize=10, color="#52514e")
    ax.set_aspect("equal")
    return wedges


def _save_figures(rows, matched_pairs, rankmat_sum, rankmat_cnt,
                  attr_corr, crystal_rows, top_n, match_thr, consensus_thr,
                  crystal_thresholds, out_dir):
    """Write every graph as its own standalone PNG under ``out_dir``."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    saved = []

    overlaps = [r["overlap"] for r in rows if r["overlap"] == r["overlap"]]
    top1 = [r["top1_dist"] for r in rows if r["top1_dist"] == r["top1_dist"]]
    cons = sum(r["consensus_sites"] for r in rows)
    fponly = sum(r["fpocket_only_sites"] for r in rows)
    pronly = sum(r["p2rank_only_sites"] for r in rows)

    # ── overlap histogram ────────────────────────────────────────────────
    if overlaps:
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.hist(overlaps, bins=np.linspace(0, 1, 11), color=_C_FP,
                edgecolor="white")
        ax.axvline(np.mean(overlaps), color="k", ls="--",
                   label=f"mean = {np.mean(overlaps):.2f}")
        ax.set_title(f"Top-{top_n} pocket-centre overlap between fpocket and p2rank")
        ax.set_xlabel(f"Fraction of top-{top_n} pockets agreeing "
                      f"(centres ≤ {match_thr:g} Å apart)")
        ax.set_ylabel("Number of receptors")
        ax.legend()
        _save_one(fig, ax, out_dir, "pocket_overlap_hist.png", saved)

    # ── overlap ECDF ─────────────────────────────────────────────────────
    if overlaps:
        x, y = _ecdf_xy(overlaps)
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.step(x, y, where="post", color=_C_ECDF, lw=2)
        ax.axvline(np.mean(overlaps), color="k", ls="--", lw=0.9,
                   label=f"mean = {np.mean(overlaps):.2f}")
        ax.set_title(f"Per-receptor top-{top_n} overlap (cumulative)")
        ax.set_xlabel("Overlap fraction")
        ax.set_ylabel("Cumulative fraction of receptors")
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02); ax.legend(fontsize=9)
        _save_one(fig, ax, out_dir, "pocket_overlap_ecdf.png", saved)

    # ── top-1 ↔ top-1 centre-distance histogram ──────────────────────────
    if top1:
        fig, ax = plt.subplots(figsize=(6.4, 4.6))
        ax.hist(top1, bins=30, color=_C_PR, edgecolor="white")
        ax.axvline(np.median(top1), color="k", ls="--",
                   label=f"median = {np.median(top1):.1f} Å")
        ax.set_title("Distance between the top-ranked fpocket and p2rank pocket")
        ax.set_xlabel("Centre-to-centre distance (Å)")
        ax.set_ylabel("Number of receptors")
        ax.legend()
        _save_one(fig, ax, out_dir, "pocket_top1_distance_hist.png", saved)

    # ── distinct-site composition — DONUT (replaces the old bar chart) ────
    tot = cons + fponly + pronly
    if tot:
        fig, ax = plt.subplots(figsize=(6.6, 5.6))
        wedges = _donut(
            ax, [cons, fponly, pronly], [_C_CONS, _C_FP, _C_PR],
            None, f"{tot:,}", "distinct\nsites",
            pctfmt=lambda pct: f"{pct:.0f}%")
        ax.legend(wedges,
                  [f"Consensus (both) — {cons:,}",
                   f"fpocket only — {fponly:,}",
                   f"p2rank only — {pronly:,}"],
                  loc="upper center", bbox_to_anchor=(0.5, -0.01),
                  frameon=False, fontsize=9)
        ax.set_title("Distinct binding sites pooled from both programs\n"
                     f"(top-{top_n} each, sites clustered ≤ {consensus_thr:g} Å)")
        _save_one(fig, ax, out_dir, "pocket_site_composition.png", saved, grid=False)

    # ── rank × rank mean centre distance (heatmap) ───────────────────────
    with np.errstate(invalid="ignore"):
        mean_d = np.where(rankmat_cnt > 0,
                          rankmat_sum / np.maximum(rankmat_cnt, 1), np.nan)
    fig, ax = plt.subplots(figsize=(5.9, 5.1))
    im = ax.imshow(mean_d, cmap="viridis_r", origin="upper")
    ax.set_title("Mean pocket-centre distance by rank pair")
    ax.set_xlabel("p2rank rank"); ax.set_ylabel("fpocket rank")
    ax.set_xticks(range(top_n)); ax.set_xticklabels(range(1, top_n + 1))
    ax.set_yticks(range(top_n)); ax.set_yticklabels(range(1, top_n + 1))
    for i in range(top_n):
        for j in range(top_n):
            if mean_d[i, j] == mean_d[i, j]:
                ax.text(j, i, f"{mean_d[i, j]:.1f}", ha="center", va="center",
                        color="white", fontsize=9)
    fig.colorbar(im, ax=ax, label="Mean centre distance (Å)")
    _save_one(fig, ax, out_dir, "pocket_rank_distance_matrix.png", saved, grid=False)

    # ── rank crossover counts among matched pairs (heatmap) ──────────────
    cross = np.zeros((top_n, top_n), dtype=int)
    for fr, prr, d in matched_pairs:
        if fr <= top_n and prr <= top_n:
            cross[fr - 1, prr - 1] += 1
    fig, ax = plt.subplots(figsize=(5.9, 5.1))
    im = ax.imshow(cross, cmap="Blues", origin="upper")
    ax.set_title(f"Matched pocket pairs by rank (agree ≤ {match_thr:g} Å)")
    ax.set_xlabel("p2rank rank"); ax.set_ylabel("fpocket rank")
    ax.set_xticks(range(top_n)); ax.set_xticklabels(range(1, top_n + 1))
    ax.set_yticks(range(top_n)); ax.set_yticklabels(range(1, top_n + 1))
    cmax = cross.max() if cross.size else 0
    for i in range(top_n):
        for j in range(top_n):
            ax.text(j, i, cross[i, j], ha="center", va="center",
                    color="white" if cross[i, j] > 0.6 * cmax else "black",
                    fontsize=9)
    fig.colorbar(im, ax=ax, label="Number of matched pairs")
    _save_one(fig, ax, out_dir, "pocket_rank_crossover.png", saved, grid=False)

    # ── rank agreement: same rank vs crossover (donut) ───────────────────
    if matched_pairs:
        same = sum(1 for fr, prr, _ in matched_pairs if fr == prr)
        crossn = len(matched_pairs) - same
        fig, ax = plt.subplots(figsize=(6.4, 5.4))
        wedges = _donut(
            ax, [same, crossn], [_C_CONS, "#CCB974"],
            None, f"{len(matched_pairs)}", "matched\npairs",
            pctfmt=lambda pct: f"{pct:.1f}%")
        ax.legend(wedges,
                  [f"Same rank (i ↔ i) — {same}",
                   f"Rank crossover — {crossn}"],
                  loc="upper center", bbox_to_anchor=(0.5, -0.01),
                  frameon=False, fontsize=9)
        ax.set_title("Do matched pockets share the same rank in both programs?")
        _save_one(fig, ax, out_dir, "pocket_rank_agreement.png", saved, grid=False)

    # ── attribute-ranking correlations (one figure per program) ──────────
    for prog, pos_col, fname in (
            ("fpocket", _C_FP, "pocket_attr_ranking_fpocket.png"),
            ("p2rank", _C_PR, "pocket_attr_ranking_p2rank.png")):
        items = attr_corr.get(prog, {}) if attr_corr else {}
        if not items:
            continue
        order = sorted(items.items(), key=lambda kv: kv[1])
        labels = [a for a, _ in order]; vals = [v for _, v in order]
        colors = [_C_NEG if v < 0 else pos_col for v in vals]
        fig, ax = plt.subplots(figsize=(6.8, 0.55 * len(labels) + 2.0))
        ax.barh(labels, vals, color=colors)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_title(f"What {prog} ranks pockets by")
        ax.set_xlabel("Spearman ρ (pocket rank vs attribute)\n"
                      "← higher value gives a better rank   |   "
                      "higher value gives a worse rank →")
        _save_one(fig, ax, out_dir, fname, saved)

    # ── crystal accuracy: distance of each program's pockets from the
    #    experimentally-validated ligand centre ───────────────────────────
    if crystal_rows:
        fp_t1 = np.array([r["fp_top1_dist"] for r in crystal_rows])
        pr_t1 = np.array([r["pr_top1_dist"] for r in crystal_rows])
        fp_best = np.array([r["fp_best_dist"] for r in crystal_rows])
        pr_best = np.array([r["pr_best_dist"] for r in crystal_rows])

        for (a, b, name, title, xlab) in (
            (fp_t1, pr_t1, "crystal_top1_distance_ecdf.png",
             "How far the top-ranked pocket sits from the crystal ligand",
             "Top-1 pocket centre → crystal ligand centroid (Å)"),
            (fp_best, pr_best, "crystal_topN_best_distance_ecdf.png",
             f"Does either program place the true site in its top {top_n}?",
             f"Nearest of the top-{top_n} pocket centres "
             f"→ crystal ligand centroid (Å)")):
            fig, ax = plt.subplots(figsize=(6.6, 4.7))
            for vals, c, lab in ((a, _C_FP, "fpocket"), (b, _C_PR, "p2rank")):
                x, y = _ecdf_xy(vals)
                hit = float((vals <= match_thr).mean())
                ax.step(x, y, where="post", color=c, lw=2,
                        label=f"{lab}  (median {np.median(vals):.1f} Å, "
                              f"≤{match_thr:g} Å: {hit:.0%})")
            ax.axvline(match_thr, color="k", ls=":", lw=1.1,
                       label=f"match threshold = {match_thr:g} Å")
            ax.set_xlim(left=0); ax.set_ylim(0, 1.02)
            ax.set_title(title)
            ax.set_xlabel(xlab)
            ax.set_ylabel("Cumulative fraction of receptors")
            ax.legend(fontsize=8)
            _save_one(fig, ax, out_dir, name, saved)

        # paired top-1 scatter: which program's top pocket is closer per receptor
        lim = float(max(fp_t1.max(), pr_t1.max())) * 1.05 if len(fp_t1) else 1.0
        fig, ax = plt.subplots(figsize=(5.9, 5.7))
        ax.scatter(fp_t1, pr_t1, s=24, alpha=0.6, color=_C_CONS,
                   edgecolor="white", linewidth=0.4)
        ax.plot([0, lim], [0, lim], color="k", ls="--", lw=0.8,
                label="equal distance")
        ax.axvline(match_thr, color=_C_FP, ls=":", lw=1)
        ax.axhline(match_thr, color=_C_PR, ls=":", lw=1)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("fpocket top-1 → crystal (Å)")
        ax.set_ylabel("p2rank top-1 → crystal (Å)")
        ax.set_title("Which program's top pocket is closer to the crystal ligand?\n"
                     "(below the diagonal ⇒ p2rank closer; above ⇒ fpocket closer)")
        ax.legend(fontsize=8)
        _save_one(fig, ax, out_dir, "crystal_paired_top1_scatter.png", saved)

        # ── distance to the crystal ligand as a function of pocket rank ──
        fp_rank, pr_rank = _rank_dist_arrays(crystal_rows, top_n)
        ranks = np.arange(1, top_n + 1)
        fig, ax = plt.subplots(figsize=(6.6, 4.7))
        with np.errstate(all="ignore"):
            for lab, c, arr in (("fpocket", _C_FP, fp_rank), ("p2rank", _C_PR, pr_rank)):
                med = np.nanmedian(arr, axis=0)
                q1 = np.nanpercentile(arr, 25, axis=0)
                q3 = np.nanpercentile(arr, 75, axis=0)
                ax.plot(ranks, med, "-o", color=c, lw=2.2, label=f"{lab} (median)")
                ax.fill_between(ranks, q1, q3, color=c, alpha=0.15)
        ax.axhline(match_thr, color="k", ls=":", lw=1,
                   label=f"match threshold = {match_thr:g} Å")
        ax.set_xticks(ranks); ax.set_ylim(bottom=0)
        ax.set_xlabel("Pocket rank (1 = top-ranked pocket)")
        ax.set_ylabel("Distance to crystal ligand centroid (Å)")
        ax.set_title("How far each ranked pocket sits from the crystal ligand\n"
                     "(line = median across receptors, band = inter-quartile range)")
        ax.legend(fontsize=8)
        _save_one(fig, ax, out_dir, "crystal_rank_distance_by_program.png", saved)

        # ── rank × threshold hit-rate matrix, one per program ────────────
        for lab, arr, fname in (
                ("fpocket", fp_rank, "crystal_hit_matrix_fpocket.png"),
                ("p2rank", pr_rank, "crystal_hit_matrix_p2rank.png")):
            rate, cnt = _rank_hit_rates(arr, crystal_thresholds)
            fig, ax = plt.subplots(
                figsize=(1.25 * len(crystal_thresholds) + 2.4, 0.75 * top_n + 2.4))
            im = ax.imshow(rate, cmap="Greens", vmin=0, vmax=1, origin="upper",
                           aspect="auto")
            ax.set_xticks(range(len(crystal_thresholds)))
            ax.set_xticklabels([f"≤ {t:g} Å" for t in crystal_thresholds])
            ax.set_yticks(range(top_n))
            ax.set_yticklabels([f"rank {r}" for r in ranks])
            ax.set_xlabel("Distance to crystal ligand centroid")
            ax.set_ylabel("Pocket rank")
            ax.set_title(f"{lab}: how often the rank-n pocket reaches the crystal\n"
                         "(fraction of receptors within the distance; count below)")
            for ri in range(top_n):
                for ti in range(len(crystal_thresholds)):
                    if rate[ri, ti] == rate[ri, ti]:
                        ax.text(ti, ri, f"{rate[ri, ti]:.0%}\n(n={cnt[ri, ti]})",
                                ha="center", va="center", fontsize=9,
                                color="white" if rate[ri, ti] > 0.55 else "black")
            fig.colorbar(im, ax=ax, label="Fraction of receptors within distance")
            _save_one(fig, ax, out_dir, fname, saved, grid=False)

        # ── NEW: DCA (centre → nearest ligand atom) ECDFs ─────────────
        # Same view as the two DCC ECDFs above but with the benchmark-standard
        # DCA metric, which credits a pocket that touches the ligand even if its
        # centroid is offset. Only drawn when DCA was computed.
        has_dca = any("fp_top1_dca" in r for r in crystal_rows)
        if has_dca:
            fp_t1a = np.array([r.get("fp_top1_dca", np.nan) for r in crystal_rows])
            pr_t1a = np.array([r.get("pr_top1_dca", np.nan) for r in crystal_rows])
            fp_ba = np.array([r.get("fp_best_dca", np.nan) for r in crystal_rows])
            pr_ba = np.array([r.get("pr_best_dca", np.nan) for r in crystal_rows])
            for (a, b, name, title, xlab) in (
                (fp_t1a, pr_t1a, "crystal_dca_top1_distance_ecdf.png",
                 "DCA: top-ranked pocket → nearest crystal-ligand atom",
                 "Top-1 pocket centre → nearest ligand atom (Å)"),
                (fp_ba, pr_ba, "crystal_dca_topN_best_distance_ecdf.png",
                 f"DCA: does either program place the true site in its top {top_n}?",
                 f"Nearest of the top-{top_n} pocket centres "
                 f"→ nearest ligand atom (Å)")):
                fig, ax = plt.subplots(figsize=(6.6, 4.7))
                for vals, c, lab in ((a, _C_FP, "fpocket"), (b, _C_PR, "p2rank")):
                    v = vals[~np.isnan(vals)]
                    if not len(v):
                        continue
                    x, y = _ecdf_xy(v)
                    hit = float((v <= match_thr).mean())
                    ax.step(x, y, where="post", color=c, lw=2,
                            label=f"{lab}  (median {np.median(v):.1f} Å, "
                                  f"≤{match_thr:g} Å: {hit:.0%})")
                ax.axvline(match_thr, color="k", ls=":", lw=1.1,
                           label=f"match threshold = {match_thr:g} Å")
                ax.set_xlim(left=0); ax.set_ylim(0, 1.02)
                ax.set_title(title); ax.set_xlabel(xlab)
                ax.set_ylabel("Cumulative fraction of receptors")
                ax.legend(fontsize=8)
                _save_one(fig, ax, out_dir, name, saved)

        # ── NEW: Top-N success curve — true site found in the top-N pockets ──
        # For each N, the fraction of receptors with ≥1 top-N pocket within
        # match_thr of the crystal. Solid = DCA, dashed = DCC; Wilson 95% CI band.
        ranks = np.arange(1, top_n + 1)
        series = ([("DCA", "dca", "-")] if has_dca else []) + [("DCC", "dist", "--")]
        fig, ax = plt.subplots(figsize=(6.9, 4.9))
        for metric_lab, metric_key, ls in series:
            fp_a, pr_a = _rank_dist_arrays(crystal_rows, top_n, metric_key)
            for prog, c, arr in (("fpocket", _C_FP, fp_a), ("p2rank", _C_PR, pr_a)):
                frac, hits, ns = _topn_success(arr, match_thr)
                ax.plot(ranks, frac, ls, color=c, marker="o", lw=2,
                        label=f"{prog} · {metric_lab}")
                los = [_wilson_ci(int(hits[i]), int(ns[i]))[0] for i in range(top_n)]
                his = [_wilson_ci(int(hits[i]), int(ns[i]))[1] for i in range(top_n)]
                ax.fill_between(ranks, los, his, color=c, alpha=0.08)
        ax.set_ylim(0, 1.02); ax.set_xticks(ranks)
        ax.set_xlabel("Number of top pockets considered (N)")
        ax.set_ylabel(f"Receptors with the true site found (≤ {match_thr:g} Å)")
        ax.set_title("Top-N success rate — is the crystal site among the top N pockets?\n"
                     "(solid = DCA nearest-atom, dashed = DCC centroid; band = Wilson 95% CI)")
        ax.legend(fontsize=8, ncol=2)
        _save_one(fig, ax, out_dir, "crystal_success_vs_topn.png", saved)

        # ── NEW: success vs distance-threshold sweep (best of top-N) ─────
        # Continuous form of the hit matrix: fraction of receptors whose best
        # top-N pocket is within a sweeping cutoff. Solid = DCA (+ Wilson band),
        # dashed = DCC; dotted vertical = the match threshold.
        thr_grid = np.linspace(0, 15, 61)
        sweep = ([("DCA", "dca", "-", True)] if has_dca else []) + \
                [("DCC", "dist", "--", False)]
        fig, ax = plt.subplots(figsize=(6.9, 4.9))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            for metric_lab, metric_key, ls, band in sweep:
                fp_a, pr_a = _rank_dist_arrays(crystal_rows, top_n, metric_key)
                for prog, c, arr in (("fpocket", _C_FP, fp_a), ("p2rank", _C_PR, pr_a)):
                    best = np.nanmin(arr, axis=1)
                    best = best[~np.isnan(best)]
                    n = len(best)
                    if not n:
                        continue
                    fr = np.array([(best <= t).mean() for t in thr_grid])
                    ax.plot(thr_grid, fr, ls, color=c, lw=2,
                            label=f"{prog} · {metric_lab}")
                    if band:
                        cis = [_wilson_ci(int((best <= t).sum()), n) for t in thr_grid]
                        ax.fill_between(thr_grid, [lo for lo, _ in cis],
                                        [hi for _, hi in cis], color=c, alpha=0.10)
        ax.axvline(match_thr, color="k", ls=":", lw=1.1,
                   label=f"match threshold = {match_thr:g} Å")
        ax.set_xlim(0, 15); ax.set_ylim(0, 1.02)
        ax.set_xlabel(f"Distance cutoff — best of top-{top_n} pocket → crystal (Å)")
        ax.set_ylabel("Fraction of receptors with the true site found")
        ax.set_title("How often the true site is recovered vs. how strict the cutoff is\n"
                     "(best of top-N; solid = DCA, dashed = DCC; band = Wilson 95% CI on DCA)")
        ax.legend(fontsize=8, ncol=2)
        _save_one(fig, ax, out_dir, "crystal_success_vs_threshold.png", saved)

    return saved


# ════════════════════════════════════════════════════════════════════════
# Driver
# ════════════════════════════════════════════════════════════════════════

def _attribute_rank_correlation(all_fp, all_pr) -> Dict[str, Dict[str, float]]:
    """Spearman rho between rank and each ranking attribute, pooled over all
    receptors. Negative rho => higher attribute -> better (lower) rank."""
    from scipy.stats import spearmanr
    out: Dict[str, Dict[str, float]] = {"fpocket": {}, "p2rank": {}}
    specs = {
        "fpocket": (all_fp, ["score", "druggability", "volume", "alpha_spheres",
                             "total_sasa", "hydrophobicity"]),
        "p2rank": (all_pr, ["score", "probability", "sas_points", "surf_atoms"]),
    }
    for prog, (pockets, attrs) in specs.items():
        ranks = np.asarray([p["rank"] for p in pockets], dtype=float)
        for a in attrs:
            vals = np.asarray([p.get(a, np.nan) for p in pockets], dtype=float)
            mask = ~np.isnan(vals)
            if mask.sum() < 3 or len(set(ranks[mask])) < 2:
                continue
            rho, _ = spearmanr(ranks[mask], vals[mask])
            if rho == rho:
                out[prog][a] = round(float(rho), 3)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-dir", default="Data/PoseBuster Benchmark Set")
    ap.add_argument("--fpocket-dir", default="pocket_results/fpocket_results")
    ap.add_argument("--p2rank-dir", default="pocket_results/p2rank_results")
    ap.add_argument("--out-dir", default="pocket_results/comparison")
    ap.add_argument("--top-n", type=int, default=3,
                    help="Top pockets per program to compare (default 3).")
    ap.add_argument("--match-thr", type=float, default=5.0,
                    help="Center distance (Å) for two pockets to 'agree' (default 5).")
    ap.add_argument("--consensus-thr", type=float, default=None,
                    help="Distance (Å) for pooled-site clustering "
                         "(default: --match-thr).")
    ap.add_argument("--crystal-thresholds", default="3,5,7",
                    help="Comma-separated distances (Å) for the rank × threshold "
                         "crystal hit-rate matrix (default: 3,5,7).")
    ap.add_argument("--ids-file", default=None,
                    help="Optional text file of '<PDBID>_<CCD>' ids (one per line, "
                         "'#' comments allowed) restricting the comparison to those "
                         "receptors (e.g. the official PoseBusters benchmark set). "
                         "Default: every on-disk benchmark receptor.")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    consensus_thr = args.match_thr if args.consensus_thr is None else args.consensus_thr
    crystal_thresholds = [float(x) for x in str(args.crystal_thresholds).split(",")
                          if x.strip()]

    bench = Path(args.benchmark_dir)
    fp_dir = Path(args.fpocket_dir)
    pr_dir = Path(args.p2rank_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    receptors = sorted(d.name for d in bench.iterdir()
                       if d.is_dir() and not d.name.startswith(("_", ".")))
    if not receptors:
        print(f"No benchmark receptors under {bench}")
        return 1
    if args.ids_file:
        ids = _load_ids(Path(args.ids_file))
        kept = [r for r in receptors if r in ids]
        print(f"Restricting to {len(kept)}/{len(receptors)} on-disk receptors listed "
              f"in {args.ids_file} ({len(ids)} ids).")
        receptors = kept
        if not receptors:
            print("No listed id matched an on-disk benchmark receptor.")
            return 1

    rows: List[dict] = []
    matched_pairs: List[Tuple[int, int, float]] = []   # (fp_rank, pr_rank, dist)
    crystal_rows: List[dict] = []                       # per-receptor crystal accuracy
    rankmat_sum = np.zeros((args.top_n, args.top_n))
    rankmat_cnt = np.zeros((args.top_n, args.top_n))
    all_fp: List[dict] = []
    all_pr: List[dict] = []
    n_missing = 0
    n_no_crystal = 0

    for name in receptors:
        stem = f"{name}_protein"
        out_fp = fp_dir / f"{stem}_out"
        csv_pr = pr_dir / f"{stem}.pdb_predictions.csv"
        fp_all = parse_fpocket(out_fp) if out_fp.exists() else []
        pr_all = parse_p2rank(csv_pr) if csv_pr.exists() else []
        if not fp_all or not pr_all:
            n_missing += 1
            continue
        all_fp.extend(fp_all[:args.top_n])
        all_pr.extend(pr_all[:args.top_n])

        r = analyze_receptor(name, fp_all, pr_all, args.top_n,
                             args.match_thr, consensus_thr)
        # rank matrices + matched pairs
        for fpk, prk, d, ok in r.pop("_matches"):
            if ok:
                matched_pairs.append((fpk["rank"], prk["rank"], d))
        D = _dist_matrix(fp_all[:args.top_n], pr_all[:args.top_n])
        for i in range(min(args.top_n, D.shape[0])):
            for j in range(min(args.top_n, D.shape[1])):
                rankmat_sum[i, j] += D[i, j]
                rankmat_cnt[i, j] += 1
        rows.append(r)

        # ── crystal (experimental) accuracy: how far each program's pockets
        #    sit from the true ligand centre ──────────────────────────────
        lig = bench / name / f"{name}_ligand.sdf"
        cc = _sdf_centroid(lig) if lig.exists() else None
        if cc is None:
            n_no_crystal += 1
            continue
        fp_top = fp_all[:args.top_n]
        pr_top = pr_all[:args.top_n]
        fp_ds = [float(np.linalg.norm(np.subtract(p["center"], cc))) for p in fp_top]
        pr_ds = [float(np.linalg.norm(np.subtract(p["center"], cc))) for p in pr_top]
        fp_bi = int(np.argmin(fp_ds))
        pr_bi = int(np.argmin(pr_ds))
        # DCA (pocket centre → nearest ligand heavy atom): the fpocket / p2rank
        # benchmark metric, computed alongside the centroid distance (DCC) above.
        # NaN everywhere if the ligand atom block cannot be read.
        atoms = _sdf_heavy_atoms(lig)
        if atoms is not None:
            fp_da = [_dca(p["center"], atoms) for p in fp_top]
            pr_da = [_dca(p["center"], atoms) for p in pr_top]
            fp_dbi = int(np.argmin(fp_da))
            pr_dbi = int(np.argmin(pr_da))
        else:
            fp_da, pr_da = [], []
            fp_dbi = pr_dbi = 0
        crow = {
            "protein": name,
            "crystal_x": round(cc[0], 3), "crystal_y": round(cc[1], 3),
            "crystal_z": round(cc[2], 3),
            "fp_top1_dist": round(fp_ds[0], 3),
            "pr_top1_dist": round(pr_ds[0], 3),
            "fp_best_dist": round(fp_ds[fp_bi], 3),
            "fp_best_rank": fp_top[fp_bi]["rank"],
            "pr_best_dist": round(pr_ds[pr_bi], 3),
            "pr_best_rank": pr_top[pr_bi]["rank"],
            "fp_top1_hit": bool(fp_ds[0] <= args.match_thr),
            "pr_top1_hit": bool(pr_ds[0] <= args.match_thr),
        }
        # DCA scalars (centre → nearest ligand atom), mirroring the DCC fields.
        if atoms is not None:
            crow.update({
                "fp_top1_dca": round(fp_da[0], 3),
                "pr_top1_dca": round(pr_da[0], 3),
                "fp_best_dca": round(fp_da[fp_dbi], 3),
                "fp_best_dca_rank": fp_top[fp_dbi]["rank"],
                "pr_best_dca": round(pr_da[pr_dbi], 3),
                "pr_best_dca_rank": pr_top[pr_dbi]["rank"],
                "fp_top1_dca_hit": bool(fp_da[0] <= args.match_thr),
                "pr_top1_dca_hit": bool(pr_da[0] <= args.match_thr),
            })
        # per-rank distances of each ranked pocket to the crystal ligand — both
        # centroid distance (DCC, *_dist) and nearest-atom distance (DCA, *_dca);
        # NaN-padded to top_n when a program returned fewer pockets.
        for k in range(args.top_n):
            crow[f"fp_rank{k+1}_dist"] = round(fp_ds[k], 3) if k < len(fp_ds) else np.nan
            crow[f"pr_rank{k+1}_dist"] = round(pr_ds[k], 3) if k < len(pr_ds) else np.nan
            crow[f"fp_rank{k+1}_dca"] = round(fp_da[k], 3) if k < len(fp_da) else np.nan
            crow[f"pr_rank{k+1}_dca"] = round(pr_da[k], 3) if k < len(pr_da) else np.nan
        crystal_rows.append(crow)

    if not rows:
        print("No receptor had both fpocket and p2rank results.")
        return 1

    import pandas as pd
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "per_receptor_pocket_comparison.csv", index=False)
    pd.DataFrame(matched_pairs, columns=["fpocket_rank", "p2rank_rank", "distance_A"]
                 ).to_csv(out_dir / "matched_pairs.csv", index=False)

    attr_corr = _attribute_rank_correlation(all_fp, all_pr)
    pd.DataFrame([{"program": prog, "attribute": a, "spearman_rank_vs_attr": v}
                  for prog, d in attr_corr.items() for a, v in d.items()]
                 ).to_csv(out_dir / "attribute_rank_correlation.csv", index=False)

    # ── Crystal-ligand accuracy CSV + summary ────────────────────────────
    crystal_summary = None
    if crystal_rows:
        cdf = pd.DataFrame(crystal_rows)
        cdf.to_csv(out_dir / "per_receptor_crystal_accuracy.csv", index=False)
        fpt1, prt1 = cdf["fp_top1_dist"], cdf["pr_top1_dist"]
        fpb, prb = cdf["fp_best_dist"], cdf["pr_best_dist"]
        crystal_summary = {
            "n_receptors": int(len(cdf)),
            "fpocket_top1_median_A": round(float(fpt1.median()), 3),
            "p2rank_top1_median_A": round(float(prt1.median()), 3),
            "fpocket_top1_mean_A": round(float(fpt1.mean()), 3),
            "p2rank_top1_mean_A": round(float(prt1.mean()), 3),
            "fpocket_top1_hit_frac": round(float((fpt1 <= args.match_thr).mean()), 3),
            "p2rank_top1_hit_frac": round(float((prt1 <= args.match_thr).mean()), 3),
            "fpocket_best_topN_median_A": round(float(fpb.median()), 3),
            "p2rank_best_topN_median_A": round(float(prb.median()), 3),
            "fpocket_best_topN_hit_frac": round(float((fpb <= args.match_thr).mean()), 3),
            "p2rank_best_topN_hit_frac": round(float((prb <= args.match_thr).mean()), 3),
        }
        # rank × threshold hit rate (the quantitative form of the matrix figures)
        fp_rank, pr_rank = _rank_dist_arrays(crystal_rows, args.top_n)
        crystal_summary["hit_matrix_thresholds_A"] = crystal_thresholds
        crystal_summary["rank_hit_rate"] = {}
        for prog, arr in (("fpocket", fp_rank), ("p2rank", pr_rank)):
            rate, _ = _rank_hit_rates(arr, crystal_thresholds)
            crystal_summary["rank_hit_rate"][prog] = {
                f"rank{ri+1}": {f"within_{t:g}A":
                                (round(float(rate[ri, ti]), 3)
                                 if rate[ri, ti] == rate[ri, ti] else None)
                                for ti, t in enumerate(crystal_thresholds)}
                for ri in range(args.top_n)}

        # ── DCA (centre → nearest ligand atom) medians/hit-fractions ─────
        # Mirrors the DCC block above with the standard pocket-detection metric.
        if "fp_top1_dca" in cdf.columns:
            fpt1d, prt1d = cdf["fp_top1_dca"].dropna(), cdf["pr_top1_dca"].dropna()
            fpbd, prbd = cdf["fp_best_dca"].dropna(), cdf["pr_best_dca"].dropna()
            crystal_summary.update({
                "fpocket_top1_dca_median_A": round(float(fpt1d.median()), 3),
                "p2rank_top1_dca_median_A": round(float(prt1d.median()), 3),
                "fpocket_top1_dca_hit_frac": round(float((fpt1d <= args.match_thr).mean()), 3),
                "p2rank_top1_dca_hit_frac": round(float((prt1d <= args.match_thr).mean()), 3),
                "fpocket_best_topN_dca_median_A": round(float(fpbd.median()), 3),
                "p2rank_best_topN_dca_median_A": round(float(prbd.median()), 3),
                "fpocket_best_topN_dca_hit_frac": round(float((fpbd <= args.match_thr).mean()), 3),
                "p2rank_best_topN_dca_hit_frac": round(float((prbd <= args.match_thr).mean()), 3),
            })

        # ── Headline Top-N success-rate table (DCC + DCA), Wilson 95% CI ──
        # For each metric × program × level (top-1, top-3, best-of-top-N), the
        # fraction of receptors with a pocket within --match-thr of the crystal.
        # Written to crystal_success_rate_table.csv and echoed into the summary.
        metric_arrs = {"DCC": _rank_dist_arrays(crystal_rows, args.top_n, "dist")}
        if "fp_top1_dca" in cdf.columns:
            metric_arrs["DCA"] = _rank_dist_arrays(crystal_rows, args.top_n, "dca")
        levels = {1: "top1"}
        if args.top_n >= 3:
            levels[3] = "top3"
        levels[args.top_n] = f"best_top{args.top_n}"
        success_table = []
        for metric, (fp_a, pr_a) in metric_arrs.items():
            for prog, arr in (("fpocket", fp_a), ("p2rank", pr_a)):
                frac, hits, ns = _topn_success(arr, args.match_thr)
                for N, lname in sorted(levels.items()):
                    if not 1 <= N <= args.top_n:
                        continue
                    k, nn = int(hits[N - 1]), int(ns[N - 1])
                    lo, hi = _wilson_ci(k, nn)
                    success_table.append({
                        "metric": metric, "program": prog, "level": lname,
                        "threshold_A": args.match_thr, "n": nn, "hits": k,
                        "success_frac": round(k / nn, 4) if nn else None,
                        "ci95_lo": round(lo, 4) if nn else None,
                        "ci95_hi": round(hi, 4) if nn else None,
                    })
        if success_table:
            pd.DataFrame(success_table).to_csv(
                out_dir / "crystal_success_rate_table.csv", index=False)
            crystal_summary["success_table"] = success_table

    # ── Aggregate stats ──────────────────────────────────────────────────
    ov = df["overlap"].dropna()
    t1 = df["top1_dist"].dropna()
    cf = df["consensus_frac"].dropna()
    cons = int(df["consensus_sites"].sum())
    fponly = int(df["fpocket_only_sites"].sum())
    pronly = int(df["p2rank_only_sites"].sum())
    tot_sites = cons + fponly + pronly
    same_rank = sum(1 for fr, prr, _ in matched_pairs if fr == prr)
    n_pairs = len(matched_pairs)

    print("=" * 74)
    print(f"fpocket vs p2rank — top-{args.top_n} comparison "
          f"(match ≤ {args.match_thr:g} Å)")
    print("=" * 74)
    print(f"  Receptors compared:     {len(rows)}  "
          f"(skipped {n_missing} missing one program's pockets)")
    print(f"  Mean top-{args.top_n} overlap:      {ov.mean():.2f}  "
          f"(median {ov.median():.2f})")
    print(f"  Receptors w/ full agreement (overlap=1): "
          f"{int((ov >= 0.999).sum())}/{len(ov)}")
    print(f"  Receptors w/ zero overlap:               "
          f"{int((ov <= 0.001).sum())}/{len(ov)}")
    print(f"  Top1↔Top1 distance:      median {t1.median():.2f} Å  "
          f"mean {t1.mean():.2f} Å")
    print(f"  Top1 agree (≤ {args.match_thr:g} Å):   "
          f"{int((t1 <= args.match_thr).sum())}/{len(t1)}")
    print(f"\n  Pooled sites (≤ {consensus_thr:g} Å clusters):  {tot_sites}")
    print(f"    consensus (both):  {cons:>5}  ({100*cons/tot_sites:.1f}%)")
    print(f"    fpocket-only:      {fponly:>5}  ({100*fponly/tot_sites:.1f}%)")
    print(f"    p2rank-only:       {pronly:>5}  ({100*pronly/tot_sites:.1f}%)")
    print(f"    mean consensus frac/receptor: {cf.mean():.2f}")
    if n_pairs:
        print(f"\n  Matched pairs:          {n_pairs}")
        print(f"    same rank (i↔i):     {same_rank} ({100*same_rank/n_pairs:.1f}%)")
        print(f"    rank crossover:       {n_pairs-same_rank} "
              f"({100*(n_pairs-same_rank)/n_pairs:.1f}%)")
    print("\n  How each program ranks (Spearman rank × attribute; "
          "negative ⇒ drives better rank):")
    for prog in ("fpocket", "p2rank"):
        items = attr_corr.get(prog, {})
        ordered = sorted(items.items(), key=lambda kv: kv[1])
        s = "  ".join(f"{a}={v:+.2f}" for a, v in ordered)
        print(f"    {prog:<8} {s}")

    if crystal_summary:
        cs, n = crystal_summary, crystal_summary["n_receptors"]
        print(f"\n  Distance from the experimentally-validated (crystal) ligand "
              f"(n={n} receptors):")
        print(f"    top-1 pocket → crystal:   fpocket median "
              f"{cs['fpocket_top1_median_A']:.2f} Å "
              f"(finds site ≤ {args.match_thr:g} Å: {cs['fpocket_top1_hit_frac']:.0%})"
              f"   p2rank median {cs['p2rank_top1_median_A']:.2f} Å "
              f"({cs['p2rank_top1_hit_frac']:.0%})")
        print(f"    best of top-{args.top_n} → crystal: fpocket median "
              f"{cs['fpocket_best_topN_median_A']:.2f} Å "
              f"({cs['fpocket_best_topN_hit_frac']:.0%})"
              f"   p2rank median {cs['p2rank_best_topN_median_A']:.2f} Å "
              f"({cs['p2rank_best_topN_hit_frac']:.0%})")
        if "fpocket_top1_dca_median_A" in cs:
            print(f"    [DCA centre→nearest atom] top-1: fpocket median "
                  f"{cs['fpocket_top1_dca_median_A']:.2f} Å "
                  f"(≤ {args.match_thr:g} Å: {cs['fpocket_top1_dca_hit_frac']:.0%})"
                  f"   p2rank median {cs['p2rank_top1_dca_median_A']:.2f} Å "
                  f"({cs['p2rank_top1_dca_hit_frac']:.0%})")
            print(f"    [DCA centre→nearest atom] best of top-{args.top_n}: fpocket median "
                  f"{cs['fpocket_best_topN_dca_median_A']:.2f} Å "
                  f"({cs['fpocket_best_topN_dca_hit_frac']:.0%})"
                  f"   p2rank median {cs['p2rank_best_topN_dca_median_A']:.2f} Å "
                  f"({cs['p2rank_best_topN_dca_hit_frac']:.0%})")
        if n_no_crystal:
            print(f"    ({n_no_crystal} compared receptors had no readable "
                  f"crystal ligand SDF)")

    summary = {
        "n_receptors": len(rows), "n_skipped_missing": n_missing,
        "top_n": args.top_n, "match_thr_A": args.match_thr,
        "consensus_thr_A": consensus_thr,
        "overlap_mean": round(float(ov.mean()), 3),
        "overlap_median": round(float(ov.median()), 3),
        "n_full_agreement": int((ov >= 0.999).sum()),
        "n_zero_overlap": int((ov <= 0.001).sum()),
        "top1_dist_median_A": round(float(t1.median()), 3),
        "top1_dist_mean_A": round(float(t1.mean()), 3),
        "top1_agree_count": int((t1 <= args.match_thr).sum()),
        "pooled_sites": tot_sites,
        "consensus_sites": cons, "fpocket_only_sites": fponly,
        "p2rank_only_sites": pronly,
        "consensus_frac_overall": round(cons / tot_sites, 3) if tot_sites else None,
        "mean_consensus_frac_per_receptor": round(float(cf.mean()), 3),
        "matched_pairs": n_pairs, "same_rank_pairs": same_rank,
        "rank_crossover_pairs": n_pairs - same_rank,
        "attribute_rank_correlation": attr_corr,
        "crystal_accuracy": crystal_summary,
    }
    (out_dir / "pocket_comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  CSVs + summary written to: {out_dir}/")

    if not args.no_plot:
        figs = _save_figures(rows, matched_pairs, rankmat_sum, rankmat_cnt,
                             attr_corr, crystal_rows, args.top_n, args.match_thr,
                             consensus_thr, crystal_thresholds, out_dir)
        for f in figs:
            print(f"  Figure: {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
