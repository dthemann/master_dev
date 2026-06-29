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

Outputs CSVs, a summary JSON, and multi-panel PNG figures under --out-dir.

Run in the analysis env (conda env ``vina`` has scipy/sklearn/matplotlib):

    python Scripts/Analysis/pocket_comparison_report.py --top-n 3 --match-thr 5.0
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import string
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

    Shared by the multi-panel figures here and in
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


def _save_figures(rows, matched_pairs, rankmat_sum, rankmat_cnt,
                  attr_corr, top_n, match_thr, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return []
    saved = []

    # ── Figure 1: overlap + distances + consensus composition ────────────
    overlaps = [r["overlap"] for r in rows if r["overlap"] == r["overlap"]]
    top1 = [r["top1_dist"] for r in rows if r["top1_dist"] == r["top1_dist"]]
    cons = sum(r["consensus_sites"] for r in rows)
    fponly = sum(r["fpocket_only_sites"] for r in rows)
    pronly = sum(r["p2rank_only_sites"] for r in rows)

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    ax[0, 0].hist(overlaps, bins=np.linspace(0, 1, 11), color="#4C72B0",
                  edgecolor="white")
    ax[0, 0].axvline(np.mean(overlaps), color="k", ls="--",
                     label=f"mean={np.mean(overlaps):.2f}")
    ax[0, 0].set_title(f"Top-{top_n} center overlap (match ≤ {match_thr:g} Å)")
    ax[0, 0].set_xlabel("fraction of top pockets agreeing")
    ax[0, 0].set_ylabel("receptors"); ax[0, 0].legend()

    srt = np.sort(overlaps)
    ecdf = np.arange(1, len(srt) + 1) / len(srt)
    ax[0, 1].step(srt, ecdf, where="post", color="#55A868", lw=2)
    ax[0, 1].axvline(np.mean(overlaps), color="k", ls="--", lw=0.8,
                     label=f"mean={np.mean(overlaps):.2f}")
    ax[0, 1].set_title("Per-receptor overlap (ECDF)")
    ax[0, 1].set_xlabel("overlap fraction")
    ax[0, 1].set_ylabel("cumulative fraction of receptors")
    ax[0, 1].set_xlim(0, 1); ax[0, 1].set_ylim(0, 1.02)
    ax[0, 1].grid(alpha=0.3); ax[0, 1].set_axisbelow(True)
    ax[0, 1].legend(fontsize=8)

    ax[1, 0].hist(top1, bins=30, color="#C44E52", edgecolor="white")
    ax[1, 0].axvline(np.median(top1), color="k", ls="--",
                     label=f"median={np.median(top1):.1f} Å")
    ax[1, 0].set_title("Top-1 fpocket ↔ top-1 p2rank distance")
    ax[1, 0].set_xlabel("center distance (Å)")
    ax[1, 0].set_ylabel("receptors"); ax[1, 0].legend()

    tot = cons + fponly + pronly
    comp = [cons / tot, fponly / tot, pronly / tot] if tot else [0, 0, 0]
    ax[1, 1].bar(["consensus\n(both)", "fpocket\nonly", "p2rank\nonly"], comp,
                 color=["#8172B3", "#4C72B0", "#CCB974"])
    ax[1, 1].set_title(f"Pooled-site composition ({tot} sites, ≤ {match_thr:g} Å)")
    ax[1, 1].set_ylabel("fraction of sites")
    for i, v in enumerate(comp):
        ax[1, 1].text(i, v + 0.01, f"{v:.0%}", ha="center")

    _label_panels(ax)
    fig.suptitle(f"fpocket vs p2rank — overlap & distance (n={len(rows)} receptors)",
                 fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = out_dir / "pocket_overlap_distance.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); saved.append(p)

    # ── Figure 2: rank×rank mean distance + rank crossover counts ────────
    with np.errstate(invalid="ignore"):
        mean_d = np.where(rankmat_cnt > 0, rankmat_sum / np.maximum(rankmat_cnt, 1),
                          np.nan)
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    im0 = ax[0].imshow(mean_d, cmap="viridis_r", origin="upper")
    ax[0].set_title("Mean center distance: fpocket rank × p2rank rank")
    ax[0].set_xlabel("p2rank rank"); ax[0].set_ylabel("fpocket rank")
    ax[0].set_xticks(range(top_n)); ax[0].set_xticklabels(range(1, top_n + 1))
    ax[0].set_yticks(range(top_n)); ax[0].set_yticklabels(range(1, top_n + 1))
    for i in range(top_n):
        for j in range(top_n):
            if mean_d[i, j] == mean_d[i, j]:
                ax[0].text(j, i, f"{mean_d[i, j]:.1f}", ha="center", va="center",
                           color="white", fontsize=9)
    fig.colorbar(im0, ax=ax[0], label="Å")

    # crossover counts among matched pairs
    cross = np.zeros((top_n, top_n), dtype=int)
    for fr, prr, d in matched_pairs:
        if fr <= top_n and prr <= top_n:
            cross[fr - 1, prr - 1] += 1
    im1 = ax[1].imshow(cross, cmap="Blues", origin="upper")
    ax[1].set_title(f"Matched pairs (≤ {match_thr:g} Å): fpocket rank × p2rank rank")
    ax[1].set_xlabel("p2rank rank"); ax[1].set_ylabel("fpocket rank")
    ax[1].set_xticks(range(top_n)); ax[1].set_xticklabels(range(1, top_n + 1))
    ax[1].set_yticks(range(top_n)); ax[1].set_yticklabels(range(1, top_n + 1))
    for i in range(top_n):
        for j in range(top_n):
            ax[1].text(j, i, cross[i, j], ha="center", va="center",
                       color="black", fontsize=9)
    fig.colorbar(im1, ax=ax[1], label="pairs")
    _label_panels(ax)
    fig.suptitle("Rank agreement vs crossover (diagonal = same rank)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "pocket_rank_distance.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); saved.append(p)

    # ── Figure 3: ranking-attribute correlations ─────────────────────────
    if attr_corr:
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
        for k, prog in enumerate(("fpocket", "p2rank")):
            items = attr_corr.get(prog, {})
            if not items:
                continue
            labels = list(items.keys())
            vals = [items[a] for a in labels]
            colors = ["#C44E52" if v < 0 else "#4C72B0" for v in vals]
            ax[k].barh(labels, vals, color=colors)
            ax[k].axvline(0, color="k", lw=0.8)
            ax[k].set_xlim(-1, 1)
            ax[k].set_title(f"{prog}: Spearman(rank, attribute)")
            ax[k].set_xlabel("← better pocket has higher value | lower value →")
        _label_panels(ax)
        fig.suptitle("What each program ranks by "
                     "(negative ⇒ attribute drives a better/lower rank)", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        p = out_dir / "pocket_attribute_ranking.png"
        fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig); saved.append(p)
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
    ap.add_argument("--ids-file", default=None,
                    help="Optional text file of '<PDBID>_<CCD>' ids (one per line, "
                         "'#' comments allowed) restricting the comparison to those "
                         "receptors (e.g. the official PoseBusters benchmark set). "
                         "Default: every on-disk benchmark receptor.")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)
    consensus_thr = args.match_thr if args.consensus_thr is None else args.consensus_thr

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
    rankmat_sum = np.zeros((args.top_n, args.top_n))
    rankmat_cnt = np.zeros((args.top_n, args.top_n))
    all_fp: List[dict] = []
    all_pr: List[dict] = []
    n_missing = 0

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
    }
    (out_dir / "pocket_comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  CSVs + summary written to: {out_dir}/")

    if not args.no_plot:
        figs = _save_figures(rows, matched_pairs, rankmat_sum, rankmat_cnt,
                             attr_corr, args.top_n, args.match_thr, out_dir)
        for f in figs:
            print(f"  Figure: {f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
