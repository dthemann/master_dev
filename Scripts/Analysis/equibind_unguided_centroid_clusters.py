#!/usr/bin/env python3
"""Cluster EquiBind *unguided* pose centroids to test pocket diversity.

EquiBind is a one-shot regressor: given a (ligand-graph, receptor-graph) pair it
predicts a single binding region, driven almost entirely by the receptor. The
pipeline diversifies "unguided" poses by feeding different RDKit conformers /
per-seed random vectors — but that mostly jitters the pose *within* EquiBind's
preferred site rather than exploring different pockets.

This script quantifies that directly. For each receptor it pulls the unguided
pose centroids from ``pipeline_summary.json`` (deduped to the raw EquiBind
centroid — clampOFF + raw variant when present), clusters them by 3D distance,
and reports how many *distinct* spatial sites the unguided poses actually cover.

Run in the analysis env (e.g. conda env ``vina``):

    python Scripts/Analysis/equibind_unguided_centroid_clusters.py \
        --summaries-dir Dockings/Benchmark_Equibind \
        --threshold 8.0 --out-dir pocket_results/unguided_clusters

A result of "most receptors collapse to 1 cluster" is evidence that extra
conformers buy little pocket diversity (see the pocket_comparison_report and the
EquiBind PoseBusters failure note).
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
# Shared (A),(B),(C)… panel labeller for multi-panel figures.
from pocket_comparison_report import _label_panels  # noqa: E402


# ── Data extraction ────────────────────────────────────────────────────────

def _raw_unguided_centroids(summary_path: Path) -> List[Tuple[float, float, float]]:
    """One raw-EquiBind centroid per unguided pose_num in a pipeline_summary.json.

    Prefers the clampOFF + raw variant (the un-clamped, un-refined EquiBind
    centroid). Falls back to any unguided variant if those tags are absent
    (legacy single-variant runs).
    """
    try:
        data = json.loads(summary_path.read_text())
    except Exception:
        return []
    results = data.get("all_results", [])
    by_pose: Dict[int, Tuple[int, Tuple[float, float, float]]] = {}
    for r in results:
        if r.get("mode") != "unguided" or not r.get("success"):
            continue
        c = r.get("pose_centroid")
        if not c:
            continue
        # priority: prefer clampOFF (raw centroid) and refine 'raw'
        clamp = r.get("clamp_variant")
        refine = r.get("refine_variant")
        prio = 0
        if clamp in (None, "clampOFF"):
            prio += 1
        if refine in (None, "raw"):
            prio += 1
        pose = r.get("pose_num", 1)
        best = by_pose.get(pose)
        if best is None or prio > best[0]:
            by_pose[pose] = (prio, (float(c[0]), float(c[1]), float(c[2])))
    return [v[1] for v in by_pose.values()]


# ── Clustering ─────────────────────────────────────────────────────────────

def _single_linkage_clusters(points: np.ndarray, threshold: float) -> np.ndarray:
    """Connected-component clustering: points within ``threshold`` (single
    linkage chains) share a cluster. Returns integer labels."""
    n = len(points)
    if n == 1:
        return np.zeros(1, dtype=int)
    from sklearn.cluster import AgglomerativeClustering
    model = AgglomerativeClustering(
        n_clusters=None, distance_threshold=threshold, linkage="single")
    return model.fit_predict(points)


def _spread(points: np.ndarray) -> Tuple[float, float]:
    """(max pairwise distance, radius of gyration) for a centroid set."""
    if len(points) < 2:
        return 0.0, 0.0
    from scipy.spatial.distance import pdist
    max_pair = float(pdist(points).max())
    rg = float(np.sqrt(((points - points.mean(axis=0)) ** 2).sum(axis=1).mean()))
    return max_pair, rg


# ── Per-receptor analysis ──────────────────────────────────────────────────

def analyze_receptor(summary_path: Path, threshold: float) -> Optional[dict]:
    cents = _raw_unguided_centroids(summary_path)
    if not cents:
        return None
    pts = np.asarray(cents, dtype=np.float64)
    labels = _single_linkage_clusters(pts, threshold)
    counts = Counter(labels.tolist())
    n_clusters = len(counts)
    largest = max(counts.values())
    max_pair, rg = _spread(pts)
    # receptor name from the directory: <id>/pipeline_summary.json
    protein = summary_path.parent.name
    return {
        "protein": protein,
        "n_poses": int(len(pts)),
        "n_clusters": int(n_clusters),
        "largest_cluster_frac": round(largest / len(pts), 4),
        "max_pairwise_dist": round(max_pair, 3),
        "radius_gyration": round(rg, 3),
        "collapsed": bool(n_clusters == 1),
    }


# ── Reporting ──────────────────────────────────────────────────────────────

def _save_plots(rows: List[dict], out_dir: Path, threshold: float) -> Optional[Path]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    n_clusters = [r["n_clusters"] for r in rows]
    largest = [r["largest_cluster_frac"] for r in rows]
    spread = [r["max_pairwise_dist"] for r in rows]

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    maxc = max(n_clusters) if n_clusters else 1
    axes[0].hist(n_clusters, bins=range(1, maxc + 2), align="left",
                 color="#4C72B0", edgecolor="white", rwidth=0.85)
    axes[0].set_xlabel(f"distinct unguided sites (<{threshold:g} Å single-linkage)")
    axes[0].set_ylabel("receptors")
    axes[0].set_title("Distinct sites per receptor")
    axes[0].set_xticks(range(1, maxc + 1))

    axes[1].hist(largest, bins=20, color="#55A868", edgecolor="white")
    axes[1].set_xlabel("largest-cluster fraction of unguided poses")
    axes[1].set_ylabel("receptors")
    axes[1].set_title("Pose concentration in the dominant site")

    axes[2].hist(spread, bins=30, color="#C44E52", edgecolor="white")
    axes[2].set_xlabel("max pairwise centroid distance (Å)")
    axes[2].set_ylabel("receptors")
    axes[2].set_title("Spatial spread of unguided poses")

    _label_panels(axes)
    fig.suptitle("EquiBind unguided pose diversity "
                 f"(n={len(rows)} receptors)", y=1.02, fontsize=13)
    fig.tight_layout()
    out = out_dir / "unguided_centroid_clusters.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summaries-dir", default="Dockings/Benchmark_Equibind",
                    help="Directory searched recursively for pipeline_summary.json "
                         "(default: Dockings/Benchmark_Equibind).")
    ap.add_argument("--glob", default="*/pipeline_summary.json",
                    help="Glob (relative to --summaries-dir) for the summaries.")
    ap.add_argument("--threshold", type=float, default=8.0,
                    help="Single-linkage distance (Å) below which two centroids "
                         "count as the same site (default: 8.0 = pocket_match_threshold).")
    ap.add_argument("--out-dir", default="pocket_results/unguided_clusters")
    ap.add_argument("--ids-file", default=None,
                    help="Optional text file of '<PDBID>_<CCD>' ids (one per line, "
                         "'#' comments allowed) restricting analysis to those "
                         "receptors (matched against each summary's parent dir name).")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    sdir = Path(args.summaries_dir)
    summaries = sorted(sdir.glob(args.glob))
    if not summaries:
        # fall back to a recursive search
        summaries = sorted(sdir.rglob("pipeline_summary.json"))
    if not summaries:
        print(f"No pipeline_summary.json found under {sdir}")
        return 1
    if args.ids_file:
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}
        kept = [s for s in summaries if s.parent.name in ids]
        print(f"Restricting to {len(kept)}/{len(summaries)} summaries listed in "
              f"{args.ids_file} ({len(ids)} ids).")
        summaries = kept
        if not summaries:
            print("No listed id matched a summary directory.")
            return 1

    rows: List[dict] = []
    for sp in summaries:
        row = analyze_receptor(sp, args.threshold)
        if row is not None:
            rows.append(row)
    if not rows:
        print("No unguided pose centroids found in any summary.")
        return 1

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    import pandas as pd
    df = pd.DataFrame(rows).sort_values("n_clusters", ascending=False)
    csv_path = out_dir / "per_receptor_unguided_clusters.csv"
    df.to_csv(csv_path, index=False)

    n = len(rows)
    dist = Counter(r["n_clusters"] for r in rows)
    collapsed = sum(1 for r in rows if r["collapsed"])
    mean_clusters = float(np.mean([r["n_clusters"] for r in rows]))
    mean_largest = float(np.mean([r["largest_cluster_frac"] for r in rows]))
    mean_spread = float(np.mean([r["max_pairwise_dist"] for r in rows]))
    mean_poses = float(np.mean([r["n_poses"] for r in rows]))

    print("=" * 72)
    print("EquiBind UNGUIDED pose diversity")
    print("=" * 72)
    print(f"  Receptors analysed:            {n}")
    print(f"  Mean unguided poses/receptor:  {mean_poses:.1f}")
    print(f"  Single-linkage threshold:      {args.threshold:g} Å")
    print(f"  Mean distinct sites/receptor:  {mean_clusters:.2f}")
    print(f"  Mean largest-cluster fraction: {mean_largest:.2f} "
          f"(1.0 = all poses in one site)")
    print(f"  Mean max pairwise spread:      {mean_spread:.2f} Å")
    print(f"  Receptors collapsing to 1 site: {collapsed}/{n} "
          f"({100*collapsed/n:.1f}%)")
    print("\n  Distribution of distinct sites per receptor:")
    for k in sorted(dist):
        bar = "█" * int(40 * dist[k] / n)
        print(f"    {k:>2} site(s): {dist[k]:>4}  {bar}")
    print("\n  Verdict:", _verdict(collapsed / n, mean_largest))
    print(f"\n  Per-receptor CSV: {csv_path}")

    if not args.no_plot:
        png = _save_plots(rows, out_dir, args.threshold)
        if png:
            print(f"  Figure:           {png}")

    summary = {
        "n_receptors": n, "mean_poses_per_receptor": round(mean_poses, 2),
        "threshold_A": args.threshold,
        "mean_distinct_sites": round(mean_clusters, 3),
        "mean_largest_cluster_frac": round(mean_largest, 3),
        "mean_max_pairwise_spread_A": round(mean_spread, 3),
        "receptors_collapsed_to_one_site": collapsed,
        "frac_collapsed": round(collapsed / n, 3),
        "distinct_site_distribution": {str(k): dist[k] for k in sorted(dist)},
    }
    (out_dir / "unguided_cluster_summary.json").write_text(json.dumps(summary, indent=2))
    return 0


def _verdict(frac_collapsed: float, mean_largest: float) -> str:
    if frac_collapsed >= 0.8 or mean_largest >= 0.9:
        return ("unguided poses overwhelmingly collapse into a single site — "
                "extra conformers/seeds buy little pocket diversity; use pocket "
                "cropping (guided) or DiffDock for multi-pocket coverage.")
    if frac_collapsed >= 0.5:
        return ("unguided poses mostly cluster in one site with occasional spread — "
                "diversity is limited; guided cropping / DiffDock cover pockets better.")
    return ("unguided poses spread across several sites — conformer/seed sampling "
            "does provide some pocket diversity here.")


if __name__ == "__main__":
    raise SystemExit(main())
