#!/usr/bin/env python3
"""Compare the computational effort of docking the PoseBusters benchmark set
across AutoDock Vina, DiffDock, and EquiBind (best configuration).

For every benchmark receptor-ligand pair this reads the per-complex timing each
tool already wrote to disk, joins it with the per-pose PoseBusters verdicts, and
reports two complementary views of cost, per method:

  1. WALL-CLOCK (time-to-result): how long you actually wait. Comparable across
     tools as elapsed time, but NOT as identical hardware work (AutoDock saturates
     many CPU threads; DiffDock/EquiBind use a GPU).
  2. RESOURCE-SECONDS (hardware work): the actual compute consumed, split by the
     hardware that did it —
       * CPU-core-seconds  = serial-equivalent CPU work (parallelism removed)
       * GPU-seconds        = time a GPU was occupied
     CPU-core-seconds and GPU-seconds are DIFFERENT currencies — a GPU-second is
     not a CPU-core-second — so they are reported side by side, never summed.

Both views are given overall and per docked pose, for *all generated* poses and
for only the poses that *survive the PoseBusters tests* (``pb_valid``).

Timing sources (per complex):
  * AutoDock : Dockings/Benchmark/<id>/<prep>/docking_summary.json
               -> overall.total_time_seconds              (CPU wall; Vina, N threads)
               CPU-core-seconds = wall × N threads (--autodock-cpu, default 32);
               GPU-seconds = 0.
  * DiffDock : Dockings/Benchmark_DiffDock/<id>/docking_summary.json
               -> overall.total_time_seconds              (single-GPU wall)
               GPU-seconds = wall; CPU side not separately recorded (GPU-bound).
  * EquiBind : Dockings/Benchmark_Equibind/<id>/pipeline_summary.json
               -> global_timing.pipeline_wall_time_s       (wall)
               -> all_results[].{prep,dock,post,refine}_time_s  (per-pose work)
               GPU-seconds = Σ dock_time_s (model inference);
               CPU-core-seconds = Σ (prep+post+refine)_time_s (conformers/UFF/smina).

EquiBind "best configuration": the benchmark run produced the full clamp×refine×
pocket fan-out, but only the **unguided + smina-refined** variant gives the best
PoseBusters validity. We therefore (a) count only that variant's poses, (b) sum
*only those poses'* per-pose resource-seconds, and (c) attribute EquiBind's
wall-clock to the best config by the fraction of per-pose compute time those poses
consumed. The unscaled full-run wall-clock is kept in the CSV for transparency.

Run in the analysis env (conda env ``vina`` — pandas + matplotlib):

    python Scripts/Analysis/docking_effort_comparison.py \
        --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
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
try:                                                # reuse the shared panel-labeller
    from pocket_comparison_report import _label_panels  # noqa: E402
except Exception:                                   # pragma: no cover - fallback
    def _label_panels(axes, fontsize=12):
        flat = np.atleast_1d(np.asarray(axes, dtype=object)).ravel()
        for i, ax in enumerate(flat):
            if ax is not None and getattr(ax, "get_visible", lambda: True)():
                ax.text(-0.08, 1.04, f"({chr(97 + i)})", transform=ax.transAxes,
                        fontsize=fontsize, fontweight="bold", va="bottom", ha="right")

# Method label -> (pretty name, bar colour, device used)
METHODS = {
    "autodock":  ("AutoDock Vina", "#4C72B0", "CPU"),
    "diffdock":  ("DiffDock",      "#55A868", "GPU"),
    "equibind":  ("EquiBind (unguided+smina)", "#C44E52", "GPU+CPU"),
}
# how each method appears in the PoseBusters per_pose_metrics ``method`` column
PER_POSE_METHOD = {
    "autodock": "autodock",
    "diffdock": "diffdock",
    "equibind": "equibind_unguided_smina_clampOFF",
}
# EquiBind per-pose timing fields, classified by the hardware that does the work
EQ_GPU_FIELD = "dock_time_s"                                   # model inference (GPU)
EQ_CPU_FIELDS = ("prep_time_s", "post_time_s", "refine_time_s")  # conformers/UFF/smina (CPU)

# resource-bar colours (by hardware, shared across methods)
CPU_COLOR, GPU_COLOR = "#8C8C8C", "#E1A100"


# ════════════════════════════════════════════════════════════════════════
# Per-complex timing readers (one per tool)
# ════════════════════════════════════════════════════════════════════════

def _autodock_time(cid: str, root: Path, prep: str) -> Optional[float]:
    p = root / cid / prep / "docking_summary.json"
    if not p.exists():
        return None
    try:
        return float(json.loads(p.read_text())["overall"]["total_time_seconds"])
    except Exception:
        return None


def _diffdock_time(cid: str, root: Path) -> Optional[float]:
    p = root / cid / "docking_summary.json"
    if not p.exists():
        return None
    try:
        return float(json.loads(p.read_text())["overall"]["total_time_seconds"])
    except Exception:
        return None


def _equibind_time(cid: str, root: Path, eq_mode: str, eq_refine: str,
                   eq_clamp: str) -> Optional[tuple]:
    """Return (wall_full_s, wall_best_s, best_gpu_s, best_cpu_core_s) for one complex.

    wall_best is the full wall-clock scaled by the share of per-pose compute time
    the best-config (unguided+smina+clampOFF) poses consumed; best_gpu_s /
    best_cpu_core_s are the summed per-pose resource-seconds of those poses.
    """
    p = root / cid / "pipeline_summary.json"
    if not p.exists():
        return None
    try:
        d = json.loads(p.read_text())
    except Exception:
        return None
    wall_full = float(d.get("global_timing", {}).get("pipeline_wall_time_s", np.nan))
    if not np.isfinite(wall_full):
        return None
    rows = d.get("all_results") or []

    def gpu_s(r):
        return float(r.get(EQ_GPU_FIELD) or 0.0)

    def cpu_s(r):
        return sum(float(r.get(k) or 0.0) for k in EQ_CPU_FIELDS)

    total_c = sum(gpu_s(r) + cpu_s(r) for r in rows)
    best = [r for r in rows if r.get("mode") == eq_mode
            and r.get("refine_variant") == eq_refine
            and r.get("clamp_variant") == eq_clamp]
    best_gpu = sum(gpu_s(r) for r in best)
    best_cpu = sum(cpu_s(r) for r in best)
    best_c = best_gpu + best_cpu
    frac = (best_c / total_c) if total_c > 0 else np.nan
    wall_best = wall_full * frac if np.isfinite(frac) else np.nan
    return wall_full, wall_best, best_gpu, best_cpu


# ════════════════════════════════════════════════════════════════════════
# Aggregation
# ════════════════════════════════════════════════════════════════════════

def load_pose_counts(csv: Path, ids: Optional[set]) -> Dict[str, pd.DataFrame]:
    """Per method: a DataFrame indexed by complex id with generated + pb_valid counts."""
    df = pd.read_csv(csv, low_memory=False)
    if "pb_valid" in df.columns:
        df["pb_valid"] = df["pb_valid"].astype(str).str.lower().isin(("true", "1", "1.0"))
    else:
        df["pb_valid"] = False
    df["cid"] = df["protein"].astype(str)
    if ids is not None:
        df = df[df["cid"].isin(ids)].copy()
    out = {}
    for key, per_pose_name in PER_POSE_METHOD.items():
        sub = df[df["method"] == per_pose_name]
        g = sub.groupby("cid").agg(generated=("pb_valid", "size"),
                                   pb_valid=("pb_valid", "sum"))
        out[key] = g
    return out


def build_table(args) -> tuple:
    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    counts = load_pose_counts(Path(args.per_pose_csv), ids)
    ad_root = Path(args.autodock_dir)
    dd_root = Path(args.diffdock_dir)
    eb_root = Path(args.equibind_dir)

    rows: List[dict] = []
    for key in METHODS:
        for cid, crow in counts[key].iterrows():
            wall = wall_full = cpu_core_s = gpu_s = None
            if key == "autodock":
                wall = _autodock_time(cid, ad_root, args.autodock_prep)
                if wall is not None:
                    wall_full, cpu_core_s, gpu_s = wall, wall * args.autodock_cpu, 0.0
            elif key == "diffdock":
                wall = _diffdock_time(cid, dd_root)
                if wall is not None:
                    wall_full, cpu_core_s, gpu_s = wall, 0.0, wall
            else:
                t = _equibind_time(cid, eb_root, args.eq_mode, args.eq_refine, args.eq_clamp)
                if t:
                    wall_full, wall, gpu_s, cpu_core_s = t[0], t[1], t[2], t[3]
            if wall is None or not np.isfinite(wall):
                continue
            rows.append({
                "method": key, "cid": cid,
                "wall_s": float(wall), "wall_full_s": float(wall_full),
                "cpu_core_s": float(cpu_core_s), "gpu_s": float(gpu_s),
                "generated": int(crow["generated"]), "pb_valid": int(crow["pb_valid"]),
            })
    pc = pd.DataFrame(rows)

    summary_rows = []
    for key, (pretty, _c, device) in METHODS.items():
        m = pc[pc["method"] == key]
        if m.empty:
            continue
        tot_wall = float(m["wall_s"].sum())
        tot_cpu = float(m["cpu_core_s"].sum())
        tot_gpu = float(m["gpu_s"].sum())
        n_gen = int(m["generated"].sum())
        n_val = int(m["pb_valid"].sum())
        summary_rows.append({
            "method": key, "name": pretty, "device": device,
            "n_complexes": int(len(m)),
            "total_wall_h": round(tot_wall / 3600.0, 4),
            "total_wall_full_h": round(float(m["wall_full_s"].sum()) / 3600.0, 4),
            "total_cpu_core_h": round(tot_cpu / 3600.0, 4),
            "total_gpu_h": round(tot_gpu / 3600.0, 4),
            "poses_generated": n_gen,
            "poses_pb_valid": n_val,
            "pb_valid_rate_pct": round(100.0 * n_val / n_gen, 2) if n_gen else np.nan,
            # wall-clock per pose
            "wall_s_per_generated": round(tot_wall / n_gen, 4) if n_gen else np.nan,
            "wall_s_per_pb_valid": round(tot_wall / n_val, 4) if n_val else np.nan,
            # resource-seconds per pose
            "cpu_core_s_per_generated": round(tot_cpu / n_gen, 4) if n_gen else np.nan,
            "gpu_s_per_generated": round(tot_gpu / n_gen, 4) if n_gen else np.nan,
            "cpu_core_s_per_pb_valid": round(tot_cpu / n_val, 4) if n_val else np.nan,
            "gpu_s_per_pb_valid": round(tot_gpu / n_val, 4) if n_val else np.nan,
        })
    summ = pd.DataFrame(summary_rows)
    return pc, summ


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════

def _bars(ax, x, ys, colors, fmt, devices=None, ymax_mult=1.25):
    ax.bar(x, ys, color=colors)
    ymax = max([v for v in ys if v == v] + [0]) or 1
    for i, (xi, y) in enumerate(zip(x, ys)):
        if y == y:
            lab = fmt.format(y) + (f"\n[{devices[i]}]" if devices else "")
            ax.text(xi, y + 0.01 * ymax, lab, ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * ymax_mult)


def make_wall_figure(pc, summ, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    colors = [METHODS[k][1] for k in order]
    devices = [METHODS[k][2] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order))

    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    _bars(ax[0, 0], x, s["total_wall_h"].to_numpy(float), colors, "{:.2f} h", devices)
    ax[0, 0].set_title("Overall computational effort\n(total wall-clock to dock the benchmark)")
    ax[0, 0].set_ylabel("Total wall-clock time (hours)")

    data = [pc[pc["method"] == k]["wall_s"].to_numpy(float) for k in order]
    bp = ax[0, 1].boxplot(data, positions=x, widths=0.6, patch_artist=True,
                          showfliers=False, medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.75)
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_title("Per-complex docking wall-clock time\n(distribution over complexes)")
    ax[0, 1].set_ylabel("Wall-clock time per complex (seconds, log scale)")

    w = 0.38
    gen = s["poses_generated"].to_numpy(float)
    val = s["poses_pb_valid"].to_numpy(float)
    ax[0, 2].bar(x - w / 2, gen, w, color=colors, alpha=0.55, label="generated")
    ax[0, 2].bar(x + w / 2, val, w, color=colors, label="survive PoseBusters")
    for i in range(len(order)):
        rate = s["pb_valid_rate_pct"].to_numpy(float)[i]
        if rate == rate:
            ax[0, 2].text(x[i] + w / 2, val[i], f"{rate:.1f}%", ha="center",
                          va="bottom", fontsize=8)
    ax[0, 2].set_title("Poses generated vs. poses that survive\nthe PoseBusters tests")
    ax[0, 2].set_ylabel("Number of poses")
    ax[0, 2].legend(fontsize=8)

    _bars(ax[1, 0], x, s["wall_s_per_generated"].to_numpy(float), colors, "{:.2f} s")
    ax[1, 0].set_title("Timing per generated pose")
    ax[1, 0].set_ylabel("Wall-clock time per generated pose (seconds)")

    _bars(ax[1, 1], x, s["wall_s_per_pb_valid"].to_numpy(float), colors, "{:.2f} s")
    ax[1, 1].set_title("Timing per pose that survives the PoseBusters tests")
    ax[1, 1].set_ylabel("Wall-clock time per valid pose (seconds)")

    _bars(ax[1, 2], x, s["pb_valid_rate_pct"].to_numpy(float), colors, "{:.1f}%", ymax_mult=1.3)
    ax[1, 2].set_title("PoseBusters validity rate\n(fraction of generated poses that survive)")
    ax[1, 2].set_ylabel("Valid poses (percent of generated)")
    ax[1, 2].set_ylim(0, min(100, ax[1, 2].get_ylim()[1]))

    for a in ax.ravel():
        a.set_xticks(x); a.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        a.grid(axis="y", alpha=0.25)
    _label_panels(ax)
    n_cx = int(pc.groupby("method")["cid"].nunique().max()) if not pc.empty else 0
    fig.suptitle("Docking WALL-CLOCK (time-to-result) — AutoDock Vina vs DiffDock vs "
                 f"EquiBind (best config), PoseBuster benchmark (≈{n_cx} complexes/method).\n"
                 "AutoDock = CPU (multi-thread), DiffDock/EquiBind = single GPU — comparable "
                 "as elapsed time, not as identical-hardware work.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "docking_effort_comparison.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


def make_resource_figure(summ, out_dir):
    """CPU-core-seconds vs GPU-seconds — the true hardware work (different currencies)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order)); w = 0.38

    def grouped(ax, cpu, gpu, fmt):
        b1 = ax.bar(x - w / 2, cpu, w, color=CPU_COLOR, label="CPU-core-seconds")
        b2 = ax.bar(x + w / 2, gpu, w, color=GPU_COLOR, label="GPU-seconds")
        ymax = max(list(cpu) + list(gpu) + [0]) or 1
        for bars, vals in ((b1, cpu), (b2, gpu)):
            for rect, v in zip(bars, vals):
                if v == v:
                    ax.text(rect.get_x() + rect.get_width() / 2, v + 0.01 * ymax,
                            fmt.format(v), ha="center", va="bottom", fontsize=8, rotation=0)
        ax.set_ylim(0, ymax * 1.25)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25); ax.legend(fontsize=8)

    fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
    grouped(ax[0], s["total_cpu_core_h"].to_numpy(float), s["total_gpu_h"].to_numpy(float), "{:.1f}")
    ax[0].set_title("Overall hardware effort\n(total resource consumed)")
    ax[0].set_ylabel("Resource used (hours): CPU-core-hours | GPU-hours")

    grouped(ax[1], s["cpu_core_s_per_generated"].to_numpy(float),
            s["gpu_s_per_generated"].to_numpy(float), "{:.2f}")
    ax[1].set_title("Hardware effort per generated pose")
    ax[1].set_ylabel("Resource per generated pose (seconds)")

    grouped(ax[2], s["cpu_core_s_per_pb_valid"].to_numpy(float),
            s["gpu_s_per_pb_valid"].to_numpy(float), "{:.2f}")
    ax[2].set_title("Hardware effort per pose that survives\nthe PoseBusters tests")
    ax[2].set_ylabel("Resource per valid pose (seconds)")

    _label_panels(ax)
    fig.suptitle("Docking HARDWARE EFFORT — CPU-core-seconds vs GPU-seconds (the actual work, "
                 "parallelism removed).\nA GPU-second is NOT a CPU-core-second — the two bars are "
                 "different currencies and are never summed. AutoDock is CPU-only; DiffDock is "
                 "GPU-only; EquiBind uses GPU for inference and CPU for UFF/smina.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = out_dir / "docking_resource_effort.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/benchmark/dock_best_equi_top5/"
                            "pose_comparison_report/per_pose_metrics.csv")
    ap.add_argument("--autodock-dir", default="Dockings/Benchmark")
    ap.add_argument("--autodock-prep", default="meeko", choices=("meeko", "mgl_tools"))
    ap.add_argument("--autodock-cpu", type=int, default=32,
                    help="CPU threads Vina used per complex (for CPU-core-seconds). "
                         "Benchmark runs used 32.")
    ap.add_argument("--diffdock-dir", default="Dockings/Benchmark_DiffDock")
    ap.add_argument("--equibind-dir", default="Dockings/Benchmark_Equibind")
    ap.add_argument("--eq-mode", default="unguided")
    ap.add_argument("--eq-refine", default="smina")
    ap.add_argument("--eq-clamp", default="clampOFF")
    ap.add_argument("--ids-file", default="Data/PoseBuster Benchmark Set/"
                                          "posebusters_pdb_ccd_ids.txt")
    ap.add_argument("--out-dir", default="posebusters_results/benchmark/docking_effort")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pc, summ = build_table(args)
    if summ.empty:
        print("No data found — check the docking dirs and per-pose CSV path.")
        return 1

    pc.to_csv(out_dir / "per_complex_effort.csv", index=False)
    summ.to_csv(out_dir / "effort_summary.csv", index=False)

    # ── console report ───────────────────────────────────────────────────
    print("\n" + "=" * 100)
    print("DOCKING COMPUTATIONAL EFFORT — benchmark set "
          f"(EquiBind best config: {args.eq_mode}+{args.eq_refine}+{args.eq_clamp})")
    print("=" * 100)
    print("WALL-CLOCK (time-to-result):")
    print(f"  {'method':<26}{'dev':>8}{'cplx':>6}{'total h':>9}{'gen':>9}"
          f"{'valid':>8}{'valid%':>8}{'s/gen':>8}{'s/valid':>9}")
    for _, r in summ.iterrows():
        print(f"  {r['name']:<26}{r['device']:>8}{r['n_complexes']:>6}{r['total_wall_h']:>9.2f}"
              f"{r['poses_generated']:>9,}{r['poses_pb_valid']:>8,}{r['pb_valid_rate_pct']:>7.1f}%"
              f"{r['wall_s_per_generated']:>8.2f}{r['wall_s_per_pb_valid']:>9.2f}")
    print("\nHARDWARE EFFORT (resource-seconds; CPU-core-s and GPU-s are different currencies):")
    print(f"  {'method':<26}{'CPU core-h':>12}{'GPU-h':>9}"
          f"{'CPUs/valid':>12}{'GPUs/valid':>12}")
    for _, r in summ.iterrows():
        print(f"  {r['name']:<26}{r['total_cpu_core_h']:>12.2f}{r['total_gpu_h']:>9.2f}"
              f"{r['cpu_core_s_per_pb_valid']:>12.2f}{r['gpu_s_per_pb_valid']:>12.2f}")
    print(f"\n  AutoDock CPU-core-seconds = wall × {args.autodock_cpu} threads. "
          "DiffDock GPU-bound (CPU prep not separately recorded).")
    eqfull = summ.loc[summ['method'] == 'equibind', 'total_wall_full_h'].values
    print(f"  EquiBind best-config wall scaled from per-pose compute share "
          f"(full-run wall = {eqfull} h).")
    print(f"\n  CSVs written to: {out_dir}/")

    if not args.no_plot:
        print(f"  Figure (wall-clock): {make_wall_figure(pc, summ, out_dir)}")
        print(f"  Figure (resources):  {make_resource_figure(summ, out_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
