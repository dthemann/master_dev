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
  * DiffDock : Dockings/Benchmark_DiffDock/<id>/docking_log.csv (docking wall, GPU)
               + optimization_log.csv (post-hoc rescoring, --diffdock-refine)
               GPU-seconds = docking wall + gnina rescoring (serial, CUDA CNN, GPU
               because the benchmark ran gnina_use_gpu true); with --diffdock-refine
               smina — or --no-diffdock-gnina-gpu — the rescoring is CPU instead.
  * EquiBind : Dockings/Benchmark_Equibind/<id>/pipeline_summary.json
               -> global_timing.pipeline_wall_time_s       (wall)
               -> all_results[].{prep,dock,post,refine}_time_s  (per-pose work)
               GPU-seconds = Σ dock_time_s (inference) + Σ refine_time_s only if gnina
               ran on the GPU (--eq-gnina-gpu). The benchmark ran EquiBind's gnina with
               --no_gpu (gnina_use_gpu false), so gnina refine is CPU-core-seconds;
               CPU-core-seconds = Σ (prep+post)_time_s (+ refine unless it was GPU gnina).

EquiBind "best configuration" (default **unguided + gnina**): the benchmark run
produced the full clamp×refine×pocket fan-out; of those, unguided + gnina gives the
best PoseBusters validity (~52% vs ~40% for smina). We therefore (a) count only that
variant's poses, (b) sum *only those poses'* per-pose resource-seconds (GPU/CPU), and
(c) estimate its wall-clock as those poses' OWN compute spread across the pipeline's
worker pool (config.n_parallel_workers) — a pure best-config-only estimate that does
NOT depend on the rest of the fan-out. The real full-run (all-variant) wall is kept in
the CSV (total_wall_full_h) for transparency. The refiner is selectable with
--eq-refine (smina is CPU; gnina is CPU too unless --eq-gnina-gpu, matching the
benchmark's gnina_use_gpu false); DiffDock's rescoring is selected in parallel with
--diffdock-refine / --diffdock-gnina-gpu (its gnina ran on the GPU).

All methods are compared on the COMMON set of complexes that have usable per-complex
timing for *every* method: complexes whose timing was lost to a skip/empty re-run
(elapsed_time_s / total_time_seconds == 0) are dropped per method, then intersected,
so total effort, pose counts and validity are like-for-like rather than spanning
different subsets.

Run in the analysis env (conda env ``vina`` — pandas + matplotlib):

    python Scripts/Analysis/docking_effort_comparison.py \
        --ids-file "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
"""
from __future__ import annotations

import argparse
import json
import os
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

# Method label -> (pretty name, bar colour, device used). The label + device are
# rewritten per run in main() to reflect the selected refine variant.
METHODS = {
    "autodock":  ("AutoDock Vina", "#4C72B0", "CPU"),
    "diffdock":  ("DiffDock (gnina)", "#55A868", "GPU"),
    "equibind":  ("EquiBind (unguided+gnina)", "#C44E52", "GPU+CPU"),
}
# EquiBind per-pose timing fields, classified by the hardware that does the work.
# Inference is GPU; conformer prep / UFF / IO are CPU. The REFINE step's hardware
# depends on the tool AND the pipeline's gnina_use_gpu flag: smina is always CPU;
# gnina is a CUDA-capable CNN but is run with --no_gpu (CPU) unless gnina_use_gpu is
# set. The benchmark ran EquiBind's gnina on CPU (gnina_use_gpu: false) and DiffDock's
# gnina on GPU (gnina_use_gpu: true) — hence the per-pipeline --*-gnina-gpu flags.
EQ_INFER_FIELD = "dock_time_s"                     # EquiBind model inference (GPU)
EQ_HOST_FIELDS = ("prep_time_s", "post_time_s")    # RDKit conformers / UFF / IO (CPU)

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
        t = float(json.loads(p.read_text())["overall"]["total_time_seconds"])
    except Exception:
        return None
    # A skip/empty re-run overwrites the original timed summary with
    # total_combinations=0 / total_time_seconds=0; treat that as "no timing"
    # (mirrors the DiffDock reader) rather than a real 0-second dock.
    return t if t > 0 else None


def _diffdock_time(cid: str, root: Path) -> Optional[float]:
    """Real per-complex DiffDock docking wall-clock from docking_log.csv.

    Use docking_log.csv's elapsed_time_s (recorded on the FIRST real dock of each
    complex, e.g. the 2026-05-30 production run). Multiple rows (appended skip-runs
    write 0.0) -> take the max, i.e. the one real dock. Returns None when no positive
    time is recorded (a skip-run whose original production timing was overwritten).

    docking_summary.json is NOT consulted: its overall.total_time_seconds is
    co-zeroed on every skip-run (skipped:1 -> 0.0), so it can never recover a
    complex the log has already lost.
    """
    log = root / cid / "docking_log.csv"
    if log.exists():
        try:
            import csv
            with open(log) as f:
                ts = [float(r.get("elapsed_time_s") or 0) for r in csv.DictReader(f)]
            mx = max(ts) if ts else 0.0
            if mx > 0:
                return mx
        except Exception:
            pass
    return None


def _diffdock_opt_time(cid: str, root: Path, tool: str) -> Optional[float]:
    """Total post-hoc {tool} rescoring/minimization time from optimization_log.csv.

    DiffDock's ``optimized_<tool>`` variants are produced by minimizing each DiffDock
    pose with smina (CPU) or gnina (CUDA CNN, GPU or CPU depending on gnina_use_gpu).
    The rescoring is serial per pose (verified: timestamp span ≈ Σ elapsed), so the
    summed per-pose elapsed_time_s is the added wall-clock. Returns None if the log
    has no rows for ``tool``.
    """
    p = root / cid / "optimization_log.csv"
    if not p.exists():
        return None
    try:
        import csv
        with open(p) as f:
            ts = [float(r.get("elapsed_time_s") or 0.0)
                  for r in csv.DictReader(f) if r.get("tool") == tool]
    except Exception:
        return None
    return sum(ts) if ts else None


def _diffdock_effort(cid: str, root: Path, refine: str,
                     gnina_gpu: bool = True) -> Optional[tuple]:
    """Return (wall_s, gpu_s, cpu_core_s) for the requested DiffDock variant.

    ``raw``   = DiffDock docking only (GPU).
    ``smina`` = docking (GPU) + smina minimization (CPU).
    ``gnina`` = docking (GPU) + gnina CNN minimization (GPU if gnina_gpu else CPU).
    None if the docking wall — or, for a rescored variant, its optimization time —
    is missing.
    """
    dock = _diffdock_time(cid, root)          # docking wall (GPU)
    if dock is None:
        return None
    if refine in ("smina", "gnina"):
        opt = _diffdock_opt_time(cid, root, refine)
        if opt is None:
            return None
        if refine == "gnina" and gnina_gpu:
            return dock + opt, dock + opt, 0.0        # docking + gnina both GPU
        return dock + opt, dock, opt                  # docking GPU; smina / --no_gpu gnina on CPU
    return dock, dock, 0.0                            # raw docking only


def _equibind_time(cid: str, root: Path, eq_mode: str, eq_refine: str,
                   eq_clamp: str, gnina_gpu: bool = False) -> Optional[tuple]:
    """Return (wall_full_s, wall_best_s, best_gpu_s, best_cpu_core_s) for one complex.

    best_gpu_s / best_cpu_core_s are the summed per-pose resource-seconds of the
    best-config (unguided + ``eq_refine`` + clampOFF) poses ONLY, split by hardware:
    inference is GPU, prep/post are CPU, refine follows the tool AND its device
    (smina -> CPU; gnina -> GPU only if gnina_gpu, else CPU via --no_gpu). wall_best
    estimates the wall to run *only* that best config — those poses' compute spread
    across config.n_parallel_workers — so it depends solely on the best-config poses,
    not the rest of the fan-out. wall_full is the real full-run (all-variant) wall,
    kept only for transparency.
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
    n_workers = float(d.get("config", {}).get("n_parallel_workers") or 32) or 32.0

    refine_on_gpu = (eq_refine == "gnina") and gnina_gpu       # gnina only, and only if run on GPU

    def gpu_s(r):
        g = float(r.get(EQ_INFER_FIELD) or 0.0)                 # inference (GPU)
        if refine_on_gpu:
            g += float(r.get("refine_time_s") or 0.0)          # gnina CNN (GPU)
        return g

    def cpu_s(r):
        c = sum(float(r.get(k) or 0.0) for k in EQ_HOST_FIELDS)  # prep/post (CPU)
        if not refine_on_gpu:
            c += float(r.get("refine_time_s") or 0.0)           # smina (CPU)
        return c

    best = [r for r in rows if r.get("mode") == eq_mode
            and r.get("refine_variant") == eq_refine
            and r.get("clamp_variant") == eq_clamp]
    best_gpu = sum(gpu_s(r) for r in best)
    best_cpu = sum(cpu_s(r) for r in best)
    best_c = best_gpu + best_cpu
    # Wall to run ONLY the best config: its OWN compute-seconds spread across the
    # pipeline's worker pool (config.n_parallel_workers). Depends solely on the
    # unguided+<refine> poses — the rest of the fan-out (other modes / raw / smina)
    # no longer enters. Idealised (assumes perfect n-worker overlap; the full run's
    # measured overlap was a touch lower, so this is a slight lower bound).
    wall_best = best_c / n_workers if n_workers > 0 else np.nan
    return wall_full, wall_best, best_gpu, best_cpu


# ════════════════════════════════════════════════════════════════════════
# Aggregation
# ════════════════════════════════════════════════════════════════════════

def _per_pose_methods(args) -> Dict[str, str]:
    """Map each method key to its ``method`` value in per_pose_metrics.csv, per the
    selected variant (EquiBind mode/refine, DiffDock refine)."""
    dd = "diffdock" if args.diffdock_refine == "raw" else f"diffdock_{args.diffdock_refine}"
    return {
        "autodock": "autodock",
        "diffdock": dd,
        "equibind": f"equibind_{args.eq_mode}_{args.eq_refine}",   # clamp-off drops the clamp token
    }


def load_pose_counts(csv: Path, ids: Optional[set],
                     per_pose_method: Dict[str, str]) -> Dict[str, pd.DataFrame]:
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
    for key, per_pose_name in per_pose_method.items():
        sub = df[df["method"] == per_pose_name]
        if per_pose_name == "diffdock" and "pose_file" in sub.columns:
            # Base DiffDock writes its top pose under TWO filenames: rankN.sdf and
            # rankN_confidence-X.sdf (byte-identical). Drop the no-confidence copies
            # so each physical pose is counted once. (The smina/gnina rescored variants
            # carry a tool suffix, so they have no such twin.)
            base = sub["pose_file"].astype(str).map(os.path.basename)
            sub = sub[~base.str.match(r"rank\d+\.sdf$")]
        g = sub.groupby("cid").agg(generated=("pb_valid", "size"),
                                   pb_valid=("pb_valid", "sum"))
        out[key] = g
    return out


def build_table(args) -> tuple:
    ids = None
    if args.ids_file and Path(args.ids_file).exists():
        ids = {ln.strip() for ln in Path(args.ids_file).read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    counts = load_pose_counts(Path(args.per_pose_csv), ids, _per_pose_methods(args))
    ad_root = Path(args.autodock_dir)
    dd_root = Path(args.diffdock_dir)
    eb_root = Path(args.equibind_dir)

    rows: List[dict] = []
    per_method_timed = {k: 0 for k in METHODS}     # complexes with usable timing
    eq_processed = eq_best = 0                      # EquiBind best-config filter sanity
    for key in METHODS:
        for cid, crow in counts[key].iterrows():
            wall = wall_full = cpu_core_s = gpu_s = None
            if key == "autodock":
                wall = _autodock_time(cid, ad_root, args.autodock_prep)
                if wall is not None:
                    wall_full, cpu_core_s, gpu_s = wall, wall * args.autodock_cpu, 0.0
            elif key == "diffdock":
                dd = _diffdock_effort(cid, dd_root, args.diffdock_refine, args.diffdock_gnina_gpu)
                if dd is not None:
                    wall, gpu_s, cpu_core_s = dd
                    wall_full = wall
            else:
                eq_processed += 1
                t = _equibind_time(cid, eb_root, args.eq_mode, args.eq_refine, args.eq_clamp,
                                   args.eq_gnina_gpu)
                if t:
                    wall_full, wall, gpu_s, cpu_core_s = t[0], t[1], t[2], t[3]
                    if (gpu_s + cpu_core_s) > 0:    # best-config filter matched real poses
                        eq_best += 1
            if wall is None or not np.isfinite(wall) or wall <= 0:
                continue
            per_method_timed[key] += 1
            rows.append({
                "method": key, "cid": cid,
                "wall_s": float(wall), "wall_full_s": float(wall_full),
                "cpu_core_s": float(cpu_core_s), "gpu_s": float(gpu_s),
                "generated": int(crow["generated"]), "pb_valid": int(crow["pb_valid"]),
            })
    if eq_processed and not eq_best:
        raise SystemExit(
            "ERROR: the EquiBind best-config filter matched no poses for any complex "
            f"(--eq-mode={args.eq_mode} --eq-refine={args.eq_refine} --eq-clamp={args.eq_clamp}).\n"
            "       Check these against the mode / refine_variant / clamp_variant values in "
            "pipeline_summary.json — clamp_variant is None for clamp-off runs, so leave "
            "--eq-clamp unset (do NOT pass clampOFF).")
    pc = pd.DataFrame(rows)

    # Compare methods on the COMMON set of complexes that have usable timing for
    # ALL methods, so total effort / pose counts / validity are like-for-like
    # (methods otherwise cover different subsets: skip-runs drop out per method).
    common = set()
    dropped_to_common = {}
    if not pc.empty:
        present = [k for k in METHODS if (pc["method"] == k).any()]
        common = set.intersection(*(set(pc.loc[pc["method"] == k, "cid"]) for k in present))
        dropped_to_common = {k: per_method_timed[k] - len(common) for k in present}
        pc = pc[pc["cid"].isin(common)].copy()
    meta = {"common_n": len(common), "per_method_timed": per_method_timed,
            "dropped_to_common": dropped_to_common}

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
    return pc, summ, meta


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


def _save_panels(panels, panel_dir, plt, figsize=(6.2, 4.8)):
    """Render each (name, draw_fn) as its own single-panel PNG (no (a)/(b) label)."""
    panel_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for nm, draw in panels:
        f, a = plt.subplots(figsize=figsize)
        draw(a)
        f.tight_layout()
        pp = panel_dir / f"{nm}.png"
        f.savefig(pp, dpi=150, bbox_inches="tight"); plt.close(f)
        written.append(pp)
    return written


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

    def _style(ax):
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25)

    def p_overall(ax):
        _bars(ax, x, s["total_wall_h"].to_numpy(float), colors, "{:.2f} h", devices)
        ax.set_title("Overall computational effort (wall-clock to dock the benchmark)\n"
                     "EquiBind = best-config-only estimate (full fan-out wall in CSV)")
        ax.set_ylabel("Total wall-clock time (hours)"); _style(ax)

    def p_distribution(ax):
        data = [pc[pc["method"] == k]["wall_s"].to_numpy(float) for k in order]
        bp = ax.boxplot(data, positions=x, widths=0.6, patch_artist=True,
                        showfliers=False, medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.75)
        ax.set_yscale("log")
        ax.set_title("Per-complex docking wall-clock time (distribution)\n"
                     "EquiBind = best-config-only estimate (compute ÷ workers)")
        ax.set_ylabel("Wall-clock time per complex (seconds, log scale)"); _style(ax)

    def p_counts(ax):
        w = 0.38
        gen = s["poses_generated"].to_numpy(float)
        val = s["poses_pb_valid"].to_numpy(float)
        ax.bar(x - w / 2, gen, w, color=colors, alpha=0.55, label="generated")
        ax.bar(x + w / 2, val, w, color=colors, label="survive PoseBusters")
        for i in range(len(order)):
            rate = s["pb_valid_rate_pct"].to_numpy(float)[i]
            if rate == rate:
                ax.text(x[i] + w / 2, val[i], f"{rate:.1f}%", ha="center",
                        va="bottom", fontsize=8)
        ax.set_title("Poses generated vs. poses that survive\nthe PoseBusters tests")
        ax.set_ylabel("Number of poses"); ax.legend(fontsize=8); _style(ax)

    def p_per_generated(ax):
        _bars(ax, x, s["wall_s_per_generated"].to_numpy(float), colors, "{:.2f} s")
        ax.set_title("Timing per generated pose")
        ax.set_ylabel("Wall-clock time per generated pose (seconds)"); _style(ax)

    def p_per_valid(ax):
        _bars(ax, x, s["wall_s_per_pb_valid"].to_numpy(float), colors, "{:.2f} s")
        ax.set_title("Timing per pose that survives the PoseBusters tests")
        ax.set_ylabel("Wall-clock time per valid pose (seconds)"); _style(ax)

    def p_rate(ax):
        _bars(ax, x, s["pb_valid_rate_pct"].to_numpy(float), colors, "{:.1f}%", ymax_mult=1.3)
        ax.set_title("PoseBusters validity rate\n(fraction of generated poses that survive)")
        ax.set_ylabel("Valid poses (percent of generated)")
        ax.set_ylim(0, min(100, ax.get_ylim()[1])); _style(ax)

    panels = [
        ("wall_overall_total", p_overall),
        ("wall_per_complex_distribution", p_distribution),
        ("wall_poses_generated_vs_valid", p_counts),
        ("wall_per_generated_pose", p_per_generated),
        ("wall_per_valid_pose", p_per_valid),
        ("wall_validity_rate", p_rate),
    ]

    # combined 2×3 grid (kept — the notebook displays this file)
    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    for a, (_n, draw) in zip(ax.ravel(), panels):
        draw(a)
    _label_panels(ax)
    fig.suptitle("Docking WALL-CLOCK (time-to-result) — AutoDock Vina vs DiffDock vs "
                 "EquiBind (best config), PoseBuster benchmark (common timed set).\n"
                 "AutoDock = CPU (multi-thread), DiffDock/EquiBind = single GPU — comparable "
                 "as elapsed time, not as identical-hardware work.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = out_dir / "docking_effort_comparison.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    return [p] + _save_panels(panels, out_dir / "panels", plt)


def make_resource_figure(summ, out_dir):
    """CPU-core-seconds vs GPU-seconds — the true hardware work (different currencies)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order)); w = 0.38

    def grouped(ax, cpu, gpu, fmt, labels=("CPU-core-seconds", "GPU-seconds")):
        b1 = ax.bar(x - w / 2, cpu, w, color=CPU_COLOR, label=labels[0])
        b2 = ax.bar(x + w / 2, gpu, w, color=GPU_COLOR, label=labels[1])
        ymax = max(list(cpu) + list(gpu) + [0]) or 1
        for bars, vals in ((b1, cpu), (b2, gpu)):
            for rect, v in zip(bars, vals):
                if v == v:
                    ax.text(rect.get_x() + rect.get_width() / 2, v + 0.01 * ymax,
                            fmt.format(v), ha="center", va="bottom", fontsize=8, rotation=0)
        ax.set_ylim(0, ymax * 1.25)
        ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
        ax.grid(axis="y", alpha=0.25); ax.legend(fontsize=8)

    def p_overall(ax):
        grouped(ax, s["total_cpu_core_h"].to_numpy(float), s["total_gpu_h"].to_numpy(float),
                "{:.1f}", labels=("CPU-core-hours", "GPU-hours"))
        ax.set_title("Overall hardware effort\n(total resource consumed)")
        ax.set_ylabel("Resource used (hours): CPU-core-hours | GPU-hours")

    def p_per_generated(ax):
        grouped(ax, s["cpu_core_s_per_generated"].to_numpy(float),
                s["gpu_s_per_generated"].to_numpy(float), "{:.2f}")
        ax.set_title("Hardware effort per generated pose")
        ax.set_ylabel("Resource per generated pose (seconds)")

    def p_per_valid(ax):
        grouped(ax, s["cpu_core_s_per_pb_valid"].to_numpy(float),
                s["gpu_s_per_pb_valid"].to_numpy(float), "{:.2f}")
        ax.set_title("Hardware effort per pose that survives\nthe PoseBusters tests")
        ax.set_ylabel("Resource per valid pose (seconds)")

    panels = [
        ("resource_overall_total", p_overall),
        ("resource_per_generated_pose", p_per_generated),
        ("resource_per_valid_pose", p_per_valid),
    ]

    # combined 1×3 grid (kept — the notebook displays this file)
    fig, ax = plt.subplots(1, 3, figsize=(18, 5.6))
    for a, (_n, draw) in zip(ax, panels):
        draw(a)
    _label_panels(ax)
    fig.suptitle("Docking HARDWARE EFFORT — CPU-core-seconds vs GPU-seconds (the actual work, "
                 "parallelism removed).\nA GPU-second is NOT a CPU-core-second — the two bars are "
                 "different currencies and are never summed. AutoDock is CPU-only; DiffDock docking + "
                 "GPU gnina rescoring are GPU; EquiBind is GPU only for inference — its gnina refine "
                 "ran on CPU (--no_gpu), so it is CPU-dominated.", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    p = out_dir / "docking_resource_effort.png"
    fig.savefig(p, dpi=130, bbox_inches="tight"); plt.close(fig)

    return [p] + _save_panels(panels, out_dir / "panels", plt)


def make_time_per_pose_figure(summ, out_dir):
    """Wall-clock time to produce ONE pose vs ONE PoseBusters-valid pose, per toolchain,
    with each bar split into the share of that time delivered by CPU vs GPU work.

    Two grouped bars per method (per generated pose, per surviving pose); each bar is
    stacked into CPU-delivered and GPU-delivered wall-seconds, apportioned by the
    method's CPU-core-second : GPU-second ratio (so the stack total is the per-pose time).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    order = [k for k in METHODS if k in set(summ["method"])]
    names = [METHODS[k][0] for k in order]
    s = summ.set_index("method").reindex(order)
    x = np.arange(len(order)); w = 0.38

    cpu_h = s["total_cpu_core_h"].to_numpy(float)
    gpu_h = s["total_gpu_h"].to_numpy(float)
    tot_h = cpu_h + gpu_h
    gpu_frac = np.divide(gpu_h, tot_h, out=np.zeros_like(gpu_h), where=tot_h > 0)
    cpu_frac = 1.0 - gpu_frac                              # share of the wall each hardware delivered

    metrics = [(s["wall_s_per_generated"].to_numpy(float), -w / 2),
               (s["wall_s_per_pb_valid"].to_numpy(float), +w / 2)]

    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    allv = [v for wall, _ in metrics for v in wall if v == v]
    ymax = (max(allv) if allv else 1) or 1
    for wall, off in metrics:
        cpu_part = cpu_frac * wall
        gpu_part = gpu_frac * wall
        ax.bar(x + off, cpu_part, w, color=CPU_COLOR)
        ax.bar(x + off, gpu_part, w, bottom=cpu_part, color=GPU_COLOR)
        for xi, wv in zip(x + off, wall):
            if wv == wv:
                ax.text(xi, wv + 0.012 * ymax, f"{wv:.2f} s",
                        ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, ymax * 1.24)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation=15, ha="right", fontsize=9)
    ax.set_ylabel("Wall-clock time per pose (seconds)")
    ax.set_title("Wall-clock time per pose, split into CPU- vs GPU-delivered work\n"
                 "left bar = per generated pose; right bar = per pose surviving PoseBusters")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(handles=[Patch(facecolor=CPU_COLOR, label="CPU-delivered"),
                       Patch(facecolor=GPU_COLOR, label="GPU-delivered")], fontsize=9)
    fig.tight_layout()
    p = out_dir / "docking_time_per_pose.png"
    fig.savefig(p, dpi=150, bbox_inches="tight"); plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv",
                    default="posebusters_results/benchmark/dock/"
                            "pose_comparison_report/per_pose_metrics.csv")
    ap.add_argument("--autodock-dir", default="Dockings/Benchmark")
    ap.add_argument("--autodock-prep", default="meeko", choices=("meeko", "mgl_tools"))
    ap.add_argument("--autodock-cpu", type=int, default=32,
                    help="CPU threads Vina used per complex (for CPU-core-seconds). "
                         "Benchmark runs used 32.")
    ap.add_argument("--diffdock-dir", default="Dockings/Benchmark_DiffDock")
    ap.add_argument("--diffdock-refine", default="gnina", choices=("raw", "smina", "gnina"),
                    help="DiffDock variant to time: raw docking, or docking + smina/gnina "
                         "post-hoc rescoring (adds optimization_log.csv time).")
    ap.add_argument("--diffdock-gnina-gpu", action=argparse.BooleanOptionalAction, default=True,
                    help="DiffDock ran gnina's CNN on the GPU (benchmark: gnina_use_gpu true). "
                         "Use --no-diffdock-gnina-gpu if it was run with --no_gpu (CPU).")
    ap.add_argument("--equibind-dir", default="Dockings/Benchmark_Equibind")
    ap.add_argument("--eq-mode", default="unguided")
    ap.add_argument("--eq-refine", default="gnina", choices=("raw", "smina", "gnina"),
                    help="EquiBind refine variant. gnina gives the best PoseBusters validity; "
                         "smina is CPU.")
    ap.add_argument("--eq-gnina-gpu", action=argparse.BooleanOptionalAction, default=False,
                    help="EquiBind ran gnina's CNN on the GPU. Benchmark used --no_gpu "
                         "(gnina_use_gpu false, GPU busy with inference), so gnina refine is "
                         "CPU work by default; pass --eq-gnina-gpu if it was run on the GPU.")
    ap.add_argument("--eq-clamp", default=None,
                    help="EquiBind clamp variant (clampON/clampOFF). Leave unset for "
                         "clamp-off runs, whose poses carry no clamp token "
                         "(clamp_variant=None in pipeline_summary.json).")
    ap.add_argument("--ids-file", default="Data/PoseBuster Benchmark Set/"
                                          "posebusters_pdb_ccd_ids.txt")
    ap.add_argument("--out-dir", default="posebusters_results/benchmark/docking_effort")
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    # Rewrite method labels + device tags to reflect the selected variants and gnina device.
    if args.diffdock_refine == "raw" or (args.diffdock_refine == "gnina" and args.diffdock_gnina_gpu):
        _dd_dev = "GPU"                                   # docking (+ GPU gnina) only
    else:
        _dd_dev = "GPU+CPU"                               # docking GPU + CPU rescoring (smina / --no_gpu gnina)
    _dd_name = "DiffDock" if args.diffdock_refine == "raw" else f"DiffDock ({args.diffdock_refine})"
    METHODS["diffdock"] = (_dd_name, METHODS["diffdock"][1], _dd_dev)
    METHODS["equibind"] = (f"EquiBind ({args.eq_mode}+{args.eq_refine})",
                           METHODS["equibind"][1], "GPU+CPU")  # inference GPU + prep/post (+CPU gnina) CPU

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pc, summ, meta = build_table(args)
    if summ.empty:
        print("No data found — check the docking dirs and per-pose CSV path.")
        return 1

    pc.to_csv(out_dir / "per_complex_effort.csv", index=False)
    summ.to_csv(out_dir / "effort_summary.csv", index=False)

    # ── console report ───────────────────────────────────────────────────
    print("\n" + "=" * 100)
    _eq_cfg = f"{args.eq_mode}+{args.eq_refine}" + (f"+{args.eq_clamp}" if args.eq_clamp else "")
    print("DOCKING COMPUTATIONAL EFFORT — benchmark set "
          f"(EquiBind best config: {_eq_cfg})")
    print("=" * 100)
    _drops = ", ".join(f"{k} {v}" for k, v in meta.get("dropped_to_common", {}).items() if v)
    print(f"Compared on the COMMON set of {meta['common_n']} complexes timed for all methods"
          + (f" (dropped to reach it — {_drops})." if _drops else "."))
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
          "DiffDock docking is GPU (CPU prep not separately recorded).")
    if args.diffdock_refine != "raw":
        _dd_dev_note = ("GPU (CUDA CNN)" if (args.diffdock_refine == "gnina" and args.diffdock_gnina_gpu)
                        else "CPU")
        print(f"  DiffDock timing = docking wall + {args.diffdock_refine} rescoring "
              f"(optimization_log.csv, serial; {_dd_dev_note}).")
    if args.eq_refine == "gnina":
        print(f"  EquiBind gnina refine ran on {'GPU' if args.eq_gnina_gpu else 'CPU (--no_gpu)'} "
              "(gnina_use_gpu " + ("true" if args.eq_gnina_gpu else "false") + ").")
    eqfull = summ.loc[summ['method'] == 'equibind', 'total_wall_full_h'].values
    eqfull_str = f"{eqfull[0]:.2f}" if len(eqfull) else "n/a"
    print(f"  EquiBind wall = best-config (unguided+{args.eq_refine}) compute ÷ worker pool "
          f"(pure best-config-only estimate; full fan-out run wall = {eqfull_str} h).")
    print(f"\n  CSVs written to: {out_dir}/")

    if not args.no_plot:
        wall_files = make_wall_figure(pc, summ, out_dir)
        res_files = make_resource_figure(summ, out_dir)
        time_fig = make_time_per_pose_figure(summ, out_dir)
        print(f"  Figure (wall-clock): {wall_files[0]}")
        print(f"  Figure (resources):  {res_files[0]}")
        print(f"  Figure (time/pose):  {time_fig}")
        panels = wall_files[1:] + res_files[1:]
        print(f"  Individual panels ({len(panels)}) in {out_dir / 'panels'}/:")
        for pp in panels:
            print(f"      - {pp.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
