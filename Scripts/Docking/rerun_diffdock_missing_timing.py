#!/usr/bin/env python3
"""Re-run DiffDock on the complexes whose docking timing was lost.

Some DiffDock complexes carry ``elapsed_time_s == 0`` in their ``docking_log.csv``.
That happens when a complex already had poses on disk during a later run: the
runner classifies it as *skipped* (``elapsed_time = 0.0``) and writes that zero as
the row's timing. ``docking_effort_comparison.py`` then DROPS every such complex
from the common set ("complexes whose timing was lost to a skip/empty re-run …
are dropped per method"), so those complexes silently fall out of the effort
comparison. This script re-docks exactly those complexes so a real wall-clock
time is recorded.

Two DiffDock datasets are affected (JKU has none):

  * benchmark  -> Dockings/Benchmark_DiffDock/<id>/docking_log.csv   (one call
                  per complex, protein+ligand from Data/PoseBuster Benchmark Set/;
                  reproduces notebook cell 9, optimization = smina+gnina)
  * orai       -> Dockings/Orai_Benchmark_DiffDock/docking_log.csv   (batched,
                  Orai receptors × benchmark ligands; reproduces notebook cell 68,
                  optimization = gnina)

Why pruning is required first
-----------------------------
``append_log_rows()`` skips any ``combo_name`` already present in the log
(``overwrite=False`` on every call site), so a bare re-run would re-dock but the
new wall-clock would never replace the stale 0. We therefore delete the affected
rows from ``docking_log.csv`` (backed up to ``*.pretiming.bak``) and drop those
combos from ``failed_docking.json`` before docking, then dock with
``overwrite_existing=True``.

Modes
-----
  (no mode)      DRY RUN — list the affected complexes and exit. Nothing changes.

  --timing-only  RECOMMENDED. Re-dock each target into a scratch directory with
                 optimization off, read the fresh docking wall-clock, and patch
                 ONLY ``elapsed_time_s`` back into the original log. Existing
                 poses / optimized_*/ dirs / PoseBusters caches are left intact.

  --redock       Full re-dock in place with overwrite_existing=True. This
                 REGENERATES the poses (DiffDock is stochastic) and re-runs the
                 configured optimizer, overwriting the existing SDFs for those
                 complexes. Any downstream PoseBusters / analysis for them must be
                 re-run. Use only if you want brand-new poses, not just timing.

Usage (from the repo root, same env as the notebook docking cells)::

    python Scripts/Docking/rerun_diffdock_missing_timing.py                 # preview
    python Scripts/Docking/rerun_diffdock_missing_timing.py --timing-only
    python Scripts/Docking/rerun_diffdock_missing_timing.py --redock
    python Scripts/Docking/rerun_diffdock_missing_timing.py --timing-only --datasets benchmark
    python Scripts/Docking/rerun_diffdock_missing_timing.py --timing-only --limit 3   # smoke test
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import pandas as pd
import yaml

_HERE = Path(__file__).resolve().parent
_REPO = _HERE.parent.parent
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

from Scripts.Docking.run_diffdock import (  # noqa: E402
    run_diffdock, get_gpu_model, ERROR_LOG_FILENAME,
)

# ── Fixed locations (match notebook cells 9 & 68) ────────────────────────────
CONFIG_PATH        = _REPO / "Scripts/Docking/diffdock_docking_config.yaml"

BENCH_DIFFDOCK_DIR = _REPO / "Dockings/Benchmark_DiffDock"
BENCH_INPUT_DIR    = _REPO / "Data/PoseBuster Benchmark Set"
BENCH_LOG_DIR      = _REPO / "Dockings/Logs/benchmark_diffdock_logs"

ORAI_DIFFDOCK_DIR  = _REPO / "Dockings/Orai_Benchmark_DiffDock"
ORAI_LIGAND_DIR    = _REPO / "Data/Ligands/PoseBuster_Benchmark_Set"
ORAI_RECEPTOR_DIR  = _REPO / "Data/Receptors"
ORAI_LOG_DIR       = _REPO / "Dockings/Logs/orai_benchmark_diffdock_logs"

SCRATCH_ROOT       = _REPO / "Dockings/_timing_rerun_scratch"


def _is_zero(v) -> bool:
    """True when a docking_log elapsed_time_s cell counts as 'no timing'."""
    if pd.isna(v) or v == "":
        return True
    try:
        return float(v) == 0.0
    except (TypeError, ValueError):
        return False


# ═══════════════════════════════════════════════════════════════════════════
# Discovery
# ═══════════════════════════════════════════════════════════════════════════

def discover_benchmark() -> List[dict]:
    """Complexes in Benchmark_DiffDock with a zero/blank elapsed_time_s row."""
    targets: List[dict] = []
    for log_path in sorted(BENCH_DIFFDOCK_DIR.glob("*/docking_log.csv")):
        pdb_id = log_path.parent.name
        try:
            df = pd.read_csv(log_path)
        except Exception as exc:
            print(f"  ! unreadable {log_path}: {exc}")
            continue
        zero = df[df["elapsed_time_s"].map(_is_zero)]
        if zero.empty:
            continue
        protein = BENCH_INPUT_DIR / pdb_id / f"{pdb_id}_protein.pdb"
        ligand = BENCH_INPUT_DIR / pdb_id / f"{pdb_id}_ligand_start_conf.sdf"
        if not (protein.exists() and ligand.exists()):
            print(f"  ! {pdb_id}: input files missing "
                  f"(protein={protein.exists()}, ligand={ligand.exists()}) — skipping")
            continue
        targets.append({
            "pdb_id": pdb_id,
            "combos": list(zero["combo_name"]),
            "log_path": log_path,
            "out_dir": log_path.parent,
            "protein": protein,
            "ligand": ligand,
        })
    return targets


def discover_orai() -> List[dict]:
    """(ligand, receptor) pairs in Orai_Benchmark_DiffDock with zero/blank timing."""
    log_path = ORAI_DIFFDOCK_DIR / "docking_log.csv"
    if not log_path.exists():
        print(f"  ! no Orai log at {log_path}")
        return []
    df = pd.read_csv(log_path)
    zero = df[df["elapsed_time_s"].map(_is_zero)]
    targets: List[dict] = []
    for _, r in zero.iterrows():
        receptor = str(r["protein_name"])                       # e.g. Orai1WT-START-Fr0
        lig_name = str(r["ligand_name"])                        # e.g. 7KC5_BJZ_start_conf
        pdb_id = lig_name.replace("_start_conf", "")
        protein = ORAI_RECEPTOR_DIR / f"{receptor}.pdb"
        ligand = ORAI_LIGAND_DIR / f"{pdb_id}_ligand_start_conf.sdf"
        if not (protein.exists() and ligand.exists()):
            print(f"  ! {r['combo_name']}: input files missing "
                  f"(protein={protein.exists()}, ligand={ligand.exists()}) — skipping")
            continue
        targets.append({
            "combo": str(r["combo_name"]),
            "receptor": receptor,
            "pdb_id": pdb_id,
            "protein": protein,
            "ligand": ligand,
        })
    return targets


# ═══════════════════════════════════════════════════════════════════════════
# Log / error-log surgery
# ═══════════════════════════════════════════════════════════════════════════

def prune_log_rows(log_path: Path, combos: set) -> int:
    """Remove rows whose combo_name is in ``combos``; back up once. Returns #removed."""
    if not log_path.exists():
        return 0
    df = pd.read_csv(log_path)
    mask = df["combo_name"].isin(combos)
    n = int(mask.sum())
    if n == 0:
        return 0
    backup = log_path.with_suffix(log_path.suffix + ".pretiming.bak")
    if not backup.exists():
        shutil.copy2(log_path, backup)
    df[~mask].to_csv(log_path, index=False)
    return n


def prune_error_log(out_dir: Path, combos: set) -> int:
    """Drop combos from failed_docking.json so they aren't skipped as 'previously failed'."""
    err_path = out_dir / ERROR_LOG_FILENAME
    if not err_path.exists():
        return 0
    try:
        data = json.loads(err_path.read_text())
    except Exception:
        return 0
    removed = [k for k in combos if k in data]
    for k in removed:
        del data[k]
    if removed:
        err_path.write_text(json.dumps(data, indent=2))
    return len(removed)


# ═══════════════════════════════════════════════════════════════════════════
# Config builders (mirror the notebook cells)
# ═══════════════════════════════════════════════════════════════════════════

def _base_cfg() -> dict:
    with open(CONFIG_PATH) as f:
        return yaml.safe_load(f)


def _read_new_timing(scratch_log: Path) -> Dict[str, float]:
    if not scratch_log.exists():
        return {}
    df = pd.read_csv(scratch_log)
    return {str(c): float(t) for c, t in zip(df["combo_name"], df["elapsed_time_s"])
            if not _is_zero(t)}


def patch_timing(log_path: Path, new_timing: Dict[str, float]) -> int:
    """Update elapsed_time_s in an existing log for the given combos. Returns #patched."""
    if not log_path.exists() or not new_timing:
        return 0
    df = pd.read_csv(log_path)
    backup = log_path.with_suffix(log_path.suffix + ".pretiming.bak")
    if not backup.exists():
        shutil.copy2(log_path, backup)
    n = 0
    for combo, secs in new_timing.items():
        m = df["combo_name"] == combo
        if m.any():
            df.loc[m, "elapsed_time_s"] = round(secs, 2)
            n += int(m.sum())
    df.to_csv(log_path, index=False)
    return n


# ═══════════════════════════════════════════════════════════════════════════
# Runners
# ═══════════════════════════════════════════════════════════════════════════

def run_benchmark(targets: List[dict], mode: str, gpu_model: str,
                  dock_timeout: int = None) -> None:
    base = _base_cfg()
    base["optimization"] = "all"          # cell 9: smina + gnina
    if dock_timeout is not None:
        base["timeout_per_complex"] = dock_timeout
    for i, t in enumerate(targets, 1):
        pid = t["pdb_id"]
        print(f"\n[bench {i}/{len(targets)}] {pid}  ({len(t['combos'])} row(s))")
        if mode == "timing-only":
            scratch = SCRATCH_ROOT / "benchmark" / pid
            if scratch.exists():
                shutil.rmtree(scratch)
            scratch.mkdir(parents=True, exist_ok=True)
            cfg = dict(base)
            cfg["optimization"] = "none"          # only need the docking wall-clock
            cfg["overwrite_existing"] = True
            cfg["overwrite_error_log"] = True
            run_diffdock([t["protein"]], [t["ligand"]], scratch, cfg, gpu_model)
            new_timing = _read_new_timing(scratch / "docking_log.csv")
            patched = patch_timing(t["log_path"], new_timing)
            print(f"    patched elapsed_time_s for {patched} row(s): {new_timing}")
            shutil.rmtree(scratch, ignore_errors=True)
        else:  # redock
            prune_log_rows(t["log_path"], set(t["combos"]))
            prune_error_log(t["out_dir"], set(t["combos"]))
            cfg = dict(base)
            cfg["output_dir"] = str(t["out_dir"])
            cfg["log_dir"] = str(BENCH_LOG_DIR)
            cfg["overwrite_existing"] = True
            cfg["overwrite_error_log"] = False
            run_diffdock([t["protein"]], [t["ligand"]], t["out_dir"], cfg, gpu_model)


def run_orai(targets: List[dict], mode: str, gpu_model: str,
             dock_timeout: int = None) -> None:
    log_path = ORAI_DIFFDOCK_DIR / "docking_log.csv"
    base = _base_cfg()
    base["optimization"] = "gnina"        # cell 68: gnina only
    base["batch_mode"] = True
    if dock_timeout is not None:
        base["timeout_per_complex"] = dock_timeout

    by_rec: Dict[str, List[dict]] = defaultdict(list)
    for t in targets:
        by_rec[t["receptor"]].append(t)

    if mode == "timing-only":
        all_new: Dict[str, float] = {}
        for rec, ts in by_rec.items():
            scratch = SCRATCH_ROOT / "orai" / rec
            if scratch.exists():
                shutil.rmtree(scratch)
            scratch.mkdir(parents=True, exist_ok=True)
            cfg = dict(base)
            cfg["optimization"] = "none"
            cfg["overwrite_existing"] = True
            cfg["overwrite_error_log"] = True
            proteins = [ts[0]["protein"]]
            ligands = [t["ligand"] for t in ts]
            print(f"\n[orai timing] {rec}: {len(ligands)} ligand(s)")
            run_diffdock(proteins, ligands, scratch, cfg, gpu_model)
            all_new.update(_read_new_timing(scratch / "docking_log.csv"))
            shutil.rmtree(scratch, ignore_errors=True)
        patched = patch_timing(log_path, all_new)
        print(f"\n    patched elapsed_time_s for {patched} Orai row(s)")
    else:  # redock
        combos = {t["combo"] for t in targets}
        prune_log_rows(log_path, combos)
        prune_error_log(ORAI_DIFFDOCK_DIR, combos)
        for rec, ts in by_rec.items():
            cfg = dict(base)
            cfg["output_dir"] = str(ORAI_DIFFDOCK_DIR)
            cfg["log_dir"] = str(ORAI_LOG_DIR)
            cfg["overwrite_existing"] = True
            cfg["overwrite_error_log"] = False
            proteins = [ts[0]["protein"]]
            ligands = [t["ligand"] for t in ts]
            print(f"\n[orai redock] {rec}: {len(ligands)} ligand(s)")
            run_diffdock(proteins, ligands, ORAI_DIFFDOCK_DIR, cfg, gpu_model)


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--timing-only", action="store_true",
                      help="Re-measure timing in a scratch dir and patch it back; "
                           "keeps existing poses (recommended).")
    mode.add_argument("--redock", action="store_true",
                      help="Full re-dock in place; OVERWRITES the existing poses.")
    ap.add_argument("--datasets", default="benchmark,orai",
                    help="Comma list: benchmark, orai (default both).")
    ap.add_argument("--limit", type=int, default=None,
                    help="Only process the first N targets per dataset (smoke test).")
    ap.add_argument("--dock-timeout", type=int, default=None,
                    help="Override timeout_per_complex (s) for the docking; use a "
                         "larger value (e.g. 1800) to re-measure complexes that "
                         "exceeded the config's 600s cap.")
    args = ap.parse_args()

    wanted = {d.strip() for d in args.datasets.split(",") if d.strip()}
    run_mode = "timing-only" if args.timing_only else "redock" if args.redock else "dry-run"

    print("Discovering DiffDock complexes with no timing recorded (elapsed_time_s == 0)\n")
    bench = discover_benchmark() if "benchmark" in wanted else []
    orai = discover_orai() if "orai" in wanted else []

    print(f"\nBenchmark_DiffDock : {len(bench)} complexes")
    for t in bench:
        print(f"    {t['pdb_id']}")
    n_orai_rec = len({t['receptor'] for t in orai})
    print(f"\nOrai_Benchmark_DiffDock : {len(orai)} (ligand,receptor) pairs "
          f"across {n_orai_rec} receptor(s)")
    for t in orai:
        print(f"    {t['combo']}")

    if run_mode == "dry-run":
        print("\n[DRY RUN] No changes made. Re-run with --timing-only (keeps poses) "
              "or --redock (regenerates poses).")
        return

    if args.limit is not None:
        bench = bench[:args.limit]
        orai = orai[:args.limit]

    if run_mode == "redock":
        print("\n⚠ --redock will OVERWRITE existing poses for the complexes above "
              "and re-run their optimizer.\n  Downstream PoseBusters / analysis for "
              "them will need to be re-run.")

    gpu_model = get_gpu_model()
    print(f"\nMode: {run_mode}   GPU: {gpu_model}")

    if args.dock_timeout is not None:
        print(f"Dock timeout override: {args.dock_timeout}s per complex")
    if bench:
        run_benchmark(bench, run_mode, gpu_model, dock_timeout=args.dock_timeout)
    if orai:
        run_orai(orai, run_mode, gpu_model, dock_timeout=args.dock_timeout)

    print("\nDone. Re-scan with a no-arg run to confirm no zero-timing rows remain.")


if __name__ == "__main__":
    main()
