#!/usr/bin/env python3
"""Refine a named list of poses with an arbitrary gnina/smina flag set.

Built for the refiner-parameter study. ``rerun_equibind_refine.py`` deliberately
exposes only the two published search modes so it cannot drift from the
pipeline; this script is its experimental sibling and takes whatever extra flags
a sweep needs (``--minimize_iters``, ``--accurate_line``, ``--seed`` ...).

It refines a POSE LIST rather than a whole tree, because the interesting
experiments run on the ~1,851 poses where ``--minimize`` and ``--local_only``
disagree, not on all 9,088.

The command it builds mirrors equibind_pipeline.refine.refine_pose exactly --
same receptor resolution, same autobox, same --num_modes 1 -- and then appends
the sweep flags, so a run with no extra flags reproduces the published arm.

USAGE
    python3 Scripts/Docking/gnina_refine_param_sweep.py \
        --pose-list poses.txt \
        --dest-root Dockings/_sweep/minimize_iters_100 \
        --search minimize --extra "--minimize_iters 100" --workers 24

``--pose-list`` is one path per line, each pointing at a committed
``*__refRAW.sdf`` (the exact input the published refinement received).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional

os.environ.setdefault("EQ_GNINA", "/home/manndo/docking_tools/gnina")
os.environ.setdefault("EQ_SMINA", "/home/manndo/anaconda3/envs/equibind/bin/smina")
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

RAW_SUFFIX = "__refRAW.sdf"
_TAG_KEYS = ("minimizedAffinity", "CNNscore", "CNNaffinity", "minimizedRMSD")


def _prepared_receptor(raw: Path) -> Optional[Path]:
    """The prepared receptor for a pose, from its tree layout.

    nested:  <root>/<complex>/<lig>__<prot>/pose.sdf  with _prep_<prot> beside the combo
    flat:    <root>/<lig>__<prot>/pose.sdf            with _prep_<prot> at the root
    """
    combo = raw.parent
    protein = combo.name.split("__", 1)[1] if "__" in combo.name else None
    if not protein:
        return None
    for parent in (combo.parent, combo.parent.parent):
        prep = parent / f"_prep_{protein}"
        if prep.is_dir():
            pdbs = [p for p in prep.glob("*.pdb") if "_crop_" not in p.name]
            if len(pdbs) == 1:
                return pdbs[0]
    return None


def _parse_tags(sdf: Path, stdout: str) -> Dict[str, float]:
    out: Dict[str, float] = {}
    try:
        lines = sdf.read_text(errors="ignore").splitlines()
    except Exception:
        return out
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("> ") and "<" in s:
            key = s[s.index("<") + 1:s.rindex(">")]
            if key in _TAG_KEYS and i + 1 < len(lines):
                try:
                    out[key] = float(lines[i + 1].strip())
                except ValueError:
                    pass
    return out


def _run_one(task: dict) -> dict:
    raw = Path(task["raw"])
    rec = Path(task["receptor"])
    dest = Path(task["dest"])
    dest.parent.mkdir(parents=True, exist_ok=True)
    row = {"pose": task["key"], "config": task["label"], "status": "", "elapsed_s": "",
           "minimizedAffinity": "", "CNNscore": "", "CNNaffinity": "", "minimizedRMSD": "",
           "error": ""}
    if dest.exists() and dest.stat().st_size > 0 and not task["force"]:
        row["status"] = "skipped_exists"
        return row
    cmd = [task["exe"], "--receptor", str(rec), "--ligand", str(raw),
           "--autobox_ligand", str(raw), "--autobox_add", str(task["autobox"]),
           "--out", str(dest), "--cpu", "1", "--seed", str(task["seed"]),
           "--num_modes", "1",
           "--minimize" if task["search"] == "minimize" else "--local_only"]
    if not task["gpu"]:
        cmd.append("--no_gpu")
    cmd += task["extra"]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=task["timeout"])
    except subprocess.TimeoutExpired:
        row.update(status="timeout", elapsed_s=f"{time.time()-t0:.3f}")
        return row
    except Exception as e:
        row.update(status="invoke_failed", error=repr(e))
        return row
    row["elapsed_s"] = f"{time.time()-t0:.3f}"
    if proc.returncode != 0 or not dest.exists() or dest.stat().st_size == 0:
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        row.update(status="failed", error=err[-1] if err else f"rc={proc.returncode}")
        return row
    row["status"] = "ok"
    for k, v in _parse_tags(dest, proc.stdout).items():
        row[k] = f"{v:.6f}"
    return row


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pose-list", required=True, type=Path,
                    help="File of *__refRAW.sdf paths, one per line.")
    ap.add_argument("--dest-root", required=True, type=Path)
    ap.add_argument("--search", default="minimize", choices=("minimize", "local_only"))
    ap.add_argument("--extra", default="", help="Extra gnina flags, e.g. '--minimize_iters 100'.")
    ap.add_argument("--label", default=None, help="Config label recorded in the manifest.")
    ap.add_argument("--tool", default="gnina", choices=("gnina", "smina"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--autobox-add", type=float, default=4.0)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    exe = os.environ["EQ_GNINA"] if args.tool == "gnina" else os.environ["EQ_SMINA"]
    if not Path(exe).exists():
        print(f"ERROR: {args.tool} not found at {exe}", file=sys.stderr)
        return 2
    label = args.label or f"{args.search}{(' ' + args.extra) if args.extra else ''}"
    extra = shlex.split(args.extra)

    poses = [Path(l.strip()) for l in args.pose_list.read_text().splitlines() if l.strip()]
    tasks, missing = [], 0
    for raw in poses:
        rec = _prepared_receptor(raw)
        if rec is None or not raw.is_file():
            missing += 1
            continue
        key = f"{raw.parent.name}/{raw.name[:-len(RAW_SUFFIX)]}"
        tasks.append({"raw": str(raw), "receptor": str(rec), "key": key, "label": label,
                      "dest": str(args.dest_root / raw.parent.name / raw.name.replace(RAW_SUFFIX, ".sdf")),
                      "search": args.search, "extra": extra, "exe": exe, "seed": args.seed,
                      "autobox": args.autobox_add, "timeout": args.timeout,
                      "gpu": args.gpu, "force": args.force})
    print(f"config      : {label}")
    print(f"tool/exe    : {args.tool} {exe}   gpu={args.gpu}  seed={args.seed}")
    print(f"poses       : {len(tasks)} runnable ({missing} unresolved)")
    if not tasks:
        return 1

    args.dest_root.mkdir(parents=True, exist_ok=True)
    man = args.dest_root / "sweep_manifest.csv"
    fields = ["pose", "config", "status", "elapsed_s", "minimizedAffinity",
              "CNNscore", "CNNaffinity", "minimizedRMSD", "error"]
    t0 = time.time(); done = 0; counts: Dict[str, int] = {}
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_run_one, t): t for t in tasks}
            for f in as_completed(futs):
                try:
                    r = f.result()
                except Exception as e:
                    r = {"pose": futs[f]["key"], "config": label, "status": "crash",
                         "error": repr(e), "elapsed_s": ""}
                w.writerow({k: r.get(k, "") for k in fields})
                counts[r["status"]] = counts.get(r["status"], 0) + 1
                done += 1
                if done % 250 == 0 or done == len(tasks):
                    el = time.time() - t0
                    print(f"  {done}/{len(tasks)}  {el/60:.1f} min  {done/el:.1f} poses/s", flush=True)
                    fh.flush()
    (args.dest_root / "sweep_summary.json").write_text(json.dumps(
        {"label": label, "search": args.search, "extra": args.extra, "tool": args.tool,
         "seed": args.seed, "gpu": args.gpu, "n": len(tasks), "status_counts": counts,
         "wall_s": round(time.time() - t0, 1)}, indent=2))
    print(f"\ndone in {(time.time()-t0)/60:.1f} min: {counts}")
    print(f"manifest: {man}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
