#!/usr/bin/env python
"""Matched-thread wall-clock re-timing for the exhaustiveness 32-vs-64 comparison.

Why this exists
---------------
The two production arms cannot be compared on their own ``elapsed_time_s``:

* the exhaustiveness-32 campaign ran on 2026-08-03 with ``--cpu 32`` (all 308 Vina
  logs record ``CPU: 32``), while
* the exhaustiveness-64 arm runs with ``--cpu 20`` (the value ``cpus_per_worker``
  holds in the config today),

so a naive ratio conflates 2x the search effort with 1.6x fewer threads, plus
whatever contention the 2026-08-03 machine was under.

Vina seeds every Monte-Carlo run deterministically from ``--seed``, so the pose
ensemble is thread-count invariant: re-running a complex at a different ``--cpu``
reproduces the pose file byte for byte. That makes a timing-only re-run safe --
it touches no scientific output and can be validated by byte-comparing its poses
against the committed ones.

This script re-docks each cohort complex at BOTH exhaustiveness values from the
*committed* exhaustiveness-32 staged inputs, on the same thread count, back to
back on an otherwise idle machine, and records the wall clock. It writes poses to
a scratch directory that is never read by the analysis.

Output: a CSV with one row per (complex, exhaustiveness) carrying elapsed_s,
best_affinity, the model count, the 1-minute load average, and -- for the
validation -- whether the pose file reproduced the committed one byte for byte.

Usage
-----
    python Scripts/Docking/time_exhaustiveness_arms.py \
        --ids-file Dockings/vina_results_full_protein_vina_scoring_exh64/exhaustiveness_sample_ids.txt \
        --out Dockings/vina_results_full_protein_vina_scoring_exh64/timing_recalibration.csv
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# The pinned build. Bare `vina` on PATH is /usr/bin/vina (v1.2.5) -- a different
# scoring engine from the v1.2.7 build the campaign used.
VINA_BIN = "/home/manndo/AutoDock-Vina/build/linux/release/vina"

TREE32 = Path("Dockings/vina_results_full_protein_vina_scoring")
TREE64 = Path("Dockings/vina_results_full_protein_vina_scoring_exh64")


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def n_models(p: Path) -> int:
    return sum(1 for line in p.read_text().splitlines() if line.startswith("MODEL"))


def best_affinity(p: Path) -> float | None:
    for line in p.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT"):
            return float(line.split()[3])
    return None


def committed_pose(tree: Path, pid: str) -> Path | None:
    d = tree / pid / "mgl_tools" / "docking"
    hits = sorted(d.glob("*_vina_out.pdbqt")) if d.is_dir() else []
    return hits[0] if len(hits) == 1 else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids-file", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--cpu", type=int, default=32, help="thread count for BOTH arms")
    ap.add_argument("--exhaustiveness", type=int, nargs="+", default=[32, 64])
    ap.add_argument("--num-modes", type=int, default=30)
    ap.add_argument("--energy-range", type=int, default=6)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--scratch", type=Path, default=None)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args(argv)

    if not os.access(VINA_BIN, os.X_OK):
        raise RuntimeError(f"Vina binary not executable: {VINA_BIN}")

    ids = [ln.strip() for ln in args.ids_file.read_text().splitlines()
           if ln.strip() and not ln.lstrip().startswith("#")]
    if args.limit:
        ids = ids[: args.limit]

    scratch = args.scratch or Path(tempfile.mkdtemp(prefix="exh_timing_"))
    scratch.mkdir(parents=True, exist_ok=True)

    version = subprocess.run([VINA_BIN, "--version"], text=True,
                             capture_output=True).stdout.strip().splitlines()[0]
    print(f"Vina: {version}")
    print(f"Threads: {args.cpu} (identical for every arm)   complexes: {len(ids)}")
    print(f"Scratch: {scratch}\n")

    # Resume: keep rows already measured.
    done: set[tuple[str, int]] = set()
    rows: list[dict] = []
    if args.out.exists():
        with args.out.open() as fh:
            for r in csv.DictReader(fh):
                rows.append(r)
                done.add((r["pdb_id"], int(r["exhaustiveness"])))
        print(f"Resuming — {len(done)} measurements already recorded\n")

    fields = ["pdb_id", "exhaustiveness", "cpu", "elapsed_s", "best_affinity",
              "n_models", "load_1min", "reproduces_committed", "committed_arm", "vina_version"]

    for i, pid in enumerate(ids, 1):
        stage = TREE32 / pid / "_staging"
        rec = stage / "receptors" / "pdbqt" / f"{pid}_protein_mgl_tools.pdbqt"
        box = stage / "receptors" / "pdbqt" / f"{pid}_protein_mgl_tools.box.txt"
        lig = stage / "ligands" / "pdbqt" / f"{pid}_ligand_start_conf.pdbqt"
        for p in (rec, box, lig):
            if not p.exists():
                raise RuntimeError(f"{pid}: missing committed staged input {p}")

        for E in args.exhaustiveness:
            if (pid, E) in done:
                continue
            out_pose = scratch / f"{pid}_e{E}.pdbqt"
            cmd = [VINA_BIN, "--receptor", str(rec), "--config", str(box), "--ligand", str(lig),
                   "--exhaustiveness", str(E), "--num_modes", str(args.num_modes),
                   "--energy_range", str(args.energy_range), "--seed", str(args.seed),
                   "--cpu", str(args.cpu), "--out", str(out_pose)]
            t0 = time.perf_counter()
            proc = subprocess.run(cmd, text=True, capture_output=True)
            elapsed = time.perf_counter() - t0
            if proc.returncode != 0 or not out_pose.exists():
                raise RuntimeError(f"{pid} E={E} failed: {proc.stderr[-800:]}")

            # Byte-identity against the committed pose of the matching arm proves
            # the production run is scientifically equivalent to this one.
            arm_tree = TREE32 if E == 32 else TREE64
            committed = committed_pose(arm_tree, pid)
            reproduces = ""
            if committed is not None:
                reproduces = str(sha256(out_pose) == sha256(committed))

            row = {"pdb_id": pid, "exhaustiveness": E, "cpu": args.cpu,
                   "elapsed_s": round(elapsed, 3), "best_affinity": best_affinity(out_pose),
                   "n_models": n_models(out_pose), "load_1min": round(os.getloadavg()[0], 2),
                   "reproduces_committed": reproduces,
                   "committed_arm": str(arm_tree) if committed is not None else "",
                   "vina_version": version}
            rows.append(row)
            print(f"[{i}/{len(ids)}] {pid} E={E:<3} {elapsed:7.1f}s  "
                  f"best={row['best_affinity']}  models={row['n_models']}  "
                  f"reproduces={reproduces or 'n/a'}")

            # Checkpoint after every measurement.
            with args.out.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=fields)
                w.writeheader()
                for r in rows:
                    w.writerow({k: r.get(k, "") for k in fields})

    print(f"\nWrote {len(rows)} measurements → {args.out}")
    if args.scratch is None:
        shutil.rmtree(scratch, ignore_errors=True)
        print(f"Removed scratch {scratch}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
