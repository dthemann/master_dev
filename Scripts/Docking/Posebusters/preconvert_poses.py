#!/usr/bin/env python3
"""Warm a PoseBusters run's PDBQT->SDF conversion cache in parallel.

run_posebusters.py converts every multi-model Vina PDBQT to per-model SDFs
inside collect_all_pose_files(), which is single-threaded: on a few thousand
poses that serial pass dominates the wall clock, because Meeko's mk_export is a
subprocess per model. The conversion is a pure function of (pdbqt, template,
converter fingerprint) and publishes a content-verified manifest, so it can be
done ahead of time from several processes and the real run then hits cache.

Run this against the SAME config the run will use, with no run in flight -- two
processes converting into one converted_pdbqt/ race on the same target paths and
the loser dies with "Converted only N/N models".
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_posebusters as R  # noqa: E402

_CTX = None


def _init(config_path: str) -> None:
    global _CTX
    R._apply_quiet_mode(True)
    R._limit_worker_threads()
    _CTX = R.load_config(config_path)


def _convert(pdbqt: str) -> tuple[str, int, str]:
    try:
        out = R._convert_pdbqt_to_sdf(pdbqt, _CTX.converted_dir, _CTX)
        return pdbqt, len(out), ""
    except Exception as exc:  # a failure here is a real data problem: report it
        return pdbqt, 0, f"{type(exc).__name__}: {exc}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", "-c", required=True)
    ap.add_argument("--workers", "-j", type=int,
                    default=max(1, (os.cpu_count() or 4) - 2))
    args = ap.parse_args()

    ctx = R.load_config(args.config)
    sources: list[str] = []
    for root in ctx.docking_directories.values():
        sources.extend(str(p) for p in sorted(Path(root).glob("**/*_vina_out.pdbqt")))
    print(f"{len(sources)} PDBQT sources -> {ctx.converted_dir}  ({args.workers} workers)",
          flush=True)

    from multiprocessing import Pool
    done = models = 0
    failures: list[tuple[str, str]] = []
    with Pool(args.workers, initializer=_init, initargs=(args.config,)) as pool:
        for pdbqt, n, err in pool.imap_unordered(_convert, sources, chunksize=1):
            done += 1
            models += n
            if err:
                failures.append((pdbqt, err))
            if done % 20 == 0:
                print(f"  {done}/{len(sources)} sources, {models:,} models, "
                      f"{len(failures)} failed", flush=True)
    print(f"DONE {done}/{len(sources)} sources, {models:,} models converted, "
          f"{len(failures)} failed")
    for pdbqt, err in failures[:20]:
        print(f"  FAIL {Path(pdbqt).name}: {err}")


if __name__ == "__main__":
    main()
