#!/usr/bin/env python3
"""Per-pose accuracy metrics for a single RAW AutoDock arm on the 308 benchmark.

Answers "how do PoseBusters validity, the <=2 A RMSD rate and the <=1 A Kabsch
RMSD rate change when the receptor is prepared with Meeko instead of MGLTools"
without running the full posebusters_pose_comparison report, which assumes a
multi-engine cohort and emits ~77 tables that a single-arm run does not need.

Every metric is computed by importing the definitions the project already uses,
so the numbers are directly comparable with per_pose_metrics.csv:

  rmsd            symmetry_rmsd()   heavy-atom RMSD, NO superposition
  pb_rmsd         posebusters_rmsd() PoseBusters check_rmsd, same frame
  pb_kabsch_rmsd  posebusters_rmsd() PoseBusters check_rmsd after superposition
  bestfit_rmsd    pb_kabsch_rmsd, falling back to rigid_body_fit() when
                  PoseBusters cannot score the pair  (posebusters_pose_comparison
                  applies exactly this preference in process_pair)
  pb_valid        AND over PB_CRITICAL_CHECKS (the canonical 22-test set)

Ranking uses autodock_rank, the immutable Vina rank; these arms carry no
gnina/smina optimizer output, so there is no optimized_rank to prefer.
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "Scripts" / "Analysis"))

from posebusters_pose_comparison import (  # noqa: E402
    PB_CRITICAL_CHECKS,
    _to_bool,
    load_first_mol,
    posebusters_rmsd,
    reassign_template,
    rigid_body_fit,
    symmetry_rmsd,
)

RDLogger.DisableLog("rdApp.*")


def _load_index(pb_csv: Path, ids: set[str] | None) -> pd.DataFrame:
    df = pd.read_csv(pb_csv, low_memory=False)
    checks = [c for c in PB_CRITICAL_CHECKS if c in df.columns]
    if not checks:
        raise SystemExit(f"{pb_csv} carries none of the canonical PoseBusters test columns")
    df["pb_valid"] = pd.DataFrame({c: _to_bool(df[c]) for c in checks}).all(axis=1)
    df = df[df["docking_method"].astype(str).str.lower() == "autodock"]
    if "optimizer" in df.columns:
        opt = df["optimizer"].astype(str).str.lower()
        df = df[opt.isin(("original", "nan", ""))]
    if ids is not None:
        df = df[df["protein"].astype(str).isin(ids)]
    return df.reset_index(drop=True)


def _score_pair(job: tuple[str, str, str, list[dict]]) -> list[dict]:
    protein, ligand, crystal_path, records = job
    crystal = load_first_mol(Path(crystal_path))
    if crystal is None:
        return []
    crystal_h = Chem.RemoveHs(crystal)
    out: list[dict] = []
    for rec in records:
        pose = load_first_mol(Path(rec["pose_file"]))
        if pose is None:
            continue
        pose = reassign_template(pose, crystal_h) or pose
        try:
            rmsd = symmetry_rmsd(pose, crystal_h)
        except Exception:
            rmsd = float("nan")
        pb_rmsd, pb_kabsch, pb_within = posebusters_rmsd(pose, crystal_h)
        try:
            _, self_bestfit = rigid_body_fit(pose, crystal_h)
        except Exception:
            self_bestfit = float("nan")
        bestfit = pb_kabsch if not math.isnan(pb_kabsch) else self_bestfit
        out.append({
            "protein": protein,
            "ligand": ligand,
            "pose_name": rec.get("pose_name"),
            "pose_file": rec["pose_file"],
            "rank": rec["rank"],
            "rmsd": rmsd,
            "pb_rmsd": pb_rmsd,
            "pb_kabsch_rmsd": pb_kabsch,
            "pb_rmsd_within_2A": pb_within,
            "bestfit_rmsd": bestfit,
            "pb_valid": bool(rec["pb_valid"]),
        })
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pb-csv", type=Path, required=True)
    ap.add_argument("--benchmark-dir", type=Path, default=Path("Data/PoseBuster Benchmark Set"))
    ap.add_argument("--ids-file", type=Path, default=None)
    ap.add_argument("--out-csv", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 4) - 2))
    args = ap.parse_args()

    ids = None
    if args.ids_file:
        ids = {ln.strip() for ln in args.ids_file.read_text().splitlines()
               if ln.strip() and not ln.lstrip().startswith("#")}

    df = _load_index(args.pb_csv, ids)
    print(f"pose rows: {len(df):,}   complexes: {df['protein'].nunique()}")

    jobs = []
    missing_crystal = []
    for (protein, ligand), sub in df.groupby(["protein", "ligand"], sort=True):
        crystal = args.benchmark_dir / str(protein) / f"{protein}_ligand.sdf"
        if not crystal.is_file():
            missing_crystal.append(protein)
            continue
        records = [{
            "pose_file": r["pose_file"],
            "pose_name": r.get("pose_name"),
            "rank": int(r["autodock_rank"]) if pd.notna(r.get("autodock_rank")) else 999,
            "pb_valid": bool(r["pb_valid"]),
        } for _, r in sub.iterrows()]
        jobs.append((str(protein), str(ligand), str(crystal), records))
    if missing_crystal:
        print(f"WARNING: no crystal ligand for {len(missing_crystal)} complexes: "
              f"{sorted(missing_crystal)[:5]}")

    from multiprocessing import Pool
    rows: list[dict] = []
    with Pool(args.workers) as pool:
        for i, res in enumerate(pool.imap_unordered(_score_pair, jobs, chunksize=1), start=1):
            rows.extend(res)
            if i % 25 == 0:
                print(f"  {i}/{len(jobs)} complexes scored", flush=True)

    out = pd.DataFrame(rows).sort_values(["protein", "rank"])
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv}  ({len(out):,} poses, {out['protein'].nunique()} complexes)")


if __name__ == "__main__":
    main()
