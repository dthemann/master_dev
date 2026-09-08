#!/usr/bin/env python
"""Intention-to-treat re-score of the five complexes dropped for missing DiffDock
output (plan v2 Phase 2.6, App. I :1074).

The five ids (7B2C_TP7, 7D6O_MTE, 7FRX_O88, 7M31_TDR, 8F4J_PHO) have no DiffDock
output and were removed from the 303-complex cohort. The thesis reports a
sensitivity check in which they are kept and scored as DiffDock failures while
the AutoDock exh128 + gnina poses are scored as produced. No registered generator
existed for that check; this script is it.

Engine: the hub's own loaders and RMSD (posebusters_pose_comparison.load_first_mol,
load_all_mols, reassign_template, symmetry_rmsd, _reference_copy_index,
_nearest_copy_index), so the instance value reproduces the hub's ``rmsd`` under
``--reference-convention instance`` and the nearest value reproduces it under
``nearest``. Poses are ranked by the CNNaffinity SDF tag, descending, which is
the gnina rank the hub writes as ``rank`` for the *_gnina AutoDock arms.

Validity: the canonical PoseBusters CSV holds no rows for the five (it is the
303-cohort file), so no PB-valid gate is applied here. The RMSD result makes the
gate immaterial for the counts (no pose is within 2 A at rank-1 or top-15 under
either convention). Pass --pb-csv to join a pb_valid column when one exists.

Outputs (both inside --out-dir, which must be NEW or empty):
  itt_dropped_poses.csv      one row per scored pose, both conventions
  itt_dropped_complexes.csv  one row per complex: rank-1, top-15 min, all-pose min
  itt_summary.json           ITT counts on the 308 basis, margins, provenance

Example:
  /home/manndo/anaconda3/envs/vina/bin/python \
      Scripts/Analysis/nearest_copy_program/itt_dropped_complexes.py \
      --out-dir posebusters_results/itt_nearest \
      --per-pose-csv posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report_nearest/per_pose_metrics.csv
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/manndo/master_dev")
ANA = ROOT / "Scripts" / "Analysis"
sys.path.insert(0, str(ANA))

from rdkit import Chem, RDLogger  # noqa: E402

RDLogger.DisableLog("rdApp.*")
import posebusters_pose_comparison as ppc  # noqa: E402

DROPPED_IDS = ["7B2C_TP7", "7D6O_MTE", "7FRX_O88", "7M31_TDR", "8F4J_PHO"]
DEFAULT_AUTODOCK_ROOT = "Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128"
DEFAULT_BENCHMARK_DIR = "Data/PoseBuster Benchmark Set"
AUTODOCK_ARM = "autodock_mgltools_exh128_gnina"
DIFFDOCK_ARM = "diffdock_smina"
# 303-basis counts (autodock_rank1, diffdock_rank1, autodock_top15, diffdock_top15)
DEFAULT_COUNTS_NEAREST = (149, 132, 213, 179)
DEFAULT_COUNTS_INSTANCE = (111, 103, 198, 167)
RANK_RE = re.compile(r"_rank(\d+)_gnina\.sdf$")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def prepare_out_dir(out_dir: Path) -> None:
    """The output directory must be new or empty. Nothing existing is touched."""
    if out_dir.exists():
        if not out_dir.is_dir():
            sys.exit(f"ERROR: --out-dir {out_dir} exists and is not a directory")
        if any(out_dir.iterdir()):
            sys.exit(f"ERROR: --out-dir {out_dir} exists and is not empty; choose a NEW directory")
    else:
        out_dir.mkdir(parents=True)


def load_copies(bench_dir: Path, pid: str):
    """Mirror of the hub's copy construction in process_pair: every record of
    <ID>_ligands.sdf with bond orders from the reference, the reference record
    REPLACED by the <ID>_ligand.sdf mol itself."""
    cdir = bench_dir / pid
    crystal = ppc.load_first_mol(cdir / f"{pid}_ligand.sdf")
    if crystal is None:
        sys.exit(f"ERROR: cannot read {cdir / f'{pid}_ligand.sdf'}")
    crystal_h = Chem.RemoveHs(crystal)
    copies_sdf = cdir / f"{pid}_ligands.sdf"
    copies_raw = ppc.load_all_mols(copies_sdf) if copies_sdf.exists() else []
    if not copies_raw:
        copies_raw = [crystal]
    copies_h = [Chem.RemoveHs(ppc.reassign_template(m, crystal_h) or m) for m in copies_raw]
    ref_idx = ppc._reference_copy_index(copies_h, crystal_h)
    if ref_idx is None:
        print(f"  [{pid}] WARNING: no record matches {pid}_ligand.sdf; prepending the reference as copy 0")
        copies_h = [crystal_h] + copies_h
        ref_idx = 0
    else:
        copies_h[ref_idx] = crystal_h
    return crystal_h, copies_h, ref_idx


def list_gnina_poses(autodock_root: Path, pid: str) -> list[Path]:
    d = autodock_root / pid / "mgl_tools" / "docking" / "optimized_gnina"
    if not d.is_dir():
        sys.exit(f"ERROR: no optimized_gnina directory for {pid}: {d}")
    files = sorted(p for p in d.glob("*_rank*_gnina.sdf") if RANK_RE.search(p.name))
    if not files:
        sys.exit(f"ERROR: no *_rank*_gnina.sdf poses under {d}")
    return files


def read_cnn_affinity(sdf: Path) -> float:
    """CNNaffinity tag from the untouched supplier record (before any template
    reassignment, which drops properties)."""
    for m in Chem.SDMolSupplier(str(sdf), removeHs=False, sanitize=False):
        if m is None:
            continue
        if m.HasProp("CNNaffinity"):
            try:
                return float(m.GetProp("CNNaffinity"))
            except (TypeError, ValueError):
                return float("nan")
        return float("nan")
    return float("nan")


def centroid(mol: Chem.Mol) -> np.ndarray:
    return Chem.RemoveHs(mol).GetConformer().GetPositions().mean(axis=0)


def score_complex(pid: str, autodock_root: Path, bench_dir: Path) -> list[dict]:
    crystal_h, copies_h, ref_idx = load_copies(bench_dir, pid)
    copy_cent = [centroid(c) for c in copies_h]
    files = list_gnina_poses(autodock_root, pid)

    rows = []
    for f in files:
        vina_rank = int(RANK_RE.search(f.name).group(1))
        cnn = read_cnn_affinity(f)
        pose = ppc.load_first_mol(f)
        if pose is None:
            rows.append(dict(protein=pid, pose_file=str(f), vina_rank=vina_rank, cnn_affinity=cnn,
                             load_failed=True))
            continue
        pose = ppc.reassign_template(pose, crystal_h) or pose
        r_all = [ppc._safe_symmetry_rmsd(pose, c) for c in copies_h]
        j_star = ppc._nearest_copy_index(r_all, ref_idx)
        pc = centroid(pose)
        cd_all = [float(np.linalg.norm(pc - c)) for c in copy_cent]
        alt = [j for j in range(len(copies_h)) if j != ref_idx]
        rows.append(dict(
            protein=pid, pose_file=str(f), vina_rank=vina_rank, cnn_affinity=cnn, load_failed=False,
            n_copies=len(copies_h), ref_copy_index=ref_idx,
            rmsd_ref_instance=r_all[ref_idx], rmsd_nearest=r_all[j_star],
            nearest_copy_index=j_star, nearest_copy_is_ref=bool(j_star == ref_idx),
            rmsd_min_alternate=(float(np.nanmin([r_all[j] for j in alt])) if alt else float("nan")),
            centroid_dist_ref_instance=cd_all[ref_idx], centroid_dist_nearest=cd_all[j_star],
        ))
    df = pd.DataFrame(rows)
    if df["load_failed"].any():
        sys.exit(f"ERROR: {int(df.load_failed.sum())} pose(s) of {pid} could not be loaded")
    if df["cnn_affinity"].isna().any():
        sys.exit(f"ERROR: {int(df.cnn_affinity.isna().sum())} pose(s) of {pid} carry no CNNaffinity tag")
    # gnina rank: CNNaffinity descending, ties broken by the Vina rank (ascending)
    df = df.sort_values(["cnn_affinity", "vina_rank"], ascending=[False, True]).reset_index(drop=True)
    df["gnina_rank"] = np.arange(1, len(df) + 1)

    # cross-check against the optimiser log when present (informational only)
    log = autodock_root / pid / "mgl_tools" / "optimization_log.csv"
    df["log_optimized_rank"] = np.nan
    if log.exists():
        try:
            lg = pd.read_csv(log, usecols=["optimized_file", "optimized_rank"])
            lg["pose_file"] = lg["optimized_file"].astype(str)
            m = df.merge(lg[["pose_file", "optimized_rank"]], on="pose_file", how="left")
            df["log_optimized_rank"] = m["optimized_rank"].to_numpy()
            mism = int((df["log_optimized_rank"].notna()
                        & (df["log_optimized_rank"] != df["gnina_rank"])).sum())
            if mism:
                print(f"  [{pid}] WARNING: {mism} pose(s) where the CNNaffinity order differs from "
                      f"optimization_log.csv optimized_rank (ties or log drift); the SDF tag order is used")
        except Exception as exc:  # noqa: BLE001
            print(f"  [{pid}] WARNING: optimization_log.csv unreadable ({exc}); no cross-check")
    return df.to_dict("records")


def per_complex_summary(poses: pd.DataFrame, top_n: int, thr: float) -> pd.DataFrame:
    out = []
    for pid, g in poses.groupby("protein", sort=True):
        g = g.sort_values("gnina_rank")
        r1 = g[g.gnina_rank == 1].iloc[0]
        top = g[g.gnina_rank <= top_n]
        best_near_all = g.loc[g.rmsd_nearest.idxmin()]
        best_inst_all = g.loc[g.rmsd_ref_instance.idxmin()]
        row = dict(
            protein=pid, n_poses=int(len(g)), n_copies=int(g.n_copies.iloc[0]),
            ref_copy_index=int(g.ref_copy_index.iloc[0]),
            rank1_rmsd_ref_instance=float(r1.rmsd_ref_instance),
            rank1_rmsd_nearest=float(r1.rmsd_nearest),
            rank1_nearest_copy_index=int(r1.nearest_copy_index),
            rank1_pb_valid=(bool(r1.pb_valid) if "pb_valid" in g.columns and pd.notna(r1.get("pb_valid")) else None),
            top15_min_rmsd_ref_instance=float(top.rmsd_ref_instance.min()),
            top15_min_rmsd_nearest=float(top.rmsd_nearest.min()),
            all_min_rmsd_ref_instance=float(g.rmsd_ref_instance.min()),
            all_min_rmsd_nearest=float(g.rmsd_nearest.min()),
            all_min_rmsd_nearest_rank=int(best_near_all.gnina_rank),
            all_min_rmsd_nearest_copy_index=int(best_near_all.nearest_copy_index),
            all_min_rmsd_ref_instance_rank=int(best_inst_all.gnina_rank),
        )
        for conv, col in (("nearest", "rmsd_nearest"), ("instance", "rmsd_ref_instance")):
            row[f"within{thr:g}_rank1_{conv}"] = bool(r1[col] <= thr)
            row[f"within{thr:g}_top{top_n}_{conv}"] = bool((top[col] <= thr).any())
            row[f"within{thr:g}_all_{conv}"] = bool((g[col] <= thr).any())
        out.append(row)
    return pd.DataFrame(out)


def effective_rank(df: pd.DataFrame) -> pd.Series:
    """Hub rank, except EquiBind whose ``rank`` is a 999 sentinel: rank by the
    refinement affinity (gnina for *_gnina, smina otherwise), ascending."""
    rank = pd.to_numeric(df["rank"], errors="coerce").astype(float)
    eq = df["method"].str.startswith("equibind")
    if eq.any():
        e = df[eq]
        key = pd.to_numeric(e.get("gnina_affinity"), errors="coerce")
        if "smina_affinity" in e.columns:
            key = key.where(e["method"].str.endswith("_gnina"),
                            pd.to_numeric(e["smina_affinity"], errors="coerce"))
        order = e.assign(_k=key).sort_values("_k", kind="mergesort").groupby(["method", "protein"]).cumcount() + 1
        rank.loc[eq] = order.reindex(e.index).astype(float)
    return rank


def counts_from_table(per_pose_csv: Path, top_n: int, thr: float) -> dict:
    """303-basis recovery counts for both arms and both conventions from a
    nearest-convention hub table (PB-valid and rmsd <= thr, existence over depth)."""
    usecols = ["method", "protein", "rank", "rmsd", "rmsd_ref_instance", "pb_valid",
               "gnina_affinity", "smina_affinity", "reference_convention"]
    header = pd.read_csv(per_pose_csv, nrows=0).columns
    missing = [c for c in ("rmsd_ref_instance", "reference_convention") if c not in header]
    if missing:
        sys.exit(f"ERROR: --per-pose-csv lacks {missing}; a nearest-convention hub table is required")
    df = pd.read_csv(per_pose_csv, usecols=[c for c in usecols if c in header], low_memory=False)
    conv = set(df["reference_convention"].dropna().unique())
    if conv != {"nearest"}:
        sys.exit(f"ERROR: --per-pose-csv reference_convention is {sorted(conv)}, expected ['nearest']")
    df = df[df["method"].isin([AUTODOCK_ARM, DIFFDOCK_ARM])].copy()
    df["eff_rank"] = effective_rank(df)
    valid = df["pb_valid"].astype(str).str.lower().isin(["true", "1", "1.0"])
    n_ids = int(df["protein"].nunique())
    out = {"n_complexes": n_ids, "source": str(per_pose_csv)}
    for conv_name, col in (("nearest", "rmsd"), ("instance", "rmsd_ref_instance")):
        vals = []
        for depth in (1, top_n):
            for arm in (AUTODOCK_ARM, DIFFDOCK_ARM):
                m = (df["method"] == arm) & (df["eff_rank"] <= depth) & valid & (df[col] <= thr)
                vals.append(int(df.loc[m, "protein"].nunique()))
        # order: autodock_rank1, diffdock_rank1, autodock_top15, diffdock_top15
        out[conv_name] = tuple(vals)
    return out


def itt_block(counts_303: tuple, summary: pd.DataFrame, conv: str, top_n: int, thr: float,
              n_total: int) -> dict:
    a1, d1, a15, d15 = counts_303
    n_303 = n_total - len(summary)
    add_r1 = int(summary[f"within{thr:g}_rank1_{conv}"].sum())
    add_t15 = int(summary[f"within{thr:g}_top{top_n}_{conv}"].sum())
    add_all = int(summary[f"within{thr:g}_all_{conv}"].sum())
    blk = {
        "counts_303": {"autodock_rank1": a1, "diffdock_rank1": d1,
                       f"autodock_top{top_n}": a15, f"diffdock_top{top_n}": d15, "n": n_303},
        "dropped_recovered": {"rank1": add_r1, f"top{top_n}": add_t15, "all_poses": add_all},
        "counts_308": {"autodock_rank1": a1 + add_r1, "diffdock_rank1": d1,
                       f"autodock_top{top_n}": a15 + add_t15, f"diffdock_top{top_n}": d15, "n": n_total},
        "margin_pp_303": {"rank1": round((a1 - d1) / n_303 * 100, 2),
                          f"top{top_n}": round((a15 - d15) / n_303 * 100, 2)},
        "margin_pp_308": {"rank1": round((a1 + add_r1 - d1) / n_total * 100, 2),
                          f"top{top_n}": round((a15 + add_t15 - d15) / n_total * 100, 2)},
        "rate_pct_308": {"autodock_rank1": round((a1 + add_r1) / n_total * 100, 1),
                         "diffdock_rank1": round(d1 / n_total * 100, 1),
                         f"autodock_top{top_n}": round((a15 + add_t15) / n_total * 100, 1),
                         f"diffdock_top{top_n}": round(d15 / n_total * 100, 1)},
        "closest_rank1_A": round(float(summary[f"rank1_rmsd_{'nearest' if conv == 'nearest' else 'ref_instance'}"].min()), 2),
        f"closest_top{top_n}_A": round(float(summary[f"top{top_n}_min_rmsd_{'nearest' if conv == 'nearest' else 'ref_instance'}"].min()), 2),
        "closest_all_A": round(float(summary[f"all_min_rmsd_{'nearest' if conv == 'nearest' else 'ref_instance'}"].min()), 2),
    }
    col = "nearest" if conv == "nearest" else "ref_instance"
    i_r1 = summary[f"rank1_rmsd_{col}"].idxmin()
    i_t = summary[f"top{top_n}_min_rmsd_{col}"].idxmin()
    i_all = summary[f"all_min_rmsd_{col}"].idxmin()
    blk["closest_rank1_id"] = str(summary.loc[i_r1, "protein"])
    blk[f"closest_top{top_n}_id"] = str(summary.loc[i_t, "protein"])
    blk["closest_all_id"] = str(summary.loc[i_all, "protein"])
    blk["closest_all_rank"] = int(summary.loc[i_all, f"all_min_rmsd_{col}_rank"])
    blk["complexes_within_thr_beyond_top_n"] = [
        str(p) for p, a, t in zip(summary.protein, summary[f"within{thr:g}_all_{conv}"],
                                  summary[f"within{thr:g}_top{top_n}_{conv}"]) if a and not t]
    return blk


# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ids", nargs="+", default=DROPPED_IDS)
    ap.add_argument("--autodock-root", default=DEFAULT_AUTODOCK_ROOT, help="read only")
    ap.add_argument("--benchmark-dir", default=DEFAULT_BENCHMARK_DIR, help="read only")
    ap.add_argument("--out-dir", required=True, help="NEW (or empty) directory for the outputs")
    ap.add_argument("--counts-303", nargs=4, type=int, default=list(DEFAULT_COUNTS_NEAREST),
                    metavar=("AD_R1", "DD_R1", "AD_T15", "DD_T15"),
                    help="303-basis counts under the NEAREST convention: autodock_rank1 diffdock_rank1 "
                         "autodock_top15 diffdock_top15 (default 149 132 213 179)")
    ap.add_argument("--counts-303-instance", nargs=4, type=int, default=list(DEFAULT_COUNTS_INSTANCE),
                    metavar=("AD_R1", "DD_R1", "AD_T15", "DD_T15"),
                    help="the same under the single-instance convention (default 111 103 198 167)")
    ap.add_argument("--per-pose-csv", default=None,
                    help="optional nearest-convention hub table; when given the 303-basis counts are "
                         "DERIVED from it for both conventions and the passed values only cross-checked")
    ap.add_argument("--pb-csv", default=None,
                    help="optional PoseBusters CSV carrying the five ids (joins pb_valid on pose_file)")
    ap.add_argument("--top-n", type=int, default=15)
    ap.add_argument("--threshold", type=float, default=2.0)
    ap.add_argument("--n-total", type=int, default=308)
    args = ap.parse_args(argv)

    os.chdir(ROOT)
    out_dir = Path(args.out_dir)
    prepare_out_dir(out_dir)
    autodock_root, bench_dir = Path(args.autodock_root), Path(args.benchmark_dir)

    rows: list[dict] = []
    for pid in args.ids:
        print(f"[{pid}] scoring")
        rows.extend(score_complex(pid, autodock_root, bench_dir))
    poses = pd.DataFrame(rows)

    validity_note = ("not applied: the canonical posebusters_filtered_results.csv is the 303-cohort file and "
                     "holds no rows for the five dropped ids")
    if args.pb_csv:
        pb = pd.read_csv(args.pb_csv, low_memory=False)
        if {"pose_file", "pb_valid"} <= set(pb.columns):
            pb = pb[["pose_file", "pb_valid"]].drop_duplicates("pose_file")
            poses = poses.merge(pb, on="pose_file", how="left")
            validity_note = f"joined from {args.pb_csv}: {int(poses.pb_valid.notna().sum())} of {len(poses)} poses matched"
        else:
            print("WARNING: --pb-csv lacks pose_file/pb_valid columns; validity not joined")

    summary = per_complex_summary(poses, args.top_n, args.threshold)

    counts_source = "passed or default"
    counts_nearest, counts_instance = tuple(args.counts_303), tuple(args.counts_303_instance)
    derived = None
    if args.per_pose_csv:
        derived = counts_from_table(Path(args.per_pose_csv), args.top_n, args.threshold)
        for conv, passed in (("nearest", counts_nearest), ("instance", counts_instance)):
            if tuple(derived[conv]) != tuple(passed):
                print(f"WARNING: {conv} 303-basis counts derived from the table {derived[conv]} differ from "
                      f"the passed/default {passed}; the DERIVED values are used")
        counts_nearest, counts_instance = tuple(derived["nearest"]), tuple(derived["instance"])
        counts_source = f"derived from {args.per_pose_csv} ({derived['n_complexes']} complexes)"

    blocks = {"nearest": itt_block(counts_nearest, summary, "nearest", args.top_n, args.threshold, args.n_total),
              "instance": itt_block(counts_instance, summary, "instance", args.top_n, args.threshold, args.n_total)}

    poses_out = out_dir / "itt_dropped_poses.csv"
    summ_out = out_dir / "itt_dropped_complexes.csv"
    json_out = out_dir / "itt_summary.json"
    poses.sort_values(["protein", "gnina_rank"]).to_csv(poses_out, index=False)
    summary.to_csv(summ_out, index=False)
    result = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(__file__).resolve()),
        "ids": list(args.ids),
        "n_poses_scored": int(len(poses)),
        "autodock_arm": AUTODOCK_ARM, "diffdock_arm": DIFFDOCK_ARM,
        "ranking": "CNNaffinity SDF tag descending (ties by Vina rank)",
        "threshold_A": args.threshold, "top_n": args.top_n, "n_total": args.n_total,
        "validity_gate": validity_note,
        "counts_303_source": counts_source,
        "counts_303_passed": {"nearest": list(args.counts_303), "instance": list(args.counts_303_instance)},
        "counts_303_derived": ({k: list(v) if isinstance(v, tuple) else v for k, v in derived.items()}
                               if derived else None),
        "autodock_root": str(autodock_root), "benchmark_dir": str(bench_dir),
        "conventions": blocks,
    }
    json_out.write_text(json.dumps(result, indent=2))

    print()
    print(summary[["protein", "n_copies", "rank1_rmsd_ref_instance", "rank1_rmsd_nearest",
                   f"top{args.top_n}_min_rmsd_ref_instance", f"top{args.top_n}_min_rmsd_nearest",
                   "all_min_rmsd_ref_instance", "all_min_rmsd_nearest", "all_min_rmsd_nearest_rank"]]
          .round(2).to_string(index=False))
    for conv in ("nearest", "instance"):
        b = blocks[conv]
        print(f"\n[{conv}] 303: {b['counts_303']}  dropped recovered: {b['dropped_recovered']}")
        print(f"[{conv}] 308: {b['counts_308']}  margins pp 308: {b['margin_pp_308']}  (303: {b['margin_pp_303']})")
        print(f"[{conv}] closest rank-1 {b['closest_rank1_A']} A ({b['closest_rank1_id']}), "
              f"top-{args.top_n} {b[f'closest_top{args.top_n}_A']} A ({b[f'closest_top{args.top_n}_id']}), "
              f"all {b['closest_all_A']} A ({b['closest_all_id']} rank {b['closest_all_rank']})")
    print(f"\nwrote {poses_out}\nwrote {summ_out}\nwrote {json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
