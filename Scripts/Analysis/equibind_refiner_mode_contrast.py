#!/usr/bin/env python3
"""Contrast the EquiBind --local_only and --minimize refiner modes on PB validity.

Answers examiner item M5: the published EquiBind arm was refined with gnina/smina
in local-search mode while AutoDock and DiffDock used energy minimisation, so RQ2
compares arms that did not receive the same operation. This script quantifies
what the matched operation would have given.

The two tables are joined on (pose_name, refine_variant). That is only valid
because the minimize tree was built by re-refining the committed __refRAW poses
in place, preserving file names exactly. The Orai arms needed a full pipeline
replay, which permutes pose numbering, and must NOT be joined this way.

TWO ESTIMATES ARE REPORTED, because the refiners fail on different poses:
  * intention-to-treat -- every pose. A failed refinement leaves the unrefined
    pose under the refined label, and raw EquiBind poses are almost never valid,
    so an arm with more failures is penalised. This is what a reader of the
    thesis table would get.
  * per-protocol -- only poses where BOTH refiners succeeded. This isolates the
    effect of the flag from the effect of one refiner crashing or timing out.
If the two disagree materially, the difference is itself a finding.
"""
from __future__ import annotations

import argparse
from math import comb
from pathlib import Path

import pandas as pd

CANONICAL = [
    "mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "protein-ligand_maximum_distance", "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors", "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters", "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors", "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
]
KEEP = ["pose_name", "protein", "refine_variant", "refine_succeeded"] + CANONICAL


def _as_bool(s: pd.Series) -> pd.Series:
    """Object-dtype columns round-trip through NaN; fill before casting."""
    return s.astype("object").where(s.notna(), False).astype(bool)


def _load(path: Path, unguided_only: bool = True) -> pd.DataFrame:
    frames = []
    for chunk in pd.read_csv(path, chunksize=200_000, low_memory=False):
        c = chunk[chunk["docking_method"] == "equibind_guided"]
        if unguided_only:
            c = c[c["pose_name"].str.contains("/unguided_", regex=False)]
        c = c[c["refine_variant"].isin(["smina", "gnina"])]
        cols = [k for k in KEEP if k in c.columns]
        if len(c):
            frames.append(c[cols])
    df = pd.concat(frames, ignore_index=True)
    tests = [t for t in CANONICAL if t in df.columns]
    df["pb_valid"] = df[tests].fillna(False).astype(bool).all(axis=1)
    return df


def _mcnemar_exact(b01: int, b10: int) -> float:
    n = b01 + b10
    if n == 0:
        return 1.0
    lo = min(b01, b10)
    return min(1.0, 2 * sum(comb(n, k) for k in range(lo + 1)) / 2 ** n)


def _contrast(a: pd.Series, b: pd.Series, label: str) -> None:
    """a = local_only, b = minimize, aligned index."""
    n = len(a)
    b01 = int((~a & b).sum())
    b10 = int((a & ~b).sum())
    p = _mcnemar_exact(b01, b10)
    print(f"  {label:<16s} n={n:6d}   "
          f"local_only {100*a.mean():5.1f}%  ->  minimize {100*b.mean():5.1f}%   "
          f"delta {100*(b.mean()-a.mean()):+5.1f} pp")
    print(f"  {'':<16s} gained {b01:5d}   lost {b10:5d}   "
          f"McNemar exact p = {p:.3g}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--published", type=Path, default=Path(
        "posebusters_results/benchmark_full_protein_vina_scoring/dock/"
        "posebusters_filtered_results.csv"))
    ap.add_argument("--minimize", type=Path, default=Path(
        "posebusters_results/benchmark_equibind_minimize/dock/"
        "posebusters_filtered_results.csv"))
    args = ap.parse_args()

    print("loading published (--local_only) ...")
    L = _load(args.published)
    print("loading minimize ...")
    M = _load(args.minimize)
    print(f"  local_only rows {len(L)}   minimize rows {len(M)}")

    key = ["pose_name", "refine_variant"]
    merged = L.merge(M, on=key, suffixes=("_loc", "_min"))
    print(f"  joined on (pose_name, refine_variant): {len(merged)} pairs")
    only_l = len(L) - len(merged)
    only_m = len(M) - len(merged)
    if only_l or only_m:
        print(f"  NOTE unmatched: {only_l} local_only-only, {only_m} minimize-only "
              f"(cohort differs: published is the 303-complex intersection, "
              f"the new run covers 308)")

    for variant in ("gnina", "smina"):
        sub = merged[merged["refine_variant"] == variant]
        if not len(sub):
            continue
        print(f"\n=== EquiBind + {variant}, unguided ===")
        _contrast(sub["pb_valid_loc"], sub["pb_valid_min"], "intention-to-treat")

        if {"refine_succeeded_loc", "refine_succeeded_min"} <= set(sub.columns):
            ok = (_as_bool(sub["refine_succeeded_loc"])
                  & _as_bool(sub["refine_succeeded_min"]))
            pp = sub[ok]
            if len(pp):
                _contrast(pp["pb_valid_loc"], pp["pb_valid_min"], "per-protocol")
                print(f"  {'':<16s} excluded {len(sub)-len(pp)} pose(s) where a "
                      f"refiner failed in either arm")

        print("  per-test FAILURE rate (%):")
        rows = []
        for t in CANONICAL:
            cl, cm = f"{t}_loc", f"{t}_min"
            if cl in sub.columns and cm in sub.columns:
                fl = 100 * (~sub[cl].fillna(False).astype(bool)).mean()
                fm = 100 * (~sub[cm].fillna(False).astype(bool)).mean()
                if fl > 0.05 or fm > 0.05:
                    rows.append((t, fl, fm, fm - fl))
        for t, fl, fm, d in sorted(rows, key=lambda r: -abs(r[3])):
            print(f"    {t:<40s} {fl:6.1f} -> {fm:6.1f}   {d:+6.1f}")

    # Complex-level: how many complexes have at least one valid pose.
    print("\n=== complex-level coverage (>=1 PB-valid pose) ===")
    for variant in ("gnina", "smina"):
        sub = merged[merged["refine_variant"] == variant]
        if not len(sub):
            continue
        g = sub.groupby("protein_loc")[["pb_valid_loc", "pb_valid_min"]].any()
        a, b = g["pb_valid_loc"], g["pb_valid_min"]
        b01 = int((~a & b).sum())
        b10 = int((a & ~b).sum())
        print(f"  {variant:6s} complexes n={len(g)}   "
              f"local_only {int(a.sum())}  ->  minimize {int(b.sum())}   "
              f"({b01} gained, {b10} lost, p = {_mcnemar_exact(b01, b10):.3g})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
