#!/usr/bin/env python3
"""Compare RAW AutoDock arms on the 308-complex PoseBusters benchmark.

Consumes the per-pose CSVs written by boxed308_raw_arm_metrics.py (one per
receptor-preparation arm) and reports, per arm and per ranking depth:

  valid     PoseBusters-valid                    (all 22 canonical tests pass)
  rmsd2     symmetry-corrected RMSD <= 2 A       (column `rmsd`, no superposition)
  kabsch1   best-fit RMSD <= 1 A                 (column `bestfit_rmsd`)
  headline  rmsd <= 2 A AND PoseBusters-valid

Arms are compared on the complexes both cover, so the Meeko-vs-MGLTools contrast
is paired; McNemar's exact test and a Newcombe interval for the difference of the
two correlated proportions accompany every depth.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "Scripts" / "Analysis"))

from stats_utils import mcnemar_exact, newcombe_paired_diff_ci, wilson_ci  # noqa: E402

DEPTHS = [1, 5, 15, 30, 0]          # 0 = every pose the arm produced
CRITERIA = ["valid", "rmsd2", "kabsch1", "headline"]


def _load(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    for col in ("rmsd", "bestfit_rmsd"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["pb_valid"] = df["pb_valid"].map(lambda x: str(x).strip().lower() in ("true", "1"))
    df["rank"] = pd.to_numeric(df["rank"], errors="coerce").fillna(999).astype(int)
    df["valid"] = df["pb_valid"]
    df["rmsd2"] = df["rmsd"] <= 2.0
    df["kabsch1"] = df["bestfit_rmsd"] <= 1.0
    df["headline"] = df["rmsd2"] & df["valid"]
    df["arm"] = label
    return df


def _per_complex(df: pd.DataFrame, depth: int) -> pd.DataFrame:
    sub = df if depth == 0 else df[df["rank"] <= depth]
    return sub.groupby("protein")[CRITERIA].any()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", action="append", nargs=2, metavar=("LABEL", "CSV"),
                    required=True, help="repeatable: --arm meeko path/to.csv")
    ap.add_argument("--reference", default=None,
                    help="LABEL of the arm every other arm is contrasted against")
    ap.add_argument("--strata-csv", type=Path, default=None,
                    help="receptor_parity.csv; adds a paired contrast restricted "
                         "to complexes whose two arms saw an identical receptor, "
                         "which separates the converter's chemistry from the "
                         "residues one converter drops")
    ap.add_argument("--drop-ids", default="",
                    help="comma-separated complex ids to exclude from every arm; "
                         "use it to remove complexes the arms did not dock under "
                         "the same conditions")
    ap.add_argument("--out-prefix", type=Path, required=True)
    args = ap.parse_args()

    dropped = {i.strip() for i in args.drop_ids.split(",") if i.strip()}

    arms = {label: _load(Path(csv), label) for label, csv in args.arm}
    if dropped:
        arms = {k: v[~v["protein"].astype(str).isin(dropped)] for k, v in arms.items()}
        print(f"excluded {len(dropped)} complexes: {', '.join(sorted(dropped))}\n")
    ref = args.reference or list(arms)[0]

    same_receptor: set[str] | None = None
    if args.strata_csv:
        parity = pd.read_csv(args.strata_csv)
        flag = parity["receptor_identical"].map(
            lambda x: str(x).strip().lower() in ("true", "1"))
        same_receptor = set(parity.loc[flag, "protein"].astype(str))
        print(f"receptor-parity stratum: {len(same_receptor)} of {len(parity)} "
              f"complexes had a byte-identical heavy-atom receptor in both arms\n")

    print("=" * 88)
    print("POSE-LEVEL COVERAGE AND VALIDITY")
    print("=" * 88)
    cov_rows = []
    for label, df in arms.items():
        n_pose, n_cplx = len(df), df["protein"].nunique()
        k = int(df["valid"].sum())
        lo, hi = wilson_ci(k, n_pose)
        print(f"{label:10s}  complexes {n_cplx:4d}   poses {n_pose:6,d}   "
              f"PB-valid {k:6,d}/{n_pose:,d} = {100*k/n_pose:5.1f}%  "
              f"[{100*lo:.1f}, {100*hi:.1f}]")
        cov_rows.append({"arm": label, "complexes": n_cplx, "poses": n_pose,
                         "valid_poses": k, "valid_rate": k / n_pose,
                         "valid_lo": lo, "valid_hi": hi})
    pd.DataFrame(cov_rows).to_csv(f"{args.out_prefix}_pose_level.csv", index=False)

    per_arm_rows, paired_rows = [], []
    for depth in DEPTHS:
        tag = "all" if depth == 0 else f"top-{depth}" if depth > 1 else "rank-1"
        tables = {label: _per_complex(df, depth) for label, df in arms.items()}
        print()
        print("=" * 88)
        print(f"DEPTH {tag}")
        print("=" * 88)
        print(f"{'arm':10s} {'n':>5s}  " + "  ".join(f"{c:>22s}" for c in CRITERIA))
        for label, tbl in tables.items():
            cells = []
            for c in CRITERIA:
                k, n = int(tbl[c].sum()), len(tbl)
                lo, hi = wilson_ci(k, n)
                cells.append(f"{k:3d}/{n:3d} {100*k/n:5.1f}% ")
                per_arm_rows.append({"arm": label, "depth": tag, "criterion": c,
                                     "k": k, "n": n, "rate": k / n,
                                     "lo": lo, "hi": hi})
            print(f"{label:10s} {len(tbl):5d}  " + "  ".join(f"{s:>22s}" for s in cells))

        common_all = sorted(set.intersection(*(set(t.index) for t in tables.values())))
        if len(arms) < 2 or not common_all:
            continue
        strata = [("all", common_all)]
        if same_receptor is not None:
            strata.append(("same-receptor",
                           [c for c in common_all if c in same_receptor]))
            strata.append(("reduced-receptor",
                           [c for c in common_all if c not in same_receptor]))
        for stratum, common in strata:
            if not common:
                continue
            print(f"\n  paired on {len(common)} complexes [{stratum}] "
                  f"(reference = {ref})")
            for label, tbl in tables.items():
                if label == ref:
                    continue
                for c in CRITERIA:
                    a = tables[ref].loc[common, c].astype(int).to_numpy()
                    b = tbl.loc[common, c].astype(int).to_numpy()
                    n10, n01, p = mcnemar_exact(a, b)
                    ci = newcombe_paired_diff_ci(a, b)
                    print(f"    {c:9s} {ref}={100*ci['p_a']:5.1f}%  {label}={100*ci['p_b']:5.1f}%  "
                          f"delta={100*ci['diff']:+5.1f} pp [{100*ci['lo']:+.1f}, {100*ci['hi']:+.1f}]  "
                          f"discordant {n10}/{n01}  p={p:.4f}")
                    paired_rows.append({"depth": tag, "stratum": stratum,
                                        "criterion": c, "reference": ref,
                                        "arm": label, "n_common": len(common),
                                        "rate_reference": ci["p_a"], "rate_arm": ci["p_b"],
                                        "diff_pp": 100 * ci["diff"],
                                        "diff_lo_pp": 100 * ci["lo"],
                                        "diff_hi_pp": 100 * ci["hi"],
                                        "n10": n10, "n01": n01, "p_mcnemar": p})

    pd.DataFrame(per_arm_rows).to_csv(f"{args.out_prefix}_per_arm.csv", index=False)
    if paired_rows:
        pd.DataFrame(paired_rows).to_csv(f"{args.out_prefix}_paired.csv", index=False)
    print(f"\nwrote {args.out_prefix}_pose_level.csv, _per_arm.csv, _paired.csv")


if __name__ == "__main__":
    main()
