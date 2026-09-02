#!/usr/bin/env python3
"""Re-bust raw AutoDock poses under ONE hydrogen convention.

WHY
---
PoseBusters' ``internal_energy`` check scores two arms of the same experiment on
different molecules. ``posebusters/tools/molecules.py`` ``add_hydrogens_with_uff_positions``
pins every atom already present in the file, ``AddHs(addCoords=True)``, then UFF-minimises
only the ADDED atoms; ``optimize_positions`` short-circuits and returns the UNRELAXED
energy when nothing was added. So:

  * Meeko-ligand arms convert via ``mk_export`` and ship FULL explicit hydrogens
    (``num_h_added == 0`` on 9,043/9,043 raw poses) -> scored on the unrelaxed UFF energy,
    including every reconstructed hydrogen.
  * MGLTools-ligand arms convert via authoritative-template reconstruction and ship
    HEAVY ATOMS ONLY (``num_h_added > 0`` on 8,991/8,991) -> PoseBusters adds and
    UFF-relaxes the entire hydrogen shell before taking the energy.

That single asymmetry is worth ~1.50 pp of the arms' 1.56 pp validity difference, so the
published PB-valid column is not comparable between preparations.

Vina uses a UNITED-ATOM model: nonpolar hydrogens are merged and never predicted. Every
hydrogen in a Meeko-arm pose SDF is therefore a reconstruction artefact, not a docking
output. Penalising their positions measures the converter, not the docking. The defensible
common convention is thus the HEAVY-ATOM one: strip hydrogens and let PoseBusters add and
relax them identically for every arm.

WHAT THIS DOES
--------------
Re-busts the EXISTING converted pose SDFs -- the very files that produced the published
numbers -- with hydrogens stripped. Nothing is re-docked and nothing is re-converted, so
hydrogen presence is the only variable that moves. Results go to a NEW output directory;
no published artefact is touched.

Built-in control: the MGLTools arms are already heavy-atom-only, so stripping is a no-op
for them and their re-busted numbers MUST reproduce the published ones. ``--mode asis``
re-busts unmodified files and must reproduce the published numbers for every arm; run it
on a sample first to validate the harness before trusting the ``strip`` numbers.

USAGE
-----
    # 1. validate the harness reproduces published numbers on unmodified input
    python rebust_hydrogen_normalised.py --mode asis --limit 400

    # 2. the real run
    python rebust_hydrogen_normalised.py --mode strip

    # 3. cascade table in posebusters_pose_comparison.py's own format
    python rebust_hydrogen_normalised.py --report
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

REPO = Path("/home/manndo/master_dev")
SRC_CSV = (REPO / "posebusters_results/benchmark_full_protein_vina_scoring/dock"
           / "posebusters_filtered_results.csv")
PER_POSE = (REPO / "posebusters_results/benchmark_full_protein_vina_scoring/dock"
            / "pose_comparison_report/per_pose_metrics.csv")
OUT_DIR = REPO / "posebusters_results/benchmark_full_protein_hnorm"

# Every raw AutoDock arm. autodock/autodock_vinardo ship hydrogens; the mgltools
# arms do not and therefore act as the no-op control. All four exhaustiveness
# settings belong here — listing only 32 (the unsuffixed key) and 64 silently
# dropped exh18/exh92 from the harness and from its --report cascade.
ARMS = ("autodock", "autodock_vinardo",
        "autodock_mgltools", "autodock_mgltools_exh18",
        "autodock_mgltools_exh64", "autodock_mgltools_exh92")

_CTX: dict = {}


# ---------------------------------------------------------------- worker


def _init(mode: str, staging: str, timeout: int) -> None:
    global _CTX
    from posebusters import PoseBusters
    from rdkit import RDLogger

    RDLogger.DisableLog("rdApp.*")
    try:
        import run_posebusters as R
        R._limit_worker_threads()
    except Exception:
        pass
    # NOTE: the buster is built ONCE per worker and bust() is called DIRECTLY, not in a
    # forked grandchild. PoseBusters @caches the 100-conformer ensemble energy by InChI
    # (`modules/energy_ratio.py: get_energies`), and that cache is the single dominant
    # cost of the run. Forking per pose throws it away every time, so the ensemble gets
    # regenerated 36,092 times instead of once per distinct ligand -- measured ~6x
    # slower. Tasks are therefore grouped by COMPLEX (below) so every pose of a ligand,
    # across all four arms, hits a warm cache.
    _CTX = {"pb": PoseBusters(config="dock"), "mode": mode,
            "staging": Path(staging), "timeout": timeout}


def _one(rec: dict) -> dict:
    """Bust a single pose. Returns the PoseBusters row plus identity columns."""
    from rdkit import Chem

    pose = Path(rec["pose_file"])
    out = {k: rec[k] for k in ("docking_method", "protein", "ligand",
                               "pose_name", "pose_file")}
    try:
        pose_path = pose
        out["n_h_in_file"] = -1
        if _CTX["mode"] == "strip":
            mol = Chem.MolFromMolFile(str(pose), removeHs=False, sanitize=True)
            if mol is None:
                out["_error"] = "load failed"
                return out
            n_h = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 1)
            out["n_h_in_file"] = n_h
            # Already heavy-atom-only (the MGLTools arms): stripping is the identity, so
            # bust the ORIGINAL file. That keeps those arms byte-identical to the
            # published run and makes them an exact reproduction control rather than an
            # SDF round-trip that could perturb coordinate precision.
            if n_h:
                mol = Chem.RemoveHs(mol)
                target = _CTX["staging"] / f"{pose.stem}__noH.sdf"
                if not target.exists():
                    w = Chem.SDWriter(str(target))
                    w.write(mol)
                    w.close()
                pose_path = target

        df = _CTX["pb"].bust(str(pose_path), None, rec["receptor_file_used"],
                             full_report=True)
        row = df.iloc[0].to_dict()
        for k, v in row.items():
            if k not in out:
                out[k] = v
        out["_error"] = ""
    except Exception as exc:
        out["_error"] = f"{type(exc).__name__}: {exc}"
    return out


def _one_complex(recs: list) -> list:
    """Bust every pose of one complex in a single worker, warm-caching its ensemble."""
    return [_one(r) for r in recs]


# ---------------------------------------------------------------- bust phase


def run_bust(mode: str, arms, workers: int, limit: int, timeout: int,
              group_by: str = "complex") -> Path:
    from multiprocessing import Pool

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    staging = OUT_DIR / "stripped_sdf"
    staging.mkdir(exist_ok=True)
    out_csv = OUT_DIR / f"rebust_{mode}.csv"

    df = pd.read_csv(SRC_CSV, low_memory=False)
    df = df[(df["optimizer"] == "original") & (df["docking_method"].isin(arms))].copy()
    need = ["pose_file", "receptor_file_used", "docking_method",
            "protein", "ligand", "pose_name"]
    df = df[need]
    if limit:
        # stratified head per arm so a calibration run covers every arm
        df = df.groupby("docking_method", group_keys=False).head(
            max(1, limit // max(1, len(arms))))

    done: set[str] = set()
    if out_csv.exists():
        prev = pd.read_csv(out_csv, low_memory=False)
        done = set(prev["pose_file"].astype(str))
        print(f"resuming: {len(done):,} poses already done", flush=True)
    recs = [r for r in df.to_dict("records") if str(r["pose_file"]) not in done]

    # Group by COMPLEX, not by pose. Every arm docks the same ligand, so all ~120 poses
    # of one (protein, ligand) share a single InChI and therefore a single cached
    # conformer-ensemble energy. Handing a whole complex to one worker turns 36,092
    # ensemble generations into ~303.
    if group_by == "pose":
        # Tail-latency escape hatch. Complex grouping serialises ~120 poses onto one
        # core, so a complex with a large receptor (7SUC_COM: 3.0 MB receptor, ~20 s
        # per bust) holds the whole run open for ~40 min on its own. Their ensembles
        # are cheap (~0.2 s), so for those the cache is worth far less than the
        # parallelism -- spread the poses across all workers instead.
        tasks = [[r] for r in recs]
    else:
        groups: dict = {}
        for r in recs:
            groups.setdefault((r["protein"], r["ligand"]), []).append(r)
        tasks = list(groups.values())
    n_poses = sum(len(t) for t in tasks)
    print(f"mode={mode}  arms={list(arms)}  poses to bust: {n_poses:,} "
          f"in {len(tasks):,} complex groups (workers={workers})", flush=True)
    if not tasks:
        return out_csv

    t0 = time.time()
    rows: list[dict] = []
    first_write = not out_csv.exists()
    seen = 0
    with Pool(workers, initializer=_init,
              initargs=(mode, str(staging), timeout), maxtasksperchild=8) as pool:
        for gi, res in enumerate(pool.imap_unordered(_one_complex, tasks, chunksize=1), 1):
            rows.extend(res)
            seen += len(res)
            if gi % 5 == 0 or gi == len(tasks):
                pd.DataFrame(rows).to_csv(out_csv, mode="a", index=False,
                                          header=first_write)
                first_write = False
                rows = []
                el = time.time() - t0
                rate = seen / max(el, 1e-9)
                print(f"  {seen:,}/{n_poses:,} poses  ({gi}/{len(tasks)} complexes)  "
                      f"{rate:.1f} pose/s  eta {(n_poses-seen)/max(rate,1e-9)/60:.1f} min",
                      flush=True)
    print(f"wrote {out_csv}", flush=True)
    return out_csv


# ---------------------------------------------------------------- report phase


def _validity(df: pd.DataFrame) -> pd.Series:
    """All-pass verdict over the canonical stock PoseBusters dock columns."""
    import run_posebusters as R

    cols = [c for c in R.CANONICAL_TEST_COLUMNS if c in df.columns]
    if not cols:
        raise SystemExit("no PoseBusters test columns in the re-bust output")
    ok = pd.Series(True, index=df.index)
    for c in cols:
        v = df[c].map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0})
        v = pd.to_numeric(v, errors="coerce")
        v = v.fillna(pd.to_numeric(df[c], errors="coerce")).fillna(1.0)
        ok &= v.astype(bool)
    return ok


def build_report(mode: str) -> None:
    reb = pd.read_csv(OUT_DIR / f"rebust_{mode}.csv", low_memory=False)
    reb = reb[reb["_error"].fillna("") == ""].copy()
    reb["pb_valid_new"] = _validity(reb)
    key = ["docking_method", "protein", "ligand", "pose_name"]
    reb["_k"] = reb[key].astype(str).agg("\x00".join, axis=1)

    pp = pd.read_csv(PER_POSE, low_memory=False)
    pp = pp[pp["method"].isin(reb["docking_method"].unique())].copy()
    pp["_k"] = (pp[["method", "protein", "ligand", "pose_name"]]
                .astype(str).agg("\x00".join, axis=1))

    m = reb.set_index("_k")["pb_valid_new"]
    hit = pp["_k"].map(m)
    print(f"joined {hit.notna().sum():,}/{len(pp):,} per-pose rows "
          f"({hit.isna().sum():,} unmatched)")
    pp["pb_valid_published"] = pp["pb_valid"]
    pp["pb_valid"] = hit.fillna(pp["pb_valid"])

    sys.path.insert(0, str(REPO / "Scripts/Analysis"))
    import posebusters_pose_comparison as P

    specs = P._all_variant_yield_specs(pp)
    casc = P.aggregate_pose_validity_cascade(pp, specs)
    table = P._render_pose_validity_cascade_table(casc)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    casc.to_csv(OUT_DIR / f"pose_validity_cascade_{mode}.csv", index=False)
    hdr = (f"Pose-level validity & accuracy cascade — hydrogen-normalised re-bust "
           f"(mode={mode})\n"
           f"PB-validity recomputed with hydrogens stripped so PoseBusters adds and\n"
           f"UFF-relaxes them identically for every arm. RMSD/Kabsch columns are the\n"
           f"published values (unaffected: the RMSD path strips H on both sides).\n")
    (OUT_DIR / f"pose_validity_cascade_{mode}.txt").write_text(hdr + "\n" + table + "\n")
    print("\n" + hdr)
    print(table)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", choices=("strip", "asis"), default="strip")
    ap.add_argument("--arms", nargs="*", default=list(ARMS))
    ap.add_argument("--workers", "-j", type=int, default=14)
    ap.add_argument("--limit", type=int, default=0,
                    help="calibration: cap poses (split across arms)")
    ap.add_argument("--timeout", type=int, default=600)
    ap.add_argument("--group", choices=("complex", "pose"), default="complex",
                    help="complex = warm ensemble cache (default); pose = max parallelism, "
                         "use for the slow-receptor tail")
    ap.add_argument("--report", action="store_true",
                    help="skip busting, just build the cascade from an existing CSV")
    args = ap.parse_args()

    if args.report:
        build_report(args.mode)
        return
    run_bust(args.mode, args.arms, args.workers, args.limit, args.timeout,
             args.group)


if __name__ == "__main__":
    main()
