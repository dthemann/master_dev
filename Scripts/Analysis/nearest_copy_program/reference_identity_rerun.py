#!/usr/bin/env python
"""Re-run the reference-dependent PoseBusters identity checks with the multi-copy ligand file as reference.

The thesis footnote (Results, primary-endpoint paragraph) reports what supplying a co-crystallised reference
to PoseBusters costs the headline endpoint: the four identity checks (molecular formula, molecular bonds,
tetrahedral chirality, double-bond stereochemistry) compare the InChI of the docked pose with the InChI of the
reference ligand (redock config, inchi_options "w"). This script recomputes those four checks per pose against
three choices of reference: the single deposited instance (<ID>_ligand.sdf == record 0 of <ID>_ligands.sdf),
the multi-copy file loaded the way PoseBusters loads it (load_all=True, one molecule with all copies as
conformers) and the nearest deposited copy of each pose (per-pose column nearest_copy_index). It then counts,
per arm and depth, how many complexes the identity requirement removes from the recovery sets under the
instance and the nearest-copy convention. Read-only on every input; writes only to --out-dir.
"""
from __future__ import annotations
import argparse, json, sys
from multiprocessing import Pool
from pathlib import Path
import pandas as pd
from rdkit import Chem, RDLogger
from posebusters.tools.loading import safe_load_mol
from posebusters.modules.identity import standardize_and_get_inchi, _compare_inchis, stereo_tetrahedral_layers

RDLogger.DisableLog("rdApp.*")
CHECKS = ["formula", "connections", "stereo_tetrahedral", "stereo_dbond"]

def _inchi(mol):
    try:
        return standardize_and_get_inchi(mol, options="w") if mol is not None else None
    except Exception:
        return None

def _compare(inchi_pred, inchi_true):
    if not inchi_pred or not inchi_true:
        return {c: False for c in CHECKS} | {"inchi_valid": False}
    r = _compare_inchis(inchi_pred, inchi_true)
    return {"formula": bool(r.get("formula", False)), "connections": bool(r.get("connections", False)),
            "stereo_tetrahedral": all(bool(r.get(k, False)) for k in stereo_tetrahedral_layers),
            "stereo_dbond": bool(r.get("stereo_dbond", False)), "inchi_valid": True}

def _protein_job(args):
    pid, bench_dir, rows = args
    d = Path(bench_dir) / pid
    inst = _inchi(safe_load_mol(d / f"{pid}_ligand.sdf"))
    comb = _inchi(safe_load_mol(d / f"{pid}_ligands.sdf", load_all=True))
    copies = [m for m in Chem.SDMolSupplier(str(d / f"{pid}_ligands.sdf"), removeHs=False)]
    copy_inchi = [_inchi(m) for m in copies]
    out = []
    for r in rows:
        pred = _inchi(safe_load_mol(Path(r["pose_file"])))
        j = r.get("nearest_copy_index")
        j = int(j) if j is not None and j == j else 0
        near = copy_inchi[j] if 0 <= j < len(copy_inchi) else inst
        rec = {"protein": pid, "method": r["method"], "pose_file": r["pose_file"], "n_copies_file": len(copies),
               "nearest_copy_index": j, "copy_inchi_all_equal": len(set(copy_inchi)) == 1}
        for tag, ref in (("instance", inst), ("combined", comb), ("nearest", near)):
            for k, v in _compare(pred, ref).items():
                rec[f"{tag}_{k}"] = v
            rec[f"{tag}_identity_pass"] = all(rec[f"{tag}_{c}"] for c in CHECKS)
        out.append(rec)
    return out

def eff_rank(df):
    r = pd.to_numeric(df["rank"], errors="coerce").astype(float)
    eq = df["method"].astype(str).str.startswith("equibind")
    if eq.any():
        e = df[eq].copy()
        g = pd.to_numeric(e.get("gnina_affinity"), errors="coerce"); s = pd.to_numeric(e.get("smina_affinity"), errors="coerce")
        e["_k"] = g.where(e["method"].str.endswith("_gnina"), s)
        o = e.sort_values(["method", "protein", "_k", "pose_file"], kind="mergesort")
        r.loc[o.index] = o.groupby(["method", "protein"]).cumcount().values + 1
    return r

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv", type=Path, required=True, help="nearest-convention per-pose table (has rmsd_ref_instance)")
    ap.add_argument("--benchmark-dir", type=Path, default=Path("Data/PoseBuster Benchmark Set"))
    ap.add_argument("--arms", nargs="+", default=["autodock_mgltools_exh128_gnina", "diffdock_smina", "equibind_unguided_gnina"])
    ap.add_argument("--depths", nargs="+", type=int, default=[1, 5, 10, 15, 30])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args()
    if a.out_dir.exists() and any(a.out_dir.iterdir()) and not a.overwrite:
        sys.exit(f"REFUSING: {a.out_dir} exists and is not empty")
    a.out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(a.per_pose_csv, low_memory=False)
    df = df[df["method"].isin(a.arms)].copy()
    df["eff_rank"] = eff_rank(df)
    if "nearest_copy_index" not in df.columns:      # instance-convention table: the nearest block collapses onto the instance
        df["nearest_copy_index"] = 0
    jobs = [(pid, str(a.benchmark_dir), g[["method", "pose_file", "nearest_copy_index"]].to_dict("records"))
            for pid, g in df.groupby("protein")]
    with Pool(a.workers) as pool:
        recs = [r for chunk in pool.imap_unordered(_protein_job, jobs, chunksize=2) for r in chunk]
    ident = pd.DataFrame(recs)
    ident.to_csv(a.out_dir / "identity_checks_per_pose.csv", index=False)
    m = df.merge(ident, on=["protein", "method", "pose_file"], how="left")
    valid = m["pb_valid"].astype(str).str.lower().eq("true")
    near_n = pd.to_numeric(m["rmsd"], errors="coerce") <= 2.0
    near_i = pd.to_numeric(m["rmsd_ref_instance"], errors="coerce") <= 2.0 if "rmsd_ref_instance" in m else near_n
    ids = sorted(df["protein"].unique())
    rows = []
    for arm in a.arms:
        for depth in a.depths:
            sel = (m["method"] == arm) & (m["eff_rank"] <= depth)
            for conv, near, idtag in (("instance", near_i, "instance"), ("nearest", near_n, "nearest"), ("nearest_combined_ref", near_n, "combined")):
                base = m[sel].assign(ok=(valid & near)[sel]).groupby("protein")["ok"].any().reindex(ids).fillna(False)
                gated = m[sel].assign(ok=(valid & near & m[f"{idtag}_identity_pass"].fillna(False).astype(bool))[sel]).groupby("protein")["ok"].any().reindex(ids).fillna(False)
                lost = sorted(base.index[base & ~gated])
                rows.append(dict(arm=arm, depth=depth, convention=conv, identity_reference=idtag, recovered=int(base.sum()),
                                 recovered_with_identity=int(gated.sum()), cost=int(base.sum() - gated.sum()), lost_complexes=";".join(lost)))
    cost = pd.DataFrame(rows); cost.to_csv(a.out_dir / "identity_cost_by_depth.csv", index=False)
    summ = {"n_complexes": len(ids), "arms": a.arms,
            "per_pose_identity_fail_rate": {tag: float((~ident[f"{tag}_identity_pass"]).mean()) for tag in ("instance", "combined", "nearest")},
            "poses_where_nearest_differs_from_instance": int((ident["instance_identity_pass"] != ident["nearest_identity_pass"]).sum()),
            "complexes_with_unequal_copy_inchi": sorted(ident.loc[~ident["copy_inchi_all_equal"], "protein"].unique().tolist()),
            "cost": {f"{r.arm}|{r.convention}|{r.depth}": r.cost for r in cost.itertuples()}}
    (a.out_dir / "identity_cost_summary.json").write_text(json.dumps(summ, indent=2))
    print(cost[cost.depth.isin([1, 15, 30])].to_string(index=False))
    print(json.dumps({k: v for k, v in summ.items() if k != "cost"}, indent=1))

if __name__ == "__main__":
    main()
