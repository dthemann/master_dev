#!/usr/bin/env python
"""Metal-adjacent stratum generator for the PoseBusters benchmark (plan D4).

Rule (decision D4 of PLAN_nearest_copy_primary_endpoint_2026-09-07_v2.md): a
complex is metal-adjacent when the minimum heavy-atom distance from the
reference ligand instance (record 0 of ``<ID>_ligands.sdf``, identical to
``<ID>_ligand.sdf``) to ANY atom of a residue that contains at least one metal
element is <= 5 A, measured on ``Data/PoseBuster Benchmark Set/<ID>/<ID>_protein.pdb``
(the PoseBusters-shipped, cofactor-retaining file, not the docked receptor).
The metal element set is the alkali, alkaline-earth, transition and
post-transition metals. On the 303 analysed ids the rule gives 73 / 78 / 81
complexes at 4 / 5 / 6 A, which is what body_appendix_short.tex:165 prints.
The atom-level rule (distance to the metal atom itself) gives 65 / 77 / 80 and
does not reproduce the print.

Membership is a complex-level property of the reference instance and does not
follow the scored copy. Two complexes disagree across copies (7JHQ_VAJ is in by
the reference only, 7TB0_UD1 is in by copy 1 only); the summary lists them.

No generator for this stratum existed before 2026-09-08; the 78-list behind
the appendix was never registered.

Outputs (all inside --out-dir, which must be a new or empty directory):
  metal_stratum.csv           one row per deposited ligand copy
  metal_stratum_summary.json  counts per cutoff, cofactor/ion split,
                              copy-dependence list, provenance
  metal_stratum_contrasts.csv per arm, depth 1/15/30, cutoff 4/5/6 A and
                              stratum: recovery under PB-valid & rmsd <= 2,
                              between-tool differences (exact McNemar,
                              Fisher on the discordant pairs) and the
                              no-qualifying-pose shares at depth 30; the same
                              block for the nearest-copy convention when
                              --nearest-csv is given

Example:
  python Scripts/Analysis/metal_stratum.py --out-dir posebusters_results/metal_stratum
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stats_utils import mcnemar_exact  # noqa: E402

# ---------------------------------------------------------------------------
# defaults
# ---------------------------------------------------------------------------
DEFAULT_BENCHMARK_DIR = "Data/PoseBuster Benchmark Set"
DEFAULT_REPORT_DIR = "posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report"
DEFAULT_IDS_FILE = f"{DEFAULT_REPORT_DIR}/analysed_cohort_ids.txt"
DEFAULT_PER_POSE_CSV = f"{DEFAULT_REPORT_DIR}/per_pose_metrics.csv"
DEFAULT_ARMS = ["autodock_mgltools_exh128_gnina", "diffdock_smina", "equibind_unguided_gnina"]
DEFAULT_CUTOFFS = [4.0, 5.0, 6.0]
PRIMARY_CUTOFF = 5.0
DEPTHS = (1, 15, 30)
MAX_DEPTH = 30
RMSD_GATE_A = 2.0
OUTPUT_FILES = ("metal_stratum.csv", "metal_stratum_summary.json", "metal_stratum_contrasts.csv")

# alkali, alkaline-earth, transition (incl. lanthanides) and post-transition metals
METALS = set("""
Li Na K Rb Cs Fr
Be Mg Ca Sr Ba Ra
Sc Ti V Cr Mn Fe Co Ni Cu Zn
Y Zr Nb Mo Tc Ru Rh Pd Ag Cd
La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu
Hf Ta W Re Os Ir Pt Au Hg
Al Ga In Sn Tl Pb Bi Po
""".split())


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------
def read_pdb_atoms(path: Path) -> pd.DataFrame:
    """ATOM/HETATM records of a PDB file with the element column resolved."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = line[76:78].strip().capitalize()
            if not el and line.startswith("HETATM"):
                # Name-derived fallback for HETATM only: an ATOM record named CA
                # is an alpha carbon, not calcium (the known PDBQT-CA trap).
                name = line[12:16].strip()
                el = "".join(c for c in name if c.isalpha())[:2].capitalize()
            rows.append(dict(
                rec=line[:6].strip(), name=line[12:16].strip(),
                resn=line[17:20].strip(), chain=line[21], resi=line[22:27].strip(),
                el=el, x=float(line[30:38]), y=float(line[38:46]), z=float(line[46:54])))
    df = pd.DataFrame(rows)
    df["reskey"] = df.rec + "|" + df.chain + "|" + df.resn + "|" + df.resi
    # The D4 rule counts HETATM metals (free ions and metal-bearing cofactors);
    # no shipped <ID>_protein.pdb has a metal in an ATOM record, and the guard
    # keeps that true if one ever does.
    df["is_metal"] = df.el.isin(METALS) & (df.rec == "HETATM")
    n_atom_metal = int((df.el.isin(METALS) & (df.rec == "ATOM")).sum())
    if n_atom_metal:
        print(f"  WARNING: {path.name}: {n_atom_metal} ATOM-record atom(s) carry a metal "
              "element and are NOT counted (the rule counts HETATM metals only).")
    return df


def heavy_xyz(mol) -> np.ndarray:
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())
                     if mol.GetAtomWithIdx(i).GetAtomicNum() > 1], dtype=float)


def _min_dist(prot_xyz: np.ndarray, lig_xyz: np.ndarray) -> tuple[float, int]:
    """(minimum distance, index of the nearest protein atom); inf when empty."""
    if prot_xyz.shape[0] == 0 or lig_xyz.shape[0] == 0:
        return float("inf"), -1
    from scipy.spatial import cKDTree
    d, idx = cKDTree(prot_xyz).query(lig_xyz)
    j = int(np.argmin(d))
    return float(d[j]), int(idx[j])


def stratum_geometry(pid: str, bench_dir: Path, cutoffs: list[float]) -> list[dict]:
    """One dict per deposited ligand copy of <pid>."""
    from rdkit import Chem
    d = bench_dir / pid
    prot = read_pdb_atoms(d / f"{pid}_protein.pdb")
    metal_reskeys = set(prot.loc[prot.is_metal, "reskey"])
    res_sel = prot.reskey.isin(metal_reskeys)
    res_atoms = prot[res_sel].reset_index(drop=True)
    metal_atoms = prot[prot.is_metal].reset_index(drop=True)
    res_natoms = prot.groupby("reskey").size()
    res_metal_els = prot[prot.is_metal].groupby("reskey").el.apply(lambda s: "/".join(sorted(set(s))))

    ref_mol = Chem.SDMolSupplier(str(d / f"{pid}_ligand.sdf"), removeHs=False, sanitize=False)[0]
    if ref_mol is None:
        raise RuntimeError(f"{pid}: cannot read {pid}_ligand.sdf")
    ref_xyz = heavy_xyz(ref_mol)
    copies = [m for m in Chem.SDMolSupplier(str(d / f"{pid}_ligands.sdf"), removeHs=False, sanitize=False)
              if m is not None]
    if not copies:
        raise RuntimeError(f"{pid}: no records in {pid}_ligands.sdf")

    out = []
    for k, m in enumerate(copies):
        xyz = heavy_xyz(m)
        is_ref = (xyz.shape == ref_xyz.shape
                  and np.allclose(np.sort(xyz, axis=0), np.sort(ref_xyz, axis=0), atol=1e-3))
        d_atom, j_atom = _min_dist(metal_atoms[["x", "y", "z"]].to_numpy(), xyz)
        d_res, j_res = _min_dist(res_atoms[["x", "y", "z"]].to_numpy(), xyz)
        row = dict(protein=pid, copy_index=k, n_copies=len(copies), is_ref=bool(is_ref),
                   n_heavy_atoms=int(xyz.shape[0]),
                   d_metal_atom=d_atom, d_metal_res=d_res,
                   nearest_metal_element="", nearest_metal_resname="",
                   nearest_metal_res_elements="", nearest_metal_res_natoms=0,
                   cofactor_or_ion="none",
                   n_metal_atoms_in_protein=int(prot.is_metal.sum()),
                   n_metal_residues_in_protein=len(metal_reskeys))
        if j_atom >= 0:
            row["nearest_metal_element"] = str(metal_atoms.iloc[j_atom].el)
        if j_res >= 0:
            rk = res_atoms.iloc[j_res].reskey
            row["nearest_metal_resname"] = str(res_atoms.iloc[j_res].resn)
            row["nearest_metal_res_elements"] = str(res_metal_els.get(rk, ""))
            row["nearest_metal_res_natoms"] = int(res_natoms[rk])
            row["cofactor_or_ion"] = "cofactor" if res_natoms[rk] > 1 else "ion"
        for c in cutoffs:
            row[f"member_{_ctag(c)}"] = bool(d_res <= c)
        out.append(row)
    return out


def _ctag(c: float) -> str:
    return f"{int(c)}A" if float(c).is_integer() else f"{c:g}A".replace(".", "p")


# ---------------------------------------------------------------------------
# recovery contrasts
# ---------------------------------------------------------------------------
def _to_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(("true", "1", "1.0", "yes"))


def effective_rank(df: pd.DataFrame) -> pd.Series:
    """Rank within (method, protein). EquiBind carries a rank==999 sentinel and
    is ordered by its refinement energy (gnina_affinity for *_gnina arms,
    smina_affinity otherwise), lowest first, NaN last, pose_file as the
    deterministic tie-breaker (a stable sort on the zero-padded unguided_NNN
    index, the same order as the hub's _gnina_affinity_rank). Every other arm
    keeps its own rank column. If the table already carries eff_rank it is used
    verbatim."""
    given = None
    if "eff_rank" in df.columns:
        # Never trust a supplied eff_rank silently: derive the rank from the
        # affinities (the hub's own ordering) and report every disagreement.
        given = pd.to_numeric(df["eff_rank"], errors="coerce")
    eff = pd.to_numeric(df["rank"], errors="coerce").copy()
    eq = df.method.str.startswith("equibind")
    if eq.any():
        e = df[eq].copy()
        key_g = pd.to_numeric(e.get("gnina_affinity"), errors="coerce") if "gnina_affinity" in e else np.nan
        key_s = pd.to_numeric(e.get("smina_affinity"), errors="coerce") if "smina_affinity" in e else np.nan
        e["_k"] = pd.Series(key_g, index=e.index).where(e.method.str.endswith("_gnina"),
                                                       pd.Series(key_s, index=e.index))
        e = e.sort_values(["method", "protein", "_k", "pose_file"], na_position="last", kind="mergesort")
        e["_eff"] = e.groupby(["method", "protein"]).cumcount() + 1
        eff.loc[e.index] = e["_eff"]
    if given is not None:
        n_dis = int((given.fillna(-1).astype(int) != eff.fillna(-1).astype(int)).sum())
        print(f"  eff_rank column present: {n_dis} of {len(df)} rows disagree with the "
              "affinity-derived rank; the derived rank is used.")
    return eff.astype(int)


def recovery_matrix(poses: pd.DataFrame, ids: list[str], arm: str, rmsd_col: str,
                    depth: int) -> pd.Series:
    """Boolean per id: any pose with eff_rank <= depth passing PB-valid & rmsd <= 2."""
    sub = poses[(poses.method == arm) & (poses.eff_rank <= depth)]
    ok = sub[_to_bool(sub.pb_valid) & (pd.to_numeric(sub[rmsd_col], errors="coerce") <= RMSD_GATE_A)]
    hit = set(ok.protein)
    return pd.Series([pid in hit for pid in ids], index=ids, dtype=bool)


def fisher_discordant(n01_free, n10_free, n01_adj, n10_adj) -> float:
    from scipy.stats import fisher_exact
    table = [[n01_free, n10_free], [n01_adj, n10_adj]]
    if sum(table[0]) == 0 or sum(table[1]) == 0:
        return float("nan")
    return float(fisher_exact(table)[1])


def contrasts_block(poses: pd.DataFrame, ids: list[str], member: dict[float, pd.Series],
                    arms: list[str], cutoffs: list[float], rmsd_col: str, convention: str) -> list[dict]:
    rows = []
    rec = {(a, d): recovery_matrix(poses, ids, a, rmsd_col, d) for a in arms for d in DEPTHS}
    for c in cutoffs:
        mem = member[c].reindex(ids).astype(bool)
        strata = {"metal_free": ~mem, "metal_adjacent": mem}
        for d in DEPTHS:
            for sname, mask in strata.items():
                n = int(mask.sum())
                for a in arms:
                    k = int(rec[(a, d)][mask].sum())
                    rows.append(dict(convention=convention, cutoff_A=c, stratum=sname, n=n, depth=d,
                                     kind="recovery", arm=a, arm_b="", k=k,
                                     rate_pct=round(100 * k / n, 6) if n else np.nan))
            # between-tool differences b - a for every pair in arm order
            for i, a in enumerate(arms):
                for b in arms[i + 1:]:
                    disc = {}
                    for sname, mask in strata.items():
                        n = int(mask.sum())
                        va, vb = rec[(a, d)][mask].astype(int), rec[(b, d)][mask].astype(int)
                        n10, n01, p = mcnemar_exact(va.to_numpy(), vb.to_numpy())  # n10 = a only
                        ka, kb = int(va.sum()), int(vb.sum())
                        disc[sname] = (n01, n10)
                        rows.append(dict(convention=convention, cutoff_A=c, stratum=sname, n=n, depth=d,
                                         kind="contrast", arm=a, arm_b=b, k_a=ka, k_b=kb,
                                         rate_a_pct=round(100 * ka / n, 6) if n else np.nan,
                                         rate_b_pct=round(100 * kb / n, 6) if n else np.nan,
                                         diff_pp=round(100 * (kb - ka) / n, 6) if n else np.nan,
                                         n_a_only=n10, n_b_only=n01, p_mcnemar=p))
                    pf = fisher_discordant(*disc["metal_free"], *disc["metal_adjacent"])
                    rows.append(dict(convention=convention, cutoff_A=c, stratum="between_strata",
                                     n=len(ids), depth=d, kind="fisher_discordant", arm=a, arm_b=b,
                                     n_b_only_free=disc["metal_free"][0], n_a_only_free=disc["metal_free"][1],
                                     n_b_only_adj=disc["metal_adjacent"][0], n_a_only_adj=disc["metal_adjacent"][1],
                                     p_fisher=pf))
        # no qualifying pose at depth MAX_DEPTH among rank-1 failures
        for sname, mask in strata.items():
            for a in arms:
                fail1 = mask & ~rec[(a, 1)]
                none30 = fail1 & ~rec[(a, MAX_DEPTH)]
                nf, nn = int(fail1.sum()), int(none30.sum())
                rows.append(dict(convention=convention, cutoff_A=c, stratum=sname, n=nf, depth=MAX_DEPTH,
                                 kind="no_qualifying_pose_among_rank1_failures", arm=a, arm_b="", k=nn,
                                 rate_pct=round(100 * nn / nf, 6) if nf else np.nan))
    return rows


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def load_ids(ids_file: Path | None, per_pose: pd.DataFrame) -> tuple[list[str], str]:
    if ids_file is not None and ids_file.is_file():
        ids = [l.strip() for l in ids_file.read_text().splitlines() if l.strip()]
        if ids:
            return sorted(set(ids)), str(ids_file)
    ids = sorted(per_pose.protein.dropna().unique())
    return ids, "derived from per-pose table protein column"


def prepare_out_dir(out_dir: Path, overwrite: bool) -> None:
    if out_dir.exists():
        if not out_dir.is_dir():
            sys.exit(f"ERROR: --out-dir {out_dir} exists and is not a directory")
        present = sorted(p.name for p in out_dir.iterdir())
        if present:
            foreign = [p for p in present if p not in OUTPUT_FILES]
            if foreign or not overwrite:
                sys.exit(f"ERROR: --out-dir {out_dir} exists and is not empty ({len(present)} entries). "
                         "Choose a NEW directory" + ("" if foreign else ", or pass --overwrite to replace "
                                                    "this script's own previous outputs") + ".")
    out_dir.mkdir(parents=True, exist_ok=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-dir", default=DEFAULT_BENCHMARK_DIR)
    ap.add_argument("--ids-file", default=DEFAULT_IDS_FILE,
                    help="one id per line; when absent the ids come from the per-pose table")
    ap.add_argument("--per-pose-csv", default=DEFAULT_PER_POSE_CSV, help="read only")
    ap.add_argument("--out-dir", required=True, help="NEW directory; refused when it exists and is non-empty")
    ap.add_argument("--arms", nargs="+", default=DEFAULT_ARMS)
    ap.add_argument("--cutoffs", nargs="+", type=float, default=DEFAULT_CUTOFFS)
    ap.add_argument("--primary-cutoff", type=float, default=PRIMARY_CUTOFF)
    ap.add_argument("--nearest-csv", default=None,
                    help="per-pose CSV carrying rmsd_nearest_copy (and optionally eff_rank) for the "
                         "nearest-copy convention block")
    ap.add_argument("--overwrite", action="store_true",
                    help="allow an --out-dir that holds only this script's own previous outputs")
    args = ap.parse_args(argv)

    bench_dir = Path(args.benchmark_dir)
    out_dir = Path(args.out_dir)
    cutoffs = sorted(set(float(c) for c in args.cutoffs) | {float(args.primary_cutoff)})
    prepare_out_dir(out_dir, args.overwrite)

    usecols = ["method", "protein", "pose_file", "rank", "rmsd", "pb_valid", "gnina_affinity", "smina_affinity"]
    header = pd.read_csv(args.per_pose_csv, nrows=0).columns
    per_pose = pd.read_csv(args.per_pose_csv, usecols=[c for c in usecols if c in header], low_memory=False)
    ids, ids_source = load_ids(Path(args.ids_file) if args.ids_file else None, per_pose)
    print(f"ids: {len(ids)} ({ids_source})")
    missing_arms = [a for a in args.arms if a not in set(per_pose.method)]
    if missing_arms:
        sys.exit(f"ERROR: arms not in per-pose table: {missing_arms}")
    poses = per_pose[per_pose.method.isin(args.arms) & per_pose.protein.isin(ids)].copy()
    poses["eff_rank"] = effective_rank(poses)

    # ---- geometry --------------------------------------------------------
    rows = []
    for i, pid in enumerate(ids, 1):
        rows.extend(stratum_geometry(pid, bench_dir, cutoffs))
        if i % 50 == 0 or i == len(ids):
            print(f"  geometry {i}/{len(ids)}")
    geo = pd.DataFrame(rows)
    ref = geo[geo.is_ref].drop_duplicates("protein").set_index("protein")
    no_ref = sorted(set(ids) - set(ref.index))
    if no_ref:
        sys.exit(f"ERROR: reference instance not found among _ligands.sdf records for {no_ref}")
    ref = ref.reindex(ids)
    rec0_is_ref = int((geo[geo.is_ref].groupby("protein").copy_index.min() == 0).sum())

    member = {c: ref[f"member_{_ctag(c)}"].astype(bool) for c in cutoffs}
    counts = {_ctag(c): int(member[c].sum()) for c in cutoffs}
    print("metal-adjacent counts:", counts)

    # cofactor / ion split per cutoff (classification of the nearest metal-bearing residue)
    split = {}
    for c in cutoffs:
        m = ref[member[c]]
        cof = m[m.cofactor_or_ion == "cofactor"]
        ion = m[m.cofactor_or_ion == "ion"]
        split[_ctag(c)] = dict(
            n_cofactor=int(len(cof)), n_ion=int(len(ion)),
            cofactor_resnames=dict(sorted(Counter(cof.nearest_metal_resname).items())),
            cofactor_complexes={r: sorted(cof.index[cof.nearest_metal_resname == r]) for r in sorted(set(cof.nearest_metal_resname))},
            ion_elements=dict(sorted(Counter(ion.nearest_metal_element).items())),
            nearest_metal_atom_elements=dict(sorted(Counter(m.nearest_metal_element).items())))

    # copy-dependence: complexes whose membership differs across deposited copies
    copy_dep = {}
    for c in cutoffs:
        col = f"member_{_ctag(c)}"
        lst = []
        for pid, g in geo.groupby("protein"):
            if g[col].nunique() > 1:
                lst.append(dict(protein=pid, n_copies=int(g.n_copies.iloc[0]),
                                reference_member=bool(g.loc[g.is_ref, col].iloc[0]),
                                per_copy=[dict(copy_index=int(r["copy_index"]), is_ref=bool(r["is_ref"]),
                                               d_metal_res=round(float(r["d_metal_res"]), 3),
                                               d_metal_atom=round(float(r["d_metal_atom"]), 3),
                                               nearest_metal_element=r["nearest_metal_element"],
                                               nearest_metal_resname=r["nearest_metal_resname"],
                                               member=bool(r[col])) for _, r in g.iterrows()]))
        copy_dep[_ctag(c)] = lst
    print("copy-dependent complexes at primary cutoff:",
          [x["protein"] for x in copy_dep[_ctag(args.primary_cutoff)]])

    # ---- contrasts -------------------------------------------------------
    contrasts = contrasts_block(poses, ids, member, args.arms, cutoffs, "rmsd", "reference")
    nearest_info = None
    if args.nearest_csv:
        near = pd.read_csv(args.nearest_csv, low_memory=False)
        need = {"method", "protein", "pb_valid", "rmsd_nearest_copy"}
        if not need <= set(near.columns):
            sys.exit(f"ERROR: --nearest-csv lacks columns {sorted(need - set(near.columns))}")
        near = near[near.method.isin(args.arms) & near.protein.isin(ids)].copy()
        missing = [a for a in args.arms if a not in set(near.method)]
        if missing:
            print(f"WARNING: --nearest-csv has no rows for {missing}; those arms count as unrecovered")
        if "eff_rank" not in near.columns and near.method.str.startswith("equibind").any():
            # EquiBind needs its refinement energy for ranking; take it from the
            # per-pose table when the nearest CSV does not carry it.
            aff_cols = [c for c in ("gnina_affinity", "smina_affinity") if c in per_pose.columns]
            if not aff_cols:
                sys.exit("ERROR: EquiBind rows need eff_rank or gnina/smina affinity for ranking")
            near = near.drop(columns=[c for c in aff_cols if c in near.columns])
            near = near.merge(per_pose[["method", "pose_file"] + aff_cols].drop_duplicates(["method", "pose_file"]),
                              on=["method", "pose_file"], how="left")
            print("  nearest CSV lacks eff_rank: EquiBind ranked from the per-pose table's affinities")
        near["eff_rank"] = effective_rank(near)
        contrasts += contrasts_block(near, ids, member, args.arms, cutoffs, "rmsd_nearest_copy", "nearest")
        nearest_info = dict(path=str(args.nearest_csv), n_rows=int(len(near)),
                            eff_rank_source="column" if "eff_rank" in pd.read_csv(args.nearest_csv, nrows=0).columns
                            else "derived")
    contrasts_df = pd.DataFrame(contrasts)

    # ---- write -----------------------------------------------------------
    geo_cols = ["protein", "copy_index", "n_copies", "is_ref", "n_heavy_atoms", "d_metal_atom", "d_metal_res",
                "nearest_metal_element", "nearest_metal_resname", "nearest_metal_res_elements",
                "nearest_metal_res_natoms", "cofactor_or_ion"] + [f"member_{_ctag(c)}" for c in cutoffs] + \
               ["n_metal_atoms_in_protein", "n_metal_residues_in_protein"]
    geo[geo_cols].to_csv(out_dir / "metal_stratum.csv", index=False, float_format="%.6f")
    contrasts_df.to_csv(out_dir / "metal_stratum_contrasts.csv", index=False)

    prim = _ctag(args.primary_cutoff)
    summary = dict(
        rule=("metal-adjacent := min heavy-atom distance from the reference ligand instance (record 0 of "
              "<ID>_ligands.sdf == <ID>_ligand.sdf) to ANY atom of a residue containing >= 1 metal element "
              f"<= {args.primary_cutoff:g} A on <ID>_protein.pdb; membership is a complex-level property of the "
              "reference instance and does not follow the scored copy"),
        metal_elements=sorted(METALS),
        inputs=dict(benchmark_dir=str(bench_dir), ids_source=ids_source, per_pose_csv=str(args.per_pose_csv),
                    nearest_csv=nearest_info, arms=args.arms, cutoffs_A=cutoffs, primary_cutoff_A=args.primary_cutoff,
                    gate=f"pb_valid & rmsd <= {RMSD_GATE_A:g} A, existence over eff_rank <= depth", depths=list(DEPTHS)),
        n_ids=len(ids),
        n_copies_total=int(len(geo)),
        n_multi_copy_ids=int((ref.n_copies > 1).sum()),
        reference_is_record0_in=rec0_is_ref,
        n_ids_with_metal_in_protein=int((ref.n_metal_atoms_in_protein > 0).sum()),
        counts_metal_adjacent=counts,
        counts_metal_free={k: len(ids) - v for k, v in counts.items()},
        primary=dict(cutoff=prim, n_metal_adjacent=counts[prim], n_metal_free=len(ids) - counts[prim],
                     metal_adjacent_ids=sorted(ref.index[member[args.primary_cutoff]])),
        cofactor_ion_split=split,
        copy_dependent_membership=copy_dep,
        atom_level_counts_for_reference={_ctag(c): int((ref.d_metal_atom <= c).sum()) for c in cutoffs},
        outputs=list(OUTPUT_FILES),
    )
    (out_dir / "metal_stratum_summary.json").write_text(json.dumps(summary, indent=1))

    # console digest of the primary-cutoff contrasts
    c5 = contrasts_df[(contrasts_df.cutoff_A == args.primary_cutoff)]
    for conv in sorted(c5.convention.unique()):
        print(f"\n[{conv}] recovery at {prim}")
        r = c5[(c5.convention == conv) & (c5.kind == "recovery")]
        print(r.pivot_table(index=["arm", "depth"], columns="stratum", values="k").to_string())
        ct = c5[(c5.convention == conv) & (c5.kind == "contrast")]
        print(ct[["depth", "stratum", "arm", "arm_b", "k_a", "k_b", "diff_pp", "n_a_only", "n_b_only", "p_mcnemar"]]
              .to_string(index=False))
    print(f"\nwrote {out_dir}/{{{', '.join(OUTPUT_FILES)}}}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
