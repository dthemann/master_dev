#!/usr/bin/env python
"""Reference-convention sensitivity table (plan v2 D9, Phase 2.7).

Reads a NEAREST-convention hub table (per_pose_metrics.csv written with
``--reference-convention nearest``) and contrasts, for the three headline arms
at depths 1, 5, 15 and 30:

  recovery   PB-valid and rmsd <= 2 A, existence over the top-d poses per complex
             (primary, nearest copy) against PB-valid and rmsd_ref_instance <= 2 A
             (sensitivity, single deposited instance), with the gain
  Form       bestfit_rmsd <= 1 A and rmsd < 1000 against the *_ref_instance twins,
             pooled over the top-d poses and per complex (any qualifying pose)
  Combined   PB-valid and rmsd <= 2 A and bestfit_rmsd <= 1 A (the triple gate),
             pooled and per complex

plus the single-copy stratum (n_copies == 1, identical under both conventions by
construction, asserted) and the multi-copy stratum, the AutoDock against DiffDock
margin under both conventions with exact McNemar and the Newcombe paired
interval, and the multi-copy counts on the 303 (from n_copies) and on the 308
(from the record count of <ID>_ligands.sdf) bases.

No re-scoring from SDFs happens here: every number is read from the hub columns.
EquiBind's ``rank`` is a 999 sentinel, so its effective rank is the ascending
gnina_affinity (smina_affinity for *_smina arms) within (method, protein).

Outputs (inside --out-dir, which must be NEW or empty):
  reference_convention_sensitivity.csv   long table, one row per arm x depth x stratum
  reference_convention_margins.csv       AutoDock against DiffDock per depth and convention
  reference_convention_sensitivity.tex   LaTeX table (house style captions)
  reference_convention_summary.json      headline numbers and provenance
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/home/manndo/master_dev")
ANA = ROOT / "Scripts" / "Analysis"
sys.path.insert(0, str(ANA))
from stats_utils import mcnemar_exact, newcombe_paired_diff_ci  # noqa: E402

DEFAULT_ARMS = ["autodock_mgltools_exh128_gnina", "diffdock_smina", "equibind_unguided_gnina"]
ARM_LABEL = {"autodock_mgltools_exh128_gnina": "AutoDock Vina + gnina",
             "diffdock_smina": "DiffDock + smina",
             "equibind_unguided_gnina": "EquiBind + gnina"}
DEFAULT_DEPTHS = [1, 5, 15, 30]
DEFAULT_BENCHMARK_DIR = "Data/PoseBuster Benchmark Set"
DEFAULT_IDS_308 = "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
REQUIRED = ["method", "protein", "rank", "rmsd", "rmsd_ref_instance", "bestfit_rmsd",
            "bestfit_rmsd_ref_instance", "pb_valid", "n_copies", "reference_convention"]
NEAR_THR, FORM_THR, EXPLODED = 2.0, 1.0, 1000.0


def prepare_out_dir(out_dir: Path) -> None:
    if out_dir.exists():
        if not out_dir.is_dir():
            sys.exit(f"ERROR: --out-dir {out_dir} exists and is not a directory")
        if any(out_dir.iterdir()):
            sys.exit(f"ERROR: --out-dir {out_dir} exists and is not empty; choose a NEW directory")
    else:
        out_dir.mkdir(parents=True)


def effective_rank(df: pd.DataFrame) -> pd.Series:
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


def load_table(path: Path, arms: list[str]) -> pd.DataFrame:
    header = list(pd.read_csv(path, nrows=0).columns)
    missing = [c for c in REQUIRED if c not in header]
    if missing:
        sys.exit(f"ERROR: {path} lacks columns {missing}; a nearest-convention hub table is required")
    usecols = REQUIRED + [c for c in ("gnina_affinity", "smina_affinity", "nearest_copy_is_ref",
                                      "nearest_copy_index", "ref_copy_index") if c in header]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    conv = sorted(set(df["reference_convention"].dropna().astype(str).unique()))
    if conv != ["nearest"]:
        sys.exit(f"ERROR: reference_convention column is {conv}, expected ['nearest']")
    ids = sorted(df["protein"].unique())
    ncop = df.groupby("protein")["n_copies"].agg(["min", "max"])
    if (ncop["min"] != ncop["max"]).any():
        sys.exit("ERROR: n_copies is not constant within a complex")
    df = df[df["method"].isin(arms)].copy()
    absent = [a for a in arms if a not in set(df["method"])]
    if absent:
        sys.exit(f"ERROR: arms not present in the table: {absent}")
    df["eff_rank"] = effective_rank(df)
    df["valid"] = df["pb_valid"].astype(str).str.lower().isin(["true", "1", "1.0"])
    for c in ("rmsd", "rmsd_ref_instance", "bestfit_rmsd", "bestfit_rmsd_ref_instance"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # gates, per pose, both conventions
    df["near_nearest"] = df["valid"] & (df["rmsd"] <= NEAR_THR)
    df["near_instance"] = df["valid"] & (df["rmsd_ref_instance"] <= NEAR_THR)
    df["form_nearest"] = (df["bestfit_rmsd"] <= FORM_THR) & (df["rmsd"] < EXPLODED)
    df["form_instance"] = (df["bestfit_rmsd_ref_instance"] <= FORM_THR) & (df["rmsd_ref_instance"] < EXPLODED)
    df["comb_nearest"] = df["near_nearest"] & (df["bestfit_rmsd"] <= FORM_THR)
    df["comb_instance"] = df["near_instance"] & (df["bestfit_rmsd_ref_instance"] <= FORM_THR)
    df.attrs["ids"] = ids
    df.attrs["n_copies"] = ncop["min"].to_dict()
    return df


def complex_flags(sub: pd.DataFrame, ids: list[str], col: str) -> pd.Series:
    """Existence over the poses of ``sub`` per complex, indexed by every id (False when absent)."""
    got = sub.groupby("protein")[col].any()
    return got.reindex(ids, fill_value=False).astype(bool)


def build_rows(df: pd.DataFrame, arms: list[str], depths: list[int]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ids = df.attrs["ids"]
    ncop = pd.Series(df.attrs["n_copies"]).reindex(ids)
    strata = {"all": pd.Series(True, index=ids),
              "single_copy": (ncop == 1),
              "multi_copy": (ncop > 1)}
    rows, flags = [], {}
    for arm in arms:
        d_arm = df[df["method"] == arm]
        for depth in depths:
            sub = d_arm[d_arm["eff_rank"] <= depth]
            f = {k: complex_flags(sub, ids, k) for k in
                 ("near_nearest", "near_instance", "form_nearest", "form_instance", "comb_nearest", "comb_instance")}
            flags[(arm, depth)] = f
            for sname, smask in strata.items():
                sub_s = sub[sub["protein"].isin(smask[smask].index)]
                n_c = int(smask.sum())
                r_near = int(f["near_nearest"][smask].sum())
                r_inst = int(f["near_instance"][smask].sum())
                if sname == "single_copy" and r_near != r_inst:
                    sys.exit(f"ERROR: single-copy stratum differs between conventions ({arm}, d={depth}); "
                             "the table is not a faithful nearest-convention table")
                lost_form = int((f["form_instance"][smask] & ~f["form_nearest"][smask]).sum())
                gained_form = int((~f["form_instance"][smask] & f["form_nearest"][smask]).sum())
                rows.append(dict(
                    arm=arm, arm_label=ARM_LABEL.get(arm, arm), depth=depth, stratum=sname,
                    n_complexes=n_c, n_poses=int(len(sub_s)),
                    recovered_nearest=r_near, recovered_instance=r_inst, gain=r_near - r_inst,
                    recovered_nearest_pct=round(r_near / n_c * 100, 1) if n_c else np.nan,
                    recovered_instance_pct=round(r_inst / n_c * 100, 1) if n_c else np.nan,
                    lost_under_nearest=int((f["near_instance"][smask] & ~f["near_nearest"][smask]).sum()),
                    form_poses_nearest=int(sub_s["form_nearest"].sum()),
                    form_poses_instance=int(sub_s["form_instance"].sum()),
                    form_complexes_nearest=int(f["form_nearest"][smask].sum()),
                    form_complexes_instance=int(f["form_instance"][smask].sum()),
                    form_complexes_lost=lost_form, form_complexes_gained=gained_form,
                    combined_poses_nearest=int(sub_s["comb_nearest"].sum()),
                    combined_poses_instance=int(sub_s["comb_instance"].sum()),
                    combined_complexes_nearest=int(f["comb_nearest"][smask].sum()),
                    combined_complexes_instance=int(f["comb_instance"][smask].sum()),
                ))
    table = pd.DataFrame(rows)

    # AutoDock against DiffDock margin, both conventions, exact McNemar and Newcombe
    mrows = []
    ad, dd = arms[0], arms[1]
    for depth in depths:
        for conv in ("nearest", "instance"):
            a = flags[(ad, depth)][f"near_{conv}"].to_numpy().astype(int)
            b = flags[(dd, depth)][f"near_{conv}"].to_numpy().astype(int)
            n10, n01, p = mcnemar_exact(a, b)          # n10 = AutoDock only, n01 = DiffDock only
            nc = newcombe_paired_diff_ci(b, a)         # diff = p_autodock - p_diffdock
            mrows.append(dict(depth=depth, convention=conv, n=len(ids),
                              autodock=int(a.sum()), diffdock=int(b.sum()),
                              margin_pp=round((a.sum() - b.sum()) / len(ids) * 100, 2),
                              autodock_only=n10, diffdock_only=n01, mcnemar_p=p,
                              newcombe_lo_pp=round(nc["lo"] * 100, 2), newcombe_hi_pp=round(nc["hi"] * 100, 2)))
    margins = pd.DataFrame(mrows)
    return table, margins, flags


def count_308(ids_file: Path, bench_dir: Path) -> dict | None:
    if not ids_file.exists() or not bench_dir.is_dir():
        return None
    ids = [l.strip() for l in ids_file.read_text().splitlines() if l.strip()]
    multi, n_read = 0, 0
    for pid in ids:
        f = bench_dir / pid / f"{pid}_ligands.sdf"
        if not f.exists():
            continue
        n_read += 1
        n_rec = sum(1 for line in f.open(errors="replace") if line.startswith("$$$$"))
        multi += int(n_rec > 1)
    return {"n_ids": len(ids), "n_read": n_read, "multi_copy": multi, "single_copy": n_read - multi}


def fmt_p(p: float) -> str:
    if p < 1e-4:
        return f"{p:.1e}".replace("e-0", "e-")
    return f"{p:.3f}" if p < 0.995 else "1.000"


def write_tex(table: pd.DataFrame, margins: pd.DataFrame, depths: list[int], n303: int,
              multi303: int, c308: dict | None, path: Path) -> None:
    allrows = table[table.stratum == "all"]
    lines = []
    lines.append("% Generated by Scripts/Analysis/nearest_copy_program/reference_convention_table.py")
    lines.append(f"% {datetime.now().isoformat(timespec='seconds')}; do not edit by hand")
    lines.append("\\begin{table}[tb]")
    lines.append("\\centering\\small")
    mc = f"{multi303} of {n303}"
    c308s = (f" and {c308['multi_copy']} of {c308['n_read']} source entries" if c308 else "")
    lines.append("\\caption{Sensitivity of the near-native endpoint to the reference convention. "
                 "Each cell counts complexes with at least one pose in the top-$d$ that passes the gate. "
                 "Recovery is PoseBusters-valid and within 2\\,\\AA{} of the nearest deposited copy "
                 "(primary) or of the single reference instance (sensitivity). Form is a best-fit RMSD "
                 "of at most 1\\,\\AA{} without a validity gate and Combined adds the form condition to "
                 "the recovery gate. The single-copy stratum is identical under both conventions by "
                 f"construction. Multi-copy complexes number {mc} analysed{c308s}.}}")
    lines.append("\\label{tab:reference_convention_sensitivity}")
    lines.append("\\begin{tabular}{llrrrrrrrr}")
    lines.append("\\toprule")
    lines.append("Tool & $d$ & \\multicolumn{3}{c}{Recovery} & \\multicolumn{2}{c}{Form} & "
                 "\\multicolumn{2}{c}{Combined} & Single copy \\\\")
    lines.append("\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}\\cmidrule(lr){10-10}")
    lines.append(" & & nearest & instance & gain & nearest & instance & nearest & instance & both \\\\")
    lines.append("\\midrule")
    single = table[table.stratum == "single_copy"].set_index(["arm", "depth"])
    for arm in allrows.arm.unique():
        sub = allrows[allrows.arm == arm]
        for i, (_, r) in enumerate(sub.iterrows()):
            label = ARM_LABEL.get(arm, arm) if i == 0 else ""
            s_rec = int(single.loc[(arm, r.depth), "recovered_nearest"])
            s_n = int(single.loc[(arm, r.depth), "n_complexes"])
            lines.append(f"{label} & {int(r.depth)} & {r.recovered_nearest} & {r.recovered_instance} & "
                         f"{r.gain:+d} & {r.form_complexes_nearest} & {r.form_complexes_instance} & "
                         f"{r.combined_complexes_nearest} & {r.combined_complexes_instance} & "
                         f"{s_rec} of {s_n} \\\\")
        lines.append("\\addlinespace")
    lines.append("\\midrule")
    lines.append("\\multicolumn{10}{l}{AutoDock Vina + gnina against DiffDock + smina, recovery gate, "
                 "exact McNemar and Newcombe 95\\,\\% interval in percentage points} \\\\")
    for _, m in margins.iterrows():
        lines.append(f"\\multicolumn{{10}}{{l}}{{$d={int(m.depth)}$, {m.convention}. "
                     f"{m.autodock} against {m.diffdock} of {m.n}, margin {m.margin_pp:+.1f}\\,pp "
                     f"({m.newcombe_lo_pp:+.1f} to {m.newcombe_hi_pp:+.1f}), discordant {m.autodock_only} "
                     f"against {m.diffdock_only}, $p$ {fmt_p(m.mcnemar_p)}}} \\\\")
    lines.append("\\bottomrule")
    lines.append("\\end{tabular}")
    lines.append("\\end{table}")
    path.write_text("\n".join(lines) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-pose-csv", required=True, help="nearest-convention hub table (read only)")
    ap.add_argument("--out-dir", required=True, help="NEW (or empty) directory for the outputs")
    ap.add_argument("--arms", nargs="+", default=DEFAULT_ARMS,
                    help="three arms; the first two form the AutoDock against DiffDock margin")
    ap.add_argument("--depths", nargs="+", type=int, default=DEFAULT_DEPTHS)
    ap.add_argument("--benchmark-dir", default=DEFAULT_BENCHMARK_DIR, help="read only, for the 308 count")
    ap.add_argument("--ids-308", default=DEFAULT_IDS_308, help="the 308 source ids, one per line")
    args = ap.parse_args(argv)
    if len(args.arms) < 2:
        sys.exit("ERROR: at least two arms are needed")

    os.chdir(ROOT)
    out_dir = Path(args.out_dir)
    prepare_out_dir(out_dir)

    df = load_table(Path(args.per_pose_csv), args.arms)
    ids = df.attrs["ids"]
    ncop = pd.Series(df.attrs["n_copies"]).reindex(ids)
    table, margins, _ = build_rows(df, args.arms, args.depths)
    c308 = count_308(Path(args.ids_308), Path(args.benchmark_dir))
    n303, multi303 = len(ids), int((ncop > 1).sum())

    csv_out = out_dir / "reference_convention_sensitivity.csv"
    m_out = out_dir / "reference_convention_margins.csv"
    tex_out = out_dir / "reference_convention_sensitivity.tex"
    json_out = out_dir / "reference_convention_summary.json"
    table.to_csv(csv_out, index=False)
    margins.to_csv(m_out, index=False)
    write_tex(table, margins, args.depths, n303, multi303, c308, tex_out)

    head = {}
    for depth in args.depths:
        sub = table[(table.stratum == "all") & (table.depth == depth)].set_index("arm")
        head[f"depth_{depth}"] = {
            "nearest": [int(sub.loc[a, "recovered_nearest"]) for a in args.arms],
            "instance": [int(sub.loc[a, "recovered_instance"]) for a in args.arms],
            "form_complexes_nearest": [int(sub.loc[a, "form_complexes_nearest"]) for a in args.arms],
            "form_complexes_instance": [int(sub.loc[a, "form_complexes_instance"]) for a in args.arms],
            "combined_complexes_nearest": [int(sub.loc[a, "combined_complexes_nearest"]) for a in args.arms],
            "combined_complexes_instance": [int(sub.loc[a, "combined_complexes_instance"]) for a in args.arms],
        }
    single = table[(table.stratum == "single_copy") & (table.depth == 1)].set_index("arm")
    summary = {
        "generated": datetime.now().isoformat(timespec="seconds"),
        "script": str(Path(__file__).resolve()),
        "per_pose_csv": str(Path(args.per_pose_csv).resolve()),
        "arms": args.arms, "depths": args.depths,
        "n_complexes": n303, "multi_copy_303": multi303, "single_copy_303": n303 - multi303,
        "count_308": c308,
        "headline": head,
        "single_copy_rank1": {a: [int(single.loc[a, "recovered_nearest"]), int(single.loc[a, "n_complexes"])]
                              for a in args.arms},
        "margins": margins.to_dict("records"),
        "gates": {"recovery": "pb_valid and rmsd <= 2", "form": "bestfit_rmsd <= 1 and rmsd < 1000",
                  "combined": "pb_valid and rmsd <= 2 and bestfit_rmsd <= 1",
                  "instance_twins": "rmsd_ref_instance, bestfit_rmsd_ref_instance"},
    }
    json_out.write_text(json.dumps(summary, indent=2, default=float))

    pd.set_option("display.width", 220)
    print(f"complexes {n303}, multi-copy {multi303}, single-copy {n303 - multi303}; 308 basis: {c308}")
    print(table[table.stratum == "all"][["arm", "depth", "recovered_nearest", "recovered_instance", "gain",
                                          "form_complexes_nearest", "form_complexes_instance",
                                          "combined_complexes_nearest", "combined_complexes_instance"]]
          .to_string(index=False))
    print(table[table.stratum != "all"][["arm", "depth", "stratum", "n_complexes", "recovered_nearest",
                                          "recovered_instance", "gain"]].to_string(index=False))
    print(margins.to_string(index=False))
    print(f"wrote {csv_out}\nwrote {m_out}\nwrote {tex_out}\nwrote {json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
