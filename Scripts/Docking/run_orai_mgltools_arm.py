#!/usr/bin/env python
"""Dock the Orai x Benchmark ADFRsuite-receptor + ADFRsuite-ligand arm.

Why this exists instead of the notebook cell or run_autodock.py's CLI:

  * Notebook cell 87 (id ac223686) re-runs run_workflow, which has no existence guard and
    re-prepares receptors into the SHARED Dockings/Orai_Benchmark/_staging. Orai1WT-START-Fr0
    hits an unguarded OpenMM minimisation that is not bit-reproducible, so a re-prep yields a
    THIRD distinct receptor and retroactively breaks the published Vinardo arm's provenance.
  * run_autodock.py main() calls run_workflow with convert_ligands_with_meeko=True regardless of
    prep_tool, and writes to output_base/<ligdir>_<converter>, not the layout used here.

This driver therefore calls run_autodock_vina() and optimize_autodock_results() ONLY, over
inputs that were staged once by hand. It never invokes a converter.

Usage:
    python Scripts/Docking/run_orai_mgltools_arm.py -c <config.yaml> [--dry-run] [--limit N]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from Scripts.Docking.run_autodock import (  # noqa: E402
    run_autodock_vina,
    optimize_autodock_results,
    generate_summary,
    get_cpu_model,
    preflight_mgl_pose_templates,
)

# The four receptors the published control docked against (Jun-22 prep).
# Orai1WT-START-Fr0 is NOT reproducible; it must be copied, never regenerated.
EXPECTED_FR0_SHA256 = "6b3ab996af727c806ca481e0f1f5ca7b0fe923df8dddf48d002e930e468206e9"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fail(msg: str) -> None:
    raise SystemExit(f"REFUSING TO RUN: {msg}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-c", "--config", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true",
                    help="run every guard and the template preflight, then exit before docking")
    ap.add_argument("--limit", type=int, default=0,
                    help="dock only the first N ligands (validation runs); 0 = all")
    ap.add_argument("--ligand-prep", choices=("adfrsuite", "meeko"), default="adfrsuite",
                    help="which ligand converter this arm's staged PDBQTs must carry. "
                         "'meeko' is for the paired twin arms whose only difference from the "
                         "ADFRsuite arm is the ligand converter; it flips Guard 1's polarity so "
                         "a wrong-converter staging is still caught, rather than unguarded.")
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    root = Path.cwd()

    out_dir = (root / cfg["output_dir"]).resolve()
    log_dir = (root / cfg["log_dir"]).resolve()
    rec_pdbqt_dir = (root / cfg["receptors_dir"] / "pdbqt").resolve()
    lig_pdbqt_dir = (root / cfg["ligand_dirs"][0] / "pdbqt").resolve()

    lig_label = "ADFRsuite" if args.ligand_prep == "adfrsuite" else "Meeko"
    print("=" * 78)
    print(f"{out_dir.parent.name} | ADFRsuite receptor + {lig_label} ligand | "
          f"Vina exh {cfg['exhaustiveness']} | num_modes {cfg['num_modes']} | "
          f"energy_range {cfg['energy_range']} | optimisation={cfg.get('optimization')}")
    print("=" * 78)
    print(f"  receptors : {rec_pdbqt_dir}")
    print(f"  ligands   : {lig_pdbqt_dir}")
    print(f"  output    : {out_dir}")
    print(f"  logs      : {log_dir}")

    # ── Guard 1: ligands must carry the converter this arm declares ──────────
    # Both polarities are guarded. Staging the wrong converter is the single most
    # dangerous mistake available here: the two arms differ ONLY in the ligand
    # converter, so a mis-staged directory silently turns the contrast into a
    # comparison of an arm against itself.
    ligs = sorted(lig_pdbqt_dir.glob("*.pdbqt"))
    if not ligs:
        fail(f"no ligand PDBQTs in {lig_pdbqt_dir}")
    meeko = [p for p in ligs if "REMARK SMILES" in p.read_text(errors="ignore")[:4000]]
    if args.ligand_prep == "adfrsuite":
        if meeko:
            fail(f"{len(meeko)} ligand PDBQTs carry REMARK SMILES (Meeko). "
                 "This arm requires ADFRsuite ligands; mk_export would silently take over.")
        print(f"\n  [ok] ligands  : {len(ligs)} ADFRsuite PDBQTs, 0 Meeko")
    else:
        if len(meeko) != len(ligs):
            fail(f"only {len(meeko)}/{len(ligs)} ligand PDBQTs carry REMARK SMILES. "
                 "This arm requires MEEKO ligands (its whole point is to be the "
                 "ligand-converter twin of the ADFRsuite arm).")
        print(f"\n  [ok] ligands  : {len(ligs)} Meeko PDBQTs (REMARK SMILES on all)")

    # ── Guard 2: receptors must be the Jun-22 set ────────────────────────────
    recs = sorted(rec_pdbqt_dir.glob("*.pdbqt"))
    boxes = sorted(rec_pdbqt_dir.glob("*.box.txt"))
    if len(recs) != 4 or len(boxes) != 4:
        fail(f"expected 4 receptor PDBQTs + 4 box files, found {len(recs)} + {len(boxes)}")
    fr0 = rec_pdbqt_dir / "Orai1WT-START-Fr0.pdbqt"
    got = sha256(fr0)
    if got != EXPECTED_FR0_SHA256:
        fail(f"Orai1WT-START-Fr0.pdbqt sha256 is {got}, expected {EXPECTED_FR0_SHA256} "
             "(the Jun-22 receptor the control docked against). A different hash means an "
             "Aug-06 or freshly re-minimised receptor, which confounds 25% of the arm.")
    print(f"  [ok] receptors: 4 PDBQT + 4 box; START-Fr0 sha256 verified (Jun-22 set)")

    # ── Guard 3: output/log dirs must not already hold a run ─────────────────
    existing = list((out_dir / "docking").glob("*_vina_out.pdbqt")) if out_dir.exists() else []
    if existing:
        fail(f"{out_dir/'docking'} already holds {len(existing)} poses. Docking would mark every "
             "combo 'skipped' and the optimiser would then rescore THOSE poses. Use a new dir.")
    if log_dir.exists() and any(log_dir.glob("*.log")):
        fail(f"{log_dir} already holds Vina logs; they would be truncated. Use a new log dir.")
    print(f"  [ok] output/log dirs are clean")

    # ── Guard 4: nesting must satisfy run_posebusters' parents[2] hop ────────
    # pose at <ROOT>/<X>/docking/<p>.pdbqt -> parents[2] == <ROOT>, which must hold _staging/.
    probe = out_dir / "docking" / "probe.pdbqt"
    resolved_root = probe.parents[2]
    prepared_probe = resolved_root / "_staging" / "ligands" / "pdbqt" / ligs[0].name
    if not prepared_probe.is_file():
        fail("output tree is not nested correctly: run_posebusters resolves the prepared ligand "
             f"at {prepared_probe}, which does not exist. Poses must sit at "
             "<ROOT>/<X>/docking/ with staging at <ROOT>/_staging/ligands/pdbqt/.")
    print(f"  [ok] nesting   : parents[2] -> {resolved_root.name}/_staging/ligands/pdbqt resolves")

    # ── Guard 5: template preflight is ALL-OR-NOTHING; do it BEFORE docking ──
    # Only meaningful for ADFRsuite ligands. Meeko poses carry REMARK SMILES, so
    # _convert_pose_to_sdf takes the mk_export branch and never consults a template —
    # which is also what the published Meeko arms did, so it is the parity path.
    names = [p.stem for p in ligs]
    if cfg.get("optimize_require_template_reconstruction"):
        templates = preflight_mgl_pose_templates(names, cfg)
        missing = [n for n in names if templates.get(n) is None]
        if missing:
            fail(f"{len(missing)} ligand templates unresolved, e.g. {missing[:5]}. With "
                 "optimize_require_template_reconstruction: true this aborts the optimiser AFTER "
                 "the multi-hour docking stage.")
        print(f"  [ok] templates : {len(templates)}/{len(names)} resolved")
    else:
        if args.ligand_prep == "adfrsuite":
            fail("ADFRsuite ligands carry no REMARK SMILES, so mk_export cannot rebuild them "
                 "and template reconstruction is mandatory — set "
                 "optimize_require_template_reconstruction: true.")
        print("  [ok] templates : not required (Meeko ligands rebuild via mk_export)")

    if args.limit:
        ligs = ligs[: args.limit]
        print(f"\n  !! --limit {args.limit}: docking only {len(ligs)} ligand(s)")

    # ── Hand-built manifest: no converter is ever invoked ────────────────────
    manifest = {
        "proteins": [
            *({"name": p.name, "path": str(p), "filetype": "pdbqt",
               "source": f"{p.stem} [pdbqt, staged Jun-22]"} for p in recs),
            *({"name": b.name, "path": str(b), "filetype": "txt",
               "source": f"{b.name} [box, staged]"} for b in boxes),
        ],
        "ligands": [
            {"name": p.name, "path": str(p), "filetype": "pdbqt",
             "source": f"{p.stem} [pdbqt, {lig_label}]"} for p in ligs
        ],
    }
    n_combo = len(recs) * len(ligs)
    print(f"\n  manifest: {len(recs)} receptors x {len(ligs)} ligands = {n_combo} combinations")

    if args.dry_run:
        print("\n[DRY RUN] all guards passed; exiting before docking.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # ── Dock ─────────────────────────────────────────────────────────────────
    _, results = run_autodock_vina(
        base_dir=out_dir,
        log_dir=log_dir,
        prepared_manifest=manifest,
        cfg=cfg,
        cpu_model=get_cpu_model(),
        protein_workflow_data=None,   # never re-prepare
        ligand_workflow_data=None,
    )
    n_ok = sum(r.status == "success" for r in results)
    n_skip = sum(r.status == "skipped" for r in results)
    n_fail = sum(r.status == "failed" for r in results)
    print(f"\nDocking: {n_ok} success | {n_fail} failed | {n_skip} skipped")
    if n_skip:
        fail(f"{n_skip} combinations were SKIPPED — poses already existed. The optimiser would "
             "rescore pre-existing poses. Investigate before continuing.")

    summary = generate_summary(results, cfg)
    (out_dir / "docking_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # ── gnina rescore, inline ────────────────────────────────────────────────
    print("\nRunning inline gnina rescore (template reconstruction)...")
    opt = optimize_autodock_results(results, cfg, out_dir)
    print(f"gnina: {opt}")

    print(f"\nDone. Results: {out_dir}")
    print("VERIFY: every optimization_log.csv row must read converter=rdkit_template_map. "
          "A single mk_export row means Meeko poses were optimised by mistake.")


if __name__ == "__main__":
    main()
