#!/usr/bin/env python
"""AutoDock Vina whole-protein docking for ONE exhaustiveness arm of the
exhaustiveness 32-vs-64 comparison.

This is a restricted, resumable clone of the authoritative whole-protein
benchmark driver (Master_Docking_AD_Full_Protein.ipynb, cell "Benchmark Set ||
Autodock (Vina Scoring)"). It differs from that cell in exactly one respect: it
docks only the complexes named in an explicit ids file instead of all 308
authoritative benchmark complexes. Every preparation, staging, docking and
gnina-rescore step is the same call into Scripts/Docking/run_autodock.py, so the
two arms are produced by identical code and differ only in the config's
``exhaustiveness`` value.

Driven by Exhaustiveness_32_vs_64.ipynb.

Usage
-----
    python Scripts/Docking/run_autodock_exhaustiveness_arm.py \
        -c Scripts/Docking/autodock_vina_docking_config_full_protein_search_vina_exh64.yaml \
        --ids-file Dockings/vina_results_full_protein_vina_scoring_exh64/exhaustiveness_sample_ids.txt

Re-running is safe: a complex whose pose file already validates is not re-docked,
and gnina optimisation fills only what is missing.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import yaml  # noqa: E402

from Scripts.Docking.run_autodock import (  # noqa: E402
    DockingResult,
    OptimizationError,
    build_prepared_manifest,
    generate_summary,
    get_cpu_model,
    get_pdbqt_dir,
    optimize_autodock_results,
    preflight_optimizers,
    resolve_optimizers,
    run_autodock_vina,
    validate_vina_output,
)
from Scripts.Utilities.inject_hetatms import inject_dropped_atom_residues  # noqa: E402
from Scripts.Utilities.prep_docking import run_workflow, _write_box_file  # noqa: E402


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True, type=Path,
                    help="AutoDock YAML config for this arm (defines exhaustiveness + output_dir)")
    ap.add_argument("--ids-file", required=True, type=Path,
                    help="Newline-separated <PDBID>_<LIG> ids to dock (the sampled cohort)")
    ap.add_argument("--benchmark-dir", type=Path, default=Path("Data/PoseBuster Benchmark Set"),
                    help="Directory holding one <PDBID>_<LIG>/ folder per benchmark complex")
    ap.add_argument("--limit", type=int, default=0,
                    help="Dock only the first N ids (0 = all). For pilot runs.")
    ap.add_argument("--expect-exhaustiveness", type=int, default=None,
                    help="Refuse to launch unless the config's exhaustiveness equals this")
    ap.add_argument("--prepared-inputs-from", type=Path, default=None,
                    help="Control-arm result tree whose ALREADY-PREPARED receptor and ligand "
                         "PDBQTs (and search box) this arm should dock verbatim, skipping "
                         "preparation entirely. Required for a controlled exhaustiveness "
                         "comparison: re-deriving the ligand with a newer Meeko can change its "
                         "atom ordering, torsion-tree root and even its rotatable-bond count, "
                         "varying the search space alongside exhaustiveness.")
    return ap.parse_args(argv)


def _pdbqt_atom_multiset(path: Path):
    """Coordinates, partial charge and AutoDock type per atom, order-insensitive.

    This is the ligand's physics. Two PDBQTs with the same multiset describe the
    same molecule in the same conformation even if Meeko wrote the atoms in a
    different order or rooted the torsion tree elsewhere.
    """
    atoms = []
    for ln in path.read_text().splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            atoms.append((ln[30:38].strip(), ln[38:46].strip(), ln[46:54].strip(),
                          ln[70:76].strip(), ln[77:79].strip()))
    return sorted(atoms)


def _torsdof(path: Path) -> int | None:
    for ln in path.read_text().splitlines():
        if ln.startswith("TORSDOF"):
            return int(ln.split()[1])
    return None


def main(argv=None) -> int:
    args = parse_args(argv)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    # ── Guard the scientific invariants of this arm ──────────────────────────
    # The comparison is whole-protein (blind) docking; a ligand box would make the
    # arms incomparable for reasons unrelated to exhaustiveness.
    if str(cfg.get("box_mode", "")).strip().lower() != "whole_protein":
        raise RuntimeError("This comparison is whole-protein (blind) docking only")

    # Post-dock optimisation is optional. `none` produces the raw Vina poses only,
    # which is what a raw-pose exhaustiveness comparison analyses; `gnina` keeps the
    # campaign's Variant-A rescore, where the empirical function drives the
    # minimisation and the CNN only scores the result; `gnina_refinement` is
    # Variant B, where CNN gradients drive the minimisation so the pose itself
    # differs. Anything else is a config slip. Each mode pins its own CNN scoring
    # value, and a mismatch there is the specific slip this guard exists to catch,
    # because run_autodock.py would otherwise be free to run a mode the arm name
    # does not describe. Both optimised modes rank on cnn_affinity so the two are
    # comparable on the same head, isolating the geometry difference.
    mode = str(cfg.get("optimization", "none")).strip().lower()
    _REQUIRED_CNN = {"gnina": "rescore", "gnina_refinement": "refinement"}
    if mode not in ("none", *_REQUIRED_CNN):
        raise RuntimeError(
            "optimization must be 'none', 'gnina' or 'gnina_refinement' for this "
            f"arm (got {mode!r})")
    if mode in _REQUIRED_CNN:
        want = _REQUIRED_CNN[mode]
        if str(cfg.get("gnina_cnn_scoring", "")).strip().lower() != want:
            raise RuntimeError(
                f"optimization: {mode} requires gnina_cnn_scoring: {want}")
        if str(cfg.get("optimize_rank_by", "")).strip().lower() != "cnn_affinity":
            raise RuntimeError(f"optimization: {mode} requires optimize_rank_by: cnn_affinity")

    exhaustiveness = int(cfg["exhaustiveness"])
    if args.expect_exhaustiveness is not None and exhaustiveness != args.expect_exhaustiveness:
        raise RuntimeError(
            f"Config exhaustiveness is {exhaustiveness}, expected "
            f"{args.expect_exhaustiveness} — refusing to write into the wrong arm")

    output_base = Path(cfg["output_dir"])
    log_dir = Path(cfg["log_dir"])
    output_base.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    overwrite = bool(cfg.get("overwrite_existing") or cfg.get("overwrite_poses"))
    overwrite_opt = bool(cfg.get("overwrite_optimization", False))

    opt_tools = resolve_optimizers(cfg)
    if opt_tools:
        preflight_optimizers(cfg)  # fail before a long run, not during it
    else:
        print("Refinement: OFF — committing raw Vina poses only")

    prep_tool = cfg.get("prep_tool", "mgltools")
    converters: list[tuple[str, str]] = []
    if prep_tool in ("mgltools", "both"):
        converters.append(("mgl_tools", "mgltools"))
    if prep_tool in ("meeko", "both"):
        converters.append(("meeko", "meeko"))

    ids = [ln.strip() for ln in args.ids_file.read_text().splitlines()
           if ln.strip() and not ln.lstrip().startswith("#")]
    if not ids:
        raise RuntimeError(f"No ids found in {args.ids_file}")
    if args.limit:
        ids = ids[: args.limit]

    missing = [pid for pid in ids if not (args.benchmark_dir / pid).is_dir()]
    if missing:
        raise RuntimeError(f"Refusing launch: ids absent from disk: {missing[:10]}")

    cpu_model = get_cpu_model()
    print(f"Arm: exhaustiveness={exhaustiveness}  →  {output_base}")
    print(f"Complexes: {len(ids)}   converters: {[c[0] for c in converters]}")
    print(f"CPU: {cpu_model}\n")

    status_path = output_base / "exhaustiveness_arm_status.json"
    issues: list[dict] = []

    def write_status(state, **extra):
        payload = {
            "schema_version": 1,
            "arm": f"exhaustiveness_{exhaustiveness}",
            "exhaustiveness": exhaustiveness,
            "config": str(args.config),
            "ids_file": str(args.ids_file),
            "n_ids": len(ids),
            "status": state,
            "updated_at": datetime.now().isoformat(),
            **extra,
        }
        tmp = status_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2, default=str) + "\n")
        tmp.replace(status_path)

    write_status("running", processed=0, last_complex=None)

    def validated_outputs(vina_out: Path):
        """The single strict committed pose file expected for one complex."""
        cands = sorted((vina_out / "docking").glob("*_vina_out.pdbqt"))
        if len(cands) != 1:
            return []
        try:
            validate_vina_output(cands[0])
        except (OSError, ValueError):
            return []
        return cands

    def optimize_with_reporting(results, vina_out: Path):
        if not opt_tools:
            return None
        try:
            summary = optimize_autodock_results(results, cfg, vina_out, overwrite=overwrite_opt)
        except OptimizationError as exc:
            summary = exc.summary
        if summary.status in {"partial", "failed"}:
            issues.append({"output_dir": str(vina_out), **summary.to_dict()})
            detail = summary.errors[0] if summary.errors else summary.status
            print(f"    [optimizer warning] {summary.failed} pose(s) failed: {detail}")
        return summary

    def optimize_existing(vina_out: Path):
        """gnina rescore over poses that are already committed for this complex."""
        if not opt_tools:
            return
        scoring = cfg.get("scoring_function", "vina")
        dock_dir = vina_out / "docking"
        rec_dir = vina_out.parent / "_staging" / "receptors" / "pdbqt"
        existing = []
        for pose in validated_outputs(vina_out):
            if "__" not in pose.stem:
                continue
            rec_stem, _, rest = pose.stem.partition("__")
            for suf in (f"_{scoring}_vina_out", "_vina_out"):
                if rest.endswith(suf):
                    rest = rest[: -len(suf)]
                    break
            rec_pdbqt = rec_dir / f"{rec_stem}.pdbqt"
            if not rec_pdbqt.exists():
                issues.append({"output_dir": str(vina_out), "status": "failed",
                               "errors": [f"receptor PDBQT missing: {rec_pdbqt}"]})
                continue
            n_models = len(validate_vina_output(pose))
            existing.append(DockingResult(
                protein_name=rec_stem, ligand_name=rest,
                protein_path=rec_pdbqt.resolve(), ligand_path=pose.resolve(),
                output_dir=dock_dir, status="success", scoring_function=scoring,
                num_poses=n_models, pose_files=[pose.resolve()]))
        if existing:
            optimize_with_reporting(existing, vina_out)

    all_results: list[DockingResult] = []
    skipped = failed_prep = 0
    t_start = time.time()

    for idx, pdb_id in enumerate(ids, 1):
        write_status("running", processed=idx - 1, last_complex=pdb_id)
        cdir = args.benchmark_dir / pdb_id
        protein_pdb = cdir / f"{pdb_id}_protein.pdb"
        ligand_start = cdir / f"{pdb_id}_ligand_start_conf.sdf"
        if not protein_pdb.exists() or not ligand_start.exists():
            raise RuntimeError(f"{pdb_id}: missing inputs — cohort must be complete")

        if not overwrite and all(validated_outputs(output_base / pdb_id / cn)
                                 for cn, _ in converters):
            print(f"[{idx}/{len(ids)}] ✓ {pdb_id} — already docked"
                  + (", running gnina optimisation" if opt_tools else ""))
            for cn, _ in converters:
                optimize_existing(output_base / pdb_id / cn)
            skipped += 1
            continue

        print(f"\n{'─' * 60}\n[{idx}/{len(ids)}] {pdb_id}  (exhaustiveness={exhaustiveness})\n{'─' * 60}")

        staging = output_base / pdb_id / "_staging"
        rec_staging, lig_staging = staging / "receptors", staging / "ligands"
        rec_staging.mkdir(parents=True, exist_ok=True)
        lig_staging.mkdir(parents=True, exist_ok=True)
        rec_link, lig_link = rec_staging / protein_pdb.name, lig_staging / ligand_start.name
        if not rec_link.exists():
            rec_link.symlink_to(protein_pdb.resolve())
        if not lig_link.exists():
            lig_link.symlink_to(ligand_start.resolve())

        for conv_name, conv_arg in converters:
            vina_out = output_base / pdb_id / conv_name
            if not overwrite and validated_outputs(vina_out):
                print(f"  ✓ [{conv_name}] already docked"
                      + (" — gnina optimisation only" if opt_tools else " — skipping"))
                optimize_existing(vina_out)
                continue

            # whole_protein mode keeps the box written during receptor prep; the
            # ligand box is never applied. _write_box_file stays imported so this
            # file mirrors the authoritative cell exactly.
            assert cfg["box_mode"] == "whole_protein", "guarded above"

            if args.prepared_inputs_from:
                # ── Dock the control arm's already-prepared inputs ────────────
                # Re-deriving the inputs would re-run PDBFixer, MGLTools and Meeko
                # with today's versions. MGLTools is deterministic so the receptor
                # reproduces byte for byte, but Meeko's atom ordering, torsion-tree
                # root and even the rotatable-bond count can change between
                # versions — an uncontrolled difference sitting right next to
                # exhaustiveness. Copying the control arm's committed files makes
                # the two arms dock byte-identical inputs by construction, and
                # leaves this tree self-contained for PoseBusters, which resolves
                # each pose's receptor from the arm's own _staging directory.
                protein_outputs = ligand_outputs = None
                src = args.prepared_inputs_from / pdb_id / "_staging"
                protein_pdbqt_dir = get_pdbqt_dir(rec_staging)
                ligand_pdbqt_dir = get_pdbqt_dir(lig_staging)
                copied = []
                for sub, dst in (("receptors", protein_pdbqt_dir), ("ligands", ligand_pdbqt_dir)):
                    src_dir = src / sub / "pdbqt"
                    if not src_dir.is_dir():
                        raise RuntimeError(f"{pdb_id}: no prepared {sub} at {src_dir}")
                    for f in sorted(src_dir.iterdir()):
                        if f.is_file():
                            target = dst / f.name
                            if not target.exists() or target.read_bytes() != f.read_bytes():
                                shutil.copyfile(f, target)
                                copied.append(f.name)
                if copied:
                    print(f"    ↳ copied {len(copied)} prepared input file(s) from the control arm")

                rec_pdbqt = protein_pdbqt_dir / f"{pdb_id}_protein_{conv_name}.pdbqt"
                box_file = protein_pdbqt_dir / f"{pdb_id}_protein_{conv_name}.box.txt"
                lig_pdbqt = ligand_pdbqt_dir / f"{pdb_id}_ligand_start_conf.pdbqt"
                for f in (rec_pdbqt, box_file, lig_pdbqt):
                    if not f.exists():
                        raise RuntimeError(f"{pdb_id}: expected prepared input missing: {f}")
                manifest = {
                    "proteins": [
                        {"name": rec_pdbqt.name, "path": str(rec_pdbqt.resolve()),
                         "filetype": "pdbqt", "source": f"{pdb_id} [pdbqt, control arm]"},
                        {"name": box_file.name, "path": str(box_file.resolve()),
                         "filetype": "txt", "source": f"{pdb_id} [box, control arm]"},
                    ],
                    "ligands": [
                        {"name": lig_pdbqt.name, "path": str(lig_pdbqt.resolve()),
                         "filetype": "pdbqt", "source": f"{pdb_id} [pdbqt, control arm]"},
                    ],
                }
            else:
                protein_pdbqt_dir = get_pdbqt_dir(rec_staging)
                protein_outputs = run_workflow(
                    input_dir=rec_staging, contains="proteins", output_dir=protein_pdbqt_dir,
                    skip_pdb_validation=cfg.get("skip_pdb_validation", False),
                    custom_postfix=f"_{conv_name}",
                    process_postfixes=cfg.get("process_postfixes", False),
                    repair_terminals=cfg.get("repair_terminals", False),
                    converter=conv_arg, convert_proteins=True, verbose=False)

                # Re-add the modified residues / ATOM-cofactors MGLTools drops, so
                # Vina docks the same complete protein PoseBusters validates against.
                for rec in (sorted(protein_pdbqt_dir.glob("*.pdbqt"))
                            + sorted(protein_pdbqt_dir.glob("*.pdb"))):
                    inj = inject_dropped_atom_residues(protein_pdb, rec)
                    if inj["atoms_injected"]:
                        print(f"    ↳ completed receptor {rec.name}: +{inj['atoms_injected']} atoms")

                ligand_pdbqt_dir = get_pdbqt_dir(lig_staging)
                ligand_outputs = run_workflow(
                    input_dir=lig_staging, contains="ligands", output_dir=ligand_pdbqt_dir,
                    process_postfixes=False, convert_ligands_with_meeko=True, verbose=False)
                manifest = build_prepared_manifest(protein_outputs, ligand_outputs)

            if not manifest["proteins"] or not manifest["ligands"]:
                print(f"  ✗ [{conv_name}] manifest empty")
                failed_prep += 1
                continue

            vina_out.mkdir(parents=True, exist_ok=True)
            _, results = run_autodock_vina(
                base_dir=vina_out, log_dir=log_dir, prepared_manifest=manifest, cfg=cfg,
                cpu_model=cpu_model, protein_workflow_data=protein_outputs,
                ligand_workflow_data=ligand_outputs)

            with open(vina_out / "docking_summary.json", "w") as fh:
                json.dump(generate_summary(results, cfg), fh, indent=2, default=str)

            optimize_with_reporting(results, vina_out)

            for r in results:
                icon = "✓" if r.status == "success" else "✗"
                aff = f"{r.best_affinity:.2f}" if r.best_affinity else "N/A"
                print(f"  {icon} [{conv_name}] {r.num_poses} poses | best: {aff} kcal/mol")
            all_results.extend(results)

    n_ok = sum(1 for r in all_results if r.status == "success")
    n_fail = sum(1 for r in all_results if r.status == "failed")
    run_failed = bool(n_fail or failed_prep or issues)
    elapsed = time.time() - t_start
    write_status("failed" if run_failed else "completed",
                 finished_at=datetime.now().isoformat(), processed=len(ids),
                 docked_this_run=n_ok, failed_docking=n_fail, skipped=skipped,
                 failed_preparation=failed_prep, optimizer_issues=issues,
                 wall_clock_s=round(elapsed, 1))

    print(f"\n{'=' * 60}")
    print("ARM FAILED" if run_failed else "ARM COMPLETE")
    print(f"  Exhaustiveness: {exhaustiveness}")
    print(f"  Docked:         {n_ok}")
    print(f"  Failed dock:    {n_fail}")
    print(f"  Skipped:        {skipped}")
    print(f"  Failed prep:    {failed_prep}")
    print(f"  Optimizer issues: {len(issues)}")
    print(f"  Wall clock:     {elapsed / 3600:.2f} h")
    print(f"  Results:        {output_base}")
    print(f"  Status:         {status_path}")
    print(f"{'=' * 60}")
    if run_failed:
        raise RuntimeError(f"Arm finished with recorded failures; see {status_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
