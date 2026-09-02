#!/usr/bin/env python
"""DiffDock docking for ONE arm of the inference-parameter sweep.

A restricted, resumable clone of the authoritative benchmark DiffDock driver
(Master_Docking.ipynb, cell "Benchmark Set || DiffDock"). It differs from that
cell in exactly two respects: it docks only the complexes named in an explicit
ids file, and it applies one arm's declared parameter change. Preparation,
sampling and smina/gnina refinement are the same call into
Scripts/Docking/run_diffdock.py, so every arm is produced by identical code.

Arms are declared in Scripts/Docking/diffdock_param_sweep_config.yaml. An arm
may only change what it declares; the runner rebuilds the effective config from
the shared base and refuses to launch if anything else drifted.

Driven by DiffDock_Parameter_Sweep.ipynb.

Usage
-----
    python Scripts/Docking/run_diffdock_arm.py \
        --arm noise_hi \
        --ids-file Dockings/DiffDock_ParamSweep/param_sweep_sample_ids.txt

Re-running is safe and is the intended way to resume. A complex whose pose set is
COMPLETE is not re-docked; one that was interrupted part-way through writing its
poses is discarded and re-docked, because a partial set would otherwise be
treated as finished forever; a complex that failed is retried; and refinement
fills only what is missing.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
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

from Scripts.Docking.run_diffdock import (  # noqa: E402
    generate_summary,
    get_gpu_model,
    resolve_inference_args,
    run_diffdock,
)

DEFAULT_CONFIG = Path("Scripts/Docking/diffdock_param_sweep_config.yaml")

# Default retries after a driver-level CUDA fault. Zero: such a fault tracks host
# load rather than the complex, so a retry usually burns the time to fail again.
# Re-running the arm on a quiet machine is the cheaper recovery, and it re-docks
# exactly the complexes left without poses.
TRANSIENT_RETRIES = 0
_TRANSIENT_CUDA = ("unspecified launch failure", "CUDA driver error",
                   "an illegal memory access", "device-side assert")


def _has_transient_cuda_fault(complex_out: Path) -> bool:
    for log in (complex_out / "_batch_logs").glob("*stderr*.log"):
        text = log.read_text(errors="ignore")
        if any(m in text for m in _TRANSIENT_CUDA):
            return True
    return False

# What an arm is allowed to declare. Anything else in an arm block is a config
# slip, not an experiment.
ARM_TOP_LEVEL_KEYS = {"num_samples", "inference_steps",
                      "optimize_minimize_iters", "optimize_force_cap",
                      "optimize_search", "optimize_autobox_add"}
ARM_META_KEYS = {"label", "targets", "rationale", "overrides", "posthoc_of", "enabled"}

# Keys that must be identical across every arm. If one of these drifts the arms
# are no longer a controlled comparison, whatever the arm block says.
SHARED_CRITICAL_KEYS = ["diffdock_dir", "diffdock_python", "device",
                        "timeout_per_complex", "batch_mode", "batch_size",
                        "group_by_receptor", "max_workers", "optimization",
                        "optimize_cpu", "optimize_seed", "gnina_use_gpu",
                        "smina_executable", "gnina_executable"]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True,
                    help="Arm key from the config's arms: registry")
    ap.add_argument("-c", "--config", type=Path, default=DEFAULT_CONFIG,
                    help="Sweep config holding the shared base and the arms registry")
    ap.add_argument("--ids-file", required=True, type=Path,
                    help="Newline-separated <PDBID>_<LIG> ids to dock (the sampled cohort)")
    ap.add_argument("--benchmark-dir", type=Path, default=Path("Data/PoseBuster Benchmark Set"),
                    help="Directory holding one <PDBID>_<LIG>/ folder per benchmark complex")
    ap.add_argument("--seed", type=int, default=None,
                    help="Override inference_seed for this run (replicate draws)")
    ap.add_argument("--limit", type=int, default=0,
                    help="Dock only the first N ids (0 = all). For pilot runs.")
    ap.add_argument("--output-root", type=Path, default=None,
                    help="Override output_root; the arm writes to <root>/<arm>/")
    ap.add_argument("--transient-retries", type=int, default=TRANSIENT_RETRIES,
                    help="Retries after a driver-level CUDA fault (default %(default)s). "
                         "0 = record the complex and move on.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Resolve and print the effective config, then exit without docking")
    return ap.parse_args(argv)


def build_arm_config(base: dict, arm_key: str, arm: dict) -> dict:
    """Effective config for one arm: the shared base plus its declared change."""
    unknown = set(arm) - ARM_TOP_LEVEL_KEYS - ARM_META_KEYS
    if unknown:
        raise RuntimeError(
            f"arm {arm_key!r} declares keys it may not set: {sorted(unknown)}. "
            f"Allowed: {sorted(ARM_TOP_LEVEL_KEYS | ARM_META_KEYS)}")

    cfg = copy.deepcopy(base)
    cfg.pop("arms", None)

    for key in ARM_TOP_LEVEL_KEYS:
        if key in arm:
            cfg[key] = arm[key]

    # The base's inference_overrides are settings EVERY arm shares — an engine
    # setting like DiffDock's inference batch_size, which changes peak GPU memory
    # and must be identical across arms or they are not comparable. The arm's own
    # overrides are its one declared delta on top. An arm may not redeclare a base
    # key: that would be a second, undeclared change hiding inside the first.
    shared = dict(base.get("inference_overrides") or {})
    own = dict(arm.get("overrides") or {})
    clash = sorted(set(shared) & set(own))
    if clash:
        raise RuntimeError(
            f"arm {arm_key!r} redeclares shared inference setting(s) {clash}. Those are held "
            f"constant across every arm; an arm that changes one is varying two things at once.")
    cfg["inference_overrides"] = {**shared, **own}
    return cfg


def _sha_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def main(argv=None) -> int:
    args = parse_args(argv)

    base = yaml.safe_load(args.config.read_text())
    arms = base.get("arms") or {}
    if args.arm not in arms:
        raise RuntimeError(f"unknown arm {args.arm!r}; config declares {sorted(arms)}")
    arm = arms[args.arm] or {}

    if not arm.get("enabled", True):
        raise RuntimeError(
            f"arm {args.arm!r} is disabled in {args.config}. Set `enabled: true` on it to "
            f"run it; the arm block records why it was retired.")

    if arm.get("posthoc_of"):
        raise RuntimeError(
            f"arm {args.arm!r} is post-hoc (posthoc_of: {arm['posthoc_of']}) and needs no "
            "docking — re-refine that arm's existing poses from the notebook instead")

    cfg = build_arm_config(base, args.arm, arm)
    if args.seed is not None:
        cfg["inference_seed"] = int(args.seed)

    # ── Guard the scientific invariants of the sweep ─────────────────────────
    # Every arm must reach DiffDock through the seeded batch path, or its poses
    # are not reproducible and a re-run silently re-samples.
    if not cfg.get("batch_mode", True):
        raise RuntimeError("batch_mode must stay true; the non-batch path is not part of this sweep")
    if cfg.get("inference_seed") is None:
        raise RuntimeError("inference_seed must be set; unseeded arms are not comparable")

    # Rebuilding from the shared base is what makes the arms controlled, so
    # confirm the base itself carries every key the comparison depends on.
    missing = [k for k in SHARED_CRITICAL_KEYS if k not in cfg]
    if missing:
        raise RuntimeError(f"config is missing shared keys the sweep controls: {missing}")

    diffdock_dir = Path(cfg["diffdock_dir"])
    samples = int(cfg["num_samples"])
    steps = int(cfg["inference_steps"])

    # DiffDock's low-temperature branch is gated on `if temp_sampling[i] != 1.0`
    # (utils/sampling.py:173,178,183). At exactly 1.0 the whole branch is skipped
    # and temp_psi / temp_sigma_data for that axis are silently ignored, so an arm
    # tuning psi against a neutral temp_sampling would measure nothing.
    _ov = cfg.get("inference_overrides") or {}
    for axis in ("tr", "rot", "tor"):
        if any(f"temp_{k}_{axis}" in _ov for k in ("psi", "sigma_data")):
            eff = _ov.get(f"temp_sampling_{axis}",
                          resolve_inference_args(samples, steps, diffdock_dir, None)[f"temp_sampling_{axis}"])
            if float(eff) == 1.0:
                raise RuntimeError(
                    f"arm {args.arm!r} sets a temp_psi/temp_sigma_data on the {axis} axis, but "
                    f"temp_sampling_{axis} is exactly 1.0, which makes DiffDock skip the whole "
                    f"low-temperature branch for that axis and ignore the setting.")


    # DiffDock sizes its sampling tensors by its own inference batch_size, so a
    # sample count that is not a multiple of it dies with a tensor-shape error
    # that is caught per complex and turned into a silent failure.
    dd_batch = int(cfg.get("inference_overrides", {}).get("batch_size", 10))
    if samples % dd_batch:
        raise RuntimeError(
            f"num_samples ({samples}) must be a multiple of DiffDock's inference "
            f"batch_size ({dd_batch}) or complexes fail with a tensor-shape error")

    # Resolve the DiffDock inference args now so a bad override key fails here,
    # before a multi-hour run, and so the arm's real parameters are recorded.
    effective_inference = resolve_inference_args(
        samples, steps, diffdock_dir, cfg.get("inference_overrides"))
    # Diff against the SHARED BASE's sampling settings, not this arm's, so an arm
    # that changes num_samples or inference_steps shows that as its departure
    # rather than reporting none.
    baseline = resolve_inference_args(
        int(base["num_samples"]), int(base["inference_steps"]), diffdock_dir,
        base.get("inference_overrides"))
    declared = {k: v for k, v in effective_inference.items() if baseline.get(k) != v}

    # actual_steps is where the reverse diffusion really stops, and it is also a
    # legitimate knob. sampling.py:98 gives the LAST executed iteration a dt equal
    # to the whole remaining schedule rather than one increment, so the shipped
    # 19-of-20 makes that step apply about 1.5x the model's estimated remaining
    # translational offset and 2.9x the torsional one. Running 20 of 20 halves
    # both. Anything BELOW steps - 1 is different in kind: it halts the diffusion
    # early at a large sigma and still writes plausible-looking SDFs, leaving no
    # trace of the truncation.
    _as, _is = effective_inference["actual_steps"], effective_inference["inference_steps"]
    if _as not in (_is - 1, _is):
        raise RuntimeError(
            f"actual_steps ({_as}) must be inference_steps ({_is}) or one less; "
            f"a lower value silently truncates the reverse diffusion")

    output_root = Path(args.output_root or cfg["output_root"])
    arm_root = output_root / args.arm
    log_dir = Path(cfg.get("log_dir", "Dockings/Logs/diffdock_param_sweep_logs"))
    arm_root.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    ids = [ln.strip() for ln in args.ids_file.read_text().splitlines()
           if ln.strip() and not ln.lstrip().startswith("#")]
    if not ids:
        raise RuntimeError(f"No ids found in {args.ids_file}")
    if args.limit:
        ids = ids[: args.limit]

    absent = [pid for pid in ids if not (args.benchmark_dir / pid).is_dir()]
    if absent:
        raise RuntimeError(f"Refusing launch: ids absent from disk: {absent[:10]}")

    gpu_model = get_gpu_model()
    print("=" * 72)
    print(f"DiffDock parameter sweep — arm {args.arm!r}  ({arm.get('label', '')})")
    print("=" * 72)
    print(f"  Complexes      : {len(ids)}")
    print(f"  Samples/steps  : {samples} / {steps} (actual {effective_inference['actual_steps']})")
    print(f"  Seed           : {cfg['inference_seed']}")
    print(f"  Refinement     : {cfg.get('optimization')}"
          + (f"  (minimize_iters={cfg.get('optimize_minimize_iters')},"
             f" force_cap={cfg.get('optimize_force_cap')})"
             if cfg.get("optimize_minimize_iters") or cfg.get("optimize_force_cap") else ""))
    print(f"  Output         : {arm_root}")
    print(f"  GPU            : {gpu_model}")
    if declared:
        print("  Departures from the shared base:")
        for k in sorted(declared):
            print(f"    {k}: {baseline.get(k)}  ->  {declared[k]}")
    else:
        print("  Departures from the shared base: none (control arm)")
    print()

    if args.dry_run:
        print(json.dumps({"effective_config": {k: v for k, v in cfg.items() if k != "arms"},
                          "effective_inference_args": effective_inference}, indent=2, default=str))
        return 0

    status_path = arm_root / "param_sweep_arm_status.json"
    ids_sha = _sha_text("\n".join(ids) + "\n")

    def write_status(state, **extra):
        payload = {
            "schema_version": 1,
            "arm": args.arm,
            "label": arm.get("label"),
            "targets": arm.get("targets"),
            "config": str(args.config),
            "ids_file": str(args.ids_file),
            "ids_sha256": ids_sha,
            "n_ids": len(ids),
            "num_samples": samples,
            "inference_steps": steps,
            "inference_seed": cfg["inference_seed"],
            "inference_overrides": cfg.get("inference_overrides") or {},
            "declared_departures": declared,
            "optimization": cfg.get("optimization"),
            "optimize_minimize_iters": cfg.get("optimize_minimize_iters"),
            "optimize_force_cap": cfg.get("optimize_force_cap"),
            "gpu_model": gpu_model,
            "status": state,
            "updated_at": datetime.now().isoformat(),
            **extra,
        }
        tmp = status_path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(payload, indent=2, default=str) + "\n")
        tmp.replace(status_path)

    write_status("running", processed=0, last_complex=None)

    docked = skipped = failed = 0
    short: list[dict] = []
    reset: list[dict] = []
    transient: list[dict] = []
    durations: list[float] = []   # docked complexes only, for the rate estimate
    cuda_failed: list[dict] = []  # lost to a driver fault, recoverable by re-running
    retries = max(0, int(args.transient_retries))
    oom_suspects: list[str] = []
    t_start = time.time()

    for idx, pdb_id in enumerate(ids, 1):
        write_status("running", processed=idx - 1, last_complex=pdb_id)

        cdir = args.benchmark_dir / pdb_id
        protein_pdb = cdir / f"{pdb_id}_protein.pdb"
        ligand_start = cdir / f"{pdb_id}_ligand_start_conf.sdf"
        if not protein_pdb.exists() or not ligand_start.exists():
            raise RuntimeError(f"{pdb_id}: missing inputs — the cohort must be complete")

        complex_out = arm_root / pdb_id

        # Resume repair. run_diffdock treats a combo directory holding ANY pose as
        # finished, so a complex interrupted while its poses were being written
        # stays half-docked forever: it is skipped on every later run, and the
        # completeness check below then fails the arm on every later run too.
        # Re-docking is deterministic under the fixed seed, so discarding the
        # partial complex costs one complex and resolves it.
        existing_raw = [p for p in complex_out.rglob("rank*_confidence*.sdf")
                        if not p.parent.name.startswith("optimized_")] if complex_out.is_dir() else []
        if existing_raw and len(existing_raw) != samples:
            print(f"  ↻ {pdb_id}: {len(existing_raw)} of {samples} poses on disk — "
                  f"discarding the partial complex and re-docking it")
            shutil.rmtree(complex_out, ignore_errors=True)
            reset.append({"pdb_id": pdb_id, "found": len(existing_raw), "expected": samples})

        # run_diffdock initialises its docking log before it prepares inputs, so
        # the per-complex directory has to exist first.
        complex_out.mkdir(parents=True, exist_ok=True)
        # run_diffdock never rewrites a pre-retry stderr log, and each retry
        # stage writes under a different name, so a RECOVERED out-of-memory
        # leaves its message on disk permanently. Clear the logs first and scan
        # only what this invocation wrote.
        shutil.rmtree(complex_out / "_batch_logs", ignore_errors=True)

        # Two lines per complex: one when it starts, so the in-flight complex is
        # visible during a multi-hour arm, and one when it finishes carrying the
        # duration. Both flush, because the notebook reads this stream live.
        n_prot = sum(1 for ln in protein_pdb.read_text().splitlines()
                     if ln.startswith("ATOM"))
        print(f"[{args.arm}] {idx:>3}/{len(ids)} {pdb_id:<12} docking "
              f"({n_prot:,} receptor atoms) ...", flush=True)
        t_complex = time.time()

        # run_diffdock prepares inputs, samples, and runs the smina/gnina
        # refinement, writing docking_log.csv and optimization_log.csv into
        # complex_out. One protein and one ligand keeps its cross product to the
        # single intended combination.
        #
        # A driver-level CUDA fault ("unspecified launch failure") is not caused by
        # the complex — 6YR2_T1C failed here repeatedly yet docked in 332 s in the
        # published campaign on this same card and parameters. It tracks host load
        # instead. Retrying therefore spends 2-3x the time to fail again while the
        # machine is busy, so the default is NOT to retry: record the complex, move
        # on, and re-run the arm when the box is quiet, which re-docks exactly the
        # complexes that have no poses. Raise --transient-retries only if you have
        # reason to believe a retry would land differently.
        for attempt in range(1, retries + 2):
            results = run_diffdock(
                proteins=[protein_pdb],
                ligands=[ligand_start],
                output_dir=complex_out,
                cfg=cfg,
                gpu_model=gpu_model,
            )
            produced = len([p for p in complex_out.rglob("rank*_confidence*.sdf")
                            if not p.parent.name.startswith("optimized_")])
            if produced or attempt > retries or not _has_transient_cuda_fault(complex_out):
                break
            print(f"  ↻ {pdb_id}: CUDA fault on attempt {attempt}; clearing and retrying")
            transient.append({"pdb_id": pdb_id, "attempt": attempt})
            shutil.rmtree(complex_out, ignore_errors=True)
            complex_out.mkdir(parents=True, exist_ok=True)
            time.sleep(20)      # let the driver settle before the next launch

        (complex_out / "docking_summary.json").write_text(
            json.dumps(generate_summary(results, cfg), indent=2, default=str) + "\n")

        # A complex that quietly produced fewer poses than the arm asked for is a
        # truncated run, not a success, and it would silently shrink this arm's
        # pool relative to every other arm.
        n_poses = len([p for p in complex_out.rglob("rank*_confidence*.sdf")
                       if not p.parent.name.startswith("optimized_")])
        lost_to_cuda = n_poses == 0 and _has_transient_cuda_fault(complex_out)
        if lost_to_cuda:
            cuda_failed.append({"pdb_id": pdb_id, "attempts": retries + 1})
            print(f"  ✗ {pdb_id}: lost to a CUDA driver fault, not retried. Re-run this arm "
                  f"on a quiet machine to pick it up.")
        elif n_poses != samples:
            short.append({"pdb_id": pdb_id, "poses": n_poses, "expected": samples})
            print(f"  ⚠ {pdb_id}: {n_poses} poses, expected {samples}")

        # The runner's OOM retry matches only two of the three strings CUDA
        # raises, and its CPU fallback keeps logging device cuda:0. Either turns
        # a lost or slow complex into a plausible-looking row.
        for log in (complex_out / "_batch_logs").glob("*stderr.log"):
            text = log.read_text(errors="ignore")
            if not n_poses and any(m in text for m in ("OutOfMemoryError", "CUDA out of memory",
                                                       "CUDA error: out of memory")):
                oom_suspects.append(str(log.relative_to(arm_root)))
                print(f"  ⚠ {pdb_id}: CUDA out-of-memory recorded in {log.name}")

        was_skipped = all(r.status == "skipped" for r in results) if results else False
        for r in results:
            if r.status == "success":
                docked += 1
            elif r.status == "skipped":
                skipped += 1
            else:
                failed += 1
                print(f"  ✗ {r.ligand_name}: {r.error_message[:160]}")

        # Closing line for this complex. Only complexes actually docked feed the
        # rate estimate — a skipped one returns instantly and would otherwise make
        # the remaining time look far shorter than it is.
        dt = time.time() - t_complex
        if not results:
            outcome = "no result"
        elif was_skipped:
            outcome = f"already done ({n_poses} poses)"
        elif n_poses == samples:
            outcome = f"{n_poses} poses"
            durations.append(dt)
        elif lost_to_cuda:
            outcome = "CUDA fault, skipped"
        else:
            outcome = f"FAILED ({n_poses} poses)"
        eta = ""
        if durations:
            remaining = sum(1 for j in range(idx, len(ids))
                            if not (arm_root / ids[j]).is_dir()
                            or len([p for p in (arm_root / ids[j]).rglob("rank*_confidence*.sdf")
                                    if not p.parent.name.startswith("optimized_")]) != samples)
            eta = (f"  eta {remaining * (sum(durations) / len(durations)) / 60:.0f} min"
                   f" for {remaining} left")
        print(f"[{args.arm}] {idx:>3}/{len(ids)} {pdb_id:<12} {outcome:<22} "
              f"{dt:>6.0f}s   elapsed {(time.time() - t_start) / 60:>5.0f} min{eta}",
              flush=True)

    elapsed = time.time() - t_start
    # A CUDA-lost complex is recoverable by re-running and must not halt a sweep
    # that still has arms to go. But an arm missing complexes is NOT "completed":
    # the notebook skips a completed arm whose cohort hash matches, so reporting
    # success here would strand the gaps permanently. Three states are needed.
    incomplete = [pid for pid in ids
                  if len([q for q in (arm_root / pid).rglob("rank*_confidence*.sdf")
                          if not q.parent.name.startswith("optimized_")]) != samples]
    run_failed = bool(short or oom_suspects or (failed - len(cuda_failed)) > 0)
    if run_failed:
        final_status = "failed"
    elif incomplete:
        final_status = "incomplete"      # recoverable: re-run picks up exactly these
    else:
        final_status = "completed"
    write_status(final_status,
                 finished_at=datetime.now().isoformat(), processed=len(ids),
                 docked_this_run=docked, skipped=skipped, failed=failed,
                 short_pose_complexes=short, oom_suspects=oom_suspects,
                 partial_complexes_redocked=reset,
                 transient_cuda_retries=transient,
                 cuda_failed_complexes=cuda_failed,
                 incomplete_complexes=incomplete,
                 wall_clock_s=round(elapsed, 1))

    print(f"\n{'=' * 72}")
    print({"failed": "ARM FAILED", "incomplete": "ARM INCOMPLETE — re-run to fill the gaps",
           "completed": "ARM COMPLETE"}[final_status])
    print(f"  Arm:        {args.arm}  ({arm.get('label', '')})")
    print(f"  Docked:     {docked}")
    print(f"  Skipped:    {skipped}  (already present)")
    print(f"  Failed:     {failed}")
    print(f"  Re-docked partials: {len(reset)}")
    print(f"  Transient CUDA retries: {len(transient)}  (retries allowed: {retries})")
    print(f"  Complexes still missing: {len(incomplete)}"
          + (f"  -> {incomplete}" if incomplete else ""))
    print(f"  Lost to CUDA faults:    {len(cuda_failed)}"
          + (f"  -> {[c['pdb_id'] for c in cuda_failed]}" if cuda_failed else ""))
    print(f"  Short pose count: {len(short)}")
    print(f"  OOM suspects:     {len(oom_suspects)}")
    print(f"  Wall clock: {elapsed / 3600:.2f} h")
    print(f"  Results:    {arm_root}")
    print(f"  Status:     {status_path}")
    print(f"{'=' * 72}")
    if run_failed:
        raise RuntimeError(f"Arm finished with recorded failures; see {status_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
