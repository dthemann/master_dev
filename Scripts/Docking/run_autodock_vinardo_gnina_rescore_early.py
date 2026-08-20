#!/usr/bin/env python3
"""Run Vinardo GNINA-rescore early while the main sequence refines Vina poses.

The authoritative four-stage sequence remains the owner of the final stage
transition.  This helper executes the exact pinned Vinardo-rescore notebook
cell in optimizer-only mode: all 308 raw Vinardo outputs are strictly validated
before launch and any attempted docking call fails closed.  Work is traversed
in descending complex order, while the later authoritative stage traverses
ascending, and a separate watchdog stops this process group at a conservative
Vina-refinement threshold.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import sys
import traceback
from typing import Any

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Scripts.Docking import run_autodock_vinardo_raw_prefill as guard
from Scripts.Docking.run_autodock import validate_vina_output


STAGE_NAME = "early_vinardo_gnina_rescore"
DEFAULT_LOCK = Path("/tmp/master_docking_autodock_vinardo_gnina_rescore_early.lock")
LAUNCHER_STATUS_NAME = "benchmark_early_vinardo_gnina_rescore_launcher_status.json"
WATCHDOG_STATUS_NAME = "benchmark_early_vinardo_gnina_rescore_watchdog_status.json"
LOG_PREFIX = "autodock_vinardo_gnina_rescore_early"
EXPECTED_COMPLEXES = 308


def _validate_complete_raw_vinardo(context: dict[str, Any]) -> dict[str, int]:
    """Prove that optimizer-only execution cannot need a docking fallback."""
    config = yaml.safe_load(Path(context["config_path"]).read_text())
    prep_tool = str(config.get("prep_tool", "mgltools")).strip().lower()
    converter_names: list[str] = []
    if prep_tool in {"mgltools", "both"}:
        converter_names.append("mgl_tools")
    if prep_tool in {"meeko", "both"}:
        converter_names.append("meeko")
    if converter_names != ["mgl_tools"]:
        raise guard.PrefillGuardError(
            "early Vinardo rescore requires exactly the pinned mgl_tools converter"
        )

    benchmark_dir = guard.REPO_ROOT / "Data/PoseBuster Benchmark Set"
    ids_path = benchmark_dir / "posebusters_pdb_ccd_ids.txt"
    authoritative_ids = sorted({
        line.strip() for line in ids_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    })
    if len(authoritative_ids) != EXPECTED_COMPLEXES:
        raise guard.PrefillGuardError(
            f"expected {EXPECTED_COMPLEXES} authoritative IDs, found {len(authoritative_ids)}"
        )

    output_dir = Path(context["output_dir"])
    pose_count = 0
    for complex_id in authoritative_ids:
        if not (benchmark_dir / complex_id).is_dir():
            raise guard.PrefillGuardError(
                f"authoritative benchmark directory is missing: {complex_id}"
            )
        docking_dir = output_dir / complex_id / "mgl_tools" / "docking"
        outputs = sorted(docking_dir.glob("*_vina_out.pdbqt"))
        if len(outputs) != 1:
            raise guard.PrefillGuardError(
                f"{complex_id}: expected one raw Vinardo output, found {len(outputs)}"
            )
        pose_file = outputs[0]
        if not pose_file.stem.endswith("_vinardo_vina_out"):
            raise guard.PrefillGuardError(
                f"{complex_id}: raw output is not explicitly Vinardo-scored"
            )
        try:
            models = validate_vina_output(pose_file)
        except (OSError, UnicodeError, ValueError) as exc:
            raise guard.PrefillGuardError(
                f"{complex_id}: invalid raw Vinardo output: {exc}"
            ) from exc
        receptor_stem, separator, _ = pose_file.stem.partition("__")
        receptor = (
            output_dir / complex_id / "_staging" / "receptors" / "pdbqt"
            / f"{receptor_stem}.pdbqt"
        )
        if not separator or not receptor.is_file():
            raise guard.PrefillGuardError(
                f"{complex_id}: staged prepared receptor is missing"
            )
        pose_count += len(models)
    return {"validated_raw_complexes": len(authoritative_ids), "validated_raw_poses": pose_count}


def _transform_rescore_source(source: str) -> str:
    optimizer_block = '''OPT_PREFLIGHT = preflight_optimizers(base_cfg)
print("Variant A: empirical minimization + GNINA CNN rescore "
      f"(rank by {base_cfg.get('optimize_rank_by', 'minimized_affinity')}; "
      f"overwrite optimization={OVERWRITE_OPTIMIZATION})")'''
    guarded_block = '''OPT_PREFLIGHT = preflight_optimizers(base_cfg)
def _early_rescore_docking_forbidden(*args, **kwargs):
    raise RuntimeError("early Vinardo GNINA rescore attempted to enter a docking path")
run_autodock_vina = _early_rescore_docking_forbidden
print("Early Vinardo GNINA rescore (raw docking is validated complete and forbidden)")'''
    source = guard._replace_once(
        source, optimizer_block, guarded_block, "optimizer-only guard",
    )
    source = guard._replace_once(
        source,
        'print(f"Found exactly {len(complex_dirs)} authoritative benchmark complexes\\n")',
        'complex_dirs = list(reversed(complex_dirs))\n'
        'print(f"Found exactly {len(complex_dirs)} authoritative benchmark complexes; "\n'
        '      "early rescore traverses them in descending order\\n")',
        "reverse traversal",
    )
    source = guard._replace_once(
        source,
        '"rank_fields": ["autodock_rank", "optimized_rank"],',
        '"rank_fields": ["autodock_rank", "optimized_rank"],\n'
        '        "early_rescore": True,\n'
        '        "traversal_order": "descending",\n'
        '        "stop_at_refinement_complexes": _EARLY_RESCORE_STOP_THRESHOLD,',
        "status execution metadata",
    )
    source = guard._replace_once(
        source,
        "for idx, cdir in enumerate(complex_dirs, 1):\n"
        "    pdb_id         = cdir.name",
        "for idx, cdir in enumerate(complex_dirs, 1):\n"
        "    _assert_early_rescore_window()\n"
        "    pdb_id         = cdir.name",
        "complex-boundary handoff",
    )
    source = guard._replace_once(
        source,
        "        for cn, _ in converters:\n"
        "            optimize_existing_poses(output_base / pdb_id / cn)",
        "        for cn, _ in converters:\n"
        "            _assert_early_rescore_window()\n"
        "            optimize_existing_poses(output_base / pdb_id / cn)",
        "completed-complex optimization handoff",
    )
    source = guard._replace_once(
        source,
        "            optimize_existing_poses(vina_out)\n"
        "            continue",
        "            _assert_early_rescore_window()\n"
        "            optimize_existing_poses(vina_out)\n"
        "            continue",
        "converter optimization handoff",
    )
    compile(source, f"early-rescore#{guard.VINARDO_CELL_ID}", "exec")
    return source


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-status", type=Path, default=guard.DEFAULT_SEQUENCE_STATUS)
    parser.add_argument("--lock-file", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--launcher-status", type=Path, default=None)
    parser.add_argument("--log-path", type=Path, default=None)
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    sequence_status_path = args.sequence_status.expanduser().resolve()
    context = guard._load_launch_context(sequence_status_path)
    raw_validation = _validate_complete_raw_vinardo(context)
    transformed_source = _transform_rescore_source(context["source"])
    transformed_sha = guard._sha256_bytes(transformed_source.encode())
    output_dir = Path(context["output_dir"])
    launcher_status_path = (
        args.launcher_status.expanduser().resolve()
        if args.launcher_status else output_dir / LAUNCHER_STATUS_NAME
    )
    rescore_status_path = output_dir / "benchmark_gnina_rescore_status.json"
    lock_path = args.lock_file.expanduser().resolve()

    plan: dict[str, Any] = {
        "status": "ready" if args.verify_only else "running",
        "stage": STAGE_NAME,
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "process_start_time": guard._process_start_time(os.getpid()),
        "started_at": guard._now(),
        "updated_at": guard._now(),
        "sequence_status": str(sequence_status_path),
        "sequence_pid": context["sequence_pid"],
        "sequence_start_time": context["sequence_start_time"],
        "sequence_process_group_id": context["sequence_process_group_id"],
        "sequence_lock_file": str(context["sequence_lock_path"]),
        "notebook": str(context["notebook_path"]),
        "source_cell_id": guard.VINARDO_CELL_ID,
        "source_cell_sha256": context["source_sha256"],
        "transformed_source_sha256": transformed_sha,
        "vinardo_config": str(context["config_path"]),
        "vinardo_config_sha256": guard._sha256_file(context["config_path"]),
        "launcher_status_path": str(launcher_status_path),
        "rescore_status_path": str(rescore_status_path),
        "refinement_status_path": str(context["refinement_status_path"]),
        "log_path": str(args.log_path.expanduser().resolve()) if args.log_path else None,
        "lock_file": str(lock_path),
        "optimization": "gnina",
        "cnn_scoring": "rescore",
        "traversal_order": "descending",
        "stop_at_refinement_complexes": guard.STOP_AT_REFINEMENT_COMPLEXES,
        "refinement_complexes_at_launch": context["completed_refinement"],
        "target_script": Path(__file__).name,
        **raw_validation,
    }
    if args.verify_only:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0
    if os.getpgrp() != os.getpid():
        raise guard.PrefillGuardError("early rescore must be launched with setsid")

    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        lock_handle.close()
        raise guard.PrefillGuardError(f"another early rescore holds {lock_path}") from exc
    lock_handle.seek(0)
    lock_handle.truncate()
    json.dump({"pid": os.getpid(), "started_at": plan["started_at"]}, lock_handle)
    lock_handle.write("\n")
    lock_handle.flush()
    guard._atomic_write_json(launcher_status_path, plan)

    def commit_terminal(state: str, **details: Any) -> None:
        plan.update({
            "status": state, "updated_at": guard._now(),
            "finished_at": guard._now(), **details,
        })
        guard._atomic_write_json(launcher_status_path, plan)

    def signal_handler(signum: int, _frame: Any) -> None:
        raise guard.PrefillSignal(signum)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    def assert_early_rescore_window() -> None:
        sequence = guard._read_json(sequence_status_path)
        live_start = guard._process_start_time(context["sequence_pid"])
        if (
            sequence.get("status") != "running"
            or sequence.get("current_stage") != guard.REQUIRED_SEQUENCE_STAGE
            or sequence.get("pid") != context["sequence_pid"]
            or live_start != context["sequence_start_time"]
        ):
            raise guard.PrefillHandoff(
                "authoritative sequence left the Vina-refinement ownership window"
            )
        refinement = guard._read_json(context["refinement_status_path"])
        try:
            completed = int(refinement.get("completed_complexes") or 0)
        except (TypeError, ValueError) as exc:
            raise guard.PrefillHandoff("Vina refinement checkpoint became invalid") from exc
        if (
            refinement.get("status") != "running"
            or completed >= guard.STOP_AT_REFINEMENT_COMPLEXES
        ):
            raise guard.PrefillHandoff(
                "conservative refinement threshold reached "
                f"(completed={completed}, threshold={guard.STOP_AT_REFINEMENT_COMPLEXES})"
            )

    namespace = {
        "__name__": "__main__",
        "__file__": f"{context['notebook_path']}#{guard.VINARDO_CELL_ID}:early-rescore",
        "__package__": None,
        "_EARLY_RESCORE_STOP_THRESHOLD": guard.STOP_AT_REFINEMENT_COMPLEXES,
        "_assert_early_rescore_window": assert_early_rescore_window,
    }
    previous_cwd = Path.cwd()
    os.chdir(guard.REPO_ROOT)
    try:
        exec(
            compile(
                transformed_source,
                f"{context['notebook_path']}#{guard.VINARDO_CELL_ID}:early-rescore",
                "exec",
            ),
            namespace,
            namespace,
        )
    except guard.PrefillHandoff as exc:
        guard._write_raw_terminal_status(rescore_status_path, "handed_off", str(exc))
        commit_terminal("handed_off", handoff_reason=str(exc))
        print(f"EARLY VINARDO RESCORE HANDED OFF: {exc}", flush=True)
        return 0
    except guard.PrefillSignal as exc:
        guard._write_raw_terminal_status(rescore_status_path, "interrupted", str(exc))
        commit_terminal("interrupted", signal=exc.signum, error=str(exc))
        print(f"EARLY VINARDO RESCORE INTERRUPTED: {exc}", file=sys.stderr, flush=True)
        return 128 + exc.signum
    except BaseException as exc:
        guard._write_raw_terminal_status(
            rescore_status_path, "failed", f"{type(exc).__name__}: {exc}",
        )
        commit_terminal(
            "failed", error_type=type(exc).__name__, error=str(exc),
            traceback=traceback.format_exc(),
        )
        print(f"EARLY VINARDO RESCORE FAILED: {type(exc).__name__}: {exc}",
              file=sys.stderr, flush=True)
        traceback.print_exc()
        return 1
    finally:
        os.chdir(previous_cwd)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()

    commit_terminal("completed")
    print("EARLY VINARDO GNINA RESCORE COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
