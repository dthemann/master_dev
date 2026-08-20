#!/usr/bin/env python3
"""Prefill raw benchmark Vinardo dockings while Vina refinement is running.

This runner deliberately executes the reviewed Vinardo notebook cell with only
its raw AutoDock/Vina work enabled.  GNINA optimization is left to the existing
four-stage sequence.  The transformation is fail-closed: it is accepted only
when the notebook cell and Vinardo YAML still match the hashes pinned by the
live sequence, and every source edit anchor occurs exactly once.

The prefill writes a separate status file and owns a separate process lock.  At
each complex boundary it verifies that the authoritative sequence is still in
the Vina-refinement stage.  If the sequence advances, the prefill hands off
cleanly instead of competing with the scheduled Vinardo stage.
"""

from __future__ import annotations

import argparse
import datetime as dt
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import sys
import tempfile
import traceback
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
VINARDO_CELL_ID = "a3b4dd6c"
VINARDO_STAGE_NAME = "vinardo_gnina_rescore"
PREFILL_STAGE_NAME = "raw_vinardo_prefill"
REQUIRED_SEQUENCE_STAGE = "gnina_refinement"
# This is the source reviewed for the raw-only transformation below.  Requiring
# both this digest and the live sequence's digest prevents a newly pinned cell
# from gaining an optimizer path that the transformer does not understand.
REVIEWED_VINARDO_CELL_SHA256 = (
    "af983718f3b00348005fa9ca7c341ab3ad17071a4f403d3e36f3f62e80bb0661"
)
STOP_AT_REFINEMENT_COMPLEXES = 250

DEFAULT_SEQUENCE_STATUS = (
    REPO_ROOT
    / "Dockings/vina_results_full_protein_vina_scoring"
    / "benchmark_gnina_vina_vinardo_sequence_status.json"
)
DEFAULT_LOCK = Path("/tmp/master_docking_autodock_vinardo_raw_prefill.lock")
DEFAULT_SEQUENCE_LOCK = Path("/tmp/master_docking_autodock_gnina_sequence.lock")
PREFILL_STATUS_NAME = "benchmark_raw_vinardo_prefill_status.json"
LAUNCHER_STATUS_NAME = "benchmark_raw_vinardo_prefill_launcher_status.json"
WATCHDOG_STATUS_NAME = "benchmark_raw_vinardo_prefill_watchdog_status.json"


class PrefillGuardError(RuntimeError):
    """The pinned sequence, notebook, config, or source transform is unsafe."""


class PrefillHandoff(RuntimeError):
    """The main sequence advanced, so raw-prefill ownership has ended."""


class PrefillSignal(BaseException):
    """Termination request that notebook ``except Exception`` blocks cannot hide."""

    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"received {signal.Signals(signum).name}")


def _now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _read_json(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except Exception as exc:
        raise PrefillGuardError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise PrefillGuardError(f"expected a JSON object in {path}")
    return payload


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        Path(temporary_name).unlink(missing_ok=True)


def _process_start_time(pid: int) -> str | None:
    """Return Linux /proc start ticks, handling command names with spaces."""
    try:
        raw = Path(f"/proc/{pid}/stat").read_text()
        close = raw.rfind(")")
        if close < 0:
            return None
        # Fields after ``comm`` begin with field 3 (state); starttime is field
        # 22, therefore offset 19 in this suffix.
        fields = raw[close + 2 :].split()
        return fields[19]
    except (OSError, IndexError, ValueError):
        return None


def _process_cmdline(pid: int) -> list[str]:
    try:
        return [
            item.decode(errors="replace")
            for item in Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
            if item
        ]
    except OSError:
        return []


def _lock_is_held(path: Path) -> bool:
    """Probe an advisory lock without disturbing its current owner."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a+") as handle:
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
        return False


def _replace_once(source: str, old: str, new: str, label: str) -> str:
    count = source.count(old)
    if count != 1:
        raise PrefillGuardError(
            f"raw-prefill transform anchor {label!r} occurs {count} times, expected 1"
        )
    return source.replace(old, new, 1)


def _transform_cell_source(source: str) -> str:
    """Disable only optimizer execution and isolate status for raw prefill."""
    optimizer_block = '''OPT_TOOLS   = resolve_optimizers(base_cfg)
if OPT_TOOLS != ["gnina"]:
    raise RuntimeError(f"Variant A resolved unexpected optimizer aliases: {OPT_TOOLS}")
OPT_PREFLIGHT = preflight_optimizers(base_cfg)
print("Variant A: empirical minimization + GNINA CNN rescore "
      f"(rank by {base_cfg.get('optimize_rank_by', 'minimized_affinity')}; "
      f"overwrite optimization={OVERWRITE_OPTIMIZATION})")'''
    raw_block = '''# Raw-only prefill: retain the exact pinned docking configuration, but leave
# every GNINA pose to the authoritative sequence and its provenance namespace.
OPT_TOOLS = ()
OPT_PREFLIGHT = None
def _raw_prefill_optimizer_forbidden(*args, **kwargs):
    raise RuntimeError("raw Vinardo prefill attempted to enter an optimizer path")
optimize_autodock_results = _raw_prefill_optimizer_forbidden
preflight_optimizers = _raw_prefill_optimizer_forbidden
print("Raw Vinardo docking prefill (GNINA deliberately deferred to sequence)")'''
    source = _replace_once(source, optimizer_block, raw_block, "optimizer block")
    source = _replace_once(
        source,
        'RESCORE_STATUS_PATH = output_base / "benchmark_gnina_rescore_status.json"',
        f'RESCORE_STATUS_PATH = output_base / "{PREFILL_STATUS_NAME}"',
        "status path",
    )
    source = _replace_once(
        source, '"stage": "gnina_rescore",',
        '"stage": "raw_vinardo_prefill",', "status stage",
    )
    source = _replace_once(
        source, '"optimization": "gnina",',
        '"optimization": "none",', "status optimization",
    )
    source = _replace_once(
        source, '"cnn_scoring": "rescore",',
        '"cnn_scoring": "none",', "status CNN mode",
    )
    source = _replace_once(
        source, '"rank_fields": ["autodock_rank", "optimized_rank"],',
        '"rank_fields": ["autodock_rank"],\n'
        '        "source_cell_id": _RAW_PREFILL_CELL_ID,\n'
        '        "source_cell_sha256": _RAW_PREFILL_SOURCE_SHA256,\n'
        '        "traversal_order": "descending",\n'
        '        "stop_at_refinement_complexes": _RAW_PREFILL_STOP_THRESHOLD,',
        "status provenance",
    )
    source = _replace_once(
        source,
        'print(f"Found exactly {len(complex_dirs)} authoritative benchmark complexes\\n")',
        'complex_dirs = list(reversed(complex_dirs))\n'
        'print(f"Found exactly {len(complex_dirs)} authoritative benchmark complexes; "\n'
        '      "raw prefill traverses them in descending order\\n")',
        "reverse traversal",
    )
    source = _replace_once(
        source,
        "for idx, cdir in enumerate(complex_dirs, 1):\n"
        "    pdb_id         = cdir.name",
        "for idx, cdir in enumerate(complex_dirs, 1):\n"
        "    _assert_raw_prefill_window()\n"
        "    pdb_id         = cdir.name",
        "complex-boundary handoff",
    )
    source = _replace_once(
        source,
        "        results_df, results = run_autodock_vina(\n",
        "        _assert_raw_prefill_window()\n"
        "        results_df, results = run_autodock_vina(\n",
        "pre-docking handoff",
    )
    compile(source, f"raw-prefill#{VINARDO_CELL_ID}", "exec")
    return source


def _load_launch_context(sequence_status_path: Path) -> dict[str, Any]:
    sequence = _read_json(sequence_status_path)
    if sequence.get("status") != "running":
        raise PrefillGuardError("authoritative AutoDock sequence is not running")
    if sequence.get("current_stage") != REQUIRED_SEQUENCE_STAGE:
        raise PrefillGuardError(
            "raw Vinardo prefill may start only while the sequence is in "
            f"{REQUIRED_SEQUENCE_STAGE!r}; current={sequence.get('current_stage')!r}"
        )

    try:
        sequence_pid = int(sequence["pid"])
    except (KeyError, TypeError, ValueError) as exc:
        raise PrefillGuardError("sequence status has no valid PID") from exc
    sequence_start_time = _process_start_time(sequence_pid)
    if sequence_start_time is None:
        raise PrefillGuardError(f"sequence PID {sequence_pid} is not alive")
    cmdline = _process_cmdline(sequence_pid)
    if not cmdline or not any(
        Path(value).name == "run_autodock_gnina_sequence.py" for value in cmdline
    ):
        raise PrefillGuardError("sequence PID is not the guarded AutoDock launcher")
    try:
        status_arg_index = cmdline.index("--status-path") + 1
        command_status_path = Path(cmdline[status_arg_index]).expanduser().resolve()
    except (ValueError, IndexError) as exc:
        raise PrefillGuardError("sequence command does not pin --status-path") from exc
    if command_status_path != sequence_status_path:
        raise PrefillGuardError("sequence PID command/status path identity mismatch")
    try:
        expected_pgid = int(sequence["process_group_id"])
        actual_pgid = os.getpgid(sequence_pid)
    except (KeyError, TypeError, ValueError, ProcessLookupError) as exc:
        raise PrefillGuardError("sequence process-group identity is invalid") from exc
    if expected_pgid != sequence_pid or actual_pgid != expected_pgid:
        raise PrefillGuardError("sequence PID/process-group identity mismatch")
    try:
        lock_arg_index = cmdline.index("--lock-file") + 1
        sequence_lock_path = Path(cmdline[lock_arg_index]).expanduser().resolve()
    except ValueError:
        sequence_lock_path = DEFAULT_SEQUENCE_LOCK.resolve()
    except IndexError as exc:
        raise PrefillGuardError("sequence command has an empty --lock-file") from exc
    if not _lock_is_held(sequence_lock_path):
        raise PrefillGuardError("authoritative sequence does not hold its campaign lock")

    stages = sequence.get("stages") or []
    matches = [s for s in stages if s.get("name") == VINARDO_STAGE_NAME]
    if len(matches) != 1 or matches[0].get("cell_id") != VINARDO_CELL_ID:
        raise PrefillGuardError("sequence does not pin the expected Vinardo stage/cell")
    stage = matches[0]
    if stage.get("status") != "pending":
        raise PrefillGuardError(
            f"scheduled Vinardo stage is not pending: {stage.get('status')!r}"
        )
    expected_source_sha = str(stage.get("source_sha256") or "")
    if len(expected_source_sha) != 64:
        raise PrefillGuardError("sequence Vinardo cell hash is missing")
    if expected_source_sha != REVIEWED_VINARDO_CELL_SHA256:
        raise PrefillGuardError("sequence Vinardo cell is not the reviewed source revision")

    notebook_path = Path(str(sequence.get("notebook") or "")).expanduser().resolve()
    config_path = Path(str(sequence.get("vinardo_config") or "")).expanduser().resolve()
    if not notebook_path.is_file() or not config_path.is_file():
        raise PrefillGuardError("pinned notebook or Vinardo config is missing")
    if _sha256_file(config_path) != str(sequence.get("vinardo_config_sha256") or ""):
        raise PrefillGuardError("Vinardo YAML changed after sequence launch")

    notebook = _read_json(notebook_path)
    cells = [c for c in notebook.get("cells") or [] if c.get("id") == VINARDO_CELL_ID]
    if len(cells) != 1 or cells[0].get("cell_type") != "code":
        raise PrefillGuardError("expected exactly one pinned Vinardo code cell")
    source = "".join(cells[0].get("source") or [])
    source_sha = _sha256_bytes(source.encode())
    if source_sha != expected_source_sha:
        raise PrefillGuardError("Vinardo notebook cell changed after sequence launch")

    config = yaml.safe_load(config_path.read_text())
    if not isinstance(config, dict):
        raise PrefillGuardError("Vinardo YAML is not a mapping")
    if str(config.get("scoring_function", "")).strip().lower() != "vinardo":
        raise PrefillGuardError("pinned config is not Vinardo scoring")
    if bool(config.get("overwrite_existing") or config.get("overwrite_poses")):
        raise PrefillGuardError("raw prefill refuses overwrite-enabled docking")
    output_dir = Path(str(config.get("output_dir") or ""))
    if not output_dir.is_absolute():
        output_dir = REPO_ROOT / output_dir
    output_dir = output_dir.resolve()

    vina_config_path = Path(str(sequence.get("config") or "")).expanduser().resolve()
    if not vina_config_path.is_file():
        raise PrefillGuardError("pinned Vina config is missing")
    if _sha256_file(vina_config_path) != str(sequence.get("config_sha256") or ""):
        raise PrefillGuardError("Vina YAML changed after sequence launch")
    vina_config = yaml.safe_load(vina_config_path.read_text())
    if not isinstance(vina_config, dict):
        raise PrefillGuardError("Vina YAML is not a mapping")
    vina_output_dir = Path(str(vina_config.get("output_dir") or ""))
    if not vina_output_dir.is_absolute():
        vina_output_dir = REPO_ROOT / vina_output_dir
    refinement_status_path = (
        vina_output_dir.resolve() / "benchmark_gnina_refinement_status.json"
    )
    refinement = _read_json(refinement_status_path)
    try:
        completed_refinement = int(refinement.get("completed_complexes") or 0)
    except (TypeError, ValueError) as exc:
        raise PrefillGuardError("Vina refinement checkpoint count is invalid") from exc
    if refinement.get("status") != "running":
        raise PrefillGuardError("Vina refinement checkpoint is not running")
    if completed_refinement >= STOP_AT_REFINEMENT_COMPLEXES:
        raise PrefillGuardError(
            "too little handoff buffer remains for raw Vinardo prefill "
            f"({completed_refinement}/{STOP_AT_REFINEMENT_COMPLEXES})"
        )

    transformed_source = _transform_cell_source(source)
    return {
        "sequence": sequence,
        "sequence_status_path": sequence_status_path,
        "sequence_pid": sequence_pid,
        "sequence_start_time": sequence_start_time,
        "sequence_process_group_id": expected_pgid,
        "sequence_lock_path": sequence_lock_path,
        "notebook_path": notebook_path,
        "config_path": config_path,
        "output_dir": output_dir,
        "refinement_status_path": refinement_status_path,
        "completed_refinement": completed_refinement,
        "source": source,
        "source_sha256": source_sha,
        "transformed_source": transformed_source,
        "transformed_source_sha256": _sha256_bytes(transformed_source.encode()),
    }


def _write_raw_terminal_status(path: Path, state: str, reason: str) -> None:
    existing: dict[str, Any] = {}
    if path.is_file():
        try:
            existing = json.loads(path.read_text())
        except Exception:
            existing = {}
    _atomic_write_json(path, {
        **existing,
        "status": state,
        "updated_at": _now(),
        "finished_at": _now(),
        "terminal_reason": reason,
    })


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-status", type=Path, default=DEFAULT_SEQUENCE_STATUS)
    parser.add_argument("--lock-file", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--launcher-status", type=Path, default=None)
    parser.add_argument("--log-path", type=Path, default=None,
                        help="Detached stdout/stderr log recorded as metadata.")
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    sequence_status_path = args.sequence_status.expanduser().resolve()
    context = _load_launch_context(sequence_status_path)
    output_dir = context["output_dir"]
    raw_status_path = output_dir / PREFILL_STATUS_NAME
    launcher_status_path = (
        args.launcher_status.expanduser().resolve()
        if args.launcher_status else output_dir / LAUNCHER_STATUS_NAME
    )
    lock_path = args.lock_file.expanduser().resolve()

    plan = {
        "status": "ready" if args.verify_only else "running",
        "stage": PREFILL_STAGE_NAME,
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "process_start_time": _process_start_time(os.getpid()),
        "started_at": _now(),
        "updated_at": _now(),
        "sequence_status": str(sequence_status_path),
        "sequence_pid": context["sequence_pid"],
        "sequence_start_time": context["sequence_start_time"],
        "sequence_process_group_id": context["sequence_process_group_id"],
        "sequence_lock_file": str(context["sequence_lock_path"]),
        "notebook": str(context["notebook_path"]),
        "source_cell_id": VINARDO_CELL_ID,
        "source_cell_sha256": context["source_sha256"],
        "transformed_source_sha256": context["transformed_source_sha256"],
        "vinardo_config": str(context["config_path"]),
        "vinardo_config_sha256": _sha256_file(context["config_path"]),
        "raw_status_path": str(raw_status_path),
        "launcher_status_path": str(launcher_status_path),
        "log_path": str(args.log_path.expanduser().resolve()) if args.log_path else None,
        "lock_file": str(lock_path),
        "optimization": "none",
        "handoff_stage": VINARDO_STAGE_NAME,
        "traversal_order": "descending",
        "refinement_status_path": str(context["refinement_status_path"]),
        "stop_at_refinement_complexes": STOP_AT_REFINEMENT_COMPLEXES,
        "refinement_complexes_at_launch": context["completed_refinement"],
    }
    if args.verify_only:
        print(json.dumps(plan, indent=2, sort_keys=True))
        return 0

    if os.getpgrp() != os.getpid():
        raise PrefillGuardError(
            "raw Vinardo prefill must be launched in its own session/process group "
            "(use setsid)"
        )

    lock_path.parent.mkdir(parents=True, exist_ok=True)
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        lock_handle.close()
        raise PrefillGuardError(f"another raw Vinardo prefill holds {lock_path}") from exc
    lock_handle.seek(0)
    lock_handle.truncate()
    json.dump({"pid": os.getpid(), "started_at": plan["started_at"]}, lock_handle)
    lock_handle.write("\n")
    lock_handle.flush()

    _atomic_write_json(launcher_status_path, plan)

    def commit_terminal(state: str, **details: Any) -> None:
        plan.update({"status": state, "updated_at": _now(), "finished_at": _now(), **details})
        _atomic_write_json(launcher_status_path, plan)

    def signal_handler(signum: int, _frame: Any) -> None:
        raise PrefillSignal(signum)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    def assert_raw_prefill_window() -> None:
        payload = _read_json(sequence_status_path)
        current_pid = payload.get("pid")
        current_stage = payload.get("current_stage")
        current_status = payload.get("status")
        live_start = _process_start_time(context["sequence_pid"])
        if (
            current_status != "running"
            or current_pid != context["sequence_pid"]
            or live_start != context["sequence_start_time"]
            or current_stage != REQUIRED_SEQUENCE_STAGE
        ):
            raise PrefillHandoff(
                "authoritative sequence left the Vina-refinement ownership window "
                f"(status={current_status!r}, stage={current_stage!r}, pid={current_pid!r})"
            )
        refinement = _read_json(context["refinement_status_path"])
        try:
            completed = int(refinement.get("completed_complexes") or 0)
        except (TypeError, ValueError) as exc:
            raise PrefillHandoff("Vina refinement checkpoint became invalid") from exc
        if refinement.get("status") != "running" or completed >= STOP_AT_REFINEMENT_COMPLEXES:
            raise PrefillHandoff(
                "conservative refinement threshold reached "
                f"(status={refinement.get('status')!r}, completed={completed}, "
                f"threshold={STOP_AT_REFINEMENT_COMPLEXES})"
            )

    namespace = {
        "__name__": "__main__",
        "__file__": f"{context['notebook_path']}#{VINARDO_CELL_ID}:raw-prefill",
        "__package__": None,
        "_RAW_PREFILL_CELL_ID": VINARDO_CELL_ID,
        "_RAW_PREFILL_SOURCE_SHA256": context["source_sha256"],
        "_RAW_PREFILL_STOP_THRESHOLD": STOP_AT_REFINEMENT_COMPLEXES,
        "_assert_raw_prefill_window": assert_raw_prefill_window,
    }
    previous_cwd = Path.cwd()
    os.chdir(REPO_ROOT)
    try:
        exec(
            compile(
                context["transformed_source"],
                f"{context['notebook_path']}#{VINARDO_CELL_ID}:raw-prefill",
                "exec",
            ),
            namespace,
            namespace,
        )
    except PrefillHandoff as exc:
        _write_raw_terminal_status(raw_status_path, "handed_off", str(exc))
        commit_terminal("handed_off", handoff_reason=str(exc))
        print(f"RAW VINARDO PREFILL HANDED OFF: {exc}", flush=True)
        return 0
    except PrefillSignal as exc:
        _write_raw_terminal_status(raw_status_path, "interrupted", str(exc))
        commit_terminal("interrupted", signal=exc.signum, error=str(exc))
        print(f"RAW VINARDO PREFILL INTERRUPTED: {exc}", file=sys.stderr, flush=True)
        return 128 + exc.signum
    except BaseException as exc:
        _write_raw_terminal_status(
            raw_status_path, "failed", f"{type(exc).__name__}: {exc}",
        )
        commit_terminal(
            "failed", error_type=type(exc).__name__, error=str(exc),
            traceback=traceback.format_exc(),
        )
        print(f"RAW VINARDO PREFILL FAILED: {type(exc).__name__}: {exc}",
              file=sys.stderr, flush=True)
        traceback.print_exc()
        return 1
    finally:
        os.chdir(previous_cwd)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()

    commit_terminal("completed")
    print("RAW VINARDO PREFILL COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
