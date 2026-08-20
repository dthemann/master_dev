#!/usr/bin/env python3
"""Stop the detached raw-Vinardo prefill well before sequence handoff.

This watchdog is intentionally a separate process group.  If the authoritative
Vina refinement reaches the conservative threshold, or the sequence leaves its
refinement stage, it terminates the *entire* prefill process group so an active
Vina child cannot be orphaned.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import signal
import tempfile
import time
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SEQUENCE_STATUS = (
    REPO_ROOT / "Dockings/vina_results_full_protein_vina_scoring"
    / "benchmark_gnina_vina_vinardo_sequence_status.json"
)
DEFAULT_REFINEMENT_STATUS = (
    REPO_ROOT / "Dockings/vina_results_full_protein_vina_scoring"
    / "benchmark_gnina_refinement_status.json"
)
DEFAULT_LAUNCHER_STATUS = (
    REPO_ROOT / "Dockings/vina_results_full_protein_vinardo_scoring"
    / "benchmark_raw_vinardo_prefill_launcher_status.json"
)
DEFAULT_STATUS = (
    REPO_ROOT / "Dockings/vina_results_full_protein_vinardo_scoring"
    / "benchmark_raw_vinardo_prefill_watchdog_status.json"
)
TERMINAL_PREFILL_STATES = {"completed", "failed", "interrupted", "handed_off"}


def _now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise RuntimeError(f"expected JSON object: {path}")
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
    try:
        raw = Path(f"/proc/{pid}/stat").read_text()
        close = raw.rfind(")")
        return raw[close + 2 :].split()[19] if close >= 0 else None
    except (OSError, IndexError, ValueError):
        return None


def _process_cmdline(pid: int) -> list[str]:
    try:
        return [
            value.decode(errors="replace")
            for value in Path(f"/proc/{pid}/cmdline").read_bytes().split(b"\0")
            if value
        ]
    except OSError:
        return []


def _same_target_process(
    pid: int, pgid: int, start_time: str, target_script: str,
) -> bool:
    if pid <= 1 or pgid != pid or _process_start_time(pid) != start_time:
        return False
    try:
        if os.getpgid(pid) != pgid:
            return False
    except ProcessLookupError:
        return False
    return any(
        Path(value).name == target_script
        for value in _process_cmdline(pid)
    )


def _same_prefill_process(pid: int, pgid: int, start_time: str) -> bool:
    """Backward-compatible raw-prefill identity helper used by focused tests."""
    return _same_target_process(
        pid, pgid, start_time, "run_autodock_vinardo_raw_prefill.py",
    )


def _live_group_identities(pgid: int) -> dict[int, str]:
    """Return non-zombie members of *pgid* keyed by immutable start ticks."""
    members: dict[int, str] = {}
    for stat_path in Path("/proc").glob("[0-9]*/stat"):
        try:
            raw = stat_path.read_text()
            close = raw.rfind(")")
            fields = raw[close + 2 :].split()
            # suffix fields: state(3), ppid(4), pgrp(5), ... starttime(22)
            if close < 0 or fields[0] == "Z" or int(fields[2]) != pgid:
                continue
            members[int(stat_path.parent.name)] = fields[19]
        except (OSError, IndexError, ValueError):
            continue
    return members


def _terminate_process_group(pgid: int, term_grace: float) -> bool:
    """TERM then KILL a validated group's surviving original members.

    Tracking every member by PID *and* start ticks ensures that leader exit does
    not hide a surviving Vina child, while PID/PGID reuse cannot redirect KILL.
    Returns whether SIGKILL escalation was required.
    """
    original = _live_group_identities(pgid)
    if not original:
        return False

    def survivors() -> dict[int, str]:
        current = _live_group_identities(pgid)
        return {
            pid: start for pid, start in current.items()
            if original.get(pid) == start
        }

    os.killpg(pgid, signal.SIGTERM)
    deadline = time.monotonic() + term_grace
    while time.monotonic() < deadline and survivors():
        time.sleep(0.25)
    remaining = survivors()
    if not remaining:
        return False
    os.killpg(pgid, signal.SIGKILL)
    kill_deadline = time.monotonic() + 10.0
    while time.monotonic() < kill_deadline and survivors():
        time.sleep(0.1)
    if survivors():
        raise RuntimeError("prefill process-group members survived SIGKILL")
    return True


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--launcher-status", type=Path, default=DEFAULT_LAUNCHER_STATUS)
    parser.add_argument("--sequence-status", type=Path, default=DEFAULT_SEQUENCE_STATUS)
    parser.add_argument("--refinement-status", type=Path, default=DEFAULT_REFINEMENT_STATUS)
    parser.add_argument("--status-path", type=Path, default=DEFAULT_STATUS)
    parser.add_argument("--stop-at", type=int, default=250)
    parser.add_argument("--poll-interval", type=float, default=15.0)
    parser.add_argument("--term-grace", type=float, default=60.0)
    parser.add_argument(
        "--target-script", default="run_autodock_vinardo_raw_prefill.py",
        help="Exact basename of the guarded process-group leader.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    if args.stop_at <= 0 or args.poll_interval <= 0 or args.term_grace <= 0:
        raise ValueError("watchdog numeric controls must be positive")
    launcher_path = args.launcher_status.expanduser().resolve()
    sequence_path = args.sequence_status.expanduser().resolve()
    refinement_path = args.refinement_status.expanduser().resolve()
    status_path = args.status_path.expanduser().resolve()

    # The launcher writes this before the watchdog is started; fail rather than
    # monitoring an inferred or stale process identity.
    launcher = _read_json(launcher_path)
    try:
        pid = int(launcher["pid"])
        pgid = int(launcher["process_group_id"])
        start_time = str(launcher["process_start_time"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("prefill launcher status lacks process identity") from exc
    target_script = str(args.target_script).strip()
    if not target_script or Path(target_script).name != target_script:
        raise ValueError("--target-script must be a plain executable basename")
    pinned_target_script = str(
        launcher.get("target_script") or "run_autodock_vinardo_raw_prefill.py"
    )
    if pinned_target_script != target_script:
        raise RuntimeError("watchdog target script disagrees with launcher metadata")
    if not _same_target_process(pid, pgid, start_time, target_script):
        raise RuntimeError("prefill process identity does not match launcher status")
    if pgid == os.getpgrp():
        raise RuntimeError("watchdog must run in a separate process group")
    try:
        pinned_sequence = Path(str(launcher["sequence_status"])).expanduser().resolve()
        pinned_refinement = Path(str(launcher["refinement_status_path"])).expanduser().resolve()
        pinned_threshold = int(launcher["stop_at_refinement_complexes"])
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError("launcher status lacks watchdog protocol metadata") from exc
    if (
        pinned_sequence != sequence_path
        or pinned_refinement != refinement_path
        or pinned_threshold != args.stop_at
    ):
        raise RuntimeError("watchdog arguments disagree with pinned launcher protocol")

    state: dict[str, Any] = {
        "status": "monitoring",
        "started_at": _now(),
        "updated_at": _now(),
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "prefill_pid": pid,
        "prefill_process_group_id": pgid,
        "prefill_process_start_time": start_time,
        "launcher_status": str(launcher_path),
        "sequence_status": str(sequence_path),
        "refinement_status": str(refinement_path),
        "stop_at_refinement_complexes": args.stop_at,
        "target_script": target_script,
    }
    _atomic_write_json(status_path, state)

    while True:
        launcher = _read_json(launcher_path)
        launcher_state = str(launcher.get("status") or "")
        if launcher_state in TERMINAL_PREFILL_STATES:
            state.update({
                "status": "completed_without_intervention",
                "updated_at": _now(), "finished_at": _now(),
                "prefill_terminal_status": launcher_state,
            })
            _atomic_write_json(status_path, state)
            return 0
        if not _same_target_process(pid, pgid, start_time, target_script):
            state.update({
                "status": "failed", "updated_at": _now(), "finished_at": _now(),
                "error": "prefill disappeared without a terminal launcher status",
            })
            _atomic_write_json(status_path, state)
            return 1

        sequence = _read_json(sequence_path)
        refinement = _read_json(refinement_path)
        try:
            completed = int(refinement.get("completed_complexes") or 0)
        except (TypeError, ValueError):
            completed = args.stop_at
        trigger: str | None = None
        if sequence.get("status") != "running" or sequence.get("current_stage") != "gnina_refinement":
            trigger = (
                "authoritative sequence left gnina_refinement "
                f"(status={sequence.get('status')!r}, stage={sequence.get('current_stage')!r})"
            )
        elif refinement.get("status") != "running" or completed >= args.stop_at:
            trigger = (
                "conservative refinement threshold reached "
                f"(status={refinement.get('status')!r}, completed={completed}, "
                f"threshold={args.stop_at})"
            )

        if trigger is not None:
            # Revalidate immediately before the destructive signal so PID reuse
            # can never redirect a process-group kill.
            if not _same_target_process(pid, pgid, start_time, target_script):
                state.update({
                    "status": "failed", "updated_at": _now(), "finished_at": _now(),
                    "error": "prefill identity changed before watchdog termination",
                })
                _atomic_write_json(status_path, state)
                return 1
            escalated = _terminate_process_group(pgid, args.term_grace)
            state.update({
                "status": "terminated_prefill", "updated_at": _now(),
                "finished_at": _now(), "reason": trigger,
                "sigkill_escalated": escalated,
            })
            _atomic_write_json(status_path, state)
            return 0

        state.update({
            "updated_at": _now(),
            "last_completed_refinement_complexes": completed,
        })
        _atomic_write_json(status_path, state)
        time.sleep(args.poll_interval)


if __name__ == "__main__":
    raise SystemExit(main())
