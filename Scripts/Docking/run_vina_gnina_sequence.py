#!/usr/bin/env python3
"""Run the benchmark Vina/GNINA notebook stages as one guarded chain.

The launcher deliberately does not daemonize.  A caller may run it in a
terminal or place it under its preferred process supervisor.  It executes the
two audited notebook cells with ordinary ``compile``/``exec`` semantics and a
fresh global namespace for each cell:

1. AutoDock Vina followed by GNINA empirical minimization/CNN rescore.
2. GNINA CNN refinement of the validated raw Vina poses.

The second stage is never entered unless the first returns normally.  Exact
source hashes and lightweight semantic checks make cell-ID drift or an
unreviewed notebook edit a fail-closed event.
"""

from __future__ import annotations

import argparse
import ast
import fcntl
import hashlib
import json
import math
import os
import signal
import sys
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import FrameType
from typing import Any, Mapping, Sequence

import yaml


RESCORE_CELL_ID = "ce11c3f1"
REFINEMENT_CELL_ID = "gnina-refinement-code"


@dataclass(frozen=True)
class CellSpec:
    """The reviewed identity and minimum semantics of one notebook stage."""

    cell_id: str
    label: str
    expected_sha256: str
    config_binding: str | None
    required_fragments: tuple[str, ...] = ()
    forbidden_fragments: tuple[str, ...] = ()


DEFAULT_CELL_SPECS: tuple[CellSpec, ...] = (
    CellSpec(
        cell_id=RESCORE_CELL_ID,
        label="autodock_vina_gnina_rescore",
        expected_sha256="ae2a60c652a1ecc6eb1ff0fcd97443e4b51d48e58b8f08a447994e79f524974b",
        config_binding="config_path",
        required_fragments=(
            '!= "gnina"',
            '!= "rescore"',
            'OPT_TOOLS != ["gnina"]',
            '"benchmark_gnina_rescore_status.json"',
            "run_autodock_vina(",
            "optimize_autodock_results(",
        ),
    ),
    CellSpec(
        cell_id=REFINEMENT_CELL_ID,
        label="autodock_vina_gnina_refinement",
        expected_sha256="5a4736ccb42885814ff246111f4653ce8e27b468ac7ec9236a2744a52aaf59c4",
        config_binding="REFINE_CONFIG_PATH",
        required_fragments=(
            'REFINE_CFG["optimization"] = "gnina_refinement"',
            'REFINE_CFG["gnina_cnn_scoring"] = "refinement"',
            '"benchmark_gnina_refinement_status.json"',
            "optimize_autodock_results(",
        ),
        forbidden_fragments=("run_autodock_vina(",),
    ),
)


class ChainGuardError(RuntimeError):
    """The notebook/configuration no longer matches the reviewed protocol."""


class ChainLockBusy(RuntimeError):
    """Another process currently owns the non-blocking chain lock."""


class CellExecutionError(RuntimeError):
    """A guarded notebook cell raised, so the chain stopped."""


class ChainSignal(BaseException):
    """Raised from a termination-signal handler without being swallowed by cells."""

    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"received {signal.Signals(signum).name}")


class ChainSignalInterrupt(KeyboardInterrupt):
    """SIGINT that remains visible to notebook ``except KeyboardInterrupt`` blocks."""

    def __init__(self, signum: int = signal.SIGINT):
        self.signum = signum
        super().__init__(f"received {signal.Signals(signum).name}")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _resolve_from_root(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path).expanduser()
    if not candidate.is_absolute():
        candidate = repo_root / candidate
    return candidate.resolve()


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Durably replace *path* without exposing a partial JSON document."""

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=False, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        try:
            directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        except (AttributeError, OSError):
            return
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


class ExclusiveChainLock:
    """A process-lifetime, fail-fast ``flock`` with inspectable owner metadata."""

    def __init__(self, path: Path, owner: Mapping[str, Any]):
        self.path = path
        self.owner = dict(owner)
        self._handle: Any | None = None

    def __enter__(self) -> "ExclusiveChainLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            existing = handle.read().strip()
            handle.close()
            suffix = f"; owner metadata: {existing}" if existing else ""
            raise ChainLockBusy(f"chain lock is already held: {self.path}{suffix}") from exc

        self._handle = handle
        handle.seek(0)
        handle.truncate()
        json.dump(self.owner, handle, indent=2, default=str)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        return self

    def set_terminal_status(self, status: str) -> None:
        if self._handle is None:
            return
        terminal = {**self.owner, "status": status, "updated_at": _utc_now()}
        self._handle.seek(0)
        self._handle.truncate()
        json.dump(terminal, self._handle, indent=2, default=str)
        self._handle.write("\n")
        self._handle.flush()
        os.fsync(self._handle.fileno())

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        if self._handle is None:
            return
        try:
            fcntl.flock(self._handle.fileno(), fcntl.LOCK_UN)
        finally:
            self._handle.close()
            self._handle = None


def _path_literal_assignment(source: str, binding: str) -> str | None:
    """Return the string inside ``binding = Path(<string>)``, if present."""

    tree = ast.parse(source)
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if not any(isinstance(target, ast.Name) and target.id == binding for target in targets):
            continue
        value = node.value
        if not (
            isinstance(value, ast.Call)
            and isinstance(value.func, ast.Name)
            and value.func.id == "Path"
            and len(value.args) == 1
            and isinstance(value.args[0], ast.Constant)
            and isinstance(value.args[0].value, str)
        ):
            return None
        return value.args[0].value
    return None


def _load_guarded_cells(
    notebook_path: Path,
    config_path: Path,
    repo_root: Path,
    specs: Sequence[CellSpec],
) -> tuple[list[dict[str, Any]], str]:
    raw_notebook = notebook_path.read_bytes()
    try:
        notebook = json.loads(raw_notebook)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ChainGuardError(f"invalid notebook JSON: {notebook_path}: {exc}") from exc

    cells = notebook.get("cells")
    if not isinstance(cells, list):
        raise ChainGuardError(f"notebook has no cells list: {notebook_path}")

    indexed: dict[str, list[tuple[int, Mapping[str, Any]]]] = {}
    for index, cell in enumerate(cells):
        if isinstance(cell, Mapping) and isinstance(cell.get("id"), str):
            indexed.setdefault(cell["id"], []).append((index, cell))

    selected: list[dict[str, Any]] = []
    previous_index = -1
    for spec in specs:
        matches = indexed.get(spec.cell_id, [])
        if len(matches) != 1:
            raise ChainGuardError(
                f"expected exactly one notebook cell id {spec.cell_id!r}, found {len(matches)}"
            )
        index, cell = matches[0]
        if index <= previous_index:
            raise ChainGuardError("guarded notebook cells are not in the required sequence")
        previous_index = index
        if cell.get("cell_type") != "code":
            raise ChainGuardError(f"guarded cell {spec.cell_id!r} is not a code cell")
        source_field = cell.get("source")
        if isinstance(source_field, list) and all(isinstance(line, str) for line in source_field):
            source = "".join(source_field)
        elif isinstance(source_field, str):
            source = source_field
        else:
            raise ChainGuardError(f"guarded cell {spec.cell_id!r} has invalid source")
        try:
            compile(source, f"{notebook_path}#{spec.cell_id}", "exec")
        except SyntaxError as exc:
            raise ChainGuardError(f"guarded cell {spec.cell_id!r} is not valid Python: {exc}") from exc

        digest = _sha256_bytes(source.encode("utf-8"))
        if not spec.expected_sha256 or digest != spec.expected_sha256:
            raise ChainGuardError(
                f"source hash mismatch for cell {spec.cell_id!r}: "
                f"expected {spec.expected_sha256}, found {digest}; review the edit and "
                "update this launcher's pinned hash before running"
            )
        missing = [fragment for fragment in spec.required_fragments if fragment not in source]
        forbidden = [fragment for fragment in spec.forbidden_fragments if fragment in source]
        if missing or forbidden:
            raise ChainGuardError(
                f"semantic guard failed for cell {spec.cell_id!r}; "
                f"missing={missing}, forbidden={forbidden}"
            )
        if spec.config_binding:
            literal = _path_literal_assignment(source, spec.config_binding)
            if literal is None:
                raise ChainGuardError(
                    f"cell {spec.cell_id!r} does not assign {spec.config_binding} = Path(<literal>)"
                )
            cell_config = _resolve_from_root(literal, repo_root)
            if cell_config != config_path:
                raise ChainGuardError(
                    f"cell {spec.cell_id!r} reads {cell_config}, but launcher was given {config_path}"
                )
        selected.append(
            {
                "spec": spec,
                "notebook_index": index,
                "source": source,
                "source_sha256": digest,
            }
        )
    return selected, _sha256_bytes(raw_notebook)


def _load_config(config_path: Path, repo_root: Path) -> tuple[dict[str, Any], Path]:
    try:
        config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise ChainGuardError(f"cannot load YAML config {config_path}: {exc}") from exc
    if not isinstance(config, dict):
        raise ChainGuardError(f"YAML config must contain a mapping: {config_path}")

    expected = {
        "scoring_function": "vina",
        "optimization": "gnina",
        "gnina_cnn_scoring": "rescore",
    }
    for key, expected_value in expected.items():
        actual = str(config.get(key, "")).strip().lower()
        if actual != expected_value:
            raise ChainGuardError(
                f"config {key!r} must be {expected_value!r}, found {config.get(key)!r}"
            )
    try:
        refinement_timeout = float(config["gnina_refinement_timeout"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ChainGuardError("config gnina_refinement_timeout must be finite and > 0") from exc
    if not math.isfinite(refinement_timeout) or refinement_timeout <= 0:
        raise ChainGuardError("config gnina_refinement_timeout must be finite and > 0")

    output_value = config.get("output_dir")
    if not isinstance(output_value, str) or not output_value.strip():
        raise ChainGuardError("config output_dir must be a non-empty path string")
    output_dir = _resolve_from_root(output_value, repo_root)
    return config, output_dir


def _error_payload(exc: BaseException, *, include_traceback: bool = True) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "type": type(exc).__name__,
        "message": str(exc),
    }
    if include_traceback:
        payload["traceback"] = "".join(
            traceback.format_exception(type(exc), exc, exc.__traceback__)
        )
    signum = getattr(exc, "signum", None)
    if signum is not None:
        payload["signal_number"] = int(signum)
        payload["signal_name"] = signal.Signals(signum).name
    return payload


class _SignalHandlers:
    def __init__(self) -> None:
        self._previous: dict[int, Any] = {}

    @staticmethod
    def _handle(signum: int, frame: FrameType | None) -> None:
        del frame
        if signum == signal.SIGINT:
            raise ChainSignalInterrupt(signum)
        raise ChainSignal(signum)

    def __enter__(self) -> "_SignalHandlers":
        for signum in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            self._previous[signum] = signal.getsignal(signum)
            signal.signal(signum, self._handle)
        return self

    def __exit__(self, exc_type: Any, exc: Any, tb: Any) -> None:
        for signum, previous in self._previous.items():
            signal.signal(signum, previous)


def _new_status(
    *,
    repo_root: Path,
    notebook_path: Path,
    notebook_sha256: str,
    config_path: Path,
    output_dir: Path,
    status_path: Path,
    lock_path: Path,
    log_path: Path | None,
    guarded_cells: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    started = _utc_now()
    return {
        "schema_version": 1,
        "chain": "autodock_vina_gnina_rescore_then_refinement",
        "status": "running",
        "started_at": started,
        "updated_at": started,
        "finished_at": None,
        "pid": os.getpid(),
        "pgid": os.getpgid(0),
        "repo_root": str(repo_root),
        "notebook": {
            "path": str(notebook_path),
            "sha256": notebook_sha256,
        },
        "config": {
            "path": str(config_path),
            "sha256": _sha256_file(config_path),
            "output_dir": str(output_dir),
        },
        "status_file": str(status_path),
        "lock_file": str(lock_path),
        "log_file": str(log_path) if log_path else None,
        "current_cell_id": None,
        "error": None,
        "cells": [
            {
                "position": position,
                "id": item["spec"].cell_id,
                "label": item["spec"].label,
                "notebook_index": item["notebook_index"],
                "source_sha256": item["source_sha256"],
                "expected_source_sha256": item["spec"].expected_sha256,
                "status": "pending",
                "started_at": None,
                "finished_at": None,
                "duration_seconds": None,
                "error": None,
            }
            for position, item in enumerate(guarded_cells, 1)
        ],
    }


def _mark_future_blocked(status: dict[str, Any], current_index: int, reason: str) -> None:
    for cell in status["cells"][current_index + 1 :]:
        if cell["status"] == "pending":
            cell["status"] = "blocked"
            cell["error"] = {"type": "UpstreamStageFailed", "message": reason}


def run_chain(
    *,
    repo_root: Path,
    notebook_path: Path,
    config_path: Path,
    status_path: Path | None = None,
    lock_path: Path | None = None,
    log_path: Path | None = None,
    cell_specs: Sequence[CellSpec] = DEFAULT_CELL_SPECS,
) -> dict[str, Any]:
    """Run both guarded stages and return the final status document.

    Exceptions are re-raised only after an atomic terminal checkpoint has been
    written.  ``cell_specs`` is injectable solely to keep focused tests small;
    the command-line interface always uses ``DEFAULT_CELL_SPECS``.
    """

    root = Path(repo_root).expanduser().resolve()
    notebook = _resolve_from_root(notebook_path, root)
    config = _resolve_from_root(config_path, root)
    _, output_dir = _load_config(config, root)
    output_dir.mkdir(parents=True, exist_ok=True)

    status_file = (
        _resolve_from_root(status_path, root)
        if status_path is not None
        else output_dir / "benchmark_vina_gnina_chain_status.json"
    )
    if not _is_within(status_file, output_dir):
        raise ChainGuardError(
            f"chain status must be under configured output_dir {output_dir}: {status_file}"
        )
    lock_file = (
        _resolve_from_root(lock_path, root)
        if lock_path is not None
        else output_dir / ".benchmark_vina_gnina_chain.lock"
    )
    log_file = _resolve_from_root(log_path, root) if log_path is not None else None

    owner = {
        "schema_version": 1,
        "pid": os.getpid(),
        "pgid": os.getpgid(0),
        "started_at": _utc_now(),
        "status": "running",
        "notebook": str(notebook),
        "config": str(config),
        "status_file": str(status_file),
        "log_file": str(log_file) if log_file else None,
    }

    previous_cwd = Path.cwd()
    os.chdir(root)
    try:
        with ExclusiveChainLock(lock_file, owner) as chain_lock:
            try:
                guarded_cells, notebook_sha256 = _load_guarded_cells(
                    notebook, config, root, cell_specs
                )
            except BaseException as exc:
                failure = {
                    **owner,
                    "chain": "autodock_vina_gnina_rescore_then_refinement",
                    "status": "failed",
                    "updated_at": _utc_now(),
                    "finished_at": _utc_now(),
                    "phase": "guard",
                    "error": _error_payload(exc),
                    "notebook_sha256": _sha256_file(notebook) if notebook.is_file() else None,
                    "config_sha256": _sha256_file(config),
                    "cells": [],
                }
                _atomic_write_json(status_file, failure)
                chain_lock.set_terminal_status("failed")
                raise

            status = _new_status(
                repo_root=root,
                notebook_path=notebook,
                notebook_sha256=notebook_sha256,
                config_path=config,
                output_dir=output_dir,
                status_path=status_file,
                lock_path=lock_file,
                log_path=log_file,
                guarded_cells=guarded_cells,
            )
            _atomic_write_json(status_file, status)
            print(f"Sequential chain status: {status_file}", flush=True)
            print(
                "Execution order: "
                + " -> ".join(item["spec"].label for item in guarded_cells),
                flush=True,
            )

            with _SignalHandlers():
                for index, guarded in enumerate(guarded_cells):
                    cell_status = status["cells"][index]
                    started_monotonic = time.monotonic()
                    cell_status["status"] = "running"
                    cell_status["started_at"] = _utc_now()
                    status["current_cell_id"] = cell_status["id"]
                    status["updated_at"] = _utc_now()
                    _atomic_write_json(status_file, status)
                    print(
                        f"Starting stage {index + 1}/{len(guarded_cells)}: "
                        f"{cell_status['label']} (cell {cell_status['id']})",
                        flush=True,
                    )

                    namespace = {
                        "__name__": "__main__",
                        "__file__": f"{notebook}#{cell_status['id']}",
                        "__builtins__": __builtins__,
                    }
                    try:
                        code = compile(
                            guarded["source"],
                            f"{notebook}#{cell_status['id']}",
                            "exec",
                        )
                        exec(code, namespace, namespace)
                    except ChainSignalInterrupt as exc:
                        cell_status["status"] = "interrupted"
                        cell_status["error"] = _error_payload(exc, include_traceback=False)
                        _mark_future_blocked(status, index, "upstream stage was interrupted")
                        status["status"] = "interrupted"
                        status["error"] = cell_status["error"]
                        raise
                    except ChainSignal as exc:
                        cell_status["status"] = "terminated"
                        cell_status["error"] = _error_payload(exc, include_traceback=False)
                        _mark_future_blocked(status, index, "upstream stage was terminated")
                        status["status"] = "terminated"
                        status["error"] = cell_status["error"]
                        raise
                    except KeyboardInterrupt as exc:
                        cell_status["status"] = "interrupted"
                        cell_status["error"] = _error_payload(exc, include_traceback=False)
                        _mark_future_blocked(status, index, "upstream stage was interrupted")
                        status["status"] = "interrupted"
                        status["error"] = cell_status["error"]
                        raise
                    except BaseException as exc:
                        cell_status["status"] = "failed"
                        cell_status["error"] = _error_payload(exc)
                        _mark_future_blocked(
                            status,
                            index,
                            f"upstream stage {cell_status['id']} failed",
                        )
                        status["status"] = "failed"
                        status["error"] = cell_status["error"]
                        raise CellExecutionError(
                            f"stage {cell_status['label']} failed; refinement was not started"
                            if index == 0
                            else f"stage {cell_status['label']} failed"
                        ) from exc
                    else:
                        cell_status["status"] = "completed"
                        print(f"Completed stage: {cell_status['label']}", flush=True)
                    finally:
                        cell_status["finished_at"] = _utc_now()
                        cell_status["duration_seconds"] = round(
                            time.monotonic() - started_monotonic, 6
                        )
                        status["updated_at"] = _utc_now()
                        if status["status"] != "running":
                            status["finished_at"] = _utc_now()
                        _atomic_write_json(status_file, status)

            status["status"] = "completed"
            status["current_cell_id"] = None
            status["finished_at"] = _utc_now()
            status["updated_at"] = status["finished_at"]
            _atomic_write_json(status_file, status)
            chain_lock.set_terminal_status("completed")
            print("Sequential Vina/GNINA chain completed.", flush=True)
            return status
    except (ChainSignalInterrupt, ChainSignal, KeyboardInterrupt):
        # The current-cell finally block already wrote the terminal checkpoint.
        raise
    finally:
        os.chdir(previous_cwd)


def _default_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _build_parser() -> argparse.ArgumentParser:
    root = _default_repo_root()
    parser = argparse.ArgumentParser(
        description=(
            "Execute the audited Vina+GNINA rescore notebook cell, then execute "
            "the GNINA-refinement cell only if rescore succeeds."
        )
    )
    parser.add_argument(
        "--notebook",
        type=Path,
        default=root / "Master_Docking_AD_Full_Protein.ipynb",
        help="Notebook containing the two pinned code-cell IDs.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=root / "Scripts/Docking/autodock_vina_docking_config_full_protein_search_vina.yaml",
        help="Vina YAML read by both guarded cells.",
    )
    parser.add_argument(
        "--status-file",
        type=Path,
        default=None,
        help="Atomic chain-status JSON (must be inside the YAML output_dir).",
    )
    parser.add_argument(
        "--lock-file",
        type=Path,
        default=None,
        help="Exclusive non-blocking chain lock (defaults inside output_dir).",
    )
    parser.add_argument(
        "--log-file",
        type=Path,
        default=None,
        help=(
            "Optional stdout/stderr log path to record in status metadata. "
            "The launcher does not redirect or daemonize itself."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        run_chain(
            repo_root=_default_repo_root(),
            notebook_path=args.notebook,
            config_path=args.config,
            status_path=args.status_file,
            lock_path=args.lock_file,
            log_path=args.log_file,
        )
    except ChainLockBusy as exc:
        print(f"Refusing duplicate launch: {exc}", file=sys.stderr, flush=True)
        return 75
    except ChainSignalInterrupt:
        print("Sequential chain interrupted by SIGINT.", file=sys.stderr, flush=True)
        return 130
    except ChainSignal as exc:
        print(f"Sequential chain terminated by {signal.Signals(exc.signum).name}.", file=sys.stderr, flush=True)
        return 128 + exc.signum
    except KeyboardInterrupt:
        print("Sequential chain interrupted.", file=sys.stderr, flush=True)
        return 130
    except BaseException:
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
