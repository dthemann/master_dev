#!/usr/bin/env python3
"""Run the Vina/Vinardo GNINA notebook stages as one fail-closed sequence.

The notebook remains the source of truth.  This launcher resolves the cells by
stable cell ID, records the exact source hashes, and executes Vina rescore,
Vina refinement, Vinardo rescore, then Vinardo refinement. Each stage starts
only when its predecessor returns normally. It is intended for a detached
process (for example ``setsid ... &``); it does not daemonize itself.
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
DEFAULT_NOTEBOOK = REPO_ROOT / "Master_Docking_AD_Full_Protein.ipynb"
DEFAULT_CONFIG = (
    REPO_ROOT
    / "Scripts/Docking/autodock_vina_docking_config_full_protein_search_vina.yaml"
)
DEFAULT_VINARDO_CONFIG = (
    REPO_ROOT
    / "Scripts/Docking/autodock_vina_docking_config_full_protein_search_vinardo.yaml"
)
DEFAULT_LOCK = Path("/tmp/master_docking_autodock_gnina_sequence.lock")

STAGE_SPECS = (
    {
        "name": "gnina_rescore",
        "cell_id": "ce11c3f1",
        "required_source": (
            "Variant A requires YAML optimization: gnina",
            "Variant A requires YAML gnina_cnn_scoring: rescore",
            "benchmark_gnina_rescore_status.json",
            "optimize_autodock_results",
        ),
        "forbidden_source": (),
    },
    {
        "name": "gnina_refinement",
        "cell_id": "gnina-refinement-code",
        "required_source": (
            'REFINE_CFG["optimization"] = "gnina_refinement"',
            'REFINE_CFG["gnina_cnn_scoring"] = "refinement"',
            "benchmark_gnina_refinement_status.json",
            "optimize_autodock_results",
        ),
        # Refinement is deliberately optimizer-only.  This guard prevents a
        # later notebook edit from quietly turning stage B into a second docking
        # pass while the detached sequence continues to trust its cell ID.
        "forbidden_source": ("run_autodock_vina(",),
    },
    {
        "name": "vinardo_gnina_rescore",
        "cell_id": "a3b4dd6c",
        "required_source": (
            "autodock_vina_docking_config_full_protein_search_vinardo.yaml",
            "Variant A requires YAML optimization: gnina",
            "Variant A requires YAML gnina_cnn_scoring: rescore",
            "benchmark_gnina_rescore_status.json",
            "optimize_autodock_results",
        ),
        "forbidden_source": (),
    },
    {
        "name": "vinardo_gnina_refinement",
        "cell_id": "vinardo-gnina-refinement-code",
        "required_source": (
            'VIN_REFINE_CFG["optimization"] = "gnina_refinement"',
            'VIN_REFINE_CFG["gnina_cnn_scoring"] = "refinement"',
            "benchmark_gnina_refinement_status.json",
            "optimize_autodock_results",
        ),
        "forbidden_source": ("run_autodock_vina(",),
    },
)


class SequenceInterrupted(RuntimeError):
    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"received signal {signum}")


def _now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    finally:
        Path(tmp_name).unlink(missing_ok=True)


def _load_stage_sources(notebook_path: Path) -> tuple[str, list[dict[str, Any]]]:
    notebook_bytes = notebook_path.read_bytes()
    notebook = json.loads(notebook_bytes)
    cells_by_id: dict[str, dict[str, Any]] = {}
    for cell in notebook.get("cells", []):
        cell_id = str(cell.get("id") or "")
        if cell_id in cells_by_id:
            raise RuntimeError(f"duplicate notebook cell ID: {cell_id}")
        cells_by_id[cell_id] = cell

    stages: list[dict[str, Any]] = []
    for spec in STAGE_SPECS:
        cell = cells_by_id.get(spec["cell_id"])
        if cell is None or cell.get("cell_type") != "code":
            raise RuntimeError(
                f"required code cell {spec['cell_id']!r} ({spec['name']}) is missing"
            )
        source = "".join(cell.get("source", []))
        for marker in spec["required_source"]:
            if marker not in source:
                raise RuntimeError(
                    f"cell {spec['cell_id']} failed semantic guard; missing {marker!r}"
                )
        for marker in spec["forbidden_source"]:
            if marker in source:
                raise RuntimeError(
                    f"cell {spec['cell_id']} failed optimizer-only guard; found {marker!r}"
                )
        # Compile both stages before any expensive work starts.
        code = compile(source, f"{notebook_path}#{spec['cell_id']}", "exec")
        stages.append({
            "name": spec["name"],
            "cell_id": spec["cell_id"],
            "source": source,
            "source_sha256": _sha256_bytes(source.encode()),
            "code": code,
            "status": "pending",
        })
    return _sha256_bytes(notebook_bytes), stages


def _public_stages(stages: list[dict[str, Any]]) -> list[dict[str, Any]]:
    keys = ("name", "cell_id", "source_sha256", "status", "started_at", "finished_at")
    return [{key: stage[key] for key in keys if key in stage} for stage in stages]


def _restore_resume_state(
    status_path: Path,
    stages: list[dict[str, Any]],
    *,
    notebook_path: Path,
    config_path: Path,
    vinardo_config_path: Path,
) -> dict[str, Any]:
    """Validate a prior sequence checkpoint and restore its completed prefix."""
    try:
        previous = json.loads(status_path.read_text())
    except FileNotFoundError as exc:
        raise RuntimeError(f"resume status does not exist: {status_path}") from exc
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"cannot read resume status {status_path}: {exc}") from exc

    if previous.get("schema_version") != 1:
        raise RuntimeError("resume status has an unsupported schema_version")
    if Path(str(previous.get("notebook") or "")).resolve() != notebook_path:
        raise RuntimeError("resume status references a different notebook")

    for key, path in (
        ("config", config_path),
        ("vinardo_config", vinardo_config_path),
    ):
        if Path(str(previous.get(key) or "")).resolve() != path:
            raise RuntimeError(f"resume status references a different {key}")
        recorded_hash = str(previous.get(f"{key}_sha256") or "")
        if not recorded_hash or recorded_hash != _sha256_bytes(path.read_bytes()):
            raise RuntimeError(f"resume status {key} hash no longer matches")

    previous_stages = previous.get("stages")
    if not isinstance(previous_stages, list) or len(previous_stages) != len(stages):
        raise RuntimeError("resume status has an unexpected stage list")

    found_incomplete = False
    for current, recorded in zip(stages, previous_stages, strict=True):
        for key in ("name", "cell_id", "source_sha256"):
            if recorded.get(key) != current.get(key):
                raise RuntimeError(
                    f"resume status stage mismatch for {current['name']}: {key}"
                )
        recorded_state = str(recorded.get("status") or "")
        if recorded_state == "completed":
            if found_incomplete:
                raise RuntimeError("resume status completed stages are not a prefix")
            current["status"] = "completed"
            for key in ("started_at", "finished_at"):
                if recorded.get(key):
                    current[key] = recorded[key]
            continue
        if recorded_state not in {"pending", "running", "interrupted", "failed"}:
            raise RuntimeError(
                f"resume status has invalid state for {current['name']}: {recorded_state!r}"
            )
        found_incomplete = True
        current["status"] = "pending"
        current.pop("started_at", None)
        current.pop("finished_at", None)

    if not found_incomplete:
        raise RuntimeError("sequence checkpoint is already complete; nothing to resume")
    return previous


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--notebook", type=Path, default=DEFAULT_NOTEBOOK)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--vinardo-config", type=Path, default=DEFAULT_VINARDO_CONFIG)
    parser.add_argument("--status-path", type=Path, default=None)
    parser.add_argument("--log-path", type=Path, default=None,
                        help="Log path recorded in status (stdout redirection is external).")
    parser.add_argument("--lock-file", type=Path, default=DEFAULT_LOCK)
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Resume the first incomplete stage from --status-path after verifying "
            "the recorded config and notebook-cell hashes."
        ),
    )
    parser.add_argument("--verify-only", action="store_true",
                        help="Validate configuration/cells and print the planned sequence.")
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    notebook_path = args.notebook.expanduser().resolve()
    config_path = args.config.expanduser().resolve()
    vinardo_config_path = args.vinardo_config.expanduser().resolve()
    cfg = yaml.safe_load(config_path.read_text())
    vinardo_cfg = yaml.safe_load(vinardo_config_path.read_text())
    if str(cfg.get("optimization", "")).strip().lower() != "gnina":
        raise RuntimeError("sequence source YAML must pin optimization: gnina")
    if str(cfg.get("gnina_cnn_scoring", "")).strip().lower() != "rescore":
        raise RuntimeError("sequence source YAML must pin gnina_cnn_scoring: rescore")
    if bool(cfg.get("overwrite_existing") or cfg.get("overwrite_poses")):
        raise RuntimeError("refusing sequence: Vina docking overwrite is enabled")
    if bool(cfg.get("overwrite_optimization", False)):
        raise RuntimeError("refusing sequence: optimizer overwrite is enabled")
    if str(cfg.get("scoring_function", "")).strip().lower() != "vina":
        raise RuntimeError("Vina sequence source YAML must pin scoring_function: vina")
    if str(cfg.get("optimize_rank_by", "")).strip().lower() != "cnn_affinity":
        raise RuntimeError("Vina sequence source YAML must rank by cnn_affinity")

    if str(vinardo_cfg.get("optimization", "")).strip().lower() != "gnina":
        raise RuntimeError("Vinardo sequence source YAML must pin optimization: gnina")
    if str(vinardo_cfg.get("gnina_cnn_scoring", "")).strip().lower() != "rescore":
        raise RuntimeError("Vinardo sequence source YAML must pin gnina_cnn_scoring: rescore")
    if str(vinardo_cfg.get("scoring_function", "")).strip().lower() != "vinardo":
        raise RuntimeError("Vinardo sequence source YAML must pin scoring_function: vinardo")
    if str(vinardo_cfg.get("optimize_rank_by", "")).strip().lower() != "cnn_affinity":
        raise RuntimeError("Vinardo sequence source YAML must rank by cnn_affinity")
    if bool(vinardo_cfg.get("overwrite_existing") or vinardo_cfg.get("overwrite_poses")):
        raise RuntimeError("refusing sequence: Vinardo docking overwrite is enabled")
    if bool(vinardo_cfg.get("overwrite_optimization", False)):
        raise RuntimeError("refusing sequence: Vinardo optimizer overwrite is enabled")

    output_dir = Path(cfg["output_dir"])
    if not output_dir.is_absolute():
        output_dir = REPO_ROOT / output_dir
    status_path = (
        args.status_path.expanduser().resolve()
        if args.status_path
        else output_dir / "benchmark_gnina_sequence_status.json"
    )
    notebook_sha256, stages = _load_stage_sources(notebook_path)

    if args.verify_only:
        previous_status = None
        if args.resume:
            previous_status = _restore_resume_state(
                status_path,
                stages,
                notebook_path=notebook_path,
                config_path=config_path,
                vinardo_config_path=vinardo_config_path,
            )
        print(json.dumps({
            "status": "ready",
            "resume": args.resume,
            "resume_count": (
                int(previous_status.get("resume_count", 0)) + 1
                if previous_status else 0
            ),
            "notebook": str(notebook_path),
            "notebook_sha256": notebook_sha256,
            "config": str(config_path),
            "vinardo_config": str(vinardo_config_path),
            "output_dir": str(output_dir),
            "status_path": str(status_path),
            "stages": _public_stages(stages),
        }, indent=2))
        return 0

    args.lock_file.parent.mkdir(parents=True, exist_ok=True)
    lock_handle = open(args.lock_file, "a+")
    try:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        lock_handle.close()
        raise RuntimeError(
            f"another AutoDock GNINA sequence holds {args.lock_file}"
        ) from exc

    previous_status = None
    if args.resume:
        previous_status = _restore_resume_state(
            status_path,
            stages,
            notebook_path=notebook_path,
            config_path=config_path,
            vinardo_config_path=vinardo_config_path,
        )

    now = _now()
    previous_logs: list[str] = []
    if previous_status:
        previous_logs.extend(
            str(value) for value in previous_status.get("previous_log_paths", []) if value
        )
        if previous_status.get("log_path"):
            previous_logs.append(str(previous_status["log_path"]))
        previous_logs = list(dict.fromkeys(previous_logs))

    status: dict[str, Any] = {
        "schema_version": 1,
        "status": "running",
        "started_at": (
            previous_status.get("started_at", now) if previous_status else now
        ),
        "updated_at": now,
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "notebook": str(notebook_path),
        "notebook_sha256": notebook_sha256,
        "config": str(config_path),
        "config_sha256": _sha256_bytes(config_path.read_bytes()),
        "vinardo_config": str(vinardo_config_path),
        "vinardo_config_sha256": _sha256_bytes(vinardo_config_path.read_bytes()),
        "log_path": str(args.log_path.expanduser().resolve()) if args.log_path else None,
        "resume_count": (
            int(previous_status.get("resume_count", 0)) + 1
            if previous_status else 0
        ),
        "resumed_at": now if previous_status else None,
        "previous_log_paths": previous_logs,
        "current_stage": None,
        "stages": _public_stages(stages),
    }

    def commit() -> None:
        status["updated_at"] = _now()
        status["stages"] = _public_stages(stages)
        _atomic_write_json(status_path, status)

    def signal_handler(signum: int, _frame: Any) -> None:
        raise SequenceInterrupted(signum)

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    commit()
    print(f"Sequential AutoDock/GNINA status: {status_path}", flush=True)
    if previous_status:
        completed = [stage["name"] for stage in stages if stage["status"] == "completed"]
        print(
            f"Resuming after {len(completed)} completed stage(s): "
            f"{', '.join(completed) if completed else 'none'}",
            flush=True,
        )

    old_cwd = Path.cwd()
    os.chdir(REPO_ROOT)
    try:
        for stage in stages:
            if stage["status"] == "completed":
                continue
            status["current_stage"] = stage["name"]
            stage["status"] = "running"
            stage["started_at"] = _now()
            commit()
            print(
                f"[{stage['started_at']}] START {stage['name']} "
                f"(cell {stage['cell_id']}, sha256={stage['source_sha256']})",
                flush=True,
            )
            namespace = {
                "__name__": "__main__",
                "__file__": f"{notebook_path}#{stage['cell_id']}",
                "__package__": None,
            }
            exec(stage["code"], namespace, namespace)
            stage["status"] = "completed"
            stage["finished_at"] = _now()
            commit()
            print(f"[{stage['finished_at']}] COMPLETE {stage['name']}", flush=True)
        status.update({
            "status": "completed",
            "current_stage": None,
            "finished_at": _now(),
        })
        commit()
        print(f"[{status['finished_at']}] SEQUENCE COMPLETE", flush=True)
        return 0
    except SequenceInterrupted as exc:
        for stage in stages:
            if stage["status"] == "running":
                stage["status"] = "interrupted"
                stage["finished_at"] = _now()
        status.update({
            "status": "interrupted",
            "finished_at": _now(),
            "signal": exc.signum,
            "error": str(exc),
        })
        commit()
        print(f"SEQUENCE INTERRUPTED: {exc}", file=sys.stderr, flush=True)
        return 128 + exc.signum
    except BaseException as exc:
        for stage in stages:
            if stage["status"] == "running":
                stage["status"] = "failed"
                stage["finished_at"] = _now()
        status.update({
            "status": "failed",
            "finished_at": _now(),
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        })
        commit()
        print(f"SEQUENCE FAILED: {type(exc).__name__}: {exc}",
              file=sys.stderr, flush=True)
        traceback.print_exc()
        return 1
    finally:
        os.chdir(old_cwd)
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)
        lock_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())
