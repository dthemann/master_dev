#!/usr/bin/env python3
"""Wait for all AutoDock/GNINA stages, then run strict PoseBusters validation.

The supervisor is deliberately separate from ``run_autodock_gnina_sequence``
so it can be attached to an already-running docking sequence. It holds the
sequence lock from artifact preflight through PoseBusters postflight, preventing
a new docking run from changing the inputs mid-validation.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
EXPECTED_STAGE_NAMES = (
    "gnina_rescore",
    "gnina_refinement",
    "vinardo_gnina_rescore",
    "vinardo_gnina_refinement",
)
EXPECTED_ARMS = (
    ("autodock", "gnina"),
    ("autodock", "gnina_refinement"),
    ("autodock_vinardo", "gnina"),
    ("autodock_vinardo", "gnina_refinement"),
)
EXPECTED_COMPLEXES = 308
EXPECTED_PROTOCOLS = {
    ("autodock", "gnina"): ("vina", "rescore"),
    ("autodock", "gnina_refinement"): ("vina", "refinement"),
    ("autodock_vinardo", "gnina"): ("vinardo", "rescore"),
    ("autodock_vinardo", "gnina_refinement"): ("vinardo", "refinement"),
}

DEFAULT_SEQUENCE_STATUS = (
    REPO_ROOT / "Dockings/vina_results_full_protein_vina_scoring/"
    "benchmark_gnina_vina_vinardo_sequence_status.json"
)
DEFAULT_CONFIG = (
    REPO_ROOT / "Scripts/Docking/Posebusters/"
    "posebusters_autodock_gnina_variants_config.yaml"
)
DEFAULT_RUNNER = REPO_ROOT / "Scripts/Docking/Posebusters/run_posebusters.py"
DEFAULT_PYTHON = Path("/home/manndo/anaconda3/envs/vina/bin/python")
DEFAULT_OUTPUT_BASE = (
    REPO_ROOT / "posebusters_results/"
    "benchmark_full_protein_autodock_gnina_variants"
)
DEFAULT_STATUS = DEFAULT_OUTPUT_BASE / "posebusters_after_autodock_status.json"
DEFAULT_MANIFEST = DEFAULT_OUTPUT_BASE / "posebusters_input_manifest.json"
DEFAULT_SEQUENCE_LOCK = Path("/tmp/master_docking_autodock_gnina_sequence.lock")
DEFAULT_WATCHER_LOCK = Path("/tmp/master_docking_posebusters_after_autodock.lock")
DEFAULT_POSEBUSTERS_LOCK = Path(
    "/tmp/master_docking_posebusters_autodock_gnina_variants.lock"
)

CHECKPOINTS = {
    "gnina_rescore": (
        REPO_ROOT / "Dockings/vina_results_full_protein_vina_scoring/"
        "benchmark_gnina_rescore_status.json",
        {"stage": "gnina_rescore", "optimization": "gnina", "cnn_scoring": "rescore"},
    ),
    "gnina_refinement": (
        REPO_ROOT / "Dockings/vina_results_full_protein_vina_scoring/"
        "benchmark_gnina_refinement_status.json",
        {"stage": "gnina_refinement", "optimization": "gnina_refinement",
         "cnn_scoring": "refinement"},
    ),
    "vinardo_gnina_rescore": (
        REPO_ROOT / "Dockings/vina_results_full_protein_vinardo_scoring/"
        "benchmark_gnina_rescore_status.json",
        {"stage": "gnina_rescore", "source_scoring": "vinardo",
         "optimization": "gnina", "cnn_scoring": "rescore"},
    ),
    "vinardo_gnina_refinement": (
        REPO_ROOT / "Dockings/vina_results_full_protein_vinardo_scoring/"
        "benchmark_gnina_refinement_status.json",
        {"stage": "gnina_refinement", "source_scoring": "vinardo",
         "optimization": "gnina_refinement", "cnn_scoring": "refinement"},
    ),
}


class SupervisorInterrupted(RuntimeError):
    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"received signal {signum}")


def _now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(value: Any) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return value


def _atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence-status", type=Path, default=DEFAULT_SEQUENCE_STATUS)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--python", type=Path, default=DEFAULT_PYTHON)
    parser.add_argument("--status-path", type=Path, default=DEFAULT_STATUS)
    parser.add_argument("--manifest-path", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--log-path", type=Path, default=None)
    parser.add_argument("--sequence-lock", type=Path, default=DEFAULT_SEQUENCE_LOCK)
    parser.add_argument("--watcher-lock", type=Path, default=DEFAULT_WATCHER_LOCK)
    parser.add_argument("--posebusters-lock", type=Path, default=DEFAULT_POSEBUSTERS_LOCK)
    parser.add_argument("--poll-interval", type=int, default=60)
    parser.add_argument("--monitor-interval", type=int, default=30)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--max-attempts", type=int, default=4)
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def _resolve(path: Path) -> Path:
    path = path.expanduser()
    return path.resolve() if path.is_absolute() else (REPO_ROOT / path).resolve()


def _validate_static_config(config_path: Path) -> dict[str, Any]:
    raw = yaml.safe_load(config_path.read_text()) or {}
    expected_dirs = {
        "autodock": "Dockings/vina_results_full_protein_vina_scoring",
        "autodock_vinardo": "Dockings/vina_results_full_protein_vinardo_scoring",
    }
    errors: list[str] = []
    if raw.get("docking_directories") != expected_dirs:
        errors.append("docking_directories must contain only the Vina and Vinardo roots")
    if raw.get("validation_receptor_directories") != expected_dirs:
        errors.append("validation_receptor_directories must match the two docking roots")
    if raw.get("config_mode") != "dock" or raw.get("receptor_policy") != "prepared_per_pose":
        errors.append("PoseBusters must use dock/prepared_per_pose mode")
    if raw.get("require_validation_receptors") is not True:
        errors.append("require_validation_receptors must be true")
    if raw.get("require_all_methods") is not True:
        errors.append("require_all_methods must be true")
    if raw.get("require_complete_results") is not True:
        errors.append("require_complete_results must be true")
    if raw.get("optional_methods") not in ([], None):
        errors.append("optional_methods must be empty")
    if raw.get("expected_common_combinations") != EXPECTED_COMPLEXES:
        errors.append(f"expected_common_combinations must be {EXPECTED_COMPLEXES}")
    variant_filter = raw.get("variant_filter") or {}
    variants = variant_filter.get("autodock") or []
    if (set(variant_filter) != {"autodock"}
            or {str(value).strip().lower() for value in variants} != {
                "gnina", "gnina_refinement",
            }):
        errors.append(
            "variant_filter must contain only autodock=gnina/gnina_refinement"
        )
    if raw.get("overwrite") is not False:
        errors.append("overwrite must be false for content-aware resume")
    try:
        configured_output = _resolve(Path(str(raw.get("output_base_dir"))))
    except (TypeError, ValueError):
        configured_output = None
    if configured_output != DEFAULT_OUTPUT_BASE.resolve():
        errors.append(
            f"output_base_dir must be the dedicated directory {DEFAULT_OUTPUT_BASE}"
        )
    pose_timeout = raw.get("pose_timeout")
    if not isinstance(pose_timeout, int) or pose_timeout <= 0:
        errors.append("pose_timeout must be a finite positive integer")
    if raw.get("monitor") is not True:
        errors.append("PoseBusters internal monitoring must be enabled")
    for key in ("stall_timeout", "max_restarts"):
        value = raw.get(key)
        if not isinstance(value, int) or value <= 0:
            errors.append(f"{key} must be a positive integer")
    if errors:
        raise RuntimeError("strict PoseBusters config failed validation: " + "; ".join(errors))
    return raw


def _probe_runtime(python: Path, runner: Path) -> dict[str, str]:
    if not python.is_file() or not os.access(python, os.X_OK):
        raise FileNotFoundError(f"PoseBusters Python is unavailable: {python}")
    if not runner.is_file():
        raise FileNotFoundError(runner)
    probe = subprocess.run(
        [str(python), "-c", (
            "import importlib.metadata,sys; import posebusters; "
            "print(sys.version.split()[0]); "
            "print(importlib.metadata.version('posebusters'))"
        )],
        cwd=REPO_ROOT, text=True, capture_output=True, timeout=30,
    )
    if probe.returncode:
        raise RuntimeError(f"PoseBusters runtime probe failed: {probe.stderr.strip()}")
    values = probe.stdout.splitlines()
    if len(values) < 2:
        raise RuntimeError(f"unexpected PoseBusters runtime probe output: {probe.stdout!r}")
    return {"python": values[0].strip(), "posebusters": values[1].strip()}


def _launch_posebusters(
    command: list[str],
    sequence_lock_handle: Any,
    posebusters_lock_handle: Any,
) -> subprocess.Popen:
    """Launch in a new group while retaining both campaign lock descriptors."""
    return subprocess.Popen(
        command,
        cwd=REPO_ROOT,
        start_new_session=True,
        stdin=subprocess.DEVNULL,
        pass_fds=(
            posebusters_lock_handle.fileno(),
            sequence_lock_handle.fileno(),
        ),
    )


def _stop_process_group(process: subprocess.Popen, timeout_s: int = 30) -> None:
    """Terminate and reap a child group; escalate only after a grace period."""
    if process.poll() is None:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        try:
            process.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait(timeout=timeout_s)
    else:
        process.wait()


def _validate_sequence_completion(payload: dict[str, Any]) -> dict[str, Any]:
    if payload.get("status") != "completed" or payload.get("current_stage") is not None:
        raise RuntimeError("sequence has not reached a clean completed state")
    stages = payload.get("stages") or []
    names = tuple(stage.get("name") for stage in stages)
    if names != EXPECTED_STAGE_NAMES:
        raise RuntimeError(f"unexpected sequence stage order: {names}")
    if any(stage.get("status") != "completed" for stage in stages):
        raise RuntimeError("not every sequence stage is completed")

    # The notebook can legitimately gain unrelated documentation or analysis
    # edits while this long sequence runs. Pin the four executed producer cells
    # below; pin the two complete scientific YAMLs here.
    for path_key, hash_key in (
        ("config", "config_sha256"),
        ("vinardo_config", "vinardo_config_sha256"),
    ):
        path = Path(str(payload.get(path_key) or ""))
        expected_hash = str(payload.get(hash_key) or "")
        if not path.is_file() or not expected_hash or _sha256(path) != expected_hash:
            raise RuntimeError(f"sequence source changed or is missing: {path_key}")

    notebook_path = Path(str(payload.get("notebook") or ""))
    if not notebook_path.is_file():
        raise RuntimeError("sequence notebook is missing")
    notebook = _read_json(notebook_path)
    cells = {str(cell.get("id") or ""): "".join(cell.get("source") or [])
             for cell in notebook.get("cells") or []}
    for stage in stages:
        source = cells.get(str(stage.get("cell_id") or ""))
        if source is None:
            raise RuntimeError(f"completed sequence cell is missing: {stage.get('cell_id')}")
        digest = hashlib.sha256(source.encode()).hexdigest()
        if digest != stage.get("source_sha256"):
            raise RuntimeError(f"completed sequence cell hash drifted: {stage.get('name')}")

    checkpoint_report: dict[str, Any] = {}
    for name in EXPECTED_STAGE_NAMES:
        path, expected_fields = CHECKPOINTS[name]
        checkpoint = _read_json(path)
        if checkpoint.get("status") != "completed":
            raise RuntimeError(f"stage checkpoint is not completed: {path}")
        for key, expected in expected_fields.items():
            if checkpoint.get(key) != expected:
                raise RuntimeError(
                    f"stage checkpoint protocol mismatch for {name}: {key}"
                )
        count_field = (
            "processed_complexes"
            if "processed_complexes" in checkpoint
            else "completed_complexes"
        )
        if checkpoint.get(count_field) != EXPECTED_COMPLEXES:
            raise RuntimeError(f"stage checkpoint does not cover 308 complexes: {name}")
        for key in ("issues", "optimizer_issues"):
            if checkpoint.get(key):
                raise RuntimeError(f"stage checkpoint records {key}: {name}")
        totals = checkpoint.get("optimizer_totals") or {}
        if int(totals.get("failed") or 0):
            raise RuntimeError(f"stage checkpoint records optimizer failures: {name}")
        for key in ("missing_inputs", "failed_preparation", "failed_docking"):
            if int(checkpoint.get(key) or 0):
                raise RuntimeError(f"stage checkpoint records {key}: {name}")
        checkpoint_report[name] = {
            "path": str(path), "sha256": _sha256(path),
            "updated_at": checkpoint.get("updated_at"),
        }
    return checkpoint_report


def _benchmark_ids() -> set[str]:
    ids_path = REPO_ROOT / "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
    values = {
        line.strip() for line in ids_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    if len(values) != EXPECTED_COMPLEXES:
        raise RuntimeError(
            f"expected {EXPECTED_COMPLEXES} authoritative benchmark IDs, found {len(values)}"
        )
    return values


def _pose_complex_id(pose: dict[str, Any], root: Path, expected_ids: set[str]) -> str:
    for field in ("pose_file", "source_pose_file"):
        value = pose.get(field)
        if not value:
            continue
        try:
            relative = Path(str(value)).resolve().relative_to(root.resolve())
        except (OSError, ValueError):
            continue
        if relative.parts and relative.parts[0] in expected_ids:
            return relative.parts[0]
    raise RuntimeError(f"pose is outside its authoritative result root: {pose.get('pose_file')}")


_IDENTITY_FIELDS = (
    "docking_method", "optimizer", "protein", "ligand",
    "pose_file", "pose_name", "file_format",
    "protein_file_used", "receptor_file_used", "posebusters_mode",
    "source_pose_file", "receptor_source_file", "receptor_stage",
    "pose_sha256", "source_pose_sha256", "receptor_sha256",
    "receptor_source_sha256", "posebusters_runtime_signature",
    "autodock_rank", "autodock_affinity", "optimized_rank", "rank_metric",
    "minimized_affinity", "cnn_score", "cnn_affinity",
    "optimizer_elapsed_time_s", "optimizer_processing_elapsed_time_s",
    "optimization_log_file", "optimizer_provenance_file",
    "optimizer_provenance_fingerprint", "source_pose_model_sha256",
    "source_scoring", "optimizer_scoring", "optimizer_search",
    "cnn_scoring", "cnn_model",
)


_INTEGER_IDENTITY_FIELDS = {"autodock_rank", "optimized_rank"}
_FLOAT_IDENTITY_FIELDS = {
    "autodock_affinity", "minimized_affinity", "cnn_score", "cnn_affinity",
    "optimizer_elapsed_time_s", "optimizer_processing_elapsed_time_s",
}


def _identity_value(field: str, value: Any) -> str | int | float:
    if value is None or (isinstance(value, str) and not value.strip()):
        return ""
    if field in _INTEGER_IDENTITY_FIELDS:
        number = float(value)
        if not math.isfinite(number) or not number.is_integer() or number < 1:
            raise ValueError(f"invalid positive integer identity value for {field}: {value!r}")
        return int(number)
    if field in _FLOAT_IDENTITY_FIELDS:
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"invalid finite numeric identity value for {field}: {value!r}")
        return number
    return str(value).strip()


def _identity_mismatches(
    actual: dict[str, Any], expected: dict[str, Any],
) -> list[str]:
    mismatched: list[str] = []
    for field in _IDENTITY_FIELDS:
        actual_value = actual[field]
        expected_value = expected[field]
        if field in _FLOAT_IDENTITY_FIELDS:
            if not math.isclose(
                float(actual_value), float(expected_value), rel_tol=1e-12, abs_tol=1e-12,
            ):
                mismatched.append(field)
        elif actual_value != expected_value:
            mismatched.append(field)
    return mismatched


def _pose_identity(pose: dict[str, Any]) -> dict[str, str | int | float]:
    source = {
        "docking_method": pose.get("method"),
        "protein_file_used": pose.get("receptor_file"),
        "receptor_file_used": pose.get("receptor_file"),
        "posebusters_mode": "dock",
        **{field: pose.get(field) for field in _IDENTITY_FIELDS if field != "docking_method"},
    }
    # Derived worker-output fields above must not be overwritten by absent pose
    # keys from the broad metadata copy.
    source["protein_file_used"] = pose.get("receptor_file")
    source["receptor_file_used"] = pose.get("receptor_file")
    source["posebusters_mode"] = "dock"
    identity = {
        field: _identity_value(field, source.get(field)) for field in _IDENTITY_FIELDS
    }
    missing = [field for field in _IDENTITY_FIELDS if identity[field] == ""]
    if missing:
        raise RuntimeError(
            "strict AutoDock pose identity is missing: " + ", ".join(missing)
        )
    return identity


def _validate_pose_scope(
    pose: dict[str, Any],
    arm: tuple[str, str],
    root: Path,
    complex_id: str,
) -> None:
    """Require the canonical raw container and exact four-arm protocol."""
    source_scoring, cnn_scoring = EXPECTED_PROTOCOLS[arm]
    pair = (str(pose.get("protein") or ""), str(pose.get("ligand") or ""))
    if pair != (complex_id, complex_id):
        raise RuntimeError(
            f"{arm} normalized pair {pair} does not match benchmark ID {complex_id}"
        )

    expected_source = (
        root.resolve() / complex_id / "mgl_tools" / "docking" /
        f"{complex_id}_protein_mgl_tools__{complex_id}_ligand_start_conf_"
        f"{source_scoring}_vina_out.pdbqt"
    )
    source = Path(str(pose.get("source_pose_file") or "")).resolve()
    if source != expected_source or not source.is_file():
        raise RuntimeError(
            f"{arm} pose is not tied to its canonical raw container: {source}"
        )

    rank = _identity_value("autodock_rank", pose.get("autodock_rank"))
    optimizer = arm[1]
    expected_pose = (
        expected_source.parent / f"optimized_{optimizer}" /
        f"{expected_source.stem}_rank{rank}_{optimizer}.sdf"
    )
    pose_file = Path(str(pose.get("pose_file") or "")).resolve()
    if pose_file != expected_pose or not pose_file.is_file():
        raise RuntimeError(f"{arm} optimized pose path is noncanonical: {pose_file}")
    expected_log = expected_source.parent.parent / "optimization_log.csv"
    if Path(str(pose.get("optimization_log_file") or "")).resolve() != expected_log:
        raise RuntimeError(f"{arm} optimization log path is noncanonical")
    expected_sidecar = Path(f"{expected_pose}.provenance.json")
    if (Path(str(pose.get("optimizer_provenance_file") or "")).resolve()
            != expected_sidecar or not expected_sidecar.is_file()):
        raise RuntimeError(f"{arm} optimizer provenance path is noncanonical")

    expected_protocol = {
        "source_scoring": source_scoring,
        "optimizer_scoring": "default",
        "optimizer_search": "minimize",
        "cnn_scoring": cnn_scoring,
        "cnn_model": "crossdock_default2018_ensemble",
        "rank_metric": "cnn_affinity",
    }
    mismatches = {
        field: (pose.get(field), expected)
        for field, expected in expected_protocol.items()
        if str(pose.get(field) or "").strip().lower() != expected
    }
    if mismatches:
        raise RuntimeError(f"{arm} optimizer protocol mismatch: {mismatches}")


def _build_input_manifest(
    config_path: Path, manifest_path: Path,
) -> tuple[dict, dict, list[str]]:
    from Scripts.Docking.Posebusters import run_posebusters as pb

    pb._apply_quiet_mode(True)
    ctx = pb.load_config(config_path)
    pose_rows, _ = pb.collect_pose_rows(ctx)
    pose_rows = pb.filter_common_combos(pose_rows, ctx.optional_methods)
    actual_common = len(pose_rows[["protein", "ligand"]].drop_duplicates())
    if actual_common != EXPECTED_COMPLEXES:
        raise RuntimeError(
            f"strict four-arm cohort has {actual_common}/{EXPECTED_COMPLEXES} complexes"
        )
    poses = pb.collect_all_pose_files(pose_rows, ctx)
    poses = pb.attach_validation_receptors(poses, ctx, {})
    verdict_columns = list(pb.expected_test_columns("dock"))
    if not verdict_columns:
        raise RuntimeError("PoseBusters exposes no canonical dock-mode verdict schema")

    expected_ids = _benchmark_ids()
    grouped: dict[tuple[str, str], dict[str, Any]] = {
        arm: {"records": {}, "pairs": set(), "complex_counts": {}}
        for arm in EXPECTED_ARMS
    }
    for pose in poses:
        arm = (str(pose.get("method") or ""), str(pose.get("optimizer") or ""))
        if arm not in grouped:
            raise RuntimeError(f"out-of-scope pose reached strict manifest: {arm}")
        signature = str(pose.get("validation_signature") or "")
        if not signature:
            raise RuntimeError("PoseBusters pose lacks a validation signature")
        group = grouped[arm]
        if signature in group["records"]:
            raise RuntimeError(f"duplicate input validation signature in {arm}: {signature}")
        pair = (str(pose.get("protein") or ""), str(pose.get("ligand") or ""))
        complex_id = _pose_complex_id(
            pose, ctx.docking_directories[arm[0]], expected_ids,
        )
        _validate_pose_scope(
            pose, arm, ctx.docking_directories[arm[0]], complex_id,
        )
        group["records"][signature] = _pose_identity(pose)
        group["pairs"].add(pair)
        group["complex_counts"][complex_id] = group["complex_counts"].get(complex_id, 0) + 1

    report: dict[str, Any] = {}
    serializable_arms: dict[str, Any] = {}
    for arm in EXPECTED_ARMS:
        key = ":".join(arm)
        group = grouped[arm]
        complex_ids = set(group["complex_counts"])
        if complex_ids != expected_ids or len(group["pairs"]) != EXPECTED_COMPLEXES:
            missing = sorted(expected_ids - complex_ids)
            extra = sorted(complex_ids - expected_ids)
            raise RuntimeError(
                f"{key} coverage mismatch: complexes={len(complex_ids)}, "
                f"pairs={len(group['pairs'])}, missing={missing[:10]}, extra={extra[:10]}"
            )
        signatures = sorted(group["records"])
        report[key] = {
            "complexes": len(complex_ids),
            "poses": len(signatures),
            "signature_sha256": _json_sha256(signatures),
        }
        serializable_arms[key] = {
            **report[key],
            "pairs": sorted([list(pair) for pair in group["pairs"]]),
            "complex_pose_counts": dict(sorted(group["complex_counts"].items())),
            "signature_identities": {
                signature: group["records"][signature] for signature in signatures
            },
        }

    for method in ("autodock", "autodock_vinardo"):
        rescore = grouped[(method, "gnina")]["complex_counts"]
        refinement = grouped[(method, "gnina_refinement")]["complex_counts"]
        if rescore != refinement:
            raise RuntimeError(
                f"{method} rescore/refinement rank coverage differs by complex"
            )

    manifest = {
        "schema_version": 2,
        "generated_at": _now(),
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "expected_complexes": EXPECTED_COMPLEXES,
        "verdict_columns": verdict_columns,
        "arms": serializable_arms,
    }
    manifest["manifest_sha256"] = _json_sha256(manifest)
    _atomic_write_json(manifest_path, manifest)
    return grouped, report, verdict_columns


def _verify_output(
    results_csv: Path,
    incomplete_file: Path,
    expected: dict[tuple[str, str], dict[str, Any]],
    verdict_columns: list[str],
) -> dict[str, Any]:
    if incomplete_file.exists():
        raise RuntimeError(f"PoseBusters reported incomplete work: {incomplete_file}")
    if not results_csv.is_file():
        raise FileNotFoundError(results_csv)

    actual = {
        arm: {"records": {}, "pairs": set(), "rows": 0}
        for arm in EXPECTED_ARMS
    }
    all_signatures: set[str] = set()
    with results_csv.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = set(_IDENTITY_FIELDS) | {"validation_signature"} | set(verdict_columns)
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            missing = sorted(required - set(reader.fieldnames or []))
            raise RuntimeError(
                "PoseBusters output is missing required provenance/verdict columns: "
                + ", ".join(missing)
            )
        for row in reader:
            arm = (str(row.get("docking_method") or ""), str(row.get("optimizer") or ""))
            if arm not in actual:
                raise RuntimeError(f"PoseBusters output contains an out-of-scope arm: {arm}")
            signature = str(row.get("validation_signature") or "")
            if not signature or signature in all_signatures:
                raise RuntimeError(f"missing/duplicate output validation signature: {signature}")
            expected_identity = expected[arm]["records"].get(signature)
            if expected_identity is None:
                raise RuntimeError(f"unexpected output validation signature in {arm}: {signature}")
            actual_identity = {
                field: _identity_value(field, row.get(field))
                for field in _IDENTITY_FIELDS
            }
            mismatched = _identity_mismatches(actual_identity, expected_identity)
            if mismatched:
                raise RuntimeError(
                    f"output identity mismatch for {signature}: {', '.join(mismatched)}"
                )
            for column in verdict_columns:
                value = str(row.get(column) or "").strip().lower()
                if value not in {"true", "false", "1", "0", "1.0", "0.0"}:
                    raise RuntimeError(
                        f"non-boolean/blank PoseBusters verdict {column!r} for {signature}"
                    )
            all_signatures.add(signature)
            actual[arm]["records"][signature] = actual_identity
            actual[arm]["pairs"].add((str(row.get("protein") or ""),
                                      str(row.get("ligand") or "")))
            actual[arm]["rows"] += 1

    report: dict[str, Any] = {}
    for arm in EXPECTED_ARMS:
        expected_signatures = set(expected[arm]["records"])
        actual_signatures = set(actual[arm]["records"])
        if actual_signatures != expected_signatures:
            missing = len(expected_signatures - actual_signatures)
            extra = len(actual_signatures - expected_signatures)
            raise RuntimeError(f"{arm} output signature mismatch: missing={missing}, extra={extra}")
        if actual[arm]["pairs"] != expected[arm]["pairs"]:
            raise RuntimeError(f"{arm} output complex-pair coverage mismatch")
        report[":".join(arm)] = {
            "rows": actual[arm]["rows"],
            "complexes": len(actual[arm]["pairs"]),
            "signature_sha256": _json_sha256(sorted(actual_signatures)),
        }
    return report


def _progress_snapshot(sequence: dict[str, Any]) -> dict[str, Any]:
    stage = str(sequence.get("current_stage") or "")
    checkpoint_path = CHECKPOINTS.get(stage, (None, None))[0]
    checkpoint = {}
    if checkpoint_path and checkpoint_path.is_file():
        try:
            checkpoint = _read_json(checkpoint_path)
        except (OSError, ValueError, json.JSONDecodeError):
            checkpoint = {}
    sequence_pid = sequence.get("pid")
    pid_alive = bool(
        isinstance(sequence_pid, int)
        and Path(f"/proc/{sequence_pid}/cmdline").is_file()
    )
    return {
        "sequence_state": sequence.get("status"),
        "sequence_stage": sequence.get("current_stage"),
        "sequence_pid": sequence_pid,
        "sequence_pid_alive": pid_alive,
        "stage_processed_complexes": checkpoint.get(
            "processed_complexes", checkpoint.get("completed_complexes"),
        ),
        "stage_last_complex": checkpoint.get("last_complex"),
        "stage_updated_at": checkpoint.get("updated_at"),
    }


def main() -> int:
    args = _parse_args()
    sequence_status = _resolve(args.sequence_status)
    config_path = _resolve(args.config)
    runner = _resolve(args.runner)
    python = _resolve(args.python)
    status_path = _resolve(args.status_path)
    manifest_path = _resolve(args.manifest_path)
    log_path = _resolve(args.log_path) if args.log_path else None
    if args.poll_interval < 10 or args.monitor_interval < 10:
        raise ValueError("poll and monitor intervals must be at least 10 seconds")
    if args.workers < 1 or args.max_attempts < 1:
        raise ValueError("workers and max-attempts must be positive")

    config = _validate_static_config(config_path)
    runtime = _probe_runtime(python, runner)
    pinned_config_sha256 = _sha256(config_path)
    pinned_runner_sha256 = _sha256(runner)
    plan = {
        "sequence_status_path": str(sequence_status),
        "config": str(config_path),
        "config_sha256": pinned_config_sha256,
        "runner": str(runner),
        "runner_sha256": pinned_runner_sha256,
        "python": str(python),
        "runtime": runtime,
        "status_path": str(status_path),
        "manifest_path": str(manifest_path),
        "output_base_dir": str((REPO_ROOT / config["output_base_dir"]).resolve()),
        "expected_stages": list(EXPECTED_STAGE_NAMES),
        "expected_arms": [list(arm) for arm in EXPECTED_ARMS],
        "sequence_lock": str(args.sequence_lock),
        "posebusters_lock": str(args.posebusters_lock),
    }
    if args.verify_only:
        print(json.dumps({"status": "ready", **plan}, indent=2))
        return 0

    args.watcher_lock.parent.mkdir(parents=True, exist_ok=True)
    watcher_lock = open(args.watcher_lock, "a+")
    try:
        fcntl.flock(watcher_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        watcher_lock.close()
        raise RuntimeError(f"another PoseBusters supervisor holds {args.watcher_lock}") from exc

    status: dict[str, Any] = {
        "schema_version": 2,
        "status": "starting",
        "pid": os.getpid(),
        "process_group_id": os.getpgrp(),
        "started_at": _now(),
        "poll_interval_s": args.poll_interval,
        "monitor_interval_s": args.monitor_interval,
        "max_attempts": args.max_attempts,
        "log_path": str(log_path) if log_path else None,
        **plan,
    }

    def commit(state: str, **details: Any) -> None:
        status.update(details)
        status["status"] = state
        status["updated_at"] = _now()
        _atomic_write_json(status_path, status)

    def revalidate_pinned_sources() -> dict[str, Any]:
        if _sha256(config_path) != pinned_config_sha256:
            raise RuntimeError("strict PoseBusters config changed after supervisor startup")
        if _sha256(runner) != pinned_runner_sha256:
            raise RuntimeError("PoseBusters runner changed after supervisor startup")
        return _validate_static_config(config_path)

    child: subprocess.Popen | None = None

    def signal_handler(signum: int, _frame: Any) -> None:
        raise SupervisorInterrupted(signum)

    def stop_child() -> None:
        """Stop and reap PoseBusters before releasing either inherited lock."""
        nonlocal child
        process = child
        if process is None:
            return
        _stop_process_group(process)
        child = None

    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)

    sequence_lock_handle = None
    posebusters_lock_handle = None
    try:
        commit("waiting_sequence")
        while sequence_lock_handle is None:
            try:
                sequence = _read_json(sequence_status)
            except (OSError, ValueError, json.JSONDecodeError) as exc:
                commit("waiting_sequence", waiting_reason=str(exc))
                time.sleep(args.poll_interval)
                continue
            snapshot = _progress_snapshot(sequence)
            if sequence.get("status") != "completed":
                commit("waiting_sequence", **snapshot)
                time.sleep(args.poll_interval)
                continue

            args.sequence_lock.parent.mkdir(parents=True, exist_ok=True)
            candidate = open(args.sequence_lock, "a+")
            try:
                fcntl.flock(candidate.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                candidate.close()
                commit("waiting_sequence_lock", **snapshot)
                time.sleep(args.poll_interval)
                continue
            sequence_lock_handle = candidate

            # Re-read under the producer lock to close the completion/preflight race.
            sequence = _read_json(sequence_status)
            try:
                checkpoints = _validate_sequence_completion(sequence)
            except Exception as exc:
                fcntl.flock(sequence_lock_handle.fileno(), fcntl.LOCK_UN)
                sequence_lock_handle.close()
                sequence_lock_handle = None
                commit("waiting_sequence", waiting_reason=str(exc), **_progress_snapshot(sequence))
                time.sleep(args.poll_interval)
                continue

        args.posebusters_lock.parent.mkdir(parents=True, exist_ok=True)
        while posebusters_lock_handle is None:
            candidate = open(args.posebusters_lock, "a+")
            try:
                fcntl.flock(candidate.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                candidate.close()
                commit(
                    "waiting_posebusters_lock",
                    waiting_reason=(
                        "an existing/orphan PoseBusters process still holds the run lock"
                    ),
                )
                time.sleep(args.poll_interval)
                continue
            posebusters_lock_handle = candidate

        config = revalidate_pinned_sources()
        completed_snapshot = _progress_snapshot(sequence)
        commit(
            "gating_artifacts",
            waiting_reason=None,
            checkpoints=checkpoints,
            sequence_completion_snapshot=completed_snapshot,
            sequence_status_sha256=_sha256(sequence_status),
            sequence_finished_at=sequence.get("finished_at"),
            **completed_snapshot,
        )
        try:
            expected, input_report, verdict_columns = _build_input_manifest(
                config_path, manifest_path,
            )
        except Exception as exc:
            commit(
                "failed", failure_stage="artifact_gate", error=str(exc),
                finished_at=_now(), posebusters_pid=None,
            )
            return 1
        commit(
            "artifacts_ready",
            input_manifest_sha256=_sha256(manifest_path),
            input_coverage=input_report,
        )

        output_base = (REPO_ROOT / config["output_base_dir"]).resolve()
        results_csv = output_base / "dock/posebusters_filtered_results.csv"
        incomplete_file = output_base / "dock/posebusters_incomplete.txt"
        failures: list[dict[str, Any]] = []
        workers = args.workers
        for attempt in range(1, args.max_attempts + 1):
            config = revalidate_pinned_sources()
            command = [
                str(python), "-u", str(runner), "--config", str(config_path),
                "--quiet", "--monitor", "--monitor-interval",
                str(args.monitor_interval), "--workers", str(workers),
            ]
            print(
                f"[{_now()}] START PoseBusters attempt {attempt}/{args.max_attempts} "
                f"with {workers} workers",
                flush=True,
            )
            child = _launch_posebusters(
                command, sequence_lock_handle, posebusters_lock_handle,
            )
            while child.poll() is None:
                output_stat = results_csv.stat() if results_csv.exists() else None
                commit(
                    "posebusters_running",
                    attempt=attempt,
                    posebusters_pid=child.pid,
                    workers=workers,
                    command=command,
                    results_csv=str(results_csv),
                    results_size_bytes=(output_stat.st_size if output_stat else 0),
                    results_modified_at=(
                        dt.datetime.fromtimestamp(output_stat.st_mtime).astimezone().isoformat()
                        if output_stat else None
                    ),
                    prior_failures=failures,
                )
                time.sleep(args.poll_interval)
            returncode = child.wait()
            child = None
            print(f"[{_now()}] PoseBusters attempt {attempt} exited {returncode}", flush=True)

            try:
                if returncode != 0:
                    raise RuntimeError(f"PoseBusters exited with return code {returncode}")
                revalidate_pinned_sources()
                output_report = _verify_output(
                    results_csv, incomplete_file, expected, verdict_columns,
                )
            except Exception as exc:
                failures.append({
                    "attempt": attempt, "workers": workers,
                    "returncode": returncode, "error": str(exc), "at": _now(),
                })
                if attempt == args.max_attempts:
                    commit(
                        "failed", failure_stage="posebusters",
                        attempts=attempt, failures=failures,
                        finished_at=_now(), posebusters_pid=None,
                    )
                    return 1
                workers = max(1, workers // 2)
                commit(
                    "posebusters_retrying", attempt=attempt,
                    next_workers=workers, failures=failures,
                    posebusters_pid=None,
                )
                time.sleep(args.poll_interval)
                continue

            commit(
                "completed",
                finished_at=_now(),
                attempts=attempt,
                posebusters_returncode=returncode,
                output_coverage=output_report,
                results_csv=str(results_csv),
                results_csv_sha256=_sha256(results_csv),
                failures=failures,
                posebusters_pid=None,
                sequence_pid_alive=False,
            )
            print(f"[{_now()}] POSEBUSTERS SUPERVISOR COMPLETE", flush=True)
            return 0
        raise AssertionError("unreachable PoseBusters attempt loop")
    except SupervisorInterrupted as exc:
        stop_child()
        commit(
            "interrupted", finished_at=_now(), signal=exc.signum,
            error=str(exc), posebusters_pid=None,
        )
        return 128 + exc.signum
    except BaseException as exc:
        stop_child()
        commit(
            "failed", finished_at=_now(), error_type=type(exc).__name__,
            error=str(exc), posebusters_pid=None,
        )
        raise
    finally:
        stop_child()
        if posebusters_lock_handle is not None:
            fcntl.flock(posebusters_lock_handle.fileno(), fcntl.LOCK_UN)
            posebusters_lock_handle.close()
        if sequence_lock_handle is not None:
            fcntl.flock(sequence_lock_handle.fileno(), fcntl.LOCK_UN)
            sequence_lock_handle.close()
        fcntl.flock(watcher_lock.fileno(), fcntl.LOCK_UN)
        watcher_lock.close()


if __name__ == "__main__":
    raise SystemExit(main())
