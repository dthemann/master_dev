from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest
import yaml

from Scripts.Docking import run_autodock_vinardo_raw_prefill as prefill
from Scripts.Docking import run_autodock_vinardo_gnina_rescore_early as early_rescore
from Scripts.Docking import watch_autodock_vinardo_raw_prefill as watchdog


def _vinardo_cell_source(repo_root: Path) -> str:
    notebook = json.loads((repo_root / "Master_Docking_AD_Full_Protein.ipynb").read_text())
    matches = [
        cell for cell in notebook["cells"]
        if cell.get("id") == prefill.VINARDO_CELL_ID
    ]
    assert len(matches) == 1
    return "".join(matches[0]["source"])


def test_transform_is_raw_only_separate_and_boundary_guarded() -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    transformed = prefill._transform_cell_source(source)

    assert "OPT_TOOLS = ()" in transformed
    assert "OPT_PREFLIGHT = None" in transformed
    assert "optimize_autodock_results = _raw_prefill_optimizer_forbidden" in transformed
    assert "preflight_optimizers(base_cfg)" not in transformed
    assert '"optimization": "none"' in transformed
    assert '"cnn_scoring": "none"' in transformed
    assert '"stage": "raw_vinardo_prefill"' in transformed
    assert prefill.PREFILL_STATUS_NAME in transformed
    assert "benchmark_gnina_rescore_status.json" not in transformed
    assert (
        "for idx, cdir in enumerate(complex_dirs, 1):\n"
        "    _assert_raw_prefill_window()"
    ) in transformed
    assert "complex_dirs = list(reversed(complex_dirs))" in transformed
    assert (
        "_assert_raw_prefill_window()\n"
        "        results_df, results = run_autodock_vina("
    ) in transformed
    # Existing-pose and fresh-pose optimizer paths are both disabled by the
    # same immutable false-y tuple; the reviewed docking calls remain present.
    assert "if not OPT_TOOLS:\n        return" in transformed
    assert "if OPT_TOOLS:\n            optimize_with_reporting(results, vina_out)" in transformed
    assert "run_autodock_vina(" in transformed
    compile(transformed, "test-raw-prefill", "exec")


def test_transform_fails_closed_on_source_drift() -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    drifted = source.replace(
        'RESCORE_STATUS_PATH = output_base / "benchmark_gnina_rescore_status.json"',
        'RESCORE_STATUS_PATH = output_base / "changed.json"',
    )
    with pytest.raises(prefill.PrefillGuardError, match="status path"):
        prefill._transform_cell_source(drifted)


def test_early_rescore_transform_is_optimizer_only_and_reverse_guarded() -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    transformed = early_rescore._transform_rescore_source(source)

    assert "OPT_TOOLS   = resolve_optimizers(base_cfg)" in transformed
    assert "OPT_PREFLIGHT = preflight_optimizers(base_cfg)" in transformed
    assert "run_autodock_vina = _early_rescore_docking_forbidden" in transformed
    assert "complex_dirs = list(reversed(complex_dirs))" in transformed
    assert transformed.count("_assert_early_rescore_window()") == 3
    assert '"optimization": "gnina"' in transformed
    assert '"cnn_scoring": "rescore"' in transformed
    assert '"early_rescore": True' in transformed
    assert "benchmark_gnina_rescore_status.json" in transformed
    compile(transformed, "test-early-rescore", "exec")


def test_early_rescore_transform_fails_closed_on_docking_guard_drift() -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    drifted = source.replace("OPT_PREFLIGHT = preflight_optimizers(base_cfg)", "OPT_PREFLIGHT = None")
    with pytest.raises(prefill.PrefillGuardError, match="optimizer-only guard"):
        early_rescore._transform_rescore_source(drifted)


def _write_launch_fixture(tmp_path: Path, source: str) -> tuple[Path, Path, Path]:
    notebook = tmp_path / "notebook.ipynb"
    notebook.write_text(json.dumps({
        "cells": [{
            "id": prefill.VINARDO_CELL_ID,
            "cell_type": "code",
            "source": source.splitlines(keepends=True),
        }]
    }))
    config = tmp_path / "vinardo.yaml"
    config.write_text(yaml.safe_dump({
        "scoring_function": "vinardo",
        "overwrite_existing": False,
        "overwrite_poses": False,
        "output_dir": "outputs",
    }))
    vina_config = tmp_path / "vina.yaml"
    vina_config.write_text(yaml.safe_dump({"output_dir": "vina_outputs"}))
    refinement_status = tmp_path / "vina_outputs" / "benchmark_gnina_refinement_status.json"
    refinement_status.parent.mkdir()
    refinement_status.write_text(json.dumps({
        "status": "running", "completed_complexes": 100,
    }))
    sequence = tmp_path / "sequence.json"
    sequence.write_text(json.dumps({
        "status": "running",
        "current_stage": prefill.REQUIRED_SEQUENCE_STAGE,
        "pid": 12345,
        "process_group_id": 12345,
        "notebook": str(notebook),
        "config": str(vina_config),
        "config_sha256": hashlib.sha256(vina_config.read_bytes()).hexdigest(),
        "vinardo_config": str(config),
        "vinardo_config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
        "stages": [{
            "name": prefill.VINARDO_STAGE_NAME,
            "cell_id": prefill.VINARDO_CELL_ID,
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
            "status": "pending",
        }],
    }))
    return sequence, notebook, config


def test_launch_context_pins_live_sequence_cell_and_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    sequence, _, _ = _write_launch_fixture(tmp_path, source)
    monkeypatch.setattr(prefill, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(prefill, "_process_start_time", lambda pid: "start-123")
    monkeypatch.setattr(
        prefill, "_process_cmdline",
        lambda pid: ["python", "run_autodock_gnina_sequence.py", "--status-path", str(sequence)],
    )
    monkeypatch.setattr(prefill.os, "getpgid", lambda pid: 12345)
    monkeypatch.setattr(prefill, "_lock_is_held", lambda path: True)

    context = prefill._load_launch_context(sequence)

    assert context["sequence_pid"] == 12345
    assert context["sequence_start_time"] == "start-123"
    assert context["source_sha256"] == hashlib.sha256(source.encode()).hexdigest()
    assert context["output_dir"] == (tmp_path / "outputs").resolve()
    assert context["completed_refinement"] == 100
    assert '"optimization": "none"' in context["transformed_source"]


@pytest.mark.parametrize("mutation,match", [
    ("config", "YAML changed"),
    ("notebook", "notebook cell changed"),
    ("stage", "not pending"),
])
def test_launch_context_rejects_pinned_input_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mutation: str, match: str,
) -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    sequence, notebook, config = _write_launch_fixture(tmp_path, source)
    monkeypatch.setattr(prefill, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(prefill, "_process_start_time", lambda pid: "start-123")
    monkeypatch.setattr(
        prefill, "_process_cmdline",
        lambda pid: ["python", "run_autodock_gnina_sequence.py", "--status-path", str(sequence)],
    )
    monkeypatch.setattr(prefill.os, "getpgid", lambda pid: 12345)
    monkeypatch.setattr(prefill, "_lock_is_held", lambda path: True)

    if mutation == "config":
        config.write_text(config.read_text() + "# drift\n")
    elif mutation == "notebook":
        payload = json.loads(notebook.read_text())
        payload["cells"][0]["source"].append("# drift\n")
        notebook.write_text(json.dumps(payload))
    else:
        payload = json.loads(sequence.read_text())
        payload["stages"][0]["status"] = "running"
        sequence.write_text(json.dumps(payload))

    with pytest.raises(prefill.PrefillGuardError, match=match):
        prefill._load_launch_context(sequence)


def test_launch_context_refuses_overwrite(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    sequence, _, config = _write_launch_fixture(tmp_path, source)
    config_payload = yaml.safe_load(config.read_text())
    config_payload["overwrite_existing"] = True
    config.write_text(yaml.safe_dump(config_payload))
    sequence_payload = json.loads(sequence.read_text())
    sequence_payload["vinardo_config_sha256"] = hashlib.sha256(config.read_bytes()).hexdigest()
    sequence.write_text(json.dumps(sequence_payload))
    monkeypatch.setattr(prefill, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(prefill, "_process_start_time", lambda pid: "start-123")
    monkeypatch.setattr(
        prefill, "_process_cmdline",
        lambda pid: ["python", "run_autodock_gnina_sequence.py", "--status-path", str(sequence)],
    )
    monkeypatch.setattr(prefill.os, "getpgid", lambda pid: 12345)
    monkeypatch.setattr(prefill, "_lock_is_held", lambda path: True)

    with pytest.raises(prefill.PrefillGuardError, match="overwrite"):
        prefill._load_launch_context(sequence)


def test_launch_context_rejects_wrong_process_identity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = _vinardo_cell_source(prefill.REPO_ROOT)
    sequence, _, _ = _write_launch_fixture(tmp_path, source)
    monkeypatch.setattr(prefill, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(prefill, "_process_start_time", lambda pid: "start-123")
    monkeypatch.setattr(prefill, "_process_cmdline", lambda pid: ["python", "unrelated.py"])

    with pytest.raises(prefill.PrefillGuardError, match="guarded AutoDock launcher"):
        prefill._load_launch_context(sequence)


def test_raw_terminal_status_replaces_stale_running_state(tmp_path: Path) -> None:
    status = tmp_path / "raw.json"
    status.write_text(json.dumps({"status": "running", "processed_complexes": 10}))

    prefill._write_raw_terminal_status(status, "interrupted", "test stop")

    payload = json.loads(status.read_text())
    assert payload["status"] == "interrupted"
    assert payload["processed_complexes"] == 10
    assert payload["terminal_reason"] == "test stop"


def test_lock_probe_detects_other_process_flock(tmp_path: Path) -> None:
    lock_path = tmp_path / "held.lock"
    code = (
        "import fcntl,sys,time; "
        "h=open(sys.argv[1],'a+'); fcntl.flock(h.fileno(),fcntl.LOCK_EX); "
        "print('held',flush=True); time.sleep(10)"
    )
    process = subprocess.Popen(
        [sys.executable, "-c", code, str(lock_path)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
    )
    try:
        assert process.stdout is not None
        assert process.stdout.readline().strip() == "held"
        assert prefill._lock_is_held(lock_path)
    finally:
        process.terminate()
        process.wait(timeout=10)
    assert not prefill._lock_is_held(lock_path)


def test_watchdog_process_identity_is_fail_closed(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(watchdog, "_process_start_time", lambda pid: "ticks")
    monkeypatch.setattr(watchdog.os, "getpgid", lambda pid: pid)
    monkeypatch.setattr(
        watchdog, "_process_cmdline",
        lambda pid: ["python", "run_autodock_vinardo_raw_prefill.py"],
    )
    assert watchdog._same_prefill_process(9000, 9000, "ticks")
    assert not watchdog._same_prefill_process(9000, 9001, "ticks")
    assert not watchdog._same_prefill_process(9000, 9000, "other")
    assert watchdog._same_target_process(
        9000, 9000, "ticks", "run_autodock_vinardo_raw_prefill.py",
    )
    assert not watchdog._same_target_process(9000, 9000, "ticks", "different.py")


def test_watchdog_kills_child_when_group_leader_exits_on_term() -> None:
    leader_code = r'''
import os
import signal
import subprocess
import sys
import time

child = subprocess.Popen([
    sys.executable, "-c",
    "import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); "
    "print('ready', flush=True); time.sleep(60)",
], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
assert child.stdout is not None and child.stdout.readline().strip() == "ready"
print(os.getpid(), child.pid, flush=True)
signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
while True:
    time.sleep(1)
'''
    leader = subprocess.Popen(
        [sys.executable, "-c", leader_code],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        start_new_session=True,
    )
    try:
        assert leader.stdout is not None
        leader_pid, child_pid = map(int, leader.stdout.readline().split())
        assert leader_pid == leader.pid
        identities = watchdog._live_group_identities(leader.pid)
        assert leader.pid in identities and child_pid in identities

        escalated = watchdog._terminate_process_group(leader.pid, term_grace=0.5)

        assert escalated
        leader.wait(timeout=10)
        assert not watchdog._live_group_identities(leader.pid)
    finally:
        if leader.poll() is None:
            os.killpg(leader.pid, signal.SIGKILL)
            leader.wait(timeout=10)
