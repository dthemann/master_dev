from __future__ import annotations

import csv
import fcntl
import hashlib
import json
import sys
from pathlib import Path

import pytest

from Scripts.Docking import run_posebusters_after_autodock_sequence as supervisor


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value))


def _identity(index: int, arm: tuple[str, str]) -> dict:
    identity = {
        field: f"{field}-{index}" for field in supervisor._IDENTITY_FIELDS
    }
    for field in supervisor._INTEGER_IDENTITY_FIELDS:
        identity[field] = 1
    for field in supervisor._FLOAT_IDENTITY_FIELDS:
        identity[field] = float(index + 1)
    source_scoring, cnn_scoring = supervisor.EXPECTED_PROTOCOLS[arm]
    identity.update({
        "docking_method": arm[0],
        "optimizer": arm[1],
        "protein": f"protein-{index}",
        "ligand": f"ligand-{index}",
        "source_scoring": source_scoring,
        "optimizer_scoring": "default",
        "optimizer_search": "minimize",
        "cnn_scoring": cnn_scoring,
        "cnn_model": "crossdock_default2018_ensemble",
        "rank_metric": "cnn_affinity",
    })
    return identity


def test_static_config_is_strict() -> None:
    raw = supervisor._validate_static_config(supervisor.DEFAULT_CONFIG)
    assert raw["require_all_methods"] is True
    assert raw["require_complete_results"] is True
    assert raw["expected_common_combinations"] == 308
    assert set(raw["variant_filter"]["autodock"]) == {"gnina", "gnina_refinement"}


def test_sequence_gate_requires_hashes_cells_and_checkpoints(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    notebook = tmp_path / "notebook.ipynb"
    config = tmp_path / "vina.yaml"
    vinardo = tmp_path / "vinardo.yaml"
    config.write_text("scoring: vina\n")
    vinardo.write_text("scoring: vinardo\n")
    cells = []
    stages = []
    for index, name in enumerate(supervisor.EXPECTED_STAGE_NAMES):
        cell_id = f"cell-{index}"
        source = f"print({name!r})\n"
        cells.append({"id": cell_id, "cell_type": "code", "source": [source]})
        stages.append({
            "name": name, "cell_id": cell_id, "status": "completed",
            "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
        })
    _write_json(notebook, {"cells": cells})

    checkpoints = {}
    for name in supervisor.EXPECTED_STAGE_NAMES:
        path = tmp_path / f"{name}.json"
        expected = dict(supervisor.CHECKPOINTS[name][1])
        count_field = (
            "completed_complexes" if name == "gnina_refinement"
            else "processed_complexes"
        )
        _write_json(path, {
            **expected, "status": "completed", count_field: 308,
            "issues": [], "optimizer_issues": [], "optimizer_totals": {"failed": 0},
        })
        checkpoints[name] = (path, supervisor.CHECKPOINTS[name][1])
    monkeypatch.setattr(supervisor, "CHECKPOINTS", checkpoints)

    payload = {
        "status": "completed", "current_stage": None, "stages": stages,
        "notebook": str(notebook), "notebook_sha256": supervisor._sha256(notebook),
        "config": str(config), "config_sha256": supervisor._sha256(config),
        "vinardo_config": str(vinardo),
        "vinardo_config_sha256": supervisor._sha256(vinardo),
    }
    # An unrelated notebook edit after launch must not invalidate the four
    # separately pinned producer-cell hashes.
    _write_json(notebook, {"cells": cells + [{
        "id": "unrelated", "cell_type": "markdown", "source": ["updated docs\n"],
    }]})
    report = supervisor._validate_sequence_completion(payload)
    assert set(report) == set(supervisor.EXPECTED_STAGE_NAMES)

    stages[-1]["status"] = "running"
    with pytest.raises(RuntimeError, match="not every sequence stage"):
        supervisor._validate_sequence_completion(payload)


def test_output_postflight_requires_exact_four_arm_signatures(tmp_path: Path) -> None:
    expected = {}
    rows = []
    verdict_columns = ["mol_pred_loaded", "sanitization"]
    for index, arm in enumerate(supervisor.EXPECTED_ARMS):
        signature = f"sig-{index}"
        pair = (f"protein-{index}", f"ligand-{index}")
        identity = _identity(index, arm)
        expected[arm] = {
            "records": {signature: identity}, "pairs": {pair},
            "complex_counts": {"ID": 1},
        }
        rows.append({
            **identity, "validation_signature": signature,
            **{column: "True" for column in verdict_columns},
        })
    output = tmp_path / "posebusters_filtered_results.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    report = supervisor._verify_output(
        output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
    )
    assert set(report) == {":".join(arm) for arm in supervisor.EXPECTED_ARMS}

    original_score = rows[0]["cnn_affinity"]
    rows[0]["cnn_affinity"] = float(original_score) + 0.5
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="output identity mismatch"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )
    rows[0]["cnn_affinity"] = original_score

    original_path = rows[0]["pose_file"]
    rows[0]["pose_file"] = f"{original_path}.corrupt"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="output identity mismatch"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )
    rows[0]["pose_file"] = original_path

    original_protocol = rows[0]["cnn_scoring"]
    rows[0]["cnn_scoring"] = "refinement"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="output identity mismatch"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )
    rows[0]["cnn_scoring"] = original_protocol

    rows.pop()
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="output signature mismatch"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )


def test_output_postflight_rejects_metadata_only_or_permuted_rows(tmp_path: Path) -> None:
    verdict_columns = ["mol_pred_loaded"]
    expected = {}
    rows = []
    for index, arm in enumerate(supervisor.EXPECTED_ARMS):
        signature = f"sig-{index}"
        identity = _identity(index, arm)
        pair = (identity["protein"], identity["ligand"])
        expected[arm] = {
            "records": {signature: identity}, "pairs": {pair},
            "complex_counts": {"ID": 1},
        }
        rows.append({**identity, "validation_signature": signature})

    output = tmp_path / "posebusters_filtered_results.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="missing required provenance/verdict"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )

    for row in rows:
        row["mol_pred_loaded"] = "True"
    rows[0]["protein"], rows[1]["protein"] = rows[1]["protein"], rows[0]["protein"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(RuntimeError, match="output identity mismatch"):
        supervisor._verify_output(
            output, tmp_path / "posebusters_incomplete.txt", expected, verdict_columns,
        )


def test_pose_scope_requires_canonical_container_and_protocol(tmp_path: Path) -> None:
    complex_id = "7ABC_XYZ"
    arm = ("autodock_vinardo", "gnina_refinement")
    docking = tmp_path / complex_id / "mgl_tools" / "docking"
    docking.mkdir(parents=True)
    source = docking / (
        f"{complex_id}_protein_mgl_tools__{complex_id}_ligand_start_conf_"
        "vinardo_vina_out.pdbqt"
    )
    source.write_text("MODEL 1\nENDMDL\n")
    optimized = docking / "optimized_gnina_refinement" / (
        f"{source.stem}_rank1_gnina_refinement.sdf"
    )
    optimized.parent.mkdir()
    optimized.write_text("pose\n$$$$\n")
    sidecar = Path(f"{optimized}.provenance.json")
    sidecar.write_text("{}\n")
    log = docking.parent / "optimization_log.csv"
    log.write_text("header\n")
    pose = {
        "protein": complex_id,
        "ligand": complex_id,
        "source_pose_file": str(source),
        "pose_file": str(optimized),
        "autodock_rank": 1,
        "optimization_log_file": str(log),
        "optimizer_provenance_file": str(sidecar),
        "source_scoring": "vinardo",
        "optimizer_scoring": "default",
        "optimizer_search": "minimize",
        "cnn_scoring": "refinement",
        "cnn_model": "crossdock_default2018_ensemble",
        "rank_metric": "cnn_affinity",
    }
    supervisor._validate_pose_scope(pose, arm, tmp_path, complex_id)

    pose["cnn_scoring"] = "rescore"
    with pytest.raises(RuntimeError, match="protocol mismatch"):
        supervisor._validate_pose_scope(pose, arm, tmp_path, complex_id)
    pose["cnn_scoring"] = "refinement"

    extra = docking / "nested" / source.name
    extra.parent.mkdir()
    extra.write_text(source.read_text())
    pose["source_pose_file"] = str(extra)
    with pytest.raises(RuntimeError, match="canonical raw container"):
        supervisor._validate_pose_scope(pose, arm, tmp_path, complex_id)


def test_posebusters_child_inherits_both_campaign_locks(tmp_path: Path) -> None:
    lock_paths = [tmp_path / "sequence.lock", tmp_path / "posebusters.lock"]
    handles = [path.open("a+") for path in lock_paths]
    for handle in handles:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    process = supervisor._launch_posebusters(
        [sys.executable, "-c", "import time; time.sleep(60)"],
        handles[0], handles[1],
    )
    try:
        # Simulate an abruptly killed supervisor: closing (without LOCK_UN)
        # leaves each open-file-description lock held by the surviving child.
        for handle in handles:
            handle.close()
        probes = [path.open("a+") for path in lock_paths]
        try:
            for probe in probes:
                with pytest.raises(BlockingIOError):
                    fcntl.flock(probe.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        finally:
            for probe in probes:
                probe.close()
    finally:
        supervisor._stop_process_group(process, timeout_s=5)

    for path in lock_paths:
        with path.open("a+") as probe:
            fcntl.flock(probe.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            fcntl.flock(probe.fileno(), fcntl.LOCK_UN)
