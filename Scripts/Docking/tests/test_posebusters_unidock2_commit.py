import hashlib
import json
import sys
import types
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
# The default test interpreter has RDKit but not the optional PoseBusters
# package. These tests exercise collection/provenance only, so a tiny import
# stub keeps them hermetic without weakening the production module.
if "posebusters" not in sys.modules:
    fake_posebusters = types.ModuleType("posebusters")
    fake_posebusters.PoseBusters = type("PoseBusters", (), {})
    fake_posebusters.__file__ = __file__
    sys.modules["posebusters"] = fake_posebusters
from Scripts.Docking import run_unidock2 as ud2
from Scripts.Docking.Posebusters import run_posebusters as pb


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _provenance(path: Path) -> dict:
    path = path.resolve()
    return {"path": str(path), "size": path.stat().st_size, "sha256": _sha(path)}


def _tool_identity(tmp_path: Path) -> dict:
    env = tmp_path / "fake_unidock2_env"
    binary = env / "bin" / "unidock2"
    binary.parent.mkdir(parents=True)
    binary.write_text("#!/bin/sh\necho 'unidock2 0.6.3'\n")
    binary.chmod(0o755)
    for name in ("tleap", "teLeap", "ambpdb"):
        executable = binary.parent / name
        executable.write_text("#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)
    engine = (env / "lib" / "python3.10" / "site-packages" /
              "unidock_engine" / "api" / "python" / "pipeline.so")
    engine.parent.mkdir(parents=True)
    engine.write_text("fixture engine library\n")
    processing = (env / "lib" / "python3.10" / "site-packages" /
                  "unidock_processing")
    processing.mkdir(parents=True)
    (processing / "__init__.py").write_text("# fixture package\n")
    (processing / "pipeline.py").write_text("VALUE = 1\n")
    records = env / "conda-meta"
    records.mkdir()
    for name in ("unidock2-0.6.3-0.json", "ambertools_stable-23.0-0.json",
        "pdbfixer-1.9-0.json"):
        record = records / name
        record.write_text("{}\n")
    return ud2.probe_tool_identity(str(binary))


def _sdf_record(affinity: float, title: str) -> str:
    return (
        f"{title}\n"
        "  pytest\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        ">  <vina_binding_free_energy>\n"
        f"{affinity}\n\n"
        "$$$$\n"
    )


def _input_sdf_record(title: str = "input ligand") -> str:
    return (
        f"{title}\n"
        "  pytest\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        "$$$$\n"
    )


def _pdb_atom(serial: int, atom_name: str, element: str, x: float) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:>4s} ALA A   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def _committed_tree(tmp_path: Path) -> dict:
    combo = "5SAK_ZRY"
    out_dir = tmp_path / combo
    out_dir.mkdir()
    sdf = out_dir / f"{combo}_unidock2_out.sdf"
    sdf.write_text(_sdf_record(-8.0, "pose1") + _sdf_record(-7.0, "pose2"))
    prepared = out_dir / f"{combo}_unidock2_receptor_prepared.pdb"
    prepared.write_text(
        _pdb_atom(1, "C", "C", 0.0) + _pdb_atom(2, "H", "H", 1.0)
    )

    inputs_dir = tmp_path / "inputs"
    inputs_dir.mkdir()
    receptor = inputs_dir / f"{combo}_protein.pdb"
    ligand = inputs_dir / f"{combo}_ligand_start_conf.sdf"
    box = inputs_dir / f"{combo}_protein_mgl_tools.box.txt"
    receptor.write_text(_pdb_atom(1, "C", "C", 0.0))
    ligand.write_text(_input_sdf_record())
    box.write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
        "size_x = 10\nsize_y = 11\nsize_z = 12\n"
    )

    effective = {
        "Advanced": {
            "exhaustiveness": 512,
            "mc_steps": 40,
            "num_pose": 30,
            "energy_range": 6.0,
            "rmsd_limit": 1.0,
            "seed": 42,
            "randomize": True,
            "opt_steps": -1,
            "refine_steps": 5,
            "use_tor_lib": False,
            "energy_decomp": False,
        },
        "Hardware": {"gpu_device_id": 0, "n_cpu": 1},
        "Settings": {
            "box_size": [10.0, 11.0, 12.0],
            "task": "screen",
            "search_mode": "free",
        },
        "Preprocessing": {
            "construct_ff": False,
            "template_docking": False,
            "compute_center": True,
            "covalent_ligand": False,
            "preserve_receptor_hydrogen": False,
            "engine_checkpoint": False,
        },
    }
    effective_json = json.dumps(
        effective, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    provenance = {
        "schema_version": 3,
        "engine": "unidock2",
        "tool": _tool_identity(tmp_path),
        "inputs": {
            "receptor": _provenance(receptor),
            "ligand": _provenance(ligand),
            "box": _provenance(box),
        },
        "center": [0.0, 0.0, 0.0],
        "size": [10.0, 11.0, 12.0],
        "effective_config": effective,
        "effective_config_sha256": hashlib.sha256(effective_json.encode()).hexdigest(),
        "driver_execution": {
            "gpu_lock_file": "/tmp/master_docking_gpu0.lock",
        },
    }
    canonical = json.dumps(
        provenance, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    runtime = ud2._requested_runtime_settings(provenance)
    manifest = {
        **provenance,
        "status": "success",
        "fingerprint": hashlib.sha256(canonical.encode()).hexdigest(),
        "generation_id": "a" * 32,
        "output_file": sdf.name,
        "output_sha256": _sha(sdf),
        "prepared_receptor_file": prepared.name,
        "prepared_receptor_sha256": _sha(prepared),
        "prepared_receptor_atoms": 2,
        "prepared_receptor_heavy_atoms": 1,
        "engine_receptor_heavy_atoms_in_box": 1,
        "receptor_prmtop_sha256": "b" * 64,
        "receptor_inpcrd_sha256": "c" * 64,
        "num_poses": 2,
        "best_affinity": -8.0,
        "requested_search_settings": runtime,
        "runtime_attestation": runtime,
        "effective_search_settings": runtime,
        "applied_center": [0.0, 0.0, 0.0],
        "applied_box_size": [10.0, 11.0, 12.0],
        "input_protein_heavy_atoms": 1,
        "completed_unix_s": 1.0,
    }
    completion = out_dir / f"{combo}_unidock2_completion.json"
    completion.write_text(json.dumps(manifest))
    return {
        "root": tmp_path,
        "out_dir": out_dir,
        "sdf": sdf,
        "prepared": prepared,
        "completion": completion,
        "manifest": manifest,
        "inputs": {"receptor": receptor, "ligand": ligand, "box": box},
    }


def _resign_manifest(tree: dict, manifest: dict) -> None:
    """Make a forged protocol internally consistent before testing rejection."""
    effective_json = json.dumps(
        manifest["effective_config"],
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    manifest["effective_config_sha256"] = hashlib.sha256(
        effective_json.encode()
    ).hexdigest()
    provenance = {key: manifest[key] for key in pb._UNIDOCK2_PROVENANCE_KEYS}
    canonical = json.dumps(
        provenance, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    manifest["fingerprint"] = hashlib.sha256(canonical.encode()).hexdigest()
    runtime = ud2._requested_runtime_settings(provenance)
    manifest["requested_search_settings"] = runtime
    manifest["runtime_attestation"] = dict(runtime)
    manifest["effective_search_settings"] = dict(runtime)
    manifest["applied_center"] = list(runtime["center"])
    manifest["applied_box_size"] = list(runtime["box_size"])
    tree["completion"].write_text(json.dumps(manifest))


def test_completion_is_authoritative_and_summary_is_optional(tmp_path):
    tree = _committed_tree(tmp_path)

    rows = pb._collect_unidock2_rows(tree["root"])

    assert len(rows) == 1
    row = rows[0]
    assert "docking_summary_file" not in row
    assert row["unidock2_completion_file"] == str(tree["completion"].resolve())
    assert row["docking_receptor_file"] == str(tree["prepared"].resolve())
    assert row["unidock2_prepared_receptor_heavy_atoms"] == 1

    # A stale resume summary is an optional index, not a second commit protocol.
    summary = tree["out_dir"] / "docking_summary.json"
    summary.write_text(json.dumps({"status": "done", "generation_id": None}))
    rows = pb._collect_unidock2_rows(tree["root"])
    assert len(rows) == 1
    assert rows[0]["docking_summary_file"] == str(summary.resolve())


def test_legacy_summary_only_output_is_rejected(tmp_path):
    tree = _committed_tree(tmp_path)
    tree["completion"].unlink()
    (tree["out_dir"] / "docking_summary.json").write_text(
        json.dumps({"status": "success", "num_poses": 2})
    )

    assert pb._collect_unidock2_rows(tree["root"]) == []


def test_self_consistent_schema2_completion_is_explicitly_rejected(tmp_path):
    tree = _committed_tree(tmp_path)
    manifest = json.loads(tree["completion"].read_text())
    manifest["schema_version"] = 2
    provenance = {
        key: manifest[key] for key in (
            "schema_version", "engine", "tool", "inputs", "center", "size",
            "effective_config", "effective_config_sha256", "driver_execution",
        )
    }
    canonical = json.dumps(
        provenance, sort_keys=True, separators=(",", ":"), allow_nan=False)
    manifest["fingerprint"] = hashlib.sha256(canonical.encode()).hexdigest()
    tree["completion"].write_text(json.dumps(manifest))

    assert pb._collect_unidock2_rows(tree["root"]) == []


@pytest.mark.parametrize(
    "tamper",
    [
        "different_exhaustiveness",
        "different_randomize",
        "different_num_pose",
        "different_preprocessing",
        "different_hardware",
        "extra_advanced_key",
        "attempt_temp_dir",
        "top_level_size_mismatch",
        "string_integer",
        "bool_integer",
        "integer_box_coordinates",
        "different_gpu_lock",
    ],
)
def test_self_consistent_nonbenchmark_protocol_is_rejected(tmp_path, tamper):
    tree = _committed_tree(tmp_path)
    manifest = json.loads(tree["completion"].read_text())
    effective = manifest["effective_config"]
    if tamper == "different_exhaustiveness":
        effective["Advanced"]["exhaustiveness"] = 256
    elif tamper == "different_randomize":
        effective["Advanced"]["randomize"] = False
    elif tamper == "different_num_pose":
        effective["Advanced"]["num_pose"] = 31
    elif tamper == "different_preprocessing":
        effective["Preprocessing"]["compute_center"] = False
    elif tamper == "different_hardware":
        effective["Hardware"]["n_cpu"] = 2
    elif tamper == "extra_advanced_key":
        effective["Advanced"]["unsupported_knob"] = True
    elif tamper == "attempt_temp_dir":
        effective["Preprocessing"]["temp_dir_name"] = "/tmp/forged-attempt"
    elif tamper == "top_level_size_mismatch":
        manifest["size"] = [99.0, 99.0, 99.0]
    elif tamper == "string_integer":
        effective["Advanced"]["exhaustiveness"] = "512"
    elif tamper == "bool_integer":
        effective["Advanced"]["opt_steps"] = True
    elif tamper == "integer_box_coordinates":
        effective["Settings"]["box_size"] = [10, 11, 12]
    elif tamper == "different_gpu_lock":
        manifest["driver_execution"]["gpu_lock_file"] = str(
            (tmp_path / "private-gpu.lock").resolve()
        )
    _resign_manifest(tree, manifest)

    assert pb._collect_unidock2_rows(tree["root"]) == []


@pytest.mark.parametrize(
    "tamper",
    [
        "input", "output_order", "fingerprint", "prepared",
        "requested_runtime", "runtime", "effective_runtime", "applied_box",
        "missing_runtime_key", "extra_runtime_key", "generation_id",
        "tool_artifact", "artifact_omission", "python_source",
        "wrong_input_basename", "resigned_box_geometry",
        "invalid_ligand_resigned",
    ],
)
def test_collector_rejects_tampered_commit_artifacts(tmp_path, tamper):
    tree = _committed_tree(tmp_path)
    if tamper == "input":
        tree["inputs"]["ligand"].write_text("changed staged ligand bytes\n")
    elif tamper == "output_order":
        tree["sdf"].write_text(_sdf_record(-7.0, "pose1") + _sdf_record(-8.0, "pose2"))
        manifest = json.loads(tree["completion"].read_text())
        manifest["output_sha256"] = _sha(tree["sdf"])
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "fingerprint":
        manifest = json.loads(tree["completion"].read_text())
        manifest["effective_config"]["Advanced"]["mc_steps"] = 41
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "prepared":
        tree["prepared"].write_text(_pdb_atom(1, "N", "N", 0.0))
    elif tamper == "requested_runtime":
        manifest = json.loads(tree["completion"].read_text())
        manifest["requested_search_settings"]["mc_steps"] = 41
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "runtime":
        manifest = json.loads(tree["completion"].read_text())
        manifest["runtime_attestation"]["mc_steps"] = 41
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "effective_runtime":
        manifest = json.loads(tree["completion"].read_text())
        manifest["effective_search_settings"]["mc_steps"] = 41
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "applied_box":
        manifest = json.loads(tree["completion"].read_text())
        manifest["applied_box_size"] = [10.0, 11.0, 13.0]
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "missing_runtime_key":
        manifest = json.loads(tree["completion"].read_text())
        del manifest["runtime_attestation"]["engine_checkpoint"]
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "extra_runtime_key":
        manifest = json.loads(tree["completion"].read_text())
        manifest["runtime_attestation"]["temp_dir_name"] = "/tmp/forged"
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "generation_id":
        manifest = json.loads(tree["completion"].read_text())
        manifest["generation_id"] = "generation-1"
        tree["completion"].write_text(json.dumps(manifest))
    elif tamper == "tool_artifact":
        manifest = json.loads(tree["completion"].read_text())
        Path(manifest["tool"]["path"]).write_text("changed executable bytes\n")
    elif tamper == "artifact_omission":
        manifest = json.loads(tree["completion"].read_text())
        manifest["tool"]["artifacts"] = [
            item for item in manifest["tool"]["artifacts"]
            if not Path(item["path"]).name.startswith("pipeline")
        ]
        _resign_manifest(tree, manifest)
    elif tamper == "python_source":
        manifest = json.loads(tree["completion"].read_text())
        source_root = Path(manifest["tool"]["python_source_tree"]["path"])
        (source_root / "pipeline.py").write_text("VALUE = 2\n")
    elif tamper == "wrong_input_basename":
        manifest = json.loads(tree["completion"].read_text())
        substitute = tmp_path / "5SB2_1K2_protein.pdb"
        substitute.write_bytes(tree["inputs"]["receptor"].read_bytes())
        manifest["inputs"]["receptor"] = _provenance(substitute)
        _resign_manifest(tree, manifest)
    elif tamper == "resigned_box_geometry":
        tree["inputs"]["box"].write_text(
            "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
            "size_x = 13\nsize_y = 11\nsize_z = 12\n"
        )
        manifest = json.loads(tree["completion"].read_text())
        manifest["inputs"]["box"] = _provenance(tree["inputs"]["box"])
        _resign_manifest(tree, manifest)
    elif tamper == "invalid_ligand_resigned":
        tree["inputs"]["ligand"].write_text("not an SDF\n")
        manifest = json.loads(tree["completion"].read_text())
        manifest["inputs"]["ligand"] = _provenance(tree["inputs"]["ligand"])
        _resign_manifest(tree, manifest)

    assert pb._collect_unidock2_rows(tree["root"]) == []


def test_expansion_and_resolver_require_sibling_prepared_receptor(tmp_path):
    tree = _committed_tree(tmp_path)
    row = pb._collect_unidock2_rows(tree["root"])[0]
    poses = pb._expand_unidock2_poses(row, tmp_path / "split", None)

    assert len(poses) == 2
    assert poses[0]["unidock2_fingerprint"] == tree["manifest"]["fingerprint"]
    assert poses[0]["unidock2_generation_id"] == "a" * 32
    receptor, source, stage = pb._resolve_unidock2_receptor(poses[0], tmp_path)
    assert receptor == source == tree["prepared"]
    assert stage == "unidock2_manifest_attested_engine_prepared_pdb"

    # Even if both a pose row and manifest are rewritten consistently after
    # collection, the resolver must not switch to a different protocol.
    original_manifest = tree["completion"].read_text()
    original_fingerprint = poses[0]["unidock2_fingerprint"]
    forged = json.loads(original_manifest)
    forged["effective_config"]["Advanced"]["randomize"] = False
    _resign_manifest(tree, forged)
    poses[0]["unidock2_fingerprint"] = forged["fingerprint"]
    with pytest.raises(ValueError, match="provenance is stale"):
        pb._resolve_unidock2_receptor(poses[0], tmp_path)
    tree["completion"].write_text(original_manifest)
    poses[0]["unidock2_fingerprint"] = original_fingerprint

    substitute = tmp_path / tree["prepared"].name
    substitute.write_bytes(tree["prepared"].read_bytes())
    poses[0]["unidock2_prepared_receptor_file"] = str(substitute)
    with pytest.raises(ValueError, match="does not point beside"):
        pb._resolve_unidock2_receptor(poses[0], tmp_path)
