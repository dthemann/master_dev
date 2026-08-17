import signal
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from Scripts.Docking import run_unidock2 as ud2


def pdb_atom(serial=1, x=0.0):
    return (
        f"ATOM  {serial:5d}  C   ALA A   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          C "
    )


def sdf_record(affinity=-7.0, title="pose"):
    affinity_block = (
        ">  <vina_binding_free_energy>\n"
        f"{affinity}\n\n"
        if affinity is not None
        else ""
    )
    return (
        f"{title}\n"
        "  pytest\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        f"{affinity_block}"
        "$$$$\n"
    )


def make_complex(base: Path, name="5SAK_ZRY", *, malformed_box=False, missing=False):
    cdir = base / name
    cdir.mkdir(parents=True, exist_ok=True)
    if missing:
        return cdir

    receptor_dir = cdir / "_staging" / "receptors"
    box_dir = receptor_dir / "pdbqt"
    ligand_dir = cdir / "_staging" / "ligands"
    box_dir.mkdir(parents=True)
    ligand_dir.mkdir(parents=True)

    receptor = receptor_dir / f"{name}_protein.pdb"
    receptor.write_text(pdb_atom() + "\n")
    ligand = ligand_dir / f"{name}_ligand_start_conf.sdf"
    ligand.write_text(sdf_record(None, title="input ligand"))
    box = box_dir / f"{name}_protein_mgl_tools.box.txt"
    if malformed_box:
        box.write_text(
            "center_x = not-a-number\ncenter_y = 0\ncenter_z = 0\n"
            "size_x = 10\nsize_y = 10\nsize_z = 10\n"
        )
    else:
        box.write_text(
            "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
            "size_x = 10\nsize_y = 11\nsize_z = 12\n"
        )
    return cdir


def make_cfg(tmp_path: Path, ids=("5SAK_ZRY",), **overrides):
    prep = tmp_path / "prep"
    prep.mkdir(exist_ok=True)
    fake_bin = tmp_path / "fake_env" / "bin" / "unidock2"
    fake_bin.parent.mkdir(parents=True, exist_ok=True)
    fake_bin.write_text("#!/bin/sh\necho 'unidock2 0.6.3'\n")
    fake_bin.chmod(0o755)
    # Keep any future dependency preflight hermetic too.
    tleap = fake_bin.parent / "tleap"
    tleap.write_text("#!/bin/sh\nexit 0\n")
    tleap.chmod(0o755)
    teleap = fake_bin.parent / "teLeap"
    teleap.write_text("#!/bin/sh\nexit 0\n")
    teleap.chmod(0o755)
    ambpdb = fake_bin.parent / "ambpdb"
    ambpdb.write_text("#!/bin/sh\nprintf '%s\\n' '" + pdb_atom() + "'\n")
    ambpdb.chmod(0o755)
    engine_library = (
        fake_bin.parent.parent / "lib" / "python3.10" / "site-packages"
        / "unidock_engine" / "api" / "python" / "pipeline_fake.so"
    )
    engine_library.parent.mkdir(parents=True, exist_ok=True)
    engine_library.write_text("fake engine library\n")
    processing_package = (
        fake_bin.parent.parent / "lib" / "python3.10" / "site-packages"
        / "unidock_processing"
    )
    processing_package.mkdir(parents=True, exist_ok=True)
    (processing_package / "__init__.py").write_text("# hermetic test package\n")
    conda_meta = fake_bin.parent.parent / "conda-meta"
    conda_meta.mkdir(exist_ok=True)
    for filename in (
        "unidock2-0.6.3-0.json",
        "ambertools_stable-22.5-0.json",
        "pdbfixer-1.12-0.json",
    ):
        (conda_meta / filename).write_text("{}\n")
    ids_file = tmp_path / "official_ids.txt"
    ids_file.write_text("".join(f"{name}\n" for name in ids))
    cfg = {
        "unidock2_bin": str(fake_bin),
        "prep_base": str(prep),
        "converter": "mgl_tools",
        "ids_file": str(ids_file),
        "output_dir": str(tmp_path / "out"),
        "log_dir": str(tmp_path / "logs"),
        "scratch_dir": str(tmp_path / "scratch"),
        "overwrite": False,
        "exhaustiveness": 512,
        "randomize": True,
        "mc_steps": 40,
        "opt_steps": -1,
        "refine_steps": 5,
        "num_pose": 30,
        "energy_range": 6.0,
        "seed": 42,
        "rmsd_limit": 1.0,
        "use_tor_lib": False,
        "energy_decomp": False,
        "construct_ff": False,
        "template_docking": False,
        "compute_center": True,
        "covalent_ligand": False,
        "preserve_receptor_hydrogen": False,
        "engine_checkpoint": False,
        "search_mode": "free",
        "task": "screen",
        "gpu_device_id": 0,
        "gpu_lock_file": str(tmp_path / "gpu0.lock"),
        "gpu_lock_timeout_s": 30,
        "n_cpu": 1,
        "timeout_s": 10,
        "terminate_grace_s": 5.0,
        "max_receptor_atoms": None,
        "exclude": [],
    }
    cfg.update(overrides)
    return cfg


def write_cfg(tmp_path: Path, cfg: dict):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return path


def tool_identity(cfg):
    binary = Path(cfg["unidock2_bin"]).resolve()
    tleap = binary.parent / "tleap"
    teleap = binary.parent / "teLeap"
    ambpdb = binary.parent / "ambpdb"
    env_prefix = binary.parent.parent
    engine_library = next(
        env_prefix.glob(
            "lib/python*/site-packages/unidock_engine/api/python/pipeline*.so"
        )
    ).resolve()
    package_records = [
        next(env_prefix.glob(pattern)).resolve()
        for pattern in (
            "conda-meta/unidock2-*.json",
            "conda-meta/ambertools_stable-*.json",
            "conda-meta/pdbfixer-*.json",
        )
    ]
    processing_package = next(
        env_prefix.glob("lib/python*/site-packages/unidock_processing")
    ).resolve()
    return {
        "path": str(binary),
        "version": "0.6.3",
        "launcher_sha256": ud2._sha256(binary),
        "tleap_path": str(tleap),
        "ambpdb_path": str(ambpdb),
        "artifacts": [
            ud2._file_provenance(path)
            for path in (
                binary, tleap, teleap, ambpdb, engine_library, *package_records
            )
        ],
        "python_source_tree": ud2._python_tree_provenance(processing_package),
    }


def output_arg(cmd):
    return Path(cmd[cmd.index("-o") + 1])


def assert_no_scratch_attempts(cfg):
    scratch_base = Path(cfg["scratch_dir"])
    assert not scratch_base.exists() or not any(scratch_base.iterdir())


def gpu_lock_is_available(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a+", encoding="utf-8") as handle:
        try:
            ud2.fcntl.flock(handle.fileno(), ud2.fcntl.LOCK_EX | ud2.fcntl.LOCK_NB)
        except BlockingIOError:
            return False
        ud2.fcntl.flock(handle.fileno(), ud2.fcntl.LOCK_UN)
        return True


def append_matching_attestation(cmd, log_path: Path, receptor_atoms=1):
    generated = yaml.safe_load(Path(cmd[cmd.index("-cf") + 1]).read_text())
    advanced = generated["Advanced"]
    python_parameters = {
        **advanced,
        **generated["Hardware"],
        **generated["Settings"],
        **generated.get("Preprocessing", {}),
    }
    center_index = cmd.index("-c") + 1
    center = [float(value) for value in cmd[center_index:center_index + 3]]
    size = [float(value) for value in generated["Settings"]["box_size"]]
    temp_dir_name = generated.get("Preprocessing", {}).get("temp_dir_name")
    if temp_dir_name:
        prep_dir = Path(temp_dir_name) / "docking_fake"
        prep_dir.mkdir(parents=True, exist_ok=True)
        (prep_dir / "receptor.prmtop").write_text("fake topology\n")
        (prep_dir / "receptor.inpcrd").write_text("fake coordinates\n")
    bounds = []
    for axis, midpoint, width in zip("xyz", center, size):
        bounds.append(f"{axis}_lo={midpoint - width / 2.0} Angstrom")
        bounds.append(f"{axis}_hi={midpoint + width / 2.0} Angstrom")
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    with open(log_path, "a", encoding="utf-8") as log:
        log.write(
            repr(python_parameters)
            + "\n"
            "DockParam:\n"
            f"seed={advanced['seed']}\n"
            f"exhaustiveness={advanced['exhaustiveness']}\n"
            f"mc_steps={advanced['mc_steps']}\n"
            f"num_pose={advanced['num_pose']}\n"
            f"opt_steps={advanced['opt_steps']}\n"
            f"refine_steps={advanced['refine_steps']}\n"
            f"rmsd_limit={advanced['rmsd_limit']} Angstrom\n"
            + "\n".join(bounds)
            + "\nGroup finished\n"
            f"Receptor has {receptor_atoms} heavy atoms in box\n"
        )


def test_strict_sdf_validation_rejects_corrupt_and_truncated_files(tmp_path):
    valid = sdf_record(-8.25)
    corrupt = {
        "delimiter_only": "$$$$\n",
        "truncated_record": valid.rsplit("$$$$", 1)[0],
        "valid_then_truncated": valid + valid.rsplit("$$$$", 1)[0],
        "missing_m_end": valid.replace("M  END", "M  BAD"),
        "missing_affinity": valid.replace(
            ">  <vina_binding_free_energy>\n-8.25\n\n", ""
        ),
        "nonfinite_affinity": valid.replace("\n-8.25\n", "\nnan\n"),
        "duplicate_affinity": valid.replace(
            "\n$$$$\n",
            "\n>  <vina_binding_free_energy>\n-1.0\n\n$$$$\n",
        ),
    }
    for label, text in corrupt.items():
        path = tmp_path / f"{label}.sdf"
        path.write_text(text)
        assert not ud2._valid_out(path), label
        with pytest.raises(ValueError):
            ud2.parse_sdf_scores(path)


def test_strict_sdf_validation_accepts_complete_scored_records(tmp_path):
    output = tmp_path / "poses.sdf"
    output.write_text(sdf_record(-8.5, "one") + sdf_record(-7.0, "two"))
    assert ud2._valid_out(output)
    assert ud2.parse_sdf_scores(output) == [-8.5, -7.0]
    assert ud2.count_poses(output) == 2
    assert ud2.best_affinity(output) == -8.5


def test_manifest_invalidates_config_inputs_and_box(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)
    calls = []

    def successful_run(cmd, _timeout, log_path, _env, **_kwargs):
        assert not gpu_lock_is_available(cfg["gpu_lock_file"])
        calls.append(list(cmd))
        output_arg(cmd).write_text(sdf_record(-7.0 - len(calls)))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", successful_run)

    first = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert first["status"] == "success"
    assert_no_scratch_attempts(cfg)
    fingerprints = {first["fingerprint"]}
    assert len(calls) == 1

    same = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert same["status"] == "done"
    assert same["fingerprint"] == first["fingerprint"]
    assert len(calls) == 1

    cfg["energy_range"] = 5.0
    changed_config = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed_config["status"] == "success"
    fingerprints.add(changed_config["fingerprint"])
    assert len(calls) == 2

    ligand = Path(ud2.plan_complex(cdir, cfg)[0]["ligand"])
    ligand.write_text(sdf_record(None, title="changed ligand bytes"))
    changed_ligand = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed_ligand["status"] == "success"
    fingerprints.add(changed_ligand["fingerprint"])
    assert len(calls) == 3

    receptor = Path(ud2.plan_complex(cdir, cfg)[0]["receptor"])
    receptor.write_text(pdb_atom() + "\n" + pdb_atom(serial=2, x=1.0) + "\n")
    changed_receptor = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed_receptor["status"] == "success"
    fingerprints.add(changed_receptor["fingerprint"])
    assert len(calls) == 4

    box = Path(ud2.plan_complex(cdir, cfg)[0]["box_file"])
    box.write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
        "size_x = 13\nsize_y = 11\nsize_z = 12\n"
    )
    changed_box = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed_box["status"] == "success"
    fingerprints.add(changed_box["fingerprint"])
    assert len(calls) == 5

    cfg["gpu_lock_file"] = str(tmp_path / "another_gpu0.lock")
    changed_lock = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed_lock["status"] == "success"
    fingerprints.add(changed_lock["fingerprint"])
    assert len(calls) == 6
    manifest = ud2._read_json(
        Path(cfg["output_dir"]) / cdir.name / f"{cdir.name}{ud2._DONE_SUFFIX}"
    )
    assert manifest["driver_execution"]["gpu_lock_file"] == str(
        Path(cfg["gpu_lock_file"]).resolve()
    )
    assert len(fingerprints) == 6


def test_completion_rejects_tampered_provenance_outputs_and_duplicate_attestation(
    tmp_path, monkeypatch
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)

    def successful_run(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", successful_run)
    assert ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    info, ok = ud2.plan_complex(cdir, cfg)
    assert ok
    protocol = ud2.effective_unidock2_config(info, cfg)
    provenance, fingerprint = ud2.build_provenance(
        info, protocol, identity, cfg["gpu_lock_file"]
    )
    out_dir = Path(cfg["output_dir"]) / cdir.name
    output = out_dir / f"{cdir.name}_unidock2_out.sdf"
    marker = out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}"
    original = ud2._read_json(marker)
    assert set(original["runtime_attestation"]) == set(ud2._STABLE_RUNTIME_KEYS)
    assert "temp_dir_name" not in repr(original["runtime_attestation"])

    cases = (
        (("generation_id",), "not-a-generation"),
        (("prepared_receptor_file",), "wrong_prepared.pdb"),
        (("num_poses",), 99),
        (("best_affinity",), -999.0),
        (("tool", "version"), "tampered"),
        (("inputs", "ligand", "sha256"), "0" * 64),
        (("requested_search_settings", "randomize"), False),
        (("runtime_attestation", "n_cpu"), 7),
        (("effective_search_settings", "energy_range"), 999.0),
        (("applied_center",), [1.0, 2.0, 3.0]),
        (("applied_box_size",), [30.0, 30.0, 30.0]),
    )
    for keys, value in cases:
        tampered = ud2.json.loads(ud2.json.dumps(original))
        location = tampered
        for key in keys[:-1]:
            location = location[key]
        location[keys[-1]] = value
        ud2._atomic_write_json(marker, tampered)
        assert not ud2._completion_matches(
            output, marker, fingerprint, provenance, cfg
        ), keys
        ud2._atomic_write_json(marker, original)

    wrong_marker = out_dir / "wrong_unidock2_completion.json"
    ud2._atomic_write_json(wrong_marker, original)
    assert not ud2._completion_matches(output, wrong_marker, fingerprint, provenance, cfg)

    ligand = Path(info["ligand"])
    old_ligand = ligand.read_bytes()
    ligand.write_text(sdf_record(None, title="changed after provenance"))
    assert not ud2._completion_matches(output, marker, fingerprint, provenance, cfg)
    ligand.write_bytes(old_ligand)
    assert ud2._completion_matches(output, marker, fingerprint, provenance, cfg)


def test_gpu_lock_releases_after_success_and_exception(tmp_path):
    lock = tmp_path / "shared_gpu.lock"
    assert gpu_lock_is_available(lock)
    with ud2._gpu_lock(lock):
        assert not gpu_lock_is_available(lock)
    assert gpu_lock_is_available(lock)

    with pytest.raises(RuntimeError, match="engine failed"):
        with ud2._gpu_lock(lock):
            assert not gpu_lock_is_available(lock)
            raise RuntimeError("engine failed")
    assert gpu_lock_is_available(lock)


def test_gpu_lock_rejects_symlink_and_has_bounded_contention_wait(tmp_path):
    if not hasattr(ud2.os, "O_NOFOLLOW"):
        pytest.skip("platform has no O_NOFOLLOW")
    target = tmp_path / "attacker-selected.lock"
    target.write_text("not a lock\n")
    symlink = tmp_path / "gpu.lock"
    symlink.symlink_to(target)
    with pytest.raises(OSError):
        with ud2._gpu_lock(symlink, timeout_s=0.01):
            pass

    real_lock = tmp_path / "real_gpu.lock"
    with ud2._gpu_lock(real_lock, timeout_s=1):
        started = ud2.time.monotonic()
        with pytest.raises(TimeoutError, match="waiting for GPU lock"):
            with ud2._gpu_lock(real_lock, timeout_s=0.02):
                pass
        assert ud2.time.monotonic() - started < 0.5


def test_gpu_lock_queues_through_temporary_contention(tmp_path, monkeypatch):
    lock = tmp_path / "queued_gpu.lock"
    real_flock = ud2.fcntl.flock
    attempts = 0
    waits = []

    def contended_twice(fd, operation):
        nonlocal attempts
        if operation & ud2.fcntl.LOCK_NB:
            attempts += 1
            if attempts <= 2:
                raise BlockingIOError(ud2.errno.EAGAIN, "GPU busy")
        return real_flock(fd, operation)

    monkeypatch.setattr(ud2.fcntl, "flock", contended_twice)
    monkeypatch.setattr(ud2.time, "sleep", lambda seconds: waits.append(seconds))
    with ud2._gpu_lock(lock, gpu_device_id=0, timeout_s=1800):
        assert attempts == 3
    assert len(waits) == 2
    assert all(0 < wait <= 0.1 for wait in waits)


def test_global_gpu_lock_is_shared_across_output_trees(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path, gpu_lock_timeout_s=0.02)
    cdir = make_complex(Path(cfg["prep_base"]))
    second = dict(cfg)
    second["output_dir"] = str(tmp_path / "other-output-tree")
    second["log_dir"] = str(tmp_path / "other-logs")
    second["scratch_dir"] = str(tmp_path / "other-scratch")
    launched = False

    def should_not_launch(*_args, **_kwargs):
        nonlocal launched
        launched = True
        raise AssertionError("contended cross-tree work must not launch")

    monkeypatch.setattr(ud2, "_run", should_not_launch)
    with ud2._gpu_lock(cfg["gpu_lock_file"], timeout_s=1):
        result = ud2.dock_complex(
            cdir, second, second["unidock2_bin"], tool_identity(second)
        )
    assert result["status"] == "locked"
    assert result["gpu_lock_timed_out"] is True
    assert not launched


def test_legacy_sdf_is_not_done_and_failed_partial_is_not_promoted_or_summarized(
    tmp_path, monkeypatch
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)
    out_dir = Path(cfg["output_dir"]) / cdir.name
    out_dir.mkdir(parents=True)
    public = out_dir / f"{cdir.name}_unidock2_out.sdf"
    public.write_text(sdf_record(-99.0, "legacy"))
    calls = []

    def fails_after_writing_valid_partial(cmd, *_args, **_kwargs):
        calls.append(list(cmd))
        output_arg(cmd).write_text(sdf_record(-4.0, "uncommitted"))
        return 17, False

    monkeypatch.setattr(ud2, "_run", fails_after_writing_valid_partial)
    result = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)

    assert len(calls) == 1, "a legacy SDF without a matching manifest must be rerun"
    assert result["status"] == "error"
    assert result.get("num_poses", 0) == 0
    assert result.get("best_affinity") is None
    assert not result.get("published_output_current", False)
    assert not (out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}").exists()
    assert_no_scratch_attempts(cfg)
    if public.exists():
        assert ud2.best_affinity(public) == -99.0


def test_failed_overwrite_invalidates_completion_and_normal_run_retries(
    tmp_path, monkeypatch
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)
    calls = []

    def succeeds(cmd, _timeout, log_path, _env, **_kwargs):
        calls.append(list(cmd))
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", succeeds)
    first = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    out_dir = Path(cfg["output_dir"]) / cdir.name
    marker = out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}"
    assert first["status"] == "success"
    assert marker.exists()

    def fails(cmd, *_args, **_kwargs):
        calls.append(list(cmd))
        output_arg(cmd).write_text(sdf_record(-1.0, "uncommitted overwrite"))
        return 17, False

    cfg["overwrite"] = True
    monkeypatch.setattr(ud2, "_run", fails)
    overwrite = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert overwrite["status"] == "error"
    assert overwrite["previous_generation_preserved"] is False
    assert not marker.exists()
    assert ud2.best_affinity(out_dir / f"{cdir.name}_unidock2_out.sdf") == -8.0

    # Completion authority was invalidated before the failed overwrite, so a
    # normal run must launch again rather than silently accepting stale output.
    cfg["overwrite"] = False
    retry = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert retry["status"] == "error"
    assert len(calls) == 3
    assert not marker.exists()

    # A scientific-setting change likewise remains retryable and cannot revive
    # the invalidated pre-overwrite completion.
    cfg["energy_range"] = 5.0
    changed = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert changed["status"] == "error"
    assert len(calls) == 4


def test_zero_exit_validation_failure_invalidates_then_recovers_on_retry(
    tmp_path, monkeypatch
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)

    def succeeds(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", succeeds)
    assert ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    out_dir = Path(cfg["output_dir"]) / cdir.name
    payload_paths = [
        out_dir / f"{cdir.name}_unidock2_out.sdf",
        out_dir / f"{cdir.name}_unidock2_receptor_prepared.pdb",
        out_dir / f"{cdir.name}_unidock2_config.yaml",
    ]
    marker = out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}"
    committed = {path: path.read_bytes() for path in payload_paths}

    def exits_zero_with_invalid_output(cmd, *_args, **_kwargs):
        output_arg(cmd).write_text("not an SDF\n")
        return 0, False

    cfg["overwrite"] = True
    monkeypatch.setattr(ud2, "_run", exits_zero_with_invalid_output)
    failed = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert failed["status"] == "error"
    assert failed["returncode"] == 0
    assert failed["previous_generation_preserved"] is False
    assert {path: path.read_bytes() for path in payload_paths} == committed
    assert not marker.exists()

    cfg["overwrite"] = False
    monkeypatch.setattr(ud2, "_run", succeeds)
    assert ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    assert marker.exists()


@pytest.mark.parametrize("destination_kind", ["output", "prepared", "config", "manifest"])
def test_atomic_promotion_failure_rolls_back_every_public_artifact(
    tmp_path, monkeypatch, destination_kind
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)

    def succeeds(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0 if not cfg["overwrite"] else -9.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", succeeds)
    assert ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    out_dir = Path(cfg["output_dir"]) / cdir.name
    destinations = {
        "output": out_dir / f"{cdir.name}_unidock2_out.sdf",
        "prepared": out_dir / f"{cdir.name}_unidock2_receptor_prepared.pdb",
        "config": out_dir / f"{cdir.name}_unidock2_config.yaml",
        "manifest": out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}",
    }
    committed_payloads = {
        name: path.read_bytes() for name, path in destinations.items()
        if name != "manifest"
    }
    real_replace = ud2.os.replace
    injected = False

    def fail_one_public_replace(source, destination):
        nonlocal injected
        if Path(destination) == destinations[destination_kind] and not injected:
            injected = True
            raise OSError(f"injected {destination_kind} promotion failure")
        return real_replace(source, destination)

    cfg["overwrite"] = True
    monkeypatch.setattr(ud2.os, "replace", fail_one_public_replace)
    failed = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert injected
    assert failed["status"] == "error"
    assert "generation publication failed" in failed["error"]
    assert failed["previous_generation_preserved"] is False
    assert {
        name: path.read_bytes() for name, path in destinations.items()
        if name != "manifest"
    } == committed_payloads
    assert not destinations["manifest"].exists()

    cfg["overwrite"] = False
    assert ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    assert destinations["manifest"].exists()


def test_generation_promotion_rolls_back_mid_commit_and_recovers(tmp_path, monkeypatch):
    candidates = [tmp_path / f"candidate-{index}" for index in range(4)]
    destinations = [tmp_path / f"public-{index}" for index in range(4)]
    for index, (candidate, destination) in enumerate(zip(candidates, destinations)):
        candidate.write_text(f"new-{index}\n")
        destination.write_text(f"old-{index}\n")

    real_replace = ud2.os.replace
    injected = False

    def fail_second_promotion(source, destination):
        nonlocal injected
        if not injected and Path(source) == candidates[1]:
            injected = True
            raise OSError("injected mid-promotion failure")
        return real_replace(source, destination)

    monkeypatch.setattr(ud2.os, "replace", fail_second_promotion)
    with pytest.raises(OSError, match="injected mid-promotion failure"):
        ud2._promote_generation(list(zip(candidates, destinations)), lambda: True)
    assert [path.read_text() for path in destinations] == [
        f"old-{index}\n" for index in range(4)
    ]

    # A subsequent attempt can publish a complete generation after rollback.
    for index, candidate in enumerate(candidates):
        candidate.write_text(f"retry-{index}\n")
    ud2._promote_generation(
        list(zip(candidates, destinations)),
        lambda: all(
            path.read_text() == f"retry-{index}\n"
            for index, path in enumerate(destinations)
        ),
    )
    assert [path.read_text() for path in destinations] == [
        f"retry-{index}\n" for index in range(4)
    ]


@pytest.mark.parametrize("limit", [0, -1])
def test_nonpositive_limit_is_rejected(tmp_path, monkeypatch, limit):
    cfg = make_cfg(tmp_path)
    make_complex(Path(cfg["prep_base"]))
    config_path = write_cfg(tmp_path, cfg)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))
    with pytest.raises(ValueError, match="limit"):
        ud2.main(config_path, plan_only=True, limit=limit)


def test_plan_excludes_completed_work_and_retains_plan_skips(tmp_path, monkeypatch):
    names = ("5SAK_ZRY", "5SB2_1K2", "5SD5_HWI")
    cfg = make_cfg(tmp_path, ids=names, max_receptor_atoms=1)
    done = make_complex(Path(cfg["prep_base"]), names[0])
    make_complex(Path(cfg["prep_base"]), names[1])
    capped = make_complex(Path(cfg["prep_base"]), names[2])
    receptor = capped / "_staging" / "receptors" / f"{names[2]}_protein.pdb"
    receptor.write_text(pdb_atom() + "\n" + pdb_atom(serial=2, x=1.0) + "\n")
    identity = tool_identity(cfg)

    def successful_run(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", successful_run)
    assert ud2.dock_complex(done, cfg, cfg["unidock2_bin"], identity)["status"] == "success"
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: identity)

    plan = ud2.main(write_cfg(tmp_path, cfg), plan_only=True)
    assert plan["already_done"] == 1
    assert plan["remaining"] == 1
    assert plan["dockable"] == 1
    assert any(item[0] == names[2] for item in plan["skipped"])


def test_main_resume_uses_manifest_result_without_overwriting_rich_success_summary(
    tmp_path, monkeypatch
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    identity = tool_identity(cfg)

    def successful_run(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    monkeypatch.setattr(ud2, "_run", successful_run)
    first = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], identity)
    assert first["status"] == "success"
    summary_path = Path(cfg["output_dir"]) / cdir.name / "docking_summary.json"
    rich_summary = {
        **first,
        "downstream_validation": {"status": "passed", "details": ["keep me"]},
    }
    ud2._atomic_write_json(summary_path, rich_summary)
    before = summary_path.read_bytes()
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: identity)

    results = ud2.main(write_cfg(tmp_path, cfg), plan_only=False)
    resumed = next(result for result in results if result["pdb_id"] == cdir.name)
    assert resumed["status"] == "done"
    assert set(resumed["runtime_attestation"]) == set(ud2._STABLE_RUNTIME_KEYS)
    assert summary_path.read_bytes() == before


@pytest.mark.parametrize("rogue", ["6ZZZ_ABC", "not_an_official_id"])
def test_ids_file_is_authoritative_and_extra_staging_dirs_are_reported_not_selected(
    tmp_path, monkeypatch, rogue
):
    cfg = make_cfg(tmp_path, ids=("5SAK_ZRY",))
    make_complex(Path(cfg["prep_base"]), "5SAK_ZRY")
    make_complex(Path(cfg["prep_base"]), rogue)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))
    plan = ud2.main(write_cfg(tmp_path, cfg), plan_only=True)
    assert plan["dockable"] == 1
    assert plan["unlisted_staged"] == [rogue]


@pytest.mark.parametrize(
    "bad_id",
    ["../5SAK_ZRY", "5SAK", "ABCD_ZRY", "5sak_ZRY", "5SAK_ZRY_extra"],
)
def test_authoritative_ids_require_exact_canonical_shape(tmp_path, bad_id):
    ids_file = tmp_path / "ids.txt"
    ids_file.write_text(bad_id + "\n")
    with pytest.raises(ValueError, match="invalid benchmark ID"):
        ud2.read_authoritative_ids(ids_file)


def test_limit_does_not_reclassify_deferred_official_ids_as_unlisted(tmp_path, monkeypatch):
    names = ("5SAK_ZRY", "5SB2_1K2", "5SD5_HWI", "5SIS_JSM")
    cfg = make_cfg(tmp_path, ids=names)
    for name in names:
        make_complex(Path(cfg["prep_base"]), name)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))
    plan = ud2.main(write_cfg(tmp_path, cfg), plan_only=True, limit=3)
    assert plan["dockable"] == 3
    assert plan["unlisted_staged"] == []
    assert plan["deferred"] == 1


def test_missing_raw_receptor_does_not_fall_back_to_converter_output(tmp_path):
    cfg = make_cfg(tmp_path)
    cdir = Path(cfg["prep_base"]) / "5SAK_ZRY"
    converted = cdir / "_staging" / "receptors" / "pdbqt"
    ligands = cdir / "_staging" / "ligands"
    converted.mkdir(parents=True)
    ligands.mkdir(parents=True)
    (converted / "5SAK_ZRY_protein_mgl_tools.pdb").write_text(pdb_atom() + "\n")
    (converted / "5SAK_ZRY_protein_mgl_tools.box.txt").write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
        "size_x = 10\nsize_y = 10\nsize_z = 10\n"
    )
    (ligands / "5SAK_ZRY_ligand_start_conf.sdf").write_text(sdf_record(None))
    info, ok = ud2.plan_complex(cdir, cfg)
    assert not ok
    assert "raw receptor" in info["plan_error"]


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("exhaustiveness", True),
        ("n_cpu", 1.5),
        ("seed", "42"),
        ("gpu_device_id", False),
        ("overwrite", "false"),
    ],
)
def test_config_rejects_bool_and_lossy_integer_coercions(tmp_path, key, value):
    cfg = make_cfg(tmp_path)
    cfg[key] = value
    with pytest.raises(ValueError, match=key):
        ud2.validate_config(cfg)


def test_config_rejects_unknown_keys_and_invalid_protocol_types(tmp_path):
    cfg = make_cfg(tmp_path)
    cfg["typo_num_poses"] = 30
    with pytest.raises(ValueError, match="unknown configuration keys: typo_num_poses"):
        ud2.validate_config(cfg)

    for key, value in (("randomize", 1), ("opt_steps", -2), ("refine_steps", True)):
        invalid = make_cfg(tmp_path)
        invalid[key] = value
        with pytest.raises(ValueError, match=key):
            ud2.validate_config(invalid)


def test_runtime_estimator_matches_anchors_and_is_monotonic():
    assert ud2._est_seconds(1) == 14.0
    assert ud2._est_seconds(1300) == 14.0
    assert ud2._est_seconds(3700) == 85.0
    assert ud2._est_seconds(41600) == pytest.approx(653.5)
    samples = [ud2._est_seconds(value) for value in (1, 1300, 2000, 3700, 10000, 41600)]
    assert samples == sorted(samples)


def test_effective_config_emits_explicit_cpu_count(tmp_path):
    cfg = make_cfg(tmp_path, n_cpu=3)
    cdir = make_complex(Path(cfg["prep_base"]))
    info, ok = ud2.plan_complex(cdir, cfg)
    assert ok
    generated = ud2.effective_unidock2_config(info, cfg)
    assert generated["Hardware"] == {"gpu_device_id": 0, "n_cpu": 3}


def test_effective_config_pins_scientific_advanced_and_preprocessing_defaults(tmp_path):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    info, ok = ud2.plan_complex(cdir, cfg)
    assert ok
    generated = ud2.effective_unidock2_config(info, cfg)
    assert {
        key: generated["Advanced"][key]
        for key in ("randomize", "opt_steps", "refine_steps", "use_tor_lib", "energy_decomp")
    } == {
        "randomize": True,
        "opt_steps": -1,
        "refine_steps": 5,
        "use_tor_lib": False,
        "energy_decomp": False,
    }
    assert generated["Preprocessing"] == {
        "construct_ff": False,
        "template_docking": False,
        "compute_center": True,
        "covalent_ligand": False,
        "preserve_receptor_hydrogen": False,
        "engine_checkpoint": False,
    }


def test_runtime_attestation_rejects_preset_override_and_default_box(tmp_path):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    info, ok = ud2.plan_complex(cdir, cfg)
    assert ok
    generated_path = tmp_path / "generated.yaml"
    generated_path.write_text(yaml.safe_dump(ud2.effective_unidock2_config(info, cfg)))
    cmd = [
        "fake", "docking", "-c", *[str(value) for value in info["center"]],
        "-cf", str(generated_path), "-o", str(tmp_path / "output.sdf"),
    ]
    log_path = tmp_path / "attestation.log"
    append_matching_attestation(cmd, log_path)
    matching = log_path.read_text()
    assert ud2._parse_runtime_attestation(matching, info, cfg)["box_size"] == info["size"]

    preset = matching.replace("exhaustiveness=512", "exhaustiveness=256").replace(
        "mc_steps=40", "mc_steps=30"
    )
    with pytest.raises(ValueError, match="runtime settings mismatch"):
        ud2._parse_runtime_attestation(preset, info, cfg)

    cpp_override = matching.replace("opt_steps=-1", "opt_steps=17")
    with pytest.raises(ValueError, match="opt_steps"):
        ud2._parse_runtime_attestation(cpp_override, info, cfg)

    default_box = matching
    for axis in "xyz":
        requested = info["size"]["xyz".index(axis)]
        midpoint = info["center"]["xyz".index(axis)]
        default_box = default_box.replace(
            f"{axis}_lo={midpoint - requested / 2.0} Angstrom",
            f"{axis}_lo={midpoint - 15.0} Angstrom",
        ).replace(
            f"{axis}_hi={midpoint + requested / 2.0} Angstrom",
            f"{axis}_hi={midpoint + 15.0} Angstrom",
        )
    with pytest.raises(ValueError, match="runtime settings mismatch"):
        ud2._parse_runtime_attestation(default_box, info, cfg)


def test_prepared_receptor_count_mismatch_cannot_publish(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))

    def mismatched_run(cmd, _timeout, log_path, _env, **_kwargs):
        output_arg(cmd).write_text(sdf_record(-8.0))
        # Fake ambpdb exports one heavy atom, while the engine claims it scored two.
        append_matching_attestation(cmd, log_path, receptor_atoms=2)
        return 0, False

    monkeypatch.setattr(ud2, "_run", mismatched_run)
    result = ud2.dock_complex(
        cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg)
    )
    out_dir = Path(cfg["output_dir"]) / cdir.name
    assert result["status"] == "error"
    assert "heavy-atom count" in result["error"]
    assert result["num_poses"] == 0
    assert result["best_affinity"] is None
    assert not (out_dir / f"{cdir.name}_unidock2_out.sdf").exists()
    assert not (out_dir / f"{cdir.name}_unidock2_receptor_prepared.pdb").exists()
    assert not (out_dir / f"{cdir.name}{ud2._DONE_SUFFIX}").exists()
    assert_no_scratch_attempts(cfg)


def test_malformed_complex_isolated_from_good_complex_in_plan(tmp_path, monkeypatch):
    names = ("5SAK_ZRY", "5SB2_1K2")
    cfg = make_cfg(tmp_path, ids=names)
    make_complex(Path(cfg["prep_base"]), names[0])
    make_complex(Path(cfg["prep_base"]), names[1], malformed_box=True)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))

    plan = ud2.main(write_cfg(tmp_path, cfg), plan_only=True)
    assert plan["dockable"] == 1
    assert any(item[0] == names[1] for item in plan["planning_errors"])


def test_missing_authoritative_complex_is_a_plan_error(tmp_path, monkeypatch):
    names = ("5SAK_ZRY", "5SB2_1K2")
    cfg = make_cfg(tmp_path, ids=names)
    make_complex(Path(cfg["prep_base"]), names[0])
    make_complex(Path(cfg["prep_base"]), names[1], missing=True)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))

    plan = ud2.main(write_cfg(tmp_path, cfg), plan_only=True)
    assert plan["dockable"] == 1
    assert any(item[0] == names[1] for item in plan["planning_errors"])
    assert ud2.results_exit_code(plan) == 1


class WaitRaisesOnce:
    def __init__(self, exception):
        self.pid = 4242
        self.exception = exception
        self.wait_calls = []

    def wait(self, timeout=None):
        self.wait_calls.append(timeout)
        if len(self.wait_calls) == 1:
            raise self.exception
        return -signal.SIGKILL

    def poll(self):
        return None


def install_fake_process(monkeypatch, exception):
    proc = WaitRaisesOnce(exception)
    popen_kwargs = {}
    killed = []
    alive = True

    def fake_popen(_cmd, **kwargs):
        popen_kwargs.update(kwargs)
        return proc

    def fake_killpg(pgid, sig):
        nonlocal alive
        if sig == 0:
            if alive:
                return
            raise ProcessLookupError
        killed.append((pgid, sig))
        # Model a cooperative process group so the driver need not spend its
        # grace period before escalating to SIGKILL in this unit test.
        alive = False

    monkeypatch.setattr(ud2.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(ud2.os, "killpg", fake_killpg)
    return proc, popen_kwargs, killed


def test_timeout_kills_and_reaps_process_group(tmp_path, monkeypatch):
    exc = subprocess.TimeoutExpired(cmd=["fake"], timeout=1)
    proc, popen_kwargs, killed = install_fake_process(monkeypatch, exc)

    assert ud2._run(["fake"], 1, tmp_path / "run.log", {}) == (None, True)
    assert popen_kwargs["start_new_session"] is True
    assert killed == [(proc.pid, signal.SIGTERM)]
    assert len(proc.wait_calls) == 2


def test_keyboard_interrupt_kills_and_reaps_process_group(tmp_path, monkeypatch):
    proc, popen_kwargs, killed = install_fake_process(monkeypatch, KeyboardInterrupt())

    with pytest.raises(KeyboardInterrupt):
        ud2._run(["fake"], 1, tmp_path / "run.log", {})
    assert popen_kwargs["start_new_session"] is True
    assert killed == [(proc.pid, signal.SIGTERM)]
    assert len(proc.wait_calls) == 2


@pytest.mark.parametrize("outcome", ["timeout", "interrupt"])
def test_dock_removes_driver_scratch_after_timeout_or_interrupt(
    tmp_path, monkeypatch, outcome
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))

    def fake_run(*_args, **_kwargs):
        if outcome == "interrupt":
            raise KeyboardInterrupt
        return None, True

    monkeypatch.setattr(ud2, "_run", fake_run)
    if outcome == "interrupt":
        with pytest.raises(KeyboardInterrupt):
            ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg))
    else:
        result = ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg))
        assert result["status"] == "timeout"
    assert gpu_lock_is_available(cfg["gpu_lock_file"])
    assert_no_scratch_attempts(cfg)


@pytest.mark.parametrize("outcome", ["success", "failure", "interrupt"])
def test_scratch_cleanup_failure_warns_without_changing_primary_outcome(
    tmp_path, monkeypatch, outcome
):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))

    def fake_run(cmd, _timeout, log_path, _env, **_kwargs):
        if outcome == "interrupt":
            raise KeyboardInterrupt("primary interrupt")
        if outcome == "failure":
            return 23, False
        output_arg(cmd).write_text(sdf_record(-8.0))
        append_matching_attestation(cmd, log_path)
        return 0, False

    def cleanup_fails(_path):
        raise OSError("scratch is busy")

    monkeypatch.setattr(ud2, "_run", fake_run)
    monkeypatch.setattr(ud2.shutil, "rmtree", cleanup_fails)
    monkeypatch.setattr(ud2.time, "sleep", lambda _seconds: None)

    if outcome == "interrupt":
        with pytest.warns(RuntimeWarning, match="scratch is busy"):
            with pytest.raises(KeyboardInterrupt, match="primary interrupt"):
                ud2.dock_complex(cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg))
        return

    with pytest.warns(RuntimeWarning, match="scratch is busy"):
        result = ud2.dock_complex(
            cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg)
        )
    assert result["status"] == outcome.replace("failure", "error")
    assert "scratch is busy" in result["scratch_cleanup_warning"]
    if outcome == "success":
        assert Path(result["completion_manifest"]).is_file()
        assert result["published_output_current"] is True
    else:
        assert result["returncode"] == 23


def test_lock_contention_is_fail_closed(tmp_path):
    lock_path = tmp_path / "complex.lock"
    with ud2._complex_lock(lock_path) as acquired:
        assert acquired
        with ud2._complex_lock(lock_path) as second:
            assert not second


def test_dock_reports_lock_contention_without_launching(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir = make_complex(Path(cfg["prep_base"]))
    launched = False

    def should_not_run(*_args, **_kwargs):
        nonlocal launched
        launched = True
        raise AssertionError("contended work must not launch")

    @contextmanager
    def contended(_path):
        yield False

    monkeypatch.setattr(ud2, "_run", should_not_run)
    monkeypatch.setattr(ud2, "_complex_lock", contended)
    result = ud2.dock_complex(
        cdir, cfg, cfg["unidock2_bin"], tool_identity(cfg)
    )
    assert result["status"] == "locked"
    assert not launched


def test_batch_lock_blocks_all_selected_ids_and_sets_failure_status(tmp_path, monkeypatch):
    names = ("5SAK_ZRY", "5SB2_1K2")
    cfg = make_cfg(tmp_path, ids=names)
    for name in names:
        make_complex(Path(cfg["prep_base"]), name)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))
    launched = []

    def should_not_dock(cdir, *_args, **_kwargs):
        launched.append(Path(cdir).name)
        raise AssertionError("a contended batch must not start any complex")

    monkeypatch.setattr(ud2, "dock_complex", should_not_dock)
    batch_lock = Path(cfg["output_dir"]) / ".unidock2.batch.lock"
    with ud2._complex_lock(batch_lock) as acquired:
        assert acquired
        results = ud2.main(write_cfg(tmp_path, cfg), plan_only=False)

    assert launched == []
    assert any(item["pdb_id"] == "__batch__" and item["status"] == "locked"
               for item in results)
    assert ud2.results_exit_code(results) == 1


def test_batch_pays_global_gpu_lock_timeout_once_then_stops_launching(
    tmp_path, monkeypatch
):
    names = ("5SAK_ZRY", "5SB2_1K2", "5SD5_HWI")
    cfg = make_cfg(tmp_path, ids=names, gpu_lock_timeout_s=0.01)
    for name in names:
        make_complex(Path(cfg["prep_base"]), name)
    monkeypatch.setattr(ud2, "probe_tool_identity", lambda _binary: tool_identity(cfg))
    launched = []

    def first_lock_timeout(cdir, *_args, **_kwargs):
        launched.append(Path(cdir).name)
        return {
            "pdb_id": Path(cdir).name,
            "status": "locked",
            "gpu_lock_timed_out": True,
            "error": "timed out waiting for shared GPU",
            "num_poses": 0,
            "best_affinity": None,
        }

    monkeypatch.setattr(ud2, "dock_complex", first_lock_timeout)
    results = ud2.main(write_cfg(tmp_path, cfg), plan_only=False)

    assert launched == [names[0]]
    by_id = {result["pdb_id"]: result for result in results}
    assert by_id[names[0]]["gpu_lock_timed_out"] is True
    for name in names[1:]:
        assert by_id[name]["status"] == "locked"
        assert by_id[name]["gpu_lock_wait_skipped"] is True
    assert ud2.results_exit_code(results) == 1


def test_cli_result_status_failure_semantics():
    assert ud2.results_exit_code([{"status": "success"}, {"status": "done"}]) == 0
    for status in ("timeout", "error", "locked"):
        assert ud2.results_exit_code([{"status": "success"}, {"status": status}]) == 1


def test_probe_identity_deduplicates_symlinked_python_library_alias(tmp_path):
    env = tmp_path / "env"
    binary = env / "bin" / "unidock2"
    binary.parent.mkdir(parents=True)
    binary.write_text("#!/bin/sh\necho 'unidock2 0.6.3'\n")
    binary.chmod(0o755)
    for executable in ("tleap", "teLeap", "ambpdb"):
        path = binary.parent / executable
        path.write_text("#!/bin/sh\nexit 0\n")
        path.chmod(0o755)

    engine = (
        env / "lib" / "python3.10" / "site-packages" / "unidock_engine"
        / "api" / "python" / "pipeline.so"
    )
    engine.parent.mkdir(parents=True)
    engine.write_bytes(b"engine")
    processing_package = (
        env / "lib" / "python3.10" / "site-packages" / "unidock_processing"
    )
    processing_package.mkdir(parents=True)
    (processing_package / "__init__.py").write_text("# hermetic test package\n")
    # Mirrors the production conda environment: glob sees two spellings, but
    # both resolve to the same inode and must represent one identity artifact.
    (env / "lib" / "python3.1").symlink_to("python3.10")
    records = env / "conda-meta"
    records.mkdir()
    for filename in (
        "unidock2-0.6.3-0.json",
        "ambertools_stable-23.0-0.json",
        "pdbfixer-1.9-0.json",
    ):
        (records / filename).write_text("{}\n")

    identity = ud2.probe_tool_identity(str(binary))
    artifact_paths = [item["path"] for item in identity["artifacts"]]
    assert artifact_paths.count(str(engine.resolve())) == 1
    assert identity["python_source_tree"]["file_count"] == 1
    source_hash = identity["python_source_tree"]["sha256"]
    (processing_package / "worker.py").write_text("VALUE = 2\n")
    changed = ud2.probe_tool_identity(str(binary))
    assert changed["python_source_tree"]["file_count"] == 2
    assert changed["python_source_tree"]["sha256"] != source_hash


def test_cli_returns_nonzero_for_plan_error_and_failed_fake_engine(tmp_path):
    names = ("5SAK_ZRY", "5SB2_1K2")
    cfg = make_cfg(tmp_path, ids=names)
    make_complex(Path(cfg["prep_base"]), names[0])
    make_complex(Path(cfg["prep_base"]), names[1], missing=True)
    config_path = write_cfg(tmp_path, cfg)
    base_cmd = [sys.executable, str(Path(ud2.__file__).resolve()), "-c", str(config_path)]

    planned = subprocess.run(
        [*base_cmd, "--plan-only"], capture_output=True, text=True, timeout=10
    )
    assert planned.returncode == 1
    assert names[1] in planned.stdout

    # Restrict the authoritative set to the complete fixture. The executable is
    # a harmless stub that exits zero without writing an SDF, which must still
    # make an execution-mode CLI invocation fail.
    Path(cfg["ids_file"]).write_text(names[0] + "\n")
    executed = subprocess.run(base_cmd, capture_output=True, text=True, timeout=10)
    assert executed.returncode == 1
    summary = ud2._read_json(
        Path(cfg["output_dir"]) / names[0] / "docking_summary.json"
    )
    assert summary["status"] == "error"
    assert summary["num_poses"] == 0
