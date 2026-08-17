from __future__ import annotations

import json
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from pathlib import Path

import pandas as pd
import pytest
import yaml
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from Scripts.Docking import run_autodock as ra


def _pdbqt_atom(serial: int, symbol: str, xyz: tuple[float, float, float], atom_type: str) -> str:
    x, y, z = xyz
    return (
        f"ATOM  {serial:5d} {symbol:>2s}   UNL     1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 {atom_type:>2s}\n"
    )


def _ethanol_pdbqt(*, wrapped: bool = False, complete: bool = True) -> str:
    body = (
        "REMARK VINA RESULT: -7.500 0.000 0.000\n"
        "REMARK SMILES CCO\n"
        "REMARK SMILES IDX 1 1 2 2 3 3\n"
        "ROOT\n"
        + _pdbqt_atom(1, "C", (0.0, 0.0, 0.1), "C")
        + _pdbqt_atom(2, "C", (1.5, 0.0, 0.1), "C")
        + _pdbqt_atom(3, "O", (3.0, 0.0, 0.1), "OA")
        + "ENDROOT\nTORSDOF 0\n"
    )
    if not wrapped:
        return body
    return "MODEL 1\n" + body + ("ENDMDL\n" if complete else "")


def _write_sdf(
    path: Path,
    smiles: str = "CCO",
    *,
    scores: dict[str, float] | None = None,
    renumber: bool = False,
) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for idx in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(idx, (idx * 1.5, 0.0, 0.1))
    mol.AddConformer(conf)
    if renumber:
        order = list(reversed(range(mol.GetNumAtoms())))
        mol = Chem.RenumberAtoms(mol, order)
    for key, value in (scores or {}).items():
        mol.SetProp(key, str(value))
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return mol


def _gnina_cfg() -> dict:
    return {
        "optimization": "gnina",
        "scoring_function": "vina",
        "optimize_search": "minimize",
        "optimize_scoring": "default",
        "optimize_rank_by": "cnn_affinity",
        "optimize_autobox_add": 4.0,
        "optimize_cpu": 1,
        "optimize_seed": 0,
        "optimize_timeout": 30,
        "gnina_cnn_scoring": "rescore",
        "gnina_cnn_model": "crossdock_default2018_ensemble",
        "gnina_cnn_model_file": "",
        "gnina_use_gpu": False,
    }


def test_split_rejects_trailing_incomplete_model_without_partial_files(tmp_path: Path) -> None:
    pose = tmp_path / "poses.pdbqt"
    pose.write_text(_ethanol_pdbqt(wrapped=True) + "MODEL 2\n" + _ethanol_pdbqt())
    split_dir = tmp_path / "split"

    with pytest.raises(ValueError, match="truncated"):
        ra._split_pdbqt_models(pose, split_dir)

    assert not list(split_dir.glob("*.pdbqt"))


@pytest.mark.parametrize("content", [
    _ethanol_pdbqt(wrapped=True, complete=False),
    _ethanol_pdbqt(wrapped=True).replace("REMARK VINA RESULT:", "REMARK NO SCORE:"),
    _ethanol_pdbqt(wrapped=True) + "trailing junk\n",
    _ethanol_pdbqt(wrapped=False),
])
def test_public_vina_validator_rejects_partial_scoreless_junk_and_unwrapped(
    tmp_path: Path, content: str,
) -> None:
    output = tmp_path / "bad_vina_out.pdbqt"
    output.write_text(content)

    with pytest.raises(ValueError):
        ra.validate_vina_output(output)


def test_public_vina_validator_returns_authoritative_pose_scores(tmp_path: Path) -> None:
    output = tmp_path / "good_vina_out.pdbqt"
    output.write_text(_ethanol_pdbqt(wrapped=True))

    assert ra.validate_vina_output(output) == [{
        "mode": 1, "affinity": -7.5, "rmsd_lb": 0.0, "rmsd_ub": 0.0,
    }]


def _raw_docking_inputs(tmp_path: Path) -> tuple[dict, dict, Path, Path]:
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("ATOM\n")
    box = tmp_path / "protein.box.txt"
    box.write_text("center_x=0\ncenter_y=0\ncenter_z=0\nsize_x=20\nsize_y=20\nsize_z=20\n")
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text(_ethanol_pdbqt())
    dock_dir = tmp_path / "docking"
    log_dir = tmp_path / "logs"
    dock_dir.mkdir()
    log_dir.mkdir()
    return (
        {"pdbqt": receptor, "box": box}, {"pdbqt": ligand}, dock_dir, log_dir,
    )


def _raw_job_kwargs(dock_dir: Path, log_dir: Path) -> dict:
    return {
        "scoring": "vina", "vina_bin": Path("/fake/vina"),
        "dock_dir": dock_dir, "log_dir": log_dir, "exhaustiveness": 1,
        "num_modes": 1, "energy_range": 3, "cpu": 1, "seed": 0,
        "timeout": 10, "overwrite": False, "overwrite_error_log": False,
        "error_log": {}, "failed_job_policy": "retry",
    }


def test_single_resume_quarantines_invalid_cache_retries_failure_and_commits_atomically(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    rec, lig, dock_dir, log_dir = _raw_docking_inputs(tmp_path)
    final = dock_dir / "protein__ligand_vina_vina_out.pdbqt"
    final.write_text("partial junk")
    kwargs = _raw_job_kwargs(dock_dir, log_dir)
    kwargs["error_log"] = {"ligand__protein": {"error": "old timeout"}}
    observed: dict[str, Path | bool] = {}

    class FakePopen:
        def __init__(self, command, **_kwargs):
            tmp_output = Path(command[command.index("--out") + 1])
            observed["tmp"] = tmp_output
            observed["final_absent_during_run"] = not final.exists()
            tmp_output.write_text(_ethanol_pdbqt(wrapped=True))

        def wait(self, timeout=None):
            return 0

    monkeypatch.setattr(ra.subprocess, "Popen", FakePopen)
    job = ra._dock_one_job(rec, lig, **kwargs)

    assert job["result"].status == "success" and job["ran_vina"]
    assert observed["tmp"] != final and observed["final_absent_during_run"] is True
    assert len(ra.validate_vina_output(final)) == 1
    quarantined = list((dock_dir / "_invalid_vina_outputs").glob("*.invalid"))
    assert len(quarantined) == 1 and quarantined[0].read_text() == "partial junk"


def test_single_timeout_never_promotes_partial_output(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rec, lig, dock_dir, log_dir = _raw_docking_inputs(tmp_path)
    final = dock_dir / "protein__ligand_vina_vina_out.pdbqt"

    class TimeoutPopen:
        def __init__(self, command, **_kwargs):
            self.output = Path(command[command.index("--out") + 1])
            self.output.write_text(_ethanol_pdbqt(wrapped=True, complete=False))
            self.first_wait = True

        def wait(self, timeout=None):
            if self.first_wait:
                self.first_wait = False
                raise subprocess.TimeoutExpired("vina", timeout)
            return -9

        def kill(self):
            return None

    monkeypatch.setattr(ra.subprocess, "Popen", TimeoutPopen)
    job = ra._dock_one_job(rec, lig, **_raw_job_kwargs(dock_dir, log_dir))

    assert job["result"].status == "failed"
    assert not final.exists()
    assert len(list((dock_dir / "_invalid_vina_outputs").glob("*.invalid"))) == 1


def test_failed_overwrite_archives_old_commit_then_normal_resume_retries(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    rec, lig, dock_dir, log_dir = _raw_docking_inputs(tmp_path)
    final = dock_dir / "protein__ligand_vina_vina_out.pdbqt"
    final.write_text(_ethanol_pdbqt(wrapped=True))
    launches = 0

    class FailThenSucceedPopen:
        def __init__(self, command, **_kwargs):
            nonlocal launches
            launches += 1
            self.returncode = 2 if launches == 1 else 0
            tmp_output = Path(command[command.index("--out") + 1])
            tmp_output.write_text(
                _ethanol_pdbqt(wrapped=True, complete=self.returncode == 0)
            )

        def wait(self, timeout=None):
            return self.returncode

    monkeypatch.setattr(ra.subprocess, "Popen", FailThenSucceedPopen)
    overwrite_kwargs = _raw_job_kwargs(dock_dir, log_dir)
    overwrite_kwargs["overwrite"] = True
    first = ra._dock_one_job(rec, lig, **overwrite_kwargs)

    assert first["result"].status == "failed" and not final.exists()
    archived = list((dock_dir / "_superseded_vina_outputs").glob("*.superseded"))
    assert len(archived) == 1 and len(ra.validate_vina_output(archived[0])) == 1

    retry_kwargs = _raw_job_kwargs(dock_dir, log_dir)
    retry_kwargs["error_log"] = {"ligand__protein": {"error": "failed overwrite"}}
    second = ra._dock_one_job(rec, lig, **retry_kwargs)

    assert launches == 2 and second["result"].status == "success"
    assert len(ra.validate_vina_output(final)) == 1


def test_batch_resume_validates_cache_and_retries_instead_of_trusting_nonempty(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    rec, lig, dock_dir, log_dir = _raw_docking_inputs(tmp_path)
    final = dock_dir / "protein__ligand_vina_vina_out.pdbqt"
    final.write_text("not a pose")

    class FakeBatchPopen:
        def __init__(self, command, **_kwargs):
            batch_dir = Path(command[command.index("--dir") + 1])
            ligand_path = Path(command[command.index("--batch") + 1])
            (batch_dir / f"{ligand_path.stem}_out.pdbqt").write_text(
                _ethanol_pdbqt(wrapped=True)
            )

        def wait(self, timeout=None):
            return 0

    monkeypatch.setattr(ra.subprocess, "Popen", FakeBatchPopen)
    kwargs = _raw_job_kwargs(dock_dir, log_dir)
    kwargs["error_log"] = {"ligand__protein": {"error": "old failure"}}
    jobs = ra._dock_batch_job(rec, [lig], **kwargs, batch_timeout=10)

    assert len(jobs) == 1 and jobs[0]["result"].status == "success"
    assert len(ra.validate_vina_output(final)) == 1
    assert len(list((dock_dir / "_invalid_vina_outputs").glob("*.invalid"))) == 1


def test_explicit_skip_policy_keeps_prior_failure_failed_not_skipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    rec, lig, dock_dir, log_dir = _raw_docking_inputs(tmp_path)
    kwargs = _raw_job_kwargs(dock_dir, log_dir)
    kwargs.update({
        "failed_job_policy": "skip",
        "error_log": {"ligand__protein": {"error": "deliberately deferred"}},
    })
    monkeypatch.setattr(
        ra.subprocess, "Popen",
        lambda *args, **kw: pytest.fail("policy=skip must not invoke Vina"),
    )

    job = ra._dock_one_job(rec, lig, **kwargs)
    assert job["result"].status == "failed"
    assert job["raw"]["status"] == "cached-failure"
    assert not job["ran_vina"]


def test_meeko_05_uses_output_flag_o_and_never_obabel_for_remark_smiles(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "pose.pdbqt"
    source.write_text(_ethanol_pdbqt())
    output = tmp_path / "pose.sdf"
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append([str(part) for part in command])
        if "--help" in command:
            return subprocess.CompletedProcess(command, 0, stdout="--output_filename", stderr="")
        assert command[0] == "/fake/mk_export.py"
        assert "-o" in command and "-s" not in command
        _write_sdf(Path(command[command.index("-o") + 1]))
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    ra._MK_EXPORT_OUTPUT_FLAG_CACHE.clear()
    ok, converter = ra._convert_pose_to_sdf(
        source, output,
        {"mk_export": "/fake/mk_export.py", "obabel": "/fake/obabel"},
    )

    assert ok and converter == "mk_export"
    assert all(call[0] != "/fake/obabel" for call in calls)


def test_installed_meeko_cli_output_flags_are_version_correct() -> None:
    meeko_05 = Path("/home/manndo/anaconda3/bin/mk_export.py")
    meeko_07 = Path("/home/manndo/anaconda3/envs/vina/bin/mk_export.py")
    if not (meeko_05.exists() and meeko_07.exists()):
        pytest.skip("both benchmark Meeko installations are not available")
    ra._MK_EXPORT_OUTPUT_FLAG_CACHE.clear()
    assert ra._mk_export_output_flag(str(meeko_05)) == "-o"
    assert ra._mk_export_output_flag(str(meeko_07)) == "-s"
    assert ra._converter_version("mk_export", str(meeko_07)) == "meeko 0.7.1"


def test_real_7uaw_meeko07_reconstruction_preserves_exact_ligand(tmp_path: Path) -> None:
    pose = Path(
        "Dockings/Benchmark/7UAW_MF6/meeko/docking/"
        "7UAW_MF6_protein_meeko__7UAW_MF6_ligand_start_conf_vina_vina_out.pdbqt"
    )
    meeko = Path("/home/manndo/anaconda3/envs/vina/bin/mk_export.py")
    if not (pose.exists() and meeko.exists()):
        pytest.skip("7UAW benchmark output or pinned Meeko 0.7.1 is unavailable")
    rank, model, _ = ra._split_pdbqt_models(pose, tmp_path / "split")[0]
    output = tmp_path / "7uaw_rank1.sdf"

    ok, converter = ra._convert_pose_to_sdf(
        model, output, {"mk_export": str(meeko), "obabel": "/usr/bin/obabel"},
    )
    valid, message, _ = ra._validate_reconstructed_sdf(model, output)
    mol = next(mol for mol in Chem.SDMolSupplier(str(output), removeHs=False) if mol)

    assert rank == 1 and ok and converter == "mk_export"
    assert valid, message
    assert mol.GetNumAtoms() == 56 and mol.GetNumHeavyAtoms() == 35
    assert len(Chem.GetMolFrags(Chem.RemoveHs(mol))) == 1


def test_explicit_smiles_hydrogen_mapping_is_complete_and_element_checked(tmp_path: Path) -> None:
    smiles = "[2H]N"
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None and mol.GetNumAtoms() == 2 and mol.GetNumHeavyAtoms() == 1
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    coordinates = []
    for idx in range(mol.GetNumAtoms()):
        xyz = (idx * 1.25, idx * 0.1, 0.2)
        coordinates.append(xyz)
        conf.SetAtomPosition(idx, xyz)
    mol.AddConformer(conf)
    sdf = tmp_path / "explicit_h.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(mol)
    writer.close()

    mapping = " ".join(f"{idx} {idx}" for idx in range(1, mol.GetNumAtoms() + 1))
    atom_types = {"H": "HD", "N": "N", "C": "C"}
    atoms = "".join(
        _pdbqt_atom(idx + 1, atom.GetSymbol(), coordinates[idx], atom_types[atom.GetSymbol()])
        for idx, atom in enumerate(mol.GetAtoms())
    )
    source = tmp_path / "explicit_h.pdbqt"
    source.write_text(f"REMARK SMILES {smiles}\nREMARK SMILES IDX {mapping}\n{atoms}")

    valid, message, identity = ra._validate_reconstructed_sdf(source, sdf)
    assert valid, message
    assert identity is not None and len(ra._remark_smiles_mapping(source.read_text())) == 2

    # Retain strictness: the same complete coordinate map must also preserve
    # element identity, including the explicitly represented hydrogen.
    source.write_text(source.read_text().replace(" 0.000 HD\n", " 0.000 C \n", 1))
    valid, message, _ = ra._validate_reconstructed_sdf(source, sdf)
    assert not valid and "element mapping mismatch" in message


def test_real_5sak_explicit_h_reconstruction_preserves_full_mapping(tmp_path: Path) -> None:
    pose = Path(
        "Dockings/vina_results_full_protein_vina_scoring/5SAK_ZRY/mgl_tools/docking/"
        "5SAK_ZRY_protein_mgl_tools__5SAK_ZRY_ligand_start_conf_vina_vina_out.pdbqt"
    )
    meeko = Path("/home/manndo/anaconda3/envs/vina/bin/mk_export.py")
    if not (pose.exists() and meeko.exists()):
        pytest.skip("5SAK production pose or pinned Meeko 0.7.1 is unavailable")
    rank, model, _ = ra._split_pdbqt_models(pose, tmp_path / "split_5sak")[0]
    output = tmp_path / "5sak_rank1.sdf"

    ok, converter = ra._convert_pose_to_sdf(
        model, output, {"mk_export": str(meeko), "obabel": ""},
    )
    valid, message, _ = ra._validate_reconstructed_sdf(model, output)
    mol = next(mol for mol in Chem.SDMolSupplier(str(output), removeHs=False) if mol)
    expected = Chem.MolFromSmiles(ra._remark_smiles(model.read_text()))

    assert rank == 1 and ok and converter == "mk_export"
    assert valid, message
    assert expected is not None and expected.GetNumAtoms() == 19
    assert expected.GetNumHeavyAtoms() == 18
    assert len(ra._remark_smiles_mapping(model.read_text())) == expected.GetNumAtoms()
    assert mol.GetNumAtoms() == 30 and mol.GetNumHeavyAtoms() == 18


def test_remark_smiles_failure_and_macrocycle_never_fall_back_to_obabel(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "remark_pose.pdbqt"
    source.write_text(_ethanol_pdbqt())
    calls: list[str] = []

    def fake_run(command, **kwargs):
        calls.append(str(command[0]))
        if "--help" in command:
            return subprocess.CompletedProcess(command, 0, stdout="--output_filename", stderr="")
        return subprocess.CompletedProcess(command, 2, stdout="", stderr="failed")

    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    ra._MK_EXPORT_OUTPUT_FLAG_CACHE.clear()
    ok, message = ra._convert_pose_to_sdf(
        source, tmp_path / "remark.sdf",
        {"mk_export": "/fake/mk", "obabel": "/fake/obabel"},
    )
    assert not ok and "fallback disabled" in message
    assert "/fake/obabel" not in calls

    macrocycle = tmp_path / "macrocycle.pdbqt"
    macrocycle.write_text(
        "REMARK VINA RESULT: -1.0 0 0\nROOT\n"
        + _pdbqt_atom(1, "C", (0.0, 0.0, 0.0), "CG0")
        + "ENDROOT\n"
    )
    ok, message = ra._convert_pose_to_sdf(
        macrocycle, tmp_path / "macrocycle.sdf",
        {"mk_export": "/fake/mk", "obabel": "/fake/obabel"},
    )
    assert not ok and "macrocycle glue" in message
    assert "/fake/obabel" not in calls


def test_optimizer_is_atomic_pins_scoring_and_accepts_valid_atom_reordering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    ligand = tmp_path / "input.sdf"
    input_mol = _write_sdf(ligand)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n")
    output = tmp_path / "optimized.sdf"
    output.write_text("stale output")
    ra._provenance_path(output).write_text("stale sidecar")
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append([str(part) for part in command])
        assert not output.exists(), "stale final output must be removed before execution"
        tmp_output = Path(command[command.index("--out") + 1])
        assert tmp_output != output
        _write_sdf(
            tmp_output,
            scores={
                "minimizedAffinity": -8.0,
                "CNNscore": 0.7,
                "CNNaffinity": 6.1,
            },
            renumber=True,
        )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    monkeypatch.setattr(
        ra, "_interprocess_gpu_lock",
        lambda *_args, **_kwargs: pytest.fail("CPU gnina must not acquire the GPU lock"),
    )
    ok, scores, _ = ra.run_optimizer_tool(
        "gnina", "/fake/gnina", receptor, ligand, output, _gnina_cfg(),
        expected_identity=ra._molecule_identity(input_mol),
    )

    assert ok and output.exists()
    assert scores == {
        "minimized_affinity": -8.0, "cnn_score": 0.7, "cnn_affinity": 6.1,
    }
    command = commands[0]
    assert command[command.index("--scoring") + 1] == "default"
    assert command[command.index("--cnn_scoring") + 1] == "rescore"
    assert command[command.index("--cnn") + 1] == "crossdock_default2018_ensemble"


def test_gpu_gnina_requires_and_holds_configured_lock_around_subprocess(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    ligand = tmp_path / "input.sdf"
    input_mol = _write_sdf(ligand)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n")
    output = tmp_path / "optimized.sdf"
    cfg = _gnina_cfg()
    cfg["gnina_use_gpu"] = True

    with pytest.raises(ValueError, match="gpu_lock_file"):
        ra._validate_optimizer_config(cfg)

    cfg["gpu_lock_file"] = str(tmp_path / "gpu.lock")
    state = {"held": False, "path": None}

    @contextmanager
    def fake_lock(path):
        state["held"] = True
        state["path"] = str(path)
        try:
            yield
        finally:
            state["held"] = False

    def fake_run(command, **kwargs):
        assert state["held"], "GPU lock must cover the actual gnina subprocess"
        assert "--no_gpu" not in command
        _write_sdf(
            Path(command[command.index("--out") + 1]),
            scores={
                "minimizedAffinity": -8.0, "CNNscore": 0.7, "CNNaffinity": 6.1,
            },
        )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(ra, "_interprocess_gpu_lock", fake_lock)
    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    ok, _, _ = ra.run_optimizer_tool(
        "gnina", "/fake/gnina", receptor, ligand, output, cfg,
        expected_identity=ra._molecule_identity(input_mol),
    )

    assert ok and state["path"] == cfg["gpu_lock_file"] and not state["held"]


def test_gpu_flock_serializes_independent_processes(tmp_path: Path) -> None:
    lock_path = tmp_path / "shared-gpu.lock"
    code = """
import sys
import time
from Scripts.Docking.run_autodock import _interprocess_gpu_lock
with _interprocess_gpu_lock(sys.argv[1]):
    print(time.monotonic(), flush=True)
    time.sleep(0.20)
    print(time.monotonic(), flush=True)
"""
    procs = [
        subprocess.Popen(
            [sys.executable, "-c", code, str(lock_path)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        for _ in range(2)
    ]
    intervals = []
    for proc in procs:
        stdout, stderr = proc.communicate(timeout=20)
        assert proc.returncode == 0, stderr
        start, end = (float(value) for value in stdout.splitlines())
        intervals.append((start, end))

    intervals.sort()
    assert intervals[1][0] >= intervals[0][1] - 0.01


def test_optimizer_rejects_missing_score_tags_and_does_not_promote_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    ligand = tmp_path / "input.sdf"
    input_mol = _write_sdf(ligand)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n")
    output = tmp_path / "optimized.sdf"

    def fake_run(command, **kwargs):
        _write_sdf(
            Path(command[command.index("--out") + 1]),
            scores={"minimizedAffinity": -8.0},
        )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    ok, _, message = ra.run_optimizer_tool(
        "gnina", "/fake/gnina", receptor, ligand, output, _gnina_cfg(),
        expected_identity=ra._molecule_identity(input_mol),
    )
    assert not ok and "missing required score" in message
    assert not output.exists()


def test_legacy_or_corrupt_cache_is_never_reused(tmp_path: Path) -> None:
    output = tmp_path / "optimized.sdf"
    input_mol = _write_sdf(
        output,
        scores={
            "minimizedAffinity": -8.0, "CNNscore": 0.7, "CNNaffinity": 6.1,
        },
    )
    identity = ra._molecule_identity(input_mol)
    expected = {"fingerprint": "expected"}

    ok, _, message, _ = ra._validate_cached_optimizer_output(
        output, expected, identity,
        {"minimized_affinity", "cnn_score", "cnn_affinity"},
    )
    assert not ok and "provenance" in message

    sidecar = ra._provenance_path(output)
    sidecar.write_text(json.dumps({
        "schema_version": ra.OPTIMIZER_CACHE_SCHEMA,
        "fingerprint": "wrong",
        "output_sha256": ra._sha256_file(output),
        "optimizer_elapsed_time_s": 3.2,
    }))
    ok, _, message, _ = ra._validate_cached_optimizer_output(
        output, expected, identity,
        {"minimized_affinity", "cnn_score", "cnn_affinity"},
    )
    assert not ok and "fingerprint mismatch" in message


def test_valid_cache_reuse_preserves_execution_time_and_detects_tampering(tmp_path: Path) -> None:
    output = tmp_path / "optimized.sdf"
    input_mol = _write_sdf(
        output,
        scores={
            "minimizedAffinity": -8.0, "CNNscore": 0.7, "CNNaffinity": 6.1,
        },
    )
    identity = ra._molecule_identity(input_mol)
    provenance = {
        "schema_version": ra.OPTIMIZER_CACHE_SCHEMA,
        "fingerprint": "matching-fingerprint",
    }
    ra._write_optimizer_provenance(output, provenance, 4.25)

    ok, scores, message, elapsed = ra._validate_cached_optimizer_output(
        output, provenance, identity,
        {"minimized_affinity", "cnn_score", "cnn_affinity"},
    )
    assert ok and message == "reused validated cache" and elapsed == 4.25
    assert scores["cnn_affinity"] == 6.1

    output.write_text(output.read_text() + "\n")
    ok, _, message, _ = ra._validate_cached_optimizer_output(
        output, provenance, identity,
        {"minimized_affinity", "cnn_score", "cnn_affinity"},
    )
    assert not ok and "checksum mismatch" in message


def test_log_scope_replacement_removes_obsolete_ranks_and_is_thread_safe(tmp_path: Path) -> None:
    log = tmp_path / "optimization_log.csv"
    old = pd.DataFrame([
        {"combo_name": "A", "tool": "gnina", "autodock_rank": rank}
        for rank in (1, 2, 3)
    ] + [{"combo_name": "KEEP", "tool": "gnina", "autodock_rank": 1}])
    old.to_csv(log, index=False)

    ra._write_optimization_log(
        log,
        [{"combo_name": "A", "tool": "gnina", "autodock_rank": 1, "status": "success"}],
        replace_scopes={("A", "gnina")},
    )

    def write_combo(name: str) -> None:
        ra._write_optimization_log(
            log,
            [{"combo_name": name, "tool": "gnina", "autodock_rank": 1}],
            replace_scopes={(name, "gnina")},
        )

    threads = [threading.Thread(target=write_combo, args=(name,)) for name in ("B", "C")]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    current = pd.read_csv(log)
    assert set(zip(current.combo_name, current.autodock_rank)) == {
        ("A", 1), ("KEEP", 1), ("B", 1), ("C", 1),
    }


def test_prune_removes_only_obsolete_pose_outputs(tmp_path: Path) -> None:
    tool_dir = tmp_path / "optimized_gnina"
    tool_dir.mkdir()
    keep = tool_dir / "poses_rank1_gnina.sdf"
    stale = tool_dir / "poses_rank2_gnina.sdf"
    unrelated = tool_dir / "other_rank2_gnina.sdf"
    for path in (keep, stale, unrelated):
        path.write_text("data")
        ra._provenance_path(path).write_text("{}")

    removed = ra._prune_obsolete_optimizer_outputs(
        tool_dir, "poses", "gnina", {keep},
    )
    assert removed == 1
    assert keep.exists() and unrelated.exists()
    assert not stale.exists() and not ra._provenance_path(stale).exists()


def test_all_pose_failure_raises_with_structured_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    pose = tmp_path / "protein__ligand_vina_vina_out.pdbqt"
    pose.write_text(_ethanol_pdbqt(wrapped=True))
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("ATOM\n")
    result = ra.DockingResult(
        protein_name="protein", ligand_name="ligand", protein_path=receptor,
        ligand_path=pose, output_dir=tmp_path, status="success", num_poses=1,
        pose_files=[pose],
    )
    monkeypatch.setattr(ra, "preflight_optimizers", lambda cfg: {
        "status": "ready",
        "requested_tools": ["gnina"],
        "optimizers": {"gnina": {"path": "/fake/gnina", "version": "1.0"}},
        "converters": {"mk_export": {"path": "/fake/mk", "version": "0.7"}},
        "settings": {"gnina": ra._optimizer_settings("gnina", cfg)},
    })
    monkeypatch.setattr(
        ra, "_convert_pose_to_sdf",
        lambda *args, **kwargs: (False, "simulated reconstruction failure"),
    )

    with pytest.raises(ra.OptimizationError) as exc_info:
        ra.optimize_autodock_results([result], _gnina_cfg(), tmp_path)

    summary = exc_info.value.summary
    assert summary.status == "failed"
    assert summary.selected_poses == 1 and summary.failed == 1
    assert summary.optimized == summary.reused == 0
    assert Path(summary.log_path).exists()


def test_public_preflight_is_cached(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = 0
    cfg = _gnina_cfg()

    def fake_preflight(_cfg):
        nonlocal calls
        calls += 1
        return (
            {"gnina": "/fake/gnina"}, {"gnina": "1.0"},
            {"mk_export": "/fake/mk", "obabel": ""},
            {"mk_export": "meeko 0.7"}, [],
        )

    monkeypatch.setattr(ra, "_optimizer_preflight", fake_preflight)
    ra._OPTIMIZER_PREFLIGHT_CACHE.clear()
    first = ra.preflight_optimizers(cfg)
    second = ra.preflight_optimizers(cfg)
    assert first == second and calls == 1


def test_dry_run_returns_nonzero_when_requested_optimizer_is_unavailable(tmp_path: Path) -> None:
    receptors = tmp_path / "receptors"
    ligands = tmp_path / "ligands"
    receptors.mkdir()
    ligands.mkdir()
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({
        "receptors_dir": str(receptors),
        "ligand_dirs": [str(ligands)],
        "output_dir": str(tmp_path / "output"),
        "log_dir": str(tmp_path / "logs"),
        "vina_bin": "/bin/true",
        "prep_tool": "meeko",
        "scoring_function": "vina",
        "optimization": "gnina",
        "gnina_executable": "/definitely/missing/gnina",
        "mk_export_executable": "/home/manndo/anaconda3/envs/vina/bin/mk_export.py",
        "optimize_scoring": "default",
        "gnina_cnn_scoring": "rescore",
        "gnina_cnn_model": "crossdock_default2018_ensemble",
    }))

    proc = subprocess.run(
        [sys.executable, str(Path(ra.__file__)), "--dry-run", "--config", str(config)],
        capture_output=True, text=True, timeout=60,
    )
    assert proc.returncode != 0
    assert "optimizer pre-flight failed" in proc.stdout
    assert "Preparing proteins" not in proc.stdout


def test_discovery_strips_vinardo_suffix_and_records_source_scoring(tmp_path: Path) -> None:
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("ATOM\n")
    pose = tmp_path / "protein__ligand_vinardo_vina_out.pdbqt"
    pose.write_text(_ethanol_pdbqt(wrapped=True))

    found = ra.discover_existing_autodock_results(tmp_path, {"protein": receptor})

    assert len(found) == 1
    assert found[0].ligand_name == "ligand"
    assert found[0].scoring_function == "vinardo"
