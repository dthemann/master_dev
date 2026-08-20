from __future__ import annotations

import csv
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


def _write_reconstruction_fixture(
    tmp_path: Path,
    source_smiles: str,
    output_smiles: str | None = None,
) -> tuple[Path, Path]:
    """Write a REMARK-mapped source pose and coordinate-aligned 3D SDF."""
    expected = Chem.MolFromSmiles(source_smiles)
    output = Chem.MolFromSmiles(output_smiles or source_smiles)
    assert expected is not None and output is not None
    assert expected.GetNumAtoms() == output.GetNumAtoms()

    coordinates = [(idx * 1.1, ((idx * idx) % 7) * 0.2, 0.1) for idx in range(output.GetNumAtoms())]
    conformer = Chem.Conformer(output.GetNumAtoms())
    conformer.Set3D(True)
    for idx, xyz in enumerate(coordinates):
        conformer.SetAtomPosition(idx, xyz)
    output.AddConformer(conformer)
    # Match Meeko's reconstruction path: it adds hydrogens before SDWriter.
    output = Chem.AddHs(output, addCoords=True)
    sdf = tmp_path / "reconstructed.sdf"
    writer = Chem.SDWriter(str(sdf))
    writer.write(output)
    writer.close()

    atom_types = {
        "C": "C", "N": "N", "O": "OA", "S": "S", "P": "P",
        "F": "F", "Cl": "CL", "Br": "BR", "I": "I",
    }
    atoms = "".join(
        _pdbqt_atom(
            idx + 1, atom.GetSymbol(), coordinates[idx], atom_types[atom.GetSymbol()],
        )
        for idx, atom in enumerate(expected.GetAtoms())
    )
    mapping = " ".join(f"{idx} {idx}" for idx in range(1, expected.GetNumAtoms() + 1))
    source = tmp_path / "source.pdbqt"
    source.write_text(
        f"REMARK SMILES {source_smiles}\nREMARK SMILES IDX {mapping}\n{atoms}"
    )
    return source, sdf


def _write_mgl_template_fixture(tmp_path: Path) -> tuple[dict, Path, Path]:
    """Write an SDF graph whose prepared PDBQT uses a different atom order."""
    template_root = tmp_path / "templates"
    prepared_root = tmp_path / "prepared"
    template_root.mkdir(parents=True)
    prepared_root.mkdir(parents=True)
    template = template_root / "ligand.sdf"
    _write_sdf(template, "CCO")

    # Template indices 0,1,2 map to PDBQT serials 2,3,1 respectively.
    prepared = prepared_root / "ligand.pdbqt"
    prepared.write_text(
        _pdbqt_atom(1, "O", (3.0, 0.0, 0.1), "OA")
        + _pdbqt_atom(2, "C", (0.0, 0.0, 0.1), "C")
        + _pdbqt_atom(3, "C", (1.5, 0.0, 0.1), "C")
    )
    pose = tmp_path / "pose.pdbqt"
    pose.write_text(
        _pdbqt_atom(1, "O", (12.0, 3.0, -1.0), "OA")
        + _pdbqt_atom(2, "C", (9.0, 3.0, -1.0), "C")
        + _pdbqt_atom(3, "C", (10.5, 3.0, -1.0), "C")
    )
    cfg = {
        "optimize_require_template_reconstruction": True,
        "optimize_ligand_template_root": str(template_root),
        "optimize_prepared_ligand_pdbqt_dir": str(prepared_root),
    }
    return cfg, pose, template


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


def _gnina_refinement_cfg() -> dict:
    cfg = _gnina_cfg()
    cfg.update({
        "optimization": "gnina_refinement",
        "gnina_cnn_scoring": "refinement",
    })
    return cfg


def test_gnina_refinement_alias_is_explicit_and_legacy_all_is_unchanged(
    tmp_path: Path,
) -> None:
    assert ra.resolve_optimizers(_gnina_cfg()) == ["gnina"]
    assert ra.resolve_optimizers(_gnina_refinement_cfg()) == ["gnina_refinement"]
    assert ra.resolve_optimizers({"optimization": "all"}) == ["smina", "gnina"]
    assert ra._optimizer_base_tool("gnina_refinement") == "gnina"
    gnina = tmp_path / "gnina"
    gnina.write_text("fixture")
    assert ra._resolve_optimizer_exe(
        "gnina_refinement", {"gnina_executable": str(gnina)},
    ) == str(gnina)

    wrong_rescore = {**_gnina_cfg(), "gnina_cnn_scoring": "refinement"}
    with pytest.raises(ValueError, match="requires gnina_cnn_scoring='rescore'"):
        ra._validate_optimizer_config(wrong_rescore)

    wrong_refinement = {
        **_gnina_refinement_cfg(), "gnina_cnn_scoring": "rescore",
    }
    with pytest.raises(ValueError, match="requires gnina_cnn_scoring='refinement'"):
        ra._validate_optimizer_config(wrong_refinement)


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


def test_reconstruction_accepts_noninvariant_sdf_stereo_encoding(tmp_path: Path) -> None:
    smiles = "Nc1cc(C(Cl)=C(Cl)Cl)c(S(N)(=O)=O)cc1S(N)(=O)=O"
    source, sdf = _write_reconstruction_fixture(tmp_path, smiles)
    expected = Chem.MolFromSmiles(smiles)
    reconstructed = next(
        mol for mol in Chem.SDMolSupplier(str(sdf), removeHs=False) if mol
    )
    expected_identity = ra._molecule_identity(expected)
    reconstructed_identity = ra._molecule_identity(reconstructed)

    assert expected_identity["canonical_isomeric_smiles"] == reconstructed_identity["canonical_isomeric_smiles"]
    assert expected_identity["bonds_by_index"] != reconstructed_identity["bonds_by_index"]
    valid, message, _ = ra._validate_reconstructed_sdf(source, sdf)
    assert valid, message


@pytest.mark.parametrize(("source_smiles", "mutated_smiles"), [
    ("CCCO", "CCOC"),
    ("F[C@H](Cl)Br", "F[C@@H](Cl)Br"),
])
def test_reconstruction_rejects_real_topology_or_stereochemistry_mutation(
    tmp_path: Path, source_smiles: str, mutated_smiles: str,
) -> None:
    source, sdf = _write_reconstruction_fixture(
        tmp_path, source_smiles, mutated_smiles,
    )
    valid, message, identity = ra._validate_reconstructed_sdf(source, sdf)
    assert not valid
    assert "does not match REMARK SMILES topology" in message
    assert identity is None


def test_mgl_template_reconstruction_copies_coordinates_without_bond_perception(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg, pose, template_path = _write_mgl_template_fixture(tmp_path)
    context = ra._resolve_mgl_pose_template("ligand", cfg)
    assert context is not None
    assert context.atom_serial_map == {0: 2, 1: 3, 2: 1}
    monkeypatch.setattr(
        ra.subprocess, "run",
        lambda *args, **kwargs: pytest.fail(
            "template reconstruction must not invoke Open Babel or Meeko"
        ),
    )

    output = tmp_path / "reconstructed_from_template.sdf"
    ok, converter = ra._convert_pose_to_sdf(
        pose, output,
        {"mk_export": "/fake/mk", "obabel": "/fake/obabel"},
        context,
    )
    valid, message, identity = ra._validate_reconstructed_sdf(
        pose, output, context,
    )
    template = next(
        mol for mol in Chem.SDMolSupplier(str(template_path), removeHs=False) if mol
    )
    reconstructed = next(
        mol for mol in Chem.SDMolSupplier(str(output), removeHs=False) if mol
    )

    assert ok and converter == "rdkit_template_map"
    assert valid, message
    assert identity == ra._molecule_identity(template)
    assert Chem.MolToSmiles(reconstructed, isomericSmiles=True) == "CCO"
    conformer = reconstructed.GetConformer()
    assert tuple(conformer.GetAtomPosition(0)) == pytest.approx((9.0, 3.0, -1.0))
    assert tuple(conformer.GetAtomPosition(1)) == pytest.approx((10.5, 3.0, -1.0))
    assert tuple(conformer.GetAtomPosition(2)) == pytest.approx((12.0, 3.0, -1.0))


def test_mgl_template_reconstruction_fails_closed_on_pose_signature_drift(
    tmp_path: Path,
) -> None:
    cfg, pose, _ = _write_mgl_template_fixture(tmp_path)
    context = ra._resolve_mgl_pose_template("ligand", cfg)
    assert context is not None
    pose.write_text(pose.read_text().replace(" OA\n", " N \n", 1))

    ok, message = ra._convert_pose_to_sdf(pose, tmp_path / "bad.sdf", {}, context)

    assert not ok
    assert "serial/type signature differs" in message
    assert not (tmp_path / "bad.sdf").exists()


def test_required_template_preflight_registers_internal_converter_without_mk_export(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg, _, _ = _write_mgl_template_fixture(tmp_path)
    cfg.update(_gnina_cfg())
    cfg["optimize_require_mk_export"] = True
    model_file = tmp_path / "cnn.model"
    model_file.write_text("fixture")
    cfg["gnina_cnn_model_file"] = str(model_file)
    monkeypatch.setattr(ra, "_resolve_optimizer_exe", lambda tool, config: "/fake/gnina")
    monkeypatch.setattr(ra, "_probe_optimizer", lambda tool, exe, config: (True, "1.0"))
    monkeypatch.setattr(ra, "_resolve_pose_converters", lambda config: {
        "mk_export": "", "obabel": "",
    })

    _, _, converters, versions, errors = ra._optimizer_preflight(cfg)

    assert errors == []
    assert converters["rdkit_template_map"] == str(Path(ra.__file__).resolve())
    assert "mgltools-template-coordinate-map-v1" in versions["rdkit_template_map"]


def test_template_inputs_are_part_of_optimizer_cache_provenance(tmp_path: Path) -> None:
    cfg, pose, _ = _write_mgl_template_fixture(tmp_path)
    context = ra._resolve_mgl_pose_template("ligand", cfg)
    assert context is not None
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n")
    provenance = ra._build_optimizer_provenance(
        tool="gnina_refinement",
        source_pose=pose,
        source_pose_container=pose,
        receptor=receptor,
        cfg=_gnina_refinement_cfg(),
        converter="rdkit_template_map",
        converter_exe=str(Path(ra.__file__).resolve()),
        converter_version="test",
        optimizer_exe="/fake/gnina",
        optimizer_version="1.0",
        input_identity=ra._molecule_identity(context.molecule),
        converter_inputs=context.converter_inputs(),
    )

    inputs = provenance["converter"]["inputs"]
    assert inputs["schema"] == "mgltools-template-coordinate-map-v1"
    assert inputs["template_sdf_sha256"] == context.template_sha256
    assert inputs["prepared_pdbqt_sha256"] == context.prepared_pdbqt_sha256
    assert inputs["atom_serial_mapping_sha256"] == context.mapping_sha256
    assert provenance["fingerprint"]


def test_optimizer_wires_required_template_converter_and_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    template_cfg, raw_pose, _ = _write_mgl_template_fixture(tmp_path)
    body = (
        "REMARK VINA RESULT: -7.500 0.000 0.000\nROOT\n"
        + raw_pose.read_text()
        + "ENDROOT\nTORSDOF 0\n"
    )
    pose = tmp_path / "protein__ligand_vinardo_vina_out.pdbqt"
    pose.write_text("MODEL 1\n" + body + "ENDMDL\n")
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("ATOM\n")
    result = ra.DockingResult(
        protein_name="protein", ligand_name="ligand", protein_path=receptor,
        ligand_path=pose, output_dir=tmp_path, status="success", num_poses=1,
        pose_files=[pose], scoring_function="vinardo",
    )
    cfg = {**_gnina_refinement_cfg(), **template_cfg}
    cfg["optimize_require_mk_export"] = False

    monkeypatch.setattr(ra, "preflight_optimizers", lambda config: {
        "status": "ready",
        "requested_tools": ["gnina_refinement"],
        "optimizers": {
            "gnina_refinement": {"path": "/fake/gnina", "version": "1.0"},
        },
        "converters": {
            "rdkit_template_map": {
                "path": str(Path(ra.__file__).resolve()),
                "version": "test",
            },
        },
        "settings": {
            "gnina_refinement": ra._optimizer_settings(
                "gnina_refinement", config,
            ),
        },
    })

    def fake_optimizer(
        tool, exe, rec, pose_ligand, out_sdf, config, expected_identity=None,
    ):
        mol = next(
            molecule for molecule in Chem.SDMolSupplier(
                str(pose_ligand), removeHs=False,
            ) if molecule
        )
        scores = {
            "minimized_affinity": -8.0,
            "cnn_score": 0.7,
            "cnn_affinity": 6.1,
        }
        mol.SetProp("minimizedAffinity", str(scores["minimized_affinity"]))
        mol.SetProp("CNNscore", str(scores["cnn_score"]))
        mol.SetProp("CNNaffinity", str(scores["cnn_affinity"]))
        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(out_sdf))
        writer.write(mol)
        writer.close()
        return True, scores, "gnina_refinement ok"

    monkeypatch.setattr(ra, "run_optimizer_tool", fake_optimizer)

    summary = ra.optimize_autodock_results([result], cfg, tmp_path)
    output = (
        tmp_path / "optimized_gnina_refinement"
        / f"{pose.stem}_rank1_gnina_refinement.sdf"
    )
    provenance = json.loads(ra._provenance_path(output).read_text())

    assert summary.status == "completed" and summary.optimized == 1
    assert provenance["converter"]["name"] == "rdkit_template_map"
    assert provenance["converter"]["inputs"]["template_sdf_sha256"]
    assert provenance["converter"]["inputs"]["prepared_pdbqt_sha256"]


@pytest.mark.parametrize(("complex_id", "representation_field"), [
    ("7PRI_7TI", "bonds_by_index"),
    ("7T0D_FPP", "bonds_by_index"),
    ("7UJ4_OQ4", "atoms_by_index"),
    ("7XQZ_FPF", "bonds_by_index"),
])
def test_real_meeko07_roundtrip_representation_changes_are_chemically_equivalent(
    tmp_path: Path, complex_id: str, representation_field: str,
) -> None:
    result_root = Path("Dockings/vina_results_full_protein_vina_scoring") / complex_id
    candidates = list(result_root.rglob("*_vina_out.pdbqt"))
    meeko = Path("/home/manndo/anaconda3/envs/vina/bin/mk_export.py")
    if len(candidates) != 1 or not meeko.exists():
        pytest.skip("production Vina output or pinned Meeko 0.7.1 is unavailable")
    _, model, _ = ra._split_pdbqt_models(candidates[0], tmp_path / "split")[0]
    output = tmp_path / f"{complex_id}.sdf"

    ok, converter = ra._convert_pose_to_sdf(
        model, output, {"mk_export": str(meeko), "obabel": ""},
    )
    expected = Chem.MolFromSmiles(ra._remark_smiles(model.read_text()))
    reconstructed = next(
        mol for mol in Chem.SDMolSupplier(str(output), removeHs=False) if mol
    )
    expected_identity = ra._molecule_identity(expected)
    reconstructed_identity = ra._molecule_identity(reconstructed)

    assert ok and converter == "mk_export"
    assert expected_identity["canonical_isomeric_smiles"] == reconstructed_identity["canonical_isomeric_smiles"]
    assert expected_identity[representation_field] != reconstructed_identity[representation_field]
    valid, message, _ = ra._validate_reconstructed_sdf(model, output)
    assert valid, message


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


def test_real_orai_vinardo_5sak_uses_authoritative_template_and_stereo_h(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    template_root = Path("Data/PoseBuster Benchmark Set")
    prepared_root = Path("Dockings/Orai_Benchmark/_staging/ligands/pdbqt")
    pose = Path(
        "Dockings/Orai_Benchmark_Vinardo/docking/"
        "Orai1WT-MDSnap-Fr300__5SAK_ZRY_ligand_start_conf_"
        "vinardo_vina_out.pdbqt"
    )
    if not (
        pose.is_file()
        and (template_root / "5SAK_ZRY/5SAK_ZRY_ligand_start_conf.sdf").is_file()
        and (prepared_root / "5SAK_ZRY_ligand_start_conf.pdbqt").is_file()
    ):
        pytest.skip("Orai Vinardo 5SAK production inputs are unavailable")
    cfg = {
        "optimize_require_template_reconstruction": True,
        "optimize_ligand_template_root": str(template_root),
        "optimize_prepared_ligand_pdbqt_dir": str(prepared_root),
    }
    context = ra._resolve_mgl_pose_template(
        "5SAK_ZRY_ligand_start_conf", cfg,
    )
    assert context is not None
    explicit_hydrogens = [
        atom for atom in context.molecule.GetAtoms() if atom.GetAtomicNum() == 1
    ]
    assert context.molecule.GetNumAtoms() == 19
    assert context.molecule.GetNumHeavyAtoms() == 18
    assert len(explicit_hydrogens) == 1

    rank, model, _ = ra._split_pdbqt_models(pose, tmp_path / "orai_split")[0]
    output = tmp_path / "orai_5sak_rank1.sdf"
    monkeypatch.setattr(
        ra.subprocess, "run",
        lambda *args, **kwargs: pytest.fail(
            "authoritative MGLTools reconstruction must not invoke Open Babel"
        ),
    )
    ok, converter = ra._convert_pose_to_sdf(
        model, output,
        {"mk_export": "/fake/mk", "obabel": "/fake/obabel"},
        context,
    )
    valid, message, identity = ra._validate_reconstructed_sdf(
        model, output, context,
    )

    assert rank == 1 and ok and converter == "rdkit_template_map"
    assert valid, message
    assert identity is not None
    assert identity["canonical_isomeric_smiles"] == (
        "[H]/N=C1/N/C(=N\\Nc2ccccc2)c2ccccc21"
    )
    assert identity["formula"] == "C14H12N4"
    assert context.converter_inputs()["mapped_atom_count"] == 19


def test_real_orai_all_authoritative_mgl_template_pairs_preflight() -> None:
    ids_path = Path(
        "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
    )
    template_root = Path("Data/PoseBuster Benchmark Set")
    prepared_root = Path("Dockings/Orai_Benchmark/_staging/ligands/pdbqt")
    if not (ids_path.is_file() and template_root.is_dir() and prepared_root.is_dir()):
        pytest.skip("Orai benchmark template inputs are unavailable")
    identifiers = sorted({
        line.strip() for line in ids_path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    })
    if len(identifiers) != 308:
        pytest.skip("the local Orai benchmark is not the authoritative 308-ligand set")
    cfg = {
        "optimize_require_template_reconstruction": True,
        "optimize_ligand_template_root": str(template_root),
        "optimize_prepared_ligand_pdbqt_dir": str(prepared_root),
    }

    contexts = ra.preflight_mgl_pose_templates(
        (f"{identifier}_ligand_start_conf" for identifier in identifiers), cfg,
    )
    retained_hydrogen_ligands = {
        name for name, context in contexts.items()
        if any(atom.GetAtomicNum() == 1 for atom in context.molecule.GetAtoms())
    }

    assert len(contexts) == 308
    assert retained_hydrogen_ligands == {
        "5SAK_ZRY_ligand_start_conf",
        "8D5D_5DK_ligand_start_conf",
    }


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


def test_gnina_refinement_alias_invokes_gnina_with_cnn_refinement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    ligand = tmp_path / "input.sdf"
    input_mol = _write_sdf(ligand)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n")
    output = tmp_path / "optimized.sdf"
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append([str(part) for part in command])
        _write_sdf(
            Path(command[command.index("--out") + 1]),
            scores={
                "minimizedAffinity": -8.0,
                "CNNscore": 0.7,
                "CNNaffinity": 6.1,
            },
        )
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(ra.subprocess, "run", fake_run)
    ok, _, _ = ra.run_optimizer_tool(
        "gnina_refinement", "/fake/gnina", receptor, ligand, output,
        _gnina_refinement_cfg(),
        expected_identity=ra._molecule_identity(input_mol),
    )

    assert ok and output.exists()
    command = commands[0]
    assert command[0] == "/fake/gnina"
    assert command[command.index("--cnn_scoring") + 1] == "refinement"


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
    ] + [
        {"combo_name": "A", "tool": "gnina_refinement", "autodock_rank": rank}
        for rank in (1, 2)
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
    assert set(zip(current.combo_name, current.tool, current.autodock_rank)) == {
        ("A", "gnina", 1),
        ("A", "gnina_refinement", 1),
        ("A", "gnina_refinement", 2),
        ("KEEP", "gnina", 1),
        ("B", "gnina", 1),
        ("C", "gnina", 1),
    }


def test_optimization_log_merge_does_not_use_pandas_csv_parser(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A valid cached log must remain resumable if pd.read_csv is unsafe."""
    log = tmp_path / "optimization_log.csv"
    log.write_text(
        "combo_name,tool,autodock_rank,status\n"
        "KEEP,gnina,1,success\n"
        "REPLACE,gnina,1,success\n"
    )

    def fail_native_parser(*_args, **_kwargs):
        raise AssertionError("optimization log must use the stdlib CSV reader")

    monkeypatch.setattr(ra.pd, "read_csv", fail_native_parser)
    ra._write_optimization_log(
        log,
        [{
            "combo_name": "REPLACE", "tool": "gnina",
            "autodock_rank": 1, "status": "success",
        }],
        replace_scopes={("REPLACE", "gnina")},
    )

    with log.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert {(row["combo_name"], row["tool"], row["autodock_rank"]) for row in rows} == {
        ("KEEP", "gnina", "1"),
        ("REPLACE", "gnina", "1"),
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


def test_gnina_rescore_and_refinement_outputs_logs_and_caches_coexist(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    pose = tmp_path / "protein__ligand_vina_vina_out.pdbqt"
    pose.write_text(_ethanol_pdbqt(wrapped=True))
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("ATOM\n")
    result = ra.DockingResult(
        protein_name="protein", ligand_name="ligand", protein_path=receptor,
        ligand_path=pose, output_dir=tmp_path, status="success", num_poses=1,
        pose_files=[pose], scoring_function="vina",
    )

    def fake_preflight(cfg):
        requested = ra.resolve_optimizers(cfg)
        return {
            "status": "ready",
            "requested_tools": requested,
            "optimizers": {
                tool: {"path": "/fake/gnina", "version": "1.0"}
                for tool in requested
            },
            "converters": {
                "mk_export": {"path": "/fake/mk", "version": "0.7"},
            },
            "settings": {
                tool: ra._optimizer_settings(tool, cfg) for tool in requested
            },
        }

    def fake_convert(_pose_file, out_sdf, _converters):
        _write_sdf(out_sdf)
        return True, "mk_export"

    calls: list[str] = []

    def fake_optimizer(
        tool, _exe, _receptor, _pose_ligand, out_sdf, _cfg,
        expected_identity=None,
    ):
        calls.append(tool)
        out_sdf.parent.mkdir(parents=True, exist_ok=True)
        scores = {
            "minimized_affinity": -8.0,
            "cnn_score": 0.7,
            "cnn_affinity": 6.1,
        }
        _write_sdf(
            out_sdf,
            scores={
                "minimizedAffinity": scores["minimized_affinity"],
                "CNNscore": scores["cnn_score"],
                "CNNaffinity": scores["cnn_affinity"],
            },
        )
        return True, scores, f"{tool} ok"

    monkeypatch.setattr(ra, "preflight_optimizers", fake_preflight)
    monkeypatch.setattr(ra, "_convert_pose_to_sdf", fake_convert)
    monkeypatch.setattr(ra, "run_optimizer_tool", fake_optimizer)

    rescore = ra.optimize_autodock_results([result], _gnina_cfg(), tmp_path)
    refinement = ra.optimize_autodock_results(
        [result], _gnina_refinement_cfg(), tmp_path,
    )

    stem = pose.stem
    rescore_out = tmp_path / "optimized_gnina" / f"{stem}_rank1_gnina.sdf"
    refinement_out = (
        tmp_path / "optimized_gnina_refinement"
        / f"{stem}_rank1_gnina_refinement.sdf"
    )
    assert rescore.status == refinement.status == "completed"
    assert calls == ["gnina", "gnina_refinement"]
    assert rescore_out.is_file() and refinement_out.is_file()

    log = pd.read_csv(tmp_path / "optimization_log.csv")
    assert set(log["tool"]) == {"gnina", "gnina_refinement"}
    assert dict(zip(log["tool"], log["cnn_scoring"])) == {
        "gnina": "rescore", "gnina_refinement": "refinement",
    }
    for tool, output in (
        ("gnina", rescore_out), ("gnina_refinement", refinement_out),
    ):
        provenance = json.loads(ra._provenance_path(output).read_text())
        assert provenance["tool"] == tool
        assert provenance["settings"]["cnn_scoring"] == (
            "refinement" if tool == "gnina_refinement" else "rescore"
        )

    # Re-running the legacy protocol must reuse its old cache without touching
    # the independently committed refinement output or log scope.
    rerun = ra.optimize_autodock_results([result], _gnina_cfg(), tmp_path)
    assert rerun.reused == 1 and rerun.optimized == 0
    assert calls == ["gnina", "gnina_refinement"]
    assert refinement_out.is_file()
    assert set(pd.read_csv(tmp_path / "optimization_log.csv")["tool"]) == {
        "gnina", "gnina_refinement",
    }


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
    receptor = tmp_path / "5SAK_ZRY_protein_mgl_tools.pdbqt"
    receptor.write_text("ATOM\n")
    pose = tmp_path / "5SAK_ZRY_protein_mgl_tools__ligand_vinardo_vina_out.pdbqt"
    pose.write_text(_ethanol_pdbqt(wrapped=True))

    found = ra.discover_existing_autodock_results(
        tmp_path, {"5SAK_ZRY_protein_mgl_tools": receptor},
    )

    assert len(found) == 1
    assert found[0].protein_name == "5SAK_ZRY_protein_mgl_tools"
    assert found[0].ligand_name == "ligand"
    assert found[0].scoring_function == "vinardo"
