from pathlib import Path
from types import SimpleNamespace
import sys

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from Scripts.Docking import run_unidock_tiled as tiled


def atom_line(x=0.0, y=0.0, z=0.0, serial=1):
    return (
        f"ATOM  {serial:5d}  C   LIG A   1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 C"
    )


def pose_file(*affinity_xyz):
    chunks = []
    for rank, (affinity, x) in enumerate(affinity_xyz, start=1):
        chunks.extend(
            [
                f"MODEL {rank}",
                f"REMARK VINA RESULT: {affinity} 0.000 0.000",
                atom_line(x=x),
                "ENDMDL",
            ]
        )
    return "\n".join(chunks) + "\n"


def make_complex(base: Path, name="C1", ligand_text=None, malformed_box=False):
    cdir = base / name
    recdir = cdir / "_staging" / "receptors" / "pdbqt"
    ligdir = cdir / "_staging" / "ligands" / "pdbqt"
    recdir.mkdir(parents=True)
    ligdir.mkdir(parents=True)
    (recdir / "protein_mgl_tools.pdbqt").write_text(atom_line() + "\n")
    if malformed_box:
        (recdir / "protein_mgl_tools.box.txt").write_text("center_x = nope\n")
    else:
        (recdir / "protein_mgl_tools.box.txt").write_text(
            "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
            "size_x = 1\nsize_y = 1\nsize_z = 1\n"
        )
    ligand = ligdir / "test_ligand_start_conf.pdbqt"
    ligand.write_text(ligand_text or atom_line() + "\n")
    return cdir, ligand


def make_cfg(tmp_path: Path, **overrides):
    cfg = {
        "unidock_bin": str(tmp_path / "unidock"),
        "prep_base": str(tmp_path / "prep"),
        "converter": "mgl_tools",
        "output_dir": str(tmp_path / "out"),
        "log_dir": str(tmp_path / "logs"),
        "overwrite": False,
        "keep_tiles": True,
        "scoring": "vina",
        "execution_mode": "auto",
        "subbox_size": 29.0,
        "spacing": 0.375,
        "overlap_slack": 1.0,
        "prefilter_cutoff": 8.0,
        "prefilter_min_atoms": 0,
        "exhaustiveness": 8,
        "num_modes": 30,
        "energy_range": 6.0,
        "seed": 42,
        "paired_batch_size": 6,
        "exclude": [],
        "max_receptor_atoms": None,
    }
    cfg.update(overrides)
    return cfg


@pytest.mark.parametrize(
    "text",
    [
        "ENDMDL\nMODEL 1\n" + atom_line() + "\n",
        "REMARK mentions MODEL and ENDMDL but is not a pose\n",
        pose_file((-7.0, 0.0)) + "MODEL 2\nREMARK VINA RESULT: -6\n" + atom_line() + "\n",
        "MODEL 1\n" + atom_line() + "\nENDMDL\n",
        "MODEL 1\nREMARK VINA RESULT: -7\nENDMDL\n",
    ],
)
def test_strict_output_validation_rejects_prior_false_positives(tmp_path, text):
    output = tmp_path / "tile.pdbqt"
    output.write_text(text)
    assert not tiled._valid_out(output)
    with pytest.raises(ValueError):
        tiled.parse_out_pdbqt(output)


def test_strict_output_validation_accepts_complete_scored_models(tmp_path):
    output = tmp_path / "tile.pdbqt"
    output.write_text(pose_file((-8.0, 0.0), (-4.0, 3.0)))
    assert tiled._valid_out(output)
    assert [pose.affinity for pose in tiled.parse_out_pdbqt(output)] == [-8.0, -4.0]


def test_unidock_120_auto_falls_back_for_six_kcal_range(tmp_path):
    cfg = make_cfg(tmp_path, energy_range=6)
    assert tiled.resolve_execution_mode(cfg, "1.2.0") == "standard"
    with pytest.raises(ValueError, match="cannot be verified"):
        tiled.resolve_execution_mode({**cfg, "execution_mode": "paired"}, "1.2.0")
    compatible = {**cfg, "energy_range": 3, "execution_mode": "paired"}
    assert tiled.resolve_execution_mode(compatible, "1.2.0") == "paired"
    assert tiled.resolve_execution_mode(cfg, "9.9.9") == "standard"


def test_grid_cap_and_overlap_validation(tmp_path):
    tiled.validate_config(make_cfg(tmp_path, subbox_size=30.0))  # 81^3, exactly allowed
    with pytest.raises(ValueError, match="grid points"):
        tiled.validate_config(make_cfg(tmp_path, subbox_size=30.1))
    with pytest.raises(ValueError, match="smaller than subbox_size"):
        tiled.validate_config(make_cfg(tmp_path, overlap_slack=29.0))


def test_planning_rejects_inputs_without_parseable_heavy_atoms(tmp_path):
    cfg = make_cfg(tmp_path)
    cdir, ligand = make_complex(Path(cfg["prep_base"]))
    ligand.write_text("REMARK no ligand atoms\n")
    info, boxes = tiled.plan_complex(cdir, cfg)
    assert boxes is None
    assert info["skip"] == "ligand has no parseable heavy atoms"

    ligand.write_text(atom_line() + "\n")
    receptor = next((cdir / "_staging" / "receptors" / "pdbqt").glob("*.pdbqt"))
    receptor.write_text("REMARK no receptor atoms\n")
    info, boxes = tiled.plan_complex(cdir, cfg)
    assert boxes is None
    assert info["skip"] == "receptor has no parseable heavy atoms"


def test_standard_mode_honours_six_kcal_range_and_publishes_atomically(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir, _ = make_complex(Path(cfg["prep_base"]))
    seen_commands = []

    def fake_run(cmd, **_kwargs):
        seen_commands.append(cmd)
        out = Path(cmd[cmd.index("--out") + 1])
        out.write_text(pose_file((-10.0, 0.0), (-5.0, 3.0)))
        return SimpleNamespace(returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(tiled.subprocess, "run", fake_run)
    result = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )

    assert result["status"] == "success"
    assert result["num_poses"] == 2
    assert result["execution_mode"] == "standard"
    assert any(
        cmd[cmd.index("--energy_range") + 1] in {"6", "6.0"} for cmd in seen_commands
    )
    merged = Path(cfg["output_dir"]) / cdir.name / f"{cdir.name}_unidock_out.pdbqt"
    assert [pose.affinity for pose in tiled.parse_out_pdbqt(merged)] == [-10.0, -5.0]
    marker = tiled._read_json(merged.parent / ".unidock_done")
    manifest = tiled._read_json(merged.parent / "run_manifest.json")
    assert marker["status"] == marker["commit_status"] == "complete"
    assert manifest["status"] == manifest["commit_status"] == "complete"
    for key in ("fingerprint", "generation_id", "output_file", "output_sha256",
                "num_poses", "n_subboxes", "commit_status"):
        assert result[key] == marker[key] == manifest[key]

    command_count = len(seen_commands)
    cached = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert cached["status"] == "done"
    assert len(seen_commands) == command_count
    for key in ("fingerprint", "generation_id", "output_file", "output_sha256",
                "num_poses", "n_subboxes", "commit_status"):
        assert cached[key] == marker[key]


def test_failed_overwrite_ignores_stale_tiles_and_removes_stale_sentinel(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path, overwrite=True)
    cdir, _ = make_complex(Path(cfg["prep_base"]))
    out_dir = Path(cfg["output_dir"]) / cdir.name
    legacy_tiles = out_dir / "_tiles"
    legacy_tiles.mkdir(parents=True)
    (legacy_tiles / f"{cdir.name}_s0_out.pdbqt").write_text(pose_file((-99.0, 0.0)))
    done = out_dir / ".unidock_done"
    done.write_text("legacy completion marker\n")

    monkeypatch.setattr(
        tiled.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(returncode=1, stdout="", stderr="failed"),
    )
    first = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert first["status"] == "incomplete"
    assert first["best_affinity"] is None
    assert first["n_subboxes_missing"] == 1
    assert first["process_returncodes"] == [1]
    assert not done.exists()

    # Turning overwrite off resumes the failed generation; it must not trust the
    # old marker or the legacy -99 tile on the next run.
    cfg["overwrite"] = False
    second = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert second["status"] == "incomplete"
    assert second["generation_id"] == first["generation_id"]
    assert not done.exists()


def test_failed_overwrite_of_current_success_cannot_reuse_current_tiles(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir, _ = make_complex(Path(cfg["prep_base"]))
    fail = False

    def fake_run(cmd, **_kwargs):
        if fail:
            return SimpleNamespace(returncode=1, stdout="", stderr="failed")
        Path(cmd[cmd.index("--out") + 1]).write_text(pose_file((-8.0, 0.0)))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(tiled.subprocess, "run", fake_run)
    first = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert first["status"] == "success"

    fail = True
    # Even a missing marker must not turn overwrite into cache repair/skip.
    (Path(cfg["output_dir"]) / cdir.name / ".unidock_done").unlink()
    cfg["overwrite"] = True
    second = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert second["status"] == "incomplete"
    assert second["generation_id"] != first["generation_id"]
    assert second["best_affinity"] is None
    assert second["n_subboxes_missing"] == 1
    assert not (Path(cfg["output_dir"]) / cdir.name / ".unidock_done").exists()


def test_nonzero_paired_return_code_cannot_commit_success(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path, energy_range=3)
    cdir, _ = make_complex(Path(cfg["prep_base"]))
    calls = 0

    def writes_output_but_fails(cmd, **_kwargs):
        nonlocal calls
        calls += 1
        index = tiled._read_json(Path(cmd[cmd.index("--ligand_index") + 1]))
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        for entry in index.values():
            ligand = Path(entry["ligand"])
            (out_dir / f"{ligand.stem}_out.pdbqt").write_text(pose_file((-8.0, 0.0)))
        return SimpleNamespace(returncode=137, stdout="", stderr="OOM")

    monkeypatch.setattr(tiled.subprocess, "run", writes_output_but_fails)
    first = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="paired",
    )
    assert first["status"] == "incomplete"
    assert first["n_subboxes_missing"] == 0
    assert not (Path(cfg["output_dir"]) / cdir.name / ".unidock_done").exists()

    # On retry, strict complete outputs from the interrupted paired job are
    # reusable. With no failing subprocess in this attempt, completion commits.
    second = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="paired",
    )
    assert second["status"] == "success"
    assert second["n_subboxes_attempted"] == 0
    assert calls == 1


def test_provenance_change_starts_clean_generation_and_refreshes_ligand(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    cdir, ligand = make_complex(Path(cfg["prep_base"]))
    should_fail = False

    def fake_run(cmd, **_kwargs):
        if should_fail:
            return SimpleNamespace(returncode=1, stdout="", stderr="failed")
        Path(cmd[cmd.index("--out") + 1]).write_text(pose_file((-8.0, 0.0)))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(tiled.subprocess, "run", fake_run)
    first = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert first["status"] == "success"

    new_ligand = atom_line(x=1.25) + "\n"
    ligand.write_text(new_ligand)
    should_fail = True
    second = tiled.dock_complex(
        cdir, cfg, cfg["unidock_bin"], unidock_version="1.2.0",
        execution_mode="standard",
    )
    assert second["status"] == "incomplete"
    assert second["generation_id"] != first["generation_id"]
    assert second["fingerprint"] != first["fingerprint"]
    copied = (
        Path(cfg["output_dir"]) / cdir.name / "_tiles" / second["generation_id"]
        / f"{cdir.name}_s0.pdbqt"
    )
    assert copied.read_text() == new_ligand
    assert not (Path(cfg["output_dir"]) / cdir.name / ".unidock_done").exists()


def test_plan_errors_are_reported_without_stopping_other_complexes(tmp_path, monkeypatch):
    cfg = make_cfg(tmp_path)
    Path(cfg["unidock_bin"]).write_text("placeholder\n")
    make_complex(Path(cfg["prep_base"]), name="GOOD")
    make_complex(Path(cfg["prep_base"]), name="BAD", malformed_box=True)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(cfg))
    monkeypatch.setattr(tiled, "probe_unidock_version", lambda _binary: "1.2.0")

    plan = tiled.main(config_path, plan_only=True)
    assert plan["dockable"] == 1
    assert len(plan["errors"]) == 1
    assert plan["errors"][0][0] == "BAD"


def test_cli_exit_code_fails_for_incomplete_or_error():
    assert tiled.results_exit_code([{"status": "success"}, {"status": "done"}]) == 0
    assert tiled.results_exit_code([{"status": "incomplete"}]) == 1
    assert tiled.results_exit_code([{"status": "error"}]) == 1
