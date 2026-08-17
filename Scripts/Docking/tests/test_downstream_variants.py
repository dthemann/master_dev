from pathlib import Path
from types import ModuleType, SimpleNamespace
import hashlib
import json
import sys

import numpy as np
import pandas as pd
import pytest


# The default lightweight test environment has RDKit/pandas but not the
# PoseBusters package. Collector/expander tests do not instantiate PoseBusters,
# so provide only the import surface required to load the module.
if "posebusters" not in sys.modules:
    posebusters_stub = ModuleType("posebusters")
    posebusters_stub.PoseBusters = object
    sys.modules["posebusters"] = posebusters_stub

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from Scripts.Analysis import docking_effort_comparison as effort
from Scripts.Analysis import pose_topn
from Scripts.Analysis import posebusters_pose_comparison as pose_comparison
from Scripts.Analysis import posebusters_validity_report as validity_report
from Scripts.Docking import run_unidock2 as ud2
from Scripts.Docking.Posebusters import run_posebusters as pb


def _atom_line(x=0.0, serial=1):
    return (
        f"ATOM  {serial:5d}  C   LIG A   1    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00     0.000 C"
    )


def _pdbqt(*affinities):
    lines = []
    for rank, affinity in enumerate(affinities, start=1):
        lines.extend([
            f"MODEL {rank}",
            f"REMARK VINA RESULT: {affinity:.3f} 0.000 0.000",
            _atom_line(float(rank), rank),
            "ENDMDL",
        ])
    return "\n".join(lines) + "\n"


def _sdf(title="pose"):
    return (
        f"{title}\n"
        "fixture\n"
        "\n"
        "  1  0  0  0  0  0  0  0  0  0  1 V2000\n"
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "M  END\n"
        "$$$$\n"
    )


def _unidock2_sdf(*affinities):
    records = []
    for rank, affinity in enumerate(affinities, start=1):
        record = _sdf(f"pose-{rank}").rsplit("$$$$", 1)[0]
        records.append(
            record
            + ">  <vina_binding_free_energy>\n"
            + f"{affinity:.3f}\n\n$$$$\n"
        )
    return "".join(records)


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _unidock2_tool_identity(tmp_path):
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


def _write_optimizer_sidecar(output, source, receptor, tool, model_digest):
    material = {
        "schema_version": 2,
        "pipeline": "run_autodock.optimize_autodock_results",
        "tool": tool,
        "source_pose_sha256": model_digest,
        "source_pose_container_sha256": _digest(source),
        "receptor_sha256": _digest(receptor),
        "source_scoring": "vina",
        "rank_metric": "minimized_affinity",
        "settings": {"search": "minimize"},
        "optimizer": {"path": "/fixture/gnina", "version": "fixture"},
        "converter": {"name": "mk_export", "path": "/fixture/mk_export",
                      "version": "fixture"},
        "input_identity": {"heavy_atom_count": 1, "bond_count": 0},
    }
    fingerprint = hashlib.sha256(
        json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    sidecar = Path(f"{output}.provenance.json")
    sidecar.write_text(json.dumps({
        **material,
        "fingerprint": fingerprint,
        "source_pose": str(Path(source).resolve()),
        "receptor": str(Path(receptor).resolve()),
        "created_at": "fixture",
        "output_sha256": _digest(output),
        "optimizer_elapsed_time_s": 2.0,
    }))
    return sidecar, fingerprint


def _autodock_fixture(tmp_path):
    root = tmp_path / "autodock"
    prep = root / "CPLX" / "mgl_tools"
    docking = prep / "docking"
    docking.mkdir(parents=True)
    raw = docking / (
        "CPLX_protein_mgl_tools__CPLX_ligand_start_conf_vina_vina_out.pdbqt"
    )
    raw.write_text(_pdbqt(-9.0, -8.0))
    receptor = root / "CPLX" / "_staging" / "receptor.pdbqt"
    receptor.parent.mkdir(parents=True)
    receptor.write_text(_atom_line() + "\n")

    gnina = docking / "optimized_gnina"
    smina = docking / "optimized_smina"
    gnina.mkdir()
    smina.mkdir()
    g1, g2, stale = (gnina / name for name in ("rank1.sdf", "rank2.sdf", "rank3.sdf"))
    failed = smina / "rank1.sdf"
    junk = gnina / "junk.sdf"
    for path in (g1, g2, stale, failed):
        path.write_text(_sdf(path.stem))
    junk.write_text("partial junk")

    provenance = {}
    model_digests = pb._autodock_model_sha256s(raw)
    for rank, path in enumerate((g1, g2, stale), start=1):
        model_digest = (model_digests[rank] if rank in model_digests
                        else hashlib.sha256(f"obsolete-model-{rank}".encode()).hexdigest())
        sidecar, fingerprint = _write_optimizer_sidecar(
            path, raw, receptor, "gnina", model_digest)
        provenance[rank] = {
            "source_pose_sha256": model_digest,
            "receptor_sha256": _digest(receptor),
            "provenance_fingerprint": fingerprint,
            "provenance_file": str(sidecar),
        }

    rows = [
        dict(tool="gnina", autodock_rank=1, optimized_rank=2,
             vina_affinity=-9.0, minimized_affinity=-10.0, cnn_score=0.4,
             cnn_affinity=-8.5, pose_file=raw.name, optimized_file=str(g1),
             status="success", rank_metric="minimized_affinity", elapsed_time_s=2.0,
             **provenance[1]),
        dict(tool="gnina", autodock_rank=2, optimized_rank=1,
             vina_affinity=-8.0, minimized_affinity=-11.0, cnn_score=0.8,
             cnn_affinity=-9.5, pose_file=raw.name, optimized_file=str(g2),
             status="success", rank_metric="minimized_affinity", elapsed_time_s=3.0,
             **provenance[2]),
        # Obsolete rank from a prior longer Vina output: must not be collected.
        dict(tool="gnina", autodock_rank=3, optimized_rank=3,
             pose_file=raw.name, optimized_file=str(stale), status="success",
             minimized_affinity=-7.0, rank_metric="minimized_affinity", elapsed_time_s=4.0,
             **provenance[3]),
        # A failed log row has no committed optimized variant.
        dict(tool="smina", autodock_rank=1, optimized_rank=1,
             pose_file=raw.name, optimized_file=str(failed), status="failed",
             rank_metric="minimized_affinity", elapsed_time_s=5.0),
        # A successful row pointing at a corrupt/truncated SDF is unavailable.
        dict(tool="gnina", autodock_rank=2, optimized_rank=1,
             pose_file="different_source.pdbqt", optimized_file=str(junk),
             status="success", rank_metric="minimized_affinity", elapsed_time_s=1.0),
    ]
    pd.DataFrame(rows).to_csv(prep / "optimization_log.csv", index=False)
    return root, raw, g1, g2


def _unidock2_fixture(tmp_path):
    root = tmp_path / "unidock2"
    out_dir = root / "CPLX"
    out_dir.mkdir(parents=True)
    output = out_dir / "CPLX_unidock2_out.sdf"
    output.write_text(_unidock2_sdf(-9.0, -8.0))
    prepared = out_dir / "CPLX_unidock2_receptor_prepared.pdb"
    prepared.write_text(_atom_line() + "\n")

    inputs_dir = tmp_path / "inputs"
    inputs_dir.mkdir()
    receptor = inputs_dir / "CPLX_protein.pdb"
    ligand = inputs_dir / "CPLX_ligand_start_conf.sdf"
    box = inputs_dir / "CPLX_protein_mgl_tools.box.txt"
    receptor.write_text(_atom_line() + "\n")
    ligand.write_text(_sdf("input"))
    box.write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\n"
        "size_x = 30\nsize_y = 30\nsize_z = 30\n"
    )

    def item(path):
        path = path.resolve()
        return {"path": str(path), "size": path.stat().st_size, "sha256": _digest(path)}

    effective = {
        "Advanced": {"energy_range": 6.0, "exhaustiveness": 512,
                     "mc_steps": 40, "num_pose": 30, "rmsd_limit": 1.0,
                     "seed": 42, "randomize": True, "opt_steps": -1,
                     "refine_steps": 5, "use_tor_lib": False,
                     "energy_decomp": False},
        "Hardware": {"gpu_device_id": 0, "n_cpu": 1},
        "Settings": {"box_size": [30.0, 30.0, 30.0],
                     "task": "screen", "search_mode": "free"},
        "Preprocessing": {
            "construct_ff": False, "template_docking": False,
            "compute_center": True, "covalent_ligand": False,
            "preserve_receptor_hydrogen": False, "engine_checkpoint": False,
        },
    }
    effective_json = json.dumps(effective, sort_keys=True, separators=(",", ":"))
    provenance = {
        "schema_version": 3,
        "engine": "unidock2",
        "tool": _unidock2_tool_identity(tmp_path),
        "inputs": {"receptor": item(receptor), "ligand": item(ligand), "box": item(box)},
        "center": [0.0, 0.0, 0.0],
        "size": [30.0, 30.0, 30.0],
        "effective_config": effective,
        "effective_config_sha256": hashlib.sha256(effective_json.encode()).hexdigest(),
        "driver_execution": {"gpu_lock_file": "/tmp/master_docking_gpu0.lock"},
    }
    canonical = json.dumps(provenance, sort_keys=True, separators=(",", ":"))
    fingerprint = hashlib.sha256(canonical.encode()).hexdigest()
    generation = "a" * 32
    runtime = ud2._requested_runtime_settings(provenance)
    completion_path = out_dir / "CPLX_unidock2_completion.json"
    completion = {
        **provenance,
        "status": "success", "fingerprint": fingerprint,
        "generation_id": generation, "output_file": output.name,
        "output_sha256": _digest(output),
        "prepared_receptor_file": prepared.name,
        "prepared_receptor_sha256": _digest(prepared),
        "prepared_receptor_atoms": 1, "prepared_receptor_heavy_atoms": 1,
        "engine_receptor_heavy_atoms_in_box": 1,
        "receptor_prmtop_sha256": "b" * 64,
        "receptor_inpcrd_sha256": "c" * 64,
        "num_poses": 2, "best_affinity": -9.0,
        "requested_search_settings": runtime,
        "runtime_attestation": runtime,
        "effective_search_settings": runtime,
        "applied_center": [0.0, 0.0, 0.0],
        "applied_box_size": [30.0, 30.0, 30.0],
        "input_protein_heavy_atoms": 1,
        "completed_unix_s": 1.0,
    }
    completion_path.write_text(json.dumps(completion))
    summary_path = out_dir / "docking_summary.json"
    summary = {
        "pdb_id": "CPLX", "status": "success", "engine": "unidock2",
        "scoring": "vina", "fingerprint": fingerprint,
        "generation_id": generation, "num_poses": 2, "best_affinity": -9.0,
        "receptor": str(receptor), "output_sdf": str(output),
        "completion_manifest": str(completion_path),
        "prepared_receptor_pdb": str(prepared),
        "published_output_current": True, "inputs_changed_during_run": False,
        "elapsed_s": 7.0,
    }
    summary_path.write_text(json.dumps(summary))
    return root, output, completion_path, prepared, receptor


def _unidock_tiled_fixture(tmp_path, root_name="unidock"):
    root = tmp_path / root_name
    complete = root / "CPLX"
    complete.mkdir(parents=True)
    merged = complete / "CPLX_unidock_out.pdbqt"
    merged.write_text(_pdbqt(-10.0, -5.0))
    receptor = tmp_path / f"{root_name}_prepared" / "CPLX.pdbqt"
    receptor.parent.mkdir()
    receptor.write_text(_atom_line() + "\n")
    ligand = receptor.parent / "CPLX_ligand.pdbqt"
    ligand.write_text(_atom_line(1.0) + "\n")
    box = receptor.parent / "CPLX.box.txt"
    box.write_text("center_x = 0\n")

    def item(path):
        path = path.resolve()
        return {"path": str(path), "size": path.stat().st_size,
                "sha256": _digest(path)}

    provenance = {
        "schema_version": 2,
        "engine": {"path": "/fixture/unidock", "version": "1.2.0"},
        "inputs": {"receptor": item(receptor), "ligand": item(ligand),
                   "box": item(box)},
        "settings": {"scoring": "vina"},
        "subboxes": [[0.0, 0.0, 0.0, 29.0, 29.0, 29.0]],
    }
    fingerprint = hashlib.sha256(
        json.dumps(provenance, sort_keys=True, separators=(",", ":"),
                   allow_nan=False).encode()
    ).hexdigest()
    generation = "fixture-generation"
    output_sha = _digest(merged)
    common = {
        "commit_status": "complete", "fingerprint": fingerprint,
        "generation_id": generation, "output_file": merged.name,
        "output_sha256": output_sha, "num_poses": 2, "n_subboxes": 1,
    }
    (complete / "run_manifest.json").write_text(json.dumps({
        "schema_version": 2, "status": "complete", "provenance": provenance,
        "best_affinity": -10.0, **common,
    }))
    (complete / ".unidock_done").write_text(json.dumps({
        "schema_version": 2, "engine": "unidock-tiled", "status": "complete",
        **common,
    }))
    (complete / "docking_summary.json").write_text(json.dumps({
        "status": "success", "scoring": "vina", "engine": "unidock-tiled",
        "receptor": str(receptor), "elapsed_s": 12.5, **common,
    }))
    return root, complete, merged, receptor


def test_autodock_collector_joins_only_current_committed_optimizer_rows(tmp_path):
    root, raw, g1, g2 = _autodock_fixture(tmp_path)

    rows = pb._collect_autodock_rows(root)
    assert len(rows) == 3
    raw_rows = [row for row in rows if row["optimizer"] == "original"]
    optimized = [row for row in rows if row["optimizer"] == "gnina"]
    assert len(raw_rows) == 1
    assert raw_rows[0]["file_path"] == str(raw)
    assert raw_rows[0]["pose_count"] == 2
    assert {(row["autodock_rank"], row["optimized_rank"]) for row in optimized} == {
        (1, 2), (2, 1),
    }
    assert {Path(row["file_path"]) for row in optimized} == {g1, g2}

    assert len(pb._collect_autodock_rows(root, {"gnina"})) == 2
    assert len(pb._collect_autodock_rows(root, {"raw"})) == 1
    assert pb._apply_variant_filter(
        [{"method": "autodock", "optimizer": "original"}], {"autodock": {"raw"}}
    )

    sidecar = Path(f"{g1}.provenance.json")
    sidecar_text = sidecar.read_text()
    sidecar.unlink()
    assert {Path(row["file_path"]) for row in pb._collect_autodock_rows(root, {"gnina"})} == {g2}
    sidecar.write_text(sidecar_text)

    pose = pb._expand_autodock_poses(optimized[0], tmp_path / "converted", None)[0]
    assert pose["source_pose_file"] == str(raw.resolve())
    assert pose["optimizer"] == "gnina"
    assert pose["autodock_rank"] == 1
    assert pose["optimized_rank"] == 2


def test_topn_keeps_raw_and_optimized_rank_axes_separate():
    frame = pd.DataFrame([
        dict(docking_method="autodock", protein="P", ligand="L", optimizer="original",
             autodock_rank=1, optimized_rank=None, pose_file="raw1.sdf", pose_name="raw1"),
        dict(docking_method="autodock", protein="P", ligand="L", optimizer="original",
             autodock_rank=2, optimized_rank=None, pose_file="raw2.sdf", pose_name="raw2"),
        dict(docking_method="autodock", protein="P", ligand="L", optimizer="gnina",
             autodock_rank=1, optimized_rank=2, pose_file="opt1.sdf", pose_name="opt1"),
        dict(docking_method="autodock", protein="P", ligand="L", optimizer="gnina",
             autodock_rank=2, optimized_rank=1, pose_file="opt2.sdf", pose_name="opt2"),
    ])
    capped = pose_topn.cap_frame(frame, top_n=1)
    assert set(capped["pose_file"]) == {"raw1.sdf", "opt2.sdf"}


def test_reports_split_and_rank_autodock_variants(tmp_path):
    csv_path = tmp_path / "posebusters.csv"
    rows = [
        dict(docking_method="autodock", protein="P", ligand="L",
             optimizer="original", autodock_rank=1, optimized_rank=None,
             pose_file="converted/raw_model1.sdf", pose_name="raw1",
             mol_pred_loaded=True),
        dict(docking_method="autodock", protein="P", ligand="L",
             optimizer="original", autodock_rank=2, optimized_rank=None,
             pose_file="converted/raw_model2.sdf", pose_name="raw2",
             mol_pred_loaded=True),
        # The source-rank order is reversed by the committed optimizer ranking.
        dict(docking_method="autodock", protein="P", ligand="L",
             optimizer="gnina", autodock_rank=1, optimized_rank=2,
             pose_file="docking/optimized_gnina/opt_source1.sdf", pose_name="opt1",
             mol_pred_loaded=True),
        dict(docking_method="autodock", protein="P", ligand="L",
             optimizer="gnina", autodock_rank=2, optimized_rank=1,
             pose_file="docking/optimized_gnina/opt_source2.sdf", pose_name="opt2",
             mol_pred_loaded=False),
    ]
    pd.DataFrame(rows).to_csv(csv_path, index=False)

    validity = validity_report.load_and_score(csv_path)
    assert set(validity["docking_method"]) == {"autodock", "autodock_gnina"}
    summary = validity_report.per_tool_summary(
        validity, validity_report._method_order(validity))
    assert summary.loc["autodock", "valid_poses"] == 2
    assert summary.loc["autodock_gnina", "valid_poses"] == 1

    index = pose_comparison._build_pose_index(csv_path)
    assert set(index["docking_method"]) == {"autodock", "autodock_gnina"}
    ranks = {
        row.pose_name: pose_comparison._pose_rank(
            row._asdict(), row.docking_method, row.pose_file)
        for row in index.itertuples(index=False)
    }
    assert ranks == {"raw1": 1, "raw2": 2, "opt1": 2, "opt2": 1}
    assert "autodock_gnina" in pose_comparison.RANKING_TOOLS


def test_complex_resume_sums_raw_and_optimizer_rows(tmp_path):
    discovered = pd.DataFrame([
        {"docking_tool": "autodock", "protein": "P", "ligand": "L", "pose_count": 2},
        {"docking_tool": "autodock", "protein": "P", "ligand": "L", "pose_count": 1},
        {"docking_tool": "autodock", "protein": "P", "ligand": "L", "pose_count": 1},
    ])
    prior = tmp_path / "results.csv"
    pd.DataFrame([
        {"docking_method": "autodock", "protein": "P", "ligand": "L"},
        {"docking_method": "autodock", "protein": "P", "ligand": "L"},
    ]).to_csv(prior, index=False)
    assert len(pb.filter_already_processed_complexes(discovered, prior)) == 3

    pd.DataFrame([
        {"docking_method": "autodock", "protein": "P", "ligand": "L"}
        for _ in range(4)
    ]).to_csv(prior, index=False)
    assert pb.filter_already_processed_complexes(discovered, prior).empty


def test_unidock_collector_requires_commit_and_preserves_merged_rank(tmp_path, monkeypatch):
    root, complete, merged, _receptor = _unidock_tiled_fixture(tmp_path)

    incomplete = root / "INCOMPLETE"
    incomplete.mkdir()
    (incomplete / "INCOMPLETE_unidock_out.pdbqt").write_text(_pdbqt(-99.0))
    (incomplete / "docking_summary.json").write_text(json.dumps({
        "status": "incomplete", "num_poses": 1,
    }))

    rows = pb._collect_unidock_rows(root)
    assert len(rows) == 1
    assert rows[0]["protein"] == rows[0]["ligand"] == "CPLX"
    assert rows[0]["pose_count"] == 2

    marker_path = complete / ".unidock_done"
    committed_marker = marker_path.read_text()
    marker_path.write_text("legacy plain completion marker\n")
    assert pb._collect_unidock_rows(root) == []
    marker_path.write_text(committed_marker)

    summary_path = complete / "docking_summary.json"
    committed_summary = json.loads(summary_path.read_text())
    stale_summary = dict(committed_summary, output_sha256="0" * 64)
    summary_path.write_text(json.dumps(stale_summary))
    assert pb._collect_unidock_rows(root) == []
    summary_path.write_text(json.dumps(committed_summary))

    manifest_path = complete / "run_manifest.json"
    committed_manifest = json.loads(manifest_path.read_text())
    stale_manifest = dict(committed_manifest, generation_id="stale-generation")
    manifest_path.write_text(json.dumps(stale_manifest))
    assert pb._collect_unidock_rows(root) == []
    manifest_path.write_text(json.dumps(committed_manifest))

    converted = tmp_path / "converted"
    converted.mkdir()
    outputs = [converted / "merged_model1.sdf", converted / "merged_model2.sdf"]
    for output in outputs:
        output.write_text(_sdf(output.stem))
    monkeypatch.setattr(pb, "_convert_pdbqt_to_sdf",
                        lambda *_args, **_kwargs: [str(path) for path in outputs])
    poses = pb._expand_unidock_poses(rows[0], converted, SimpleNamespace())
    assert [(pose["unidock_rank"], pose["unidock_affinity"]) for pose in poses] == [
        (1, -10.0), (2, -5.0),
    ]


def test_unidock_parser_rejects_truncated_or_scoreless_models(tmp_path):
    truncated = tmp_path / "truncated.pdbqt"
    truncated.write_text(_pdbqt(-8.0) + "MODEL 2\nREMARK VINA RESULT: -7\n")
    scoreless = tmp_path / "scoreless.pdbqt"
    scoreless.write_text("MODEL 1\n" + _atom_line() + "\nENDMDL\n")
    assert pb._closed_scored_pdbqt_models(truncated) is None
    assert pb._closed_scored_pdbqt_models(scoreless) is None


def test_unidock2_requires_authoritative_current_schema3_commit(tmp_path):
    root, output, completion_path, prepared, receptor = _unidock2_fixture(tmp_path)
    rows = pb._collect_unidock2_rows(root)
    assert len(rows) == 1
    assert rows[0]["pose_count"] == 2
    assert rows[0]["unidock2_output_sha256"] == _digest(output)
    assert rows[0]["docking_receptor_file"] == str(prepared.resolve())

    expanded = pb._expand_unidock2_poses(
        rows[0], tmp_path / "unidock2_split", SimpleNamespace())
    assert [(pose["unidock2_rank"], pose["unidock2_affinity"])
            for pose in expanded] == [(1, -9.0), (2, -8.0)]
    resolved, prepared_used, stage = pb._resolve_unidock2_receptor(
        expanded[0], tmp_path / "unused")
    assert resolved == prepared_used == prepared.resolve()
    assert stage == "unidock2_manifest_attested_engine_prepared_pdb"

    committed = completion_path.read_text()
    completion_path.unlink()
    assert pb._collect_unidock2_rows(root) == []
    completion_path.write_text(json.dumps({"schema_version": 1, "status": "success"}))
    assert pb._collect_unidock2_rows(root) == []
    completion_path.write_text(committed)

    output_text = output.read_text()
    output.write_text(output_text.replace("-9.000", "-7.000", 1))
    assert pb._collect_unidock2_rows(root) == []
    output.write_text(output_text)

    receptor_text = receptor.read_text()
    receptor.write_text(receptor_text + "REMARK tampered\n")
    assert pb._collect_unidock2_rows(root) == []
    receptor.write_text(receptor_text)

    prepared_text = prepared.read_text()
    prepared.write_text(prepared_text + "REMARK tampered\n")
    assert pb._collect_unidock2_rows(root) == []
    with pytest.raises(ValueError, match="missing, malformed, or stale"):
        pb._resolve_unidock2_receptor(expanded[0], tmp_path / "unused")
    prepared.write_text(prepared_text)


def test_autodock_effort_adds_optimizer_wall_and_device_work(tmp_path):
    methods = effort._per_pose_methods(SimpleNamespace(
        autodock_refine="gnina", diffdock_refine="raw",
        eq_mode="unguided", eq_refine="raw"))
    assert methods["autodock"] == "autodock_gnina"

    prep = tmp_path / "CPLX" / "meeko"
    prep.mkdir(parents=True)
    (prep / "docking_summary.json").write_text(json.dumps({
        "overall": {"total_time_seconds": 10.0},
    }))
    pd.DataFrame([
        {"tool": "gnina", "elapsed_time_s": 2.0, "status": "success"},
        {"tool": "gnina", "elapsed_time_s": 3.0, "status": "success"},
    ]).to_csv(prep / "optimization_log.csv", index=False)

    cpu = effort._autodock_effort(
        "CPLX", tmp_path, "meeko", "gnina", docking_cpu=32,
        optimizer_cpu=4, gnina_gpu=False,
    )
    gpu = effort._autodock_effort(
        "CPLX", tmp_path, "meeko", "gnina", docking_cpu=32,
        optimizer_cpu=4, gnina_gpu=True,
    )
    assert cpu == (15.0, 0.0, 340.0)
    assert gpu == (15.0, 5.0, 320.0)

    pd.DataFrame([
        {"tool": "gnina", "elapsed_time_s": 4.0, "status": "failed"},
    ]).to_csv(prep / "optimization_log.csv", index=False)
    assert effort._autodock_opt_time("CPLX", tmp_path, "meeko", "gnina") is None

    # Once the variant has a committed-success row, retain time spent on its
    # failed attempt instead of making that consumed compute disappear.
    pd.DataFrame([
        {"tool": "gnina", "elapsed_time_s": 2.0, "status": "success"},
        {"tool": "gnina", "elapsed_time_s": 4.0, "status": "failed"},
    ]).to_csv(prep / "optimization_log.csv", index=False)
    assert effort._autodock_opt_time("CPLX", tmp_path, "meeko", "gnina") == 6.0

    mixed = pd.DataFrame({
        "pose_file": ["converted/raw_model1.sdf", "docking/optimized_gnina/rank1.sdf"],
    })
    selected = effort._filter_optimizer_variant(
        mixed, "autodock", {"autodock": "gnina"})
    assert selected["pose_file"].tolist() == ["docking/optimized_gnina/rank1.sdf"]


def test_unidock_effort_requires_matching_commit_and_keeps_engines_separate(tmp_path):
    tiled_root, tiled, _merged, _receptor = _unidock_tiled_fixture(
        tmp_path, root_name="tiled")
    assert effort._unidock_effort("CPLX", tiled_root, "unidock") == (
        12.5, 12.5, 0.0,
    )

    single_root, _output, _completion, _prepared, _input = _unidock2_fixture(tmp_path)
    single = single_root / "CPLX"
    wall, gpu_s, cpu_s = effort._unidock_effort(
        "CPLX", single_root, "unidock2")
    assert wall == 7.0
    assert np.isnan(gpu_s) and np.isnan(cpu_s)

    completion_path = single / "CPLX_unidock2_completion.json"
    summary_path = single / "docking_summary.json"
    committed_completion = json.loads(completion_path.read_text())
    committed_summary = json.loads(summary_path.read_text())
    legacy = dict(committed_completion, schema_version=2)
    provenance = {
        key: legacy[key] for key in (
            "schema_version", "engine", "tool", "inputs", "center", "size",
            "effective_config", "effective_config_sha256", "driver_execution",
        )
    }
    canonical = json.dumps(provenance, sort_keys=True, separators=(",", ":"))
    legacy["fingerprint"] = hashlib.sha256(canonical.encode()).hexdigest()
    completion_path.write_text(json.dumps(legacy))
    summary_path.write_text(json.dumps(
        dict(committed_summary, fingerprint=legacy["fingerprint"])))
    assert effort._unidock_effort("CPLX", single_root, "unidock2") is None
    completion_path.write_text(json.dumps(committed_completion))
    summary_path.write_text(json.dumps(committed_summary))

    # A stale summary must not be attributed to the current completion marker.
    summary = json.loads(summary_path.read_text())
    summary["fingerprint"] = "stale"
    summary_path.write_text(json.dumps(summary))
    assert effort._unidock_effort("CPLX", single_root, "unidock2") is None
