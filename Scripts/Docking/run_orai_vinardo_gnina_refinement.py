#!/usr/bin/env python3
"""Refine completed Orai x benchmark AutoDock Vina/Vinardo poses with GNINA.

This is intentionally an optimizer-only driver.  It reconstructs
``DockingResult`` objects from committed PDBQT files and never calls the
AutoDock docking entry point.  Vina poses use their embedded Meeko SMILES atom
map; legacy Vinardo poses are rebuilt from the authoritative benchmark SDF
topology using a validated prepared-PDBQT atom serial map.  Each
receptor-ligand complex is checkpointed independently so an interrupted run
can safely reuse validated SDF/provenance caches and the merged optimization
log.
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
import re
import signal
import sys
import tempfile
from collections import Counter
from pathlib import Path
from types import FrameType
from typing import Any, Mapping

import yaml
from rdkit import Chem


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
DEFAULT_OUTPUT_DIR = REPO_ROOT / "Dockings/Orai_Benchmark_Vinardo"
DEFAULT_VINA_OUTPUT_DIR = REPO_ROOT / "Dockings/Orai_Benchmark"
DEFAULT_RECEPTOR_DIR = (
    REPO_ROOT / "Dockings/Orai_Benchmark/_staging/receptors/pdbqt"
)
DEFAULT_VINA_RECEPTOR_DIR = (
    REPO_ROOT / "Dockings/Orai_Benchmark_old/_staging/receptors/pdbqt"
)
DEFAULT_IDS_FILE = (
    REPO_ROOT / "Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"
)
DEFAULT_CONFIG = REPO_ROOT / "Scripts/Docking/orai_benchmark_autodock_config.yaml"
DEFAULT_SOURCE_LOG = (
    REPO_ROOT
    / "Dockings/Logs/orai_benchmark_vinardo_logs"
    / "notebook_cell_vinardo_tmux_20260806-225202.log"
)
DEFAULT_STATUS = (
    DEFAULT_OUTPUT_DIR / "orai_benchmark_vinardo_gnina_refinement_status.json"
)
DEFAULT_VINA_STATUS = (
    DEFAULT_VINA_OUTPUT_DIR / "orai_benchmark_vina_gnina_refinement_status.json"
)
DEFAULT_LOCK = Path(
    "/tmp/master_docking_orai_benchmark_vinardo_gnina_refinement.lock"
)
DEFAULT_VINA_LOCK = Path(
    "/tmp/master_docking_orai_benchmark_vina_gnina_refinement.lock"
)
DEFAULT_GPU_LOCK = Path("/tmp/master_docking_gpu0.lock")
# gnina CNN modes this driver can commit, mapped to the run_autodock optimizer key
# they resolve to. Rescoring and refinement are separate variants that must never
# share a pose directory, log or provenance cache: rescoring minimises against the
# empirical function and uses the CNN only to rank, whereas refinement minimises
# against the CNN itself, so the two produce different geometries from the same
# input pose. Keys double as the ``--mode`` choices.
GNINA_MODES = {"refinement": "gnina_refinement", "rescore": "gnina"}

EXPECTED_RECEPTORS = 4
EXPECTED_LIGANDS = 308
EXPECTED_COMPLEXES = EXPECTED_RECEPTORS * EXPECTED_LIGANDS
EXPECTED_POSES_PER_COMPLEX = 10


class RefinementInterrupted(KeyboardInterrupt):
    """A termination signal translated into a status-aware interruption."""

    def __init__(self, signum: int):
        self.signum = signum
        super().__init__(f"received {signal.Signals(signum).name}")


def _now() -> str:
    return dt.datetime.now().astimezone().isoformat()


def _resolve(path: Path) -> Path:
    path = path.expanduser()
    if not path.is_absolute():
        path = REPO_ROOT / path
    return path.resolve()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, default=str)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    finally:
        Path(temporary_name).unlink(missing_ok=True)


class JobLock:
    """Process-lifetime, fail-fast lock with inspectable owner metadata."""

    def __init__(self, path: Path):
        self.path = path
        self.handle: Any | None = None

    def __enter__(self) -> "JobLock":
        self.path.parent.mkdir(parents=True, exist_ok=True)
        handle = self.path.open("a+", encoding="utf-8")
        try:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            handle.seek(0)
            owner = handle.read().strip()
            handle.close()
            suffix = f"; owner metadata: {owner}" if owner else ""
            raise RuntimeError(f"refinement lock is already held: {self.path}{suffix}") from exc
        self.handle = handle
        handle.seek(0)
        handle.truncate()
        json.dump(
            {
                "pid": os.getpid(),
                "process_group_id": os.getpgrp(),
                "started_at": _now(),
                "status": "running",
            },
            handle,
            indent=2,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        if self.handle is None:
            return
        try:
            fcntl.flock(self.handle.fileno(), fcntl.LOCK_UN)
        finally:
            self.handle.close()
            self.handle = None


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-scoring", choices=("vina", "vinardo"), default="vinardo"
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--receptor-dir", type=Path, default=None)
    parser.add_argument("--ids-file", type=Path, default=DEFAULT_IDS_FILE)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--source-log", type=Path, default=None)
    parser.add_argument("--status-path", type=Path, default=None)
    parser.add_argument("--log-path", type=Path, default=None)
    parser.add_argument("--lock-file", type=Path, default=None)
    parser.add_argument("--gpu-lock-file", type=Path, default=DEFAULT_GPU_LOCK)
    parser.add_argument(
        "--mode", choices=tuple(GNINA_MODES), default="refinement",
        help="gnina CNN mode. 'refinement' minimises the pose against the learned "
             "score (the original behaviour of this driver). 'rescore' minimises "
             "against the empirical function and applies the CNN only to rank, which "
             "is the mode used by the benchmark AutoDock arm and by the Orai x JKU "
             "experimental panel, so it is what the Orai x Benchmark control panel "
             "needs to be comparable with them. The two write to separate pose "
             "directories, logs and status files and never share a cache.",
    )
    parser.add_argument("--verify-only", action="store_true")
    return parser.parse_args()


def _build_refinement_config(
    config_path: Path,
    output_dir: Path,
    gpu_lock: Path,
    source_scoring: str,
    mode: str = "refinement",
) -> dict:
    with config_path.open(encoding="utf-8") as handle:
        source_cfg = yaml.safe_load(handle)
    if not isinstance(source_cfg, dict):
        raise RuntimeError(f"configuration is not a mapping: {config_path}")
    if str(source_cfg.get("scoring_function", "")).strip().lower() != "vina":
        raise RuntimeError("Orai source YAML must retain scoring_function: vina")
    if str(source_cfg.get("optimization", "")).strip().lower() != "none":
        raise RuntimeError("Orai source YAML must retain optimization: none")
    if bool(source_cfg.get("overwrite_existing") or source_cfg.get("overwrite_poses")):
        raise RuntimeError("refusing optimizer-only launch while docking overwrite is enabled")

    cfg = copy.deepcopy(source_cfg)
    use_template_reconstruction = source_scoring == "vinardo"
    cfg.update(
        {
            "output_dir": str(output_dir),
            "scoring_function": source_scoring,
            "optimization": GNINA_MODES[mode],
            "gnina_cnn_scoring": mode,
            "optimize_rank_by": "cnn_affinity",
            "optimize_search": "minimize",
            "optimize_scoring": "default",
            "optimize_autobox_add": 4.0,
            "optimize_top_n": 0,
            "optimize_cpu": 8,
            "optimize_seed": 42,
            "optimize_timeout": 1800,
            "overwrite_optimization": False,
            "gnina_use_gpu": True,
            "gnina_gpu_preflight": True,
            "gpu_lock_file": str(gpu_lock),
            "gnina_cnn_model": "crossdock_default2018_ensemble",
            "gnina_cnn_model_file": "",
            "gnina_executable": "/home/manndo/docking_tools/gnina",
            "mk_export_executable": "/home/manndo/anaconda3/envs/vina/bin/mk_export.py",
            "obabel_executable": "/usr/bin/obabel",
            "optimize_require_mk_export": not use_template_reconstruction,
            "optimize_require_template_reconstruction": use_template_reconstruction,
            "optimize_ligand_template_root": str(
                REPO_ROOT / "Data/PoseBuster Benchmark Set"
            ) if use_template_reconstruction else "",
            "optimize_prepared_ligand_pdbqt_dir": str(
                REPO_ROOT / "Dockings/Orai_Benchmark/_staging/ligands/pdbqt"
            ) if use_template_reconstruction else "",
        }
    )
    return cfg


def _read_authoritative_ids(ids_file: Path) -> list[str]:
    identifiers = sorted(
        {
            line.strip()
            for line in ids_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
    )
    if len(identifiers) != EXPECTED_LIGANDS:
        raise RuntimeError(
            f"expected {EXPECTED_LIGANDS} authoritative ligand IDs, found {len(identifiers)}"
        )
    return identifiers


def _validate_source_run(
    output_dir: Path,
    source_log: Path | None,
    source_scoring: str,
) -> None:
    summary_path = output_dir / "docking_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    overall = summary.get("overall") or {}
    configuration = summary.get("configuration") or {}
    if configuration.get("scoring_function") != source_scoring:
        raise RuntimeError(
            "docking summary scoring_function disagrees with requested source: "
            f"{configuration.get('scoring_function')!r}"
        )
    if overall.get("total_combinations") != EXPECTED_COMPLEXES:
        raise RuntimeError(
            "docking summary does not cover the expected 1232 combinations"
        )

    if source_scoring == "vinardo":
        if source_log is None:
            raise RuntimeError("Vinardo source validation requires --source-log")
        log_text = source_log.read_text(encoding="utf-8", errors="replace")
        required_markers = (
            "ORAI x BENCHMARK VINARDO DOCKING COMPLETE".replace("x", "×"),
            "Success: 1232 | Cached/skipped: 0 | Failed: 0",
            "Completed combinations: 1232/1232",
            "New Vinardo cell finished successfully",
        )
        missing = [marker for marker in required_markers if marker not in log_text]
        if missing:
            raise RuntimeError(
                f"source log lacks successful completion markers: {missing}"
            )
        expected_summary = {
            "successful": EXPECTED_COMPLEXES,
            "failed": 0,
            "skipped": 0,
            "total_poses": EXPECTED_COMPLEXES * EXPECTED_POSES_PER_COMPLEX,
        }
        mismatches = {
            key: (overall.get(key), value)
            for key, value in expected_summary.items()
            if overall.get(key) != value
        }
        if mismatches:
            raise RuntimeError(
                f"docking summary disagrees with expected completed run: {mismatches}"
            )

    suffix = "_Vinardo" if source_scoring == "vinardo" else ""
    failure_path = output_dir / f"failed_docking_Orai_Benchmark{suffix}.json"
    failures = json.loads(failure_path.read_text(encoding="utf-8"))
    if failures:
        raise RuntimeError(f"source run still records unresolved failures: {failure_path}")


def _validate_docking_csv(output_dir: Path, source_scoring: str) -> set[str]:
    suffix = "_Vinardo" if source_scoring == "vinardo" else ""
    csv_path = output_dir / f"docking_log_Orai_Benchmark{suffix}.csv"
    with csv_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != EXPECTED_COMPLEXES:
        raise RuntimeError(f"expected {EXPECTED_COMPLEXES} docking-log rows, found {len(rows)}")

    combinations: set[str] = set()
    errors: list[str] = []
    for row in rows:
        combo = str(row.get("combo_name") or "")
        if not combo or combo in combinations:
            errors.append(f"missing or duplicate combo_name {combo!r}")
        combinations.add(combo)
        if row.get("status") != "success":
            errors.append(f"{combo}: status={row.get('status')!r}")
        if str(row.get("scoring_function") or "").lower() != source_scoring:
            errors.append(f"{combo}: scoring_function={row.get('scoring_function')!r}")
        try:
            pose_count = int(row.get("num_poses") or 0)
        except ValueError:
            pose_count = -1
        if pose_count != EXPECTED_POSES_PER_COMPLEX:
            errors.append(f"{combo}: num_poses={row.get('num_poses')!r}")
    if errors:
        raise RuntimeError("invalid docking log: " + "; ".join(errors[:10]))
    return combinations


def _validate_embedded_smiles_metadata(pose_files: list[Path]) -> int:
    """Require a complete Meeko SMILES/atom map in every raw Vina model."""
    from Scripts.Docking.run_autodock import (
        _pdbqt_atom_coordinates,
        _pdbqt_atom_elements,
        _remark_smiles,
        _remark_smiles_mapping,
    )

    validated = 0
    model_pattern = re.compile(
        r"^MODEL\s+\d+\s*$\n(.*?)^ENDMDL\s*$", re.MULTILINE | re.DOTALL
    )
    for pose_file in pose_files:
        text = pose_file.read_text(encoding="utf-8")
        models = model_pattern.findall(text)
        if len(models) != EXPECTED_POSES_PER_COMPLEX:
            raise RuntimeError(
                f"{pose_file}: expected {EXPECTED_POSES_PER_COMPLEX} embedded models, "
                f"found {len(models)}"
            )
        for rank, model_text in enumerate(models, 1):
            smiles = _remark_smiles(model_text)
            molecule = Chem.MolFromSmiles(smiles) if smiles else None
            if molecule is None:
                raise RuntimeError(
                    f"{pose_file} rank {rank}: missing or invalid REMARK SMILES"
                )
            mapping = _remark_smiles_mapping(model_text)
            expected_atoms = molecule.GetNumAtoms()
            smiles_indices = [pair[0] for pair in mapping]
            pdbqt_indices = [pair[1] for pair in mapping]
            if (
                len(mapping) != expected_atoms
                or len(set(smiles_indices)) != expected_atoms
                or set(smiles_indices) != set(range(1, expected_atoms + 1))
                or len(set(pdbqt_indices)) != expected_atoms
            ):
                raise RuntimeError(
                    f"{pose_file} rank {rank}: incomplete/non-unique REMARK SMILES IDX map"
                )
            elements = _pdbqt_atom_elements(model_text)
            coordinates = _pdbqt_atom_coordinates(model_text)
            for smiles_index, pdbqt_index in mapping:
                expected_element = molecule.GetAtomWithIdx(smiles_index - 1).GetSymbol()
                if elements.get(pdbqt_index) != expected_element:
                    raise RuntimeError(
                        f"{pose_file} rank {rank}: element mismatch for mapped atom "
                        f"{smiles_index}->{pdbqt_index}"
                    )
                coordinate = coordinates.get(pdbqt_index)
                if coordinate is None or not all(
                    value == value and abs(value) != float("inf") for value in coordinate
                ):
                    raise RuntimeError(
                        f"{pose_file} rank {rank}: missing/non-finite mapped coordinate"
                    )
            validated += 1
    return validated


def _discover_and_validate_targets(
    output_dir: Path,
    receptor_dir: Path,
    authoritative_ids: list[str],
    docking_combinations: set[str],
    source_scoring: str,
) -> tuple[list[Any], dict[str, Any]]:
    from Scripts.Docking.run_autodock import discover_existing_autodock_results

    pose_suffix = (
        "_vinardo_vina_out.pdbqt"
        if source_scoring == "vinardo"
        else "_vina_vina_out.pdbqt"
    )
    pose_files = sorted((output_dir / "docking").glob(f"*{pose_suffix}"))
    if len(pose_files) != EXPECTED_COMPLEXES:
        raise RuntimeError(
            f"expected {EXPECTED_COMPLEXES} committed {source_scoring} PDBQTs, "
            f"found {len(pose_files)}"
        )
    embedded_smiles_pose_models = None
    if source_scoring == "vina":
        embedded_smiles_pose_models = _validate_embedded_smiles_metadata(pose_files)

    receptor_names = sorted({pose.stem.partition("__")[0] for pose in pose_files})
    if len(receptor_names) != EXPECTED_RECEPTORS or any(not name for name in receptor_names):
        raise RuntimeError(
            f"expected {EXPECTED_RECEPTORS} receptor identities, found {receptor_names}"
        )
    receptors = {name: receptor_dir / f"{name}.pdbqt" for name in receptor_names}
    missing_receptors = [str(path) for path in receptors.values() if not path.is_file()]
    if missing_receptors:
        raise RuntimeError(f"prepared receptor PDBQTs are missing: {missing_receptors}")

    results = discover_existing_autodock_results(output_dir / "docking", receptors)
    if len(results) != EXPECTED_COMPLEXES:
        raise RuntimeError(
            f"reconstructed {len(results)}/{EXPECTED_COMPLEXES} optimizer targets"
        )

    allowed_ids = set(authoritative_ids)
    receptor_counts: Counter[str] = Counter()
    ligand_counts: Counter[str] = Counter()
    result_combinations: set[str] = set()
    total_poses = 0
    errors: list[str] = []
    for result in results:
        receptor_counts[result.protein_name] += 1
        suffix = "_ligand_start_conf"
        ligand_id = (
            result.ligand_name[: -len(suffix)]
            if result.ligand_name.endswith(suffix)
            else ""
        )
        ligand_counts[ligand_id] += 1
        combo = f"{result.ligand_name}__{result.protein_name}"
        result_combinations.add(combo)
        total_poses += int(result.num_poses)
        if result.scoring_function != source_scoring:
            errors.append(f"{combo}: source scoring is {result.scoring_function!r}")
        if result.num_poses != EXPECTED_POSES_PER_COMPLEX:
            errors.append(f"{combo}: reconstructed {result.num_poses} poses")
        if ligand_id not in allowed_ids:
            errors.append(f"{combo}: ligand is outside the authoritative ID set")

    if result_combinations != docking_combinations:
        errors.append(
            "committed-output combinations differ from docking-log combinations "
            f"(outputs-only={len(result_combinations - docking_combinations)}, "
            f"log-only={len(docking_combinations - result_combinations)})"
        )
    if any(count != EXPECTED_LIGANDS for count in receptor_counts.values()):
        errors.append(f"per-receptor counts are not all {EXPECTED_LIGANDS}: {receptor_counts}")
    if set(ligand_counts) != allowed_ids or any(
        count != EXPECTED_RECEPTORS for count in ligand_counts.values()
    ):
        errors.append("authoritative ligands do not each occur once per receptor")
    if total_poses != EXPECTED_COMPLEXES * EXPECTED_POSES_PER_COMPLEX:
        errors.append(f"expected 12320 total poses, reconstructed {total_poses}")
    if errors:
        raise RuntimeError(
            f"raw {source_scoring} validation failed: " + "; ".join(errors[:10])
        )

    results.sort(key=lambda result: (result.protein_name, result.ligand_name))
    validation = {
        "receptors": receptor_names,
        "receptor_count": len(receptor_names),
        "receptor_sources": {
            name: {
                "path": str(receptors[name].resolve()),
                "sha256": hashlib.sha256(receptors[name].read_bytes()).hexdigest(),
            }
            for name in receptor_names
        },
        "authoritative_ligand_count": len(authoritative_ids),
        "validated_complexes": len(results),
        "validated_poses": total_poses,
    }
    if source_scoring == "vina":
        validation["receptor_provenance_note"] = (
            "The exact Orai1WT-START-Fr0 receptor PDBQT used for the June 23 "
            "Vina docking was overwritten by a later nondeterministic OpenMM "
            "preparation and is unavailable. The receptor path and SHA256 recorded "
            "above are the declared June 22 same-workflow surrogate used for "
            "refinement; the other three receptor PDBQTs are byte-identical between "
            "the June archive and current staging."
        )
    if embedded_smiles_pose_models is not None:
        validation["embedded_smiles_pose_models"] = embedded_smiles_pose_models
    return results, validation


def _protocol(cfg: Mapping[str, Any]) -> dict[str, Any]:
    source_scoring = str(cfg["scoring_function"])
    use_template_reconstruction = bool(
        cfg["optimize_require_template_reconstruction"]
    )
    protocol = {
        "source_scoring": cfg["scoring_function"],
        "optimization": cfg["optimization"],
        "cnn_scoring": cfg["gnina_cnn_scoring"],
        "cnn_model": cfg["gnina_cnn_model"],
        "rank_by": cfg["optimize_rank_by"],
        "rank_fields": ["autodock_rank", "optimized_rank"],
        "search": cfg["optimize_search"],
        "empirical_scoring": cfg["optimize_scoring"],
        "autobox_add": cfg["optimize_autobox_add"],
        "top_n": cfg["optimize_top_n"],
        "cpu": cfg["optimize_cpu"],
        "seed": cfg["optimize_seed"],
        "timeout_s": cfg["optimize_timeout"],
        "use_gpu": cfg["gnina_use_gpu"],
        "gpu_lock_file": cfg["gpu_lock_file"],
        "overwrite_optimization": cfg["overwrite_optimization"],
    }
    if use_template_reconstruction:
        protocol["pose_reconstruction"] = {
            "method": "rdkit_template_map",
            "authoritative_sdf_root": cfg["optimize_ligand_template_root"],
            "prepared_pdbqt_root": cfg["optimize_prepared_ligand_pdbqt_dir"],
            "required": True,
            "chemistry_note": (
                "GNINA uses the authoritative benchmark SDF topology and docked "
                "heavy-atom coordinates; MGLTools-added hydrogens absent from the "
                "authoritative ligand are discarded. autodock_rank remains the "
                "original Vinardo ordering."
            ),
        }
    else:
        protocol["pose_reconstruction"] = {
            "method": "mk_export_embedded_smiles",
            "required": True,
            "chemistry_note": (
                "GNINA input chemistry is reconstructed from the SMILES, atom-index, "
                "and hydrogen-parent metadata embedded in every AutoDock Vina pose. "
                f"autodock_rank remains the original {source_scoring} ordering."
            ),
        }
    return protocol


def _signal_handler(signum: int, frame: FrameType | None) -> None:
    del frame
    raise RefinementInterrupted(signum)


def main() -> int:
    args = _parse_args()
    source_scoring = str(args.source_scoring)
    mode = str(args.mode)
    optimizer_key = GNINA_MODES[mode]
    default_output_dir = (
        DEFAULT_OUTPUT_DIR if source_scoring == "vinardo" else DEFAULT_VINA_OUTPUT_DIR
    )
    default_status = DEFAULT_STATUS if source_scoring == "vinardo" else DEFAULT_VINA_STATUS
    default_lock = DEFAULT_LOCK if source_scoring == "vinardo" else DEFAULT_VINA_LOCK
    # Rescoring is a separate variant from refinement, so it gets its own status file
    # and its own process lock. Sharing either would let a rescore run resume from a
    # refinement checkpoint, or block on a refinement run that is doing different work.
    if mode != "refinement":
        default_status = default_status.with_name(
            default_status.name.replace("gnina_refinement", f"gnina_{mode}")
        )
        default_lock = default_lock.with_name(
            default_lock.name.replace("gnina_refinement", f"gnina_{mode}")
        )
    output_dir = _resolve(args.output_dir or default_output_dir)
    default_receptor_dir = (
        DEFAULT_RECEPTOR_DIR
        if source_scoring == "vinardo"
        else DEFAULT_VINA_RECEPTOR_DIR
    )
    receptor_dir = _resolve(args.receptor_dir or default_receptor_dir)
    ids_file = _resolve(args.ids_file)
    config_path = _resolve(args.config)
    source_log = (
        _resolve(args.source_log)
        if args.source_log
        else (DEFAULT_SOURCE_LOG.resolve() if source_scoring == "vinardo" else None)
    )
    status_path = _resolve(args.status_path or default_status)
    lock_file = (args.lock_file or default_lock).expanduser().resolve()
    gpu_lock = args.gpu_lock_file.expanduser().resolve()
    log_path = _resolve(args.log_path) if args.log_path else None

    cfg = _build_refinement_config(
        config_path, output_dir, gpu_lock, source_scoring, mode
    )
    from Scripts.Docking.run_autodock import (
        OptimizationError,
        optimize_autodock_results,
        preflight_optimizers,
        preflight_mgl_pose_templates,
        resolve_optimizers,
    )

    if resolve_optimizers(cfg) != [optimizer_key]:
        raise RuntimeError(
            f"protocol did not resolve exclusively to {optimizer_key}"
        )

    def validate() -> tuple[list[Any], dict[str, Any], dict[str, Any]]:
        _validate_source_run(output_dir, source_log, source_scoring)
        authoritative_ids = _read_authoritative_ids(ids_file)
        docking_combinations = _validate_docking_csv(output_dir, source_scoring)
        targets, validation = _discover_and_validate_targets(
            output_dir,
            receptor_dir,
            authoritative_ids,
            docking_combinations,
            source_scoring,
        )
        if cfg["optimize_require_template_reconstruction"]:
            templates = preflight_mgl_pose_templates(
                (target.ligand_name for target in targets), cfg,
            )
            if len(templates) != EXPECTED_LIGANDS:
                raise RuntimeError(
                    f"validated {len(templates)}/{EXPECTED_LIGANDS} authoritative "
                    "ligand template mappings"
                )
            retained_stereo_h = sorted(
                ligand_name
                for ligand_name, context in templates.items()
                if any(
                    atom.GetAtomicNum() == 1
                    for atom in context.molecule.GetAtoms()
                )
            )
            validation["pose_reconstruction"] = {
                "method": "rdkit_template_map",
                "validated_ligand_templates": len(templates),
                "retained_stereo_hydrogen_ligands": retained_stereo_h,
                "required": True,
            }
        else:
            embedded_pose_models = int(
                validation.get("embedded_smiles_pose_models") or 0
            )
            if embedded_pose_models != EXPECTED_COMPLEXES * EXPECTED_POSES_PER_COMPLEX:
                raise RuntimeError(
                    "embedded Meeko metadata validation did not cover all 12320 poses"
                )
            validation["pose_reconstruction"] = {
                "method": "mk_export_embedded_smiles",
                "validated_pose_models": embedded_pose_models,
                "required": True,
            }
        preflight = preflight_optimizers(cfg)
        return targets, validation, preflight

    if args.verify_only:
        targets, validation, preflight = validate()
        print(
            json.dumps(
                {
                    "status": "ready",
                    "output_dir": str(output_dir),
                    "source_log": str(source_log) if source_log else None,
                    "status_path": str(status_path),
                    "protocol": _protocol(cfg),
                    "validation": validation,
                    "preflight": preflight,
                    "target_complexes": len(targets),
                },
                indent=2,
                default=str,
            )
        )
        return 0

    started_at = _now()
    totals = {
        "selected_poses": 0,
        "attempted": 0,
        "optimized": 0,
        "reused": 0,
        "failed": 0,
        "pruned_outputs": 0,
    }
    issues: list[dict[str, Any]] = []
    processed = 0
    validation: dict[str, Any] = {}

    def write_status(state: str, *, last_complex: str | None = None) -> None:
        _atomic_write_json(
            status_path,
            {
                "schema_version": 1,
                "stage": f"orai_benchmark_{source_scoring}_{optimizer_key}",
                "status": state,
                "started_at": started_at,
                "updated_at": _now(),
                "pid": os.getpid(),
                "process_group_id": os.getpgrp(),
                "log_path": str(log_path) if log_path else None,
                "source_log": str(source_log) if source_log else None,
                "output_dir": str(output_dir),
                "optimization_log": str(output_dir / "optimization_log.csv"),
                "optimized_pose_dir": str(
                    output_dir / f"docking/optimized_{optimizer_key}"
                ),
                "protocol": _protocol(cfg),
                "validation": validation,
                "authoritative_complexes": EXPECTED_COMPLEXES,
                "processed_complexes": processed,
                "last_complex": last_complex,
                "optimizer_totals": dict(totals),
                "issues": list(issues),
            },
        )

    with JobLock(lock_file):
        try:
            targets, validation, preflight = validate()
            validation["preflight"] = preflight
        except Exception as exc:
            issues.append({"stage": "validation_or_preflight", "error": str(exc)})
            write_status("failed")
            raise

        print(
            f"Validated {validation['validated_complexes']} {source_scoring} complexes / "
            f"{validation['validated_poses']} poses; starting GNINA CNN {mode}.",
            flush=True,
        )
        print(
            "Original autodock_rank is retained; optimized_rank uses cnn_affinity. "
            f"GPU work is serialized through {gpu_lock}.",
            flush=True,
        )
        write_status("running")

        previous_handlers = {
            signum: signal.getsignal(signum) for signum in (signal.SIGINT, signal.SIGTERM)
        }
        for signum in previous_handlers:
            signal.signal(signum, _signal_handler)
        try:
            for index, result in enumerate(targets, 1):
                combo = f"{result.protein_name}__{result.ligand_name}"
                try:
                    summary = optimize_autodock_results(
                        [result], cfg, output_dir, overwrite=False
                    )
                except OptimizationError as exc:
                    summary = exc.summary
                except Exception as exc:
                    processed += 1
                    issues.append(
                        {"complex": combo, "status": "failed", "errors": [str(exc)]}
                    )
                    write_status("running", last_complex=combo)
                    print(f"[{index}/{EXPECTED_COMPLEXES}] FAILED {combo}: {exc}", flush=True)
                    continue

                for key in totals:
                    totals[key] += int(getattr(summary, key, 0) or 0)
                if summary.status != "completed":
                    issues.append({"complex": combo, **summary.to_dict()})
                processed += 1
                write_status("running", last_complex=combo)
                print(
                    f"[{index}/{EXPECTED_COMPLEXES}] {combo}: {summary.status}; "
                    f"optimized={summary.optimized}, reused={summary.reused}, "
                    f"failed={summary.failed}",
                    flush=True,
                )
        except RefinementInterrupted as exc:
            issues.append({"stage": "driver", "error": str(exc)})
            write_status("interrupted")
            raise
        except BaseException as exc:
            issues.append({"stage": "driver", "error": str(exc)})
            write_status("failed")
            raise
        finally:
            for signum, handler in previous_handlers.items():
                signal.signal(signum, handler)

        failed = bool(issues or processed != EXPECTED_COMPLEXES or totals["failed"])
        write_status("failed" if failed else "completed")
        print(f"Refinement status: {status_path}", flush=True)
        if failed:
            raise RuntimeError(
                f"GNINA refinement finished with {len(issues)} recorded issue(s); "
                f"see {status_path}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
