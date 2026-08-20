#!/usr/bin/env python3
"""Whole-protein Uni-Dock2 (GPU) docking of the Benchmark set — one shot per complex.

Uni-Dock2 (DP Technology, v0.6.3) is a ground-up rewrite of Uni-Dock. Unlike the
Uni-Dock-1.x engine — whose GPU affinity grid is capped at 128 pts/axis and 531441
pts total, forcing the whole-protein boxes here to be covered by overlapping 0.375 Å
sub-boxes (see ``run_unidock_tiled.py``) — Uni-Dock2 manages its grid internally and
docks a whole-protein box in a SINGLE GPU call. The benchmark's longest-axis case
(7M31_TDR, 135x111x192 Å) peaked at only ~140 MiB of GPU memory; 8F4J_PHO is
larger by volume. This driver therefore does one ``unidock2 docking`` call per complex.

Interface differences vs Uni-Dock 1.x:
  * Inputs are a raw PDB receptor + an SDF ligand (not PDBQT). Uni-Dock2 runs its own
    receptor preparation (pdbfixer + AmberTools) internally, so we feed it the clean
    original protein PDB and ligand start-conformer SDF that the AutoDock Vina
    full-protein cell staged, preserving the benchmark IDs and boxes. Uni-Dock2 may
    normalize/remove receptor atoms; the exact prepared receptor scored is exported.
  * The box is given as a centre (-c) plus a ``box_size: [x, y, z]`` in the YAML config.
  * Output is an SDF; the Vina affinity of each pose is stored in the SDF property
    ``<vina_binding_free_energy>`` (lower = better).

Cost model (RTX 4070 Laptop): the GPU docking itself is ~2 s regardless of box size;
per-complex WALL-TIME is dominated by CPU-side receptor preparation, which scales with
receptor size — ~14 s at ~1.3k input protein heavy atoms, ~85 s near the ~3.6k median,
and >10 min for the largest inputs (up to ~41.6k heavy ATOM records). Two guards keep
a 308-complex batch robust:
  * ``timeout_s`` — a per-complex wall-clock limit; a complex that exceeds it is killed
    (whole process group) and recorded as ``timeout`` so the batch keeps going.
  * ``max_receptor_atoms`` — optionally pre-skip by raw protein heavy-ATOM count;
    this is a pre-preparation cost proxy, not the final engine atom count.

Prep/dock split (keeps the idle GPU fed): because that CPU prep is single-threaded and
GPU-free, it is decoupled from the GPU docking. Each receptor is prepared once with the
``unidock2 protein_prep`` sub-command — pdbfixer + AmberTools ``tleap`` — into a
parameterized ``.dms`` receptor; the fast ``unidock2 docking`` shot then consumes that
``.dms`` (its ``-r`` accepts PDB *or* DMS) and skips its own internal prep entirely. The
two produce a byte-identical scored receptor and pose SDF, so this is purely a scheduling
change, not a scientific one. ``prep_workers`` (config) receptors are prepared CONCURRENTLY
on the CPU while a single GPU-locked consumer docks the already-prepared ones in submission
order — so the GPU runs its ~2 s shots back-to-back instead of once every prep-time. The
exported prepared-receptor PDB is materialized (via ``ambpdb``) from the prep step's
``receptor.prmtop``/``receptor.inpcrd``, because ``docking`` from a ``.dms`` no longer
regenerates them; the completion manifest, fingerprint, and runtime attestation are
otherwise unchanged, so existing completions still resume and re-docks stay bit-identical.

Runs are RESUMABLE: only a strictly validated, scored SDF with a matching JSON
completion sentinel (input/config/tool/output hashes) is trusted. Docking writes
to a generation-specific temporary SDF and atomically publishes it only after a
successful process exit and full validation.

All knobs live in ``unidock2_docking_config_full_protein.yaml``.
"""
from __future__ import annotations

import ast
import contextlib
import errno
import fcntl
import hashlib
import json
import math
import os
import re
import shlex
import signal
import shutil
import stat
import subprocess
import tempfile
import threading
import time
import uuid
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import yaml

# Uni-Dock2 stores each pose's Vina affinity in this SDF property (lower = better).
_AFFINITY_TAG = "vina_binding_free_energy"
_AFFINITY_RE = re.compile(r"^>\s*<%s>" % _AFFINITY_TAG)
_SENTINEL_SCHEMA = 3
_ID_RE = re.compile(r"^[0-9][A-Z0-9]{3}_[A-Z0-9]{3}$")
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_DONE_SUFFIX = "_unidock2_completion.json"
# CPU receptor-prep concurrency (GPU docking is always serial). The default keeps
# well under core count so the few large-receptor preps do not starve the machine.
_DEFAULT_PREP_WORKERS = 8
_MAX_PREP_WORKERS = 64
# Measured per-complex GPU docking cost once the receptor is pre-prepared (a ~2 s
# GPU shot plus process startup). Used only for the wall-time estimate.
_EST_DOCK_SECONDS = 5.0
_CONFIG_KEYS = {
    "unidock2_bin", "prep_base", "converter", "ids_file", "output_dir", "log_dir",
    "overwrite", "exhaustiveness", "mc_steps", "num_pose", "energy_range", "seed",
    "rmsd_limit", "randomize", "opt_steps", "refine_steps", "use_tor_lib",
    "energy_decomp", "construct_ff", "template_docking", "compute_center",
    "covalent_ligand", "preserve_receptor_hydrogen", "engine_checkpoint",
    "search_mode", "task", "gpu_device_id", "gpu_lock_file", "gpu_lock_timeout_s",
    "n_cpu", "prep_workers", "timeout_s", "terminate_grace_s", "scratch_dir",
    "max_receptor_atoms", "exclude", "path_base", "_config_path", "_path_base",
}


# ────────────────────────────────────────────────────────────────────────────
# Small parsers
# ────────────────────────────────────────────────────────────────────────────
def read_box(box_file: Path):
    """Parse an AutoDock '<key> = <value>' box.txt -> (center, size) 3-tuples."""
    v = {}
    for line_number, line in enumerate(Path(box_file).read_text().splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if "=" not in line:
            raise ValueError(f"box file {box_file}:{line_number} is not '<key> = <value>'")
        k, val = line.split("=", 1)
        key = k.strip()
        if key in v:
            raise ValueError(f"box file {box_file} repeats {key!r}")
        try:
            v[key] = float(val.strip())
        except ValueError as exc:
            raise ValueError(f"box file {box_file}:{line_number} has a non-numeric value") from exc
    required = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
    missing = [key for key in required if key not in v]
    if missing:
        raise ValueError(f"box file {box_file} is missing: {', '.join(missing)}")
    values = [v[key] for key in required]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"box file {box_file} contains a non-finite value")
    if any(v[key] <= 0 for key in ("size_x", "size_y", "size_z")):
        raise ValueError(f"box file {box_file} contains a non-positive size")
    return (v["center_x"], v["center_y"], v["center_z"]), (v["size_x"], v["size_y"], v["size_z"])


def receptor_heavy_atoms(pdb: Path) -> int:
    """Count heavy ATOM records in the staged raw protein PDB.

    This is deliberately named and reported as a *raw* count. Uni-Dock2 can
    remove/normalize atoms during preparation; its in-box count is separately
    attested from the engine log after a run.
    """
    n = 0
    for line in Path(pdb).read_text(errors="ignore").splitlines():
        if line.startswith("ATOM"):
            atom_name = re.sub(r"^[0-9]+", "", line[12:16].strip()).upper()
            element = (line[76:78].strip() or atom_name[:1]).upper()
            if element not in {"H", "D", "T"}:
                n += 1
    return n


def _validate_sdf_record(
    lines, path: Path, record_number: int, require_affinity: bool = True
) -> float | None:
    if len(lines) < 5:
        raise ValueError(f"{path}: SDF record {record_number} is too short")
    mend_indices = [i for i, line in enumerate(lines) if line.strip() == "M  END"]
    if len(mend_indices) != 1:
        raise ValueError(
            f"{path}: SDF record {record_number} has {len(mend_indices)} M  END lines; expected 1"
        )
    mend = mend_indices[0]

    # Validate a non-empty V2000 or V3000 structure without requiring RDKit in
    # the notebook environment.
    counts = lines[3] if len(lines) > 3 else ""
    if "V2000" in counts:
        try:
            atom_count = int(counts[:3])
            bond_count = int(counts[3:6])
        except ValueError as exc:
            raise ValueError(f"{path}: record {record_number} has an invalid V2000 counts line") from exc
        if atom_count < 1 or bond_count < 0 or mend < 4 + atom_count + bond_count:
            raise ValueError(f"{path}: record {record_number} has a truncated V2000 mol block")
        for atom_line in lines[4:4 + atom_count]:
            try:
                xyz = (float(atom_line[0:10]), float(atom_line[10:20]), float(atom_line[20:30]))
            except ValueError as exc:
                raise ValueError(f"{path}: record {record_number} has invalid atom coordinates") from exc
            if not all(math.isfinite(value) for value in xyz) or not atom_line[31:34].strip():
                raise ValueError(f"{path}: record {record_number} has an invalid atom line")
        for bond_line in lines[4 + atom_count:4 + atom_count + bond_count]:
            try:
                first = int(bond_line[0:3])
                second = int(bond_line[3:6])
                order = int(bond_line[6:9])
            except ValueError as exc:
                raise ValueError(f"{path}: record {record_number} has an invalid bond line") from exc
            if not (1 <= first <= atom_count and 1 <= second <= atom_count and order >= 1):
                raise ValueError(f"{path}: record {record_number} has an invalid bond reference")
    elif "V3000" in counts:
        count_line = next((line for line in lines if "M  V30 COUNTS" in line), None)
        if count_line is None:
            raise ValueError(f"{path}: record {record_number} has no V3000 COUNTS line")
        try:
            count_fields = count_line.split("COUNTS", 1)[1].split()
            atom_count = int(count_fields[0])
            bond_count = int(count_fields[1])
        except (IndexError, ValueError) as exc:
            raise ValueError(f"{path}: record {record_number} has invalid V3000 counts") from exc
        begin_atom = next((i for i, line in enumerate(lines) if "M  V30 BEGIN ATOM" in line), None)
        end_atom = next((i for i, line in enumerate(lines) if "M  V30 END ATOM" in line), None)
        if begin_atom is None or end_atom is None or end_atom <= begin_atom:
            raise ValueError(f"{path}: record {record_number} has a truncated V3000 atom block")
        atom_lines = [line for line in lines[begin_atom + 1:end_atom] if "M  V30" in line]
        if atom_count < 1 or bond_count < 0 or len(atom_lines) != atom_count:
            raise ValueError(f"{path}: record {record_number} has inconsistent V3000 atom counts")
        for atom_line in atom_lines:
            fields = atom_line.split()
            try:
                xyz = tuple(float(value) for value in fields[4:7])
            except (IndexError, ValueError) as exc:
                raise ValueError(f"{path}: record {record_number} has an invalid V3000 atom") from exc
            if len(fields) < 7 or len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
                raise ValueError(f"{path}: record {record_number} has an invalid V3000 atom")
        begin_bond = next((i for i, line in enumerate(lines) if "M  V30 BEGIN BOND" in line), None)
        end_bond = next((i for i, line in enumerate(lines) if "M  V30 END BOND" in line), None)
        if bond_count:
            if begin_bond is None or end_bond is None or end_bond <= begin_bond:
                raise ValueError(f"{path}: record {record_number} has a truncated V3000 bond block")
            bond_lines = [line for line in lines[begin_bond + 1:end_bond] if "M  V30" in line]
            if len(bond_lines) != bond_count:
                raise ValueError(f"{path}: record {record_number} has inconsistent V3000 bond counts")
    else:
        raise ValueError(f"{path}: record {record_number} has no recognized molfile counts line")
    if atom_count < 1:
        raise ValueError(f"{path}: record {record_number} has no atoms")
    if mend <= 3:
        raise ValueError(f"{path}: record {record_number} has a truncated mol block")

    values = []
    for i, line in enumerate(lines):
        if _AFFINITY_RE.match(line):
            if i + 1 >= len(lines) or not lines[i + 1].strip():
                raise ValueError(f"{path}: record {record_number} has an empty affinity")
            try:
                value = float(lines[i + 1].strip())
            except ValueError as exc:
                raise ValueError(f"{path}: record {record_number} has an invalid affinity") from exc
            if not math.isfinite(value):
                raise ValueError(f"{path}: record {record_number} has a non-finite affinity")
            values.append(value)
    if require_affinity and len(values) != 1:
        raise ValueError(
            f"{path}: record {record_number} has {len(values)} {_AFFINITY_TAG} fields; expected 1"
        )
    if not require_affinity and values:
        raise ValueError(f"{path}: input ligand unexpectedly contains output affinity fields")
    return values[0] if values else None


def parse_sdf_scores(sdf: Path):
    """Return one finite Vina affinity per complete, structurally valid SDF record."""
    try:
        lines = Path(sdf).read_text(errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read output {sdf}: {exc}") from exc
    records, current = [], []
    for line in lines:
        if line.strip() == "$$$$":
            if not any(part.strip() for part in current):
                raise ValueError(f"{sdf}: empty SDF record")
            records.append(current)
            current = []
        else:
            current.append(line)
    if any(line.strip() for line in current):
        raise ValueError(f"{sdf}: truncated SDF record without $$$$")
    if not records:
        raise ValueError(f"{sdf}: no complete SDF records")
    return [_validate_sdf_record(record, Path(sdf), i) for i, record in enumerate(records, start=1)]


def validate_sdf(sdf: Path, max_poses: int | None = None):
    """Strictly validate every SDF record and its finite Vina affinity."""
    scores = parse_sdf_scores(sdf)
    if max_poses is not None and len(scores) > max_poses:
        raise ValueError(f"{sdf}: {len(scores)} poses exceeds requested maximum {max_poses}")
    return {"num_poses": len(scores), "best_affinity": min(scores), "affinities": scores}


def validate_single_input_ligand(sdf: Path):
    """Validate the staged ligand as exactly one V2000/V3000 molecule."""
    try:
        lines = Path(sdf).read_text(errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read input ligand {sdf}: {exc}") from exc
    delimiter_indices = [i for i, line in enumerate(lines) if line.strip() == "$$$$"]
    if len(delimiter_indices) > 1:
        raise ValueError(f"{sdf}: expected one input molecule, found multiple delimiters")
    if delimiter_indices:
        delimiter = delimiter_indices[0]
        if any(line.strip() for line in lines[delimiter + 1:]):
            raise ValueError(f"{sdf}: data follows the sole molecule delimiter")
        lines = lines[:delimiter]
    _validate_sdf_record(lines, Path(sdf), 1, require_affinity=False)
    return True


def _validate_docking_output(sdf: Path, cfg: dict):
    validation = validate_sdf(sdf, max_poses=int(cfg["num_pose"]))
    scores = validation["affinities"]
    tolerance = 0.002  # Uni-Dock2 writes affinities rounded to three decimals.
    if any(later < earlier - tolerance for earlier, later in zip(scores, scores[1:])):
        raise ValueError(f"{sdf}: output affinities are not ordered best-to-worst")
    if scores[-1] > scores[0] + float(cfg["energy_range"]) + tolerance:
        raise ValueError(f"{sdf}: output exceeds the configured energy_range")
    return validation


def best_affinity(sdf: Path):
    """Return the minimum validated <vina_binding_free_energy>."""
    return min(parse_sdf_scores(sdf))


def count_poses(sdf: Path) -> int:
    """Number of complete, structurally valid, scored pose records."""
    try:
        return len(parse_sdf_scores(sdf))
    except (OSError, ValueError):
        return 0


def _valid_out(sdf: Path, max_poses: int | None = None) -> bool:
    try:
        validate_sdf(sdf, max_poses=max_poses)
        return True
    except (OSError, ValueError):
        return False


def _atomic_write_text(path: Path, text: str):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _atomic_write_json(path: Path, value):
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _python_tree_provenance(package_root: Path):
    """Deterministically fingerprint installed Python sources imported by the CLI."""
    package_root = Path(package_root).resolve()
    sources = sorted(
        path for path in package_root.rglob("*.py")
        if "__pycache__" not in path.parts and path.is_file()
    )
    if not sources:
        raise RuntimeError(f"no Python sources found under {package_root}")
    digest = hashlib.sha256()
    total_size = 0
    for source in sources:
        relative = source.relative_to(package_root).as_posix().encode("utf-8")
        content_hash = _sha256(source)
        size = source.stat().st_size
        total_size += size
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(size.to_bytes(8, "big"))
        digest.update(bytes.fromhex(content_hash))
    return {
        "path": str(package_root), "file_count": len(sources),
        "size": total_size, "sha256": digest.hexdigest(),
    }


def _read_json(path: Path):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def validate_config(cfg: dict) -> dict:
    if not isinstance(cfg, dict):
        raise ValueError("configuration must be a YAML mapping")
    unknown = sorted(set(cfg) - _CONFIG_KEYS)
    if unknown:
        raise ValueError(f"unknown configuration keys: {', '.join(unknown)}")
    required = (
        "unidock2_bin", "prep_base", "converter", "ids_file", "output_dir", "log_dir",
        "exhaustiveness", "mc_steps", "num_pose", "energy_range", "seed",
        "randomize", "opt_steps", "refine_steps", "use_tor_lib", "energy_decomp",
        "construct_ff", "template_docking", "compute_center", "covalent_ligand",
        "preserve_receptor_hydrogen", "engine_checkpoint", "search_mode", "task",
        "gpu_device_id", "gpu_lock_file", "n_cpu", "timeout_s",
    )
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"missing configuration keys: {', '.join(missing)}")

    def require_int(key, minimum=None):
        value = cfg[key]
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"{key} must be an integer")
        if minimum is not None and value < minimum:
            raise ValueError(f"{key} must be at least {minimum}")
        return value

    for key in ("exhaustiveness", "mc_steps", "num_pose", "timeout_s", "n_cpu"):
        require_int(key, 1)
    require_int("seed")
    require_int("gpu_device_id", 0)
    require_int("opt_steps", -1)
    require_int("refine_steps", 0)
    # Concurrent CPU receptor-prep workers (GPU docking stays serial on the lock).
    # Optional; absent means the split-pipeline default.
    workers = cfg.get("prep_workers", _DEFAULT_PREP_WORKERS)
    if isinstance(workers, bool) or not isinstance(workers, int) or workers < 1:
        raise ValueError("prep_workers must be an integer of at least 1")
    if workers > _MAX_PREP_WORKERS:
        raise ValueError(f"prep_workers must not exceed {_MAX_PREP_WORKERS}")
    cfg["prep_workers"] = workers
    for key in (
        "randomize", "use_tor_lib", "energy_decomp", "construct_ff",
        "template_docking", "compute_center", "covalent_ligand",
        "preserve_receptor_hydrogen", "engine_checkpoint",
    ):
        if not isinstance(cfg[key], bool):
            raise ValueError(f"{key} must be true or false")
    if cfg.get("max_receptor_atoms") is not None:
        if isinstance(cfg["max_receptor_atoms"], bool) or not isinstance(
            cfg["max_receptor_atoms"], int
        ) or cfg["max_receptor_atoms"] < 1:
            raise ValueError("max_receptor_atoms must be null or an integer of at least 1")
    for key in ("energy_range",):
        value = cfg[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{key} must be numeric")
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be a finite non-negative number")
    if cfg.get("rmsd_limit") is not None:
        value = cfg["rmsd_limit"]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("rmsd_limit must be null or numeric")
        value = float(value)
        if not math.isfinite(value) or value < 0:
            raise ValueError("rmsd_limit must be null or a finite non-negative number")
    if not isinstance(cfg.get("overwrite", False), bool):
        raise ValueError("overwrite must be true or false")
    exclude = cfg.get("exclude") or []
    if not isinstance(exclude, list):
        raise ValueError("exclude must be a YAML list")
    invalid_exclude = [value for value in exclude if not isinstance(value, str) or not _ID_RE.fullmatch(value)]
    if invalid_exclude:
        raise ValueError(f"exclude contains invalid benchmark IDs: {invalid_exclude!r}")
    if len(exclude) != len(set(exclude)):
        raise ValueError("exclude contains duplicate benchmark IDs")
    for key in ("search_mode", "task", "converter"):
        if not isinstance(cfg[key], str) or not cfg[key].strip():
            raise ValueError(f"{key} must be a non-empty string")
    if cfg["search_mode"] != "free":
        raise ValueError(
            "search_mode must be 'free' so configured exhaustiveness/mc_steps are honored"
        )
    if cfg["task"] != "screen":
        raise ValueError("task must be 'screen' for this benchmark driver")
    for key in ("unidock2_bin", "prep_base", "output_dir", "log_dir", "gpu_lock_file"):
        if not isinstance(cfg[key], str) or not cfg[key].strip():
            raise ValueError(f"{key} must be a non-empty path string")
    if not isinstance(cfg["ids_file"], str) or not cfg["ids_file"].strip():
        raise ValueError("ids_file must be a non-empty path string")
    if cfg.get("scratch_dir") is not None and not re.fullmatch(
        r"[A-Za-z0-9_./-]+", cfg["scratch_dir"]
    ):
        raise ValueError(
            "scratch_dir may contain only letters, digits, '_', '-', '.', and '/' "
            "because Uni-Dock2/AmberTools does not shell-quote it"
        )
    grace = cfg.get("terminate_grace_s", 5.0)
    if isinstance(grace, bool) or not isinstance(grace, (int, float)):
        raise ValueError("terminate_grace_s must be numeric")
    if not math.isfinite(float(grace)) or float(grace) <= 0:
        raise ValueError("terminate_grace_s must be finite and positive")
    lock_timeout = cfg.get(
        "gpu_lock_timeout_s", float(cfg["timeout_s"]) + float(grace) + 60.0
    )
    if isinstance(lock_timeout, bool) or not isinstance(lock_timeout, (int, float)):
        raise ValueError("gpu_lock_timeout_s must be numeric")
    if not math.isfinite(float(lock_timeout)) or float(lock_timeout) <= 0:
        raise ValueError("gpu_lock_timeout_s must be finite and positive")
    cfg["gpu_lock_timeout_s"] = float(lock_timeout)
    return cfg


def probe_tool_identity(unidock2_bin: str):
    binary = Path(unidock2_bin).resolve()
    env_prefix = binary.parent.parent
    if not binary.is_file() or not os.access(binary, os.X_OK):
        raise FileNotFoundError(f"Uni-Dock2 binary is missing or not executable: {binary}")
    tleap = binary.parent / "tleap"
    if not tleap.is_file() or not os.access(tleap, os.X_OK):
        raise FileNotFoundError(
            f"AmberTools tleap is missing or not executable beside Uni-Dock2: {tleap}"
        )
    ambpdb = binary.parent / "ambpdb"
    if not ambpdb.is_file() or not os.access(ambpdb, os.X_OK):
        raise FileNotFoundError(
            f"AmberTools ambpdb is missing or not executable beside Uni-Dock2: {ambpdb}"
        )
    teleap = binary.parent / "teLeap"
    if not teleap.is_file() or not os.access(teleap, os.X_OK):
        raise FileNotFoundError(f"AmberTools teLeap binary is missing or not executable: {teleap}")

    def one_match(pattern: str, label: str):
        matches = sorted(env_prefix.glob(pattern))
        unique = {}
        for path in matches:
            stat = path.stat()
            unique.setdefault((stat.st_dev, stat.st_ino), path.resolve())
        matches = sorted(unique.values())
        if len(matches) != 1:
            raise RuntimeError(
                f"expected exactly one {label} matching {pattern!r}; found {len(matches)}"
            )
        return matches[0]

    engine_library = one_match(
        "lib/python*/site-packages/unidock_engine/api/python/pipeline*.so",
        "Uni-Dock2 engine library",
    )
    processing_package = one_match(
        "lib/python*/site-packages/unidock_processing",
        "unidock_processing Python package",
    )
    package_records = [
        one_match("conda-meta/unidock2-*.json", "Uni-Dock2 conda record"),
        one_match("conda-meta/ambertools_stable-*.json", "AmberTools conda record"),
        one_match("conda-meta/pdbfixer-*.json", "pdbfixer conda record"),
    ]
    try:
        proc = subprocess.run(
            [str(binary), "--version"], capture_output=True, text=True,
            timeout=15, env=_subprocess_env(str(binary)),
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError("Uni-Dock2 --version timed out after 15 seconds") from exc
    if proc.returncode != 0:
        raise RuntimeError(
            f"could not determine Uni-Dock2 version (rc={proc.returncode}): "
            f"{proc.stderr.strip() or proc.stdout.strip()}"
        )
    match = re.search(r"unidock2\s+([^\s]+)", f"{proc.stdout}\n{proc.stderr}", re.IGNORECASE)
    if not match:
        raise RuntimeError("could not parse Uni-Dock2 --version output")
    artifacts = [binary, tleap, teleap, ambpdb, engine_library, *package_records]
    return {
        "path": str(binary), "version": match.group(1),
        "launcher_sha256": _sha256(binary), "tleap_path": str(tleap),
        "ambpdb_path": str(ambpdb),
        "artifacts": [_file_provenance(path) for path in artifacts],
        "python_source_tree": _python_tree_provenance(processing_package),
    }


# ────────────────────────────────────────────────────────────────────────────
# Per-complex docking (one Uni-Dock2 call, timeout-guarded, resumable)
# ────────────────────────────────────────────────────────────────────────────
def _subprocess_env(unidock2_bin: str) -> dict:
    """Environment for the Uni-Dock2 subprocess. Its receptor prep shells out to
    AmberTools `tleap`, which lives in the unidock2 env's bin — so that bin must be
    on PATH even when this driver runs from a DIFFERENT conda env (the notebook
    kernel). Prepend it and set AMBERHOME/CONDA_PREFIX so the run matches a normal
    `conda activate unidock2`. (Shared libs already load via RPATH.)"""
    env_bin = Path(unidock2_bin).resolve().parent
    env_prefix = env_bin.parent
    e = os.environ.copy()
    e["PATH"] = f"{env_bin}{os.pathsep}{e.get('PATH', '')}"
    e["AMBERHOME"] = str(env_prefix)
    e["CONDA_PREFIX"] = str(env_prefix)
    return e


# Prep subprocesses run on a background thread pool, so KeyboardInterrupt is
# delivered to the main thread and cannot unwind their proc.wait(). Track every
# live prep process group here so the main thread can reap them all on shutdown.
_ACTIVE_PREP_LOCK = threading.Lock()
_ACTIVE_PREP: dict[int, subprocess.Popen] = {}


def _register_prep(proc: subprocess.Popen):
    with _ACTIVE_PREP_LOCK:
        _ACTIVE_PREP[id(proc)] = proc


def _deregister_prep(proc: subprocess.Popen):
    with _ACTIVE_PREP_LOCK:
        _ACTIVE_PREP.pop(id(proc), None)


def _kill_active_prep(grace_s: float = 2.0):
    """Reap every still-running prep process group (used on batch shutdown)."""
    with _ACTIVE_PREP_LOCK:
        procs = list(_ACTIVE_PREP.values())
    for proc in procs:
        try:
            _terminate_process_group(proc, grace_s=grace_s)
        except Exception:
            pass


def _run(cmd, timeout_s, log_path: Path, env: dict, terminate_grace_s: float = 5.0,
         on_spawn=None, on_reap=None):
    """Run in a new process group and always reap it after timeout/interruption.

    The log is append-only so every retry remains auditable. The return contract
    stays ``(returncode_or_None, timed_out)`` for notebook/test compatibility.
    ``on_spawn``/``on_reap`` (optional) receive the Popen so a caller running on a
    worker thread can register/deregister it for main-thread shutdown reaping.
    """
    log_path.parent.mkdir(parents=True, exist_ok=True)
    grace_s = float(terminate_grace_s)
    with open(log_path, "a", encoding="utf-8") as log:
        log.write(
            f"\n===== Uni-Dock2 attempt {time.strftime('%Y-%m-%dT%H:%M:%S%z')} "
            f"pid={os.getpid()} =====\n$ {shlex.join([str(c) for c in cmd])}\n\n"
        )
        log.flush()
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT,
                                 start_new_session=True, env=env)
        if on_spawn is not None:
            on_spawn(proc)
        try:
            return proc.wait(timeout=timeout_s), False
        except subprocess.TimeoutExpired:
            log.write(f"\n[driver] timeout after {timeout_s}s; terminating process group\n")
            log.flush()
            _terminate_process_group(proc, grace_s=grace_s)
            return None, True
        except BaseException:
            log.write("\n[driver] interrupted; terminating process group\n")
            log.flush()
            _terminate_process_group(proc, grace_s=grace_s)
            raise
        finally:
            if on_reap is not None:
                on_reap(proc)


def _process_group_exists(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
    except ProcessLookupError:
        return False
    return True


def _terminate_process_group(proc: subprocess.Popen, grace_s: float = 5.0):
    """Gracefully TERM, then KILL, a detached group and reap its leader."""
    pgid = proc.pid  # start_new_session=True makes the child its process-group leader.
    if _process_group_exists(pgid):
        try:
            os.killpg(pgid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + grace_s
    while _process_group_exists(pgid) and time.monotonic() < deadline:
        proc.poll()
        time.sleep(0.05)
    if _process_group_exists(pgid):
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    try:
        proc.wait(timeout=max(1.0, grace_s))
    except subprocess.TimeoutExpired:
        try:
            os.killpg(pgid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


@contextlib.contextmanager
def _complex_lock(lock_path: Path):
    """Yield whether an advisory interprocess lock was acquired, without waiting."""
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            if exc.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            yield False
            return
        lock.seek(0)
        lock.truncate()
        lock.write(json.dumps({"pid": os.getpid(), "acquired_unix_s": time.time()}) + "\n")
        lock.flush()
        try:
            yield True
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _open_secure_lock_file(lock_path: Path):
    """Open a user-owned regular lock file without following a final symlink."""
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(lock_path, flags, 0o600)
    try:
        metadata = os.fstat(fd)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"GPU lock is not a regular file: {lock_path}")
        if hasattr(os, "getuid") and metadata.st_uid != os.getuid():
            raise PermissionError(f"GPU lock is not owned by the current user: {lock_path}")
        os.fchmod(fd, 0o600)
        return os.fdopen(fd, "r+", encoding="utf-8")
    except BaseException:
        os.close(fd)
        raise


@contextlib.contextmanager
def _gpu_lock(
    lock_path: Path, gpu_device_id: int | None = None, timeout_s: float = 3600.0
):
    """Serialize use of one GPU across independent docking pipelines.

    Another valid GPU job is normal contention, but the wait is bounded so a
    wedged live holder cannot stall a batch forever. POSIX releases the advisory
    lock if this process dies.
    """
    lock_path = Path(lock_path)
    timeout_s = float(timeout_s)
    deadline = time.monotonic() + timeout_s
    with _open_secure_lock_file(lock_path) as lock:
        while True:
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except OSError as exc:
                if exc.errno not in (errno.EACCES, errno.EAGAIN):
                    raise
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise TimeoutError(
                        f"timed out after {timeout_s:g}s waiting for GPU lock {lock_path}"
                    ) from exc
                time.sleep(min(0.1, remaining))
        try:
            lock.seek(0)
            lock.truncate()
            lock.write(json.dumps({
                "pid": os.getpid(), "engine": "unidock2",
                "gpu_device_id": gpu_device_id, "acquired_unix_s": time.time(),
            }) + "\n")
            lock.flush()
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def plan_complex(cdir: Path, cfg: dict):
    """Resolve inputs for one complex. Returns (info, ok) where ok is False if the
    complex should be skipped (info['skip'] says why)."""
    pdb_id = cdir.name
    conv = cfg["converter"]
    stg = cdir / "_staging"
    # Require the exact raw receptor and start conformer staged by the Vina cell.
    # Silently falling back to a converter-produced receptor changes the experiment.
    rec = stg / "receptors" / f"{pdb_id}_protein.pdb"
    lig = stg / "ligands" / f"{pdb_id}_ligand_start_conf.sdf"
    box_file = stg / "receptors" / "pdbqt" / f"{pdb_id}_protein_{conv}.box.txt"
    info = {
        "pdb_id": pdb_id, "receptor": str(rec), "ligand": str(lig),
        "box_file": str(box_file),
    }

    missing = [label for label, path in (
        ("raw receptor PDB", rec), ("start-conformer SDF", lig), ("box file", box_file)
    ) if not path.is_file()]
    if missing:
        info["plan_error"] = (
            f"missing {', '.join(missing)} (run the AutoDock Vina staging cell first)"
        )
        return info, False

    validate_single_input_ligand(lig)
    n_atoms = receptor_heavy_atoms(rec)
    if n_atoms < 1:
        raise ValueError(f"{rec}: raw protein has no heavy ATOM records")
    info["input_protein_heavy_atoms"] = n_atoms
    mra = cfg.get("max_receptor_atoms")
    if mra is not None and n_atoms > int(mra):
        info["skip"] = (
            f"input protein {n_atoms} heavy ATOM records > max_receptor_atoms {mra}"
        )
        return info, False

    (cx, cy, cz), (sx, sy, sz) = read_box(box_file)
    info.update(center=[cx, cy, cz], size=[sx, sy, sz])
    return info, True


def effective_unidock2_config(info: dict, cfg: dict, temp_dir_name: Path | None = None):
    """Exact nested YAML consumed by Uni-Dock2 (also fingerprinted)."""
    # Per-complex Uni-Dock2 config. NOTE: Uni-Dock2's YAML parser (io/yaml.py) only
    # reads keys NESTED under the group headers Required/Advanced/Hardware/Settings/
    # Preprocessing — flat keys are silently ignored (fall back to defaults). The box
    # extent is `box_size` under Settings (NOT `size`); there is no CLI flag for it, so
    # it MUST come from here or docking silently runs in the default 30 Å box.
    # receptor/ligand/center/output are supplied on the command line (they override YAML).
    advanced = {
        "exhaustiveness": int(cfg["exhaustiveness"]),
        "mc_steps": int(cfg["mc_steps"]),
        "num_pose": int(cfg["num_pose"]),
        "energy_range": float(cfg["energy_range"]),
        "seed": int(cfg["seed"]),
        "randomize": cfg["randomize"],
        "opt_steps": int(cfg["opt_steps"]),
        "refine_steps": int(cfg["refine_steps"]),
        "use_tor_lib": cfg["use_tor_lib"],
        "energy_decomp": cfg["energy_decomp"],
    }
    if cfg.get("rmsd_limit") is not None:
        advanced["rmsd_limit"] = float(cfg["rmsd_limit"])
    preprocessing = {
        "construct_ff": cfg["construct_ff"],
        "template_docking": cfg["template_docking"],
        "compute_center": cfg["compute_center"],
        "covalent_ligand": cfg["covalent_ligand"],
        "preserve_receptor_hydrogen": cfg["preserve_receptor_hydrogen"],
        "engine_checkpoint": cfg["engine_checkpoint"],
    }
    if temp_dir_name is not None:
        preprocessing["temp_dir_name"] = str(Path(temp_dir_name).resolve())
    result = {
        "Advanced": advanced,
        "Hardware": {
            "gpu_device_id": int(cfg["gpu_device_id"]),
            "n_cpu": int(cfg["n_cpu"]),
        },
        "Settings": {
            "box_size": [float(s) for s in info["size"]],
            "task": cfg["task"],
            "search_mode": cfg["search_mode"],
        },
        "Preprocessing": preprocessing,
    }
    return result


def _file_provenance(path: Path):
    path = Path(path).resolve()
    return {"path": str(path), "size": path.stat().st_size, "sha256": _sha256(path)}


def build_provenance(
    info: dict, effective_cfg: dict, tool_identity: dict, gpu_lock_file: str
):
    effective_json = json.dumps(
        effective_cfg, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    payload = {
        "schema_version": _SENTINEL_SCHEMA,
        "engine": "unidock2",
        "tool": tool_identity,
        "inputs": {
            "receptor": _file_provenance(Path(info["receptor"])),
            "ligand": _file_provenance(Path(info["ligand"])),
            "box": _file_provenance(Path(info["box_file"])),
        },
        "center": [float(value) for value in info["center"]],
        "size": [float(value) for value in info["size"]],
        "effective_config": effective_cfg,
        "effective_config_sha256": hashlib.sha256(effective_json.encode("utf-8")).hexdigest(),
        "driver_execution": {
            "gpu_lock_file": str(Path(gpu_lock_file).expanduser().resolve()),
        },
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return payload, hashlib.sha256(canonical.encode("utf-8")).hexdigest()


_STABLE_RUNTIME_KEYS = (
    "search_mode", "task", "energy_range", "n_cpu", "gpu_device_id",
    "seed", "exhaustiveness", "mc_steps", "num_pose", "rmsd_limit",
    "randomize", "opt_steps", "refine_steps", "use_tor_lib", "energy_decomp",
    "construct_ff", "template_docking", "compute_center", "covalent_ligand",
    "preserve_receptor_hydrogen", "engine_checkpoint",
    "center", "box_size",
)


def _requested_runtime_settings(expected_provenance: dict):
    protocol = expected_provenance["effective_config"]
    advanced = protocol["Advanced"]
    hardware = protocol["Hardware"]
    settings = protocol["Settings"]
    preprocessing = protocol["Preprocessing"]
    return {
        "search_mode": settings["search_mode"],
        "task": settings["task"],
        "energy_range": float(advanced["energy_range"]),
        "n_cpu": int(hardware["n_cpu"]),
        "gpu_device_id": int(hardware["gpu_device_id"]),
        "seed": int(advanced["seed"]),
        "exhaustiveness": int(advanced["exhaustiveness"]),
        "mc_steps": int(advanced["mc_steps"]),
        "num_pose": int(advanced["num_pose"]),
        # Uni-Dock2 0.6.3 defaults to 1.0 when the optional YAML key is absent.
        "rmsd_limit": float(advanced.get("rmsd_limit", 1.0)),
        "randomize": advanced["randomize"],
        "opt_steps": int(advanced["opt_steps"]),
        "refine_steps": int(advanced["refine_steps"]),
        "use_tor_lib": advanced["use_tor_lib"],
        "energy_decomp": advanced["energy_decomp"],
        "construct_ff": preprocessing["construct_ff"],
        "template_docking": preprocessing["template_docking"],
        "compute_center": preprocessing["compute_center"],
        "covalent_ligand": preprocessing["covalent_ligand"],
        "preserve_receptor_hydrogen": preprocessing["preserve_receptor_hydrogen"],
        "engine_checkpoint": preprocessing["engine_checkpoint"],
        "center": [float(value) for value in expected_provenance["center"]],
        "box_size": [float(value) for value in settings["box_size"]],
    }


def _json_equal(first, second) -> bool:
    try:
        return json.dumps(first, sort_keys=True, separators=(",", ":"), allow_nan=False) == json.dumps(
            second, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
    except (TypeError, ValueError):
        return False


def _runtime_attestation_matches(runtime: dict, expected: dict) -> bool:
    if not isinstance(runtime, dict) or set(runtime) != set(_STABLE_RUNTIME_KEYS):
        return False
    for key in _STABLE_RUNTIME_KEYS:
        if key in {"center", "box_size"}:
            observed, wanted = runtime[key], expected[key]
            if not isinstance(observed, list) or len(observed) != 3:
                return False
            try:
                if any(
                    not math.isfinite(float(value))
                    or not math.isclose(float(value), float(target), rel_tol=0.0, abs_tol=0.005)
                    for value, target in zip(observed, wanted)
                ):
                    return False
            except (TypeError, ValueError):
                return False
        elif not _json_equal(runtime[key], expected[key]):
            return False
    return True


def _completion_matches(
    out_sdf: Path, done_marker: Path, fingerprint: str,
    expected_provenance: dict, cfg: dict,
) -> bool:
    done = _read_json(done_marker)
    out_sdf = Path(out_sdf)
    suffix = "_unidock2_out.sdf"
    if not out_sdf.name.endswith(suffix):
        return False
    pdb_id = out_sdf.name[:-len(suffix)]
    expected_prepared_name = f"{pdb_id}_unidock2_receptor_prepared.pdb"
    expected_requested = _requested_runtime_settings(expected_provenance)
    if not (
        done
        and done.get("schema_version") == _SENTINEL_SCHEMA
        and done.get("engine") == "unidock2"
        and done.get("status") == "success"
        and done.get("fingerprint") == fingerprint
        and isinstance(done.get("generation_id"), str)
        and re.fullmatch(r"[0-9a-f]{32}", done["generation_id"])
        and done.get("output_file") == out_sdf.name
        and Path(done_marker).name == f"{pdb_id}{_DONE_SUFFIX}"
        and all(_json_equal(done.get(key), value) for key, value in expected_provenance.items())
    ):
        return False
    try:
        for expected_input in expected_provenance["inputs"].values():
            if not _json_equal(
                _file_provenance(Path(expected_input["path"])), expected_input
            ):
                return False
        output_validation = _validate_docking_output(out_sdf, cfg)
        if (
            done.get("num_poses") != output_validation["num_poses"]
            or not isinstance(done.get("best_affinity"), (int, float))
            or isinstance(done.get("best_affinity"), bool)
            or not math.isclose(
                float(done["best_affinity"]), output_validation["best_affinity"],
                rel_tol=0.0, abs_tol=1e-12,
            )
        ):
            return False
        if done.get("output_sha256") != _sha256(out_sdf):
            return False
        prepared_name = done.get("prepared_receptor_file")
        if prepared_name != expected_prepared_name:
            return False
        prepared = out_sdf.parent / prepared_name
        prepared_validation = _validate_prepared_receptor(prepared)
        runtime = done.get("runtime_attestation")
        return (
            done.get("prepared_receptor_sha256") == _sha256(prepared)
            and done.get("prepared_receptor_atoms")
            == prepared_validation["prepared_receptor_atoms"]
            and done.get("prepared_receptor_heavy_atoms")
            == prepared_validation["prepared_receptor_heavy_atoms"]
            and done.get("engine_receptor_heavy_atoms_in_box")
            == prepared_validation["prepared_receptor_heavy_atoms"]
            and done.get("input_protein_heavy_atoms")
            == receptor_heavy_atoms(Path(expected_provenance["inputs"]["receptor"]["path"]))
            and _json_equal(done.get("requested_search_settings"), expected_requested)
            and _runtime_attestation_matches(runtime, expected_requested)
            and _json_equal(done.get("effective_search_settings"), runtime)
            and _json_equal(done.get("applied_center"), runtime.get("center"))
            and _json_equal(done.get("applied_box_size"), runtime.get("box_size"))
        )
    except (KeyError, OSError, TypeError, UnicodeError, ValueError):
        return False


_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def _parse_runtime_attestation(
    attempt_log: str, info: dict, cfg: dict, runtime_cfg: dict | None = None
):
    """Parse and verify the settings the GPU engine says it actually used."""
    parameter_maps = []
    for line in attempt_log.splitlines():
        stripped = line.strip()
        if not (stripped.startswith("{") and stripped.endswith("}")):
            continue
        try:
            parsed = ast.literal_eval(stripped)
        except (SyntaxError, ValueError):
            continue
        if isinstance(parsed, dict) and "box_size" in parsed and "search_mode" in parsed:
            parameter_maps.append(parsed)
    if len(parameter_maps) != 1:
        raise ValueError(
            f"expected one parsed Uni-Dock2 parameter mapping; found {len(parameter_maps)}"
        )
    python_params = parameter_maps[0]

    blocks = re.findall(
        r"DockParam:\s*(.*?)(?=\n(?:\[[^\n]+\]\s+\[info\]\s+)?Group\s|\Z)",
        attempt_log,
        flags=re.DOTALL,
    )
    if len(blocks) != 1:
        raise ValueError(f"expected one DockParam block in attempt log; found {len(blocks)}")
    block = blocks[0]

    def integer(name):
        match = re.search(rf"\b{re.escape(name)}=(-?\d+)\b", block)
        if not match:
            raise ValueError(f"runtime attestation is missing {name}")
        return int(match.group(1))

    actual = {
        "seed": integer("seed"),
        "exhaustiveness": integer("exhaustiveness"),
        "mc_steps": integer("mc_steps"),
        "num_pose": integer("num_pose"),
        "opt_steps": integer("opt_steps"),
        "refine_steps": integer("refine_steps"),
    }
    expected = {
        "seed": int(cfg["seed"]),
        "exhaustiveness": int(cfg["exhaustiveness"]),
        "mc_steps": int(cfg["mc_steps"]),
        "num_pose": int(cfg["num_pose"]),
        "opt_steps": int(cfg["opt_steps"]),
        "refine_steps": int(cfg["refine_steps"]),
    }
    mismatches = [f"{key}: requested {expected[key]}, engine used {actual[key]}"
                  for key in expected if actual[key] != expected[key]]
    rmsd_match = re.search(rf"\brmsd_limit=({_FLOAT_PATTERN})\s+Angstrom", block)
    if not rmsd_match:
        raise ValueError("runtime attestation is missing rmsd_limit")
    actual["rmsd_limit"] = float(rmsd_match.group(1))
    requested_rmsd = float(cfg["rmsd_limit"] if cfg.get("rmsd_limit") is not None else 1.0)
    if not math.isclose(actual["rmsd_limit"], requested_rmsd, rel_tol=0.0, abs_tol=1e-6):
        mismatches.append(
            f"rmsd_limit: requested {requested_rmsd}, engine used {actual['rmsd_limit']}"
        )

    expected_python = {
        "search_mode": cfg["search_mode"], "task": cfg["task"],
        "energy_range": float(cfg["energy_range"]),
        "box_size": [float(value) for value in info["size"]],
        "n_cpu": int(cfg["n_cpu"]), "gpu_device_id": int(cfg["gpu_device_id"]),
        "rmsd_limit": requested_rmsd,
        "randomize": cfg["randomize"],
        "opt_steps": int(cfg["opt_steps"]),
        "refine_steps": int(cfg["refine_steps"]),
        "use_tor_lib": cfg["use_tor_lib"],
        "energy_decomp": cfg["energy_decomp"],
        "construct_ff": cfg["construct_ff"],
        "template_docking": cfg["template_docking"],
        "compute_center": cfg["compute_center"],
        "covalent_ligand": cfg["covalent_ligand"],
        "preserve_receptor_hydrogen": cfg["preserve_receptor_hydrogen"],
        "engine_checkpoint": cfg["engine_checkpoint"],
    }
    if runtime_cfg is not None:
        expected_python["temp_dir_name"] = runtime_cfg["Preprocessing"]["temp_dir_name"]
    for key, wanted in expected_python.items():
        if key not in python_params:
            mismatches.append(f"Python parameter mapping omitted {key}")
            continue
        observed = python_params[key]
        if isinstance(wanted, list):
            if not isinstance(observed, (list, tuple)) or len(observed) != len(wanted) or any(
                not math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1e-6)
                for a, b in zip(wanted, observed)
            ):
                mismatches.append(f"{key}: requested {wanted}, parser used {observed}")
        elif isinstance(wanted, float):
            try:
                same = math.isclose(wanted, float(observed), rel_tol=0.0, abs_tol=1e-9)
            except (TypeError, ValueError):
                same = False
            if not same:
                mismatches.append(f"{key}: requested {wanted}, parser used {observed}")
        elif observed != wanted:
            mismatches.append(f"{key}: requested {wanted}, parser used {observed}")

    bounds = {}
    for axis in "xyz":
        for side in ("lo", "hi"):
            key = f"{axis}_{side}"
            match = re.search(rf"\b{key}=({_FLOAT_PATTERN})\s+Angstrom", block)
            if not match:
                raise ValueError(f"runtime attestation is missing box bound {key}")
            bounds[key] = float(match.group(1))
    actual_center = [(bounds[f"{a}_lo"] + bounds[f"{a}_hi"]) / 2.0 for a in "xyz"]
    actual_size = [bounds[f"{a}_hi"] - bounds[f"{a}_lo"] for a in "xyz"]
    for label, wanted, observed in (
        ("center", info["center"], actual_center), ("box_size", info["size"], actual_size)
    ):
        for axis, (requested, used) in zip("xyz", zip(wanted, observed)):
            if not math.isclose(float(requested), used, rel_tol=0.0, abs_tol=0.005):
                mismatches.append(f"{label}.{axis}: requested {requested}, engine used {used}")

    atom_matches = re.findall(r"Receptor has\s+(\d+)\s+heavy atoms in box", attempt_log)
    if len(atom_matches) != 1:
        raise ValueError(
            "runtime attestation must contain exactly one engine receptor-heavy-atom count"
        )
    if int(atom_matches[0]) < 1:
        raise ValueError("engine reported no receptor heavy atoms in the docking box")
    if mismatches:
        raise ValueError("runtime settings mismatch: " + "; ".join(mismatches))
    actual.update({
        "search_mode": python_params["search_mode"],
        "task": python_params["task"],
        "energy_range": float(python_params["energy_range"]),
        "n_cpu": int(python_params["n_cpu"]),
        "gpu_device_id": int(python_params["gpu_device_id"]),
        "randomize": python_params["randomize"],
        "use_tor_lib": python_params["use_tor_lib"],
        "energy_decomp": python_params["energy_decomp"],
        "construct_ff": python_params["construct_ff"],
        "template_docking": python_params["template_docking"],
        "compute_center": python_params["compute_center"],
        "covalent_ligand": python_params["covalent_ligand"],
        "preserve_receptor_hydrogen": python_params["preserve_receptor_hydrogen"],
        "engine_checkpoint": python_params["engine_checkpoint"],
    })
    return {
        **actual,
        "center": actual_center,
        "box_size": actual_size,
        "engine_receptor_heavy_atoms_in_box": int(atom_matches[0]),
        "python_parameters": {key: python_params.get(key) for key in expected_python},
    }


def _stable_runtime_attestation(attestation: dict):
    """Strip attempt-specific parser details, especially the temporary path."""
    return {key: attestation[key] for key in _STABLE_RUNTIME_KEYS}


def _read_log_segment(log_path: Path, start_offset: int) -> str:
    with open(log_path, "rb") as fh:
        fh.seek(start_offset)
        return fh.read().decode("utf-8", errors="replace")


def _quarantine_artifact(path: Path, reason: str) -> Path | None:
    """Move one explicit artifact aside; never delete it or overwrite a prior copy."""
    path = Path(path)
    if not path.exists() and not path.is_symlink():
        return None
    safe_reason = re.sub(r"[^A-Za-z0-9_.-]+", "-", reason).strip("-")[:40] or "stale"
    quarantine = path.parent / "quarantine"
    quarantine.mkdir(parents=True, exist_ok=True)
    target = quarantine / f"{path.name}.{time.time_ns()}.{safe_reason}"
    os.replace(path, target)
    return target


def _backup_artifact(path: Path, reason: str) -> Path | None:
    """Copy a public artifact aside before replacement without disturbing it."""
    path = Path(path)
    if not path.exists() and not path.is_symlink():
        return None
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"public artifact is not a regular file: {path}")
    safe_reason = re.sub(r"[^A-Za-z0-9_.-]+", "-", reason).strip("-")[:40] or "backup"
    quarantine = path.parent / "quarantine"
    quarantine.mkdir(parents=True, exist_ok=True)
    target = quarantine / f"{path.name}.{time.time_ns()}.{safe_reason}"
    with open(path, "rb") as source, open(target, "xb") as destination:
        shutil.copyfileobj(source, destination)
        destination.flush()
        os.fsync(destination.fileno())
    shutil.copystat(path, target, follow_symlinks=False)
    return target


def _restore_backup(backup: Path, destination: Path):
    """Atomically restore a backup while retaining the audit copy."""
    backup = Path(backup)
    destination = Path(destination)
    restore = destination.with_name(f".{destination.name}.{uuid.uuid4().hex}.restore")
    try:
        with open(backup, "rb") as source, open(restore, "xb") as output:
            shutil.copyfileobj(source, output)
            output.flush()
            os.fsync(output.fileno())
        shutil.copystat(backup, restore, follow_symlinks=False)
        os.replace(restore, destination)
    finally:
        try:
            restore.unlink()
        except FileNotFoundError:
            pass


def _promote_generation(promotions, verify):
    """Publish one multi-file generation, rolling back every target on failure.

    POSIX rename is atomic per file rather than across files. Keeping independent
    audit backups lets us restore the complete prior generation if any rename,
    manifest write, or post-publish verification fails.
    """
    backups = {}
    promoted = []
    try:
        for _candidate, destination in promotions:
            backups[destination] = _backup_artifact(destination, "pre-promotion-backup")
        for candidate, destination in promotions:
            os.replace(candidate, destination)
            promoted.append(destination)
        if not verify():
            raise RuntimeError("published generation failed completion verification")
        return backups
    except BaseException as original:
        rollback_errors = []
        for _candidate, destination in reversed(promotions):
            backup = backups.get(destination)
            try:
                if backup is not None:
                    _restore_backup(backup, destination)
                elif destination in promoted and (destination.exists() or destination.is_symlink()):
                    _quarantine_artifact(destination, "rolled-back-new-generation")
            except Exception as exc:
                rollback_errors.append(f"{destination}: {exc}")
        if rollback_errors:
            raise RuntimeError(
                f"generation promotion failed ({original}); rollback also failed: "
                + "; ".join(rollback_errors)
            ) from original
        raise


def _emit_cleanup_warning(message: str):
    """Report cleanup trouble without warning filters changing run semantics."""
    try:
        warnings.warn(message, RuntimeWarning, stacklevel=3)
    except Exception:
        pass


def _cleanup_scratch(scratch_root: Path | None, scratch_parent: Path | None):
    messages = []
    if scratch_root is not None:
        last_error = None
        for _ in range(3):
            try:
                shutil.rmtree(scratch_root)
                last_error = None
                break
            except FileNotFoundError:
                last_error = None
                break
            except OSError as exc:
                last_error = exc
                time.sleep(0.1)
        if last_error is not None:
            messages.append(f"could not remove driver scratch {scratch_root}: {last_error}")
    if scratch_parent is not None:
        try:
            scratch_parent.rmdir()
        except OSError:
            pass
    if messages:
        message = "; ".join(messages)
        _emit_cleanup_warning(message)
        return message
    return None


def _ensure_private_directory(path: Path):
    """Create/verify a non-symlink, current-user-only scratch directory."""
    path = Path(path)
    path.mkdir(parents=True, mode=0o700, exist_ok=True)
    metadata = path.lstat()
    if not stat.S_ISDIR(metadata.st_mode) or path.is_symlink():
        raise ValueError(f"scratch path is not a real directory: {path}")
    if hasattr(os, "getuid") and metadata.st_uid != os.getuid():
        raise PermissionError(f"scratch directory is not owned by the current user: {path}")
    path.chmod(0o700)
    return path


def _validate_prepared_receptor(pdb: Path):
    atom_count = heavy_count = 0
    for line_number, line in enumerate(Path(pdb).read_text(errors="strict").splitlines(), 1):
        if not line.startswith("ATOM"):
            continue
        atom_count += 1
        try:
            xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        except ValueError as exc:
            raise ValueError(f"{pdb}:{line_number} has invalid PDB coordinates") from exc
        if not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"{pdb}:{line_number} has non-finite PDB coordinates")
        atom_name = re.sub(r"^[0-9]+", "", line[12:16].strip()).upper()
        element = (line[76:78].strip() or atom_name[:1]).upper()
        if element not in {"H", "D", "T"}:
            heavy_count += 1
    if atom_count < 1 or heavy_count < 1:
        raise ValueError(f"{pdb}: prepared receptor has no valid heavy ATOM records")
    return {"prepared_receptor_atoms": atom_count, "prepared_receptor_heavy_atoms": heavy_count}


def _materialize_prepared_receptor(
    scratch_root: Path, tool_identity: dict, destination: Path, env: dict
):
    topologies = sorted(Path(scratch_root).rglob("receptor.prmtop"))
    candidates = [(top, top.with_name("receptor.inpcrd")) for top in topologies
                  if top.with_name("receptor.inpcrd").is_file()]
    if len(candidates) != 1:
        raise ValueError(
            f"expected one receptor.prmtop/inpcrd pair in scratch; found {len(candidates)}"
        )
    prmtop, inpcrd = candidates[0]
    stderr = b""
    with open(destination, "xb") as output:
        proc = subprocess.Popen(
            [tool_identity["ambpdb_path"], "-p", str(prmtop), "-c", str(inpcrd)],
            stdout=output, stderr=subprocess.PIPE, start_new_session=True, env=env,
        )
        try:
            _, stderr = proc.communicate(timeout=120)
        except subprocess.TimeoutExpired:
            _terminate_process_group(proc, grace_s=5.0)
            raise RuntimeError("ambpdb conversion timed out after 120 seconds")
        except BaseException:
            _terminate_process_group(proc, grace_s=5.0)
            raise
        output.flush()
        os.fsync(output.fileno())
    if proc.returncode != 0:
        raise RuntimeError(
            f"ambpdb failed with rc={proc.returncode}: "
            f"{stderr.decode('utf-8', errors='replace').strip()}"
        )
    validation = _validate_prepared_receptor(destination)
    return {
        **validation,
        "prepared_receptor_sha256": _sha256(destination),
        "receptor_prmtop_sha256": _sha256(prmtop),
        "receptor_inpcrd_sha256": _sha256(inpcrd),
    }


def _rmtree_quiet(path: Path):
    """Remove a driver-owned scratch tree, warning (not raising) on trouble."""
    try:
        shutil.rmtree(path)
    except FileNotFoundError:
        pass
    except OSError as exc:
        _emit_cleanup_warning(f"could not remove prep scratch {path}: {exc}")


def _cleanup_prep_bundle(bundle):
    """Drop the transient prep scratch once its docking stage has consumed it."""
    scratch = bundle.get("prep_scratch") if isinstance(bundle, dict) else None
    if scratch:
        _rmtree_quiet(Path(scratch))


def _prep_config(cfg: dict, temp_dir: Path) -> dict:
    """Nested YAML for ``unidock2 protein_prep`` (reads only the Preprocessing group).

    ``preserve_receptor_hydrogen`` MUST match the docking config so the parameterized
    receptor is identical to the one the all-in-one ``docking`` path would have built.
    """
    return {
        "Preprocessing": {
            "preserve_receptor_hydrogen": cfg["preserve_receptor_hydrogen"],
            "temp_dir_name": str(Path(temp_dir).resolve()),
        }
    }


def prepare_complex(cdir: Path, info: dict, cfg: dict, unidock2_bin: str,
                    tool_identity: dict, prep_parent: Path):
    """CPU receptor preparation for one complex — no GPU.

    Runs ``unidock2 protein_prep`` (pdbfixer + AmberTools tleap) on the raw protein
    PDB, producing a parameterized ``.dms`` receptor plus the exported prepared-receptor
    PDB (materialized from tleap's prmtop/inpcrd, exactly as the all-in-one path exports
    it). Executes on a worker thread; returns a bundle the serial GPU docking stage
    consumes. Never touches the output tree, so a failure or interrupt leaves nothing to
    resume — the complex is simply re-prepared next run.
    """
    t0 = time.time()
    pdb_id = Path(cdir).name
    rec = Path(info["receptor"]).resolve()
    prep_log = Path(cfg["log_dir"]) / f"{pdb_id}.prep.log"
    prep_scratch = Path(tempfile.mkdtemp(prefix=f"prep-{pdb_id}-", dir=prep_parent))

    def _fail(status, error):
        _rmtree_quiet(prep_scratch)
        return {"pdb_id": pdb_id, "status": status, "error": error,
                "elapsed_s": round(time.time() - t0, 1)}

    try:
        prep_tmp = _ensure_private_directory(prep_scratch / "tmp")
        dms = prep_scratch / f"{pdb_id}.dms"
        prepared_pdb = prep_scratch / f"{pdb_id}_prepared.pdb"
        prep_cfg_path = prep_scratch / f"{pdb_id}_prep_config.yaml"
        _atomic_write_text(
            prep_cfg_path, yaml.safe_dump(_prep_config(cfg, prep_tmp), sort_keys=False)
        )
        cmd = [
            unidock2_bin, "protein_prep", "-r", str(rec),
            "-o", str(dms.resolve()), "-cf", str(prep_cfg_path.resolve()),
        ]
        rc, timed_out = _run(
            cmd, int(cfg["timeout_s"]), prep_log, _subprocess_env(unidock2_bin),
            terminate_grace_s=float(cfg.get("terminate_grace_s", 5.0)),
            on_spawn=_register_prep, on_reap=_deregister_prep,
        )
        if timed_out:
            return _fail("timeout", f"protein_prep exceeded timeout_s={cfg['timeout_s']}")
        if rc != 0:
            return _fail("error", f"protein_prep exited with return code {rc}")
        if not dms.is_file() or dms.stat().st_size == 0:
            return _fail("error", "protein_prep produced no .dms receptor")
        # Export the prepared receptor here (docking from a .dms no longer emits
        # prmtop/inpcrd); the box-coverage parity check happens in dock_complex.
        prepared_validation = _materialize_prepared_receptor(
            prep_tmp, tool_identity, prepared_pdb, _subprocess_env(unidock2_bin)
        )
        return {
            "pdb_id": pdb_id, "status": "prepared",
            "dms": str(dms), "prepared_pdb": str(prepared_pdb),
            "prepared_validation": prepared_validation,
            "prep_scratch": str(prep_scratch),
            "elapsed_s": round(time.time() - t0, 1),
        }
    except (KeyboardInterrupt, SystemExit):
        _rmtree_quiet(prep_scratch)
        raise
    except BaseException as exc:  # noqa: BLE001 - report, do not crash the pool
        return _fail("error", f"receptor preparation failed: {type(exc).__name__}: {exc}")


def _install_prepared_receptor(prepared: dict, destination: Path) -> dict:
    """Copy the prep-stage prepared-receptor PDB into the docking generation and
    re-validate it, returning the same shape ``_materialize_prepared_receptor`` does.
    Fails loudly if the receptor changed between preparation and docking."""
    source = Path(prepared["prepared_pdb"])
    destination = Path(destination)
    with open(source, "rb") as src, open(destination, "xb") as dst:
        shutil.copyfileobj(src, dst)
        dst.flush()
        os.fsync(dst.fileno())
    validation = _validate_prepared_receptor(destination)
    pv = prepared["prepared_validation"]
    if (
        validation["prepared_receptor_atoms"] != pv["prepared_receptor_atoms"]
        or validation["prepared_receptor_heavy_atoms"] != pv["prepared_receptor_heavy_atoms"]
        or _sha256(destination) != pv["prepared_receptor_sha256"]
    ):
        raise RuntimeError("prepared receptor changed between preparation and docking")
    return {
        **validation,
        "prepared_receptor_sha256": pv["prepared_receptor_sha256"],
        "receptor_prmtop_sha256": pv["receptor_prmtop_sha256"],
        "receptor_inpcrd_sha256": pv["receptor_inpcrd_sha256"],
    }


def _result_from_manifest(
    info: dict, manifest: dict, out_sdf: Path, done_marker: Path,
    *, status: str = "done", elapsed_s: float = 0.0,
):
    """Return the same rich, canonical scientific result for run and resume."""
    prepared = Path(out_sdf).parent / manifest["prepared_receptor_file"]
    return {
        **info,
        "status": status,
        "num_poses": manifest["num_poses"],
        "best_affinity": manifest["best_affinity"],
        "scoring": "vina",
        "engine": "unidock2",
        "tool_version": manifest["tool"]["version"],
        "fingerprint": manifest["fingerprint"],
        "generation_id": manifest.get("generation_id"),
        "requested_search_settings": manifest["requested_search_settings"],
        "runtime_attestation": manifest["runtime_attestation"],
        "effective_search_settings": manifest["effective_search_settings"],
        "applied_center": manifest["applied_center"],
        "applied_box_size": manifest["applied_box_size"],
        "input_protein_heavy_atoms": manifest["input_protein_heavy_atoms"],
        "engine_receptor_heavy_atoms_in_box": manifest[
            "engine_receptor_heavy_atoms_in_box"
        ],
        "prepared_receptor_atoms": manifest["prepared_receptor_atoms"],
        "prepared_receptor_heavy_atoms": manifest["prepared_receptor_heavy_atoms"],
        "output_sha256": manifest["output_sha256"],
        "prepared_receptor_sha256": manifest["prepared_receptor_sha256"],
        "published_output_current": True,
        "output_sdf": str(out_sdf),
        "completion_manifest": str(done_marker),
        "prepared_receptor_pdb": str(prepared),
        "elapsed_s": round(float(elapsed_s), 1),
    }


def dock_complex(cdir: Path, cfg: dict, unidock2_bin: str, tool_identity: dict | None = None,
                 prepared: dict | None = None):
    """Dock one complex with locking, runtime attestation, and atomic publication.

    ``prepared`` (from :func:`prepare_complex`) supplies a pre-built ``.dms`` receptor
    and the already-exported prepared-receptor PDB: docking then consumes the ``.dms``
    (skipping its internal prep) and the export is installed rather than re-materialized.
    When ``prepared`` is None the receptor is the raw PDB and Uni-Dock2 preps it inline.
    """
    t0 = time.time()
    validate_config(cfg)
    identity = tool_identity or probe_tool_identity(unidock2_bin)
    pdb_id = Path(cdir).name
    out_dir = Path(cfg["output_dir"]) / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = out_dir / f"{pdb_id}_unidock2_out.sdf"
    prepared_pdb = out_dir / f"{pdb_id}_unidock2_receptor_prepared.pdb"
    done_marker = out_dir / f"{pdb_id}{_DONE_SUFFIX}"
    stable_cfg_path = out_dir / f"{pdb_id}_unidock2_config.yaml"
    log = Path(cfg["log_dir"]) / f"{pdb_id}.log"

    with _complex_lock(out_dir / ".unidock2.lock") as acquired:
        if not acquired:
            return {
                "pdb_id": pdb_id, "status": "locked",
                "error": "another process is operating on this complex",
                "num_poses": 0, "best_affinity": None,
                "elapsed_s": round(time.time() - t0, 1),
            }

        info, ok = plan_complex(Path(cdir), cfg)
        if not ok:
            status = "error" if info.get("plan_error") else "skipped"
            return {
                **info, "status": status,
                "error": info.get("plan_error"),
                "elapsed_s": round(time.time() - t0, 1),
            }

        protocol_cfg = effective_unidock2_config(info, cfg)
        provenance, fingerprint = build_provenance(
            info, protocol_cfg, identity, cfg["gpu_lock_file"]
        )
        prior_completion_valid = _completion_matches(
            out_sdf, done_marker, fingerprint, provenance, cfg
        )
        if not cfg.get("overwrite", False) and prior_completion_valid:
            done = _read_json(done_marker)
            return _result_from_manifest(
                info, done, out_sdf, done_marker,
                status="done", elapsed_s=time.time() - t0,
            )

        # A required attempt supersedes completion authority immediately. Keep
        # prior payload bytes in place until marker-last promotion, but remove
        # the canonical marker so a failed overwrite followed by overwrite=false
        # retries instead of silently accepting the pre-overwrite generation.
        _quarantine_artifact(done_marker, "invalidated-before-attempt")
        _quarantine_artifact(
            out_dir / ".unidock2_done", "invalidated-legacy-manifest"
        )

        scratch_base = (
            Path(cfg["scratch_dir"]) if cfg.get("scratch_dir")
            else Path("/tmp") / f"unidock2_driver_{os.getuid()}"
        )
        output_token = hashlib.sha256(
            str(Path(cfg["output_dir"])).encode("utf-8")
        ).hexdigest()[:16]
        scratch_parent = scratch_base / output_token
        _ensure_private_directory(scratch_base)

        # With both batch and complex locks held, same-output-tree leftovers are
        # abandoned. Quarantine diagnostic files and remove driver-owned scratch.
        for abandoned in sorted(out_dir.glob(f".{pdb_id}_unidock2_out.*.partial.sdf")):
            _quarantine_artifact(abandoned, "abandoned-partial")
        for pattern, reason in (
            (f".{pdb_id}_unidock2_config.*.partial.yaml", "abandoned-config"),
            (f".{pdb_id}_unidock2_receptor_prepared.*.partial.pdb", "abandoned-receptor"),
            (f".{pdb_id}_unidock2_completion.*.partial.json", "abandoned-manifest"),
        ):
            for abandoned in sorted(out_dir.glob(pattern)):
                _quarantine_artifact(abandoned, reason)
        for parent in (scratch_parent, out_dir / ".scratch"):
            if parent.is_dir():
                for abandoned in sorted(parent.glob(f"{pdb_id}-*")):
                    if abandoned.is_dir() and abandoned.parent.resolve() == parent.resolve():
                        shutil.rmtree(abandoned)

        generation_id = uuid.uuid4().hex
        tmp_sdf = out_dir / f".{pdb_id}_unidock2_out.{generation_id}.partial.sdf"
        tmp_prepared = out_dir / f".{pdb_id}_unidock2_receptor_prepared.{generation_id}.partial.pdb"
        attempt_cfg = out_dir / f".{pdb_id}_unidock2_config.{generation_id}.partial.yaml"
        attempt_manifest = out_dir / f".{pdb_id}_unidock2_completion.{generation_id}.partial.json"
        _ensure_private_directory(scratch_parent)
        scratch_root = Path(tempfile.mkdtemp(prefix=f"{pdb_id}-{generation_id}-", dir=scratch_parent))
        runtime_cfg = effective_unidock2_config(info, cfg, temp_dir_name=scratch_root)
        try:
            _atomic_write_text(attempt_cfg, yaml.safe_dump(runtime_cfg, sort_keys=False))
        except BaseException:
            _cleanup_scratch(scratch_root, scratch_parent)
            raise

        lig = Path(info["ligand"]).resolve()
        cx, cy, cz = info["center"]
        # Fingerprint/provenance always key on the RAW receptor (info["receptor"]); the
        # .dms is a deterministic derivative, so passing it here does not shift identity.
        if prepared is not None:
            if prepared.get("pdb_id") != pdb_id or prepared.get("status") != "prepared":
                raise ValueError("prepared receptor bundle does not match this complex")
            receptor_arg = Path(prepared["dms"]).resolve()
        else:
            receptor_arg = Path(info["receptor"]).resolve()
        cmd = [
            unidock2_bin, "docking", "-r", str(receptor_arg), "-l", str(lig),
            "-c", str(cx), str(cy), str(cz), "-o", str(tmp_sdf.resolve()),
            "-cf", str(attempt_cfg.resolve()),
        ]
        log.parent.mkdir(parents=True, exist_ok=True)
        log_start = log.stat().st_size if log.exists() else 0
        rc, timed_out = None, False
        validation = None
        attestation = None
        prepared_validation = None
        committed_manifest = None
        failure = None
        gpu_lock_timed_out = False
        inputs_changed_during_run = False
        cleanup_warning = None
        try:
            try:
                with _gpu_lock(
                    Path(cfg["gpu_lock_file"]), gpu_device_id=int(cfg["gpu_device_id"]),
                    timeout_s=float(cfg["gpu_lock_timeout_s"]),
                ):
                    rc, timed_out = _run(
                        cmd, int(cfg["timeout_s"]), log, _subprocess_env(unidock2_bin),
                        terminate_grace_s=float(cfg.get("terminate_grace_s", 5.0)),
                    )
            except TimeoutError as exc:
                gpu_lock_timed_out = True
                failure = str(exc)
            except BaseException:
                _quarantine_artifact(tmp_sdf, "interrupted-partial")
                _quarantine_artifact(tmp_prepared, "interrupted-prepared-receptor")
                _quarantine_artifact(attempt_cfg, "interrupted-config")
                _quarantine_artifact(attempt_manifest, "interrupted-manifest")
                raise

            if rc == 0 and not timed_out and not gpu_lock_timed_out:
                try:
                    validation = _validate_docking_output(tmp_sdf, cfg)
                    attempt_log = _read_log_segment(log, log_start)
                    attestation = _parse_runtime_attestation(
                        attempt_log, info, cfg, runtime_cfg=runtime_cfg
                    )
                    if prepared is not None:
                        # docking from a .dms does not regenerate prmtop/inpcrd, so the
                        # export comes from the prep stage; install and re-verify it.
                        prepared_validation = _install_prepared_receptor(
                            prepared, tmp_prepared
                        )
                    else:
                        prepared_validation = _materialize_prepared_receptor(
                            scratch_root, identity, tmp_prepared, _subprocess_env(unidock2_bin)
                        )
                    if (
                        prepared_validation["prepared_receptor_heavy_atoms"]
                        != attestation["engine_receptor_heavy_atoms_in_box"]
                    ):
                        raise ValueError(
                            "prepared receptor heavy-atom count does not equal the engine "
                            "in-box count; the configured box does not cover the scored receptor"
                        )
                    _, final_fingerprint = build_provenance(
                        info, protocol_cfg, identity, cfg["gpu_lock_file"]
                    )
                    inputs_changed_during_run = final_fingerprint != fingerprint
                    if inputs_changed_during_run:
                        raise ValueError("inputs or box changed while docking was running")
                except (OSError, ValueError, RuntimeError) as exc:
                    failure = str(exc)

            if rc == 0 and not timed_out and validation and attestation and not failure:
                try:
                    candidate_prepared = _validate_prepared_receptor(tmp_prepared)
                    requested_search = _requested_runtime_settings(provenance)
                    runtime_attestation = _stable_runtime_attestation(attestation)
                    manifest = {
                        **provenance,
                        "status": "success", "fingerprint": fingerprint,
                        "generation_id": generation_id, "output_file": out_sdf.name,
                        "output_sha256": _sha256(tmp_sdf),
                        "prepared_receptor_file": prepared_pdb.name,
                        "prepared_receptor_sha256": _sha256(tmp_prepared),
                        **candidate_prepared,
                        "receptor_prmtop_sha256": prepared_validation[
                            "receptor_prmtop_sha256"
                        ],
                        "receptor_inpcrd_sha256": prepared_validation[
                            "receptor_inpcrd_sha256"
                        ],
                        "num_poses": validation["num_poses"],
                        "best_affinity": validation["best_affinity"],
                        "requested_search_settings": requested_search,
                        "runtime_attestation": runtime_attestation,
                        "effective_search_settings": runtime_attestation,
                        "applied_center": runtime_attestation["center"],
                        "applied_box_size": runtime_attestation["box_size"],
                        "input_protein_heavy_atoms": info["input_protein_heavy_atoms"],
                        "engine_receptor_heavy_atoms_in_box": attestation[
                            "engine_receptor_heavy_atoms_in_box"
                        ],
                        "completed_unix_s": time.time(),
                    }
                    _atomic_write_json(attempt_manifest, manifest)
                    promotions = [
                        (tmp_sdf, out_sdf),
                        (tmp_prepared, prepared_pdb),
                        (attempt_cfg, stable_cfg_path),
                        (attempt_manifest, done_marker),
                    ]
                    _promote_generation(
                        promotions,
                        lambda: _completion_matches(
                            out_sdf, done_marker, fingerprint, provenance, cfg
                        ),
                    )
                    committed_manifest = _read_json(done_marker)
                    validation = _validate_docking_output(out_sdf, cfg)
                    try:
                        _quarantine_artifact(
                            out_dir / ".unidock2_done", "superseded-legacy-manifest"
                        )
                    except OSError as exc:
                        _emit_cleanup_warning(
                            f"could not quarantine legacy completion marker: {exc}"
                        )
                except (OSError, ValueError, RuntimeError) as exc:
                    failure = f"generation publication failed: {exc}"
        finally:
            cleanup_warning = _cleanup_scratch(scratch_root, scratch_parent)

        success = (
            rc == 0 and not timed_out and validation is not None and attestation is not None
            and prepared_validation is not None
            and not failure and _completion_matches(
                out_sdf, done_marker, fingerprint, provenance, cfg
            )
        )
        if not success:
            partial_copy = _quarantine_artifact(
                tmp_sdf, "timeout-partial" if timed_out else "failed-partial"
            )
            config_copy = _quarantine_artifact(attempt_cfg, "failed-config")
            prepared_copy = _quarantine_artifact(tmp_prepared, "failed-prepared-receptor")
            manifest_copy = _quarantine_artifact(attempt_manifest, "failed-manifest")
            if gpu_lock_timed_out:
                status = "locked"
            elif timed_out:
                status = "timeout"
                failure = failure or f"exceeded timeout_s={cfg['timeout_s']}"
            else:
                status = "error"
                failure = failure or f"Uni-Dock2 exited with return code {rc}"
        else:
            status = "success"
            partial_copy = config_copy = prepared_copy = manifest_copy = None

        previous_generation_preserved = bool(
            not success and prior_completion_valid
            and _completion_matches(out_sdf, done_marker, fingerprint, provenance, cfg)
        )
        operational = {
            "returncode": rc,
            "timed_out": timed_out,
            "gpu_lock_timed_out": gpu_lock_timed_out,
            "error": failure,
            "inputs_changed_during_run": inputs_changed_during_run,
            "previous_generation_preserved": previous_generation_preserved,
            "scratch_cleanup_warning": cleanup_warning,
            "quarantined_partial": str(partial_copy) if partial_copy else None,
            "quarantined_attempt_config": str(config_copy) if config_copy else None,
            "quarantined_prepared_receptor": str(prepared_copy) if prepared_copy else None,
            "quarantined_attempt_manifest": str(manifest_copy) if manifest_copy else None,
            "log": str(log),
        }
        if success:
            result = _result_from_manifest(
                info, committed_manifest, out_sdf, done_marker,
                status="success", elapsed_s=time.time() - t0,
            )
            result.update(operational)
            return result

        return {
            **info, "status": status, "returncode": rc, "timed_out": timed_out,
            "gpu_lock_timed_out": gpu_lock_timed_out,
            "error": failure,
            "num_poses": 0,
            "best_affinity": None,
            "scoring": "vina", "engine": "unidock2",
            "tool_version": identity["version"], "fingerprint": fingerprint,
            "generation_id": generation_id,
            "requested_search_settings": _requested_runtime_settings(provenance),
            "runtime_attestation": (
                _stable_runtime_attestation(attestation) if attestation else None
            ),
            "effective_search_settings": (
                _stable_runtime_attestation(attestation)
                if attestation else None
            ),
            "engine_receptor_heavy_atoms_in_box": (
                attestation["engine_receptor_heavy_atoms_in_box"] if attestation else None
            ),
            "inputs_changed_during_run": inputs_changed_during_run,
            "previous_generation_preserved": previous_generation_preserved,
            "published_output_current": False,
            "output_sdf": None,
            "completion_manifest": None,
            "prepared_receptor_pdb": None,
            "quarantined_partial": str(partial_copy) if partial_copy else None,
            "quarantined_attempt_config": str(config_copy) if config_copy else None,
            "quarantined_prepared_receptor": str(prepared_copy) if prepared_copy else None,
            "quarantined_attempt_manifest": str(manifest_copy) if manifest_copy else None,
            "scratch_cleanup_warning": cleanup_warning,
            "log": str(log), "elapsed_s": round(time.time() - t0, 1),
        }


# ────────────────────────────────────────────────────────────────────────────
# Driver
# ────────────────────────────────────────────────────────────────────────────
def load_config(path) -> dict:
    """Load config with deterministic paths, independent of notebook CWD.

    Repository configs default to the repository root because their paths are
    repository-relative. External configs default to their own directory. An
    explicit ``path_base`` (relative to the config file) overrides that rule.
    """
    requested = Path(path).expanduser()
    if requested.is_absolute():
        config_path = requested.resolve()
    else:
        cwd_candidate = (Path.cwd() / requested).resolve()
        project_candidate = (_PROJECT_ROOT / requested).resolve()
        if cwd_candidate.is_file():
            config_path = cwd_candidate
        elif project_candidate.is_file():
            config_path = project_candidate
        else:
            raise FileNotFoundError(f"configuration file not found: {path}")
    try:
        raw = yaml.safe_load(config_path.read_text())
    except yaml.YAMLError as exc:
        raise ValueError(f"invalid YAML in {config_path}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ValueError(f"configuration {config_path} must contain a YAML mapping")

    if raw.get("path_base") is not None:
        configured_base = Path(str(raw["path_base"])).expanduser()
        base = (
            configured_base.resolve() if configured_base.is_absolute()
            else (config_path.parent / configured_base).resolve()
        )
    else:
        try:
            config_path.relative_to(_PROJECT_ROOT)
            base = _PROJECT_ROOT
        except ValueError:
            base = config_path.parent

    cfg = dict(raw)
    for key in (
        "unidock2_bin", "prep_base", "output_dir", "log_dir", "ids_file",
        "scratch_dir", "gpu_lock_file",
    ):
        value = cfg.get(key)
        if value is None:
            continue
        value_path = Path(str(value)).expanduser()
        cfg[key] = str(
            value_path.resolve() if value_path.is_absolute() else (base / value_path).resolve()
        )
    cfg["_config_path"] = str(config_path)
    cfg["_path_base"] = str(base)
    return validate_config(cfg)


def read_authoritative_ids(ids_file: Path):
    ids, seen = [], set()
    for line_number, raw in enumerate(Path(ids_file).read_text().splitlines(), 1):
        pdb_id = raw.strip()
        if not pdb_id or pdb_id.startswith("#"):
            continue
        if not _ID_RE.fullmatch(pdb_id):
            raise ValueError(
                f"invalid benchmark ID in {ids_file}:{line_number}: {pdb_id!r}; "
                "expected NXXX_XXX uppercase form"
            )
        if pdb_id in seen:
            raise ValueError(f"duplicate benchmark ID in {ids_file}: {pdb_id}")
        ids.append(pdb_id)
        seen.add(pdb_id)
    if not ids:
        raise ValueError(f"authoritative IDs file is empty: {ids_file}")
    return ids


def discover_complexes(cfg: dict):
    """Select authoritative IDs; unrelated staged directories are reported by main."""
    prep = Path(cfg["prep_base"])
    if not prep.is_dir():
        raise FileNotFoundError(f"prepared-input directory not found: {prep}")
    exclude = set(cfg.get("exclude") or [])
    staged = [
        path for path in sorted(prep.iterdir())
        if path.is_dir() and not path.name.startswith(("_", ".")) and path.name != "Logs"
    ]
    if cfg.get("ids_file"):
        ids = read_authoritative_ids(Path(cfg["ids_file"]))
        authoritative = set(ids)
        unknown_excludes = sorted(exclude - authoritative)
        if unknown_excludes:
            raise ValueError(f"exclude contains IDs absent from ids_file: {unknown_excludes!r}")
        return [prep / pdb_id for pdb_id in ids if pdb_id not in exclude]
    invalid = [path.name for path in staged if not _ID_RE.fullmatch(path.name)]
    if invalid:
        raise ValueError(f"prepared-input directory contains invalid complex IDs: {invalid!r}")
    return [path for path in staged if path.name not in exclude]


def _est_seconds(n_atoms: int) -> float:
    """Piecewise prep-bound heuristic matching the measured anchors and long tail."""
    n_atoms = int(n_atoms)
    if n_atoms <= 1300:
        return 14.0
    if n_atoms <= 3700:
        return 14.0 + (71.0 / 2400.0) * (n_atoms - 1300)
    return 85.0 + 0.015 * (n_atoms - 3700)


def main(config_path, plan_only=False, limit=None):
    """Plan or run the authoritative benchmark set.

    ``limit`` selects a positive number of *remaining* complexes, not the first
    IDs before completed/skipped entries are removed.
    """
    if not isinstance(plan_only, bool):
        raise ValueError("plan_only must be true or false")
    if limit is not None and (isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0):
        raise ValueError("limit must be None or a positive integer")
    cfg = load_config(config_path)
    unidock2_bin = cfg["unidock2_bin"]
    tool_identity = probe_tool_identity(unidock2_bin)
    Path(cfg["output_dir"]).mkdir(parents=True, exist_ok=True)
    Path(cfg["log_dir"]).mkdir(parents=True, exist_ok=True)

    complexes = discover_complexes(cfg)
    prep = Path(cfg["prep_base"])
    authoritative = (
        set(read_authoritative_ids(Path(cfg["ids_file"])))
        if cfg.get("ids_file") else {cdir.name for cdir in complexes}
    )
    unlisted_staged = sorted(
        path.name for path in prep.iterdir()
        if path.is_dir() and not path.name.startswith(("_", "."))
        and path.name != "Logs" and path.name not in authoritative
    )

    remaining, completed, skipped, plan_errors, giants = [], [], [], [], []
    for cdir in complexes:
        try:
            info, ok = plan_complex(cdir, cfg)
        except Exception as exc:
            plan_errors.append((cdir.name, f"plan error: {exc}"))
            continue
        if not ok:
            if info.get("plan_error"):
                plan_errors.append((cdir.name, info["plan_error"]))
            else:
                skipped.append((cdir.name, info.get("skip", "preflight skip")))
            continue
        try:
            protocol_cfg = effective_unidock2_config(info, cfg)
            provenance, fingerprint = build_provenance(
                info, protocol_cfg, tool_identity, cfg["gpu_lock_file"]
            )
            out_dir = Path(cfg["output_dir"]) / cdir.name
            out_sdf = out_dir / f"{cdir.name}_unidock2_out.sdf"
            done_marker = out_dir / f"{cdir.name}{_DONE_SUFFIX}"
            if (not cfg.get("overwrite", False)
                    and _completion_matches(
                        out_sdf, done_marker, fingerprint, provenance, cfg,
                    )):
                validation = _validate_docking_output(out_sdf, cfg)
                completed.append((cdir, info, fingerprint, validation))
                continue
        except Exception as exc:
            plan_errors.append((cdir.name, f"provenance error: {exc}"))
            continue
        remaining.append((cdir, info, fingerprint))
        if info["input_protein_heavy_atoms"] > 15000:
            giants.append((info["input_protein_heavy_atoms"], cdir.name))

    selected = remaining if limit is None else remaining[:limit]
    deferred = remaining[len(selected):]
    prep_workers = int(cfg["prep_workers"])
    est_prep_s = sum(_est_seconds(info["input_protein_heavy_atoms"]) for _, info, _ in selected)
    # Wall time ≈ the (parallelized) prep sum plus the serial ~few-second GPU shots.
    est_dock_s = len(selected) * _EST_DOCK_SECONDS
    est_wall_s = est_prep_s / max(prep_workers, 1) + est_dock_s
    excluded = sorted(set(cfg.get("exclude") or []))
    print(f"Uni-Dock2 {tool_identity['version']} whole-protein (parallel prep -> serial GPU dock) "
          f"| vina scoring | requested exhaustiveness {cfg['exhaustiveness']} "
          f"| mc_steps {cfg['mc_steps']} | num_pose {cfg['num_pose']} "
          f"| search_mode {cfg['search_mode']} | prep_workers {prep_workers}")
    print(f"  complexes: {len(remaining)} remaining ({len(selected)} selected), "
          f"{len(completed)} already-done, {len(skipped)} skipped, "
          f"{len(plan_errors)} errors, {len(excluded)} excluded "
          f"({', '.join(excluded) or 'none'})")
    if deferred:
        print(f"  limit defers {len(deferred)} otherwise-dockable complexes")
    if unlisted_staged:
        print(f"  ignored {len(unlisted_staged)} staged directories absent from ids_file")
    print(f"  rough estimate: ~{est_wall_s / 3600:.1f} h wall with {prep_workers} parallel "
          f"prep workers (~{est_prep_s / 3600:.1f} h serial receptor prep + "
          f"~{est_dock_s / 60:.0f} min serial GPU docking). "
          f"Per-complex timeout: {cfg['timeout_s']} s")
    if cfg.get("max_receptor_atoms"):
        print(f"  receptor-size cap: skip input protein PDBs with > "
              f"{cfg['max_receptor_atoms']} heavy ATOM records (pre-preparation proxy)")
    if giants:
        print(f"  {len(giants)} large input proteins (>15k heavy ATOM records) may take minutes "
              f"or hit the timeout: {', '.join(name for _, name in sorted(giants, reverse=True))}")
    for name, why in skipped[:20]:
        print(f"    skip {name}: {why}")
    for name, why in plan_errors[:20]:
        print(f"    error {name}: {why}")
    if plan_only:
        return {
            "dockable": len(selected), "remaining": len(remaining),
            "already_done": len(completed), "completed": len(completed),
            "skipped": skipped, "planning_errors": plan_errors,
            "errors": plan_errors, "deferred": len(deferred),
            "unlisted_staged": unlisted_staged,
            "prep_workers": prep_workers,
            "est_hours": round(est_wall_s / 3600, 2),
            "est_prep_hours": round(est_prep_s / 3600, 2),
            "giants": giants,
        }

    results = []
    for cdir, info, fingerprint, validation in completed:
        manifest_path = Path(cfg["output_dir"]) / cdir.name / f"{cdir.name}{_DONE_SUFFIX}"
        out_sdf = Path(cfg["output_dir"]) / cdir.name / f"{cdir.name}_unidock2_out.sdf"
        manifest = _read_json(manifest_path)
        results.append(_result_from_manifest(
            info, manifest, out_sdf, manifest_path, status="done", elapsed_s=0.0,
        ))
    results.extend({"pdb_id": name, "status": "skipped", "skip": why}
                   for name, why in skipped)
    results.extend({"pdb_id": name, "status": "error", "phase": "planning", "error": why}
                   for name, why in plan_errors)
    results.extend({"pdb_id": cdir.name, "status": "deferred",
                    "reason": "outside requested limit"}
                   for cdir, _info, _fingerprint in deferred)

    with _complex_lock(Path(cfg["output_dir"]) / ".unidock2.batch.lock") as batch_acquired:
        if not batch_acquired:
            results.append({
                "pdb_id": "__batch__", "status": "locked",
                "error": "another Uni-Dock2 batch owns this output tree",
            })
        else:
            for result in results:
                if result["status"] != "deferred":
                    summary_path = (
                        Path(cfg["output_dir"]) / result["pdb_id"] / "docking_summary.json"
                    )
                    if result["status"] != "done" or not summary_path.is_file():
                        _atomic_write_json(summary_path, result)

            # Parallel CPU receptor prep feeds the serial GPU docking loop. Prep
            # bundles land in a private driver-scratch tree, never the output tree,
            # so an interrupt or crash simply re-preps the remaining complexes.
            scratch_base = (
                Path(cfg["scratch_dir"]) if cfg.get("scratch_dir")
                else Path("/tmp") / f"unidock2_driver_{os.getuid()}"
            )
            output_token = hashlib.sha256(
                str(Path(cfg["output_dir"])).encode("utf-8")
            ).hexdigest()[:16]
            prep_parent = scratch_base / f"prep_{output_token}"
            _ensure_private_directory(scratch_base)
            _ensure_private_directory(prep_parent)

            executor = (
                ThreadPoolExecutor(max_workers=prep_workers, thread_name_prefix="u2prep")
                if selected else None
            )
            try:
                prep_futures = [
                    executor.submit(
                        prepare_complex, cdir, info, cfg, unidock2_bin,
                        tool_identity, prep_parent,
                    )
                    for cdir, info, _fingerprint in selected
                ] if executor is not None else []
                for idx, (cdir, info, _fingerprint) in enumerate(selected, 1):
                    bundle = None
                    try:
                        bundle = prep_futures[idx - 1].result()
                        if bundle.get("status") != "prepared":
                            result = {
                                "pdb_id": cdir.name,
                                "status": "timeout" if bundle.get("status") == "timeout" else "error",
                                "phase": "prep",
                                "error": bundle.get("error", "receptor preparation failed"),
                                "input_protein_heavy_atoms": info.get("input_protein_heavy_atoms"),
                                "num_poses": 0, "best_affinity": None,
                                "elapsed_s": bundle.get("elapsed_s", 0),
                            }
                        else:
                            result = dock_complex(
                                cdir, cfg, unidock2_bin,
                                tool_identity=tool_identity, prepared=bundle,
                            )
                    except Exception as exc:
                        result = {
                            "pdb_id": cdir.name, "status": "error", "phase": "docking",
                            "error": f"{type(exc).__name__}: {exc}",
                            "num_poses": 0, "best_affinity": None,
                        }
                    finally:
                        if bundle is not None:
                            _cleanup_prep_bundle(bundle)
                    results.append(result)
                    _atomic_write_json(
                        Path(cfg["output_dir"]) / cdir.name / "docking_summary.json", result
                    )
                    icon = {"success": "✓", "timeout": "⏱", "skipped": "–",
                            "error": "✗", "done": "≡", "locked": "🔒"}.get(
                                result["status"], "?"
                            )
                    extra = ""
                    if result["status"] == "success" and result.get("best_affinity") is not None:
                        extra = (f" | {result['num_poses']} poses, best "
                                 f"{result['best_affinity']:.2f} | "
                                 f"{result.get('elapsed_s', 0):.0f}s")
                    elif result["status"] in ("timeout", "error"):
                        extra = (f" | {result.get('input_protein_heavy_atoms', '?')} input "
                                 f"protein atoms | {result.get('elapsed_s', 0):.0f}s")
                    print(f"[{idx}/{len(selected)}] {icon} {cdir.name}{extra}")
                    if result.get("gpu_lock_timed_out"):
                        reason = (
                            "batch stopped after the shared GPU lock wait timed out; "
                            "no additional complexes were launched"
                        )
                        for blocked_cdir, _blocked_info, _blocked_fingerprint in selected[idx:]:
                            blocked = {
                                "pdb_id": blocked_cdir.name,
                                "status": "locked",
                                "phase": "gpu-lock-batch-stop",
                                "gpu_lock_timed_out": True,
                                "gpu_lock_wait_skipped": True,
                                "error": reason,
                                "num_poses": 0,
                                "best_affinity": None,
                            }
                            results.append(blocked)
                            _atomic_write_json(
                                Path(cfg["output_dir"]) / blocked_cdir.name /
                                "docking_summary.json",
                                blocked,
                            )
                        break
            finally:
                if executor is not None:
                    _kill_active_prep()
                    executor.shutdown(wait=True, cancel_futures=True)
                _rmtree_quiet(prep_parent)

    counts = {status: sum(1 for result in results if result["status"] == status)
              for status in ("success", "done", "timeout", "skipped", "error", "locked", "deferred")}
    print(f"\n{'=' * 60}\nUNI-DOCK2 COMPLETE\n  success: {counts['success']}  "
          f"already-done: {counts['done']}  timeout: {counts['timeout']}  "
          f"skipped: {counts['skipped']}  error: {counts['error']}  "
          f"locked: {counts['locked']}  deferred: {counts['deferred']}")
    if counts["timeout"] or counts["error"] or counts["locked"]:
        print("  (re-run to retry timeout/error/locked complexes; completed manifests are resumed)")
    print(f"  results: {cfg['output_dir']}\n{'=' * 60}")
    return results


def results_exit_code(results) -> int:
    if isinstance(results, list):
        if any(result.get("status") in {"timeout", "error", "locked"} for result in results):
            return 1
    elif isinstance(results, dict) and (results.get("planning_errors") or results.get("errors")):
        return 1
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Whole-protein Uni-Dock2 docking (Benchmark set)")
    ap.add_argument("-c", "--config", required=True)
    ap.add_argument("--plan-only", action="store_true", help="print the workload + estimate, do not dock")
    ap.add_argument("--limit", type=int, default=None, help="dock only the first N complexes (testing)")
    a = ap.parse_args()
    outcome = main(a.config, plan_only=a.plan_only, limit=a.limit)
    raise SystemExit(results_exit_code(outcome))
