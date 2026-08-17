"""PoseBusters Parallel Pose Validation Pipeline.

Validates docked poses from one or more docking methods using PoseBusters.
Supports AutoDock Vina and Vinardo scoring (raw and gnina/smina-optimized
variants), tiled Uni-Dock, single-shot Uni-Dock2, DiffDock, and EquiBind
outputs."""



from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import re
import shutil
import pickle
import select
import signal
import struct
import subprocess
import sys
import threading
import time
import warnings
from collections import Counter
from dataclasses import dataclass, field
from functools import lru_cache
from multiprocessing import Pool, TimeoutError as MPTimeoutError, cpu_count, get_context
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from matplotlib.colors import LinearSegmentedColormap
from posebusters import PoseBusters

# Legacy shared-receptor workflows may reuse this curated cofactor/metal
# allowlist (see Scripts/Utilities/inject_hetatms.py). It is not provenance and
# is therefore forbidden by the benchmark's ``prepared_per_pose`` policy.
# Import is best-effort for backwards-compatible non-benchmark configurations.
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
try:
    from Scripts.Utilities.inject_hetatms import (
        DEFAULT_KEEP_COFACTORS, DEFAULT_KEEP_METALS,
    )
    _DEFAULT_KEEP_HETATM: frozenset[str] = frozenset(DEFAULT_KEEP_COFACTORS) | frozenset(DEFAULT_KEEP_METALS)
except Exception:  # pragma: no cover - optional dependency on repo layout
    _DEFAULT_KEEP_HETATM = frozenset()


# ============================================================================
# RDKIT LOG NOISE
# ============================================================================

# PoseBusters loads every pose through RDKit, which emits one warning per
# molecule for SDFs tagged 2D but carrying real (non-zero Z) coordinates:
#   "Warning: molecule is tagged as 2D, but at least one Z coordinate is not
#    zero. Marking the mol as 3D."
# It is harmless (RDKit just promotes the mol to 3D, which is what we want) but
# drowns the real progress output. Route RDKit's C++ logs through Python stderr
# and drop only this one message — every other RDKit warning/error (unparseable
# elements, sanitization failures, ...) still gets through.
_RDKIT_SUPPRESS_SUBSTR = "is tagged as 2D, but at least one Z coordinate is not zero"


class _RDKitStderrFilter:
    """stderr wrapper that swallows only the harmless 2D/3D RDKit warning."""

    def __init__(self, stream):
        self._stream = stream

    def write(self, msg):
        if _RDKIT_SUPPRESS_SUBSTR not in msg:
            self._stream.write(msg)
        return len(msg)

    def flush(self):
        self._stream.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)


def _silence_rdkit_2d3d_warning() -> None:
    """Install the 2D/3D-warning filter on this process's stderr (idempotent).

    Must run in every process that loads molecules: the main process (PDBQT→SDF
    conversion) and each Pool worker (where the ``bust`` calls live).
    """
    from rdkit import rdBase

    rdBase.LogToPythonStderr()
    if not isinstance(sys.stderr, _RDKitStderrFilter):
        sys.stderr = _RDKitStderrFilter(sys.stderr)


_silence_rdkit_2d3d_warning()


# ============================================================================
# QUIET MODE  (suppress warning chatter) + WORKER THREAD LIMITS
# ============================================================================

# Module-global set once at startup (main process) and re-applied in each Pool
# worker via _init_worker. Gates the pipeline's own verbose per-file listings.
_QUIET = False


def _apply_quiet_mode(quiet: bool) -> None:
    """Silence Python warnings and *all* RDKit logging when *quiet* is set.

    Complements _silence_rdkit_2d3d_warning (which drops only the one harmless
    2D/3D message): quiet mode drops every Python warning and every RDKit
    warning/error log line. Runs in the main process and in every Pool worker.
    Idempotent; a no-op when *quiet* is False.
    """
    global _QUIET
    _QUIET = bool(quiet)
    if not quiet:
        return
    warnings.filterwarnings("ignore")
    os.environ.setdefault("PYTHONWARNINGS", "ignore")
    try:
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
    except Exception:  # pragma: no cover - RDKit always present in this env
        pass


def _limit_worker_threads() -> None:
    """Constrain per-worker thread fan-out so N Pool workers don't oversubscribe.

    Two independent sources of hidden intra-worker parallelism blow up CPU
    contention when the pipeline already runs one worker per core:

    1. NumPy/pandas/BLAS honour OMP/MKL/OPENBLAS/NUMEXPR thread-count env vars.
    2. PoseBusters' ``internal_energy`` check embeds + UFF-minimises an ensemble
       of 50 conformers per pose (posebusters/modules/energy_ratio.py) with
       ``num_threads=0`` — RDKit's "use every core" — so each of N workers spawns
       an all-core RDKit thread pool. That is N x cores threads fighting over
       cores, which degrades throughput catastrophically under load.

    We already parallelise at the pose level, so each worker should stay
    single-threaded. Env caps handle (1); a guarded monkeypatch forcing RDKit's
    conformer generation to one thread handles (2). Both fail safe.
    """
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.setdefault(var, "1")

    try:
        from posebusters.modules import energy_ratio as _er

        _orig_new_conformation = _er.new_conformation

        def _single_thread_new_conformation(mol, n_confs=1, num_threads=0,
                                            energy_minimization=True):
            # Force one RDKit thread regardless of what PoseBusters requests.
            return _orig_new_conformation(mol, n_confs, 1, energy_minimization)

        # get_energies() looks up new_conformation as a module global, so
        # replacing the attribute is picked up on the next call.
        if getattr(_er.new_conformation, "__name__", "") != "_single_thread_new_conformation":
            _er.new_conformation = _single_thread_new_conformation
    except Exception:  # pragma: no cover - posebusters internals may change
        pass


# ============================================================================
# CONFIG
# ============================================================================

@dataclass
class PipelineConfig:
    """Resolved pipeline configuration."""

    receptors_dir: Path                 # folder of receptor PDB files
    docking_directories: dict[str, Path]  # method_key -> ligand/poses folder
    work_dir: Path                      # working dir (used as ligand-template search root)
    output_base_dir: Path               # base output folder
    config_mode: str = "dock"           # "dock" or "mol"
    save_interval: int = 100
    overwrite: bool = False
    num_workers: int | None = None      # None -> all cores
    ligand_template_dirs: list[Path] = field(default_factory=list)
    make_plots: bool = True
    copy_proved_poses: bool = True
    # Receptor selection policy. ``prepared_per_pose`` resolves the receptor
    # artifact actually used by each docking method/variant and records its
    # content digest on every result. ``legacy_shared`` retains the historical
    # protein-name lookup under ``receptors_dir`` for non-benchmark workflows.
    receptor_policy: str = "legacy_shared"
    validation_receptor_directories: dict[str, Path] = field(default_factory=dict)
    require_all_methods: bool = False
    require_validation_receptors: bool = True
    require_complete_results: bool = False
    expected_common_combinations: int | None = None
    # Methods that are collected + busted but do NOT constrain the common-combo
    # intersection. The shared benchmark cohort is the intersection over the
    # REQUIRED (non-optional) methods; an optional method contributes its poses
    # for whichever cohort complexes it happens to have docked, so a partially
    # docked tool (or a brand-new one still being run) never collapses the whole
    # cohort down to its own subset. See filter_common_combos.
    optional_methods: set[str] = field(default_factory=set)
    # HETATM residue names to strip from each receptor PDB before validation
    # (e.g. crystallisation artefacts that the docker never saw but PoseBusters
    # otherwise checks against). Stored as a set of upper-cased 3-letter codes.
    strip_hetatm_residues: set[str] = field(default_factory=set)
    # Auto-derive mode: keep ONLY these HETATM residues and strip every other
    # HETATM from each receptor PDB before validation. This is retained only for
    # legacy shared-receptor workflows; residue names alone cannot prove docking
    # provenance. ``None`` disables keep-mode and falls back to the explicit
    # strip list. Takes precedence when set.
    keep_hetatm_residues: set[str] | None = None

    # Suppress warning chatter (Python warnings + all RDKit logs) and the
    # verbose per-file listings (HETATM strip table, protein file dump,
    # per-model conversion notes). Progress / checkpoint / summary lines stay.
    quiet: bool = False
    # Recycle each Pool worker after this many poses so RDKit conformer-cache /
    # RSS growth over a long run can't accumulate into swap thrashing — the
    # classic "stalls after a few thousand poses". None -> workers live forever.
    worker_maxtasks: int | None = 200
    # Per-pose wall-clock cap (seconds) for a single bust() call. A pathological
    # ligand can otherwise send the 50-conformer UFF ensemble into an effectively
    # unbounded compute and stall a worker. 0/None disables the guard.
    pose_timeout: int | None = 600
    # Cap intra-worker thread fan-out (BLAS + RDKit conformer generation) so N
    # Pool workers don't each try to use every core. See _limit_worker_threads.
    limit_worker_threads: bool = True
    # Live resource monitor: a background thread prints core utilisation per
    # process + memory/swap every ``monitor_interval`` seconds and flags
    # oversubscription / swap pressure / stalls. Off by default (diagnostic).
    monitor: bool = False
    monitor_interval: int = 30
    # Stall recovery. A dead worker (e.g. OOM-killed under memory pressure) makes
    # Pool.imap_unordered hang forever waiting for a result that never comes — the
    # "stuck after a few thousand busts" with no error. If no pose completes for
    # ``stall_timeout`` seconds the pool is declared stuck, terminated, and a
    # fresh pool is restarted on the not-yet-done poses (up to ``max_restarts``
    # times). None -> auto (max(pose_timeout, 300) x 3). 0 disables detection.
    stall_timeout: int | None = None
    max_restarts: int = 5

    # Variant restriction: keep only the listed variant(s) per tool when
    # collecting poses. AutoDock and DiffDock use optimizer
    # (original/smina/gnina), EquiBind uses refine_variant (raw/smina/gnina),
    # and tiled Uni-Dock uses variant=tiled. Empty -> keep every variant
    # (default). See _apply_variant_filter / config key `variant_filter`.
    variant_filter: dict[str, set[str]] = field(default_factory=dict)

    # Derived
    output_dir: Path = field(init=False)
    converted_dir: Path = field(init=False)
    validation_receptors_dir: Path = field(init=False)

    def __post_init__(self):
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = self.output_base_dir / self.config_mode
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.converted_dir = self.output_dir / "converted_pdbqt"
        self.converted_dir.mkdir(parents=True, exist_ok=True)
        self.validation_receptors_dir = self.output_dir / "validation_receptors"


# Provenance field carrying each tool's variant label (used by variant_filter).
_VARIANT_FIELD = {
    "autodock": "optimizer",
    "unidock": "variant",
    "diffdock": "optimizer",
    "equibind": "refine_variant",
}


def _variant_base_method(key: str) -> str:
    """Normalise a method/tool key to its base tool name. Every EquiBind flavour
    (equibind_guided / equibind_exclusion / equibind_docked_poses / …) collapses
    to 'equibind'; every tiled Uni-Dock alias collapses to 'unidock'; Uni-Dock2
    stays 'unidock2'; every AutoDock scoring flavour (autodock / autodock_vinardo)
    collapses to 'autodock' so they share the receptor resolver + optimizer
    variant axis; diffdock maps to itself.

    Order matters: the unidock2 test precedes the unidock test because
    ``"unidock2".startswith("unidock")`` is also true."""
    k = str(key).strip().lower()
    if k.startswith("equibind"):
        return "equibind"
    if k.startswith("unidock2") or k in ("uni-dock2", "uni-dock-2"):
        return "unidock2"
    if k.startswith("unidock") or k == "uni-dock":
        return "unidock"
    if k.startswith("autodock"):
        return "autodock"
    return k


def load_config(config_path: str | Path) -> PipelineConfig:
    """Load and validate a YAML config file into a PipelineConfig."""
    config_path = Path(config_path)
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path) as f:
        raw: dict[str, Any] = yaml.safe_load(f) or {}

    work_dir = Path(raw.get("work_dir", ".")).expanduser().resolve()

    def _resolve(p: str | Path) -> Path:
        p = Path(p).expanduser()
        return p if p.is_absolute() else (work_dir / p).resolve()

    if "receptors_dir" not in raw:
        raise KeyError("Config must define 'receptors_dir'")
    receptors_dir = _resolve(raw["receptors_dir"])

    docking_raw = raw.get("docking_directories") or {}
    if not docking_raw:
        raise KeyError("Config must define at least one entry under 'docking_directories'")
    docking_directories = {k: _resolve(v) for k, v in docking_raw.items()}

    receptor_policy = str(raw.get("receptor_policy", "legacy_shared")).strip().lower()
    if receptor_policy not in {"legacy_shared", "prepared_per_pose"}:
        raise ValueError(
            "receptor_policy must be 'legacy_shared' or 'prepared_per_pose' "
            f"(got {receptor_policy!r})"
        )
    validation_receptor_directories = {
        str(k): _resolve(v)
        for k, v in (raw.get("validation_receptor_directories") or {}).items()
    }
    if receptor_policy == "prepared_per_pose":
        missing_roots = sorted(set(docking_directories) - set(validation_receptor_directories))
        if missing_roots:
            raise KeyError(
                "prepared_per_pose requires validation_receptor_directories entries for: "
                + ", ".join(missing_roots)
            )

    output_base_dir = _resolve(raw.get("output_base_dir", "posebusters_results"))
    config_mode = str(raw.get("config_mode", "dock")).strip().lower()
    if config_mode not in {"dock", "mol"}:
        raise ValueError(f"config_mode must be 'dock' or 'mol' (got {config_mode!r})")
    require_validation_receptors = bool(raw.get("require_validation_receptors", True))
    if config_mode == "dock" and not require_validation_receptors:
        raise ValueError(
            "dock mode requires require_validation_receptors: true; ligand-only "
            "fallback is intentionally unsupported"
        )

    template_dirs = [_resolve(p) for p in (raw.get("ligand_template_dirs") or [])]

    strip_raw = raw.get("strip_hetatm_residues") or []
    strip_set = {str(r).strip().upper() for r in strip_raw if str(r).strip()}

    # keep_hetatm_residues: None/false -> off (use strip list);
    # True/"auto" -> the injected cofactor+metal set (mirrors the docking
    # config's `retain_hetatm_residues: true`); explicit list -> those codes.
    keep_raw = raw.get("keep_hetatm_residues", None)
    if keep_raw in (None, False):
        keep_set: set[str] | None = None
    elif keep_raw is True or (isinstance(keep_raw, str) and keep_raw.strip().lower() == "auto"):
        if not _DEFAULT_KEEP_HETATM:
            raise RuntimeError(
                "keep_hetatm_residues: auto requested but the injected cofactor "
                "set could not be imported from Scripts/Utilities/inject_hetatms.py. "
                "Provide an explicit list instead."
            )
        keep_set = set(_DEFAULT_KEEP_HETATM)
    else:
        keep_set = {str(r).strip().upper() for r in keep_raw if str(r).strip()}
    if receptor_policy == "prepared_per_pose" and (strip_set or keep_set is not None):
        raise ValueError(
            "strip_hetatm_residues/keep_hetatm_residues cannot be combined with "
            "prepared_per_pose; receptor membership must come from provenance"
        )

    # variant_filter: restrict which per-tool variant is validated. Keys are the
    # docking_directories method keys (or base tool names); values are the variant
    # label(s) to KEEP (a scalar or a list). AutoDock/DiffDock use optimizer
    # (original/smina/gnina), EquiBind uses refine_variant (raw/smina/gnina), and
    # Uni-Dock uses variant=tiled. A tool not listed keeps every variant.
    # See _apply_variant_filter.
    vf_raw = raw.get("variant_filter") or {}
    variant_filter: dict[str, set[str]] = {}
    for _k, _v in vf_raw.items():
        base = _variant_base_method(_k)
        if isinstance(_v, (list, tuple, set)):
            vals = {str(x).strip().lower() for x in _v if str(x).strip()}
        elif _v in (None, "", False):
            vals = set()
        else:
            vals = {str(_v).strip().lower()}
        if vals:
            variant_filter[base] = variant_filter.get(base, set()) | vals

    return PipelineConfig(
        receptors_dir=receptors_dir,
        docking_directories=docking_directories,
        work_dir=work_dir,
        output_base_dir=output_base_dir,
        config_mode=config_mode,
        save_interval=int(raw.get("save_interval", 100)),
        overwrite=bool(raw.get("overwrite", False)),
        num_workers=raw.get("num_workers"),
        ligand_template_dirs=template_dirs,
        make_plots=bool(raw.get("make_plots", True)),
        copy_proved_poses=bool(raw.get("copy_proved_poses", True)),
        receptor_policy=receptor_policy,
        validation_receptor_directories=validation_receptor_directories,
        require_all_methods=bool(raw.get("require_all_methods", False)),
        require_validation_receptors=require_validation_receptors,
        require_complete_results=bool(raw.get("require_complete_results", False)),
        expected_common_combinations=(
            None if raw.get("expected_common_combinations") is None
            else int(raw["expected_common_combinations"])
        ),
        optional_methods={
            str(m).strip() for m in (raw.get("optional_methods") or []) if str(m).strip()
        },
        strip_hetatm_residues=strip_set,
        keep_hetatm_residues=keep_set,
        quiet=bool(raw.get("quiet", False)),
        worker_maxtasks=(None if raw.get("worker_maxtasks", 200) in (None, 0, False)
                         else int(raw.get("worker_maxtasks", 200))),
        pose_timeout=(None if raw.get("pose_timeout", 600) in (None, 0, False)
                      else int(raw.get("pose_timeout", 600))),
        limit_worker_threads=bool(raw.get("limit_worker_threads", True)),
        monitor=bool(raw.get("monitor", False)),
        monitor_interval=int(raw.get("monitor_interval", 30)),
        stall_timeout=(None if raw.get("stall_timeout", None) is None
                       else int(raw.get("stall_timeout"))),
        max_restarts=int(raw.get("max_restarts", 5)),
        variant_filter=variant_filter,
    )


# ============================================================================
# PROTEIN DISCOVERY
# ============================================================================

_NAME_SUFFIXES = ("_protein", "_clean", "_cleaned", "_receptor")


def strip_hetatm_residues_from_pdb(
    src: Path, dst: Path,
    drop_resnames: set[str] | None = None,
    keep_resnames: set[str] | None = None,
) -> tuple[int, int]:
    """Write *dst* as a copy of *src* with HETATM records filtered.

    With *keep_resnames*: drop every HETATM whose 3-letter residue name is NOT
    in the set (legacy allowlist mode).
    Otherwise with *drop_resnames*: drop HETATM whose residue name IS in the set
    (legacy explicit-strip mode). ATOM and other records pass through unchanged.
    Returns (kept_lines, dropped_lines).
    """
    kept = dropped = 0
    out_lines: list[str] = []
    for line in src.read_text().splitlines(keepends=True):
        if line.startswith("HETATM"):
            resn = line[17:20].strip().upper()
            drop = (resn not in keep_resnames) if keep_resnames is not None \
                else (resn in (drop_resnames or set()))
            if drop:
                dropped += 1
                continue
            kept += 1
        elif line.startswith("ATOM"):
            kept += 1
        out_lines.append(line)
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text("".join(out_lines))
    return kept, dropped


def prepare_cleaned_receptor_dir(ctx: PipelineConfig) -> Path:
    """Write HETATM-filtered copies of every PDB under ``ctx.receptors_dir`` to
    ``ctx.output_dir / hetatm_cleaned`` and return that directory.

    Keep-mode (``ctx.keep_hetatm_residues`` set) takes precedence: only those
    residues survive and every other HETATM is stripped. This legacy heuristic
    does not prove what a docker saw and is never used by ``prepared_per_pose``.
    Otherwise the explicit strip list (``ctx.strip_hetatm_residues``) is applied.
    If neither is configured, the original directory is returned unchanged.
    """
    keep = ctx.keep_hetatm_residues
    drop = ctx.strip_hetatm_residues
    if keep is None and not drop:
        return ctx.receptors_dir

    cleaned_dir = ctx.output_dir / "hetatm_cleaned"
    cleaned_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 80)
    if keep is not None:
        print(f"HETATM keep-mode: keeping ONLY {sorted(keep)}; stripping all other HETATM")
    else:
        print(f"Stripping HETATM residues: {sorted(drop)}")
    print("=" * 80)

    total_kept = total_dropped = 0
    n_files = 0
    for src in sorted(Path(ctx.receptors_dir).glob("**/*.pdb")):
        dst = cleaned_dir / src.relative_to(ctx.receptors_dir)
        # Always rewrite — the keep/strip set may have changed since last run
        if keep is not None:
            kept, dropped = strip_hetatm_residues_from_pdb(
                src.resolve(), dst, keep_resnames=keep)
        else:
            kept, dropped = strip_hetatm_residues_from_pdb(
                src.resolve(), dst, drop_resnames=drop)
        total_kept += kept
        total_dropped += dropped
        n_files += 1
        if dropped and not _QUIET:
            print(f"  {src.name}: kept {kept}, dropped {dropped} HETATM atoms")

    print(f"\n  Receptors cleaned: {n_files}  "
          f"(total atoms kept {total_kept}, dropped {total_dropped})")
    print(f"  Cleaned receptors → {cleaned_dir}")
    return cleaned_dir


def discover_proteins(receptors_dir: Path) -> tuple[dict[str, str], dict[str, str]]:
    """Build {stem -> path} and a normalized variant map from receptor PDBs."""
    pdbs = sorted(Path(receptors_dir).glob("**/*.pdb"))
    file_map = {p.stem: str(p) for p in pdbs}
    normalized_map = {stem.replace("-", "_").lower(): path for stem, path in file_map.items()}

    if not _QUIET:
        print("=" * 80)
        print("AVAILABLE PROTEIN PDB FILES")
        print("=" * 80)
        for stem, path in sorted(file_map.items()):
            print(f"  {stem}  →  {path}")
    print(f"\nTotal: {len(file_map)} PDB files found in {receptors_dir}")
    return file_map, normalized_map


def find_protein_file(
    protein_name: str,
    file_map: dict[str, str],
    normalized_map: dict[str, str],
) -> str | None:
    """Resolve a PDB file path for *protein_name* using flexible matching."""
    # 1. Exact match
    if protein_name in file_map:
        return file_map[protein_name]

    # 2. Common suffix variants (add and strip)
    for suffix in _NAME_SUFFIXES:
        candidate = protein_name + suffix
        if candidate in file_map:
            return file_map[candidate]
        if protein_name.endswith(suffix):
            stripped = protein_name[: -len(suffix)]
            if stripped in file_map:
                return file_map[stripped]

    # 3. Substring containment (shortest stem wins)
    name_lower = protein_name.lower()
    matches = [(s, p) for s, p in file_map.items() if name_lower in s.lower()]
    if matches:
        return min(matches, key=lambda t: len(t[0]))[1]

    # 4. Normalized
    norm = protein_name.replace("-", "_").lower()
    if norm in normalized_map:
        return normalized_map[norm]
    for nstem, path in normalized_map.items():
        if norm in nstem:
            return path

    return None


# ============================================================================
# TEST-COLUMN DETECTION
# ============================================================================

_METADATA_COLS = {
    "file_path", "filepath", "file", "path", "sdf_file", "sdf_path",
    "method", "docking_method", "protein", "ligand", "pose_rank", "rank",
    "molecule", "mol_name", "name", "complex", "protein_path", "ligand_path",
    "mol_pred", "mol_true", "mol_cond",
    "receptor_file", "receptor_file_used", "receptor_source_file",
    "source_pose_file", "source_pose_sha256", "receptor_source_sha256",
    "receptor_stage", "receptor_sha256", "pose_sha256",
    "validation_signature", "posebusters_runtime_signature",
    # EquiBind variant provenance (never a pass/fail test column).
    "pocket_source", "pocket_id", "clamp_variant", "refine_variant",
    "smina_affinity", "gnina_affinity",
    # AutoDock Vina native rank + affinity.
    "autodock_rank", "autodock_affinity", "optimized_rank", "rank_metric",
    "minimized_affinity", "cnn_score", "cnn_affinity", "optimizer_elapsed_time_s",
    "optimizer_processing_elapsed_time_s",
    # AutoDock/DiffDock post-pose optimizer provenance (original / smina / gnina).
    "optimizer", "optimization_log_file", "optimizer_provenance_file",
    "optimizer_provenance_fingerprint", "source_pose_model_sha256",
    # Tiled Uni-Dock native rank/affinity and provenance.
    "unidock_rank", "unidock_affinity", "variant", "scoring", "engine",
    "docking_receptor_file", "docking_summary_file",
    "unidock_manifest_file", "unidock_fingerprint", "unidock_generation_id",
    # Uni-Dock2 native rank + affinity (vina_binding_free_energy).
    "unidock2_rank", "unidock2_affinity", "unidock2_completion_file",
    "unidock2_fingerprint", "unidock2_generation_id", "unidock2_output_sha256",
    "unidock2_prepared_receptor_file", "unidock2_prepared_receptor_sha256",
    "unidock2_prepared_receptor_atoms",
    "unidock2_prepared_receptor_heavy_atoms",
    "unidock2_engine_receptor_heavy_atoms_in_box",
}

_EXCLUDE_COLS = {
    "mol_true_loaded",
    "number_short_outlier_bonds", "number_long_outlier_bonds",
    "number_outlier_angles", "number_clashes",
    "number_non-aromatic_rings_pass", "number_aromatic_rings_pass",
    "number_non-aromatic_rings_checked", "number_aromatic_rings_checked",
    "number_double_bonds_checked", "number_double_bonds_pass",
    "number_valid_bonds", "number_valid_angles", "number_valid_noncov_pairs",
    "number_noncov_pairs", "number_bonds", "number_angles",
    "num_h_added",
    "not_too_far_away_organic_cofactors",
    "not_too_far_away_inorganic_cofactors",
    "not_too_far_away_waters",
}

_BOOL_LIKE_VALUES = frozenset({True, False, 1, 0, 1.0, 0.0,
                               "True", "False", "true", "false"})
_BOOL_STR_MAP = {"True": True, "true": True, "False": False, "false": False}

# The canonical PoseBusters pass/fail tests — the ONLY columns that define a
# pose's validity verdict. With ``full_report=True`` the result table also
# carries many diagnostic columns (e.g. ``most_extreme_clash_protein``,
# ``volume_overlap_<group>``, ``*_maximum_distance_from_plane``); several are
# boolean and would otherwise be mistaken for tests by the heuristic below —
# notably ``most_extreme_clash_protein`` (always False) which silently fails
# every pose. Restricting to this allowlist keeps the verdict correct while
# the diagnostic columns remain in the CSV for inspection.
# These are the complete stock PoseBusters verdicts for v0.6.x. Project-specific
# concepts such as allosteric-site acceptance must be reported as separate
# metrics rather than silently changing PoseBusters' all-pass result.
MOLECULE_TEST_COLUMNS = (
    "mol_pred_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
)

CANONICAL_TEST_COLUMNS = (
    "mol_pred_loaded", "mol_cond_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "protein-ligand_maximum_distance",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
)


def expected_test_columns(config_mode: str) -> list[str]:
    """Return the exact stock PoseBusters pass/fail columns for a mode."""
    mode = str(config_mode).strip().lower()
    if mode == "mol":
        return list(MOLECULE_TEST_COLUMNS)
    if mode == "dock":
        return list(CANONICAL_TEST_COLUMNS)
    raise ValueError(f"Unsupported PoseBusters config mode: {config_mode!r}")


def require_test_columns(df: pd.DataFrame, config_mode: str) -> list[str]:
    """Require every stock test column; partial reports must never pass."""
    expected = expected_test_columns(config_mode)
    missing = [column for column in expected if column not in df.columns]
    if missing:
        raise ValueError(
            f"PoseBusters {config_mode!r} results are missing {len(missing)} "
            "required stock test column(s): " + ", ".join(missing)
        )
    return expected


def identify_test_columns(df: pd.DataFrame) -> list[str]:
    """Return the list of boolean PoseBusters test columns in *df*.

    Prefer the canonical PoseBusters test set (those present in *df*). Only if
    none are found — e.g. a CSV from a non-PoseBusters tool — fall back to the
    older boolean-detection heuristic.
    """
    canonical = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
    if canonical:
        return canonical

    test_cols = []
    for col in df.columns:
        cl = col.lower().strip()
        if cl in _METADATA_COLS or col in _EXCLUDE_COLS or cl in _EXCLUDE_COLS:
            continue
        if cl.startswith("number_") or cl.startswith("num_"):
            continue
        unique_vals = set(df[col].dropna().unique())
        if unique_vals.issubset(_BOOL_LIKE_VALUES):
            test_cols.append(col)
    if not test_cols:
        test_cols = [c for c in df.columns if df[c].dtype == bool]
    return test_cols


def _strict_bool(value: Any) -> bool:
    """Coerce one PoseBusters verdict without making nulls truthy."""
    if value is None or (not isinstance(value, (str, bytes)) and pd.isna(value)):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        if value in (0, 0.0):
            return False
        if value in (1, 1.0):
            return True
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"true", "1"}:
            return True
        if normalized in {"false", "0", "", "nan", "none", "null"}:
            return False
    raise ValueError(f"Unrecognized PoseBusters boolean value: {value!r}")


def coerce_test_cols_to_bool(df: pd.DataFrame, test_cols: list[str]) -> None:
    """Strictly convert stock verdicts in-place; missing values are failures."""
    for tc in test_cols:
        df[tc] = df[tc].map(_strict_bool).astype(bool)


# ============================================================================
# POSE ROW COLLECTORS  (one per docking method)
# Each returns rows with: docking_tool, protein, ligand, file_path, pose_count
# ============================================================================

# Method-specific suffixes that get appended to receptor / ligand stems by the
# upstream docking pipelines. We strip them so the (protein, ligand) tuple
# emitted by every collector aligns and `filter_common_combos` can intersect.
_PROTEIN_SUFFIXES = ("_cleaned", "_clean", "_protein", "_receptor")
_LIGAND_SUFFIXES = (
    "_ligand_start_conf_vina",   # AutoDock (post-meeko/mgltools naming)
    "_ligand_start_conf",        # EquiBind
    "_start_conf_vina",
    "_start_conf",               # DiffDock
    "_ligand",
    "_vina",
)


def _strip_suffixes(name: str, suffixes: tuple[str, ...]) -> str:
    """Iteratively strip any of *suffixes* from the end of *name*."""
    changed = True
    while changed:
        changed = False
        for s in suffixes:
            if name.endswith(s) and len(name) > len(s):
                name = name[: -len(s)]
                changed = True
                break
    return name


def _normalize_protein(name: str) -> str:
    # Direct Vina artifacts carry the preparation backend in the receptor stem,
    # e.g. ``7SUC_COM_protein_mgl_tools`` or ``*_protein_meeko_retry1``.
    name = re.sub(r"_(?:mgl_tools|meeko(?:_retry\d+)?)$", "", name)
    return _strip_suffixes(name, _PROTEIN_SUFFIXES)


def _normalize_ligand(name: str) -> str:
    return _strip_suffixes(name, _LIGAND_SUFFIXES)


def _diffdock_optimizer_of(sdf_file: Path) -> str:
    """A DiffDock pose's optimizer, read from its path: 'original' for a top-level
    rank*.sdf, else the tool of its optimized_<tool>/ subfolder (smina / gnina).
    Mirrors the tagging in _expand_diffdock_poses so counts and busting agree."""
    for part in sdf_file.parts:
        if part.startswith("optimized_"):
            return part[len("optimized_"):] or "original"
    return "original"


# DiffDock writes its top pose twice: a bare ``rankN.sdf`` that is byte-identical
# to the ``rankN_confidence-*.sdf`` beside it. Only confidence-named poses are
# re-optimised, so no bare copy exists in optimized_<tool>/ subfolders. The
# recursive ``**/*.sdf`` collectors must skip the confidence-less bare copy, or
# every DiffDock complex gains a duplicate top-ranked (rank-1) pose. Mirrors the
# drop the analysis scripts already apply (docking_effort_comparison,
# posebusters_pose_comparison, pose_cluster_crystal_pocket_report).
_DIFFDOCK_BARE_RANK = re.compile(r"^rank\d+\.sdf$", re.IGNORECASE)


def _is_diffdock_bare_rank(sdf_file: Path) -> bool:
    """True for a confidence-less bare ``rankN.sdf`` (byte-identical duplicate of
    its ``rankN_confidence-*.sdf`` sibling)."""
    return bool(_DIFFDOCK_BARE_RANK.match(sdf_file.name))


def _equibind_refine_of(sdf_file: Path) -> str | None:
    """A EquiBind pose's refine_variant (raw / smina / gnina): the SDF tag wins,
    else the filename suffix (__refSMINA / __refRAW). None when neither is present.
    Mirrors _expand_equibind_poses so counts and busting agree."""
    tag = _read_sdf_tags(sdf_file, ("refine_variant",)).get("refine_variant")
    return tag or _equibind_refine_variant(sdf_file.name)


def _variant_ok(value: str | None, keep: set[str] | None) -> bool:
    """True if a pose's variant *value* is in the *keep* set (case-insensitive).
    A falsy *keep* keeps every pose; a missing/empty *value* is dropped when a
    non-empty keep set is given."""
    if not keep:
        return True
    v = str(value).strip().lower() if value not in (None, "") else None
    return v in keep


def _autodock_variant_ok(value: str, keep: set[str] | None) -> bool:
    """Variant predicate with friendly aliases for raw AutoDock output."""
    if not keep:
        return True
    aliases = {"raw": "original", "native": "original", "none": "original"}
    wanted = {aliases.get(str(v).strip().lower(), str(v).strip().lower()) for v in keep}
    value = aliases.get(str(value).strip().lower(), str(value).strip().lower())
    return value in wanted


def _closed_scored_pdbqt_models(
    path: Path,
    *,
    allow_unwrapped_single: bool = False,
) -> list[tuple[int, float]] | None:
    """Return closed ``(MODEL, affinity)`` records, or ``None`` if malformed.

    A completion check must validate the sequence, not merely search for the two
    marker substrings: ENDMDL-before-MODEL and a complete first model followed by
    a truncated second model are both invalid. Uni-Dock merged output also needs
    one score per model so a scoreless file can never become a completed run.
    """
    try:
        lines = Path(path).read_text(errors="ignore").splitlines()
    except OSError:
        return None
    models: list[tuple[int, float]] = []
    in_model = False
    model_number: int | None = None
    affinity: float | None = None
    atom_count = 0
    saw_marker = False
    for line in lines:
        if line.startswith("MODEL"):
            saw_marker = True
            if in_model:
                return None
            fields = line.split()
            try:
                model_number = int(fields[1])
            except (IndexError, ValueError):
                return None
            if model_number <= 0:
                return None
            in_model, affinity, atom_count = True, None, 0
            continue
        if line.startswith("ENDMDL"):
            saw_marker = True
            if not in_model or model_number is None or affinity is None or atom_count <= 0:
                return None
            models.append((model_number, affinity))
            in_model, model_number, affinity = False, None, None
            continue
        if in_model and line.startswith("REMARK VINA RESULT:") and affinity is None:
            match = re.search(r"REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)", line)
            if match:
                affinity = float(match.group(1))
        if in_model and line.startswith(("ATOM", "HETATM")):
            atom_count += 1
    if in_model or (saw_marker and not models):
        return None
    if models:
        ranks = [rank for rank, _ in models]
        return models if len(ranks) == len(set(ranks)) else None
    if not allow_unwrapped_single:
        return None
    text = "\n".join(lines)
    match = re.search(r"REMARK VINA RESULT:\s*(-?\d+(?:\.\d+)?)", text)
    has_atom = any(line.startswith(("ATOM", "HETATM")) for line in lines)
    return [(1, float(match.group(1)))] if match and has_atom else None


def _autodock_model_sha256s(path: Path) -> dict[int, str]:
    """Hashes of the exact wrapper-free model bytes optimized upstream."""
    try:
        text = Path(path).read_text(errors="strict")
    except (OSError, UnicodeError):
        return {}
    lines = text.splitlines(keepends=True)
    if not any(line.startswith(("MODEL", "ENDMDL")) for line in lines):
        return {1: hashlib.sha256(text.encode()).hexdigest()}
    output: dict[int, str] = {}
    rank: int | None = None
    body: list[str] | None = None
    for line in lines:
        if line.startswith("MODEL"):
            fields = line.split()
            try:
                rank = int(fields[1])
            except (IndexError, ValueError):
                return {}
            body = []
        elif line.startswith("ENDMDL"):
            if rank is None or body is None:
                return {}
            output[rank] = hashlib.sha256("".join(body).encode()).hexdigest()
            rank, body = None, None
        elif body is not None:
            body.append(line)
    return {} if body is not None else output


def _valid_sdf_pose(path: Path) -> bool:
    """Cheap fail-closed structural check for one optimizer-produced SDF."""
    try:
        text = Path(path).read_text(errors="ignore")
    except OSError:
        return False
    if not text.strip().endswith("$$$$") or text.count("$$$$") != 1:
        return False
    record = text.split("$$$$", 1)[0]
    lines = record.splitlines()
    if len(lines) < 4:
        return False
    counts = lines[3]
    if "V3000" in counts:
        return ("M  V30 BEGIN ATOM" in record and "M  V30 END ATOM" in record
                and "M  V30 END CTAB" in record)
    try:
        return int(counts[:3]) > 0 and "M  END" in record
    except ValueError:
        return False


def _positive_int(value: Any) -> int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if np.isfinite(number) and number > 0 and number.is_integer() else None


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if np.isfinite(number) else None


def _find_autodock_optimization_log(pdbqt_file: Path, poses_dir: Path) -> Path | None:
    """Nearest optimization log inside the configured AutoDock tree."""
    root = Path(poses_dir).resolve()
    current = pdbqt_file.resolve().parent
    while current == root or current.is_relative_to(root):
        candidate = current / "optimization_log.csv"
        if candidate.is_file():
            return candidate
        if current == root:
            break
        current = current.parent
    return None


def _resolve_optimized_sdf(raw_value: Any, log_path: Path, pdbqt_file: Path,
                           tool: str) -> Path | None:
    value = str(raw_value or "").strip()
    if not value or value.lower() == "nan":
        return None
    given = Path(value).expanduser()
    candidates = ([given] if given.is_absolute() else
                  [given, log_path.parent / given, pdbqt_file.parent / given])
    expected = (pdbqt_file.parent / f"optimized_{tool}").resolve()
    for candidate in candidates:
        path = candidate.resolve()
        if path.is_file() and path.is_relative_to(expected) and _valid_sdf_pose(path):
            return path
    return None


def _resolve_logged_file(raw_value: Any, log_path: Path) -> Path | None:
    value = str(raw_value or "").strip()
    if not value or value.lower() == "nan":
        return None
    given = Path(value).expanduser()
    candidates = [given] if given.is_absolute() else [given, log_path.parent / given]
    for candidate in candidates:
        path = candidate.resolve()
        if path.is_file():
            return path
    return None


def _validated_optimizer_provenance(
    record: dict,
    log_path: Path,
    pdbqt_file: Path,
    optimized_file: Path,
    tool: str,
    expected_model_digest: str,
) -> tuple[Path, dict] | None:
    """Validate the mandatory v2 optimizer sidecar against current artifacts."""
    sidecar = _resolve_logged_file(record.get("provenance_file"), log_path)
    expected_sidecar = Path(f"{optimized_file}.provenance.json").resolve()
    if sidecar is None or sidecar != expected_sidecar:
        return None
    try:
        provenance = json.loads(sidecar.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(provenance, dict) or provenance.get("schema_version") != 2:
        return None
    fingerprint = str(record.get("provenance_fingerprint") or "").strip().lower()
    if not fingerprint or fingerprint == "nan" or provenance.get("fingerprint") != fingerprint:
        return None

    # Recompute the fingerprint exactly as run_autodock does: it hashes the
    # immutable material before adding dynamic paths/timestamps/output digest.
    material = {
        key: value for key, value in provenance.items()
        if key not in {
            "fingerprint", "source_pose", "receptor", "created_at", "output_sha256",
            "optimizer_elapsed_time_s",
        }
    }
    encoded = json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(encoded).hexdigest() != fingerprint:
        return None
    if provenance.get("tool") != tool:
        return None
    if any(key not in provenance for key in (
        "settings", "optimizer", "converter", "input_identity", "source_scoring",
        "source_pose_sha256", "source_pose_container_sha256", "receptor_sha256",
        "output_sha256", "optimizer_elapsed_time_s",
    )):
        return None

    source_model_digest = str(record.get("source_pose_sha256") or "").strip().lower()
    receptor_digest = str(record.get("receptor_sha256") or "").strip().lower()
    if (not source_model_digest or source_model_digest == "nan"
            or source_model_digest != expected_model_digest
            or provenance.get("source_pose_sha256") != source_model_digest):
        return None
    if (not receptor_digest or receptor_digest == "nan"
            or provenance.get("receptor_sha256") != receptor_digest):
        return None
    optimizer_elapsed = _finite_float(provenance.get("optimizer_elapsed_time_s"))
    logged_elapsed = _finite_float(record.get("elapsed_time_s"))
    if optimizer_elapsed is None or optimizer_elapsed < 0 or logged_elapsed is None or logged_elapsed < 0:
        return None
    try:
        source_container = Path(str(provenance["source_pose"])).expanduser().resolve()
        receptor = Path(str(provenance["receptor"])).expanduser().resolve()
    except (KeyError, TypeError):
        return None
    if source_container != pdbqt_file.resolve() or not receptor.is_file():
        return None
    try:
        if provenance["source_pose_container_sha256"] != _sha256_file(pdbqt_file):
            return None
        if provenance["receptor_sha256"] != _sha256_file(receptor):
            return None
        if provenance["output_sha256"] != _sha256_file(optimized_file):
            return None
    except OSError:
        return None
    return sidecar, provenance


def _autodock_optimized_rows(
    pdbqt_file: Path,
    poses_dir: Path,
    protein: str,
    ligand: str,
    current_models: list[tuple[int, float]],
    keep_variants: set[str] | None,
) -> list[dict]:
    """Join current Vina models to successful, loadable optimizer log rows.

    The log is the commit record. Merely finding a non-empty file under an
    ``optimized_*`` directory is insufficient because it may be stale, partial,
    or left by a failed/shorter rerun.
    """
    log_path = _find_autodock_optimization_log(pdbqt_file, poses_dir)
    if log_path is None:
        return []
    try:
        log = pd.read_csv(log_path, low_memory=False)
    except Exception:
        return []
    required = {"tool", "autodock_rank", "optimized_rank", "optimized_file", "status"}
    if not required.issubset(log.columns):
        return []
    affinities = dict(current_models)
    model_digests = _autodock_model_sha256s(pdbqt_file)
    if set(model_digests) != set(affinities):
        return []
    committed: dict[tuple[str, int], dict] = {}
    for record in log.to_dict("records"):
        tool = str(record.get("tool") or "").strip().lower()
        if tool not in {"smina", "gnina"} or not _autodock_variant_ok(tool, keep_variants):
            continue
        if str(record.get("status") or "").strip().lower() != "success":
            continue
        source_name = str(record.get("pose_file") or "").strip()
        if source_name and source_name.lower() != "nan" and Path(source_name).name != pdbqt_file.name:
            continue
        autodock_rank = _positive_int(record.get("autodock_rank"))
        optimized_rank = _positive_int(record.get("optimized_rank"))
        if autodock_rank not in affinities or optimized_rank is None:
            continue
        rank_metric = str(record.get("rank_metric") or "").strip().lower()
        if rank_metric not in {"minimized_affinity", "cnn_affinity", "cnn_score"}:
            continue
        if _finite_float(record.get(rank_metric)) is None:
            continue
        optimized_file = _resolve_optimized_sdf(
            record.get("optimized_file"), log_path, pdbqt_file, tool)
        if optimized_file is None:
            continue
        validated = _validated_optimizer_provenance(
            record, log_path, pdbqt_file, optimized_file, tool,
            model_digests[autodock_rank])
        if validated is None:
            continue
        provenance_file, provenance = validated
        row = {
            "docking_tool": "autodock",
            "protein": protein,
            "ligand": ligand,
            "file_path": str(optimized_file),
            "source_pdbqt": str(pdbqt_file.resolve()),
            "pose_count": 1,
            "optimizer": tool,
            "autodock_rank": autodock_rank,
            "optimized_rank": optimized_rank,
            "autodock_affinity": (
                _finite_float(record.get("vina_affinity"))
                if _finite_float(record.get("vina_affinity")) is not None
                else affinities[autodock_rank]
            ),
            "rank_metric": rank_metric,
            "optimization_log_file": str(log_path.resolve()),
            "optimizer_provenance_file": str(provenance_file),
            "optimizer_provenance_fingerprint": provenance["fingerprint"],
            "source_pose_model_sha256": provenance["source_pose_sha256"],
        }
        for source, target in (
            ("minimized_affinity", "minimized_affinity"),
            ("cnn_score", "cnn_score"),
            ("cnn_affinity", "cnn_affinity"),
            ("elapsed_time_s", "optimizer_elapsed_time_s"),
            ("processing_elapsed_time_s", "optimizer_processing_elapsed_time_s"),
        ):
            value = _finite_float(record.get(source))
            if value is not None:
                row[target] = value
        committed[(tool, autodock_rank)] = row
    output: list[dict] = []
    for tool in ("smina", "gnina"):
        tool_rows = [row for (row_tool, _), row in committed.items() if row_tool == tool]
        optimized_ranks = [row["optimized_rank"] for row in tool_rows]
        rank_metrics = {row["rank_metric"] for row in tool_rows}
        if (len(optimized_ranks) != len(set(optimized_ranks))
                or sorted(optimized_ranks) != list(range(1, len(optimized_ranks) + 1))
                or len(rank_metrics) > 1):
            # An internally inconsistent rank table cannot safely define a
            # variant; fail the whole tool closed instead of choosing rows.
            continue
        output.extend(sorted(tool_rows, key=lambda row: row["autodock_rank"]))
    return output


def _collect_autodock_rows(poses_dir: Path, keep_variants: set[str] | None = None) -> list[dict]:
    rows = []
    # Accept both the historical flat staging directory and the actual nested
    # Vina result tree. Only final docking outputs match this suffix.
    for pdbqt_file in sorted(Path(poses_dir).glob("**/*_vina_out.pdbqt")):
        stem = pdbqt_file.stem.replace("_vina_out", "")
        parts = stem.split("__")
        if len(parts) == 2:
            protein, ligand = parts
            models = _closed_scored_pdbqt_models(pdbqt_file, allow_unwrapped_single=True)
            if not models:
                continue
            protein = _normalize_protein(protein)
            ligand = _normalize_ligand(ligand)
            if _autodock_variant_ok("original", keep_variants):
                rows.append({
                    "docking_tool": "autodock",
                    "protein": protein,
                    "ligand": ligand,
                    "file_path": str(pdbqt_file),
                    "pose_count": len(models),
                    "optimizer": "original",
                })
            rows.extend(_autodock_optimized_rows(
                pdbqt_file, Path(poses_dir), protein, ligand, models, keep_variants))
    return rows


def _validated_unidock_commit(
    pdbqt_file: Path,
    summary_path: Path,
    done_marker: Path,
    manifest_path: Path,
    models: list[tuple[int, float]],
) -> tuple[dict, dict, dict] | None:
    """Validate the v2 manifest + JSON sentinel + current input/output hashes."""
    try:
        summary = json.loads(summary_path.read_text())
        marker = json.loads(done_marker.read_text())
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not all(isinstance(value, dict) for value in (summary, marker, manifest)):
        return None
    if ([rank for rank, _ in models] != list(range(1, len(models) + 1))
            or [affinity for _, affinity in models]
            != sorted(affinity for _, affinity in models)):
        return None
    if (marker.get("schema_version") != 2 or marker.get("engine") != "unidock-tiled"
            or marker.get("status") != "complete"
            or marker.get("output_file") != pdbqt_file.name):
        return None
    if manifest.get("schema_version") != 2 or manifest.get("status") != "complete":
        return None
    if str(summary.get("status") or "").strip().lower() not in {"success", "done"}:
        return None
    for record in (summary, marker, manifest):
        if record.get("commit_status") != "complete":
            return None
    for key in (
        "fingerprint", "generation_id", "output_file", "output_sha256",
        "num_poses", "n_subboxes",
    ):
        if marker.get(key) in (None, ""):
            return None
        if not (marker.get(key) == manifest.get(key) == summary.get(key)):
            return None
    if (marker.get("num_poses") != len(models)
            or _positive_int(marker.get("n_subboxes")) is None):
        return None
    try:
        if _sha256_file(pdbqt_file) != marker["output_sha256"]:
            return None
        provenance = manifest["provenance"]
        if (provenance.get("schema_version") != 2
                or not isinstance(provenance.get("engine"), dict)
                or marker.get("n_subboxes") != len(provenance.get("subboxes") or [])):
            return None
        canonical = json.dumps(
            provenance, sort_keys=True, separators=(",", ":"), allow_nan=False)
        if hashlib.sha256(canonical.encode()).hexdigest() != marker["fingerprint"]:
            return None
        inputs = provenance["inputs"]
        for name in ("receptor", "ligand", "box"):
            item = inputs[name]
            path = Path(str(item["path"])).expanduser().resolve()
            if (not path.is_file() or int(item["size"]) != path.stat().st_size
                    or item["sha256"] != _sha256_file(path)):
                return None
    except (KeyError, TypeError, ValueError, OSError):
        return None
    return summary, marker, manifest


def _collect_unidock_rows(poses_dir: Path, keep_variants: set[str] | None = None) -> list[dict]:
    """Collect only transactionally complete tiled Uni-Dock merged outputs."""
    if not _variant_ok("tiled", keep_variants):
        return []
    rows: list[dict] = []
    for pdbqt_file in sorted(Path(poses_dir).glob("**/*_unidock_out.pdbqt")):
        out_dir = pdbqt_file.parent
        done_marker = out_dir / ".unidock_done"
        summary_path = out_dir / "docking_summary.json"
        manifest_path = out_dir / "run_manifest.json"
        if not done_marker.is_file() or not summary_path.is_file() or not manifest_path.is_file():
            continue
        models = _closed_scored_pdbqt_models(pdbqt_file)
        if not models:
            continue
        committed = _validated_unidock_commit(
            pdbqt_file, summary_path, done_marker, manifest_path, models)
        if committed is None:
            continue
        summary, marker, manifest = committed
        combo = pdbqt_file.stem[: -len("_unidock_out")]
        if "__" in combo:
            protein, ligand = combo.split("__", 1)
        else:
            # Benchmark tiling uses one <PDB>_<CCD> id for both sides, matching
            # the normalized AutoDock/DiffDock/EquiBind intersection key.
            protein = ligand = combo
        row = {
            "docking_tool": "unidock",
            "protein": _normalize_protein(protein),
            "ligand": _normalize_ligand(ligand),
            "file_path": str(pdbqt_file),
            "pose_count": len(models),
            "variant": "tiled",
            "scoring": str(manifest["provenance"].get("settings", {}).get("scoring") or "vina"),
            "engine": str(marker["engine"]),
            "docking_summary_file": str(summary_path.resolve()),
            "unidock_manifest_file": str(manifest_path.resolve()),
            "unidock_fingerprint": marker["fingerprint"],
            "unidock_generation_id": marker["generation_id"],
        }
        receptor = str(
            manifest.get("provenance", {}).get("inputs", {}).get("receptor", {}).get("path")
            or ""
        ).strip()
        if receptor:
            receptor_path = Path(receptor).expanduser()
            if not receptor_path.is_absolute():
                receptor_path = (summary_path.parent / receptor_path).resolve()
            row["docking_receptor_file"] = str(receptor_path)
        rows.append(row)
    return rows


_UNIDOCK2_AFFINITY_HEADER = re.compile(
    r"^\s*>\s*<vina_binding_free_energy>", re.IGNORECASE)


def _strict_unidock2_sdf_scores(sdf_path: Path) -> list[float] | None:
    """Validate every Uni-Dock2 SDF record and return its ordered affinities.

    This mirrors the transactional driver's fail-closed output contract: every
    record is ``$$$$``-terminated, has one non-empty V2000/V3000 mol block, and
    carries exactly one finite ``vina_binding_free_energy`` value. Merely counting
    terminators is insufficient because a complete-looking corrupt/unscored SDF
    must never become a benchmark pose.
    """
    try:
        lines = Path(sdf_path).read_text(errors="strict").splitlines()
    except (OSError, UnicodeError):
        return None
    records: list[list[str]] = []
    current: list[str] = []
    for line in lines:
        if line.strip() == "$$$$":
            if not any(part.strip() for part in current):
                return None
            records.append(current)
            current = []
        else:
            current.append(line)
    if any(part.strip() for part in current) or not records:
        return None

    scores: list[float] = []
    try:
        for record in records:
            if len(record) < 5:
                return None
            mend = [i for i, line in enumerate(record) if line.strip() == "M  END"]
            if len(mend) != 1:
                return None
            counts = record[3]
            if "V2000" in counts:
                atom_count = int(counts[:3])
                bond_count = int(counts[3:6])
                if atom_count < 1 or bond_count < 0 or mend[0] < 4 + atom_count + bond_count:
                    return None
                for atom_line in record[4:4 + atom_count]:
                    xyz = tuple(float(atom_line[start:start + 10]) for start in (0, 10, 20))
                    if not all(np.isfinite(value) for value in xyz) or not atom_line[31:34].strip():
                        return None
                for bond_line in record[4 + atom_count:4 + atom_count + bond_count]:
                    first, second, order = (int(bond_line[start:start + 3])
                                            for start in (0, 3, 6))
                    if not (1 <= first <= atom_count and 1 <= second <= atom_count and order >= 1):
                        return None
            elif "V3000" in counts:
                count_line = next((line for line in record if "M  V30 COUNTS" in line), None)
                if count_line is None:
                    return None
                fields = count_line.split("COUNTS", 1)[1].split()
                if len(fields) < 2:
                    return None
                atom_count, bond_count = int(fields[0]), int(fields[1])
                begin_atom = next((i for i, line in enumerate(record)
                                   if "M  V30 BEGIN ATOM" in line), None)
                end_atom = next((i for i, line in enumerate(record)
                                 if "M  V30 END ATOM" in line), None)
                if atom_count < 1 or bond_count < 0 or begin_atom is None or end_atom is None:
                    return None
                atom_lines = [line for line in record[begin_atom + 1:end_atom]
                              if "M  V30" in line]
                if len(atom_lines) != atom_count:
                    return None
                for atom_line in atom_lines:
                    fields = atom_line.split()
                    if len(fields) < 7:
                        return None
                    xyz = tuple(float(value) for value in fields[4:7])
                    if not all(np.isfinite(value) for value in xyz):
                        return None
                if bond_count:
                    begin_bond = next((i for i, line in enumerate(record)
                                       if "M  V30 BEGIN BOND" in line), None)
                    end_bond = next((i for i, line in enumerate(record)
                                     if "M  V30 END BOND" in line), None)
                    if begin_bond is None or end_bond is None or end_bond <= begin_bond:
                        return None
                    bond_lines = [line for line in record[begin_bond + 1:end_bond]
                                  if "M  V30" in line]
                    if len(bond_lines) != bond_count:
                        return None
            else:
                return None

            headers = [i for i, line in enumerate(record)
                       if _UNIDOCK2_AFFINITY_HEADER.match(line)]
            if len(headers) != 1 or headers[0] + 1 >= len(record):
                return None
            value_text = record[headers[0] + 1].strip()
            if not value_text:
                return None
            affinity = float(value_text)
            if not np.isfinite(affinity):
                return None
            scores.append(affinity)
    except (IndexError, TypeError, ValueError):
        return None

    tolerance = 0.002
    if any(later < earlier - tolerance for earlier, later in zip(scores, scores[1:])):
        return None
    return scores


def _pdb_atom_counts(path: Path) -> tuple[int, int] | None:
    """Strict ``(ATOM records, heavy ATOM records)`` for an attested PDB."""
    atoms = heavy = 0
    try:
        lines = Path(path).read_text(errors="strict").splitlines()
    except (OSError, UnicodeError):
        return None
    for line in lines:
        if not line.startswith("ATOM"):
            continue
        atoms += 1
        try:
            xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
        except ValueError:
            return None
        if not all(np.isfinite(value) for value in xyz):
            return None
        atom_name = re.sub(r"^[0-9]+", "", line[12:16].strip()).upper()
        element = (line[76:78].strip() if len(line) >= 78 else "") or atom_name[:1]
        if element.upper() not in {"H", "D", "T"}:
            heavy += 1
    return (atoms, heavy) if atoms > 0 and heavy > 0 else None


_UNIDOCK2_SENTINEL_SCHEMA = 3
_UNIDOCK2_PROVENANCE_KEYS = (
    "schema_version", "engine", "tool", "inputs", "center", "size",
    "effective_config", "effective_config_sha256", "driver_execution",
)
_UNIDOCK2_ADVANCED_PROTOCOL = {
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
}
_UNIDOCK2_HARDWARE_PROTOCOL = {"gpu_device_id": 0, "n_cpu": 1}
_UNIDOCK2_PREPROCESSING_PROTOCOL = {
    "construct_ff": False,
    "template_docking": False,
    "compute_center": True,
    "covalent_ligand": False,
    "preserve_receptor_hydrogen": False,
    "engine_checkpoint": False,
}
_UNIDOCK2_GPU_LOCK_FILE = "/tmp/master_docking_gpu0.lock"


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _strict_typed_mapping(value: object, expected: dict) -> bool:
    """Exact mapping equality without Python's ``True == 1`` coercion."""
    return (
        isinstance(value, dict)
        and set(value) == set(expected)
        and all(type(value[key]) is type(wanted) and value[key] == wanted
                for key, wanted in expected.items())
    )


def _strict_float_vector(value: object, *, positive: bool = False) -> bool:
    return (
        isinstance(value, list)
        and len(value) == 3
        and all(type(item) is float and np.isfinite(item)
                and (not positive or item > 0) for item in value)
    )


def _unidock2_benchmark_protocol_matches(completion: dict) -> bool:
    """Fail closed unless the manifest is for this exact benchmark protocol.

    Schema 3 makes every Uni-Dock2 0.6.3 parser default explicit. PoseBusters is
    wired to the one benchmark output tree, so accepting a self-consistent but
    different protocol would silently mix experiments under the same method.
    Only the per-complex center/box values may vary.
    """
    try:
        effective = completion["effective_config"]
        if not isinstance(effective, dict) or set(effective) != {
            "Advanced", "Hardware", "Settings", "Preprocessing",
        }:
            return False
        if (not _strict_typed_mapping(
                effective["Advanced"], _UNIDOCK2_ADVANCED_PROTOCOL)
                or not _strict_typed_mapping(
                    effective["Hardware"], _UNIDOCK2_HARDWARE_PROTOCOL)
                or not _strict_typed_mapping(
                    effective["Preprocessing"], _UNIDOCK2_PREPROCESSING_PROTOCOL)):
            return False
        settings = effective["Settings"]
        if (not isinstance(settings, dict)
                or set(settings) != {"box_size", "task", "search_mode"}
                or type(settings["task"]) is not str or settings["task"] != "screen"
                or type(settings["search_mode"]) is not str
                or settings["search_mode"] != "free"):
            return False
        center = completion["center"]
        size = completion["size"]
        box_size = settings["box_size"]
        if (not _strict_float_vector(center)
                or not _strict_float_vector(size, positive=True)
                or not _strict_float_vector(box_size, positive=True)
                or size != box_size):
            return False
        driver_execution = completion["driver_execution"]
        return (
            isinstance(driver_execution, dict)
            and set(driver_execution) == {"gpu_lock_file"}
            and type(driver_execution["gpu_lock_file"]) is str
            and driver_execution["gpu_lock_file"] == _UNIDOCK2_GPU_LOCK_FILE
        )
    except (KeyError, TypeError):
        return False


def _unidock2_manifest_fingerprint_matches(completion: dict) -> bool:
    """Recompute both nested-config and whole-provenance schema-3 digests."""
    try:
        effective_json = json.dumps(
            completion["effective_config"],
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        if completion.get("effective_config_sha256") != hashlib.sha256(
            effective_json.encode("utf-8")
        ).hexdigest():
            return False
        provenance = {key: completion[key] for key in _UNIDOCK2_PROVENANCE_KEYS}
        canonical = json.dumps(
            provenance, sort_keys=True, separators=(",", ":"), allow_nan=False
        )
        return (
            _is_sha256(completion.get("fingerprint"))
            and completion["fingerprint"]
            == hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        )
    except (KeyError, TypeError, ValueError):
        return False


def _unidock2_runtime_fields_match(completion: dict) -> bool:
    """Check the exact stable runtime-attestation contract used by v3."""
    try:
        from Scripts.Docking import run_unidock2 as unidock2
        provenance = {key: completion[key] for key in _UNIDOCK2_PROVENANCE_KEYS}
        requested = unidock2._requested_runtime_settings(provenance)
        runtime = completion["runtime_attestation"]
        return (
            unidock2._json_equal(completion.get("requested_search_settings"), requested)
            and unidock2._runtime_attestation_matches(runtime, requested)
            and unidock2._json_equal(completion.get("effective_search_settings"), runtime)
            and unidock2._json_equal(completion.get("applied_center"), runtime.get("center"))
            and unidock2._json_equal(completion.get("applied_box_size"), runtime.get("box_size"))
        )
    except (ImportError, KeyError, TypeError, ValueError):
        return False


@lru_cache(maxsize=8)
def _live_unidock2_tool_identity(binary_path: str) -> dict:
    """Driver-owned identity probe, cached by canonical launcher path."""
    from Scripts.Docking import run_unidock2 as unidock2
    return unidock2.probe_tool_identity(binary_path)


def _current_unidock2_tool_identity(tool: object) -> bool:
    """Validate every file/source hash in the schema-3 Uni-Dock2 identity."""
    expected_keys = {
        "path", "version", "launcher_sha256", "tleap_path", "ambpdb_path",
        "artifacts", "python_source_tree",
    }
    if not isinstance(tool, dict) or set(tool) != expected_keys:
        return False
    if tool.get("version") != "0.6.3":
        return False
    try:
        artifacts = tool["artifacts"]
        if not isinstance(artifacts, list) or not artifacts:
            return False
        artifact_paths: set[str] = set()
        for item in artifacts:
            if not isinstance(item, dict) or set(item) != {"path", "size", "sha256"}:
                return False
            raw = item["path"]
            path = Path(raw).expanduser()
            if (not isinstance(raw, str) or not path.is_absolute()
                    or str(path.resolve()) != raw or not path.is_file()
                    or isinstance(item["size"], bool) or not isinstance(item["size"], int)
                    or item["size"] != path.stat().st_size
                    or not _is_sha256(item["sha256"])
                    or item["sha256"] != _sha256_file(path)
                    or raw in artifact_paths):
                return False
            artifact_paths.add(raw)

        for key in ("path", "tleap_path", "ambpdb_path"):
            raw = tool[key]
            path = Path(raw).expanduser()
            if (not isinstance(raw, str) or not path.is_absolute()
                    or str(path.resolve()) != raw or raw not in artifact_paths):
                return False
        if (not _is_sha256(tool["launcher_sha256"])
                or tool["launcher_sha256"] != _sha256_file(tool["path"])):
            return False

        source_tree = tool["python_source_tree"]
        if (not isinstance(source_tree, dict)
                or set(source_tree) != {"path", "file_count", "size", "sha256"}):
            return False
        source_root = Path(source_tree["path"]).expanduser()
        if (not source_root.is_absolute()
                or str(source_root.resolve()) != source_tree["path"]
                or not source_root.is_dir()):
            return False
        from Scripts.Docking import run_unidock2 as unidock2
        current_tree = unidock2._python_tree_provenance(source_root)
        if not unidock2._json_equal(source_tree, current_tree):
            return False
        # Schema shape alone is insufficient: a re-signed manifest could omit
        # the engine library, teLeap, or one of the conda records while retaining
        # the three named launcher paths above. The driver's live probe is the
        # single source of truth for the complete ordered artifact inventory.
        live_identity = _live_unidock2_tool_identity(str(Path(tool["path"]).resolve()))
        return unidock2._json_equal(tool, live_identity)
    except (ImportError, KeyError, OSError, RuntimeError, TypeError, ValueError):
        return False


def _validated_unidock2_commit(
    sdf_file: Path,
    completion_path: Path,
    scores: list[float],
) -> dict | None:
    """Validate a schema-3 Uni-Dock2 completion against all current artifacts.

    The completion file is the commit record. ``docking_summary.json`` is
    intentionally absent from this contract: it is a replaceable progress/index
    view and may be missing or stale without invalidating an otherwise complete
    transaction.
    """
    try:
        completion = json.loads(completion_path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(completion, dict):
        return None
    suffix = "_unidock2_out"
    if not sdf_file.stem.endswith(suffix):
        return None
    combo = sdf_file.stem[: -len(suffix)]
    if completion_path != sdf_file.parent / f"{combo}_unidock2_completion.json":
        return None
    if (completion.get("schema_version") != _UNIDOCK2_SENTINEL_SCHEMA
            or completion.get("engine") != "unidock2"
            or completion.get("status") != "success"
            or completion.get("output_file") != sdf_file.name):
        return None
    fingerprint = completion.get("fingerprint")
    generation = completion.get("generation_id")
    output_sha = completion.get("output_sha256")
    if (not _is_sha256(fingerprint) or not _is_sha256(output_sha)
            or not isinstance(generation, str)
            or re.fullmatch(r"[0-9a-f]{32}", generation) is None):
        return None
    if (_positive_int(completion.get("num_poses")) != len(scores)
            or output_sha != _sha256_file(sdf_file)):
        return None
    best = _finite_float(completion.get("best_affinity"))
    if best is None or abs(best - min(scores)) > 0.002:
        return None

    try:
        provenance = {key: completion[key] for key in _UNIDOCK2_PROVENANCE_KEYS}
        effective = completion["effective_config"]
        if not _unidock2_manifest_fingerprint_matches(completion):
            return None
        if not _unidock2_benchmark_protocol_matches(completion):
            return None
        if not _current_unidock2_tool_identity(completion["tool"]):
            return None
        driver_execution = completion["driver_execution"]
        if not isinstance(driver_execution, dict) or set(driver_execution) != {"gpu_lock_file"}:
            return None
        gpu_lock_file = driver_execution.get("gpu_lock_file")
        if (not isinstance(gpu_lock_file, str) or not gpu_lock_file
                or not Path(gpu_lock_file).is_absolute()
                or str(Path(gpu_lock_file).expanduser().resolve()) != gpu_lock_file):
            return None
        inputs = completion["inputs"]
        if not isinstance(inputs, dict) or set(inputs) != {"receptor", "ligand", "box"}:
            return None
        expected_input_names = {
            "receptor": f"{combo}_protein.pdb",
            "ligand": f"{combo}_ligand_start_conf.sdf",
            "box": f"{combo}_protein_mgl_tools.box.txt",
        }
        input_paths: dict[str, Path] = {}
        for name in ("receptor", "ligand", "box"):
            item = inputs[name]
            if (not isinstance(item, dict)
                    or set(item) != {"path", "size", "sha256"}
                    or not isinstance(item.get("path"), str)
                    or not _is_sha256(item.get("sha256"))):
                return None
            raw_path = item["path"]
            current = Path(raw_path).expanduser()
            if (not current.is_absolute() or str(current.resolve()) != raw_path
                    or not current.is_file()
                    or current.name != expected_input_names[name]
                    or isinstance(item.get("size"), bool)
                    or not isinstance(item.get("size"), int)
                    or item["size"] != current.stat().st_size
                    or item["sha256"] != _sha256_file(current)):
                return None
            input_paths[name] = current
        from Scripts.Docking import run_unidock2 as unidock2
        unidock2.validate_single_input_ligand(input_paths["ligand"])
        parsed_center, parsed_size = unidock2.read_box(input_paths["box"])
        if (not unidock2._json_equal(list(parsed_center), completion["center"])
                or not unidock2._json_equal(list(parsed_size), completion["size"])):
            return None
        advanced = effective["Advanced"]
        requested_num_pose = _positive_int(advanced["num_pose"])
        energy_range = _finite_float(advanced["energy_range"])
        if (requested_num_pose is None or energy_range is None or energy_range < 0):
            return None
        if not unidock2._completion_matches(
            sdf_file, completion_path, fingerprint, provenance,
            {"num_pose": requested_num_pose, "energy_range": energy_range},
        ):
            return None
        if (not _is_sha256(completion.get("receptor_prmtop_sha256"))
                or not _is_sha256(completion.get("receptor_inpcrd_sha256"))
                or _positive_int(completion.get("input_protein_heavy_atoms")) is None
                or (_finite_float(completion.get("completed_unix_s")) or 0) <= 0):
            return None
    except (ImportError, KeyError, TypeError, ValueError, OSError):
        return None
    return completion


def _collect_unidock2_rows(poses_dir: Path, keep_variants: set[str] | None = None) -> list[dict]:
    """Collect only schema-3 committed Uni-Dock2 outputs.

    ``docking_summary.json`` is never sufficient on its own. The authoritative
    ``<id>_unidock2_completion.json`` must match the current scored SDF, all
    fingerprinted inputs, and the materialized engine receptor byte-for-byte.
    """
    if not _variant_ok("single", keep_variants):
        return []
    rows: list[dict] = []
    for sdf_file in sorted(Path(poses_dir).glob("**/*_unidock2_out.sdf")):
        out_dir = sdf_file.parent
        summary_path = out_dir / "docking_summary.json"
        combo = sdf_file.stem[: -len("_unidock2_out")]
        completion_path = out_dir / f"{combo}_unidock2_completion.json"
        if not completion_path.is_file():
            continue
        scores = _strict_unidock2_sdf_scores(sdf_file)
        if not scores:
            continue
        completion = _validated_unidock2_commit(sdf_file, completion_path, scores)
        if completion is None:
            continue
        if "__" in combo:
            protein, ligand = combo.split("__", 1)
        else:
            # Whole-protein docking uses one <PDB>_<CCD> id for both sides, the
            # same normalized intersection key the other collectors emit.
            protein = ligand = combo
        row = {
            "docking_tool": "unidock2",
            "protein": _normalize_protein(protein),
            "ligand": _normalize_ligand(ligand),
            "file_path": str(sdf_file),
            "pose_count": len(scores),
            "variant": "single",
            "scoring": "vina",
            "engine": "unidock2",
            "unidock2_completion_file": str(completion_path.resolve()),
            "unidock2_fingerprint": completion["fingerprint"],
            "unidock2_generation_id": completion["generation_id"],
            "unidock2_output_sha256": completion["output_sha256"],
        }
        prepared = sdf_file.parent / completion["prepared_receptor_file"]
        row["docking_receptor_file"] = str(prepared.resolve())
        row["unidock2_prepared_receptor_file"] = str(prepared.resolve())
        row["unidock2_prepared_receptor_sha256"] = completion[
            "prepared_receptor_sha256"]
        row["unidock2_prepared_receptor_atoms"] = completion[
            "prepared_receptor_atoms"]
        row["unidock2_prepared_receptor_heavy_atoms"] = completion[
            "prepared_receptor_heavy_atoms"]
        row["unidock2_engine_receptor_heavy_atoms_in_box"] = completion[
            "engine_receptor_heavy_atoms_in_box"]
        # The summary is only a human-readable progress/index view. Preserve its
        # path when it is parseable, but never let absence, staleness, or status
        # wording override the schema-3 commit above.
        try:
            summary = json.loads(summary_path.read_text())
        except (OSError, json.JSONDecodeError):
            summary = None
        if isinstance(summary, dict):
            row["docking_summary_file"] = str(summary_path.resolve())
        rows.append(row)
    return rows


def _collect_diffdock_rows(poses_dir: Path, keep_variants: set[str] | None = None) -> list[dict]:
    rows = []
    skip_names = {"prepared_proteins", "converted_pdbqt", "prepared_ligands", "Orai"}
    for subdir in Path(poses_dir).iterdir():
        # NB: do NOT skip on "_ligand__" — the benchmark staging labels every
        # pair <pdb>_ligand__<pdb>_protein (same convention the equibind
        # collector accepts), so that filter would drop every staged pair.
        # Real artefact dirs are handled by skip_names above.
        if not subdir.is_dir() or subdir.name in skip_names:
            continue
        parts = subdir.name.split("__")
        if len(parts) == 2:
            ligand, protein = parts
            sdfs = [f for f in subdir.glob("**/*.sdf")
                    if not _is_diffdock_bare_rank(f)]   # drop confidence-less duplicate
            if keep_variants:                       # restrict to allowed optimizer(s)
                sdfs = [f for f in sdfs
                        if _variant_ok(_diffdock_optimizer_of(f), keep_variants)]
            pose_count = len(sdfs)
            if pose_count > 0:
                rows.append({
                    "docking_tool": "diffdock",
                    "protein": _normalize_protein(protein),
                    "ligand": _normalize_ligand(ligand),
                    "file_path": str(subdir), "pose_count": pose_count,
                })
    return rows


def _collect_equibind_rows(poses_dir: Path, keep_variants: set[str] | None = None) -> list[dict]:
    rows = []
    for subdir in Path(poses_dir).iterdir():
        if not subdir.is_dir() or "_pymol" in subdir.name:
            continue
        dir_clean = (
            subdir.name[: -len("_spatial_sites")]
            if subdir.name.endswith("_spatial_sites")
            else subdir.name
        )
        parts = dir_clean.split("__")
        if len(parts) == 2:
            ligand, protein = parts
            # The finalized-results record is authoritative. A non-empty SDF in
            # an unfinished combo is not proof that post-processing completed.
            sdfs = _equibind_finalized_sdfs(subdir)
            if keep_variants:                       # restrict to allowed refine_variant(s)
                sdfs = [f for f in sdfs
                        if _variant_ok(_equibind_refine_of(f), keep_variants)]
            pose_count = len(sdfs)
            if pose_count > 0:
                rows.append({
                    "docking_tool": "equibind",
                    "protein": _normalize_protein(protein),
                    "ligand": _normalize_ligand(ligand),
                    "file_path": str(subdir), "pose_count": pose_count,
                })
    return rows


ROW_COLLECTORS = {
    "autodock": _collect_autodock_rows,
    "autodock_vinardo": _collect_autodock_rows,   # same Vina-output tree, Vinardo scoring
    "unidock": _collect_unidock_rows,
    "unidock_tiled": _collect_unidock_rows,
    "unidock2": _collect_unidock2_rows,
    "diffdock": _collect_diffdock_rows,
    "equibind_guided": _collect_equibind_rows,
    "equibind_exclusion": _collect_equibind_rows,
    "equibind_docked_poses": _collect_equibind_rows,
}


# ============================================================================
# POSE-FILE EXPANDERS  (row -> per-pose dicts)
# ============================================================================

def _autodock_affinities_from_pdbqt(pdbqt_path: Path) -> dict[int, float]:
    """{model_number: Vina affinity (kcal/mol)} from a Vina output PDBQT.

    Vina writes one ``REMARK VINA RESULT:`` line per ``MODEL`` (more-negative =
    better; models are already sorted best-first, so the model number IS the Vina
    rank). This is the AutoDock score that never made it into the results CSV —
    captured here so future runs carry it natively.
    """
    out: dict[int, float] = {}
    cur: int | None = None
    try:
        for ln in Path(pdbqt_path).read_text(errors="ignore").splitlines():
            if ln.startswith("MODEL"):
                parts = ln.split()
                cur = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else None
            elif ln.startswith("REMARK VINA RESULT:") and cur is not None:
                m = re.search(r"(-?\d+\.\d+)", ln)
                if m:
                    out[cur] = float(m.group(1))
    except (OSError, ValueError):
        pass
    return out


def _expand_autodock_poses(row: dict, conv_dir: Path, ctx: PipelineConfig) -> list[dict]:
    file_path = Path(row["file_path"])
    if not file_path.exists():
        return []
    if file_path.suffix.lower() == ".sdf":
        if not _valid_sdf_pose(file_path):
            return []
        source_pdbqt = Path(str(row.get("source_pdbqt") or ""))
        if not source_pdbqt.is_file():
            return []
        optimizer = str(row.get("optimizer") or "").strip().lower()
        if optimizer not in {"smina", "gnina"}:
            return []
        pose = {
            "method": row["docking_tool"],
            "protein": row["protein"],
            "ligand": row["ligand"],
            "pose_file": str(file_path),
            "source_pose_file": str(source_pdbqt.resolve()),
            "pose_name": (
                f"{row['protein']}__{row['ligand']}/optimized_{optimizer}/"
                f"{file_path.name}"
            ),
            "file_format": "sdf",
            "optimizer": optimizer,
        }
        for key in (
            "autodock_rank", "optimized_rank", "autodock_affinity", "rank_metric",
            "minimized_affinity", "cnn_score", "cnn_affinity",
            "optimizer_elapsed_time_s", "optimizer_processing_elapsed_time_s",
            "optimization_log_file",
            "optimizer_provenance_file", "optimizer_provenance_fingerprint",
            "source_pose_model_sha256",
        ):
            if key in row and row[key] not in (None, ""):
                pose[key] = row[key]
        return [pose]
    if file_path.suffix.lower() != ".pdbqt":
        return []
    if not _QUIET:
        print(f"  Converting {file_path.name} to SDF...")
    sdf_files = _convert_pdbqt_to_sdf(str(file_path), conv_dir, ctx)
    affinities = _autodock_affinities_from_pdbqt(file_path)
    poses = []
    for sdf_file in sdf_files:
        sdf_path = Path(sdf_file)
        model_num = sdf_path.stem.split("_model")[-1] if "_model" in sdf_path.stem else "1"
        pose = {
            "method": row["docking_tool"], "protein": row["protein"], "ligand": row["ligand"],
            "pose_file": sdf_file,
            "source_pose_file": str(file_path.resolve()),
            "pose_name": f"{row['protein']}__{row['ligand']}_pose{model_num}",
            "file_format": "sdf",
            "optimizer": str(row.get("optimizer") or "original"),
        }
        try:
            mn = int(model_num)
        except ValueError:
            mn = None
        if mn is not None:
            pose["autodock_rank"] = mn               # Vina rank (1 = best)
            if mn in affinities:
                pose["autodock_affinity"] = affinities[mn]   # kcal/mol
        poses.append(pose)
    return poses


def _expand_unidock_poses(row: dict, conv_dir: Path, ctx: PipelineConfig) -> list[dict]:
    """Convert the committed merged PDBQT and preserve its global score order."""
    file_path = Path(row["file_path"])
    models = _closed_scored_pdbqt_models(file_path)
    if not models:
        return []
    if not _QUIET:
        print(f"  Converting {file_path.name} to SDF...")
    sdf_files = _convert_pdbqt_to_sdf(str(file_path), conv_dir, ctx)
    affinities = dict(models)
    poses: list[dict] = []
    for position, sdf_file in enumerate(sdf_files, start=1):
        sdf_path = Path(sdf_file)
        model_text = sdf_path.stem.split("_model")[-1] if "_model" in sdf_path.stem else str(position)
        try:
            model_number = int(model_text)
        except ValueError:
            model_number = position
        if model_number not in affinities:
            continue
        pose = {
            "method": row["docking_tool"],
            "protein": row["protein"],
            "ligand": row["ligand"],
            "pose_file": str(sdf_path),
            "source_pose_file": str(file_path.resolve()),
            "pose_name": f"{row['protein']}__{row['ligand']}_pose{model_number}",
            "file_format": "sdf",
            "unidock_rank": model_number,
            "unidock_affinity": affinities[model_number],
            "variant": str(row.get("variant") or "tiled"),
            "scoring": str(row.get("scoring") or "vina"),
            "engine": str(row.get("engine") or "unidock-tiled"),
        }
        for key in (
            "docking_receptor_file", "docking_summary_file", "unidock_manifest_file",
            "unidock_fingerprint", "unidock_generation_id",
        ):
            if row.get(key):
                pose[key] = row[key]
        poses.append(pose)
    return poses


def _sdf_record_affinity(record_lines: list[str]) -> float | None:
    """Read Uni-Dock2's ``vina_binding_free_energy`` tag from one SDF record."""
    for i, line in enumerate(record_lines):
        if line.lstrip().startswith(">") and "vina_binding_free_energy" in line:
            for follow in record_lines[i + 1:]:
                value = follow.strip()
                if value:
                    return _finite_float(value)
            break
    return None


def _split_unidock2_records(sdf_path: Path, conv_dir: Path) -> list[tuple[Path, float | None]]:
    """Split a multi-record SDF into single-molecule ``<stem>_model<N>.sdf`` files.

    Records are written back byte-for-byte (no RDKit re-perception), so
    PoseBusters validates the exact docked geometry Uni-Dock2 produced. Returns
    ``(path, affinity)`` in the file's native best-first order; the 1-based index
    is therefore the Uni-Dock2 affinity rank, and the ``_model<N>`` naming is the
    token :func:`parse_rank` already understands in the analysis scripts."""
    scores = _strict_unidock2_sdf_scores(Path(sdf_path))
    if not scores:
        return []
    try:
        text = Path(sdf_path).read_text(errors="strict")
    except (OSError, UnicodeError):
        return []
    conv_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(sdf_path).stem
    results: list[tuple[Path, float | None]] = []
    record_lines: list[str] = []
    idx = 0
    for line in text.splitlines():
        if line.strip() == "$$$$":
            block = "\n".join(record_lines).strip("\n")
            if block:
                idx += 1
                out = conv_dir / f"{stem}_model{idx}.sdf"
                out.write_text(block + "\n$$$$\n")
                affinity = _sdf_record_affinity(record_lines)
                if affinity is None or abs(affinity - scores[idx - 1]) > 1e-12:
                    return []
                results.append((out, affinity))
            record_lines = []
        else:
            record_lines.append(line)
    return results if len(results) == len(scores) else []


def _expand_unidock2_poses(row: dict, conv_dir: Path, ctx: PipelineConfig) -> list[dict]:
    """Expand a Uni-Dock2 multi-pose SDF into one pose dict per record.

    The SDF is best-first affinity-ordered, so the 1-based record index is the
    Uni-Dock2 rank; the geometry is preserved verbatim (see
    :func:`_split_unidock2_records`)."""
    file_path = Path(row["file_path"])
    if not _QUIET:
        print(f"  Splitting {file_path.name} into per-pose SDFs...")
    records = _split_unidock2_records(file_path, conv_dir)
    if not records:
        return []
    poses: list[dict] = []
    for rank, (sdf_path, affinity) in enumerate(records, start=1):
        pose = {
            "method": row["docking_tool"],
            "protein": row["protein"],
            "ligand": row["ligand"],
            "pose_file": str(sdf_path),
            "source_pose_file": str(file_path.resolve()),
            "pose_name": f"{row['protein']}__{row['ligand']}_pose{rank}",
            "file_format": "sdf",
            "unidock2_rank": rank,
            "variant": str(row.get("variant") or "single"),
            "scoring": str(row.get("scoring") or "vina"),
            "engine": str(row.get("engine") or "unidock2"),
        }
        if affinity is not None:
            pose["unidock2_affinity"] = affinity
        for key in (
            "docking_receptor_file", "docking_summary_file",
            "unidock2_completion_file", "unidock2_fingerprint",
            "unidock2_generation_id", "unidock2_output_sha256",
            "unidock2_prepared_receptor_file", "unidock2_prepared_receptor_sha256",
            "unidock2_prepared_receptor_atoms",
            "unidock2_prepared_receptor_heavy_atoms",
            "unidock2_engine_receptor_heavy_atoms_in_box",
        ):
            if row.get(key):
                pose[key] = row[key]
        poses.append(pose)
    return poses


def _expand_diffdock_poses(row: dict, _conv_dir: Path, _ctx: PipelineConfig) -> list[dict]:
    file_path = Path(row["file_path"])
    if not file_path.is_dir():
        return []
    poses = []
    for sdf_file in sorted(file_path.glob("**/*.sdf")):
        if _is_diffdock_bare_rank(sdf_file):
            continue                                 # confidence-less duplicate of rankN_confidence-*.sdf
        confidence = None
        if "confidence" in sdf_file.stem.lower():
            m = re.search(r"confidence[_-]?([\d.]+)", sdf_file.stem, re.IGNORECASE)
            if m:
                try:
                    confidence = float(m.group(1))
                except ValueError:
                    pass
        try:
            rel = str(sdf_file.relative_to(file_path))
        except ValueError:
            rel = sdf_file.name
        # Optimizer provenance: smina/gnina re-optimised poses are written to
        # optimized_<tool>/ subfolders next to the original rank*.sdf (see
        # run_diffdock.optimize_results). Tag each pose with its optimizer so the
        # downstream reports keep optimised vs original DiffDock poses separable
        # (mirrors EquiBind's refine_variant) instead of silently pooling them.
        optimizer = "original"
        for part in sdf_file.parts:
            if part.startswith("optimized_"):
                optimizer = part[len("optimized_"):] or "original"
                break
        info = {
            "method": row["docking_tool"], "protein": row["protein"], "ligand": row["ligand"],
            "pose_file": str(sdf_file),
            "source_pose_file": str(sdf_file.resolve()),
            "pose_name": f"{row['ligand']}__{row['protein']}/{rel}",
            "file_format": "sdf",
            "optimizer": optimizer,
        }
        if confidence is not None:
            info["confidence"] = confidence
        poses.append(info)
    return poses


def _equibind_pocket_source(filename: str) -> str:
    """Derive the EquiBind pose's pocket source from its filename.

    EquiBind emits three kinds of pose per complex, encoded in the file stem:
      unguided_NNN.sdf            -> "unguided"  (blind EquiBind, no pocket)
      fpocket_raw_pXXX_poseYY.sdf -> "fpocket"   (guided by an fpocket pocket)
      p2rank_raw_pXXX_poseYY.sdf  -> "p2rank"    (guided by a p2rank pocket)
    """
    n = filename.lower()
    if n.startswith("unguided"):
        return "unguided"
    if n.startswith("p2rank"):
        return "p2rank"
    if n.startswith("fpocket"):
        return "fpocket"
    return "unknown"


# EquiBind embeds per-pose provenance both as a filename suffix
# (__clampON / __refSMINA / …) and as SDF tags written by
# equibind_pipeline.pipeline._tag_pose_source. We prefer the SDF tags (the
# authoritative record) and fall back to the filename so older poses — which
# carry neither tags nor suffixes (clamp_mode='on', refine_mode='off') — still
# resolve their pocket source from the filename prefix.
_EQ_PROVENANCE_TAGS = (
    "pocket_source", "pocket_id", "clamp_variant", "refine_variant",
    "smina_affinity", "gnina_affinity",
    # New EquiBind outputs persist the exact final-stage receptor provenance.
    "validation_receptor", "receptor_stage", "receptor_sha256",
    "refine_succeeded", "uff_context",
)


def _read_sdf_tags(sdf_path: Path, tags: tuple[str, ...]) -> dict[str, str]:
    """Read selected ``>  <tag>`` SDF properties at the text level (no RDKit).

    Returns {tag: value} for those *tags* present with a non-empty value.
    """
    want = set(tags)
    out: dict[str, str] = {}
    try:
        lines = Path(sdf_path).read_text(errors="ignore").splitlines()
    except OSError:
        return out
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith(">") and "<" in s and s.endswith(">"):
            tag = s[s.index("<") + 1: s.rindex(">")]
            if tag in want and i + 1 < len(lines):
                val = lines[i + 1].strip()
                if val:
                    out[tag] = val
    return out


def _equibind_clamp_variant(filename: str) -> str | None:
    """Centroid-clamp variant from an EquiBind pose filename (clamp_mode='both')."""
    n = filename.lower()
    if "__clampoff" in n:
        return "clampOFF"
    if "__clampon" in n:
        return "clampON"
    return None


def _equibind_refine_variant(filename: str) -> str | None:
    """Docking re-search variant from an EquiBind pose filename (refine_mode='both')."""
    n = filename.lower()
    if "__refsmina" in n:
        return "smina"
    if "__refgnina" in n:
        return "gnina"
    if "__refraw" in n:
        return "raw"
    return None


def _equibind_result_metadata(combo_dir: Path) -> dict[str, dict[str, Any]]:
    """Return finalized-pose metadata keyed by resolved path and basename.

    Existing benchmark outputs predate receptor tags in their SDFs, but their
    ``guided_results.json`` records whether UFF/refinement succeeded. That is
    required to distinguish a successful smina/GNINA pose (full receptor) from
    a failed refinement that merely copied the raw EquiBind geometry.
    """
    path = combo_dir / "guided_results.json"
    try:
        rows = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    out: dict[str, dict[str, Any]] = {}
    combo = combo_dir.resolve()
    for row in rows if isinstance(rows, list) else []:
        pose_path = row.get("sdf_path")
        if not pose_path:
            continue
        p = Path(str(pose_path)).expanduser()
        if not p.is_absolute():
            p = combo / p
        key = str(p.resolve())
        if key in out:
            raise ValueError(f"Duplicate EquiBind result metadata path: {key}")
        out[key] = row
    return out


def _equibind_finalized_sdfs(combo_dir: Path) -> list[Path]:
    """Return only EquiBind poses committed by its finalized-results record.

    Current benchmark artifacts use ``guided_results.json``. Future outputs may
    be accepted without that legacy file only when every SDF carries an exact
    ``validation_receptor`` tag. This deliberately excludes scratch/partial
    runs such as the unfinished 7WCF directory.
    """
    combo = combo_dir.resolve()
    results_path = combo / "guided_results.json"
    if results_path.is_file():
        try:
            records = json.loads(results_path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"Invalid EquiBind results manifest {results_path}: {exc}") from exc
        if not isinstance(records, list):
            raise ValueError(f"EquiBind results manifest is not a list: {results_path}")

        finalized: list[Path] = []
        seen: set[str] = set()
        for record in records:
            if not isinstance(record, dict) or record.get("success") is not True:
                continue
            raw_path = record.get("sdf_path")
            if not raw_path:
                continue
            path = Path(str(raw_path)).expanduser()
            if not path.is_absolute():
                path = combo / path
            path = path.resolve()
            if path.parent != combo:
                raise ValueError(
                    f"EquiBind manifest pose escapes its combo directory: {path}"
                )
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(
                    f"Finalized EquiBind pose is missing or empty: {path}"
                )
            key = str(path)
            if key in seen:
                raise ValueError(f"Duplicate EquiBind finalized pose: {path}")
            seen.add(key)
            finalized.append(path)
        return sorted(finalized)

    tagged: list[Path] = []
    uncommitted = [
        path for path in sorted(combo.glob("*.sdf"))
        if path.is_file() and path.stat().st_size > 0
    ]
    for path in uncommitted:
        tags = _read_sdf_tags(path, ("validation_receptor",))
        if tags.get("validation_receptor"):
            tagged.append(path.resolve())
    if uncommitted and len(tagged) != len(uncommitted):
        print(
            f"  WARNING: excluding incomplete EquiBind combo {combo.name}: "
            f"{len(uncommitted)} SDF(s), {len(tagged)} exact receptor tag(s), "
            "and no guided_results.json"
        )
        return []
    return tagged


def _expand_equibind_poses(row: dict, _conv_dir: Path, _ctx: PipelineConfig) -> list[dict]:
    file_path = Path(row["file_path"])
    if not file_path.is_dir():
        return []
    poses = []
    result_meta = _equibind_result_metadata(file_path.resolve())
    for sdf_file in _equibind_finalized_sdfs(file_path):
        try:
            rel = str(sdf_file.relative_to(file_path))
        except ValueError:
            rel = sdf_file.name
        fname = sdf_file.name
        # Full per-pose provenance so every EquiBind variant axis is analysable
        # downstream without re-parsing filenames: pocket source (fpocket /
        # p2rank / unguided), centroid-clamp and re-search variants, the matched
        # pocket id and smina/gnina affinity. SDF tags win; the filename is the fallback.
        tags = _read_sdf_tags(sdf_file, _EQ_PROVENANCE_TAGS)
        legacy_meta = result_meta.get(str(sdf_file.resolve())) or {}
        pose = {
            "method": row["docking_tool"], "protein": row["protein"], "ligand": row["ligand"],
            "pose_file": str(sdf_file),
            "source_pose_file": str(sdf_file.resolve()),
            "pose_name": f"{row['ligand']}__{row['protein']}/{rel}",
            "file_format": "sdf",
            "pocket_source": tags.get("pocket_source") or _equibind_pocket_source(fname),
        }
        clamp = tags.get("clamp_variant") or _equibind_clamp_variant(fname)
        refine = tags.get("refine_variant") or _equibind_refine_variant(fname)
        if clamp:
            pose["clamp_variant"] = clamp
        if refine:
            pose["refine_variant"] = refine
        if tags.get("pocket_id"):
            pose["pocket_id"] = tags["pocket_id"]
        if "uff_minimized" in legacy_meta:
            pose["uff_minimized"] = bool(legacy_meta.get("uff_minimized"))
        if legacy_meta.get("error"):
            pose["refine_error"] = str(legacy_meta["error"])
        if "refine_succeeded" in tags:
            pose["refine_succeeded"] = tags["refine_succeeded"].strip().lower() == "true"
        elif (refine in {"smina", "gnina"} and legacy_meta
              and legacy_meta.get("success") is True and "error" in legacy_meta):
            pose["refine_succeeded"] = not str(legacy_meta.get("error", "")).startswith("refine_failed")
        for key in ("validation_receptor", "receptor_stage", "receptor_sha256", "uff_context"):
            if tags.get(key):
                pose[key] = tags[key]
        for _aff_tag in ("smina_affinity", "gnina_affinity"):
            if tags.get(_aff_tag):
                try:
                    pose[_aff_tag] = float(tags[_aff_tag])
                except ValueError:
                    pass
        poses.append(pose)
    return poses


POSE_EXPANDERS = {
    "autodock": _expand_autodock_poses,
    "autodock_vinardo": _expand_autodock_poses,   # same Vina-output tree, Vinardo scoring
    "unidock": _expand_unidock_poses,
    "unidock_tiled": _expand_unidock_poses,
    "unidock2": _expand_unidock2_poses,
    "diffdock": _expand_diffdock_poses,
    "equibind_guided": _expand_equibind_poses,
    "equibind_exclusion": _expand_equibind_poses,
    "equibind_docked_poses": _expand_equibind_poses,
}


# ============================================================================
# PER-POSE VALIDATION RECEPTOR PROVENANCE
# ============================================================================

_VALIDATION_SCHEMA = "prepared-receptor-stock-posebusters-v3"
_EQ_RECEPTOR_MATERIALIZER_SCHEMA = "equibind-effective-graph-v1"
_FILE_HASH_CACHE: dict[tuple[str, int, int, int], str] = {}


def _sha256_file(path: str | Path) -> str:
    """Content SHA-256 with a process-local stat-keyed cache."""
    p = Path(path).resolve()
    st = p.stat()
    # ctime closes the same-size/restored-mtime hole while still allowing shared
    # multi-pose receptor artifacts to be hashed only once per process.
    key = (str(p), int(st.st_size), int(st.st_mtime_ns), int(st.st_ctime_ns))
    cached = _FILE_HASH_CACHE.get(key)
    if cached:
        return cached
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    _FILE_HASH_CACHE[key] = value
    return value


@lru_cache(maxsize=4)
def _posebusters_runtime_signature(config_mode: str) -> str:
    """Fingerprint PoseBusters/RDKit code and config, including local patches."""
    import rdkit
    import posebusters as pb_package

    package_root = Path(pb_package.__file__).resolve().parent
    digest = hashlib.sha256()
    digest.update(_VALIDATION_SCHEMA.encode())
    digest.update(str(config_mode).encode())
    digest.update(importlib.metadata.version("posebusters").encode())
    digest.update(str(rdkit.__version__).encode())
    for path in sorted(package_root.rglob("*")):
        if path.is_file() and path.suffix.lower() in {".py", ".yml", ".yaml"}:
            digest.update(str(path.relative_to(package_root)).encode())
            digest.update(bytes.fromhex(_sha256_file(path)))
    return digest.hexdigest()


def _validation_root(ctx: PipelineConfig, method: str) -> Path | None:
    root = ctx.validation_receptor_directories.get(method)
    if root is None:
        root = ctx.validation_receptor_directories.get(_variant_base_method(method))
    return Path(root).resolve() if root is not None else None


def _pdb_element(line: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    atom_name = "".join(ch for ch in line[12:16] if ch.isalpha()).upper()
    return atom_name[:1]


def _pdb_residue_key(line: str) -> tuple[str, str, str, str, str]:
    return (line[:6], line[17:20], line[21:22], line[22:26], line[26:27])


def _crop_restore_offset(lines: list[str]) -> tuple[float, float, float]:
    for line in lines:
        if "To restore original coords: add offset" not in line:
            continue
        match = re.search(
            r"add offset\s*\(\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*\)",
            line,
        )
        if match:
            return tuple(float(match.group(i)) for i in (1, 2, 3))
    return 0.0, 0.0, 0.0


@lru_cache(maxsize=2048)
def _materialize_equibind_receptor(
    source: Path,
    output_dir: Path,
    mode: str,
) -> Path:
    """Write a PB-ready receptor matching EquiBind's effective atom membership.

    ``graph`` retains exactly residues accepted by EquiBind inference (non-water
    residues containing N, CA and C) and restores translated crop coordinates.
    Dynamic UFF receptor subsets are not reconstructed heuristically: those must
    be tagged by the producing pipeline with their exact validation artifact.
    """
    source = source.resolve()
    source_hash = _sha256_file(source)
    target = (
        output_dir / _normalize_protein(source.stem)
        / f"{source_hash[:20]}_{_EQ_RECEPTOR_MATERIALIZER_SCHEMA}_{mode}.pdb"
    )
    lines = source.read_text(errors="replace").splitlines()
    atoms = [line for line in lines if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms:
        raise ValueError(f"No PDB atoms found in EquiBind receptor {source}")

    keep_keys: set[tuple[str, str, str, str, str]] = set()
    if mode == "graph":
        names_by_residue: dict[tuple[str, str, str, str, str], set[str]] = {}
        for line in atoms:
            key = _pdb_residue_key(line)
            names_by_residue.setdefault(key, set()).add(line[12:16].strip())
        keep_keys = {
            key for key, names in names_by_residue.items()
            if key[1].strip().upper() != "HOH" and {"N", "CA", "C"} <= names
        }
    else:
        raise ValueError(f"Unknown EquiBind receptor materialization mode: {mode}")

    ox, oy, oz = _crop_restore_offset(lines)
    output_lines: list[str] = [
        f"REMARK   VALIDATION RECEPTOR source_sha256={source_hash}\n",
        f"REMARK   VALIDATION RECEPTOR mode={mode}\n",
    ]
    serial = 1
    for line in atoms:
        if _pdb_residue_key(line) not in keep_keys:
            continue
        try:
            x = float(line[30:38]) + ox
            y = float(line[38:46]) + oy
            z = float(line[46:54]) + oz
        except ValueError:
            continue
        padded = line.ljust(80)
        output_lines.append(
            f"{padded[:6]}{serial:5d}{padded[11:30]}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}{padded[54:].rstrip()}\n"
        )
        serial += 1
    output_lines.append("END\n")
    if serial == 1:
        raise ValueError(f"EquiBind {mode} receptor contains no effective atoms: {source}")

    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(".pdb.tmp")
    temporary.write_text("".join(output_lines))
    temporary.replace(target)
    return target.resolve()


def _single_existing(candidates: list[Path], description: str) -> Path:
    existing = sorted({p.resolve() for p in candidates if p.is_file()})
    if len(existing) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {description}; found {len(existing)}: "
            + ", ".join(str(p) for p in existing[:5])
        )
    return existing[0]


_AD4_ELEMENT = {
    "H": "H", "HD": "H", "HS": "H",
    "C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
    "P": "P", "F": "F", "Cl": "CL", "Br": "BR", "I": "I",
    "Si": "SI", "B": "B", "Zn": "ZN", "Fe": "FE", "Mg": "MG",
    "Mn": "MN", "Ca": "CA", "Na": "NA", "K": "K", "Co": "CO",
    "Cu": "CU", "Ni": "NI", "Cd": "CD", "Hg": "HG",
}


def _prepared_heavy_atoms(path: Path, pdbqt: bool) -> Counter:
    """Parse PB-relevant heavy-atom membership from a prepared receptor."""
    atoms: Counter = Counter()
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if pdbqt:
            fields = line.split()
            atom_type = fields[-1] if fields else ""
            element = _AD4_ELEMENT.get(atom_type)
            if element is None:
                raise ValueError(f"Unknown AutoDock atom type {atom_type!r} in {path}")
        else:
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                raise ValueError(f"Missing PDB element field in {path}: {line!r}")
        if element == "H":
            continue
        try:
            coordinates = tuple(
                int(round(float(line[start:end]) * 1000))
                for start, end in ((30, 38), (38, 46), (46, 54))
            )
        except ValueError as exc:
            raise ValueError(f"Invalid receptor coordinates in {path}: {line!r}") from exc
        atoms[(line[:6].strip(), line[17:20].strip(), element, *coordinates)] += 1
    if not atoms:
        raise ValueError(f"No heavy receptor atoms found in {path}")
    return atoms


@lru_cache(maxsize=1024)
def _assert_vina_receptor_parity(pdb_path: str, pdbqt_path: str) -> None:
    """Fail unless the PB PDB exactly represents the heavy atoms Vina saw."""
    pdb_atoms = _prepared_heavy_atoms(Path(pdb_path), pdbqt=False)
    pdbqt_atoms = _prepared_heavy_atoms(Path(pdbqt_path), pdbqt=True)
    if pdb_atoms == pdbqt_atoms:
        return
    pdb_only = pdb_atoms - pdbqt_atoms
    pdbqt_only = pdbqt_atoms - pdb_atoms
    raise ValueError(
        "Prepared Vina PDB/PDBQT heavy-atom mismatch: "
        f"PDB-only={sum(pdb_only.values())}, PDBQT-only={sum(pdbqt_only.values())}; "
        f"PDB-only sample={list(pdb_only.elements())[:3]!r}; "
        f"PDBQT-only sample={list(pdbqt_only.elements())[:3]!r}"
    )


def _resolve_autodock_receptor(pose: dict, root: Path) -> tuple[Path, Path, str]:
    protein = pose["protein"]
    source_pose = Path(str(pose.get("source_pose_file") or "")).resolve()
    if not source_pose.is_file() or not source_pose.name.endswith("_vina_out.pdbqt"):
        raise FileNotFoundError(f"Missing authoritative Vina output provenance: {source_pose}")
    output_stem = source_pose.stem[: -len("_vina_out")]
    parts = output_stem.split("__", 1)
    if len(parts) != 2:
        raise ValueError(f"Unexpected Vina output name: {source_pose.name}")
    receptor_stem = parts[0]
    if source_pose.parent.name != "docking":
        raise ValueError(f"Vina output is not under a docking directory: {source_pose}")
    complex_dir = source_pose.parent.parent.parent
    if complex_dir.name != protein or not complex_dir.is_relative_to(root):
        raise ValueError(
            f"Vina provenance/root mismatch for {protein}: {source_pose} (root {root})"
        )
    receptor = _single_existing(
        [complex_dir / "_staging" / "receptors" / "pdbqt" / f"{receptor_stem}.pdb"],
        f"prepared Vina receptor for {protein}",
    )
    prepared = _single_existing(
        [complex_dir / "_staging" / "receptors" / "pdbqt" / f"{receptor_stem}.pdbqt"],
        f"prepared Vina PDBQT receptor for {protein}",
    )
    _assert_vina_receptor_parity(str(receptor), str(prepared))
    return receptor, prepared, "vina_prepared_pdb_exact_to_pdbqt"


def _resolve_unidock_receptor(pose: dict, root: Path) -> tuple[Path, Path, str]:
    """Resolve the exact prepared receptor recorded by tiled Uni-Dock."""
    raw = str(pose.get("docking_receptor_file") or "").strip()
    if not raw:
        raise FileNotFoundError(
            f"Uni-Dock summary did not record its receptor for {pose['pose_file']}"
        )
    prepared = Path(raw).expanduser().resolve()
    if not prepared.is_file() or prepared.suffix.lower() != ".pdbqt":
        raise FileNotFoundError(f"Missing prepared Uni-Dock receptor: {prepared}")
    if not prepared.is_relative_to(root.resolve()):
        raise ValueError(
            f"Uni-Dock receptor is outside its configured preparation root: "
            f"{prepared} (root {root})"
        )
    receptor = prepared.with_suffix(".pdb")
    if not receptor.is_file():
        raise FileNotFoundError(
            f"Missing PDB companion for prepared Uni-Dock receptor: {receptor}"
        )
    _assert_vina_receptor_parity(str(receptor), str(prepared))
    return receptor, prepared, "unidock_tiled_vina_prepared_pdb_exact_to_pdbqt"


def _resolve_unidock2_receptor(pose: dict, root: Path) -> tuple[Path, Path, str]:
    """Use the exact engine-prepared receptor attested by the schema-3 sentinel.

    The raw staged PDB is only an *input* to Uni-Dock2's pdbfixer/AmberTools
    preparation. Validating clashes against it would not reproduce the atoms the
    engine scored. The transactional driver therefore materializes and hashes its
    prepared receptor beside the output; the collector carries that identity here.
    """
    del root  # authoritative receptor is self-contained beside the committed output
    source_raw = str(pose.get("source_pose_file") or "").strip()
    manifest_raw = str(pose.get("unidock2_completion_file") or "").strip()
    receptor_raw = str(pose.get("unidock2_prepared_receptor_file") or "").strip()
    expected_sha = str(pose.get("unidock2_prepared_receptor_sha256") or "").strip()
    if not source_raw or not manifest_raw or not receptor_raw or not expected_sha:
        raise ValueError(
            f"Uni-Dock2 pose lacks schema-3 prepared-receptor provenance: {pose['pose_file']}"
        )

    source = Path(source_raw).expanduser().resolve()
    suffix = "_unidock2_out"
    if not source.is_file() or not source.stem.endswith(suffix):
        raise ValueError(f"Missing authoritative Uni-Dock2 source SDF: {source}")
    combo = source.stem[: -len(suffix)]
    manifest = source.parent / f"{combo}_unidock2_completion.json"
    receptor = source.parent / f"{combo}_unidock2_receptor_prepared.pdb"
    if (Path(manifest_raw).expanduser().resolve() != manifest.resolve()
            or Path(receptor_raw).expanduser().resolve() != receptor.resolve()
            or str(pose.get("docking_receptor_file") or "").strip()
            and Path(str(pose["docking_receptor_file"])).expanduser().resolve()
            != receptor.resolve()):
        raise ValueError(
            f"Uni-Dock2 provenance does not point beside its source SDF: {source}"
        )
    if receptor.is_symlink():
        raise ValueError(f"Uni-Dock2 prepared receptor must be a retained regular file: {receptor}")
    try:
        completion = json.loads(manifest.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"Invalid Uni-Dock2 completion manifest: {manifest}") from exc
    if (not isinstance(completion, dict)
            or completion.get("schema_version") != _UNIDOCK2_SENTINEL_SCHEMA
            or completion.get("engine") != "unidock2"
            or completion.get("status") != "success"
            or completion.get("output_file") != source.name
            or completion.get("output_sha256") != pose.get("unidock2_output_sha256")
            or completion.get("output_sha256") != _sha256_file(source)
            or completion.get("fingerprint") != pose.get("unidock2_fingerprint")
            or completion.get("generation_id") != pose.get("unidock2_generation_id")
            or not isinstance(completion.get("generation_id"), str)
            or re.fullmatch(r"[0-9a-f]{32}", completion["generation_id"]) is None
            or completion.get("prepared_receptor_file") != receptor.name
            or completion.get("prepared_receptor_sha256") != expected_sha
            or not _unidock2_manifest_fingerprint_matches(completion)
            or not _unidock2_benchmark_protocol_matches(completion)
            or not _unidock2_runtime_fields_match(completion)):
        raise ValueError(f"Uni-Dock2 completion provenance is stale: {manifest}")

    counts = _pdb_atom_counts(receptor)
    if counts is None:
        raise ValueError(f"Uni-Dock2 prepared receptor is malformed: {receptor}")
    atoms, heavy = counts
    if (not receptor.is_file()
            or _sha256_file(receptor) != expected_sha
            or _positive_int(completion.get("prepared_receptor_atoms")) != atoms
            or _positive_int(completion.get("prepared_receptor_heavy_atoms")) != heavy
            or _positive_int(completion.get("engine_receptor_heavy_atoms_in_box")) != heavy
            or _positive_int(pose.get("unidock2_prepared_receptor_atoms")) != atoms
            or _positive_int(pose.get("unidock2_prepared_receptor_heavy_atoms")) != heavy
            or _positive_int(pose.get("unidock2_engine_receptor_heavy_atoms_in_box")) != heavy):
        raise ValueError(
            f"Uni-Dock2 prepared receptor is missing, malformed, or stale: {receptor}"
        )
    return receptor, receptor, "unidock2_manifest_attested_engine_prepared_pdb"


def _resolve_diffdock_receptor(pose: dict, root: Path) -> tuple[Path, Path, str]:
    protein = pose["protein"]
    source_pose = Path(str(pose.get("source_pose_file") or pose["pose_file"])).resolve()
    protein_dir = (root / protein).resolve()
    if not source_pose.is_file() or not source_pose.is_relative_to(protein_dir):
        raise ValueError(
            f"DiffDock pose provenance/root mismatch for {protein}: {source_pose}"
        )
    relative = source_pose.relative_to(protein_dir)
    if len(relative.parts) < 2 or "__" not in relative.parts[0]:
        raise ValueError(f"Unexpected DiffDock artifact path: {source_pose}")
    combo_protein = _normalize_protein(relative.parts[0].split("__", 1)[1])
    if combo_protein != protein:
        raise ValueError(
            f"DiffDock combo/receptor mismatch: {relative.parts[0]} vs {protein}"
        )
    exact = root / protein / "prepared_proteins" / f"{protein}_protein_prepared.pdb"
    candidates = [exact]
    prep_dir = root / protein / "prepared_proteins"
    if prep_dir.is_dir():
        candidates.extend(prep_dir.glob("*.pdb"))
    receptor = _single_existing(candidates, f"prepared DiffDock receptor for {protein}")
    return receptor, receptor, "diffdock_prepared_pdb"


def _resolve_tagged_receptor(pose: dict) -> tuple[Path, Path, str] | None:
    tagged = pose.get("validation_receptor")
    if not tagged:
        return None
    tagged_path = Path(str(tagged)).expanduser()
    if not tagged_path.is_absolute():
        tagged_path = Path(pose["pose_file"]).resolve().parent / tagged_path
    tagged_path = tagged_path.resolve()
    if not tagged_path.is_file():
        raise FileNotFoundError(f"Tagged validation receptor is missing: {tagged_path}")
    expected = str(pose.get("receptor_sha256", "")).strip().lower()
    if expected and _sha256_file(tagged_path) != expected:
        raise ValueError(f"Tagged receptor SHA-256 mismatch: {tagged_path}")
    return tagged_path, tagged_path, str(pose.get("receptor_stage") or "tagged")


def _equibind_pipeline_config(root: Path, protein: str) -> dict[str, Any]:
    path = root / protein / "pipeline_summary.json"
    try:
        data = json.loads(path.read_text())
        return data.get("config", {}) if isinstance(data, dict) else {}
    except (OSError, json.JSONDecodeError):
        return {}


def _resolve_equibind_receptor(
    pose: dict,
    root: Path,
    ctx: PipelineConfig,
) -> tuple[Path, Path, str]:
    protein = pose["protein"]
    source_pose = Path(str(pose.get("source_pose_file") or pose["pose_file"])).resolve()
    protein_dir = (root / protein).resolve()
    if not source_pose.is_file() or not source_pose.is_relative_to(protein_dir):
        raise ValueError(
            f"EquiBind pose provenance/root mismatch for {protein}: {source_pose}"
        )
    relative = source_pose.relative_to(protein_dir)
    if len(relative.parts) != 2 or "__" not in relative.parts[0]:
        raise ValueError(f"Unexpected finalized EquiBind artifact path: {source_pose}")
    combo_protein = _normalize_protein(relative.parts[0].split("__", 1)[1])
    if combo_protein != protein:
        raise ValueError(
            f"EquiBind combo/receptor mismatch: {relative.parts[0]} vs {protein}"
        )
    tagged = _resolve_tagged_receptor(pose)
    if tagged is not None:
        return tagged

    prep_dir = root / protein / f"_prep_{protein}_protein"
    full_exact = prep_dir / f"{protein}_protein_protein.pdb"
    full_candidates = [full_exact]
    if prep_dir.is_dir():
        full_candidates.extend(prep_dir.glob("*_protein.pdb"))
    full = _single_existing(full_candidates, f"full prepared EquiBind receptor for {protein}")

    refine = str(pose.get("refine_variant") or "raw").lower()
    refine_succeeded = pose.get("refine_succeeded")
    if refine in {"smina", "gnina"}:
        if refine_succeeded is None:
            raise ValueError(
                f"Cannot determine whether EquiBind {refine} refinement succeeded for "
                f"{pose['pose_file']}"
            )
        if bool(refine_succeeded):
            return full, full, f"equibind_{refine}_full_prepared"

    # A failed refinement copied the raw geometry. Determine whether the raw pose
    # was UFF-minimized; legacy summaries explicitly record this per pose.
    uff_minimized = pose.get("uff_minimized")
    if uff_minimized is None:
        pipeline_cfg = _equibind_pipeline_config(root, protein)
        if pipeline_cfg.get("uff_minimize") is False:
            uff_minimized = False
        else:
            raise ValueError(f"Missing EquiBind UFF provenance for {pose['pose_file']}")
    if bool(uff_minimized):
        raise ValueError(
            "EquiBind UFF changed the receptor context, but this legacy pose has "
            f"no exact validation_receptor tag: {pose['pose_file']}"
        )

    pocket_source = str(pose.get("pocket_source") or "unknown").lower()
    if pocket_source == "unguided":
        search_source = full
        stage = "equibind_graph_full"
    else:
        pocket_id = str(pose.get("pocket_id") or "").strip()
        if not pocket_id:
            raise ValueError(f"Missing EquiBind pocket_id for guided pose {pose['pose_file']}")
        crop = prep_dir / f"{full.stem}_crop_{pocket_id}.pdb"
        search_source = _single_existing(
            [crop], f"EquiBind crop {pocket_id} for {protein}")
        stage = "equibind_graph_crop"
    receptor = _materialize_equibind_receptor(
        search_source, ctx.validation_receptors_dir / "equibind", "graph")
    return receptor, search_source, stage


def _resolve_prepared_receptor(pose: dict, ctx: PipelineConfig) -> tuple[Path, Path, str]:
    method = str(pose["method"])
    root = _validation_root(ctx, method)
    if root is None:
        raise FileNotFoundError(f"No validation receptor directory configured for {method}")
    base = _variant_base_method(method)
    if base == "autodock":
        return _resolve_autodock_receptor(pose, root)
    if base == "unidock":
        return _resolve_unidock_receptor(pose, root)
    if base == "unidock2":
        return _resolve_unidock2_receptor(pose, root)
    if base == "diffdock":
        return _resolve_diffdock_receptor(pose, root)
    if base == "equibind":
        return _resolve_equibind_receptor(pose, root, ctx)
    raise ValueError(f"No prepared-receptor resolver for method {method!r}")


def _fingerprint_pose(pose: dict, config_mode: str) -> None:
    pose_path = Path(pose["pose_file"]).resolve()
    source_pose = (
        Path(pose["source_pose_file"]).resolve()
        if pose.get("source_pose_file") else pose_path
    )
    receptor_path = Path(pose["receptor_file"]).resolve() if pose.get("receptor_file") else None
    receptor_source = (
        Path(pose["receptor_source_file"]).resolve()
        if pose.get("receptor_source_file") else None
    )
    pose["pose_sha256"] = _sha256_file(pose_path)
    pose["source_pose_sha256"] = _sha256_file(source_pose)
    pose["receptor_sha256"] = _sha256_file(receptor_path) if receptor_path else "none"
    pose["receptor_source_sha256"] = (
        _sha256_file(receptor_source) if receptor_source else "none"
    )
    runtime = _posebusters_runtime_signature(config_mode)
    pose["posebusters_runtime_signature"] = runtime
    metadata = {
        key: value for key, value in pose.items()
        if key not in {"validation_signature"}
        and isinstance(value, (str, int, float, bool, type(None)))
    }
    payload = {
        "schema": _VALIDATION_SCHEMA,
        "config_mode": config_mode,
        "runtime": runtime,
        "metadata": metadata,
    }
    pose["validation_signature"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def attach_validation_receptors(
    poses: list[dict],
    ctx: PipelineConfig,
    legacy_cache: dict[str, str | None] | None = None,
) -> list[dict]:
    """Attach exact receptor paths/hashes and fail before busting on any gap."""
    failures: list[str] = []
    for pose in poses:
        try:
            if ctx.config_mode == "mol":
                pose["receptor_file"] = None
                pose["receptor_source_file"] = None
                pose["receptor_stage"] = "none"
            elif ctx.receptor_policy == "prepared_per_pose":
                receptor, source, stage = _resolve_prepared_receptor(pose, ctx)
                pose["receptor_file"] = str(receptor)
                pose["receptor_source_file"] = str(source)
                pose["receptor_stage"] = stage
            else:
                receptor = (legacy_cache or {}).get(str(pose["protein"]))
                if not receptor or not Path(receptor).is_file():
                    raise FileNotFoundError(
                        f"No shared receptor resolved for {pose['protein']}"
                    )
                pose["receptor_file"] = str(Path(receptor).resolve())
                pose["receptor_source_file"] = str(Path(receptor).resolve())
                pose["receptor_stage"] = "legacy_shared"
            _fingerprint_pose(pose, ctx.config_mode)
        except Exception as exc:
            failures.append(f"{pose.get('method')} {pose.get('pose_file')}: {exc}")

    if failures and ctx.require_validation_receptors:
        preview = "\n".join(f"  - {msg}" for msg in failures[:20])
        extra = f"\n  ... and {len(failures) - 20} more" if len(failures) > 20 else ""
        raise RuntimeError(
            f"Validation receptor preflight failed for {len(failures)} pose(s):\n"
            f"{preview}{extra}"
        )
    return poses


# ============================================================================
# PDBQT -> SDF CONVERSION (AutoDock only)
# ============================================================================

def _split_pdbqt_models(pdbqt_file: str) -> list[str]:
    with open(pdbqt_file) as f:
        content = f.read()
    models, current, in_model = [], [], False
    for line in content.split("\n"):
        if line.startswith("MODEL"):
            in_model = True
            current = [line]
        elif line.startswith("ENDMDL"):
            current.append(line)
            models.append("\n".join(current))
            current, in_model = [], False
        elif in_model:
            current.append(line)
    if not models and content.strip():
        models = [content]
    return models


def _find_template_mol(ligand_name: str, ctx: PipelineConfig):
    """Return ``(molecule, path)`` for the deterministic best template match."""
    from rdkit import Chem

    search_dirs: list[Path] = list(ctx.ligand_template_dirs)
    # Also try a few common defaults relative to work_dir
    for d in ("ligands", "Ligands", "docking_ready_mgltools",
              "docking_ready_mgltools/ligands", ""):
        candidate = ctx.work_dir / d
        if candidate not in search_dirs:
            search_dirs.append(candidate)

    # Template SDFs may live in per-complex subdirectories (PoseBuster Benchmark
    # layout), so search recursively. Prefer exact-stem matches over fuzzy globs.
    for sdir in search_dirs:
        if not sdir.exists():
            continue
        for pat in (f"{ligand_name}.sdf", f"{ligand_name}_*.sdf", f"*{ligand_name}*.sdf"):
            for match in sorted(sdir.rglob(pat)):
                try:
                    mol = Chem.MolFromMolFile(str(match), removeHs=True, sanitize=True)
                    if mol is not None:
                        if not _QUIET:
                            print(f"    Using bond-order template: {match.relative_to(sdir)}")
                        return mol, match.resolve()
                except Exception:
                    pass
    return None, None


@lru_cache(maxsize=1)
def _conversion_tool_fingerprint() -> str:
    """Fingerprint conversion semantics and installed chemistry tool versions."""
    from rdkit import rdBase

    versions: dict[str, str] = {
        "schema": "pdbqt-to-sdf-manifest-v5",
        "rdkit": rdBase.rdkitVersion,
    }
    for package in ("meeko", "openbabel-wheel", "openbabel"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            pass
    selected_executables = {
        "mk_export": _mk_export_available(),
        "obabel": _obabel_available(),
    }
    for executable, resolved in selected_executables.items():
        if resolved and Path(resolved).is_file():
            versions[f"executable:{executable}"] = _sha256_file(resolved)
    return hashlib.sha256(
        json.dumps(versions, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


_MK_EXPORT_BIN: str | None | bool = False  # False = not yet looked up
_OBABEL_BIN: str | None | bool = False


def _resolve_environment_executable(*names: str) -> str | None:
    """Prefer a tool beside the running Python, then fall back to ``PATH``."""
    environment_bin = Path(sys.executable).resolve().parent
    for name in names:
        candidate = environment_bin / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    for name in names:
        resolved = shutil.which(name)
        if resolved:
            return resolved
    return None


def _mk_export_available() -> str | None:
    """Locate the Meeko ``mk_export.py`` CLI once and cache the result."""
    global _MK_EXPORT_BIN
    if _MK_EXPORT_BIN is False:
        _MK_EXPORT_BIN = _resolve_environment_executable("mk_export.py", "mk_export")
    return _MK_EXPORT_BIN


def _obabel_available() -> str | None:
    """Locate Open Babel in the same environment as the running Python."""
    global _OBABEL_BIN
    if _OBABEL_BIN is False:
        _OBABEL_BIN = _resolve_environment_executable("obabel")
    return _OBABEL_BIN


def _converted_model_is_valid(model_pdbqt: Path, output_sdf: Path) -> bool:
    """Require sane chemistry and preservation of every docked heavy coordinate."""
    from rdkit import Chem

    try:
        mol = next(iter(Chem.SDMolSupplier(
            str(output_sdf), removeHs=False, sanitize=True)), None)
    except Exception:
        mol = None
    if mol is None or any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        return False

    pdbqt_atoms: Counter = Counter()
    ligand_types: dict[str, str | None] = dict(_AD4_ELEMENT)
    ligand_types.update({"CG0": "C", "G0": None})
    try:
        for line in model_pdbqt.read_text(errors="replace").splitlines():
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_type = line.split()[-1]
            if atom_type not in ligand_types:
                return False
            element = ligand_types[atom_type]
            if element in {None, "H"}:  # macrocycle glue dummy or hydrogen
                continue
            xyz = tuple(
                int(round(float(line[start:end]) * 1000))
                for start, end in ((30, 38), (38, 46), (46, 54))
            )
            pdbqt_atoms[(element, *xyz)] += 1
    except (OSError, ValueError):
        return False

    conformer = mol.GetConformer()
    sdf_atoms: Counter = Counter()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        point = conformer.GetAtomPosition(atom.GetIdx())
        xyz = tuple(int(round(value * 1000)) for value in (point.x, point.y, point.z))
        sdf_atoms[(atom.GetSymbol().upper(), *xyz)] += 1
    return bool(pdbqt_atoms) and pdbqt_atoms == sdf_atoms


def _mk_export_model(model_pdbqt: Path, output_sdf: Path) -> bool:
    """Convert a single-model PDBQT to *output_sdf* with Meeko ``mk_export``.

    Returns True only if a chemically sane, radical-free molecule was written;
    otherwise removes any partial output and returns False so the caller can
    fall back to the template/obabel strategies.
    """
    mk_bin = _mk_export_available()
    if not mk_bin:
        return False

    try:
        res = subprocess.run(
            # The Meeko installed in the configured Vina environment exposes
            # ``-s/--write_sdf`` for the output file and ``--suffix`` separately.
            # Use the long option so it cannot be confused with suffix handling.
            [mk_bin, str(model_pdbqt), "--write_sdf", str(output_sdf)],
            capture_output=True, text=True,
        )
    except OSError:
        return False
    if res.returncode != 0 or not output_sdf.exists():
        output_sdf.unlink(missing_ok=True)
        return False

    # Guard chemistry and coordinates so a converter cannot silently publish a
    # different ligand or crystal-template conformer as the docked pose.
    if not _converted_model_is_valid(model_pdbqt, output_sdf):
        output_sdf.unlink(missing_ok=True)
        return False
    return True


def _convert_pdbqt_to_sdf(pdbqt_file: str, out_dir: Path, ctx: PipelineConfig) -> list[str]:
    """Convert a multi-pose PDBQT with a content-verified conversion cache."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pdbqt_path = Path(pdbqt_file)
    pdbqt_path = pdbqt_path.resolve()
    base_name = pdbqt_path.stem
    models = _split_pdbqt_models(pdbqt_file)
    if not models:
        print(f"  Warning: No models found in {pdbqt_file}")
        return []

    stem_clean = base_name.replace("_vina_out", "")
    parts = stem_clean.split("__")
    needs_template = any("REMARK SMILES" not in model for model in models)
    template_mol, template_path = (
        _find_template_mol(_normalize_ligand(parts[1]), ctx)
        if len(parts) == 2 and needs_template else (None, None)
    )

    source_sha = _sha256_file(pdbqt_path)
    template_sha = _sha256_file(template_path) if template_path else "none"
    source_key = hashlib.sha256(
        f"{pdbqt_path}\0{source_sha}".encode()
    ).hexdigest()[:12]
    expected_outputs = [
        out_dir / f"{base_name}_{source_key}_model{i}.sdf"
        for i in range(1, len(models) + 1)
    ]
    manifest_path = out_dir / f"{base_name}_{source_key}_conversion.json"
    cache_inputs = {
        "schema": "pdbqt-to-sdf-manifest-v5",
        "source_path": str(pdbqt_path),
        "source_sha256": source_sha,
        "template_path": str(template_path) if template_path else None,
        "template_sha256": template_sha,
        "converter_fingerprint": _conversion_tool_fingerprint(),
        "model_count": len(models),
    }

    if not ctx.overwrite and manifest_path.is_file():
        try:
            manifest = json.loads(manifest_path.read_text())
            if not isinstance(manifest, dict):
                raise TypeError("conversion manifest is not a mapping")
            outputs = manifest.get("outputs", [])
            if (not isinstance(outputs, list)
                    or any(not isinstance(record, dict) for record in outputs)):
                raise TypeError("conversion manifest outputs are not mappings")
            cache_valid = all(manifest.get(k) == v for k, v in cache_inputs.items())
            cache_valid = cache_valid and len(outputs) == len(expected_outputs)
            if cache_valid:
                for expected, record in zip(expected_outputs, outputs):
                    if (Path(record.get("path", "")).resolve() != expected.resolve()
                            or not expected.is_file()
                            or _sha256_file(expected) != record.get("sha256")):
                        cache_valid = False
                        break
            if cache_valid:
                return [str(path) for path in expected_outputs]
        except (OSError, json.JSONDecodeError, TypeError):
            pass

    # A legacy/existence-only or stale cache is never trusted. Remove only the
    # exact generated targets for this source key; a valid manifest is published
    # after every expected model has converted successfully.
    manifest_path.unlink(missing_ok=True)
    for output in expected_outputs:
        output.unlink(missing_ok=True)

    converted: list[str] = []
    for i, (model_content, output_sdf) in enumerate(
        zip(models, expected_outputs), start=1
    ):

        temp_pdbqt = out_dir / f"{base_name}_model{i}.pdbqt"
        temp_pdb = out_dir / f"{base_name}_model{i}.pdb"
        with open(temp_pdbqt, "w") as f:
            f.write(model_content)

        mol_final = None

        # Strategy 0 (preferred): Meeko mk_export. Ligand prep embeds the input
        # SMILES (`REMARK SMILES`) in the PDBQT, so Meeko rebuilds the exact
        # bond orders and formal charges instead of guessing. This avoids the
        # radical / valence artefacts obabel produces when it has to infer bond
        # orders from coordinates alone (see _mk_export_model). Embedded-SMILES
        # poses fail closed if Meeko cannot reconstruct them; only genuinely
        # non-Meeko PDBQTs without that header use the fallback strategies.
        if "REMARK SMILES" in model_content:
            if _mk_export_model(temp_pdbqt, output_sdf):
                converted.append(str(output_sdf))
                temp_pdbqt.unlink(missing_ok=True)
                continue
            temp_pdbqt.unlink(missing_ok=True)
            manifest_path.unlink(missing_ok=True)
            for partial_output in expected_outputs:
                partial_output.unlink(missing_ok=True)
            raise RuntimeError(
                "Meeko failed to reconstruct an embedded-SMILES Vina pose "
                f"without changing its docked heavy coordinates: {pdbqt_path} model {i}"
            )

        # Strategy 1: obabel PDBQT -> PDB -> RDKit + template
        if template_mol is not None:
            try:
                obabel_bin = _obabel_available()
                if not obabel_bin:
                    raise FileNotFoundError("obabel is unavailable")
                res = subprocess.run(
                    [obabel_bin, str(temp_pdbqt), "-O", str(temp_pdb)],
                    capture_output=True, text=True,
                )
                if res.returncode == 0 and temp_pdb.exists():
                    raw = Chem.MolFromPDBFile(str(temp_pdb), removeHs=True, sanitize=False)
                    if raw is not None:
                        try:
                            mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw)
                            Chem.SanitizeMol(mol_final)
                        except Exception as e:
                            print(f"    Template assignment failed for model {i}: {e}")
                            mol_final = None
                if temp_pdb.exists():
                    temp_pdb.unlink()
            except FileNotFoundError:
                pass

        # Strategy 2: obabel PDBQT -> SDF directly
        if mol_final is None:
            try:
                obabel_bin = _obabel_available()
                if not obabel_bin:
                    raise FileNotFoundError("obabel is unavailable")
                res = subprocess.run(
                    [obabel_bin, str(temp_pdbqt), "-O", str(output_sdf)],
                    capture_output=True, text=True,
                )
                if res.returncode == 0 and output_sdf.exists():
                    if template_mol is not None:
                        try:
                            raw = Chem.MolFromMolFile(str(output_sdf), removeHs=True, sanitize=False)
                            if raw is not None:
                                mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw)
                                Chem.SanitizeMol(mol_final)
                        except Exception:
                            if _converted_model_is_valid(temp_pdbqt, output_sdf):
                                converted.append(str(output_sdf))
                                temp_pdbqt.unlink(missing_ok=True)
                                continue
                            output_sdf.unlink(missing_ok=True)
                    else:
                        if _converted_model_is_valid(temp_pdbqt, output_sdf):
                            converted.append(str(output_sdf))
                            temp_pdbqt.unlink(missing_ok=True)
                            continue
                        output_sdf.unlink(missing_ok=True)
            except FileNotFoundError:
                try:
                    raw = Chem.MolFromPDBFile(str(temp_pdbqt), removeHs=True, sanitize=False)
                    if raw is not None and template_mol is not None:
                        mol_final = AllChem.AssignBondOrdersFromTemplate(template_mol, raw)
                        Chem.SanitizeMol(mol_final)
                    elif raw is not None:
                        try:
                            Chem.SanitizeMol(raw)
                        except Exception:
                            pass
                        mol_final = raw
                except Exception as e:
                    print(f"    RDKit fallback failed for model {i}: {e}")

        if mol_final is not None:
            writer = Chem.SDWriter(str(output_sdf))
            writer.write(mol_final)
            writer.close()
            if (output_sdf.exists()
                    and _converted_model_is_valid(temp_pdbqt, output_sdf)):
                converted.append(str(output_sdf))
            else:
                output_sdf.unlink(missing_ok=True)
        elif not output_sdf.exists():
            print(f"    Warning: Could not convert model {i} from {pdbqt_path.name}")

        temp_pdbqt.unlink(missing_ok=True)

    if len(converted) != len(expected_outputs) or any(not p.is_file() for p in expected_outputs):
        manifest_path.unlink(missing_ok=True)
        for partial_output in expected_outputs:
            partial_output.unlink(missing_ok=True)
        raise RuntimeError(
            f"Converted only {len(converted)}/{len(expected_outputs)} models from {pdbqt_path}"
        )

    manifest = dict(cache_inputs)
    manifest["outputs"] = [
        {"model": i, "path": str(path.resolve()), "sha256": _sha256_file(path)}
        for i, path in enumerate(expected_outputs, start=1)
    ]
    manifest_tmp = manifest_path.with_suffix(".json.tmp")
    manifest_tmp.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    manifest_tmp.replace(manifest_path)
    return [str(path) for path in expected_outputs]


def _apply_variant_filter(poses: list[dict], variant_filter: dict[str, set[str]]) -> list[dict]:
    """Keep only poses whose tool-variant is allowed by *variant_filter*.

    *variant_filter* maps a base tool name to the set of labels to keep.
    AutoDock/DiffDock use ``optimizer`` (original/smina/gnina), EquiBind uses
    ``refine_variant`` (raw/smina/gnina), and Uni-Dock uses ``variant=tiled``.
    Poses whose variant label is missing or not allowed are dropped. Returns
    *poses* unchanged when the filter is empty.
    """
    if not variant_filter:
        return poses
    kept: list[dict] = []
    for p in poses:
        base = _variant_base_method(p.get("method", ""))
        allowed = variant_filter.get(base)
        field_name = _VARIANT_FIELD.get(base)
        if not allowed or field_name is None:
            kept.append(p)                       # no filter for this tool / no variant axis
            continue
        val = p.get(field_name)
        val = str(val).strip().lower() if val not in (None, "") else None
        normalized_allowed = allowed
        if base == "autodock":
            aliases = {"raw": "original", "native": "original", "none": "original"}
            val = aliases.get(val, val)
            normalized_allowed = {aliases.get(item, item) for item in allowed}
        if val in normalized_allowed:
            kept.append(p)
    return kept


def collect_all_pose_files(filtered_df: pd.DataFrame, ctx: PipelineConfig) -> list[dict]:
    """Expand every row in *filtered_df* into individual pose-file dicts."""
    all_poses: list[dict] = []
    for _, row in filtered_df.iterrows():
        expander = POSE_EXPANDERS.get(row["docking_tool"])
        if expander is None:
            print(f"  WARNING: No pose expander for '{row['docking_tool']}', skipping.")
            continue
        all_poses.extend(expander(row.to_dict(), ctx.converted_dir, ctx))
    if ctx.variant_filter:
        before = len(all_poses)
        all_poses = _apply_variant_filter(all_poses, ctx.variant_filter)
        desc = ", ".join(f"{k}={'/'.join(sorted(v))}"
                         for k, v in sorted(ctx.variant_filter.items()))
        print(f"   variant_filter [{desc}]: kept {len(all_poses)}/{before} poses")
    return all_poses


# ============================================================================
# PARALLEL WORKER FUNCTIONS  (module-level so they're picklable)
# ============================================================================

_worker_buster_dock: PoseBusters | None = None
_worker_buster_mol: PoseBusters | None = None
_worker_protein_cache: dict[str, str | None] = {}
_worker_config_mode: str = "mol"
_worker_pose_timeout: int | None = None
# Shared across workers (a Manager dict, keyed by ligand): once one pose of a
# ligand times out, every remaining pose of that same molecule is skipped
# instead of re-hanging for another full timeout — EquiBind emits up to 540
# poses per ligand, so one poison molecule would otherwise cost 540 x timeout.
_worker_poison: Any = None


def _init_worker(
    config_mode: str,
    protein_cache: dict[str, str | None],
    quiet: bool = False,
    pose_timeout: int | None = None,
    limit_threads: bool = True,
    poison=None,
) -> None:
    global _worker_buster_dock, _worker_buster_mol, _worker_protein_cache
    global _worker_config_mode, _worker_pose_timeout, _worker_poison
    from posebusters import PoseBusters as _PB

    _silence_rdkit_2d3d_warning()  # workers run bust(); install the filter here too
    _apply_quiet_mode(quiet)       # and drop all warning chatter if requested
    if limit_threads:
        _limit_worker_threads()    # one BLAS/RDKit thread per worker
    if _LIBC is not None:          # die with the main process, so an interrupted
        try:                        # run (killed/crashed main) leaves no orphaned
            _LIBC.prctl(_PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)  # workers behind
        except Exception:
            pass

    _worker_config_mode = config_mode
    _worker_protein_cache = protein_cache or {}
    _worker_pose_timeout = pose_timeout
    _worker_poison = poison
    if config_mode == "dock":
        _worker_buster_dock = _PB(config="dock")
        _worker_buster_mol = None
    else:
        _worker_buster_mol = _PB(config="mol")


class _PoseTimeout(Exception):
    """Raised when a single bust() call exceeds the per-pose wall-clock cap."""


# A single bust() can hang inside RDKit's C++ conformer generation (the
# 50-conformer internal-energy ensemble): all workers pin a core at 100% and
# complete nothing. A SIGALRM cannot break that — the Python handler only runs
# when control returns to the interpreter, which a tight C++ loop never does.
# So we run each bust() in a FORKED grandchild and hard-kill it on timeout:
# SIGTERM/SIGKILL to a child is actioned by the kernel no matter what C++ code it
# is stuck in, so the wedged pose is dropped and the worker moves to the next.
# The grandchild inherits the constructed PoseBusters objects via fork (no
# re-construction) and, being a separate process, still uses exactly one core.
# NB: multiprocessing.Process can't be used here — Pool workers are daemonic and
# "daemonic processes are not allowed to have children". The raw os.fork() syscall
# has no such restriction, so we fork by hand and talk over an os.pipe().
_PR_SET_PDEATHSIG = 1
try:
    import ctypes as _ctypes
    _LIBC = _ctypes.CDLL("libc.so.6", use_errno=True)
except Exception:  # pragma: no cover - non-glibc platform
    _LIBC = None


def _run_bust(use_dock, pose_file, protein_file):
    """Direct bust() in this process — used when no per-pose timeout is set."""
    buster = _worker_buster_dock if use_dock else _worker_buster_mol
    return buster.bust(pose_file, None, protein_file, full_report=True)


def _write_all(fd, data: bytes) -> None:
    while data:
        data = data[os.write(fd, data):]


def _read_exact(fd, n: int, deadline: float) -> bytes | None:
    """Read exactly *n* bytes from *fd* before *deadline* (monotonic); None if not."""
    buf = b""
    while len(buf) < n:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            return None
        if not select.select([fd], [], [], remaining)[0]:
            return None
        chunk = os.read(fd, n - len(buf))
        if not chunk:                      # EOF: child died without full payload
            return None
        buf += chunk
    return buf


def _reap(pid: int, killed: bool) -> None:
    """Reap the forked child *pid*.

    On the success path (*killed* False) the child has already sent its result
    and is calling os._exit, so a blocking waitpid returns in microseconds — no
    polling/sleeping (the old WNOHANG+sleep loop idled the worker ~50 ms per pose
    on a busy box, throttling throughput). On timeout, force-kill it.
    """
    try:
        if not killed:
            os.waitpid(pid, 0)
            return
    except (ChildProcessError, OSError):
        return
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.kill(pid, sig)
        except OSError:
            return
        for _ in range(30):                # up to ~1.5s per signal
            time.sleep(0.05)
            try:
                if os.waitpid(pid, os.WNOHANG)[0]:
                    return
            except (ChildProcessError, OSError):
                return


def _bust_with_timeout(use_dock, pose_file, protein_file, timeout):
    """Run bust() in a hand-forked child, hard-killed if it exceeds *timeout* (s).

    Returns the result DataFrame, or raises _PoseTimeout if the child had to be
    killed (an effectively-infinite RDKit computation). The child inherits the
    constructed PoseBusters objects via fork, so there is no re-construction cost.
    """
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:                           # ---- child ----
        try:
            os.close(r)
            if _LIBC is not None:          # die if the worker (parent) dies
                try:
                    _LIBC.prctl(_PR_SET_PDEATHSIG, signal.SIGKILL, 0, 0, 0)
                except Exception:
                    pass
            try:
                buster = _worker_buster_dock if use_dock else _worker_buster_mol
                df = buster.bust(pose_file, None, protein_file, full_report=True)
                blob = pickle.dumps(("ok", df), protocol=pickle.HIGHEST_PROTOCOL)
            except Exception as e:
                blob = pickle.dumps(("err", f"{type(e).__name__}: {e}"))
            _write_all(w, struct.pack(">Q", len(blob)))
            _write_all(w, blob)
        except Exception:
            pass
        finally:
            os._exit(0)                    # no atexit — parent already has the result

    # ---- parent ----
    os.close(w)
    deadline = time.monotonic() + timeout
    got = False
    try:
        header = _read_exact(r, 8, deadline)
        if header is None:
            raise _PoseTimeout()
        (n,) = struct.unpack(">Q", header)
        body = _read_exact(r, n, deadline)
        if body is None:
            raise _PoseTimeout()
        status, payload = pickle.loads(body)
        got = True
    finally:
        os.close(r)
        _reap(pid, killed=not got)
    if status == "err":
        raise RuntimeError(payload)
    return payload


def _process_single_pose(pose_info: dict):
    """Validate one pose; returns (task_id, result_df | None, error_msg | None).

    The pose_file is always returned (even on error) so the caller can track
    exactly which poses a pool pass consumed — needed to compute the remaining
    set when a stalled pool is terminated and restarted.
    """
    pose_file = pose_info["pose_file"]
    task_id = str(pose_info.get("validation_signature") or pose_file)
    protein_name = pose_info["protein"]
    lig_key = f"{protein_name}||{pose_info.get('ligand', '')}"

    if not os.path.exists(pose_file):
        return task_id, None, f"File not found: {pose_file}"

    # Skip poses of a ligand already known to hang bust() (a sibling pose timed
    # out) — re-running would just burn another full timeout on the same molecule.
    if _worker_poison is not None:
        try:
            poisoned = lig_key in _worker_poison
        except Exception:
            poisoned = False
        if poisoned:
            return task_id, None, f"POISON-SKIP ({lig_key}): {pose_file}"

    try:
        protein_file = None
        used_mode = _worker_config_mode
        timeout = _worker_pose_timeout
        guarded = bool(timeout and timeout > 0)

        def _bust(use_dock, prot):
            return (_bust_with_timeout(use_dock, pose_file, prot, timeout)
                    if guarded else _run_bust(use_dock, pose_file, prot))

        if _worker_config_mode == "dock":
            protein_file = (pose_info.get("receptor_file")
                            or _worker_protein_cache.get(protein_name))
            if protein_file and Path(protein_file).exists():
                df = _bust(True, protein_file)
                used_mode = "dock"
            else:
                return task_id, None, (
                    f"RECEPTOR-ERROR: required receptor is missing for {pose_file}: "
                    f"{protein_file or 'unresolved'}"
                )
        else:
            df = _bust(False, None)
            used_mode = "mol"

        df["docking_method"] = pose_info["method"]
        df["protein"] = protein_name
        df["ligand"] = pose_info["ligand"]
        df["pose_file"] = pose_info["pose_file"]
        df["pose_name"] = pose_info["pose_name"]
        df["file_format"] = pose_info.get("file_format", "sdf")
        df["protein_file_used"] = protein_file or "none"
        df["receptor_file_used"] = protein_file or "none"
        df["posebusters_mode"] = used_mode
        for key in (
            "source_pose_file", "receptor_source_file", "receptor_stage",
            "receptor_sha256", "receptor_source_sha256", "pose_sha256",
            "source_pose_sha256",
            "validation_signature",
            "posebusters_runtime_signature",
        ):
            if key in pose_info:
                df[key] = pose_info[key]
        if "confidence" in pose_info:
            df["diffdock_confidence"] = pose_info["confidence"]
        # AutoDock/DiffDock optimizer provenance (original / smina / gnina) so
        # optimized poses stay distinguishable from their raw source downstream.
        if "optimizer" in pose_info:
            df["optimizer"] = pose_info["optimizer"]
        # AutoDock Vina rank + affinity (kcal/mol). Vina produces them natively but
        # they were never propagated to the CSV; carry them through so AutoDock is
        # rankable from the results table like DiffDock/EquiBind.
        for key in (
            "autodock_rank", "autodock_affinity", "optimized_rank", "rank_metric",
            "minimized_affinity", "cnn_score", "cnn_affinity",
            "optimizer_elapsed_time_s", "optimizer_processing_elapsed_time_s",
            "optimization_log_file",
            "optimizer_provenance_file", "optimizer_provenance_fingerprint",
            "source_pose_model_sha256",
        ):
            if key in pose_info:
                df[key] = pose_info[key]
        # Tiled Uni-Dock's merged MODEL order is its global affinity ranking;
        # Uni-Dock2's best-first SDF record order is its single-shot ranking.
        for key in (
            "unidock_rank", "unidock_affinity", "variant", "scoring", "engine",
            "docking_receptor_file", "docking_summary_file",
            "unidock_manifest_file", "unidock_fingerprint", "unidock_generation_id",
            "unidock2_rank", "unidock2_affinity", "unidock2_completion_file",
            "unidock2_fingerprint", "unidock2_generation_id", "unidock2_output_sha256",
            "unidock2_prepared_receptor_file", "unidock2_prepared_receptor_sha256",
            "unidock2_prepared_receptor_atoms",
            "unidock2_prepared_receptor_heavy_atoms",
            "unidock2_engine_receptor_heavy_atoms_in_box",
        ):
            if key in pose_info:
                df[key] = pose_info[key]
        # Carry EquiBind variant provenance through to the results CSV so every
        # pose axis (pocket source, clamp, re-search, pocket id, smina affinity)
        # is analysable without re-parsing filenames.
        for key in _EQ_PROVENANCE_TAGS:
            if key in pose_info:
                df[key] = pose_info[key]
        return task_id, df, None

    except _PoseTimeout:
        # Mark the whole ligand poison so its remaining poses are skipped, then
        # report. Sentinel prefix "TIMEOUT" lets the main loop tally these
        # separately and record the pose path; the pose is NOT written to the
        # results CSV, so a resume run retries it unless the caller excludes it.
        if _worker_poison is not None:
            try:
                _worker_poison[lig_key] = True
            except Exception:
                pass
        return task_id, None, f"TIMEOUT (>{_worker_pose_timeout}s): {pose_file}"
    except Exception as e:
        return task_id, None, f"Error processing {Path(pose_file).name}: {e}"


# ============================================================================
# RESOURCE MONITOR  (core utilisation + bottleneck detection)
# ============================================================================

def _fmt_gb(nbytes: float) -> str:
    return f"{nbytes / 1024 ** 3:.1f}G"


class _ResourceMonitor:
    """Background sampler: per-process core utilisation + bottleneck flags.

    Runs as a daemon thread in the main process while the Pool loop blocks on
    results, printing one snapshot every *interval* seconds: system CPU / load /
    RAM / swap, then the worker pool's aggregate CPU%, thread count and RSS, the
    derived "cores busy", and pose throughput. Heuristics flag the failure modes
    that stall a long run — CPU oversubscription, swap / memory pressure, core
    under-utilisation, and busy-but-no-progress (the slow / dead-worker
    signature). Uses psutil for per-process detail; without it, still reports
    load average + system memory / swap, which already expose the two main
    bottlenecks. Never raises into the run — a sampling error just prints a note.
    """

    _CPU_BUSY_PCT = 20.0   # a worker above this counts as "busy" this tick

    def __init__(self, interval, ncores, progress_getter=None):
        self.interval = max(2, int(interval))
        self.ncores = ncores or (os.cpu_count() or 1)
        self._get_progress = progress_getter or (lambda: None)
        self._stop = threading.Event()
        self._thread = None
        self._t0 = 0.0
        self._main = None
        self._proc_cache = {}          # pid -> psutil.Process (reused for cpu deltas)
        self._prev_done = 0
        self._prev_swap = None
        self._peak_swap = 0
        self._peak_load = 0.0
        self._min_rate = None
        self._stall_ticks = 0
        try:
            import psutil
            self._psutil = psutil
        except Exception:
            self._psutil = None

    def start(self):
        self._t0 = time.monotonic()
        if self._psutil is not None:
            self._psutil.cpu_percent(interval=None)            # prime system CPU
            try:
                self._main = self._psutil.Process()
                self._main.cpu_percent(interval=None)          # prime main CPU
                self._proc_cache[self._main.pid] = self._main
            except Exception:
                self._main = None
        self._thread = threading.Thread(target=self._run, name="pb-monitor", daemon=True)
        self._thread.start()
        detail = "per-process" if self._psutil is not None else "system-only (no psutil)"
        print(f"  [MON] resource monitor on — every {self.interval}s, "
              f"{self.ncores} logical cores, {detail}")

    def stop(self):
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=self.interval + 2)
        self._final_summary()

    # --- internals ----------------------------------------------------------
    def _run(self):
        # wait() -> True when stop is set (exit); False on timeout (sample)
        while not self._stop.wait(self.interval):
            try:
                self._sample()
            except Exception as e:  # a monitor bug must never take down the run
                print(f"  [MON] sample failed: {e}")

    def _elapsed(self):
        return int(time.monotonic() - self._t0)

    def _sample(self):
        done = self._get_progress()
        t = self._elapsed()
        if self._psutil is None:
            self._sample_basic(t, done)
            return

        ps = self._psutil
        vm = ps.virtual_memory()
        sw = ps.swap_memory()
        sys_cpu = ps.cpu_percent(interval=None)
        try:
            load1 = os.getloadavg()[0]
        except OSError:
            load1 = float("nan")

        main = self._main or ps.Process()
        try:
            children = main.children(recursive=True)
        except Exception:
            children = []

        live = {main.pid} | {c.pid for c in children}
        for pid in list(self._proc_cache):
            if pid not in live:
                self._proc_cache.pop(pid, None)

        try:
            main_cpu = main.cpu_percent(interval=None)
        except Exception:
            main_cpu = 0.0

        worker_cpu, worker_thr, worker_rss, busy = [], 0, 0, 0
        for c in children:
            p = self._proc_cache.get(c.pid)
            first = p is None
            if first:                              # first time we see this worker
                self._proc_cache[c.pid] = c
                p = c
            try:
                # cpu_percent needs a prior call to form a delta; threads/RSS
                # are instantaneous, so collect those even on a worker's 1st tick.
                cpu = None if first else p.cpu_percent(interval=None)
                if first:
                    p.cpu_percent(interval=None)   # prime for the next tick
                worker_thr += p.num_threads()
                worker_rss += p.memory_info().rss
            except Exception:
                continue
            if cpu is not None:
                worker_cpu.append(cpu)
                if cpu > self._CPU_BUSY_PCT:
                    busy += 1

        n_procs = len(children)
        # Cores busy from SYSTEM CPU, not per-process summation: each bust runs in
        # a short-lived forked grandchild that the per-process sampler can't catch
        # a CPU delta for, which made this read ~0 even at full tilt. sys_cpu is
        # robust to that churn.
        cores_busy = sys_cpu / 100.0 * self.ncores

        delta = (done - self._prev_done) if done is not None else None
        thr_str = ""
        if done is not None:
            rate = delta / self.interval
            thr_str = f"  |  poses {done} (+{delta}, {rate:.1f}/s)"
            self._prev_done = done
            if delta > 0:
                self._min_rate = rate if self._min_rate is None else min(self._min_rate, rate)

        self._peak_swap = max(self._peak_swap, sw.used)
        if load1 == load1:
            self._peak_load = max(self._peak_load, load1)

        print(
            f"  [MON t+{t}s] cores~{cores_busy:.1f}/{self.ncores} (sys-CPU {sys_cpu:.0f}%)  "
            f"load {load1:.1f}/{self.ncores}  "
            f"RAM {vm.percent:.0f}% ({_fmt_gb(vm.used)}/{_fmt_gb(vm.total)})  "
            f"swap {_fmt_gb(sw.used)}\n"
            f"             procs {n_procs}  threads {worker_thr}  "
            f"RSS {_fmt_gb(worker_rss)}  main-CPU {main_cpu:.0f}%{thr_str}"
        )
        self._flag(load1, vm, sw, cores_busy, delta)

    def _flag(self, load1, vm, sw, cores_busy, delta):
        out = []
        if load1 == load1 and load1 > self.ncores * 1.5:
            out.append(f"HIGH LOAD {load1:.0f} vs {self.ncores} cores — CPU oversubscribed"
                       " (use fewer workers, or check limit_worker_threads)")

        if self._prev_swap is not None and sw.used > self._prev_swap + 256 * 1024 ** 2:
            out.append(f"SWAP GROWING +{_fmt_gb(sw.used - self._prev_swap)} this tick"
                       " — memory-driven stall risk; lower num_workers / worker_maxtasks")
        elif vm.percent >= 92:
            out.append(f"MEMORY PRESSURE: {vm.percent:.0f}% RAM used — recycle workers sooner")
        self._prev_swap = sw.used

        if delta is not None and delta == 0 and cores_busy > self.ncores * 0.25:
            self._stall_ticks += 1
            out.append(f"BUSY BUT 0 POSES DONE for {self._stall_ticks * self.interval}s"
                       " — a pose is hanging; it will be killed at pose_timeout")
        else:
            self._stall_ticks = 0
            if delta is not None and cores_busy < self.ncores * 0.4:
                out.append(f"UNDER-UTILISED: only ~{cores_busy:.1f}/{self.ncores} cores busy"
                           " — too few workers (num_workers 0/blank = all cores), a hang"
                           " holding workers, or main-process/I/O serialisation")

        for f in out:
            print(f"             ⚠ {f}")

    def _sample_basic(self, t, done):
        try:
            load1 = os.getloadavg()[0]
        except OSError:
            load1 = float("nan")
        thr_str = ""
        if done is not None:
            delta = done - self._prev_done
            thr_str = f"  |  poses {done} (+{delta}, {delta / self.interval:.1f}/s)"
            self._prev_done = done
        print(f"  [MON t+{t}s] load {load1:.1f}/{self.ncores}  {self._read_meminfo()}{thr_str}")
        if load1 == load1 and load1 > self.ncores * 1.5:
            print(f"             ⚠ HIGH LOAD: {load1:.0f} vs {self.ncores} cores — oversubscribed")

    @staticmethod
    def _read_meminfo():
        try:
            info = {}
            for line in Path("/proc/meminfo").read_text().splitlines():
                k, _, rest = line.partition(":")
                info[k] = int(rest.strip().split()[0]) * 1024   # kB -> bytes
            used = info["MemTotal"] - info.get("MemAvailable", info["MemTotal"])
            swap_used = info.get("SwapTotal", 0) - info.get("SwapFree", 0)
            return f"RAM {_fmt_gb(used)}/{_fmt_gb(info['MemTotal'])}  swap {_fmt_gb(swap_used)}"
        except Exception:
            return "mem n/a"

    def _final_summary(self):
        bits = []
        if self._peak_load:
            bits.append(f"peak load {self._peak_load:.0f}/{self.ncores}")
        if self._peak_swap:
            bits.append(f"peak swap {_fmt_gb(self._peak_swap)}")
        if self._min_rate is not None:
            bits.append(f"slowest {self._min_rate:.1f} poses/s")
        if bits:
            print(f"  [MON] monitor off — {' | '.join(bits)}")


# ============================================================================
# CHECKPOINT + ANALYSIS
# ============================================================================

def _save_checkpoint(frames, output_file, n_new, n_pending, final=False):
    df_out = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    temporary = output_file.with_suffix(output_file.suffix + ".tmp")
    df_out.to_csv(temporary, index=False)
    temporary.replace(output_file)
    tag = "FINAL" if final else "CHECKPOINT"
    print(f"    [{tag}] {n_new}/{n_pending} new  |  total rows: {len(df_out)}  →  {output_file.name}")


def _reconcile_previous_results(
    previous: pd.DataFrame,
    poses_list: list[dict],
) -> tuple[pd.DataFrame, set[str]]:
    """Keep only prior rows whose content-aware validation signature is current."""
    signatures = [str(p.get("validation_signature") or "") for p in poses_list]
    if any(not value for value in signatures):
        raise ValueError("Every pose must have a validation_signature before resume reconciliation")
    duplicates = [sig for sig, count in Counter(signatures).items() if count > 1]
    if duplicates:
        raise ValueError(
            f"Current pose set contains {len(duplicates)} duplicate validation signatures"
        )
    if previous.empty:
        return previous, set()
    if "validation_signature" not in previous.columns:
        print(
            "  Existing CSV predates content-aware receptor provenance; "
            "all rows are stale and will be recomputed."
        )
        return previous.iloc[0:0].copy(), set()
    current = set(signatures)
    valid = previous[
        previous["validation_signature"].fillna("").astype(str).isin(current)
    ].copy()
    valid = valid.drop_duplicates(subset=["validation_signature"], keep="last")
    done = set(valid["validation_signature"].astype(str))
    n_stale = len(previous) - len(valid)
    if n_stale:
        print(f"  Discarding {n_stale} stale/out-of-scope previous result row(s)")
    return valid, done


def analyze_poses_with_posebusters(
    poses_list: list[dict],
    protein_cache: dict[str, str | None],
    config: str = "mol",
    output_file: Path | None = None,
    save_interval: int = 100,
    overwrite: bool = False,
    num_workers: int | None = None,
    quiet: bool = False,
    worker_maxtasks: int | None = 200,
    pose_timeout: int | None = 600,
    limit_threads: bool = True,
    monitor: bool = False,
    monitor_interval: int = 30,
    stall_timeout: int | None = None,
    max_restarts: int = 5,
    require_complete: bool = False,
) -> pd.DataFrame:
    """Run PoseBusters on *poses_list* using multiprocessing.Pool.

    *worker_maxtasks* recycles each worker after that many poses (bounds RDKit
    memory growth); *pose_timeout* caps a single bust() call; *limit_threads*
    keeps each worker single-threaded so N workers don't oversubscribe the CPU.

    Stall recovery: if no pose completes for *stall_timeout* seconds the pool is
    assumed wedged (a dead / OOM-killed worker makes imap_unordered block
    forever), so it is terminated and restarted on the remaining poses — halving
    the worker count each time to relieve memory pressure — up to *max_restarts*
    times. *monitor* starts a background resource sampler for the run.
    """
    if not poses_list:
        # An empty current input universe must not resurrect an unrelated CSV
        # left by an older method/variant selection.
        print("No poses to analyze!")
        return pd.DataFrame()

    previous_results = pd.DataFrame()
    already_done: set[str] = set()
    if output_file is not None and output_file.exists():
        if overwrite:
            print(f"  OVERWRITE — deleting {output_file.name}")
            output_file.unlink()
        else:
            previous_results = pd.read_csv(output_file)
            previous_row_count = len(previous_results)
            previous_results, already_done = _reconcile_previous_results(
                previous_results, poses_list)
            if len(previous_results) != previous_row_count:
                # Publish the reconciled universe immediately. If every pending
                # pose later errors, an old stale CSV must not survive and look
                # like the result of the current inputs.
                temporary = output_file.with_suffix(output_file.suffix + ".tmp")
                previous_results.to_csv(temporary, index=False)
                temporary.replace(output_file)
            print(f"  Loaded {len(previous_results)} previous results; {len(already_done)} poses done")

    pending = [
        p for p in poses_list
        if str(p["validation_signature"]) not in already_done
    ]
    n_skipped = len(poses_list) - len(pending)
    if n_skipped:
        print(f"  Skipping {n_skipped} already-inspected poses")
    if not pending:
        print("  All poses already inspected.")
        if output_file is not None:
            output_file.with_name("posebusters_incomplete.txt").unlink(missing_ok=True)
        return previous_results

    # num_workers None/0/negative all mean "use every core" (0 is the common
    # 'auto' spelling; before, 0 silently ran the whole job on ONE worker).
    n_workers = num_workers if (num_workers and num_workers > 0) else cpu_count()
    n_workers = min(n_workers, len(pending))
    print(f"\n  Workers: {n_workers} ({cpu_count()} CPUs available)")

    new_frames: list[pd.DataFrame] = []
    timeout_msgs: list[str] = []
    failure_msgs: list[str] = []
    total = len(pending)
    consumed: set[str] = set()          # validation signatures consumed across passes
    progress = {"n": 0}                  # live counter the monitor reads
    n_new = n_errors = n_timeout = n_poison = n_dock = 0
    print(f"  {total} poses to process (checkpoint every {save_interval})...\n")

    def _checkpoint(final: bool = False) -> None:
        if not output_file:
            return
        frames = ([previous_results] if not previous_results.empty else []) + new_frames
        if frames:
            _save_checkpoint(frames, output_file, n_new, total, final=final)
            # Collapse the growing list of 1-row frames into one so the next
            # checkpoint and the final concat stay O(N) instead of O(N^2) — a
            # thousands-of-tiny-frames concat every interval is a real slowdown.
            if not final and len(new_frames) > 1:
                new_frames[:] = [pd.concat(new_frames, ignore_index=True)]

    maxtasks = worker_maxtasks if (worker_maxtasks and worker_maxtasks > 0) else None
    if maxtasks:
        print(f"  Recycling each worker every {maxtasks} poses"
              f"{f'; per-pose timeout {pose_timeout}s' if pose_timeout else ''}")

    # No-result-for-this-long => the pool is wedged (dead/OOM-killed worker makes
    # imap_unordered block forever). Default to comfortably longer than a single
    # legit slow pose so only a real stall trips it. 0 disables detection.
    if stall_timeout is None:
        stall_secs = max(pose_timeout or 0, 300) * 3
    elif stall_timeout <= 0:
        stall_secs = 0
    else:
        stall_secs = stall_timeout
    if stall_secs:
        print(f"  Stall watchdog: restart the pool if no pose completes for {stall_secs}s "
              f"(up to {max_restarts} restarts)")

    # Shared poison-ligand registry so a molecule that hangs bust() is validated
    # at most once, not once per (up to 540) sibling pose. Best-effort: if a
    # Manager can't start, the skip is simply disabled.
    poison_mgr = poison = None
    if pose_timeout and pose_timeout > 0:
        try:
            poison_mgr = get_context("fork").Manager()
            poison = poison_mgr.dict()
        except Exception as e:
            print(f"  (poison-ligand skip unavailable: {e})")
            poison_mgr = poison = None

    def _handle(task_id, result_df, error_msg) -> None:
        nonlocal n_new, n_errors, n_timeout, n_poison, n_dock
        consumed.add(str(task_id))
        progress["n"] = len(consumed)
        if error_msg:
            es = str(error_msg)
            failure_msgs.append(es)
            if es.startswith("TIMEOUT"):
                n_timeout += 1
                timeout_msgs.append(es)
                print(f"  {es}")
            elif es.startswith("POISON-SKIP"):
                n_poison += 1          # many per poison ligand — tally silently
            else:
                n_errors += 1
                print(f"  {es}")
            return
        if result_df is not None:
            mode_val = result_df["posebusters_mode"].iloc[0] if "posebusters_mode" in result_df.columns else ""
            if mode_val == "dock":
                n_dock += 1
            new_frames.append(result_df)
            n_new += 1
            if n_new % 10 == 0 or n_new == 1:
                print(f"  Processed {n_new}/{total}  (overall {len(already_done) + n_new}/{len(poses_list)})")
            if n_new % save_interval == 0:
                _checkpoint()

    def _run_one_pass(work: list[dict]) -> bool:
        """Run one Pool over *work*; return True if it stalled and needs a restart."""
        procs = max(1, min(n_workers, len(work)))
        stalled = False
        mon = (_ResourceMonitor(monitor_interval, cpu_count(), lambda: progress["n"])
               if monitor else None)
        pool = Pool(processes=procs, maxtasksperchild=maxtasks,
                    initializer=_init_worker,
                    initargs=(config, dict(protein_cache), quiet, pose_timeout,
                              limit_threads, poison))
        if mon:
            mon.start()
        try:
            it = pool.imap_unordered(_process_single_pose, work)
            seen = 0
            while seen < len(work):
                try:
                    task_id, result_df, error_msg = (
                        it.next(stall_secs) if stall_secs else next(it))
                except StopIteration:
                    break
                except MPTimeoutError:
                    stalled = True
                    print(f"\n  [STALL] no pose completed in {stall_secs}s — a worker likely "
                          f"died (OOM?) or is wedged in an uninterruptible call. "
                          f"Terminating the pool and restarting on the remainder.")
                    break
                seen += 1
                _handle(task_id, result_df, error_msg)
        finally:
            if stalled:
                pool.terminate()
            else:
                pool.close()
            pool.join()
            if mon:
                mon.stop()
        return stalled

    remaining = list(pending)
    restart = 0
    while remaining:
        stalled = _run_one_pass(remaining)
        remaining = [
            p for p in pending
            if str(p["validation_signature"]) not in consumed
        ]
        if not stalled:
            break
        restart += 1
        if restart > max_restarts:
            print(f"  [STALL] gave up after {max_restarts} restarts; {len(remaining)} pose(s) "
                  f"left unprocessed — they are not in the CSV, so the next run resumes them.")
            break
        old = n_workers
        n_workers = max(1, n_workers // 2)
        _checkpoint()  # persist what we have before retrying
        print(f"  [STALL] restart {restart}/{max_restarts} on {len(remaining)} remaining pose(s)"
              + (f"; workers {old}->{n_workers} to ease memory pressure" if n_workers != old else ""))

    _checkpoint(final=True)

    poison_ligands = sorted(poison.keys()) if poison is not None else []
    if poison_mgr is not None:
        try:
            poison_mgr.shutdown()
        except Exception:
            pass

    print(f"\n  Done:  new={n_new}  prev={n_skipped}  total={len(already_done)+n_new}"
          f"  errors={n_errors}  timeouts={n_timeout}"
          + (f"  poison_skipped={n_poison}" if n_poison else "")
          + (f"  stalls={restart}" if restart else ""))
    if config == "dock":
        print(f"    dock mode: {n_dock}  |  ligand-only fallback: disabled")
    if poison_ligands:
        # A molecule whose conformer generation hangs bust(); every pose of it was
        # dropped after the first timeout. Log so they can be inspected/excluded.
        print(f"\n  {len(poison_ligands)} ligand(s) hung bust() (conformer generation) and had "
              f"all their poses skipped after the first {pose_timeout}s timeout:")
        for lk in poison_ligands:
            print(f"    ✗ {lk}")
        if output_file is not None:
            pf = output_file.with_name("posebusters_poison_ligands.txt")
            pf.write_text("\n".join(poison_ligands) + "\n")
            print(f"  Poison ligands logged to: {pf}")
    if timeout_msgs:
        # Timed-out poses were skipped (not in the CSV) so a plain resume RETRIES
        # them next run. Record them so a permanently-pathological pose can be
        # inspected or excluded instead of retried forever.
        print(f"\n  {n_timeout} pose(s) exceeded the {pose_timeout}s cap and were "
              f"skipped — NOT written to the CSV, so the next run will retry them.")
        if output_file is not None:
            timeout_file = output_file.with_name("posebusters_timeouts.txt")
            existing = (timeout_file.read_text().splitlines()
                        if timeout_file.exists() else [])
            merged = sorted(set(existing) | set(timeout_msgs))
            timeout_file.write_text("\n".join(merged) + "\n")
            print(f"  Timed-out poses logged to: {timeout_file}")

    incomplete = len(failure_msgs) + len(remaining)
    if incomplete and output_file is not None:
        incomplete_file = output_file.with_name("posebusters_incomplete.txt")
        details = list(failure_msgs)
        details.extend(
            f"UNPROCESSED: {pose.get('method')} {pose.get('pose_file')}"
            for pose in remaining
        )
        incomplete_file.write_text("\n".join(details) + "\n")
        print(f"  Incomplete validation details: {incomplete_file}")
    if incomplete and require_complete:
        raise RuntimeError(
            f"PoseBusters validation incomplete: {len(failure_msgs)} failed/timed-out "
            f"and {len(remaining)} unprocessed pose(s). Partial results were checkpointed; "
            "fix the reported inputs and resume before using benchmark verdicts."
        )
    if not incomplete and output_file is not None:
        output_file.with_name("posebusters_incomplete.txt").unlink(missing_ok=True)

    all_frames = ([previous_results] if not previous_results.empty else []) + new_frames
    return pd.concat(all_frames, ignore_index=True) if all_frames else pd.DataFrame()


# ============================================================================
# STAGE: pose discovery + overview
# ============================================================================

def collect_pose_rows(ctx: PipelineConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collect rows from every configured docking method.

    Returns (filtered_poses_df, overview_df).
    """
    print("=" * 100)
    print("DOCKING POSES OVERVIEW")
    print("=" * 100)
    print("\nCollecting pose counts from all docking methods...\n")

    all_rows: list[dict] = []
    method_dfs: dict[str, pd.DataFrame] = {}
    collected_methods: set[str] = set()

    for method_key, poses_dir in ctx.docking_directories.items():
        collector = ROW_COLLECTORS.get(method_key)
        if collector is None:
            if ctx.require_all_methods:
                raise KeyError(f"No row collector for required method {method_key!r}")
            print(f"  WARNING: No row collector for '{method_key}', skipping.")
            continue
        if not poses_dir.exists():
            if ctx.require_all_methods:
                raise FileNotFoundError(
                    f"Required docking directory not found for {method_key}: {poses_dir}"
                )
            print(f"  WARNING: Directory not found for '{method_key}': {poses_dir}, skipping.")
            continue

        # variant_filter (if set) restricts pose_count to the kept variant(s) so
        # the overview, common-combo intersection and per-complex resume all agree
        # with the busted set (see _apply_variant_filter for the pose-level pass).
        keep_variants = (ctx.variant_filter.get(_variant_base_method(method_key))
                         if ctx.variant_filter else None)
        method_rows = collector(poses_dir, keep_variants)
        for r in method_rows:
            r["docking_tool"] = method_key
        all_rows.extend(method_rows)

        df_m = pd.DataFrame(method_rows)
        if not df_m.empty:
            collected_methods.add(method_key)
            col_name = method_key.replace(" ", "_").title()
            summary = df_m.groupby(["protein", "ligand"])["pose_count"].sum().reset_index()
            summary.columns = ["Protein", "Ligand", col_name]
            method_dfs[col_name] = summary
        print(f"  {method_key}: {len(method_rows)} protein-ligand combinations")

    if ctx.require_all_methods:
        missing = sorted(set(ctx.docking_directories) - collected_methods)
        if missing:
            raise RuntimeError(
                "Required docking method(s) produced no pose rows: " + ", ".join(missing)
            )

    df_combined = pd.DataFrame()
    for col_name, df_m in method_dfs.items():
        df_combined = df_m if df_combined.empty else df_combined.merge(
            df_m, on=["Protein", "Ligand"], how="outer"
        )

    if not df_combined.empty:
        method_cols = [c for c in df_combined.columns if c not in ("Protein", "Ligand")]
        df_combined = df_combined.fillna(0)
        for col in method_cols:
            df_combined[col] = df_combined[col].astype(int)
        df_combined["Total"] = df_combined[method_cols].sum(axis=1)
        df_combined = df_combined.sort_values(["Protein", "Ligand"])

        # The full per-complex table (one row per protein-ligand) is hundreds of
        # lines on the benchmark set; under --quiet keep only the summary totals
        # below (the whole table is still written to pose_counts_overview.csv).
        if not _QUIET:
            print("\n" + "=" * 100)
            print("POSE COUNTS BY PROTEIN-LIGAND COMBINATION")
            print("=" * 100)
            print(df_combined.to_string(index=False))

        print("\n" + "=" * 100)
        print("SUMMARY STATISTICS")
        print("=" * 100)
        print(f"\nTotal unique protein-ligand combinations: {len(df_combined)}")
        print("\nTotal poses by method:")
        for col in method_cols:
            print(f"  {col}: {df_combined[col].sum():,} poses")
        print(f"  Combined Total: {df_combined['Total'].sum():,} poses")
        print(f"\nUnique proteins: {df_combined['Protein'].nunique()}")
        print(f"Unique ligands: {df_combined['Ligand'].nunique()}")

        if not _QUIET:
            print("\n" + "-" * 100 + "\nPOSES BY LIGAND\n" + "-" * 100)
            print(df_combined.groupby("Ligand")[method_cols + ["Total"]].sum().to_string())
            print("\n" + "-" * 100 + "\nPOSES BY PROTEIN\n" + "-" * 100)
            print(df_combined.groupby("Protein")[method_cols + ["Total"]].sum().to_string())

        out_csv = ctx.output_base_dir / "pose_counts_overview.csv"
        df_combined.to_csv(out_csv, index=False)
        print(f"\n{'=' * 100}\nTable saved to: {out_csv}\n{'=' * 100}")
    else:
        print("\nNo docking results found!")

    return pd.DataFrame(all_rows), df_combined


def filter_common_combos(
    filtered_poses_df: pd.DataFrame,
    optional_methods: set[str] | frozenset[str] = frozenset(),
) -> pd.DataFrame:
    """Keep only protein-ligand combinations present across the REQUIRED methods.

    The shared cohort is the intersection over the required (non-``optional``)
    methods; ``optional_methods`` are kept for whichever cohort combos they
    happen to cover but never shrink the cohort, so a partially docked or
    brand-new tool cannot collapse the benchmark down to its own subset. With no
    optional methods this is the historical "present in ALL methods" filter."""
    print(f"Total entries: {len(filtered_poses_df)}")
    if filtered_poses_df.empty:
        print("filtered_poses_df is empty. Check docking_directories paths.")
        return filtered_poses_df

    optional = {str(m).strip() for m in (optional_methods or set())}
    print(f"Methods found: {filtered_poses_df['docking_tool'].unique().tolist()}")
    combos_by_method = {
        m: set(zip(g["protein"], g["ligand"]))
        for m, g in filtered_poses_df.groupby("docking_tool")
    }
    # The cohort is defined by the required methods only; optional methods that
    # produced rows are still busted, but they don't constrain the intersection.
    required_by_method = {m: c for m, c in combos_by_method.items() if m not in optional}
    basis = required_by_method or combos_by_method   # all-optional → fall back to all

    print("=" * 100)
    print("FILTERING: Keep Only Protein-Ligand Combinations in ALL REQUIRED Methods")
    print("=" * 100)
    for method, combos in combos_by_method.items():
        tag = "  (optional — does not constrain the cohort)" if method in optional else ""
        print(f"  {method}: {len(combos)} combinations{tag}")

    common_combos = set.intersection(*basis.values()) if basis else set()
    print(f"\nProtein-ligand combinations in ALL required methods: {len(common_combos)}")

    if not common_combos:
        print("\nWARNING: No combinations found in all methods!")
        return filtered_poses_df.iloc[0:0]

    if not _QUIET:
        for protein, ligand in sorted(common_combos):
            print(f"  {protein} + {ligand}")

    combo_index = pd.MultiIndex.from_arrays(
        [filtered_poses_df["protein"], filtered_poses_df["ligand"]]
    )
    filtered_poses_df = filtered_poses_df[combo_index.isin(common_combos)].copy()

    print(f"\n{'-' * 100}\nFILTERED SUBSET\n{'-' * 100}")
    per_method = filtered_poses_df.groupby("docking_tool")["pose_count"].agg(["count", "sum"])
    for method, (n_combos, n_poses) in per_method.iterrows():
        print(f"  {method}: {n_combos} combinations, {n_poses:,} poses")

    methods_in_df = sorted(filtered_poses_df["docking_tool"].unique())
    pivot = (
        filtered_poses_df
        .pivot_table(index=["protein", "ligand"], columns="docking_tool",
                     values="pose_count", aggfunc="sum", fill_value=0)
        .reindex(columns=methods_in_df, fill_value=0)
        .sort_index()
    )
    pivot["__total__"] = pivot.sum(axis=1)

    header = f"{'Protein':<35} {'Ligand':<30}" + "".join(f" {m:>15}" for m in methods_in_df) + f" {'Total':>10}"
    if not _QUIET:                     # per-complex rows: hundreds of lines
        print(f"\n{header}\n" + "-" * len(header))
        for (protein, ligand), row in pivot.iterrows():
            cells = "".join(f" {int(row[m]):>15,}" for m in methods_in_df)
            print(f"{protein:<35} {ligand:<30}{cells} {int(row['__total__']):>10,}")
        print("-" * len(header))
    grand_totals = pivot[methods_in_df].sum().astype(int)
    grand_all = int(pivot["__total__"].sum())
    totals_str = f"{'TOTAL':<35} {'':<30}" + "".join(f" {grand_totals[m]:>15,}" for m in methods_in_df) + f" {grand_all:>10,}"
    print(totals_str)
    print(f"\nFinal filtered_poses_df shape: {filtered_poses_df.shape}")
    return filtered_poses_df


def filter_already_processed_complexes(
    filtered_poses_df: pd.DataFrame,
    results_csv_path: Path,
    overwrite: bool = False,
) -> pd.DataFrame:
    """Drop complex rows already fully validated in a previous run.

    A "complex" is a (docking_tool, protein, ligand) triple. The per-pose resume
    in ``analyze_poses_with_posebusters`` only skips the (cheap) ``bust()`` call
    for poses already in the results CSV — the costly pose-file collection and
    PDBQT->SDF conversion in ``collect_all_pose_files`` still runs for every
    complex. Removing complexes that are already complete here lets the pipeline
    skip that work entirely so only never-validated complexes are processed.

    A complex counts as done when the number of its poses already present in
    *results_csv_path* is >= the ``pose_count`` discovered for it, so a complex
    that was only partially validated (e.g. an interrupted run) is still
    reprocessed. With *overwrite* the CSV is about to be rebuilt, so nothing is
    skipped.
    """
    if overwrite or filtered_poses_df.empty or not results_csv_path.exists():
        return filtered_poses_df

    key_cols = ["docking_method", "protein", "ligand"]
    try:
        prev = pd.read_csv(results_csv_path, usecols=lambda c: c in set(key_cols))
    except (ValueError, pd.errors.EmptyDataError):
        return filtered_poses_df  # no key columns / empty file -> nothing to skip
    if prev.empty or not set(key_cols) <= set(prev.columns):
        return filtered_poses_df

    # done_counts maps (docking_tool, protein, ligand) -> #poses already in CSV.
    done_counts = (
        prev.dropna(subset=key_cols)
            .astype({c: str for c in key_cols})
            .groupby(key_cols).size().to_dict()
    )

    # A collector may emit several rows for one complex (AutoDock raw PDBQT plus
    # one committed SDF row per optimized pose). Compare the previous result
    # count with the SUM across that whole tool/complex, then keep/drop all of its
    # rows together. Per-row comparisons would let two already-busted raw poses
    # incorrectly mark every one-pose optimized row as complete.
    expected_counts = (
        filtered_poses_df.assign(
            docking_tool=filtered_poses_df["docking_tool"].astype(str),
            protein=filtered_poses_df["protein"].astype(str),
            ligand=filtered_poses_df["ligand"].astype(str),
        )
        .groupby(["docking_tool", "protein", "ligand"])["pose_count"]
        .sum().to_dict()
    )
    complete_keys = {
        key for key, expected in expected_counts.items()
        if done_counts.get(key, 0) >= expected
    }
    row_keys = list(zip(
        filtered_poses_df["docking_tool"].astype(str),
        filtered_poses_df["protein"].astype(str),
        filtered_poses_df["ligand"].astype(str),
    ))
    done_mask = pd.Series([key in complete_keys for key in row_keys],
                          index=filtered_poses_df.index)
    n_done = int(done_mask.sum())

    print("\n" + "=" * 80)
    print("SKIPPING ALREADY-VALIDATED COMPLEXES")
    print("=" * 80)
    if n_done:
        done_keys = sorted(
            {(t, p, l) for t, p, l in zip(
                filtered_poses_df.loc[done_mask, "docking_tool"],
                filtered_poses_df.loc[done_mask, "protein"],
                filtered_poses_df.loc[done_mask, "ligand"])}
        )
        for tool, protein, ligand in done_keys:
            print(f"  ✓ done: {tool:>22}  {protein} + {ligand}")
    print(f"\n  {n_done}/{len(filtered_poses_df)} complex rows already complete "
          f"→ {len(filtered_poses_df) - n_done} remaining to process")
    return filtered_poses_df.loc[~done_mask].copy()


def resolve_protein_files(
    filtered_poses_df: pd.DataFrame,
    file_map: dict[str, str],
    normalized_map: dict[str, str],
) -> dict[str, str | None]:
    """Resolve a PDB file for every unique protein in the filtered set."""
    cache: dict[str, str | None] = {}
    print("-" * 80 + "\nPROTEIN FILE RESOLUTION\n" + "-" * 80)
    if filtered_poses_df.empty:
        return cache

    for pname in sorted(filtered_poses_df["protein"].unique()):
        pfile = find_protein_file(pname, file_map, normalized_map)
        cache[pname] = pfile
        status = f"✓ {Path(pfile).name}" if pfile else "✗ NOT FOUND"
        print(f"  {pname}: {status}")

    n_found = sum(1 for v in cache.values() if v)
    n_missing = sum(1 for v in cache.values() if not v)
    print(f"\n  Resolved: {n_found}/{len(cache)} proteins")
    if n_missing:
        print(f"  WARNING: {n_missing} proteins have no PDB — dock-mode preflight will reject them")
    return cache


# ============================================================================
# STAGE: copy proved poses
# ============================================================================

def _infer_method(filepath: str, docking_directories: dict[str, Path]) -> str:
    fp = str(filepath)
    for mk, dp in docking_directories.items():
        dp_str = str(dp)
        if dp_str in fp or Path(dp_str).name in fp:
            return mk
    fl = fp.lower()
    if "unidock2" in fl or "uni-dock2" in fl:
        return "unidock2"
    if "unidock" in fl or "uni-dock" in fl:
        return next((k for k in docking_directories if k.startswith("unidock")), "unidock")
    if any(x in fl for x in ("autodock", "vina", "_vina_", "converted_pdbqt")):
        return "autodock"
    if any(x in fl for x in ("diffdock", "diff_dock")):
        return "diffdock"
    if any(x in fl for x in ("equibind", "equi_bind")):
        if "exclusion" in fl or "spatial" in fl:
            return next((k for k in docking_directories if "exclusion" in k), "equibind_exclusion")
        return next((k for k in docking_directories if "guided" in k), "equibind_guided")
    return "unknown"


def copy_proved_poses(ctx: PipelineConfig) -> None:
    """Copy poses that passed all PoseBusters tests into organised folders."""
    results_csv = ctx.output_dir / "posebusters_filtered_results.csv"
    if not results_csv.exists():
        csvs = list(ctx.output_dir.glob("*.csv"))
        for c in csvs:
            if "result" in c.name.lower() or "posebust" in c.name.lower():
                results_csv = c
                break
        if not results_csv.exists() and csvs:
            results_csv = csvs[0]
        if not results_csv.exists():
            raise FileNotFoundError(f"No PoseBusters results CSV in {ctx.output_dir}")

    print("=" * 100)
    print(f"Loading PoseBusters results from: {results_csv}")
    print("=" * 100)

    df_pb = pd.read_csv(results_csv)
    print(f"\nTotal entries: {len(df_pb)}")

    test_cols = require_test_columns(df_pb, ctx.config_mode)
    coerce_test_cols_to_bool(df_pb, test_cols)

    # Scope copied artifacts to the exact results file. This prevents stale
    # passes and the historical ``_dupN`` accumulation from contaminating a new
    # benchmark run while keeping each generated set reproducible.
    result_key = _sha256_file(results_csv)[:16]
    output_base = ctx.output_dir / "proved_poses" / result_key
    folders = {key: output_base / key for key in ctx.docking_directories}

    print(f"\nIdentified {len(test_cols)} PoseBusters test columns:")
    for tc in test_cols:
        n_pass = df_pb[tc].sum()
        print(f"  {tc}: {n_pass}/{len(df_pb)} passed ({100*n_pass/len(df_pb):.1f}%)")

    excluded_found = [
        c for c in df_pb.columns
        if c in _EXCLUDE_COLS or c.lower() in _EXCLUDE_COLS
        or c.lower().startswith(("number_", "num_"))
    ]
    if excluded_found:
        print(f"\nExcluded {len(excluded_found)} non-test columns:")
        for ec in excluded_found:
            print(f"  ✗ {ec}")

    df_pb["all_passed"] = df_pb[test_cols].all(axis=1)
    df_passed = df_pb[df_pb["all_passed"]].copy()

    print(f"\n{'=' * 100}")
    print(f"PASSED ALL {len(test_cols)} TESTS: {len(df_passed)} / {len(df_pb)} "
          f"({100*len(df_passed)/len(df_pb):.1f}%)")
    print(f"{'=' * 100}")

    print("\nBottleneck tests (lowest pass rates):")
    for tc, rate in sorted(
        ((tc, df_pb[tc].mean() * 100) for tc in test_cols), key=lambda x: x[1]
    ):
        if rate < 99.0:
            print(f"  {tc}: {rate:.1f}%")

    file_col = next(
        (c for c in ("file_path", "filepath", "sdf_file", "sdf_path", "path",
                     "file", "mol_pred", "File", "FILE_PATH", "pose_file")
         if c in df_passed.columns), None,
    )
    if file_col is None:
        file_col = next(
            (c for c in df_passed.columns if "path" in c.lower() or "file" in c.lower()),
            None,
        )
    if file_col is None:
        raise ValueError("Cannot find file path column in results.")

    method_col = next(
        (c for c in ("method", "docking_method", "Method", "DOCKING_METHOD")
         if c in df_passed.columns), None,
    )
    print(f"\nFile column: '{file_col}'")
    print(f"Method column: '{method_col}'" if method_col else "Method column: NOT FOUND (inferring)")

    df_passed["_method"] = (
        df_passed[method_col]
        if method_col
        else df_passed[file_col].apply(lambda fp: _infer_method(fp, ctx.docking_directories))
    )
    print("\nPassed poses by method:")
    for method, count in df_passed["_method"].value_counts().items():
        print(f"  {method}: {count}")

    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=True)

    copied_count = {key: 0 for key in ctx.docking_directories}
    copied_count["unknown"] = 0
    skipped_count = {"not_found": 0, "unknown_method": 0}

    print(f"\n{'=' * 100}\nCOPYING POSEBUSTERS-PROVED POSES...\n{'=' * 100}")

    for _, row in df_passed.iterrows():
        src_path = Path(str(row[file_col]))
        method = row["_method"]
        protein_name = str(row["protein"]) if "protein" in row.index and pd.notna(row.get("protein")) else ""
        ligand_name = str(row["ligand"]) if "ligand" in row.index and pd.notna(row.get("ligand")) else ""

        if not src_path.is_absolute():
            src_path = ctx.work_dir / src_path
        if not src_path.exists():
            for alt in (ctx.output_dir / src_path.name, ctx.converted_dir / src_path.name):
                if alt.exists():
                    src_path = alt
                    break
            else:
                skipped_count["not_found"] += 1
                continue

        if method not in folders:
            skipped_count["unknown_method"] += 1
            print(f"  WARNING: Unknown method '{method}' for {src_path.name}")
            continue

        dest_dir = folders[method]
        pose_key = str(row.get("validation_signature") or "").strip()[:16]
        if not pose_key:
            pose_key = hashlib.sha256(str(src_path.resolve()).encode()).hexdigest()[:16]
        prefix = f"{protein_name}__{ligand_name}__" if protein_name and ligand_name else ""
        dest_name = f"{prefix}{pose_key}__{src_path.name}"
        dest_file = dest_dir / dest_name

        shutil.copy2(str(src_path), str(dest_file))
        copied_count[method] = copied_count.get(method, 0) + 1

    print(f"\n{'=' * 100}\nCOPY COMPLETE\n{'=' * 100}")
    print(f"\n  Output: {output_base}")
    total_copied = 0
    for mk, cnt in copied_count.items():
        if cnt > 0:
            print(f"    {mk:>25}: {cnt:>4} files  →  {folders.get(mk, 'N/A')}")
            total_copied += cnt
    print(f"\n  Total copied: {total_copied}")
    if skipped_count["not_found"]:
        print(f"  Skipped (not found): {skipped_count['not_found']}")
    if skipped_count["unknown_method"]:
        print(f"  Skipped (unknown method): {skipped_count['unknown_method']}")

    for mn, fp in folders.items():
        print(f"    {fp.relative_to(ctx.work_dir)}: {len(list(fp.glob('*')))} files")

    manifest_cols = [file_col, "_method"]
    for extra in ("protein", "ligand"):
        if extra in df_passed.columns:
            manifest_cols.append(extra)
    manifest_cols += test_cols
    manifest = df_passed[manifest_cols].copy()
    manifest.rename(columns={"_method": "docking_method"}, inplace=True)
    manifest_path = output_base / "posebuster_proved_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"\n  Manifest: {manifest_path}\n{'=' * 100}")


# ============================================================================
# STAGE: plots
# ============================================================================

TEST_DISPLAY = {
    "mol_pred_loaded": "Molecule Loaded",
    "sanitization": "Sanitization",
    "inchi_convertible": "InChI Convertible",
    "all_atoms_connected": "All Atoms Connected",
    "no_radicals": "No Radicals",
    "bond_lengths": "Bond Lengths",
    "bond_angles": "Bond Angles",
    "internal_steric_clash": "No Steric Clash",
    "aromatic_ring_flatness": "Aromatic Flatness",
    "non-aromatic_ring_non-flatness": "Non-Arom. Ring Shape",
    "double_bond_flatness": "Double Bond Flatness",
    "internal_energy": "Internal Energy",
    "passes_valence_checks": "Valence Checks",
    "passes_kekulization": "Kekulization",
    "no_radicals_before_sanitization": "No Pre-Sanit. Radicals",
}

_COLOR_PALETTE = [
    "#3498db", "#e74c3c", "#2ecc71", "#9b59b6", "#f39c12",
    "#1abc9c", "#e67e22", "#34495e", "#d35400", "#8e44ad",
]


def make_plots(ctx: PipelineConfig) -> None:
    """Generate all summary figures from the PoseBusters results CSV."""
    pb_csv = ctx.output_dir / "posebusters_filtered_results.csv"
    plot_output_dir = ctx.output_dir
    df = pd.read_csv(pb_csv)

    test_cols = require_test_columns(df, ctx.config_mode)
    coerce_test_cols_to_bool(df, test_cols)
    df["all_passed"] = df[test_cols].all(axis=1)

    methods = sorted(df["docking_method"].unique())
    method_colors = {m: _COLOR_PALETTE[i % len(_COLOR_PALETTE)] for i, m in enumerate(methods)}

    n_total = len(df)
    n_passed = int(df["all_passed"].sum())
    print(f"PoseBusters Results: {n_total} poses, {n_passed} passed all ({n_passed/n_total*100:.1f}%)")
    for m in methods:
        g = df[df["docking_method"] == m]
        print(f"  {m}: {g['all_passed'].sum()}/{len(g)} ({g['all_passed'].mean()*100:.1f}%)")

    # FIGURE 1: heatmap
    variable_tests = [t for t in test_cols if df[t].mean() < 1.0]
    trivial_tests = [t for t in test_cols if t not in variable_tests]
    pass_rates = df.groupby("docking_method")[test_cols].mean() * 100
    pass_rates_var = pass_rates[variable_tests] if variable_tests else pass_rates

    fig1, ax1 = plt.subplots(figsize=(max(10, len(test_cols) * 0.8), max(4, len(methods) * 0.8)))
    cmap = LinearSegmentedColormap.from_list("ryg", ["#e74c3c", "#f39c12", "#27ae60"])
    im = ax1.imshow(pass_rates_var.values, cmap=cmap, aspect="auto", vmin=50, vmax=100)
    ax1.set_yticks(range(len(pass_rates_var.index)))
    ax1.set_yticklabels(pass_rates_var.index, fontsize=11, fontweight="bold")
    ax1.set_xticks(range(len(pass_rates_var.columns)))
    ax1.set_xticklabels(
        [TEST_DISPLAY.get(c, c.replace("_", " ").title()) for c in pass_rates_var.columns],
        rotation=45, ha="right", fontsize=10,
    )
    for i in range(pass_rates_var.shape[0]):
        for j in range(pass_rates_var.shape[1]):
            val = pass_rates_var.values[i, j]
            ax1.text(j, i, f"{val:.1f}%", ha="center", va="center", fontsize=10,
                     fontweight="bold", color="white" if val < 75 else "black")
    cbar = plt.colorbar(im, ax=ax1, shrink=0.8, pad=0.02)
    cbar.set_label("Pass Rate (%)", fontsize=11)
    if trivial_tests:
        ax1.set_xlabel(
            f"Tests at 100% (not shown): {', '.join(TEST_DISPLAY.get(t, t) for t in trivial_tests)}",
            fontsize=8, style="italic",
        )
    ax1.set_title("PoseBusters Test Pass Rates by Docking Method", fontsize=14, fontweight="bold", pad=12)
    fig1.tight_layout()
    fig1.savefig(plot_output_dir / "pb_passrate_heatmap.png", dpi=200, bbox_inches="tight")
    plt.close(fig1)

    # FIGURE 2: stacked bars
    fig2, ax2 = plt.subplots(figsize=(max(7, len(methods) * 2), 5))
    counts = df.groupby("docking_method")["all_passed"].value_counts().unstack(fill_value=0)
    counts = counts.reindex(columns=[True, False], fill_value=0).rename(
        columns={True: "Passed All", False: "Failed >=1"}
    )

    ax2.bar(range(len(counts)), counts["Passed All"],
            color=[method_colors.get(m, "#888") for m in counts.index],
            edgecolor="white", linewidth=1.5, label="Passed All Tests")
    ax2.bar(range(len(counts)), counts["Failed >=1"], bottom=counts["Passed All"],
            color=[method_colors.get(m, "#888") for m in counts.index],
            alpha=0.3, edgecolor="white", linewidth=1.5, hatch="///", label="Failed >=1 Test")
    for i, (m, row) in enumerate(counts.iterrows()):
        total = row.sum()
        pct = row["Passed All"] / total * 100
        ax2.text(i, total + 5, f'{int(row["Passed All"])}/{int(total)}\n({pct:.1f}%)',
                 ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax2.set_xticks(range(len(counts)))
    ax2.set_xticklabels(counts.index, fontsize=12, fontweight="bold")
    ax2.set_ylabel("Number of Poses", fontsize=12)
    ax2.set_title("PoseBusters Validation Summary by Docking Method", fontsize=14, fontweight="bold")
    ax2.legend(loc="upper right", fontsize=10)
    ax2.set_ylim(0, counts.sum(axis=1).max() * 1.25)
    ax2.grid(axis="y", alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(plot_output_dir / "pb_pass_fail_bars.png", dpi=200, bbox_inches="tight")
    plt.close(fig2)

    # FIGURE 3: per-test failure rates
    if variable_tests:
        fail_rates = (1 - df.groupby("docking_method")[variable_tests].mean()) * 100
        fig3, ax3 = plt.subplots(figsize=(max(10, len(variable_tests) * 1.5), 5))
        x = np.arange(len(variable_tests))
        width = 0.8 / len(methods)
        for i, m in enumerate(methods):
            vals = fail_rates.loc[m].values
            ax3.bar(x + i * width - 0.4 + width / 2, vals, width,
                    label=m, color=method_colors.get(m, "#888"), edgecolor="white")
            for j, v in enumerate(vals):
                if v > 2:
                    ax3.text(x[j] + i * width - 0.4 + width / 2, v + 0.5,
                             f"{v:.1f}%", ha="center", va="bottom", fontsize=8, fontweight="bold")
        ax3.set_xticks(x)
        ax3.set_xticklabels(
            [TEST_DISPLAY.get(c, c.replace("_", " ").title()) for c in variable_tests],
            rotation=40, ha="right", fontsize=10,
        )
        ax3.set_ylabel("Failure Rate (%)", fontsize=12)
        ax3.set_title("PoseBusters Failure Rates by Test and Method", fontsize=14, fontweight="bold")
        ax3.legend(fontsize=10)
        ax3.grid(axis="y", alpha=0.3)
        fig3.tight_layout()
        fig3.savefig(plot_output_dir / "pb_failure_rates.png", dpi=200, bbox_inches="tight")
        plt.close(fig3)

    # FIGURE 4: per protein-ligand pass rate
    n_methods = len(methods)
    fig4, axes4 = plt.subplots(1, n_methods, figsize=(6 * n_methods, 5), sharey=True, squeeze=False)
    for ax, m in zip(axes4.flatten(), methods):
        sub = df[df["docking_method"] == m]
        combo_pass = sub.groupby(["protein", "ligand"])["all_passed"].mean() * 100
        combo_pass = combo_pass.sort_values(ascending=True)
        labels = [f"{p}\n{l}" for p, l in combo_pass.index]
        colors = [("#27ae60" if v >= 75 else "#f39c12" if v >= 50 else "#e74c3c") for v in combo_pass.values]
        ax.barh(range(len(combo_pass)), combo_pass.values, color=colors, edgecolor="white")
        ax.set_yticks(range(len(combo_pass)))
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Pass Rate (%)", fontsize=11)
        ax.set_title(m, fontsize=13, fontweight="bold", color=method_colors.get(m, "black"))
        ax.set_xlim(0, 105)
        ax.axvline(75, color="gray", linestyle="--", alpha=0.5)
        ax.grid(axis="x", alpha=0.3)
        for i, v in enumerate(combo_pass.values):
            ax.text(v + 1, i, f"{v:.0f}%", va="center", fontsize=9)
    fig4.suptitle("PoseBusters Pass Rate by Protein-Ligand Combination", fontsize=14, fontweight="bold")
    fig4.tight_layout()
    fig4.savefig(plot_output_dir / "pb_pass_by_combo.png", dpi=200, bbox_inches="tight")
    plt.close(fig4)

    # FIGURE 5: distribution of #tests failed
    df["n_failed"] = (~df[test_cols]).sum(axis=1)
    fig5, ax5 = plt.subplots(figsize=(max(8, len(methods) * 2.5), 5))
    for i, m in enumerate(methods):
        vals = df[df["docking_method"] == m]["n_failed"]
        parts = ax5.violinplot([vals], positions=[i], showmedians=True, widths=0.7)
        for pc in parts["bodies"]:
            pc.set_facecolor(method_colors.get(m, "#888"))
            pc.set_alpha(0.4)
        for key in ("cbars", "cmins", "cmaxes", "cmedians"):
            parts[key].set_color(method_colors.get(m, "#888"))
        jitter = np.random.normal(0, 0.08, len(vals))
        ax5.scatter(np.full(len(vals), i) + jitter, vals, s=10, alpha=0.3,
                    color=method_colors.get(m, "#888"))
    ax5.set_xticks(range(len(methods)))
    ax5.set_xticklabels(methods, fontsize=12, fontweight="bold")
    ax5.set_ylabel("Number of Tests Failed", fontsize=12)
    ax5.set_title("Distribution of Failed Tests per Pose", fontsize=14, fontweight="bold")
    ax5.set_ylim(-0.5, df["n_failed"].max() + 1)
    ax5.grid(axis="y", alpha=0.3)
    for i, m in enumerate(methods):
        med = df[df["docking_method"] == m]["n_failed"].median()
        ax5.text(i, med + 0.3, f"median={med:.0f}", ha="center", fontsize=9, fontweight="bold")
    fig5.tight_layout()
    fig5.savefig(plot_output_dir / "pb_nfailed_violin.png", dpi=200, bbox_inches="tight")
    plt.close(fig5)

    print(f"\nAll figures saved to: {plot_output_dir}")
    for name, label in (
        ("pb_passrate_heatmap.png", "Test pass rate heatmap"),
        ("pb_pass_fail_bars.png",   "Pass/fail stacked bars"),
        ("pb_failure_rates.png",    "Per-test failure rates"),
        ("pb_pass_by_combo.png",    "Pass rate by protein-ligand"),
        ("pb_nfailed_violin.png",   "Distribution of #tests failed"),
    ):
        print(f"  {name}  — {label}")


# ============================================================================
# DIAGNOSTIC
# ============================================================================

def print_bottleneck_diagnostics(ctx: PipelineConfig) -> None:
    csv_path = ctx.output_dir / "posebusters_filtered_results.csv"
    if not csv_path.exists():
        return
    df_pb = pd.read_csv(csv_path, low_memory=False)
    for test in ("no_radicals", "non-aromatic_ring_non-flatness", "internal_steric_clash"):
        if test in df_pb.columns:
            print(f"\n{test} pass rate by method:")
            series = pd.to_numeric(df_pb[test], errors="coerce")
            print(series.groupby(df_pb["docking_method"]).mean().round(3) * 100)


# ============================================================================
# TOP-LEVEL ENTRY POINT
# ============================================================================

def run_pipeline(
    config_path: str | Path,
    overwrite: bool | None = None,
    quiet: bool | None = None,
    num_workers: int | None = None,
    pose_timeout: int | None = None,
    worker_maxtasks: int | None = None,
    monitor: bool | None = None,
    monitor_interval: int | None = None,
) -> pd.DataFrame:
    """Run the full PoseBusters validation pipeline driven by a YAML config.

    Any of *overwrite*, *quiet*, *num_workers*, *pose_timeout*, *worker_maxtasks*,
    *monitor*, *monitor_interval* override the corresponding config value when
    provided (used by the CLI).
    """
    ctx = load_config(config_path)
    if overwrite is not None:
        ctx.overwrite = overwrite
    if quiet is not None:
        ctx.quiet = quiet
    if num_workers is not None:
        ctx.num_workers = num_workers
    if pose_timeout is not None:
        ctx.pose_timeout = pose_timeout or None
    if worker_maxtasks is not None:
        ctx.worker_maxtasks = worker_maxtasks or None
    if monitor is not None:
        ctx.monitor = monitor
    if monitor_interval is not None:
        ctx.monitor_interval = monitor_interval

    # Apply quiet mode + thread caps in the main process now, before the Pool is
    # forked, so workers inherit them and the noisy per-file listings below are
    # already suppressed.
    _apply_quiet_mode(ctx.quiet)
    if ctx.limit_worker_threads:
        _limit_worker_threads()

    # Prepared-per-pose mode bypasses the historical global HETATM cleaner. Each
    # expanded pose is instead paired with its method/variant-specific receptor.
    protein_cache: dict[str, str | None] = {}
    if ctx.receptor_policy == "legacy_shared":
        ctx.receptors_dir = prepare_cleaned_receptor_dir(ctx)
        file_map, normalized_map = discover_proteins(ctx.receptors_dir)
    else:
        file_map, normalized_map = {}, {}

    filtered_poses_df, _ = collect_pose_rows(ctx)
    filtered_poses_df = filter_common_combos(filtered_poses_df, ctx.optional_methods)
    if ctx.expected_common_combinations is not None:
        actual_common = len(
            filtered_poses_df[["protein", "ligand"]].drop_duplicates()
        ) if not filtered_poses_df.empty else 0
        if actual_common != ctx.expected_common_combinations:
            raise RuntimeError(
                "Common-complex cohort changed: expected "
                f"{ctx.expected_common_combinations}, found {actual_common}. "
                "Review incomplete/missing docking artifacts before changing the benchmark lock."
            )
        print(f"  Cohort lock verified: {actual_common} common complexes")
    if ctx.receptor_policy == "legacy_shared":
        protein_cache = resolve_protein_files(filtered_poses_df, file_map, normalized_map)

    print("\n" + "=" * 80)
    print("PoseBusters Pose Validation — FILTERED SUBSET (PARALLEL)")
    print("=" * 80)
    results_csv_path = ctx.output_dir / "posebusters_filtered_results.csv"
    print("\n  *** OVERWRITE mode ON ***" if ctx.overwrite
          else f"\n  Resume mode: will skip poses already in {results_csv_path.name}")

    if filtered_poses_df.empty:
        print("\nERROR: No filtered poses available!")
        return pd.DataFrame()

    print(f"\nFiltered subset: {len(filtered_poses_df)} combinations, "
          f"~{filtered_poses_df['pose_count'].sum() if not filtered_poses_df.empty else 0:,} poses")

    print("\n" + "-" * 80 + "\n1. Collecting individual pose files...\n" + "-" * 80)
    all_poses = collect_all_pose_files(filtered_poses_df, ctx)
    all_poses = attach_validation_receptors(all_poses, ctx, protein_cache)
    print(f"\n   TOTAL: {len(all_poses)} individual poses to validate")
    for m, c in sorted(Counter(p["method"] for p in all_poses).items()):
        print(f"      {m}: {c}")

    print("\n" + "-" * 80)
    print(f"2. Running PoseBusters (config='{ctx.config_mode}', "
          f"workers={ctx.num_workers or 'all CPUs'}, "
          f"interval={ctx.save_interval}, overwrite={ctx.overwrite})")
    print("-" * 80)

    results_df = analyze_poses_with_posebusters(
        all_poses,
        protein_cache=protein_cache,
        config=ctx.config_mode,
        output_file=results_csv_path,
        save_interval=ctx.save_interval,
        overwrite=ctx.overwrite,
        num_workers=ctx.num_workers,
        quiet=ctx.quiet,
        worker_maxtasks=ctx.worker_maxtasks,
        pose_timeout=ctx.pose_timeout,
        limit_threads=ctx.limit_worker_threads,
        monitor=ctx.monitor,
        monitor_interval=ctx.monitor_interval,
        stall_timeout=ctx.stall_timeout,
        max_restarts=ctx.max_restarts,
        require_complete=ctx.require_complete_results,
    )

    if not results_df.empty:
        print(f"\n   Results saved to: {results_csv_path}")
        print(f"   Total rows: {len(results_df)}")
        if "posebusters_mode" in results_df.columns:
            print("\n   Mode breakdown:")
            print(results_df["posebusters_mode"].value_counts().to_string(header=False))
    else:
        print("\n   No results!")

    print("\n" + "=" * 80 + "\nDone!\n" + "=" * 80)

    print_bottleneck_diagnostics(ctx)

    if ctx.copy_proved_poses and not results_df.empty:
        copy_proved_poses(ctx)

    if ctx.make_plots and not results_df.empty:
        make_plots(ctx)

    return results_df


def _parse_cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", "-c", default="config.yaml",
        help="Path to the YAML config file (default: config.yaml)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", default=None,
        help="Overwrite existing PoseBusters results instead of resuming",
    )
    parser.add_argument(
        "--quiet", "-q", action="store_true", default=None,
        help="Suppress warning chatter (Python + RDKit warnings) and the verbose "
             "per-file listings; keep progress/checkpoint/summary output",
    )
    parser.add_argument(
        "--workers", "-j", type=int, default=None, metavar="N",
        help="Number of parallel worker processes (default: all CPUs)",
    )
    parser.add_argument(
        "--pose-timeout", type=int, default=None, metavar="SECONDS",
        help="Per-pose wall-clock cap for a single bust() call; 0 disables "
             "(overrides config; default 600s)",
    )
    parser.add_argument(
        "--max-tasks-per-child", type=int, default=None, metavar="N",
        help="Recycle each worker after N poses to bound memory growth; 0 keeps "
             "workers for the whole run (overrides config; default 200)",
    )
    parser.add_argument(
        "--monitor", action="store_true", default=None,
        help="Print a live resource snapshot (per-process core use, RAM, swap) "
             "and flag bottlenecks during the bust loop",
    )
    parser.add_argument(
        "--monitor-interval", type=int, default=None, metavar="SECONDS",
        help="Seconds between resource-monitor snapshots (default 30)",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_cli()
    run_pipeline(
        args.config,
        overwrite=args.overwrite,
        quiet=args.quiet,
        num_workers=args.workers,
        pose_timeout=args.pose_timeout,
        worker_maxtasks=args.max_tasks_per_child,
        monitor=args.monitor,
        monitor_interval=args.monitor_interval,
    )
