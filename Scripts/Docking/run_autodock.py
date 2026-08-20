#!/usr/bin/env python3
"""
AutoDock Vina batch docking CLI.

Usage:
    python run_autodock.py                          # uses default docking_config.yaml
    python run_autodock.py -c my_config.yaml        # custom config
    python run_autodock.py --dry-run                 # validate setup without docking
"""
from __future__ import annotations

import argparse
import csv
import contextlib
import hashlib
import itertools
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

try:  # Linux/macOS: protects optimization_log.csv across notebook processes.
    import fcntl
except ImportError:  # pragma: no cover - Windows fallback still has atomic replace
    fcntl = None

import numpy as np
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

# ── Ensure project root is importable ────────────────────────────────────────
_SCRIPT_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SCRIPT_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from Scripts.Utilities.prep_docking import run_workflow  # noqa: E402
from Scripts.Utilities.inject_hetatms import (  # noqa: E402
    inject_hetatms, DEFAULT_KEEP_COFACTORS, DEFAULT_KEEP_METALS,
)


# ═══════════════════════════════════════════════════════════════════════════════
# Configuration
# ═══════════════════════════════════════════════════════════════════════════════

FAILED_JOB_POLICIES = {"retry", "skip"}


def _failed_job_policy(cfg: dict) -> str:
    """Return the explicit policy for entries in the persistent failure log.

    ``retry`` is deliberately the default: a transient timeout or launch failure
    must not become a permanent cache entry that silently removes a complex from
    every later run. ``skip`` remains available for deliberate triage, but such
    entries are reported as failed (not as successfully skipped work).
    """
    policy = str(cfg.get("failed_job_policy", "retry")).strip().lower()
    if policy not in FAILED_JOB_POLICIES:
        raise ValueError(
            "failed_job_policy must be one of "
            f"{', '.join(sorted(FAILED_JOB_POLICIES))} (got {policy!r})"
        )
    return policy


def load_config(config_path: Path) -> dict:
    """Load and validate a YAML configuration file."""
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    # ── Validate required keys ──
    required = ["receptors_dir", "ligand_dirs", "output_dir", "vina_bin"]
    for key in required:
        if key not in cfg or cfg[key] is None:
            raise ValueError(f"Missing required config key: '{key}'")

    # ── Normalise prep_tool ──
    prep_tool = cfg.get("prep_tool", "both").lower().strip()
    if prep_tool not in ("mgltools", "meeko", "both"):
        raise ValueError(f"Invalid prep_tool '{prep_tool}'. Must be 'mgltools', 'meeko', or 'both'.")
    cfg["prep_tool"] = prep_tool

    scoring = cfg.get("scoring_function", "vina").lower()
    if scoring not in ("vina", "vinardo", "ad4"):
        raise ValueError(f"Invalid scoring_function '{scoring}'. Must be 'vina', 'vinardo', or 'ad4'.")
    cfg["scoring_function"] = scoring

    cfg["failed_job_policy"] = _failed_job_policy(cfg)

    # ── Validate post-pose optimisation mode (including explicit CNN refinement) ──
    # Raises early on a typo rather than silently docking without re-ranking.
    resolve_optimizers(cfg)
    _validate_optimizer_config(cfg)

    # retain_hetatm_residues: None/false = off, True = default KEEP list (cofactors+metals),
    # explicit list = use those 3-letter residue codes only.
    raw_keep = cfg.get("retain_hetatm_residues", None)
    if raw_keep in (None, False):
        cfg["retain_hetatm_residues"] = None
    elif raw_keep is True:
        cfg["retain_hetatm_residues"] = sorted(DEFAULT_KEEP_COFACTORS | DEFAULT_KEEP_METALS)
    else:
        cfg["retain_hetatm_residues"] = sorted(
            {str(r).strip().upper() for r in raw_keep if str(r).strip()}
        )

    return cfg


def print_config(cfg: dict) -> None:
    """Print a formatted summary of the configuration."""
    print("=" * 80)
    print("AutoDock Vina Batch Docking Configuration")
    print("=" * 80)
    print(f"Vina executable:     {cfg['vina_bin']}")
    print(f"Scoring function:    {cfg['scoring_function']}")
    print(f"Prep tool(s):        {cfg['prep_tool']}")

    batch_mode = cfg.get("batch_mode")
    batch_str = {None: "auto", True: "forced ON", False: "disabled"}.get(batch_mode, str(batch_mode))
    print(f"Batch mode:          {batch_str}")
    print(f"Batch size:          {cfg.get('batch_size', 0) or 'all at once'}")
    bt = cfg.get("batch_timeout", 0)
    print(f"Batch timeout:       {bt}s" if bt else "Batch timeout:       derived from per-complex timeout")

    if cfg["scoring_function"] == "ad4":
        print(f"AutoGrid4:           {cfg.get('autogrid_bin', 'N/A')}")
    print(f"Exhaustiveness:      {cfg.get('exhaustiveness', 32)}")
    print(f"Num modes:           {cfg.get('num_modes', 10)}")
    print(f"Energy range:        {cfg.get('energy_range', 3)} kcal/mol")
    tpc = cfg.get("timeout_per_complex", 0)
    print(f"Timeout per complex: {tpc}s ({tpc / 60:.0f} min)" if tpc else "Timeout per complex: unlimited")
    print(f"Workers:             {cfg.get('max_workers', 1)}  (CPUs/worker: {cfg.get('cpus_per_worker', 4)})")
    keep = cfg.get("retain_hetatm_residues")
    if keep:
        n = len(keep)
        sample = ", ".join(keep[:8]) + (f", … (+{n-8} more)" if n > 8 else "")
        print(f"Retain HETATM:       YES  ({n} residue types: {sample})")
    else:
        print(f"Retain HETATM:       no")
    print(f"Overwrite docking:   {cfg.get('overwrite_existing', False)}")
    print(f"Overwrite poses:     {cfg.get('overwrite_poses', False)}")
    print(f"Prior failures:      {_failed_job_policy(cfg)}")
    opt_tools = resolve_optimizers(cfg)
    if opt_tools:
        print(f"Pose optimization:   {', '.join(opt_tools)}  "
              f"(search={cfg.get('optimize_search', 'minimize')}, "
              f"autobox_add={cfg.get('optimize_autobox_add', 4.0)} Å)")
        print(f"  Re-rank metric:    {cfg.get('optimize_rank_by', 'minimized_affinity')}")
        print(f"  Empirical score:   {cfg.get('optimize_scoring', 'default')}")
        top_n = int(cfg.get("optimize_top_n", 0) or 0)
        print(f"  Re-rank poses:     top {top_n}" if top_n else "  Re-rank poses:     all poses")
        if any(_is_gnina_optimizer(tool) for tool in opt_tools):
            print(f"  gnina GPU:         {'on' if cfg.get('gnina_use_gpu', False) else 'off (--no_gpu)'}")
            if cfg.get("gnina_use_gpu", False):
                print(f"  GPU lock:          {cfg.get('gpu_lock_file')}")
            print(f"  CNN scoring/model: {cfg.get('gnina_cnn_scoring', 'rescore')} / "
                  f"{cfg.get('gnina_cnn_model', 'crossdock_default2018_ensemble')}")
    else:
        print(f"Pose optimization:   off")
    print(f"Output dir:          {cfg['output_dir']}")
    print(f"Log dir:             {cfg.get('log_dir', 'logs')}")
    print(f"Receptor dir:        {cfg['receptors_dir']}")
    print(f"Ligand dirs ({len(cfg['ligand_dirs'])}):")
    for d in cfg["ligand_dirs"]:
        print(f"  • {d}")
    print()


# ═══════════════════════════════════════════════════════════════════════════════
# Data classes & helpers
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class DockingResult:
    """Result of docking a protein-ligand combination."""
    protein_name: str
    ligand_name: str
    protein_path: Path
    ligand_path: Path
    output_dir: Path
    status: str  # "success", "failed", "skipped"
    scoring_function: str = "vina"
    num_poses: int = 0
    pose_files: List[Path] = field(default_factory=list)
    best_affinity: Optional[float] = None
    error_message: str = ""
    elapsed_time: float = 0.0

    def to_dict(self) -> dict:
        return {
            "protein_name": self.protein_name,
            "ligand_name": self.ligand_name,
            "protein_path": str(self.protein_path),
            "ligand_path": str(self.ligand_path),
            "output_dir": str(self.output_dir),
            "status": self.status,
            "scoring_function": self.scoring_function,
            "num_poses": self.num_poses,
            "pose_files": [str(p) for p in self.pose_files],
            "best_affinity": self.best_affinity,
            "error_message": self.error_message,
            "elapsed_time": self.elapsed_time,
        }


def get_file_stem(path: Path) -> str:
    return path.stem.replace("_ligand", "").replace("_protein", "")


def _inject_hetatms_into_protein_outputs(
    protein_outputs: dict,
    keep_residues: list[str] | None,
    label: str = "",
) -> None:
    """For every (cleaned_pdb -> pdbqt) entry produced by `run_workflow`, append
    KEEP-list HETATM atoms (cofactors + metals) from the *original* source PDB
    into the receptor PDBQT in place.

    The dict key in `protein_outputs` points to the PDBFixer-cleaned PDB, which
    has already had heterogens stripped. We therefore look up the matching
    original PDB in ``protein_outputs["protein_files"]`` by stem (the cleaner
    appends a converter postfix like ``_meeko`` or ``_mgl_tools``) and extract
    HETATMs from there. No-op when keep_residues is falsy.

    Idempotent only when called once per pdbqt; calling twice would double-inject.
    """
    if not keep_residues:
        return
    keep_set = set(keep_residues)

    # Build {original_stem -> original_pdb_path} from the workflow's
    # original input list. `protein_files` is a list/set of pre-PDBFixer paths.
    originals = list(protein_outputs.get("protein_files") or [])
    orig_by_stem: dict[str, Path] = {Path(p).stem: Path(p) for p in originals}

    # Common converter postfixes that the cleaning step appends to the stem.
    _POSTFIXES = ("_meeko", "_mgl_tools", "_pymol", "_openbabel", "_mko", "_mgl")

    def _resolve_original(cleaned_pdb: Path) -> Path | None:
        """Map a cleaned PDB path back to its original source PDB."""
        stem = cleaned_pdb.stem
        for sfx in _POSTFIXES:
            if stem.endswith(sfx):
                cand = stem[: -len(sfx)]
                if cand in orig_by_stem:
                    return orig_by_stem[cand]
        # Fallback: exact stem match (no postfix applied)
        if stem in orig_by_stem:
            return orig_by_stem[stem]
        # Fallback: any original whose stem is a prefix of cleaned stem
        for s, p in orig_by_stem.items():
            if stem.startswith(s):
                return p
        return None

    for key in ("meeko_outputs", "converter_protein_outputs"):
        produced = protein_outputs.get(key) or {}
        injected = 0
        atoms_added = 0
        for cleaned_pdb_str, pdbqt in produced.items():
            if not pdbqt:
                continue
            cleaned_path = Path(cleaned_pdb_str)
            pdbqt_path = Path(pdbqt)
            if not pdbqt_path.exists():
                continue
            src_path = _resolve_original(cleaned_path)
            if src_path is None or not src_path.exists():
                # As a last resort, try the dict key (works for callers that
                # passed the original directly, with no PDBFixer cleaning).
                src_path = cleaned_path if cleaned_path.exists() else None
            if src_path is None:
                continue
            try:
                result = inject_hetatms(src_path, pdbqt_path, keep_residues=keep_set)
            except Exception as exc:
                print(f"  ⚠ HETATM injection failed for {pdbqt_path.name}: {exc}")
                continue
            if result["atoms_injected"] > 0:
                injected += 1
                atoms_added += result["atoms_injected"]
                print(f"  ✓ {pdbqt_path.name}: +{result['atoms_injected']} atoms "
                      f"({result['by_resname']})  ← {src_path.name}")
        if produced:
            print(f"  [{label or key}] HETATM injection: {injected}/{len(produced)} "
                  f"receptors augmented, +{atoms_added} atoms total")


def collect_files(root: Path, extensions: List[str]) -> List[Path]:
    files = []
    for ext in extensions:
        files.extend(sorted(p for p in root.glob(f"*{ext}") if p.is_file()))
    return sorted(set(files))


def get_pdbqt_dir(base_dir: Path) -> Path:
    d = base_dir / "pdbqt"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ═══════════════════════════════════════════════════════════════════════════════
# Vina log / PDBQT parsing
# ═══════════════════════════════════════════════════════════════════════════════

def parse_vina_affinities(log_text: str) -> List[dict]:
    modes: List[dict] = []
    mode_pattern = r'^\s*(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)\s*$'
    for m in re.finditer(mode_pattern, log_text, re.MULTILINE):
        modes.append({
            'mode': int(m.group(1)),
            'affinity': float(m.group(2)),
            'rmsd_lb': float(m.group(3)),
            'rmsd_ub': float(m.group(4)),
        })
    return modes


def parse_pdbqt_affinities(pdbqt_text: str) -> List[dict]:
    modes: List[dict] = []
    remark_pattern = r'^REMARK VINA RESULT:\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)'
    for idx, m in enumerate(re.finditer(remark_pattern, pdbqt_text, re.MULTILINE), 1):
        modes.append({
            'mode': idx,
            'affinity': float(m.group(1)),
            'rmsd_lb': float(m.group(2)),
            'rmsd_ub': float(m.group(3)),
        })
    return modes


_VINA_RESULT_LINE = re.compile(
    r"^REMARK VINA RESULT:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s+"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)\s*$"
)


def validate_vina_output(pdbqt_path: Path | str) -> List[dict]:
    """Validate a committed raw Vina output and return its scored poses.

    A resumable result must be a completely written multi-model PDBQT: every
    model has a unique consecutive rank, a structurally closed ligand torsion
    tree, exactly one finite ``REMARK VINA RESULT`` record, and a matching
    ``ENDMDL``. Any non-whitespace outside the model envelopes is rejected. This
    intentionally does *not* treat an arbitrary non-empty or partially written
    file as completed work.

    Raises ``ValueError``/``OSError`` on invalid input. Returned dictionaries use
    the same ``mode``/``affinity``/RMSD fields as the other public parsers.
    """
    path = Path(pdbqt_path)
    if not path.is_file():
        raise ValueError(f"Vina output is not a regular file: {path}")
    text = path.read_text(errors="strict")
    if not text.strip():
        raise ValueError(f"empty Vina output: {path}")
    if "\x00" in text:
        raise ValueError(f"Vina output contains NUL bytes: {path}")

    models: List[dict] = []
    current_rank: Optional[int] = None
    current_lines: List[str] = []
    seen_ranks: Set[int] = set()

    for line_number, line in enumerate(text.splitlines(), 1):
        if line.startswith("MODEL"):
            if current_rank is not None:
                raise ValueError(
                    f"nested MODEL at line {line_number}; MODEL {current_rank} is truncated"
                )
            match = re.fullmatch(r"MODEL\s+(\d+)\s*", line)
            if not match:
                raise ValueError(f"malformed MODEL record at line {line_number}")
            current_rank = int(match.group(1))
            if current_rank < 1 or current_rank in seen_ranks:
                raise ValueError(f"invalid or duplicate MODEL rank {current_rank}")
            current_lines = []
            continue

        if line.startswith("ENDMDL"):
            if not re.fullmatch(r"ENDMDL\s*", line):
                raise ValueError(f"malformed ENDMDL record at line {line_number}")
            if current_rank is None:
                raise ValueError(f"ENDMDL without open MODEL at line {line_number}")
            body = "\n".join(current_lines) + "\n"
            _validate_pdbqt_pose_body(body, f"MODEL {current_rank}")
            score_lines = [
                _VINA_RESULT_LINE.fullmatch(record)
                for record in current_lines
                if record.startswith("REMARK VINA RESULT:")
            ]
            if len(score_lines) != 1 or score_lines[0] is None:
                raise ValueError(
                    f"MODEL {current_rank} must contain exactly one valid Vina result"
                )
            affinity, rmsd_lb, rmsd_ub = (
                float(value) for value in score_lines[0].groups()
            )
            if not all(math.isfinite(value) for value in (affinity, rmsd_lb, rmsd_ub)):
                raise ValueError(f"MODEL {current_rank} has non-finite score values")
            models.append({
                "mode": current_rank,
                "affinity": affinity,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub,
            })
            seen_ranks.add(current_rank)
            current_rank = None
            current_lines = []
            continue

        if current_rank is None:
            if line.strip():
                raise ValueError(
                    f"content outside MODEL…ENDMDL at line {line_number}: {line[:60]!r}"
                )
        else:
            current_lines.append(line)

    if current_rank is not None:
        raise ValueError(f"MODEL {current_rank} is truncated (missing ENDMDL)")
    if not models:
        raise ValueError(f"no complete scored MODEL…ENDMDL poses in {path}")
    expected_ranks = list(range(1, len(models) + 1))
    ranks = [model["mode"] for model in models]
    if ranks != expected_ranks:
        raise ValueError(
            f"MODEL ranks must be consecutive and ordered from 1 (got {ranks})"
        )
    return models


def _quarantine_vina_output(path: Path, reason: str) -> Optional[Path]:
    """Atomically move an invalid/partial pose away from resumable results."""
    path = Path(path)
    if not path.exists():
        return None
    quarantine_dir = path.parent / "_invalid_vina_outputs"
    quarantine_dir.mkdir(parents=True, exist_ok=True)
    token = f"{time.time_ns()}-{os.getpid()}-{threading.get_ident()}"
    quarantined = quarantine_dir / f"{path.name}.{token}.invalid"
    os.replace(path, quarantined)
    try:
        _atomic_write_json(quarantined.with_name(f"{quarantined.name}.json"), {
            "schema_version": 1,
            "original_path": str(path),
            "quarantined_path": str(quarantined),
            "reason": str(reason),
            "timestamp": datetime.now().isoformat(),
        })
    except Exception as exc:  # the pose is already safe; retain it even if audit metadata fails
        print(f"  ⚠ Could not write quarantine metadata for {quarantined}: {exc}")
    return quarantined


def _archive_superseded_vina_output(path: Path, reason: str) -> Optional[Path]:
    """Recoverably remove a valid prior commit before an overwrite attempt.

    The old file cannot remain at its canonical path: if its replacement fails,
    a later overwrite=false resume would otherwise accept the old generation and
    silently collapse the failed attempt back to "done".
    """
    path = Path(path)
    if not path.exists():
        return None
    archive_dir = path.parent / "_superseded_vina_outputs"
    archive_dir.mkdir(parents=True, exist_ok=True)
    token = f"{time.time_ns()}-{os.getpid()}-{threading.get_ident()}"
    archived = archive_dir / f"{path.name}.{token}.superseded"
    os.replace(path, archived)
    try:
        _atomic_write_json(archived.with_name(f"{archived.name}.json"), {
            "schema_version": 1,
            "original_path": str(path),
            "archived_path": str(archived),
            "reason": str(reason),
            "timestamp": datetime.now().isoformat(),
        })
    except Exception as exc:
        print(f"  ⚠ Could not write archive metadata for {archived}: {exc}")
    return archived


# ═══════════════════════════════════════════════════════════════════════════════
# Error / docking logs
# ═══════════════════════════════════════════════════════════════════════════════

def _error_log_filename(dir_name: str) -> str:
    return f"failed_docking_{dir_name}.json"


def _docking_log_filename(dir_name: str) -> str:
    return f"docking_log_{dir_name}.csv"


def load_error_log(output_dir: Path) -> Dict[str, dict]:
    path = output_dir / _error_log_filename(output_dir.name)
    if path.exists():
        try:
            with open(path, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, Exception) as e:
            print(f"  ⚠ Could not read error log {path}: {e}")
    return {}


def save_error_log(output_dir: Path, error_log: Dict[str, dict]):
    path = output_dir / _error_log_filename(output_dir.name)
    _atomic_write_json(path, error_log)


def record_failure(output_dir: Path, combo_name: str, result: DockingResult):
    error_log = load_error_log(output_dir)
    error_log[combo_name] = {
        "error": result.error_message,
        "timestamp": datetime.now().isoformat(),
        "protein": result.protein_name,
        "ligand": result.ligand_name,
        "elapsed_time": round(result.elapsed_time, 2),
    }
    save_error_log(output_dir, error_log)


# ═══════════════════════════════════════════════════════════════════════════════
# Molecular-property helpers
# ═══════════════════════════════════════════════════════════════════════════════

DOCKING_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name",
    "status", "scoring_function", "error_reason", "num_poses", "elapsed_time_s",
    "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula",
    "prot_num_residues", "prot_num_atoms", "prot_num_chains",
    "exhaustiveness", "num_modes", "energy_range", "cpu_model",
    "best_affinity_kcal",
]


def get_cpu_model() -> str:
    try:
        with open("/proc/cpuinfo", "r") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    try:
        import platform
        return platform.processor() or "N/A"
    except Exception:
        return "N/A"


def get_ligand_properties_from_sdf(sdf_path: Path) -> Dict:
    props = {
        "lig_molecular_weight": None, "lig_heavy_atoms": None,
        "lig_total_atoms": None, "lig_rotatable_bonds": None,
        "lig_num_rings": None, "lig_aromatic_rings": None,
        "lig_hbd": None, "lig_hba": None, "lig_tpsa": None,
        "lig_logp": None, "lig_formula": None,
    }
    try:
        supplier = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        mol = next(iter(supplier), None)
        if mol is None:
            return props
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        props["lig_molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        props["lig_heavy_atoms"] = mol.GetNumHeavyAtoms()
        props["lig_total_atoms"] = mol.GetNumAtoms()
        props["lig_rotatable_bonds"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        ri = mol.GetRingInfo()
        props["lig_num_rings"] = ri.NumRings()
        props["lig_aromatic_rings"] = rdMolDescriptors.CalcNumAromaticRings(mol)
        props["lig_hbd"] = Lipinski.NumHDonors(mol)
        props["lig_hba"] = Lipinski.NumHAcceptors(mol)
        props["lig_tpsa"] = round(Descriptors.TPSA(mol), 2)
        props["lig_logp"] = round(Descriptors.MolLogP(mol), 2)
        props["lig_formula"] = rdMolDescriptors.CalcMolFormula(mol)
    except Exception as e:
        print(f"  ⚠ Could not compute ligand properties for {sdf_path.name}: {e}")
    return props


def get_ligand_properties_from_pdbqt(pdbqt_path: Path) -> Dict:
    props = {
        "lig_molecular_weight": None, "lig_heavy_atoms": None,
        "lig_total_atoms": None, "lig_rotatable_bonds": None,
        "lig_num_rings": None, "lig_aromatic_rings": None,
        "lig_hbd": None, "lig_hba": None, "lig_tpsa": None,
        "lig_logp": None, "lig_formula": None,
    }
    try:
        atom_count = 0
        heavy_count = 0
        rotatable = 0
        with open(pdbqt_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_count += 1
                    element = line[77:79].strip() if len(line) > 77 else ""
                    if element and element not in ("H", "HD"):
                        heavy_count += 1
                elif line.startswith("BRANCH"):
                    rotatable += 1
        props["lig_total_atoms"] = atom_count
        props["lig_heavy_atoms"] = heavy_count if heavy_count > 0 else None
        # Each rotatable bond is one BRANCH/ENDBRANCH pair, but only "BRANCH"
        # lines are counted above ("ENDBRANCH".startswith("BRANCH") is False),
        # so the raw count already equals the number of active torsions.
        props["lig_rotatable_bonds"] = rotatable if rotatable > 0 else None
    except Exception as e:
        print(f"  ⚠ Could not parse PDBQT ligand {pdbqt_path.name}: {e}")
    return props


def get_ligand_properties(lig_path: Path) -> Dict:
    if lig_path.suffix.lower() in (".sdf", ".mol2"):
        return get_ligand_properties_from_sdf(lig_path)
    return get_ligand_properties_from_pdbqt(lig_path)


def get_protein_properties(pdb_path: Path) -> Dict:
    props = {"prot_num_residues": None, "prot_num_atoms": None, "prot_num_chains": None}
    try:
        residues = set()
        chains = set()
        atom_count = 0
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1
                    chain_id = line[21]
                    res_seq = line[22:27].strip()
                    residues.add((chain_id, res_seq))
                    chains.add(chain_id)
        props["prot_num_residues"] = len(residues)
        props["prot_num_atoms"] = atom_count
        props["prot_num_chains"] = len(chains)
    except Exception as e:
        print(f"  ⚠ Could not parse protein {pdb_path.name}: {e}")
    return props


def precompute_properties(
    proteins: List[Path],
    ligands: List[Path],
) -> Tuple[Dict[Path, Dict], Dict[Path, Dict]]:
    lig_props_cache: Dict[Path, Dict] = {}
    prot_props_cache: Dict[Path, Dict] = {}
    n_workers = min(os.cpu_count() or 4, 8)

    print("Pre-computing ligand properties...")
    with ProcessPoolExecutor(max_workers=min(n_workers, len(ligands) or 1)) as pool:
        futures = {pool.submit(get_ligand_properties, lig): lig for lig in ligands}
        for fut in as_completed(futures):
            lig = futures[fut]
            key = lig.resolve()  # match the resolved paths used at lookup time
            try:
                lig_props_cache[key] = fut.result()
            except Exception as e:
                print(f"  ⚠ Failed computing props for {lig.name}: {e}")
                lig_props_cache[key] = get_ligand_properties_from_pdbqt(lig)
    print(f"  {len(lig_props_cache)} ligands analysed")

    print("Pre-computing protein properties...")
    with ProcessPoolExecutor(max_workers=min(n_workers, len(proteins) or 1)) as pool:
        futures = {pool.submit(get_protein_properties, prot): prot for prot in proteins}
        for fut in as_completed(futures):
            prot = futures[fut]
            key = prot.resolve()  # match the resolved paths used at lookup time
            try:
                prot_props_cache[key] = fut.result()
            except Exception as e:
                print(f"  ⚠ Failed computing props for {prot.name}: {e}")
                prot_props_cache[key] = {"prot_num_residues": None, "prot_num_atoms": None, "prot_num_chains": None}
    print(f"  {len(prot_props_cache)} proteins analysed")

    return lig_props_cache, prot_props_cache


def _build_log_row(
    result: DockingResult,
    lig_props: Dict,
    prot_props: Dict,
    cfg: dict,
    cpu_model: str,
) -> Dict:
    row = {
        "timestamp": datetime.now().isoformat(),
        "combo_name": f"{result.ligand_name}__{result.protein_name}",
        "protein_name": result.protein_name,
        "ligand_name": result.ligand_name,
        "status": result.status,
        "scoring_function": result.scoring_function,
        "error_reason": result.error_message if result.status == "failed" else "",
        "num_poses": result.num_poses,
        "elapsed_time_s": round(result.elapsed_time, 2),
        "best_affinity_kcal": result.best_affinity,
        "exhaustiveness": cfg.get("exhaustiveness", 32),
        "num_modes": cfg.get("num_modes", 10),
        "energy_range": cfg.get("energy_range", 3),
        "cpu_model": cpu_model,
    }
    row.update(lig_props)
    row.update(prot_props)
    return row


def init_docking_log(output_dir: Path) -> Tuple[Path, set]:
    log_path = output_dir / _docking_log_filename(output_dir.name)
    existing_combos: set = set()
    if log_path.exists():
        try:
            df = pd.read_csv(log_path)
            existing_combos = set(df["combo_name"])
        except Exception:
            pass
    else:
        pd.DataFrame(columns=DOCKING_LOG_COLUMNS).to_csv(log_path, index=False)
    return log_path, existing_combos


def append_log_rows(
    results: List[DockingResult],
    log_path: Path,
    existing_combos: set,
    lig_props_cache: Dict[Path, Dict],
    prot_props_cache: Dict[Path, Dict],
    cfg: dict,
    cpu_model: str,
    overwrite: bool = False,
) -> None:
    """Append log rows for `results` to the docking-log CSV.

    Always appends (O(1) amortised). Stale rows for re-docked combos (overwrite
    mode) are left in place and superseded by the freshly-appended rows, then
    collapsed in a single pass by compact_docking_log() at the end of the run.
    This replaces the former per-result full-file read+filter+rewrite, which was
    O(N^2) across a run and ran while holding the log lock.
    """
    rows = []
    for r in results:
        combo = f"{r.ligand_name}__{r.protein_name}"
        if combo in existing_combos and not overwrite:
            continue
        lig_props = lig_props_cache.get(r.ligand_path)
        if lig_props is None:
            lig_props = get_ligand_properties(r.ligand_path)
        prot_props = prot_props_cache.get(r.protein_path)
        if prot_props is None:
            prot_props = get_protein_properties(r.protein_path)
        rows.append(_build_log_row(r, lig_props, prot_props, cfg, cpu_model))
        existing_combos.add(combo)

    if not rows:
        return

    # The CSV header is written once by init_docking_log(); always append.
    df = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
    df.to_csv(log_path, mode="a", header=False, index=False)


def compact_docking_log(log_path: Path) -> None:
    """Collapse the docking-log CSV to one row per combo, keeping the most
    recent. Run once at the end of a docking run so stale/duplicate rows left
    by append-only writes and overwrite re-docks are removed in a single O(N)
    pass instead of the former per-result O(N^2) rewrite.

    "Keep last" is correct for both modes: on overwrite the new row is appended
    after the stale one; without overwrite an already-logged combo is never
    re-appended, so its single row is preserved unchanged.
    """
    if not log_path.exists():
        return
    try:
        df = pd.read_csv(log_path)
    except Exception as e:
        print(f"  ⚠ Could not compact docking log {log_path.name}: {e}")
        return
    if df.empty or "combo_name" not in df.columns:
        return
    before = len(df)
    df = df.drop_duplicates(subset="combo_name", keep="last")
    if len(df) != before:
        df.to_csv(log_path, index=False)
        print(f"  📝 Compacted {log_path.name}: {before} → {len(df)} rows")


# ═══════════════════════════════════════════════════════════════════════════════
# Summary & reporting
# ═══════════════════════════════════════════════════════════════════════════════

def generate_summary(results: List[DockingResult], cfg: dict) -> Dict:
    summary = {
        "timestamp": datetime.now().isoformat(),
        "configuration": {
            "scoring_function": cfg.get("scoring_function", "vina"),
            "exhaustiveness": cfg.get("exhaustiveness", 32),
            "num_modes": cfg.get("num_modes", 10),
            "energy_range": cfg.get("energy_range", 3),
            "vina_bin": str(cfg["vina_bin"]),
        },
        "overall": {
            "total_combinations": len(results),
            "successful": sum(1 for r in results if r.status == "success"),
            "failed": sum(1 for r in results if r.status == "failed"),
            "skipped": sum(1 for r in results if r.status == "skipped"),
            "total_poses": sum(r.num_poses for r in results),
            "total_time_seconds": sum(r.elapsed_time for r in results),
        },
        "by_protein": defaultdict(lambda: {"combinations": 0, "poses": 0, "success": 0}),
        "by_ligand": defaultdict(lambda: {"combinations": 0, "poses": 0, "success": 0}),
        "combinations": [],
    }
    for r in results:
        summary["by_protein"][r.protein_name]["combinations"] += 1
        summary["by_protein"][r.protein_name]["poses"] += r.num_poses
        if r.status == "success":
            summary["by_protein"][r.protein_name]["success"] += 1
        summary["by_ligand"][r.ligand_name]["combinations"] += 1
        summary["by_ligand"][r.ligand_name]["poses"] += r.num_poses
        if r.status == "success":
            summary["by_ligand"][r.ligand_name]["success"] += 1
        summary["combinations"].append({
            "protein": r.protein_name,
            "ligand": r.ligand_name,
            "status": r.status,
            "num_poses": r.num_poses,
            "time_seconds": round(r.elapsed_time, 2),
            "best_affinity": r.best_affinity,
            "error": r.error_message if r.error_message else None,
        })
    summary["by_protein"] = dict(summary["by_protein"])
    summary["by_ligand"] = dict(summary["by_ligand"])
    return summary


def print_summary(summary: Dict):
    print("\n" + "=" * 80)
    print("AUTODOCK VINA DOCKING SUMMARY")
    print("=" * 80)
    overall = summary["overall"]
    print(f"\nOverall Statistics:")
    print(f"  Total combinations: {overall['total_combinations']}")
    print(f"  Successful: {overall['successful']}")
    print(f"  Failed: {overall['failed']}")
    print(f"  Skipped: {overall['skipped']}")
    print(f"  Total poses generated: {overall['total_poses']}")
    print(f"  Total time: {overall['total_time_seconds']:.1f}s ({overall['total_time_seconds']/60:.1f} min)")
    print("\n" + "-" * 80)
    print("Results by Protein:")
    print("-" * 80)
    for protein, stats in sorted(summary["by_protein"].items()):
        print(f"  {protein:40s} | Combos: {stats['combinations']:3d} | Poses: {stats['poses']:4d}")
    print("\n" + "-" * 80)
    print("Results by Ligand:")
    print("-" * 80)
    for ligand, stats in sorted(summary["by_ligand"].items()):
        print(f"  {ligand:40s} | Combos: {stats['combinations']:3d} | Poses: {stats['poses']:4d}")


# ═══════════════════════════════════════════════════════════════════════════════
# Receptor / ligand collection helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _normalize_key(path: Path) -> str:
    name = path.name
    if name.endswith(".box.txt"):
        return name[: -len(".box.txt")]
    if name.endswith(".pdbqt"):
        return name[: -len(".pdbqt")]
    return path.stem


def _collect_receptors(
    base_dir: Path,
    prepared_manifest: dict,
    workflow_data: dict | None,
    receptor_dirs: "tuple[Path, ...] | list[Path]" = (),
) -> list[dict]:
    groups: defaultdict[str, dict] = defaultdict(dict)

    for entry in prepared_manifest.get("proteins", []):
        entry_path = Path(entry.get("path", "")).resolve()
        if not entry_path.exists():
            continue
        key = _normalize_key(entry_path)
        if entry_path.suffix.lower() == ".pdbqt":
            groups[key]["pdbqt"] = entry_path
        elif entry_path.name.lower().endswith(".box.txt"):
            groups[key]["box"] = entry_path

    if workflow_data:
        pdbqt_map = {
            **(workflow_data.get("meeko_outputs") or {}),
            **(workflow_data.get("converter_protein_outputs") or {}),
        }
        box_map = {
            **(workflow_data.get("meeko_box_files") or {}),
            **(workflow_data.get("converter_box_files") or {}),
        }
        for original, pdbqt in pdbqt_map.items():
            box = box_map.get(original)
            if not pdbqt or not box:
                continue
            pdbqt_path = Path(pdbqt).resolve()
            box_path = Path(box).resolve()
            if not (pdbqt_path.exists() and box_path.exists()):
                continue
            key = _normalize_key(pdbqt_path)
            groups[key].setdefault("pdbqt", pdbqt_path)
            groups[key].setdefault("box", box_path)

    receptors = []
    for key, data in groups.items():
        pdbqt = data.get("pdbqt")
        box = data.get("box")
        if pdbqt and box:
            receptors.append({"name": key, "pdbqt": pdbqt, "box": box})

    if not receptors:
        for pdbqt_path in sorted(base_dir.glob("*.pdbqt")):
            box_path = pdbqt_path.with_suffix(".box.txt")
            if box_path.exists():
                receptors.append({
                    "name": _normalize_key(pdbqt_path),
                    "pdbqt": pdbqt_path.resolve(),
                    "box": box_path.resolve(),
                })

    if not receptors and receptor_dirs:
        # Receptor preparation produced nothing this run — e.g. mk_prepare_receptor
        # failing against a newer Meeko CLI — but the receptors that produced the
        # existing poses are still on disk. Recover exactly the stems those poses
        # were docked against, so a re-optimisation pass (``optimization: gnina``)
        # can run without re-preparing or re-docking anything. Deriving the stems
        # from the existing output filenames is what keeps receptor-to-pose parity:
        # it can never introduce a receptor variant the poses were not docked with.
        stems = {p.name.split("__", 1)[0]
                 for p in (base_dir / "docking").glob("*__*_vina_vina_out.pdbqt")
                 if "__" in p.name}
        for rdir in receptor_dirs:
            rdir = Path(rdir)
            if not rdir.is_dir():
                continue
            for stem in sorted(stems):
                pdbqt_path = rdir / f"{stem}.pdbqt"
                box_path = rdir / f"{stem}.box.txt"
                if pdbqt_path.exists() and box_path.exists():
                    receptors.append({
                        "name": _normalize_key(pdbqt_path),
                        "pdbqt": pdbqt_path.resolve(),
                        "box": box_path.resolve(),
                    })
        if receptors:
            print(f"  (receptor prep produced none; reusing {len(receptors)} already-prepared "
                  f"receptor(s) matching the existing poses)")

    receptors.sort(key=lambda r: r["name"])
    return receptors


def _collect_ligands(
    base_dir: Path,
    prepared_manifest: dict,
    workflow_data: dict | None,
) -> list[dict]:
    ligands: dict[str, Path] = {}

    for entry in prepared_manifest.get("ligands", []):
        entry_path = Path(entry.get("path", "")).resolve()
        if entry_path.exists() and entry_path.suffix.lower() == ".pdbqt":
            ligands.setdefault(_normalize_key(entry_path), entry_path)

    if workflow_data:
        lig_map = {
            **(workflow_data.get("meeko_ligand_outputs") or {}),
            **(workflow_data.get("converter_ligand_outputs") or {}),
        }
        for original, pdbqt in lig_map.items():
            if not pdbqt:
                continue
            pdbqt_path = Path(pdbqt).resolve()
            if pdbqt_path.exists():
                ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path)

    ligand_dir = base_dir / "ligands"
    if not ligands and ligand_dir.exists():
        for pdbqt_path in sorted(ligand_dir.glob("*.pdbqt")):
            ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path.resolve())

    if not ligands:
        for pdbqt_path in sorted(base_dir.glob("*.pdbqt")):
            if pdbqt_path.with_suffix(".box.txt").exists():
                continue
            ligands.setdefault(_normalize_key(pdbqt_path), pdbqt_path.resolve())

    return [{"name": n, "pdbqt": p} for n, p in sorted(ligands.items())]


# ═══════════════════════════════════════════════════════════════════════════════
# AutoGrid4 (ad4 scoring)
# ═══════════════════════════════════════════════════════════════════════════════

def _parse_box_txt(box_path: Path) -> dict:
    params = {}
    with open(box_path, "r") as f:
        for line in f:
            line = line.strip()
            if "=" in line:
                key, _, val = line.partition("=")
                key = key.strip()
                val = val.strip()
                try:
                    params[key] = float(val)
                except ValueError:
                    params[key] = val
    return params


def _extract_receptor_atom_types(pdbqt_path: Path) -> list[str]:
    atom_types: set[str] = set()
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 78:
                atype = line[77:79].strip()
                if atype:
                    atom_types.add(atype)
    return sorted(atom_types)


def _generate_ad4_maps(
    rec_pdbqt: Path,
    box_path: Path,
    output_dir: Path,
    autogrid_bin: Path,
    spacing: float = 0.375,
) -> str:
    output_dir.mkdir(parents=True, exist_ok=True)
    maps_prefix = output_dir / rec_pdbqt.stem
    fld_path = Path(f"{maps_prefix}.maps.fld")

    if fld_path.exists() and fld_path.stat().st_mtime > rec_pdbqt.stat().st_mtime:
        print(f"  ✓ ad4 maps cached for {rec_pdbqt.name}")
        return str(maps_prefix)

    box = _parse_box_txt(box_path)
    center_x = box.get("center_x", 0.0)
    center_y = box.get("center_y", 0.0)
    center_z = box.get("center_z", 0.0)
    size_x = box.get("size_x", 20.0)
    size_y = box.get("size_y", 20.0)
    size_z = box.get("size_z", 20.0)
    npts_x = int(size_x / spacing)
    npts_y = int(size_y / spacing)
    npts_z = int(size_z / spacing)

    ligand_types = ["A", "C", "HD", "H", "N", "NA", "OA", "F", "Cl", "Br", "I", "S", "SA"]
    rec_types = _extract_receptor_atom_types(rec_pdbqt)
    all_types = sorted(set(ligand_types) | set(rec_types))

    gpf_path = Path(f"{maps_prefix}.gpf")
    with open(gpf_path, "w") as f:
        f.write(f"npts {npts_x} {npts_y} {npts_z}\n")
        f.write(f"gridfld {maps_prefix.name}.maps.fld\n")
        f.write(f"spacing {spacing}\n")
        f.write(f"receptor_types {' '.join(rec_types)}\n")
        f.write(f"ligand_types {' '.join(ligand_types)}\n")
        f.write(f"receptor {rec_pdbqt.name}\n")
        f.write(f"gridcenter {center_x} {center_y} {center_z}\n")
        f.write("smooth 0.5\n")
        for atype in all_types:
            f.write(f"map {maps_prefix.name}.{atype}.map\n")
        f.write(f"elecmap {maps_prefix.name}.e.map\n")
        f.write(f"dsolvmap {maps_prefix.name}.d.map\n")
        f.write("dielectric -0.1465\n")

    glg_path = Path(f"{maps_prefix}.glg")
    print(f"  Running autogrid4 for {rec_pdbqt.name}...")
    result = subprocess.run(
        [str(autogrid_bin), "-p", str(gpf_path), "-l", str(glg_path)],
        capture_output=True, text=True, cwd=str(output_dir), timeout=300,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"autogrid4 failed (exit {result.returncode}) for {rec_pdbqt.name}:\n"
            f"{(result.stdout + result.stderr)[:500]}"
        )
    if not fld_path.exists():
        raise FileNotFoundError(f"autogrid4 did not produce {fld_path}")

    print(f"  ✓ ad4 maps generated: {maps_prefix.name}.maps.fld")
    return str(maps_prefix)


# ═══════════════════════════════════════════════════════════════════════════════
# Batch-mode docking
# ═══════════════════════════════════════════════════════════════════════════════

def _dock_batch_job(
    rec: dict,
    ligands: list[dict],
    *,
    scoring: str,
    vina_bin: Path,
    dock_dir: Path,
    log_dir: Path,
    exhaustiveness: int,
    num_modes: int | None,
    energy_range: int | None,
    cpu: int | None,
    seed: int | None,
    timeout: int,
    overwrite: bool,
    overwrite_error_log: bool,
    error_log: dict,
    failed_job_policy: str = "retry",
    batch_timeout: int = 0,
) -> list[dict]:
    rec_path = rec["pdbqt"]
    box_path = rec["box"]
    rec_name = _normalize_key(rec_path)

    results_out: list[dict] = []
    ligands_to_dock: list[dict] = []

    for lig in ligands:
        lig_path = lig["pdbqt"]
        lig_name = _normalize_key(lig_path)
        combo_name = f"{lig_name}__{rec_name}"
        stem = f"{rec_path.stem}__{lig_path.stem}_{scoring}_vina_out"
        out_path = dock_dir / f"{stem}.pdbqt"
        vina_log_path = log_dir / f"{stem}.log"

        if out_path.exists() and not overwrite:
            try:
                modes = validate_vina_output(out_path)
            except (OSError, UnicodeError, ValueError) as exc:
                quarantined = _quarantine_vina_output(
                    out_path, f"invalid committed output during resume: {exc}",
                )
                print(f"    ⚠ Invalid cached output quarantined: {quarantined}")
            else:
                r = DockingResult(
                    protein_name=rec_name, ligand_name=lig_name,
                    protein_path=rec_path, ligand_path=lig_path,
                    output_dir=dock_dir, status="skipped",
                    scoring_function=scoring,
                    num_poses=len(modes), best_affinity=modes[0]["affinity"],
                    pose_files=[out_path],
                )
                results_out.append({
                    "result": r,
                    "raw": {"receptor": str(rec_path), "ligand": str(lig_path),
                            "output": str(out_path), "log": str(vina_log_path),
                            "status": "skip-valid"},
                    "combo_name": combo_name, "ran_vina": False,
                })
                continue
        elif out_path.exists():
            try:
                validate_vina_output(out_path)
            except (OSError, UnicodeError, ValueError) as exc:
                _quarantine_vina_output(
                    out_path, f"invalid committed output before overwrite: {exc}",
                )
            else:
                _archive_superseded_vina_output(
                    out_path, "valid prior commit archived before overwrite attempt",
                )

        if (
            not overwrite_error_log
            and combo_name in error_log
            and failed_job_policy == "skip"
        ):
            previous_error = error_log[combo_name].get("error", "unknown")
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="failed",
                scoring_function=scoring,
                error_message=f"Cached prior failure (policy=skip): {previous_error[:180]}",
            )
            results_out.append({
                "result": r,
                "raw": {"receptor": str(rec_path), "ligand": str(lig_path),
                        "output": str(out_path), "log": str(vina_log_path),
                        "status": "cached-failure"},
                "combo_name": combo_name, "ran_vina": False,
            })
            continue

        ligands_to_dock.append(lig)

    if not ligands_to_dock:
        return results_out

    batch_out_dir = dock_dir / f"_batch_{rec_name}_{scoring}"
    batch_out_dir.mkdir(parents=True, exist_ok=True)
    batch_log_path = log_dir / f"batch_{rec_name}_{scoring}.log"

    # Batch scratch survives interrupted processes. Remove it from consideration
    # before launch so stale files can never masquerade as this invocation's
    # completed ligands.
    for lig in ligands_to_dock:
        stale_batch_output = batch_out_dir / f"{lig['pdbqt'].stem}_out.pdbqt"
        if stale_batch_output.exists():
            _quarantine_vina_output(
                stale_batch_output, "stale batch scratch before a new Vina invocation",
            )

    cmd: list[str] = [str(vina_bin)]

    if scoring == "ad4":
        maps_prefix = rec.get("maps_prefix", "")
        cmd += ["--maps", str(maps_prefix), "--scoring", "ad4"]
    else:
        cmd += ["--receptor", str(rec_path), "--config", str(box_path)]
        if scoring == "vinardo":
            cmd += ["--scoring", "vinardo"]

    for lig in ligands_to_dock:
        cmd += ["--batch", str(lig["pdbqt"])]

    cmd += ["--dir", str(batch_out_dir)]
    cmd += ["--exhaustiveness", str(exhaustiveness)]
    if num_modes is not None:
        cmd += ["--num_modes", str(num_modes)]
    if energy_range is not None:
        cmd += ["--energy_range", str(energy_range)]
    if cpu is not None:
        cmd += ["--cpu", str(cpu)]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    effective_timeout = batch_timeout if batch_timeout > 0 else (
        (timeout * len(ligands_to_dock)) if timeout > 0 else None
    )
    t0 = time.time()
    timed_out = False
    retcode: int | None = None  # stays None if open()/Popen() raises before wait()
    try:
        with open(batch_log_path, "w") as log_f:
            proc = subprocess.Popen(cmd, stdout=log_f, stderr=subprocess.STDOUT)
            try:
                retcode = proc.wait(timeout=effective_timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                retcode = -1
                timed_out = True
        batch_output = batch_log_path.read_text()
        batch_ok = retcode == 0
    except Exception as exc:
        batch_ok = False
        batch_output = str(exc)
        timed_out = False
    total_elapsed = time.time() - t0

    # On timeout: identify undocked ligands and split-retry
    if timed_out and len(ligands_to_dock) > 1:
        print(f"    ⏱ Batch timed out after {total_elapsed:.0f}s — splitting into sub-batches")

        completed_modes: Dict[str, List[dict]] = {}
        for lig in ligands_to_dock:
            lig_path = lig["pdbqt"]
            batch_out_file = batch_out_dir / f"{lig_path.stem}_out.pdbqt"
            if not batch_out_file.exists():
                continue
            try:
                completed_modes[lig_path.stem] = validate_vina_output(batch_out_file)
            except (OSError, UnicodeError, ValueError) as exc:
                _quarantine_vina_output(
                    batch_out_file, f"invalid/partial batch output after timeout: {exc}",
                )

        finished_ligs = [l for l in ligands_to_dock if l["pdbqt"].stem in completed_modes]
        remaining_ligs = [l for l in ligands_to_dock if l["pdbqt"].stem not in completed_modes]
        per_lig_time = total_elapsed / max(len(ligands_to_dock), 1)

        for lig in finished_ligs:
            lig_path = lig["pdbqt"]
            lig_name = _normalize_key(lig_path)
            combo_name = f"{lig_name}__{rec_name}"
            batch_out_file = batch_out_dir / f"{lig_path.stem}_out.pdbqt"
            stem = f"{rec_path.stem}__{lig_path.stem}_{scoring}_vina_out"
            final_out = dock_dir / f"{stem}.pdbqt"
            per_lig_log = log_dir / f"{stem}.log"

            os.replace(batch_out_file, final_out)
            per_lig_log.write_text(
                f"# Extracted from batch run (partial, timeout): {batch_log_path.name}\n"
                f"# Receptor: {rec_name}  Ligand: {lig_name}  Scoring: {scoring}\n\n"
            )
            modes = completed_modes[lig_path.stem]
            n_poses = len(modes)
            best_aff = modes[0]["affinity"] if modes else None

            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="success",
                scoring_function=scoring,
                num_poses=n_poses, best_affinity=best_aff,
                pose_files=[final_out],
                elapsed_time=per_lig_time,
            )
            raw = {"receptor": str(rec_path), "ligand": str(lig_path),
                   "output": str(final_out), "log": str(per_lig_log), "status": "ok"}
            results_out.append({"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": True})

        if remaining_ligs:
            n_sub = min(4, len(remaining_ligs))
            sub_size = max(1, len(remaining_ligs) // n_sub)
            sub_batches = [remaining_ligs[i:i + sub_size] for i in range(0, len(remaining_ligs), sub_size)]
            # We only reach this split-retry after a real timeout fired, so
            # effective_timeout is guaranteed to be a positive number. Base the
            # sub-batch timeout on it and never fall back to 0 — a 0/None timeout
            # would make the recursive proc.wait() block forever on a stuck ligand.
            sub_timeout = max(effective_timeout, timeout) if timeout > 0 else effective_timeout
            print(f"    Retrying {len(remaining_ligs)} remaining ligands in {len(sub_batches)} sub-batches")

            common_kwargs = dict(
                scoring=scoring, vina_bin=vina_bin, dock_dir=dock_dir,
                log_dir=log_dir, exhaustiveness=exhaustiveness,
                num_modes=num_modes, energy_range=energy_range, cpu=cpu,
                seed=seed, timeout=timeout, overwrite=overwrite,
                overwrite_error_log=overwrite_error_log, error_log=error_log,
                failed_job_policy=failed_job_policy,
                batch_timeout=sub_timeout,
            )

            for sb_idx, sub_batch in enumerate(sub_batches, 1):
                print(f"    Sub-batch {sb_idx}/{len(sub_batches)}: {len(sub_batch)} ligands")
                sub_results = _dock_batch_job(rec, sub_batch, **common_kwargs)
                results_out.extend(sub_results)

        return results_out
    total_elapsed = time.time() - t0

    per_lig_time = total_elapsed / max(len(ligands_to_dock), 1)
    for lig in ligands_to_dock:
        lig_path = lig["pdbqt"]
        lig_name = _normalize_key(lig_path)
        combo_name = f"{lig_name}__{rec_name}"

        batch_out_file = batch_out_dir / f"{lig_path.stem}_out.pdbqt"
        stem = f"{rec_path.stem}__{lig_path.stem}_{scoring}_vina_out"
        final_out = dock_dir / f"{stem}.pdbqt"
        per_lig_log = log_dir / f"{stem}.log"

        modes: List[dict] = []
        validation_error = ""
        if batch_out_file.exists():
            try:
                modes = validate_vina_output(batch_out_file)
            except (OSError, UnicodeError, ValueError) as exc:
                validation_error = str(exc)
                _quarantine_vina_output(
                    batch_out_file, f"invalid/partial batch output: {exc}",
                )

        if modes:
            os.replace(batch_out_file, final_out)
            per_lig_log.write_text(
                f"# Extracted from batch run: {batch_log_path.name}\n"
                f"# Receptor: {rec_name}  Ligand: {lig_name}  Scoring: {scoring}\n\n"
            )
            # Do NOT fall back to parsing batch_output here: it is the shared
            # batch log covering every ligand, so parse_vina_affinities would
            # return all ligands' modes and attribute them to this one. The
            # per-ligand pose file's REMARK VINA RESULT records are authoritative.
            n_poses = len(modes)
            best_aff = modes[0]["affinity"] if modes else None

            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="success",
                scoring_function=scoring,
                num_poses=n_poses, best_affinity=best_aff,
                pose_files=[final_out],
                elapsed_time=per_lig_time,
            )
            raw = {"receptor": str(rec_path), "ligand": str(lig_path),
                   "output": str(final_out), "log": str(per_lig_log), "status": "ok"}
        else:
            if validation_error:
                err_msg = f"Invalid/partial Vina output: {validation_error}"
            elif batch_ok:
                err_msg = "No output produced"
            elif retcode is not None:
                err_msg = f"Batch docking failed (rc={retcode})"
            else:
                err_msg = f"Batch docking failed to launch: {batch_output[:200]}"
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="failed",
                scoring_function=scoring,
                error_message=err_msg[:300],
                elapsed_time=per_lig_time,
            )
            raw = {"receptor": str(rec_path), "ligand": str(lig_path),
                   "output": str(final_out), "log": str(per_lig_log), "status": "fail"}

        results_out.append({"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": True})

    return results_out


# ═══════════════════════════════════════════════════════════════════════════════
# Single-job docking
# ═══════════════════════════════════════════════════════════════════════════════

def _dock_one_job(
    rec: dict,
    lig: dict,
    *,
    scoring: str = "vina",
    vina_bin: Path,
    dock_dir: Path,
    log_dir: Path,
    exhaustiveness: int,
    num_modes: int | None,
    energy_range: int | None,
    cpu: int | None,
    seed: int | None,
    timeout: int,
    overwrite: bool,
    overwrite_error_log: bool,
    error_log: dict,
    failed_job_policy: str = "retry",
) -> dict:
    rec_path = rec["pdbqt"]
    lig_path = lig["pdbqt"]
    box_path = rec["box"]
    rec_name = _normalize_key(rec_path)
    lig_name = _normalize_key(lig_path)
    combo_name = f"{lig_name}__{rec_name}"
    stem = f"{rec_path.stem}__{lig_path.stem}_{scoring}_vina_out"
    out_path = dock_dir / f"{stem}.pdbqt"
    vina_log_path = log_dir / f"{stem}.log"

    if out_path.exists() and not overwrite:
        try:
            modes = validate_vina_output(out_path)
        except (OSError, UnicodeError, ValueError) as exc:
            quarantined = _quarantine_vina_output(
                out_path, f"invalid committed output during resume: {exc}",
            )
            print(f"    ⚠ Invalid cached output quarantined: {quarantined}")
        else:
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="skipped",
                scoring_function=scoring,
                num_poses=len(modes), best_affinity=modes[0]["affinity"],
                pose_files=[out_path],
            )
            raw = {
                "receptor": str(rec_path), "ligand": str(lig_path),
                "output": str(out_path), "log": str(vina_log_path),
                "status": "skip-valid",
            }
            return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": False}
    elif out_path.exists():
        try:
            validate_vina_output(out_path)
        except (OSError, UnicodeError, ValueError) as exc:
            _quarantine_vina_output(
                out_path, f"invalid committed output before overwrite: {exc}",
            )
        else:
            _archive_superseded_vina_output(
                out_path, "valid prior commit archived before overwrite attempt",
            )

    if (
        not overwrite_error_log
        and combo_name in error_log
        and failed_job_policy == "skip"
    ):
        previous_error = error_log[combo_name].get("error", "unknown")
        r = DockingResult(
            protein_name=rec_name, ligand_name=lig_name,
            protein_path=rec_path, ligand_path=lig_path,
            output_dir=dock_dir, status="failed",
            scoring_function=scoring,
            error_message=f"Cached prior failure (policy=skip): {previous_error[:180]}",
        )
        raw = {
            "receptor": str(rec_path), "ligand": str(lig_path),
            "output": str(out_path), "log": str(vina_log_path),
            "status": "cached-failure",
        }
        return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": False}

    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{out_path.stem}.", suffix=".pdbqt", dir=str(dock_dir),
    )
    os.close(fd)
    tmp_out = Path(tmp_name)
    tmp_out.unlink(missing_ok=True)

    cmd = [str(vina_bin)]

    if scoring == "ad4":
        maps_prefix = rec.get("maps_prefix", "")
        cmd += ["--maps", str(maps_prefix), "--scoring", "ad4"]
    else:
        cmd += ["--receptor", str(rec_path), "--config", str(box_path)]
        if scoring == "vinardo":
            cmd += ["--scoring", "vinardo"]

    cmd += [
        "--ligand", str(lig_path),
        "--exhaustiveness", str(exhaustiveness),
        "--out", str(tmp_out),
    ]
    if num_modes is not None:
        cmd += ["--num_modes", str(num_modes)]
    if energy_range is not None:
        cmd += ["--energy_range", str(energy_range)]
    if cpu is not None:
        cmd += ["--cpu", str(cpu)]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    effective_timeout = timeout if timeout > 0 else None
    t0 = time.time()
    try:
        with open(vina_log_path, "w") as log_f:
            proc = subprocess.Popen(cmd, stdout=log_f, stderr=subprocess.STDOUT)
            try:
                retcode = proc.wait(timeout=effective_timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()
                if tmp_out.exists():
                    _quarantine_vina_output(
                        tmp_out, f"Vina subprocess timed out after {timeout}s",
                    )
                elapsed = time.time() - t0
                r = DockingResult(
                    protein_name=rec_name, ligand_name=lig_name,
                    protein_path=rec_path, ligand_path=lig_path,
                    output_dir=dock_dir, status="failed",
                    scoring_function=scoring,
                    error_message=f"Timeout ({timeout}s)",
                    elapsed_time=elapsed,
                )
                raw = {
                    "receptor": str(rec_path), "ligand": str(lig_path),
                    "output": str(out_path), "log": str(vina_log_path),
                    "status": "fail",
                }
                return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": True}
        output = vina_log_path.read_text()
        ok = retcode == 0
        if not ok:
            output = output.strip()
    except Exception as exc:
        ok = False
        output = str(exc)

    elapsed = time.time() - t0
    modes: List[dict] = []
    validation_error = ""
    if ok:
        try:
            modes = validate_vina_output(tmp_out)
        except (OSError, UnicodeError, ValueError) as exc:
            validation_error = str(exc)
            ok = False
    if ok and modes:
        os.replace(tmp_out, out_path)
        status_str = "ok"
    else:
        status_str = "fail"
        if tmp_out.exists():
            reason = (
                f"invalid/partial Vina output: {validation_error}"
                if validation_error else f"Vina invocation failed: {(output or 'unknown error')[:200]}"
            )
            _quarantine_vina_output(tmp_out, reason)

    n_poses = len(modes)
    best_aff = modes[0]["affinity"] if modes else None
    error_message = ""
    if status_str != "ok":
        error_message = (
            f"Invalid/partial Vina output: {validation_error}"
            if validation_error else (output or "Unknown error")[:300]
        )

    r = DockingResult(
        protein_name=rec_name, ligand_name=lig_name,
        protein_path=rec_path, ligand_path=lig_path,
        output_dir=dock_dir,
        status="success" if status_str == "ok" else "failed",
        scoring_function=scoring,
        num_poses=n_poses,
        best_affinity=best_aff,
        pose_files=[out_path] if status_str == "ok" else [],
        error_message=error_message,
        elapsed_time=elapsed,
    )
    raw = {
        "receptor": str(rec_path), "ligand": str(lig_path),
        "output": str(out_path), "log": str(vina_log_path),
        "status": status_str,
    }
    return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": True}


# ═══════════════════════════════════════════════════════════════════════════════
# Main parallel docking orchestrator
# ═══════════════════════════════════════════════════════════════════════════════

def run_autodock_vina(
    base_dir: Path | str,
    log_dir: Path | str,
    prepared_manifest: dict,
    cfg: dict,
    cpu_model: str,
    protein_workflow_data: dict | None = None,
    ligand_workflow_data: dict | None = None,
    *,
    lig_props_cache: Dict[Path, Dict] | None = None,
    prot_props_cache: Dict[Path, Dict] | None = None,
    progress_counters: dict | None = None,
) -> Tuple[pd.DataFrame, List[DockingResult]]:
    """Run AutoDock Vina for all receptor-ligand combinations."""
    base_dir = Path(base_dir).resolve()
    log_dir = Path(log_dir).resolve()
    vina_bin = Path(cfg["vina_bin"]).expanduser().resolve()

    scoring = cfg.get("scoring_function", "vina")
    batch_mode = cfg.get("batch_mode")
    batch_size = cfg.get("batch_size", 0)
    batch_timeout = cfg.get("batch_timeout", 0)
    exhaustiveness = cfg.get("exhaustiveness", 32)
    num_modes = cfg.get("num_modes", 10)
    energy_range = cfg.get("energy_range", 3)
    cpu = cfg.get("cpus_per_worker")
    seed = cfg.get("seed")
    overwrite = cfg.get("overwrite_existing", False) or cfg.get("overwrite_poses", False)
    overwrite_error_log = cfg.get("overwrite_error_log", False)
    timeout = cfg.get("timeout_per_complex", 0)
    max_workers = cfg.get("max_workers", 1)
    failed_job_policy = _failed_job_policy(cfg)

    dock_dir = base_dir / "docking"
    dock_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    _rec_root = cfg.get("receptors_dir")
    receptors = _collect_receptors(
        base_dir, prepared_manifest, protein_workflow_data,
        receptor_dirs=([Path(_rec_root).expanduser() / "pdbqt"] if _rec_root else []))
    ligands = _collect_ligands(base_dir, prepared_manifest, ligand_workflow_data)
    print(f"Docking receptors: {len(receptors)} | ligands: {len(ligands)} | scoring: {scoring}")

    vina_ok = vina_bin.is_file() and os.access(vina_bin, os.X_OK)
    if not vina_ok:
        print(f"WARNING: Vina executable not found or not executable ('{vina_bin}').")
    if not receptors:
        print("WARNING: No receptor PDBQT + box.txt pairs found; aborting docking step.")
    if not ligands:
        print("WARNING: No ligand PDBQT files found; aborting docking step.")

    raw_results: list[dict] = []
    docking_results: List[DockingResult] = []
    if not vina_ok or not receptors or not ligands:
        return pd.DataFrame(raw_results), docking_results

    # ad4 scoring: generate maps
    if scoring == "ad4":
        autogrid_bin = Path(cfg.get("autogrid_bin", "/usr/local/bin/autogrid4")).expanduser().resolve()
        if not (autogrid_bin.is_file() and os.access(autogrid_bin, os.X_OK)):
            print(f"ERROR: autogrid4 not found or not executable ('{autogrid_bin}').")
            return pd.DataFrame(raw_results), docking_results
        maps_dir = dock_dir / "ad4_maps"
        for rec in receptors:
            try:
                maps_prefix = _generate_ad4_maps(
                    rec_pdbqt=rec["pdbqt"], box_path=rec["box"],
                    output_dir=maps_dir, autogrid_bin=autogrid_bin,
                )
                rec["maps_prefix"] = maps_prefix
            except Exception as exc:
                print(f"  ✗ Failed to generate ad4 maps for {rec['name']}: {exc}")
                return pd.DataFrame(raw_results), docking_results

    use_batch = batch_mode
    if use_batch is None:
        use_batch = len(ligands) > 1
    if use_batch:
        print(f"Batch mode: ON ({len(receptors)} receptor(s), {len(ligands)} ligands per Vina invocation)")
    else:
        print("Batch mode: OFF (per-ligand docking)")

    effective_output_dir = base_dir
    error_log = load_error_log(effective_output_dir)
    if overwrite_error_log:
        error_log = {}
        save_error_log(effective_output_dir, error_log)

    csv_log_path, existing_combos = init_docking_log(effective_output_dir)
    print(f"Docking log: {csv_log_path} ({len(existing_combos)} existing entries)")

    jobs = list(itertools.product(receptors, ligands))
    n_total = len(jobs)
    print(f"Total docking jobs: {n_total}  |  Workers: {max_workers}")

    if progress_counters is not None:
        with progress_counters["lock"]:
            progress_counters["total"] = n_total

    start_all = time.time()
    log_lock = threading.Lock()

    def _process_result(job_out, job_idx):
        r = job_out["result"]
        combo_name = job_out["combo_name"]

        if r.status == "failed" and job_out["ran_vina"]:
            with log_lock:
                # error_log (in-memory) is authoritative; update it and persist
                # write-only. Avoids record_failure()'s per-failure full re-read
                # + rewrite (O(N^2) over a run) and the duplicate dict update.
                error_log[combo_name] = {
                    "error": r.error_message,
                    "timestamp": datetime.now().isoformat(),
                    "protein": r.protein_name,
                    "ligand": r.ligand_name,
                    "elapsed_time": round(r.elapsed_time, 2),
                }
                save_error_log(effective_output_dir, error_log)
        elif r.status in {"success", "skipped"} and combo_name in error_log:
            # A newly validated commit (including a valid cached pose) resolves
            # an older failure record. Leaving it behind would make a later
            # policy=skip run incorrectly suppress good work.
            with log_lock:
                error_log.pop(combo_name, None)
                save_error_log(effective_output_dir, error_log)

        icon = "✓" if r.status == "success" else ("⊘" if r.status == "skipped" else "✗")
        aff_str = f"{r.best_affinity:.2f} kcal/mol" if r.best_affinity is not None else "N/A"
        elapsed_str = f"{r.elapsed_time:.1f}s" if r.elapsed_time > 0 else "cached"
        print(f"  [{job_idx}/{n_total}] {icon} {r.protein_name} x {r.ligand_name} -> {r.status} "
              f"| {r.num_poses} poses | best: {aff_str} | {elapsed_str}")

        with log_lock:
            append_log_rows(
                [r], csv_log_path, existing_combos,
                lig_props_cache or {}, prot_props_cache or {},
                cfg, cpu_model,
                overwrite=overwrite or job_out["ran_vina"],
            )

        if progress_counters is not None:
            with progress_counters["lock"]:
                if r.status == "success":
                    progress_counters["done"] += 1
                elif r.status == "skipped":
                    progress_counters["skipped"] += 1
                else:
                    progress_counters["failed"] += 1
                progress_counters["running"] = max(0, progress_counters["running"] - 1)

        return r, job_out["raw"]

    worker_kwargs = dict(
        scoring=scoring,
        vina_bin=vina_bin,
        dock_dir=dock_dir,
        log_dir=log_dir,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        energy_range=energy_range,
        cpu=cpu,
        seed=seed,
        timeout=timeout,
        overwrite=overwrite,
        overwrite_error_log=overwrite_error_log,
        error_log=error_log,
        failed_job_policy=failed_job_policy,
    )

    if use_batch:
        for rec_idx, rec in enumerate(receptors, 1):
            if batch_size > 0:
                chunks = [ligands[i:i + batch_size] for i in range(0, len(ligands), batch_size)]
            else:
                chunks = [ligands]
            n_chunks = len(chunks)
            ligs_label = (f"{len(ligands)} ligands in {n_chunks} batch(es) of ≤{batch_size}"
                          if n_chunks > 1 else f"{len(ligands)} ligands")
            print(f"  Receptor {rec_idx}/{len(receptors)}: {rec['name']} — batch docking {ligs_label}")
            if progress_counters is not None:
                with progress_counters["lock"]:
                    progress_counters["running"] += len(ligands)
            for chunk_idx, chunk in enumerate(chunks, 1):
                if n_chunks > 1:
                    print(f"    Batch {chunk_idx}/{n_chunks}: {len(chunk)} ligands")
                try:
                    batch_results = _dock_batch_job(rec, chunk, **worker_kwargs, batch_timeout=batch_timeout)
                except Exception as exc:
                    print(f"    ✗ BATCH EXCEPTION for {rec['name']} (batch {chunk_idx}): {exc}")
                    continue
                for i, job_out in enumerate(batch_results):
                    job_idx = len(docking_results) + 1
                    r, raw = _process_result(job_out, job_idx)
                    docking_results.append(r)
                    raw_results.append(raw)

    elif max_workers > 1:
        futures = {}
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            for idx, (rec, lig) in enumerate(jobs, 1):
                if progress_counters is not None:
                    with progress_counters["lock"]:
                        progress_counters["running"] += 1
                fut = executor.submit(_dock_one_job, rec, lig, **worker_kwargs)
                futures[fut] = idx
            for fut in as_completed(futures):
                job_idx = futures[fut]
                try:
                    job_out = fut.result()
                except Exception as exc:
                    print(f"  [{job_idx}/{n_total}] ✗ WORKER EXCEPTION: {exc}")
                    continue
                r, raw = _process_result(job_out, job_idx)
                docking_results.append(r)
                raw_results.append(raw)
    else:
        for idx, (rec, lig) in enumerate(jobs, 1):
            if progress_counters is not None:
                with progress_counters["lock"]:
                    progress_counters["running"] += 1
            job_out = _dock_one_job(rec, lig, **worker_kwargs)
            r, raw = _process_result(job_out, idx)
            docking_results.append(r)
            raw_results.append(raw)

    total_time = time.time() - start_all
    n_success = sum(1 for r in docking_results if r.status == "success")
    n_failed = sum(1 for r in docking_results if r.status == "failed")
    n_skipped = sum(1 for r in docking_results if r.status == "skipped")

    print(f"\nDocking complete in {total_time:.1f}s ({total_time/60:.1f} min)")
    print(f"  Success: {n_success}  |  Failed: {n_failed}  |  Skipped: {n_skipped}")

    save_error_log(effective_output_dir, error_log)
    compact_docking_log(csv_log_path)

    results_df = pd.DataFrame(raw_results)
    return results_df, docking_results


# ═══════════════════════════════════════════════════════════════════════════════
# Post-pose optimization / gnina re-ranking (smina / gnina)
# ═══════════════════════════════════════════════════════════════════════════════
#
# AutoDock Vina ranks poses by its own empirical score (REMARK VINA RESULT; model
# 1 = best). This optional step keeps that ranking intact and, in parallel, locally
# minimises every Vina pose against the receptor with smina or gnina — boxed around
# the pose, so it never leaves the original pocket — and attaches a physics/CNN
# score. gnina additionally runs a CNN rescorer (CNNscore / CNNaffinity).
#
# Two rankings are kept side by side for every pose (satisfying "keep the original
# order but also re-rank on gnina"):
#   • autodock_rank  — Vina's original rank (1 = best), never reordered.
#   • optimized_rank — rank after sorting the same poses on the optimiser's
#                      minimised affinity (1 = most negative / best). This is the
#                      metric diffdock_gnina_rerank_analysis.py re-ranks on.
#
# Each pose is first rebuilt from its Vina PDBQT into a clean SDF (Meeko mk_export
# via embedded REMARK SMILES, an authoritative SDF+serial map for MGLTools, or a
# guarded Open Babel fallback) because smina/gnina cannot parse Meeko's macrocycle
# glue pseudo-atoms (CG0/G0) — see the macrocycle ring-opening fix. The original
# Vina PDBQT is never modified. Optimised poses are
# written to a per-complex optimized_<tool>/ subfolder (one SDF per pose, rank
# encoded in the filename) and both rankings + scores are recorded in
# optimization_log.csv — mirroring run_diffdock.optimize_results so the downstream
# reports treat AutoDock+gnina like the other tools' optimised variants.

# ``gnina`` is deliberately the historical/default CNN-rescore protocol.  Its
# name, output directory, log identity, and cache provenance must stay stable:
# existing benchmark trees contain many validated ``optimized_gnina`` caches.
# ``gnina_refinement`` is a distinct protocol identity that uses the same gnina
# executable but lets the CNN participate in geometry refinement.  It therefore
# gets a separate output/cache/log namespace.
OPTIMIZER_TOOLS = ("smina", "gnina", "gnina_refinement")
_LEGACY_ALL_OPTIMIZERS = ("smina", "gnina")
OPTIMIZER_CACHE_SCHEMA = 2
_TEMPLATE_RECONSTRUCTION_CONVERTER = "rdkit_template_map"
_TEMPLATE_RECONSTRUCTION_SCHEMA = "mgltools-template-coordinate-map-v1"
_TEMPLATE_MAPPING_TOLERANCE_A = 0.002


@dataclass(frozen=True)
class MGLPoseTemplate:
    """Authoritative ligand graph plus the MGLTools serial-coordinate map.

    MGLTools PDBQT files do not carry bond orders.  The atom serials do remain
    stable through Vina, however, so the prepared ligand's unchanged starting
    coordinates provide a deterministic bridge from the authoritative SDF atom
    indices to every docked pose.
    """

    ligand_name: str
    template_sdf: Path
    prepared_pdbqt: Path
    molecule: Chem.Mol
    atom_serial_map: Dict[int, int]
    prepared_atom_signature: Tuple[Tuple[int, str], ...]
    template_sha256: str
    prepared_pdbqt_sha256: str
    mapping_sha256: str

    def converter_inputs(self) -> dict:
        """Stable conversion inputs included in optimizer cache provenance."""
        return {
            "schema": _TEMPLATE_RECONSTRUCTION_SCHEMA,
            "template_sdf": str(self.template_sdf),
            "template_sdf_sha256": self.template_sha256,
            "prepared_pdbqt": str(self.prepared_pdbqt),
            "prepared_pdbqt_sha256": self.prepared_pdbqt_sha256,
            "atom_serial_mapping_sha256": self.mapping_sha256,
            "mapped_atom_count": len(self.atom_serial_map),
        }


@dataclass
class OptimizationSummary:
    """Machine-readable outcome returned by :func:`optimize_autodock_results`."""

    requested_tools: List[str] = field(default_factory=list)
    available_tools: List[str] = field(default_factory=list)
    target_complexes: int = 0
    selected_poses: int = 0
    attempted: int = 0
    optimized: int = 0
    reused: int = 0
    failed: int = 0
    pruned_outputs: int = 0
    status: str = "disabled"
    log_path: str = ""
    errors: List[str] = field(default_factory=list)
    time_by_tool_s: Dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            "requested_tools": list(self.requested_tools),
            "available_tools": list(self.available_tools),
            "target_complexes": self.target_complexes,
            "selected_poses": self.selected_poses,
            "attempted": self.attempted,
            "optimized": self.optimized,
            "reused": self.reused,
            "failed": self.failed,
            "pruned_outputs": self.pruned_outputs,
            "status": self.status,
            "log_path": self.log_path,
            "errors": list(self.errors),
            "time_by_tool_s": dict(self.time_by_tool_s),
        }


class OptimizationError(RuntimeError):
    """Fatal optimizer setup/run error with its partial structured summary."""

    def __init__(self, message: str, summary: OptimizationSummary):
        super().__init__(message)
        self.summary = summary

OPTIMIZATION_LOG_COLUMNS = [
    "timestamp", "combo_name", "protein_name", "ligand_name", "tool",
    "autodock_rank", "optimized_rank", "rank_metric", "vina_affinity",
    "minimized_affinity", "cnn_score", "cnn_affinity",
    "pose_file", "optimized_file", "status", "message", "elapsed_time_s",
    "converter", "converter_version", "optimizer_version",
    "source_scoring", "optimizer_scoring", "cnn_scoring", "cnn_model",
    "source_pose_sha256", "receptor_sha256", "provenance_fingerprint",
    "provenance_file", "processing_elapsed_time_s",
]

# How the optimised ranking is derived → (score column, higher_is_better). The
# default minimised affinity matches diffdock_gnina_rerank_analysis.py and is the
# only score smina produces; gnina additionally offers its CNN metrics.
_RERANK_METRICS = {
    "minimized_affinity": ("minimized_affinity", False),  # kcal/mol, more negative = better
    "cnn_affinity": ("cnn_affinity", True),               # predicted pK, higher = better
    "cnn_score": ("cnn_score", True),                     # pose-quality 0–1, higher = better
}

# SDF data tags written by smina/gnina, mapped to our log column names.
_SDF_TAG_RE = re.compile(r">\s*<([^>]+)>[^\r\n]*\r?\n([^\r\n]*)")
_SCORE_TAGS = {
    "minimizedAffinity": "minimized_affinity",
    "CNNscore": "cnn_score",
    "CNNaffinity": "cnn_affinity",
}

_ALLOWED_OPTIMIZER_SCORING = {
    "ad4_scoring", "default", "dkoes_fast", "dkoes_scoring",
    "dkoes_scoring_old", "vina", "vinardo",
}
_ALLOWED_CNN_SCORING = {
    "none", "rescore", "refinement", "metrorescore", "metrorefine", "all",
}
_OPTIMIZATION_LOG_LOCK = threading.Lock()
_FILE_HASH_CACHE: Dict[Tuple[str, int, int], str] = {}
_FILE_HASH_CACHE_LOCK = threading.Lock()
_EXECUTABLE_IDENTITY_CACHE: Dict[Tuple[str, str, int, int], dict] = {}
_MK_EXPORT_OUTPUT_FLAG_CACHE: Dict[str, str] = {}
_OPTIMIZER_PREFLIGHT_CACHE: Dict[str, dict] = {}


def _optimizer_base_tool(tool: str) -> str:
    """Executable family for a public optimizer protocol identity."""
    normalized = str(tool).strip().lower()
    if normalized == "gnina_refinement":
        return "gnina"
    if normalized in {"smina", "gnina"}:
        return normalized
    raise ValueError(f"unknown optimizer protocol {tool!r}")


def _is_gnina_optimizer(tool: str) -> bool:
    return _optimizer_base_tool(tool) == "gnina"


def _required_cnn_scoring(tool: str) -> Optional[str]:
    """Pinned CNN mode for a gnina protocol, else ``None`` for smina."""
    normalized = str(tool).strip().lower()
    if normalized == "gnina":
        return "rescore"
    if normalized == "gnina_refinement":
        return "refinement"
    _optimizer_base_tool(normalized)  # validate before returning
    return None


def resolve_optimizers(cfg: dict) -> List[str]:
    """Map the ``optimization`` config value to a list of tools to run."""
    mode = str(cfg.get("optimization", "none")).strip().lower()
    if mode in ("none", "off", "false", ""):
        return []
    if mode == "all":
        # Backward compatibility: ``all`` historically meant smina + gnina's
        # default rescore protocol.  CNN refinement is opt-in via the explicit
        # ``gnina_refinement`` value and is commonly run in a separate cell.
        return list(_LEGACY_ALL_OPTIMIZERS)
    if mode in OPTIMIZER_TOOLS:
        return [mode]
    raise ValueError(
        "optimization must be one of "
        f"none/smina/gnina/gnina_refinement/all (got {mode!r})"
    )


def _validate_optimizer_config(cfg: dict) -> None:
    """Validate optimizer controls early; never silently reinterpret typos."""
    tools = resolve_optimizers(cfg)
    if not tools:
        return
    search = str(cfg.get("optimize_search", "minimize")).strip().lower()
    if search not in {"minimize", "local_only"}:
        raise ValueError("optimize_search must be 'minimize' or 'local_only'")
    scoring = str(cfg.get("optimize_scoring", "default")).strip().lower()
    if scoring not in _ALLOWED_OPTIMIZER_SCORING:
        raise ValueError(
            "optimize_scoring must be a built-in smina/gnina scoring function "
            f"({', '.join(sorted(_ALLOWED_OPTIMIZER_SCORING))})"
        )
    rank_by = str(cfg.get("optimize_rank_by", "minimized_affinity")).strip().lower()
    if rank_by not in _RERANK_METRICS:
        raise ValueError(
            f"optimize_rank_by must be one of {', '.join(_RERANK_METRICS)}"
        )
    cnn_scoring = str(cfg.get("gnina_cnn_scoring", "rescore")).strip().lower()
    if cnn_scoring not in _ALLOWED_CNN_SCORING:
        raise ValueError(
            f"gnina_cnn_scoring must be one of {', '.join(sorted(_ALLOWED_CNN_SCORING))}"
        )
    for tool in tools:
        required_cnn_scoring = _required_cnn_scoring(tool)
        if required_cnn_scoring is not None and cnn_scoring != required_cnn_scoring:
            raise ValueError(
                f"optimization={tool!r} requires "
                f"gnina_cnn_scoring={required_cnn_scoring!r} "
                f"(got {cnn_scoring!r})"
            )
    if (
        any(_is_gnina_optimizer(tool) for tool in tools)
        and rank_by.startswith("cnn_")
        and cnn_scoring == "none"
    ):
        raise ValueError(
            f"optimize_rank_by={rank_by!r} requires gnina_cnn_scoring != 'none'"
        )
    if any(_is_gnina_optimizer(tool) for tool in tools) and bool(cfg.get("gnina_use_gpu", False)):
        lock_value = cfg.get("gpu_lock_file")
        if not isinstance(lock_value, (str, os.PathLike)) or not str(lock_value).strip():
            raise ValueError(
                "gnina_use_gpu=true requires a non-empty gpu_lock_file so "
                "concurrent GPU docking engines cannot overlap"
            )
        lock_path = Path(lock_value).expanduser()
        if lock_path.exists() and not lock_path.is_file():
            raise ValueError(f"gpu_lock_file is not a regular file: {lock_path}")
    try:
        top_n = int(cfg.get("optimize_top_n", 0) or 0)
    except (TypeError, ValueError) as e:
        raise ValueError("optimize_top_n must be a non-negative integer") from e
    if top_n < 0:
        raise ValueError("optimize_top_n must be a non-negative integer")
    try:
        autobox_add = float(cfg.get("optimize_autobox_add", 4.0))
        cpu = int(cfg.get("optimize_cpu", 1))
        timeout_s = float(cfg.get("optimize_timeout", 300))
        int(cfg.get("optimize_seed", cfg.get("seed", 0) or 0))
    except (TypeError, ValueError) as e:
        raise ValueError("optimizer numeric controls contain an invalid value") from e
    if not math.isfinite(autobox_add) or autobox_add <= 0:
        raise ValueError("optimize_autobox_add must be finite and > 0")
    if cpu < 1:
        raise ValueError("optimize_cpu must be >= 1")
    if not math.isfinite(timeout_s) or timeout_s <= 0:
        raise ValueError("optimize_timeout must be finite and > 0")


def _optimizer_settings(tool: str, cfg: dict) -> dict:
    """Normalized, explicit settings used in commands and provenance."""
    settings = {
        "search": str(cfg.get("optimize_search", "minimize")).strip().lower(),
        "empirical_scoring": str(cfg.get("optimize_scoring", "default")).strip().lower(),
        "autobox_add": float(cfg.get("optimize_autobox_add", 4.0)),
        "cpu": int(cfg.get("optimize_cpu", 1)),
        "seed": int(cfg.get("optimize_seed", cfg.get("seed", 0) or 0)),
        "timeout_s": float(cfg.get("optimize_timeout", 300)),
        "num_modes": 1,
    }
    if _is_gnina_optimizer(tool):
        settings.update({
            "cnn_scoring": str(cfg.get("gnina_cnn_scoring", "rescore")).strip().lower(),
            "cnn_model": str(
                cfg.get("gnina_cnn_model", "crossdock_default2018_ensemble")
            ).strip(),
            "cnn_model_file": str(cfg.get("gnina_cnn_model_file", "") or "").strip(),
            "use_gpu": bool(cfg.get("gnina_use_gpu", False)),
        })
    return settings


def _resolve_optimizer_exe(tool: str, cfg: dict) -> str:
    """Locate a tool's executable: explicit config path → PATH → ''."""
    base_tool = _optimizer_base_tool(tool)
    explicit = cfg.get(f"{base_tool}_executable")
    if explicit:
        path = os.path.expanduser(str(explicit))
        if Path(path).exists():
            return path
        return shutil.which(path) or ""
    return shutil.which(base_tool) or ""


def _gnina_subprocess_env(cfg: dict) -> dict:
    """Environment for the gnina subprocess, with CUDA/cuDNN libs on the path.

    The prebuilt gnina CUDA binary is dynamically linked against the CUDA runtime
    + cuDNN 9 shared objects (libcudart.so.12, libcudnn.so.9, …). The project's
    ``/home/manndo/docking_tools/gnina`` wrapper already seeds these from the
    diffdock env, so this is normally a no-op, but honour an explicit
    ``gnina_lib_dirs`` list and/or a ``gnina_python`` (or ``diffdock_python``)
    interpreter whose env root supplies the ``nvidia/*/lib`` dirs — matching
    run_diffdock so a bare ``gnina.bin`` also resolves its libs.
    """
    env = os.environ.copy()
    lib_dirs: List[str] = []

    for d in cfg.get("gnina_lib_dirs") or []:
        p = os.path.expanduser(str(d))
        if Path(p).is_dir():
            lib_dirs.append(p)

    python_hint = cfg.get("gnina_python") or cfg.get("diffdock_python")
    if python_hint:
        env_root = Path(os.path.expanduser(python_hint)).resolve().parent.parent
        lib_dirs.extend(
            str(p) for p in sorted(env_root.glob("lib/python*/site-packages/nvidia/*/lib"))
            if p.is_dir()
        )

    if lib_dirs:
        existing = env.get("LD_LIBRARY_PATH", "")
        env["LD_LIBRARY_PATH"] = ":".join(lib_dirs + ([existing] if existing else []))
    return env


@contextlib.contextmanager
def _interprocess_gpu_lock(lock_file: Path | str):
    """Serialize GPU engine subprocesses with a kernel-released advisory lock."""
    if fcntl is None:  # pragma: no cover - benchmark host is Linux
        raise RuntimeError("gpu_lock_file requires fcntl/flock support")
    lock_path = Path(lock_file).expanduser().resolve()
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a+") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            yield lock_path
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def _probe_optimizer(tool: str, exe: str, cfg: dict) -> Tuple[bool, str]:
    """Run ``<exe> --version`` to confirm the binary actually loads (used by
    --dry-run to catch missing shared libraries before a real run)."""
    env = _gnina_subprocess_env(cfg) if _is_gnina_optimizer(tool) else None
    try:
        proc = subprocess.run([exe, "--version"], capture_output=True, text=True,
                              timeout=30, env=env)
    except Exception as e:
        return False, str(e)
    lines = (proc.stdout or proc.stderr or "").strip().splitlines()
    if proc.returncode != 0:
        return False, lines[-1] if lines else f"rc={proc.returncode}"
    return True, lines[0] if lines else "ok"


def _sha256_file(path: Path) -> str:
    resolved = Path(path).resolve()
    stat = resolved.stat()
    cache_key = (str(resolved), stat.st_size, stat.st_mtime_ns)
    with _FILE_HASH_CACHE_LOCK:
        cached = _FILE_HASH_CACHE.get(cache_key)
    if cached is not None:
        return cached
    digest = hashlib.sha256()
    with open(resolved, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    value = digest.hexdigest()
    with _FILE_HASH_CACHE_LOCK:
        # Drop stale identities for the same path when size/mtime changes.
        for key in [key for key in _FILE_HASH_CACHE if key[0] == str(resolved) and key != cache_key]:
            _FILE_HASH_CACHE.pop(key, None)
        _FILE_HASH_CACHE[cache_key] = value
    return value


def _atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent),
    )
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, path)
    finally:
        with contextlib.suppress(FileNotFoundError):
            Path(tmp_name).unlink()


def _converter_version(name: str, exe: str) -> str:
    """Return a stable converter version string for audit/cache provenance."""
    if name == "mk_export":
        # Query the interpreter in the executable's shebang. Importing Meeko in
        # this process would report the notebook/base environment, which can be
        # different from an explicitly pinned mk_export executable.
        with contextlib.suppress(OSError, IndexError):
            with open(exe, "r", errors="ignore") as handle:
                first_line = handle.readline().rstrip("\r\n")
            if first_line.startswith("#!"):
                interpreter = first_line[2:].strip().split()[0]
                proc = subprocess.run(
                    [
                        interpreter, "-c",
                        "import importlib.metadata as m; print('meeko ' + m.version('meeko'))",
                    ],
                    capture_output=True, text=True, timeout=30,
                )
                if proc.returncode == 0 and proc.stdout.strip():
                    return proc.stdout.strip().splitlines()[0]
        args = [exe, "--help"]
    else:
        args = [exe, "-V"]
    try:
        proc = subprocess.run(args, capture_output=True, text=True, timeout=30)
        lines = (proc.stdout or proc.stderr or "").strip().splitlines()
        if lines:
            return lines[0].strip()
    except Exception:
        pass
    return "unknown"


def _mk_export_output_flag(exe: str) -> str:
    """Resolve Meeko's version-specific explicit SDF output option.

    Meeko 0.5 uses ``-o/--output_filename``; 0.7 uses
    ``-s/--write_sdf``. This deliberately inspects help instead of assuming
    that the incompatible short option has one meaning across releases.
    """
    resolved = str(Path(exe).resolve())
    cached = _MK_EXPORT_OUTPUT_FLAG_CACHE.get(resolved)
    if cached:
        return cached
    proc = subprocess.run([exe, "--help"], capture_output=True, text=True, timeout=30)
    help_text = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(f"mk_export --help failed with rc={proc.returncode}")
    if "--output_filename" in help_text:
        flag = "-o"
    elif "--write_sdf" in help_text:
        flag = "-s"
    else:
        raise RuntimeError("unsupported mk_export CLI: no explicit SDF output option")
    _MK_EXPORT_OUTPUT_FLAG_CACHE[resolved] = flag
    return flag


def _executable_identity(exe: str, version: str) -> dict:
    path = Path(exe).resolve()
    try:
        stat = path.stat()
        cache_key = (str(path), version, stat.st_size, stat.st_mtime_ns)
    except OSError:
        return {"path": str(path), "version": version}
    cached = _EXECUTABLE_IDENTITY_CACHE.get(cache_key)
    if cached is None:
        cached = {
            "path": str(path), "version": version,
            "size": stat.st_size, "mtime_ns": stat.st_mtime_ns,
        }
        _EXECUTABLE_IDENTITY_CACHE[cache_key] = cached
    return dict(cached)


def _parse_sdf_scores(sdf_path: Path) -> Dict[str, float]:
    scores: Dict[str, float] = {}
    try:
        text = Path(sdf_path).read_text()
    except Exception:
        return scores
    for m in _SDF_TAG_RE.finditer(text):
        tag = m.group(1).strip()
        col = _SCORE_TAGS.get(tag)
        if col and col not in scores:
            try:
                value = float(m.group(2).strip())
                if math.isfinite(value):
                    scores[col] = value
            except ValueError:
                pass
    return scores


def _read_single_sdf_molecule(sdf_path: Path) -> Tuple[Optional[Chem.Mol], str]:
    """Read exactly one valid molecule from an SDF, rejecting junk/partials."""
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        return None, "missing or empty SDF"
    try:
        supplier = Chem.SDMolSupplier(
            str(sdf_path), removeHs=False, sanitize=True, strictParsing=True,
        )
        records = list(supplier)
    except Exception as e:
        return None, f"unreadable SDF: {e}"
    if len(records) != 1 or records[0] is None:
        return None, f"expected exactly one valid SDF molecule, found {len(records)} record(s)"
    mol = records[0]
    if any(atom.GetAtomicNum() == 0 for atom in mol.GetAtoms()):
        return None, "SDF contains dummy/pseudo atoms"
    try:
        heavy = Chem.RemoveHs(mol)
        if len(Chem.GetMolFrags(heavy)) != 1:
            return None, "SDF contains multiple disconnected fragments"
    except Exception as e:
        return None, f"invalid molecular graph: {e}"
    if heavy.GetNumAtoms() == 0:
        return None, "SDF contains no heavy atoms"
    if mol.GetNumConformers() != 1 or not mol.GetConformer().Is3D():
        return None, "SDF does not contain one 3D conformer"
    return mol, "ok"


def _molecule_identity(mol: Chem.Mol) -> dict:
    """Canonical identity plus index-stable topology for atom-map validation."""
    heavy = Chem.RemoveHs(mol)
    Chem.SanitizeMol(heavy)
    atoms = [
        [
            atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetIsotope(),
            bool(atom.GetIsAromatic()), int(atom.GetChiralTag()),
        ]
        for atom in heavy.GetAtoms()
    ]
    bonds = sorted(
        [
            min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            str(bond.GetBondType()), bool(bond.GetIsAromatic()), int(bond.GetStereo()),
        ]
        for bond in heavy.GetBonds()
    )
    return {
        "canonical_isomeric_smiles": Chem.MolToSmiles(
            heavy, canonical=True, isomericSmiles=True,
        ),
        "formula": rdMolDescriptors.CalcMolFormula(heavy),
        "heavy_atom_count": heavy.GetNumAtoms(),
        "atoms_by_index": atoms,
        "bonds_by_index": bonds,
    }


def _required_score_columns(tool: str, cfg: dict) -> Set[str]:
    required = {"minimized_affinity"}
    if _is_gnina_optimizer(tool) and _optimizer_settings(tool, cfg)["cnn_scoring"] != "none":
        required.update({"cnn_score", "cnn_affinity"})
    return required


def _validate_optimizer_sdf(
    sdf_path: Path,
    expected_identity: dict,
    required_scores: Set[str],
) -> Tuple[bool, Dict[str, float], str]:
    mol, message = _read_single_sdf_molecule(sdf_path)
    if mol is None:
        return False, {}, message
    try:
        identity = _molecule_identity(mol)
    except Exception as e:
        return False, {}, f"cannot determine output molecule identity: {e}"
    invariant_fields = ("canonical_isomeric_smiles", "formula", "heavy_atom_count")
    if any(identity.get(field) != expected_identity.get(field) for field in invariant_fields):
        return False, {}, "optimizer changed molecule identity/topology"
    # GNINA legitimately reorders atoms in its output SDF. Validate an explicit
    # full graph+stereochemistry isomorphism instead of requiring equal indexes.
    expected_mol = Chem.MolFromSmiles(expected_identity["canonical_isomeric_smiles"])
    output_heavy = Chem.RemoveHs(mol)
    if expected_mol is None:
        return False, {}, "cannot reconstruct expected molecule for atom mapping"
    mapping = output_heavy.GetSubstructMatch(expected_mol, useChirality=True)
    reverse_mapping = expected_mol.GetSubstructMatch(output_heavy, useChirality=True)
    if len(mapping) != expected_mol.GetNumAtoms() or len(reverse_mapping) != output_heavy.GetNumAtoms():
        return False, {}, "optimizer output has no complete topology-preserving atom mapping"
    scores = _parse_sdf_scores(sdf_path)
    missing = sorted(required_scores - scores.keys())
    if missing:
        return False, scores, f"optimizer SDF missing required score tag(s): {', '.join(missing)}"
    return True, scores, "ok"


def _validate_pdbqt_pose_body(body: str, label: str) -> None:
    """Reject structurally truncated ligand torsion trees."""
    lines = body.splitlines()
    if "\x00" in body:
        raise ValueError(f"{label} contains NUL bytes")
    if not re.search(r"^(?:ATOM|HETATM)\s", body, flags=re.MULTILINE):
        raise ValueError(f"{label} contains no atoms")
    root_positions = [idx for idx, line in enumerate(lines) if line.startswith("ROOT")]
    endroot_positions = [idx for idx, line in enumerate(lines) if line.startswith("ENDROOT")]
    if len(root_positions) != 1 or len(endroot_positions) != 1:
        raise ValueError(f"{label} has an incomplete ROOT…ENDROOT block")
    if root_positions[0] >= endroot_positions[0]:
        raise ValueError(f"{label} has ENDROOT before ROOT")
    branch_stack: List[Tuple[int, int]] = []
    for line in lines:
        if line.startswith("BRANCH"):
            fields = line.split()
            if len(fields) < 3:
                raise ValueError(f"{label} has a malformed BRANCH record")
            branch_stack.append((int(fields[1]), int(fields[2])))
        elif line.startswith("ENDBRANCH"):
            fields = line.split()
            if len(fields) < 3 or not branch_stack:
                raise ValueError(f"{label} has an unmatched ENDBRANCH record")
            pair = (int(fields[1]), int(fields[2]))
            if branch_stack.pop() != pair:
                raise ValueError(f"{label} has a mismatched ENDBRANCH record")
    if branch_stack:
        raise ValueError(f"{label} has an unclosed BRANCH record")
    torsdof = [line for line in lines if line.startswith("TORSDOF")]
    if len(torsdof) != 1 or not re.fullmatch(r"TORSDOF\s+\d+\s*", torsdof[0]):
        raise ValueError(f"{label} has no single valid TORSDOF terminator")


def _split_pdbqt_models(pdbqt_path: Path, out_dir: Path) -> List[Tuple[int, Path, Optional[float]]]:
    """Split a Vina output PDBQT into one single-model PDBQT per pose.

    Returns ``[(vina_rank, pose_pdbqt_path, vina_affinity), …]``. Vina writes its
    poses best-first, so the MODEL number is the Vina rank (1 = best). Each written
    file drops the MODEL/ENDMDL wrapper but keeps the ROOT/BRANCH torsion tree and
    the ``REMARK VINA RESULT`` line, so smina/gnina read it as a single ligand.
    """
    text = Path(pdbqt_path).read_text(errors="strict")
    if not text.strip():
        raise ValueError(f"empty PDBQT output: {pdbqt_path}")
    stem = Path(pdbqt_path).stem
    parsed: List[Tuple[int, List[str], float]] = []
    lines = text.splitlines(keepends=True)
    has_wrappers = any(line.startswith(("MODEL", "ENDMDL")) for line in lines)

    if has_wrappers:
        cur_lines: Optional[List[str]] = None
        cur_rank: Optional[int] = None
        seen_ranks: Set[int] = set()
        for line_no, line in enumerate(lines, 1):
            if line.startswith("MODEL"):
                if cur_lines is not None:
                    raise ValueError(
                        f"nested MODEL at line {line_no}; previous model is truncated"
                    )
                match = re.fullmatch(r"MODEL\s+(\d+)\s*\r?\n?", line)
                if not match:
                    raise ValueError(f"malformed MODEL record at line {line_no}")
                cur_rank = int(match.group(1))
                if cur_rank in seen_ranks:
                    raise ValueError(f"duplicate MODEL rank {cur_rank}")
                cur_lines = []
                continue
            if line.startswith("ENDMDL"):
                if cur_lines is None or cur_rank is None:
                    raise ValueError(f"ENDMDL without open MODEL at line {line_no}")
                body = "".join(cur_lines)
                _validate_pdbqt_pose_body(body, f"MODEL {cur_rank}")
                affinity_match = re.search(
                    r"^REMARK VINA RESULT:\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
                    body, flags=re.MULTILINE,
                )
                if not affinity_match:
                    raise ValueError(f"MODEL {cur_rank} has no valid Vina affinity")
                affinity = float(affinity_match.group(1))
                if not math.isfinite(affinity):
                    raise ValueError(f"MODEL {cur_rank} has non-finite Vina affinity")
                parsed.append((cur_rank, cur_lines, affinity))
                seen_ranks.add(cur_rank)
                cur_lines = None
                cur_rank = None
                continue
            if cur_lines is not None:
                cur_lines.append(line)
        if cur_lines is not None:
            raise ValueError(f"MODEL {cur_rank} is truncated (missing ENDMDL)")
        if not parsed:
            raise ValueError(f"no complete MODEL…ENDMDL poses in {pdbqt_path}")
    else:
        # Legacy single-pose files are accepted only when they contain atoms and
        # a finite Vina result. A partially written multi-model file is never
        # treated as this form because the leading MODEL is detected above.
        _validate_pdbqt_pose_body(text, "single-pose PDBQT")
        affinity_match = re.search(
            r"^REMARK VINA RESULT:\s*([-+]?\d+(?:\.\d+)?(?:[Ee][-+]?\d+)?)",
            text, flags=re.MULTILINE,
        )
        if not affinity_match:
            raise ValueError(f"single-pose PDBQT has no valid Vina affinity: {pdbqt_path}")
        affinity = float(affinity_match.group(1))
        if not math.isfinite(affinity):
            raise ValueError(f"single-pose PDBQT has non-finite Vina affinity: {pdbqt_path}")
        parsed.append((1, lines, affinity))

    # Validate the entire input before writing any split files, so a truncated
    # later model cannot leave a plausible partial pose set behind.
    out_dir.mkdir(parents=True, exist_ok=True)
    poses: List[Tuple[int, Path, Optional[float]]] = []
    for rank, body_lines, affinity in parsed:
        pose_path = out_dir / f"{stem}_rank{rank}.pdbqt"
        pose_path.write_text("".join(body_lines))
        poses.append((rank, pose_path, affinity))
    return poses


def _resolve_pose_converters(cfg: dict) -> Dict[str, str]:
    """Locate the PDBQT→SDF converters (Meeko mk_export, Open Babel) once.

    Explicit config path → PATH → ''. mk_export rebuilds the exact molecule from
    the Vina pose's embedded ``REMARK SMILES``; obabel is the coordinate-only
    fallback for non-Meeko PDBQTs.
    """
    def _find(explicit, *names) -> str:
        if explicit:
            p = os.path.expanduser(str(explicit))
            return p if Path(p).exists() else (shutil.which(p) or "")
        for n in names:
            hit = shutil.which(n)
            if hit:
                return hit
        return ""

    return {
        "mk_export": _find(cfg.get("mk_export_executable"), "mk_export.py", "mk_export"),
        "obabel": _find(cfg.get("obabel_executable"), "obabel"),
    }


def _remark_smiles(text: str) -> Optional[str]:
    match = re.search(r"^REMARK SMILES\s+(?!IDX\b)(\S.*)$", text, flags=re.MULTILINE)
    return match.group(1).strip() if match else None


def _remark_smiles_mapping(text: str) -> List[Tuple[int, int]]:
    values: List[int] = []
    for match in re.finditer(r"^REMARK SMILES IDX\s+(.+)$", text, flags=re.MULTILINE):
        try:
            values.extend(int(value) for value in match.group(1).split())
        except ValueError:
            return []
    if len(values) % 2:
        return []
    return list(zip(values[0::2], values[1::2]))


def _pdbqt_atom_coordinates(text: str) -> Dict[int, Tuple[float, float, float]]:
    coordinates: Dict[int, Tuple[float, float, float]] = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            coordinates[int(line[6:11])] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]),
            )
        except (ValueError, IndexError):
            continue
    return coordinates


# Elements AutoDock Vina has no atom type for. Meeko writes them into the PDBQT with
# a carbon type (so the force field can score them) but keeps the real element in the
# ATOM-NAME column. Boron is the case that matters here: the JKU 2-APB analogue is a
# boronate, and reading its type column as chemistry made every one of its poses
# unreconstructible.
_AUTODOCK_UNTYPEABLE_ELEMENTS = {"B", "Si", "Se", "Li", "Be", "Na", "K", "Ca", "Al"}


def _pdbqt_atom_elements(text: str) -> Dict[int, str]:
    """Return PDBQT serial→element, interpreting AutoDock atom types safely."""
    autodock_elements = {
        "A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N",
        "O": "O", "OA": "O", "OS": "O", "S": "S", "SA": "S",
        "H": "H", "HD": "H", "F": "F", "P": "P", "I": "I",
        "CL": "Cl", "BR": "Br",
    }
    periodic_table = Chem.GetPeriodicTable()
    elements: Dict[int, str] = {}
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except ValueError:
            continue
        atom_type = line[77:].strip().split()[0] if line[77:].strip() else ""
        letters = re.sub(r"[^A-Za-z]", "", atom_type)
        upper = letters.upper()
        if upper.startswith("CG"):
            upper = "C"  # macrocycle carbon/glue decoration (mapped carbons only)
        elif upper.startswith("G"):
            continue  # Meeko macrocycle glue pseudo-atoms have no element
        element = autodock_elements.get(upper, "")
        if not element and letters:
            # Preserve genuine two-letter elements (Zn, Mg, ...), then fall
            # back to the first letter for decorated AutoDock types.
            for candidate in (letters[0].upper() + letters[1:2].lower(), letters[0].upper()):
                try:
                    if periodic_table.GetAtomicNumber(candidate) > 0:
                        element = candidate
                        break
                except RuntimeError:
                    continue
        # AutoDock's type set cannot express every element. Meeko types boron (and
        # other untypeable atoms) as carbon for the force field while keeping the
        # true element in the ATOM-NAME column, so a type-derived element would
        # read B as C and the SMILES-to-PDBQT element check would reject a pose
        # that is in fact mapped correctly. Prefer the atom name when it names a
        # real element that AutoDock has no type for. The per-atom coordinate
        # identity check below is untouched and remains the guarantee of pose
        # fidelity — this only stops the type column from being read as chemistry
        # it was never able to carry.
        if element == "C":
            name_letters = re.sub(r"[^A-Za-z]", "", line[12:16])
            if name_letters:
                candidate = name_letters[0].upper() + name_letters[1:2].lower()
                for sym in (candidate, name_letters[0].upper()):
                    if sym in _AUTODOCK_UNTYPEABLE_ELEMENTS:
                        element = sym
                        break
        if element:
            elements[serial] = element
    return elements


def _template_reconstruction_settings(
    cfg: dict,
) -> Tuple[bool, bool, Optional[Path], Optional[Path]]:
    """Resolve and validate the optional authoritative-template configuration."""
    required = bool(cfg.get("optimize_require_template_reconstruction", False))
    template_value = cfg.get("optimize_ligand_template_root")
    prepared_value = cfg.get("optimize_prepared_ligand_pdbqt_dir")
    enabled = bool(required or template_value or prepared_value)
    if not enabled:
        return False, False, None, None
    if not template_value or not prepared_value:
        raise ValueError(
            "template reconstruction requires both optimize_ligand_template_root "
            "and optimize_prepared_ligand_pdbqt_dir"
        )
    template_root = Path(str(template_value)).expanduser().resolve()
    prepared_root = Path(str(prepared_value)).expanduser().resolve()
    if not template_root.is_dir():
        raise ValueError(f"ligand template root is not a directory: {template_root}")
    if not prepared_root.is_dir():
        raise ValueError(
            f"prepared ligand PDBQT directory is not a directory: {prepared_root}"
        )
    return True, required, template_root, prepared_root


def _pdbqt_atom_type_signature(text: str) -> Tuple[Tuple[int, str], ...]:
    """Return the ordered ``(serial, AutoDock type)`` ligand atom signature."""
    signature: List[Tuple[int, str]] = []
    seen: Set[int] = set()
    for line in text.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            serial = int(line[6:11])
        except (ValueError, IndexError) as e:
            raise ValueError("PDBQT contains a malformed atom serial") from e
        if serial in seen:
            raise ValueError(f"PDBQT contains duplicate atom serial {serial}")
        tail = line[77:].strip() if len(line) > 77 else ""
        atom_type = tail.split()[0] if tail else ""
        if not atom_type:
            fields = line.split()
            atom_type = fields[-1] if fields else ""
        if not atom_type:
            raise ValueError(f"PDBQT atom {serial} has no AutoDock type")
        signature.append((serial, atom_type.upper()))
        seen.add(serial)
    if not signature:
        raise ValueError("PDBQT contains no ligand atoms")
    return tuple(signature)


def _load_mgl_pose_template(
    ligand_name: str,
    template_sdf: Path,
    prepared_pdbqt: Path,
) -> MGLPoseTemplate:
    """Build a strict SDF-atom to prepared-PDBQT-serial coordinate map."""
    template_sdf = template_sdf.resolve()
    prepared_pdbqt = prepared_pdbqt.resolve()
    full_mol, message = _read_single_sdf_molecule(template_sdf)
    if full_mol is None:
        raise ValueError(f"invalid authoritative ligand SDF {template_sdf}: {message}")

    # Remove ordinary explicit hydrogens, but deliberately retain hydrogens
    # needed to encode double-bond stereochemistry (RDKit's default RemoveHs
    # behaviour). RemoveAllHs would erase the directional imine identity in
    # benchmark ligands such as 5SAK_ZRY and 8D5D_5DK.
    template_mol = Chem.RemoveHs(full_mol)
    if template_mol.GetNumConformers() != 1 or not template_mol.GetConformer().Is3D():
        raise ValueError(f"authoritative ligand template is not a single 3D conformer: {template_sdf}")
    try:
        _molecule_identity(template_mol)
    except Exception as e:
        raise ValueError(f"invalid authoritative ligand chemistry in {template_sdf}: {e}") from e

    try:
        prepared_text = prepared_pdbqt.read_text(errors="strict")
    except OSError as e:
        raise ValueError(f"cannot read prepared ligand PDBQT {prepared_pdbqt}: {e}") from e
    signature = _pdbqt_atom_type_signature(prepared_text)
    coordinates = _pdbqt_atom_coordinates(prepared_text)
    elements = _pdbqt_atom_elements(prepared_text)
    signature_serials = {serial for serial, _ in signature}
    if set(coordinates) != signature_serials:
        raise ValueError(f"prepared ligand PDBQT has invalid atom coordinates: {prepared_pdbqt}")
    if any(
        not all(math.isfinite(value) for value in xyz)
        for xyz in coordinates.values()
    ):
        raise ValueError(f"prepared ligand PDBQT has non-finite coordinates: {prepared_pdbqt}")
    if set(elements) != signature_serials:
        missing = sorted(signature_serials - set(elements))
        raise ValueError(
            f"prepared ligand PDBQT has unsupported atom types at serials {missing}: "
            f"{prepared_pdbqt}"
        )

    tolerance_sq = _TEMPLATE_MAPPING_TOLERANCE_A ** 2
    conformer = template_mol.GetConformer()
    atom_serial_map: Dict[int, int] = {}
    for atom in template_mol.GetAtoms():
        atom_idx = atom.GetIdx()
        point = conformer.GetAtomPosition(atom_idx)
        candidates: List[int] = []
        for serial, element in elements.items():
            if element != atom.GetSymbol():
                continue
            xyz = coordinates[serial]
            delta_sq = sum(
                (actual - expected) ** 2
                for actual, expected in zip((point.x, point.y, point.z), xyz)
            )
            if delta_sq <= tolerance_sq:
                candidates.append(serial)
        if len(candidates) != 1:
            raise ValueError(
                f"{ligand_name}: template atom {atom_idx + 1} ({atom.GetSymbol()}) "
                f"has {len(candidates)} prepared-PDBQT coordinate matches"
            )
        atom_serial_map[atom_idx] = candidates[0]

    mapped_serials = list(atom_serial_map.values())
    if len(set(mapped_serials)) != len(mapped_serials):
        raise ValueError(f"{ligand_name}: template-to-PDBQT atom map is not one-to-one")
    mapped_heavy = {
        serial for atom_idx, serial in atom_serial_map.items()
        if template_mol.GetAtomWithIdx(atom_idx).GetAtomicNum() != 1
    }
    prepared_heavy = {serial for serial, element in elements.items() if element != "H"}
    if mapped_heavy != prepared_heavy:
        raise ValueError(
            f"{ligand_name}: authoritative template and prepared PDBQT heavy atoms differ "
            f"(template-only={sorted(mapped_heavy - prepared_heavy)}, "
            f"PDBQT-only={sorted(prepared_heavy - mapped_heavy)})"
        )

    mapping_material = json.dumps(
        sorted(atom_serial_map.items()), separators=(",", ":"),
    ).encode()
    return MGLPoseTemplate(
        ligand_name=ligand_name,
        template_sdf=template_sdf,
        prepared_pdbqt=prepared_pdbqt,
        molecule=template_mol,
        atom_serial_map=atom_serial_map,
        prepared_atom_signature=signature,
        template_sha256=_sha256_file(template_sdf),
        prepared_pdbqt_sha256=_sha256_file(prepared_pdbqt),
        mapping_sha256=hashlib.sha256(mapping_material).hexdigest(),
    )


def _resolve_mgl_pose_template(ligand_name: str, cfg: dict) -> Optional[MGLPoseTemplate]:
    """Resolve one result's authoritative SDF and prepared MGLTools PDBQT."""
    enabled, required, template_root, prepared_root = (
        _template_reconstruction_settings(cfg)
    )
    if not enabled or template_root is None or prepared_root is None:
        return None
    if Path(ligand_name).name != ligand_name:
        raise ValueError(f"unsafe ligand name for template lookup: {ligand_name!r}")

    candidates: List[Path] = []
    suffix = "_ligand_start_conf"
    if ligand_name.endswith(suffix) and len(ligand_name) > len(suffix):
        complex_id = ligand_name[: -len(suffix)]
        candidates.append(template_root / complex_id / f"{complex_id}{suffix}.sdf")
    candidates.append(template_root / f"{ligand_name}.sdf")
    template_sdf = next((path for path in candidates if path.is_file()), None)
    prepared_pdbqt = prepared_root / f"{ligand_name}.pdbqt"
    if template_sdf is None or not prepared_pdbqt.is_file():
        if required:
            raise ValueError(
                f"{ligand_name}: required template pair is missing "
                f"(SDF candidates={[str(path) for path in candidates]}, "
                f"prepared PDBQT={prepared_pdbqt})"
            )
        return None
    return _load_mgl_pose_template(ligand_name, template_sdf, prepared_pdbqt)


def preflight_mgl_pose_templates(
    ligand_names: Iterable[str],
    cfg: dict,
) -> Dict[str, MGLPoseTemplate]:
    """Load every requested authoritative template pair before optimizer work.

    Thin benchmark drivers can call this during ``--verify-only`` so a missing,
    ambiguous, or atom-incompatible ligand fails before the long-running job is
    launched.
    """
    enabled, required, _, _ = _template_reconstruction_settings(cfg)
    if not enabled:
        return {}
    resolved: Dict[str, MGLPoseTemplate] = {}
    errors: List[str] = []
    for ligand_name in sorted({str(name) for name in ligand_names}):
        try:
            context = _resolve_mgl_pose_template(ligand_name, cfg)
        except Exception as e:
            errors.append(f"{ligand_name}: {e}")
            continue
        if context is None:
            if required:
                errors.append(f"{ligand_name}: required template context was not resolved")
            continue
        resolved[ligand_name] = context
    if errors:
        preview = "; ".join(errors[:10])
        if len(errors) > 10:
            preview += f"; ... {len(errors) - 10} more"
        raise ValueError(f"template reconstruction preflight failed: {preview}")
    return resolved


def _validate_reconstructed_sdf(
    model_pdbqt: Path,
    out_sdf: Path,
    pose_template: Optional[MGLPoseTemplate] = None,
) -> Tuple[bool, str, Optional[dict]]:
    """Verify chemistry and PDBQT→SDF atom-coordinate mapping."""
    try:
        source_text = model_pdbqt.read_text(errors="strict")
    except OSError as e:
        return False, f"unreadable source pose: {e}", None
    mol, message = _read_single_sdf_molecule(out_sdf)
    if mol is None:
        return False, message, None

    smiles = _remark_smiles(source_text)
    try:
        identity = _molecule_identity(mol)
    except Exception as e:
        return False, f"cannot determine reconstructed molecule identity: {e}", None

    if pose_template is not None:
        try:
            pose_signature = _pdbqt_atom_type_signature(source_text)
        except ValueError as e:
            return False, f"invalid template-mapped source pose: {e}", None
        if pose_signature != pose_template.prepared_atom_signature:
            return False, (
                "docked pose atom serial/type signature differs from the prepared "
                "MGLTools ligand"
            ), None
        coordinates = _pdbqt_atom_coordinates(source_text)
        elements = _pdbqt_atom_elements(source_text)
        signature_serials = {serial for serial, _ in pose_signature}
        if set(coordinates) != signature_serials or set(elements) != signature_serials:
            return False, "docked pose has invalid coordinates or AutoDock atom types", None
        if any(
            not all(math.isfinite(value) for value in xyz)
            for xyz in coordinates.values()
        ):
            return False, "docked pose has non-finite atom coordinates", None

        expected = pose_template.molecule
        try:
            expected_identity = _molecule_identity(expected)
            forward_mapping = mol.GetSubstructMatch(expected, useChirality=True)
            reverse_mapping = expected.GetSubstructMatch(mol, useChirality=True)
        except Exception as e:
            return False, f"cannot compare reconstructed template topology: {e}", None
        invariant_fields = (
            "canonical_isomeric_smiles", "formula", "heavy_atom_count",
        )
        if any(
            identity.get(field) != expected_identity.get(field)
            for field in invariant_fields
        ):
            return False, (
                "reconstructed molecule does not match authoritative SDF topology"
            ), None
        if (
            mol.GetNumAtoms() != expected.GetNumAtoms()
            or len(forward_mapping) != expected.GetNumAtoms()
            or len(reverse_mapping) != mol.GetNumAtoms()
        ):
            return False, (
                "reconstructed molecule has no complete chirality-preserving "
                "authoritative-template mapping"
            ), None

        conformer = mol.GetConformer()
        tolerance_sq = _TEMPLATE_MAPPING_TOLERANCE_A ** 2
        for atom_idx, serial in pose_template.atom_serial_map.items():
            output_atom = mol.GetAtomWithIdx(atom_idx)
            expected_atom = expected.GetAtomWithIdx(atom_idx)
            if output_atom.GetAtomicNum() != expected_atom.GetAtomicNum():
                return False, (
                    f"reconstructed SDF changed template atom {atom_idx + 1} identity"
                ), None
            if elements.get(serial) != expected_atom.GetSymbol():
                return False, (
                    f"element mapping mismatch for template atom {atom_idx + 1} "
                    f"and PDBQT atom {serial}"
                ), None
            position = conformer.GetAtomPosition(atom_idx)
            pose_position = coordinates[serial]
            delta_sq = sum(
                (actual - expected_coordinate) ** 2
                for actual, expected_coordinate in zip(
                    (position.x, position.y, position.z), pose_position,
                )
            )
            if delta_sq > tolerance_sq:
                return False, (
                    f"coordinate mapping mismatch for template atom {atom_idx + 1} "
                    f"and PDBQT atom {serial}"
                ), None
    elif smiles:
        expected = Chem.MolFromSmiles(smiles)
        if expected is None:
            return False, "invalid REMARK SMILES in source PDBQT", None
        expected_identity = _molecule_identity(expected)
        # Atom chiral tags (CW/CCW) encode parity relative to the molecule's
        # internal bond ordering, and SDF can represent an unspecified,
        # non-stereogenic double bond as STEREOANY instead of STEREONONE.
        # Consequently, atoms_by_index/bonds_by_index are useful diagnostics
        # but are not serialization-stable identity invariants.  Compare the
        # canonical chemistry first, then require a complete chirality-aware
        # graph match in both directions.  The explicit REMARK mapping checks
        # below still enforce index, element and coordinate preservation.
        invariant_fields = (
            "canonical_isomeric_smiles", "formula", "heavy_atom_count",
        )
        if any(
            identity.get(field) != expected_identity.get(field)
            for field in invariant_fields
        ):
            return False, "reconstructed molecule does not match REMARK SMILES topology", None
        try:
            expected_heavy = Chem.RemoveHs(expected)
            reconstructed_heavy = Chem.RemoveHs(mol)
            forward_mapping = reconstructed_heavy.GetSubstructMatch(
                expected_heavy, useChirality=True,
            )
            reverse_mapping = expected_heavy.GetSubstructMatch(
                reconstructed_heavy, useChirality=True,
            )
        except Exception as e:
            return False, f"cannot compare reconstructed molecule topology: {e}", None
        if (
            len(forward_mapping) != expected_heavy.GetNumAtoms()
            or len(reverse_mapping) != reconstructed_heavy.GetNumAtoms()
        ):
            return False, (
                "reconstructed molecule has no complete chirality-preserving "
                "REMARK SMILES topology mapping"
            ), None

        mapping = _remark_smiles_mapping(source_text)
        # Meeko's SMILES IDX maps every atom explicitly represented by REMARK
        # SMILES, not only heavy atoms. This includes isotope/directional
        # hydrogens such as the leading [H] in 5SAK_ZRY. Hydrogens absent from
        # the SMILES are represented separately by REMARK H PARENT.
        expected_atoms = expected.GetNumAtoms()
        if len(mapping) != expected_atoms:
            return False, (
                "incomplete REMARK SMILES IDX mapping "
                f"({len(mapping)} entries for {expected_atoms} explicit SMILES atoms)"
            ), None
        smiles_indices = [pair[0] for pair in mapping]
        pdbqt_indices = [pair[1] for pair in mapping]
        if len(set(smiles_indices)) != len(mapping) or len(set(pdbqt_indices)) != len(mapping):
            return False, "REMARK SMILES IDX mapping is not one-to-one", None
        if set(smiles_indices) != set(range(1, expected_atoms + 1)):
            return False, "REMARK SMILES IDX mapping does not cover every SMILES atom exactly once", None
        coordinates = _pdbqt_atom_coordinates(source_text)
        pdbqt_elements = _pdbqt_atom_elements(source_text)
        if mol.GetNumAtoms() < expected_atoms:
            return False, "reconstructed SDF lacks explicitly mapped SMILES atoms", None
        conformer = mol.GetConformer()
        for smiles_idx, pdbqt_idx in mapping:
            if pdbqt_idx not in coordinates or pdbqt_idx not in pdbqt_elements:
                return False, "REMARK SMILES IDX mapping references a missing atom", None
            expected_atom = expected.GetAtomWithIdx(smiles_idx - 1)
            output_atom = mol.GetAtomWithIdx(smiles_idx - 1)
            expected_element = expected_atom.GetSymbol()
            if output_atom.GetAtomicNum() != expected_atom.GetAtomicNum():
                return False, f"reconstructed SDF changed SMILES atom {smiles_idx} identity", None
            if pdbqt_elements[pdbqt_idx] != expected_element:
                return False, (
                    f"element mapping mismatch for SMILES atom {smiles_idx} "
                    f"({expected_element}) and PDBQT atom {pdbqt_idx} "
                    f"({pdbqt_elements[pdbqt_idx]})"
                ), None
            sdf_pos = conformer.GetAtomPosition(smiles_idx - 1)
            pdbqt_pos = coordinates[pdbqt_idx]
            delta_sq = sum(
                (a - b) ** 2 for a, b in zip((sdf_pos.x, sdf_pos.y, sdf_pos.z), pdbqt_pos)
            )
            if delta_sq > 0.05 ** 2:
                return False, (
                    f"coordinate mapping mismatch for SMILES atom {smiles_idx} "
                    f"and PDBQT atom {pdbqt_idx}"
                ), None
    else:
        # Plain legacy PDBQT has no authoritative bond-order template. Open
        # Babel is allowed only for non-macrocycles, and must preserve the
        # heavy-atom sequence and coordinate mapping exactly.
        source_atoms: List[Tuple[str, Tuple[float, float, float]]] = []
        for line in source_text.splitlines():
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_type = line[77:].strip().split()[0] if line[77:].strip() else ""
            if atom_type.upper().startswith("H"):
                continue
            element = re.sub(r"[^A-Za-z]", "", atom_type)
            if element in {"A", "C"}:
                element = "C"
            elif element in {"N", "NA", "NS"}:
                element = "N"
            elif element in {"O", "OA", "OS"}:
                element = "O"
            elif element in {"S", "SA"}:
                element = "S"
            elif element:
                element = element[0].upper() + element[1:].lower()
            try:
                coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                return False, "invalid PDBQT atom coordinates", None
            source_atoms.append((element, coords))
        heavy = Chem.RemoveHs(mol)
        if len(source_atoms) != heavy.GetNumAtoms():
            return False, "Open Babel conversion changed the heavy-atom count", None
        conformer = heavy.GetConformer()
        for idx, ((element, coords), atom) in enumerate(zip(source_atoms, heavy.GetAtoms())):
            if element and atom.GetSymbol() != element:
                return False, f"Open Babel changed atom identity at index {idx + 1}", None
            pos = conformer.GetAtomPosition(idx)
            if sum((a - b) ** 2 for a, b in zip((pos.x, pos.y, pos.z), coords)) > 0.05 ** 2:
                return False, f"Open Babel changed atom mapping at index {idx + 1}", None

    return True, "validated", identity


def _reconstruct_template_pose_sdf(
    model_pdbqt: Path,
    out_sdf: Path,
    pose_template: MGLPoseTemplate,
) -> Tuple[bool, str]:
    """Copy docked serial coordinates onto an authoritative RDKit SDF graph."""
    try:
        source_text = model_pdbqt.read_text(errors="strict")
        pose_signature = _pdbqt_atom_type_signature(source_text)
    except (OSError, ValueError) as e:
        return False, f"cannot read template-mapped pose: {e}"
    if pose_signature != pose_template.prepared_atom_signature:
        return False, (
            "docked pose atom serial/type signature differs from the prepared "
            "MGLTools ligand"
        )
    coordinates = _pdbqt_atom_coordinates(source_text)
    elements = _pdbqt_atom_elements(source_text)
    signature_serials = {serial for serial, _ in pose_signature}
    if set(coordinates) != signature_serials or set(elements) != signature_serials:
        return False, "docked pose has invalid coordinates or AutoDock atom types"
    if any(
        not all(math.isfinite(value) for value in xyz)
        for xyz in coordinates.values()
    ):
        return False, "docked pose has non-finite atom coordinates"

    reconstructed = Chem.Mol(pose_template.molecule)
    conformer = reconstructed.GetConformer()
    for atom_idx, serial in pose_template.atom_serial_map.items():
        expected_element = reconstructed.GetAtomWithIdx(atom_idx).GetSymbol()
        if elements.get(serial) != expected_element:
            return False, (
                f"element mapping mismatch for template atom {atom_idx + 1} "
                f"and PDBQT atom {serial}"
            )
        conformer.SetAtomPosition(atom_idx, coordinates[serial])

    fd, temporary_name = tempfile.mkstemp(
        prefix=f".{out_sdf.stem}.", suffix=".sdf", dir=str(out_sdf.parent),
    )
    os.close(fd)
    temporary = Path(temporary_name)
    temporary.unlink(missing_ok=True)
    try:
        writer = Chem.SDWriter(str(temporary))
        if writer is None:
            return False, "RDKit could not create the template-reconstructed SDF"
        try:
            writer.write(reconstructed)
        finally:
            writer.close()
        ok, message, _ = _validate_reconstructed_sdf(
            model_pdbqt, temporary, pose_template,
        )
        if not ok:
            return False, f"template reconstruction validation failed: {message}"
        os.replace(temporary, out_sdf)
        return True, _TEMPLATE_RECONSTRUCTION_CONVERTER
    except Exception as e:
        return False, f"template reconstruction failed: {e}"
    finally:
        temporary.unlink(missing_ok=True)


def _convert_pose_to_sdf(
    model_pdbqt: Path,
    out_sdf: Path,
    converters: Dict[str, str],
    pose_template: Optional[MGLPoseTemplate] = None,
) -> Tuple[bool, str]:
    """Reconstruct a real-molecule SDF from a single Vina pose PDBQT.

    smina/gnina cannot parse Meeko's macrocycle glue pseudo-atom types (CG0/G0),
    so every pose is first rebuilt into a clean SDF. Meeko ``mk_export`` uses the
    embedded ``REMARK SMILES`` to restore exact bond orders (closing macrocycle
    rings and dropping glue atoms). For MGLTools poses, an optional authoritative
    SDF/prepared-PDBQT context transfers docked coordinates without asking Open
    Babel to perceive any bonds. Open Babel remains the legacy fallback when no
    authoritative topology is configured. Returns ``(ok, converter_or_error)``.
    """
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    out_sdf.unlink(missing_ok=True)
    try:
        text = Path(model_pdbqt).read_text(errors="ignore")
    except OSError as e:
        return False, f"unreadable pose: {e}"

    has_remark_smiles = _remark_smiles(text) is not None
    has_macrocycle_glue = bool(re.search(
        r"^(?:ATOM|HETATM).*(?:\sCG\d*|\sG\d*)\s*$", text, flags=re.MULTILINE,
    ))
    mk = converters.get("mk_export")
    if has_remark_smiles:
        if not mk:
            return False, (
                "pose contains REMARK SMILES and requires Meeko mk_export; "
                "unsafe Open Babel fallback is disabled"
            )
        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{out_sdf.stem}.", suffix=".sdf", dir=str(out_sdf.parent),
        )
        os.close(fd)
        Path(tmp_name).unlink(missing_ok=True)
        try:
            # Meeko 0.5 uses -o; Meeko 0.7 changed its explicit SDF output to
            # -s. Detect the long-option semantics so -s is never mistaken for
            # the old filename-suffix option that caused the original bug.
            output_flag = _mk_export_output_flag(mk)
            r = subprocess.run([mk, str(model_pdbqt), output_flag, tmp_name],
                               capture_output=True, text=True, timeout=120)
            if r.returncode != 0:
                detail = (r.stderr or r.stdout or "").strip().splitlines()
                return False, (
                    f"mk_export rc={r.returncode}: "
                    f"{detail[-1] if detail else 'unknown error'}; "
                    "unsafe Open Babel fallback disabled"
                )
            ok, message, _ = _validate_reconstructed_sdf(model_pdbqt, Path(tmp_name))
            if not ok:
                return False, f"mk_export validation failed: {message}"
            os.replace(tmp_name, out_sdf)
            return True, "mk_export"
        except Exception as e:
            return False, f"mk_export failed: {e}; unsafe Open Babel fallback disabled"
        finally:
            Path(tmp_name).unlink(missing_ok=True)

    if pose_template is not None:
        return _reconstruct_template_pose_sdf(
            model_pdbqt, out_sdf, pose_template,
        )

    if has_macrocycle_glue:
        return False, (
            "pose contains macrocycle glue atoms but no usable REMARK SMILES; "
            "unsafe Open Babel fallback is disabled"
        )

    ob = converters.get("obabel")
    if ob:
        fd, tmp_name = tempfile.mkstemp(
            prefix=f".{out_sdf.stem}.", suffix=".sdf", dir=str(out_sdf.parent),
        )
        os.close(fd)
        Path(tmp_name).unlink(missing_ok=True)
        try:
            r = subprocess.run([ob, str(model_pdbqt), "-O", tmp_name],
                               capture_output=True, text=True, timeout=120)
            if r.returncode == 0:
                ok, message, _ = _validate_reconstructed_sdf(model_pdbqt, Path(tmp_name))
                if not ok:
                    return False, f"Open Babel validation failed: {message}"
                os.replace(tmp_name, out_sdf)
                return True, "obabel"
            detail = (r.stderr or r.stdout or "").strip().splitlines()
            return False, f"obabel rc={r.returncode}: {detail[-1] if detail else 'unknown error'}"
        except Exception as e:
            return False, f"obabel failed: {e}"
        finally:
            Path(tmp_name).unlink(missing_ok=True)

    return False, "no working pdbqt→sdf converter (need Meeko mk_export or Open Babel)"


def run_optimizer_tool(
    tool: str,
    exe: str,
    receptor: Path,
    pose_ligand: Path,
    out_sdf: Path,
    cfg: dict,
    expected_identity: Optional[dict] = None,
) -> Tuple[bool, Dict[str, float], str]:
    """Locally re-optimise a single (reconstructed-SDF) Vina pose with smina/gnina.

    Returns ``(success, scores, message)``. On failure the caller keeps the
    original Vina pose and its rank, so a pose is never lost.
    """
    _validate_optimizer_config({**cfg, "optimization": tool})
    settings = _optimizer_settings(tool, cfg)
    search = settings["search"]
    timeout_s = settings["timeout_s"]

    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    out_sdf.unlink(missing_ok=True)
    _provenance_path(out_sdf).unlink(missing_ok=True)

    if expected_identity is None:
        input_mol, input_message = _read_single_sdf_molecule(pose_ligand)
        if input_mol is None:
            return False, {}, f"invalid optimizer input: {input_message}"
        expected_identity = _molecule_identity(input_mol)

    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{out_sdf.stem}.", suffix=".sdf", dir=str(out_sdf.parent),
    )
    os.close(fd)
    tmp_out = Path(tmp_name)
    tmp_out.unlink(missing_ok=True)

    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(pose_ligand),
        "--autobox_ligand", str(pose_ligand),
        "--autobox_add", str(settings["autobox_add"]),
        "--out", str(tmp_out),
        "--cpu", str(settings["cpu"]),
        "--seed", str(settings["seed"]),
        "--num_modes", "1",
        "--scoring", settings["empirical_scoring"],
    ]
    # --minimize (pure energy min) and --local_only (local MC) both keep the
    # ligand at its input pose and never run a global search, so the optimised
    # pose stays in the same pocket — preserving the pose identity we re-rank.
    cmd.append("--minimize" if search == "minimize" else "--local_only")
    if _is_gnina_optimizer(tool):
        cmd.extend(["--cnn_scoring", settings["cnn_scoring"]])
        if settings["cnn_model_file"]:
            cmd.extend(["--cnn_model", settings["cnn_model_file"]])
        elif settings["cnn_model"]:
            cmd.extend(["--cnn", settings["cnn_model"]])
        if not settings["use_gpu"]:
            cmd.append("--no_gpu")

    proc_env = _gnina_subprocess_env(cfg) if _is_gnina_optimizer(tool) else None

    gpu_lock = (
        _interprocess_gpu_lock(cfg["gpu_lock_file"])
        if _is_gnina_optimizer(tool) and settings["use_gpu"]
        else contextlib.nullcontext()
    )
    try:
        # The blocking lock wait happens before subprocess.run starts, so it does
        # not consume the optimizer's per-pose timeout. flock is released by the
        # context manager on normal exit, timeout, exception, or process death.
        with gpu_lock:
            proc = subprocess.run(
                cmd, capture_output=True, text=True, timeout=timeout_s, env=proc_env,
            )
    except subprocess.TimeoutExpired:
        tmp_out.unlink(missing_ok=True)
        return False, {}, f"{tool} timed out after {timeout_s}s"
    except FileNotFoundError:
        tmp_out.unlink(missing_ok=True)
        return False, {}, f"{tool} executable not found: {exe}"
    except Exception as e:  # pragma: no cover - defensive
        tmp_out.unlink(missing_ok=True)
        return False, {}, f"{tool} invocation failed: {e}"

    if proc.returncode != 0:
        tmp_out.unlink(missing_ok=True)
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        return False, {}, f"{tool} rc={proc.returncode}: {err[-1] if err else 'unknown error'}"

    if not tmp_out.exists() or tmp_out.stat().st_size == 0:
        tmp_out.unlink(missing_ok=True)
        return False, {}, f"{tool} produced no output pose"

    ok, scores, message = _validate_optimizer_sdf(
        tmp_out, expected_identity, _required_score_columns(tool, cfg),
    )
    if not ok:
        tmp_out.unlink(missing_ok=True)
        return False, scores, f"invalid {tool} output: {message}"
    os.replace(tmp_out, out_sdf)
    return True, scores, f"{tool} {search} ok"


def _provenance_path(out_sdf: Path) -> Path:
    return Path(f"{out_sdf}.provenance.json")


def _build_optimizer_provenance(
    *,
    tool: str,
    source_pose: Path,
    source_pose_container: Path,
    receptor: Path,
    cfg: dict,
    converter: str,
    converter_exe: str,
    converter_version: str,
    optimizer_exe: str,
    optimizer_version: str,
    input_identity: dict,
    converter_inputs: Optional[dict] = None,
) -> dict:
    converter_material = {
        "name": converter,
        **_executable_identity(converter_exe, converter_version),
    }
    if converter_inputs:
        converter_material["inputs"] = dict(converter_inputs)
    material = {
        "schema_version": OPTIMIZER_CACHE_SCHEMA,
        "pipeline": "run_autodock.optimize_autodock_results",
        "tool": tool,
        "source_pose_sha256": _sha256_file(source_pose),
        "source_pose_container_sha256": _sha256_file(source_pose_container),
        "receptor_sha256": _sha256_file(receptor),
        "source_scoring": str(cfg.get("scoring_function", "vina")).strip().lower(),
        "rank_metric": str(cfg.get("optimize_rank_by", "minimized_affinity")).strip().lower(),
        "settings": _optimizer_settings(tool, cfg),
        "optimizer": _executable_identity(optimizer_exe, optimizer_version),
        "converter": converter_material,
        "input_identity": input_identity,
    }
    encoded = json.dumps(material, sort_keys=True, separators=(",", ":")).encode()
    return {
        **material,
        "fingerprint": hashlib.sha256(encoded).hexdigest(),
        "source_pose": str(source_pose_container.resolve()),
        "receptor": str(receptor.resolve()),
        "created_at": datetime.now().isoformat(),
    }


def _write_optimizer_provenance(
    out_sdf: Path,
    provenance: dict,
    optimizer_elapsed_time_s: float,
) -> None:
    payload = dict(provenance)
    payload["output_sha256"] = _sha256_file(out_sdf)
    payload["optimizer_elapsed_time_s"] = round(float(optimizer_elapsed_time_s), 6)
    _atomic_write_json(_provenance_path(out_sdf), payload)


def _validate_cached_optimizer_output(
    out_sdf: Path,
    expected_provenance: dict,
    expected_identity: dict,
    required_scores: Set[str],
) -> Tuple[bool, Dict[str, float], str, float]:
    sidecar = _provenance_path(out_sdf)
    if not out_sdf.exists() or not sidecar.exists():
        return False, {}, "cache missing output or provenance sidecar", 0.0
    try:
        saved = json.loads(sidecar.read_text())
    except Exception as e:
        return False, {}, f"invalid provenance sidecar: {e}", 0.0
    if saved.get("schema_version") != OPTIMIZER_CACHE_SCHEMA:
        return False, {}, "legacy/incompatible provenance schema", 0.0
    if saved.get("fingerprint") != expected_provenance.get("fingerprint"):
        return False, {}, "provenance fingerprint mismatch", 0.0
    for path_field in ("source_pose", "receptor"):
        if saved.get(path_field) != expected_provenance.get(path_field):
            return False, {}, f"provenance {path_field} path mismatch", 0.0
    try:
        if saved.get("output_sha256") != _sha256_file(out_sdf):
            return False, {}, "optimized SDF checksum mismatch", 0.0
    except OSError as e:
        return False, {}, f"cannot hash optimized SDF: {e}", 0.0
    ok, scores, message = _validate_optimizer_sdf(
        out_sdf, expected_identity, required_scores,
    )
    if not ok:
        return False, scores, message, 0.0
    try:
        optimizer_elapsed = float(saved["optimizer_elapsed_time_s"])
    except (KeyError, TypeError, ValueError):
        return False, {}, "provenance lacks durable optimizer elapsed time", 0.0
    if not math.isfinite(optimizer_elapsed) or optimizer_elapsed < 0:
        return False, {}, "invalid optimizer elapsed time in provenance", 0.0
    return True, scores, "reused validated cache", optimizer_elapsed


def _discard_optimizer_cache(out_sdf: Path) -> None:
    out_sdf.unlink(missing_ok=True)
    _provenance_path(out_sdf).unlink(missing_ok=True)


def _prune_obsolete_optimizer_outputs(
    tool_dir: Path,
    pose_stem: str,
    tool: str,
    desired_outputs: Set[Path],
) -> int:
    """Remove only outputs belonging to this pose set but no longer selected."""
    removed = 0
    desired = {path.resolve() for path in desired_outputs}
    if not tool_dir.exists():
        return removed
    for sdf_path in tool_dir.glob(f"{pose_stem}_rank*_{tool}.sdf"):
        if sdf_path.resolve() not in desired:
            _discard_optimizer_cache(sdf_path)
            removed += 1
    for sidecar in tool_dir.glob(f"{pose_stem}_rank*_{tool}.sdf.provenance.json"):
        sdf_path = Path(str(sidecar)[: -len(".provenance.json")])
        if not sdf_path.exists() or sdf_path.resolve() not in desired:
            sidecar.unlink(missing_ok=True)
    return removed


def _optimizer_preflight(
    cfg: dict,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str], Dict[str, str], List[str]]:
    """Resolve and load-test all explicitly requested optimizer dependencies."""
    tools = resolve_optimizers(cfg)
    _validate_optimizer_config(cfg)
    exes: Dict[str, str] = {}
    optimizer_versions: Dict[str, str] = {}
    errors: List[str] = []
    for tool in tools:
        exe = _resolve_optimizer_exe(tool, cfg)
        if not exe:
            errors.append(f"{tool} executable not found")
            continue
        ok, version = _probe_optimizer(tool, exe, cfg)
        if not ok:
            errors.append(f"{tool} failed to load: {version}")
            continue
        exes[tool] = exe
        optimizer_versions[tool] = version
        settings = _optimizer_settings(tool, cfg)
        if _is_gnina_optimizer(tool):
            model_file = settings.get("cnn_model_file", "")
            model = settings.get("cnn_model", "")
            if model_file and not Path(model_file).is_file():
                errors.append(f"gnina CNN model file not found: {model_file}")
            elif settings.get("cnn_scoring") != "none" and not model_file:
                try:
                    help_proc = subprocess.run(
                        [exe, "--help"], capture_output=True, text=True,
                        timeout=30, env=_gnina_subprocess_env(cfg),
                    )
                    help_text = (help_proc.stdout or "") + (help_proc.stderr or "")
                    model_prefix = model.removesuffix("_ensemble")
                    if help_proc.returncode != 0 or not model_prefix or model_prefix not in help_text:
                        errors.append(f"gnina built-in CNN model is not advertised by this binary: {model}")
                except Exception as e:
                    errors.append(f"cannot validate gnina CNN model {model!r}: {e}")
            if settings.get("use_gpu"):
                lock_path = Path(cfg["gpu_lock_file"]).expanduser().resolve()
                try:
                    lock_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(lock_path, "a+"):
                        pass
                except OSError as e:
                    errors.append(f"gnina GPU lock is not usable ({lock_path}): {e}")
                if cfg.get("gnina_gpu_preflight", True):
                    nvidia_smi = shutil.which("nvidia-smi")
                    if not nvidia_smi:
                        errors.append("gnina GPU requested but nvidia-smi is unavailable")
                    else:
                        try:
                            gpu_probe = subprocess.run(
                                [nvidia_smi, "-L"], capture_output=True, text=True, timeout=30,
                            )
                            if gpu_probe.returncode != 0 or "GPU " not in gpu_probe.stdout:
                                errors.append("gnina GPU requested but no CUDA GPU is visible")
                        except Exception as e:
                            errors.append(f"gnina GPU preflight failed: {e}")

    converters = _resolve_pose_converters(cfg)
    converter_versions: Dict[str, str] = {}
    for name, exe in converters.items():
        if exe:
            probe_args = [exe, "--help"] if name == "mk_export" else [exe, "-V"]
            try:
                probe = subprocess.run(
                    probe_args, capture_output=True, text=True, timeout=30,
                )
            except Exception as e:
                errors.append(f"{name} failed to load: {e}")
                converters[name] = ""
                continue
            if probe.returncode != 0:
                detail = (probe.stderr or probe.stdout or "").strip().splitlines()
                errors.append(
                    f"{name} failed to load: "
                    f"{detail[-1] if detail else f'rc={probe.returncode}'}"
                )
                converters[name] = ""
                continue
            converter_versions[name] = _converter_version(name, exe)

    template_enabled = False
    template_required = bool(
        cfg.get("optimize_require_template_reconstruction", False)
    )
    try:
        template_enabled, template_required, _, _ = (
            _template_reconstruction_settings(cfg)
        )
    except ValueError as e:
        errors.append(f"template reconstruction preflight failed: {e}")
    if template_enabled:
        converters[_TEMPLATE_RECONSTRUCTION_CONVERTER] = str(
            Path(__file__).resolve()
        )
        converter_versions[_TEMPLATE_RECONSTRUCTION_CONVERTER] = (
            f"rdkit {Chem.rdBase.rdkitVersion}; {_TEMPLATE_RECONSTRUCTION_SCHEMA}"
        )

    require_mk_export = bool(cfg.get("optimize_require_mk_export", True))
    if (
        require_mk_export and not template_required
        and not converters.get("mk_export")
    ):
        errors.append(
            "Meeko mk_export not found (required for REMARK SMILES/macrocycle-safe reconstruction)"
        )
    elif not any(converters.values()):
        errors.append("no PDBQT→SDF converter found")
    return exes, optimizer_versions, converters, converter_versions, errors


def preflight_optimizers(cfg: dict, *, refresh: bool = False) -> dict:
    """Public fail-fast preflight for scripts and notebook thin drivers.

    Returns a JSON-serializable dictionary with ``status``, requested tools,
    resolved optimizer/converter paths and versions, and normalized settings.
    On any requested dependency/model failure it raises
    :class:`OptimizationError`; ``exc.summary.to_dict()`` is the structured
    failure result. Successful results are cached for the process so a notebook
    optimizing hundreds of per-complex folders does not repeatedly probe tools.
    """
    requested = resolve_optimizers(cfg)
    if not requested:
        return {
            "status": "disabled", "requested_tools": [],
            "optimizers": {}, "converters": {}, "settings": {},
        }
    _validate_optimizer_config(cfg)
    relevant_keys = [
        "optimization", "smina_executable", "gnina_executable",
        "mk_export_executable", "obabel_executable", "gnina_lib_dirs",
        "gnina_python", "diffdock_python", "optimize_require_mk_export",
        "optimize_require_template_reconstruction",
        "optimize_ligand_template_root", "optimize_prepared_ligand_pdbqt_dir",
        "optimize_search", "optimize_scoring", "optimize_autobox_add",
        "optimize_cpu", "optimize_seed", "gnina_cnn_scoring",
        "gnina_cnn_model", "gnina_cnn_model_file", "gnina_use_gpu",
        "gnina_gpu_preflight", "gpu_lock_file",
    ]
    cache_material = {
        key: cfg.get(key) for key in relevant_keys
    }
    cache_material.update({"PATH": os.environ.get("PATH", ""), "python": sys.executable})
    resolved_dependencies = {
        **{tool: _resolve_optimizer_exe(tool, cfg) for tool in requested},
        **{
            f"converter:{name}": path
            for name, path in _resolve_pose_converters(cfg).items()
        },
    }
    cache_material["resolved_dependencies"] = {}
    for name, path in resolved_dependencies.items():
        identity: Dict[str, Any] = {"path": path}
        if path:
            with contextlib.suppress(OSError):
                stat = Path(path).resolve().stat()
                identity.update({"size": stat.st_size, "mtime_ns": stat.st_mtime_ns})
        cache_material["resolved_dependencies"][name] = identity
    cache_key = hashlib.sha256(
        json.dumps(cache_material, sort_keys=True, default=str).encode()
    ).hexdigest()
    if not refresh and cache_key in _OPTIMIZER_PREFLIGHT_CACHE:
        return json.loads(json.dumps(_OPTIMIZER_PREFLIGHT_CACHE[cache_key]))

    exes, versions, converters, converter_versions, errors = _optimizer_preflight(cfg)
    if errors:
        summary = OptimizationSummary(
            requested_tools=list(requested),
            available_tools=[tool for tool in requested if tool in exes],
            status="preflight_failed",
            errors=list(errors),
        )
        raise OptimizationError("optimizer preflight failed: " + "; ".join(errors), summary)
    payload = {
        "status": "ready",
        "requested_tools": list(requested),
        "optimizers": {
            tool: {"path": exes[tool], "version": versions[tool]}
            for tool in requested
        },
        "converters": {
            name: {"path": exe, "version": converter_versions.get(name, "unknown")}
            for name, exe in converters.items() if exe
        },
        "settings": {tool: _optimizer_settings(tool, cfg) for tool in requested},
    }
    _OPTIMIZER_PREFLIGHT_CACHE[cache_key] = payload
    return json.loads(json.dumps(payload))


def _write_optimization_log(
    log_path: Path,
    rows: List[dict],
    replace_scopes: Optional[Set[Tuple[str, str]]] = None,
) -> None:
    """Atomically replace current combo/tool scopes under an inter-process lock.

    Replacing a complete scope, instead of only de-duplicating ranks, removes
    obsolete rows when the pose set or ``optimize_top_n`` shrinks.
    """
    if not rows and not replace_scopes:
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = Path(f"{log_path}.lock")
    with _OPTIMIZATION_LOG_LOCK, open(lock_path, "a+") as lock_handle:
        if fcntl is not None:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            if log_path.exists():
                try:
                    # Pandas' native CSV parser can abort the entire process on
                    # otherwise valid optimizer logs (observed as SIGSEGV in
                    # pandas_parser.so).  The stdlib reader is fast enough for
                    # these small per-complex logs and, importantly, turns bad
                    # input into a catchable Python exception.
                    with log_path.open("r", encoding="utf-8-sig", newline="") as handle:
                        reader = csv.DictReader(handle)
                        if not reader.fieldnames:
                            raise ValueError("optimization log has no header")
                        records = list(reader)
                        if any(None in row for row in records):
                            raise ValueError("optimization log has excess CSV fields")
                    df_old = pd.DataFrame.from_records(
                        records, columns=reader.fieldnames,
                    )
                except Exception:
                    # Never merge against a malformed/truncated concurrent log.
                    df_old = pd.DataFrame(columns=OPTIMIZATION_LOG_COLUMNS)
            else:
                df_old = pd.DataFrame(columns=OPTIMIZATION_LOG_COLUMNS)

            if replace_scopes and not df_old.empty and {"combo_name", "tool"}.issubset(df_old):
                scope_mask = pd.Series(False, index=df_old.index)
                for combo_name, tool in replace_scopes:
                    scope_mask |= (
                        (df_old["combo_name"].astype(str) == str(combo_name))
                        & (df_old["tool"].astype(str) == str(tool))
                    )
                df_old = df_old.loc[~scope_mask]

            df_new = pd.DataFrame(rows)
            frames = [frame for frame in (df_old, df_new) if not frame.empty]
            df_all = (
                pd.concat(frames, ignore_index=True, sort=False)
                if frames else pd.DataFrame(columns=OPTIMIZATION_LOG_COLUMNS)
            )
            key_cols = ["combo_name", "tool", "autodock_rank"]
            if not df_all.empty and set(key_cols).issubset(df_all.columns):
                # Existing CSV values are strings while freshly generated ranks
                # are integers. Compare normalized keys without changing the
                # values written to the merged log.
                normalized_keys = df_all[key_cols].astype(str)
                df_all = df_all.loc[~normalized_keys.duplicated(keep="last")]
            extra_columns = [
                column for column in df_all.columns
                if column not in OPTIMIZATION_LOG_COLUMNS
            ]
            df_all = df_all.reindex(columns=OPTIMIZATION_LOG_COLUMNS + extra_columns)

            fd, tmp_name = tempfile.mkstemp(
                prefix=f".{log_path.name}.", suffix=".tmp", dir=str(log_path.parent),
            )
            os.close(fd)
            try:
                df_all.to_csv(tmp_name, index=False)
                with open(tmp_name, "rb") as tmp_handle:
                    os.fsync(tmp_handle.fileno())
                os.replace(tmp_name, log_path)
            finally:
                Path(tmp_name).unlink(missing_ok=True)
        finally:
            if fcntl is not None:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def discover_existing_autodock_results(
    output_tree: Path,
    receptors: Iterable[Path] | Dict[str, Path],
) -> List[DockingResult]:
    """Build optimizer-ready results from existing ``*_vina_out.pdbqt`` files.

    ``receptors`` may be a mapping from the output filename's protein prefix to
    the exact prepared receptor, or an iterable whose stems are matched against
    that prefix. Ambiguous/unmatched receptors and malformed pose files raise;
    the helper never guesses a receptor or silently accepts truncated output.
    """
    if isinstance(receptors, dict):
        receptor_lookup = {str(key): Path(value) for key, value in receptors.items()}
    else:
        receptor_lookup = {}
        for receptor in receptors:
            path = Path(receptor)
            receptor_lookup[path.stem] = path
            receptor_lookup[get_file_stem(path)] = path

    discovered: List[DockingResult] = []
    errors: List[str] = []
    for pose_path in sorted(Path(output_tree).rglob("*_vina_out.pdbqt")):
        if any(part.startswith("optimized_") for part in pose_path.parts):
            continue
        protein_prefix, separator, ligand_part = pose_path.stem.partition("__")
        if not separator:
            errors.append(f"cannot parse protein/ligand names from {pose_path}")
            continue
        exact = receptor_lookup.get(protein_prefix)
        candidates = {
            path.resolve(): path for key, path in receptor_lookup.items()
            if protein_prefix == key or protein_prefix.startswith(f"{key}_")
        }
        receptor = exact or (next(iter(candidates.values())) if len(candidates) == 1 else None)
        if receptor is None or not receptor.exists():
            errors.append(f"no unique prepared receptor for {pose_path.name}")
            continue
        try:
            with tempfile.TemporaryDirectory() as tmp:
                models = _split_pdbqt_models(pose_path, Path(tmp))
        except Exception as e:
            errors.append(f"invalid poses in {pose_path}: {e}")
            continue
        affinities = [affinity for _, _, affinity in models if affinity is not None]
        scoring_match = re.search(r"_(vina|vinardo|ad4)_vina_out$", ligand_part)
        if not scoring_match:
            errors.append(f"cannot derive source scoring from {pose_path.name}")
            continue
        source_scoring = scoring_match.group(1)
        ligand_name = ligand_part[: scoring_match.start()]
        discovered.append(DockingResult(
            # The output prefix is the exact receptor identity used by the live
            # docking path.  Do not pass the prepared receptor through
            # ``get_file_stem`` here: that helper strips the semantic token
            # ``_protein`` and would give optimizer-only/resumed runs a
            # different combo_name from freshly docked runs.
            protein_name=protein_prefix,
            ligand_name=ligand_name,
            protein_path=receptor,
            ligand_path=pose_path,
            output_dir=pose_path.parent,
            status="skipped",
            scoring_function=source_scoring,
            num_poses=len(models),
            pose_files=[pose_path],
            best_affinity=min(affinities) if affinities else None,
        ))
    if errors:
        preview = "; ".join(errors[:5])
        if len(errors) > 5:
            preview += f"; … {len(errors) - 5} more"
        raise ValueError(f"existing AutoDock output discovery failed: {preview}")
    return discovered


def optimize_autodock_results(
    results: List[DockingResult],
    cfg: dict,
    output_dir: Path,
    overwrite: bool = False,
) -> OptimizationSummary:
    """Run the configured optimiser(s) over every docked pose and re-rank them.

    Operates on both freshly docked and previously docked (skipped) complexes, so
    enabling ``optimization`` and re-running re-optimises existing poses without
    re-docking. A cache is reused only after provenance, checksum, molecular
    identity, topology, atom mapping and score-tag validation. Explicit setup
    failures and all-pose failures raise :class:`OptimizationError`; its
    ``summary`` attribute contains the partial structured outcome.
    """
    requested_tools = resolve_optimizers(cfg)
    summary = OptimizationSummary(requested_tools=list(requested_tools))
    if not requested_tools:
        return summary

    output_dir = Path(output_dir)
    summary.log_path = str(output_dir / "optimization_log.csv")
    try:
        preflight = preflight_optimizers(cfg)
    except OptimizationError as e:
        e.summary.log_path = summary.log_path
        raise
    exes = {
        tool: details["path"] for tool, details in preflight["optimizers"].items()
    }
    optimizer_versions = {
        tool: details["version"] for tool, details in preflight["optimizers"].items()
    }
    converters = {
        name: details["path"] for name, details in preflight["converters"].items()
    }
    # Keep stable converter keys for selection logic.
    converters.setdefault("mk_export", "")
    converters.setdefault("obabel", "")
    converter_versions = {
        name: details["version"] for name, details in preflight["converters"].items()
    }
    summary.available_tools = [tool for tool in requested_tools if tool in exes]
    tools = list(requested_tools)

    targets = [
        r for r in results
        if r.status in ("success", "skipped") and r.pose_files
        and Path(r.pose_files[0]).exists()
    ]
    summary.target_complexes = len(targets)

    pose_templates: Dict[int, Optional[MGLPoseTemplate]] = {
        id(result): None for result in targets
    }
    try:
        template_enabled, _, _, _ = _template_reconstruction_settings(cfg)
    except ValueError as e:
        summary.status = "preflight_failed"
        summary.errors.append(f"template reconstruction preflight failed: {e}")
        raise OptimizationError(summary.errors[-1], summary) from e
    if template_enabled:
        try:
            templates_by_ligand = preflight_mgl_pose_templates(
                (result.ligand_name for result in targets), cfg,
            )
        except ValueError as e:
            summary.status = "preflight_failed"
            summary.errors.append(str(e))
            raise OptimizationError(str(e), summary) from e
        for result in targets:
            pose_templates[id(result)] = templates_by_ligand.get(result.ligand_name)

    print("\n" + "=" * 80)
    print(f"Pose optimization + re-ranking: {', '.join(tools)}")
    print("=" * 80)
    print(f"Optimizing {len(targets)} complex(es) × {len(tools)} tool(s) "
          f"(search={cfg.get('optimize_search', 'minimize')})")

    log_path = Path(summary.log_path)
    top_n = int(cfg.get("optimize_top_n", 0) or 0)
    rank_by = str(cfg.get("optimize_rank_by", "minimized_affinity")).strip().lower()

    rows: List[dict] = []
    time_by_tool: Dict[str, float] = defaultdict(float)
    runs_by_tool: Dict[str, int] = defaultdict(int)
    replace_scopes: Set[Tuple[str, str]] = set()

    for idx, r in enumerate(targets, 1):
        receptor = Path(r.protein_path)
        pose_pdbqt = Path(r.pose_files[0])
        pose_template = pose_templates[id(r)]
        combo_name = f"{r.ligand_name}__{r.protein_name}"
        pose_cfg = {
            **cfg,
            "scoring_function": r.scoring_function or cfg.get("scoring_function", "vina"),
        }

        for tool in tools:
            replace_scopes.add((combo_name, tool))
            tool_dir = pose_pdbqt.parent / f"optimized_{tool}"
            pose_records: List[dict] = []

            with tempfile.TemporaryDirectory() as tmp:
                try:
                    split = _split_pdbqt_models(pose_pdbqt, Path(tmp))
                except Exception as e:
                    summary.failed += 1
                    summary.errors.append(f"{combo_name}/{tool}: invalid pose file: {e}")
                    continue
                if top_n > 0:
                    split = [p for p in split if p[0] <= top_n]
                desired_outputs = {
                    tool_dir / f"{pose_pdbqt.stem}_rank{rank}_{tool}.sdf"
                    for rank, _, _ in split
                }
                summary.pruned_outputs += _prune_obsolete_optimizer_outputs(
                    tool_dir, pose_pdbqt.stem, tool, desired_outputs,
                )

                for rank, pose_file, vina_aff in split:
                    summary.selected_poses += 1
                    out_sdf = tool_dir / f"{pose_pdbqt.stem}_rank{rank}_{tool}.sdf"
                    if overwrite:
                        _discard_optimizer_cache(out_sdf)
                    t0 = time.time()
                    pose_sdf = Path(tmp) / f"{pose_file.stem}_in.sdf"
                    if pose_template is None:
                        conv_ok, conv_msg = _convert_pose_to_sdf(
                            pose_file, pose_sdf, converters,
                        )
                    else:
                        conv_ok, conv_msg = _convert_pose_to_sdf(
                            pose_file, pose_sdf, converters, pose_template,
                        )
                    converter = conv_msg if conv_ok else ""
                    active_pose_template = (
                        pose_template
                        if converter == _TEMPLATE_RECONSTRUCTION_CONVERTER
                        else None
                    )
                    converter_version = converter_versions.get(converter, "")
                    provenance: Optional[dict] = None
                    expected_identity: Optional[dict] = None
                    executed = False
                    execution_elapsed = 0.0
                    if not conv_ok:
                        _discard_optimizer_cache(out_sdf)
                        ok, scores = False, {}
                        msg = f"pdbqt→sdf failed: {conv_msg}"
                    else:
                        if active_pose_template is None:
                            valid_input, input_message, expected_identity = (
                                _validate_reconstructed_sdf(pose_file, pose_sdf)
                            )
                        else:
                            valid_input, input_message, expected_identity = (
                                _validate_reconstructed_sdf(
                                    pose_file, pose_sdf, active_pose_template,
                                )
                            )
                        if not valid_input or expected_identity is None:
                            _discard_optimizer_cache(out_sdf)
                            ok, scores = False, {}
                            msg = f"pdbqt→sdf validation failed: {input_message}"
                        else:
                            try:
                                provenance = _build_optimizer_provenance(
                                    tool=tool,
                                    source_pose=pose_file,
                                    source_pose_container=pose_pdbqt,
                                    receptor=receptor,
                                    cfg=pose_cfg,
                                    converter=converter,
                                    converter_exe=converters[converter],
                                    converter_version=converter_version,
                                    optimizer_exe=exes[tool],
                                    optimizer_version=optimizer_versions[tool],
                                    input_identity=expected_identity,
                                    converter_inputs=(
                                        active_pose_template.converter_inputs()
                                        if active_pose_template is not None else None
                                    ),
                                )
                            except Exception as e:
                                _discard_optimizer_cache(out_sdf)
                                ok, scores = False, {}
                                msg = f"cannot build provenance: {e}"
                            else:
                                cache_ok = False
                                if not overwrite:
                                    cache_ok, scores, msg, execution_elapsed = (
                                        _validate_cached_optimizer_output(
                                        out_sdf, provenance, expected_identity,
                                        _required_score_columns(tool, pose_cfg),
                                    )
                                    )
                                if cache_ok:
                                    ok = True
                                    summary.reused += 1
                                else:
                                    _discard_optimizer_cache(out_sdf)
                                    summary.attempted += 1
                                    runs_by_tool[tool] += 1
                                    executed = True
                                    execution_started = time.time()
                                    ok, scores, msg = run_optimizer_tool(
                                        tool, exes[tool], receptor, pose_sdf,
                                        out_sdf, pose_cfg, expected_identity=expected_identity,
                                    )
                                    execution_elapsed = time.time() - execution_started
                                    if ok:
                                        try:
                                            _write_optimizer_provenance(
                                                out_sdf, provenance, execution_elapsed,
                                            )
                                        except Exception as e:
                                            _discard_optimizer_cache(out_sdf)
                                            ok, scores = False, {}
                                            msg = f"cannot commit provenance: {e}"
                                        else:
                                            summary.optimized += 1
                    elapsed = time.time() - t0
                    if executed:
                        time_by_tool[tool] += execution_elapsed
                    if not ok:
                        summary.failed += 1
                        summary.errors.append(f"{combo_name}/{tool}/rank{rank}: {msg}")
                    pose_records.append({
                        "rank": rank, "vina_aff": vina_aff, "ok": ok,
                        "scores": scores, "msg": msg, "elapsed": elapsed,
                        "out_sdf": out_sdf, "converter": converter,
                        "converter_version": converter_version,
                        "provenance": provenance,
                        "execution_elapsed": execution_elapsed,
                    })

            # Re-rank on the chosen optimiser score (ties broken by the original
            # Vina rank for a stable order). smina has no CNN metrics, so a CNN
            # request falls back to its minimised affinity. Poses whose
            # optimisation failed carry no score and get no optimized_rank — they
            # keep only their autodock_rank.
            metric_col, higher_better = _RERANK_METRICS[rank_by]
            if tool == "smina" and metric_col.startswith("cnn"):
                metric_col, higher_better = "minimized_affinity", False
            ranked = [
                p for p in pose_records
                if p["ok"] and p["scores"].get(metric_col) is not None
            ]
            ranked.sort(key=lambda p: (
                -p["scores"][metric_col] if higher_better else p["scores"][metric_col],
                p["rank"],
            ))
            opt_rank_by_orig = {p["rank"]: i + 1 for i, p in enumerate(ranked)}

            for p in pose_records:
                sc = p["scores"]
                rows.append({
                    "timestamp": datetime.now().isoformat(),
                    "combo_name": combo_name,
                    "protein_name": r.protein_name,
                    "ligand_name": r.ligand_name,
                    "tool": tool,
                    "autodock_rank": p["rank"],
                    "optimized_rank": opt_rank_by_orig.get(p["rank"]),
                    "rank_metric": metric_col if p["ok"] else None,
                    "vina_affinity": p["vina_aff"],
                    "minimized_affinity": sc.get("minimized_affinity"),
                    "cnn_score": sc.get("cnn_score"),
                    "cnn_affinity": sc.get("cnn_affinity"),
                    "pose_file": pose_pdbqt.name,
                    "optimized_file": str(p["out_sdf"]) if p["ok"] else "",
                    "status": "success" if p["ok"] else "failed",
                    "message": p["msg"],
                    # Durable execution effort is restored from provenance on
                    # cache reuse; reruns therefore cannot erase runtime cost.
                    "elapsed_time_s": round(p["execution_elapsed"], 2),
                    "converter": p["converter"],
                    "converter_version": p["converter_version"],
                    "optimizer_version": optimizer_versions[tool],
                    "source_scoring": str(pose_cfg.get("scoring_function", "vina")).lower(),
                    "optimizer_scoring": _optimizer_settings(tool, pose_cfg)["empirical_scoring"],
                    "cnn_scoring": _optimizer_settings(tool, pose_cfg).get("cnn_scoring", "none"),
                    "cnn_model": _optimizer_settings(tool, pose_cfg).get("cnn_model", ""),
                    "source_pose_sha256": (
                        p["provenance"].get("source_pose_sha256", "")
                        if p["provenance"] else ""
                    ),
                    "receptor_sha256": (
                        p["provenance"].get("receptor_sha256", "")
                        if p["provenance"] else ""
                    ),
                    "provenance_fingerprint": (
                        p["provenance"].get("fingerprint", "")
                        if p["provenance"] else ""
                    ),
                    "provenance_file": (
                        str(_provenance_path(p["out_sdf"])) if p["ok"] else ""
                    ),
                    "processing_elapsed_time_s": round(p["elapsed"], 2),
                })

        print(f"  [{idx}/{len(targets)}] {combo_name}: "
              f"{r.num_poses or len(r.pose_files)} pose(s) × {len(tools)} tool(s)")

    _write_optimization_log(log_path, rows, replace_scopes=replace_scopes)
    summary.time_by_tool_s = {
        tool: round(time_by_tool.get(tool, 0.0), 6) for tool in tools
    }
    print(
        f"\nOptimization complete: {summary.optimized} optimized, "
        f"{summary.reused} reused, {summary.failed} failed"
    )
    if any(runs_by_tool.values()):
        print("Time by tool (this run, excludes reused):")
        for tool in tools:
            secs = time_by_tool.get(tool, 0.0)
            n = runs_by_tool.get(tool, 0)
            avg = secs / n if n else 0.0
            print(f"  {tool}: {secs:.1f}s total ({secs / 60:.1f} min) "
                  f"over {n} pose(s) ({avg:.1f}s/pose)")
    print(f"Optimization log: {log_path}")
    successful = summary.optimized + summary.reused
    if not targets:
        summary.status = "no_targets"
    elif successful == 0:
        summary.status = "failed"
        message = "all selected optimizer poses failed"
        if summary.errors:
            message += f": {summary.errors[0]}"
        raise OptimizationError(message, summary)
    elif summary.failed:
        summary.status = "partial"
    else:
        summary.status = "completed"
    return summary


# ═══════════════════════════════════════════════════════════════════════════════
# Manifest builder
# ═══════════════════════════════════════════════════════════════════════════════

def build_prepared_manifest(
    protein_data: dict | None,
    ligand_data: dict | None,
) -> dict:
    manifest = {"proteins": [], "ligands": []}

    def record(section: str, source: str, file_path: str | None):
        if not file_path:
            return
        path_obj = Path(file_path).resolve()
        if not path_obj.exists():
            return
        manifest[section].append({
            "name": path_obj.name,
            "path": str(path_obj),
            "filetype": path_obj.suffix.lstrip(".").lower() or "unknown",
            "source": source,
        })

    if protein_data:
        for original, prepared in (protein_data.get("protein_results") or {}).items():
            record("proteins", original, prepared)
        for key in ("meeko_outputs", "converter_protein_outputs"):
            for original, pdbqt in (protein_data.get(key) or {}).items():
                record("proteins", f"{original} [pdbqt]", pdbqt)
        for key in ("meeko_box_files", "converter_box_files"):
            for original, box_file in (protein_data.get(key) or {}).items():
                record("proteins", f"{original} [box]", box_file)

    if ligand_data:
        for original, formats in (ligand_data.get("converted_ligands") or {}).items():
            for fmt, path in formats.items():
                record("ligands", f"{original} [{fmt}]", path)
        for key in ("meeko_ligand_outputs", "converter_ligand_outputs"):
            for original, pdbqt in (ligand_data.get(key) or {}).items():
                record("ligands", f"{original} [pdbqt]", pdbqt)

    return manifest


# ═══════════════════════════════════════════════════════════════════════════════
# Progress monitor
# ═══════════════════════════════════════════════════════════════════════════════

def _monitor_progress(counters: dict, stop_event: threading.Event, poll_interval: float = 10.0):
    start_time = time.time()
    while not stop_event.is_set():
        stop_event.wait(poll_interval)
        if stop_event.is_set():
            break
        elapsed = time.time() - start_time
        mins, secs = divmod(int(elapsed), 60)
        hrs, mins = divmod(mins, 60)
        time_str = f"{hrs}h{mins:02d}m{secs:02d}s" if hrs else f"{mins}m{secs:02d}s"
        with counters["lock"]:
            done = counters["done"]
            failed = counters["failed"]
            skipped = counters["skipped"]
            running = counters["running"]
            total = counters["total"]
        processed = done + failed + skipped
        pct = f"{processed / total * 100:.0f}%" if total else "0%"
        print(f"  ⏱ {time_str} | {pct} ({processed}/{total}) | "
              f"✓ {done} ok | ✗ {failed} fail | ⊘ {skipped} skip | ⚡ {running} active")


# ═══════════════════════════════════════════════════════════════════════════════
# Main entry point
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="AutoDock Vina batch docking pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-c", "--config",
        default=str(_SCRIPT_DIR / "docking_config.yaml"),
        help="Path to YAML configuration file (default: docking_config.yaml next to this script)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Validate configuration and inputs without running docking",
    )
    parser.add_argument(
        "--optimize",
        choices=["none", "smina", "gnina", "gnina_refinement", "all"],
        default=None,
        help="Post-dock pose optimisation + re-ranking tool (overrides the config "
             "'optimization' key). Poses are locally minimised and re-ranked on the "
             "optimiser's affinity while the original Vina rank is preserved. "
             "'gnina' preserves the historical CNN-rescore protocol; "
             "'gnina_refinement' requires gnina_cnn_scoring=refinement.",
    )
    args = parser.parse_args()

    # ── Load config ──
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)
    cfg = load_config(config_path)
    if args.optimize is not None:
        cfg["optimization"] = args.optimize
        _validate_optimizer_config(cfg)  # re-validate after the CLI override

    # ── Resolve paths relative to CWD ──
    receptors_dir = Path(cfg["receptors_dir"])
    output_base = Path(cfg["output_dir"])
    log_dir = Path(cfg.get("log_dir", "logs"))
    ligand_dirs = [Path(d) for d in cfg["ligand_dirs"]]

    output_base.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    # ── Validate input directories BEFORE creating any pdbqt subdirs. get_pdbqt_dir()
    # does mkdir(parents=True), which would silently materialise a mistyped
    # receptors_dir and make the existence check below unreachable. ──
    for d in ligand_dirs:
        if not d.exists():
            print(f"ERROR: Ligand directory missing: {d}")
            sys.exit(1)
    if not receptors_dir.exists():
        print(f"ERROR: Receptor directory missing: {receptors_dir}")
        sys.exit(1)

    protein_pdbqt_dir = get_pdbqt_dir(receptors_dir)

    print_config(cfg)

    cpu_model = get_cpu_model()
    print(f"CPU: {cpu_model}")

    # Optimizer dependencies and the explicitly pinned scoring/CNN model must
    # load before receptor preparation or a potentially multi-day docking run.
    opt_tools = resolve_optimizers(cfg)
    if opt_tools:
        print("\n── Optimizer pre-flight ──")
        try:
            optimizer_preflight = preflight_optimizers(cfg)
        except OptimizationError as e:
            print("ERROR: optimizer pre-flight failed:")
            for error in e.summary.errors:
                print(f"  ✗ {error}")
            sys.exit(1)
        for tool, details in optimizer_preflight["optimizers"].items():
            print(f"  ✓ {tool}: {details['path']} — {details['version']}")
        for converter, details in optimizer_preflight["converters"].items():
            print(f"  ✓ {converter}: {details['path']} — {details['version']}")

    # ── Determine converters ──
    prep_tool = cfg["prep_tool"]
    converters: list[tuple[str, str]] = []  # (name, converter_arg)
    if prep_tool in ("mgltools", "both"):
        converters.append(("mgl_tools", "mgltools"))
    if prep_tool in ("meeko", "both"):
        converters.append(("meeko", "meeko"))

    # ── Prepare proteins ──
    all_protein_outputs: list[tuple[str, dict]] = []
    for conv_name, conv_arg in converters:
        print(f"\n{'─' * 80}")
        print(f"Preparing proteins with: {conv_name}")
        print(f"{'─' * 80}")
        protein_outputs = run_workflow(
            input_dir=receptors_dir,
            contains="proteins",
            output_dir=protein_pdbqt_dir,
            skip_pdb_validation=cfg.get("skip_pdb_validation", False),
            custom_postfix=f"_{conv_name}",
            process_postfixes=cfg.get("process_postfixes", False),
            repair_terminals=cfg.get("repair_terminals", False),
            converter=conv_arg,
            convert_proteins=True,
            log_file=log_dir / f"protein_conversion_overview_{conv_name}.txt",
        )

        # ── A1: append KEEP-list HETATMs (cofactors + metals) so AutoDock
        # avoids clashing with the cofactor that PoseBusters later validates
        # against. No-op when retain_hetatm_residues is None/false.
        if cfg.get("retain_hetatm_residues"):
            print(f"\n  ── HETATM injection ({conv_name}) ──")
            _inject_hetatms_into_protein_outputs(
                protein_outputs, cfg["retain_hetatm_residues"], label=conv_name,
            )

        all_protein_outputs.append((conv_name, protein_outputs))

    if args.dry_run:
        print("\n[DRY RUN] Configuration validated. Protein preparation done. Exiting before docking.")
        sys.exit(0)

    # ── Dock per ligand directory × converter ──
    all_vina_results: List[DockingResult] = []

    for drugs_dir in ligand_dirs:
        for converter_name, protein_outputs in all_protein_outputs:
            # Use ligand directory name as subfolder for results
            result_subfolder = f"{drugs_dir.name}_{converter_name}"
            vina_output_dir = output_base / result_subfolder
            vina_output_dir.mkdir(parents=True, exist_ok=True)

            ligand_pdbqt_dir = get_pdbqt_dir(drugs_dir)

            print(f"\n{'=' * 80}")
            print(f"Processing: {drugs_dir.name}  |  Converter: {converter_name}")
            print(f"Ligand PDBQT dir: {ligand_pdbqt_dir}")
            print(f"Output: {vina_output_dir}")
            print(f"{'=' * 80}")

            # ── Check for existing ligand PDBQTs or prepare them ──
            existing_lig_pdbqts = sorted(ligand_pdbqt_dir.glob("*_ligand*.pdbqt"))
            if (ligand_pdbqt_dir / "ligands").exists():
                existing_lig_pdbqts += sorted((ligand_pdbqt_dir / "ligands").glob("*.pdbqt"))
            if not existing_lig_pdbqts:
                existing_lig_pdbqts = sorted(
                    p for p in ligand_pdbqt_dir.glob("*.pdbqt")
                    if not p.with_suffix(".box.txt").exists()
                    and not Path(str(p).replace(".pdbqt", ".box.txt")).exists()
                )

            if existing_lig_pdbqts:
                print(f"✓ Found {len(existing_lig_pdbqts)} existing ligand PDBQT files — skipping preparation.")
                ligand_outputs = {
                    "meeko_ligand_outputs": {str(p): str(p) for p in existing_lig_pdbqts},
                    "converted_ligands": {},
                    "ligand_files": set(),
                    "protein_files": set(),
                    "protein_results": {},
                    "fixed_proteins": [],
                    "log_entries": [],
                }
                for p in existing_lig_pdbqts[:10]:
                    print(f"  • {p.name}")
                if len(existing_lig_pdbqts) > 10:
                    print(f"  ... and {len(existing_lig_pdbqts) - 10} more")
            else:
                print(f"No existing ligand PDBQTs — running Meeko preparation...")
                ligand_outputs = run_workflow(
                    input_dir=drugs_dir,
                    contains="ligands",
                    output_dir=ligand_pdbqt_dir,
                    process_postfixes=False,
                    convert_ligands_with_meeko=True,
                    log_file=log_dir / f"ligand_conversion_{drugs_dir.name}.txt",
                )

            # ── Build manifest ──
            prepared_manifest = build_prepared_manifest(protein_outputs, ligand_outputs)
            print(f"Proteins prepared: {len(prepared_manifest['proteins'])}")
            print(f"Ligands prepared: {len(prepared_manifest['ligands'])}")

            # ── Pre-flight check ──
            vina_bin = Path(cfg["vina_bin"])
            errors = []
            if not vina_bin.exists():
                errors.append(f"✗ Vina executable not found: {vina_bin}")
            elif not os.access(vina_bin, os.X_OK):
                errors.append(f"✗ Vina executable not executable: {vina_bin}")

            lig_files = collect_files(ligand_pdbqt_dir, [".pdbqt"])
            print(f"  Ligands: {len(lig_files)} PDBQT files in {ligand_pdbqt_dir}")
            if not lig_files:
                errors.append(f"✗ No ligand files in {ligand_pdbqt_dir}")

            rec_files = collect_files(protein_pdbqt_dir, [".pdb", ".pdbqt"])
            print(f"  Receptors: {len(rec_files)} files in {protein_pdbqt_dir}")
            if not rec_files:
                errors.append(f"✗ No receptor files in {protein_pdbqt_dir}")

            if errors:
                print("\n⚠ ISSUES FOUND:")
                for e in errors:
                    print(f"  {e}")
                print("Skipping this directory.")
                continue

            print("✓ All checks passed — ready to dock.")

            # ── Precompute properties ──
            lig_props_cache, prot_props_cache = precompute_properties(rec_files, lig_files)

            # ── Monitor thread ──
            progress_counters = {
                "done": 0, "failed": 0, "skipped": 0, "running": 0, "total": 0,
                "lock": threading.Lock(),
            }
            stop_monitor = threading.Event()
            monitor_thread = threading.Thread(
                target=_monitor_progress,
                args=(progress_counters, stop_monitor, 15),
                daemon=True,
            )
            monitor_thread.start()

            # ── Run docking ──
            try:
                results_df, docking_results = run_autodock_vina(
                    base_dir=vina_output_dir,
                    log_dir=log_dir,
                    prepared_manifest=prepared_manifest,
                    cfg=cfg,
                    cpu_model=cpu_model,
                    protein_workflow_data=protein_outputs,
                    ligand_workflow_data=ligand_outputs,
                    lig_props_cache=lig_props_cache,
                    prot_props_cache=prot_props_cache,
                    progress_counters=progress_counters,
                )
            finally:
                stop_monitor.set()
                monitor_thread.join(timeout=5)

            all_vina_results.extend(docking_results)

            # ── Optional gnina/smina re-optimisation + re-ranking of the Vina
            # poses (no-op unless 'optimization' is set). Runs per output
            # subfolder so its optimization_log.csv sits beside the poses. ──
            overwrite_opt = cfg.get("overwrite_existing", False) or cfg.get("overwrite_poses", False)
            optimization_summary = optimize_autodock_results(
                docking_results, cfg, vina_output_dir, overwrite=overwrite_opt,
            )

            summary = generate_summary(docking_results, cfg)
            summary["optimization"] = optimization_summary.to_dict()
            print_summary(summary)

            summary_path = vina_output_dir / "docking_summary.json"
            with open(summary_path, "w") as f:
                json.dump(summary, f, indent=2, default=str)
            print(f"\nSummary saved: {summary_path}")

    print(f"\n{'=' * 80}")
    print(f"ALL DONE — {len(all_vina_results)} total results across "
          f"{len(ligand_dirs)} ligand directories × {len(all_protein_outputs)} converters")
    print(f"{'=' * 80}")


if __name__ == "__main__":
    main()
