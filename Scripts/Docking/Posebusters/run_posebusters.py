"""PoseBusters Parallel Pose Validation Pipeline.

Validates docked poses from one or more docking methods using PoseBusters.
Supports AutoDock Vina, DiffDock, and EquiBind outputs.

Run as a script:
    python pose_busters_para_refactored.py --config config.yaml

Use as a library:
    from pose_busters_para_refactored import run_pipeline
    run_pipeline("config.yaml")
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from matplotlib.colors import LinearSegmentedColormap
from posebusters import PoseBusters

# The receptor-prep step injects this curated cofactor/metal set into the
# docking PDBQT (see Scripts/Utilities/inject_hetatms.py). Importing it lets
# the validation strip step auto-derive "keep exactly what the docker saw".
# Import is best-effort: if the project isn't importable (e.g. run standalone),
# auto keep-mode is simply unavailable and an explicit keep list is required.
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
    # HETATM residue names to strip from each receptor PDB before validation
    # (e.g. crystallisation artefacts that the docker never saw but PoseBusters
    # otherwise checks against). Stored as a set of upper-cased 3-letter codes.
    strip_hetatm_residues: set[str] = field(default_factory=set)
    # Auto-derive mode: keep ONLY these HETATM residues and strip every other
    # HETATM from each receptor PDB before validation. Set this to exactly the
    # residues the docker saw (the injected cofactors/metals) so PoseBusters can
    # only ever flag clashes against atoms Vina was actually given — eliminating
    # false clashes from buffers/glycans/cognate-ligand remnants the docker was
    # blind to. ``None`` disables keep-mode and falls back to
    # ``strip_hetatm_residues`` (explicit strip list). Takes precedence when set.
    keep_hetatm_residues: set[str] | None = None

    # Derived
    output_dir: Path = field(init=False)
    converted_dir: Path = field(init=False)

    def __post_init__(self):
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        self.output_dir = self.output_base_dir / self.config_mode
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.converted_dir = self.output_dir / "converted_pdbqt"
        self.converted_dir.mkdir(parents=True, exist_ok=True)


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

    output_base_dir = _resolve(raw.get("output_base_dir", "posebusters_results"))

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

    return PipelineConfig(
        receptors_dir=receptors_dir,
        docking_directories=docking_directories,
        work_dir=work_dir,
        output_base_dir=output_base_dir,
        config_mode=raw.get("config_mode", "dock"),
        save_interval=int(raw.get("save_interval", 100)),
        overwrite=bool(raw.get("overwrite", False)),
        num_workers=raw.get("num_workers"),
        ligand_template_dirs=template_dirs,
        make_plots=bool(raw.get("make_plots", True)),
        copy_proved_poses=bool(raw.get("copy_proved_poses", True)),
        strip_hetatm_residues=strip_set,
        keep_hetatm_residues=keep_set,
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
    in the set (keep-mode — mirror the docker's injected receptor).
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
    residues survive, every other HETATM is stripped — so PoseBusters validates
    against exactly the atoms the docker saw. Otherwise the legacy explicit
    strip list (``ctx.strip_hetatm_residues``) is applied. If neither is
    configured, the original directory is returned unchanged.
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
        if dropped:
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
    # EquiBind variant provenance (never a pass/fail test column).
    "pocket_source", "pocket_id", "clamp_variant", "refine_variant",
    "smina_affinity",
    # AutoDock Vina native rank + affinity.
    "autodock_rank", "autodock_affinity",
    # DiffDock post-pose optimizer provenance (original / smina / gnina).
    "optimizer",
}

_EXCLUDE_COLS = {
    "mol_true_loaded", "mol_cond_loaded",
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
# Protein/cofactor/water columns are only present in 'dock' mode; the
# intersection with the actual columns handles 'mol' mode automatically.
CANONICAL_TEST_COLUMNS = (
    "mol_pred_loaded", "sanitization", "inchi_convertible",
    "all_atoms_connected", "no_radicals",
    "bond_lengths", "bond_angles", "internal_steric_clash",
    "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
    "double_bond_flatness", "internal_energy",
    "minimum_distance_to_protein",
    "minimum_distance_to_organic_cofactors",
    "minimum_distance_to_inorganic_cofactors",
    "minimum_distance_to_waters",
    "volume_overlap_with_protein",
    "volume_overlap_with_organic_cofactors",
    "volume_overlap_with_inorganic_cofactors",
    "volume_overlap_with_waters",
)


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


def coerce_test_cols_to_bool(df: pd.DataFrame, test_cols: list[str]) -> None:
    """Convert string-encoded booleans to actual bool dtype in-place."""
    for tc in test_cols:
        if df[tc].dtype == object:
            df[tc] = df[tc].map(_BOOL_STR_MAP)
        df[tc] = df[tc].astype(bool)


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
    return _strip_suffixes(name, _PROTEIN_SUFFIXES)


def _normalize_ligand(name: str) -> str:
    return _strip_suffixes(name, _LIGAND_SUFFIXES)


def _collect_autodock_rows(poses_dir: Path) -> list[dict]:
    rows = []
    for pdbqt_file in Path(poses_dir).glob("*_vina_out.pdbqt"):
        stem = pdbqt_file.stem.replace("_vina_out", "")
        parts = stem.split("__")
        if len(parts) == 2:
            protein, ligand = parts
            with open(pdbqt_file) as f:
                pose_count = f.read().count("MODEL") or 1
            rows.append({
                "docking_tool": "autodock",
                "protein": _normalize_protein(protein),
                "ligand": _normalize_ligand(ligand),
                "file_path": str(pdbqt_file), "pose_count": pose_count,
            })
    return rows


def _collect_diffdock_rows(poses_dir: Path) -> list[dict]:
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
            pose_count = len(list(subdir.glob("**/*.sdf")))
            if pose_count > 0:
                rows.append({
                    "docking_tool": "diffdock",
                    "protein": _normalize_protein(protein),
                    "ligand": _normalize_ligand(ligand),
                    "file_path": str(subdir), "pose_count": pose_count,
                })
    return rows


def _collect_equibind_rows(poses_dir: Path) -> list[dict]:
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
            pose_count = len([f for f in subdir.glob("**/*.sdf") if "prep" not in f.parts])
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
    if not file_path.exists() or file_path.suffix != ".pdbqt":
        return []
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
            "pose_name": f"{row['protein']}__{row['ligand']}_pose{model_num}",
            "file_format": "sdf",
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


def _expand_diffdock_poses(row: dict, _conv_dir: Path, _ctx: PipelineConfig) -> list[dict]:
    file_path = Path(row["file_path"])
    if not file_path.is_dir():
        return []
    poses = []
    for sdf_file in sorted(file_path.glob("**/*.sdf")):
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
_EQ_PROVENANCE_TAGS = ("pocket_source", "pocket_id", "clamp_variant",
                       "refine_variant", "smina_affinity")


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
    if "__refraw" in n:
        return "raw"
    return None


def _expand_equibind_poses(row: dict, _conv_dir: Path, _ctx: PipelineConfig) -> list[dict]:
    file_path = Path(row["file_path"])
    if not file_path.is_dir():
        return []
    poses = []
    for sdf_file in sorted(f for f in file_path.glob("**/*.sdf") if "prep" not in f.parts):
        try:
            rel = str(sdf_file.relative_to(file_path))
        except ValueError:
            rel = sdf_file.name
        fname = sdf_file.name
        # Full per-pose provenance so every EquiBind variant axis is analysable
        # downstream without re-parsing filenames: pocket source (fpocket /
        # p2rank / unguided), centroid-clamp and re-search variants, the matched
        # pocket id and smina affinity. SDF tags win; the filename is the fallback.
        tags = _read_sdf_tags(sdf_file, _EQ_PROVENANCE_TAGS)
        pose = {
            "method": row["docking_tool"], "protein": row["protein"], "ligand": row["ligand"],
            "pose_file": str(sdf_file),
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
        if tags.get("smina_affinity"):
            try:
                pose["smina_affinity"] = float(tags["smina_affinity"])
            except ValueError:
                pass
        poses.append(pose)
    return poses


POSE_EXPANDERS = {
    "autodock": _expand_autodock_poses,
    "diffdock": _expand_diffdock_poses,
    "equibind_guided": _expand_equibind_poses,
    "equibind_exclusion": _expand_equibind_poses,
    "equibind_docked_poses": _expand_equibind_poses,
}


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
    """Search for a template SDF for *ligand_name* in configured directories."""
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
            for match in sdir.rglob(pat):
                try:
                    mol = Chem.MolFromMolFile(str(match), removeHs=True, sanitize=True)
                    if mol is not None:
                        print(f"    Using bond-order template: {match.relative_to(sdir)}")
                        return mol
                except Exception:
                    pass
    return None


_MK_EXPORT_BIN: str | None | bool = False  # False = not yet looked up


def _mk_export_available() -> str | None:
    """Locate the Meeko ``mk_export.py`` CLI once and cache the result."""
    global _MK_EXPORT_BIN
    if _MK_EXPORT_BIN is False:
        _MK_EXPORT_BIN = shutil.which("mk_export.py") or shutil.which("mk_export")
    return _MK_EXPORT_BIN


def _mk_export_model(model_pdbqt: Path, output_sdf: Path) -> bool:
    """Convert a single-model PDBQT to *output_sdf* with Meeko ``mk_export``.

    Returns True only if a chemically sane, radical-free molecule was written;
    otherwise removes any partial output and returns False so the caller can
    fall back to the template/obabel strategies.
    """
    mk_bin = _mk_export_available()
    if not mk_bin:
        return False
    from rdkit import Chem

    try:
        res = subprocess.run(
            [mk_bin, str(model_pdbqt), "-s", str(output_sdf)],
            capture_output=True, text=True,
        )
    except OSError:
        return False
    if res.returncode != 0 or not output_sdf.exists():
        output_sdf.unlink(missing_ok=True)
        return False

    # Meeko reconstructs bond orders from the embedded SMILES, but guard anyway:
    # reject anything that won't sanitize or carries radicals (the exact failure
    # mode obabel exhibits) so a bad export can't masquerade as a valid pose.
    mol = next(iter(Chem.SDMolSupplier(str(output_sdf), removeHs=False, sanitize=True)), None)
    if mol is None or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        output_sdf.unlink(missing_ok=True)
        return False
    return True


def _convert_pdbqt_to_sdf(pdbqt_file: str, out_dir: Path, ctx: PipelineConfig) -> list[str]:
    """Convert a multi-pose PDBQT to individual SDF files."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    pdbqt_path = Path(pdbqt_file)
    base_name = pdbqt_path.stem
    models = _split_pdbqt_models(pdbqt_file)
    if not models:
        print(f"  Warning: No models found in {pdbqt_file}")
        return []

    stem_clean = base_name.replace("_vina_out", "")
    parts = stem_clean.split("__")
    template_mol = _find_template_mol(parts[1], ctx) if len(parts) == 2 else None

    converted: list[str] = []
    for i, model_content in enumerate(models, start=1):
        output_sdf = out_dir / f"{base_name}_model{i}.sdf"
        if output_sdf.exists():
            converted.append(str(output_sdf))
            continue

        temp_pdbqt = out_dir / f"{base_name}_model{i}.pdbqt"
        temp_pdb = out_dir / f"{base_name}_model{i}.pdb"
        with open(temp_pdbqt, "w") as f:
            f.write(model_content)

        mol_final = None

        # Strategy 0 (preferred): Meeko mk_export. Ligand prep embeds the input
        # SMILES (`REMARK SMILES`) in the PDBQT, so Meeko rebuilds the exact
        # bond orders and formal charges instead of guessing. This avoids the
        # radical / valence artefacts obabel produces when it has to infer bond
        # orders from coordinates alone (see _mk_export_model). Falls through to
        # the template/obabel strategies when the SMILES header is absent (e.g.
        # a non-Meeko-prepped ligand) or mk_export is unavailable.
        if "REMARK SMILES" in model_content and _mk_export_model(temp_pdbqt, output_sdf):
            converted.append(str(output_sdf))
            temp_pdbqt.unlink(missing_ok=True)
            continue

        # Strategy 1: obabel PDBQT -> PDB -> RDKit + template
        if template_mol is not None:
            try:
                res = subprocess.run(
                    ["obabel", str(temp_pdbqt), "-O", str(temp_pdb)],
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
                res = subprocess.run(
                    ["obabel", str(temp_pdbqt), "-O", str(output_sdf)],
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
                            converted.append(str(output_sdf))
                            temp_pdbqt.unlink(missing_ok=True)
                            continue
                    else:
                        converted.append(str(output_sdf))
                        temp_pdbqt.unlink(missing_ok=True)
                        continue
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
            if output_sdf.exists():
                converted.append(str(output_sdf))
        elif not output_sdf.exists():
            print(f"    Warning: Could not convert model {i} from {pdbqt_path.name}")

        temp_pdbqt.unlink(missing_ok=True)

    return converted


def collect_all_pose_files(filtered_df: pd.DataFrame, ctx: PipelineConfig) -> list[dict]:
    """Expand every row in *filtered_df* into individual pose-file dicts."""
    all_poses: list[dict] = []
    for _, row in filtered_df.iterrows():
        expander = POSE_EXPANDERS.get(row["docking_tool"])
        if expander is None:
            print(f"  WARNING: No pose expander for '{row['docking_tool']}', skipping.")
            continue
        all_poses.extend(expander(row.to_dict(), ctx.converted_dir, ctx))
    return all_poses


# ============================================================================
# PARALLEL WORKER FUNCTIONS  (module-level so they're picklable)
# ============================================================================

_worker_buster_dock: PoseBusters | None = None
_worker_buster_mol: PoseBusters | None = None
_worker_protein_cache: dict[str, str | None] = {}
_worker_config_mode: str = "mol"


def _init_worker(config_mode: str, protein_cache: dict[str, str | None]) -> None:
    global _worker_buster_dock, _worker_buster_mol, _worker_protein_cache, _worker_config_mode
    from posebusters import PoseBusters as _PB

    _silence_rdkit_2d3d_warning()  # workers run bust(); install the filter here too

    _worker_config_mode = config_mode
    _worker_protein_cache = protein_cache or {}
    if config_mode == "dock":
        _worker_buster_dock = _PB(config="dock")
    _worker_buster_mol = _PB(config="mol")


def _process_single_pose(pose_info: dict):
    """Validate one pose; returns (result_df | None, error_msg | None)."""
    pose_file = pose_info["pose_file"]
    protein_name = pose_info["protein"]

    if not os.path.exists(pose_file):
        return None, f"File not found: {pose_file}"

    try:
        protein_file = None
        used_mode = _worker_config_mode

        if _worker_config_mode == "dock":
            protein_file = _worker_protein_cache.get(protein_name)
            if protein_file and Path(protein_file).exists():
                df = _worker_buster_dock.bust(pose_file, None, protein_file, full_report=True)
                used_mode = "dock"
            else:
                df = _worker_buster_mol.bust(pose_file, None, None, full_report=True)
                used_mode = "mol (fallback)"
        else:
            df = _worker_buster_mol.bust(pose_file, None, None, full_report=True)
            used_mode = "mol"

        df["docking_method"] = pose_info["method"]
        df["protein"] = protein_name
        df["ligand"] = pose_info["ligand"]
        df["pose_file"] = pose_info["pose_file"]
        df["pose_name"] = pose_info["pose_name"]
        df["file_format"] = pose_info.get("file_format", "sdf")
        df["protein_file_used"] = protein_file or "none"
        df["posebusters_mode"] = used_mode
        if "confidence" in pose_info:
            df["diffdock_confidence"] = pose_info["confidence"]
        # DiffDock optimizer provenance (original / smina / gnina) so optimised
        # poses stay distinguishable from the raw DiffDock output downstream.
        if "optimizer" in pose_info:
            df["optimizer"] = pose_info["optimizer"]
        # AutoDock Vina rank + affinity (kcal/mol). Vina produces them natively but
        # they were never propagated to the CSV; carry them through so AutoDock is
        # rankable from the results table like DiffDock/EquiBind.
        for key in ("autodock_rank", "autodock_affinity"):
            if key in pose_info:
                df[key] = pose_info[key]
        # Carry EquiBind variant provenance through to the results CSV so every
        # pose axis (pocket source, clamp, re-search, pocket id, smina affinity)
        # is analysable without re-parsing filenames.
        for key in _EQ_PROVENANCE_TAGS:
            if key in pose_info:
                df[key] = pose_info[key]
        return df, None

    except Exception as e:
        return None, f"Error processing {Path(pose_file).name}: {e}"


# ============================================================================
# CHECKPOINT + ANALYSIS
# ============================================================================

def _save_checkpoint(frames, output_file, n_new, n_pending, final=False):
    df_out = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    df_out.to_csv(output_file, index=False)
    tag = "FINAL" if final else "CHECKPOINT"
    print(f"    [{tag}] {n_new}/{n_pending} new  |  total rows: {len(df_out)}  →  {output_file.name}")


def analyze_poses_with_posebusters(
    poses_list: list[dict],
    protein_cache: dict[str, str | None],
    config: str = "mol",
    output_file: Path | None = None,
    save_interval: int = 100,
    overwrite: bool = False,
    num_workers: int | None = None,
) -> pd.DataFrame:
    """Run PoseBusters on *poses_list* using multiprocessing.Pool."""
    if not poses_list:
        # No new poses — either nothing was discovered or every complex was
        # already validated and skipped upstream. Return the existing results so
        # downstream steps (copy-proved-poses, plots) still run on the full set.
        if output_file is not None and output_file.exists() and not overwrite:
            print("No new poses to analyze — all complexes already validated.")
            return pd.read_csv(output_file)
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
            if "pose_file" in previous_results.columns:
                already_done = set(previous_results["pose_file"].dropna().astype(str))
            print(f"  Loaded {len(previous_results)} previous results; {len(already_done)} poses done")

    pending = [p for p in poses_list if str(p["pose_file"]) not in already_done]
    n_skipped = len(poses_list) - len(pending)
    if n_skipped:
        print(f"  Skipping {n_skipped} already-inspected poses")
    if not pending:
        print("  All poses already inspected.")
        return previous_results

    n_workers = num_workers if num_workers is not None else cpu_count()
    n_workers = min(n_workers, len(pending))
    print(f"\n  Workers: {n_workers} ({cpu_count()} CPUs available)")

    new_frames: list[pd.DataFrame] = []
    total = len(pending)
    n_new = n_errors = n_dock = n_mol_fb = 0
    print(f"  {total} poses to process (checkpoint every {save_interval})...\n")

    def _checkpoint(final: bool = False) -> None:
        if not output_file:
            return
        frames = ([previous_results] if not previous_results.empty else []) + new_frames
        if frames:
            _save_checkpoint(frames, output_file, n_new, total, final=final)

    with Pool(processes=n_workers, initializer=_init_worker,
              initargs=(config, dict(protein_cache))) as pool:
        for result_df, error_msg in pool.imap_unordered(_process_single_pose, pending):
            if error_msg:
                n_errors += 1
                print(f"  {error_msg}")
                continue
            if result_df is not None:
                mode_val = result_df["posebusters_mode"].iloc[0] if "posebusters_mode" in result_df.columns else ""
                if mode_val == "dock":
                    n_dock += 1
                elif "fallback" in str(mode_val):
                    n_mol_fb += 1
                new_frames.append(result_df)
                n_new += 1

            if n_new % 10 == 0 or n_new == 1:
                print(f"  Processed {n_new}/{total}  (overall {len(already_done) + n_new}/{len(poses_list)})")
            if n_new > 0 and n_new % save_interval == 0:
                _checkpoint()

    _checkpoint(final=True)

    print(f"\n  Done:  new={n_new}  prev={n_skipped}  total={len(already_done)+n_new}  errors={n_errors}")
    if config == "dock":
        print(f"    dock mode: {n_dock}  |  mol fallback: {n_mol_fb}")

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

    for method_key, poses_dir in ctx.docking_directories.items():
        collector = ROW_COLLECTORS.get(method_key)
        if collector is None:
            print(f"  WARNING: No row collector for '{method_key}', skipping.")
            continue
        if not poses_dir.exists():
            print(f"  WARNING: Directory not found for '{method_key}': {poses_dir}, skipping.")
            continue

        method_rows = collector(poses_dir)
        for r in method_rows:
            r["docking_tool"] = method_key
        all_rows.extend(method_rows)

        df_m = pd.DataFrame(method_rows)
        if not df_m.empty:
            col_name = method_key.replace(" ", "_").title()
            summary = df_m.groupby(["protein", "ligand"])["pose_count"].sum().reset_index()
            summary.columns = ["Protein", "Ligand", col_name]
            method_dfs[col_name] = summary
        print(f"  {method_key}: {len(method_rows)} protein-ligand combinations")

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

        print("\n" + "=" * 100)
        print("POSE COUNTS BY PROTEIN-LIGAND COMBINATION")
        print("=" * 100)
        print(df_combined.to_string(index=False))

        print("\n" + "=" * 100)
        print("SUMMARY STATISTICS")
        print("=" * 100)
        print(f"\nTotal unique protein-ligand combinations: {len(df_combined)}")
        print(f"\nTotal poses by method:")
        for col in method_cols:
            print(f"  {col}: {df_combined[col].sum():,} poses")
        print(f"  Combined Total: {df_combined['Total'].sum():,} poses")
        print(f"\nUnique proteins: {df_combined['Protein'].nunique()}")
        print(f"Unique ligands: {df_combined['Ligand'].nunique()}")

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


def filter_common_combos(filtered_poses_df: pd.DataFrame) -> pd.DataFrame:
    """Keep only protein-ligand combinations present across ALL methods."""
    print(f"Total entries: {len(filtered_poses_df)}")
    if filtered_poses_df.empty:
        print("filtered_poses_df is empty. Check docking_directories paths.")
        return filtered_poses_df

    print(f"Methods found: {filtered_poses_df['docking_tool'].unique().tolist()}")
    combos_by_method = {
        m: set(zip(g["protein"], g["ligand"]))
        for m, g in filtered_poses_df.groupby("docking_tool")
    }

    print("=" * 100)
    print("FILTERING: Keep Only Protein-Ligand Combinations in ALL Methods")
    print("=" * 100)
    for method, combos in combos_by_method.items():
        print(f"  {method}: {len(combos)} combinations")

    common_combos = set.intersection(*combos_by_method.values()) if combos_by_method else set()
    print(f"\nProtein-ligand combinations in ALL methods: {len(common_combos)}")

    if not common_combos:
        print("\nWARNING: No combinations found in all methods!")
        return filtered_poses_df.iloc[0:0]

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

    def _is_done(row) -> bool:
        key = (str(row["docking_tool"]), str(row["protein"]), str(row["ligand"]))
        return done_counts.get(key, 0) >= row["pose_count"]

    done_mask = filtered_poses_df.apply(_is_done, axis=1)
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
        print(f"  WARNING: {n_missing} proteins have no PDB — will fall back to 'mol' mode")
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
    output_base = ctx.work_dir / "posebuster_proved" / ctx.config_mode
    folders = {key: output_base / key for key in ctx.docking_directories}

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

    test_cols = identify_test_columns(df_pb)
    if not test_cols:
        raise ValueError("Cannot identify PoseBusters test columns.")
    coerce_test_cols_to_bool(df_pb, test_cols)

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

    print(f"\nBottleneck tests (lowest pass rates):")
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
    print(f"\nPassed poses by method:")
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
        dest_name = (
            f"{protein_name}__{ligand_name}__{src_path.name}"
            if (protein_name and ligand_name and "__" not in src_path.name)
            else src_path.name
        )
        dest_file = dest_dir / dest_name

        if dest_file.exists():
            stem, suffix = dest_file.stem, dest_file.suffix
            counter = 1
            while dest_file.exists():
                dest_file = dest_dir / f"{stem}_dup{counter}{suffix}"
                counter += 1

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

    test_cols = identify_test_columns(df)
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

def run_pipeline(config_path: str | Path, overwrite: bool | None = None) -> pd.DataFrame:
    """Run the full PoseBusters validation pipeline driven by a YAML config.

    *overwrite* overrides the value from the config file when provided.
    """
    ctx = load_config(config_path)
    if overwrite is not None:
        ctx.overwrite = overwrite

    # Optional: write cleaned receptor PDBs (HETATM artefacts removed) and
    # validate against those instead of the originals.
    ctx.receptors_dir = prepare_cleaned_receptor_dir(ctx)

    file_map, normalized_map = discover_proteins(ctx.receptors_dir)

    filtered_poses_df, _ = collect_pose_rows(ctx)
    filtered_poses_df = filter_common_combos(filtered_poses_df)
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

    # Skip complexes already fully validated so the expensive pose collection /
    # PDBQT->SDF conversion below only touches complexes not yet run.
    filtered_poses_df = filter_already_processed_complexes(
        filtered_poses_df, results_csv_path, overwrite=ctx.overwrite)

    print(f"\nFiltered subset: {len(filtered_poses_df)} combinations, "
          f"~{filtered_poses_df['pose_count'].sum() if not filtered_poses_df.empty else 0:,} poses")

    print("\n" + "-" * 80 + "\n1. Collecting individual pose files...\n" + "-" * 80)
    all_poses = collect_all_pose_files(filtered_poses_df, ctx)
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
    )

    if not results_df.empty:
        print(f"\n   Results saved to: {results_csv_path}")
        print(f"   Total rows: {len(results_df)}")
        if "posebusters_mode" in results_df.columns:
            print(f"\n   Mode breakdown:")
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
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_cli()
    run_pipeline(args.config, overwrite=args.overwrite)
