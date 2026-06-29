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
import itertools
import json
import os
import re
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(error_log, f, indent=2)


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
        props["lig_rotatable_bonds"] = rotatable // 2 if rotatable > 0 else None
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
            try:
                lig_props_cache[lig] = fut.result()
            except Exception as e:
                print(f"  ⚠ Failed computing props for {lig.name}: {e}")
                lig_props_cache[lig] = get_ligand_properties_from_pdbqt(lig)
    print(f"  {len(lig_props_cache)} ligands analysed")

    print("Pre-computing protein properties...")
    with ProcessPoolExecutor(max_workers=min(n_workers, len(proteins) or 1)) as pool:
        futures = {pool.submit(get_protein_properties, prot): prot for prot in proteins}
        for fut in as_completed(futures):
            prot = futures[fut]
            try:
                prot_props_cache[prot] = fut.result()
            except Exception as e:
                print(f"  ⚠ Failed computing props for {prot.name}: {e}")
                prot_props_cache[prot] = {"prot_num_residues": None, "prot_num_atoms": None, "prot_num_chains": None}
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
    rows = []
    replaced = 0
    for r in results:
        combo = f"{r.ligand_name}__{r.protein_name}"
        if combo in existing_combos and not overwrite:
            continue
        lig_props = lig_props_cache.get(r.ligand_path, get_ligand_properties(r.ligand_path))
        prot_props = prot_props_cache.get(r.protein_path, get_protein_properties(r.protein_path))
        rows.append(_build_log_row(r, lig_props, prot_props, cfg, cpu_model))
        if combo in existing_combos:
            replaced += 1
        existing_combos.add(combo)

    if not rows:
        return

    if replaced > 0 and log_path.exists():
        df_existing = pd.read_csv(log_path)
        replace_combos = {row["combo_name"] for row in rows}
        df_existing = df_existing[~df_existing["combo_name"].isin(replace_combos)]
        df_new = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
        frames = [df for df in (df_existing, df_new) if not df.empty]
        df_all = pd.concat(frames, ignore_index=True) if frames else df_new
        df_all.to_csv(log_path, index=False)
        print(f"  📝 Replaced {replaced} + appended {len(rows) - replaced} entries in {log_path.name}")
    else:
        df = pd.DataFrame(rows, columns=DOCKING_LOG_COLUMNS)
        df.to_csv(log_path, mode="a", header=False, index=False)
        print(f"  📝 Appended {len(rows)} entries to {log_path.name}")


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

        if not overwrite_error_log and combo_name in error_log:
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="skipped",
                scoring_function=scoring,
                error_message=f"Previously failed: {error_log[combo_name].get('error', 'unknown')[:100]}",
            )
            results_out.append({
                "result": r,
                "raw": {"receptor": str(rec_path), "ligand": str(lig_path),
                        "output": str(out_path), "log": str(vina_log_path),
                        "status": "skip-failed"},
                "combo_name": combo_name, "ran_vina": False,
            })
            continue

        if out_path.exists() and not overwrite:
            n_poses = 0
            best_aff = None
            if vina_log_path.exists():
                modes = parse_vina_affinities(vina_log_path.read_text())
                n_poses = len(modes)
                if modes:
                    best_aff = modes[0]["affinity"]
            r = DockingResult(
                protein_name=rec_name, ligand_name=lig_name,
                protein_path=rec_path, ligand_path=lig_path,
                output_dir=dock_dir, status="skipped",
                scoring_function=scoring,
                num_poses=n_poses, best_affinity=best_aff,
                pose_files=[out_path],
            )
            results_out.append({
                "result": r,
                "raw": {"receptor": str(rec_path), "ligand": str(lig_path),
                        "output": str(out_path), "log": str(vina_log_path),
                        "status": "skip-exists"},
                "combo_name": combo_name, "ran_vina": False,
            })
            continue

        ligands_to_dock.append(lig)

    if not ligands_to_dock:
        return results_out

    batch_out_dir = dock_dir / f"_batch_{rec_name}_{scoring}"
    batch_out_dir.mkdir(parents=True, exist_ok=True)
    batch_log_path = log_dir / f"batch_{rec_name}_{scoring}.log"

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

        docked_stems: set[str] = set()
        for lig in ligands_to_dock:
            lig_path = lig["pdbqt"]
            batch_out_file = batch_out_dir / f"{lig_path.stem}_out.pdbqt"
            if batch_out_file.exists() and batch_out_file.stat().st_size > 0:
                docked_stems.add(lig_path.stem)

        finished_ligs = [l for l in ligands_to_dock if l["pdbqt"].stem in docked_stems]
        remaining_ligs = [l for l in ligands_to_dock if l["pdbqt"].stem not in docked_stems]
        per_lig_time = total_elapsed / max(len(ligands_to_dock), 1)

        for lig in finished_ligs:
            lig_path = lig["pdbqt"]
            lig_name = _normalize_key(lig_path)
            combo_name = f"{lig_name}__{rec_name}"
            batch_out_file = batch_out_dir / f"{lig_path.stem}_out.pdbqt"
            stem = f"{rec_path.stem}__{lig_path.stem}_{scoring}_vina_out"
            final_out = dock_dir / f"{stem}.pdbqt"
            per_lig_log = log_dir / f"{stem}.log"

            batch_out_file.rename(final_out)
            per_lig_log.write_text(
                f"# Extracted from batch run (partial, timeout): {batch_log_path.name}\n"
                f"# Receptor: {rec_name}  Ligand: {lig_name}  Scoring: {scoring}\n\n"
            )
            modes = []
            try:
                modes = parse_pdbqt_affinities(final_out.read_text())
            except Exception:
                pass
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
            sub_timeout = max(effective_timeout or (timeout * len(remaining_ligs)), timeout) if timeout > 0 else 0
            print(f"    Retrying {len(remaining_ligs)} remaining ligands in {len(sub_batches)} sub-batches")

            common_kwargs = dict(
                scoring=scoring, vina_bin=vina_bin, dock_dir=dock_dir,
                log_dir=log_dir, exhaustiveness=exhaustiveness,
                num_modes=num_modes, energy_range=energy_range, cpu=cpu,
                seed=seed, timeout=timeout, overwrite=overwrite,
                overwrite_error_log=overwrite_error_log, error_log=error_log,
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

        if batch_out_file.exists() and batch_out_file.stat().st_size > 0:
            batch_out_file.rename(final_out)
            per_lig_log.write_text(
                f"# Extracted from batch run: {batch_log_path.name}\n"
                f"# Receptor: {rec_name}  Ligand: {lig_name}  Scoring: {scoring}\n\n"
            )
            modes = []
            try:
                pose_text = final_out.read_text()
                modes = parse_pdbqt_affinities(pose_text)
            except Exception:
                pass
            if not modes:
                modes = parse_vina_affinities(batch_output)
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
            err_msg = f"Batch docking failed (rc={retcode})" if not batch_ok else "No output produced"
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

    if not overwrite_error_log and combo_name in error_log:
        r = DockingResult(
            protein_name=rec_name, ligand_name=lig_name,
            protein_path=rec_path, ligand_path=lig_path,
            output_dir=dock_dir, status="skipped",
            scoring_function=scoring,
            error_message=f"Previously failed: {error_log[combo_name].get('error', 'unknown')[:100]}",
        )
        raw = {
            "receptor": str(rec_path), "ligand": str(lig_path),
            "output": str(out_path), "log": str(vina_log_path),
            "status": "skip-failed",
        }
        return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": False}

    if out_path.exists() and not overwrite:
        n_poses = 0
        best_aff = None
        if vina_log_path.exists():
            modes = parse_vina_affinities(vina_log_path.read_text())
            n_poses = len(modes)
            if modes:
                best_aff = modes[0]["affinity"]
        r = DockingResult(
            protein_name=rec_name, ligand_name=lig_name,
            protein_path=rec_path, ligand_path=lig_path,
            output_dir=dock_dir, status="skipped",
            scoring_function=scoring,
            num_poses=n_poses, best_affinity=best_aff,
            pose_files=[out_path],
        )
        raw = {
            "receptor": str(rec_path), "ligand": str(lig_path),
            "output": str(out_path), "log": str(vina_log_path),
            "status": "skip-exists",
        }
        return {"result": r, "raw": raw, "combo_name": combo_name, "ran_vina": False}

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
        "--out", str(out_path),
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
    status_str = "ok" if ok and out_path.exists() and out_path.stat().st_size > 0 else "fail"

    if status_str != "ok" and out_path.exists():
        try:
            out_path.unlink()
        except Exception:
            pass

    modes = parse_vina_affinities(output or "")
    n_poses = len(modes)
    best_aff = modes[0]["affinity"] if modes else None

    r = DockingResult(
        protein_name=rec_name, ligand_name=lig_name,
        protein_path=rec_path, ligand_path=lig_path,
        output_dir=dock_dir,
        status="success" if status_str == "ok" else "failed",
        scoring_function=scoring,
        num_poses=n_poses,
        best_affinity=best_aff,
        pose_files=[out_path] if status_str == "ok" else [],
        error_message="" if status_str == "ok" else (output or "Unknown error")[:300],
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

    dock_dir = base_dir / "docking"
    dock_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    receptors = _collect_receptors(base_dir, prepared_manifest, protein_workflow_data)
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

        if r.status == "failed":
            with log_lock:
                record_failure(effective_output_dir, combo_name, r)
                error_log[combo_name] = {
                    "error": r.error_message,
                    "timestamp": datetime.now().isoformat(),
                    "protein": r.protein_name,
                    "ligand": r.ligand_name,
                    "elapsed_time": round(r.elapsed_time, 2),
                }

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
                overwrite=overwrite,
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

    results_df = pd.DataFrame(raw_results)
    return results_df, docking_results


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
    args = parser.parse_args()

    # ── Load config ──
    config_path = Path(args.config).resolve()
    if not config_path.exists():
        print(f"ERROR: Config file not found: {config_path}")
        sys.exit(1)
    cfg = load_config(config_path)

    # ── Resolve paths relative to CWD ──
    receptors_dir = Path(cfg["receptors_dir"])
    output_base = Path(cfg["output_dir"])
    log_dir = Path(cfg.get("log_dir", "logs"))
    ligand_dirs = [Path(d) for d in cfg["ligand_dirs"]]

    output_base.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    protein_pdbqt_dir = get_pdbqt_dir(receptors_dir)

    # ── Validate directories ──
    for d in ligand_dirs:
        if not d.exists():
            print(f"ERROR: Ligand directory missing: {d}")
            sys.exit(1)
    if not receptors_dir.exists():
        print(f"ERROR: Receptor directory missing: {receptors_dir}")
        sys.exit(1)

    print_config(cfg)

    cpu_model = get_cpu_model()
    print(f"CPU: {cpu_model}")

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

            summary = generate_summary(docking_results, cfg)
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
