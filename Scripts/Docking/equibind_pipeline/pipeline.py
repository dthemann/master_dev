"""3-phase EquiBind docking pipeline orchestration.

Phase 1 (CPU, parallel): prep proteins, ligands, conformers, build DockJobs.
Phase 2 (GPU, batched):  EquiBind forward pass.
Phase 3 (CPU, parallel): corrections, pocket enforcement, UFF, write SDFs.
"""
from __future__ import annotations

import itertools
import json
import math
import multiprocessing as _mp
import shutil
import subprocess
import threading
import time
from collections import defaultdict
from concurrent.futures import (
    ProcessPoolExecutor, ThreadPoolExecutor, as_completed,
)
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Geometry import Point3D

# Suppress RDKit deprecation warnings (e.g. GetValence)
RDLogger.DisableLog("rdApp.*")

from .config import CFG
from .conformers import ensure_ligand_sdf, generate_conformer_sdfs
from .cropping import (
    clamp_pose_to_pocket, get_cropped_protein_for_pocket,
    translate_pose_back_to_original_frame,
)
from .equibind_model import (
    get_args, get_model, get_receptor_graph_cached, receptor_graph_count,
    run_corrections_inproc, run_gpu_batch,
    get_lig_graph_revised, get_geometry_graph,
)
from .logging_utils import (
    DOCKING_LOG_COLUMNS, append_log_row, init_docking_log, record_failure,
)
from .monitor import monitor
from .pockets import PocketInfo, parse_fpocket_results, parse_p2rank_results
from .properties import (
    get_ligand_properties, get_protein_properties, precompute_properties,
)
from .types import DockJob, GuidedPoseResult, P3Intermediate, TimingRecord
from .refine import refine_pose
from .uff import prewarm_protein_cache, uff_minimize_pose

_BAR = "\u2500" * 80


# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# Pose-variant fan-out  (centroid-clamp axis  x  re-search axis)
# \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
# Each EquiBind pose can fan out into several finalized variants so the effect
# of centroid clamping and of the docking re-search can be measured directly.
# When both axes are single-valued (clamp_mode='on', refine_mode='off') exactly
# one file is written with the legacy name \u2014 fully back-compatible.

def _clamp_variants() -> List[Tuple[str, str]]:
    """(filename_suffix, behaviour) per clamp variant for ``CFG.clamp_mode``.

    behaviour \u2208 {'default' (legacy force_pocket + clamp_pose_to_pocket; may
    reject far poses), 'force_clamp' (always clamp inward, never reject),
    'no_clamp' (keep the raw EquiBind centroid)}.
    """
    m = CFG.clamp_mode
    if m == "off":
        return [("", "no_clamp")]
    if m == "both":
        return [("__clampON", "force_clamp"), ("__clampOFF", "no_clamp")]
    return [("", "default")]  # "on"


def _refine_variants() -> List[Tuple[str, str]]:
    """(filename_suffix, kind) per re-search variant for ``CFG.refine_mode``.

    ``kind`` ∈ {'raw', 'smina', 'gnina'}. With one active tool the refined file
    keeps the legacy unsuffixed name ('on') / '__refRAW'+'__refSMINA' ('both');
    with ``refine_tool='both'`` each tool gets its own '__ref<TOOL>' suffix so a
    pose fans out into one refined variant per backend.
    """
    m = CFG.refine_mode
    if m == "off":
        return [("", "raw")]
    tools = CFG.refine_tool_list()

    def _tool_suffix(tool: str) -> str:
        return f"__ref{tool.upper()}" if len(tools) > 1 else ""

    if m == "on":
        return [(_tool_suffix(t), t) for t in tools]
    # "both": the un-refined pose plus one refined pose per active tool.
    return [("__refRAW", "raw")] + [(f"__ref{t.upper()}", t) for t in tools]


def _clamp_label(suffix: str) -> Optional[str]:
    return suffix.lstrip("_") or None


def _variant_specs() -> List[Tuple[str, Optional[str], Optional[str]]]:
    """(filename_suffix, clamp_variant_label, refine_kind_label) per variant."""
    refine_active = CFG.refine_mode != "off"
    specs: List[Tuple[str, Optional[str], Optional[str]]] = []
    for cs, _ in _clamp_variants():
        for rs, rk in _refine_variants():
            specs.append((cs + rs, _clamp_label(cs), rk if refine_active else None))
    return specs


def _variant_count() -> int:
    return len(_variant_specs())


def _variant_final_paths(base_sdf: Path) -> List[Path]:
    """Every output SDF a single pose fans out into (clamp \u00d7 re-search)."""
    d, base = base_sdf.parent, base_sdf.stem
    return [d / f"{base}{suf}.sdf" for suf, _, _ in _variant_specs()]


# ════════════════════════════════════════════════════════════════════════
# Small helpers
# ════════════════════════════════════════════════════════════════════════

def _sdf_centroid(sdf_path: Path) -> Optional[Tuple[float, float, float]]:
    """Read an SDF coord block directly without RDKit parse."""
    try:
        with open(sdf_path) as f:
            lines = f.readlines()
        if len(lines) < 5:
            return None
        n_atoms = int(lines[3].split()[0])
        coords = []
        for i in range(4, min(4 + n_atoms, len(lines))):
            p = lines[i].split()
            if len(p) >= 3:
                coords.append((float(p[0]), float(p[1]), float(p[2])))
        if not coords:
            return None
        arr = np.asarray(coords)
        c = arr.mean(axis=0)
        return float(c[0]), float(c[1]), float(c[2])
    except Exception:
        return None


def _translate_sdf_to_pocket(sdf_path: Path,
                              pocket_center: Tuple[float, float, float],
                              output_path: Path) -> bool:
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            mol = Chem.MolFromMol2File(str(sdf_path), removeHs=False)
        if mol is None:
            return False
        mol = Chem.RWMol(mol)
        conf = mol.GetConformer()
        positions = conf.GetPositions()
        centroid = positions.mean(axis=0)
        sx, sy, sz = (float(pocket_center[i] - centroid[i]) for i in range(3))
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x + sx, p.y + sy, p.z + sz))
        w = Chem.SDWriter(str(output_path)); w.write(mol); w.close()
        return output_path.exists() and output_path.stat().st_size > 0
    except Exception as e:
        monitor.warning(f"Could not translate ligand to pocket: {e}")
        return False


def _prepare_protein(protein_pdb: Path, prep_dir: Path) -> Path:
    """Run ``reduce`` (if available) on protein; otherwise copy."""
    prepared = prep_dir / f"{protein_pdb.stem}_protein.pdb"
    if prepared.exists():
        return prepared
    if CFG.reduce_executable:
        try:
            res = subprocess.run(
                [CFG.reduce_executable, "-Quiet", str(protein_pdb)],
                capture_output=True, text=True, timeout=60,
            )
            if res.returncode == 0 and res.stdout.strip():
                prepared.write_text(res.stdout)
                return prepared
        except Exception:
            pass
    shutil.copy2(protein_pdb, prepared)
    return prepared


def _collect_input_files(directory: Path, extensions: List[str],
                         name_filter: str = "") -> List[Path]:
    files: List[Path] = []
    for ext in extensions:
        files.extend(
            p for p in directory.glob(f"*{ext}")
            if p.is_file() and (not name_filter or name_filter in p.stem)
        )
    return sorted(set(files))


# ════════════════════════════════════════════════════════════════════════
# DockJob construction
# ════════════════════════════════════════════════════════════════════════

def _prepare_dock_job(
    conf_path: Path, seed: int, pocket: Optional[PocketInfo],
    prepared_protein: Path, protein_key: str,
    combo_dir: Path, prep_dir: Path,
    protein_name: str, ligand_name: str, lig_sdf: Path,
    mode: str, pose_num: int, final_sdf: Path,
    cropped_protein: Optional[Path], crop_offset: Optional[np.ndarray],
) -> Optional[DockJob]:
    t0 = time.time()
    combo_name = f"{ligand_name}__{protein_name}"

    if cropped_protein is not None:
        dock_protein = cropped_protein
        dock_key = str(cropped_protein)
        target_center: Optional[Tuple[float, float, float]] = (0.0, 0.0, 0.0)
    else:
        dock_protein = prepared_protein
        dock_key = protein_key
        target_center = pocket.center if pocket else None

    job_id = (f"{ligand_name}__{protein_name}__{mode}_"
              f"{pocket.unique_id if pocket else 'unguided'}_p{pose_num:02d}")
    translated = prep_dir / f"lig_{job_id}.sdf"

    if target_center is not None:
        ok = _translate_sdf_to_pocket(conf_path, target_center, translated)
        if not ok:
            ok = _translate_sdf_to_pocket(lig_sdf, target_center, translated)
        if not ok:
            return None
        lig_file = translated
    else:
        lig_file = conf_path

    try:
        suppl = Chem.SDMolSupplier(str(lig_file), sanitize=False, removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            return None
        Chem.SanitizeMol(mol)
        if not mol.HasProp("_Name"):
            mol.SetProp("_Name", ligand_name)
    except Exception:
        return None

    args = get_args()
    dp = args.dataset_params
    try:
        lig_graph = get_lig_graph_revised(
            mol, ligand_name,
            max_neighbors=dp["lig_max_neighbors"],
            use_rdkit_coords=False,
            radius=dp["lig_graph_radius"],
        )
        lig_graph.ndata["new_x"] = lig_graph.ndata["x"]
    except Exception:
        return None

    try:
        geometry_graph = get_geometry_graph(mol) if dp.get("geometry_regularization", True) else None
    except Exception:
        geometry_graph = None
    lig_coord = lig_graph.ndata["new_x"].clone()

    return DockJob(
        job_id=job_id, combo_name=combo_name,
        protein_pdb=dock_protein, protein_key=dock_key,
        prepared_protein=prepared_protein,
        ligand_sdf=lig_file, ligand_mol=mol,
        lig_graph=lig_graph, lig_coord=lig_coord, geometry_graph=geometry_graph,
        seed=seed, pocket=pocket,
        pocket_center=pocket.center if pocket else None,
        crop_offset=crop_offset,
        combo_dir=combo_dir, final_sdf=final_sdf, pose_num=pose_num,
        protein_name=protein_name, ligand_name=ligand_name,
        mode=pocket.source if pocket else "unguided",
        prep_time=time.time() - t0,
    )


# ════════════════════════════════════════════════════════════════════════
# Phase 3 helpers — Stage 1 (pre-UFF) + Stage 3 (finalize)
# ════════════════════════════════════════════════════════════════════════

def _result_failed(job: DockJob, error: str, post_time: float = 0.0,
                   clamp_variant: Optional[str] = None) -> GuidedPoseResult:
    return GuidedPoseResult(
        protein_name=job.protein_name, ligand_name=job.ligand_name,
        mode=job.mode,
        pocket_id=job.pocket.pocket_id if job.pocket else None,
        pocket_unique_id=job.pocket.unique_id if job.pocket else None,
        pose_num=job.pose_num,
        pocket_center=job.pocket_center, pose_centroid=None,
        sdf_path=None, success=False, error=error, clamp_variant=clamp_variant,
        prep_time_s=job.prep_time, dock_time_s=job.dock_time, post_time_s=post_time,
    )


def _result_failed_inter(inter: P3Intermediate, error: str,
                         post_time: float = 0.0) -> GuidedPoseResult:
    return GuidedPoseResult(
        protein_name=inter.protein_name, ligand_name=inter.ligand_name,
        mode=inter.mode, pocket_id=inter.pocket_id,
        pocket_unique_id=inter.pocket_unique_id, pose_num=inter.pose_num,
        pocket_center=inter.pocket_center, pose_centroid=None,
        sdf_path=None, success=False, error=error, clamp_variant=inter.clamp_variant,
        prep_time_s=inter.prep_time, dock_time_s=inter.dock_time, post_time_s=post_time,
    )


def _skipped_variant_result(job: DockJob, path: Path,
                            clamp_variant: Optional[str],
                            refine_variant: Optional[str]) -> GuidedPoseResult:
    return GuidedPoseResult(
        protein_name=job.protein_name, ligand_name=job.ligand_name,
        mode=job.mode,
        pocket_id=job.pocket.pocket_id if job.pocket else None,
        pocket_unique_id=job.pocket.unique_id if job.pocket else None,
        pose_num=job.pose_num,
        pocket_center=job.pocket_center, pose_centroid=_sdf_centroid(path),
        sdf_path=path, success=True,
        clamp_variant=clamp_variant, refine_variant=refine_variant,
    )


def _tag_pose_source(sdf_path: Path, mode: Optional[str],
                     pocket_unique_id: Optional[str] = None,
                     clamp_variant: Optional[str] = None,
                     refine_variant: Optional[str] = None,
                     refine_affinity: Optional[float] = None,
                     cnn_score: Optional[float] = None,
                     cnn_affinity: Optional[float] = None) -> None:
    """Embed pose provenance as SDF properties on a finalized pose.

    Records ``<pocket_source>`` (fpocket | p2rank | unguided), ``<pocket_id>``
    (guided), and — when the corresponding axis is active — ``<clamp_variant>``,
    ``<refine_variant>``, the per-tool affinity (``<smina_affinity>`` /
    ``<gnina_affinity>``) and, for gnina, ``<cnn_score>`` / ``<cnn_affinity>``.
    Done at the text level so the molecule block stays byte-for-byte intact;
    best-effort — a tagging failure must never lose the pose.
    """
    try:
        text = sdf_path.read_text()
        block = f">  <pocket_source>\n{mode or 'unguided'}\n\n"
        if pocket_unique_id:
            block += f">  <pocket_id>\n{pocket_unique_id}\n\n"
        if clamp_variant:
            block += f">  <clamp_variant>\n{clamp_variant}\n\n"
        if refine_variant:
            block += f">  <refine_variant>\n{refine_variant}\n\n"
        if refine_affinity is not None:
            # Per-tool tag: <smina_affinity> stays byte-identical for smina poses.
            aff_tag = (f"{refine_variant}_affinity"
                       if refine_variant in ("smina", "gnina") else "refine_affinity")
            block += f">  <{aff_tag}>\n{refine_affinity:.3f}\n\n"
        if cnn_score is not None:
            block += f">  <cnn_score>\n{cnn_score:.3f}\n\n"
        if cnn_affinity is not None:
            block += f">  <cnn_affinity>\n{cnn_affinity:.3f}\n\n"
        if "$$$$" in text:
            i = text.rindex("$$$$")
            text = text[:i] + block + text[i:]
        else:
            text = text.rstrip("\n") + "\n" + block + "$$$$\n"
        sdf_path.write_text(text)
    except Exception:
        pass


def _postprocess_stage1(job: DockJob, clamp_suffix: str,
                        clamp_behavior: str) -> Tuple[Optional[GuidedPoseResult],
                                                      Optional[P3Intermediate]]:
    """Stage 1 for one clamp variant: corrections + SDF write + centroid clamp.

    Returns ``(failed_result | None, P3Intermediate | None)``. Never finalizes \u2014
    Stage 3 (:func:`_emit_variants`) writes the clamp \u00d7 re-search output files.
    """
    t0 = time.time()
    clamp_variant = _clamp_label(clamp_suffix)

    if not job.gpu_success or job.predicted_coords is None:
        return _result_failed(job, job.error, clamp_variant=clamp_variant), None

    try:
        optimized_mol = run_corrections_inproc(
            job.ligand_mol, job.lig_coord, job.predicted_coords)
    except Exception:
        optimized_mol = deepcopy(job.ligand_mol)
        conf = optimized_mol.GetConformer()
        coords_np = job.predicted_coords.detach().cpu().numpy()
        for i in range(optimized_mol.GetNumAtoms()):
            conf.SetAtomPosition(i, Point3D(*coords_np[i].tolist()))

    pose_dir = job.combo_dir / f"{job.job_id}{clamp_suffix}_run"
    pose_dir.mkdir(parents=True, exist_ok=True)
    out_sdf = pose_dir / "output.sdf"
    try:
        w = Chem.SDWriter(str(out_sdf)); w.write(optimized_mol); w.close()
    except Exception as e:
        return _result_failed(job, f"SDF write failed: {e}",
                              post_time=time.time() - t0,
                              clamp_variant=clamp_variant), None

    if job.crop_offset is not None and np.any(job.crop_offset != 0):
        restored = pose_dir / "output_restored.sdf"
        if translate_pose_back_to_original_frame(out_sdf, job.crop_offset, restored):
            out_sdf = restored

    centroid = _sdf_centroid(out_sdf)

    # \u2500\u2500 Centroid clamping \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
    # 'no_clamp'    \u2192 keep the raw EquiBind centroid (no translate, no reject).
    # 'force_clamp' \u2192 always pull a far centroid back inside the pocket.
    # 'default'     \u2192 legacy: honour CFG.force_pocket + CFG.clamp_pose_to_pocket
    #                 (clamp far poses, or reject them if clamping is disabled).
    if clamp_behavior != "no_clamp" and centroid is not None and job.pocket_center is not None:
        dist = math.dist(centroid, job.pocket_center)
        if dist > CFG.pocket_match_threshold:
            do_clamp = (clamp_behavior == "force_clamp"
                        or (CFG.force_pocket and CFG.clamp_pose_to_pocket))
            do_reject = (clamp_behavior == "default"
                         and CFG.force_pocket and not CFG.clamp_pose_to_pocket)
            if do_clamp:
                clamped = pose_dir / "output_clamped.sdf"
                ok, _, _ = clamp_pose_to_pocket(
                    sdf_path=out_sdf, pocket_center=job.pocket_center,
                    max_distance=CFG.pocket_match_threshold, output_path=clamped,
                )
                if ok:
                    out_sdf = clamped
                elif clamp_behavior == "default":
                    shutil.rmtree(pose_dir, ignore_errors=True)
                    return _result_failed(
                        job,
                        f"Pose rejected: {dist:.1f}\u00c5 from pocket (clamping failed)",
                        post_time=time.time() - t0,
                        clamp_variant=clamp_variant), None
                # force_clamp + clamp failure \u2192 keep the unclamped pose.
            elif do_reject:
                shutil.rmtree(pose_dir, ignore_errors=True)
                return _result_failed(
                    job, f"Pose rejected: {dist:.1f}\u00c5 from pocket",
                    post_time=time.time() - t0,
                    clamp_variant=clamp_variant), None

    return None, P3Intermediate(
        out_sdf=out_sdf, pose_dir=pose_dir,
        minimized_sdf=pose_dir / "output_uff_minimized.sdf",
        final_sdf=job.final_sdf,
        protein_pdb=job.prepared_protein,
        protein_name=job.protein_name, ligand_name=job.ligand_name,
        mode=job.mode,
        pocket_id=job.pocket.pocket_id if job.pocket else None,
        pocket_unique_id=job.pocket.unique_id if job.pocket else None,
        pose_num=job.pose_num,
        pocket_center=job.pocket_center,
        prep_time=job.prep_time, dock_time=job.dock_time,
        pre_uff_time=time.time() - t0,
        clamp_suffix=clamp_suffix, clamp_variant=clamp_variant,
    )


def _uff_minimize_timed(docked_sdf: Path, protein_pdb: Path,
                         output_sdf: Path) -> tuple:
    """Pickle-safe wrapper used by the ProcessPool stage."""
    t0 = time.time()
    ok, e_before, e_after, msg = uff_minimize_pose(
        docked_sdf=docked_sdf, protein_pdb=protein_pdb, output_sdf=output_sdf)
    return ok, e_before, e_after, msg, time.time() - t0


def _emit_variants(inter: P3Intermediate, geom_sdf: Path, post_time: float,
                   uff_minimized: bool, uff_e_before: Optional[float],
                   uff_e_after: Optional[float]) -> List[GuidedPoseResult]:
    """Stage 3: write the re-search variants for one (already-clamped) pose.

    ``geom_sdf`` is the finished geometry (UFF-minimized if available, else the
    pre-UFF pose). For each re-search variant: 'raw' just copies it; 'smina' /
    'gnina' run a local docking re-search and fall back to the raw pose on
    failure so a pose is never lost. Cleans up the per-pose scratch dir when done.
    """
    base_sdf = inter.final_sdf
    d, base = base_sdf.parent, base_sdf.stem
    d.mkdir(parents=True, exist_ok=True)
    refine_active = CFG.refine_mode != "off"
    results: List[GuidedPoseResult] = []

    for rsuf, rkind in _refine_variants():
        target = d / f"{base}{inter.clamp_suffix}{rsuf}.sdf"
        affinity: Optional[float] = None
        cnn_score: Optional[float] = None
        cnn_affinity: Optional[float] = None
        err = ""
        # Time the re-search per refine variant so 'raw' (a plain copy, ~0 s) and
        # 'smina'/'gnina' (a local docking subprocess) are directly comparable.
        # post_time (corrections + clamp + UFF) is the shared geometry work and is
        # kept separate; refine_time_s isolates the re-search cost on top of it.
        r0 = time.time()
        if rkind in ("smina", "gnina"):
            ok, scores, msg = refine_pose(
                geom_sdf, inter.protein_pdb, target, tool=rkind)
            if not ok or not (target.exists() and target.stat().st_size > 0):
                monitor.warning(
                    f"{rkind} refine failed for {target.name}: {msg}; keeping pre-refine pose")
                shutil.copy2(geom_sdf, target)
                err = f"refine_failed: {msg}"
            else:
                affinity = scores.get("minimized_affinity")
                cnn_score = scores.get("cnn_score")
                cnn_affinity = scores.get("cnn_affinity")
        else:
            shutil.copy2(geom_sdf, target)
        refine_time = time.time() - r0

        refine_variant = rkind if refine_active else None
        _tag_pose_source(target, inter.mode, inter.pocket_unique_id,
                         clamp_variant=inter.clamp_variant,
                         refine_variant=refine_variant, refine_affinity=affinity,
                         cnn_score=cnn_score, cnn_affinity=cnn_affinity)
        results.append(GuidedPoseResult(
            protein_name=inter.protein_name, ligand_name=inter.ligand_name,
            mode=inter.mode, pocket_id=inter.pocket_id,
            pocket_unique_id=inter.pocket_unique_id, pose_num=inter.pose_num,
            pocket_center=inter.pocket_center, pose_centroid=_sdf_centroid(target),
            sdf_path=target, success=True, error=err,
            uff_minimized=uff_minimized,
            uff_energy_before=uff_e_before, uff_energy_after=uff_e_after,
            clamp_variant=inter.clamp_variant, refine_variant=refine_variant,
            refine_affinity=affinity, cnn_score=cnn_score, cnn_affinity=cnn_affinity,
            prep_time_s=inter.prep_time, dock_time_s=inter.dock_time,
            post_time_s=post_time, refine_time_s=refine_time,
        ))

    shutil.rmtree(inter.pose_dir, ignore_errors=True)
    return results


# ════════════════════════════════════════════════════════════════════════
# Pipeline entry point
# ════════════════════════════════════════════════════════════════════════

def run_pipeline() -> dict:  # noqa: C901  (long but linear)
    """Run the complete 3-phase EquiBind pipeline.

    Returns a summary dict (also written to ``pipeline_summary.json``).
    """
    CFG.validate()
    CFG.ensure_output_dirs()
    CFG.banner()

    out_dir = CFG.pocket_guided_output_dir

    ligand_files = _collect_input_files(CFG.drugs_dir, CFG.ligand_extensions)
    receptor_files = _collect_input_files(CFG.receptors_dir, [".pdb"],
                                          name_filter=CFG.receptor_name_filter)
    if not ligand_files or not receptor_files:
        raise RuntimeError(
            f"No input files found (ligands={len(ligand_files)}, "
            f"receptors={len(receptor_files)})"
        )

    n_per_combo = (CFG.n_top_pockets * CFG.poses_per_pocket * 2) + CFG.n_unguided_poses
    total_combos = len(receptor_files) * len(ligand_files)

    # Trigger model load (so receptor graph builds work)
    model, args, device = get_model()

    monitor.header("GLOBALLY DECOUPLED 3-PHASE EQUIBIND PIPELINE")
    print(f"  EquiBind device: {device}")
    print(f"  GPU batch size: {CFG.gpu_batch_size}")
    print(f"  Proteins: {len(receptor_files)}")
    print(f"  Ligands: {len(ligand_files)}")
    print(f"  Total combos: {total_combos}  (\u2248 {n_per_combo} poses/combo)")
    print(f"  Workers: {CFG.n_parallel_workers}")
    print()

    # ── Pocket discovery ─────────────────────────────────────────────────
    protein_fpocket: Dict[str, List[PocketInfo]] = {}
    protein_p2rank: Dict[str, List[PocketInfo]] = {}
    for prot in receptor_files:
        pname = prot.stem
        fp = parse_fpocket_results(CFG.fpocket_results_folder, pname, CFG.n_top_pockets)
        pr = parse_p2rank_results(CFG.p2rank_folder, pname, CFG.n_top_pockets)
        protein_fpocket[pname] = fp
        protein_p2rank[pname] = pr
        monitor.info(f"{pname}: {len(fp)} fpocket + {len(pr)} p2rank pockets")
    print()

    lig_props_cache, prot_props_cache = precompute_properties(receptor_files, ligand_files)
    log_path, existing_combos = init_docking_log(out_dir)
    print(f"Docking log: {log_path} ({len(existing_combos)} existing entries)")

    pipeline_start = time.time()

    # ════════════════════════════════════════════════════════════════════
    # PHASE 1 — CPU prep
    # ════════════════════════════════════════════════════════════════════
    phase1_start = time.time()
    monitor.header("PHASE 1: CPU PREPARATION \u2014 ALL COMBOS")

    all_jobs: List[DockJob] = []
    cached_results: List[GuidedPoseResult] = []
    combo_phase1_times: Dict[str, float] = {}

    # Prepared proteins + receptor graphs per unique protein
    prepared_proteins: Dict[str, Path] = {}
    protein_keys: Dict[str, str] = {}
    for protein in receptor_files:
        pname = protein.stem
        prep_dir = out_dir / f"_prep_{pname}"
        prep_dir.mkdir(parents=True, exist_ok=True)
        prepared = _prepare_protein(protein, prep_dir)
        prepared_proteins[pname] = prepared
        protein_keys[pname] = str(prepared)
        monitor.info(f"Building receptor graph for {prepared.name}")
        get_receptor_graph_cached(prepared)

    # Cropped proteins per pocket (parallel)
    all_pockets_flat: List[Tuple[str, PocketInfo]] = []
    for prot in receptor_files:
        pname = prot.stem
        for p in protein_fpocket.get(pname, []) + protein_p2rank.get(pname, []):
            all_pockets_flat.append((pname, p))

    cropped_map: Dict[str, Dict[str, Tuple[Optional[Path], np.ndarray]]] = defaultdict(dict)

    def _prepare_cropped_for_pocket(pn: str, pk: PocketInfo):
        if not CFG.use_protein_cropping:
            return pn, pk.unique_id, None, np.zeros(3)
        prep_dir = out_dir / f"_prep_{pn}"
        cropped, offset = get_cropped_protein_for_pocket(
            protein_pdb=prepared_proteins[pn], pocket=pk,
            prep_dir=prep_dir, crop_radius=CFG.pocket_crop_radius,
        )
        if cropped is not None:
            get_receptor_graph_cached(cropped)
        return pn, pk.unique_id, cropped, offset

    if all_pockets_flat:
        n_workers = min(CFG.n_parallel_workers, len(all_pockets_flat))
        monitor.info(f"Preparing cropped proteins for {len(all_pockets_flat)} pockets "
                     f"({n_workers} workers)")
        with ThreadPoolExecutor(max_workers=n_workers) as ex:
            futures = {ex.submit(_prepare_cropped_for_pocket, pn, pk): (pn, pk)
                       for pn, pk in all_pockets_flat}
            for fut in as_completed(futures):
                try:
                    pname, uid, cropped, offset = fut.result()
                    cropped_map[pname][uid] = (cropped, offset)
                except Exception as e:
                    pn, pk = futures[fut]
                    monitor.warning(
                        f"Cropped protein prep failed for {pn}/{pk.unique_id}: {e}")
                    cropped_map[pn][pk.unique_id] = (None, np.zeros(3))

    # Conformers per ligand (parallel)
    ligand_conformers: Dict[str, List[Tuple[Path, int]]] = {}
    ligand_sdfs: Dict[str, Path] = {}

    def _prepare_ligand_conformers(ligand: Path):
        lname = ligand.stem
        lig_prep_dir = out_dir / f"_prep_lig_{lname}"
        lig_prep_dir.mkdir(parents=True, exist_ok=True)
        lig_sdf = ensure_ligand_sdf(ligand, lig_prep_dir)

        conf_dir = lig_prep_dir / "conformers"
        existing = sorted(conf_dir.glob("seed*_c*.sdf")) if conf_dir.exists() else []
        if existing:
            files: List[Tuple[Path, int]] = []
            for cp in existing:
                if cp.stat().st_size == 0:
                    continue
                try:
                    s = int(cp.stem.split("_c")[0].replace("seed", ""))
                except (ValueError, IndexError):
                    s = 42
                files.append((cp, s))
            if files:
                monitor.info(f"Ligand {lname}: {len(files)} existing conformers reused")
                return lname, lig_sdf, files
        files = generate_conformer_sdfs(ligand, conf_dir)
        if not files:
            files = [(lig_sdf, 42)]
        monitor.info(f"Ligand {lname}: {len(files)} conformers")
        return lname, lig_sdf, files

    monitor.info(f"Generating conformers for {len(ligand_files)} ligands "
                 f"({CFG.n_parallel_workers} workers)")
    with ThreadPoolExecutor(max_workers=CFG.n_parallel_workers) as ex:
        futures = {ex.submit(_prepare_ligand_conformers, lig): lig for lig in ligand_files}
        for i, fut in enumerate(as_completed(futures)):
            try:
                lname, lig_sdf, files = fut.result()
                ligand_sdfs[lname] = lig_sdf
                ligand_conformers[lname] = files
            except Exception as e:
                monitor.warning(f"Conformer gen failed for {futures[fut].stem}: {e}")
            if (i + 1) % 50 == 0:
                monitor.info(f"[HEARTBEAT] Conformer gen: {i+1}/{len(ligand_files)} done")

    # Build DockJobs per combo (parallel)
    all_combos = list(itertools.product(receptor_files, ligand_files))

    def _process_combo(protein: Path, ligand: Path):
        pname, lname = protein.stem, ligand.stem
        combo_name = f"{lname}__{pname}"
        t0 = time.time()
        combo_dir = out_dir / combo_name
        combo_dir.mkdir(parents=True, exist_ok=True)
        prep_dir = combo_dir / "prep"
        prep_dir.mkdir(parents=True, exist_ok=True)

        prepared_protein = prepared_proteins[pname]
        pkey = protein_keys[pname]
        lig_sdf = ligand_sdfs[lname]
        confs = ligand_conformers[lname]
        fp_pockets = protein_fpocket.get(pname, [])
        pr_pockets = protein_p2rank.get(pname, [])
        all_pockets = fp_pockets + pr_pockets

        local_jobs: List[DockJob] = []
        local_cached: List[GuidedPoseResult] = []

        # Combo-level cache
        results_file = combo_dir / "guided_results.json"
        if CFG.skip_existing and results_file.exists():
            try:
                cached = json.loads(results_file.read_text())
                expected = ((len(all_pockets) * CFG.poses_per_pocket
                             + CFG.n_unguided_poses) * _variant_count())
                if len(cached) >= expected:
                    for r in cached:
                        local_cached.append(GuidedPoseResult(
                            protein_name=r["protein_name"], ligand_name=r["ligand_name"],
                            mode=r["mode"], pocket_id=r.get("pocket_id"),
                            pocket_unique_id=r.get("pocket_unique_id"),
                            pose_num=r.get("pose_num", 1),
                            pocket_center=tuple(r["pocket_center"]) if r.get("pocket_center") else None,
                            pose_centroid=tuple(r["pose_centroid"]) if r.get("pose_centroid") else None,
                            sdf_path=Path(r["sdf_path"]) if r.get("sdf_path") else None,
                            success=r["success"], error=r.get("error", ""),
                        ))
                    monitor.info(f"[SKIP] {combo_name}: {len(cached)} cached results")
                    return combo_name, local_jobs, local_cached, time.time() - t0
            except Exception:
                pass

        for pocket in all_pockets:
            uid = pocket.unique_id
            cropped, offset = cropped_map.get(pname, {}).get(uid, (None, np.zeros(3)))
            conf_idx = 0
            for pose_num in range(1, CFG.poses_per_pocket + 1):
                final_sdf = combo_dir / f"{uid}_pose{pose_num:02d}.sdf"
                vpaths = _variant_final_paths(final_sdf)
                if CFG.skip_existing and all(vp.exists() for vp in vpaths):
                    for vp, (_, clab, rk) in zip(vpaths, _variant_specs()):
                        local_cached.append(GuidedPoseResult(
                            protein_name=pname, ligand_name=lname,
                            mode=pocket.source, pocket_id=pocket.pocket_id,
                            pocket_unique_id=pocket.unique_id, pose_num=pose_num,
                            pocket_center=pocket.center,
                            pose_centroid=_sdf_centroid(vp),
                            sdf_path=vp, success=True,
                            clamp_variant=clab, refine_variant=rk,
                        ))
                    continue
                if conf_idx >= len(confs):
                    break
                conf_path, seed = confs[conf_idx]
                conf_idx += 1
                job = _prepare_dock_job(
                    conf_path=conf_path, seed=seed, pocket=pocket,
                    prepared_protein=prepared_protein, protein_key=pkey,
                    combo_dir=combo_dir, prep_dir=prep_dir,
                    protein_name=pname, ligand_name=lname,
                    lig_sdf=lig_sdf, mode=pocket.source, pose_num=pose_num,
                    final_sdf=final_sdf, cropped_protein=cropped, crop_offset=offset,
                )
                if job is not None:
                    local_jobs.append(job)

        unguided_count = 0
        for conf_path, seed in confs:
            if unguided_count >= CFG.n_unguided_poses:
                break
            pose_num = unguided_count + 1
            final_sdf = combo_dir / f"unguided_{pose_num:03d}.sdf"
            vpaths = _variant_final_paths(final_sdf)
            if CFG.skip_existing and all(vp.exists() for vp in vpaths):
                for vp, (_, clab, rk) in zip(vpaths, _variant_specs()):
                    local_cached.append(GuidedPoseResult(
                        protein_name=pname, ligand_name=lname,
                        mode="unguided", pocket_id=None, pocket_unique_id=None,
                        pose_num=pose_num, pocket_center=None,
                        pose_centroid=_sdf_centroid(vp),
                        sdf_path=vp, success=True,
                        clamp_variant=clab, refine_variant=rk,
                    ))
                unguided_count += 1
                continue
            job = _prepare_dock_job(
                conf_path=conf_path, seed=seed, pocket=None,
                prepared_protein=prepared_protein, protein_key=pkey,
                combo_dir=combo_dir, prep_dir=prep_dir,
                protein_name=pname, ligand_name=lname,
                lig_sdf=lig_sdf, mode="unguided", pose_num=pose_num,
                final_sdf=final_sdf, cropped_protein=None, crop_offset=None,
            )
            if job is not None:
                local_jobs.append(job)
                unguided_count += 1
        return combo_name, local_jobs, local_cached, time.time() - t0

    monitor.info(f"Building DockJobs for {len(all_combos)} combos "
                 f"({CFG.n_parallel_workers} workers)")
    with ThreadPoolExecutor(max_workers=CFG.n_parallel_workers) as ex:
        futures = {ex.submit(_process_combo, prot, lig): (prot, lig)
                   for prot, lig in all_combos}
        for i, fut in enumerate(as_completed(futures)):
            try:
                combo_name, jobs, cached, elapsed = fut.result()
                all_jobs.extend(jobs)
                cached_results.extend(cached)
                combo_phase1_times[combo_name] = elapsed
            except Exception as e:
                prot, lig = futures[fut]
                monitor.warning(f"Combo {lig.stem}__{prot.stem} failed: {e}")
            if (i + 1) % 100 == 0:
                monitor.info(f"[HEARTBEAT] DockJob building: {i+1}/{len(all_combos)} done")

    phase1_time = time.time() - phase1_start
    n_jobs = len(all_jobs)
    print(f"\n{_BAR}")
    print(f"  PHASE 1 COMPLETE in {phase1_time:.1f}s ({phase1_time/60:.1f} min)")
    print(f"  Jobs prepared:    {n_jobs}")
    print(f"  Cached results:   {len(cached_results)}")
    print(f"  Unique proteins:  {len({j.protein_key for j in all_jobs})}")
    print(f"  Unique ligands:   {len({j.ligand_name for j in all_jobs})}")
    print(f"  Receptor graphs:  {receptor_graph_count()}")
    print(f"{_BAR}")

    # ════════════════════════════════════════════════════════════════════
    # PHASE 2 — GPU inference
    # ════════════════════════════════════════════════════════════════════
    phase2_start = time.time()
    monitor.header(f"PHASE 2: GPU INFERENCE \u2014 {n_jobs} TOTAL JOBS")

    if n_jobs > 0:
        jobs_by_protein: Dict[str, List[DockJob]] = defaultdict(list)
        for job in all_jobs:
            jobs_by_protein[job.protein_key].append(job)

        total_batches = 0
        bs = CFG.gpu_batch_size
        for pkey, pjobs in jobs_by_protein.items():
            n_batches = (len(pjobs) + bs - 1) // bs
            total_batches += n_batches
            monitor.info(f"  {Path(pkey).name}: {len(pjobs)} jobs \u2192 {n_batches} batches")
            for i in range(0, len(pjobs), bs):
                run_gpu_batch(pjobs[i:i + bs])
        phase2_time = time.time() - phase2_start

        gpu_ok = sum(1 for j in all_jobs if j.gpu_success)
        gpu_fail = n_jobs - gpu_ok
        print(f"\n{_BAR}")
        print(f"  PHASE 2 COMPLETE in {phase2_time:.1f}s ({phase2_time/60:.1f} min)")
        print(f"  Total batches: {total_batches}")
        print(f"  Success: {gpu_ok}  |  Failed: {gpu_fail}")
        if gpu_ok > 0:
            avg = sum(j.dock_time for j in all_jobs if j.gpu_success) / gpu_ok
            print(f"  Avg inference: {avg*1000:.1f}ms per ligand")
        print(f"{_BAR}")
    else:
        phase2_time = 0.0
        print("\n  No jobs to dock (all cached). Phase 2 skipped.")

    # ════════════════════════════════════════════════════════════════════
    # PHASE 3 — CPU post-processing
    # ════════════════════════════════════════════════════════════════════
    phase3_start = time.time()
    monitor.header(f"PHASE 3: CPU POST-PROCESSING \u2014 {n_jobs} TOTAL JOBS")

    n_variants = _variant_count()
    if n_variants > 1:
        monitor.info(f"Pose-variant fan-out active: clamp_mode={CFG.clamp_mode}, "
                     f"refine_mode={CFG.refine_mode} \u2192 {n_variants} output(s)/pose")

    new_results: List[GuidedPoseResult] = []
    p3_done_skipped: List[GuidedPoseResult] = []
    p3_todo: List[DockJob] = []
    variant_specs = _variant_specs()
    for job in all_jobs:
        vpaths = _variant_final_paths(job.final_sdf)
        if all(p.exists() and p.stat().st_size > 0 for p in vpaths):
            for p, (_, clab, rk) in zip(vpaths, variant_specs):
                p3_done_skipped.append(_skipped_variant_result(job, p, clab, rk))
        else:
            p3_todo.append(job)
    if p3_done_skipped:
        monitor.info(f"Skipping {len(p3_done_skipped)} already-postprocessed pose(s)")
    new_results.extend(p3_done_skipped)
    n_todo = len(p3_todo)
    monitor.info(f"Post-processing {n_todo} new jobs ({len(p3_done_skipped)} skipped)")

    if n_todo > 0:
        # Each job fans out into one Stage-1 unit per clamp variant (clamping
        # happens before UFF, so each variant needs its own relaxation).
        units: List[Tuple[DockJob, str, str]] = [
            (job, csuf, cbehav)
            for job in p3_todo for (csuf, cbehav) in _clamp_variants()
        ]
        n_units = len(units)
        n_workers = min(CFG.n_parallel_workers, n_units)

        if CFG.uff_minimize:
            prewarm_protein_cache({j.prepared_protein for j in p3_todo})

        # Stage 1 \u2014 corrections + centroid clamp (ThreadPool)
        s1_start = time.time()
        s1_results: List[Optional[Tuple[Optional[GuidedPoseResult],
                                        Optional[P3Intermediate]]]] = [None] * n_units

        def _safe_stage1(idx: int, job: DockJob, csuf: str, cbehav: str):
            try:
                return idx, _postprocess_stage1(job, csuf, cbehav)
            except Exception as e:
                return idx, (_result_failed(job, f"Stage1 exception: {e}",
                                            clamp_variant=_clamp_label(csuf)), None)

        if n_workers > 1:
            with ThreadPoolExecutor(max_workers=n_workers) as ex:
                futures = [ex.submit(_safe_stage1, i, j, cs, cb)
                           for i, (j, cs, cb) in enumerate(units)]
                for fut in as_completed(futures):
                    idx, res = fut.result()
                    s1_results[idx] = res
        else:
            for i, (j, cs, cb) in enumerate(units):
                _, s1_results[i] = _safe_stage1(i, j, cs, cb)

        s1_time = time.time() - s1_start

        intermediates: List[P3Intermediate] = []
        for result, intermediate in s1_results:
            if result is not None:
                new_results.append(result)
            elif intermediate is not None:
                intermediates.append(intermediate)
        n_inter = len(intermediates)
        monitor.info(f"Stage 1 (corrections+clamp): {s1_time:.1f}s \u2014 "
                     f"{n_inter} pose(s) to finalize")

        # Stage 2 \u2014 UFF in ProcessPool (one relaxation per clamp variant)
        uff_results: List[Optional[tuple]] = [None] * n_inter
        if CFG.uff_minimize and intermediates:
            s2_start = time.time()
            n_uff_workers = min(n_workers, n_inter)
            monitor.info(f"Stage 2: UFF \u2014 {n_inter} poses, {n_uff_workers} process workers")
            mp_ctx = _mp.get_context("fork")
            with ProcessPoolExecutor(max_workers=n_uff_workers, mp_context=mp_ctx) as pool:
                futures = {
                    pool.submit(_uff_minimize_timed,
                                docked_sdf=inter.out_sdf,
                                protein_pdb=inter.protein_pdb,
                                output_sdf=inter.minimized_sdf): j
                    for j, inter in enumerate(intermediates)
                }
                for fut in as_completed(futures):
                    j = futures[fut]
                    try:
                        uff_results[j] = fut.result()
                    except Exception as e:
                        uff_results[j] = (False, 0.0, 0.0, f"UFF process error: {e}", 0.0)
            s2_time = time.time() - s2_start
            monitor.info(f"Stage 2 (UFF): {s2_time:.1f}s "
                         f"({n_inter / s2_time if s2_time > 0 else 0:.1f} poses/s)")

        # Stage 3 \u2014 finalize + re-search variants (ThreadPool; smina subprocesses)
        s3_start = time.time()

        def _safe_emit(idx: int) -> List[GuidedPoseResult]:
            inter = intermediates[idx]
            if CFG.uff_minimize:
                ok, e_before, e_after, _msg, uff_time = (
                    uff_results[idx] or (False, None, None, "", 0.0))
                uff_ok = bool(ok) and inter.minimized_sdf.exists()
                geom = inter.minimized_sdf if uff_ok else inter.out_sdf
                e_b = e_before if uff_ok else None
                e_a = e_after if uff_ok else None
                post = inter.pre_uff_time + uff_time
            else:
                geom, uff_ok, e_b, e_a, post = (
                    inter.out_sdf, False, None, None, inter.pre_uff_time)
            try:
                return _emit_variants(inter, geom, post, uff_ok, e_b, e_a)
            except Exception as e:
                return [_result_failed_inter(inter, f"emit failed: {e}", post)]

        if n_inter:
            if n_workers > 1:
                with ThreadPoolExecutor(max_workers=n_workers) as ex:
                    futures = [ex.submit(_safe_emit, i) for i in range(n_inter)]
                    for fut in as_completed(futures):
                        new_results.extend(fut.result())
            else:
                for i in range(n_inter):
                    new_results.extend(_safe_emit(i))
        s3_time = time.time() - s3_start
        if CFG.refine_mode != "off":
            monitor.info(f"Stage 3 (finalize+re-search): {s3_time:.1f}s")

        _save_checkpoint(out_dir, cached_results, new_results)

    phase3_time = time.time() - phase3_start
    n_new_ok = sum(1 for r in new_results if r.success)
    n_new_fail = sum(1 for r in new_results if not r.success)
    print(f"\n{_BAR}")
    print(f"  PHASE 3 COMPLETE in {phase3_time:.1f}s ({phase3_time/60:.1f} min)")
    print(f"  Already done: {len(p3_done_skipped)}")
    print(f"  Newly processed: {n_todo}")
    print(f"  Success: {n_new_ok}  |  Failed: {n_new_fail}")
    print(f"{_BAR}")

    # ════════════════════════════════════════════════════════════════════
    # Results, logging, summary
    # ════════════════════════════════════════════════════════════════════
    pipeline_elapsed = time.time() - pipeline_start
    all_results = cached_results + new_results

    summary = _finalize(
        out_dir=out_dir, log_path=log_path, existing_combos=existing_combos,
        receptor_files=receptor_files, ligand_files=ligand_files,
        all_jobs=all_jobs, new_results=new_results, cached_results=cached_results,
        all_results=all_results,
        protein_fpocket=protein_fpocket, protein_p2rank=protein_p2rank,
        lig_props_cache=lig_props_cache, prot_props_cache=prot_props_cache,
        combo_phase1_times=combo_phase1_times,
        phase1_time=phase1_time, phase2_time=phase2_time, phase3_time=phase3_time,
        pipeline_elapsed=pipeline_elapsed,
        n_cached=len(cached_results), n_new_ok=n_new_ok, n_new_fail=n_new_fail,
        n_jobs=n_jobs, total_combos=total_combos,
    )

    monitor.print_summary()
    return summary


# ════════════════════════════════════════════════════════════════════════
# Result writing & timing report
# ════════════════════════════════════════════════════════════════════════

def _save_checkpoint(out_dir: Path, cached: List[GuidedPoseResult],
                     new: List[GuidedPoseResult]) -> None:
    by_combo: Dict[str, List[GuidedPoseResult]] = defaultdict(list)
    for r in cached + new:
        by_combo[f"{r.ligand_name}__{r.protein_name}"].append(r)
    saved = 0
    for combo, results in by_combo.items():
        d = out_dir / combo
        d.mkdir(parents=True, exist_ok=True)
        try:
            with open(d / "guided_results.json", "w") as f:
                json.dump([r.to_dict() for r in results], f, indent=2)
            saved += 1
        except Exception:
            pass
    monitor.info(f"Checkpoint: saved guided_results.json for {saved} combos")


def _finalize(*, out_dir: Path, log_path: Path, existing_combos: set,
              receptor_files, ligand_files,
              all_jobs, new_results, cached_results, all_results,
              protein_fpocket, protein_p2rank,
              lig_props_cache, prot_props_cache,
              combo_phase1_times,
              phase1_time, phase2_time, phase3_time, pipeline_elapsed,
              n_cached, n_new_ok, n_new_fail, n_jobs, total_combos) -> dict:

    timing_records: List[TimingRecord] = []
    results_by_combo: Dict[str, List[GuidedPoseResult]] = defaultdict(list)
    for r in all_results:
        results_by_combo[f"{r.ligand_name}__{r.protein_name}"].append(r)

    combo_dock_times: Dict[str, float] = defaultdict(float)
    combo_post_times: Dict[str, float] = defaultdict(float)
    combo_n_attempted: Dict[str, int] = defaultdict(int)
    for job in all_jobs:
        combo_dock_times[job.combo_name] += job.dock_time
        combo_n_attempted[job.combo_name] += 1
    for r in new_results:
        combo_post_times[f"{r.ligand_name}__{r.protein_name}"] += (
            r.post_time_s + r.refine_time_s)

    for protein, ligand in itertools.product(receptor_files, ligand_files):
        pname, lname = protein.stem, ligand.stem
        combo_name = f"{lname}__{pname}"
        results = results_by_combo.get(combo_name, [])
        p1 = combo_phase1_times.get(combo_name, 0.0)
        p2 = combo_dock_times.get(combo_name, 0.0)
        p3 = combo_post_times.get(combo_name, 0.0)
        timing_records.append(TimingRecord(
            protein=pname, ligand=lname,
            phase1_prep_s=p1, phase2_dock_s=p2, phase3_post_s=p3, total_s=p1 + p2 + p3,
            n_poses_attempted=combo_n_attempted.get(combo_name, 0),
            n_poses_success=sum(1 for r in results if r.success),
        ))

        # Final per-combo JSON
        d = out_dir / combo_name
        d.mkdir(parents=True, exist_ok=True)
        with open(d / "guided_results.json", "w") as f:
            json.dump([r.to_dict() for r in results], f, indent=2)

        n_ok = sum(1 for r in results if r.success)
        n_fail = sum(1 for r in results if not r.success)
        status = "success" if n_ok > 0 else "failed"
        error = f"{n_fail} pose(s) failed" if n_fail > 0 else ""
        elapsed = p1 + p2 + p3

        lig_props = lig_props_cache.get(ligand) or get_ligand_properties(ligand)
        prot_props = prot_props_cache.get(protein) or get_protein_properties(protein)
        append_log_row(
            log_path=log_path, existing_combos=existing_combos,
            combo_name=combo_name, protein_name=pname, ligand_name=lname,
            status=status, error_reason=error, elapsed_time=elapsed,
            combo_results=results, lig_props=lig_props, prot_props=prot_props,
        )
        if status == "failed":
            record_failure(out_dir, combo_name, error, pname, lname, elapsed)

    # Unguided pose statistics
    monitor.header("STATISTICS: UNGUIDED POSES vs POCKETS")
    pockets_by_protein: Dict[str, Dict[str, List[PocketInfo]]] = defaultdict(
        lambda: {"fpocket": [], "p2rank": []})
    for pname, pks in protein_fpocket.items():
        pockets_by_protein[pname]["fpocket"] = pks
    for pname, pks in protein_p2rank.items():
        pockets_by_protein[pname]["p2rank"] = pks

    unguided = [r for r in all_results
                if r.mode == "unguided" and r.success and r.pose_centroid]
    n_in_fp = n_in_pr = n_in_either = 0
    for r in unguided:
        fp_pockets = pockets_by_protein[r.protein_name]["fpocket"]
        pr_pockets = pockets_by_protein[r.protein_name]["p2rank"]
        in_fp = any(math.dist(r.pose_centroid, p.center) <= CFG.pocket_match_threshold for p in fp_pockets)
        in_pr = any(math.dist(r.pose_centroid, p.center) <= CFG.pocket_match_threshold for p in pr_pockets)
        n_in_fp += int(in_fp)
        n_in_pr += int(in_pr)
        n_in_either += int(in_fp or in_pr)
    print(f"  Unguided poses: {len(unguided)}")
    print(f"  In fpocket: {n_in_fp}  In p2rank: {n_in_pr}  In either: {n_in_either}")

    # Timing report
    monitor.header("TIMING REPORT")
    total_phases = phase1_time + phase2_time + phase3_time
    print("\n  GLOBAL PHASE TIMING:")
    print(f"  {'Phase':<25} {'Time (s)':<12} {'Time (min)':<12} {'%':>6}")
    print(f"  {'-'*25} {'-'*12} {'-'*12} {'-'*6}")
    for label, t in [("Phase 1: CPU Prep", phase1_time),
                     ("Phase 2: GPU Dock", phase2_time),
                     ("Phase 3: CPU Post", phase3_time)]:
        pct = (t / total_phases * 100) if total_phases > 0 else 0
        print(f"  {label:<25} {t:<12.1f} {t/60:<12.1f} {pct:>5.1f}%")
    print(f"  {'TOTAL (phases)':<25} {total_phases:<12.1f} {total_phases/60:<12.1f} {'100.0%':>6}")
    print(f"  {'Pipeline wall time':<25} {pipeline_elapsed:<12.1f} {pipeline_elapsed/60:<12.1f}")

    # Timing CSV
    timing_csv = out_dir / "pipeline_timing.csv"
    pd.DataFrame([{
        "protein": t.protein, "ligand": t.ligand,
        "phase1_prep_s": round(t.phase1_prep_s, 2),
        "phase2_dock_s": round(t.phase2_dock_s, 2),
        "phase3_post_s": round(t.phase3_post_s, 2),
        "total_s": round(t.total_s, 2),
        "n_attempted": t.n_poses_attempted,
        "n_success": t.n_poses_success,
    } for t in timing_records]).to_csv(timing_csv, index=False)
    print(f"\n  Timing saved to: {timing_csv}")

    # Per-variant timing CSV — one row per finalized pose variant so the cost of
    # each axis (clamp on/off, raw vs smina/gnina re-search) is directly comparable
    # without digging through the per-combo guided_results.json files. prep/dock
    # are shared across a pose's variants; post is per clamp variant; refine is
    # per refine variant.
    variant_csv = out_dir / "pose_variant_timing.csv"
    pd.DataFrame([{
        "combo_name": f"{r.ligand_name}__{r.protein_name}",
        "protein": r.protein_name, "ligand": r.ligand_name,
        "pocket_source": r.mode, "pocket_id": r.pocket_unique_id,
        "pose_num": r.pose_num,
        "clamp_variant": r.clamp_variant, "refine_variant": r.refine_variant,
        "success": r.success,
        "prep_time_s": round(r.prep_time_s, 4),
        "dock_time_s": round(r.dock_time_s, 4),
        "post_time_s": round(r.post_time_s, 4),
        "refine_time_s": round(r.refine_time_s, 4),
        "total_time_s": round(
            r.prep_time_s + r.dock_time_s + r.post_time_s + r.refine_time_s, 4),
        "uff_minimized": r.uff_minimized,
        "refine_affinity": r.refine_affinity,
        "cnn_score": r.cnn_score,
        "cnn_affinity": r.cnn_affinity,
        "sdf_path": str(r.sdf_path) if r.sdf_path else None,
    } for r in all_results]).to_csv(variant_csv, index=False)
    print(f"  Per-variant timing saved to: {variant_csv}")

    # Concise per-axis timing summary (freshly computed, successful poses only —
    # cached/skipped poses carry no timing). Directly answers "how much does
    # clamping / the smina/gnina re-search cost per pose".
    def _avg(rows: List[GuidedPoseResult], attr: str) -> float:
        return sum(getattr(x, attr) for x in rows) / len(rows) if rows else 0.0

    timed = [r for r in new_results if r.success]
    if timed and _variant_count() > 1:
        print("\n  PER-VARIANT TIMING (mean per pose, this run):")
        by_clamp: Dict[str, List[GuidedPoseResult]] = defaultdict(list)
        for r in timed:
            by_clamp[r.clamp_variant or "(n/a)"].append(r)
        for cv, rows in sorted(by_clamp.items()):
            print(f"    clamp={cv:<9} n={len(rows):<5} "
                  f"post(correct+clamp+UFF)={_avg(rows, 'post_time_s'):.3f}s  "
                  f"refine={_avg(rows, 'refine_time_s'):.3f}s")
        by_refine: Dict[str, List[GuidedPoseResult]] = defaultdict(list)
        for r in timed:
            by_refine[r.refine_variant or "(n/a)"].append(r)
        for rv, rows in sorted(by_refine.items()):
            print(f"    refine={rv:<8} n={len(rows):<5} "
                  f"re-search={_avg(rows, 'refine_time_s'):.3f}s")

    summary = {
        "timestamp": datetime.now().isoformat(),
        "architecture": "globally decoupled 3-phase pipeline",
        "config": {
            "gpu_batch_size": CFG.gpu_batch_size,
            "n_top_pockets": CFG.n_top_pockets,
            "poses_per_pocket": CFG.poses_per_pocket,
            "n_unguided_poses": CFG.n_unguided_poses,
            "pocket_match_threshold_A": CFG.pocket_match_threshold,
            "force_pocket": CFG.force_pocket,
            "uff_minimize": CFG.uff_minimize,
            "n_parallel_workers": CFG.n_parallel_workers,
            "use_protein_cropping": CFG.use_protein_cropping,
            "clamp_mode": CFG.clamp_mode,
            "refine_mode": CFG.refine_mode,
            "refine_tool": CFG.refine_tool,
            "refine_search": CFG.smina_search,
        },
        "global_timing": {
            "pipeline_wall_time_s": round(pipeline_elapsed, 1),
            "phase1_prep_s": round(phase1_time, 1),
            "phase2_dock_s": round(phase2_time, 1),
            "phase3_post_s": round(phase3_time, 1),
        },
        "totals": {
            "combos": total_combos,
            "total_jobs_prepared": n_jobs,
            "total_cached": n_cached,
            "poses_success": n_new_ok + n_cached,
            "poses_failed": n_new_fail,
            "unguided_in_pocket": n_in_either,
        },
        "monitor_counters": {
            "equibind_calls": monitor.call_count,
            "equibind_success": monitor.success_count,
            "equibind_fail": monitor.fail_count,
        },
        "per_combo_timing": [
            {"protein": t.protein, "ligand": t.ligand,
             "prep_s": round(t.phase1_prep_s, 2),
             "dock_s": round(t.phase2_dock_s, 2),
             "post_s": round(t.phase3_post_s, 2),
             "total_s": round(t.total_s, 2)}
            for t in timing_records
        ],
        "all_results": [r.to_dict() for r in all_results],
    }
    summary_json = out_dir / "pipeline_summary.json"
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"  Summary saved to: {summary_json}")
    return summary
