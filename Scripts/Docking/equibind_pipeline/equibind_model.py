"""EquiBind model loader, receptor-graph cache, and GPU batch runner.

Importing this module:
    1. Patches ``dgl.function.copy_edge`` -> ``copy_e`` for newer DGL.
    2. Inserts ``CFG.equibind_dir`` into ``sys.path``.
    3. Lazily loads the EquiBind model on first call to :func:`get_model`.
"""
from __future__ import annotations

import sys
import threading
import time
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import yaml

import dgl  # noqa: F401
import dgl.function as _dgl_fn

from rdkit import Chem
from rdkit.Geometry import Point3D

from .config import CFG
from .monitor import monitor
from .types import DockJob


# ── DGL compatibility shim ──────────────────────────────────────────────
if not hasattr(_dgl_fn, "copy_edge"):
    _dgl_fn.copy_edge = _dgl_fn.copy_e
    print("  \u2713 Patched dgl.function.copy_edge \u2192 copy_e (DGL compat shim)")

# ── Make EquiBind importable ────────────────────────────────────────────
_eb_dir = str(CFG.equibind_dir)
if _eb_dir not in sys.path:
    sys.path.insert(0, _eb_dir)

from models.equibind import EquiBind as EquiBindModel  # noqa: E402
from commons.process_mols import (  # noqa: E402
    get_receptor_inference,
    get_rec_graph,
    get_lig_graph_revised,
    get_geometry_graph,
)
from commons.geometry_utils import (  # noqa: E402
    rigid_transform_Kabsch_3D,
    get_torsions,
    get_dihedral_vonMises,
    apply_changes,
)
from commons.utils import seed_all  # noqa: E402


# ── Lazy singletons ──────────────────────────────────────────────────────
_model = None
_args: Optional[SimpleNamespace] = None
_device: Optional[torch.device] = None
_load_lock = threading.Lock()

_rec_graph_cache: Dict[str, object] = {}
_rec_graph_lock = threading.Lock()


def _load_model_inner() -> Tuple[object, SimpleNamespace, torch.device]:
    t0 = time.time()
    ckpt_path = CFG.equibind_dir / "runs" / "flexible_self_docking" / "best_checkpoint.pt"
    train_args_path = ckpt_path.parent / "train_arguments.yaml"

    with open(train_args_path) as f:
        train_args = yaml.safe_load(f)
    train_args["model_parameters"]["noise_initial"] = 0

    args = SimpleNamespace(**train_args)
    args.checkpoint = str(ckpt_path)
    args.use_rdkit_coords = args.dataset_params.get("use_rdkit_coords", True)
    args.device = CFG.equibind_device

    dev = torch.device("cuda:0" if torch.cuda.is_available() and CFG.equibind_device == "cuda" else "cpu")
    checkpoint = torch.load(str(ckpt_path), map_location=dev)

    model = EquiBindModel(
        device=dev,
        lig_input_edge_feats_dim=15,
        rec_input_edge_feats_dim=27,
        **args.model_parameters,
    )
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(dev)
    model.eval()
    seed_all(args.seed)

    monitor.info(f"EquiBind model loaded from {ckpt_path}")
    monitor.info(f"  Device: {dev}  |  Load time: {time.time() - t0:.2f}s")
    return model, args, dev


def get_model() -> Tuple[object, SimpleNamespace, torch.device]:
    """Lazy-load the EquiBind model (thread-safe)."""
    global _model, _args, _device
    if _model is None:
        with _load_lock:
            if _model is None:
                _model, _args, _device = _load_model_inner()
    return _model, _args, _device


def get_args() -> SimpleNamespace:
    return get_model()[1]


def _build_receptor_graph(protein_pdb: Path, args: SimpleNamespace):
    dp = args.dataset_params
    rec, rec_coords, c_alpha_coords, n_coords, c_coords = get_receptor_inference(str(protein_pdb))
    return get_rec_graph(
        rec, rec_coords, c_alpha_coords, n_coords, c_coords,
        use_rec_atoms=dp["use_rec_atoms"],
        rec_radius=dp["rec_graph_radius"],
        surface_max_neighbors=dp["surface_max_neighbors"],
        surface_graph_cutoff=dp["surface_graph_cutoff"],
        surface_mesh_cutoff=dp["surface_mesh_cutoff"],
        c_alpha_max_neighbors=dp["c_alpha_max_neighbors"],
    )


def get_receptor_graph_cached(protein_pdb: Path) -> object:
    """Build (or fetch) the EquiBind receptor graph for a protein."""
    key = str(protein_pdb)
    if key in _rec_graph_cache:
        return _rec_graph_cache[key]
    with _rec_graph_lock:
        if key not in _rec_graph_cache:
            _, args, _ = get_model()
            _rec_graph_cache[key] = _build_receptor_graph(protein_pdb, args)
    return _rec_graph_cache[key]


def receptor_graph_count() -> int:
    return len(_rec_graph_cache)


# ── Post-EquiBind corrections (torsion fitting + Kabsch) ────────────────

def run_corrections_inproc(lig: Chem.Mol, lig_coord: torch.Tensor,
                           predicted_coords: torch.Tensor) -> Chem.Mol:
    input_coords = lig_coord.detach().cpu().numpy()
    prediction = predicted_coords.detach().cpu().numpy()

    lig_input = deepcopy(lig)
    conf = lig_input.GetConformer()
    for i in range(lig_input.GetNumAtoms()):
        x, y, z = input_coords[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    lig_equibind = deepcopy(lig)
    conf = lig_equibind.GetConformer()
    for i in range(lig_equibind.GetNumAtoms()):
        x, y, z = prediction[i]
        conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))

    coords_pred = lig_equibind.GetConformer().GetPositions()
    Z_pt_cloud = coords_pred
    rotable_bonds = get_torsions([lig_input])
    new_dihedrals = np.zeros(len(rotable_bonds))
    for idx_t, r in enumerate(rotable_bonds):
        new_dihedrals[idx_t] = get_dihedral_vonMises(
            lig_input, lig_input.GetConformer(), r, Z_pt_cloud)
    optimized_mol = apply_changes(lig_input, new_dihedrals, rotable_bonds)
    optimized_conf = optimized_mol.GetConformer()
    coords_pred_optimized = optimized_conf.GetPositions()
    R, t = rigid_transform_Kabsch_3D(coords_pred_optimized.T, coords_pred.T)
    coords_pred_optimized = (R @ coords_pred_optimized.T).T + t.squeeze()
    for i in range(optimized_mol.GetNumAtoms()):
        x, y, z = coords_pred_optimized[i]
        optimized_conf.SetAtomPosition(i, Point3D(float(x), float(y), float(z)))
    return optimized_mol


# ── GPU batch inference ─────────────────────────────────────────────────

def run_gpu_batch(jobs: List[DockJob]) -> List[DockJob]:
    """Run EquiBind inference on a batch of DockJobs sharing a receptor graph."""
    if not jobs:
        return jobs

    model, _, device = get_model()
    rec_graph = _rec_graph_cache[jobs[0].protein_key]

    t_batch = time.time()
    try:
        rec_graph_dev = rec_graph.to(device)
        lig_graphs = [j.lig_graph.to(device) for j in jobs]
        geom_graphs = [g.to(device) if g is not None else None for g in (j.geometry_graph for j in jobs)]

        with torch.no_grad():
            for i, job in enumerate(jobs):
                jt0 = time.time()
                try:
                    seed_all(job.seed)
                    predictions = model(lig_graphs[i], rec_graph_dev, geom_graphs[i])
                    job.predicted_coords = predictions[0][0]
                    job.gpu_success = True
                except Exception as e:
                    job.error = f"GPU inference failed: {e}"
                    job.gpu_success = False
                job.dock_time = time.time() - jt0
    except Exception as e:
        elapsed = time.time() - t_batch
        for job in jobs:
            job.error = f"Batch inference failed: {e}"
            job.gpu_success = False
            job.dock_time = elapsed / len(jobs)

    return jobs
