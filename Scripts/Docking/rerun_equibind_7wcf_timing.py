#!/usr/bin/env python3
"""Re-dock 7WCF_ACP with EquiBind to recover its lost wall-clock record.

Why this exists. The 2026-07-10 benchmark EquiBind campaign was interrupted
(SIGINT) partway through 7WCF_ACP's phase-3 post-processing, so that complex
never wrote a ``pipeline_summary.json``. On restart the launcher's resume test
(cell b1a9baa1) skipped it, because that test only asks whether ``__refGNINA``
pose SDFs exist on disk -- 75 of them did -- and not whether the pipeline ever
finished. The result is the one complex of 308 with no timing, which is why
docking_effort_comparison.py's cross-tool intersection reports 302 rather than
303 complexes.

This replays cell b1a9baa1's environment for that single complex, reading
``equibind_benchmark_config.yaml`` (the split file validated key-by-key against
the benchmark's own pipeline_summary records) instead of the deprecated shared
YAML.

It writes to a SEPARATE output root. The canonical
``Dockings/Benchmark_Equibind/7WCF_ACP`` holds the truncated 248-pose set that
the published accuracy numbers were computed from, and a clean re-dock produces
the full 270, so overwriting it in place would silently change results that are
already reported.

Takes a complex id so a control complex -- one whose July wall time IS on record
-- can be replayed under the same conditions to check that today's machine
reproduces the campaign before the recovered number is trusted. ``gnina_use_gpu``
is overridable for the same reason: the runner never logged it, so the benchmark
YAML's value is an inference rather than a recorded fact.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import yaml

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--pdb-id", default="7WCF_ACP")
ap.add_argument("--gnina-gpu", choices=["true", "false"], default=None,
                help="override gnina_use_gpu (default: the YAML's value)")
ap.add_argument("--tag", default=None,
                help="output subdirectory suffix, to keep conditions apart")
cli = ap.parse_args()

PDB_ID = cli.pdb_id
REPO = Path(__file__).resolve().parents[2]
BENCHMARK_DIR = REPO / "Data/PoseBuster Benchmark Set"
CONFIG_PATH = REPO / "Scripts/Docking/equibind_benchmark_config.yaml"
RUNNER = REPO / "Scripts/Docking/run_equibind.py"
OUTPUT_BASE = REPO / "Dockings/Benchmark_Equibind_timing_rerun"
LOG_DIR = REPO / "Dockings/Logs/benchmark_equibind_timing_rerun_logs"
EQUIBIND_PYTHON = Path("/home/manndo/anaconda3/envs/equibind/bin/python")

cfg = yaml.safe_load(CONFIG_PATH.read_text())

complex_dir = BENCHMARK_DIR / PDB_ID
protein_pdb = complex_dir / f"{PDB_ID}_protein.pdb"
ligand_sdf = complex_dir / f"{PDB_ID}_ligand_start_conf.sdf"
for p in (protein_pdb, ligand_sdf, RUNNER, EQUIBIND_PYTHON):
    if not p.exists():
        sys.exit(f"missing required input: {p}")

run_name = PDB_ID if not cli.tag else f"{PDB_ID}__{cli.tag}"
complex_out = OUTPUT_BASE / run_name
staging = OUTPUT_BASE / "_staging" / run_name
rec_staging, lig_staging = staging / "receptors", staging / "ligands"
for d in (complex_out, rec_staging, lig_staging, LOG_DIR):
    d.mkdir(parents=True, exist_ok=True)

for link, target in ((rec_staging / protein_pdb.name, protein_pdb),
                     (lig_staging / ligand_sdf.name, ligand_sdf)):
    if not link.exists():
        link.symlink_to(target.resolve())

# refine_mode / refine_tool are hardcoded to "both" in the launcher cell rather
# than read from the YAML; the benchmark config carries the same pair, so read
# them here and assert the match instead of silently diverging.
refine_mode, refine_tool = "both", "both"
assert cfg["refine_mode"] == refine_mode and cfg["refine_tool"] == refine_tool

env = os.environ.copy()
env.update({
    "EQ_RECEPTORS_DIR": str(rec_staging.resolve()),
    "EQ_DRUGS_DIR": str(lig_staging.resolve()),
    "EQ_OUTPUT_DIR": str(complex_out.resolve()),
    "EQ_RECEPTOR_FILTER": "",
    "EQ_EQUIBIND_DIR": os.path.expanduser(cfg["equibind_dir"]),
    "EQ_DEVICE": cfg["device"],
    "EQ_GPU_BATCH_SIZE": str(cfg["gpu_batch_size"]),
    "EQ_FPOCKET_DIR": str((REPO / cfg["fpocket_results_dir"]).resolve()),
    "EQ_P2RANK_DIR": str((REPO / cfg["p2rank_results_dir"]).resolve()),
    "EQ_CLAMP_MODE": str(cfg["clamp_mode"]),
    "EQ_REFINE_MODE": refine_mode,
    "EQ_REFINE_TOOL": refine_tool,
    "EQ_SMINA_SEARCH": str(cfg["smina_search"]),
    "EQ_SMINA_AUTOBOX_ADD": str(cfg["smina_autobox_add"]),
    "EQ_SMINA_CPU": str(cfg["smina_cpu"]),
    "EQ_SMINA_SEED": str(cfg["smina_seed"]),
    "EQ_SMINA_TIMEOUT": str(cfg["smina_timeout_s"]),
    "EQ_SMINA_RECEPTOR_EXT": str(cfg["smina_receptor_ext"]),
    "EQ_GNINA_USE_GPU": (cli.gnina_gpu.capitalize() if cli.gnina_gpu
                         else str(cfg["gnina_use_gpu"])),
    "EQ_UFF_VDW_THRESH": str(cfg["uff_vdw_thresh"]),
    "EQ_N_TOP_POCKETS": str(cfg["n_top_pockets"]),
    "EQ_POSES_PER_POCKET": str(cfg["poses_per_pocket"]),
    "EQ_N_UNGUIDED_POSES": str(cfg["n_unguided_poses"]),
    "EQ_POCKET_MATCH_THRESHOLD": str(cfg["pocket_match_threshold"]),
    "EQ_FORCE_POCKET": str(cfg["force_pocket"]),
    "EQ_CLAMP_POSE": str(cfg["clamp_pose_to_pocket"]),
    "EQ_USE_CROPPING": str(cfg["use_protein_cropping"]),
    "EQ_CROP_RADIUS": str(cfg["pocket_crop_radius"]),
    "EQ_CROP_BUFFER": str(cfg["pocket_crop_buffer"]),
    "EQ_CONFORMERS_PER_SEED": str(cfg["conformers_per_seed"]),
    "EQ_RDKIT_SEEDS": ",".join(str(s) for s in cfg["rdkit_seeds"]),
    "EQ_UFF_MINIMIZE": str(cfg["uff_minimize"]),
    "EQ_UFF_MAX_ITERS": str(cfg["uff_max_iters"]),
    "EQ_UFF_ENERGY_TOL": str(cfg["uff_energy_tol"]),
    "EQ_UFF_FORCE_TOL": str(cfg["uff_force_tol"]),
    "EQ_UFF_PROXIMITY_RADIUS": str(cfg["uff_proximity_radius"]),
    "EQ_UFF_ADD_H": str(cfg["uff_add_hydrogens"]),
    "EQ_UFF_STRIP_NONSTANDARD": str(cfg["uff_strip_nonstandard"]),
    "EQ_WORKERS": str(cfg["n_parallel_workers"]),
    "EQ_PARALLEL_CONFORMER_GEN": str(cfg["parallel_conformer_gen"]),
    "EQ_SKIP_EXISTING": str(cfg["skip_existing"]),
})
if str(cfg.get("smina_executable", "")).strip():
    env["EQ_SMINA"] = str(cfg["smina_executable"]).strip()
env["EQ_GNINA"] = (str(cfg.get("gnina_executable", "")).strip()
                   or "/home/manndo/docking_tools/gnina")

log_file = LOG_DIR / f"{run_name}.log"
print(f"re-docking {PDB_ID} -> {complex_out}")
print(f"log: {log_file}")
with open(log_file, "w") as lf:
    proc = subprocess.run([str(EQUIBIND_PYTHON), str(RUNNER)],
                          env=env, cwd=str(REPO), stdout=lf,
                          stderr=subprocess.STDOUT)

summary_file = complex_out / "pipeline_summary.json"
if proc.returncode != 0 or not summary_file.exists():
    sys.exit(f"FAILED rc={proc.returncode}; see {log_file}")

s = json.loads(summary_file.read_text())
gt, tot = s.get("global_timing", {}), s.get("totals", {})
print(f"  poses ok={tot.get('poses_success')} failed={tot.get('poses_failed')}")
print(f"  wall={gt.get('pipeline_wall_time_s')}s  "
      f"prep={gt.get('phase1_prep_s')}s dock={gt.get('phase2_dock_s')}s "
      f"post={gt.get('phase3_post_s')}s")
