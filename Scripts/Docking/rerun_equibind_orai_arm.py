#!/usr/bin/env python3
"""Replay an Orai EquiBind arm into a fresh output tree, retaining the raw poses.

WHY THIS EXISTS
Both Orai EquiBind arms ran with ``refine_mode: "on"``, which REPLACES each pose
with its refined version. No ``__refRAW.sdf`` survives in either tree, so unlike
the calibration-benchmark arm they cannot be re-refined in place: the pre-refine
input no longer exists. Reproducing the arm with ``refine_mode=both`` regenerates
the raw poses alongside the refined ones, after which
``rerun_equibind_refine.py`` can produce the matched second refiner mode from
identical inputs.

WHAT IT REPRODUCES
The launcher notebook cells are the real source of truth (``run_equibind.py``
reads no YAML; it is driven entirely by ``EQ_*`` environment variables). This
script rebuilds that environment from the same YAML the cell reads, and applies
the overrides the cell hardcodes. It never writes to the published tree.

THE PROTONATION TRAP
``config.reduce_executable`` is a bare ``shutil.which("reduce")`` with no
sibling-of-``sys.executable`` fallback, and ``_prepare_protein`` gates receptor
protonation on it. Both original runs were launched by absolute interpreter path
from a shell whose PATH lacked the equibind env, so ``reduce`` was not found and
the receptors were NOT protonated (verified on disk: 10362 atoms, 0 hydrogens).
Re-running from an activated env would find ``reduce`` and silently protonate,
changing a second axis. Two independent guards here:
  1. ``--copy-prep`` pre-copies the published ``_prep_*`` directories into the
     destination, so ``_prepare_protein`` hits its cache and never re-preps.
  2. The child PATH is stripped of the equibind env's bin directory, and the
     script refuses to launch if ``reduce`` is still resolvable.

USAGE
    python3 Scripts/Docking/rerun_equibind_orai_arm.py --arm orai_jku \
        --dest-root Dockings/equibind_results_uffoff_refineboth \
        --refine-mode both --smina-search local_only
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict

import yaml

REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "Scripts" / "Docking" / "run_equibind.py"
EQUIBIND_PYTHON = Path("/home/manndo/anaconda3/envs/equibind/bin/python")

ARMS = {
    # Orai1 receptors x PoseBusters benchmark ligands (the control arm).
    # Launcher: Master_Docking_AD_Full_Protein.ipynb cell 97 (id 6da72bd0).
    # That cell reads every value through a strict cfg_req(), so its YAML is
    # authoritative and nothing needs overriding here.
    "orai_benchmark": {
        "config": "Scripts/Docking/equibind_orai_benchmark_config.yaml",
        "published_root": "Dockings/Orai_Benchmark_Equibind",
        # The cell stages inputs; the staging dirs still exist and are reused
        # read-only rather than rebuilt, so the input set cannot drift.
        "receptors_dir": "Dockings/Orai_Benchmark_Equibind/_staging/receptors",
        "drugs_dir": "Dockings/Orai_Benchmark_Equibind/_staging/ligands",
        "cell_overrides": {},
        "thread_pins": True,   # cell 97 pins OMP/MKL/BLAS/NUMEXPR to 1
    },
    # Orai1 receptors x JKU experimental ligands (UFF-off run of 2026-08-17).
    # Launcher: cell 69 (id eqbnd-run-001), which HARDCODES five variables
    # after env.update(), so the YAML alone would not reproduce the run.
    "orai_jku": {
        "config": "Scripts/Docking/equibind_orai_jku_config.yaml",
        "published_root": "Dockings/equibind_results_uffoff/pdbqt",
        "receptors_dir": "Data/Receptors",
        "drugs_dir": "Data/Ligands/JKU/pdbqt",
        "cell_overrides": {
            "EQ_N_TOP_POCKETS": "0",
            "EQ_CLAMP_MODE": "off",
            "EQ_REFINE_MODE": "on",
            "EQ_REFINE_TOOL": "gnina",
            "EQ_GNINA": "/home/manndo/docking_tools/gnina",
        },
        "thread_pins": False,  # cell 69 sets no thread pins
    },
}


def _req(cfg: dict, key: str, path: Path):
    if key not in cfg:
        raise KeyError(f"'{key}' missing from {path}")
    return cfg[key]


def build_env(arm: dict, cfg: dict, cfg_path: Path, dest_root: Path) -> Dict[str, str]:
    """Rebuild the launcher cell's EQ_* environment for this arm."""
    env = os.environ.copy()

    # gnina's CNN rescore grabs all cores via OpenMP regardless of --cpu 1;
    # without these pins the worker pool oversubscribes and each refine crawls.
    if arm["thread_pins"]:
        for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[var] = "1"

    def r(key):
        return _req(cfg, key, cfg_path)

    env.update({
        "EQ_RECEPTORS_DIR": str((REPO / arm["receptors_dir"]).resolve()),
        "EQ_DRUGS_DIR": str((REPO / arm["drugs_dir"]).resolve()),
        "EQ_OUTPUT_DIR": str(dest_root.resolve()),
        "EQ_RECEPTOR_FILTER": str(r("receptor_name_filter")),
        "EQ_EQUIBIND_DIR": os.path.expanduser(str(r("equibind_dir"))),
        "EQ_DEVICE": str(r("device")),
        "EQ_GPU_BATCH_SIZE": str(r("gpu_batch_size")),
        "EQ_FPOCKET_DIR": str((REPO / str(r("fpocket_results_dir"))).resolve()),
        "EQ_P2RANK_DIR": str((REPO / str(r("p2rank_results_dir"))).resolve()),
        "EQ_CLAMP_MODE": str(r("clamp_mode")),
        "EQ_REFINE_MODE": str(r("refine_mode")),
        "EQ_REFINE_TOOL": str(r("refine_tool")),
        "EQ_SMINA_SEARCH": str(r("smina_search")),
        "EQ_SMINA_AUTOBOX_ADD": str(r("smina_autobox_add")),
        "EQ_SMINA_CPU": str(r("smina_cpu")),
        "EQ_SMINA_SEED": str(r("smina_seed")),
        "EQ_SMINA_TIMEOUT": str(r("smina_timeout_s")),
        "EQ_SMINA_RECEPTOR_EXT": str(r("smina_receptor_ext")),
        "EQ_GNINA_USE_GPU": str(r("gnina_use_gpu")),
        "EQ_UFF_VDW_THRESH": str(r("uff_vdw_thresh")),
        "EQ_LIGAND_EXTENSIONS": ",".join(r("ligand_extensions")),
        "EQ_N_TOP_POCKETS": str(r("n_top_pockets")),
        "EQ_POSES_PER_POCKET": str(r("poses_per_pocket")),
        "EQ_N_UNGUIDED_POSES": str(r("n_unguided_poses")),
        "EQ_POCKET_MATCH_THRESHOLD": str(r("pocket_match_threshold")),
        "EQ_FORCE_POCKET": str(r("force_pocket")),
        "EQ_CLAMP_POSE": str(r("clamp_pose_to_pocket")),
        "EQ_USE_CROPPING": str(r("use_protein_cropping")),
        "EQ_CROP_RADIUS": str(r("pocket_crop_radius")),
        "EQ_CROP_BUFFER": str(r("pocket_crop_buffer")),
        "EQ_CONFORMERS_PER_SEED": str(r("conformers_per_seed")),
        "EQ_RDKIT_SEEDS": ",".join(str(s) for s in r("rdkit_seeds")),
        "EQ_UFF_MINIMIZE": str(r("uff_minimize")),
        "EQ_UFF_MAX_ITERS": str(r("uff_max_iters")),
        "EQ_UFF_ENERGY_TOL": str(r("uff_energy_tol")),
        "EQ_UFF_FORCE_TOL": str(r("uff_force_tol")),
        "EQ_UFF_PROXIMITY_RADIUS": str(r("uff_proximity_radius")),
        "EQ_UFF_ADD_H": str(r("uff_add_hydrogens")),
        "EQ_UFF_STRIP_NONSTANDARD": str(r("uff_strip_nonstandard")),
        "EQ_WORKERS": str(r("n_parallel_workers") or os.cpu_count() or 8),
        "EQ_PARALLEL_CONFORMER_GEN": str(r("parallel_conformer_gen")),
        "EQ_SKIP_EXISTING": str(r("skip_existing")),
    })
    smina_exe = str(r("smina_executable")).strip()
    gnina_exe = str(r("gnina_executable")).strip()
    if smina_exe:
        env["EQ_SMINA"] = smina_exe
    if gnina_exe:
        env["EQ_GNINA"] = gnina_exe

    # Values the launcher cell sets in code, after the YAML forward.
    env.update(arm["cell_overrides"])
    return env


def strip_env_bin_from_path(env: Dict[str, str]) -> str:
    """Drop the equibind env's bin/ from PATH so `reduce` stays unresolvable."""
    env_bin = str(EQUIBIND_PYTHON.parent)
    parts = [p for p in env.get("PATH", "").split(os.pathsep)
             if p and Path(p).resolve() != Path(env_bin).resolve()]
    return os.pathsep.join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--arm", required=True, choices=sorted(ARMS))
    ap.add_argument("--dest-root", required=True, type=Path)
    ap.add_argument("--refine-mode", default="both", choices=("on", "off", "both"))
    ap.add_argument("--smina-search", default="local_only",
                    choices=("local_only", "minimize"))
    ap.add_argument("--copy-prep", action="store_true", default=True,
                    help="Pre-copy the published _prep_* dirs (default: on).")
    ap.add_argument("--no-copy-prep", dest="copy_prep", action="store_false")
    ap.add_argument("--workers", type=int, default=None,
                    help="Override EQ_WORKERS (the YAML value otherwise). Use to leave\n                         cores free for concurrent analysis.")
    ap.add_argument("--log-file", type=Path, default=None)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    arm = ARMS[args.arm]
    cfg_path = REPO / arm["config"]
    cfg = yaml.safe_load(cfg_path.read_text())
    published = REPO / arm["published_root"]
    dest_root = (args.dest_root if args.dest_root.is_absolute()
                 else REPO / args.dest_root).resolve()

    if dest_root == published.resolve():
        print("ERROR: --dest-root must not be the published tree", file=sys.stderr)
        return 2

    env = build_env(arm, cfg, cfg_path, dest_root)
    if args.workers:
        env["EQ_WORKERS"] = str(args.workers)
    env["PATH"] = strip_env_bin_from_path(env)
    # Both CLI flags override the environment inside run_equibind.py.
    argv = [str(EQUIBIND_PYTHON), str(RUNNER),
            "--refine-mode", args.refine_mode,
            "--smina-search", args.smina_search]

    # Guard the protonation axis: refuse to run if reduce became resolvable.
    resolved_reduce = shutil.which("reduce", path=env["PATH"])
    print(f"arm            : {args.arm}")
    print(f"config         : {cfg_path.relative_to(REPO)}")
    print(f"published tree : {published.relative_to(REPO)}  (READ-ONLY)")
    print(f"dest tree      : {dest_root}")
    print(f"refine_mode    : {args.refine_mode}   smina_search: {args.smina_search}")
    print(f"refine_tool    : {env['EQ_REFINE_TOOL']}   gnina_gpu: {env['EQ_GNINA_USE_GPU']}")
    print(f"workers        : {env['EQ_WORKERS']}   uff_minimize: {env['EQ_UFF_MINIMIZE']}")
    print(f"receptors      : {env['EQ_RECEPTORS_DIR']}")
    print(f"ligands        : {env['EQ_DRUGS_DIR']}")
    print(f"reduce on PATH : {resolved_reduce or '(none - receptors stay unprotonated)'}")
    if resolved_reduce:
        print("ERROR: `reduce` is resolvable in the child PATH. The original runs did "
              "NOT protonate their receptors; continuing would change a second axis.",
              file=sys.stderr)
        return 3

    # Reuse the published prepared receptors so _prepare_protein hits its cache.
    if args.copy_prep:
        preps = sorted(published.glob("_prep_*"))
        dest_root.mkdir(parents=True, exist_ok=True)
        copied = 0
        for p in preps:
            target = dest_root / p.name
            if not target.exists():
                shutil.copytree(p, target)
                copied += 1
        print(f"prep dirs      : {len(preps)} published, {copied} copied into dest "
              f"(receptor prep will hit cache)")

    if args.dry_run:
        print("\n[dry-run] argv:", " ".join(argv))
        return 0

    log_file = args.log_file or (dest_root / f"rerun_{args.arm}_{args.smina_search}.log")
    log_file.parent.mkdir(parents=True, exist_ok=True)
    print(f"\nlogging to     : {log_file}\nlaunching...", flush=True)

    with open(log_file, "w") as lf:
        proc = subprocess.run(argv, env=env, cwd=str(REPO), stdout=lf,
                              stderr=subprocess.STDOUT)

    print(f"\nreturn code: {proc.returncode}")
    summary = dest_root / "pipeline_summary.json"
    if summary.exists():
        s = json.loads(summary.read_text())
        print(f"  config : {json.dumps(s.get('config'))}")
        print(f"  totals : {json.dumps(s.get('totals'))}")
        print(f"  timing : {json.dumps(s.get('global_timing'))}")
    else:
        print(f"  no pipeline_summary.json at {summary} - check {log_file}")
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main())
