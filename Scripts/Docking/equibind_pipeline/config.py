"""Centralised configuration for the EquiBind docking pipeline.

All tunables live here. Override via environment variables (``EQ_*``) or by
mutating ``CFG`` before calling :func:`equibind_pipeline.pipeline.run_pipeline`.
"""
from __future__ import annotations

import os
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List


def _env_int(name: str, default: int) -> int:
    v = os.environ.get(name)
    return int(v) if v else default


def _env_float(name: str, default: float) -> float:
    v = os.environ.get(name)
    return float(v) if v else default


def _env_bool(name: str, default: bool) -> bool:
    v = os.environ.get(name)
    if v is None:
        return default
    return v.strip().lower() in ("1", "true", "yes", "on")


def _env_path(name: str, default: Path) -> Path:
    v = os.environ.get(name)
    return Path(v) if v else default


_TRISTATE = ("on", "off", "both")


def _env_tristate(name: str, default: str) -> str:
    """Parse an enable/disable/both switch. Accepts on|off|both (+ a few aliases)."""
    v = os.environ.get(name)
    if not v:
        return default
    s = v.strip().lower()
    aliases = {"true": "on", "yes": "on", "enable": "on", "enabled": "on", "1": "on",
               "false": "off", "no": "off", "disable": "off", "disabled": "off", "0": "off"}
    s = aliases.get(s, s)
    if s not in _TRISTATE:
        raise ValueError(f"{name} must be one of {_TRISTATE} (got {v!r})")
    return s


def _env_ext_list(name: str, default: List[str]) -> List[str]:
    """Parse a comma/space-separated extension list; normalize to lower-case '.ext'."""
    v = os.environ.get(name)
    if not v:
        return [e if e.startswith(".") else f".{e}" for e in default]
    parts = [p.strip().lower() for p in v.replace(",", " ").split() if p.strip()]
    return [p if p.startswith(".") else f".{p}" for p in parts]


def _env_int_list(name: str, default: List[int]) -> List[int]:
    """Parse a comma/space-separated integer list (e.g. RDKit seeds)."""
    v = os.environ.get(name)
    if not v:
        return list(default)
    parts = [p.strip() for p in v.replace(",", " ").split() if p.strip()]
    try:
        return [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(f"{name} must be a list of integers (got {v!r})") from e


_REFINE_TOOLS = ("smina", "gnina", "both")


def _env_refine_tool(name: str, default: str) -> str:
    """Parse the re-search backend selector: smina | gnina | both (alias: all)."""
    v = os.environ.get(name)
    s = (v or default).strip().lower()
    if s == "all":
        s = "both"
    if s not in _REFINE_TOOLS:
        raise ValueError(f"{name} must be one of {_REFINE_TOOLS} (got {v!r})")
    return s


def _find_exe(name: str, env_var: str) -> str:
    """Locate a re-search executable, in priority order:

    1. ``<env_var>`` env var (explicit override).
    2. ``name`` on ``PATH``.
    3. ``name`` next to the running Python interpreter (``sys.executable``'s
       dir = the active conda env's ``bin/``). This covers the common launch
       pattern where the env's Python is invoked by absolute path without
       activating the env, so ``PATH`` lacks ``…/envs/<env>/bin`` even though
       the tool was installed into that very env.

    Returns "" when the tool cannot be found (validate() reports it only if a
    re-search that needs it is actually enabled).
    """
    explicit = os.environ.get(env_var)
    if explicit:
        return explicit
    on_path = shutil.which(name)
    if on_path:
        return on_path
    sibling = Path(sys.executable).resolve().parent / name
    if sibling.exists():
        return str(sibling)
    return ""


def _find_smina() -> str:
    return _find_exe("smina", "EQ_SMINA")


def _find_gnina() -> str:
    return _find_exe("gnina", "EQ_GNINA")


@dataclass
class Config:
    # ── Workspace / inputs ────────────────────────────────────────────────
    workspace_root: Path = field(default_factory=Path.cwd)
    drugs_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_DRUGS_DIR", Path("storage/ligands_sdf_small_symmetric_approved")))
    receptors_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_RECEPTORS_DIR", Path.cwd() / "Orai"))
    receptor_name_filter: str = os.environ.get("EQ_RECEPTOR_FILTER", "_cleaned")
    ligand_extensions: List[str] = field(default_factory=lambda: _env_ext_list(
        "EQ_LIGAND_EXTENSIONS", [".sdf", ".pdb"]))

    # ── Pocket result dirs ────────────────────────────────────────────────
    fpocket_results_folder: Path = field(default_factory=lambda: _env_path(
        "EQ_FPOCKET_DIR", Path.cwd() / "pocket_results" / "fpocket_results"))
    p2rank_folder: Path = field(default_factory=lambda: _env_path(
        "EQ_P2RANK_DIR", Path.cwd() / "pocket_results" / "p2rank_results"))

    # ── EquiBind ──────────────────────────────────────────────────────────
    equibind_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_EQUIBIND_DIR", Path.home() / "tools" / "EquiBind"))
    equibind_device: str = os.environ.get("EQ_DEVICE", "cuda")
    gpu_batch_size: int = _env_int("EQ_GPU_BATCH_SIZE", 8)

    # ── Output ────────────────────────────────────────────────────────────
    pocket_guided_output_dir: Path = field(default_factory=lambda: _env_path(
        "EQ_OUTPUT_DIR", Path.cwd() / "small_equibind_pocket_guided"))

    # ── Pocket-guided settings ────────────────────────────────────────────
    n_top_pockets: int = _env_int("EQ_N_TOP_POCKETS", 3)
    poses_per_pocket: int = _env_int("EQ_POSES_PER_POCKET", 5)
    n_unguided_poses: int = _env_int("EQ_N_UNGUIDED_POSES", 5)
    pocket_match_threshold: float = _env_float("EQ_POCKET_MATCH_THRESHOLD", 8.0)

    # ── Pocket enforcement ────────────────────────────────────────────────
    force_pocket: bool = _env_bool("EQ_FORCE_POCKET", True)
    clamp_pose_to_pocket: bool = _env_bool("EQ_CLAMP_POSE", True)
    # Centroid-clamping mode (enable / disable / both):
    #   "on"   → legacy behaviour: honour force_pocket + clamp_pose_to_pocket.
    #   "off"  → never clamp and never reject; keep the raw EquiBind centroid.
    #   "both" → emit two pose variants per job, one clamped ("__clampON") and
    #            one left at the predicted centroid ("__clampOFF"), so the
    #            effect of clamping on PB-validity can be measured directly.
    clamp_mode: str = _env_tristate("EQ_CLAMP_MODE", "on")

    # ── Docking / force-field re-search (post-pose refinement) ─────────────
    # A local docking re-search that can move the ligand as a body out of
    # protein clashes (UFF alone cannot un-bury a deeply overlapping pose).
    #   "off"  → no re-search (default; back-compat).
    #   "on"   → replace each pose with its re-searched ("refined") pose.
    #   "both" → emit both the un-refined ("__refRAW") and refined
    #            ("__refSMINA"/"__refGNINA") pose so the gain can be measured.
    refine_mode: str = _env_tristate("EQ_REFINE_MODE", "off")
    # Re-search backend: "smina" (physics), "gnina" (physics + CNN rescore), or
    # "both" (emit one refined variant per tool, e.g. __refSMINA + __refGNINA).
    refine_tool: str = _env_refine_tool("EQ_REFINE_TOOL", "smina")
    # Local optimisation mode (shared by both backends): "local_only" (local MC +
    # scoring; recommended) or "minimize" (pure energy minimisation, even more
    # local). Named smina_* for back-compat; applies to gnina identically.
    smina_search: str = os.environ.get("EQ_SMINA_SEARCH", "local_only").strip().lower()
    smina_autobox_add: float = _env_float("EQ_SMINA_AUTOBOX_ADD", 4.0)
    smina_cpu: int = _env_int("EQ_SMINA_CPU", 1)
    smina_seed: int = _env_int("EQ_SMINA_SEED", 0)
    smina_timeout_s: int = _env_int("EQ_SMINA_TIMEOUT", 300)
    # Receptor file for the re-search: the prepared protein PDB by default. smina
    # and gnina both read PDB/PDBQT directly; point this at a pdbqt via
    # EQ_SMINA_RECEPTOR_EXT if you prefer the meeko/mgltools receptors.
    smina_receptor_ext: str = os.environ.get("EQ_SMINA_RECEPTOR_EXT", ".pdb").strip().lower()
    smina_executable: str = field(default_factory=_find_smina)
    # gnina re-search backend. gnina runs a CNN by default; keep it off the GPU
    # (--no_gpu) unless gnina_use_gpu is set, since the GPU is usually busy with
    # EquiBind inference. Only needed when refine_tool is "gnina" or "both".
    gnina_executable: str = field(default_factory=_find_gnina)
    gnina_use_gpu: bool = _env_bool("EQ_GNINA_USE_GPU", False)

    # ── Protein cropping ──────────────────────────────────────────────────
    use_protein_cropping: bool = _env_bool("EQ_USE_CROPPING", True)
    pocket_crop_radius: float = _env_float("EQ_CROP_RADIUS", 8.0)
    pocket_crop_buffer: float = _env_float("EQ_CROP_BUFFER", 3.0)
    include_full_residues: bool = True

    # ── RDKit conformer generation ────────────────────────────────────────
    rdkit_seeds: List[int] = field(default_factory=lambda: _env_int_list(
        "EQ_RDKIT_SEEDS", [42, 123, 456, 789, 1001, 2022, 3141, 5926, 8675, 9999]))
    conformers_per_seed: int = _env_int("EQ_CONFORMERS_PER_SEED", 1)

    # ── UFF minimization ──────────────────────────────────────────────────
    uff_minimize: bool = _env_bool("EQ_UFF_MINIMIZE", True)
    uff_max_iters: int = _env_int("EQ_UFF_MAX_ITERS", 200)
    uff_energy_tol: float = _env_float("EQ_UFF_ENERGY_TOL", 1e-4)
    uff_force_tol: float = _env_float("EQ_UFF_FORCE_TOL", 1e-3)
    uff_proximity_radius: float = _env_float("EQ_UFF_PROXIMITY_RADIUS", 6.0)
    uff_add_hydrogens: bool = _env_bool("EQ_UFF_ADD_H", True)
    # vdW cutoff for UFF: pairs farther than vdwThresh x sum-of-vdw-radii are
    # dropped. The old default 0.1 excluded essentially all vdW interactions
    # (even active clashes), so the minimizer could not feel — or relieve —
    # protein overlaps. 10.0 is RDKit's default and includes the clash range.
    uff_vdw_thresh: float = _env_float("EQ_UFF_VDW_THRESH", 10.0)
    uff_strip_nonstandard: bool = _env_bool("EQ_UFF_STRIP_NONSTANDARD", True)

    # ── Parallelism ───────────────────────────────────────────────────────
    n_parallel_workers: int = _env_int(
        "EQ_WORKERS", min(os.cpu_count() or 22, 26))
    parallel_conformer_gen: bool = _env_bool("EQ_PARALLEL_CONFORMER_GEN", True)

    # ── Skip / overwrite ──────────────────────────────────────────────────
    skip_existing: bool = _env_bool("EQ_SKIP_EXISTING", True)

    # ── Filenames ─────────────────────────────────────────────────────────
    docking_log_filename: str = "docking_log.csv"
    error_log_filename: str = "failed_docking.json"

    # ── Derived / discovered ──────────────────────────────────────────────
    reduce_executable: str = field(default_factory=lambda: shutil.which("reduce") or "")

    def refine_tool_list(self) -> List[str]:
        """Re-search backends to run: ['smina'] | ['gnina'] | ['smina','gnina']."""
        return ["smina", "gnina"] if self.refine_tool == "both" else [self.refine_tool]

    def _refine_executable(self, tool: str) -> str:
        return self.gnina_executable if tool == "gnina" else self.smina_executable

    def validate(self) -> None:
        if not self.drugs_dir.exists():
            raise FileNotFoundError(f"Ligand directory missing: {self.drugs_dir}")
        if not self.receptors_dir.exists():
            raise FileNotFoundError(f"Receptor directory missing: {self.receptors_dir}")
        if self.clamp_mode not in _TRISTATE:
            raise ValueError(f"clamp_mode must be one of {_TRISTATE} (got {self.clamp_mode!r})")
        if self.refine_mode not in _TRISTATE:
            raise ValueError(f"refine_mode must be one of {_TRISTATE} (got {self.refine_mode!r})")
        if self.refine_mode != "off":
            if self.refine_tool not in _REFINE_TOOLS:
                raise ValueError(
                    f"refine_tool must be one of {_REFINE_TOOLS} (got {self.refine_tool!r})")
            if self.smina_search not in ("local_only", "minimize"):
                raise ValueError(
                    f"smina_search must be 'local_only' or 'minimize' (got {self.smina_search!r})")
            _install_hint = {
                "smina": "`conda install -c conda-forge smina` or download the static binary from "
                         "https://sourceforge.net/projects/smina/, then put it on PATH or set EQ_SMINA",
                "gnina": "download the static binary from https://github.com/gnina/gnina/releases "
                         "(or build it), then put it on PATH or set EQ_GNINA",
            }
            for tool in self.refine_tool_list():
                exe = self._refine_executable(tool)
                if not exe or not (shutil.which(exe) or Path(exe).exists()):
                    raise FileNotFoundError(
                        f"refine_mode is enabled with refine_tool={self.refine_tool!r} but the "
                        f"{tool!r} executable was not found. Install it: {_install_hint[tool]}"
                        f"=/path/to/{tool}.")

    def ensure_output_dirs(self) -> None:
        self.pocket_guided_output_dir.mkdir(parents=True, exist_ok=True)

    def banner(self) -> None:
        print("=" * 80)
        print("EquiBind Batch Docking Configuration")
        print("=" * 80)
        print(f"  Workspace:        {self.workspace_root}")
        print(f"  Ligand dir:       {self.drugs_dir}")
        print(f"  Ligand exts:      {', '.join(self.ligand_extensions)}")
        print(f"  Receptor dir:     {self.receptors_dir}")
        print(f"  EquiBind dir:     {self.equibind_dir}")
        print(f"  EquiBind device:  {self.equibind_device}")
        print(f"  Output:           {self.pocket_guided_output_dir}")
        print(f"  Pocket results:   fpocket={self.fpocket_results_folder}")
        print(f"                    p2rank={self.p2rank_folder}")
        print(f"  Top pockets:      {self.n_top_pockets} per method")
        print(f"  Poses/pocket:     {self.poses_per_pocket}")
        print(f"  Unguided poses:   {self.n_unguided_poses}")
        print(f"  Workers:          {self.n_parallel_workers}")
        print(f"  GPU batch size:   {self.gpu_batch_size}")
        print(f"  UFF minimize:     {self.uff_minimize}")
        print(f"  Skip existing:    {self.skip_existing}")
        print(f"  Force pocket:     {self.force_pocket}")
        print(f"  Clamp mode:       {self.clamp_mode}  "
              f"(clamp_pose_to_pocket={self.clamp_pose_to_pocket})")
        if self.refine_mode == "off":
            print(f"  Re-search:        off")
        else:
            tools = "+".join(self.refine_tool_list())
            print(f"  Re-search:        {self.refine_mode}  "
                  f"(tool={tools}/{self.smina_search}, "
                  f"autobox_add={self.smina_autobox_add}\u00c5)")
            for tool in self.refine_tool_list():
                line = f"                    {tool}={self._refine_executable(tool) or 'NOT FOUND'}"
                if tool == "gnina":
                    line += f"  (GPU={'on' if self.gnina_use_gpu else 'off'})"
                print(line)
        print(f"  Protein cropping: {self.use_protein_cropping} "
              f"(radius={self.pocket_crop_radius}\u00c5 + buffer={self.pocket_crop_buffer}\u00c5)")
        print("=" * 80)


# Singleton — import this everywhere
CFG = Config()
