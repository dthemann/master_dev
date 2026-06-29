"""Post-pose docking re-search (force-field / CNN re-optimisation).

EquiBind regresses ligand coordinates with no excluded-volume term, so its
poses routinely interpenetrate the protein (median ~45% volume overlap in the
benchmark). UFF with the protein held fixed is a *local* minimiser and cannot
translate a deeply buried ligand back out. A docking re-search can: it samples
rigid-body + torsional moves of the ligand against the receptor and relaxes it
out of the clash while staying in the same pocket.

This module implements an in-place local optimisation with either **smina** or
**gnina** (``--local_only`` or ``--minimize``) seeded with the EquiBind pose and
boxed around it (``--autobox_ligand``). Both share the AutoDock-Vina CLI, so a
single wrapper drives either; gnina additionally runs a CNN rescorer and writes
``<CNNscore>`` / ``<CNNaffinity>`` on top of smina's ``<minimizedAffinity>``.
Each reads/writes SDF directly, so no PDBQT round-trip is needed for the ligand.
"""
from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Dict, Optional, Tuple

from .config import CFG
from .monitor import monitor

# Supported re-search backends. Both speak the AutoDock-Vina CLI.
REFINE_TOOLS = ("smina", "gnina")

# Data tags smina/gnina write into the output SDF, mapped to our score keys.
_SDF_TAG_RE = re.compile(r">\s*<([^>]+)>\s*\n([^\n]*)")
_SCORE_TAGS = {
    "minimizedAffinity": "minimized_affinity",  # smina + gnina (kcal/mol)
    "CNNscore": "cnn_score",                     # gnina only (0..1 pose quality)
    "CNNaffinity": "cnn_affinity",               # gnina only (predicted pK)
}
# Fallback: parse the affinity from the stdout result table ("   1   -7.3 ...").
_STDOUT_RE = re.compile(r"^\s*1\s+(-?\d+\.\d+)", re.MULTILINE)


def _resolve_receptor(protein_pdb: Path) -> Path:
    """Pick the receptor file the re-search should dock against.

    Defaults to the prepared protein PDB. If ``EQ_SMINA_RECEPTOR_EXT`` asks for
    a different extension (e.g. a ``.pdbqt`` sibling produced by meeko/mgltools)
    and that file exists next to the PDB, use it; otherwise fall back to the PDB.
    smina and gnina both read PDB/PDBQT, so this applies to either backend.
    """
    ext = CFG.smina_receptor_ext
    if ext and ext != protein_pdb.suffix.lower():
        sibling = protein_pdb.with_suffix(ext)
        if sibling.exists():
            return sibling
        monitor.warning(
            f"refine receptor {sibling.name} not found; using {protein_pdb.name}")
    return protein_pdb


def _refine_executable(tool: str) -> str:
    """Resolve the configured executable for ``tool`` (falls back to its name)."""
    exe = CFG.gnina_executable if tool == "gnina" else CFG.smina_executable
    return exe or tool


def _parse_scores(out_sdf: Path, stdout: str) -> Dict[str, float]:
    """Extract minimizedAffinity (+ gnina CNN scores) from the output SDF."""
    scores: Dict[str, float] = {}
    try:
        if out_sdf.exists():
            for m in _SDF_TAG_RE.finditer(out_sdf.read_text()):
                key = _SCORE_TAGS.get(m.group(1).strip())
                if key and key not in scores:
                    try:
                        scores[key] = float(m.group(2).strip())
                    except ValueError:
                        pass
    except Exception:
        pass
    if "minimized_affinity" not in scores:
        m = _STDOUT_RE.search(stdout or "")
        if m:
            scores["minimized_affinity"] = float(m.group(1))
    return scores


def refine_pose(
    pose_sdf: Path,
    protein_pdb: Path,
    output_sdf: Path,
    tool: str = "smina",
    search: Optional[str] = None,
    autobox_add: Optional[float] = None,
    cpu: Optional[int] = None,
    seed: Optional[int] = None,
    timeout_s: Optional[int] = None,
) -> Tuple[bool, Dict[str, float], str]:
    """Locally re-optimise ``pose_sdf`` against the receptor with smina/gnina.

    Returns ``(success, scores, message)`` where ``scores`` maps
    ``minimized_affinity`` (+ ``cnn_score`` / ``cnn_affinity`` for gnina) to
    floats. ``success`` means the tool produced a non-empty output SDF; on
    failure the caller should fall back to the input pose so it is never lost.
    """
    tool = (tool or "smina").strip().lower()
    if tool not in REFINE_TOOLS:
        return False, {}, f"unknown refine tool {tool!r} (expected one of {REFINE_TOOLS})"

    search = (search or CFG.smina_search)
    autobox_add = CFG.smina_autobox_add if autobox_add is None else autobox_add
    cpu = CFG.smina_cpu if cpu is None else cpu
    seed = CFG.smina_seed if seed is None else seed
    timeout_s = CFG.smina_timeout_s if timeout_s is None else timeout_s

    exe = _refine_executable(tool)
    receptor = _resolve_receptor(protein_pdb)
    output_sdf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        exe,
        "--receptor", str(receptor),
        "--ligand", str(pose_sdf),
        "--autobox_ligand", str(pose_sdf),
        "--autobox_add", str(autobox_add),
        "--out", str(output_sdf),
        "--cpu", str(cpu),
        "--seed", str(seed),
        "--num_modes", "1",
    ]
    # --local_only / --minimize both keep the ligand near its input pose; they
    # never run a global search, so the refined pose stays in the same pocket.
    cmd.append("--minimize" if search == "minimize" else "--local_only")
    # gnina runs a CNN by default; keep it off the GPU unless asked (the GPU is
    # usually busy with EquiBind inference). smina has no --no_gpu flag.
    if tool == "gnina" and not CFG.gnina_use_gpu:
        cmd.append("--no_gpu")

    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired:
        return False, {}, f"{tool} timed out after {timeout_s}s"
    except FileNotFoundError:
        return False, {}, f"{tool} executable not found: {exe}"
    except Exception as e:  # pragma: no cover - defensive
        return False, {}, f"{tool} invocation failed: {e}"

    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        return False, {}, f"{tool} rc={proc.returncode}: {err[-1] if err else 'unknown error'}"

    if not output_sdf.exists() or output_sdf.stat().st_size == 0:
        return False, {}, f"{tool} produced no output pose"

    scores = _parse_scores(output_sdf, proc.stdout)
    aff = scores.get("minimized_affinity")
    aff_str = f"{aff:.2f} kcal/mol" if aff is not None else "n/a"
    cnn_str = ""
    if "cnn_score" in scores:
        cnn_str = f", CNNscore={scores['cnn_score']:.3f}"
    return True, scores, f"{tool} {search} ok (affinity={aff_str}{cnn_str})"


def smina_refine_pose(
    pose_sdf: Path,
    protein_pdb: Path,
    output_sdf: Path,
    **kwargs,
) -> Tuple[bool, Optional[float], str]:
    """Back-compat shim for the original smina-only entry point.

    Returns ``(success, minimized_affinity, message)`` — the legacy 3-tuple.
    New code should call :func:`refine_pose`, which also exposes gnina + the
    full score dict.
    """
    ok, scores, msg = refine_pose(pose_sdf, protein_pdb, output_sdf, tool="smina", **kwargs)
    return ok, scores.get("minimized_affinity"), msg
