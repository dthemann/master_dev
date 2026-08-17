#!/usr/bin/env python3
"""Whole-protein Uni-Dock (GPU) via overlapping 0.375 Å sub-box tiling.

Uni-Dock is a GPU re-implementation of AutoDock Vina, but its GPU affinity grid
is capped at 128 points/axis AND 531 441 total points (81^3). The whole-protein
search boxes in this benchmark are far too large to run at AutoDock's 0.375 Å
spacing in a single box. Instead we tile each whole-protein box with OVERLAPPING
sub-boxes, each <= ``subbox_size`` Å (29 Å -> 79^3 = 493 039 grid points at
0.375 Å, safely under the caps), dock the ligand in every sub-box on the GPU
(paired-batch only when it honours the requested settings; otherwise standard
per-tile commands), then merge the poses. This reproduces AutoDock's 0.375 Å
grid resolution across the whole protein.

Design points (each validated; see git history / the design workflow):

* Overlap = ligand DIAMETER + slack, not the compact start-conformer footprint.
  Docking freely reorients (and can extend) the ligand, so the longest atom-atom
  distance can align to any single box axis; the neighbour overlap must cover it
  or a seam-straddling pose is missed. (Review finding M1.)
* Cross-sub-box affinities are directly comparable: a fully-contained pose scores
  identically regardless of which box wall is nearby, because Uni-Dock's grid
  sums receptor atoms within the interaction cutoff of each grid point, including
  atoms outside the box. Verified empirically (score invariant to box size), so
  pooling poses and ranking by affinity is valid. (Review finding M5.)
* Empty sub-boxes (little/no receptor within cutoff) are skipped before docking.
  (Review finding M3.)
* Runs are resumable at the sub-box level within a provenance-fingerprinted,
  generation-specific scratch directory. A crash/OOM is fixed by re-running the
  cell, which re-docks only missing or strictly invalid sub-boxes. A complex is
  only marked ``success`` once every expected scored sub-box is present and the
  merged output has been atomically published. (Review finding M2.)
* Exhaustiveness barely affects runtime (per-box time is grid-precompute-bound),
  so it is set low (8) by default: running exhaustiveness 32 in each of dozens to
  thousands of sub-boxes would be a huge, heterogeneous over-search versus
  single-box AutoDock and would misrepresent search-effort parity. Total search
  effort (Sum of exhaustiveness over sub-boxes) is recorded per complex.
  (Review finding M4.)

Inputs are REUSED from the AutoDock Vina full-protein cell (no re-prep -> the
same receptor/ligand PDBQT and whole-protein box), so Uni-Dock docks exactly the
same complexes. All parameters live in a YAML config; see
``unidock_docking_config_full_protein.yaml``.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
import uuid
from pathlib import Path

import numpy as np
import yaml

# AutoDock atom types that are hydrogens (everything else is a heavy atom).
_H_TYPES = {"H", "HD"}
_VINA_RE = re.compile(
    r"^\s*REMARK\s+VINA RESULT:\s*"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?)"
)

# Uni-Dock 1.2.0's paired worker does not consume the corresponding CLI
# options: these values are compiled into vina_cuda_worker.h.  Do not add a
# version here without checking that version's source/behaviour first.
_PAIRED_FIXED_SETTINGS = {
    "1.2.0": {"scoring": "vina", "spacing": 0.375, "energy_range": 3.0},
}
_GRID_AXIS_CAP = 128
_GRID_TOTAL_CAP = 531_441
_MANIFEST_SCHEMA = 2


# ────────────────────────────────────────────────────────────────────────────
# Small parsers
# ────────────────────────────────────────────────────────────────────────────
def read_box(box_file: Path):
    """Parse an AutoDock '<key> = <value>' box.txt -> (center, size) 3-tuples."""
    v = {}
    for line in Path(box_file).read_text().splitlines():
        if "=" in line:
            k, val = line.split("=", 1)
            v[k.strip()] = float(val.strip())
    required = ("center_x", "center_y", "center_z", "size_x", "size_y", "size_z")
    missing = [key for key in required if key not in v]
    if missing:
        raise ValueError(f"box file {box_file} is missing: {', '.join(missing)}")
    values = [v[key] for key in required]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"box file {box_file} contains a non-finite value")
    if any(v[key] <= 0 for key in ("size_x", "size_y", "size_z")):
        raise ValueError(f"box file {box_file} contains a non-positive size")
    return (v["center_x"], v["center_y"], v["center_z"]), (v["size_x"], v["size_y"], v["size_z"])


def _atom_type(line: str) -> str:
    return line[77:].strip() or line.split()[-1]


def heavy_coords(pdbqt: Path) -> np.ndarray:
    """Heavy-atom coordinates (N,3) from a PDBQT (ATOM/HETATM, excluding H/HD)."""
    pts = []
    for line in Path(pdbqt).read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and _atom_type(line) not in _H_TYPES:
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    return np.asarray(pts, dtype=np.float64)


def ligand_diameter(pdbqt: Path) -> float:
    """Max pairwise heavy-atom distance — a rotation-safe extent bound (M1)."""
    xyz = heavy_coords(pdbqt)
    if len(xyz) < 2:
        return 0.0
    # O(n^2) but ligands are small; avoids a scipy dependency.
    d2 = np.sum((xyz[:, None, :] - xyz[None, :, :]) ** 2, axis=-1)
    return float(math.sqrt(d2.max()))


# ────────────────────────────────────────────────────────────────────────────
# Tiling geometry (overlap = ligand diameter + slack)
# ────────────────────────────────────────────────────────────────────────────
def _axis_tiles(c: float, L: float, W: float, overlap: float):
    """Tile one axis of length L (centre c) into windows of width <= W that
    overlap by >= ``overlap`` after shrinking by the ligand, so the
    fully-contained regions cover [c-L/2, c+L/2] with no seam gap."""
    usable = W - overlap
    lo = c - L / 2.0
    if L <= usable:
        return [(c, min(W, L + overlap))]
    n = int(math.ceil(L / usable))
    first = lo + usable / 2.0
    last = lo + L - usable / 2.0
    step = (last - first) / (n - 1)
    return [(first + i * step, W) for i in range(n)]


def subboxes(center, size, overlap: float, W: float = 29.0):
    """Overlapping sub-boxes covering a whole-protein box. Returns list of
    (cx, cy, cz, sx, sy, sz). Raises ValueError if the ligand is too large to be
    fully contained in any W-wide box (overlap >= W)."""
    if W - overlap <= 0:
        raise ValueError(f"overlap {overlap:.2f} >= subbox_size {W:.2f}: ligand too "
                         f"large to tile at this spacing")
    xs = _axis_tiles(center[0], size[0], W, overlap)
    ys = _axis_tiles(center[1], size[1], W, overlap)
    zs = _axis_tiles(center[2], size[2], W, overlap)
    return [(cx, cy, cz, sx, sy, sz) for cx, sx in xs for cy, sy in ys for cz, sz in zs]


def prefilter_boxes(boxes, rec_xyz: np.ndarray, cutoff: float, min_atoms: int):
    """Drop sub-boxes with <= min_atoms receptor heavy atoms within box+cutoff (M3)."""
    if len(rec_xyz) == 0:
        return list(boxes)
    kept = []
    for b in boxes:
        cx, cy, cz, sx, sy, sz = b
        lo = np.array([cx - sx / 2 - cutoff, cy - sy / 2 - cutoff, cz - sz / 2 - cutoff])
        hi = np.array([cx + sx / 2 + cutoff, cy + sy / 2 + cutoff, cz + sz / 2 + cutoff])
        inside = int(np.all((rec_xyz >= lo) & (rec_xyz <= hi), axis=1).sum())
        if inside > min_atoms:
            kept.append(b)
    return kept


# ────────────────────────────────────────────────────────────────────────────
# Pose parsing / merge  (numpy-only; same-atom-order RMSD)
# ────────────────────────────────────────────────────────────────────────────
class _Pose:
    __slots__ = ("affinity", "lines", "heavy_xyz")

    def __init__(self, affinity, lines, heavy_xyz):
        self.affinity = affinity
        self.lines = lines
        self.heavy_xyz = heavy_xyz


def parse_out_pdbqt(path: Path):
    """Strictly parse a complete Uni-Dock/Vina multi-model PDBQT.

    Every MODEL must have a matching ENDMDL, exactly one finite Vina affinity,
    and at least one parseable heavy atom.  Any malformed block invalidates the
    *whole* tile so a completed first model cannot hide a truncated later one.
    """
    poses, cur, affinities, heavy = [], None, None, None
    atom_records = 0
    try:
        lines = Path(path).read_text(errors="strict").splitlines()
    except (OSError, UnicodeError) as exc:
        raise ValueError(f"cannot read output {path}: {exc}") from exc

    for lineno, raw in enumerate(lines, start=1):
        record = raw[:6].strip()
        if record == "MODEL":
            if cur is not None:
                raise ValueError(f"{path}:{lineno}: nested MODEL")
            cur, affinities, heavy, atom_records = [], [], [], 0
            continue
        if record == "ENDMDL":
            if cur is None:
                raise ValueError(f"{path}:{lineno}: ENDMDL without MODEL")
            if len(affinities) != 1:
                raise ValueError(
                    f"{path}:{lineno}: MODEL has {len(affinities)} affinity records; expected 1"
                )
            if not math.isfinite(affinities[0]):
                raise ValueError(f"{path}:{lineno}: non-finite affinity")
            if atom_records == 0 or not heavy:
                raise ValueError(f"{path}:{lineno}: MODEL has no parseable heavy atoms")
            poses.append(_Pose(affinities[0], cur, np.asarray(heavy, dtype=np.float64)))
            cur = affinities = heavy = None
            continue
        if cur is None:
            continue

        cur.append(raw)
        match = _VINA_RE.match(raw)
        if match:
            try:
                affinities.append(float(match.group(1)))
            except ValueError as exc:  # defensive: the regex already constrains it
                raise ValueError(f"{path}:{lineno}: invalid affinity") from exc
        if record in ("ATOM", "HETATM"):
            atom_records += 1
            if _atom_type(raw) not in _H_TYPES:
                try:
                    xyz = (float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))
                except ValueError as exc:
                    raise ValueError(f"{path}:{lineno}: malformed atom coordinates") from exc
                if not all(math.isfinite(value) for value in xyz):
                    raise ValueError(f"{path}:{lineno}: non-finite atom coordinates")
                heavy.append(xyz)

    if cur is not None:
        raise ValueError(f"{path}: truncated MODEL without ENDMDL")
    if not poses:
        raise ValueError(f"{path}: no complete scored MODEL blocks")
    heavy_counts = {len(p.heavy_xyz) for p in poses}
    if len(heavy_counts) != 1:
        raise ValueError(f"{path}: inconsistent ligand heavy-atom counts across MODEL blocks")
    return poses


def _rmsd(a: np.ndarray, b: np.ndarray) -> float:
    d = a - b
    return math.sqrt(float(np.mean(np.sum(d * d, axis=1))))


def merge_poses(pose_lists, num_modes: int, energy_range: float, rmsd_threshold: float = 2.0):
    """Pool poses from all sub-boxes, drop overlap re-finds (best-first, RMSD dedup),
    keep those within ``energy_range`` of the global best, cap at ``num_modes``."""
    pool = [p for sub in pose_lists for p in sub if p.affinity is not None]  # m8 guard
    pool.sort(key=lambda p: p.affinity)
    if not pool:
        return []
    best = pool[0].affinity
    kept = []
    for p in pool:
        if p.affinity > best + energy_range:      # m7: AutoDock energy-range contract
            break
        if any(p.heavy_xyz.shape == q.heavy_xyz.shape
               and _rmsd(p.heavy_xyz, q.heavy_xyz) < rmsd_threshold for q in kept):
            continue
        kept.append(p)
        if len(kept) >= num_modes:
            break
    return kept


def write_merged_pdbqt(poses, out_path: Path):
    """Write kept poses as one multi-MODEL PDBQT, MODEL renumbered 1..k."""
    chunks = []
    for i, p in enumerate(poses, start=1):
        chunks.extend((f"MODEL {i}\n", "\n".join(p.lines) + "\n", "ENDMDL\n"))
    _atomic_write_text(out_path, "".join(chunks))


# ────────────────────────────────────────────────────────────────────────────
# Atomic state, configuration, and execution-mode helpers
# ────────────────────────────────────────────────────────────────────────────
def _atomic_write_text(path: Path, text: str):
    """Replace ``path`` atomically with UTF-8 text (same-filesystem rename)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex}.tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _atomic_write_json(path: Path, value):
    _atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True, default=str) + "\n")


def _atomic_copy(src: Path, dst: Path):
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f".{dst.name}.{uuid.uuid4().hex}.tmp")
    try:
        shutil.copyfile(src, tmp)
        os.replace(tmp, dst)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path):
    try:
        value = json.loads(Path(path).read_text())
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return value if isinstance(value, dict) else None


def _grid_points(size: float, spacing: float) -> int:
    """Conservative grid-point count used for Uni-Dock GPU cap validation."""
    return int(math.ceil(size / spacing)) + 1


def validate_config(cfg: dict) -> dict:
    """Fail early on unsafe or scientifically ambiguous configuration."""
    if not isinstance(cfg, dict):
        raise ValueError("configuration must be a YAML mapping")
    required = (
        "unidock_bin", "prep_base", "converter", "output_dir", "log_dir",
        "scoring", "subbox_size", "spacing", "overlap_slack",
        "prefilter_cutoff", "prefilter_min_atoms", "exhaustiveness",
        "num_modes", "energy_range", "seed", "paired_batch_size",
    )
    missing = [key for key in required if key not in cfg]
    if missing:
        raise ValueError(f"missing configuration keys: {', '.join(missing)}")

    if cfg["scoring"] not in {"vina", "vinardo", "ad4"}:
        raise ValueError("scoring must be one of: vina, vinardo, ad4")
    mode = cfg.get("execution_mode", "auto")
    if mode not in {"auto", "paired", "standard"}:
        raise ValueError("execution_mode must be auto, paired, or standard")

    spacing = float(cfg["spacing"])
    width = float(cfg["subbox_size"])
    slack = float(cfg["overlap_slack"])
    if not math.isfinite(spacing) or spacing <= 0:
        raise ValueError("spacing must be a finite positive number")
    if not math.isfinite(width) or width <= 0:
        raise ValueError("subbox_size must be a finite positive number")
    if not math.isfinite(slack) or slack < 0:
        raise ValueError("overlap_slack must be a finite non-negative number")
    if slack >= width:
        raise ValueError("overlap_slack must be smaller than subbox_size")

    npoints = _grid_points(width, spacing)
    if npoints > _GRID_AXIS_CAP or npoints ** 3 > _GRID_TOTAL_CAP:
        raise ValueError(
            f"subbox_size={width:g} at spacing={spacing:g} requires approximately "
            f"{npoints}^3 grid points, exceeding Uni-Dock caps "
            f"({_GRID_AXIS_CAP}/axis and {_GRID_TOTAL_CAP} total)"
        )

    numeric_nonnegative = ("prefilter_cutoff", "energy_range")
    for key in numeric_nonnegative:
        value = float(cfg[key])
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{key} must be a finite non-negative number")
    for key in ("prefilter_min_atoms",):
        if int(cfg[key]) < 0:
            raise ValueError(f"{key} must be non-negative")
    for key in ("exhaustiveness", "num_modes", "paired_batch_size"):
        if int(cfg[key]) < 1:
            raise ValueError(f"{key} must be at least 1")
    if cfg.get("max_receptor_atoms") is not None and int(cfg["max_receptor_atoms"]) < 1:
        raise ValueError("max_receptor_atoms must be null or at least 1")
    if cfg.get("exclude") is not None and not isinstance(cfg["exclude"], list):
        raise ValueError("exclude must be a YAML list")
    return cfg


def probe_unidock_version(unidock_bin: str) -> str:
    proc = subprocess.run([unidock_bin, "--version"], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"could not determine Uni-Dock version (rc={proc.returncode}): "
            f"{proc.stderr.strip() or proc.stdout.strip()}"
        )
    match = re.search(r"Uni-Dock\s+v?([^\s]+)", f"{proc.stdout}\n{proc.stderr}")
    if not match:
        raise RuntimeError("could not parse Uni-Dock --version output")
    return match.group(1)


def resolve_execution_mode(cfg: dict, version: str) -> str:
    """Choose a mode that really honours the requested scoring parameters.

    Uni-Dock 1.2.0 paired mode silently hard-codes three settings.  ``auto``
    therefore selects standard per-tile invocations when any requested value
    differs.  Explicit ``paired`` fails closed instead of silently changing the
    experiment.  Unknown versions are treated as unverified and use standard.
    """
    requested = cfg.get("execution_mode", "auto")
    fixed = _PAIRED_FIXED_SETTINGS.get(version)
    compatible = bool(
        fixed
        and cfg["scoring"] == fixed["scoring"]
        and math.isclose(float(cfg["spacing"]), fixed["spacing"], rel_tol=0, abs_tol=1e-12)
        and math.isclose(float(cfg["energy_range"]), fixed["energy_range"], rel_tol=0, abs_tol=1e-12)
    )
    if requested == "standard":
        return "standard"
    if requested == "paired":
        if not compatible:
            known = fixed or {"scoring": "unverified", "spacing": "unverified", "energy_range": "unverified"}
            raise ValueError(
                f"paired mode in Uni-Dock {version} cannot be verified to honour the requested "
                f"scoring/spacing/energy_range; known paired values are {known}. "
                "Use execution_mode: auto or standard."
            )
        return "paired"
    return "paired" if compatible else "standard"


# ────────────────────────────────────────────────────────────────────────────
# Per-complex docking (transactional and resumable)
# ────────────────────────────────────────────────────────────────────────────
def _valid_out(path: Path) -> bool:
    """A tile is complete only when strict scored-model parsing succeeds."""
    try:
        parse_out_pdbqt(path)
    except (OSError, ValueError):
        return False
    return True


def plan_complex(cdir: Path, cfg: dict):
    """Resolve inputs and compute the (prefiltered) sub-box list for one complex.
    Returns (info, boxes) or (info, None) if it should be skipped."""
    pdb_id = cdir.name
    conv = cfg["converter"]
    rec = next(iter((cdir / "_staging" / "receptors" / "pdbqt").glob(f"*_{conv}.pdbqt")), None)
    box_file = next(iter((cdir / "_staging" / "receptors" / "pdbqt").glob(f"*_{conv}.box.txt")), None)
    lig = next(iter((cdir / "_staging" / "ligands" / "pdbqt").glob("*_ligand_start_conf.pdbqt")), None)
    info = {"pdb_id": pdb_id, "receptor": str(rec) if rec else None,
            "ligand": str(lig) if lig else None,
            "box_file": str(box_file) if box_file else None}

    if not (rec and box_file and lig):
        info["skip"] = "missing prepared inputs (run the AutoDock Vina cell first)"
        return info, None

    rec_xyz = heavy_coords(rec)
    info["receptor_atoms"] = len(rec_xyz)
    if len(rec_xyz) == 0:
        info["skip"] = "receptor has no parseable heavy atoms"
        return info, None
    mra = cfg.get("max_receptor_atoms")
    if mra is not None and len(rec_xyz) > int(mra):
        info["skip"] = f"receptor {len(rec_xyz)} atoms > max_receptor_atoms {mra}"
        return info, None

    ligand_xyz = heavy_coords(lig)
    if len(ligand_xyz) == 0:
        info["skip"] = "ligand has no parseable heavy atoms"
        return info, None
    center, size = read_box(box_file)
    diam = ligand_diameter(lig)
    overlap = diam + float(cfg["overlap_slack"])
    W = float(cfg["subbox_size"])
    info["ligand_diameter"] = round(diam, 2)
    try:
        boxes = subboxes(center, size, overlap, W)
    except ValueError as e:
        info["skip"] = f"cannot tile: {e}"
        return info, None
    boxes = prefilter_boxes(boxes, rec_xyz, float(cfg["prefilter_cutoff"]),
                            int(cfg["prefilter_min_atoms"]))
    if not boxes:
        info["skip"] = "no receptor-containing sub-boxes after prefilter"
        return info, None
    info["n_subboxes"] = len(boxes)
    return info, boxes


_PROVENANCE_CONFIG_KEYS = (
    "converter", "scoring", "spacing", "subbox_size", "overlap_slack",
    "prefilter_cutoff", "prefilter_min_atoms", "max_receptor_atoms",
    "exhaustiveness", "num_modes", "energy_range", "seed",
    "paired_batch_size", "execution_mode",
)


def _file_provenance(path: Path):
    path = Path(path).resolve()
    return {"path": str(path), "size": path.stat().st_size, "sha256": _sha256(path)}


def _build_provenance(info: dict, boxes, cfg: dict, unidock_bin: str,
                      unidock_version: str, execution_mode: str):
    settings = {key: cfg.get(key) for key in _PROVENANCE_CONFIG_KEYS}
    settings["execution_mode_resolved"] = execution_mode
    payload = {
        "schema_version": _MANIFEST_SCHEMA,
        "engine": {"path": str(Path(unidock_bin).resolve()), "version": unidock_version},
        "inputs": {
            "receptor": _file_provenance(Path(info["receptor"])),
            "ligand": _file_provenance(Path(info["ligand"])),
            "box": _file_provenance(Path(info["box_file"])),
        },
        "settings": settings,
        "subboxes": [[float(value) for value in box] for box in boxes],
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return payload, hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _output_matches_manifest(out_pose: Path, manifest: dict) -> bool:
    expected = manifest.get("output_sha256")
    expected_poses = manifest.get("num_poses")
    if not expected or isinstance(expected_poses, bool) or not isinstance(expected_poses, int):
        return False
    try:
        poses = parse_out_pdbqt(out_pose)
        return bool(poses) and len(poses) == expected_poses and _sha256(out_pose) == expected
    except (OSError, ValueError):
        return False


def _done_marker_matches(done_marker: Path, manifest: dict) -> bool:
    done = _read_json(done_marker)
    return bool(
        done
        and done.get("schema_version") == _MANIFEST_SCHEMA
        and done.get("engine") == "unidock-tiled"
        and done.get("status") == "complete"
        and done.get("commit_status") == manifest.get("commit_status") == "complete"
        and done.get("fingerprint") == manifest.get("fingerprint")
        and done.get("generation_id") == manifest.get("generation_id")
        and done.get("output_file") == manifest.get("output_file")
        and done.get("output_sha256") == manifest.get("output_sha256")
        and done.get("num_poses") == manifest.get("num_poses")
        and done.get("n_subboxes") == manifest.get("n_subboxes")
    )


def _new_manifest(provenance: dict, fingerprint: str, overwrite: bool):
    return {
        "schema_version": _MANIFEST_SCHEMA,
        "fingerprint": fingerprint,
        "generation_id": uuid.uuid4().hex,
        "status": "in_progress",
        "started_with_overwrite": overwrite,
        "created_unix_s": time.time(),
        "provenance": provenance,
        "attempts": [],
    }


def _prepare_tile_inputs(tiles_dir: Path, pdb_id: str, ligand: Path, boxes, indices):
    paired_cfg = {}
    for i in indices:
        cx, cy, cz, sx, sy, sz = boxes[i]
        lig_copy = tiles_dir / f"{pdb_id}_s{i}.pdbqt"
        box_cfg = tiles_dir / f"{pdb_id}_s{i}.box.txt"
        # Always refresh copies.  A prior interrupted attempt must not preserve
        # an obsolete or partially-written ligand/config file.
        _atomic_copy(ligand, lig_copy)
        _atomic_write_text(
            box_cfg,
            f"center_x = {cx}\ncenter_y = {cy}\ncenter_z = {cz}\n"
            f"size_x = {sx}\nsize_y = {sy}\nsize_z = {sz}\n",
        )
        paired_cfg[f"s{i}"] = {
            "protein": None,  # filled by the paired runner
            "ligand": str(lig_copy.resolve()),
            "ligand_config": str(box_cfg.resolve()),
        }
    return paired_cfg


def _run_paired(unidock_bin: str, rec: Path, tiles_dir: Path, pdb_id: str,
                boxes, missing, cfg: dict):
    paired_cfg = _prepare_tile_inputs(tiles_dir, pdb_id, Path(cfg["_ligand"]), boxes, missing)
    for value in paired_cfg.values():
        value["protein"] = str(rec)
    paired = tiles_dir / "paired.json"
    _atomic_write_json(paired, paired_cfg)
    cmd = [
        unidock_bin, "--paired_batch_size", str(int(cfg["paired_batch_size"])),
        "--ligand_index", str(paired.resolve()), "--dir", str(tiles_dir.resolve()),
        "--scoring", cfg["scoring"], "--spacing", str(cfg["spacing"]),
        "--size_x", str(cfg["subbox_size"]), "--size_y", str(cfg["subbox_size"]),
        "--size_z", str(cfg["subbox_size"]),
        "--exhaustiveness", str(int(cfg["exhaustiveness"])),
        "--num_modes", str(int(cfg["num_modes"])),
        "--energy_range", str(cfg["energy_range"]), "--seed", str(int(cfg["seed"])),
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        return [proc.returncode], [
            f"$ {' '.join(cmd)}\nrc={proc.returncode}\n"
            f"=== STDOUT ===\n{proc.stdout}\n=== STDERR ===\n{proc.stderr}\n"
        ]
    except OSError as exc:
        return [None], [f"$ {' '.join(cmd)}\nlaunch error: {exc}\n"]


def _run_standard(unidock_bin: str, rec: Path, tiles_dir: Path, pdb_id: str,
                  boxes, missing, cfg: dict):
    _prepare_tile_inputs(tiles_dir, pdb_id, Path(cfg["_ligand"]), boxes, missing)
    returncodes, logs = [], []
    for i in missing:
        cx, cy, cz, sx, sy, sz = boxes[i]
        lig_copy = tiles_dir / f"{pdb_id}_s{i}.pdbqt"
        final_out = tiles_dir / f"{pdb_id}_s{i}_out.pdbqt"
        tmp_out = tiles_dir / f".{pdb_id}_s{i}_out.{uuid.uuid4().hex}.tmp.pdbqt"
        cmd = [
            unidock_bin, "--receptor", str(rec), "--ligand", str(lig_copy.resolve()),
            "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
            "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
            "--out", str(tmp_out.resolve()), "--scoring", cfg["scoring"],
            "--spacing", str(cfg["spacing"]),
            "--exhaustiveness", str(int(cfg["exhaustiveness"])),
            "--num_modes", str(int(cfg["num_modes"])),
            "--energy_range", str(cfg["energy_range"]), "--seed", str(int(cfg["seed"])),
        ]
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
            rc, stdout, stderr = proc.returncode, proc.stdout, proc.stderr
        except OSError as exc:
            rc, stdout, stderr = None, "", f"launch error: {exc}"
        returncodes.append(rc)
        logs.append(
            f"$ {' '.join(cmd)}\nrc={rc}\n=== STDOUT ===\n{stdout}\n=== STDERR ===\n{stderr}\n"
        )
        if rc == 0 and _valid_out(tmp_out):
            os.replace(tmp_out, final_out)
        else:
            try:
                tmp_out.unlink()
            except FileNotFoundError:
                pass
    return returncodes, logs


def dock_complex(cdir: Path, cfg: dict, unidock_bin: str,
                 unidock_version: str | None = None, execution_mode: str | None = None):
    """Tile, dock, and transactionally publish one complex.

    A generation-specific scratch directory prevents legacy or previous
    overwrite tiles from ever satisfying completeness.  The merged output and
    done sentinel are published only after every strict tile validation passes.
    """
    t0 = time.time()
    validate_config(cfg)
    version = unidock_version or probe_unidock_version(unidock_bin)
    mode = execution_mode or resolve_execution_mode(cfg, version)
    if mode not in {"paired", "standard"}:
        raise ValueError(f"invalid resolved execution mode: {mode}")
    if mode == "paired":
        resolve_execution_mode({**cfg, "execution_mode": "paired"}, version)

    info, boxes = plan_complex(cdir, cfg)
    if boxes is None:
        return {**info, "status": "skipped", "elapsed_s": round(time.time() - t0, 1)}

    pdb_id = info["pdb_id"]
    rec = Path(info["receptor"]).resolve()
    lig = Path(info["ligand"]).resolve()
    out_dir = Path(cfg["output_dir"]) / pdb_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pose = out_dir / f"{pdb_id}_unidock_out.pdbqt"
    done_marker = out_dir / ".unidock_done"
    manifest_path = out_dir / "run_manifest.json"
    provenance, fingerprint = _build_provenance(info, boxes, cfg, unidock_bin, version, mode)
    current = _read_json(manifest_path)
    overwrite = bool(cfg.get("overwrite", False))

    same_run = bool(
        current
        and current.get("schema_version") == _MANIFEST_SCHEMA
        and current.get("fingerprint") == fingerprint
    )
    if (not overwrite and same_run and current.get("status") == "complete"
            and current.get("commit_status") == "complete"
            and _output_matches_manifest(out_pose, current)):
        marker_matches = _done_marker_matches(done_marker, current)
        if not marker_matches:
            _atomic_write_json(done_marker, {
                "schema_version": _MANIFEST_SCHEMA,
                "engine": "unidock-tiled",
                "status": "complete", "commit_status": "complete",
                "fingerprint": fingerprint,
                "generation_id": current["generation_id"],
                "num_poses": current.get("num_poses"),
                "n_subboxes": current.get("n_subboxes"),
                "output_file": current.get("output_file"),
                "output_sha256": current.get("output_sha256"),
            })
        return {
            **info, "status": "done", "fingerprint": fingerprint,
            "generation_id": current["generation_id"], "execution_mode": mode,
            "commit_status": "complete", "output_file": current["output_file"],
            "output_sha256": current["output_sha256"],
            "num_poses": current["num_poses"],
            "best_affinity": current.get("best_affinity"),
            "n_subboxes": current["n_subboxes"], "n_subboxes_missing": 0,
            "unidock_version": version, "elapsed_s": round(time.time() - t0, 1),
        }

    resumable = bool(
        same_run
        and current.get("status") in {"in_progress", "incomplete"}
        and (not overwrite or current.get("started_with_overwrite") is True)
        and current.get("generation_id")
    )
    if resumable:
        manifest = current
    else:
        manifest = _new_manifest(provenance, fingerprint, overwrite)
        _atomic_write_json(manifest_path, manifest)
    # A marker from an old generation is never a valid completion signal.
    try:
        done_marker.unlink()
    except FileNotFoundError:
        pass

    generation_id = manifest["generation_id"]
    tiles_dir = out_dir / "_tiles" / generation_id
    tiles_dir.mkdir(parents=True, exist_ok=True)

    def tile_out(i):
        return tiles_dir / f"{pdb_id}_s{i}_out.pdbqt"

    missing = [i for i in range(len(boxes)) if not _valid_out(tile_out(i))]
    for i in missing:
        # Remove malformed current-generation files before trying again.
        try:
            tile_out(i).unlink()
        except FileNotFoundError:
            pass

    returncodes, log_chunks = [], []
    if missing:
        run_cfg = dict(cfg)
        run_cfg["_ligand"] = str(lig)
        if mode == "paired":
            returncodes, log_chunks = _run_paired(
                unidock_bin, rec, tiles_dir, pdb_id, boxes, missing, run_cfg
            )
        else:
            returncodes, log_chunks = _run_standard(
                unidock_bin, rec, tiles_dir, pdb_id, boxes, missing, run_cfg
            )
        _atomic_write_text(
            Path(cfg["log_dir"]) / f"{pdb_id}.log",
            f"generation={generation_id}\nfingerprint={fingerprint}\nmode={mode}\n"
            f"attempted {len(missing)} of {len(boxes)} sub-boxes\n\n" + "\n".join(log_chunks),
        )

    attempted_failed = any(rc != 0 for rc in returncodes)
    # Fail closed if an input changes while a long tiled run is in flight; the
    # next invocation will compute a new fingerprint and use a new generation.
    _, final_fingerprint = _build_provenance(info, boxes, cfg, unidock_bin, version, mode)
    inputs_changed_during_run = final_fingerprint != fingerprint
    attempted_failed = attempted_failed or inputs_changed_during_run
    present = [i for i in range(len(boxes)) if _valid_out(tile_out(i))]
    n_missing = len(boxes) - len(present)
    merged = merge_poses(
        (parse_out_pdbqt(tile_out(i)) for i in present),
        num_modes=int(cfg["num_modes"]), energy_range=float(cfg["energy_range"]),
    )
    status = "success" if n_missing == 0 and not attempted_failed and merged else "incomplete"
    best = merged[0].affinity if merged else None
    n_docked_now = sum(1 for i in missing if _valid_out(tile_out(i)))

    attempt = {
        "unix_s": time.time(), "execution_mode": mode, "attempted_tiles": len(missing),
        "valid_new_tiles": n_docked_now, "returncodes": returncodes,
        "remaining_tiles": n_missing,
        "inputs_changed_during_run": inputs_changed_during_run,
    }
    manifest.setdefault("attempts", []).append(attempt)
    manifest["attempts"] = manifest["attempts"][-50:]
    manifest["status"] = "complete" if status == "success" else "incomplete"
    manifest["updated_unix_s"] = time.time()

    if status == "success":
        write_merged_pdbqt(merged, out_pose)
        # Validate our just-published deliverable before committing completion.
        published = parse_out_pdbqt(out_pose)
        if len(published) != len(merged):
            raise RuntimeError("atomic merged-output validation changed the pose count")
        manifest["output_sha256"] = _sha256(out_pose)
        manifest["num_poses"] = len(merged)
        manifest["best_affinity"] = best
        manifest["n_subboxes"] = len(boxes)
        manifest["output_file"] = out_pose.name
        manifest["commit_status"] = "complete"
        _atomic_write_json(manifest_path, manifest)
        _atomic_write_json(done_marker, {
            "schema_version": _MANIFEST_SCHEMA, "engine": "unidock-tiled",
            "status": "complete", "commit_status": "complete",
            "fingerprint": fingerprint, "generation_id": generation_id,
            "num_poses": len(merged), "n_subboxes": len(boxes),
            "output_file": out_pose.name,
            "output_sha256": manifest["output_sha256"],
        })
        if not cfg.get("keep_tiles", False):
            shutil.rmtree(tiles_dir, ignore_errors=True)
    else:
        # Do not replace the public merged output with a partial result.
        _atomic_write_json(manifest_path, manifest)

    return {
        **info, "status": status, "num_poses": len(merged), "best_affinity": best,
        "n_subboxes": len(boxes), "n_subboxes_missing": n_missing,
        "n_subboxes_attempted": len(missing), "n_docked_this_run": n_docked_now,
        "process_returncodes": returncodes, "fingerprint": fingerprint,
        "inputs_changed_during_run": inputs_changed_during_run,
        "generation_id": generation_id,
        "total_search_effort": len(boxes) * int(cfg["exhaustiveness"]),
        "scoring": cfg["scoring"], "spacing": cfg["spacing"],
        "subbox_size": cfg["subbox_size"], "energy_range": cfg["energy_range"],
        "exhaustiveness": int(cfg["exhaustiveness"]), "engine": "unidock-tiled",
        "execution_mode": mode, "unidock_version": version,
        "commit_status": "complete" if status == "success" else None,
        "output_file": out_pose.name if status == "success" else None,
        "output_sha256": manifest.get("output_sha256") if status == "success" else None,
        "published_output_current": status == "success",
        "elapsed_s": round(time.time() - t0, 1),
    }


# ────────────────────────────────────────────────────────────────────────────
# Driver
# ────────────────────────────────────────────────────────────────────────────
def load_config(path) -> dict:
    with open(path) as f:
        return validate_config(yaml.safe_load(f))


def discover_complexes(cfg: dict):
    prep = Path(cfg["prep_base"])
    exclude = set(cfg.get("exclude") or [])
    return [d for d in sorted(prep.iterdir())
            if d.is_dir() and not d.name.startswith(("_", ".")) and d.name != "Logs"
            and d.name not in exclude]


def main(config_path, plan_only=False, limit=None):
    cfg = load_config(config_path)
    unidock_bin = cfg["unidock_bin"]
    if not Path(unidock_bin).is_file():
        raise FileNotFoundError(f"Uni-Dock binary not found: {unidock_bin}")
    unidock_version = probe_unidock_version(unidock_bin)
    execution_mode = resolve_execution_mode(cfg, unidock_version)
    Path(cfg["output_dir"]).mkdir(parents=True, exist_ok=True)
    Path(cfg["log_dir"]).mkdir(parents=True, exist_ok=True)

    complexes = discover_complexes(cfg)
    if limit:
        complexes = complexes[:limit]

    # ── Plan: tile every complex, report the workload + time estimate ────────
    RATE_S = 4.7   # measured s/sub-box, paired-batch, exhaustiveness-independent
    total_boxes, dockable, skipped, plan_errors = 0, [], [], []
    for cdir in complexes:
        try:                                    # a malformed box.txt/PDBQT must not halt planning
            info, boxes = plan_complex(cdir, cfg)
        except Exception as e:
            plan_errors.append((cdir.name, f"plan error: {e}"))
            continue
        if boxes is None:
            skipped.append((cdir.name, info.get("skip")))
        else:
            dockable.append(cdir)
            total_boxes += info["n_subboxes"]
    excluded = sorted(set(cfg.get("exclude") or []))
    print(f"Uni-Dock TILED whole-protein | spacing {cfg['spacing']} Å | subbox {cfg['subbox_size']} Å "
          f"| exhaustiveness {cfg['exhaustiveness']} | mode {execution_mode} "
          f"(Uni-Dock {unidock_version})")
    print(f"  complexes: {len(dockable)} to dock, {len(skipped)} skipped, {len(plan_errors)} errors, "
          f"{len(excluded)} excluded ({', '.join(excluded) or 'none'})")
    if execution_mode == "paired":
        print(f"  sub-boxes total: {total_boxes}  →  est. {total_boxes * RATE_S / 3600:.1f} h "
              f"(~{RATE_S:.1f}s/box measured in paired mode)")
    else:
        print(f"  sub-boxes total: {total_boxes}; no calibrated runtime estimate for "
              "standard per-tile mode (the paired 4.7 s/box measurement is not applicable)")
    if cfg.get("max_receptor_atoms"):
        print(f"  receptor-size cap: skip receptors > {cfg['max_receptor_atoms']} heavy atoms")
    for name, why in skipped[:20]:
        print(f"    skip {name}: {why}")
    for name, why in plan_errors[:20]:
        print(f"    error {name}: {why}")
    if plan_only:
        return {"dockable": len(dockable), "skipped": skipped, "errors": plan_errors,
                "total_subboxes": total_boxes, "execution_mode": execution_mode,
                "unidock_version": unidock_version}

    # ── Dock ─────────────────────────────────────────────────────────────────
    results = [
        {"pdb_id": name, "status": "skipped", "skip": why} for name, why in skipped
    ] + [
        {"pdb_id": name, "status": "error", "error": why} for name, why in plan_errors
    ]
    for idx, cdir in enumerate(dockable, 1):
        try:
            r = dock_complex(
                cdir, cfg, unidock_bin, unidock_version=unidock_version,
                execution_mode=execution_mode,
            )
        except Exception as e:                                  # never let one complex halt the run
            r = {"pdb_id": cdir.name, "status": "error", "error": str(e)}
        results.append(r)
        if r["status"] != "done":   # don't clobber a finished complex's summary on rerun
            cdir_out = Path(cfg["output_dir"]) / cdir.name
            _atomic_write_json(cdir_out / "docking_summary.json", r)
        icon = {"success": "✓", "incomplete": "◐", "skipped": "–",
                "error": "✗", "done": "≡"}.get(r["status"], "?")
        extra = ""
        if r["status"] in ("success", "incomplete"):
            extra = (f" | {r['num_poses']} poses, best {r['best_affinity']:.2f} "
                     f"| {r['n_subboxes']} boxes ({r['n_docked_this_run']} docked now"
                     f"{', %d MISSING' % r['n_subboxes_missing'] if r['n_subboxes_missing'] else ''})"
                     f" | {r['elapsed_s']:.0f}s") if r["best_affinity"] is not None else \
                    f" | 0 poses | {r['n_subboxes']} boxes | {r['elapsed_s']:.0f}s"
        print(f"[{idx}/{len(dockable)}] {icon} {cdir.name}{extra}")

    ok = sum(1 for r in results if r["status"] == "success")
    inc = sum(1 for r in results if r["status"] == "incomplete")
    err = sum(1 for r in results if r["status"] == "error")
    done = sum(1 for r in results if r["status"] == "done")
    skip = sum(1 for r in results if r["status"] == "skipped")
    print(f"\n{'='*60}\nUNI-DOCK TILED COMPLETE\n  success: {ok}  already-done: {done}  "
          f"incomplete: {inc}  skipped: {skip}  error: {err}")
    if inc:
        print("  (re-run the cell to finish 'incomplete' complexes — only missing sub-boxes re-dock)")
    print(f"  results: {cfg['output_dir']}\n{'='*60}")
    return results


def results_exit_code(results) -> int:
    """CLI failure if any selected complex errored or remains incomplete."""
    if isinstance(results, list) and any(
        result.get("status") in {"incomplete", "error"} for result in results
    ):
        return 1
    return 0


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Tiled whole-protein Uni-Dock docking")
    ap.add_argument("-c", "--config", required=True)
    ap.add_argument("--plan-only", action="store_true", help="print the workload + time estimate, do not dock")
    ap.add_argument("--limit", type=int, default=None, help="dock only the first N complexes (testing)")
    a = ap.parse_args()
    outcome = main(a.config, plan_only=a.plan_only, limit=a.limit)
    raise SystemExit(0 if a.plan_only else results_exit_code(outcome))
