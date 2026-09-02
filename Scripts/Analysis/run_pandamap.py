"""PandaMap interaction-fingerprint generator for docked poses.

A config-driven, parallel re-write of ``panda_maps.ipynb`` that mirrors the
structure of ``Scripts/Docking/Posebusters/run_posebusters.py``. For every
selected docked pose it builds a combined protein-ligand PDB, runs PandaMap
(``HybridProtLigMapper``) to detect the 16 protein-ligand interaction types, and
records the result at **residue level** (which residue, which interaction type,
which ligand atom, distance) — not just per-type counts. It also maps the
crystal ligand of each pair to give a *native* interaction fingerprint, so the
companion report (``pandamap_interaction_report.py``) can measure how faithfully
each docking method reproduces the native interactions.

Input modes (choose one in the config)
--------------------------------------
* **Mode A — PoseBusters CSV** (``pb_csv``): reuse the per-pose table written by
  run_posebusters. It already carries ``pose_file / pose_name / protein / ligand
  / docking_method / protein_file_used`` plus the variant provenance columns
  (``pocket_source / clamp_variant / refine_variant``). ``pb_valid_only`` keeps
  only poses that passed every canonical PoseBusters test.
* **Mode B — directory scan** (``docking_directories``): reuse run_posebusters'
  ``ROW_COLLECTORS`` / ``POSE_EXPANDERS`` / ``find_protein_file`` so the fragile
  path parsing is never duplicated.

Pose selection: **top-N per (method, protein, ligand)** (``poses_per_combo``).

Outputs (under ``output_dir``)
------------------------------
* ``pandamap_interactions.csv`` — long/residue-level: one row per
  (pose, interaction_type, residue).
* ``pandamap_pose_summary.csv`` — one row per pose: per-type counts + total.
* ``crystal_interactions.csv`` — native fingerprint per pair (residue-level).
* ``maps/<method>/<protein>__<ligand>/*.png`` — 2D interaction maps (``render_images``).
* ``pandamap_errors.csv``.

Must run under the **vina** conda env (the only one with ``pandamap``):
    /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py \
        --config Scripts/Analysis/pandamap_config.yaml
Helpful flags: ``--limit-pairs N`` (debug), ``--workers K``, ``--poses-per-combo N``,
``--render/--no-render``, ``--overwrite``, ``--out-dir DIR``.
"""
from __future__ import annotations

import argparse
import contextlib
import io
import os
import re
import sys
import tempfile
import warnings
from dataclasses import dataclass, field
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

warnings.filterwarnings("ignore")

# run_posebusters lives one package over — add it to the path so we can reuse its
# pose collectors / protein resolution / canonical test set rather than
# duplicating them (single source of truth, same as the other analysis scripts).
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_PB_DIR = _PROJECT_ROOT / "Scripts" / "Docking" / "Posebusters"
for _p in (str(_PROJECT_ROOT), str(_PB_DIR)):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# The 16 PandaMap interaction types (order = pandamap/core.py self.interactions).
INTERACTION_TYPES = [
    "hydrogen_bonds", "carbon_pi", "pi_pi_stacking", "donor_pi", "amide_pi",
    "hydrophobic", "ionic", "halogen_bonds", "cation_pi", "metal_coordination",
    "salt_bridge", "covalent", "alkyl_pi", "attractive_charge", "pi_cation",
    "repulsion",
]

# Sodium's element symbol collides with pandas' default missing-value tokens: it
# is the literal string "NA". Reading pandamap_interactions.csv back with the
# defaults parses those cells as NaN, and the resume path below writes the merged
# frame straight back out as empty fields, so every resumed run erased sodium
# from its own fingerprints — 275 metal-coordination rows on the run that exposed
# this, with nothing recomputed to replace them. Read the frames that carry
# lig_atom_element with the default token list switched off and only the empty
# field treated as missing: "NA" survives the round trip, while a field PandaMap
# genuinely could not fill still arrives as NaN rather than "".
_CSV_NA_VALUES = [""]


def read_pandamap_csv(path, **kwargs) -> pd.DataFrame:
    """``read_csv`` that keeps element symbols literal (see ``_CSV_NA_VALUES``)."""
    return pd.read_csv(path, keep_default_na=False, na_values=_CSV_NA_VALUES, **kwargs)


# ──────────────────────────────────────────────────────────────────────────
# Config
# ──────────────────────────────────────────────────────────────────────────

@dataclass
class PandaMapConfig:
    work_dir: Path
    output_dir: Path
    # Mode A
    pb_csv: Path | None = None
    pb_valid_only: bool = False
    # Mode B
    docking_directories: dict[str, Path] = field(default_factory=dict)
    receptors_dir: Path | None = None
    # crystal reference (benchmark only)
    benchmark_dir: Path | None = None
    # behaviour
    poses_per_combo: int = 3
    render_images: bool = True
    split_equibind: bool = True
    split_diffdock: bool = True
    use_dssp: bool = False
    num_workers: int | None = None
    overwrite: bool = False
    # best-EquiBind / best-DiffDock filter (map only the single best variant)
    best_equibind_only: bool = False
    best_diffdock_only: bool = False
    # pin DiffDock to a named optimizer variant ('raw' | 'smina' | 'gnina') instead
    # of oracle-ranking it; takes precedence over best_diffdock_only when set.
    diffdock_variant: str | None = None
    # pin EquiBind to a variant by SPEC tokens (e.g. 'gnina', 'unguided_gnina') instead
    # of oracle-ranking it; takes precedence over best_equibind_only when set.
    equibind_variant: str | None = None
    # pin AutoDock to one variant (e.g. 'gnina', 'autodock_gnina', 'raw'); drops every
    # other autodock* method, Vinardo family included. None → map them all.
    autodock_variant: str | None = None
    oracle_summary: Path | None = None
    # per-pose metrics table (posebusters) — source of gnina_affinity used to rank
    # EquiBind's gnina-optimised poses (the pb_csv carries no gnina column).
    per_pose_metrics: Path | None = None
    # restrict to an explicit '<PDBID>_<LIG>' id allow-list (one per line)
    ids_file: Path | None = None
    # profile EVERY pose AND the crystal reference of a complex against ONE receptor,
    # '<canonical_receptor_dir>/<complex_id>_protein.pdb', instead of the per-pose
    # receptor recorded in protein_file_used. See _canonical_receptor.
    canonical_receptor_dir: Path | None = None


def load_config(path: str | Path) -> PandaMapConfig:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    with open(path) as f:
        raw: dict[str, Any] = yaml.safe_load(f) or {}

    work_dir = Path(raw.get("work_dir", ".")).expanduser().resolve()

    def _resolve(p):
        p = Path(p).expanduser()
        return p if p.is_absolute() else (work_dir / p).resolve()

    pb_csv = _resolve(raw["pb_csv"]) if raw.get("pb_csv") else None
    docking = {k: _resolve(v) for k, v in (raw.get("docking_directories") or {}).items()}
    if not pb_csv and not docking:
        raise KeyError("Config must define either 'pb_csv' (Mode A) or 'docking_directories' (Mode B)")

    return PandaMapConfig(
        work_dir=work_dir,
        output_dir=_resolve(raw.get("output_dir", "pandamap_results")),
        pb_csv=pb_csv,
        pb_valid_only=bool(raw.get("pb_valid_only", False)),
        docking_directories=docking,
        receptors_dir=_resolve(raw["receptors_dir"]) if raw.get("receptors_dir") else None,
        benchmark_dir=_resolve(raw["benchmark_dir"]) if raw.get("benchmark_dir") else None,
        poses_per_combo=int(raw.get("poses_per_combo", 3)),
        render_images=bool(raw.get("render_images", True)),
        split_equibind=bool(raw.get("split_equibind", True)),
        split_diffdock=bool(raw.get("split_diffdock", True)),
        use_dssp=bool(raw.get("use_dssp", False)),
        num_workers=raw.get("num_workers"),
        overwrite=bool(raw.get("overwrite", False)),
        best_equibind_only=bool(raw.get("best_equibind_only", False)),
        best_diffdock_only=bool(raw.get("best_diffdock_only", False)),
        diffdock_variant=(raw.get("diffdock_variant") or None),
        equibind_variant=(raw.get("equibind_variant") or None),
        autodock_variant=(raw.get("autodock_variant") or None),
        oracle_summary=_resolve(raw["oracle_summary"]) if raw.get("oracle_summary") else None,
        per_pose_metrics=_resolve(raw["per_pose_metrics"]) if raw.get("per_pose_metrics") else None,
        ids_file=_resolve(raw["ids_file"]) if raw.get("ids_file") else None,
        canonical_receptor_dir=(_resolve(raw["canonical_receptor_dir"])
                                if raw.get("canonical_receptor_dir") else None),
    )


# ──────────────────────────────────────────────────────────────────────────
# EquiBind variant split + rank parsing
# (mirrors posebusters_pose_comparison._classify_equibind / parse_rank — kept
#  local to avoid importing that module, whose prolif/MDAnalysis import order is
#  segfault-sensitive in a forked worker.)
# ──────────────────────────────────────────────────────────────────────────

def _col(row, key):
    v = row.get(key) if hasattr(row, "get") else None
    if v is None:
        return None
    s = str(v).strip()
    return s if s and s.lower() not in ("nan", "none", "unknown", "<na>") else None


def classify_equibind(row) -> tuple[str, str | None, str | None]:
    """(pocket_source, refine, clamp) from CSV provenance cols or the filename."""
    name = Path(str(row.get("pose_name", ""))).name.lower()
    pocket = _col(row, "pocket_source")
    if pocket not in ("unguided", "fpocket", "p2rank"):
        pocket = ("unguided" if name.startswith("unguided")
                  else "p2rank" if name.startswith("p2rank")
                  else "fpocket" if name.startswith("fpocket") else "guided")
    refine = _col(row, "refine_variant")
    if refine not in ("smina", "raw", "gnina"):
        refine = ("smina" if "__refsmina" in name else
                  "gnina" if "__refgnina" in name else
                  "raw" if "__refraw" in name else None)
    clamp = _col(row, "clamp_variant")
    if clamp not in ("clampON", "clampOFF"):
        clamp = ("clampOFF" if "__clampoff" in name else "clampON" if "__clampon" in name else None)
    return pocket, refine, clamp


def equibind_label(row) -> str:
    pocket, refine, clamp = classify_equibind(row)
    lab = f"equibind_{pocket}"
    if refine:
        lab += f"_{refine}"
    if clamp:
        lab += f"_{clamp}"
    return lab


def classify_diffdock(row) -> str | None:
    """Optimizer backend for a DiffDock pose: 'smina' | 'gnina' | None (raw).

    Prefers the ``optimizer`` provenance column (written by run_posebusters from
    each pose's ``optimized_<tool>/`` subfolder); falls back to the path in
    ``pose_name``.
    """
    opt = _col(row, "optimizer")
    if opt in ("smina", "gnina"):
        return opt
    if opt == "original":
        return None
    name = str(row.get("pose_name", "")).lower()
    if "optimized_smina" in name:
        return "smina"
    if "optimized_gnina" in name:
        return "gnina"
    return None


def diffdock_label(row) -> str:
    """``diffdock`` (raw) | ``diffdock_smina`` | ``diffdock_gnina``."""
    o = classify_diffdock(row)
    return f"diffdock_{o}" if o else "diffdock"


# Suffixes that mark an AutoDock method key as a post-dock optimizer variant
# (matches posebusters_pose_comparison._AUTODOCK_OPTIMIZER_SUFFIXES).
_AUTODOCK_OPTIMIZER_SUFFIXES = ("_gnina_refinement", "_smina", "_gnina")


def classify_autodock(row) -> str | None:
    """Optimizer backend for an AutoDock pose: 'smina' | 'gnina' |
    'smina_refinement' | 'gnina_refinement' | None (raw Vina).

    Prefers the ``optimizer`` provenance column (written by run_posebusters from
    each pose's ``optimized_<tool>/`` subfolder); falls back to the path in
    ``pose_name``. Mirrors :func:`classify_diffdock`.

    Rescoring and refinement must stay distinct: rescoring only re-ranks the Vina
    geometry while refinement minimises it, so they are different poses carrying
    independent ``optimized_rank`` sequences. Both the column values and the
    ``optimized_<tool>_refinement/`` folder names are prefixed by the plain
    rescoring token, so the refinement variants are tested first — matching on the
    prefix instead would pool the two under one key and collide their ranks.
    """
    opt = _col(row, "optimizer")
    if opt in ("smina_refinement", "gnina_refinement", "smina", "gnina"):
        return opt
    if opt in ("original", "raw", "native", "none"):
        return None
    name = str(row.get("pose_name", "")).lower()
    for tok in ("optimized_smina_refinement", "optimized_gnina_refinement",
                "optimized_smina", "optimized_gnina"):
        if tok in name:
            return tok[len("optimized_"):]
    return None


def autodock_label(row) -> str:
    """Stable AutoDock method key, preserving the run tree and appending the
    optimizer: ``autodock`` | ``autodock_gnina`` | ``autodock_vinardo`` |
    ``autodock_mgltools_exh128`` | ``autodock_mgltools_exh128_gnina`` | …

    Raw and optimized poses MUST be separate methods — they are different
    geometries with different ranking axes, so pooling them under one 'autodock'
    key would let raw rank-1 and optimized rank-1 collide in select_top_n.

    The tree base is kept verbatim rather than folded to 'autodock'/'autodock_vinardo'.
    Folding was safe while the only AutoDock trees were the Meeko Vina and Vinardo
    pair, but the whole-protein campaign now also carries the ADFRsuite/MGLTools
    exhaustiveness ladder (``autodock_mgltools[_exh18|_exh64|_exh92|_exh128]``), which
    differs in ligand preparation and search effort. Folding those onto 'autodock'
    silently pools poses from up to six independent docking runs under one label and
    lets select_top_n mix them, and it makes every ladder arm unpinnable by
    ``autodock_variant``. Trees that were already 'autodock' / 'autodock_vinardo' keep
    exactly the keys they had, so published runs are unaffected.
    """
    base = str(row.get("docking_method") or row.get("method") or "autodock").strip().lower()
    for suf in _AUTODOCK_OPTIMIZER_SUFFIXES:
        if base.endswith(suf):
            base = base[: -len(suf)]
            break
    o = classify_autodock(row)
    return f"{base}_{o}" if o else base


def _positive_int(v) -> int | None:
    """Positive int from a CSV/dict value, else None (NaN / '' / <=0 → None)."""
    try:
        n = int(float(v))
    except (TypeError, ValueError):
        return None
    return n if n > 0 else None


def _pose_rank(row, method: str, pose_name: str) -> int:
    """Pose rank for select_top_n (1 = best).

    Optimized AutoDock (``autodock*_gnina`` / ``*_smina``) is ranked by
    ``optimized_rank`` — the gnina/smina re-ranking, which genuinely reorders the
    poses — NOT the original Vina rank encoded in the ``_rank<N>_`` filename token.
    Raw AutoDock uses its native ``autodock_rank`` column when present. Everything
    else (and any missing column) falls back to the filename via parse_rank.
    ``row`` may be a pandas Series (CSV mode) or a pose dict (dir mode).
    """
    if method.startswith("autodock"):
        col = "optimized_rank" if method.endswith(_AUTODOCK_OPTIMIZER_SUFFIXES) else "autodock_rank"
        r = _positive_int(row.get(col) if hasattr(row, "get") else None)
        if r is not None:
            return r
    return parse_rank(method, pose_name)


def parse_rank(method: str, pose_name: str) -> int:
    """Pose rank from the filename (1 = top). Unranked tools → 999."""
    name = Path(pose_name).name
    if method == "autodock":
        # Vina writes models in ascending affinity order (model N = native rank N).
        # The pose NAME uses 'poseN', the pose FILE 'modelN' — match either, so the
        # native rank is used instead of falling through to the 999 sentinel (which
        # made select_top_n keep the lexicographic, not rank-ordered, top poses).
        m = re.search(r"_?(?:model|pose)(\d+)", name)
        return int(m.group(1)) if m else 999
    if method == "diffdock":
        m = re.search(r"rank(\d+)", name)
        return int(m.group(1)) if m else 999
    m = re.search(r"(?:pose|rank|model)(\d+)", name)
    return int(m.group(1)) if m else 999


# ──────────────────────────────────────────────────────────────────────────
# Pose loading
# ──────────────────────────────────────────────────────────────────────────

_PROVENANCE = ("pocket_source", "clamp_variant", "refine_variant", "smina_affinity", "pocket_id", "optimizer", "pb_valid")


def _load_allowed_ids(path: Path) -> set[str]:
    """Load '<PDBID>_<LIG>' complex ids (one per line; '#' comments ignored)."""
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def _canonical_receptor(cfg: PandaMapConfig, protein: str) -> str | None:
    """Receptor every pose of complex *protein* should be profiled against, or None.

    A typed interaction fingerprint is keyed on (resname, resnum, chain), so poses are
    only comparable across methods — and against the crystal reference — when they are
    profiled against receptor files that share one residue-numbering scheme. They do
    not by default: each toolchain records its own prepared receptor in
    ``protein_file_used``, and some of those are renumbered copies. The MGLTools
    receptor used by the whole-protein AutoDock/Uni-Dock runs renumbers residues
    sequentially from 1 and the Uni-Dock2 receptor drops chain IDs entirely, so their
    fingerprint keys cannot match a crystal reference that carries author numbering,
    which silently collapses their measured native-interaction recovery towards zero.

    Pointing this at a directory of cleaned, author-numbered receptors
    (``<dir>/<complex_id>_protein.pdb``) puts every method and the crystal on one
    structure. Only correct when those files are coordinate-compatible with what each
    tool actually docked into — verify before enabling.
    """
    if not cfg.canonical_receptor_dir:
        return None
    p = cfg.canonical_receptor_dir / f"{protein}_protein.pdb"
    return str(p) if p.exists() else None


def load_poses_from_csv(cfg: PandaMapConfig) -> list[dict]:
    from run_posebusters import CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool  # noqa: E402
    df = pd.read_csv(cfg.pb_csv, low_memory=False)
    need = {"pose_file", "protein", "ligand", "docking_method"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"PB CSV missing columns: {missing}")
    df["docking_method"] = df["docking_method"].astype(str).str.lower()

    # Compute PoseBusters validity per pose ALWAYS, so it can be carried into the
    # fingerprint output as a 'pb_valid' flag (a pose passes iff every canonical PB
    # test is True). Only DROP invalid poses when pb_valid_only is on; otherwise every
    # pose is fingerprinted and merely tagged — giving full tool×complex coverage while
    # letting downstream still select valid-only. (Tools with catastrophic geometry,
    # e.g. EquiBind, otherwise vanish from complexes where all their top poses fail PB.)
    df = df.copy()
    checks = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
    if checks:
        sub = df[checks].copy()
        coerce_test_cols_to_bool(sub, checks)
        df["_pb_valid"] = sub.all(axis=1).to_numpy()
    else:
        df["_pb_valid"] = pd.NA
    if cfg.pb_valid_only and checks:
        df = df[df["_pb_valid"]].copy()

    poses = []
    for _, r in df.iterrows():
        method = r["docking_method"]
        if cfg.split_equibind and method.startswith("equibind"):
            method = equibind_label(r)
        elif cfg.split_diffdock and method.startswith("diffdock"):
            method = diffdock_label(r)
        elif method.startswith("autodock"):
            # Always separate raw Vina/Vinardo from their gnina/smina-optimized
            # variants — they are distinct geometries with distinct ranking axes.
            method = autodock_label(r)
        pose_file = str(r["pose_file"])
        pose_name = str(r.get("pose_name", Path(pose_file).stem))
        rec = {
            "method": method,
            "protein": str(r["protein"]),
            "ligand": str(r["ligand"]),
            "pose_file": pose_file,
            "pose_name": pose_name,
            "pose_rank": _pose_rank(r, method, pose_name),
            "protein_file": (_canonical_receptor(cfg, str(r["protein"]))
                             or (str(r["protein_file_used"])
                                 if _col(r, "protein_file_used") else None)),
            "pb_valid": (bool(r["_pb_valid"]) if pd.notna(r["_pb_valid"]) else None),
        }
        for c in _PROVENANCE:
            if c in df.columns:
                rec[c] = _col(r, c)
        poses.append(rec)
    return poses


def load_poses_from_dirs(cfg: PandaMapConfig) -> list[dict]:
    from run_posebusters import (  # noqa: E402
        ROW_COLLECTORS, POSE_EXPANDERS, discover_proteins, find_protein_file,
    )

    class _Ctx:  # POSE_EXPANDERS only touch .converted_dir / .work_dir / .ligand_template_dirs
        converted_dir = cfg.output_dir / "_converted_pdbqt"
        work_dir = cfg.work_dir
        ligand_template_dirs: list = []
    _Ctx.converted_dir.mkdir(parents=True, exist_ok=True)

    file_map, norm_map = ({}, {})
    if cfg.receptors_dir and cfg.receptors_dir.exists():
        file_map, norm_map = discover_proteins(cfg.receptors_dir)

    poses = []
    for method_key, d in cfg.docking_directories.items():
        collector = ROW_COLLECTORS.get(method_key)
        expander = POSE_EXPANDERS.get(method_key)
        if collector is None or expander is None or not Path(d).exists():
            print(f"  [skip] {method_key}: no collector/expander or dir missing ({d})")
            continue
        for row in collector(Path(d)):
            row["docking_tool"] = method_key
            for p in expander(row, _Ctx.converted_dir, _Ctx):
                method = p["method"]
                if cfg.split_equibind and method.startswith("equibind"):
                    method = equibind_label(p)
                elif cfg.split_diffdock and method.startswith("diffdock"):
                    method = diffdock_label(p)
                elif method.startswith("autodock"):
                    method = autodock_label(p)
                poses.append({
                    "method": method, "protein": p["protein"], "ligand": p["ligand"],
                    "pose_file": p["pose_file"], "pose_name": p["pose_name"],
                    "pose_rank": _pose_rank(p, method, p["pose_name"]),
                    "protein_file": find_protein_file(p["protein"], file_map, norm_map),
                    "pocket_source": p.get("pocket_source"),
                    "clamp_variant": p.get("clamp_variant"),
                    "refine_variant": p.get("refine_variant"),
                    "optimizer": p.get("optimizer"),
                })
    return poses


# Default location of the oracle summary written by posebusters_pose_comparison.py
# (the only place the per-variant success metrics exist). Used to rank EquiBind/DiffDock
# variants for --best-*-only (by PB-Valid AND RMSD ≤ 2 Å) when no explicit path is configured.
DEFAULT_ORACLE_SUMMARY = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/oracle_summary.csv")


def _resolve_variant_oracle(oracle_csv: Path | None) -> Path | None:
    """Prefer the ``*_all_variants.csv`` sibling of an oracle summary for per-variant ranking.

    posebusters_pose_comparison.py writes oracle_summary.csv with each tool COLLAPSED to a
    single row (``diffdock`` / ``equibind_unguided_gnina`` — its own best variant, relabelled
    to the bare tool name). PandaMap ranks per-variant labels (``diffdock_smina`` /
    ``diffdock_gnina`` / …); against the collapsed file only the bare ``diffdock`` row matches,
    so the other variants are dropped (``m in scores.index`` is False) and raw DiffDock "wins"
    as the sole candidate.

    The comparison also writes an ``oracle_summary_all_variants.csv`` sibling with one row per
    variant; when it exists we rank against that instead. Falls back to the given path (older
    runs / crystal-free sets lacking the sibling); a path already ending in ``_all_variants.csv``
    is returned unchanged.
    """
    if not oracle_csv:
        return oracle_csv
    p = Path(oracle_csv)
    if p.name.endswith("_all_variants.csv"):
        return p
    sibling = p.with_name(f"{p.stem}_all_variants{p.suffix}")
    return sibling if sibling.exists() else p


# Metric the best-variant filters rank on: the combined docking-success criterion
# PB-Valid AND RMSD ≤ 2 Å (a pose must be BOTH near-native AND physically valid), written
# per variant by posebusters_pose_comparison.py. Older summaries that predate the combined
# column fall back to the RMSD-only rate.
BEST_VARIANT_METRIC = "oracle_pb_valid_and_rmsd2_%"
BEST_VARIANT_METRIC_FALLBACK = "oracle_rmsd_le_2.0A_%"


def _variant_ranking_scores(oracle_csv: Path, tag: str) -> tuple[pd.Series | None, str | None]:
    """Per-variant ranking scores and the metric column used, read from *oracle_csv*.

    Ranks on ``oracle_pb_valid_and_rmsd2_%`` (PB-Valid AND RMSD ≤ 2 Å — the combined
    docking-success criterion), falling back to ``oracle_rmsd_le_2.0A_%`` when a summary
    predates the combined column. Returns (scores, metric) or (None, None) — with a ``[tag]``
    note — when the file carries neither column. The caller has already resolved *oracle_csv*
    to its ``*_all_variants.csv`` sibling (see _resolve_variant_oracle)."""
    osum = pd.read_csv(oracle_csv, index_col=0)
    col = (BEST_VARIANT_METRIC if BEST_VARIANT_METRIC in osum.columns
           else BEST_VARIANT_METRIC_FALLBACK if BEST_VARIANT_METRIC_FALLBACK in osum.columns
           else None)
    if col is None:
        print(f"  [{tag}] neither '{BEST_VARIANT_METRIC}' nor "
              f"'{BEST_VARIANT_METRIC_FALLBACK}' in {oracle_csv} — keeping all variants.")
        return None, None
    return pd.to_numeric(osum[col], errors="coerce"), col


def filter_best_equibind(poses: list[dict],
                         oracle_csv: Path | None) -> tuple[list[dict], str | None]:
    """Keep non-EquiBind poses + only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest PB-Valid AND RMSD ≤ 2 Å rate
    (``oracle_pb_valid_and_rmsd2_%``, falling back to ``oracle_rmsd_le_2.0A_%`` for older
    summaries — see _variant_ranking_scores) in *oracle_csv* — the oracle_summary.csv written
    by posebusters_pose_comparison.py (PandaMap carries no RMSD of its own). The chosen variant
    is selected among the EquiBind variants actually present in *poses*; AutoDock/DiffDock poses
    are always kept.

    Returns (kept_poses, best_variant_key). ``best`` is ``None`` and *poses* is
    returned unchanged when filtering can't be applied: no EquiBind variants
    present, missing/unreadable oracle summary, or no overlap between the present
    variants and the ranking metric.
    """
    eq_present = sorted({str(p["method"]) for p in poses
                         if str(p["method"]).startswith("equibind")})
    if not eq_present:
        return poses, None
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-equibind-only] oracle summary not found at {oracle_csv} — "
              "keeping all EquiBind variants.")
        return poses, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-equibind-only")
    if scores is None:
        return poses, None
    cand = {m: float(scores[m]) for m in eq_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-equibind-only] none of the present EquiBind variants have a "
              f"score in {oracle_csv} — keeping all EquiBind variants.")
        return poses, None

    best = max(cand, key=cand.get)
    kept = [p for p in poses
            if not str(p["method"]).startswith("equibind") or str(p["method"]) == best]
    return kept, best


def filter_equibind_variant(poses: list[dict],
                            spec: str) -> tuple[list[dict], str | None]:
    """Keep non-EquiBind poses + only the EquiBind variant(s) matching *spec*.

    Explicit, oracle-free counterpart to filter_best_equibind (mirrors
    posebusters_validity_report.select_equibind_variant). *spec* tokens (split on
    '_' or '/') are matched against the label's pocket/refine/clamp tokens, e.g.
    'gnina' or 'unguided_gnina'; a token subset constrains only the axes it names.
    AutoDock/DiffDock poses are always kept. Returns (kept_poses, kept_label) —
    kept_label is the retained variant when *spec* resolves to exactly one (for the
    report's 'EquiBind*' relabel), else ``None`` (poses unchanged when nothing
    matches, still filtered when several variants match).
    """
    tokens = [t for t in str(spec).strip().lower().replace("/", "_").split("_") if t]
    eq_present = sorted({str(p["method"]) for p in poses
                         if str(p["method"]).startswith("equibind")})
    if not tokens or not eq_present:
        return poses, None
    matched = [m for m in eq_present
               if all(t in set(m.lower().split("_")[1:]) for t in tokens)]
    if not matched:
        print(f"  [equibind-variant] no EquiBind variant matches '{spec}' "
              f"(present: {eq_present}) — keeping all EquiBind variants.")
        return poses, None
    kept = [p for p in poses
            if not str(p["method"]).startswith("equibind") or str(p["method"]) in matched]
    return kept, (matched[0] if len(matched) == 1 else None)


def filter_autodock_variant(poses: list[dict],
                            spec: str) -> tuple[list[dict], str | None]:
    """Keep non-AutoDock poses + only the single AutoDock variant named by *spec*.

    Explicit, oracle-free counterpart to filter_best_diffdock's ``pin`` mode. The
    AutoDock family spans two scoring functions (``autodock`` / ``autodock_vinardo``)
    crossed with the optimizer suffixes produced by :func:`autodock_label`, so a run
    that does not pin one maps every variant and multiplies the PandaMap cost.

    *spec* is matched against the full method key first (e.g. ``autodock_gnina``), then
    as an optimizer suffix on the Vina family (``gnina`` → ``autodock_gnina``), with
    ``raw``/``original`` resolving to the unoptimised ``autodock``. Every other
    ``autodock*`` method — including the whole Vinardo family — is dropped, so pinning
    is also how a run is restricted to one scoring function. DiffDock/EquiBind poses are
    always kept. Returns (kept_poses, kept_label); ``None`` when nothing matches.
    """
    spec = str(spec).strip().lower()
    ad_present = sorted({str(p["method"]) for p in poses
                         if str(p["method"]).startswith("autodock")})
    if not spec or not ad_present:
        return poses, None
    if spec in ad_present:
        target = spec
    elif spec in ("raw", "original") and "autodock" in ad_present:
        target = "autodock"
    elif f"autodock_{spec}" in ad_present:
        target = f"autodock_{spec}"
    else:
        print(f"  [autodock-variant] no AutoDock variant matches '{spec}' "
              f"(present: {ad_present}) — keeping all AutoDock variants.")
        return poses, None
    kept = [p for p in poses
            if not str(p["method"]).startswith("autodock") or str(p["method"]) == target]
    return kept, target


def filter_best_diffdock(poses: list[dict],
                         oracle_csv: Path | None,
                         pin: str | None = None) -> tuple[list[dict], str | None]:
    """Keep non-DiffDock poses + only one DiffDock optimizer variant.

    When *pin* is given ('raw' | 'smina' | 'gnina') the named variant is kept
    directly (``diffdock`` for 'raw', else ``diffdock_<pin>``), bypassing the
    oracle entirely — use this to force gnina regardless of the benchmark ranking.

    Otherwise "best" = the DiffDock variant (diffdock / diffdock_smina / diffdock_gnina)
    with the highest PB-Valid AND RMSD ≤ 2 Å rate (``oracle_pb_valid_and_rmsd2_%``, falling
    back to ``oracle_rmsd_le_2.0A_%`` for older summaries — see _variant_ranking_scores) in
    *oracle_csv*. Mirrors filter_best_equibind; AutoDock/EquiBind poses are always kept. Returns
    (kept_poses, kept_variant_key); ``None`` when no DiffDock variant is present, a pinned
    variant is absent, the summary is missing/unreadable, or no present variant has a score.
    """
    dd_present = sorted({str(p["method"]) for p in poses
                         if str(p["method"]).startswith("diffdock")})
    if not dd_present:
        return poses, None
    if pin:
        target = "diffdock" if pin == "raw" else f"diffdock_{pin}"
        if target not in dd_present:
            print(f"  [diffdock-variant={pin}] '{target}' not among present DiffDock "
                  f"variants {dd_present} — keeping all DiffDock variants.")
            return poses, None
        kept = [p for p in poses
                if not str(p["method"]).startswith("diffdock") or str(p["method"]) == target]
        return kept, target
    oracle_csv = _resolve_variant_oracle(oracle_csv)
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-diffdock-only] oracle summary not found at {oracle_csv} — "
              "keeping all DiffDock variants.")
        return poses, None

    scores, _ = _variant_ranking_scores(oracle_csv, "best-diffdock-only")
    if scores is None:
        return poses, None
    cand = {m: float(scores[m]) for m in dd_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-diffdock-only] none of the present DiffDock variants have a "
              f"score in {oracle_csv} — keeping all DiffDock variants.")
        return poses, None

    best = max(cand, key=cand.get)
    kept = [p for p in poses
            if not str(p["method"]).startswith("diffdock") or str(p["method"]) == best]
    return kept, best


def select_top_n(poses: list[dict], n: int) -> list[dict]:
    """Keep the top-n ranked poses per (method, protein, ligand)."""
    by_combo: dict[tuple, list[dict]] = {}
    for p in poses:
        by_combo.setdefault((p["method"], p["protein"], p["ligand"]), []).append(p)
    out = []
    for recs in by_combo.values():
        recs.sort(key=lambda r: (r["pose_rank"], r["pose_name"]))
        out.extend(recs[:n])
    return out


# Default source of gnina_affinity for EquiBind ranking (posebusters per-pose table;
# the pb_csv has no gnina column). Overridable via config 'per_pose_metrics'.
DEFAULT_PER_POSE_METRICS = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv")


def apply_gnina_ranks(poses: list[dict], per_pose_csv: Path | None,
                      work_dir: Path | None = None) -> list[dict]:
    """Rank EquiBind poses by gnina affinity (most-negative = rank 1), overriding pose_rank.

    EquiBind emits no native confidence score; its gnina-optimised poses carry the ranking
    signal only in their gnina affinity, which lives in the posebusters per_pose_metrics.csv
    (the pb_csv has no gnina column). This imports that affinity, keyed on
    (method, protein, ligand, pose_name), and assigns a per-complex 1..N rank so
    :func:`select_top_n` keeps the true top-N by gnina affinity rather than the arbitrary
    generation order. Non-EquiBind poses are untouched; an EquiBind complex whose poses have
    no gnina affinity (raw/smina variants) keeps its existing pose_rank. No-ops if the table
    is absent so the pipeline still runs without it.
    """
    eq = [p for p in poses if str(p["method"]).startswith("equibind")]
    if not eq:
        return poses
    p = Path(per_pose_csv) if per_pose_csv else DEFAULT_PER_POSE_METRICS
    if not p.is_absolute() and work_dir:
        p = work_dir / p
    if not p.exists():
        print(f"  [gnina-rank] per_pose_metrics not found at {p} — EquiBind keeps filename rank.")
        return poses
    ppm = pd.read_csv(p, low_memory=False)
    need = {"method", "protein", "ligand", "pose_name", "gnina_affinity"}
    if need - set(ppm.columns):
        print(f"  [gnina-rank] {p} missing {need - set(ppm.columns)} — EquiBind keeps filename rank.")
        return poses
    gaff = pd.to_numeric(ppm["gnina_affinity"], errors="coerce")
    lookup = {(str(m), str(pr), str(lg), str(nm)): a
              for m, pr, lg, nm, a in zip(ppm["method"], ppm["protein"],
                                          ppm["ligand"], ppm["pose_name"], gaff)}
    groups: dict[tuple, list[dict]] = {}
    for pose in eq:
        pose["_gaff"] = lookup.get((str(pose["method"]), str(pose["protein"]),
                                    str(pose["ligand"]), str(pose["pose_name"])))
        groups.setdefault((pose["method"], pose["protein"], pose["ligand"]), []).append(pose)
    n_ranked = 0
    for grp in groups.values():
        if not any(pd.notna(x["_gaff"]) for x in grp):
            continue                       # raw/smina variant — no gnina signal, leave as-is
        grp.sort(key=lambda x: (pd.isna(x["_gaff"]),
                                float(x["_gaff"]) if pd.notna(x["_gaff"]) else 0.0,
                                str(x["pose_name"])))
        for i, x in enumerate(grp, 1):
            x["pose_rank"] = i
            n_ranked += 1
    for pose in eq:
        pose.pop("_gaff", None)
    if n_ranked:
        print(f"  [gnina-rank] assigned gnina-affinity ranks to {n_ranked} EquiBind poses "
              f"(from {p.name}).")
    return poses


# ──────────────────────────────────────────────────────────────────────────
# Complex assembly + PandaMap worker
# ──────────────────────────────────────────────────────────────────────────

def combine_protein_ligand_to_pdb(protein_pdb, ligand_file, out_pdb, resname="LIG") -> bool:
    """Write protein + docked ligand into one PDB for PandaMap.

    The ligand is forced to ``HETATM`` records with residue name *resname* and
    chain ``X`` (PandaMap identifies the ligand by HETATM + resname). Crystal
    waters are dropped; other receptor HETATMs (kept cofactors/metals) pass
    through. Ported and hardened from panda_maps.ipynb cell 6.
    """
    from rdkit import Chem
    try:
        prot_lines = []
        for ln in open(protein_pdb):
            if ln.startswith(("ATOM", "TER")):
                prot_lines.append(ln)
            elif ln.startswith("HETATM") and ln[17:20].strip().upper() not in ("HOH", "WAT", "DOD"):
                prot_lines.append(ln)

        ext = Path(ligand_file).suffix.lower()
        if ext == ".sdf":
            mol = Chem.MolFromMolFile(str(ligand_file), removeHs=False, sanitize=False)
        elif ext in (".pdb", ".pdbqt"):
            mol = Chem.MolFromPDBFile(str(ligand_file), removeHs=False, sanitize=False)
        else:
            return False
        if mol is None:
            return False

        lig_lines = []
        for ln in Chem.MolToPDBBlock(mol).split("\n"):
            if ln.startswith(("HETATM", "ATOM")):
                ln = ln.ljust(80)
                lig_lines.append(("HETATM" + ln[6:17] + f"{resname:>3}" + " X" + ln[22:]).rstrip())
        if not lig_lines:
            return False

        with open(out_pdb, "w") as f:
            f.writelines(prot_lines)
            f.write("\n" + "\n".join(lig_lines) + "\nEND\n")
        return True
    except Exception:
        return False


def _residue_rows(interactions: dict) -> list[dict]:
    """Flatten mapper.interactions into residue-level records."""
    rows = []
    for itype, recs in interactions.items():
        if not isinstance(recs, list):
            continue
        for rec in recs:
            pr = rec.get("protein_residue") if isinstance(rec, dict) else None
            la = rec.get("ligand_atom") if isinstance(rec, dict) else None
            resname = getattr(pr, "resname", None)
            rid = getattr(pr, "id", (None, None, None))
            resnum = rid[1] if (pr is not None and len(rid) > 1) else None
            chain = getattr(getattr(pr, "parent", None), "id", None)
            rows.append({
                "interaction_type": itype,
                "resname": resname,
                "resnum": resnum,
                "chain": chain,
                "lig_atom_element": getattr(la, "element", None),
                "distance": round(float(rec["distance"]), 3) if isinstance(rec, dict) and rec.get("distance") is not None else None,
            })
    return rows


# worker globals
_W_RENDER = False
_W_USE_DSSP = False
_W_OUTDIR: Path | None = None


def _init_worker(render: bool, use_dssp: bool, out_dir: str):
    global _W_RENDER, _W_USE_DSSP, _W_OUTDIR
    import matplotlib
    matplotlib.use("Agg")
    _W_RENDER = render
    _W_USE_DSSP = use_dssp
    _W_OUTDIR = Path(out_dir)
    warnings.filterwarnings("ignore")


def _process_pose(pose: dict):
    """Map one pose. Returns (summary_row, interaction_rows, error_row|None)."""
    from pandamap import HybridProtLigMapper
    import matplotlib.pyplot as plt

    method, protein, ligand = pose["method"], pose["protein"], pose["ligand"]
    pose_file = pose["pose_file"]
    pose_name = pose["pose_name"]
    protein_file = pose.get("protein_file")

    base = {"method": method, "protein": protein, "ligand": ligand,
            "pose_name": pose_name, "pose_file": pose_file,
            "pose_rank": pose.get("pose_rank", 999)}
    for c in _PROVENANCE:
        if c in pose:
            base[c] = pose[c]

    if not protein_file or not Path(protein_file).exists():
        return None, [], {**base, "error": f"protein file not found: {protein_file}"}
    if not Path(pose_file).exists():
        return None, [], {**base, "error": f"pose file not found: {pose_file}"}

    safe = re.sub(r"[^A-Za-z0-9._-]", "_", pose_name)
    tmp_pdb = Path(tempfile.gettempdir()) / f"pm_{os.getpid()}_{safe[:60]}.pdb"
    try:
        if not combine_protein_ligand_to_pdb(protein_file, pose_file, str(tmp_pdb)):
            return None, [], {**base, "error": "combine_protein_ligand_to_pdb failed"}

        png = None
        if _W_RENDER and _W_OUTDIR is not None:
            img_dir = _W_OUTDIR / "maps" / method / f"{protein}__{ligand}"
            img_dir.mkdir(parents=True, exist_ok=True)
            png = img_dir / f"{safe[:80]}.png"

        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            mapper = HybridProtLigMapper(str(tmp_pdb), ligand_resname="LIG")
            if png is not None:
                # run_analysis detects + renders; same fingerprint as detect_interactions.
                mapper.run_analysis(output_file=str(png), use_dssp=_W_USE_DSSP)
            else:
                mapper.detect_interactions()
            plt.close("all")

        interactions = getattr(mapper, "interactions", {}) or {}
        rows = _residue_rows(interactions)
        for r in rows:
            r.update(base)
        counts = {t: 0 for t in INTERACTION_TYPES}
        for t, lst in interactions.items():
            if isinstance(lst, list):
                counts[t] = len(lst)
        summary = {**base, **counts,
                   "total_interactions": sum(counts.values()),
                   "n_residues": len({(r["resname"], r["resnum"], r["chain"]) for r in rows}),
                   "image": str(png) if png else ""}
        return summary, rows, None
    except Exception as e:
        return None, [], {**base, "error": f"{type(e).__name__}: {e}"}
    finally:
        try:
            tmp_pdb.unlink(missing_ok=True)
        except OSError:
            pass


# ──────────────────────────────────────────────────────────────────────────
# Crystal reference tasks
# ──────────────────────────────────────────────────────────────────────────

def build_crystal_tasks(selected: list[dict], cfg: PandaMapConfig) -> list[dict]:
    """One crystal-ligand task per (protein, ligand) pair (benchmark only)."""
    if not cfg.benchmark_dir:
        return []
    # protein_file per pair: borrow the receptor used by any docked pose of the pair.
    prot_for_pair, seen = {}, set()
    for p in selected:
        key = (p["protein"], p["ligand"])
        if key not in prot_for_pair and p.get("protein_file"):
            prot_for_pair[key] = p["protein_file"]
    tasks = []
    for (protein, ligand) in {(p["protein"], p["ligand"]) for p in selected}:
        if (protein, ligand) in seen:
            continue
        seen.add((protein, ligand))
        crystal = cfg.benchmark_dir / protein / f"{protein}_ligand.sdf"
        pfile = _canonical_receptor(cfg, protein) or prot_for_pair.get((protein, ligand))
        if not pfile:
            pfile = str(cfg.benchmark_dir / protein / f"{protein}_protein.pdb")
        if crystal.exists() and Path(pfile).exists():
            tasks.append({
                "method": "crystal", "protein": protein, "ligand": ligand,
                "pose_file": str(crystal), "pose_name": f"{protein}__{ligand}/crystal",
                "pose_rank": 0, "protein_file": pfile,
            })
    return tasks


# ──────────────────────────────────────────────────────────────────────────
# Driver
# ──────────────────────────────────────────────────────────────────────────

def _run_pool(tasks, cfg, label):
    n_workers = cfg.num_workers or max(1, (cpu_count() or 4) - 1)
    n_workers = max(1, min(n_workers, len(tasks)))
    print(f"\n{label}: {len(tasks)} poses on {n_workers} workers "
          f"(render={cfg.render_images})")
    summaries, inter_rows, errors = [], [], []
    with Pool(n_workers, initializer=_init_worker,
              initargs=(cfg.render_images, cfg.use_dssp, str(cfg.output_dir))) as pool:
        for i, (summary, rows, err) in enumerate(pool.imap_unordered(_process_pose, tasks, chunksize=4), 1):
            if err:
                errors.append(err)
            else:
                summaries.append(summary)
                inter_rows.extend(rows)
            if i % 50 == 0 or i == len(tasks):
                print(f"  [{i}/{len(tasks)}] done  (ok={len(summaries)} err={len(errors)})")
    return summaries, inter_rows, errors


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", "-c", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, default=None)
    ap.add_argument("--limit-pairs", type=int, default=0, help="process only first N pairs (debug)")
    ap.add_argument("--poses-per-combo", type=int, default=None)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--render", action=argparse.BooleanOptionalAction, default=None,
                    help="render 2D PNG maps (slow). Default from config.")
    ap.add_argument("--overwrite", action="store_true", default=None)
    ap.add_argument("--pb-valid-only", action=argparse.BooleanOptionalAction, default=None,
                    help="Fingerprint only PoseBusters-valid poses (config 'pb_valid_only'). "
                         "Use --no-pb-valid-only to fingerprint ALL top-N poses regardless of "
                         "validity (every pose is still tagged pb_valid) — gives full "
                         "tool×complex coverage when a tool's top poses all fail PB "
                         "(e.g. EquiBind's distorted geometry).")
    ap.add_argument("--best-equibind-only", action="store_true", default=None,
                    help="Map only the single best EquiBind variant (highest PB-Valid AND "
                         "RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%, from --oracle-summary) "
                         "alongside AutoDock/DiffDock. Overrides config 'best_equibind_only'.")
    ap.add_argument("--equibind-variant", default=None, metavar="SPEC",
                    help="Pin EquiBind to the variant matching SPEC (tokens split on '_' or "
                         "'/', e.g. 'gnina' or 'unguided_gnina') instead of oracle-ranking it. "
                         "Takes precedence over --best-equibind-only. Overrides config "
                         "'equibind_variant'.")
    ap.add_argument("--best-diffdock-only", action="store_true", default=None,
                    help="Map only the single best DiffDock optimizer variant "
                         "(highest PB-Valid AND RMSD ≤ 2 Å = oracle_pb_valid_and_rmsd2_%%, "
                         "from --oracle-summary) alongside AutoDock/EquiBind. Overrides config "
                         "'best_diffdock_only'.")
    ap.add_argument("--diffdock-variant", choices=("raw", "smina", "gnina"), default=None,
                    help="Pin DiffDock to this optimizer variant (e.g. 'gnina') instead of "
                         "oracle-ranking it. Takes precedence over --best-diffdock-only. "
                         "Overrides config 'diffdock_variant'.")
    ap.add_argument("--autodock-variant", default=None, metavar="SPEC",
                    help="Pin AutoDock to one variant and drop every other autodock* method, "
                         "Vinardo included (e.g. 'gnina' → autodock_gnina, 'raw' → autodock, "
                         "or a full key such as 'autodock_vinardo_gnina'). Overrides config "
                         "'autodock_variant'.")
    ap.add_argument("--oracle-summary", type=Path, default=None,
                    help="oracle_summary.csv (from posebusters_pose_comparison.py) "
                         "used to rank EquiBind variants for --best-equibind-only. "
                         "Overrides config 'oracle_summary'.")
    ap.add_argument("--ids-file", type=Path, default=None,
                    help="Restrict analysis to the '<PDBID>_<LIG>' complex ids "
                         "listed in this file (one per line; '#' comments ok). "
                         "Overrides config 'ids_file'.")
    args = ap.parse_args()

    try:
        import pandamap  # noqa: F401
    except Exception:
        sys.exit("ERROR: 'pandamap' not importable. Run with the vina env python:\n"
                 "  /home/manndo/anaconda3/envs/vina/bin/python Scripts/Analysis/run_pandamap.py ...")

    cfg = load_config(args.config)
    if args.out_dir:
        cfg.output_dir = args.out_dir.resolve()
    if args.poses_per_combo is not None:
        cfg.poses_per_combo = args.poses_per_combo
    if args.workers is not None:
        cfg.num_workers = args.workers
    if args.render is not None:
        cfg.render_images = args.render
    if args.overwrite is not None:
        cfg.overwrite = args.overwrite
    if args.pb_valid_only is not None:
        cfg.pb_valid_only = args.pb_valid_only
    if args.best_equibind_only is not None:
        cfg.best_equibind_only = args.best_equibind_only
    if args.equibind_variant is not None:
        cfg.equibind_variant = args.equibind_variant
    if args.best_diffdock_only is not None:
        cfg.best_diffdock_only = args.best_diffdock_only
    if args.diffdock_variant is not None:
        cfg.diffdock_variant = args.diffdock_variant
    if args.autodock_variant is not None:
        cfg.autodock_variant = args.autodock_variant
    if args.oracle_summary is not None:
        cfg.oracle_summary = args.oracle_summary.resolve()
    if args.ids_file is not None:
        cfg.ids_file = args.ids_file.resolve()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("PandaMap interaction generator")
    print(f"  input   : {'PB CSV ' + str(cfg.pb_csv) if cfg.pb_csv else 'dir scan'}")
    print(f"  output  : {cfg.output_dir}")
    print(f"  top-N   : {cfg.poses_per_combo} per method×pair   pb_valid_only={cfg.pb_valid_only}")
    print(f"  crystal : {cfg.benchmark_dir or '(none)'}")
    print("=" * 70)

    poses = load_poses_from_csv(cfg) if cfg.pb_csv else load_poses_from_dirs(cfg)
    # Restrict to the official benchmark-set ids (complex id == 'protein' field,
    # which equals '<PDBID>_<LIG>' for the benchmark staging).
    if cfg.ids_file:
        allowed = _load_allowed_ids(cfg.ids_file)
        p0 = len({(p["protein"], p["ligand"]) for p in poses})
        poses = [p for p in poses if p["protein"] in allowed]
        p1 = len({(p["protein"], p["ligand"]) for p in poses})
        print(f"Restricted to {len(allowed)} ids from {cfg.ids_file.name}: "
              f"{p0} → {p1} pairs ({len(poses)} poses)")
    print(f"Collected {len(poses)} poses; methods: {sorted({p['method'] for p in poses})}")

    if cfg.equibind_variant:
        before = len(poses)
        poses, best_eq = filter_equibind_variant(poses, cfg.equibind_variant)
        if best_eq:
            print(f"equibind-variant={cfg.equibind_variant}: keeping '{best_eq}'; "
                  f"{before} → {len(poses)} poses. "
                  "Run the report with the same --equibind-variant to label it 'EquiBind*'.")
    elif cfg.best_equibind_only:
        oracle_csv = cfg.oracle_summary or (cfg.work_dir / DEFAULT_ORACLE_SUMMARY)
        before = len(poses)
        poses, best_eq = filter_best_equibind(poses, oracle_csv)
        if best_eq:
            print(f"best-equibind-only: keeping '{best_eq}' (top EquiBind by "
                  f"PB-Valid AND RMSD ≤ 2 Å); {before} → {len(poses)} poses. "
                  "Run the report with --best-equibind-only to label it 'EquiBind*'.")

    if cfg.autodock_variant:
        before = len(poses)
        poses, best_ad = filter_autodock_variant(poses, cfg.autodock_variant)
        if best_ad:
            print(f"autodock-variant={cfg.autodock_variant}: keeping '{best_ad}'; "
                  f"{before} → {len(poses)} poses.")

    if cfg.diffdock_variant:
        before = len(poses)
        poses, best_dd = filter_best_diffdock(poses, None, pin=cfg.diffdock_variant)
        if best_dd:
            print(f"diffdock-variant={cfg.diffdock_variant}: keeping '{best_dd}'; "
                  f"{before} → {len(poses)} poses. "
                  "Run the report with the same --diffdock-variant to label it 'DiffDock*'.")
    elif cfg.best_diffdock_only:
        oracle_csv = cfg.oracle_summary or (cfg.work_dir / DEFAULT_ORACLE_SUMMARY)
        before = len(poses)
        poses, best_dd = filter_best_diffdock(poses, oracle_csv)
        if best_dd:
            print(f"best-diffdock-only: keeping '{best_dd}' (top DiffDock by "
                  f"PB-Valid AND RMSD ≤ 2 Å); {before} → {len(poses)} poses. "
                  "Run the report with --best-diffdock-only to label it 'DiffDock*'.")

    # EquiBind has no native rank; rank its gnina-optimised poses by gnina affinity
    # (imported from posebusters per_pose_metrics.csv) so select_top_n keeps the true
    # top-N rather than the arbitrary generation order.
    poses = apply_gnina_ranks(poses, cfg.per_pose_metrics, cfg.work_dir)

    selected = select_top_n(poses, cfg.poses_per_combo)
    if args.limit_pairs:
        keep_pairs = sorted({(p["protein"], p["ligand"]) for p in selected})[:args.limit_pairs]
        keep = set(keep_pairs)
        selected = [p for p in selected if (p["protein"], p["ligand"]) in keep]
    n_pairs = len({(p["protein"], p["ligand"]) for p in selected})
    print(f"Selected {len(selected)} poses across {n_pairs} pairs "
          f"({len({p['method'] for p in selected})} methods)")

    # resume: skip poses already in an existing summary unless --overwrite
    summary_csv = cfg.output_dir / "pandamap_pose_summary.csv"
    inter_csv = cfg.output_dir / "pandamap_interactions.csv"
    prev_summary = prev_inter = None
    if summary_csv.exists() and not cfg.overwrite:
        prev_summary = read_pandamap_csv(summary_csv)
        done = set(prev_summary["pose_file"].astype(str)) if "pose_file" in prev_summary else set()
        before = len(selected)
        selected = [p for p in selected if p["pose_file"] not in done]
        if inter_csv.exists():
            prev_inter = read_pandamap_csv(inter_csv)
        print(f"Resume: {before - len(selected)} already done, {len(selected)} to do")

    summaries, inter_rows, errors = ([], [], [])
    if selected:
        summaries, inter_rows, errors = _run_pool(selected, cfg, "Docked poses")

    # crystal references (always recomputed fresh unless present)
    crystal_csv = cfg.output_dir / "crystal_interactions.csv"
    if cfg.benchmark_dir and (cfg.overwrite or not crystal_csv.exists()):
        ctasks = build_crystal_tasks(
            poses if not args.limit_pairs else selected + [p for p in poses
                if (p["protein"], p["ligand"]) in {(s["protein"], s["ligand"]) for s in summaries}],
            cfg)
        if args.limit_pairs:
            keep = {(p["protein"], p["ligand"]) for p in (summaries or selected)}
            ctasks = [t for t in ctasks if (t["protein"], t["ligand"]) in keep]
        if ctasks:
            csum, crows, cerr = _run_pool(ctasks, cfg, "Crystal references")
            pd.DataFrame(crows).to_csv(crystal_csv, index=False)
            pd.DataFrame(csum).to_csv(cfg.output_dir / "crystal_pose_summary.csv", index=False)
            print(f"  crystal: {len(csum)} mapped, {len(cerr)} failed → {crystal_csv.name}")

    # write docked outputs (merge with resume data)
    sum_df = pd.DataFrame(summaries)
    if prev_summary is not None:
        sum_df = pd.concat([prev_summary, sum_df], ignore_index=True)
    sum_df.to_csv(summary_csv, index=False)

    inter_df = pd.DataFrame(inter_rows)
    if prev_inter is not None:
        inter_df = pd.concat([prev_inter, inter_df], ignore_index=True)
    inter_df.to_csv(inter_csv, index=False)

    if errors:
        pd.DataFrame(errors).to_csv(cfg.output_dir / "pandamap_errors.csv", index=False)

    print("\n" + "=" * 70)
    print(f"DONE  ok={len(summaries)}  errors={len(errors)}")
    print(f"  {summary_csv}  ({len(sum_df)} poses)")
    print(f"  {inter_csv}  ({len(inter_df)} residue-level interaction rows)")
    if not sum_df.empty:
        print("\nMean total interactions per pose by method:")
        print(sum_df.groupby("method")["total_interactions"].mean().round(2).sort_values(ascending=False).to_string())
    print("=" * 70)


if __name__ == "__main__":
    main()
