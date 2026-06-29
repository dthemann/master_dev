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
    use_dssp: bool = False
    num_workers: int | None = None
    overwrite: bool = False
    # best-EquiBind filter (map only the single best EquiBind variant)
    best_equibind_only: bool = False
    oracle_summary: Path | None = None
    # restrict to an explicit '<PDBID>_<LIG>' id allow-list (one per line)
    ids_file: Path | None = None


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
        use_dssp=bool(raw.get("use_dssp", False)),
        num_workers=raw.get("num_workers"),
        overwrite=bool(raw.get("overwrite", False)),
        best_equibind_only=bool(raw.get("best_equibind_only", False)),
        oracle_summary=_resolve(raw["oracle_summary"]) if raw.get("oracle_summary") else None,
        ids_file=_resolve(raw["ids_file"]) if raw.get("ids_file") else None,
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
    if refine not in ("smina", "raw"):
        refine = ("smina" if "__refsmina" in name else "raw" if "__refraw" in name else None)
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


def parse_rank(method: str, pose_name: str) -> int:
    """Pose rank from the filename (1 = top). Unranked tools → 999."""
    name = Path(pose_name).name
    if method == "autodock":
        m = re.search(r"_?model(\d+)", name)
        return int(m.group(1)) if m else 999
    if method == "diffdock":
        m = re.search(r"rank(\d+)", name)
        return int(m.group(1)) if m else 999
    m = re.search(r"(?:pose|rank|model)(\d+)", name)
    return int(m.group(1)) if m else 999


# ──────────────────────────────────────────────────────────────────────────
# Pose loading
# ──────────────────────────────────────────────────────────────────────────

_PROVENANCE = ("pocket_source", "clamp_variant", "refine_variant", "smina_affinity", "pocket_id")


def _load_allowed_ids(path: Path) -> set[str]:
    """Load '<PDBID>_<LIG>' complex ids (one per line; '#' comments ignored)."""
    return {ln.strip() for ln in Path(path).read_text().splitlines()
            if ln.strip() and not ln.lstrip().startswith("#")}


def load_poses_from_csv(cfg: PandaMapConfig) -> list[dict]:
    from run_posebusters import CANONICAL_TEST_COLUMNS, coerce_test_cols_to_bool  # noqa: E402
    df = pd.read_csv(cfg.pb_csv, low_memory=False)
    need = {"pose_file", "protein", "ligand", "docking_method"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"PB CSV missing columns: {missing}")
    df["docking_method"] = df["docking_method"].astype(str).str.lower()

    if cfg.pb_valid_only:
        checks = [c for c in CANONICAL_TEST_COLUMNS if c in df.columns]
        if checks:
            sub = df[checks].copy()
            coerce_test_cols_to_bool(sub, checks)
            df = df[sub.all(axis=1)].copy()

    poses = []
    for _, r in df.iterrows():
        method = r["docking_method"]
        if cfg.split_equibind and method.startswith("equibind"):
            method = equibind_label(r)
        pose_file = str(r["pose_file"])
        pose_name = str(r.get("pose_name", Path(pose_file).stem))
        rec = {
            "method": method,
            "protein": str(r["protein"]),
            "ligand": str(r["ligand"]),
            "pose_file": pose_file,
            "pose_name": pose_name,
            "pose_rank": parse_rank(r["docking_method"], pose_name),
            "protein_file": str(r["protein_file_used"]) if _col(r, "protein_file_used") else None,
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
                poses.append({
                    "method": method, "protein": p["protein"], "ligand": p["ligand"],
                    "pose_file": p["pose_file"], "pose_name": p["pose_name"],
                    "pose_rank": parse_rank(p["method"], p["pose_name"]),
                    "protein_file": find_protein_file(p["protein"], file_map, norm_map),
                    "pocket_source": p.get("pocket_source"),
                    "clamp_variant": p.get("clamp_variant"),
                    "refine_variant": p.get("refine_variant"),
                })
    return poses


# Default location of the oracle summary written by posebusters_pose_comparison.py
# (the only place the oracle_rmsd_le_2.0A_% metric exists). Used to rank EquiBind
# variants for --best-equibind-only when no explicit path is configured.
DEFAULT_ORACLE_SUMMARY = Path(
    "posebusters_results/benchmark/dock/pose_comparison_report/oracle_summary.csv")


def filter_best_equibind(poses: list[dict],
                         oracle_csv: Path | None) -> tuple[list[dict], str | None]:
    """Keep non-EquiBind poses + only the single best EquiBind variant.

    "Best" = the EquiBind variant with the highest oracle RMSD ≤ 2 Å success rate
    (``oracle_rmsd_le_2.0A_%``) in *oracle_csv* — the oracle_summary.csv written by
    posebusters_pose_comparison.py (PandaMap carries no RMSD of its own). The
    chosen variant is selected among the EquiBind variants actually present in
    *poses*; AutoDock/DiffDock poses are always kept.

    Returns (kept_poses, best_variant_key). ``best`` is ``None`` and *poses* is
    returned unchanged when filtering can't be applied: no EquiBind variants
    present, missing/unreadable oracle summary, or no overlap between the present
    variants and the oracle metric.
    """
    eq_present = sorted({str(p["method"]) for p in poses
                         if str(p["method"]).startswith("equibind")})
    if not eq_present:
        return poses, None
    if not oracle_csv or not Path(oracle_csv).exists():
        print(f"  [best-equibind-only] oracle summary not found at {oracle_csv} — "
              "keeping all EquiBind variants.")
        return poses, None

    col = "oracle_rmsd_le_2.0A_%"
    osum = pd.read_csv(oracle_csv, index_col=0)
    if col not in osum.columns:
        print(f"  [best-equibind-only] '{col}' missing from {oracle_csv} — "
              "keeping all EquiBind variants.")
        return poses, None

    scores = pd.to_numeric(osum[col], errors="coerce")
    cand = {m: float(scores[m]) for m in eq_present
            if m in scores.index and pd.notna(scores[m])}
    if not cand:
        print("  [best-equibind-only] none of the present EquiBind variants have an "
              f"oracle score in {oracle_csv} — keeping all EquiBind variants.")
        return poses, None

    best = max(cand, key=cand.get)
    kept = [p for p in poses
            if not str(p["method"]).startswith("equibind") or str(p["method"]) == best]
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
        pfile = prot_for_pair.get((protein, ligand))
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
    ap.add_argument("--best-equibind-only", action="store_true", default=None,
                    help="Map only the single best EquiBind variant (highest "
                         "oracle_rmsd_le_2.0A_%%, from --oracle-summary) alongside "
                         "AutoDock/DiffDock. Overrides config 'best_equibind_only'.")
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
    if args.best_equibind_only is not None:
        cfg.best_equibind_only = args.best_equibind_only
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

    if cfg.best_equibind_only:
        oracle_csv = cfg.oracle_summary or (cfg.work_dir / DEFAULT_ORACLE_SUMMARY)
        before = len(poses)
        poses, best_eq = filter_best_equibind(poses, oracle_csv)
        if best_eq:
            print(f"best-equibind-only: keeping '{best_eq}' (top EquiBind by "
                  f"oracle_rmsd_le_2.0A_%); {before} → {len(poses)} poses. "
                  "Run the report with --best-equibind-only to label it 'EquiBind*'.")

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
        prev_summary = pd.read_csv(summary_csv)
        done = set(prev_summary["pose_file"].astype(str)) if "pose_file" in prev_summary else set()
        before = len(selected)
        selected = [p for p in selected if p["pose_file"] not in done]
        if inter_csv.exists():
            prev_inter = pd.read_csv(inter_csv)
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
