#!/usr/bin/env python3
"""Re-run the EquiBind post-docking refinement under a different local-optimisation mode.

WHY THIS EXISTS
The EquiBind arms were refined with gnina/smina in local-search mode
(``--local_only``). The AutoDock Vina and DiffDock arms were refined in
energy-minimisation mode (``--minimize``), which gnina's own help text
recommends ("you probably want to use --minimize"). RQ2 compares the three
arms' response to post-docking optimisation, so the arms should receive the
same operation. This script produces the matched ``--minimize`` EquiBind arm
in a separate output tree, leaving the published ``--local_only`` tree
untouched.

HOW IT STAYS LIKE-FOR-LIKE
``pipeline._emit_variants`` writes the ``__refRAW`` variant with
``shutil.copy2(geom_sdf, target)`` and then appends provenance tags at the text
level, leaving the molecule block byte-for-byte intact. The committed
``*__refRAW.sdf`` is therefore exactly the input that ``refine_pose`` received
when it produced the sibling ``*__refGNINA.sdf``. Re-refining from
``*__refRAW.sdf`` through the pipeline's own ``refine_pose`` changes one thing
only: the search flag.

The receptor, autobox padding, seed, cpu count, timeout and device are all read
from the same config surface the original run used, and the failure fallback
(keep the pre-refine pose under the refined label) is reproduced so pose counts
match the source tree exactly.

USAGE
    python3 Scripts/Docking/rerun_equibind_refine.py \
        --source-root Dockings/Benchmark_Equibind \
        --dest-root   Dockings/Benchmark_Equibind_minimize \
        --search minimize --tools gnina,smina --workers 30
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Executable + device pinning ──────────────────────────────────────────────
# Set before importing equibind_pipeline: Config resolves these at class-body
# evaluation time. gnina is not on PATH in any env, so the original runs used
# the docking_tools wrapper; smina came from the equibind env's bin/.
os.environ.setdefault("EQ_GNINA", "/home/manndo/docking_tools/gnina")
os.environ.setdefault("EQ_SMINA", "/home/manndo/anaconda3/envs/equibind/bin/smina")
os.environ.setdefault("EQ_GNINA_USE_GPU", "false")

# Pin math-library threads to 1, as the notebook launcher cell does. gnina's CNN
# rescore grabs ~all cores per process via OpenMP and IGNORES its --cpu 1 flag
# (measured: 33 threads per process, 2 once pinned). Set before the
# equibind_pipeline import so the pool workers and refine_pose's subprocesses
# inherit it. Override externally to experiment; these are setdefault.
#
# THE PIN ONLY PAYS OFF WITH ENOUGH WORKERS. Measured on this box (32 cores,
# gnina on GPU, gnina+smina per pose):
#     unpinned, 14 workers  1.24 poses/s
#     pinned,   14 workers  0.59 poses/s   <- pin alone is a REGRESSION
#     pinned,   28 workers  2.14 poses/s   <- use this
# One thread per worker across only 14 workers leaves half the machine idle. The
# notebook's "~8x faster" note describes its own CPU-mode arm with ~31 workers,
# where 31 x 33 threads genuinely thrashed; it does not transfer to a GPU run at
# low worker counts. Keep --workers near the core count when pinning.
for _var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from equibind_pipeline.config import CFG          # noqa: E402
from equibind_pipeline.refine import refine_pose  # noqa: E402

RAW_SUFFIX = "__refRAW.sdf"
# Provenance tags carried over from the raw pose onto the refined pose. The
# pipeline re-derives these from the in-memory job; here they are copied from
# the raw file, which already carries them.
_CARRY_TAGS = ("pocket_source", "pocket_id", "clamp_variant")


# ── SDF tag helpers (text level, molecule block untouched) ───────────────────
def _read_tags(sdf: Path) -> Dict[str, str]:
    """Parse ``>  <key>\\nvalue`` blocks out of an SDF."""
    tags: Dict[str, str] = {}
    try:
        lines = sdf.read_text().splitlines()
    except Exception:
        return tags
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("> ") and "<" in s and ">" in s[2:]:
            key = s[s.index("<") + 1:s.rindex(">")]
            if i + 1 < len(lines) and key not in tags:
                tags[key] = lines[i + 1].strip()
    return tags


def _tag_pose(sdf_path: Path, carried: Dict[str, str], refine_variant: str,
              affinity: Optional[float], cnn_score: Optional[float],
              cnn_affinity: Optional[float]) -> None:
    """Mirror pipeline._tag_pose_source onto a refined pose."""
    try:
        text = sdf_path.read_text()
        block = f">  <pocket_source>\n{carried.get('pocket_source', 'unguided')}\n\n"
        if carried.get("pocket_id"):
            block += f">  <pocket_id>\n{carried['pocket_id']}\n\n"
        if carried.get("clamp_variant"):
            block += f">  <clamp_variant>\n{carried['clamp_variant']}\n\n"
        block += f">  <refine_variant>\n{refine_variant}\n\n"
        if affinity is not None:
            block += f">  <{refine_variant}_affinity>\n{affinity:.3f}\n\n"
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
        pass  # best effort: a tagging failure must never lose the pose


# ── Discovery ────────────────────────────────────────────────────────────────
def _prepared_receptor(prep_parent: Path, protein_name: str) -> Optional[Path]:
    """The full prepared receptor PDB the refinement docks against.

    pipeline._emit_variants refines against ``job.prepared_protein`` -- the
    whole prepared protein, never a pocket crop -- so the ``_crop_`` siblings
    are excluded here.
    """
    prep = prep_parent / f"_prep_{protein_name}"
    if not prep.is_dir():
        return None
    pdbs = [p for p in prep.glob("*.pdb") if "_crop_" not in p.name]
    return pdbs[0] if len(pdbs) == 1 else None


def _detect_layout(source_root: Path) -> str:
    """'flat' when combo dirs sit at the root, 'nested' when under a complex dir.

    The calibration-benchmark tree is nested (one directory per complex, each
    holding its own ``<lig>__<prot>/`` and ``_prep_*/``). Both Orai trees are
    flat: every ``<lig>__<prot>/`` combo and every ``_prep_*/`` is a direct
    child of the run root.
    """
    for d in source_root.iterdir():
        if d.is_dir() and "__" in d.name and not d.name.startswith("_"):
            return "flat"
    return "nested"


def _pose_mode(filename: str) -> str:
    """Site mode from the pose file name: fpocket | p2rank | unguided."""
    for mode in ("fpocket", "p2rank", "unguided"):
        if filename.startswith(mode):
            return mode
    return "other"


def discover(source_root: Path, layout: str = "auto",
             modes: Optional[set] = None) -> Tuple[List[dict], List[str]]:
    """Every raw pose in ``source_root`` paired with its prepared receptor.

    ``modes`` restricts to given site modes. The calibration-benchmark results
    table reports the smina/gnina refined variants for the UNGUIDED poses only
    (pool 9,090 = 303 x 30); its fpocket and P2Rank rows are raw-only. So a
    re-refinement aimed at the RQ2 answer needs `--modes unguided`, which is a
    third of the tree.
    """
    if layout == "auto":
        layout = _detect_layout(source_root)
    jobs: List[dict] = []
    problems: List[str] = []

    if layout == "flat":
        groups = [(source_root, source_root)]
    else:
        groups = [(d, d) for d in sorted(p for p in source_root.iterdir() if p.is_dir())]

    for container, prep_parent in groups:
        pose_dirs = [d for d in sorted(container.iterdir())
                     if d.is_dir() and "__" in d.name and not d.name.startswith("_")]
        for pose_dir in pose_dirs:
            protein_name = pose_dir.name.split("__", 1)[1]
            receptor = _prepared_receptor(prep_parent, protein_name)
            if receptor is None:
                problems.append(f"no unique prepared receptor for {pose_dir}")
                continue
            raws = sorted(pose_dir.glob(f"*{RAW_SUFFIX}"))
            if not raws:
                problems.append(f"no {RAW_SUFFIX} poses in {pose_dir}")
                continue
            if modes:
                raws = [r for r in raws if _pose_mode(r.name) in modes]
            for raw in raws:
                jobs.append({
                    "raw": str(raw),
                    "receptor": str(receptor),
                    "rel_dir": str(pose_dir.relative_to(source_root)),
                    "base": raw.name[:-len(RAW_SUFFIX)],
                    "complex": (pose_dir.name if layout == "flat"
                                else container.name),
                })
    return jobs, problems


# ── Worker ───────────────────────────────────────────────────────────────────
def _run_one(task: dict) -> List[dict]:
    raw = Path(task["raw"])
    receptor = Path(task["receptor"])
    dest_dir = Path(task["dest_dir"])
    search = task["search"]
    tools = task["tools"]
    force = task["force"]

    dest_dir.mkdir(parents=True, exist_ok=True)
    carried = _read_tags(raw)
    rows: List[dict] = []

    # Mirror the raw pose so the destination tree is self-contained.
    raw_copy = dest_dir / raw.name
    if force or not (raw_copy.exists() and raw_copy.stat().st_size > 0):
        shutil.copy2(raw, raw_copy)

    for tool in tools:
        target = dest_dir / f"{task['base']}__ref{tool.upper()}.sdf"
        row = {
            "complex": task["complex"], "rel_dir": task["rel_dir"],
            "pose": task["base"], "tool": tool, "search": search,
            "status": "", "minimized_affinity": "", "cnn_score": "",
            "cnn_affinity": "", "elapsed_s": "", "error": "",
        }
        if not force and target.exists() and target.stat().st_size > 0:
            row["status"] = "skipped_exists"
            rows.append(row)
            continue

        t0 = time.time()
        ok, scores, msg = refine_pose(raw, receptor, target, tool=tool, search=search)
        elapsed = time.time() - t0
        row["elapsed_s"] = f"{elapsed:.3f}"

        if not ok or not (target.exists() and target.stat().st_size > 0):
            # Same fallback as pipeline._emit_variants: never lose a pose, so
            # the refined pose count matches the source arm exactly.
            shutil.copy2(raw, target)
            row["status"] = "refine_failed_kept_raw"
            row["error"] = msg
            _tag_pose(target, carried, tool, None, None, None)
        else:
            aff = scores.get("minimized_affinity")
            cnn_s = scores.get("cnn_score")
            cnn_a = scores.get("cnn_affinity")
            row["status"] = "ok"
            row["minimized_affinity"] = "" if aff is None else f"{aff:.4f}"
            row["cnn_score"] = "" if cnn_s is None else f"{cnn_s:.4f}"
            row["cnn_affinity"] = "" if cnn_a is None else f"{cnn_a:.4f}"
            _tag_pose(target, carried, tool, aff, cnn_s, cnn_a)
        rows.append(row)
    return rows


# ── Main ─────────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-root", required=True, type=Path,
                    help="Existing EquiBind output tree carrying *__refRAW.sdf poses.")
    ap.add_argument("--dest-root", required=True, type=Path,
                    help="New tree for the re-refined poses (never the source).")
    ap.add_argument("--search", default="minimize", choices=("minimize", "local_only"))
    ap.add_argument("--tools", default="gnina",
                    help="Comma-separated: gnina, smina, or gnina,smina.")
    # Near the core count on purpose: with the thread pins above, each worker is
    # single-threaded, so a low worker count leaves most of the machine idle.
    ap.add_argument("--workers", type=int, default=28)
    ap.add_argument("--layout", default="auto", choices=("auto", "flat", "nested"),
                    help="Tree shape; 'auto' detects it from the source root.")
    ap.add_argument("--modes", default="",
                    help="Comma-separated site modes to refine (fpocket, p2rank, "
                         "unguided). Empty = all. The RQ2 EquiBind figures are "
                         "unguided-only, so '--modes unguided' is the RQ2 subset.")
    ap.add_argument("--limit-complexes", type=int, default=0,
                    help="Process only the first N complexes (smoke test).")
    ap.add_argument("--force", action="store_true",
                    help="Recompute poses whose output already exists.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    source_root = args.source_root.resolve()
    dest_root = args.dest_root.resolve()
    if dest_root == source_root:
        print("ERROR: --dest-root must differ from --source-root", file=sys.stderr)
        return 2
    tools = [t.strip().lower() for t in args.tools.split(",") if t.strip()]
    if any(t not in ("gnina", "smina") for t in tools):
        print(f"ERROR: --tools must be gnina and/or smina (got {tools})", file=sys.stderr)
        return 2

    exes = {"gnina": CFG.gnina_executable, "smina": CFG.smina_executable}
    for t in tools:
        if not exes[t] or not Path(exes[t]).exists():
            print(f"ERROR: {t} executable not resolved (got {exes[t]!r})", file=sys.stderr)
            return 2

    print(f"source      : {source_root}")
    print(f"dest        : {dest_root}")
    print(f"search mode : --{args.search}")
    print(f"tools       : {', '.join(tools)}  ({', '.join(exes[t] for t in tools)})")
    print(f"receptor    : prepared protein PDB (ext {CFG.smina_receptor_ext})")
    print(f"params      : autobox_add={CFG.smina_autobox_add} cpu={CFG.smina_cpu} "
          f"seed={CFG.smina_seed} timeout={CFG.smina_timeout_s}s gpu={CFG.gnina_use_gpu}")

    layout = _detect_layout(source_root) if args.layout == "auto" else args.layout
    modes = {m.strip().lower() for m in args.modes.split(",") if m.strip()}
    print(f"layout      : {layout}")
    print(f"site modes  : {', '.join(sorted(modes)) if modes else 'all'}")
    jobs, problems = discover(source_root, layout, modes or None)
    if args.limit_complexes:
        keep = sorted({j["complex"] for j in jobs})[:args.limit_complexes]
        jobs = [j for j in jobs if j["complex"] in set(keep)]
    n_complexes = len({j["complex"] for j in jobs})
    print(f"discovered  : {len(jobs)} raw poses across {n_complexes} complexes "
          f"-> {len(jobs) * len(tools)} refinements")
    if problems:
        print(f"WARNING     : {len(problems)} discovery problems; first 5:")
        for p in problems[:5]:
            print(f"              {p}")
    if args.dry_run or not jobs:
        return 0

    for j in jobs:
        j["dest_dir"] = str(dest_root / j["rel_dir"])
        j["search"] = args.search
        j["tools"] = tools
        j["force"] = args.force

    dest_root.mkdir(parents=True, exist_ok=True)
    manifest_path = dest_root / f"refine_manifest_{args.search}.csv"
    fields = ["complex", "rel_dir", "pose", "tool", "search", "status",
              "minimized_affinity", "cnn_score", "cnn_affinity", "elapsed_s", "error"]

    t_start = time.time()
    counts: Dict[str, int] = {}
    done = 0
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(_run_one, j): j for j in jobs}
            for fut in as_completed(futures):
                try:
                    rows = fut.result()
                except Exception as e:  # pragma: no cover - defensive
                    j = futures[fut]
                    rows = [{"complex": j["complex"], "rel_dir": j["rel_dir"],
                             "pose": j["base"], "tool": t, "search": args.search,
                             "status": "worker_crash", "minimized_affinity": "",
                             "cnn_score": "", "cnn_affinity": "", "elapsed_s": "",
                             "error": repr(e)} for t in tools]
                for r in rows:
                    writer.writerow(r)
                    counts[r["status"]] = counts.get(r["status"], 0) + 1
                done += 1
                if done % 250 == 0 or done == len(jobs):
                    el = time.time() - t_start
                    rate = done / el if el else 0
                    eta = (len(jobs) - done) / rate if rate else 0
                    print(f"  {done}/{len(jobs)} poses  {el/60:.1f} min elapsed  "
                          f"{rate:.1f} poses/s  ETA {eta/60:.1f} min", flush=True)
                    fh.flush()

    elapsed = time.time() - t_start
    summary = {
        "source_root": str(source_root),
        "dest_root": str(dest_root),
        "search": args.search,
        "tools": tools,
        "executables": {t: exes[t] for t in tools},
        "refine_params": {
            "autobox_add": CFG.smina_autobox_add,
            "cpu": CFG.smina_cpu,
            "seed": CFG.smina_seed,
            "timeout_s": CFG.smina_timeout_s,
            "receptor_ext": CFG.smina_receptor_ext,
            "gnina_use_gpu": CFG.gnina_use_gpu,
        },
        "n_complexes": n_complexes,
        "n_raw_poses": len(jobs),
        "n_refinements": len(jobs) * len(tools),
        "status_counts": counts,
        "discovery_problems": problems,
        "wall_time_s": round(elapsed, 1),
        "manifest": str(manifest_path),
    }
    summary_path = dest_root / f"rerun_summary_{args.search}.json"
    summary_path.write_text(json.dumps(summary, indent=2))

    print(f"\ndone in {elapsed/60:.1f} min")
    for k, v in sorted(counts.items()):
        print(f"  {k:26s} {v}")
    print(f"manifest: {manifest_path}")
    print(f"summary : {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
