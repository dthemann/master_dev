#!/usr/bin/env python3
"""Prepare the EquiBind --minimize tree so run_posebusters.py can consume it.

The re-refinement harness writes poses only. The PoseBusters runner additionally
needs three things that live outside the pose files, and it fails in three
different ways without them:

  1. ``_prep_<protein>/`` beside each complex. ``_resolve_equibind_receptor``
     builds both the full prepared receptor and, for raw guided poses, the
     ``_crop_<pocket_id>.pdb`` variant from this directory. Missing ->
     "Expected exactly one full prepared EquiBind receptor ...; found 0".
     Symlinked, not copied, so the receptor sha256 matches the published arm
     byte for byte and the two runs stay comparable.

  2. ``guided_results.json`` per combo. ``_equibind_finalized_sdfs`` reads only
     this manifest; without it the combo is skipped as "incomplete" and the run
     collects zero poses. It is DERIVED from the published manifest rather than
     synthesized, so pose identity, pocket provenance (``mode``, ``pocket_id``,
     ``pocket_unique_id``, ``pocket_center``) and UFF provenance are inherited
     verbatim. Only what genuinely changed is rewritten: the pose path, the
     recomputed centroid, the refinement scores and the refinement outcome.
     A synthesized manifest that dropped ``pocket_id`` would silently break the
     raw variant's crop-receptor resolution.

  3. A flattened staging farm. ``_collect_equibind_rows`` splits each top-level
     directory name on "__" and requires exactly two parts, so it cannot see the
     real tree's ``<COMPLEX>/<combo>/`` nesting. The published arm is staged by a
     notebook cell whose source root is a hardcoded literal, so it cannot stage
     this tree; this script reproduces the same naming
     (``<CPLX>_ligand__<CPLX>_protein``) so ``pose_name`` matches the published
     CSV and the two tables join directly.

Refinement success is read from the pose file itself (the ``<tool>_affinity`` tag,
which the harness writes only when the refiner returned a score) rather than from
the manifest CSV, whose ``skipped_exists`` status is ambiguous across resumed passes.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

SUFFIX_TO_VARIANT = {
    "__refRAW.sdf": "raw",
    "__refSMINA.sdf": "smina",
    "__refGNINA.sdf": "gnina",
}


def _read_tags(sdf: Path) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    lines = sdf.read_text(errors="ignore").splitlines()
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith("> ") and "<" in s and ">" in s[2:]:
            key = s[s.index("<") + 1:s.rindex(">")]
            if i + 1 < len(lines):
                tags[key] = lines[i + 1].strip()  # last wins, matching the runner
    return tags


def _centroid(sdf: Path) -> Optional[List[float]]:
    try:
        lines = sdf.read_text(errors="ignore").splitlines()
        n = int(lines[3][0:3])
        pts = [(float(l[0:10]), float(l[10:20]), float(l[20:30]))
               for l in lines[4:4 + n]]
        if not pts:
            return None
        return [sum(p[i] for p in pts) / len(pts) for i in range(3)]
    except Exception:
        return None


def _variant_of(name: str) -> Optional[str]:
    for suf, var in SUFFIX_TO_VARIANT.items():
        if name.endswith(suf):
            return var
    return None


def _detect_layout(root: Path) -> str:
    """'flat' when combo dirs sit at the root, 'nested' when under a complex dir.

    The calibration-benchmark tree is nested; both Orai trees are flat.
    """
    for d in root.iterdir():
        if d.is_dir() and "__" in d.name and not d.name.startswith("_"):
            return "flat"
    return "nested"


def _combo_pairs(published: Path, minimize: Path, layout: str) -> List[Tuple[Path, Path]]:
    """(published_combo, minimize_combo) for every combo in the minimize tree."""
    pairs: List[Tuple[Path, Path]] = []
    if layout == "flat":
        for combo in sorted(d for d in minimize.iterdir()
                            if d.is_dir() and "__" in d.name and not d.name.startswith("_")):
            pairs.append((published / combo.name, combo))
    else:
        for cplx in sorted(p for p in minimize.iterdir() if p.is_dir()):
            for combo in sorted(d for d in cplx.iterdir()
                                if d.is_dir() and "__" in d.name and not d.name.startswith("_")):
                pairs.append((published / cplx.name / combo.name, combo))
    return pairs


def _pose_mode(filename: str) -> str:
    for mode in ("fpocket", "p2rank", "unguided"):
        if filename.startswith(mode):
            return mode
    return "other"


def build_manifests(published: Path, minimize: Path, layout: str,
                    modes: Optional[set] = None) -> Tuple[int, int, int, List[str]]:
    """Derive each combo's guided_results.json for the minimize tree.

    ``modes`` restricts the manifest to given site modes. This matters because
    the manifest is the ONLY thing the runner collects poses through: a tree that
    was refined for one site mode will still contain stray files from an earlier
    partial pass, and listing those would bust a mode with incomplete complex
    coverage. Restricting the manifest excludes them without deleting anything.
    """
    combos = n_rec = n_fail = 0
    problems: List[str] = []

    for src_combo, combo in _combo_pairs(published, minimize, layout):
            src_manifest = src_combo / "guided_results.json"
            if not src_manifest.is_file():
                problems.append(f"no source manifest at {src_manifest}")
                continue
            records = json.loads(src_manifest.read_text())
            out: List[dict] = []
            for rec in records:
                old = Path(rec["sdf_path"])
                if modes and _pose_mode(old.name) not in modes:
                    continue
                new = combo / old.name
                if not new.is_file() or new.stat().st_size == 0:
                    problems.append(f"missing pose {new}")
                    continue
                variant = _variant_of(old.name) or rec.get("refine_variant")
                tags = _read_tags(new)

                r = dict(rec)
                r["sdf_path"] = str(new.resolve())
                r["refine_variant"] = variant
                cen = _centroid(new)
                if cen is not None:
                    r["pose_centroid"] = cen

                if variant == "raw":
                    r["refine_affinity"] = None
                    r["cnn_score"] = None
                    r["cnn_affinity"] = None
                    r["success"] = True
                    r["error"] = ""
                else:
                    aff = tags.get(f"{variant}_affinity")
                    ok = aff is not None
                    r["refine_affinity"] = float(aff) if ok else None
                    r["cnn_score"] = float(tags["cnn_score"]) if "cnn_score" in tags else None
                    r["cnn_affinity"] = (float(tags["cnn_affinity"])
                                         if "cnn_affinity" in tags else None)
                    # The harness keeps the pre-refine pose under the refined
                    # label when a refiner fails, exactly as the pipeline does,
                    # so the pose is committed either way; only `error` differs.
                    r["success"] = True
                    r["error"] = "" if ok else "refine_failed: no affinity tag"
                    if not ok:
                        n_fail += 1
                out.append(r)
                n_rec += 1
            (combo / "guided_results.json").write_text(json.dumps(out, indent=1))
            combos += 1
    return combos, n_rec, n_fail, problems


def link_prep_dirs(published: Path, minimize: Path) -> Tuple[int, List[str]]:
    linked = 0
    problems: List[str] = []
    for cplx in sorted(p for p in minimize.iterdir() if p.is_dir()):
        src_cplx = published / cplx.name
        if not src_cplx.is_dir():
            problems.append(f"no published complex dir for {cplx.name}")
            continue
        preps = [d for d in src_cplx.iterdir() if d.is_dir() and d.name.startswith("_prep_")]
        if not preps:
            problems.append(f"no _prep_* under {src_cplx}")
            continue
        for prep in preps:
            target = cplx / prep.name
            if target.is_symlink() or target.exists():
                continue
            target.symlink_to(prep.resolve(), target_is_directory=True)
            linked += 1
    return linked, problems


def build_staging(minimize: Path, staging: Path) -> Tuple[int, List[str]]:
    staging.mkdir(parents=True, exist_ok=True)
    made = 0
    problems: List[str] = []
    for cplx in sorted(p for p in minimize.iterdir() if p.is_dir()):
        combos = [d for d in cplx.iterdir()
                  if d.is_dir() and "__" in d.name and not d.name.startswith("_")]
        if len(combos) != 1:
            problems.append(f"{cplx.name}: expected 1 combo dir, found {len(combos)}")
            continue
        # Reproduce the published staged name so pose_name matches that CSV.
        link = staging / f"{cplx.name}_ligand__{cplx.name}_protein"
        if link.is_symlink() or link.exists():
            continue
        link.symlink_to(combos[0].resolve(), target_is_directory=True)
        made += 1
    return made, problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--published", type=Path, default=Path("Dockings/Benchmark_Equibind"))
    ap.add_argument("--minimize", type=Path, default=Path("Dockings/Benchmark_Equibind_minimize"))
    ap.add_argument("--staging", type=Path,
                    default=Path("posebusters_results/_equibind_minimize_staging/equibind"))
    ap.add_argument("--expect-complexes", type=int, default=308)
    ap.add_argument("--layout", default="auto", choices=("auto", "flat", "nested"))
    # The Orai configs resolve receptors from a staged receptors_dir and read the
    # pose tree directly, so neither _prep_ links nor a staging farm applies there.
    ap.add_argument("--manifests-only", action="store_true",
                    help="Only derive guided_results.json (use for the flat Orai trees).")
    ap.add_argument("--modes", default="",
                    help="Comma-separated site modes to list in the manifests "
                         "(fpocket, p2rank, unguided). Empty = all.")
    args = ap.parse_args()

    published, minimize, staging = args.published, args.minimize, args.staging
    if not minimize.is_dir():
        print(f"ERROR: {minimize} does not exist", file=sys.stderr)
        return 2

    layout = _detect_layout(minimize) if args.layout == "auto" else args.layout
    n_units = len([p for p in minimize.iterdir() if p.is_dir()])
    print(f"minimize tree : {minimize}")
    print(f"layout        : {layout}  ({n_units} top-level dirs)")
    if layout == "nested" and n_units != args.expect_complexes:
        print(f"ERROR: expected {args.expect_complexes} complexes, found {n_units}. "
              f"Is the refinement still running?", file=sys.stderr)
        return 3

    p1: List[str] = []
    if not args.manifests_only:
        linked, p1 = link_prep_dirs(published, minimize)
        print(f"_prep_ links  : {linked} created")

    modes = {m.strip().lower() for m in args.modes.split(",") if m.strip()}
    print(f"site modes    : {', '.join(sorted(modes)) if modes else 'all'}")
    combos, recs, fails, p2 = build_manifests(published, minimize, layout,
                                              modes or None)
    print(f"manifests     : {combos} combos, {recs} pose records, {fails} refine failures")

    p3: List[str] = []
    if not args.manifests_only:
        made, p3 = build_staging(minimize, staging)
        print(f"staging links : {made} created under {staging}")

    problems = p1 + p2 + p3
    if problems:
        print(f"\nPROBLEMS ({len(problems)}):")
        for p in problems[:30]:
            print(f"  {p}")
        return 1
    print("\nok - ready for run_posebusters.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
