#!/usr/bin/env python3
"""Stage the crystal-boxed AutoDock arms of the published 308 for PoseBusters.

Two problems have to be solved before those poses can be busted at all.

1. Arm collision. ``_normalize_protein`` strips the ``_meeko`` / ``_mgl_tools``
   suffix (run_posebusters.py:800), so a tree holding both converter arms emits
   two rows per complex under one protein key. Each arm therefore gets its own
   staging root and its own output_base_dir.

2. No parity-correct validation receptor. ``_resolve_autodock_receptor`` demands
   ``<stem>.pdb`` beside ``<stem>.pdbqt`` and asserts the two hold the same heavy
   atoms. In this campaign that holds for only 142/308 mgl_tools and 70/297
   meeko receptors: the cofactor/metal HETATM block was injected into the PDBQT
   after the PDB was written, and the 151 complexes docked against a
   ``_meeko_retry1`` receptor have no matching PDB at all -- only the fixed PDB
   fed to the failed first attempt.

   The receptor is rebuilt here so that it represents exactly the heavy atoms
   Vina searched: every heavy atom of the PDBQT is emitted once, reusing the
   fixed PDB's own record for it wherever the two agree (which keeps that file's
   hydrogens and formatting, matching the convention of the whole-protein arm's
   receptors) and rendering the rest -- the injected cofactors and metals --
   from the PDBQT line. Hydrogens are carried over for any residue that keeps at
   least one heavy atom. Nothing is written unless the result passes the
   runner's own ``_assert_vina_receptor_parity``.

The pose files are hard-linked, not symlinked: the runner resolves a pose path
before deriving the complex directory from it, so a symlink would point
provenance back at Dockings/Benchmark and fail the root check. Dockings/Benchmark
itself is never modified.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_posebusters as R  # noqa: E402

CONVERTERS = ("meeko", "mgl_tools")


def _pdbqt_key(line: str) -> tuple:
    fields = line.split()
    element = R._AD4_ELEMENT.get(fields[-1] if fields else "")
    if element is None:
        raise ValueError(f"unknown AutoDock atom type in {line!r}")
    coords = tuple(int(round(float(line[a:b]) * 1000)) for a, b in ((30, 38), (38, 46), (46, 54)))
    return (line[:6].strip(), line[17:20].strip(), element, *coords)


def _pdb_key(line: str) -> tuple:
    element = line[76:78].strip().upper()
    if not element:
        element = "".join(c for c in line[12:16] if c.isalpha()).upper()[:1]
    coords = tuple(int(round(float(line[a:b]) * 1000)) for a, b in ((30, 38), (38, 46), (46, 54)))
    return (line[:6].strip(), line[17:20].strip(), element, *coords)


def _is_h_pdbqt(line: str) -> bool:
    fields = line.split()
    return R._AD4_ELEMENT.get(fields[-1] if fields else "") == "H"


def _is_h_pdb(line: str) -> bool:
    element = line[76:78].strip().upper()
    if not element:
        element = "".join(c for c in line[12:16] if c.isalpha()).upper()[:1]
    return element == "H"


def _residue(line: str) -> tuple[str, str, str]:
    return (line[21], line[22:27], line[17:20])


def _pdbqt_line_to_pdb(line: str) -> str:
    """Render a PDBQT atom record as a PDB record with a correct element field."""
    fields = line.split()
    element = R._AD4_ELEMENT[fields[-1]]
    body = line[:66].ljust(66)
    return f"{body}{'':10}{element:>2}  "


def build_receptor_pdb(pdbqt: Path, fixed_pdb: Path | None) -> str:
    """Return PDB text whose heavy atoms are exactly those of *pdbqt*.

    Record order follows the fixed PDB, so the result is that file with the
    residues the converter dropped removed; the cofactors and metals injected
    into the PDBQT after it was written are appended, which is where the
    injection put them in the PDBQT too.
    """
    wanted: Counter = Counter()
    pdbqt_lines: dict[tuple, str] = {}
    for line in pdbqt.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")) or _is_h_pdbqt(line):
            continue
        key = _pdbqt_key(line)
        wanted[key] += 1
        pdbqt_lines.setdefault(key, line)

    source = (fixed_pdb.read_text(errors="replace").splitlines()
              if fixed_pdb is not None and fixed_pdb.is_file() else [])

    # Pass 1: which residues keep at least one heavy atom, so their hydrogens
    # travel with them.
    remaining = Counter(wanted)
    kept_residues: set[tuple[str, str, str]] = set()
    for line in source:
        if not line.startswith(("ATOM", "HETATM")) or _is_h_pdb(line):
            continue
        key = _pdb_key(line)
        if remaining.get(key, 0) > 0:
            remaining[key] -= 1
            kept_residues.add(_residue(line))

    # Pass 2: re-emit in the source order.
    remaining = Counter(wanted)
    kept_lines: list[str] = []
    for line in source:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if _is_h_pdb(line):
            if _residue(line) in kept_residues:
                kept_lines.append(line)
            continue
        key = _pdb_key(line)
        if remaining.get(key, 0) > 0:
            remaining[key] -= 1
            kept_lines.append(line)

    for key, count in remaining.items():
        if count > 0:
            kept_lines.extend([_pdbqt_line_to_pdb(pdbqt_lines[key])] * count)

    return "\n".join(kept_lines) + "\nEND\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", type=Path, default=Path("Dockings/Benchmark"))
    ap.add_argument("--ids-file", type=Path,
                    default=Path("Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"))
    ap.add_argument("--out", type=Path,
                    default=Path("posebusters_results/_boxed308_staging"))
    ap.add_argument("--fresh", action="store_true", help="delete the staging root first")
    args = ap.parse_args()

    ids = [ln.strip() for ln in args.ids_file.read_text().splitlines() if ln.strip()]
    if args.fresh and args.out.exists():
        shutil.rmtree(args.out)

    failures: list[tuple[str, str, str]] = []
    for converter in CONVERTERS:
        staged = reused = 0
        for pdb_id in ids:
            src_dir = args.tree / pdb_id / converter / "docking"
            poses = sorted(src_dir.glob("*_vina_out.pdbqt")) if src_dir.is_dir() else []
            if not poses:
                continue
            pose = poses[0]
            stem = pose.stem[: -len("_vina_out")].split("__", 1)[0]
            stage_dir = args.out / converter / pdb_id / converter / "docking"
            rec_dir = args.out / converter / pdb_id / "_staging" / "receptors" / "pdbqt"
            stage_dir.mkdir(parents=True, exist_ok=True)
            rec_dir.mkdir(parents=True, exist_ok=True)

            target = stage_dir / pose.name
            if not target.exists():
                os.link(pose, target)

            src_pdbqt = args.tree / pdb_id / "_staging" / "receptors" / "pdbqt" / f"{stem}.pdbqt"
            if not src_pdbqt.is_file():
                failures.append((converter, pdb_id, f"no prepared PDBQT {src_pdbqt.name}"))
                continue
            dst_pdbqt = rec_dir / f"{stem}.pdbqt"
            if not dst_pdbqt.exists():
                os.link(src_pdbqt, dst_pdbqt)

            dst_pdb = rec_dir / f"{stem}.pdb"
            if dst_pdb.is_file():
                reused += 1
                continue
            src_pdb = args.tree / pdb_id / "_staging" / "receptors" / "pdbqt" / f"{stem}.pdb"
            if not src_pdb.is_file():
                # the _meeko_retry1 arm: fall back to the fixed PDB of the first attempt
                src_pdb = (args.tree / pdb_id / "_staging" / "receptors" / "pdbqt"
                           / f"{stem.replace('_retry1', '')}.pdb")
            try:
                text = build_receptor_pdb(dst_pdbqt, src_pdb if src_pdb.is_file() else None)
                dst_pdb.write_text(text)
                R._assert_vina_receptor_parity.__wrapped__(str(dst_pdb), str(dst_pdbqt))
            except Exception as exc:
                dst_pdb.unlink(missing_ok=True)
                failures.append((converter, pdb_id, f"{type(exc).__name__}: {exc}"))
                continue
            staged += 1
        print(f"{converter:10s} receptors rebuilt {staged:4d}, already present {reused:4d}")

    print(f"failures: {len(failures)}")
    for converter, pdb_id, why in failures[:20]:
        print(f"  FAIL {converter} {pdb_id}: {why}")
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
