#!/usr/bin/env python3
"""Write the missing ``*_protein_meeko_retry1.pdb`` validation receptors.

151 of the 297 crystal-boxed Meeko complexes were docked against a retried
receptor preparation, ``<id>_protein_meeko_retry1.pdbqt``. The staging directory
kept that PDBQT but not a matching ``.pdb``: the only PDB beside it is
``<id>_protein_meeko.pdb``, the PDBFixer output that was fed to the FIRST,
failed Meeko attempt. PoseBusters resolves an AutoDock pose's validation
receptor from the receptor stem in the pose filename
(run_posebusters.py:_resolve_autodock_receptor) and then asserts that the PDB
and the PDBQT hold exactly the same heavy atoms, so those 151 complexes cannot
be busted at all as things stand.

The retry differs from the first attempt only by dropping whole residues -- not
one residue that survives loses an atom, and no retained heavy atom moves. So
the validation receptor is recovered exactly by filtering the fixed PDB down to
the residues the retry PDBQT kept, which preserves that file's protonation and
record layout (the same convention the mgl_tools arm's ``.pdb`` already uses)
while making its heavy atoms identical to the receptor Vina actually searched.

Writes nothing unless the reconstruction passes the runner's own parity check.
"""
from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import run_posebusters as R  # noqa: E402

_HYDROGEN = {"H", "HD", "HS"}


def _element(line: str, pdbqt: bool) -> str:
    field = (line[77:79] if pdbqt else line[76:78]).strip().upper()
    if field:
        return field
    return "".join(c for c in line[12:16] if c.isalpha()).upper()[:1]


def _residue(line: str) -> tuple[str, str, str]:
    return (line[21], line[22:27], line[17:20])


def _kept_residues(pdbqt: Path) -> set[tuple[str, str, str]]:
    keep = set()
    for line in pdbqt.read_text(errors="replace").splitlines():
        if line.startswith(("ATOM", "HETATM")) and _element(line, True) not in _HYDROGEN:
            keep.add(_residue(line))
    return keep


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", type=Path, default=Path("Dockings/Benchmark"))
    ap.add_argument("--ids-file", type=Path,
                    default=Path("Data/PoseBuster Benchmark Set/posebusters_pdb_ccd_ids.txt"))
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    ids = [ln.strip() for ln in args.ids_file.read_text().splitlines() if ln.strip()]
    written = skipped = 0
    failures: list[tuple[str, str]] = []

    for pdb_id in ids:
        stage = args.tree / pdb_id / "_staging" / "receptors" / "pdbqt"
        pdbqt = stage / f"{pdb_id}_protein_meeko_retry1.pdbqt"
        target = stage / f"{pdb_id}_protein_meeko_retry1.pdb"
        source = stage / f"{pdb_id}_protein_meeko.pdb"
        if not pdbqt.is_file():
            continue
        if target.is_file():
            skipped += 1
            continue
        if not source.is_file():
            failures.append((pdb_id, "no first-attempt PDB to filter"))
            continue

        keep = _kept_residues(pdbqt)
        out_lines, seen = [], collections.Counter()
        for line in source.read_text(errors="replace").splitlines():
            if line.startswith(("ATOM", "HETATM")):
                if _residue(line) not in keep:
                    continue
                seen[_residue(line)] += 1
            elif line.startswith(("TER", "ANISOU")):
                continue
            out_lines.append(line)
        if "END" not in {ln[:3] for ln in out_lines[-3:]}:
            out_lines.append("END")

        missing = keep - set(seen)
        if missing:
            failures.append((pdb_id, f"{len(missing)} retry residues absent from the fixed PDB"))
            continue

        text = "\n".join(out_lines) + "\n"
        if args.dry_run:
            tmp = target.with_suffix(".pdb.check")
            tmp.write_text(text)
            probe = tmp
        else:
            target.write_text(text)
            probe = target
        try:
            R._assert_vina_receptor_parity.__wrapped__(str(probe), str(pdbqt))
        except Exception as exc:
            failures.append((pdb_id, f"parity: {exc}"))
            probe.unlink(missing_ok=True)
            continue
        if args.dry_run:
            probe.unlink(missing_ok=True)
        written += 1

    print(f"reconstructed {written} validation receptors "
          f"({'dry run, nothing kept' if args.dry_run else 'written'}), "
          f"{skipped} already present, {len(failures)} failed")
    for pdb_id, why in failures[:20]:
        print(f"  FAIL {pdb_id}: {why}")
    if failures:
        sys.exit(1)


if __name__ == "__main__":
    main()
