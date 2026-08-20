"""Inject biologically relevant HETATM atoms from a source PDB into a
prepared receptor PDBQT.

Meeko (``mk_prepare_receptor.py``) and MGLTools (``prepare_receptor4.py``)
strip every HETATM during receptor preparation. For docking against
proteins whose binding site is occupied by a real cofactor (HEM, FAD,
NAD/NAP, FMN, SAM, ATP, ANP, ...), this leaves the docker blind to those
atoms and produces poses that overlap with them.

This module appends the kept HETATM atoms (cofactors + metal ions) to
the receptor PDBQT after the protein atoms, using a heuristic element ->
AD4 atom-type mapping. Partial charges are set to 0.0 — adequate for
steric / shape-driven docking with Vina scoring (the dominant use case);
not adequate for AD4 scoring of cofactor binding affinities.

The receptor PDBQT remains a single, valid file that Vina can read.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Curated lists
# ---------------------------------------------------------------------------

DEFAULT_KEEP_COFACTORS: frozenset[str] = frozenset({
    # haem & related
    "HEM", "HEC", "HEA", "HEB", "HAS", "HDD", "HEO",
    # nicotinamides
    "NAD", "NAI", "NAP", "NDP", "NHD",
    # flavins
    "FMN", "FAD", "FAS", "FED", "F42",
    # nucleotides & analogues
    "ATP", "ADP", "AMP", "ANP", "GTP", "GDP", "GMP",
    "TTP", "UTP", "CTP", "CDP", "UDP",
    # CoA family
    "COA", "ACO", "COD", "HSC", "MCA", "MYA",
    # methyl-donors / vitamins / others
    "SAM", "SAH", "PLP", "PMP", "TPP", "TDP", "THF", "MTE", "MGD",
    "OXY", "ORO", "MOH",
    # photosynthesis / chlorophylls / carotenoids
    "BCL", "BCR", "BCT", "CHL", "CLA", "PHO",
    # lipids / fatty acids commonly bound
    "PLM", "OLA", "PAM", "STE", "LMG", "SQD", "DGD", "LHG",
    # iron-sulfur clusters & related metalloclusters
    "SF4", "FES", "F3S", "CFM", "CFN", "NFR", "NFU",
    # misc
    "PCY", "XO7", "SIN", "FPS", "FON", "7DG",
    # active-site cofactors found in the PoseBuster Benchmark set that the
    # docker must respect (added after a HETATM audit of the benchmark PDBs):
    "B12",   # cobalamin / vitamin B12
    "H4B",   # tetrahydrobiopterin
    "GSH",   # glutathione (reduced)
    "COM",   # coenzyme M
    "F43",   # coenzyme F430
    "TYD",   # thymidine-5'-diphosphate
    "TPW",   # thiamine diphosphate analogue
    "S3P",   # shikimate-3-phosphate
    "DTP",   # 2'-deoxyadenosine-5'-triphosphate
    "CP",    # carbamoyl phosphate
})

DEFAULT_KEEP_METALS: frozenset[str] = frozenset({
    "MG", "ZN", "CA", "FE", "FE2", "FE3", "MN", "MN3",
    "CO", "CU", "CU1", "NI", "NA", "K", "LI",
    "CD", "HG", "PT", "AU", "AG",
})

# The 20 standard amino acids (+ common alt names). Residues on ATOM records that
# are NOT in this set are modified/covalent residues (phospho-Ser/Thr/Tyr, seleno-
# Met, carboxy-Lys, PLP-Lys, oxidised Cys, methyl-Lys, ...) or cofactors deposited
# as ATOM records (e.g. UDP, V97). MGLTools' prepare_receptor4.py silently deletes
# these, leaving holes in the docking pocket that the docker then places ligands
# into — producing PoseBusters clashes against atoms that are genuinely part of
# the deposited protein. inject_dropped_atom_residues() re-adds exactly those.
STANDARD_AA: frozenset[str] = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # protonation / tautomer variants some tools emit
    "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN", "ARN", "TYM",
})

# Modified residues that MGLTools drops but which we also must NOT inject, because
# they carry an element AutoDock/AD4 (and the PoseBusters receptor-parity check)
# cannot type. MSE (selenomethionine) is a crystallographic phasing surrogate for
# methionine; its selenium has no AD4 atom type, so we leave it dropped like mgl
# does rather than inject an unscorable Se atom.
_DROP_ONLY_RESIDUES: frozenset[str] = frozenset({
    "MSE",   # selenomethionine (Se)
})

# Element symbol -> AD4 atom type (uppercase keys).
# Defaults assume "everything is a polar acceptor" since most cofactor
# N/O/S atoms in receptors are H-bond acceptors. For carbon, default to "A"
# (aromatic) is too aggressive for chained cofactors; plain "C" is safer.
_ELEMENT_TO_AD4: dict[str, str] = {
    "H": "HD",   # treat as donor H (typical receptor convention)
    "D": "HD",
    "C": "C",
    "N": "NA",   # most cofactor N atoms are H-bond acceptors
    "O": "OA",   # most cofactor O atoms are acceptors
    "S": "SA",
    "P": "P",
    "F": "F",
    "CL": "Cl",
    "BR": "Br",
    "I":  "I",
    # metals (AD4 types use mixed case)
    "MG": "Mg", "ZN": "Zn", "CA": "Ca",
    "FE": "Fe", "MN": "Mn", "CO": "Co", "CU": "Cu", "NI": "Ni",
    "NA": "Na", "K":  "K",
    "CD": "Cd", "HG": "Hg",
}


def element_to_ad4(element: str) -> str:
    """Map an element symbol (1-2 chars, any case) to an AD4 atom type."""
    if not element:
        return "C"
    return _ELEMENT_TO_AD4.get(element.strip().upper(), element.strip().capitalize())


# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

def _infer_element(atom_name: str, fallback_resname: str = "") -> str:
    """Best-effort guess at the element symbol from a 4-char PDB atom name."""
    name = atom_name.strip()
    if not name:
        return "X"
    # PDB convention: cols 77-78 hold the element if filled; otherwise the
    # first non-digit character(s) of the atom name. Handle a few special
    # cases (e.g. CA the alpha-carbon vs CA the calcium ion).
    if fallback_resname.strip().upper() == "CA" and name == "CA":
        return "CA"
    if name[0].isdigit():
        name = name[1:]
    # 2-letter elements (Mg, Fe, etc.) only when the atom name starts with them
    two = name[:2].upper()
    if two in {"MG", "ZN", "CA", "FE", "MN", "CO", "CU", "NI",
               "NA", "CL", "BR", "CD", "HG", "PT", "AU", "AG"}:
        return two
    return name[0].upper()


def extract_hetatm_lines(pdb_path: Path, keep: Iterable[str]) -> list[dict]:
    """Return parsed HETATM atoms from *pdb_path* whose resname is in *keep*."""
    keep_set = {r.strip().upper() for r in keep}
    atoms: list[dict] = []
    for line in pdb_path.read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        resn = line[17:20].strip().upper()
        if resn not in keep_set:
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        # Element: cols 77-78 if present, else infer from atom name
        elem = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not elem:
            elem = _infer_element(line[12:16], resn)
        # Drop hydrogens — receptor HETATMs without explicit polar-H knowledge
        # are safer treated as heavy-atom-only (Vina re-adds polar-H sterics
        # implicitly via atom radii).
        if elem == "H":
            continue
        atoms.append({
            "name":    line[12:16],
            "altloc":  line[16:17] if len(line) > 16 else " ",
            "resname": line[17:20],
            "chain":   line[21:22] if len(line) > 21 else " ",
            "resnum":  line[22:26],
            "icode":   line[26:27] if len(line) > 26 else " ",
            "x": x, "y": y, "z": z,
            "occupancy": _safe_float(line[54:60], 1.00),
            "bfactor":   _safe_float(line[60:66], 0.00),
            "element":   elem,
        })
    return atoms


def _safe_float(s: str, default: float) -> float:
    try:
        return float(s)
    except (ValueError, IndexError):
        return default


# ---------------------------------------------------------------------------
# PDBQT formatting
# ---------------------------------------------------------------------------

def format_pdbqt_atom(atom: dict, serial: int) -> str:
    """Format one parsed HETATM dict as a PDBQT atom line (80 chars + nl)."""
    ad4_type = element_to_ad4(atom["element"])
    # PDBQT atom-name format: 4-character right-pad if 1-2 char, left-aligned
    # to col 14 for 3-4 char names. Keep the source atom name as-is — it's
    # already 4-char padded by the PDB parser.
    name = atom["name"]
    if len(name) < 4:
        name = name.ljust(4)
    line = (
        f"HETATM"
        f"{serial:>5} "
        f"{name:<4}"
        f"{atom['altloc']:1}"
        f"{atom['resname']:>3} "
        f"{atom['chain']:1}"
        f"{atom['resnum']:>4}"
        f"{atom['icode']:1}   "
        f"{atom['x']:>8.3f}"
        f"{atom['y']:>8.3f}"
        f"{atom['z']:>8.3f}"
        f"{atom['occupancy']:>6.2f}"
        f"{atom['bfactor']:>6.2f}"
        f"    "                  # cols 67-70
        f"{0.000:>6.3f} "        # cols 71-77: partial charge (Vina ignores for sterics)
        f"{ad4_type:<2}"         # cols 78-79: AD4 atom type
    )
    return line[:80].ljust(80) + "\n"


def format_pdb_atom(atom: dict, serial: int) -> str:
    """Format one parsed HETATM dict as a plain-PDB atom line (element in 77-78).

    Companion to :func:`format_pdbqt_atom` for the ``.pdb`` twin of a prepared
    receptor. The element goes in columns 77-78 so a heavy-atom-parity check that
    derives the PDBQT element from its AD4 type matches this PDB element exactly
    (the round-trip ``AD4_ELEMENT[element_to_ad4(e)] == e`` holds for the
    C/N/O/S/P/metal atoms that ever get injected here).
    """
    name = atom["name"]
    if len(name) < 4:
        name = name.ljust(4)
    elem = atom["element"].strip().upper()
    line = (
        f"HETATM"
        f"{serial:>5} "
        f"{name:<4}"
        f"{atom['altloc']:1}"
        f"{atom['resname']:>3} "
        f"{atom['chain']:1}"
        f"{atom['resnum']:>4}"
        f"{atom['icode']:1}   "
        f"{atom['x']:>8.3f}"
        f"{atom['y']:>8.3f}"
        f"{atom['z']:>8.3f}"
        f"{atom['occupancy']:>6.2f}"
        f"{atom['bfactor']:>6.2f}"
        f"          "            # cols 67-76
        f"{elem:>2}"             # cols 77-78: element symbol
    )
    return line[:80].ljust(80) + "\n"


def _next_serial(pdbqt_lines: list[str]) -> int:
    """Return the next free atom serial in *pdbqt_lines* (1 if none)."""
    last = 0
    for L in pdbqt_lines:
        if L.startswith(("ATOM", "HETATM")):
            try:
                last = max(last, int(L[6:11]))
            except ValueError:
                pass
    return last + 1


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------

def inject_hetatms(
    source_pdb: Path,
    receptor_pdbqt: Path,
    output_pdbqt: Path | None = None,
    keep_residues: Iterable[str] | None = None,
) -> dict:
    """Inject KEEP HETATMs from *source_pdb* into *receptor_pdbqt*.

    Writes to *output_pdbqt* (defaults to overwriting *receptor_pdbqt* in place).
    Returns a dict summarising the operation.
    """
    if keep_residues is None:
        keep_residues = DEFAULT_KEEP_COFACTORS | DEFAULT_KEEP_METALS

    src = Path(source_pdb)
    rec = Path(receptor_pdbqt)
    out = Path(output_pdbqt) if output_pdbqt else rec

    pdbqt_text = rec.read_text()
    lines = pdbqt_text.splitlines(keepends=True)

    # Strip any END / TER at the very end; we'll re-add END after append.
    while lines and lines[-1].strip() in ("END", "ENDMDL"):
        lines.pop()

    atoms = extract_hetatm_lines(src, keep_residues)
    if not atoms:
        out.write_text("".join(lines) + "END\n")
        return {"source": str(src), "receptor": str(rec), "output": str(out),
                "atoms_injected": 0, "residues_injected": [], "by_resname": {}}

    serial = _next_serial(lines)
    by_resname: dict[str, int] = {}
    new_lines: list[str] = []
    for a in atoms:
        new_lines.append(format_pdbqt_atom(a, serial))
        serial += 1
        rn = a["resname"].strip().upper()
        by_resname[rn] = by_resname.get(rn, 0) + 1

    out.write_text("".join(lines) + "".join(new_lines) + "END\n")
    return {
        "source": str(src),
        "receptor": str(rec),
        "output": str(out),
        "atoms_injected": len(atoms),
        "residues_injected": sorted(by_resname),
        "by_resname": by_resname,
    }


# ---------------------------------------------------------------------------
# Dropped modified-residue / ATOM-cofactor injection
# ---------------------------------------------------------------------------

def _coord_index(pdbqt_path: Path, cell: float = 1.0) -> dict[tuple[int, int, int], list[tuple[float, float, float]]]:
    """Spatial hash of the prepared receptor's heavy-atom coordinates.

    Bins heavy ATOM/HETATM coordinates into *cell*-Å cubes so a source atom's
    presence can be checked against only its own bin + 26 neighbours.
    """
    idx: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
    for line in pdbqt_path.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        # skip hydrogens (prepared receptors carry polar-H; match on heavy atoms)
        name = line[12:16].strip()
        elem = line[76:78].strip().upper() if len(line) >= 78 else ""
        if elem == "H" or (not elem and name[:1] == "H"):
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        key = (int(x // cell), int(y // cell), int(z // cell))
        idx.setdefault(key, []).append((x, y, z))
    return idx


def _present(idx: dict, pt: tuple[float, float, float], tol: float, cell: float = 1.0) -> bool:
    """True if any indexed atom lies within *tol* Å of *pt*."""
    cx, cy, cz = int(pt[0] // cell), int(pt[1] // cell), int(pt[2] // cell)
    t2 = tol * tol
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for (x, y, z) in idx.get((cx + dx, cy + dy, cz + dz), ()):  # noqa
                    if (x - pt[0]) ** 2 + (y - pt[1]) ** 2 + (z - pt[2]) ** 2 <= t2:
                        return True
    return False


def extract_dropped_atom_residues(
    source_pdb: Path,
    prepared_pdbqt: Path,
    min_res_atoms: int = 3,
    tol: float = 0.5,
) -> list[dict]:
    """Return heavy ATOM records from *source_pdb* that MGLTools deleted.

    For every ATOM-record residue whose name is non-standard (see
    :data:`STANDARD_AA`) and which has at least *min_res_atoms* heavy atoms, any
    heavy atom NOT found within *tol* Å in *prepared_pdbqt* is returned. This
    covers both cases MGLTools produces:

    * a wholly deleted residue / ATOM-cofactor (e.g. UDP, V97) → all its atoms
      are missing and re-added;
    * a modified residue kept under its renamed parent but with the modifying
      atoms stripped (e.g. CSD→CYS drops the sulfinic oxygens, LLP→LYS drops the
      pyridoxal-phosphate, SEP→SER drops the phosphate) → only the missing
      modifying atoms are re-added.

    Atoms MGLTools kept (coordinates preserved) are never re-added, so the
    operation cannot double-inject. HETATM records are intentionally ignored;
    those are handled by the cofactor/metal keep-list in :func:`inject_hetatms`.
    """
    idx = _coord_index(prepared_pdbqt)

    # group source ATOM heavy atoms by residue
    residues: dict[tuple, list[dict]] = {}
    order: list[tuple] = []
    for line in source_pdb.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        resn = line[17:20].strip().upper()
        if resn in STANDARD_AA or resn in _DROP_ONLY_RESIDUES:
            continue
        try:
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        except (ValueError, IndexError):
            continue
        elem = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not elem:
            elem = _infer_element(line[12:16], resn)
        if elem == "H":
            continue
        key = (line[21:22], line[22:26], line[26:27], resn)  # chain, resnum, icode, resname
        if key not in residues:
            residues[key] = []
            order.append(key)
        residues[key].append({
            "name":    line[12:16],
            "altloc":  line[16:17] if len(line) > 16 else " ",
            "resname": line[17:20],
            "chain":   line[21:22] if len(line) > 21 else " ",
            "resnum":  line[22:26],
            "icode":   line[26:27] if len(line) > 26 else " ",
            "x": x, "y": y, "z": z,
            "occupancy": _safe_float(line[54:60], 1.00),
            "bfactor":   _safe_float(line[60:66], 0.00),
            "element":   elem,
        })

    dropped: list[dict] = []
    for key in order:
        atoms = residues[key]
        if len(atoms) < min_res_atoms:
            continue
        # re-add every heavy atom MGLTools left out — whether the whole residue
        # was deleted or only its modifying atoms were stripped from a renamed parent
        dropped.extend(a for a in atoms if not _present(idx, (a["x"], a["y"], a["z"]), tol))
    return dropped


def inject_dropped_atom_residues(
    source_pdb: Path,
    receptor_pdbqt: Path,
    output_pdbqt: Path | None = None,
    min_res_atoms: int = 3,
    tol: float = 0.5,
) -> dict:
    """Append MGLTools-deleted modified residues / ATOM-cofactors to a PDBQT.

    Mirrors :func:`inject_hetatms` but for ATOM-record species. Idempotent per
    call in the sense that atoms already present in the receptor are never added
    (coordinate match), so re-running against an already-augmented receptor is a
    no-op for the previously injected atoms.
    """
    src = Path(source_pdb)
    rec = Path(receptor_pdbqt)
    out = Path(output_pdbqt) if output_pdbqt else rec

    atoms = extract_dropped_atom_residues(src, rec, min_res_atoms=min_res_atoms, tol=tol)

    lines = rec.read_text().splitlines(keepends=True)
    while lines and lines[-1].strip() in ("END", "ENDMDL"):
        lines.pop()

    if not atoms:
        out.write_text("".join(lines) + "END\n")
        return {"source": str(src), "receptor": str(rec), "output": str(out),
                "atoms_injected": 0, "residues_injected": [], "by_resname": {}}

    # Write in the format matching the target file: PDBQT (AD4 type) for the file
    # Vina docks against, plain PDB (element in cols 77-78) for its .pdb twin.
    # Both must carry the SAME heavy atoms or the PB receptor-parity check fails.
    fmt = format_pdb_atom if rec.suffix.lower() == ".pdb" else format_pdbqt_atom
    serial = _next_serial(lines)
    by_resname: dict[str, int] = {}
    new_lines: list[str] = []
    for a in atoms:
        new_lines.append(fmt(a, serial))
        serial += 1
        rn = a["resname"].strip().upper()
        by_resname[rn] = by_resname.get(rn, 0) + 1

    out.write_text("".join(lines) + "".join(new_lines) + "END\n")
    return {
        "source": str(src),
        "receptor": str(rec),
        "output": str(out),
        "atoms_injected": len(atoms),
        "residues_injected": sorted(by_resname),
        "by_resname": by_resname,
    }


def _cli() -> None:
    p = argparse.ArgumentParser(description="Inject KEEP-list HETATMs into a receptor PDBQT.")
    p.add_argument("--source-pdb", required=True, type=Path,
                   help="Original PDB containing the HETATM records to inject")
    p.add_argument("--receptor-pdbqt", required=True, type=Path,
                   help="Existing Meeko/MGLTools-prepared receptor PDBQT")
    p.add_argument("--output-pdbqt", type=Path, default=None,
                   help="Destination (defaults to overwriting --receptor-pdbqt)")
    p.add_argument("--keep", nargs="*", default=None,
                   help="Override the default KEEP residue list (3-letter codes)")
    args = p.parse_args()

    result = inject_hetatms(
        args.source_pdb, args.receptor_pdbqt, args.output_pdbqt,
        keep_residues=args.keep,
    )
    print(f"Injected {result['atoms_injected']} atoms "
          f"({', '.join(f'{r}:{n}' for r, n in result['by_resname'].items())}) "
          f"into {result['output']}")


if __name__ == "__main__":
    _cli()
