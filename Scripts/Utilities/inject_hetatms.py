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
