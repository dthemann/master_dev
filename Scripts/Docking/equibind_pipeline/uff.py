"""UFF post-docking minimization (ligand relaxes inside fixed protein pocket)."""
from __future__ import annotations

import threading
from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdForceFieldHelpers
from scipy.spatial import cKDTree

from .config import CFG
from .monitor import monitor


# Maximum expected valence per element at formal charge 0.
_MAX_VALENCE = {
    "C": 4, "N": 4, "O": 3, "S": 6, "P": 5,
    "H": 1, "F": 1, "Cl": 1, "Br": 1, "I": 1,
    "B": 3, "Si": 4, "Se": 2,
}

_NORMAL_VALENCE = {"N": 3, "O": 2}

_NON_STANDARD_ELEMENTS = {"ZN", "FE", "MG", "CA", "MN", "CO", "CU", "NI", "NA",
                          "K", "CD", "HG", "PT", "AU", "AG"}

_STANDARD_RESIDUES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "HOH", "WAT",
}


# ── Shared mol cleanup helpers ─────────────────────────────────────────────

def _fix_spurious_bonds(mol: Chem.Mol) -> Tuple[Chem.Mol, int]:
    """Iteratively remove the longest bond from over-valent atoms."""
    emol = Chem.RWMol(mol)
    removed = 0
    for _ in range(20):
        found_bad = False
        conf = emol.GetConformer()
        for atom in emol.GetAtoms():
            max_v = _MAX_VALENCE.get(atom.GetSymbol())
            if max_v is None:
                continue
            explicit_v = int(round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())))
            if explicit_v <= max_v:
                continue
            worst_bond = None
            worst_len = -1.0
            for bond in atom.GetBonds():
                p1 = conf.GetAtomPosition(bond.GetBeginAtomIdx())
                p2 = conf.GetAtomPosition(bond.GetEndAtomIdx())
                d = p1.Distance(p2)
                if d > worst_len:
                    worst_len = d
                    worst_bond = bond
            if worst_bond is not None:
                emol.RemoveBond(worst_bond.GetBeginAtomIdx(), worst_bond.GetEndAtomIdx())
                removed += 1
                found_bad = True
        if not found_bad:
            break
    return emol.GetMol(), removed


def _fix_ligand_valence(mol: Chem.Mol) -> int:
    """Set missing formal charges (e.g. quaternary N+, trivalent O+)."""
    n_fixed = 0
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() != 0:
            continue
        normal_v = _NORMAL_VALENCE.get(atom.GetSymbol())
        if normal_v is None:
            continue
        explicit_v = int(round(sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())))
        if explicit_v > normal_v:
            atom.SetFormalCharge(explicit_v - normal_v)
            n_fixed += 1
    return n_fixed


def _sanitize_protein(mol: Chem.Mol) -> None:
    """Sanitize a protein mol, falling back to partial sanitize on failure."""
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            Chem.SanitizeMol(
                mol,
                Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
            )
        except Exception:
            pass
    Chem.FastFindRings(mol)


# ── Cleaned-protein cache (load + strip + fix bonds + sanitize ONCE) ─────

_cleaned_protein_cache: Dict[str, Optional[Chem.Mol]] = {}
_cleaned_protein_lock = threading.Lock()


def _strip_nonstandard_residues(mol: Chem.Mol) -> Chem.Mol:
    remove_idx = set()
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetSymbol().upper() in _NON_STANDARD_ELEMENTS:
            remove_idx.add(i)
            continue
        info = atom.GetPDBResidueInfo()
        if info and info.GetResidueName().strip().upper() not in _STANDARD_RESIDUES:
            remove_idx.add(i)
    if not remove_idx:
        return mol
    emol = Chem.RWMol(mol)
    for idx in sorted(remove_idx, reverse=True):
        emol.RemoveAtom(int(idx))
    return emol.GetMol()


def get_cleaned_protein(
    protein_pdb: Path,
    strip_nonstandard: Optional[bool] = None,
) -> Optional[Chem.Mol]:
    """Cached cleaned-protein producer."""
    strip = CFG.uff_strip_nonstandard if strip_nonstandard is None else strip_nonstandard
    key = str(protein_pdb)
    with _cleaned_protein_lock:
        if key in _cleaned_protein_cache:
            return _cleaned_protein_cache[key]

    try:
        mol = Chem.MolFromPDBFile(str(protein_pdb), removeHs=True, sanitize=False)
        if mol is None:
            with _cleaned_protein_lock:
                _cleaned_protein_cache[key] = None
            return None

        if strip:
            mol = _strip_nonstandard_residues(mol)

        mol, n_removed = _fix_spurious_bonds(mol)
        if n_removed > 0:
            monitor.info(f"UFF: removed {n_removed} spurious PDB proximity bonds (cached)")

        _sanitize_protein(mol)

        if mol.GetNumConformers() == 0:
            with _cleaned_protein_lock:
                _cleaned_protein_cache[key] = None
            return None

        with _cleaned_protein_lock:
            _cleaned_protein_cache[key] = mol
        return mol
    except Exception as e:
        monitor.warning(f"Failed to load/clean protein {protein_pdb}: {e}")
        with _cleaned_protein_lock:
            _cleaned_protein_cache[key] = None
        return None


def prewarm_protein_cache(protein_paths: Iterable[Path],
                          strip_nonstandard: Optional[bool] = None) -> None:
    """Pre-load proteins so forked workers inherit the cache via copy-on-write."""
    unique = {str(p) for p in protein_paths}
    for p in unique:
        get_cleaned_protein(Path(p), strip_nonstandard)
    monitor.info(f"UFF: pre-warmed protein cache for {len(unique)} unique protein(s)")


# ── Nearby-protein extraction ────────────────────────────────────────────

def _extract_nearby_protein_atoms(
    protein_pdb: Path, ligand_coords: np.ndarray,
    radius: float, strip_nonstandard: bool,
) -> Optional[Chem.Mol]:
    try:
        protein_mol = get_cleaned_protein(protein_pdb, strip_nonstandard)
        if protein_mol is None:
            return None

        prot_coords = protein_mol.GetConformer().GetPositions()
        prot_tree = cKDTree(prot_coords)
        near_atom_indices: set = set()
        for lc in ligand_coords:
            near_atom_indices.update(prot_tree.query_ball_point(lc, radius))
        if not near_atom_indices:
            return None

        # Expand to full residues
        residue_keys = set()
        for i in near_atom_indices:
            info = protein_mol.GetAtomWithIdx(int(i)).GetPDBResidueInfo()
            if info:
                residue_keys.add(
                    (info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode()))
        if residue_keys:
            full = set()
            for i in range(protein_mol.GetNumAtoms()):
                info = protein_mol.GetAtomWithIdx(i).GetPDBResidueInfo()
                if info:
                    key = (info.GetChainId(), info.GetResidueNumber(), info.GetInsertionCode())
                    if key in residue_keys:
                        full.add(int(i))
            near_atom_indices = full

        emol = Chem.RWMol(protein_mol)
        all_atom_idx = set(range(protein_mol.GetNumAtoms()))
        for idx in sorted(all_atom_idx - near_atom_indices, reverse=True):
            emol.RemoveAtom(int(idx))
        nearby = emol.GetMol()
        nearby, _ = _fix_spurious_bonds(nearby)
        _sanitize_protein(nearby)
        return nearby
    except Exception as e:
        monitor.warning(f"Failed to extract nearby protein atoms: {e}")
        return None


# ── Main entry point ─────────────────────────────────────────────────────

def _write_sdf(mol: Chem.Mol, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(path))
    w.write(mol)
    w.close()


def _ligand_only_minimize(lig_mol: Chem.Mol, output_sdf: Path, add_hs: bool,
                          max_iters: int, energy_tol: float, force_tol: float,
                          msg_suffix: str = "") -> Tuple[bool, float, float, str]:
    if not rdForceFieldHelpers.UFFHasAllMoleculeParams(lig_mol):
        return False, 0.0, 0.0, "UFF missing parameters for ligand atoms"
    ff = rdForceFieldHelpers.UFFGetMoleculeForceField(lig_mol)
    if ff is None:
        return False, 0.0, 0.0, "Could not create UFF force field for ligand"
    e_before = ff.CalcEnergy()
    ff.Minimize(maxIts=max_iters, energyTol=energy_tol, forceTol=force_tol)
    e_after = ff.CalcEnergy()
    if add_hs:
        lig_mol = Chem.RemoveHs(lig_mol)
    _write_sdf(lig_mol, output_sdf)
    return True, e_before, e_after, msg_suffix or ""


def uff_minimize_pose(
    docked_sdf: Path, protein_pdb: Path, output_sdf: Path,
    max_iters: Optional[int] = None,
    energy_tol: Optional[float] = None,
    force_tol: Optional[float] = None,
    proximity_radius: Optional[float] = None,
    add_hs: Optional[bool] = None,
    strip_nonstandard: Optional[bool] = None,
) -> Tuple[bool, float, float, str]:
    """Minimize a docked ligand using UFF with the protein held fixed.

    Returns ``(success, energy_before, energy_after, message)``.
    """
    max_iters = CFG.uff_max_iters if max_iters is None else max_iters
    energy_tol = CFG.uff_energy_tol if energy_tol is None else energy_tol
    force_tol = CFG.uff_force_tol if force_tol is None else force_tol
    proximity_radius = CFG.uff_proximity_radius if proximity_radius is None else proximity_radius
    add_hs = CFG.uff_add_hydrogens if add_hs is None else add_hs
    strip_nonstandard = CFG.uff_strip_nonstandard if strip_nonstandard is None else strip_nonstandard

    try:
        suppl = Chem.SDMolSupplier(str(docked_sdf), removeHs=False, sanitize=False)
        lig_mol = next(iter(suppl), None)
        if lig_mol is None:
            return False, 0.0, 0.0, "Could not load docked ligand SDF"

        _fix_ligand_valence(lig_mol)
        try:
            Chem.SanitizeMol(lig_mol)
        except Exception:
            try:
                Chem.SanitizeMol(
                    lig_mol,
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                    | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                )
            except Exception as e:
                return False, 0.0, 0.0, f"Ligand sanitization failed: {e}"
        try:
            Chem.FastFindRings(lig_mol)
        except Exception:
            pass

        lig_coords = lig_mol.GetConformer().GetPositions()
        n_lig_atoms = lig_mol.GetNumAtoms()

        prot_nearby = _extract_nearby_protein_atoms(
            protein_pdb, lig_coords, proximity_radius, strip_nonstandard)

        # ── Path A: no protein context — minimize ligand alone ───────────
        if prot_nearby is None or prot_nearby.GetNumAtoms() == 0:
            monitor.info("UFF: No nearby protein atoms found, minimizing ligand in vacuum")
            if add_hs:
                lig_mol = Chem.AddHs(lig_mol, addCoords=True)
            return _ligand_only_minimize(
                lig_mol, output_sdf, add_hs, max_iters, energy_tol, force_tol)

        # ── Add hydrogens (best-effort) ──────────────────────────────────
        if add_hs:
            try:
                lig_mol = Chem.AddHs(lig_mol, addCoords=True)
            except Exception:
                pass
            try:
                prot_nearby = Chem.AddHs(prot_nearby, addCoords=True)
            except Exception:
                pass

        n_lig_with_h = lig_mol.GetNumAtoms()
        n_prot_with_h = prot_nearby.GetNumAtoms()

        # ── Pre-check: ensure UFF can handle each part ────────────────────
        if not rdForceFieldHelpers.UFFHasAllMoleculeParams(prot_nearby):
            bad_keys = set()
            for i in range(prot_nearby.GetNumAtoms()):
                atom = prot_nearby.GetAtomWithIdx(i)
                if str(atom.GetHybridization()) == "UNSPECIFIED":
                    info = atom.GetPDBResidueInfo()
                    if info:
                        bad_keys.add((info.GetChainId(), info.GetResidueNumber(),
                                      info.GetInsertionCode()))
                    else:
                        bad_keys.add(("?", i, ""))
            if bad_keys:
                keep = set()
                for i in range(prot_nearby.GetNumAtoms()):
                    info = prot_nearby.GetAtomWithIdx(i).GetPDBResidueInfo()
                    key = (info.GetChainId(), info.GetResidueNumber(),
                           info.GetInsertionCode()) if info else ("?", i, "")
                    if key not in bad_keys:
                        keep.add(i)
                if not keep:
                    monitor.info("UFF: all nearby protein residues lack UFF params, ligand-only")
                    return _ligand_only_minimize(
                        lig_mol, output_sdf, add_hs, max_iters, energy_tol, force_tol,
                        msg_suffix="ligand-only (all protein residues unsupported)")
                emol = Chem.RWMol(prot_nearby)
                for idx in sorted(set(range(prot_nearby.GetNumAtoms())) - keep, reverse=True):
                    emol.RemoveAtom(int(idx))
                prot_nearby = emol.GetMol()
                prot_nearby, _ = _fix_spurious_bonds(prot_nearby)
                _sanitize_protein(prot_nearby)
                n_prot_with_h = prot_nearby.GetNumAtoms()
                monitor.info(f"UFF: stripped {len(bad_keys)} bad residue(s), "
                             f"{n_prot_with_h} protein atoms remain")

        # ── Combine ligand + protein ─────────────────────────────────────
        combined = Chem.CombineMols(lig_mol, prot_nearby)
        if combined.GetNumConformers() == 0:
            return False, 0.0, 0.0, "Combined molecule has no conformer"
        _fix_ligand_valence(combined)
        _sanitize_protein(combined)

        if not rdForceFieldHelpers.UFFHasAllMoleculeParams(combined):
            blame = []
            if not rdForceFieldHelpers.UFFHasAllMoleculeParams(lig_mol):
                blame.append("ligand")
            if not rdForceFieldHelpers.UFFHasAllMoleculeParams(prot_nearby):
                blame.append("protein")
            if not blame:
                blame.append("combined-only (parts OK individually)")
            monitor.warning(f"UFF params missing in: {', '.join(blame)} \u2014 ligand-only")
            return _ligand_only_minimize(
                lig_mol, output_sdf, add_hs, max_iters, energy_tol, force_tol,
                msg_suffix=f"ligand-only ({', '.join(blame)})")

        try:
            ff = rdForceFieldHelpers.UFFGetMoleculeForceField(
                combined, vdwThresh=CFG.uff_vdw_thresh)
        except Exception as e:
            monitor.warning(f"UFF FF construction failed on combined mol: {e}")
            return _ligand_only_minimize(
                lig_mol, output_sdf, add_hs, max_iters, energy_tol, force_tol,
                msg_suffix="ligand-only (combined FF failed)")
        if ff is None:
            return False, 0.0, 0.0, "Could not create UFF force field for combined system"

        # Constrain protein atoms (after the ligand atoms)
        for i in range(n_lig_with_h, n_lig_with_h + n_prot_with_h):
            ff.AddFixedPoint(i)

        e_before = ff.CalcEnergy()
        converged = ff.Minimize(maxIts=max_iters, energyTol=energy_tol, forceTol=force_tol)
        e_after = ff.CalcEnergy()
        msg = "converged" if converged == 0 else f"not converged (code {converged})"

        # Extract minimized ligand coordinates
        combined_conf = combined.GetConformer()
        out_lig = deepcopy(lig_mol)
        out_conf = out_lig.GetConformer()
        for i in range(n_lig_with_h):
            out_conf.SetAtomPosition(i, combined_conf.GetAtomPosition(i))
        if add_hs:
            out_lig = Chem.RemoveHs(out_lig)
        _write_sdf(out_lig, output_sdf)

        monitor.info(f"UFF minimization: E {e_before:.1f} \u2192 {e_after:.1f} kcal/mol "
                     f"(\u0394E={e_after - e_before:.1f}), {msg}, "
                     f"{n_lig_with_h} ligand + {n_prot_with_h} protein atoms")
        return True, e_before, e_after, msg
    except Exception as e:
        monitor.warning(f"UFF minimization failed: {e}")
        return False, 0.0, 0.0, str(e)
