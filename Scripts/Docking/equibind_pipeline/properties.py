"""Ligand & protein property extraction (precomputed once per pipeline run)."""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors


_LIG_KEYS = (
    "lig_molecular_weight", "lig_heavy_atoms", "lig_total_atoms",
    "lig_rotatable_bonds", "lig_num_rings", "lig_aromatic_rings",
    "lig_hbd", "lig_hba", "lig_tpsa", "lig_logp", "lig_formula",
)


def get_gpu_model() -> str:
    """Return the GPU model name via nvidia-smi, or 'N/A' if unavailable."""
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=gpu_name", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=5,
        )
        names = [ln.strip() for ln in result.stdout.strip().splitlines() if ln.strip()]
        return names[0] if names else "N/A"
    except Exception:
        return "N/A"


def get_ligand_properties(sdf_path: Path) -> Dict:
    """Extract molecular properties from a ligand file using RDKit."""
    props: Dict = {k: None for k in _LIG_KEYS}
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
        mol = next(iter(suppl), None)
        if mol is None:
            mol = Chem.MolFromPDBFile(str(sdf_path), removeHs=False, sanitize=False)
        if mol is None:
            return props
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

        ri = mol.GetRingInfo()
        props.update(
            lig_molecular_weight=round(Descriptors.MolWt(mol), 2),
            lig_heavy_atoms=mol.GetNumHeavyAtoms(),
            lig_total_atoms=mol.GetNumAtoms(),
            lig_rotatable_bonds=rdMolDescriptors.CalcNumRotatableBonds(mol),
            lig_num_rings=ri.NumRings(),
            lig_aromatic_rings=rdMolDescriptors.CalcNumAromaticRings(mol),
            lig_hbd=Lipinski.NumHDonors(mol),
            lig_hba=Lipinski.NumHAcceptors(mol),
            lig_tpsa=round(Descriptors.TPSA(mol), 2),
            lig_logp=round(Descriptors.MolLogP(mol), 2),
            lig_formula=rdMolDescriptors.CalcMolFormula(mol),
        )
    except Exception as e:
        print(f"  \u26a0 Could not compute ligand properties for {sdf_path.name}: {e}")
    return props


def get_protein_properties(pdb_path: Path) -> Dict:
    """Extract basic protein properties by parsing PDB ATOM records."""
    props: Dict = {"prot_num_residues": None, "prot_num_atoms": None,
                   "prot_num_chains": None}
    try:
        residues, chains, atom_count = set(), set(), 0
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_count += 1
                    chain_id = line[21]
                    res_seq = line[22:27].strip()
                    residues.add((chain_id, res_seq))
                    chains.add(chain_id)
        props.update(
            prot_num_residues=len(residues),
            prot_num_atoms=atom_count,
            prot_num_chains=len(chains),
        )
    except Exception as e:
        print(f"  \u26a0 Could not parse protein {pdb_path.name}: {e}")
    return props


def precompute_properties(
    proteins: List[Path], ligands: List[Path],
) -> Tuple[Dict[Path, Dict], Dict[Path, Dict]]:
    """Precompute properties for all proteins and ligands."""
    print("Pre-computing ligand properties...")
    lig_cache = {lig: get_ligand_properties(lig) for lig in ligands}
    print(f"  {len(lig_cache)} ligands analysed")

    print("Pre-computing protein properties...")
    prot_cache = {prot: get_protein_properties(prot) for prot in proteins}
    print(f"  {len(prot_cache)} proteins analysed")

    return lig_cache, prot_cache
