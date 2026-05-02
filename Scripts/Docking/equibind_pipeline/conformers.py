"""RDKit conformer generation utilities."""
from __future__ import annotations

import shutil
from concurrent.futures import ProcessPoolExecutor, TimeoutError as _FT, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .config import CFG
from .monitor import monitor


# Atomic numbers (H,B,C,N,O,F,P,S,Cl,Se,Br,I) — molecules with metals or other
# elements often hang in RDKit embedding/MMFF.
_ORGANIC_ATOMIC_NUMS = frozenset({1, 5, 6, 7, 8, 9, 15, 16, 17, 34, 35, 53})


def _generate_conformers_for_seed(
    mol_block: str, seed: int, per_seed: int, output_dir: Path,
) -> List[Tuple[Path, int]]:
    """Worker (multiprocess-safe): build ``per_seed`` conformers for one seed."""
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
    if mol is None:
        return []
    if not all(a.GetAtomicNum() in _ORGANIC_ATOMIC_NUMS for a in mol.GetAtoms()):
        return []

    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1  # avoid contention with outer pool
    params.useRandomCoords = True
    params.maxIterations = 500

    try:
        cids = AllChem.EmbedMultipleConfs(mol, numConfs=per_seed, params=params)
    except Exception:
        cids = []

    can_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    out: List[Tuple[Path, int]] = []
    for ci, cid in enumerate(cids):
        if can_mmff:
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
            except Exception:
                pass
        cpath = output_dir / f"seed{seed}_c{ci}.sdf"
        try:
            w = Chem.SDWriter(str(cpath))
            w.write(mol, confId=cid)
            w.close()
            if cpath.exists() and cpath.stat().st_size > 0:
                out.append((cpath, seed))
        except Exception:
            pass
    return out


def _load_ligand(ligand_file: Path) -> Optional[Chem.Mol]:
    suffix = ligand_file.suffix.lower()
    try:
        if suffix == ".sdf":
            return next(iter(Chem.SDMolSupplier(str(ligand_file), removeHs=False)), None)
        if suffix == ".mol2":
            return Chem.MolFromMol2File(str(ligand_file), removeHs=False)
        if suffix == ".pdb":
            return Chem.MolFromPDBFile(str(ligand_file), removeHs=False)
    except Exception:
        pass
    return None


def generate_conformer_sdfs(
    ligand_file: Path, output_dir: Path,
    seeds: Optional[List[int]] = None,
    per_seed: Optional[int] = None,
) -> List[Tuple[Path, int]]:
    """Generate SDF files for all (seed, conformer) pairs."""
    seeds = list(CFG.rdkit_seeds) if seeds is None else seeds
    per_seed = CFG.conformers_per_seed if per_seed is None else per_seed
    output_dir.mkdir(parents=True, exist_ok=True)

    mol = _load_ligand(ligand_file)
    if mol is None:
        return []
    mol = Chem.AddHs(mol)
    mol_block = Chem.MolToMolBlock(mol)

    files: List[Tuple[Path, int]] = []
    if CFG.parallel_conformer_gen and len(seeds) > 1:
        n_workers = min(CFG.n_parallel_workers, len(seeds))
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            futures = {
                ex.submit(_generate_conformers_for_seed, mol_block, s, per_seed, output_dir): s
                for s in seeds
            }
            for fut in as_completed(futures):
                try:
                    files.extend(fut.result(timeout=60))
                except _FT:
                    monitor.warning(f"[TIMEOUT] Conformer gen timed out for seed {futures[fut]}")
                except Exception:
                    pass
    else:
        for seed in seeds:
            files.extend(_generate_conformers_for_seed(mol_block, seed, per_seed, output_dir))
    return files


def ensure_ligand_sdf(ligand_file: Path, prep_dir: Path) -> Path:
    """Copy/convert input ligand into a normalized SDF."""
    lig_sdf = prep_dir / f"{ligand_file.stem}.sdf"
    if lig_sdf.exists():
        return lig_sdf
    suffix = ligand_file.suffix.lower()
    if suffix == ".sdf":
        shutil.copy2(ligand_file, lig_sdf)
        return lig_sdf
    if suffix in (".mol2", ".pdb"):
        loader = Chem.MolFromMol2File if suffix == ".mol2" else Chem.MolFromPDBFile
        mol = loader(str(ligand_file), removeHs=False)
        if mol is not None:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            w = Chem.SDWriter(str(lig_sdf)); w.write(mol); w.close()
            return lig_sdf
    shutil.copy2(ligand_file, lig_sdf)
    return lig_sdf
