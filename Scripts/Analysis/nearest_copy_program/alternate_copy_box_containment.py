#!/usr/bin/env python
"""How many alternate deposited ligand copies would a PoseBusters-style 25 A cube on the reference instance contain?

For every benchmark entry with more than one record in <ID>_ligands.sdf, the reference instance is record 0
(== <ID>_ligand.sdf). A cube of edge 25 A is centred on the reference's heavy-atom centroid; an alternate copy
counts as inside when its own heavy-atom centroid lies inside the cube (the rule the appendix sentence states).
Two further rules are reported for transparency: any heavy atom of the copy inside the cube, and centroid within a
12.5 A sphere. Read-only on the benchmark; writes one JSON to --out.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

def heavy_xyz(mol):
    conf = mol.GetConformer()
    return np.array([list(conf.GetAtomPosition(a.GetIdx())) for a in mol.GetAtoms() if a.GetAtomicNum() > 1])

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--benchmark-dir", type=Path, default=Path("Data/PoseBuster Benchmark Set"))
    ap.add_argument("--edge", type=float, default=25.0)
    ap.add_argument("--ids-file", type=Path, default=None, help="restrict to these entry ids (one per line); default all folders")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    half = a.edge / 2
    n_ids = n_multi = n_alt = in_centroid = in_anyatom = in_sphere = 0
    per_id = {}
    keep = set(x.strip() for x in a.ids_file.read_text().split() if x.strip()) if a.ids_file else None
    for d in sorted(p for p in a.benchmark_dir.iterdir() if p.is_dir()):
        if keep is not None and d.name not in keep:
            continue
        f = d / f"{d.name}_ligands.sdf"
        if not f.exists():
            continue
        mols = [m for m in Chem.SDMolSupplier(str(f), removeHs=False, sanitize=False) if m is not None]
        n_ids += 1
        if len(mols) < 2:
            continue
        n_multi += 1
        ref = heavy_xyz(mols[0]); c0 = ref.mean(axis=0)
        rec = []
        for m in mols[1:]:
            xyz = heavy_xyz(m); c = xyz.mean(axis=0)
            cen_in = bool(np.all(np.abs(c - c0) <= half))
            any_in = bool(np.any(np.all(np.abs(xyz - c0) <= half, axis=1)))
            sph_in = bool(np.linalg.norm(c - c0) <= half)
            n_alt += 1; in_centroid += cen_in; in_anyatom += any_in; in_sphere += sph_in
            rec.append({"centroid_dist_A": round(float(np.linalg.norm(c - c0)), 2), "centroid_in_cube": cen_in, "any_atom_in_cube": any_in, "centroid_in_sphere": sph_in})
        per_id[d.name] = rec
    out = {"ids_file": str(a.ids_file) if a.ids_file else None, "edge_A": a.edge, "n_entries": n_ids, "n_multi_copy_entries": n_multi, "n_alternate_copies": n_alt,
           "alternate_centroid_inside_cube": in_centroid, "alternate_any_heavy_atom_inside_cube": in_anyatom,
           "alternate_centroid_inside_sphere_r_half_edge": in_sphere, "per_entry": per_id}
    a.out.parent.mkdir(parents=True, exist_ok=True); a.out.write_text(json.dumps(out, indent=1))
    print({k: v for k, v in out.items() if k != "per_entry"})

if __name__ == "__main__":
    main()
