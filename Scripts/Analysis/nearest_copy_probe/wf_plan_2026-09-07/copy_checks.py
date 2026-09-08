import os, sys, math, json
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors
RDLogger.DisableLog("rdApp.*")
base = "/home/manndo/master_dev/Data/PoseBuster Benchmark Set"
out = "/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/wf"

def heavy(m):
    return [a for a in m.GetAtoms() if a.GetAtomicNum() != 1]
def elem_seq(m):
    return tuple(a.GetSymbol() for a in heavy(m))
def formula(m):
    mm = Chem.Mol(m); mm.UpdatePropertyCache(strict=False)
    return rdMolDescriptors.CalcMolFormula(Chem.RemoveHs(mm, sanitize=False))
def centroid(m):
    pos = m.GetConformer().GetPositions()
    idx = [a.GetIdx() for a in heavy(m)]
    return pos[idx].mean(axis=0)
def can_smiles(m):
    try:
        mm = Chem.RemoveHs(Chem.Mol(m), sanitize=False); Chem.SanitizeMol(mm)
        return Chem.MolToSmiles(mm)
    except Exception:
        return None
def charges(m):
    return sum(a.GetFormalCharge() for a in m.GetAtoms())

rows = []
for d in sorted(os.listdir(base)):
    p = os.path.join(base, d)
    if not os.path.isdir(p): continue
    ref_f = os.path.join(p, f"{d}_ligand.sdf"); all_f = os.path.join(p, f"{d}_ligands.sdf")
    if not (os.path.exists(ref_f) and os.path.exists(all_f)): 
        rows.append(dict(id=d, note="missing file")); continue
    ref = next(iter(Chem.SDMolSupplier(ref_f, removeHs=False, sanitize=False)))
    copies = [m for m in Chem.SDMolSupplier(all_f, removeHs=False, sanitize=False)]
    n_none = sum(m is None for m in copies)
    copies = [m for m in copies if m is not None]
    rc = centroid(ref); rseq = elem_seq(ref); rform = formula(ref); rsmi = can_smiles(ref); rchg = charges(ref)
    ref_h = Chem.RemoveHs(ref, sanitize=False)
    for k, m in enumerate(copies):
        cc = centroid(m)
        dist = float(np.linalg.norm(cc - rc))
        m_h = Chem.RemoveHs(m, sanitize=False)
        # substructure both directions on heavy atoms (bond orders as read)
        try: 
            m_h.UpdatePropertyCache(strict=False); ref_h.UpdatePropertyCache(strict=False)
            match_ref_in_copy = bool(m_h.GetSubstructMatches(ref_h, uniquify=False, useChirality=False))
            match_copy_in_ref = bool(ref_h.GetSubstructMatches(m_h, uniquify=False, useChirality=False))
        except Exception as e:
            match_ref_in_copy = match_copy_in_ref = None
        rows.append(dict(id=d, idx=k, name=m.GetProp("_Name") if m.HasProp("_Name") else "",
                         n_records=len(copies), n_none=n_none,
                         heavy_ref=len(heavy(ref)), heavy_copy=len(heavy(m)),
                         same_heavy=len(heavy(ref))==len(heavy(m)),
                         same_elem_seq=(elem_seq(m)==rseq), formula_ref=rform, formula_copy=formula(m),
                         same_formula=(formula(m)==rform), smiles_same=(can_smiles(m)==rsmi) if rsmi else None,
                         charge_ref=rchg, charge_copy=charges(m), n_atoms_all_ref=ref.GetNumAtoms(), n_atoms_all_copy=m.GetNumAtoms(),
                         cd_to_ref=dist, is_ref=dist<1e-3,
                         match_ref_in_copy=match_ref_in_copy, match_copy_in_ref=match_copy_in_ref))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(out, "copies_inventory.csv"), index=False)
multi = df[df.n_records > 1]
ids_multi = sorted(multi.id.unique())
print("ids with folder:", df.id.nunique(), " multi-copy ids:", len(ids_multi))
print("records with None (unparseable):", int(df.n_none.max()))
# reference index
refidx = df[df.is_ref].groupby("id").idx.min()
print("reference record index distribution:", refidx.value_counts().to_dict())
print("ids where no record matches the reference centroid:", sorted(set(df.id) - set(refidx.index)))
alt = multi[~multi.is_ref]
print("alternate copies total:", len(alt))
print("alt same heavy count:", int(alt.same_heavy.sum()), "/", len(alt))
print("alt same formula:", int(alt.same_formula.sum()), "/", len(alt))
print("alt same element sequence (atom order):", int(alt.same_elem_seq.sum()), "/", len(alt))
print("alt same canonical smiles:", alt.smiles_same.value_counts(dropna=False).to_dict())
print("alt same total charge:", int((alt.charge_copy==alt.charge_ref).sum()), "/", len(alt))
print("alt ref-substructure-in-copy match:", alt.match_ref_in_copy.value_counts(dropna=False).to_dict())
print("alt copy-substructure-in-ref match:", alt.match_copy_in_ref.value_counts(dropna=False).to_dict())
bad = alt[(~alt.same_heavy) | (~alt.same_formula) | (~alt.same_elem_seq) | (alt.smiles_same==False) | (alt.match_ref_in_copy==False)]
print("\nALT COPIES DIFFERING FROM REFERENCE:")
print(bad[["id","idx","name","heavy_ref","heavy_copy","formula_ref","formula_copy","same_elem_seq","smiles_same","match_ref_in_copy","match_copy_in_ref","cd_to_ref"]].to_string())
# centroid distance bands (per alternate copy and per id nearest alternate)
print("\nalt copies by centroid distance to reference: <=2:", int((alt.cd_to_ref<=2).sum()), " <=4:", int((alt.cd_to_ref<=4).sum()), " <=8:", int((alt.cd_to_ref<=8).sum()), " >8:", int((alt.cd_to_ref>8).sum()))
nearest = alt.groupby("id").cd_to_ref.min()
print("ids whose NEAREST alternate is <=2:", int((nearest<=2).sum()), " <=4:", int((nearest<=4).sum()), " <=8:", int((nearest<=8).sum()), " min:", round(nearest.min(),2), " median:", round(nearest.median(),2))
print("ids with any alternate <=8 A:", sorted(nearest[nearest<=8].index.tolist()))
print(alt[alt.cd_to_ref<=8][["id","idx","name","cd_to_ref"]].to_string())
# names -> chain parse
def chain_of(n):
    parts = n.split("_")
    return parts[2] if len(parts)>=4 else ""
multi = multi.copy(); multi["chain"] = multi.name.map(chain_of)
same_chain = multi.groupby("id").apply(lambda g: g.chain.nunique() < len(g))
print("\nids where two copies share a chain id (same-chain second site or altloc):", int(same_chain.sum()), sorted(same_chain[same_chain].index.tolist()))
# distribution of copy counts
print("copy count distribution:", multi.groupby("id").n_records.first().value_counts().sort_index().to_dict())
# dropped complexes
dropped = ["7B2C_TP7","7D6O_MTE","7FRX_O88","7M31_TDR","8F4J_PHO"]
print("dropped ids copy counts:", df[df.id.isin(dropped)].groupby("id").n_records.first().to_dict())
