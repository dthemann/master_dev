import os, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
W="/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/wf"
base="/home/manndo/master_dev/Data/PoseBuster Benchmark Set"
ppm="/home/manndo/master_dev/posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv"
prot = pd.read_csv(ppm, usecols=["protein"]).protein.unique().tolist()
dropped=["7B2C_TP7","7D6O_MTE","7FRX_O88","7M31_TDR","8F4J_PHO"]
ids308 = sorted(set(prot)|set(dropped)); print("n analysed:",len(prot)," n308:",len(ids308))
df = pd.read_csv(f"{W}/copies_inventory.csv")
df = df[df.id.isin(ids308)]
multi = df[df.n_records>1]; alt = multi[~multi.is_ref]
print("multi-copy ids in 308:", multi.id.nunique(), " in 303:", multi[multi.id.isin(prot)].id.nunique())
print("copy count distribution (308):", multi.groupby("id").n_records.first().value_counts().sort_index().to_dict())
print("alt copies (308):", len(alt), " same heavy:", int(alt.same_heavy.sum()), " same formula:", int(alt.same_formula.sum()), " same elem order:", int(alt.same_elem_seq.sum()), " same charge:", int((alt.charge_copy==alt.charge_ref).sum()), " substruct match both ways:", int((alt.match_ref_in_copy & alt.match_copy_in_ref).sum()))
print("alt smiles differs (308):", alt[alt.smiles_same==False][["id","name","cd_to_ref"]].to_string())
# stereo-free smiles for those
def smi(m, iso):
    mm=Chem.RemoveHs(Chem.Mol(m),sanitize=False); Chem.SanitizeMol(mm); return Chem.MolToSmiles(mm, isomericSmiles=iso)
for i in alt[alt.smiles_same==False].id.unique():
    ref=next(iter(Chem.SDMolSupplier(f"{base}/{i}/{i}_ligand.sdf",removeHs=False,sanitize=False)))
    cps=[m for m in Chem.SDMolSupplier(f"{base}/{i}/{i}_ligands.sdf",removeHs=False,sanitize=False)]
    print(i, "nonstereo-identical:", all(smi(c,False)==smi(ref,False) for c in cps), " stereo-identical:", all(smi(c,True)==smi(ref,True) for c in cps))
print("\nalt copies by centroid distance to ref (308): <=2:", int((alt.cd_to_ref<=2).sum())," <=4:", int((alt.cd_to_ref<=4).sum())," <=8:", int((alt.cd_to_ref<=8).sum())," <=12:", int((alt.cd_to_ref<=12).sum()), " >12:", int((alt.cd_to_ref>12).sum()), " total:", len(alt))
near = alt.groupby("id").cd_to_ref.min()
print("ids nearest alt <=2/4/8/12:", int((near<=2).sum()), int((near<=4).sum()), int((near<=8).sum()), int((near<=12).sum()), " median", round(near.median(),2), " min", round(near.min(),2), " max", round(near.max(),2))
print(alt[alt.cd_to_ref<=12][["id","name","heavy_ref","cd_to_ref"]].sort_values("cd_to_ref").to_string())
# ---- four low-Jaccard ids: alternate site characterisation
def prot_atoms(pdb):
    rows=[]
    for ln in open(pdb):
        if ln.startswith(("ATOM","HETATM")):
            el = ln[76:78].strip() or ln[12:14].strip()
            if el.upper()=="H": continue
            rows.append((ln[21], ln[17:20].strip(), int(ln[22:26]), ln[0:6].strip(), float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
    return rows
def env(m, atoms, cut=4.5):
    xyz=np.array([[a[4],a[5],a[6]] for a in atoms]); pos=m.GetConformer().GetPositions()[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()!=1]]
    d=np.linalg.norm(xyz[:,None,:]-pos[None,:,:],axis=2); close=(d.min(axis=1)<=cut)
    res=set((atoms[i][0],atoms[i][1],atoms[i][2],atoms[i][3]) for i in np.where(close)[0])
    return res, int((d.min(axis=1)<=4.0).sum())
for i in ["7A9E_R4W","7TUO_KL9","7VKZ_NOJ","7Z1Q_NIO"]:
    ref=next(iter(Chem.SDMolSupplier(f"{base}/{i}/{i}_ligand.sdf",removeHs=False,sanitize=False)))
    cps=[m for m in Chem.SDMolSupplier(f"{base}/{i}/{i}_ligands.sdf",removeHs=False,sanitize=False)]
    atoms=prot_atoms(f"{base}/{i}/{i}_protein.pdb")
    chains=sorted(set(a[0] for a in atoms)); hets=sorted(set((a[1]) for a in atoms if a[3]=="HETATM"))
    print(f"\n== {i}: formula {df[(df.id==i)&(df.idx==0)].formula_ref.iloc[0]}, heavy {ref.GetNumHeavyAtoms()}, protein chains {chains}, HETATM residues {hets}")
    rref, nref = env(ref, atoms); 
    def fmt(s): 
        return sorted(f"{c}:{r}{n}" for c,r,n,k in s)
    print("  REF", cps[0].GetProp("_Name"), "| protein heavy atoms within 4A:", nref, "| residues(4.5A):", fmt(rref))
    refpos=ref.GetConformer().GetPositions()
    for k,c in enumerate(cps[1:],1):
        rc, nc = env(c, atoms); jac=len(rc&rref)/len(rc|rref) if (rc|rref) else float('nan')
        cpos=c.GetConformer().GetPositions(); dmin=np.linalg.norm(refpos[:,None,:]-cpos[None,:,:],axis=2).min()
        cd=np.linalg.norm(refpos.mean(0)-cpos.mean(0))
        print(f"  ALT{k} {c.GetProp('_Name')} | centroid sep {cd:.2f} A | min heavy-atom sep to ref {dmin:.2f} A | prot atoms within 4A: {nc} | Jaccard {jac:.2f} | residues(4.5A): {fmt(rc)}")
        print("      shared:", fmt(rc&rref))
