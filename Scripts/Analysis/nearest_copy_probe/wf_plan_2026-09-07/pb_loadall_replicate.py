"""Replicate exactly what installed PoseBusters (redock.yml loading params, load_all=True on mol_true) returns
for the rank-1 poses of the three headline arms on the 138 multi-copy complexes, and compare with per_pose_metrics
(single reference) and with the probe's nearest-copy symmetry_rmsd."""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0,'/home/manndo/master_dev/Scripts/Analysis')
from rdkit import Chem, RDLogger; RDLogger.DisableLog('rdApp.*')
import posebusters_pose_comparison as ppc
from posebusters.tools.loading import safe_load_mol
from posebusters.modules.rmsd import check_rmsd
BENCH=Path('/home/manndo/master_dev/Data/PoseBuster Benchmark Set')
alld=pd.read_csv('/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv', low_memory=False)
arms=['autodock_mgltools_exh128_gnina','diffdock_smina','equibind_unguided_gnina']
d=alld[(alld.method.isin(arms))&(alld.eff_rank==1)&(alld.n_copies>1)].copy()
print(len(d),'rank-1 poses on multi-copy complexes')
LP=dict(cleanup=False,sanitize=False,add_hs=False,assign_stereo=False,load_all=True)
rows=[]
for prot,g in d.groupby('protein'):
    cdir=BENCH/prot
    crystal_h=Chem.RemoveHs(ppc.load_first_mol(cdir/f'{prot}_ligand.sdf'))
    mt_all=safe_load_mol(cdir/f'{prot}_ligands.sdf', **LP)         # what PoseBusters redock does
    mt_one=safe_load_mol(cdir/f'{prot}_ligand.sdf', **{**LP,'load_all':False})
    for _,r in g.iterrows():
        pose=ppc.load_first_mol(Path(r.pose_file)); 
        if pose is None: continue
        pose=ppc.reassign_template(pose,crystal_h) or pose
        ph=Chem.RemoveHs(pose)
        try: ra=check_rmsd(ph, Chem.RemoveHs(mt_all,sanitize=False))['results']
        except Exception as e: ra={'rmsd':np.nan,'kabsch_rmsd':np.nan,'centroid_distance':np.nan}
        try: ro=check_rmsd(ph, Chem.RemoveHs(mt_one,sanitize=False))['results']
        except Exception as e: ro={'rmsd':np.nan,'kabsch_rmsd':np.nan,'centroid_distance':np.nan}
        rows.append(dict(method=r.method,protein=prot,n_conf=mt_all.GetNumConformers(),
            csv_pb_rmsd=r.pb_rmsd,csv_bestfit=r.bestfit_rmsd,csv_cd=r.centroid_dist,csv_rmsd=r.rmsd,pb_valid=r.pb_valid,
            probe_nearest=r.rmsd_nearest_copy,
            pb1_rmsd=ro['rmsd'],pb1_kabsch=ro['kabsch_rmsd'],pb1_cd=ro['centroid_distance'],
            pbN_rmsd=ra['rmsd'],pbN_kabsch=ra['kabsch_rmsd'],pbN_cd=ra['centroid_distance']))
o=pd.DataFrame(rows); o.to_csv('pb_loadall_rank1.csv',index=False)
print('single-ref check_rmsd reproduces csv pb_rmsd: max|diff| =', np.nanmax(np.abs(o.pb1_rmsd-o.csv_pb_rmsd)), '; kabsch vs bestfit max|diff| =', np.nanmax(np.abs(o.pb1_kabsch-o.csv_bestfit)))
print('load_all rmsd vs probe symmetry nearest: max|diff| =', np.nanmax(np.abs(o.pbN_rmsd-o.probe_nearest)), '; n differing >0.05 A:', int((np.abs(o.pbN_rmsd-o.probe_nearest)>0.05).sum()))
print(o[np.abs(o.pbN_rmsd-o.probe_nearest)>0.05][['method','protein','probe_nearest','pbN_rmsd','pb1_rmsd','csv_rmsd']].to_string())
for m,g in o.groupby('method'):
    v=g[g.pb_valid==True]
    print(f'\n{m}: n={len(g)} valid={len(v)}')
    print('  near-native (<=2) single:', int((v.pb1_rmsd<=2).sum()), ' load_all:', int((v.pbN_rmsd<=2).sum()))
    print('  form gate kabsch<=1 on valid&near(single):', int(((v.pb1_rmsd<=2)&(v.pb1_kabsch<=1)).sum()), ' | on valid&near(load_all) with kabsch from argmin copy:', int(((v.pbN_rmsd<=2)&(v.pbN_kabsch<=1)).sum()), ' | valid&near(load_all) but kabsch vs single ref <=1:', int(((v.pbN_rmsd<=2)&(v.pb1_kabsch<=1)).sum()))
    flip=v[(v.pbN_rmsd<=2)&(v.pb1_rmsd>2)]
    print('  flips:', len(flip), ' kabsch change on flips (argmin copy - ref): median', round(float((flip.pbN_kabsch-flip.pb1_kabsch).median()),3), ' max|.|', round(float((flip.pbN_kabsch-flip.pb1_kabsch).abs().max()),3))
    print('  form verdict changes among flips (kabsch<=1):', int(((flip.pbN_kabsch<=1)!=(flip.pb1_kabsch<=1)).sum()))
    print('  4A centroid reach single:', int((g.pb1_cd<=4).sum()), ' load_all (argmin-rmsd copy centroid):', int((g.pbN_cd<=4).sum()))
    # kabsch differs between copies even for non-flips (conformer differences between copies)
    nn=g[(g.pbN_rmsd<=2)]
    print('  among all near(load_all): kabsch(argmin copy) vs kabsch(ref) max|diff|', round(float((nn.pbN_kabsch-nn.pb1_kabsch).abs().max()),3))
