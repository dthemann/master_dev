"""Re-score benchmark poses against EVERY crystallographic copy in <ID>_ligands.sdf.

Uses the thesis's own loaders and RMSD engine (imported from
posebusters_pose_comparison.py) so the reference-copy value must reproduce the
`rmsd` column exactly; the nearest-copy value is the sensitivity quantity.
"""
import sys, os, argparse
from pathlib import Path
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis')
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import posebusters_pose_comparison as ppc

BENCH = Path('/home/manndo/master_dev/Data/PoseBuster Benchmark Set')
PPM_DEFAULT = '/home/manndo/master_dev/posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv'

def load_all_mols(sdf):
    out = []
    for m in Chem.SDMolSupplier(str(sdf), removeHs=True, sanitize=False):
        if m is None: continue
        try: Chem.SanitizeMol(m)
        except Exception:
            try: m.UpdatePropertyCache(strict=False); Chem.GetSymmSSSR(m)
            except Exception: continue
        out.append(Chem.RemoveHs(m))
    return out

def work(args):
    prot, recs = args
    cdir = BENCH / prot
    crystal = ppc.load_first_mol(cdir / f'{prot}_ligand.sdf')
    crystal_h = Chem.RemoveHs(crystal)
    copies = load_all_mols(cdir / f'{prot}_ligands.sdf')
    ref_c = crystal_h.GetConformer().GetPositions().mean(0)
    copy_c = [c.GetConformer().GetPositions().mean(0) for c in copies]
    ref_idx = int(np.argmin([np.linalg.norm(c - ref_c) for c in copy_c]))
    out = []
    for rec in recs:
        pose = ppc.load_first_mol(Path(rec['pose_file']))
        if pose is None:
            out.append({**rec, 'rmsd_ref_re': np.nan}); continue
        pose = ppc.reassign_template(pose, crystal_h) or pose
        r = [ppc.symmetry_rmsd(pose, c) for c in copies]
        pc = Chem.RemoveHs(pose).GetConformer().GetPositions().mean(0)
        cd = [float(np.linalg.norm(pc - c)) for c in copy_c]
        alt = [i for i in range(len(copies)) if i != ref_idx]
        out.append({**rec,
            'n_copies': len(copies),
            'rmsd_ref_re': ppc.symmetry_rmsd(pose, crystal_h),
            'rmsd_refcopy': r[ref_idx],
            'rmsd_nearest_copy': float(np.nanmin(r)),
            'nearest_copy_is_ref': bool(int(np.nanargmin(r)) == ref_idx),
            'cd_ref': cd[ref_idx],
            'cd_min_alt': min([cd[i] for i in alt]) if alt else np.nan,
            'rmsd_min_alt': float(np.nanmin([r[i] for i in alt])) if alt else np.nan,
        })
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--arms', nargs='+', required=True)
    ap.add_argument('--rank1-only', action='store_true')
    ap.add_argument('--out', required=True)
    ap.add_argument('--procs', type=int, default=8)
    ap.add_argument('--ppm', default=PPM_DEFAULT)
    a = ap.parse_args()
    df = pd.read_csv(a.ppm, low_memory=False)
    df = df[df.method.isin(a.arms)].copy()
    # EquiBind carries a rank==999 sentinel; its rank-1 is the lowest gnina/smina affinity.
    eq = df.method.str.startswith('equibind')
    if eq.any():
        e = df[eq].copy()
        key = pd.to_numeric(e.gnina_affinity, errors='coerce').where(e.method.str.endswith('_gnina'),
              pd.to_numeric(e.smina_affinity, errors='coerce'))
        e['_k'] = key
        e['eff_rank'] = e.sort_values('_k').groupby(['method','protein']).cumcount() + 1
        df.loc[eq, 'eff_rank'] = e['eff_rank']
    df.loc[~eq, 'eff_rank'] = df.loc[~eq, 'rank']
    df['eff_rank'] = df['eff_rank'].astype(int)
    if a.rank1_only:
        df = df[df.eff_rank == 1]
    cols = ['method','protein','pose_file','rank','eff_rank','rmsd','pb_rmsd','pb_valid','centroid_dist','bestfit_rmsd']
    groups = [(p, g[cols].to_dict('records')) for p, g in df.groupby('protein')]
    print(f'{len(df)} poses, {len(groups)} complexes, {a.procs} procs', flush=True)
    with Pool(a.procs) as pool:
        res = [r for chunk in pool.imap_unordered(work, groups, chunksize=2) for r in chunk]
    out = pd.DataFrame(res)
    out.to_csv(a.out, index=False)
    print('wrote', a.out, out.shape)
