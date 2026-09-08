"""Kabsch (best-fit) RMSD of each pose against the NEAREST crystallographic copy,
for every pose whose nearest copy (by in-place symmetry RMSD) is not the reference.
Uses ppc.posebusters_rmsd (PoseBusters check_rmsd) exactly as process_pair does."""
import sys, numpy as np, pandas as pd
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis')
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis/nearest_copy_probe')
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import posebusters_pose_comparison as ppc
from nearest_copy_rescoring import load_all_mols, BENCH
S = '/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/'

def work(args):
    prot, recs = args
    cdir = BENCH / prot
    crystal_h = Chem.RemoveHs(ppc.load_first_mol(cdir / f'{prot}_ligand.sdf'))
    copies = load_all_mols(cdir / f'{prot}_ligands.sdf')
    ref_c = crystal_h.GetConformer().GetPositions().mean(0)
    copy_c = [c.GetConformer().GetPositions().mean(0) for c in copies]
    ref_idx = int(np.argmin([np.linalg.norm(c - ref_c) for c in copy_c]))
    out = []
    for rec in recs:
        pose = ppc.load_first_mol(Path(rec['pose_file']))
        if pose is None:
            continue
        pose = ppc.reassign_template(pose, crystal_h) or pose
        r = [ppc.symmetry_rmsd(pose, c) for c in copies]
        k = int(np.nanargmin(r))
        pbr_n, kab_n, _ = ppc.posebusters_rmsd(pose, copies[k])
        pbr_r, kab_r, _ = ppc.posebusters_rmsd(pose, crystal_h)
        out.append({'method': rec['method'], 'pose_file': rec['pose_file'],
                    'nearest_idx': k, 'ref_idx': ref_idx, 'n_copies': len(copies),
                    'rmsd_nearest_re': float(r[k]), 'pb_rmsd_nearest': pbr_n,
                    'kabsch_nearest': kab_n, 'kabsch_ref_re': kab_r})
    return out

if __name__ == '__main__':
    df = pd.concat([pd.read_csv(S + 'nearest_copy_alldepths.csv', low_memory=False),
                    pd.read_csv(S + 'wf/nearest_copy_extra_arms.csv', low_memory=False),
                    pd.read_csv(S + 'wf/nearest_copy_guided_equibind_arms.csv', low_memory=False)])
    sub = df[~df.nearest_copy_is_ref.astype(bool)]
    print(len(df), 'poses total;', len(sub), 'with nearest copy != reference', flush=True)
    groups = [(p, g[['method', 'pose_file']].to_dict('records')) for p, g in sub.groupby('protein')]
    with Pool(12) as pool:
        res = [r for chunk in pool.imap_unordered(work, groups, chunksize=1) for r in chunk]
    out = pd.DataFrame(res)
    out.to_csv(S + 'wf/kabsch_vs_nearest_copy.csv', index=False)
    print('wrote', out.shape, flush=True)
