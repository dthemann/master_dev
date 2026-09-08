"""Score the AutoDock exh128+gnina poses of the five DiffDock-dropped ids against
the reference copy and against every copy in <ID>_ligands.sdf (thesis engine)."""
import sys, numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis')
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis/nearest_copy_probe')
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import posebusters_pose_comparison as ppc
from nearest_copy_rescoring import load_all_mols, BENCH
IDS = ['7B2C_TP7', '7D6O_MTE', '7FRX_O88', '7M31_TDR', '8F4J_PHO']
ROOT = Path('/home/manndo/master_dev/Dockings/vina_results_full_protein_vina_scoring_mgltools_exh128')
rows = []
for pid in IDS:
    cdir = BENCH / pid
    crystal_h = Chem.RemoveHs(ppc.load_first_mol(cdir / f'{pid}_ligand.sdf'))
    copies = load_all_mols(cdir / f'{pid}_ligands.sdf')
    ref_c = crystal_h.GetConformer().GetPositions().mean(0)
    copy_c = [c.GetConformer().GetPositions().mean(0) for c in copies]
    ref_idx = int(np.argmin([np.linalg.norm(c - ref_c) for c in copy_c]))
    log = pd.read_csv(ROOT / pid / 'mgl_tools' / 'optimization_log.csv')
    for _, r in log.iterrows():
        pose = ppc.load_first_mol(Path(r['optimized_file']))
        pose = ppc.reassign_template(pose, crystal_h) or pose
        rr = [ppc.symmetry_rmsd(pose, c) for c in copies]
        pc = Chem.RemoveHs(pose).GetConformer().GetPositions().mean(0)
        cd = [float(np.linalg.norm(pc - c)) for c in copy_c]
        rows.append({'protein': pid, 'n_copies': len(copies), 'autodock_rank': int(r['autodock_rank']),
                     'optimized_rank': int(r['optimized_rank']), 'cnn_affinity': r['cnn_affinity'],
                     'status': r['status'], 'rmsd_ref': ppc.symmetry_rmsd(pose, crystal_h),
                     'rmsd_refcopy': rr[ref_idx], 'rmsd_nearest_copy': float(np.nanmin(rr)),
                     'nearest_idx': int(np.nanargmin(rr)), 'ref_idx': ref_idx,
                     'cd_ref': cd[ref_idx], 'cd_nearest': float(min(cd)),
                     'pose_file': r['optimized_file']})
out = pd.DataFrame(rows)
out.to_csv('/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/wf/itt_dropped_ids_nearest_copy.csv', index=False)
print(out.groupby('protein').size())
for d in (1, 15, 30):
    s = out[out.optimized_rank <= d]
    g = s.groupby('protein').agg(min_ref=('rmsd_ref', 'min'), min_nearest=('rmsd_nearest_copy', 'min'),
                                 min_cd_ref=('cd_ref', 'min'), min_cd_nearest=('cd_nearest', 'min'))
    print(f'--- depth {d}'); print(g.round(2).to_string())
    print('  n within 2 A of ref:', int((g.min_ref <= 2).sum()), ' of nearest copy:', int((g.min_nearest <= 2).sum()))
print('closest pose to ref in top-15 pool over the five:', round(out[out.optimized_rank <= 15].rmsd_ref.min(), 2))
print('closest pose to nearest copy in top-15 pool:', round(out[out.optimized_rank <= 15].rmsd_nearest_copy.min(), 2))
