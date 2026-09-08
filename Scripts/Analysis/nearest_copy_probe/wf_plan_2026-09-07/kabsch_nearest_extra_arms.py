"""Kabsch (superposed, symmetry-corrected, heavy-atom) RMSD of every benchmark pose in a
multi-copy complex against EVERY crystallographic copy in <ID>_ligands.sdf, using the
thesis loaders (load_first_mol, reassign_template), the thesis in-place engine
(symmetry_rmsd) and PoseBusters check_rmsd (via ppc.posebusters_rmsd) per copy.
Also runs PoseBusters' own multi-conformer path (one mol_true with all copies as
conformers, choose_by='rmsd') to show which conformer the kabsch value follows.
Read-only on the repository; writes only under scratchpad/wf/.
"""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis')
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import posebusters_pose_comparison as ppc
from posebusters.modules.rmsd import check_rmsd as pb_check_rmsd

BENCH = Path('/home/manndo/master_dev/Data/PoseBuster Benchmark Set')
WF = Path('/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/wf')
SRC = str(WF / 'nearest_copy_extra_arms.csv')

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

def same_topology(a, b):
    if a.GetNumAtoms() != b.GetNumAtoms(): return False
    if [x.GetAtomicNum() for x in a.GetAtoms()] != [x.GetAtomicNum() for x in b.GetAtoms()]: return False
    ea = sorted((min(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()), max(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx())) for bd in a.GetBonds())
    eb = sorted((min(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()), max(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx())) for bd in b.GetBonds())
    return ea == eb

def work(args):
    prot, recs = args
    cdir = BENCH / prot
    crystal = ppc.load_first_mol(cdir / f'{prot}_ligand.sdf')
    crystal_h = Chem.RemoveHs(crystal)
    copies = load_all_mols(cdir / f'{prot}_ligands.sdf')
    ref_c = crystal_h.GetConformer().GetPositions().mean(0)
    copy_c = [c.GetConformer().GetPositions().mean(0) for c in copies]
    ref_idx = int(np.argmin([np.linalg.norm(c - ref_c) for c in copy_c]))
    same_order = all(same_topology(copies[0], c) for c in copies[1:])
    # PoseBusters-style combined multi-conformer mol_true (AddConformer, no atom-order check)
    combined = None
    if len(copies) > 1:
        try:
            combined = Chem.Mol(copies[0])
            for c in copies[1:]:
                combined.AddConformer(c.GetConformer(), assignId=True)
        except Exception:
            combined = None
    out = []
    for rec in recs:
        pose = ppc.load_first_mol(Path(rec['pose_file']))
        row = {k: rec[k] for k in ('method','protein','pose_file','rank','eff_rank','pb_valid','rmsd','pb_rmsd',
                                   'bestfit_rmsd','centroid_dist','n_copies','rmsd_nearest_copy','nearest_copy_is_ref')}
        row.update(ref_idx=ref_idx, n_copies_loaded=len(copies), copies_same_atom_order=same_order)
        if pose is None:
            out.append(row); continue
        pose = ppc.reassign_template(pose, crystal_h) or pose
        pose_h = Chem.RemoveHs(pose)
        pc = pose_h.GetConformer().GetPositions().mean(0)
        r_sym, r_pb, k_pb, cd = [], [], [], []
        for c in copies:
            try: r_sym.append(ppc.symmetry_rmsd(pose, c))
            except Exception: r_sym.append(np.nan)
            a, b, _ = ppc.posebusters_rmsd(pose, c)
            r_pb.append(a); k_pb.append(b)
            cd.append(float(np.linalg.norm(pc - c.GetConformer().GetPositions().mean(0))))
        r_sym = np.array(r_sym, float); r_pb = np.array(r_pb, float); k_pb = np.array(k_pb, float)
        ni = int(np.nanargmin(np.nan_to_num(r_sym, nan=np.inf)))
        ki = int(np.nanargmin(np.nan_to_num(k_pb, nan=np.inf))) if np.isfinite(k_pb).any() else -1
        pbi = int(np.nanargmin(np.nan_to_num(r_pb, nan=np.inf))) if np.isfinite(r_pb).any() else -1
        row.update(
            sym_rmsd_ref=float(r_sym[ref_idx]), pb_rmsd_ref=float(r_pb[ref_idx]), kabsch_ref=float(k_pb[ref_idx]),
            nearest_idx=ni, sym_rmsd_nearest=float(r_sym[ni]), pb_rmsd_nearest=float(r_pb[ni]),
            kabsch_nearest=float(k_pb[ni]), cd_nearest=float(cd[ni]),
            kabsch_min_all=float(k_pb[ki]) if ki >= 0 else np.nan, kabsch_min_idx=ki,
            pb_rmsd_min_all=float(r_pb[pbi]) if pbi >= 0 else np.nan, pb_rmsd_min_idx=pbi,
            cd_min_all=float(min(cd)),
            kabsch_all=';'.join(f'{v:.4f}' for v in k_pb), sym_all=';'.join(f'{v:.4f}' for v in r_sym),
        )
        if combined is not None:
            try:
                res = pb_check_rmsd(pose_h, combined, rmsd_threshold=2.0, heavy_only=True)['results']
                row.update(pbmulti_rmsd=float(res['rmsd']), pbmulti_kabsch=float(res['kabsch_rmsd']),
                           pbmulti_centroid=float(res['centroid_distance']))
            except Exception as e:
                row.update(pbmulti_rmsd=np.nan, pbmulti_kabsch=np.nan, pbmulti_centroid=np.nan)
        out.append(row)
    return out

if __name__ == '__main__':
    df = pd.read_csv(SRC, low_memory=False)
    multi = df[df.n_copies > 1]
    single = df.iloc[0:0]
    sel = pd.concat([multi, single])
    print('rows to score:', len(sel), 'multi', len(multi), 'single-check', len(single), flush=True)
    tasks = [(p, g.to_dict('records')) for p, g in sel.groupby('protein')]
    tasks.sort(key=lambda t: -len(t[1]))
    rows = []
    with Pool(12) as pool:
        for i, res in enumerate(pool.imap_unordered(work, tasks)):
            rows.extend(res)
            if i % 10 == 0: print(f'  {i+1}/{len(tasks)} proteins, {len(rows)} rows', flush=True)
    out = pd.DataFrame(rows)
    out.to_csv(WF / 'nearest_copy_kabsch_extra_arms.csv', index=False)
    print('wrote', WF / 'nearest_copy_kabsch_extra_arms.csv', out.shape, flush=True)
