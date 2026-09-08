import numpy as np, pandas as pd, itertools
from pathlib import Path
WF = Path('/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/wf')
SRC = '/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv'
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)
a = pd.read_csv(SRC, low_memory=False)
k = pd.read_csv(WF / 'nearest_copy_kabsch.csv', low_memory=False)
print('kabsch rows', k.shape, 'missing kabsch_nearest', k.kabsch_nearest.isna().sum(), 'missing kabsch_ref', k.kabsch_ref.isna().sum())

# ---------- (2) verification ----------
print('\n== VERIFICATION ==')
print('nearest idx agrees with probe (sym_rmsd_nearest vs rmsd_nearest_copy):', np.nanmax(np.abs(k.sym_rmsd_nearest - k.rmsd_nearest_copy)))
print('sym_rmsd_ref vs rmsd column max abs:', np.nanmax(np.abs(k.sym_rmsd_ref - k.rmsd)))
print('pb_rmsd_ref vs pb_rmsd column max abs:', np.nanmax(np.abs(k.pb_rmsd_ref - k.pb_rmsd)))
d = (k.kabsch_ref - k.bestfit_rmsd).abs()
print('kabsch_ref vs bestfit_rmsd: max', d.max(), ' n>1e-6:', (d > 1e-6).sum(), 'of', len(k))
single = k[k.n_copies == 1]
print('50 single-copy check: max |kabsch_ref-bestfit|', (single.kabsch_ref - single.bestfit_rmsd).abs().max(), 'n', len(single))
ref_rows = k[(k.n_copies > 1) & (k.nearest_copy_is_ref.astype(bool))]
print('multi-copy, nearest==ref: n', len(ref_rows), 'max |kabsch_nearest-bestfit|', (ref_rows.kabsch_nearest - ref_rows.bestfit_rmsd).abs().max())
flip = k[(k.n_copies > 1) & (~k.nearest_copy_is_ref.astype(bool))]
print('flipped rows:', len(flip))
print('copies_same_atom_order by protein:', k.groupby('protein').copies_same_atom_order.first().value_counts().to_dict())
bad = k[~k.copies_same_atom_order.astype(bool)].protein.unique()
print('proteins with differing atom order across copies:', list(bad))

# ---------- (1) PoseBusters multi-conformer semantics ----------
print('\n== POSEBUSTERS MULTI-CONFORMER PATH (choose_by=rmsd) ==')
m = k[(k.n_copies > 1) & k.pbmulti_kabsch.notna()].copy()
def kab_at(row, idx):
    vals = [float(v) for v in str(row.kabsch_all).split(';')]
    return vals[int(idx)] if 0 <= int(idx) < len(vals) else np.nan
m['kabsch_at_pbrmsd_argmin'] = [kab_at(r, r.pb_rmsd_min_idx) for r in m.itertuples()]
print('n with pbmulti:', len(m))
print('pbmulti_rmsd == pb_rmsd_min_all  max abs:', (m.pbmulti_rmsd - m.pb_rmsd_min_all).abs().max())
print('pbmulti_kabsch == kabsch AT the in-place-argmin conformer  max abs:', (m.pbmulti_kabsch - m.kabsch_at_pbrmsd_argmin).abs().max())
print('pbmulti_kabsch == kabsch_min_all (independent min)  max abs:', (m.pbmulti_kabsch - m.kabsch_min_all).abs().max(),
      ' n differ >1e-6:', ((m.pbmulti_kabsch - m.kabsch_min_all).abs() > 1e-6).sum())
print('n where kabsch_min_idx != nearest_idx (independent min picks another copy):', (m.kabsch_min_idx != m.nearest_idx).sum())
print('n where independent min < follow value by >0.05 A:', ((m.kabsch_nearest - m.kabsch_min_all) > 0.05).sum())

# ---------- build the full 53,651-row frame under both conventions ----------
cols = ['method','pose_file','kabsch_nearest','kabsch_min_all','pb_rmsd_nearest','cd_nearest','pb_rmsd_min_all','cd_min_all','nearest_idx','kabsch_min_idx']
f = a.merge(k[cols], on=['method','pose_file'], how='left')
one = f.n_copies == 1
f.loc[one, 'kabsch_nearest'] = f.loc[one, 'bestfit_rmsd']
f.loc[one, 'kabsch_min_all'] = f.loc[one, 'bestfit_rmsd']
f.loc[one, 'pb_rmsd_nearest'] = f.loc[one, 'pb_rmsd']
f.loc[one, 'pb_rmsd_min_all'] = f.loc[one, 'pb_rmsd']
f.loc[one, 'cd_nearest'] = f.loc[one, 'centroid_dist']
f.loc[one, 'cd_min_all'] = f.loc[one, 'centroid_dist']
assert f.kabsch_nearest.notna().all(), f.kabsch_nearest.isna().sum()
f['pb_valid'] = f.pb_valid.astype(bool)
f.to_csv(WF / 'per_pose_two_conventions.csv', index=False)
print('\nwrote per_pose_two_conventions.csv', f.shape)

CONV = {
  'reference': dict(rmsd='rmsd', kab='bestfit_rmsd', pbr='pb_rmsd', cd='centroid_dist'),
  'nearest_follow': dict(rmsd='rmsd_nearest_copy', kab='kabsch_nearest', pbr='pb_rmsd_nearest', cd='cd_nearest'),
  'nearest_independent': dict(rmsd='rmsd_nearest_copy', kab='kabsch_min_all', pbr='pb_rmsd_min_all', cd='cd_min_all'),
}
HEAD = {'AutoDock': 'autodock_mgltools_exh128_gnina', 'DiffDock': 'diffdock_smina', 'EquiBind': 'equibind_unguided_gnina'}
universe = sorted(f[f.method == 'diffdock_smina'].protein.unique()); N = len(universe); print('N complexes', N)

# ---------- (3a) Table 1 form and triple columns ----------
rows = []
for conv, c in CONV.items():
    for meth, g in f.groupby('method'):
        valid = g.pb_valid; near = g[c['rmsd']] <= 2.0
        form = (g[c['kab']] <= 1.0) & (g[c['rmsd']] < 1000)
        triple = valid & near & (g[c['kab']] <= 1.0)
        nc = lambda m_: g[m_].protein.nunique()
        rows.append(dict(convention=conv, method=meth, poses=len(g),
            valid_complexes=nc(valid), valid_poses=int(valid.sum()), valid_pct=round(100*valid.mean(),1),
            near_complexes=nc(near), near_poses=int(near.sum()), near_pct=round(100*near.mean(),1),
            form_complexes=nc(form), form_poses=int(form.sum()), form_pct=round(100*form.mean(),1),
            triple_complexes=nc(triple), triple_poses=int(triple.sum()), triple_pct=round(100*triple.mean(),1)))
t1 = pd.DataFrame(rows); t1.to_csv(WF / 'table1_form_two_conventions.csv', index=False)
print('\n== TABLE 1 blocks under three conventions ==\n', t1.to_string(index=False))

# ---------- (3b) triple gate at top-15 (selection quantity) & double gate ----------
def exist(g, mask, depth):
    return g[(g.eff_rank <= depth) & mask].protein.nunique()
rows = []
for conv, c in CONV.items():
    for meth, g in f.groupby('method'):
        for depth in (1, 15, 30):
            trip = g.pb_valid & (g[c['rmsd']] <= 2.0) & (g[c['kab']] <= 1.0) & (g[c['rmsd']] < 1000)
            dbl = g.pb_valid & (g[c['rmsd']] <= 2.0)
            k1 = g.pb_valid & (g[c['kab']] <= 1.0)
            rows.append(dict(convention=conv, method=meth, depth=depth,
                             triple_complexes=exist(g, trip, depth), triple_poses=int(((g.eff_rank <= depth) & trip).sum()),
                             double_complexes=exist(g, dbl, depth), kabsch1_valid_complexes=exist(g, k1, depth)))
tg = pd.DataFrame(rows); tg.to_csv(WF / 'triple_gate_depths_two_conventions.csv', index=False)
print('\n== TRIPLE GATE (pb_valid & rmsd<=2 & kabsch<=1) existence by depth ==\n', tg[tg.method.isin(HEAD.values())].to_string(index=False))

# ---------- (3c) Kabsch recovery table (Table 26 Kabsch half) ----------
THR = [0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5]
rows = []
for conv, c in CONV.items():
    for lab, meth in HEAD.items():
        g = f[f.method == meth]
        for depth in (1, 15, 30):
            s = g[g.eff_rank <= depth]
            r = dict(convention=conv, tool=lab, depth=depth)
            for t in THR:
                hit = s[s.pb_valid & (s[c['kab']] <= t)].protein.nunique()
                r[f't{t}'] = round(100.0 * hit / N, 1); r[f'n{t}'] = hit
            rows.append(r)
kr = pd.DataFrame(rows); kr.to_csv(WF / 'kabsch_recovery_depths_two_conventions.csv', index=False)
print('\n== KABSCH RECOVERY (pb_valid & kabsch<=t, existence) ==\n', kr[['convention','tool','depth'] + [f't{t}' for t in THR]].to_string(index=False))
print(kr[['convention','tool','depth','n1','n2']].to_string(index=False))

# ---------- (3c') form McNemar AutoDock vs DiffDock at 1 A ----------
from scipy.stats import binomtest
rows = []
for conv, c in CONV.items():
    for depth in (1, 15, 30):
        hits = {}
        for lab, meth in HEAD.items():
            g = f[(f.method == meth) & (f.eff_rank <= depth)]
            hits[lab] = set(g[g.pb_valid & (g[c['kab']] <= 1.0)].protein)
        for x, y in (('AutoDock','DiffDock'), ('AutoDock','EquiBind'), ('DiffDock','EquiBind')):
            b = len(hits[x] - hits[y]); cc = len(hits[y] - hits[x])
            p = binomtest(b, b + cc, 0.5).pvalue if b + cc else 1.0
            rows.append(dict(convention=conv, depth=depth, pair=f'{x} vs {y}', only_first=b, only_second=cc, p_exact_uncorrected=p,
                             n_first=len(hits[x]), n_second=len(hits[y])))
mc = pd.DataFrame(rows); mc.to_csv(WF / 'form_mcnemar_two_conventions.csv', index=False)
print('\n== FORM McNEMAR (1 A, pb_valid) ==\n', mc.to_string(index=False))

# ---------- (3d) filmstrip / Table 27 per_tool_depth under conventions ----------
def mech(r):
    return pd.cut(r, [-0.01, 1/3, 2/3, 1.01], labels=['placement-limited','mixed','form-limited'])
rows = []; clip = []
for conv, c in CONV.items():
    for depth in (1, 5, 15):
        tot = 0; clipped = 0
        for lab, meth in HEAD.items():
            g = f[(f.method == meth) & (f.eff_rank <= depth) & f.pb_valid].copy()
            g['form'] = g[c['kab']]; g['inplace'] = g[c['pbr']].where(g[c['pbr']].notna(), g[c['rmsd']])
            g = g.dropna(subset=['form','inplace'])
            g = g[g[c['cd']] <= 8.0]
            r = (g.form**2 / g.inplace**2).clip(0, 1); mm = mech(r).astype(str)
            tot += len(g); clipped += int(((g.inplace > 5) | (g.form > 5)).sum())
            rows.append(dict(convention=conv, tool=lab, depth=depth, n_poses=len(g), n_valid_complexes=g.protein.nunique(),
                coverage_pct=round(100*g.protein.nunique()/N, 1),
                inplace_median=round(g.inplace.median(),3), inplace_q1=round(g.inplace.quantile(.25),3), inplace_q3=round(g.inplace.quantile(.75),3),
                form_median=round(g.form.median(),3), form_q1=round(g.form.quantile(.25),3), form_q3=round(g.form.quantile(.75),3),
                median_r=round(r.median(),3),
                pct_placement_limited=round(100*(mm=='placement-limited').mean(),1), pct_mixed=round(100*(mm=='mixed').mean(),1),
                pct_form_limited=round(100*(mm=='form-limited').mean(),1)))
        clip.append(dict(convention=conv, depth=depth, cohort_poses=tot, poses_beyond_5A_either_axis=clipped))
fs = pd.DataFrame(rows); fs.to_csv(WF / 'filmstrip_table27_two_conventions.csv', index=False)
print('\n== TABLE 27 (filmstrip per_tool_depth) ==\n', fs.to_string(index=False))
print(pd.DataFrame(clip).to_string(index=False))

# ---------- (3e) H.4 near-native pose counts (rmsd<=2, pose level) ----------
rows = []
for conv, c in CONV.items():
    for lab, meth in HEAD.items():
        g = f[f.method == meth]
        nn = g[g[c['rmsd']] <= 2.0]
        rows.append(dict(convention=conv, tool=lab, near_native_poses=len(nn), near_native_valid=int(nn.pb_valid.sum()),
                         near_native_validity_pct=round(100*nn.pb_valid.mean(),2), all_validity_pct=round(100*g.pb_valid.mean(),2)))
h4 = pd.DataFrame(rows); h4.to_csv(WF / 'h4_near_native_pose_counts_two_conventions.csv', index=False)
print('\n== H.4 near-native pose strata ==\n', h4.to_string(index=False))

# ---------- (4) follow vs independent on flipped & all multi-copy ----------
mm_ = f[f.n_copies > 1]
print('\n== FOLLOW vs INDEPENDENT ==')
print('multi-copy poses:', len(mm_), ' independent picks other copy than in-place nearest:', (mm_.kabsch_min_idx != mm_.nearest_idx).sum())
print('form gate (<=1) verdict differs follow vs independent:', ((mm_.kabsch_nearest <= 1) != (mm_.kabsch_min_all <= 1)).sum())
print('form gate verdict differs reference vs follow:', ((mm_.bestfit_rmsd <= 1) != (mm_.kabsch_nearest <= 1)).sum(),
      ' of which gained:', ((mm_.bestfit_rmsd > 1) & (mm_.kabsch_nearest <= 1)).sum(), ' lost:', ((mm_.bestfit_rmsd <= 1) & (mm_.kabsch_nearest > 1)).sum())
print('form gate verdict differs reference vs independent:', ((mm_.bestfit_rmsd <= 1) != (mm_.kabsch_min_all <= 1)).sum(),
      ' lost:', ((mm_.bestfit_rmsd <= 1) & (mm_.kabsch_min_all > 1)).sum())
fl = f[(f.n_copies > 1) & (~f.nearest_copy_is_ref.astype(bool))]
print('flipped poses:', len(fl), ' |kabsch_nearest - bestfit| median/q90/max:', fl.eval('abs(kabsch_nearest-bestfit_rmsd)').quantile([.5,.9,1.0]).round(3).tolist())
print('flipped: form verdict changed ref->follow:', ((fl.bestfit_rmsd <= 1) != (fl.kabsch_nearest <= 1)).sum(), ' gained', ((fl.bestfit_rmsd > 1) & (fl.kabsch_nearest <= 1)).sum(), ' lost', ((fl.bestfit_rmsd <= 1) & (fl.kabsch_nearest > 1)).sum())
# per-copy conformer heterogeneity: kabsch between copies themselves is not computed here; proxy = spread of kabsch_all per pose
sp = k[k.n_copies > 1].kabsch_all.str.split(';').apply(lambda v: np.ptp([float(x) for x in v]))
print('per-pose spread (max-min) of Kabsch across copies: median', round(sp.median(),3), 'q90', round(sp.quantile(.9),3), 'max', round(sp.max(),3))
