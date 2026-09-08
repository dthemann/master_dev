"""Numerical consequences of the nearest-copy convention for variant selection,
the exhaustiveness ladder, Table 1, the depth tables, the McNemar/TOST block
and the appendix top-k tables. Both conventions side by side."""
import sys, re, math, numpy as np, pandas as pd
sys.path.insert(0, '/home/manndo/master_dev/Scripts/Analysis')
from stats_utils import (mcnemar_exact, newcombe_paired_diff_ci, mcnemar_power,
                         tost_paired_proportions, holm, wilson_ci)
from scipy.stats import norm
S = '/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/'
W = S + 'wf/'
PPM = '/home/manndo/master_dev/posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/per_pose_metrics.csv'
REP = '/home/manndo/master_dev/posebusters_results/benchmark_matched_equibind/dock/pose_comparison_report/'
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40)

df = pd.concat([pd.read_csv(S + 'nearest_copy_alldepths.csv', low_memory=False),
                pd.read_csv(W + 'nearest_copy_extra_arms.csv', low_memory=False),
                pd.read_csv(W + 'nearest_copy_guided_equibind_arms.csv', low_memory=False)],
               ignore_index=True)
ppm = pd.read_csv(PPM, low_memory=False, usecols=['method', 'pose_file', 'ligand', 'pose_name', 'autodock_rank',
                                                  'optimized_rank', 'cnn_score', 'minimized_affinity',
                                                  'smina_affinity', 'gnina_affinity'])
assert not ppm.duplicated(['method', 'pose_file']).any()
df = df.merge(ppm, on=['method', 'pose_file'], how='left', validate='1:1')
kab = pd.read_csv(W + 'kabsch_vs_nearest_copy.csv', usecols=['method', 'pose_file', 'kabsch_nearest'])
df = df.merge(kab, on=['method', 'pose_file'], how='left', validate='1:1')
df['kab_ref'] = pd.to_numeric(df.bestfit_rmsd, errors='coerce')
df['kab_near'] = df.kabsch_nearest.fillna(df.kab_ref)
df['pb_valid'] = df.pb_valid.astype(bool)
# raw EquiBind arms: rank = generation index (the thesis's 'generation' kind), not the probe's NaN sort
EQ_RE = re.compile(r"_(\d{2,3})__ref", re.I)
raw_eq = df.method.str.startswith('equibind') & df.method.str.endswith('_raw')
gi = df.loc[raw_eq, 'pose_name'].astype(str).map(lambda s: int(EQ_RE.search(s).group(1)) if EQ_RE.search(s) else 10 ** 6)
df.loc[raw_eq, 'eff_rank'] = gi.groupby([df.loc[raw_eq, 'method'], df.loc[raw_eq, 'protein']]).rank(method='first').astype(int)
df['eff_rank'] = df.eff_rank.astype(int)
ALL = sorted(df.protein.unique()); N = len(ALL); assert N == 303
print('rows', len(df), 'arms', df.method.nunique(), 'complexes', N)
print('bare diffdock rank\\d+.sdf rows in diffdock:', int(df[df.method == 'diffdock'].pose_name.astype(str).str.contains(r'rank\d+\.sdf$').sum()))

RMSD = {'reference': 'rmsd', 'nearest': 'rmsd_nearest_copy'}
KAB = {'reference': 'kab_ref', 'nearest': 'kab_near'}
CONVS = ['reference', 'nearest']

def flags(sub, conv):
    r = sub[RMSD[conv]]; k = sub[KAB[conv]]
    near = r <= 2.0
    return pd.DataFrame({'near': near, 'valid_near': near & sub.pb_valid,
                         'triple': near & sub.pb_valid & (k < 1.0),
                         'form': k < 1.0, 'valid': sub.pb_valid}, index=sub.index)

def exist(sub, d, conv, gate, rank_col='eff_rank'):
    s = sub[sub[rank_col] <= d]
    f = flags(s, conv)[gate]
    return f.groupby(s.protein).any().reindex(ALL, fill_value=False)

def arm(m): return df[df.method == m]

def mcn(a, b):
    n10, n01, p = mcnemar_exact(a.astype(int).values, b.astype(int).values)
    return n10, n01, p

# ---------------------------------------------------------------- (1) ladder
LADDER = [('autodock_mgltools_exh18', 'Exhaustiveness 18', 18, 'raw'),
          ('autodock_mgltools', 'Exhaustiveness 32', 32, 'raw'),
          ('autodock_mgltools_gnina', '  + gnina rescoring (32)', 32, 'gnina'),
          ('autodock_mgltools_exh64', 'Exhaustiveness 64', 64, 'raw'),
          ('autodock_mgltools_exh64_gnina', '  + gnina rescoring (64)', 64, 'gnina'),
          ('autodock_mgltools_exh92', 'Exhaustiveness 92', 92, 'raw'),
          ('autodock_mgltools_exh128', 'Exhaustiveness 128', 128, 'raw'),
          ('autodock_mgltools_exh128_gnina', '  + gnina, selected (128)', 128, 'gnina')]
SEL = 'autodock_mgltools_exh128_gnina'
rows = []
for m, lab, e, kind in LADDER:
    for conv in CONVS:
        for gate in ('valid_near', 'near', 'triple'):
            rec = {'method': m, 'label': lab, 'exh': e, 'kind': kind, 'convention': conv, 'gate': gate}
            for d in (1, 5, 15, 30):
                rec[f'd{d}'] = int(exist(arm(m), d, conv, gate).sum())
            rows.append(rec)
ladder = pd.DataFrame(rows)
ladder.to_csv(W + 'ladder_both_conventions.csv', index=False)
print('\n=== (1) LADDER, PB-valid & <=2 A, existence gate (Table 8 order) ===')
print(ladder[ladder.gate == 'valid_near'].pivot_table(index=['exh', 'kind', 'label'], columns='convention', values=['d1', 'd5', 'd15', 'd30'], sort=False).to_string())
print('\n--- triple gate (PB-valid & <=2 A & Kabsch<1 A), Kabsch to the same copy ---')
print(ladder[ladder.gate == 'triple'].pivot_table(index=['exh', 'kind', 'label'], columns='convention', values=['d1', 'd5', 'd15', 'd30'], sort=False).to_string())
print('\n--- near-native only ---')
print(ladder[ladder.gate == 'near'].pivot_table(index=['exh', 'kind', 'label'], columns='convention', values=['d1', 'd5', 'd15', 'd30'], sort=False).to_string())

# selection rule: selected vs each of 7 at top-15, exact McNemar, Holm over 7; both gates; both conventions
sel_rows = []
for conv in CONVS:
    for gate in ('valid_near', 'triple'):
        S15 = exist(arm(SEL), 15, conv, gate)
        ps, recs = [], []
        for m, lab, e, kind in LADDER:
            if m == SEL: continue
            O = exist(arm(m), 15, conv, gate)
            n10, n01, p = mcn(S15, O)
            ps.append(p); recs.append({'convention': conv, 'gate': gate, 'other': m, 'selected_k15': int(S15.sum()),
                                       'other_k15': int(O.sum()), 'sel_only': n10, 'other_only': n01, 'p_raw': p})
        ph = holm(ps)
        for r, q in zip(recs, ph): r['p_holm'] = q
        sel_rows += recs
        # rank-1 family
        S1 = exist(arm(SEL), 1, conv, gate); ps1 = []
        for m, lab, e, kind in LADDER:
            if m == SEL: continue
            ps1.append(mcn(S1, exist(arm(m), 1, conv, gate))[2])
        n_surv = int((holm(ps1) < 0.05).sum())
        # depth sweep: depths where selected strictly leads all seven
        lead = []
        for d in range(1, 31):
            s = int(exist(arm(SEL), d, conv, gate).sum())
            others = [int(exist(arm(m), d, conv, gate).sum()) for m, *_ in LADDER if m != SEL]
            lead.append(s > max(others))
        lead_d = [d for d, l in zip(range(1, 31), lead) if l]
        print(f'\n[{conv} | {gate}] selected@15 = {int(S15.sum())}; rank-1 contrasts surviving Holm: {n_surv}/7; '
              f'strictly leads all seven at depths {lead_d[0] if lead_d else None}..{lead_d[-1] if lead_d else None} '
              f'(contiguous={lead_d == list(range(lead_d[0], lead_d[-1] + 1)) if lead_d else None}; all depths: {lead_d})')
sel = pd.DataFrame(sel_rows); sel.to_csv(W + 'ladder_selection_mcnemar_top15.csv', index=False)
print(sel.to_string())

# prose numbers of Appendix C
g128 = arm(SEL).copy()
def rank_by(sub, col, ascending):
    s = sub.copy(); s['_r'] = s.groupby('protein')[col].rank(method='first', ascending=ascending).astype(int); return s
prose = []
for conv in CONVS:
    for lab, col, asc in (('CNNaffinity (optimized_rank)', None, None), ('CNNscore', 'cnn_score', False),
                          ('minimised Vina energy', 'minimized_affinity', True), ('parent Vina order', 'autodock_rank', True)):
        s = g128 if col is None else rank_by(g128, col, asc)
        rc = 'eff_rank' if col is None else '_r'
        vals = {d: int(exist(s, d, conv, 'valid_near', rc).sum()) for d in (1, 15, 30)}
        prose.append({'convention': conv, 'item': f'exh128+gnina ordered by {lab}: PB-valid&<=2A', **{f'd{d}': v for d, v in vals.items()}})
    a = exist(g128, 1, conv, 'valid_near'); b = exist(rank_by(g128, 'cnn_score', False), 1, conv, 'valid_near', '_r')
    n10, n01, p = mcn(a, b)
    prose.append({'convention': conv, 'item': 'rank-1 CNNaffinity-only / CNNscore-only / McNemar p', 'd1': n10, 'd15': n01, 'd30': round(p, 3)})
    raw128 = arm('autodock_mgltools_exh128')
    R15, G15 = exist(raw128, 15, conv, 'valid_near'), exist(g128, 15, conv, 'valid_near')
    n10, n01, p = mcn(G15, R15)
    prose.append({'convention': conv, 'item': 'exh128 raw->gnina at top-15: gained (gnina only) / lost (raw only) / p', 'd1': n10, 'd15': n01, 'd30': round(p, 4)})
    prose.append({'convention': conv, 'item': 'rank-1 near-native (no validity): exh64 raw / exh92 raw / exh128 raw',
                  'd1': int(exist(arm('autodock_mgltools_exh64'), 1, conv, 'near').sum()), 'd15': int(exist(arm('autodock_mgltools_exh92'), 1, conv, 'near').sum()), 'd30': int(exist(raw128, 1, conv, 'near').sum())})
    prose.append({'convention': conv, 'item': 'rank-1 near-native (no validity): exh64+gnina / exh128+gnina',
                  'd1': int(exist(arm('autodock_mgltools_exh64_gnina'), 1, conv, 'near').sum()), 'd15': int(exist(g128, 1, conv, 'near').sum())})
    # triple over the whole pool per rung + marginal return per hour
    cost = pd.read_csv('/home/manndo/master_dev/posebusters_results/autodock_exhaustiveness_returns/exh_cost_per_arm.csv').set_index('exhaustiveness').total_wall_h
    pool = {e: int(exist(arm(m), 30, conv, 'triple').sum()) for m, lab, e, kind in LADDER if kind == 'raw'}
    prose.append({'convention': conv, 'item': 'triple over whole pool, raw rungs 18/32/64/92/128', 'd1': str(list(pool.values()))})
    steps = [(18, 32), (32, 64), (64, 92), (92, 128)]
    rates = [(pool[b] - pool[a]) / (cost[b] - cost[a]) for a, b in steps]
    prose.append({'convention': conv, 'item': 'triple-pool complexes per extra wall-hour, steps 18-32/32-64/64-92/92-128 (point est., bootstrap CI = recompute)', 'd1': str([round(r, 1) for r in rates])})
    prose.append({'convention': conv, 'item': 'first-vs-last measurable step contrast (point est.)', 'd1': round(rates[0] - rates[2], 1)})
    prose.append({'convention': conv, 'item': 'exh92->128: triple-pool delta / rank-1 valid_near delta', 'd1': pool[128] - pool[92],
                  'd15': int(exist(raw128, 1, conv, 'valid_near').sum()) - int(exist(arm('autodock_mgltools_exh92'), 1, conv, 'valid_near').sum())})
prose = pd.DataFrame(prose); prose.to_csv(W + 'appendix_c_prose_numbers.csv', index=False)
print('\n=== Appendix C prose numbers ==='); print(prose.to_string())

# ---------------------------------------------------------------- (2)/(3) oracle selection metric
def oracle_metric(sub, conv):
    r = sub[RMSD[conv]]
    idx = r.groupby([sub.protein, sub.ligand]).idxmin()
    o = sub.loc[idx]
    n = len(o)
    return (100 * float(((o[RMSD[conv]] <= 2.0) & o.pb_valid).mean()), 100 * float((o[RMSD[conv]] <= 2.0).mean()),
            int(((o[RMSD[conv]] <= 2.0) & o.pb_valid).sum()), n)
osum = pd.read_csv(REP + 'oracle_summary_all_variants.csv', index_col=0)
orows = []
for m in sorted(df.method.unique()):
    rec = {'method': m}
    for conv in CONVS:
        v, r2, k, n = oracle_metric(arm(m), conv)
        rec[f'oracle_pb_valid_and_rmsd2_%_{conv}'] = round(v, 4); rec[f'oracle_rmsd_le_2A_%_{conv}'] = round(r2, 4); rec[f'oracle_valid_near_complexes_{conv}'] = k
    rec['published_oracle_pb_valid_and_rmsd2_%'] = float(osum.loc[m, 'oracle_pb_valid_and_rmsd2_%']) if m in osum.index else np.nan
    orows.append(rec)
orc = pd.DataFrame(orows); orc.to_csv(W + 'oracle_selection_metric_both_conventions.csv', index=False)
print('\n=== (2)/(3) BEST_VARIANT_METRIC oracle_pb_valid_and_rmsd2_% (min-RMSD pose, then validity) ===')
print(orc.to_string())
for fam in ('diffdock', 'equibind'):
    for conv in CONVS:
        f = orc[orc.method.str.startswith(fam)].set_index('method')[f'oracle_pb_valid_and_rmsd2_%_{conv}']
        print(f'  {fam} winner [{conv}]: {f.idxmax()} = {f.max():.2f}  (ranking: {f.sort_values(ascending=False).round(2).to_dict()})')

# DiffDock stated selection rule (body :197): triple at top-15, complexes and qualifying poses
print('\n=== (2) DiffDock smina vs gnina under the stated rule and at rank-1 ===')
dd_rows = []
for conv in CONVS:
    A, B = arm('diffdock_smina'), arm('diffdock_gnina')
    for gate in ('triple', 'valid_near', 'near'):
        a15, b15 = exist(A, 15, conv, gate), exist(B, 15, conv, gate)
        n10, n01, p = mcn(a15, b15)
        pa = int(flags(A[A.eff_rank <= 15], conv)[gate].sum()); pb = int(flags(B[B.eff_rank <= 15], conv)[gate].sum())
        dd_rows.append({'convention': conv, 'gate': gate, 'depth': 15, 'smina_complexes': int(a15.sum()), 'gnina_complexes': int(b15.sum()),
                        'smina_poses': pa, 'gnina_poses': pb, 'smina_only': n10, 'gnina_only': n01, 'p': round(p, 4)})
    for gate in ('valid_near', 'near'):
        a1, b1 = exist(A, 1, conv, gate), exist(B, 1, conv, gate); n10, n01, p = mcn(a1, b1)
        dd_rows.append({'convention': conv, 'gate': gate, 'depth': 1, 'smina_complexes': int(a1.sum()), 'gnina_complexes': int(b1.sum()),
                        'smina_only': n10, 'gnina_only': n01, 'p': round(p, 4)})
    for gate in ('triple', 'near', 'form'):
        a30, b30 = exist(A, 30, conv, gate), exist(B, 30, conv, gate); n10, n01, p = mcn(a30, b30)
        dd_rows.append({'convention': conv, 'gate': gate, 'depth': 30, 'smina_complexes': int(a30.sum()), 'gnina_complexes': int(b30.sum()),
                        'smina_only': n10, 'gnina_only': n01, 'p': round(p, 4), 'agree_both': int((a30 & b30).sum())})
    # footnote :197 re-ranking by refiner score
    for m, col in (('diffdock_smina', 'smina_affinity'), ('diffdock_gnina', 'gnina_affinity')):
        s = rank_by(arm(m), col, True)
        dd_rows.append({'convention': conv, 'gate': f'near rank-1, re-ranked by {col}', 'depth': 1,
                        'smina_complexes' if 'smina' in m else 'gnina_complexes': int(exist(s, 1, conv, 'near', '_r').sum())})
    dd_rows.append({'convention': conv, 'gate': 'near rank-1, confidence order: diffdock raw / diffdock_smina / diffdock_gnina', 'depth': 1,
                    'smina_complexes': int(exist(arm('diffdock_smina'), 1, conv, 'near').sum()), 'gnina_complexes': int(exist(arm('diffdock_gnina'), 1, conv, 'near').sum()),
                    'smina_poses': int(exist(arm('diffdock'), 1, conv, 'near').sum())})
dd = pd.DataFrame(dd_rows); dd.to_csv(W + 'diffdock_refiner_selection_both_conventions.csv', index=False); print(dd.to_string())

# EquiBind 9 arms: stated rule (triple@15), near@30, valid anywhere
print('\n=== (3) EquiBind nine configurations ===')
eq_rows = []
for m in sorted(df.method.unique()):
    if not m.startswith('equibind'): continue
    rec = {'method': m, 'valid_anywhere': int(exist(arm(m), 30, 'reference', 'valid').sum())}
    for conv in CONVS:
        rec[f'triple_top15_{conv}'] = int(exist(arm(m), 15, conv, 'triple').sum())
        rec[f'near_top30_{conv}'] = int(exist(arm(m), 30, conv, 'near').sum())
        rec[f'valid_near_top15_{conv}'] = int(exist(arm(m), 15, conv, 'valid_near').sum())
        rec[f'valid_near_rank1_{conv}'] = int(exist(arm(m), 1, conv, 'valid_near').sum())
    eq_rows.append(rec)
eq = pd.DataFrame(eq_rows); eq.to_csv(W + 'equibind_nine_arms_both_conventions.csv', index=False); print(eq.to_string())
for conv in CONVS:
    G, Sm = arm('equibind_unguided_gnina'), arm('equibind_unguided_smina')
    n10, n01, p = mcn(exist(G, 15, conv, 'near'), exist(Sm, 15, conv, 'near'))
    print(f'  [{conv}] near@15 gnina vs smina: {int(exist(G,15,conv,"near").sum())} vs {int(exist(Sm,15,conv,"near").sum())}, disc {n10}/{n01}, p={p:.7f}')
    n10, n01, p = mcn(exist(G, 15, conv, 'triple'), exist(Sm, 15, conv, 'triple'))
    print(f'  [{conv}] triple@15 gnina vs smina: {int(exist(G,15,conv,"triple").sum())} vs {int(exist(Sm,15,conv,"triple").sum())}, disc {n10}/{n01}, p={p:.4f}')

# ---------------------------------------------------------------- (4) Table 1 pose-level
print('\n=== (4) Table 1 pose-level columns ===')
T1 = ['autodock_mgltools_exh128', 'autodock_mgltools_exh128_gnina', 'diffdock', 'diffdock_smina', 'diffdock_gnina',
      'equibind_unguided_raw', 'equibind_unguided_smina', 'equibind_unguided_gnina', 'equibind_fpocket_raw', 'equibind_p2rank_raw']
t1 = []
for m in T1:
    sub = arm(m); n = len(sub)
    for conv in CONVS:
        f = flags(sub, conv)
        rec = {'method': m, 'convention': conv, 'poses': n,
               'valid_complexes': int(sub.protein[f.valid].nunique()), 'valid_poses': int(f.valid.sum()), 'valid_pct': round(100 * f.valid.mean(), 1),
               'near_complexes': int(sub.protein[f.near].nunique()), 'near_poses': int(f.near.sum()), 'near_pct': round(100 * f.near.mean(), 1),
               'form_complexes': int(sub.protein[f.form].nunique()), 'form_poses': int(f.form.sum()), 'form_pct': round(100 * f.form.mean(), 1),
               'triple_complexes': int(sub.protein[f.triple].nunique()), 'triple_poses': int(f.triple.sum()), 'triple_pct': round(100 * f.triple.mean(), 1)}
        t1.append(rec)
t1 = pd.DataFrame(t1); t1.to_csv(W + 'table1_pose_level_both_conventions.csv', index=False); print(t1.to_string())

# ---------------------------------------------------------------- Table 2 + depth gain + cross-tool
print('\n=== Table 2 depth recovery, depth-gain decomposition, between-tool McNemar ===')
HEAD = [('AutoDock', SEL), ('DiffDock', 'diffdock_smina'), ('EquiBind', 'equibind_unguided_gnina')]
t2 = []
for conv in CONVS:
    E = {lab: {d: exist(arm(m), d, conv, 'valid_near') for d in (1, 15, 30)} for lab, m in HEAD}
    Nn = {lab: {d: exist(arm(m), d, conv, 'near') for d in (1, 15)} for lab, m in HEAD}
    p115 = [mcn(E[l][1], E[l][15])[2] for l, _ in HEAD]; p1530 = [mcn(E[l][15], E[l][30])[2] for l, _ in HEAD]
    h115, h1530 = holm(p115), holm(p1530)
    for i, (lab, m) in enumerate(HEAD):
        gained_near = Nn[lab][15] & ~Nn[lab][1]
        rec = {'convention': conv, 'tool': lab, 'rank1': int(E[lab][1].sum()), 'top15': int(E[lab][15].sum()), 'top30': int(E[lab][30].sum()),
               'gain_1_15': int(E[lab][15].sum() - E[lab][1].sum()), 'p_holm_1_15': h115[i],
               'gain_15_30': int(E[lab][30].sum() - E[lab][15].sum()), 'disc_15_30': mcn(E[lab][15], E[lab][30])[1], 'p_holm_15_30': h1530[i],
               'dg_gained_near': int(gained_near.sum()), 'dg_remain_invalid': -int((gained_near & ~E[lab][15]).sum()),
               'dg_valid_rescued': int((Nn[lab][1] & ~E[lab][1] & E[lab][15]).sum())}
        t2.append(rec)
    for (la, ma), (lb, mb) in [(HEAD[0], HEAD[1]), (HEAD[0], HEAD[2]), (HEAD[1], HEAD[2])]:
        for d in (1, 15, 30):
            n10, n01, p = mcn(E[la][d], E[lb][d])
            t2.append({'convention': conv, 'tool': f'{la} vs {lb} @ {d}', 'rank1': int(E[la][d].sum()), 'top15': int(E[lb][d].sum()),
                       'disc_15_30': f'{n10}/{n01}', 'p_holm_1_15': p})
t2 = pd.DataFrame(t2); t2.to_csv(W + 'table2_depth_recovery_both_conventions.csv', index=False); print(t2.to_string())

# appendix top-k Wilson (Table 9) and refinement McNemar (Table 10) and cross-tool (Table 11)
print('\n=== Appendix top-k (near% / validity-aware% with Wilson CI) ===')
T9 = [('AutoDock Vina (raw)', 'autodock_mgltools_exh128'), ('AutoDock Vina + gnina', SEL), ('DiffDock (raw)', 'diffdock'),
      ('DiffDock + smina', 'diffdock_smina'), ('DiffDock + gnina', 'diffdock_gnina'), ('EquiBind (raw)', 'equibind_unguided_raw'),
      ('EquiBind + gnina', 'equibind_unguided_gnina')]
t9 = []
for lab, m in T9:
    for k in (1, 5, 10, 15):
        rec = {'variant': lab, 'k': k}
        for conv in CONVS:
            for gate, tag in (('near', 'near'), ('valid_near', 'valid')):
                c = int(exist(arm(m), k, conv, gate).sum()); lo, hi = wilson_ci(c, N)
                rec[f'{tag}_{conv}'] = f'{100*c/N:.1f} [{100*lo:.1f}, {100*hi:.1f}] ({c})'
        t9.append(rec)
t9 = pd.DataFrame(t9); t9.to_csv(W + 'appendix_topk_wilson_both_conventions.csv', index=False); print(t9.to_string())

print('\n=== Appendix within-tool refinement McNemar (raw -> opt, near-native) ===')
T10 = [('AutoDock Vina + gnina', 'autodock_mgltools_exh128', SEL), ('DiffDock + smina', 'diffdock', 'diffdock_smina'),
       ('DiffDock + gnina', 'diffdock', 'diffdock_gnina'), ('EquiBind + gnina', 'equibind_unguided_raw', 'equibind_unguided_gnina')]
t10 = []
for conv in CONVS:
    recs, ps = [], []
    for lab, raw, opt in T10:
        for k in (1, 5, 10, 15):
            R, O = exist(arm(raw), k, conv, 'near'), exist(arm(opt), k, conv, 'near'); n10, n01, p = mcn(R, O)
            recs.append({'convention': conv, 'tool': lab, 'k': k, 'raw_pct': round(100 * R.sum() / N, 1), 'opt_pct': round(100 * O.sum() / N, 1),
                         'raw_only': n10, 'opt_only': n01, 'p_raw': p}); ps.append(p)
    for r, q in zip(recs, holm(ps)): r['p_holm'] = q
    t10 += recs
t10 = pd.DataFrame(t10); t10.to_csv(W + 'appendix_refinement_mcnemar_both_conventions.csv', index=False); print(t10.to_string())

print('\n=== Appendix cross-tool near-native McNemar family of sixteen ===')
T11 = [('AutoDock Vina + gnina vs DiffDock (raw)', SEL, 'diffdock'), ('AutoDock Vina + gnina vs DiffDock (gnina-opt)', SEL, 'diffdock_gnina'),
       ('DiffDock (gnina-opt) vs EquiBind (gnina-opt)', 'diffdock_gnina', 'equibind_unguided_gnina'), ('AutoDock Vina + gnina vs EquiBind (gnina-opt)', SEL, 'equibind_unguided_gnina')]
t11 = []
for conv in CONVS:
    recs, ps = [], []
    for lab, a, b in T11:
        for k in (1, 5, 10, 15):
            A, B = exist(arm(a), k, conv, 'near'), exist(arm(b), k, conv, 'near'); n10, n01, p = mcn(A, B)
            recs.append({'convention': conv, 'contrast': lab, 'k': k, 'a_pct': round(100 * A.sum() / N, 1), 'b_pct': round(100 * B.sum() / N, 1), 'a_only': n10, 'b_only': n01, 'p_raw': p}); ps.append(p)
    for r, q in zip(recs, holm(ps)): r['p_holm'] = q
    t11 += recs
t11 = pd.DataFrame(t11); t11.to_csv(W + 'appendix_crosstool_mcnemar_both_conventions.csv', index=False); print(t11.to_string())

# ---------------------------------------------------------------- (5)/(6) two-pipeline contrast, power, TOST
print('\n=== (5)/(6) AutoDock* vs DiffDock*: Newcombe, McNemar, power, TOST ===')
bc = []
for conv in CONVS:
    for d in (1, 5, 15, 30):
        A = exist(arm(SEL), d, conv, 'valid_near').astype(int).values; B = exist(arm('diffdock_smina'), d, conv, 'valid_near').astype(int).values
        # thesis orientation: AutoDock minus DiffDock  -> pass a=DiffDock (reference arm), b=AutoDock
        ci = newcombe_paired_diff_ci(B, A); ci90 = newcombe_paired_diff_ci(B, A, z=norm.ppf(0.95)); pw = mcnemar_power(B, A)
        _, _, p = mcnemar_exact(B, A)
        za, zb = norm.ppf(0.975), norm.ppf(0.80); psi = ci['n_discordant'] / N
        rec = {'convention': conv, 'k': d, 'AutoDock': int(A.sum()), 'DiffDock': int(B.sum()), 'both': ci['n11'], 'neither': ci['n00'],
               'AD_only': ci['n01'], 'DD_only': ci['n10'], 'discordant': ci['n_discordant'], 'discordant_pct': round(100 * psi, 1),
               'diff_pp': round(100 * ci['diff'], 2), 'ci95_lo': round(100 * ci['lo'], 2), 'ci95_hi': round(100 * ci['hi'], 2),
               'ci90_lo': round(100 * ci90['lo'], 2), 'ci90_hi': round(100 * ci90['hi'], 2), 'mcnemar_p': round(p, 4),
               'mde_pp_80pct': round(100 * pw['mde'], 2), 'observed_power': round(pw['observed_power'], 3),
               'n_needed': int(round(pw['n_needed'])) if not math.isnan(pw['n_needed']) else None,
               'mde_pp_at_428': round(100 * (za + zb) * math.sqrt(psi / 428), 2),
               'smallest_passing_margin_pp': math.ceil(max(abs(ci90['lo']), abs(ci90['hi'])) * 1000) / 10}
        for mg in (5, 10, 12, 15):
            rec[f'equiv_{mg}pp'] = tost_paired_proportions(B, A, margin=mg / 100)['equivalent']
        bc.append(rec)
bc = pd.DataFrame(bc); bc.to_csv(W + 'two_pipeline_contrast_power_tost_both_conventions.csv', index=False); print(bc.to_string())
print('\nhand check of the 9.8 pp: (z_.975+z_.80)*sqrt((112/303)/303) =', round(100 * (norm.ppf(0.975) + norm.ppf(0.8)) * math.sqrt((112 / 303) / 303), 3))
print('done')
