import pandas as pd
import numpy as np

CSV = '/home/manndo/master_dev/posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv'
df = pd.read_csv(CSV)
df['pb_valid'] = df['pb_valid'].astype(bool)

ad = df[df['method'] == 'autodock'].copy()


def best_of_top_d_repr(sub, d):
    """Generator logic: pick single lowest-RMSD pose among rank<=d, then check ITS validity."""
    s = sub.copy()
    s['_rk'] = pd.to_numeric(s['rank'], errors='coerce')
    cand = s[s['_rk'] <= d].dropna(subset=['rmsd'])
    if not len(cand):
        return None
    rep = cand.loc[cand.groupby(['protein', 'ligand'])['rmsd'].idxmin()]
    return rep


def any_gate(sub, d):
    """My gate: any pose in rank<=d with rmsd<=2 (near) / rmsd<=2 & valid (valid)."""
    s = sub.copy()
    s['_rk'] = pd.to_numeric(s['rank'], errors='coerce')
    cand = s[s['_rk'] <= d]
    near = np.any(cand['rmsd'].to_numpy() <= 2.0)
    valid = np.any((cand['rmsd'].to_numpy() <= 2.0) & cand['pb_valid'].to_numpy())
    return near, valid


for d in (15, 30):
    rep = best_of_top_d_repr(ad, d)
    rep_near = (rep['rmsd'] <= 2.0)
    rep_valid = rep_near & rep['pb_valid']
    print(f'=== autodock d={d} (GENERATOR: single lowest-RMSD representative) ===')
    print(f'  near  = {int(rep_near.sum())}')
    print(f'  valid = {int(rep_valid.sum())}')

    # my "any" gate per complex
    near_c = 0
    valid_c = 0
    per = {}
    for key, sub in ad.groupby(['protein', 'ligand']):
        n, v = any_gate(sub, d)
        per[key] = (n, v)
        near_c += int(n)
        valid_c += int(v)
    print(f'  ANY-gate near  = {near_c}')
    print(f'  ANY-gate valid = {valid_c}')

    # Diff: complexes where ANY-gate valid=True but representative valid=False
    rep_valid_keys = set(zip(rep.loc[rep_valid, 'protein'], rep.loc[rep_valid, 'ligand']))
    any_valid_keys = {k for k, (n, v) in per.items() if v}
    only_in_any = any_valid_keys - rep_valid_keys
    only_in_rep = rep_valid_keys - any_valid_keys
    print(f'  complexes valid under ANY but NOT representative: {sorted(only_in_any)}')
    print(f'  complexes valid under representative but NOT ANY: {sorted(only_in_rep)}')

    # Detail for the differing complexes
    for key in sorted(only_in_any):
        prot, lig = key
        sub = ad[(ad['protein'] == prot) & (ad['ligand'] == lig)].copy()
        sub['_rk'] = pd.to_numeric(sub['rank'], errors='coerce')
        top = sub[sub['_rk'] <= d].sort_values('_rk')
        print(f'\n  --- {prot}/{lig}  (top-{d}) ---')
        show = top[['rank', 'pose_name', 'rmsd', 'pb_valid']].copy()
        show['rmsd'] = show['rmsd'].round(3)
        # highlight lowest rmsd pose
        minidx = top['rmsd'].idxmin()
        print(show.to_string(index=False))
        r = top.loc[minidx]
        print(f'    lowest-RMSD pose in top-{d}: rank={int(r["rank"])} rmsd={r["rmsd"]:.3f} pb_valid={r["pb_valid"]}  <-- representative')
        # the near+valid pose(s)
        nv = top[(top['rmsd'] <= 2.0) & (top['pb_valid'])]
        if len(nv):
            rr = nv.sort_values('rmsd').iloc[0]
            print(f'    best near+VALID pose in top-{d}: rank={int(rr["rank"])} rmsd={rr["rmsd"]:.3f} pb_valid={rr["pb_valid"]}')
