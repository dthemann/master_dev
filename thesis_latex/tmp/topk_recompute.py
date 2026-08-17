import pandas as pd
import numpy as np

CSV = '/home/manndo/master_dev/posebusters_results/benchmark/dock/pose_comparison_report/per_pose_metrics.csv'
df = pd.read_csv(CSV)

# Normalize pb_valid to bool
df['pb_valid'] = df['pb_valid'].astype(bool)

METHODS = [
    'autodock',
    'diffdock',
    'diffdock_smina',
    'diffdock_gnina',
    'equibind_unguided_raw',
    'equibind_unguided_smina',
    'equibind_unguided_gnina',
]
KS = [1, 5, 10, 15, 30]

# universe of complexes = 303 unique protein+ligand
all_complexes = df.groupby(['protein', 'ligand']).ngroups
print('total complexes in file:', all_complexes)


def rank_method(sub, method):
    sub = sub.copy()
    if method in ('equibind_unguided_gnina',):
        # sort by gnina_affinity ascending (most negative = rank1), tie-break pose_name
        sub = sub.sort_values(['gnina_affinity', 'pose_name'], ascending=[True, True],
                              kind='mergesort', na_position='last')
    elif method in ('equibind_unguided_smina',):
        sub = sub.sort_values(['smina_affinity', 'pose_name'], ascending=[True, True],
                              kind='mergesort', na_position='last')
    else:
        # autodock, diffdock*, equibind_unguided_raw: sort by rank ascending, tie-break pose_name
        sub = sub.sort_values(['rank', 'pose_name'], ascending=[True, True],
                              kind='mergesort')
    return sub


rows = []
for method in METHODS:
    mdf = df[df['method'] == method]
    n_complexes_present = mdf.groupby(['protein', 'ligand']).ngroups
    # per complex
    per_complex = {}
    for (prot, lig), sub in mdf.groupby(['protein', 'ligand']):
        ranked = rank_method(sub, method)
        rmsd = ranked['rmsd'].to_numpy()
        valid = ranked['pb_valid'].to_numpy()
        per_complex[(prot, lig)] = (rmsd, valid)

    for k in KS:
        near_count = 0
        valid_count = 0
        for key, (rmsd, valid) in per_complex.items():
            topk_rmsd = rmsd[:k]
            topk_valid = valid[:k]
            near = np.any(topk_rmsd <= 2.0)
            vg = np.any((topk_rmsd <= 2.0) & topk_valid)
            if near:
                near_count += 1
            if vg:
                valid_count += 1
        rows.append({
            'method': method,
            'k': k,
            'near_count': near_count,
            'near_pct': round(100.0 * near_count / all_complexes, 2),
            'valid_count': valid_count,
            'valid_pct': round(100.0 * valid_count / all_complexes, 2),
            'n_present': n_complexes_present,
        })

res = pd.DataFrame(rows)
pd.set_option('display.width', 200)
pd.set_option('display.max_rows', 200)
print(res.to_string(index=False))

# Validation gate for autodock
print('\n=== AUTODOCK VALIDATION ===')
adv = res[res['method'] == 'autodock']
for k in KS:
    r = adv[adv['k'] == k].iloc[0]
    print(f"k={k}: near {r['near_count']}/{all_complexes} ({r['near_pct']}%), valid {r['valid_count']}/{all_complexes} ({r['valid_pct']}%)")
