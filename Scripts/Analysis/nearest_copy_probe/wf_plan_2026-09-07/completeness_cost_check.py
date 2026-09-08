# Completeness auditor check: cost-per-qualifying-pose medians under nearest-copy.
# Reproduces the printed 137.2 / 49.4 / 2.4 (IQR 75.6-167.4 / 14.1-127.3 / 1.0-7.3), 271/479/38 and 186/478
# on the reference basis, then swaps near2 for pb_valid & rmsd_nearest_copy<=2.
import pandas as pd, numpy as np
from scipy import stats
eff=pd.read_csv('/home/manndo/master_dev/posebusters_results/benchmark_matched_equibind/docking_effort_gnina_v2_charged/per_complex_effort.csv')
eff['charged']=eff.cpu_core_s/32+eff.gpu_s
nc=pd.read_csv('/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv')
mp={'autodock':'autodock_mgltools_exh128_gnina','diffdock':'diffdock_smina','equibind':'equibind_unguided_gnina'}
for m,k in mp.items():
    e=eff[eff.method==m].set_index('cid'); d=nc[nc.method==k]
    q=d.assign(q=(d.pb_valid.astype(bool))&(d.rmsd_nearest_copy<=2)).groupby('protein').q.sum()
    e=e.join(q.rename('q_new')); s=e[e.q_new>0]; per=s.charged/s.q_new
    print(m,len(s),round(per.median(),1),round(per.quantile(.25),1),round(per.quantile(.75),1),s.q_new.median(),round(e.charged.sum()/len(s)),round(e.gpu_s.sum()/len(s)))
