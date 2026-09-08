import sys, pandas as pd, numpy as np
sys.path.insert(0,'/home/manndo/master_dev/Scripts/Analysis')
from stats_utils import newcombe_paired_diff_ci, mcnemar_power
df=pd.read_csv("/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv")
df['pb_valid']=df['pb_valid'].astype(str).str.lower().eq('true')
AD='autodock_mgltools_exh128_gnina'; DD='diffdock_smina'
prots=sorted(df.protein.unique())
def vec(arm,k,col):
    d=df[(df.method==arm)&(df.eff_rank<=k)&(df.pb_valid)&(df[col]<=2)]
    s=set(d.protein); return np.array([p in s for p in prots])
for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
    print("==",lab)
    for k in [1,5,15,30]:
        a=vec(DD,k,col); b=vec(AD,k,col)   # diff = AD - DD
        r95=newcombe_paired_diff_ci(a,b,1.96); r90=newcombe_paired_diff_ci(a,b,1.645)
        pw=mcnemar_power(a,b) if k==1 else None
        print(f"k={k}: diff {100*r95['diff']:+.1f} pp, 95% [{100*r95['lo']:+.1f}, {100*r95['hi']:+.1f}], 90% [{100*r90['lo']:+.1f}, {100*r90['hi']:+.1f}], smallest passing margin {100*max(abs(r90['lo']),abs(r90['hi'])):.1f}; n11={r95['n11']} n00={r95['n00']} disc={r95['n_discordant']}")
        if pw: print("   power:",{k_:(round(v,3) if isinstance(v,float) else v) for k_,v in pw.items()})
