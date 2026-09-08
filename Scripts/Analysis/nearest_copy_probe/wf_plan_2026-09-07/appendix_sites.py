import pandas as pd, numpy as np
from scipy.stats import binomtest
from math import sqrt
from scipy.stats import norm
def proportion_confint(x,n,method="wilson"):
    z=norm.ppf(0.975); p=x/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den; return c-h,c+h
CSV="/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv"
df=pd.read_csv(CSV)
df['pb_valid']=df['pb_valid'].astype(str).str.lower().eq('true')
N=303
AD='autodock_mgltools_exh128_gnina'; ADR='autodock_mgltools_exh128'; DD='diffdock_smina'; DDR='diffdock'; EQ='equibind_unguided_gnina'; UD='unidock2'
names={AD:'AutoDock+gnina',ADR:'AutoDock raw',DD:'DiffDock+smina',DDR:'DiffDock raw',EQ:'EquiBind+gnina',UD:'Uni-Dock2'}
def rec(arm,k,col,valid,thr=2.0):
    d=df[(df.method==arm)&(df.eff_rank<=k)]
    m=d[col]<=thr
    if valid: m&=d.pb_valid
    s=d[m].groupby('protein').size()
    return set(s.index)
def wilson(x,n):
    lo,hi=proportion_confint(x,n,method='wilson'); return f"{100*x/n:.1f} [{100*lo:.1f}, {100*hi:.1f}]"
def mcn(a,b):
    A=set(a);B=set(b); ao=len(A-B); bo=len(B-A); n=ao+bo
    p=binomtest(min(ao,bo),n,0.5).pvalue if n>0 else float('nan')
    return ao,bo,p
print("=== A. Top-k table (:891-918) near-native% / validity-aware% ref -> nearest ===")
for arm in [ADR,AD,DDR,DD,EQ,UD]:
    for k in [1,5,10,15,30]:
        nn_r=len(rec(arm,k,'rmsd',False)); nn_n=len(rec(arm,k,'rmsd_nearest_copy',False))
        va_r=len(rec(arm,k,'rmsd',True)); va_n=len(rec(arm,k,'rmsd_nearest_copy',True))
        print(f"{names[arm]:16s} k={k:2d} NN {nn_r:3d}->{nn_n:3d} {wilson(nn_r,N)} -> {wilson(nn_n,N)} | VA {va_r:3d}->{va_n:3d} {wilson(va_r,N)} -> {wilson(va_n,N)} | acc-but-invalid {100*(nn_r-va_r)/N:.1f}->{100*(nn_n-va_n)/N:.1f} ({nn_r-va_r}->{nn_n-va_n})")
print("\n=== D. Within-tool raw->opt near-native (:945-952) ===")
for raw,opt in [(ADR,AD),(DDR,DD)]:
    for k in [1,5,10,15]:
        for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
            a=rec(raw,k,col,False); b=rec(opt,k,col,False); ao,bo,p=mcn(a,b)
            print(f"{names[raw]}->{names[opt]} k={k} {lab}: {100*len(a)/N:.1f}->{100*len(b)/N:.1f} discordant {ao}/{bo} p={p:.4g}")
print("\n=== Q. Cross-tool NN (:1024-1029) AD+gnina vs DD raw, AD+gnina vs EQ ===")
for other in [DDR,EQ,DD]:
    for k in [1,5,10,15,30]:
        for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
            a=rec(AD,k,col,False); b=rec(other,k,col,False); ao,bo,p=mcn(a,b)
            print(f"AD+gnina vs {names[other]} k={k} {lab}: {100*len(a)/N:.1f} vs {100*len(b)/N:.1f} AD-only/other-only {ao}/{bo} p={p:.4g}")
print("\n=== H. :1258 rank-1 VA AD vs DD both/neither/discordant ===")
for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
    a=rec(AD,1,col,True); b=rec(DD,1,col,True)
    print(f"{lab}: both {len(a&b)} neither {N-len(a|b)} discordant {len(a^b)} ({100*len(a^b)/N:.1f}%) AD-only {len(a-b)} DD-only {len(b-a)}; diff {100*(len(a)-len(b))/N:.1f} pp")
    for k in [5,15,30]:
        a=rec(AD,k,col,True); b=rec(DD,k,col,True); ao,bo,p=mcn(a,b)
        print(f"   k={k} {lab}: diff {100*(len(a)-len(b))/N:+.1f} pp, discordant {ao}/{bo}, exact McNemar p={p:.3g}")
print("\n=== E. Depth-gain decomposition (:1245-1248) rank-1 -> top-15 ===")
for arm in [AD,DD,EQ]:
    for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
        nn1=rec(arm,1,col,False); nn15=rec(arm,15,col,False); va1=rec(arm,1,col,True); va15=rec(arm,15,col,True)
        gain_nn=len(nn15-nn1); still_invalid=len((nn15-nn1)-va15); rescued=len((nn1-va1)&va15); net=len(va15)-len(va1)
        print(f"{names[arm]} {lab}: gaining NN +{gain_nn}, remain invalid -{still_invalid}, rescued +{rescued}, net +{net} (check {gain_nn-still_invalid+rescued})")
print("\n=== F. Validity cost NN-VA (:1194) ===")
for arm in [AD,DD,EQ,DDR]:
    for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
        costs=[len(rec(arm,k,col,False))-len(rec(arm,k,col,True)) for k in [1,5,10,15,30]]
        print(f"{names[arm]} {lab}: cost at k=1,5,10,15,30 = {costs}  pp {[round(100*c/N,2) for c in costs]}")
print("\n=== G. PB decomposition near-native poses (:1170-1177): count & valid% ===")
for arm in [ADR,AD,DDR,DD,EQ]:
    for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
        d=df[(df.method==arm)&(df[col]<=2)]
        print(f"{names[arm]} {lab}: NN poses {len(d)}, valid {100*d.pb_valid.mean():.1f}%, invalid n={int((~d.pb_valid).sum())}")
print("\n=== I. Threshold-resolved (:1301-1318) pb_valid & rmsd<=thr ===")
thr=[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5]
for arm in [AD,DD,EQ]:
    for k in [1,15,30]:
        for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
            row=[f"{100*len(rec(arm,k,col,True,t))/N:.1f}" for t in thr]
            print(f"{names[arm]} d={k} {lab}: "+" ".join(row))
print("\n=== J. Interaction pool (:543): PB-valid top-5 of three arms ===")
d=df[(df.method.isin([AD,DD,EQ]))&(df.eff_rank<=5)&(df.pb_valid)]
print(f"poses {len(d)}, complexes {d.protein.nunique()}; median rmsd ref {d.rmsd.median():.2f} nearest {d.rmsd_nearest_copy.median():.2f}; within 2A ref {100*(d.rmsd<=2).mean():.1f}% nearest {100*(d.rmsd_nearest_copy<=2).mean():.1f}%")
for arm in [AD,DD,EQ]:
    e=d[d.method==arm]
    cd_n=np.fmin(e.cd_ref, e.cd_min_alt.fillna(np.inf))
    print(f"  {names[arm]}: n={len(e)}; complexes within 10A (rmsd) ref {e[e.rmsd<=10].protein.nunique()} nearest {e[e.rmsd_nearest_copy<=10].protein.nunique()}; never-reached ref {N-e[e.rmsd<=10].protein.nunique()} nearest {N-e[e.rmsd_nearest_copy<=10].protein.nunique()}")
e10r=d[d.rmsd<=10]; e10n=d[d.rmsd_nearest_copy<=10]
print(f"  within-10A pool: share <=2A ref {100*(e10r.rmsd<=2).mean():.1f}% nearest {100*(e10n.rmsd_nearest_copy<=2).mean():.1f}%")
print("\n=== L. Placement/form cohort (:1375-1383): pb_valid & centroid<=8 ===")
for arm in [AD,DD,EQ]:
    for k in [1,5,15]:
        e=df[(df.method==arm)&(df.eff_rank<=k)&(df.pb_valid)]
        cdn=np.fmin(e.cd_ref, e.cd_min_alt.fillna(np.inf))
        r=e[e.centroid_dist<=8]; n_=e[cdn<=8]
        print(f"{names[arm]} ({k}): n poses ref {len(r)} (complexes {r.protein.nunique()}, cov {100*r.protein.nunique()/N:.1f}%) -> nearest {len(n_)} (complexes {n_.protein.nunique()}, cov {100*n_.protein.nunique()/N:.1f}%); in-place median pb_rmsd ref {r.pb_rmsd.median():.3f}")
print("\n=== R. 4A centroid reach over pool & rank-1 (:1392) ===")
for k in [1,30]:
    sets={}
    for arm in [AD,DD,EQ]:
        e=df[(df.method==arm)&(df.eff_rank<=k)]
        cdn=np.fmin(e.cd_ref, e.cd_min_alt.fillna(np.inf))
        sets[(arm,'ref')]=set(e[e.centroid_dist<=4].protein); sets[(arm,'nearest')]=set(e[cdn<=4].protein)
    for lab in ['ref','nearest']:
        ao,bo,p=mcn(sets[(AD,lab)],sets[(DD,lab)])
        un=sets[(AD,lab)]|sets[(DD,lab)]|sets[(EQ,lab)]
        print(f"k={k} {lab}: reach AD {len(sets[(AD,lab)])} DD {len(sets[(DD,lab)])} EQ {len(sets[(EQ,lab)])}; AD vs DD discordant {ao}/{bo} p={p:.3g}; pooled oracle {len(un)} ({100*len(un)/N:.1f}%)")
print("\n=== O. exh128 raw -> gnina top-15 VA discordant (:184 'gaining 17 losing 4') ===")
for k in [1,15,30]:
    for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
        a=rec(ADR,k,col,True); b=rec(AD,k,col,True); ao,bo,p=mcn(a,b)
        print(f"k={k} {lab}: raw {len(a)} gnina {len(b)} lost {ao} gained {bo} p={p:.3g}")
print("\n=== P. :246 raw/rescored rank-1 NN without validity ===")
for arm in [ADR,AD]:
    print(f"{names[arm]}: ref {len(rec(arm,1,'rmsd',False))} nearest {len(rec(arm,1,'rmsd_nearest_copy',False))}")
print("\n=== M. :267 EquiBind gnina unguided NN over full pool & VA top-15 (three-part needs bestfit) ===")
print(f"EQ NN all30 ref {len(rec(EQ,30,'rmsd',False))} nearest {len(rec(EQ,30,'rmsd_nearest_copy',False))}; VA top15 ref {len(rec(EQ,15,'rmsd',True))} nearest {len(rec(EQ,15,'rmsd_nearest_copy',True))}")
print("\n=== Fixed-rank band NN rates (:990-995) DD raw->smina per-rank pooled ===")
for band in [(1,10),(11,20),(21,30)]:
    for col,lab in [('rmsd','ref'),('rmsd_nearest_copy','nearest')]:
        r=df[(df.method==DDR)&(df.eff_rank.between(*band))]; o=df[(df.method==DD)&(df.eff_rank.between(*band))]
        print(f"DiffDock band {band} {lab}: NN raw {100*(r[col]<=2).mean():.1f} -> smina {100*(o[col]<=2).mean():.1f} (+{100*((o[col]<=2).mean()-(r[col]<=2).mean()):.1f}); NN&valid +{100*(((o[col]<=2)&o.pb_valid).mean()-((r[col]<=2)&r.pb_valid).mean()):.1f}")
