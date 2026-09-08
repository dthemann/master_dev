import pandas as pd, numpy as np
from scipy import stats
from math import sqrt
import itertools
p='/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv'
df=pd.read_csv(p)
AD='autodock_mgltools_exh128_gnina'; DD='diffdock_smina'; EB='equibind_unguided_gnina'; ADR='autodock_mgltools_exh128'; DDR='diffdock'
arms=[AD,DD,EB]
prots=sorted(set(df[df.method==DD].protein)&set(df[df.method==AD].protein)&set(df[df.method==EB].protein))
print('n complexes',len(prots))
print('multi-copy among 303:', (df[df.method==AD].groupby('protein').n_copies.first()>1).sum())
df['cd_min']=df[['cd_ref','cd_min_alt']].min(axis=1)

def recov(arm, d, col, gate=True):
    s=df[(df.method==arm)&(df.eff_rank<=d)]
    ok=(s[col]<=2.0)
    if gate: ok=ok&s.pb_valid
    hit=s[ok].groupby('protein').size()
    return pd.Series([1 if pr in hit.index else 0 for pr in prots], index=prots)

def mcnemar_exact(b,c):
    n=b+c
    if n==0: return 1.0
    k=min(b,c)
    return min(1.0, 2*stats.binom.cdf(k,n,0.5)) if b!=c else 1.0

def wilson(k,n,z=1.959964):
    ph=k/n; den=1+z*z/n; cen=(ph+z*z/(2*n))/den; half=z*sqrt(ph*(1-ph)/n+z*z/(4*n*n))/den
    return cen-half, cen+half

def newcombe_paired(a,b,c,d,z=1.959964):
    # Newcombe 1998 method 10 for difference of paired proportions p1-p2; a=both yes,b=only1,c=only2,d=neither
    n=a+b+c+d
    p1=(a+b)/n; p2=(a+c)/n
    l1,u1=wilson(a+b,n,z); l2,u2=wilson(a+c,n,z)
    # phi correction
    A=(a+b)*(a+c)*(c+d)*(b+d)
    if A==0: phi=0.0
    else:
        phi=(a*d-b*c)/sqrt(A)
        if a*d-b*c>0: phi=max(0,(a*d-b*c-n/2))/sqrt(A) if False else phi
    # Newcombe: phi estimate with continuity? standard method 10 uses phi as above (no cc) unless ad-bc>0 then (ad-bc - n/2)... implement both
    dlt=p1-p2
    lo=dlt-sqrt((p1-l1)**2-2*phi*(p1-l1)*(u2-p2)+(u2-p2)**2)
    hi=dlt+sqrt((u1-p1)**2-2*phi*(u1-p1)*(p2-l2)+(p2-l2)**2)
    return dlt,lo,hi

def newcombe10(a,b,c,d,z=1.959964):
    n=a+b+c+d
    p1=(a+b)/n; p2=(a+c)/n
    l1,u1=wilson(a+b,n,z); l2,u2=wilson(a+c,n,z)
    A=(a+b)*(a+c)*(c+d)*(b+d)
    if A==0: phi=0.0
    else:
        num=a*d-b*c
        if num>0: num=num-n/2  # Newcombe's continuity-adjusted phi
        phi=num/sqrt(A)
    dlt=p1-p2
    lo=dlt-sqrt((p1-l1)**2-2*phi*(p1-l1)*(u2-p2)+(u2-p2)**2)
    hi=dlt+sqrt((u1-p1)**2-2*phi*(u1-p1)*(p2-l2)+(p2-l2)**2)
    return dlt,lo,hi

def holm(ps):
    idx=np.argsort(ps); m=len(ps); adj=np.empty(m); run=0
    for r,i in enumerate(idx):
        v=min(1.0,(m-r)*ps[i]); run=max(run,v); adj[i]=run
    return adj

def cochranQ(mat):
    # mat: n x k binary
    k=mat.shape[1]; T=mat.sum(axis=1); G=mat.sum(axis=0); N=mat.sum()
    Q=(k-1)*(k*(G**2).sum()-N**2)/(k*N-(T**2).sum())
    return Q, stats.chi2.sf(Q,k-1)

for col,label in [('rmsd','REFERENCE'),('rmsd_nearest_copy','NEAREST')]:
    print('\n=====',label,'=====')
    R={}
    for d in [1,5,10,15,30]:
        for arm in arms: R[(arm,d)]=recov(arm,d,col)
        for arm in [ADR,DDR]: R[(arm,d)]=recov(arm,d,col)
    for arm in arms+[ADR,DDR]:
        print(arm, {d:(int(R[(arm,d)].sum()), round(100*R[(arm,d)].mean(),1)) for d in [1,5,10,15,30]})
    # Table 5 gain column
    for arm in arms:
        g1=R[(arm,15)]-R[(arm,1)]; g2=R[(arm,30)]-R[(arm,15)]
        b1=int((g1==1).sum()); c1=int((g1==-1).sum()); b2=int((g2==1).sum()); c2=int((g2==-1).sum())
        lo,hi=wilson(b1-c1,303)
        print(f'  {arm} 1->15 +{b1-c1} (+{100*(b1-c1)/303:.1f}%) p={mcnemar_exact(b1,c1):.3g} Wilson[{100*lo:.1f},{100*hi:.1f}] | 15->30 +{b2-c2} (+{100*(b2-c2)/303:.1f}%) b={b2} c={c2} p={mcnemar_exact(b2,c2):.3g}')
    # pairwise per depth
    for d in [1,5,15,30]:
        ps=[]; rows=[]
        for x,y in [(AD,DD),(AD,EB),(DD,EB)]:
            a=int(((R[(x,d)]==1)&(R[(y,d)]==1)).sum()); b=int(((R[(x,d)]==1)&(R[(y,d)]==0)).sum())
            c=int(((R[(x,d)]==0)&(R[(y,d)]==1)).sum()); dd=303-a-b-c
            pval=mcnemar_exact(b,c); ps.append(pval)
            dl,lo,hi=newcombe10(a,b,c,dd)
            rows.append((x,y,b,c,pval,100*dl,100*lo,100*hi))
        adj=holm(ps)
        mat=np.column_stack([R[(arm,d)].values for arm in arms]); Q,pq=cochranQ(mat)
        print(f' depth {d}: CochranQ={Q:.1f} p={pq:.2e}')
        for r,pa in zip(rows,adj):
            print(f'   {r[0][:8]} vs {r[1][:8]}: disc {r[2]}/{r[3]} exact p={r[4]:.3g} holm={pa:.3g} diff={r[5]:+.1f} CI[{r[6]:+.1f},{r[7]:+.1f}]')
        # 90% interval for AD vs DD
        a=int(((R[(AD,d)]==1)&(R[(DD,d)]==1)).sum()); b=int(((R[(AD,d)]==1)&(R[(DD,d)]==0)).sum()); c=int(((R[(AD,d)]==0)&(R[(DD,d)]==1)).sum())
        dl,lo,hi=newcombe10(a,b,c,303-a-b-c,z=1.644854)
        print(f'   AD-DD 90% Newcombe: {100*dl:+.1f} [{100*lo:+.1f},{100*hi:+.1f}]')
    # validity-gate removal: NN-only minus NN&valid
    for d in [1,5,15,30]:
        out=[]
        for arm in arms:
            nn=recov(arm,d,col,gate=False).sum(); out.append(int(nn-R[(arm,d)].sum()))
        print(f' validity gate removes at depth {d}:',out, 'NN-only:',[int(recov(arm,d,col,gate=False).sum()) for arm in arms])
    # rank-1 failures containing no qualifying pose anywhere (top-30)
    for arm in arms:
        f1=(R[(arm,1)]==0); none=(R[(arm,30)]==0)
        print(f' {arm}: rank-1 failures {int(f1.sum())}, of which none anywhere {int((f1&none).sum())} = {100*(f1&none).sum()/f1.sum():.1f}%')
    # MDD for rank-1 AD vs DD at 80% power (McNemar, discordant n)
    b=int(((R[(AD,1)]==1)&(R[(DD,1)]==0)).sum()); c=int(((R[(AD,1)]==0)&(R[(DD,1)]==1)).sum()); nd=b+c
    # asymptotic: detectable |b-c|/303 such that z=(b-c)/sqrt(nd) => (1.96+0.84)*sqrt(nd)/303
    print(f' MDD rank-1 (nd={nd}): {100*(1.959964+0.841621)*sqrt(nd)/303:.1f} pp')
    # Table 1 near-nativeness block (all poses, no validity), complexes and poses
    for arm in [ADR,AD,DDR,DD,EB]:
        s=df[df.method==arm]; nn=(s[col]<=2.0)
        both=nn&s.pb_valid
        print(f' T1 {arm}: poses {len(s)}, NN comp {s[nn].protein.nunique()}, NN poses {int(nn.sum())} ({100*nn.mean():.1f}%), PBvalid&NN poses {int(both.sum())}, PBvalid&NN comp {s[both].protein.nunique()}')
    # qualifying poses restricted to 303 analysed
    for arm in [AD,DD,EB,ADR]:
        s=df[(df.method==arm)&(df.protein.isin(prots))]; q=(s[col]<=2.0)&s.pb_valid
        print(f' cost {arm}: qualifying poses {int(q.sum())} over {s[q].protein.nunique()} complexes; median qualifying per succeeding complex {s[q].groupby("protein").size().median()}')
    # 4A reach pose-level any pose, per tool and pooled
    cdcol='cd_ref' if col=='rmsd' else 'cd_min'
    reach={}
    for arm in arms:
        s=df[(df.method==arm)&(df.protein.isin(prots))]
        m=s.groupby('protein')[cdcol].min().reindex(prots)
        reach[arm]=m
        r1=s[s.eff_rank==1].groupby('protein')[cdcol].min().reindex(prots)
        print(f' 4A reach {arm}: any-pose {int((m<=4).sum())} ({100*(m<=4).mean():.1f}%) median best {m[m<=4].median():.2f}; rank-1 {int((r1<=4).sum())} ({100*(r1<=4).mean():.1f}%)')
    pooled=pd.concat(reach.values(),axis=1).min(axis=1)
    print(f' pooled 3-tool 4A reach {int((pooled<=4).sum())} ({100*(pooled<=4).mean():.1f}%) median {pooled[pooled<=4].median():.2f}')
    pooled2=pd.concat([reach[AD],reach[DD]],axis=1).min(axis=1)
    print(f' AD+DD 4A reach {int((pooled2<=4).sum())} ({100*(pooled2<=4).mean():.1f}%) median {pooled2[pooled2<=4].median():.2f}')
    # rank-1 near-site 8A trim for AutoDock
    s=df[(df.method==AD)&(df.eff_rank==1)&(df.protein.isin(prots))]
    print(f' AD rank-1 within 8A of site: {int((s[cdcol]<=8).sum())}/{len(s)} = {100*(s[cdcol]<=8).mean():.1f}%')
    # 6,169 PB-valid poses within 8A at top-15 (all three tools)?
    s=df[(df.method.isin(arms))&(df.eff_rank<=15)&(df.protein.isin(prots))&(df.pb_valid)]
    print(f' PB-valid top-15 poses within 8A (3 tools): {int((s[cdcol]<=8).sum())}, of which beyond 5A on either axis(inplace>5 or bestfit>5): {int(((s[cdcol]<=8)&((s[col]>5)|(s.bestfit_rmsd>5))).sum())}')
    # DiffDock raw vs smina fixed-rank NN rates (no gate) at ranks 1..5? line 843: fixed-rank rates
    for k in [1,2,3,5]:
        r=[]
        for arm in [DDR,DD]:
            s=df[(df.method==arm)&(df.eff_rank==k)&(df.protein.isin(prots))]
            r.append(100*(s[col]<=2).mean())
        print(f' DiffDock fixed rank {k} NN rate raw {r[0]:.1f} smina {r[1]:.1f} delta {r[1]-r[0]:+.1f}')
    # AutoDock raw vs gnina rank-1 & top-15 (valid&NN) McNemar
    for d in [1,5,15]:
        x=R[(ADR,d)]; y=R[(AD,d)]
        b=int(((y==1)&(x==0)).sum()); c=int(((y==0)&(x==1)).sum())
        print(f' AD gnina vs raw depth {d}: {int(x.sum())}->{int(y.sum())} disc {b}/{c} p={mcnemar_exact(b,c):.3g}')
