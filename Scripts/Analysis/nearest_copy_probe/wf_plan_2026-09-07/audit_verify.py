import sys, numpy as np, pandas as pd
sys.path.insert(0, "/home/manndo/master_dev/Scripts/Analysis")
from stats_utils import wilson_ci, newcombe_paired_diff_ci, mcnemar_power
from scipy.stats import binomtest
P="/tmp/claude-1000/-home-manndo-master-dev/55d7a717-f5e7-4353-8e73-45c4a49f6d87/scratchpad/nearest_copy_alldepths.csv"
df=pd.read_csv(P)
print(df.shape, df.method.unique())
df["r"]=df.eff_rank
df["cd_near"]=np.minimum(df.cd_ref, df.cd_min_alt.fillna(np.inf))
prots=sorted(df.protein.unique()); N=len(prots); print("N",N)
def sets(m, col, d, valid):
    s=df[(df.method==m)&(df.r<=d)&(df[col]<=2)]
    if valid: s=s[s.pb_valid==True] if s.pb_valid.dtype==bool else s[s.pb_valid.astype(str).str.lower()=="true"]
    return set(s.protein)
print("pb_valid dtype", df.pb_valid.dtype, df.pb_valid.unique()[:5])
def mcn(a,b):
    # a,b sets; returns won(a only), lost(b only), p
    w=len(a-b); l=len(b-a); p=binomtest(min(w,l), w+l, 0.5).pvalue if w+l>0 else 1.0
    return w,l,p
def pct(k): return 100*k/N
def wl(k):
    lo,hi=wilson_ci(k,N)[:2] if isinstance(wilson_ci(k,N),tuple) else (None,None)
    return f"{pct(k):.1f} [{100*lo:.1f}, {100*hi:.1f}]"
try:
    print("wilson test", wilson_ci(111,303))
except Exception as e: print("wilson err",e)
AD="autodock_mgltools_exh128_gnina"; ADR="autodock_mgltools_exh128"; DS="diffdock_smina"; DR="diffdock"; EQ="equibind_unguided_gnina"; UD="unidock2"
print("\n=== A. baselines ref vs nearest, VA and NN, d=1,5,10,15,30")
for m in [ADR,AD,DR,DS,EQ,UD]:
    for d in [1,5,10,15,30]:
        print(m,d,"VA ref",len(sets(m,"rmsd",d,1)),"near",len(sets(m,"rmsd_nearest_copy",d,1)),"| NN ref",len(sets(m,"rmsd",d,0)),"near",len(sets(m,"rmsd_nearest_copy",d,0)))
print("\n=== B. Wilson table nearest")
for m in [ADR,AD,DR,DS,EQ]:
    for d in [1,5,10,15]:
        print(m,d,"NN",wl(len(sets(m,"rmsd_nearest_copy",d,0))),"VA",wl(len(sets(m,"rmsd_nearest_copy",d,1))))
print("\n=== C. optimisation McNemar NN nearest (raw vs opt): won by opt / lost")
for raw,opt in [(ADR,AD),(DR,DS)]:
    for d in [1,5,10,15]:
        a=sets(raw,"rmsd_nearest_copy",d,0); b=sets(opt,"rmsd_nearest_copy",d,0)
        w,l,p=mcn(b,a); print(raw,"->",opt,d,f"{pct(len(a)):.1f}->{pct(len(b)):.1f}", "lost/gained",l,w,"p",f"{p:.3g}")
print("\n=== D. cross-tool NN nearest")
for x,y in [(AD,DR),(AD,EQ),(DS,EQ)]:
    for d in [1,5,10,15]:
        a=sets(x,"rmsd_nearest_copy",d,0); b=sets(y,"rmsd_nearest_copy",d,0); w,l,p=mcn(a,b)
        print(x,"vs",y,d,"won/lost",w,l,"p",f"{p:.3g}")
print("ref check AD vs DR d=1:", mcn(sets(AD,"rmsd",1,0),sets(DR,"rmsd",1,0)))
print("\n=== E. accurate-but-invalid (NN-VA complexes) ref -> nearest")
for m in [ADR,AD,DR,DS,EQ]:
    print(m,[ (len(sets(m,"rmsd",d,0))-len(sets(m,"rmsd",d,1)), len(sets(m,"rmsd_nearest_copy",d,0))-len(sets(m,"rmsd_nearest_copy",d,1))) for d in [1,5,10,15,30]])
    print("   nearest pct:",[wl(len(sets(m,"rmsd_nearest_copy",d,0))-len(sets(m,"rmsd_nearest_copy",d,1))) for d in [1,5,10,15]])
print("\n=== F. near-native pose counts / valid% ref -> nearest (:1170)")
for m in [ADR,AD,DR,DS,EQ]:
    for col in ["rmsd","rmsd_nearest_copy"]:
        s=df[(df.method==m)&(df[col]<=2)]; v=(s.pb_valid==True).sum() if s.pb_valid.dtype==bool else (s.pb_valid.astype(str).str.lower()=="true").sum()
        print(m,col,len(s),f"{100*v/len(s):.1f}","invalid",len(s)-v)
print("\n=== H. depth gain decomposition ref -> nearest (:1245)")
for m in [AD,DS,EQ]:
    for col in ["rmsd","rmsd_nearest_copy"]:
        NN1,NN15,VA1,VA15=sets(m,col,1,0),sets(m,col,15,0),sets(m,col,1,1),sets(m,col,15,1)
        g=NN15-NN1; inv=g-VA15; resc=(VA15-VA1)&NN1
        print(m,col,"+",len(g),"-",len(inv),"+",len(resc),"net",len(VA15)-len(VA1))
print("\n=== I/J. Newcombe & power (:1256/:1258)")
for col in ["rmsd","rmsd_nearest_copy"]:
    for d in [1,5,15,30]:
        a=np.array([p in sets(AD,col,d,1) for p in prots]); b=np.array([p in sets(DS,col,d,1) for p in prots])
        r95=newcombe_paired_diff_ci(a,b); r90=newcombe_paired_diff_ci(a,b,z=1.645)
        both=(a&b).sum(); neither=(~a&~b).sum(); disc=(a!=b).sum()
        w,l,p=mcn(sets(AD,col,d,1),sets(DS,col,d,1))
        print(col,d,"diff95",r95,"90",r90,"both",both,"neither",neither,"disc",disc,"p",f"{p:.3g}")
        if d==1:
            try: print("   power:",mcnemar_power(a,b))
            except Exception as e: print("   power err",e)
print("\n=== K. threshold-resolved nearest (:1301)")
thr=[0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5]
for m in [AD,DS,EQ]:
    for d in [1,15,30]:
        for col in ["rmsd","rmsd_nearest_copy"]:
            s=df[(df.method==m)&(df.r<=d)&(df.pb_valid==True)]
            row=[f"{100*len(set(s[s[col]<=t].protein))/N:.1f}" for t in thr]; print(m,d,col," ".join(row))
print("\n=== L. placement/form cohort (:1375) pb_valid & centroid<=8")
for m in [AD,DS,EQ]:
    for d in [1,5,15]:
        s=df[(df.method==m)&(df.r<=d)&(df.pb_valid==True)]
        for cdc in ["cd_ref","cd_near"]:
            t=s[s[cdc]<=8]; print(m,d,cdc,len(t),t.protein.nunique(),f"{100*t.protein.nunique()/N:.1f}")
print("\n=== M. ladder exh128 raw vs gnina top-15 VA nearest McNemar")
a=sets(ADR,"rmsd_nearest_copy",15,1); b=sets(AD,"rmsd_nearest_copy",15,1); print("raw",len(a),"sel",len(b),"lost/gained",len(a-b),len(b-a),"p",mcn(b,a)[2])
a=sets(ADR,"rmsd",15,1); b=sets(AD,"rmsd",15,1); print("REF raw",len(a),"sel",len(b),"lost/gained",len(a-b),len(b-a),"p",mcn(b,a)[2])
for d in range(1,31):
    a=sets(ADR,"rmsd_nearest_copy",d,1); b=sets(AD,"rmsd_nearest_copy",d,1); w,l,p=mcn(b,a)
    print(" d",d,"raw",len(a),"sel",len(b),"p",f"{p:.3g}", end=";")
print()
print("\n=== Q. :543 reconstruction PB-valid top-5 of 3 arms")
s=df[(df.method.isin([AD,DS,EQ]))&(df.r<=5)&(df.pb_valid==True)]
print("n",len(s),"median ref",s.rmsd.median(),"near",s.rmsd_nearest_copy.median(),"<=2 ref",f"{100*(s.rmsd<=2).mean():.1f}","near",f"{100*(s.rmsd_nearest_copy<=2).mean():.1f}")
for m in [AD,DS,EQ]:
    t=s[s.method==m]; print(m,"never within 10 ref",N-t[t.rmsd<=10].protein.nunique(),"near",N-t[t.rmsd_nearest_copy<=10].protein.nunique(), "complexes with any valid top5",t.protein.nunique())
print("\n=== R. reach 4 A centroid (all poses, rank-1)")
for m in [AD,DS,EQ]:
    for cdc in ["cd_ref","cd_near"]:
        t=df[df.method==m]; print(m,cdc,"all",t[t[cdc]<=4].protein.nunique(),"rank1",t[(t.r==1)&(t[cdc]<=4)].protein.nunique())
for cdc in ["cd_ref","cd_near"]:
    a=set(df[(df.method==AD)&(df[cdc]<=4)].protein); b=set(df[(df.method==DS)&(df[cdc]<=4)].protein); e=set(df[(df.method==EQ)&(df[cdc]<=4)].protein)
    print(cdc,"AD vs DS discordant",len(a-b),len(b-a),"pooled oracle",f"{100*len(a|b|e)/N:.1f}")
print("\n=== S. band table DiffDock smina minus raw, NN and NN&valid, pose-weighted within bands")
for col in ["rmsd","rmsd_nearest_copy"]:
    out=[]
    for lo,hi in [(1,10),(11,20),(21,30)]:
        r=df[(df.method==DR)&(df.r.between(lo,hi))]; s2=df[(df.method==DS)&(df.r.between(lo,hi))]
        nn=100*((s2[col]<=2).mean()-(r[col]<=2).mean()); nv=100*(((s2[col]<=2)&(s2.pb_valid==True)).mean()-((r[col]<=2)&(r.pb_valid==True)).mean())
        out.append((round(nn,1),round(nv,1)))
    print(col,out)
print("\n=== U. dropped ids in probe?", df[df.protein.isin(["7B2C_TP7","7D6O_MTE","7FRX_O88","7M31_TDR","8F4J_PHO"])].shape)
print("\n=== :267 equibind gnina NN pool ref/near", len(sets(EQ,"rmsd",30,0)), len(sets(EQ,"rmsd_nearest_copy",30,0)))
print("=== :246 rank-1 NN no validity ADR/AD ref->near", len(sets(ADR,"rmsd",1,0)),len(sets(ADR,"rmsd_nearest_copy",1,0)),len(sets(AD,"rmsd",1,0)),len(sets(AD,"rmsd_nearest_copy",1,0)))
print("=== flips in (1.5,2] at rank1 per arm:")
for m in [AD,DS,EQ]:
    t=df[(df.method==m)&(df.r==1)]; f=t[(t.rmsd>2)&(t.rmsd_nearest_copy<=2)]; print(m,len(f),((f.rmsd_nearest_copy>1.5)).sum())
