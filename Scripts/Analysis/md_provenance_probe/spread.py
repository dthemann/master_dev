"""How much of the conformational spread present in bar3 do Fr300/Fr400/Fr499 actually cover?"""
import struct, json
import numpy as np
from pathlib import Path
BASE = Path("/mnt/c/Users/domin/OneDrive/00_Endless_Learning/FH Technikum/AI Engineering"
            "/Master Project/Orai1-Docking-MScProjekt-DMann/Orai/MDFiles")
NAT=182692; FRAMESZ=56+3*(8+4*NAT); BLK=8+4*NAT; NPROT=21026
ca=[]
with open(BASE/"Orai1WT-MDSnap-Fr300.pdb") as fh:
    i=0
    for line in fh:
        if line.startswith("ATOM"):
            if line[12:16].strip()=="CA": ca.append(i)
            i+=1
ca=np.array(ca)
p=BASE/"Orai1_WT_bar3.dcd"; size=p.stat().st_size
hdr=[h for h in range(4096) if (size-h)%FRAMESZ==0 and (size-h)//FRAMESZ>0][0]
nfr=(size-hdr)//FRAMESZ
frames=list(range(0,nfr,5))
X=[]
with open(p,"rb") as fh:
    for f in frames:
        base=hdr+f*FRAMESZ+56; xyz=np.empty((NPROT,3))
        for k in range(3):
            fh.seek(base+k*BLK+4); xyz[:,k]=np.frombuffer(fh.read(4*NPROT),dtype="<f4")
        X.append(xyz[ca])
X=np.array(X)
def krm(P,Q):
    Pc=P-P.mean(0); Qc=Q-Q.mean(0)
    V,S,Wt=np.linalg.svd(Pc.T@Qc); d=np.sign(np.linalg.det(V@Wt))
    R=V@np.diag([1,1,d])@Wt
    return float(np.sqrt((((Pc@R)-Qc)**2).sum()/len(P)))
n=len(frames); M=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n): M[i,j]=M[j,i]=krm(X[i],X[j])
iu=np.triu_indices(n,1)
sel=[frames.index(300),frames.index(400)]
# 499 is not a multiple of 5; use 495 as its proxy and note it
sel.append(frames.index(495))
selpairs=[M[a,b] for k,a in enumerate(sel) for b in sel[k+1:]]
# after the plateau begins (frame >= 250)
late=[i for i,f in enumerate(frames) if f>=250]
Ml=M[np.ix_(late,late)]; ilu=np.triu_indices(len(late),1)
out={
 "n_frames_sampled": n, "stride": 5,
 "all_pairs_ca_rmsd": {"mean": round(float(M[iu].mean()),3), "max": round(float(M[iu].max()),3),
                       "p95": round(float(np.percentile(M[iu],95)),3)},
 "plateau_only_pairs_frame_ge_250": {"n": len(late), "mean": round(float(Ml[ilu].mean()),3),
                       "max": round(float(Ml[ilu].max()),3), "p95": round(float(np.percentile(Ml[ilu],95)),3)},
 "chosen_three_pairs_using_495_for_499": [round(x,3) for x in selpairs],
 "chosen_three_max": round(float(max(selpairs)),3),
}
out["coverage_of_full_file_spread_pct"]=round(100*out["chosen_three_max"]/out["all_pairs_ca_rmsd"]["max"],1)
out["coverage_of_plateau_spread_pct"]=round(100*out["chosen_three_max"]/out["plateau_only_pairs_frame_ge_250"]["max"],1)
print(json.dumps(out,indent=1))
Path("/tmp/claude-1000/-home-manndo-master-dev/2aa6624d-3bb5-4674-81d2-2059c33af20b/scratchpad/spread.json").write_text(json.dumps(out,indent=1))
