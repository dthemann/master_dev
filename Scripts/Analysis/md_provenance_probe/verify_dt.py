"""Verification of the frame-interval estimate.
Guards against: (1) start-of-file artefact, (2) atom-block mis-mapping, (3) minimisation vs dynamics,
(4) reliance on a single tracer species."""
import struct, json
import numpy as np
from pathlib import Path

BASE = Path("/mnt/c/Users/domin/OneDrive/00_Endless_Learning/FH Technikum/AI Engineering"
            "/Master Project/Orai1-Docking-MScProjekt-DMann/Orai/MDFiles")
NAT = 182692; FRAMESZ = 56 + 3*(8+4*NAT); BLK = 8+4*NAT
def hdrinfo(p):
    s = p.stat().st_size
    for h in range(4096):
        if (s-h) % FRAMESZ == 0 and (s-h)//FRAMESZ > 0: return h, (s-h)//FRAMESZ
def rd(fh, hdr, fr, a0, n):
    off = hdr + fr*FRAMESZ
    fh.seek(off+4); uc = struct.unpack("<6d", fh.read(48))
    base = off+56; out = np.empty((n,3))
    for k in range(3):
        fh.seek(base + k*BLK + 4 + 4*a0)
        out[:,k] = np.frombuffer(fh.read(4*n), dtype="<f4")
    return out, np.array([uc[0],uc[2],uc[5]])
def mi(d,L): return d - L*np.round(d/L)

# ion indices from PSF
ions = {"POT":[], "CLA":[]}
with open(BASE/"hORAI1_WT_ARN.psf") as fh:
    for i,line in enumerate(fh):
        if i<7: continue
        f=line.split()
        if len(f)>=8 and f[0].isdigit() and f[1]=="IONS":
            ions.setdefault(f[3],[]).append(int(f[0])-1)
ION0 = 182433
out = {"ion_counts": {k:len(v) for k,v in ions.items()}}

def msd_series(p, hdr, start, n, a0, nat_read, stride3=True):
    pos, boxes = [], []
    with open(p,"rb") as fh:
        for f in range(start, start+n):
            xyz, L = rd(fh, hdr, f, a0, nat_read)
            pos.append(xyz[::3] if stride3 else xyz); boxes.append(L)
    pos=np.array(pos); L=np.array(boxes).mean(0)
    unw=np.empty_like(pos); unw[0]=pos[0]
    steps=[]
    for f in range(1,len(pos)):
        d=mi(pos[f]-pos[f-1],L); steps.append(float((d**2).sum(-1).mean())); unw[f]=unw[f-1]+d
    m={}
    for lag in [1,2,5,10,15,20]:
        if lag<len(unw):
            d=unw[lag:]-unw[:-lag]; m[lag]=round(float((d**2).sum(-1).mean()),3)
    return m, [round(s,3) for s in steps]

p = BASE/"Orai1_WT_bar3.dcd"; hdr,nfr = hdrinfo(p)
# (1) water, three disjoint windows and two disjoint water blocks
for tag,(start,a0) in {"win0_blockA":(0,62064), "win400_blockA":(400,62064),
                       "win200_blockB":(200,62064+30000*3), "win0_blockB":(0,62064+30000*3)}.items():
    m,st = msd_series(p,hdr,start,21,a0,3000*3)
    out[f"bar3_water_{tag}"] = {"msd_by_lag":m, "per_step_msd_first8":st[:8]}
# (2) ions as an independent tracer (K+ and Cl-, contiguous IONS block)
ionpos=[]; boxes=[]
with open(p,"rb") as fh:
    for f in range(0,21):
        xyz,L = rd(fh,hdr,f,ION0,259); ionpos.append(xyz); boxes.append(L)
ionpos=np.array(ionpos); L=np.array(boxes).mean(0)
unw=np.empty_like(ionpos); unw[0]=ionpos[0]
for f in range(1,len(ionpos)): unw[f]=unw[f-1]+mi(ionpos[f]-ionpos[f-1],L)
im={}
for lag in [1,2,5,10,20]:
    d=unw[lag:]-unw[:-lag]; im[lag]=round(float((d**2).sum(-1).mean()),3)
out["bar3_ion_msd_by_lag"]=im
out["ion_over_water_ratio_lag20"]= round(im[20]/out["bar3_water_win0_blockA"]["msd_by_lag"][20],3)

# (3) minimisation vs dynamics discriminator on fix1: is per-frame displacement constant or decaying?
pf = BASE/"Orai1_WT_fix1.dcd"; hdrf,nfrf = hdrinfo(pf)
mf, stf = msd_series(pf,hdrf,0,21,62064,3000*3)
out["fix1_water_msd_by_lag"]=mf; out["fix1_per_step_msd_first20"]=stf
mf2, stf2 = msd_series(pf,hdrf,400,21,62064,3000*3)
out["fix1_water_msd_by_lag_win400"]=mf2; out["fix1_per_step_msd_win400_first8"]=stf2[:8]
# fix1 long-lag saturation
pos=[];boxes=[]
with open(pf,"rb") as fh:
    for f in [0,50,100,200,300,400,500]:
        xyz,L=rd(fh,hdrf,f,62064,3000*3); pos.append(xyz[::3]); boxes.append(L)
pos=np.array(pos); L=np.array(boxes).mean(0)
out["fix1_msd_vs_frame0"]={fr: round(float((mi(pos[i]-pos[0],L)**2).sum(-1).mean()),3)
                           for i,fr in enumerate([0,50,100,200,300,400,500])}
# (4) bar3 volume autocorrelation - barostat coupling sanity
vols=[]
with open(p,"rb") as fh:
    for f in range(nfr):
        fh.seek(hdr+f*FRAMESZ+4); uc=struct.unpack("<6d",fh.read(48)); vols.append(uc[0]*uc[2]*uc[5])
v=np.array(vols); v=v-v.mean()
ac={lag: round(float((v[:-lag]*v[lag:]).mean()/ (v*v).mean()),3) for lag in [1,2,3,5,10,20]}
out["bar3_volume_autocorr"]=ac
out["bar3_volume_sd_nm3"]=round(float(np.std(np.array(vols)/1000)),3)

dest=Path("/tmp/claude-1000/-home-manndo-master-dev/2aa6624d-3bb5-4674-81d2-2059c33af20b/scratchpad/verify_dt.json")
dest.write_text(json.dumps(out,indent=1)); print(json.dumps(out,indent=1))
