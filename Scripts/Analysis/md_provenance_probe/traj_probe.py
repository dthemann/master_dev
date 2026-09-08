"""Raw DCD probe: box dimensions + protein CA RMSD per frame, without loading full frames.
Protein occupies the first 21026 atoms of the PSF/DCD (verified against the snapshot PDBs)."""
import struct, json, sys
import numpy as np
from pathlib import Path

BASE = Path("/mnt/c/Users/domin/OneDrive/00_Endless_Learning/FH Technikum/AI Engineering"
            "/Master Project/Orai1-Docking-MScProjekt-DMann/Orai/MDFiles")
NATOMS = 182692
NPROT  = 21026
FRAMESZ = 56 + 3 * (8 + 4 * NATOMS)

# CA indices from a snapshot PDB (same atom order as the DCD's protein block)
ca_idx = []
with open(BASE / "Orai1WT-MDSnap-Fr300.pdb") as fh:
    i = 0
    for line in fh:
        if line.startswith("ATOM"):
            if line[12:16].strip() == "CA":
                ca_idx.append(i)
            i += 1
ca_idx = np.array(ca_idx)
print(f"protein atoms in PDB: {i}, CA atoms: {len(ca_idx)}", flush=True)

def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Pr = Pc @ R
    return float(np.sqrt(((Pr - Qc) ** 2).sum() / len(P)))

def probe(name):
    p = BASE / name
    size = p.stat().st_size
    nfr = None
    for hdr in range(0, 4096):
        if (size - hdr) % FRAMESZ == 0 and (size - hdr) // FRAMESZ > 0:
            nfr = (size - hdr) // FRAMESZ
            header = hdr
            break
    print(f"{name}: size={size} header={header} n_frames={nfr}", flush=True)
    boxes, cas = [], []
    with open(p, "rb") as fh:
        for f in range(nfr):
            off = header + f * FRAMESZ
            fh.seek(off + 4)
            uc = struct.unpack("<6d", fh.read(48))
            boxes.append(uc)
            base = off + 56
            xyz = np.empty((NPROT, 3), dtype=np.float32)
            for k in range(3):
                fh.seek(base + k * (8 + 4 * NATOMS) + 4)
                xyz[:, k] = np.frombuffer(fh.read(4 * NPROT), dtype="<f4")
            cas.append(xyz[ca_idx].astype(np.float64))
            if f % 100 == 0:
                print(f"  {name} frame {f}", flush=True)
    return np.array(boxes), cas

out = {}
data = {}
for name in ["Orai1_WT_fix1.dcd", "Orai1_WT_bar3.dcd"]:
    boxes, cas = probe(name)
    data[name] = (boxes, cas)
    vols = boxes[:, 0] * boxes[:, 2] * boxes[:, 5]
    ref = cas[0]
    rmsd_vs_first = [kabsch_rmsd(c, ref) for c in cas]
    out[name] = {
        "n_frames": len(cas),
        "unitcell_first": list(boxes[0]),
        "unitcell_last": list(boxes[-1]),
        "a_trace": [round(float(x), 3) for x in boxes[:, 0]],
        "c_trace": [round(float(x), 3) for x in boxes[:, 5]],
        "volume_nm3_trace": [round(float(v) / 1000.0, 2) for v in vols],
        "ca_rmsd_vs_frame0": [round(x, 3) for x in rmsd_vs_first],
    }

# cross-file continuity
b = data["Orai1_WT_bar3.dcd"][1]; f1 = data["Orai1_WT_fix1.dcd"][1]
out["continuity"] = {
    "fix1_last_vs_bar3_first": round(kabsch_rmsd(f1[-1], b[0]), 4),
    "fix1_last_vs_bar3_first_norot": round(float(np.sqrt(((f1[-1]-b[0])**2).sum()/len(b[0]))), 4),
    "fix1_first_vs_bar3_first": round(kabsch_rmsd(f1[0], b[0]), 4),
    "bar3_299_vs_300": round(kabsch_rmsd(b[299], b[300]), 4),
}
# pairwise among the four used frames
sel = {"Fr0_fix1_0": f1[0], "Fr300": b[300], "Fr400": b[400], "Fr499": b[499]}
pw = {}
ks = list(sel)
for i in range(len(ks)):
    for j in range(i+1, len(ks)):
        pw[f"{ks[i]}|{ks[j]}"] = round(kabsch_rmsd(sel[ks[i]], sel[ks[j]]), 3)
out["four_frame_pairwise_ca_rmsd"] = pw

dest = Path("/tmp/claude-1000/-home-manndo-master-dev/2aa6624d-3bb5-4674-81d2-2059c33af20b/scratchpad/traj_probe.json")
dest.write_text(json.dumps(out, indent=1))
print("WROTE", dest, flush=True)
