"""Independent estimate of the DCD inter-frame interval from solvent and lipid displacement.
The DCD headers carry VMD placeholders, so the save frequency is estimated physically:
  water:  MSD(lag) = 6 D dt lag,  D(CHARMM TIP3) ~ 5.5e-5 cm2/s = 5.5e11 A2/s
  lipid:  MSD_2D(lag) = 4 D dt lag, D(lipid, 50% cholesterol) ~ 1e-8 .. 5e-8 cm2/s
"""
import struct, json
import numpy as np
from pathlib import Path

BASE = Path("/mnt/c/Users/domin/OneDrive/00_Endless_Learning/FH Technikum/AI Engineering"
            "/Master Project/Orai1-Docking-MScProjekt-DMann/Orai/MDFiles")
NAT = 182692
FRAMESZ = 56 + 3 * (8 + 4 * NAT)
BLK = 8 + 4 * NAT
MEMB0, MEMB1 = 21026, 62064          # 0-based, [start, end)
WAT0 = 62064                          # first water OH2, 0-based
NWAT_SAMPLE = 4000                    # waters sampled (contiguous block)

def header_bytes(p):
    size = p.stat().st_size
    for h in range(4096):
        if (size - h) % FRAMESZ == 0 and (size - h) // FRAMESZ > 0:
            return h, (size - h) // FRAMESZ
    raise RuntimeError("no header size found")

def read_block(fh, hdr, frame, a0, n):
    """Return (n,3) float64 for atoms [a0, a0+n) of the given frame, plus the unit cell."""
    off = hdr + frame * FRAMESZ
    fh.seek(off + 4); uc = struct.unpack("<6d", fh.read(48))
    base = off + 56
    out = np.empty((n, 3))
    for k in range(3):
        fh.seek(base + k * BLK + 4 + 4 * a0)
        out[:, k] = np.frombuffer(fh.read(4 * n), dtype="<f4")
    return out, np.array([uc[0], uc[2], uc[5]])

def minimage(d, L):
    return d - L * np.round(d / L)

res = {}
for name, nread in [("Orai1_WT_bar3.dcd", 26), ("Orai1_WT_fix1.dcd", 26)]:
    p = BASE / name
    hdr, nfr = header_bytes(p)
    # ---- water ----
    wat, boxes = [], []
    with open(p, "rb") as fh:
        for f in range(nread):
            xyz, L = read_block(fh, hdr, f, WAT0, NWAT_SAMPLE * 3)
            wat.append(xyz[::3])          # OH2 oxygens
            boxes.append(L)
    wat = np.array(wat); L = np.array(boxes).mean(0)
    lo, hi = wat[0].min(0), wat[0].max(0)
    # unwrap step by step
    unw = np.empty_like(wat); unw[0] = wat[0]
    step_frac = []
    for f in range(1, len(wat)):
        d = minimage(wat[f] - wat[f - 1], L)
        step_frac.append(float((np.abs(d) > L / 4).mean()))
        unw[f] = unw[f - 1] + d
    msd = {}
    for lag in [1, 2, 3, 5, 8, 12, 20, 25]:
        if lag >= len(unw): continue
        d = unw[lag:] - unw[:-lag]
        msd[lag] = float((d ** 2).sum(-1).mean())
    D_TIP3 = 5.5e11        # A^2/s
    dt_est = {lag: m / (6 * D_TIP3 * lag) for lag, m in msd.items()}
    res[name] = {
        "n_frames": nfr, "box_mean_A": [round(x, 2) for x in L],
        "water_coord_min": [round(x, 1) for x in lo], "water_coord_max": [round(x, 1) for x in hi],
        "water_msd_A2_by_lag": {k: round(v, 2) for k, v in msd.items()},
        "frac_step_components_over_L4": [round(x, 4) for x in step_frac[:6]],
        "dt_per_frame_s_from_water": {k: f"{v:.3e}" for k, v in dt_est.items()},
        "dt_per_frame_ps_from_water_D5.5": {k: round(v * 1e12, 2) for k, v in dt_est.items()},
        "dt_per_frame_ps_from_water_D4.5": {k: round(v * 1e12 * 5.5 / 4.5, 2) for k, v in dt_est.items()},
        "dt_per_frame_ps_from_water_D6.5": {k: round(v * 1e12 * 5.5 / 6.5, 2) for k, v in dt_est.items()},
    }
    print(name, "water done", flush=True)

# ---- lipid lateral MSD over the whole production run ----
p = BASE / "Orai1_WT_bar3.dcd"
hdr, nfr = header_bytes(p)
# phosphorus indices within MEMB, from the PSF
pidx = []
with open(BASE / "hORAI1_WT_ARN.psf") as fh:
    for i, line in enumerate(fh):
        if i < 7: continue
        f = line.split()
        if len(f) >= 8 and f[0].isdigit() and f[1] == "MEMB" and f[4] == "P":
            pidx.append(int(f[0]) - 1 - MEMB0)
pidx = np.array(pidx)
print("lipid P atoms:", len(pidx), flush=True)

stride = 10
frames = list(range(0, nfr, stride))
pos, boxes = [], []
with open(p, "rb") as fh:
    for f in frames:
        xyz, L = read_block(fh, hdr, f, MEMB0, MEMB1 - MEMB0)
        pos.append(xyz[pidx]); boxes.append(L)
pos = np.array(pos); L = np.array(boxes).mean(0)
zmid = pos[0][:, 2].mean()
upper = pos[0][:, 2] > zmid
unw = np.empty_like(pos); unw[0] = pos[0]
for f in range(1, len(pos)):
    unw[f] = unw[f - 1] + minimage(pos[f] - pos[f - 1], L)
# remove per-leaflet lateral centre-of-mass drift
for mask in (upper, ~upper):
    unw[:, mask, :2] -= unw[:, mask, :2].mean(axis=1, keepdims=True)
lip = {}
for lag in [1, 2, 5, 10, 20, 40, 49]:
    if lag >= len(unw): continue
    d = unw[lag:, :, :2] - unw[:-lag, :, :2]
    lip[lag * stride] = round(float((d ** 2).sum(-1).mean()), 3)
res["lipid_lateral_msd_A2_by_frame_lag"] = lip
res["lipid_note"] = ("MSD_2D = 4*D*t. For D in 1e-8..5e-8 cm2/s (= 1e9..5e9 A2/s), "
                     "t = MSD/(4D). Divide by the frame lag for per-frame dt.")
for lag, m in lip.items():
    res.setdefault("lipid_dt_per_frame_ps", {})[lag] = {
        "D_1e-8": round(m / (4 * 1e9) / lag * 1e12, 1),
        "D_3e-8": round(m / (4 * 3e9) / lag * 1e12, 1),
        "D_5e-8": round(m / (4 * 5e9) / lag * 1e12, 1),
    }

dest = Path("/tmp/claude-1000/-home-manndo-master-dev/2aa6624d-3bb5-4674-81d2-2059c33af20b/scratchpad/msd_probe.json")
dest.write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
print("WROTE", dest, flush=True)
