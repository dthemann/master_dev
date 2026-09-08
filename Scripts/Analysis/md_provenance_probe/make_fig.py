import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = "/tmp/claude-1000/-home-manndo-master-dev/2aa6624d-3bb5-4674-81d2-2059c33af20b/scratchpad/"
tr = json.load(open(S + "traj_probe.json"))
ms = json.load(open(S + "msd_probe.json"))

PROD = "#1f77b4"   # Orai1_WT_bar3
EQUI = "#ff7f0e"   # Orai1_WT_fix1
GRID = dict(color="#d9d9d9", lw=0.6, alpha=0.9)
DT_PROD = 0.66     # ps per frame, from bulk-water displacement (>45 A from the bilayer)
DT_EQUI = 0.075

def label_panels(axes):
    for ax, L in zip(axes, "ABCDEF"):
        ax.annotate(f"({L})", xy=(0, 1), xycoords="axes fraction",
                    xytext=(-30, 8), textcoords="offset points",
                    ha="left", va="bottom", fontsize=12, fontweight="bold")

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))
for ax in axes:
    ax.grid(True, **GRID); ax.set_axisbelow(True)
    for s in ("top", "right"): ax.spines[s].set_visible(False)

# (A) protein relaxation in the production-stage file
r = tr["Orai1_WT_bar3.dcd"]["ca_rmsd_vs_frame0"]
ax = axes[0]
ax.plot(range(len(r)), r, color=PROD, lw=1.6)
for f, name in [(300, "Fr300"), (400, "Fr400"), (499, "Fr499")]:
    ax.plot([f], [r[f]], "o", ms=8, color=PROD, mec="white", mew=1.5, zorder=5)
    ax.annotate(name, xy=(f, r[f]), xytext=(0, 11), textcoords="offset points",
                ha="center", fontsize=9, color="#333333")
ax.axvspan(250, 499, color=PROD, alpha=0.07, lw=0)
ax.annotate("plateau", xy=(375, 0.35), ha="center", fontsize=9, color="#52514e")
ax.set_xlabel("Frame Number in Orai1_WT_bar3")
ax.set_ylabel("Alpha-Carbon Root-Mean-Square Deviation\nfrom Frame 0 (Angstrom)")
ax.set_title("Protein relaxes, then holds", fontsize=10.5, color="#333333", pad=36)
ax.set_ylim(0, 1.85)
sec = ax.secondary_xaxis("top", functions=(lambda x: x * DT_PROD, lambda t: t / DT_PROD))
sec.set_xlabel("Estimated Elapsed Simulation Time (picoseconds)", fontsize=9)

# (B) box volume: constant volume in the restrained run, barostat in the later run
ax = axes[1]
for key, c, nm in [("Orai1_WT_fix1.dcd", EQUI, "Orai1_WT_fix1"),
                   ("Orai1_WT_bar3.dcd", PROD, "Orai1_WT_bar3")]:
    v = tr[key]["volume_nm3_trace"]
    ax.plot(range(len(v)), v, color=c, lw=1.6, label=nm)
    sd = float(np.std(v))
    ax.annotate(f"{nm}\nstandard deviation {sd:.2f}", xy=(len(v) - 1, v[-1]),
                xytext=(-4, -34 if c == EQUI else 14),
                textcoords="offset points", ha="right", fontsize=9, color=c, fontweight="bold")
ax.set_xlabel("Frame Number")
ax.set_ylabel("Simulation Box Volume (cubic nanometres)")
ax.set_title("Constant volume, then a fluctuating barostat", fontsize=10.5, color="#333333")
ax.set_ylim(1650, 2360)

# (C) water displacement sets the interval between saved frames
ax = axes[2]
for key, c, nm, dt in [("Orai1_WT_fix1.dcd", EQUI, "Orai1_WT_fix1", DT_EQUI),
                       ("Orai1_WT_bar3.dcd", PROD, "Orai1_WT_bar3", DT_PROD)]:
    d = ms[key]["water_msd_A2_by_lag"]
    lags = sorted(int(k) for k in d)
    vals = [d[str(k)] for k in lags]
    ax.plot(lags, vals, "o-", color=c, lw=1.8, ms=7, mec="white", mew=1.2, label=nm)
    j = len(lags) - 3 if c == PROD else len(lags) - 1
    ax.annotate(f"{nm}\n~{dt*1000:.0f} fs per frame", xy=(lags[j], vals[j]),
                xytext=(-4, -40) if c == PROD else (-6, -34), textcoords="offset points",
                ha="right", fontsize=9, color=c, fontweight="bold")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("Separation Between Frames (number of frames)")
ax.text(0.5, -0.29, "curves show all sampled water; the quoted interval uses bulk water only",
        transform=ax.transAxes, ha="center", fontsize=8, color="#52514e")
ax.set_ylabel("Mean Squared Displacement of\nWater Oxygen Atoms (square Angstrom)")
ax.set_title("Water motion dates the frame interval", fontsize=10.5, color="#333333")
ax.set_ylim(0.08, 90)

label_panels(axes)
fig.suptitle("Evidence for the stage and the time resolution of the supplied Orai1 trajectories",
             fontsize=13, y=1.005)
fig.tight_layout()
out = S + "orai_md_frame_provenance.png"
fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
print("wrote", out)
