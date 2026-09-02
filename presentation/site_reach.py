import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

depth  = [1, 5, 10, 15]
series = [
    ("AutoDock*", [59, 79, 81, 83], "#1f77b4", "o", (10,  9)),
    ("DiffDock*", [54, 76, 79, 81], "#ff7f0e", "s", (10, -11)),
    ("EquiBind*", [43, 46, 47, 47], "#2ca02c", "^", (10,  0)),
]

fig, ax = plt.subplots(figsize=(6.8, 4.4), dpi=220)

for y, xt, txt, dash in ((95, 0.45, "pooled 3-tool oracle   95%", (0, (6, 4))),
                         (56.8, 6.2, "best blind selection rule   57%", (0, (2, 3)))):
    ax.axhline(y, color="#777777", ls=dash, lw=1.4, zorder=1)
    ax.text(xt, y + 2.0, txt, va="bottom", ha="left", fontsize=8.5,
            color="#555555", zorder=4,
            bbox=dict(fc="white", ec="none", pad=1.0))

for label, vals, colour, marker, off in series:
    ax.plot(depth, vals, "-", color=colour, marker=marker, lw=2.2,
            ms=8, mec="white", mew=1.6, zorder=3, label=label)
    ax.annotate(f"{vals[-1]}%", (depth[-1], vals[-1]),
                textcoords="offset points", xytext=off,
                va="center", fontsize=10.5, fontweight="bold", color=colour)

ax.set_xlim(0.2, 17.4)
ax.set_ylim(0, 104)
ax.set_xticks(depth)
ax.set_xticklabels([f"top-{d}" for d in depth])
ax.set_xlabel("Ranking depth inspected", fontsize=10.5)
ax.set_ylabel("% of 303 complexes reaching the crystal pocket", fontsize=10)
ax.set_title("Reaching the right pocket (cluster centroid within 4 Å)",
             fontsize=12.5, fontweight="bold", pad=10)
ax.grid(axis="y", color="#DDDDDD", lw=0.8)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color("#AAAAAA")
ax.tick_params(colors="#444444", labelsize=9.5)
ax.legend(loc="lower right", frameon=True, framealpha=0.95, fontsize=9.5,
          edgecolor="#CCCCCC", bbox_to_anchor=(1.0, 0.03))
fig.tight_layout()
fig.savefig("site_reach.png", bbox_inches="tight", facecolor="white")
print("ok")
