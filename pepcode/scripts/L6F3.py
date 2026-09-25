"""Lecture 18, Fig. 3: P-wave ray paths in PREM for a surface source, showing the shadow zone
produced by the drop in P-wave speed at the core-mantle boundary."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem, sphray
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})

fig, ax = plt.subplots(figsize=(7.4, 6.0))
for r, col in [(sphray.R_TOP, "0.85"), (prem.R_CMB, "0.72"), (prem.R_ICB, "0.6")]:
    ax.add_patch(plt.Circle((0, 0), r, color=col, zorder=0))
for r in [sphray.R_TOP, prem.R_CMB, prem.R_ICB]:
    ax.plot(r * np.sin(np.linspace(0, 2 * np.pi, 400)), r * np.cos(np.linspace(0, 2 * np.pi, 400)), "k", lw=0.7)

def draw(phase, ps, colour):
    for p in ps:
        path = sphray.ray_path(p, phase, n_per_layer=12)
        if path is None:
            continue
        D, r = path[:, 0], path[:, 1]
        ax.plot(r * np.sin(D), r * np.cos(D), color=colour, lw=0.7)

draw("P", np.linspace(258, 770, 22), "tab:blue")
draw("PKP", np.linspace(120, 252, 9), "tab:red")
draw("PKIKP", np.linspace(5, 108, 7), "tab:green")

# shadow zone: between the largest P distance and the smallest PKP distance
dP = sphray.travel_time_curve("P", n=300)[0].max()
dK = sphray.travel_time_curve("PKP", n=300)[0].min()
th = np.radians(np.linspace(dP, dK, 100))
R = sphray.R_TOP * 1.04
ax.plot(R * np.sin(th), R * np.cos(th), "k", lw=3, solid_capstyle="butt")
mid = np.radians(0.5 * (dP + dK))
ax.text(1.09 * sphray.R_TOP * np.sin(mid), 1.09 * sphray.R_TOP * np.cos(mid), "P-wave shadow zone",
        ha="left", va="center")
for d in [0, 30, 60, 90, 150, 180]:
    t = np.radians(d)
    ax.text(1.06 * sphray.R_TOP * np.sin(t), 1.06 * sphray.R_TOP * np.cos(t), r"$%d^\circ$" % d, ha="center", va="center", fontsize=9)
ax.plot([0], [sphray.R_TOP], "k*", ms=10, zorder=5)
ax.text(-300, sphray.R_TOP + 250, "source", ha="right", va="center", fontsize=9)
ax.text(5050, -2050, "P", color="tab:blue"); ax.text(-1400, -800, "PKP", color="tab:red"); ax.text(-600, -100, "PKIKP", color="tab:green")
ax.set_aspect("equal"); ax.axis("off")
ax.set_xlim(-1.15 * sphray.R_TOP, 1.75 * sphray.R_TOP); ax.set_ylim(-1.2 * sphray.R_TOP, 1.2 * sphray.R_TOP)
fig.tight_layout()
fig.savefig(str(FIG / "L6F3.pdf"))
