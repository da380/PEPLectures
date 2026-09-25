"""Lecture 18, Fig. 4: travel-time curves of the main body-wave phases in PREM (surface focus)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import sphray

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
P_PHASES = ["P", "PP", "PcP", "PKP", "PKiKP", "PKIKP"]
S_PHASES = ["S", "SS", "ScS", "SKS", "SKKS", "SKIKS"]
from pepseis import prem  # noqa: E402
from pepseis.paths import FIG
P_MOHO = {"P": 6151.0 / prem.vp(6151.0, "lower"), "S": 6151.0 / prem.vs(6151.0, "lower")}   # rays turning below 220 km depth
# label positions: (epicentral angle at which to label, x offset, y offset in points)
LABELS = {"P": (70, 8, -12), "PcP": (30, 2, 7), "PP": (130, -12, 8), "PKP": (160, 0, 9), "PKiKP": (70, 0, 7),
          "PKIKP": (140, 0, -13), "S": (75, 12, -12), "ScS": (18, -14, 6), "SS": (120, -14, 7), "SKS": (112, 10, -13),
          "SKKS": (150, 0, 8), "SKIKS": (165, 0, -13)}

fig, ax = plt.subplots(figsize=(6.4, 7.2))
for phases, col in [(P_PHASES, "tab:red"), (S_PHASES, "k")]:
    for ph in phases:
        # rays turning below 220 km only: above this r/v is almost constant in isotropic PREM and
        # the ducted rays produce long, sparse branches that only clutter the plot
        d, t, p = sphray.travel_time_curve(ph, n=1500, p_max=P_MOHO[sphray.PHASES[ph][0][0]])
        if d.size == 0:
            continue
        keep = d <= 180                           # discard the branches beyond 180 degrees
        d, t = d[keep], t[keep]
        # break the curve where a branch jump makes it discontinuous
        jumps = np.where((np.abs(np.diff(d)) > 3) | (np.abs(np.diff(t)) > 60))[0]
        for seg in np.split(np.arange(d.size), jumps + 1):
            ax.plot(d[seg], t[seg] / 60, color=col, lw=1.0)
        dl, dx, dy = LABELS[ph]
        i = np.argmin(np.abs(d - dl))
        ax.annotate(ph, (d[i], t[i] / 60), xytext=(dx, dy), textcoords="offset points", color=col, fontsize=10, ha="center")
ax.set_xlim(0, 180); ax.set_ylim(0, 40)
ax.set_xlabel(r"epicentral angle $\Delta$ / degrees"); ax.set_ylabel("travel time / minutes")
ax.grid(True, ls=":", lw=0.5)
fig.tight_layout()
fig.savefig(str(FIG / "L6F4.pdf")); fig.savefig(str(FIG / "L6F4.png"), dpi=200)
