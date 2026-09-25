"""Lecture 23, Fig. 1: ten days of the vertical-component acceleration recorded at the Black
Forest Observatory (II.BFO) after the 2011 Tohoku earthquake (data from EarthScope via obspy,
instrument response removed; see scripts/fetch_modes_data.py)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.seismograms import tohoku, hours
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
tr = tohoku("II_BFO", "Z")[0]
bp = tr.copy().filter("bandpass", freqmin=0.3e-3, freqmax=1e-3, corners=4, zerophase=True)
fig, axes = plt.subplots(2, 1, figsize=(8.0, 5.2), constrained_layout=True)
ax = axes[0]
ax.plot(hours(tr), 1e6 * tr.data, "k", lw=0.3)
ax.set_xlim(0, 12); ax.set_xlabel("time after the earthquake / hours")
ax.text(0.99, 0.94, "first twelve hours, unfiltered", transform=ax.transAxes, ha="right", va="top")
ax = axes[1]
d = hours(bp) / 24; sel = d >= 0.25
ax.plot(d[sel], 1e9 * bp.data[sel], "k", lw=0.3)
ax.set_xlim(0, 10); ax.set_xlabel("time after the earthquake / days")
ax.text(0.99, 0.94, "from six hours to ten days, band-pass filtered between 0.3 and 1 mHz", transform=ax.transAxes, ha="right", va="top")
axes[0].set_ylabel(r"acceleration / $\mu$m s$^{-2}$"); axes[1].set_ylabel(r"acceleration / nm s$^{-2}$")
for ax in axes:
    ax.grid(True, ls=":", lw=0.5)
fig.savefig(str(FIG / "L11F1.pdf")); fig.savefig(str(FIG / "L11F1.png"), dpi=200)
