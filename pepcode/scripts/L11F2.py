"""Lecture 23, Fig. 2: amplitude spectra of the II.BFO record of Fig. 1 (vertical component, ten
days; transverse component, three days, as the horizontal noise is larger at low frequency) with the fundamental spheroidal and toroidal mode frequencies of PREM
(dispersion-corrected PREM values, see pepseis.modes) marked."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.seismograms import tohoku, amplitude_spectrum
from pepseis.modes import fundamental_frequencies_corrected
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
st = tohoku("II_BFO", "ZT")
fmin, fmax = 0.25, 3.0
# PREM fundamental-mode frequencies corrected for physical dispersion (the elastic catalogue is 0.3-0.5% high here)
lS, fS, lT, fT = fundamental_frequencies_corrected(20)
fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.6), constrained_layout=True, sharex=True)
for ax, chan, name, (lc, fc), sym, col in [(axes[0], "Z", "vertical", (lS, fS), "S", "tab:red"),
                                           (axes[1], "T", "transverse", (lT, fT), "T", "tab:blue")]:
    tr = st.select(channel="LH" + chan)[0]
    f, sp = amplitude_spectrum(tr, 0.0, 240.0 if chan == "Z" else 72.0)
    band = (f >= fmin) & (f <= fmax)
    sp = sp[band] / sp[band].max()
    ax.plot(f[band], sp, "k", lw=0.6)
    ax.set_xlim(fmin, fmax); ax.set_ylim(0, 1.45); ax.grid(True, ls=":", lw=0.5)
    ax.set_ylabel("amplitude, " + name)
    for L, fL in zip(lc, fc):
        fm = 1e3 * fL
        if fmin < fm < fmax:
            ax.axvline(fm, color=col, lw=0.5, alpha=0.5, zorder=0)
            if L <= 20:
                # the radial mode 0S0 lies just below 0S5, so its label is shifted to the left
                ax.text(fm - 0.035 if L == 0 else fm, 1.42, r"${}_0%s_{%d}$" % (sym, L), rotation=90, fontsize=9,
                        ha="center", va="top", color=col)
axes[1].set_xlabel("frequency / mHz")
fig.savefig(str(FIG / "L11F2.pdf"))
