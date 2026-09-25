"""Lecture 23, Fig. 2: amplitude spectra of the II.BFO record of Fig. 1 (vertical component, ten
days; transverse component, three days, as the horizontal noise is larger at low frequency) with the fundamental spheroidal and toroidal mode frequencies of PREM
(mode_lab_2 catalogue) marked."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.seismograms import tohoku, amplitude_spectrum
from pepseis.modes import catalogue
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
st = tohoku("II_BFO", "ZT")
fmin, fmax = 0.25, 3.0
ls, ns, fs = catalogue("spheroidal"); lt, nt, ft = catalogue("toroidal")
fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.6), constrained_layout=True, sharex=True)
for ax, chan, name, (lc, nc, fc), sym, col in [(axes[0], "Z", "vertical", (ls, ns, fs), "S", "tab:red"),
                                               (axes[1], "T", "transverse", (lt, nt, ft), "T", "tab:blue")]:
    tr = st.select(channel="LH" + chan)[0]
    f, sp = amplitude_spectrum(tr, 0.0, 240.0 if chan == "Z" else 72.0)
    band = (f >= fmin) & (f <= fmax)
    sp = sp[band] / sp[band].max()
    ax.plot(f[band], sp, "k", lw=0.6)
    ax.set_xlim(fmin, fmax); ax.set_ylim(0, 1.3); ax.grid(True, ls=":", lw=0.5)
    ax.set_ylabel("amplitude, " + name)
    for L in range(0 if sym == "S" else 2, 60):
        sel = (lc == L) & (nc == 0)
        if not sel.any():
            continue
        fm = 1e3 * fc[sel][0]
        if fmin < fm < fmax:
            ax.axvline(fm, color=col, lw=0.5, alpha=0.5, zorder=0)
            if L <= 20:
                ax.text(fm, 1.05 if L == 0 else 1.28, r"${}_0%s_{%d}$" % (sym, L), rotation=90, fontsize=9,
                        ha="center", va="top", color=col)
axes[1].set_xlabel("frequency / mHz")
fig.savefig(str(FIG / "L11F2.pdf")); fig.savefig(str(FIG / "L11F2.png"), dpi=200)
