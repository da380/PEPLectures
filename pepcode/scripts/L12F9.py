"""Lecture 24: splitting of 0S2 and 0S3 seen in twenty-day vertical records of the 2011 Tohoku
earthquake (EarthScope data, scripts/fetch_splitting_data.py).  Hann (cosine-bell) window from
three hours to twenty days after the origin time, zero-padded spectra; the PREM degenerate
frequencies (corrected for physical dispersion, see pepseis.modes.modal_q_and_dispersion) are marked."""
import os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from obspy import read
from pepseis.seismograms import DATA, TOHOKU
from pepseis.modes import catalogue, PremModes, modal_q_and_dispersion
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
# PREM degenerate frequencies corrected for physical dispersion (the catalogue is elastic at 1 s)
pm = PremModes(order=5, element_size_km=50.0)
f02 = 1e3 * modal_q_and_dispersion(pm, 2, 0)[2]; f03 = 1e3 * modal_q_and_dispersion(pm, 3, 0)[2]
stations = [("IU", "CTAO"), ("II", "NNA"), ("II", "BFO")]
fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6), constrained_layout=True)
t0 = TOHOKU["time"]
for k, (net, sta) in enumerate(stations):
    tr = read(os.path.join(DATA, "tohoku20_%s_%s.mseed" % (net, sta)))[0]
    tr.trim(t0 + 3 * 3600, t0 + 20 * 86400)
    x = tr.data - tr.data.mean(); N = len(x); dt = tr.stats.delta
    X = np.abs(np.fft.rfft(x * np.hanning(N), n=8 * N)); fr = np.fft.rfftfreq(8 * N, dt) * 1e3
    for ax, (a, b) in zip(axes, [(0.293, 0.325), (0.455, 0.481)]):
        sel = (fr > a) & (fr < b); y = X[sel] / X[sel].max()
        ax.plot(fr[sel], y + 1.1 * (2 - k), "k", lw=0.9)
        ax.text(a + 0.0005, 1.1 * (2 - k) + 0.9, sta, fontsize=9, va="top")
for ax, fm, name in [(axes[0], f02, r"${}_0S_2$"), (axes[1], f03, r"${}_0S_3$")]:
    ax.axvline(fm, color="tab:red", lw=0.8, ls="--")
    ax.text(fm, 3.45, name + " (PREM)", color="tab:red", ha="center", va="bottom", fontsize=9)
    ax.set_ylim(0, 3.6); ax.set_yticks([]); ax.set_xlabel("frequency / mHz"); ax.grid(True, axis="x", ls=":", lw=0.5)
axes[0].set_xlim(0.293, 0.325); axes[1].set_xlim(0.455, 0.481)
axes[0].set_ylabel("normalised amplitude (offset by station)")
fig.savefig(str(FIG / "L12F9.pdf"))
print("PREM 0S2 %.4f mHz, 0S3 %.4f mHz" % (f02, f03))
