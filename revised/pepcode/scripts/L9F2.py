"""Lecture 21, Fig. 2: synthetic vertical-component acceleration seismograms in PREM computed
with and without self-gravitation (yspec; see run_yspec.py), for an Mw 8 source at 20 km depth
recorded at an epicentral angle of 90 degrees.  Top: the full band of the calculation
(0.3-50 mHz); below: the same seismogram band-passed in successively lower frequency bands.  A cosine taper at the start of each trace, two periods
of its lowest frequency long, removes a transient at the origin time; the record is six hours long, without attenuation."""
import os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "yspec")
k = 2                                                   # receiver 2: Delta = 90 degrees
t, z0 = np.loadtxt(os.path.join(data, f"q0.{k}"), usecols=(0, 1), unpack=True)
_, z2 = np.loadtxt(os.path.join(data, f"q2.{k}"), usecols=(0, 1), unpack=True)
dt = t[1] - t[0]


def bandpass(x, dt, f1, f2, taper=0.15):
    """Zero-phase band-pass in the frequency domain with cosine tapers of relative width
    `taper` at each edge (no start-up transient, unlike a recursive filter)."""
    n = len(x)
    # taper the final twenty minutes (the record does not end at rest) and zero-pad to avoid
    # circular wraparound of the end of the record onto its start
    x = x.copy(); m = int(1200 / dt); x[-m:] *= 0.5 * (1 + np.cos(np.pi * np.arange(m) / m))
    x = np.concatenate([x, np.zeros(n)])
    X = np.fft.rfft(x); f = np.fft.rfftfreq(2 * n, dt)
    lo, hi = f1 * (1 - taper), f2 * (1 + taper)
    w = np.zeros_like(f)
    w[(f >= f1) & (f <= f2)] = 1.0
    m = (f >= lo) & (f < f1); w[m] = 0.5 * (1 - np.cos(np.pi * (f[m] - lo) / (f1 - lo)))
    m = (f > f2) & (f <= hi); w[m] = 0.5 * (1 + np.cos(np.pi * (f[m] - f2) / (hi - f2)))
    return np.fft.irfft(X * w, 2 * n)[:n]


bands = [None, (20e-3, 45e-3), (10e-3, 20e-3), (5e-3, 10e-3), (2e-3, 5e-3), (1e-3, 2e-3)]
labels = ["0.3–50 mHz", "20–45 mHz", "10–20 mHz", "5–10 mHz", "2–5 mHz", "1–2 mHz"]
fig, axes = plt.subplots(len(bands), 1, figsize=(9.6, 9.6), sharex=True, constrained_layout=True)
for ax, band, lab in zip(axes, bands, labels):
    if band is None:
        y0, y2 = z0, z2
    else:
        y0, y2 = bandpass(z0, dt, *band), bandpass(z2, dt, *band)
    # a cosine taper at the start, two periods of the lowest frequency in the band long (two
    # minutes for the full band), removes a transient at the origin time that dominates the
    # lowest bands
    T = 120.0 if band is None else min(2.0 / band[0], 1800.0)
    ramp = np.ones_like(t); m = t < T; ramp[m] = 0.5 * (1 - np.cos(np.pi * t[m] / T))
    y0, y2 = y0 * ramp, y2 * ramp
    scale = np.abs(y2).max()
    ax.plot(t / 60, y2 / scale, "k", lw=0.7, label="with self-gravitation")
    ax.plot(t / 60, y0 / scale, "tab:red", lw=0.7, alpha=0.85, label="without gravitation")
    ax.set_ylim(-1.15, 1.15); ax.set_yticks([-1, 0, 1]); ax.grid(True, ls=":", lw=0.5)
    ax.text(0.01, 0.9, lab, transform=ax.transAxes, fontsize=9, va="top")
axes[0].legend(loc="upper right", fontsize=8, frameon=False)
axes[-1].set_xlabel("time after the earthquake / minutes"); axes[-1].set_xlim(0, 360)
fig.savefig(str(FIG / "L9F2.pdf")); fig.savefig(str(FIG / "L9F2.png"), dpi=200)
