"""Figure for Example 4 of examples11.tex.

Amplitude spectrum of the sum of two equal-amplitude cosines, at the
frequencies of two adjacent singlets of 0S2, observed through a boxcar
window of length T = 24 h and T = 120 h.  Saved to
examples/figures_ex/ex11_window.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "examples", "figures_ex", "ex11_window.png")

f0 = 0.3094e-3          # Hz, degenerate frequency of 0S2
df = 4.6e-6             # Hz, spacing of adjacent singlets
f1, f2 = f0 - df / 2, f0 + df / 2
records = [(24.0 * 3600.0, "T = 24 h", "#2a78d6", "-"),
           (120.0 * 3600.0, "T = 120 h", "#eb6834", "--")]


def htilde(omega, T):
    """Transform of the boxcar on [0,T] with the lecture's convention."""
    x = omega * T / 2.0
    return T * np.exp(-1j * x) * np.sinc(x / np.pi)


def spectrum(f, T):
    w = 2 * np.pi * f
    tot = 0.0
    for fk in (f1, f2):
        wk = 2 * np.pi * fk
        tot = tot + 0.5 * (htilde(w - wk, T) + htilde(w + wk, T))
    return np.abs(tot)


f = np.linspace(0.290e-3, 0.330e-3, 4001)
fig, ax = plt.subplots(figsize=(6.4, 3.6), facecolor="white")
for T, lab, col, ls in records:
    s = spectrum(f, T)
    ax.plot(f * 1e3, s / s.max(), color=col, ls=ls, lw=1.8, label=lab)
for fk in (f1, f2):
    ax.axvline(fk * 1e3, color="0.6", lw=0.8, ls=":")
ax.set_xlabel("Frequency (mHz)")
ax.set_ylabel("Normalised amplitude")
ax.set_xlim(0.290, 0.330)
ax.set_ylim(0, 1.05)
ax.grid(True, color="0.9", lw=0.6)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(OUT, dpi=200, facecolor="white")
print("saved", os.path.abspath(OUT))
