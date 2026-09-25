"""Lecture 24, Fig. 6: radial eigenfunctions of some toroidal and spheroidal modes of PREM."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem
from pepseis.modes import PremModes
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
M = PremModes()
fig, axes = plt.subplots(2, 4, figsize=(10.5, 7.0), sharey=True, constrained_layout=True)
tor = [(2, 0), (2, 1), (10, 0), (40, 0)]
for ax, (l, n) in zip(axes[0], tor):
    f, r, W = M.toroidal(l, n)
    w = W[:, n] / np.abs(W[:, n]).max()
    ax.plot(w, r, "k", lw=1.3)
    ax.set_title(r"${}_{%d}T_{%d}$, %.3f mHz" % (n, l, 1e3 * f[n]), fontsize=10)
sph = [(2, 0), (2, 1), (10, 0), (0, 0)]
for ax, (l, n) in zip(axes[1], sph):
    f, r, U, V, P = M.spheroidal_mode(l, n)
    zeta = np.sqrt(l * (l + 1)) if l > 0 else 0.0
    scale = max(np.abs(U).max(), np.abs(zeta * V).max())
    ax.plot(U / scale, r, "k", lw=1.3, label=r"$U$")
    if l > 0:
        ax.plot(zeta * V / scale, r, "tab:red", lw=1.3, ls="--", label=r"$\zeta V$")
    ax.set_title(r"${}_{%d}S_{%d}$, %.3f mHz" % (n, l, 1e3 * f), fontsize=10)
axes[1, 0].legend(loc="lower right", fontsize=9, frameon=False)
for ax in axes.ravel():
    for rb in (prem.R_CMB, prem.R_ICB):
        ax.axhline(rb, color="0.6", lw=0.6, ls="--")
    ax.axvline(0, color="0.8", lw=0.6); ax.set_ylim(0, 6371); ax.set_xlim(-1.1, 1.1)
    ax.grid(True, ls=":", lw=0.4)
for ax in axes[:, 0]:
    ax.set_ylabel("radius / km")
for ax in axes[1]:
    ax.set_xlabel("normalised eigenfunction")
fig.savefig(str(FIG / "L12F6.pdf")); fig.savefig(str(FIG / "L12F6.png"), dpi=200)
