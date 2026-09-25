"""Lecture 22, Fig. 2: PREM density and the reciprocal of the hydrostatic ellipticity from
Clairaut's equation."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem
from pepseis.clairaut import ellipticity
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
r, eps = ellipticity()
rr, dd = [], []
for a, b in zip(prem.BOUNDARIES[:-1], prem.BOUNDARIES[1:]):
    x = np.linspace(a, b, 200); rr += list(x); dd += list(prem.density(np.clip(x, a + 1e-6, b - 1e-6)))
fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.2), sharey=True, constrained_layout=True)
axes[0].plot(dd, rr, "k", lw=1.3); axes[0].set_xlabel(r"density / kg m$^{-3}$"); axes[0].set_ylabel("radius / km"); axes[0].set_xlim(0, 14000); axes[0].set_xticks([0, 4000, 8000, 12000])
axes[1].plot(1 / eps, r, "k", lw=1.3); axes[1].set_xlabel(r"$1/\epsilon$"); axes[1].set_xlim(280, 430)
for ax in axes:
    for rb in (prem.R_CMB, prem.R_ICB):
        ax.axhline(rb, color="0.6", lw=0.6, ls="--")
    ax.set_ylim(0, 6400); ax.grid(True, ls=":", lw=0.5)
axes[0].text(13600, prem.R_CMB + 60, "CMB", ha="right", fontsize=9, color="0.4"); axes[0].text(13600, prem.R_ICB + 60, "ICB", ha="right", fontsize=9, color="0.4")
fig.savefig(str(FIG / "L10F2.pdf"))
print("1/eps surface %.1f, centre %.1f" % (1 / eps[-1], 1 / eps[0]))
