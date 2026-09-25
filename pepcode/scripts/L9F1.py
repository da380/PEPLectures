"""Lecture 21, Fig. 1: gravitational acceleration and pressure in PREM as functions of radius,
computed from the density (Fig. 2 of Lecture 24 shows the density itself)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})

r = np.linspace(1.0, prem.R_EARTH, 400)
g = np.array([prem.gravity(ri) for ri in r])
p = np.array([prem.pressure(ri) for ri in r]) / 1e9
fig, axes = plt.subplots(1, 2, figsize=(8.0, 4.2), sharey=True, constrained_layout=True)
axes[0].plot(g, r, "k", lw=1.3); axes[0].set_xlabel(r"gravitational acceleration / m s$^{-2}$"); axes[0].set_ylabel("radius / km")
axes[1].plot(p, r, "k", lw=1.3); axes[1].set_xlabel("pressure / GPa")
for ax in axes:
    for rb, lab in [(prem.R_CMB, "CMB"), (prem.R_ICB, "ICB")]:
        ax.axhline(rb, color="0.6", lw=0.6, ls="--")
    ax.set_ylim(0, 6400); ax.grid(True, ls=":", lw=0.5)
axes[0].text(0.3, prem.R_CMB + 70, "CMB", fontsize=9, color="0.4"); axes[0].text(0.3, prem.R_ICB + 70, "ICB", fontsize=9, color="0.4")
axes[0].set_xlim(0, 11.5); axes[1].set_xlim(0, 380)
fig.savefig(str(FIG / "L9F1.pdf"))
print("g surface %.2f, g CMB %.2f, p CMB %.0f GPa, p centre %.0f GPa" % (g[-1], prem.gravity(prem.R_CMB), prem.pressure(prem.R_CMB) / 1e9, p[0]))
