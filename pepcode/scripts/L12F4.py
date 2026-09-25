"""Lecture 24, Fig. 4: density against radius in PREM."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
rr, dd = [], []
for a, b in zip(prem.BOUNDARIES[:-1], prem.BOUNDARIES[1:]):
    r = np.linspace(a, b, 200); rr += list(r); dd += list(prem.density(np.clip(r, a + 1e-6, b - 1e-6)) / 1e3)
fig, ax = plt.subplots(figsize=(6.0, 4.4), constrained_layout=True)
ax.plot(dd, rr, "k", lw=1.3)
for rb, lab in [(prem.R_CMB, "CMB"), (prem.R_ICB, "ICB")]:
    ax.axhline(rb, color="0.6", lw=0.6, ls="--"); ax.text(13.6, rb + 60, lab, ha="right", fontsize=9, color="0.4")
ax.set_xlim(0, 14); ax.set_ylim(0, 6400); ax.set_xlabel(r"density / g cm$^{-3}$"); ax.set_ylabel("radius / km"); ax.grid(True, ls=":", lw=0.5)
fig.savefig(str(FIG / "L12F4.pdf"))
