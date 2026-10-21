"""Slides, Lecture 18: P- and S-wave speeds against radius in the Jeffreys-Bullen model of 1940 (dashed) and in
PREM (solid, as in L6F5).  The J-B speeds are those distributed with obspy's TauP (jb.nd); its inner-core S speed
was added later and is not drawn (Jeffreys and Bullen had none).  Written as PDF to ../figures/slides.
Usage:  poetry run python scripts/slides_jb_prem.py"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import obspy
from pepseis import prem
from pepseis.paths import FIG
from L6F5 import profile

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 13})

rows = []
for line in open(os.path.join(os.path.dirname(obspy.__file__), "taup", "data", "jb.nd")):
    parts = line.split()
    if len(parts) >= 3:
        rows.append([float(v) for v in parts[:3]])
depth, vp, vs = np.array(rows).T
r_jb = 6371.0 - depth
mantle = vs > 0
mantle &= r_jb > prem.R_CMB - 1

fig, ax = plt.subplots(figsize=(7.0, 5.0))
r, a = profile(prem.vp); ax.plot(a, r, "k", lw=1.4)
r, b = profile(prem.vs); ax.plot(b, r, "tab:red", lw=1.4)
ax.plot(vp, r_jb, "k--", lw=1.2)
ax.plot(vs[mantle], r_jb[mantle], "--", color="tab:red", lw=1.2)
for rb, lab in [(prem.R_CMB, "CMB"), (prem.R_ICB, "ICB")]:
    ax.axhline(rb, color="0.6", lw=0.6, ls=":"); ax.text(0.3, rb + 70, lab, ha="left", fontsize=11, color="0.4")
ax.plot([], [], "k", lw=1.4, label="PREM (1981)"); ax.plot([], [], "k--", lw=1.2, label="Jeffreys–Bullen (1940)")
ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.66), frameon=False, fontsize=12)
ax.text(12.9, 5000, r"$\alpha$", fontsize=15); ax.text(6.9, 4900, r"$\beta$", fontsize=15, color="tab:red")
ax.set_xlim(0, 14); ax.set_ylim(0, 6400)
ax.set_xlabel(r"wave speed / km s$^{-1}$"); ax.set_ylabel("radius / km")
ax.grid(True, ls=":", lw=0.5)
fig.tight_layout()
fig.savefig(str(FIG / "slides" / "jb_prem.pdf"))
