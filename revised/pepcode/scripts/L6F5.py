"""Lecture 18, Fig. 5: P- and S-wave speeds against radius in PREM."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis import prem
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})


def profile(func):
    """Radius and value arrays, layer by layer, so that a discontinuity is drawn as a
    horizontal step at the boundary radius."""
    rr, vv = [], []
    for a, b in zip(prem.BOUNDARIES[:-1], prem.BOUNDARIES[1:]):
        r = np.linspace(a, b, 300)
        v = func(np.clip(r, a + 1e-6, b - 1e-6))
        rr += list(r); vv += list(v)
    return np.array(rr), np.array(vv)


if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    r, a = profile(prem.vp); ax.plot(a, r, "k", lw=1.2, solid_capstyle="round")
    r, b = profile(prem.vs); ax.plot(b, r, "tab:red", lw=1.2, solid_capstyle="round")
    for rb, lab in [(prem.R_CMB, "CMB"), (prem.R_ICB, "ICB")]:
        ax.axhline(rb, color="0.6", lw=0.6, ls="--"); ax.text(13.9, rb + 60, lab, ha="right", fontsize=9, color="0.4")
    ax.text(12.2, 4700, r"$\alpha$", fontsize=13); ax.text(6.9, 4900, r"$\beta$", fontsize=13, color="tab:red")
    ax.set_xlim(0, 14); ax.set_ylim(0, 6400)
    ax.set_xlabel(r"wave speed / km s$^{-1}$"); ax.set_ylabel("radius / km")
    ax.grid(True, ls=":", lw=0.5)
    fig.tight_layout()
    fig.savefig(str(FIG / "L6F5.pdf")); fig.savefig(str(FIG / "L6F5.png"), dpi=200)
