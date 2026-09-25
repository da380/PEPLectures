"""Lecture 24, Figs. 1 and 2: toroidal and spheroidal eigenfrequencies of PREM against degree,
up to 50 mHz and degree 500.  Toroidal modes come from the reviewed mode_lab_2 catalogue (Myhill);
spheroidal modes from the reviewed catalogue for degrees up to 200 and from pepseis.modes.
spheroidal_extension (same solver and acceptance tests, cached in data/modes/) above that.
Radial modes (l = 0) are plotted with the spheroidal modes without distinction."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.modes import catalogue, spheroidal_extension
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})
FMAX, LMAX = 50.0, 500
for family, name in [("toroidal", "L12F1"), ("spheroidal", "L12F2")]:
    l, n, f = catalogue(family)
    if family == "spheroidal":
        le, fe = spheroidal_extension(201, LMAX, 1e-3 * FMAX)
        n0 = n[(l <= 200)]
        l, f = np.r_[l[l <= 200], le], np.r_[f[l <= 200], fe]
    keep = (l <= LMAX) & (f <= 1e-3 * FMAX)
    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    ax.scatter(l[keep], 1e3 * f[keep], s=1.2, c="k", lw=0, rasterized=True)
    ax.set_xlim(0, LMAX); ax.set_ylim(0, FMAX)
    ax.set_xlabel(r"degree $l$"); ax.set_ylabel("frequency / mHz")
    ax.grid(True, ls=":", lw=0.5)
    fig.savefig(str(FIG / f"{name}.pdf"), dpi=300)
    print(name, "modes plotted:", keep.sum())
