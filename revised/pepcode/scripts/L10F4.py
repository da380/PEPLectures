"""Lecture 22, Fig. 4: the dwarf planet Haumea (semi-axes from the 2017 stellar occultation,
Ortiz et al. 2017) compared with the Jacobi ellipsoid of the same volume and the same ratio
of equatorial axes: equatorial and meridional sections drawn to scale."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pepseis.ellipsoids import jacobi_axes, omega2
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
haumea = np.array([2322.0, 1704.0, 1138.0]) / 2.0                     # semi-axes, km
abar = haumea.prod() ** (1.0 / 3.0)
jac = jacobi_axes(haumea[1] / haumea[0]) * abar
fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.6), constrained_layout=True)
for ax, (i, j, name) in zip(axes, [(0, 1, "equatorial section"), (0, 2, "meridional section")]):
    ax.add_patch(Ellipse((0, 0), 2 * haumea[i], 2 * haumea[j], facecolor="0.88", edgecolor="k", lw=1.4, label="Haumea"))
    ax.add_patch(Ellipse((0, 0), 2 * jac[i], 2 * jac[j], facecolor="none", edgecolor="tab:red", lw=1.4, ls="--", label="Jacobi ellipsoid"))
    ax.set_xlim(-1350, 1350); ax.set_ylim(-1350, 1350); ax.set_aspect("equal")
    ax.set_xlabel(r"$a_{1}$ direction / km"); ax.set_ylabel((r"$a_{2}$" if j == 1 else r"$a_{3}$ (rotation axis)") + " direction / km")
    ax.set_title(name, fontsize=10); ax.grid(True, ls=":", lw=0.5)
axes[0].legend(loc="upper right", fontsize=8, frameon=False)
txt = ("Haumea: $2322\\times1704\\times1138$ km\nJacobi, same $a_2/a_1$: $%d\\times%d\\times%d$ km" % tuple(np.round(2 * jac)))
axes[1].text(0.03, 0.03, txt, transform=axes[1].transAxes, fontsize=8, va="bottom")
fig.savefig(str(FIG / "L10F4.pdf")); fig.savefig(str(FIG / "L10F4.png"), dpi=200)
