"""Lecture 22, Fig. 3: the Maclaurin and Jacobi sequences of equilibrium figures of a
homogeneous rotating fluid, as rotation rate against angular momentum (both nondimensional),
with the shapes of a few members drawn (equatorial section above, meridional section below)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pepseis.ellipsoids import maclaurin_axes, jacobi_axes, omega2, angular_momentum
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 9})

es = np.concatenate([np.linspace(0.02, 0.8, 60), np.linspace(0.8, 0.995, 80)])
mac = np.array([[angular_momentum(maclaurin_axes(e)), omega2(maclaurin_axes(e))] for e in es])
ratios = np.concatenate([np.linspace(0.9999, 0.9, 12), np.linspace(0.9, 0.2, 40)])
jac = np.array([[angular_momentum(jacobi_axes(r)), omega2(jacobi_axes(r))] for r in ratios])
bif = jacobi_axes(0.9999); pear = jacobi_axes(0.43218)

# drawn at its printed width (\textwidth) so that the fonts appear at their nominal size
fig = plt.figure(figsize=(6.5, 4.4), constrained_layout=True)
gs = fig.add_gridspec(2, 7, height_ratios=[3.0, 1.5])
ax = fig.add_subplot(gs[0, :])
ax.plot(mac[:, 0], mac[:, 1], "k", lw=1.4, label="Maclaurin spheroids")
ax.plot(jac[:, 0], jac[:, 1], "tab:red", lw=1.4, label="Jacobi ellipsoids")
ax.plot(angular_momentum(bif), omega2(bif), "ko", ms=6); ax.annotate("bifurcation", (angular_momentum(bif), omega2(bif)), xytext=(-58, 8), textcoords="offset points", fontsize=8)
Lp, Op = angular_momentum(pear), omega2(pear)
# the pear-shaped sequence is known only approximately: drawn schematically, leaving the
# Jacobi sequence tangentially and continuing towards larger angular momentum
s_ = np.linspace(0, 1, 50)
jac_slope = np.gradient(jac[:, 1], jac[:, 0])[np.argmin(np.abs(jac[:, 0] - Lp))]
Lpear = Lp + 0.14 * s_; Opear = Op + jac_slope * (Lpear - Lp) - 0.06 * s_ ** 2
ax.plot(Lpear, Opear, "--", color="tab:blue", lw=1.4, label="pear-shaped figures (schematic)")
ax.plot(Lp, Op, "s", color="tab:red", ms=6); ax.annotate("pear-shaped figures bifurcate", (Lp, Op), xytext=(6, 6), textcoords="offset points", fontsize=8)
imax = np.argmax(mac[:, 1]); ax.annotate("maximum rotation rate", (mac[imax, 0], mac[imax, 1]), xytext=(8, 3), textcoords="offset points", fontsize=8)
ax.set_xlim(0, 0.75); ax.set_ylim(0, 0.5)
ax.set_xlabel(r"angular momentum $L/\sqrt{GM^{3}\bar{a}}$"); ax.set_ylabel(r"$\Omega^{2}/\pi G\rho$")
ax.legend(loc="lower right", fontsize=8, frameon=False); ax.grid(True, ls=":", lw=0.5)
# shapes
members = [("sphere", np.array([1.0, 1.0, 1.0])), (r"Maclaurin," "\n" r"$e=0.6$", maclaurin_axes(0.6)), ("Maclaurin at\nbifurcation", bif),
           (r"Jacobi," "\n" r"$a_2/a_1=0.6$", jacobi_axes(0.6)), ("Jacobi at\nthe pear point", pear)]
for k, (name, a) in enumerate(members):
    axs = fig.add_subplot(gs[1, k])
    axs.add_patch(Ellipse((0, 0.9), 2 * a[0], 2 * a[1], facecolor="0.85", edgecolor="k", lw=1))   # equatorial section
    axs.add_patch(Ellipse((0, -1.2), 2 * a[0], 2 * a[2], facecolor="0.85", edgecolor="k", lw=1))  # meridional section
    axs.set_xlim(-2.1, 2.1); axs.set_ylim(-2.3, 2.3); axs.set_aspect("equal"); axs.axis("off")
    axs.set_title(name, fontsize=7.5)
# schematic pear-shaped and fissioning figures (no closed form exists; these are sketches)
def outline(ax, centre, a, b, f, ls="-"):
    x = np.linspace(-a, a, 400); half = b * np.sqrt(np.clip(1 - (x / a) ** 2, 0, None)) * f(x / a)
    ax.fill_between(x + centre[0], centre[1] - half, centre[1] + half, facecolor="0.85", edgecolor="k", lw=1, ls=ls)
for k, (name, f, ls) in enumerate([("pear-shaped\n(schematic)", lambda s: (1 + 0.25 * s) * (1 - 0.15 * np.exp(-((s + 0.3) / 0.35) ** 2)), "-"),
                                    ("fission?\n(schematic)", lambda s: (1 + 0.1 * s) * (1 - 0.6 * np.exp(-(s / 0.3) ** 2)), "--")]):
    axs = fig.add_subplot(gs[1, 5 + k]); a = pear
    outline(axs, (0, 0.9), a[0], a[1], f, ls); outline(axs, (0, -1.2), a[0], a[2], f, ls)
    axs.set_xlim(-2.1, 2.1); axs.set_ylim(-2.3, 2.3); axs.set_aspect("equal"); axs.axis("off"); axs.set_title(name, fontsize=7.5)
fig.savefig(str(FIG / "L10F3.pdf"))
print("max Omega^2 on Maclaurin sequence %.4f at e=%.4f" % (mac[imax, 1], es[imax]))
