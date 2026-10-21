"""Slides, Lecture 20: a Taylor test of the 2D delay-time kernel of L8F3/L8F4.  A fixed perturbation of the P-wave
speed (a Gaussian bump off the ray, within the first Fresnel zone) is scaled by amplitudes s; the change in the
cross-correlation delay time found by rerunning the forward problem is compared with the first-order prediction
s * integral(K_alpha dalpha) from the one adjoint calculation.  They agree to first order in s.
Written as PDF to ../figures/slides.
Usage:  poetry run python scripts/slides_taylor.py"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.kernel2d import KernelExperiment
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 13})
CENTRE, WIDTH = (1.0, 0.6), 0.06
AMPLITUDES = [-0.2, -0.15, -0.1, -0.05, 0.05, 0.1, 0.15, 0.2]

def shape(X, Z):
    """The perturbation, dalpha / alpha, at unit amplitude."""
    return np.exp(-((X - CENTRE[0]) ** 2 + (Z - CENTRE[1]) ** 2) / (2 * WIDTH ** 2))

E = KernelExperiment(); E.forward(); E.adjoint(); K = E.kernel()
X, Z = np.meshgrid(E.x, E.z, indexing="ij")
slope = np.sum(K * E.alpha0 * shape(X, Z)) * E.h ** 2           # d tau_bar / d s
exact = []
for s in AMPLITUDES:
    Ep = KernelExperiment(alpha_pert=lambda X, Z, s=s: s * shape(X, Z)); Ep.forward()
    exact.append(Ep.delay_time(E.s_rec))       # perturbed synthetic against the reference, as in check_kernel.py
    print(f"s = {s:+.2f}: exact {exact[-1]:+.5f}, first order {s * slope:+.5f}")

fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.25, 1]})
a0.pcolormesh(E.x[E.plot_x], E.z[E.plot_z], E.trim(shape(X, Z)).T, cmap="Reds", vmin=0, vmax=1, shading="auto",
              rasterized=True)
v = 4.0                                        # the kernel's first Fresnel zone, for orientation (as in L8F4)
a0.contour(E.x[E.plot_x], E.z[E.plot_z], E.trim(K).T, levels=[-0.25 * v], colors="0.4", linewidths=0.8,
           linestyles="--")
a0.plot(*E.src, "k*", ms=11); a0.plot(*E.rec, "kv", ms=9)
a0.set_aspect("equal"); a0.set_xticks([]); a0.set_yticks([])
a0.set_title(r"perturbation $\delta\alpha$", fontsize=13)
s_line = np.linspace(-0.22, 0.22, 2)
a1.plot(s_line, slope * s_line, color="tab:blue", lw=1.6, label=r"first order: $s\int K_{\alpha}\,\delta\alpha\,\mathrm{d}^2\mathbf{x}$")
a1.plot(AMPLITUDES, exact, "o", color="k", ms=6, label="exact: forward runs")
a1.axhline(0, color="0.7", lw=0.6); a1.axvline(0, color="0.7", lw=0.6)
a1.set_xlabel(r"amplitude $s$ of $\delta\alpha/\alpha$"); a1.set_ylabel(r"change in delay time $\bar{\tau}$")
a1.legend(loc="best", frameon=False, fontsize=12)
a1.grid(True, ls=":", lw=0.5)
fig.tight_layout()
fig.savefig(str(FIG / "slides" / "taylor_test.pdf"))
