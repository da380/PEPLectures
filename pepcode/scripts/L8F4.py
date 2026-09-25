"""Lecture 20, Fig. 4: the finished 2D delay-time kernel for the P-wave speed, with a
cross-section perpendicular to the ray (there is no 'doughnut hole' in two dimensions)."""
import numpy as np
from pepseis.kernel2d import KernelExperiment
from pepseis.paths import FIG
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})

E = KernelExperiment(); E.forward(); E.adjoint(); K = E.kernel()
fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.6), constrained_layout=True, gridspec_kw={"width_ratios": [2.2, 1]})
ax = axes[0]; v = 4.0                       # clipped: the kernel is singular at the source and receiver
im = ax.pcolormesh(E.x[E.plot_x], E.z[E.plot_z], E.trim(K).T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto", rasterized=True)
ax.plot([E.src[0], E.rec[0]], [E.src[1], E.rec[1]], "k--", lw=0.8)
ax.plot(*E.src, "k*", ms=10); ax.plot(*E.rec, "kv", ms=8)
xm = 0.5 * (E.src[0] + E.rec[0]); ax.axvline(xm, color="0.3", lw=0.6, ls=":")
ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(r"kernel $K_{\alpha}$", fontsize=10)
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02); cb.set_label(r"$K_{\alpha}$")
ax = axes[1]
ix = int(np.argmin(np.abs(E.x - xm)))
ax.plot(K[ix, E.plot_z], E.z[E.plot_z], "k", lw=1.2)
ax.axvline(0, color="0.6", lw=0.6); ax.axhline(E.src[1], color="0.6", lw=0.6, ls="--")
ax.set_xlabel(r"$K_{\alpha}$"); ax.set_ylabel(r"$z$"); ax.set_title("cross-section at the mid-point", fontsize=10)
ax.set_ylim(0, 1)
fig.savefig(str(FIG / "L8F4.pdf"), dpi=200)
print("kernel on the ray at mid-point: %.3g (min %.3g)" % (K[ix, int(np.argmin(np.abs(E.z - E.src[1])))], K.min()))
