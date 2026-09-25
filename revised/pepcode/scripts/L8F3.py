"""Lecture 20, Fig. 3: construction of a cross-correlation delay-time kernel by the adjoint
method in two dimensions.  Rows: three instants; columns: the forward field, the adjoint
field, and the kernel accumulated up to that instant."""
import numpy as np
from pepseis.kernel2d import KernelExperiment
from pepseis.paths import FIG
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})

E = KernelExperiment(); E.forward(); E.adjoint(); K_final = E.kernel()
times = [0.65, 0.95, 1.25]
idx = [int(np.argmin(np.abs(E.snap_t - tt))) for tt in times]
fig, axes = plt.subplots(3, 3, figsize=(10.5, 5.6), constrained_layout=True)
vf = 0.4 * max(np.abs(E.snap_div[i]).max() for i in idx)
va = 0.4 * max(np.abs(E.snap_div_adj[i]).max() for i in idx)
vk = 4.0                                   # the kernel is singular at the source and receiver; clip there
for r, i in enumerate(idx):
    fields = [(E.snap_div[i], vf, r"forward field $\nabla\cdot\mathbf{u}$"),
              (E.snap_div_adj[i], va, r"adjoint field $\nabla\cdot\mathbf{u}'$"),
              (E.kernel(upto=i + 1), vk, r"kernel $K_{\alpha}$ accumulated so far")]
    for c, (f, v, title) in enumerate(fields):
        ax = axes[r, c]
        ax.pcolormesh(E.x[E.plot_x], E.z[E.plot_z], E.trim(f).T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto", rasterized=True)
        ax.plot(*E.src, "k*", ms=9); ax.plot(*E.rec, "kv", ms=7)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if r == 0:
            ax.set_title(title, fontsize=10)
    axes[r, 0].set_ylabel(r"$t=%.2f$" % E.snap_t[i], fontsize=10)
fig.savefig(str(FIG / "L8F3.pdf"), dpi=200); fig.savefig(str(FIG / "L8F3.png"), dpi=200)
