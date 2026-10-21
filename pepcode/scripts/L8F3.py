"""Lecture 20, Fig. 3: construction of a cross-correlation delay-time kernel by the adjoint
method in two dimensions.  Rows: three instants; columns: the forward field, the adjoint
field, their product (the integrand of the kernel), and the kernel accumulated up to that instant.
With --webm, also a movie of the four columns through time, for the slides: ../figures/anim/L8F3.webm."""
import sys
import numpy as np
from pepseis.kernel2d import KernelExperiment
from pepseis.paths import FIG
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 14})

E = KernelExperiment(); E.forward(); E.adjoint(); K_final = E.kernel()
times = [0.65, 0.95, 1.25]
idx = [int(np.argmin(np.abs(E.snap_t - tt))) for tt in times]
fig, axes = plt.subplots(3, 4, figsize=(12.0, 4.9), constrained_layout=True)
vf = 0.4 * max(np.abs(E.snap_div[i]).max() for i in idx)
va = 0.4 * max(np.abs(E.snap_div_adj[i]).max() for i in idx)
vk = 4.0                                   # the kernel is singular at the source and receiver; clip there
vi = 0.4 * max(np.abs(E.integrand(i)).max() for i in idx)
for r, i in enumerate(idx):
    fields = [(E.snap_div[i], vf, r"forward field $\nabla\cdot\mathbf{u}(\mathbf{x},t)$"),
              (E.snap_div_adj[i], va, r"adjoint field $\nabla\cdot\mathbf{u}^{\dagger}(\mathbf{x},T-t)$"),
              (E.integrand(i), vi, r"integrand $-2\rho\alpha\,\nabla\cdot\mathbf{u}\,\nabla\cdot\mathbf{u}^{\dagger}$"),
              (E.kernel(upto=i + 1), vk, r"kernel $K_{\alpha}$ accumulated up to $t$")]
    for c, (f, v, title) in enumerate(fields):
        ax = axes[r, c]
        ax.pcolormesh(E.x[E.plot_x], E.z[E.plot_z], E.trim(f).T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto", rasterized=True)
        ax.plot(*E.src, "k*", ms=9); ax.plot(*E.rec, "kv", ms=7)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        if r == 0:
            ax.set_title(title, fontsize=13)
    axes[r, 0].set_ylabel(r"$t=%.2f$" % E.snap_t[i], fontsize=13)
fig.savefig(str(FIG / "L8F3.pdf"), dpi=200)

if "--webm" in sys.argv:
    import imageio_ffmpeg
    from matplotlib.animation import FFMpegWriter
    plt.rcParams.update({"font.size": 11, "animation.ffmpeg_path": imageio_ffmpeg.get_ffmpeg_exe()})
    fig, axes = plt.subplots(1, 4, figsize=(12.0, 2.6), constrained_layout=True)
    titles = [r"forward $\nabla\cdot\mathbf{u}(\mathbf{x},t)$", r"adjoint $\nabla\cdot\mathbf{u}^{\dagger}(\mathbf{x},T-t)$",
              r"integrand", r"kernel $K_{\alpha}$ accumulated up to $t$"]
    K = np.zeros((len(E.x), len(E.z)))
    writer = FFMpegWriter(fps=15, codec="libvpx-vp9", extra_args=["-crf", "34", "-b:v", "0", "-pix_fmt", "yuv420p"])
    with writer.saving(fig, str(FIG / "anim" / "L8F3.webm"), 100):
        for i in range(len(E.snap_t)):
            if E.snap_div_adj[i] is not None:
                K += E.snap_div[i] * E.snap_div_adj[i]
            if i % 2 or not 0.1 < E.snap_t[i] < E.t_end - 0.05:
                continue
            kernel = -2.0 * E.rho * E.alpha * K * E.dt * E.snap_every
            fields = [(E.snap_div[i], vf), (E.snap_div_adj[i], va), (E.integrand(i), vi), (kernel, vk)]
            for ax, (fld, v), title in zip(axes, fields, titles):
                ax.clear()
                ax.pcolormesh(E.x[E.plot_x], E.z[E.plot_z], E.trim(fld).T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto")
                ax.plot(*E.src, "k*", ms=8); ax.plot(*E.rec, "kv", ms=6)
                ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=11)
            axes[0].set_ylabel(r"$t=%.2f$" % E.snap_t[i])
            writer.grab_frame()
