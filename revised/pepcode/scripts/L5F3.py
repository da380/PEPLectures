"""Lecture 17, Fig. 3: rays from a plane wavefront passing through a smooth low-velocity
anomaly.  Behind the anomaly the rays cross and a caustic is formed."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.raytrace2d import Medium, plane_wave_fan, gaussian_lens
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

x = np.linspace(0, 5, 401)
z = np.linspace(0, 3, 241)
c = gaussian_lens(x, z, c0=1.0, dc=-0.22, centre=(1.8, 1.5), width=0.35)
med = Medium(x, z, c)
rays = plane_wave_fan(med, (0.05, 2.95), 0.0, 45, t_max=6.0, bounds=(0, 5, 0, 3))

fig, ax = plt.subplots(figsize=(7.2, 4.3))
im = ax.pcolormesh(x, z, c.T, cmap="RdBu", vmin=0.75, vmax=1.25, shading="auto", rasterized=True)
for t, rx, rz in rays:
    ax.plot(rx, rz, "k", lw=0.6)
ax.set_xlim(0, 5); ax.set_ylim(0, 3); ax.set_aspect("equal")
ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$z$")
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cb.set_label("wave speed")
fig.tight_layout()
fig.savefig(str(FIG / "L5F3.pdf"), dpi=200)
fig.savefig(str(FIG / "L5F3.png"), dpi=200)
