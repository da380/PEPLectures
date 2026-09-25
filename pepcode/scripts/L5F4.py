"""Lecture 17, Fig. 4: rays from a plane wavefront passing through a smooth random medium
with 5 per cent rms velocity fluctuations.  Caustics form generically."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.raytrace2d import Medium, plane_wave_fan, random_medium
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

x = np.linspace(0, 10, 801)
z = np.linspace(0, 4, 321)
c = random_medium(x, z, c0=1.0, rms=0.05, corr_len=0.5, seed=7, x_onset=(0.6, 1.6))
med = Medium(x, z, c)
rays = plane_wave_fan(med, (0.05, 3.95), 0.0, 80, t_max=12.0, bounds=(0, 10, 0, 4))

fig, ax = plt.subplots(figsize=(7.2, 3.4))
im = ax.pcolormesh(x, z, c.T, cmap="RdBu", vmin=0.85, vmax=1.15, shading="auto", rasterized=True)
for t, rx, rz in rays:
    ax.plot(rx, rz, "k", lw=0.5)
ax.set_xlim(0, 10); ax.set_ylim(0, 4); ax.set_aspect("equal")
ax.set_xlabel(r"$x$"); ax.set_ylabel(r"$z$")
cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cb.set_label("wave speed")
fig.tight_layout()
fig.savefig(str(FIG / "L5F4.pdf"), dpi=200)
