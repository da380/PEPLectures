"""Fan of P-wave rays in a medium with a linear velocity gradient.

alpha(z) = alpha0 + g z, with z measured downwards. A ray leaving the origin at
take-off angle theta0 from the downward vertical is a circular arc of radius
R = 1/(g p_x), p_x = sin(theta0)/alpha0, centred at depth z = -alpha0/g.
The figure accompanies Example 1 of the worked examples for Lecture 16.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

alpha0 = 5.0     # km/s
g = 0.05         # 1/s   (velocity gradient, km/s per km)
angles_deg = np.arange(20.0, 90.0, 10.0)

fig, ax = plt.subplots(figsize=(7.0, 3.4), facecolor="white")
cmap = plt.get_cmap("Blues")
for k, th0_deg in enumerate(angles_deg):
    th0 = np.deg2rad(th0_deg)
    px = np.sin(th0) / alpha0
    R = 1.0 / (g * px)
    th = np.linspace(th0, np.pi - th0, 400)      # down to the return to the surface
    x = R * (np.cos(th0) - np.cos(th))
    z = R * (np.sin(th) - np.sin(th0))
    ax.plot(x, z, lw=1.4, color=cmap(0.45 + 0.55 * k / (len(angles_deg) - 1)),
            label=r"$\theta_0=%d^\circ$" % th0_deg)

ax.axhline(0.0, color="0.3", lw=0.8)
ax.set_xlim(0.0, 600.0)
ax.set_ylim(215.0, -8.0)          # depth increases downwards
ax.set_xlabel("Horizontal distance $x$ (km)")
ax.set_ylabel("Depth $z$ (km)")
ax.set_aspect("equal")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8, ncol=2, loc="lower right")
fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                   "examples", "figures_ex", "ex5_rays.png")
fig.savefig(out, dpi=200, facecolor="white")
print("wrote", os.path.abspath(out))
