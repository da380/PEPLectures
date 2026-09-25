"""Figure for Example 3 of examples8.tex: steepest descent with exact line search on
J(m) = (1/2) m^T H m with H = diag(1, 9), starting from m_0 = (9, 1).

Saves examples/figures_ex/ex8_steepest.png.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
outdir = os.path.join(here, "..", "examples", "figures_ex")
os.makedirs(outdir, exist_ok=True)

H = np.diag([1.0, 9.0])
m = np.array([9.0, 1.0])
path = [m.copy()]
for _ in range(12):
    g = H @ m
    lam = (g @ g) / (g @ H @ g)
    m = m - lam * g
    path.append(m.copy())
path = np.array(path)

x = np.linspace(-2, 10, 400)
y = np.linspace(-2.5, 2.5, 400)
X, Y = np.meshgrid(x, y)
J = 0.5 * (X ** 2 + 9 * Y ** 2)

fig, ax = plt.subplots(figsize=(7.5, 3.6), facecolor="white")
levels = 45 * 0.64 ** np.arange(0, 12)
ax.contour(X, Y, J, levels=np.sort(levels), colors="0.6", linewidths=0.8)
ax.plot(path[:, 0], path[:, 1], "-o", color="C3", ms=3.5, lw=1.2, label="steepest descent iterates")
ax.plot(0, 0, "k*", ms=9, label=r"minimum $\overline{\mathbf{m}}$")
ax.set_xlabel(r"$m_1$ (dimensionless)")
ax.set_ylabel(r"$m_2$ (dimensionless)")
ax.set_aspect("equal")
ax.set_xlim(-2, 10)
ax.set_ylim(-2.5, 2.5)
ax.legend(loc="upper right", frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(outdir, "ex8_steepest.png"), dpi=200, facecolor="white")
print("saved", os.path.join(outdir, "ex8_steepest.png"))
