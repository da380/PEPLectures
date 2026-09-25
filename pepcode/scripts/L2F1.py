"""Lecture 14, Fig. 1: geometric meaning of the polar decomposition F = RU = VR.
A unit circle is mapped by U (stretch along the eigenvectors of U), by R (rotation)
and by F; the dashed lines follow the eigenvectors of U through each map."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 12})

# --- a concrete deformation gradient ---------------------------------------
def rot(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s], [s, c]])

theta_U = np.deg2rad(25)            # orientation of the eigenvectors of U
lam = np.array([1.7, 0.7])          # principal stretches
Q = rot(theta_U)
U = Q @ np.diag(lam) @ Q.T
R = rot(np.deg2rad(40))
F = R @ U
V = R @ U @ R.T                     # so that F = V R as well

# --- drawing helpers ---------------------------------------------------------
t = np.linspace(0, 2 * np.pi, 400)
circle = np.vstack([np.cos(t), np.sin(t)])
axes_pts = 1.35 * np.array([[1, -1, 0, 0], [0, 0, 1, -1]], float)   # dashed eigen-axes (in U's frame)
axes_pts = Q @ axes_pts

def draw(ax, A, centre, label=None):
    """Draw the image under A of the unit circle and of the dashed eigen-axes of U."""
    c = np.array(centre)[:, None]
    p = A @ circle + c
    ax.plot(p[0], p[1], "k", lw=1.4)
    a = A @ axes_pts + c
    ax.plot(a[0, :2], a[1, :2], "k--", lw=0.7)
    ax.plot(a[0, 2:], a[1, 2:], "k--", lw=0.7)

def arrow(ax, p, q, text, side):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=16, lw=1.4, color="k",
                                 shrinkA=0, shrinkB=0))
    m = (np.array(p) + np.array(q)) / 2 + np.array(side)
    ax.text(m[0], m[1], text, ha="center", va="center", fontsize=15)

fig, ax = plt.subplots(figsize=(7.2, 5.2))
I = np.eye(2)
cT, cL, cR, cB = (0, 0), (-4.6, 0), (4.6, 0), (0, -4.0)
draw(ax, I, cT)
draw(ax, R, cL)
draw(ax, U, cR)
draw(ax, F, cB)
arrow(ax, (-1.5, 0), (-3.0, 0), r"$\mathbf{R}$", (0, 0.4))
arrow(ax, (1.5, 0), (2.85, 0), r"$\mathbf{U}$", (0, 0.4))
arrow(ax, (0, -1.5), (0, -2.5), r"$\mathbf{F}$", (-0.5, 0))
arrow(ax, (-3.6, -1.2), (-1.6, -3.0), r"$\mathbf{V}$", (-0.5, -0.3))
arrow(ax, (3.6, -1.2), (1.6, -3.0), r"$\mathbf{R}$", (0.5, -0.3))
ax.set_aspect("equal")
ax.set_xlim(-6.6, 6.6)
ax.set_ylim(-5.9, 1.9)
ax.axis("off")
fig.tight_layout(pad=0.2)
fig.savefig(str(FIG / "L2F1.pdf"))
