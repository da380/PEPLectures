"""Radial eigenfunctions W(r) of the toroidal modes 0T2, 1T2, 2T2 of a
homogeneous sphere (Example 1 of the worked examples for Lecture 23).

W(r) = j_l(k r), with k b a root of  x j_l'(x) = j_l(x)  (traction-free surface).
Writes examples/figures_ex/ex12_toroidal.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import spherical_jn
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "examples", "figures_ex", "ex12_toroidal.png")

b = 6371.0  # km
l = 2


def freq_eq(x):
    return x * spherical_jn(l, x, derivative=True) - spherical_jn(l, x)


def roots(nmax, xmax=40.0):
    xs = np.linspace(0.05, xmax, 40001)
    f = freq_eq(xs)
    out = []
    for i in range(len(xs) - 1):
        if f[i] * f[i + 1] < 0:
            out.append(brentq(freq_eq, xs[i], xs[i + 1], xtol=1e-13))
            if len(out) == nmax:
                break
    return out


x_n = roots(3)
r = np.linspace(0.0, b, 1001)

# Categorical palette, fixed order (blue, orange, aqua).
colours = ["#2a78d6", "#eb6834", "#1baf7a"]

fig, ax = plt.subplots(figsize=(6.4, 4.0), facecolor="white")
ax.set_facecolor("white")
ax.axhline(0.0, color="#b0b0b0", lw=0.8, zorder=0)
for n, (x, c) in enumerate(zip(x_n, colours)):
    W = spherical_jn(l, x * r / b)
    W = W / np.max(np.abs(W))
    ax.plot(r, W, color=c, lw=1.8, label=r"${}_{%d}T_{2}$" % n)
ax.set_xlim(0.0, b)
ax.set_ylim(-1.05, 1.05)
ax.set_xlabel("Radius $r$ (km)")
ax.set_ylabel("Radial eigenfunction $W(r)$ (normalised)")
ax.grid(True, color="#e6e6e6", lw=0.6)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(OUT, dpi=200, facecolor="white")
print("roots kb:", ["%.4f" % x for x in x_n])
print("wrote", os.path.abspath(OUT))
