"""Ray fan for the two-layer earth model of Example 5 (examples/examples6.tex).

Homogeneous mantle (alpha_m = 13 km/s) above a homogeneous core (radius 3480 km,
alpha_c = 8 km/s). All ray paths are straight-line segments. Saves
examples/figures_ex/ex6_shadow_zone.png.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

b, c, am, ac = 6371.0, 3480.0, 13.0, 8.0
q1 = c / am

MANTLE = "#2a78d6"   # categorical slot 1 (blue)
CORE = "#eb6834"     # categorical slot 2 (orange)
INK = "#0b0b0b"
GREY = "#52514e"


def xy(r, delta):
    """Cartesian position (source at the top of the circle, delta measured clockwise)."""
    return r * np.sin(delta), r * np.cos(delta)


def mantle_ray(q):
    D = 2 * np.arccos(am * q / b)
    return [xy(b, 0.0), xy(b, D)]


def core_ray(q):
    Dm = np.arccos(am * q / b) - np.arccos(am * q / c)
    Dc = 2 * np.arccos(ac * q / c)
    return [xy(b, 0.0), xy(c, Dm), xy(c, Dm + Dc), xy(b, 2 * Dm + Dc)]


fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=200)
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

th = np.linspace(0, 2 * np.pi, 721)
ax.plot(b * np.sin(th), b * np.cos(th), color=INK, lw=1.0)
ax.plot(c * np.sin(th), c * np.cos(th), color=GREY, lw=0.8, ls="--")

# Mantle rays: equally spaced in epicentral angle up to the grazing ray.
for D in np.radians(np.arange(10, 114, 9)):
    q = (b / am) * np.cos(D / 2)
    P = mantle_ray(q)
    ax.plot(*zip(*P), color=MANTLE, lw=0.9)
Pg = mantle_ray(q1 * (1 - 1e-9))
ax.plot(*zip(*Pg), color=MANTLE, lw=1.6)

# Core rays: from just below the grazing parameter down to the vertical ray.
for q in np.concatenate([q1 * np.array([0.999, 0.99, 0.97, 0.94, 0.9]),
                         np.linspace(0.82 * q1, 0.0, 9)]):
    P = core_ray(q)
    ax.plot(*zip(*P), color=CORE, lw=0.9)

# Shadow-zone arcs on the surface, either side of the source.
D1 = 2 * np.arccos(c / b)
dmin = 2 * np.pi - (D1 + 2 * np.arccos(ac / am))
for sgn in (+1, -1):
    ths = sgn * np.linspace(D1, dmin, 200)
    ax.plot(1.015 * b * np.sin(ths), 1.015 * b * np.cos(ths), color=INK, lw=4.0,
            solid_capstyle="butt")

# Annotations.
ax.plot(0, b, marker="*", color=INK, ms=12, zorder=5)
ax.text(0, 1.06 * b, "source", ha="center", va="bottom", fontsize=9, color=INK)
ax.text(0, 0.45 * c, "core", ha="center", va="center", fontsize=9, color=GREY)
ax.text(-0.83 * b, 0.55 * b, "mantle", ha="center", va="center", fontsize=9, color=GREY)
for sgn, txt in ((+1, r"$113.8^\circ$"), (-1, r"$113.8^\circ$")):
    x, y = xy(1.15 * b, sgn * D1)
    ax.text(x, y, txt, ha="center", va="center", fontsize=8, color=INK)
for sgn in (+1, -1):
    x, y = xy(1.15 * b, sgn * dmin)
    ax.text(x, y, r"$142.2^\circ$", ha="center", va="center", fontsize=8, color=INK)
x, y = xy(1.26 * b, np.radians(128))
ax.text(x, y, "shadow\nzone", ha="center", va="center", fontsize=9, color=INK)
x, y = xy(1.26 * b, -np.radians(128))
ax.text(x, y, "shadow\nzone", ha="center", va="center", fontsize=9, color=INK)

ax.plot([], [], color=MANTLE, lw=1.2, label=r"mantle rays ($q>c/\alpha_{\mathrm{m}}$)")
ax.plot([], [], color=CORE, lw=1.2, label=r"core-refracted rays ($q<c/\alpha_{\mathrm{m}}$)")
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=8, ncol=1)

ax.set_aspect("equal")
ax.set_xlim(-1.4 * b, 1.4 * b)
ax.set_ylim(-1.3 * b, 1.18 * b)
ax.set_xlabel("distance (km)")
ax.set_ylabel("distance (km)")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "examples", "figures_ex")
os.makedirs(out, exist_ok=True)
fig.savefig(os.path.join(out, "ex6_shadow_zone.png"), dpi=200, bbox_inches="tight", facecolor="white")
print("saved", os.path.join(out, "ex6_shadow_zone.png"))
