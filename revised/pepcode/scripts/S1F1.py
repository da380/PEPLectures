"""Problem Set 1, solution to Q10: a section through the slowness surface of a transversely
isotropic medium and the corresponding wave surface (the energy velocity as a function of
propagation direction), for the quasi-SV sheet, which has cusps in the wave surface."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pepseis.elastic import voigt_to_tensor, ti_voigt
from pepseis.paths import FIG

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 11})

rho = 3300.0
A = voigt_to_tensor(ti_voigt(240.0, 200.0, 30.0, 60.0, 72.0))   # GPa, x3 symmetry axis; F chosen small so that the qSV wave surface has cusps
e1, e3 = np.eye(3)[0], np.eye(3)[2]
phi = np.linspace(0, 2 * np.pi, 1441)
slow = np.zeros((3, len(phi), 2)); vel = np.zeros((3, len(phi), 2))
for m, ph in enumerate(phi):
    n = np.cos(ph) * e1 + np.sin(ph) * e3
    G = np.einsum("ijkl,j,l->ik", A, n, n) / rho
    c2, a = np.linalg.eigh(G)
    for k in range(3):
        c = np.sqrt(c2[k] * 1e9) / 1e3                    # km/s
        p = n / c                                          # s/km
        ak = a[:, k]
        v = np.einsum("ijkl,i,k,l->j", A, ak, ak, p) / rho * 1e9 / 1e6   # (GPa/(kg m^-3)) * s/km -> km/s
        slow[k, m] = p[[0, 2]]; vel[k, m] = v[[0, 2]]

fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.9))
cols = ["tab:green", "tab:red", "tab:blue"]
for k, col in enumerate(cols):
    axes[0].plot(slow[k, :, 0], slow[k, :, 1], color=col, lw=1.2)
    axes[1].plot(vel[k, :, 0], vel[k, :, 1], color=col, lw=1.2)
axes[0].set_title("slowness surface"); axes[0].set_xlabel(r"$p_1$ / s km$^{-1}$"); axes[0].set_ylabel(r"$p_3$ / s km$^{-1}$")
axes[1].set_title("wave surface"); axes[1].set_xlabel(r"$v_1$ / km s$^{-1}$"); axes[1].set_ylabel(r"$v_3$ / km s$^{-1}$")
for ax in axes:
    ax.set_aspect("equal"); ax.grid(True, ls=":", lw=0.5)
fig.tight_layout()
fig.savefig(str(FIG / "S1F1.pdf")); fig.savefig(str(FIG / "S1F1.png"), dpi=200)
