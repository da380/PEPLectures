"""Lecture 24, Fig. 7: sensitivity kernels of toroidal and spheroidal modes of PREM.  Toroidal
kernels (shear modulus, density) follow analytically from the toroidal eigenvalue problem of the
lecture (pepseis.modes.PremModes.toroidal_kernels); spheroidal kernels (bulk modulus, shear
modulus, density, including self-gravitation) are obtained from mode_lab_2 by element-wise
re-assembly and Rayleigh's principle (pepseis.modes.spheroidal_kernels).  All are plotted as
fractional kernels: delta omega/omega = int [K_kappa (d kappa/kappa) + K_mu (d mu/mu) + K_rho (d rho/rho)] dr."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from pepseis import prem
from pepseis.modes import PremModes, spheroidal_kernels
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})
M = PremModes(order=5, element_size_km=50.0)
fig, axes = plt.subplots(2, 4, figsize=(10.5, 7.4), sharey=False, constrained_layout=True)
for ax, (l, n) in zip(axes[0], [(2, 0), (2, 1), (10, 0), (40, 0)]):
    f, r, w, Kmu, Krho = M.toroidal_kernels(l, n)
    rc = np.clip(r, prem.R_CMB + 1e-6, 6368 - 1e-6)
    kap, mu = prem.moduli(rc); rho = prem.density(rc); om = 2 * np.pi * f
    ax.plot(Kmu * mu / om * 1e3, r, "k", lw=1.3, label=r"$\mu$")
    ax.plot(Krho * rho / om * 1e3, r, "tab:red", lw=1.3, label=r"$\rho$")
    ax.set_title(r"${}_{%d}T_{%d}$, %.3f mHz" % (n, l, 1e3 * f), fontsize=10)
    ax.set_ylim(3400, 6371)
for ax, (l, n) in zip(axes[1], [(2, 0), (2, 1), (10, 0), (0, 0)]):
    f0, rc, Kk, Km, Kr, _ = spheroidal_kernels(M, l, n)
    rcc = np.clip(rc, 1e-3, 6368 - 1e-6)
    kap, mu = prem.moduli(rcc); rho = prem.density(rcc); om = 2 * np.pi * f0
    ax.plot(Kk * kap / om * 1e3, rc, "tab:blue", lw=1.3, label=r"$\kappa$")
    ax.plot(Km * mu / om * 1e3, rc, "k", lw=1.3, label=r"$\mu$")
    ax.plot(Kr * rho / om * 1e3, rc, "tab:red", lw=1.3, label=r"$\rho$")
    ax.set_title(r"${}_{%d}S_{%d}$, %.3f mHz" % (n, l, 1e3 * f0), fontsize=10)
    ax.set_ylim(0, 6371)
    print("%dS%d" % (n, l), "int(Kk kap + Km mu + Kr rho) dr / omega =",
          np.nansum((Kk * kap + np.nan_to_num(Km) * mu + Kr * rho) * np.gradient(rc) * 1e3) / om)
for ax in axes.ravel():
    ax.axvline(0, color="0.8", lw=0.6)
    for rb in (prem.R_CMB, prem.R_ICB):
        ax.axhline(rb, color="0.6", lw=0.6, ls="--")
    ax.grid(True, ls=":", lw=0.4); ax.xaxis.set_major_locator(MaxNLocator(3))
for ax in axes[1]:
    ax.set_xlabel(r"fractional kernel / km$^{-1}$")
for ax in axes[:, 0]:
    ax.set_ylabel("radius / km")
axes[0, 0].legend(loc="lower right", fontsize=9, frameon=False); axes[1, 0].legend(loc="lower right", fontsize=9, frameon=False)
fig.savefig(str(FIG / "L12F7.pdf")); fig.savefig(str(FIG / "L12F7.png"), dpi=200)
