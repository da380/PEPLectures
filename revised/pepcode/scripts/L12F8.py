"""Lecture 24: radial eigenfunctions of some less typical spheroidal modes of PREM: the Slichter
mode 1S1 (translation of the inner core), a core-mantle-boundary Stoneley mode 2S16, an inner-core
boundary Stoneley mode 4S8, and an inner-core mode 11S2.  Modes were identified by the partition of
their kinetic energy between the inner core, outer core and mantle."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mode_lab_2 as ml
from mode_lab_2.spheroidal import physical_ocean_fraction
from pepseis import prem
from pepseis.modes import PremModes, reference_frequency, catalogue
lc, nc, fc = catalogue("spheroidal")
def cat_f(l, n):
    return fc[(lc == l) & (nc == n)][0]
from pepseis.paths import FIG
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm", "font.size": 10})

def state(pm, l, f_target_hz, f_max_hz, pick=None):
    f, vec, pen = pm.spheroidal(l, f_max_hz)
    k = int(np.argmin(np.abs(f - f_target_hz))) if pick is None else pick(f, vec, pen)
    r = np.linspace(1.0, 6370.999, 1600)
    U, V, P = ml.sample_spheroidal_state(pen, vec[:, k], r * 1e3)
    return f[k], r, U, V

def ocean_pick(f, vec, pen):
    cand = [(physical_ocean_fraction(pen, vec[:, k]), k) for k in range(len(f)) if 0.5e-3 < f[k] < 1.5e-3]
    return max(cand)[1]

M50 = PremModes(order=5, element_size_km=50.0); M100 = PremModes(order=5, element_size_km=100.0)
panels = [(M100, 1, cat_f(1, 1), 0.3e-3, None, r"${}_1S_1$" "\nSlichter mode"),
          (M50, 16, cat_f(16, 2), 4e-3, None, r"${}_2S_{16}$" "\nCMB Stoneley mode"),
          (M50, 8, cat_f(8, 4), 4e-3, None, r"${}_4S_{8}$" "\nICB Stoneley mode"),
          (M50, 2, cat_f(2, 11), 5e-3, None, r"${}_{11}S_{2}$" "\ninner-core mode")]
# drawn at its printed width (\textwidth) so that the fonts appear at their nominal size
fig, axes = plt.subplots(1, 4, figsize=(6.5, 3.6), sharey=True, constrained_layout=True)
for ax, (pm, l, ft, fmax, pick, name) in zip(axes, panels):
    f, r, U, V = state(pm, l, ft, fmax, pick)
    zeta = np.sqrt(l * (l + 1)); scale = max(np.abs(U).max(), np.abs(zeta * V).max())
    ax.plot(U / scale, r, "k", lw=1.3, label=r"$U$"); ax.plot(zeta * V / scale, r, "tab:red", lw=1.3, ls="--", label=r"$\zeta V$")
    per = "%.1f h" % (1 / f / 3600) if f < 1e-4 else "%.1f min" % (1 / f / 60)
    ax.set_title(name + "\n%.4f mHz, %s" % (1e3 * f, per), fontsize=8)
    for rb in (prem.R_CMB, prem.R_ICB):
        ax.axhline(rb, color="0.6", lw=0.6, ls="--")
    ax.axvline(0, color="0.8", lw=0.6); ax.set_xlim(-1.1, 1.1); ax.grid(True, ls=":", lw=0.4)
    ax.set_ylim(0, 6371)
    print(name, "f = %.5f mHz" % (1e3 * f))
fig.supxlabel("normalised eigenfunction")
axes[0].set_ylabel("radius / km"); axes[0].legend(loc="lower right", fontsize=9, frameon=False)
fig.savefig(str(FIG / "L12F8.pdf")); fig.savefig(str(FIG / "L12F8.png"), dpi=200)
