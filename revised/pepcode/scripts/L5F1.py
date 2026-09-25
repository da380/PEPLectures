"""Lecture 17, Fig. 1: plane P-wave scattered by a heterogeneity of size comparable to the
wavelength (finite-difference simulation, see fd2d_psv.py)."""
import numpy as np
import sys
from pepseis.fd2d_psv import PSVSolver, gaussian_anomaly, snapshot_figure, animate
from pepseis.paths import FIG

h = 1.0 / 300
x = np.arange(0, 2.0, h); z = np.arange(0, 1.0, h)
alpha0, beta0 = 1.0, 0.58
width = 0.025                                  # pulse width; wavelength ~ 2 pi width
anom = gaussian_anomaly(x, z, centre=(0.9, 0.5), width=0.06, frac=-0.35)
alpha = alpha0 * (1 + anom); beta = beta0 * (1 + anom); rho = np.ones_like(alpha)
s = PSVSolver(x, z, alpha, beta, rho)
s.plane_p_wave(x0=0.35, width=width, alpha0=alpha0)
if "--animate" in sys.argv:          # python3 L5F1.py --animate  ->  figures/anim/L5F1.gif
    animate(s, anom, str(FIG / "anim/L5F1"), t_end=1.5, n_frames=60)
else:
    s.run(t_end=1.05)
    snapshot_figure(s, anom, str(FIG / "L5F1"), vmax=0.5)
