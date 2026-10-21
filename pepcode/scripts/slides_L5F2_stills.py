"""Slides, Lecture 17: stills of the Fig. 2 simulation (same set-up as L5F2.py), standing in for the
animation: t = 0 (plane wave in the homogeneous region) and t = 0.9 (wavefront distorted inside the anomaly).
Written as PDF to ../figures/slides.
Usage:  poetry run python scripts/slides_L5F2_stills.py"""
import numpy as np
from pepseis.fd2d_psv import PSVSolver, gaussian_anomaly, snapshot_figure
from pepseis.paths import FIG

h = 1.0 / 300
x = np.arange(0, 2.0, h); z = np.arange(0, 1.0, h)
alpha0, beta0 = 1.0, 0.58
width = 0.025
anom = gaussian_anomaly(x, z, centre=(0.9, 0.5), width=0.3, frac=-0.2)
alpha = alpha0 * (1 + anom); beta = beta0 * (1 + anom); rho = np.ones_like(alpha)
s = PSVSolver(x, z, alpha, beta, rho)
s.plane_p_wave(x0=0.3, width=width, alpha0=alpha0)
snapshot_figure(s, anom, str(FIG / "slides/L5F2_t00"), vmax=1.0)
s.run(t_end=0.9)
snapshot_figure(s, anom, str(FIG / "slides/L5F2_t09"), vmax=1.0)
