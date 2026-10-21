"""Slides, Lecture 20: wavefront healing.  A plane P-wave passes a small low-velocity heterogeneity (comparable in
size to the pulse width); the delay imprinted on the wavefront is filled in as the wave moves on, while a weak
scattered wave spreads out (finite-difference simulation, see fd2d_psv.py, drawn as L5F1/L5F2).  Three snapshots.
Written as PDF to ../figures/slides; with --webm, a movie to ../figures/anim/healing.webm instead.
Usage:  poetry run python scripts/slides_healing.py [--webm]"""
import sys
import numpy as np
from pepseis.fd2d_psv import PSVSolver, gaussian_anomaly, snapshot_figure, movie
from pepseis.paths import FIG

h = 1.0 / 300
x = np.arange(0, 2.0, h); z = np.arange(0, 1.0, h)
alpha0, beta0 = 1.0, 0.58
width = 0.025
anom = gaussian_anomaly(x, z, centre=(0.7, 0.5), width=0.04, frac=-0.3)
alpha = alpha0 * (1 + anom); beta = beta0 * (1 + anom); rho = np.ones_like(alpha)
s = PSVSolver(x, z, alpha, beta, rho)
s.plane_p_wave(x0=0.35, width=width, alpha0=alpha0)
if "--webm" in sys.argv:
    movie(s, anom, str(FIG / "anim" / "healing"), np.linspace(0.2, 1.3, 90), vmax=0.4)
    sys.exit()
for k, t_end in enumerate([0.37, 0.62, 1.05], start=1):
    s.run(t_end=t_end)
    name = FIG / "slides" / f"healing_t{k}"
    snapshot_figure(s, anom, str(name), vmax=0.4)
    (FIG / "slides" / f"healing_t{k}.png").unlink()
