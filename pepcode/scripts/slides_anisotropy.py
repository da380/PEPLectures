"""Slides, Lecture 16: waves from a vertical point force in an isotropic medium and in a strongly anisotropic one
(transversely isotropic, axis vertical; Thomsen alpha0 = 1, beta0 = 0.5, epsilon = 0.3, delta = -0.1), at the same
time (finite-difference simulation, see fd2d_psv.py, drawn as L5F1/L5F2).  In the isotropic medium the P and S
wavefronts are circles; in the anisotropic one the quasi-P and quasi-S fronts are not, the particle motion is
neither along nor across them, and the quasi-S front has cusps.  The grid is larger than the window drawn, so no
wave reaches its absorbing edges.  Written as PDF to ../figures/slides; with --webm, also movies to
../figures/anim/anisotropy_*.webm (each scaled as its snapshot).
Usage:  poetry run python scripts/slides_anisotropy.py [--webm]"""
import subprocess
import sys
from types import SimpleNamespace
import warnings
import numpy as np
from pepseis.fd2d_psv import TIPSVSolver, snapshot_figure, movie
from pepseis.paths import FIG

warnings.filterwarnings("ignore", "No contour levels")   # nothing to outline: there is no anomaly
h = 1.0 / 300
x = np.arange(0, 3.2, h); z = np.arange(0, 3.2, h)
ones = np.ones((x.size, z.size))
view = (x >= 0.4) & (x < 2.8), (z >= 0.4) & (z < 2.8)
xs = zs = 1.6
ix, iz = int(round(xs / h)), int(round(zs / h))
sig, t0, t_end = 0.018, 0.08, 0.85

def force(t):     # a Gaussian pulse differentiated, so that no static displacement is left
    return -(t - t0) / sig ** 2 * np.exp(-(t - t0) ** 2 / (2 * sig ** 2))

alpha0, beta0, eps, delta = 1.0, 0.5, 0.3, -0.1
c33, c55 = alpha0 ** 2, beta0 ** 2
MEDIA = {
    "anisotropy_iso": (c33, c33 - 2 * c55, c33, c55),
    "anisotropy_ti": (c33 * (1 + 2 * eps), np.sqrt(2 * delta * c33 * (c33 - c55) + (c33 - c55) ** 2) - c55, c33, c55),
}
X, Z = np.meshgrid(x[view[0]], z[view[1]], indexing="ij")
far = np.hypot(X - xs, Z - zs) > 0.15            # scale away from the source, where the near field is large
for name, (c11, c13, c33_, c55_) in MEDIA.items():
    s = TIPSVSolver(x, z, c11 * ones, c13 * ones, c33_ * ones, c55_ * ones, ones, sponge_z=True)
    while s.t < t_end - 1e-12:
        s.step(force=(ix, iz, 0.0, force(s.t)))
    ux, uz = s.ux[np.ix_(*view)], s.uz[np.ix_(*view)]
    scale = np.hypot(ux, uz)[far].max()
    shown = SimpleNamespace(x=x[view[0]], z=z[view[1]], ux=ux / scale, uz=uz / scale)
    snapshot_figure(shown, np.zeros_like(ux), str(FIG / "slides" / name), vmax=0.35)
    (FIG / "slides" / (name + ".png")).unlink()
    pdf = str(FIG / "slides" / (name + ".pdf"))                 # the figure is 2:1 and the view square:
    subprocess.run(["pdfcrop", "--margins", "4", pdf, pdf], check=True, capture_output=True)   # trim the sides
    if "--webm" in sys.argv:
        s = TIPSVSolver(x, z, c11 * ones, c13 * ones, c33_ * ones, c55_ * ones, ones, sponge_z=True)
        def advance(solver, t):
            while solver.t < t - 1e-12:
                solver.step(force=(ix, iz, 0.0, force(solver.t)))
        movie(s, np.zeros_like(ones), str(FIG / "anim" / name), np.linspace(0.05, 1.0, 96), advance=advance,
              view=view, unit=scale, vmax=0.35)
