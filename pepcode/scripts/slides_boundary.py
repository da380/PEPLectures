"""Slides, Lecture 16: a P-wave from an explosive point source meeting a horizontal boundary between two
homogeneous media (finite-difference simulation, see fd2d_psv.py), drawn as in L5F1/L5F2.  Two snapshots:
the incident wave reaching the boundary, and later the reflected and transmitted P- and S-waves.
The grid is larger than the window drawn, so that up to the second snapshot no wave reaches the absorbing
layers at its edges and nothing returns from them: the boundary is the only one the waves meet.
Written as PDF to ../figures/slides; with --webm, a movie to ../figures/anim/boundary.webm instead.
Usage:  poetry run python scripts/slides_boundary.py [--webm]"""
import sys
from types import SimpleNamespace
import numpy as np
from pepseis.fd2d_psv import PSVSolver, snapshot_figure, movie
from pepseis.paths import FIG

h = 1.0 / 300
x = np.arange(0, 2.6, h); z = np.arange(0, 1.8, h)
X, Z = np.meshgrid(x, z, indexing="ij")
zb = 0.7                                        # boundary; the lower medium is faster and denser
lower = Z < zb
alpha = np.where(lower, 1.4, 1.0); beta = np.where(lower, 0.81, 0.58); rho = np.where(lower, 1.2, 1.0)
s = PSVSolver(x, z, alpha, beta, rho, sponge_z=True)
view = (x >= 0.3) & (x < 2.3), (z >= 0.3) & (z < 1.3)    # the 2 x 1 window drawn

# explosion with a Gaussian moment function, so that no static field is left at the source
xs, zs = 1.3, 0.97
ix, iz = int(round(xs / h)), int(round(zs / h))
sig, t0 = 0.018, 0.08
amp = 2.5e4                                     # scales the incident pulse to displacements of order one
def moment_rate(t):
    return -amp * (t - t0) / sig ** 2 * np.exp(-(t - t0) ** 2 / (2 * sig ** 2))

boundary = -lower.astype(float)                 # snapshot_figure draws its half-minimum contour (dashed)
def advance(solver, t_end):
    while solver.t < t_end - 1e-12:
        solver.step(explosion=(ix, iz, moment_rate(solver.t)))
if "--webm" in sys.argv:
    movie(s, boundary, str(FIG / "anim" / "boundary"), np.linspace(0.05, 0.75, 85), advance=advance, view=view,
          vmax=0.6)
    sys.exit()
for name, t_end in [("boundary_t1", 0.36), ("boundary_t2", 0.62)]:
    advance(s, t_end)
    shown = SimpleNamespace(x=x[view[0]], z=z[view[1]], ux=s.ux[np.ix_(*view)], uz=s.uz[np.ix_(*view)])
    snapshot_figure(shown, boundary[np.ix_(*view)], str(FIG / "slides" / name), vmax=0.6)
    (FIG / "slides" / (name + ".png")).unlink()
