# pepseis

Code behind the figures and problem sets of the PEP seismology lectures.

* `pepseis/prem.py` - the PREM earth model (polynomials, gravity, pressure).
* `pepseis/sphray.py` - ray tracing and travel times in a spherically symmetric earth model.
* `pepseis/raytrace2d.py` - 2D Hamiltonian ray tracing through a gridded wave-speed field.
* `pepseis/fd2d_psv.py` - a small 2D P-SV finite-difference solver (staggered grid, sponge boundaries).
* `pepseis/elastic.py` - anisotropic elastic tensors, Christoffel equation, slowness sections.
* `pepseis/toytomo.py` - a toy straight-ray delay-time tomography problem.
* `scripts/L<n>F<m>.py` - the script for figure m of lecture file n (`lecture<n>.tex`); `S1F1.py` for the
  problem-set solutions. Figures are written to `../figures`.

Install and run with poetry:

    poetry install
    poetry run python scripts/make_figures.py           # all Python-made figures
    poetry run python scripts/L6F3.py                   # one figure
    poetry run python scripts/L5F1.py --animate         # animation for the finite-difference snapshots

The TikZ figure sources are in `../figsrc/` and are compiled with pdflatex.
