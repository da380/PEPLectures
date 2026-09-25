# pepseis

Code behind the figures and problem sets of the PEP seismology lectures.

* `pepseis/prem.py` - the PREM earth model (polynomials, gravity, pressure, moduli).
* `pepseis/sphray.py` - ray tracing and travel times in a spherically symmetric earth model.
* `pepseis/raytrace2d.py` - 2D Hamiltonian ray tracing through a gridded wave-speed field.
* `pepseis/fd2d_psv.py` - a small 2D P-SV finite-difference solver (staggered grid, sponge boundaries).
* `pepseis/elastic.py` - anisotropic elastic tensors, Christoffel equation, slowness sections.
* `pepseis/toytomo.py` - a toy straight-ray delay-time tomography problem.
* `pepseis/kernel2d.py` - 2D adjoint construction of a cross-correlation delay-time kernel.
* `pepseis/clairaut.py` - Clairaut's equation for the hydrostatic ellipticity of PREM.
* `pepseis/ellipsoids.py` - Maclaurin and Jacobi ellipsoids (index symbols, bifurcation).
* `pepseis/modes.py` - normal modes of PREM via Alex Myhill's `mode_lab_2` (installed editable from
  `~/dev/mode_lab_2`): catalogues, eigenfunctions, toroidal and spheroidal sensitivity kernels, modal Q and
  dispersion corrections, extension of the spheroidal catalogue to high degree.
* `pepseis/seismograms.py` - access to the cached long-period records of the 2011 Tohoku earthquake
  (`data/seismograms/`, downloaded by `scripts/fetch_modes_data.py` and `scripts/fetch_splitting_data.py`).
* `scripts/L<n>F<m>.py` - the script for figure m of lecture file n (`lecture<n>.tex`, i.e. Lecture n+12);
  `S1F1.py` for the problem-set solutions. Figures are written as PDF to `../figures`.
* `scripts/run_yspec.py` - runs David's yspec code for the synthetic seismograms of Lecture 21 (outputs cached
  in `data/yspec/`); `scripts/check_kernel.py` - validation of the 2D kernel against direct perturbation.
* `data/modes/` - cached mode computations (high-degree spheroidal catalogue extension, dispersion-corrected
  fundamental-mode frequencies).

Install and run with poetry:

    poetry install
    poetry run python scripts/make_figures.py           # all Python-made figures
    poetry run python scripts/L6F3.py                   # one figure
    poetry run python scripts/L5F1.py --animate         # animation for the finite-difference snapshots

The TikZ figure sources are in `../figsrc/` and are compiled with pdflatex.
