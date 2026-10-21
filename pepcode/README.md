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
* `scripts/slides_seismicity.py` - global seismicity map for the slides, written as PNG to `../figures/slides`
  (USGS catalogue and NASA Blue Marble cached in `data/seismicity/`).
* Other slide figures, written as PDF to `../figures/slides`: `slides_boundary.py` (Lecture 16, P-wave meeting a
  boundary), `slides_anisotropy.py` (Lecture 16, point force in isotropic and TI media, with
  `fd2d_psv.TIPSVSolver`), `slides_L5F2_stills.py` (Lecture 17, stills standing in for the L5F2 animation),
  `slides_taup.py` (Lecture 18, PKKKKP, PcP/PcS and Pdiff/Sdiff ray paths and travel times with `obspy.taup`),
  `slides_healing.py` (Lecture 20, wavefront healing), `slides_bdkernels.py` (Lecture 20, banana-doughnut kernels
  of P, S and SS traveltimes at 60 degrees in PREM without its ocean, by the adjoint method with YSpec's Python package `yspec-seismo`;
  kernels and traces cached in `data/kernels/`, 4 threads by default), `slides_jb_prem.py` (Lecture 18,
  Jeffreys–Bullen against PREM), `slides_taylor.py` (Lecture 20, a Taylor test of the 2D kernel). The FD figures use `fd2d_psv` and are drawn like L5F1/L5F2; `L5F1.py`, `L5F2.py`,
  `slides_boundary.py`, `slides_anisotropy.py`, `slides_healing.py` and `L8F3.py` take `--webm` to write the movies in
  `../figures/anim/` for the slides (`fd2d_psv.movie`, with `imageio-ffmpeg`).
* `data/modes/` - cached mode computations (high-degree spheroidal catalogue extension, dispersion-corrected
  fundamental-mode frequencies).

Install and run with poetry:

    poetry install
    poetry run python scripts/make_figures.py           # all Python-made figures
    poetry run python scripts/L6F3.py                   # one figure
    poetry run python scripts/L5F1.py --animate         # animation for the finite-difference snapshots

The TikZ figure sources are in `../figsrc/` and are compiled with pdflatex.
