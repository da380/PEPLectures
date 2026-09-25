# Plan: regenerating lecture figures in Python

The 28 figures in `figures/` fall into four groups. Only the first two are worth regenerating;
the others are photographs, hand drawings or figures reproduced from papers.

## A. Simple computed figures — regenerate now (each a short matplotlib script)

| Figure | Used in | What it shows | Script / status |
|---|---|---|---|
| `L2F1.png` | Lecture 13 | Polar decomposition: circle → ellipse under U, R, V, F | `ex1_polar_decomposition.py` (written for the Lecture 12 worked examples; a generic version with a random or user-supplied F is a five-line change). |
| `L4F1.png` | Lecture 15 | Plane sections through the slowness surface of two anisotropic media | `slowness_surface.py`: build A_ijkl (e.g. cubic or TI with chosen constants), loop over directions in a plane, `eigh` of the Christoffel matrix, plot 1/c_k. Trivial once A is set up; the TI tensor is already built in `verify_examples2.py`. |
| `L6F5.png` | Lecture 17 | PREM α(r), β(r) | `prem.py`: PREM is defined by polynomials in r on 13 layers (Dziewonski & Anderson 1981, Table I); embed the coefficients (about 60 numbers) and evaluate ρ, α, β, κ, μ, g, p. The same module then feeds L10F1 and L12F4. |
| `L10F1.png` (upper) | Lecture 21 | PREM density | from `prem.py`. |
| `L10F1.png` (lower) | Lecture 21 | 1/ε from Clairaut's equation in PREM | `clairaut.py`: integrate ε'' + 8πGρ₀g₀⁻¹(ε'+ε/r) − 6ε/r² = 0 outward from r=0 with ε(0)=1, ε'(0)=0 (linear ODE, so one shot suffices), scale to satisfy the surface condition ε'(b) = [5Ω²b³/(2GM) − 2ε(b)]/b. Needs g₀(r) from `prem.py`; density discontinuities are handled automatically by the ODE solver if ρ₀ is piecewise (use `solve_ivp` layer by layer with continuity of ε and ε'). Expected surface value 1/ε ≈ 299 (the hydrostatic figure), central value ≈ 415, matching the existing plot. |
| `L12F4.png` | Lecture 23 | PREM density vs radius | from `prem.py`. |
| `L5F3.png` | Lecture 16 | Ray fan through a low-velocity anomaly forming a caustic | `ray_caustic.py`: 2D α(x,z) = α₀ − A exp(−|x−x₀|²/2σ²), integrate Hamilton's equations with `solve_ivp` for a fan of initial slownesses on a line, plot rays; the fold caustic appears behind the lens. The ray-fan script for the Lecture 16 examples (`ex5_ray_fan.py`) is the starting point. |
| `L6F3.png` | Lecture 17 | Ray paths in the Earth incl. the shadow zone | `earth_rays.py`: shoot rays in PREM (ray parameter loop) using the spherical ray equations dΔ/dr, dT/dr of the lecture with layer-by-layer integration and Snell's law at the CMB/ICB. Moderate effort (turning-point handling, integrable singularity: substitute r = r_t cosh u or use `quad` with the weight). A first version restricted to P and PKP is enough for the shadow-zone picture. |

## B. Computed figures needing a small numerical model — plan only

| Figure | What is needed |
|---|---|
| `L12F1.png`, `L12F2.png` (PREM toroidal/spheroidal eigenfrequencies vs degree) | Toroidal modes in PREM: a radial ODE eigenvalue problem per l (the Sturm–Liouville problem of Lecture 23), solvable by shooting or by a finite-element discretisation of the weak form (both ~100 lines with scipy; the FE route is cleaner for the discontinuities). Spheroidal modes require the coupled U, V, P (gravity) system with fluid outer core — a proper mode code (e.g. a Python port of MINEOS or use of an existing package); not worth writing for one figure. Suggest: regenerate the toroidal plot in Python, keep the spheroidal plot. |
| `L5F1.png`, `L5F2.png` (snapshots of a plane wave scattered by / passing through a heterogeneity) | 2D elastic finite-difference simulation (staggered grid, ~200 lines; or `devito`/`SPECFEM2D` output). Feasible but a project in itself; low priority since the current figures are fine. |

## C. Observational figures — keep as they are

`L6F2.png` (seismogram), `L6F4.png` (travel-time curves, likely from the ISC/ak135 tables),
`L7F1.png` (Woodhouse & Dziewonski 1989 model), `L7F2.png` (2D tomography toy example — could be
regenerated with the toy-tomography code of the Lecture 18 examples if desired), `L8F1.png` (Tape et
al. 2010), `L8F2.png` (Bozdağ et al. 2016), `L8F3.png` (Liu & Tromp 2008), `L11F1.png`, `L11F2.png`
(Bolivia 1994 seismogram and spectrum — regenerable from IRIS data with obspy if wanted),
`L12F3.png` (Benioff et al. 1961), `L12F5.png` (tomographic model compilation).

## D. Schematics — keep (or redraw in TikZ, not Python)

`L1F1.png` (hand-drawn reference/physical body map), `L3F1.png`, `L3F2.png`, `L3F3.png` (fault
diagram, elastic rebound cartoon, fault photograph), `L6F1.png` (hypocentre/epicentre sketch).

## Conventions for regenerated figures

* One script per figure in `python/`, writing PNG (dpi 200) and PDF to `figures/` with the
  *same file name* as the figure it replaces, so no `.tex` changes are needed; keep the originals
  under `figures/original/` until the replacements are approved.
* Matplotlib defaults with `font.size = 11`, serif text to sit with Palatino (`mathpazo`) in the
  notes: `rcParams['font.family'] = 'serif'`, `rcParams['mathtext.fontset'] = 'cm'`.
* Axis labels with units in the notes' style ("Radius / km", "Velocity / km s$^{-1}$").
* No titles inside figures (captions carry the description); colour-blind-safe two-colour scheme
  (e.g. black and `tab:red`) for P/S pairs, matching the existing plots.
* A `Makefile` target `make figures` running all scripts in order (`prem.py` first).
