# 2026 revision: running notes (started 25 September 2026)

## Decisions taken
* Lectures renumbered 13–24 to match the 2026 syllabus; lecture6 retitled "Spherical Earth structure".
* `lectures.tex` is the single document built (`./build.sh`); the `lectureN.tex` files are bodies only
  (no preamble, no title) and are `\input` by the master. Titles are set in `lectures.tex`.
  Author and term appear on the overall title page only.
* `revised/figures/` holds links to the untouched originals plus regenerated figures; sources in `revised/figsrc/`.
* `../make_diffs.sh` wraps each body in last year's preamble so `latexdiff` can compare with the committed original.

## Problem sets (being revised alongside the lectures)
* Problem Set 1, Q2 split (25 Sept): Q2 is now energy conservation for finite elasticity (identity, the flux
  $s_j=-T_{ij}v_i$, its meaning, conservation of total energy); Q3 is the short derivation of the linearised
  equations and natural boundary conditions from the given Lagrangian. This discharges the promises in
  Lectures 13 and 15. Set 1 now has nine questions.

## Per-lecture log
* Lecture 13 (lecture1.tex): renumbered; Fig. 1 redrawn in TikZ; footnote added noting that reflections
  are excluded by $J>0$. No other changes.
* Lecture 14 (lecture2.tex): renumbered; "1B dynamics course" → "Part IB dynamics course"; Fig. 1
  (polar decomposition) regenerated from an explicit $\mathbf{F}$ (`figsrc/L2F1.py`). $U$ kept for both
  the stretch and the auxiliary strain energy function (David's decision).
* Lecture 14: repeated sentence "Let $\bphi^0$ be an equilibrium configuration of a body." cut. Transversely
  isotropic tensor rewritten in Love's $A,C,F,L,N$ notation (verified numerically: components, TI invariance,
  isotropic limit). Problem Set 1 Q4 and its solution changed to match; the solution now also checks the
  first-order qP speed against the exact axial value $C/\rho$. Old coefficients map as $\lambda=A-2N$,
  $\mu=N$, $4\xi=F-A+2N$, $\zeta=N-L$, $8\gamma=A+C-2F-4L$.
* To do at Lecture 24: its footnote on the five moduli of a spherically symmetric model could cite Love's
  $A,C,F,L,N$ from Lecture 14.
* Lecture 15 (lecture3.tex): "To rectify things..." sentence reworded; Fig. 1 (2D fault cross-section, replacing
  the Pearson figure) and Fig. 2 (elastic rebound, three panels) redrawn in TikZ (`figsrc/L3F1.tex`,
  `figsrc/L3F2.tex`); Figs. 1 and 3 now cited in the text; energy conservation now explicitly deferred to
  the first problem set.
* Problem Set 1: new Q4 (stress glut of a planar fault in an isotropic body; Green-function representation
  reduced to a single integral over the stress glut; Burridge--Knopoff form for a fault; moment tensor as a
  remark in the solution). Later questions renumbered 5–8 automatically. Adapted from `examples/examples3.tex`
  Example 1(a) as David directed.
* Lecture 16 (lecture4.tex): Fig. 1 regenerated (`figsrc/L4F1.py`): slowness sections for a transversely
  isotropic medium (section through the axis) and single-crystal olivine (Abramson et al. 1997 constants);
  caption rewritten. No text changes.
* Problem Set 1, Q6: now (a) quasi-P and (b) quasi-S speeds by degenerate perturbation theory (David: the old
  verbal caveat about not examining degenerate calculations is dropped). Solution verified numerically;
  in the (SV, SH) basis the 2x2 matrix is diagonal; SH result is exact.
* Lecture 17 (lecture5.tex): caustic paragraph reworded (envelope of the ray family; amplitude infinite on
  the caustic); $H_k$ used consistently; footnote on Euler's theorem now points to a new non-examinable
  appendix with the proof; Fig. 3 now cited; new Fig. 4 (rays through a random medium) with a sentence on the
  generic occurrence of caustics. Figures regenerated: L5F1/L5F2 from a new 2D P-SV finite-difference code
  (`figsrc/fd2d_psv.py`, 4th-order staggered grid, sponge boundaries; `--animate` writes GIFs to
  `figures/anim/`), L5F3/L5F4 from a new 2D Hamiltonian ray tracer (`figsrc/raytrace2d.py`).
* Lecture 17 figures (David's choices): white-background "Blues" colour map with black arrows, arrows halved
  in size; the FD initial condition is now sampled consistently on the staggered grid (no residual at the
  start position); Fig. 1 colour scale saturated at 0.5 so the scattered waves are visible.
* Problem Set 1 checked in full (25 Sept): all solutions re-derived and correct. Wording: Q1(b) now asks for
  the condition under which angular momentum changes only through the boundary torque; Q9 "with depth";
  Q5 states that $f_j$, $t_j$ and $\hat{\mathbf{n}}$ are evaluated at $\mathbf{x}'$, and the fault result is
  written with $\hat n_k(\mathbf{x}')$ (holds for non-planar faults).
* Problem Set 1, new Q10: energy density/flux of a plane wave in an anisotropic body (equipartition, energy
  velocity $v_j=A_{ijkl}a_ia_kp_l/\rho$, $\mathbf v\cdot\hat{\mathbf p}=c$, isotropic check), Hellmann--Feynman
  style proof that $\mathbf v=\partial H/\partial\mathbf p$ (normal to the slowness surface, ray direction,
  group velocity), and the wave surface of a point source with polar reciprocity and cusps. Solution figure
  `figsrc/S1F1.py` (slowness and wave curves of a TI medium). Set 1 now has ten questions.
* `figsrc/elastic.py` holds the shared anisotropic-tensor helpers (used by L4F1.py and S1F1.py).

## Problem Set 2: David's wishes (25 Sept), to act on when we reach it
* Ten questions in total.
* Add a gravitational binding energy question (crib from `examples/all_examples.pdf`).
* Add a question on the Coriolis splitting equations (or similar).
* Extend the moment-tensor/point-source question: have them argue the point-source approximation from the
  length scale of the eigenfunction strains being large compared with the fault, given a frequency limit.
* The adjoint question is fine as is, but consider replacing it by an analogue problem (acoustic waves, or
  another physics application) with the same core ideas and less index gymnastics.
* Instead of a tomography question: a guided Herglotz--Wiechert question (agreed 25 Sept). Steps:
  (1) tau(q) = T - q Delta = 2 int sqrt(eta^2-q^2) dr/r, d tau/dq = -Delta, hence dT/dDelta = q;
  (2) with eta = r/v monotone, Delta(q) = 2 int_q^{eta_b} q/sqrt(eta^2-q^2) (d ln r/d eta) d eta;
  (3) Abel inversion using the supplied identity int_{eta_1}^{eta} q dq / sqrt((eta^2-q^2)(q^2-eta_1^2)) = pi/2
      to get ln(b/r_1) = (1/pi) int_{eta_1}^{eta_b} Delta(q)/sqrt(q^2-eta_1^2) dq = (1/pi) int_0^{Delta_1} arccosh(q/eta_1) dDelta;
  (4) interpretation: v(r) from a travel-time curve; low-velocity zones (shadow) and velocity jumps (triplications).
  Solution to include a numerical check on the PREM mantle using figsrc/sphray.py.
* Lecture 18 (lecture6.tex): Fig. 1 redrawn in TikZ (cross-section: fault, hypocentre, epicentre, wavefronts,
  seismometer); Fig. 2 (Loma Prieta at KEV) kept; Figs. 3–5 regenerated from new code: `figsrc/prem.py`
  (PREM polynomials, checked: mass, g, p, interface values) and `figsrc/sphray.py` (spherical ray tracing by
  the lecture's Delta/T integrals; validated: PcP 8m29s, ScS 15m34s, PKIKP(180) 20m09s, P to 98 deg, PKP
  from 145 deg). Fig. 3 is now cited in the Oldham sentence; captions of Figs. 3 and 4 rewritten.
  Travel-time curves restricted to rays turning below 220 km (ducted rays in the near-constant r/v layer
  clutter the plot) and to Delta <= 180.
* Lecture 19 (lecture7.tex): Fig. 2 regenerated (`figsrc/L7F2.py`: straight-ray toy tomography, 20 sources x
  14 receivers, 48x48 cells, norm damping with lambda set by chi^2/n = 1); caption rewritten. Fig. 1
  (Woodhouse & Dziewonski 1989) kept. No text changes. David mentioned his `pygeoinf` package (PyPI) for
  function-space toy tomography; possibly useful for the Set 2 resolution question.
* Lecture 19: four further figures from the toy problem (`figsrc/toytomo.py` module; scripts L7F3–L7F6):
  null-space model (Fig. 3), damping parameter and trade-off curve (Fig. 4), checkerboard resolution test
  (Fig. 5, with a new paragraph on synthetic resolution tests pointing to the resolution matrix in Set 2),
  posterior standard deviation (Fig. 6). Each cited by a sentence or two of new text at the relevant point.
* Lecture 18, Fig. 4 caption now states that only rays turning below 220 km are shown.
* Lecture 19: closing paragraph of the Bayesian section rewritten (David: less hard on Bayesian methods;
  regularisation is an implicit prior; Bayesian formulation is honest about the prior and gives uncertainties),
  with a footnote on function-space formulations.
* Lecture 19 (later the same day): new paragraph on finding the null space via the eigendecomposition of
  A^T A, leading to the SVD (two new equations) with Fig. 4 (singular values and two singular vectors); the
  damped solution written in the singular basis (new equation) at the end of the regularisation section;
  the discrepancy principle named as the way to fix lambda (David: never mention L-curve corners); the
  Bayesian figure now uses a correlated Gaussian prior (corr. length 0.15, std 0.003) and shows input,
  posterior mean and posterior std; checkerboard test replaced by spike tests (3x3 block in a well-covered and a poorly covered cell; recovered/input amplitude shown; links to the resolution matrix). Figures renumbered 1–7,
  equation labels renumbered positionally (43 equations).
* Elastic tensor with a source present (student's point from last year): Lecture 14 gains a footnote after
  Hooke's law, Lecture 15 a sentence after eq. (15) and Lecture 21 a sentence after eq. (41), all stating
  that A_ijkl is the second derivative of the time-independent elastic part U(x,C) of W, the stress-glut term
  being first order in s and contributing only force terms.
* Python code reorganised into a poetry package `revised/pepcode/` (package `pepseis`: prem, sphray, raytrace2d,
  fd2d_psv, elastic, toytomo, paths; figure scripts in `scripts/`, `make_figures.py` runs them all; local
  `.venv` via `poetry install`; pygeoinf is a dependency). TikZ sources stay in `revised/figsrc/`.
* Lecture 19: toy problem coarsened to a 24x24 pixel basis (seed 6 for even coverage); text now names the
  pixel basis and gives the explicit sparse form of A for a block parameterisation; new non-examinable
  subsection "Inversion in function spaces" with Fig. 8 from pygeoinf (Tutorial 10, plane geometry) and a
  Stuart (2010) footnote; spike test uses a single block.
* Lecture 20 (lecture8.tex): Figs. 1 and 2 now cited; summary item (ii) no longer mentions being given other
  misfits; two new figures (L8F4 construction of a 2D cross-correlation delay-time kernel, L8F5 the finished
  kernel with a cross-section) from `pepseis/kernel2d.py` (adjoint run in reversed time, K_alpha =
  -2 rho alpha int div u div u' dt, verified against direct perturbations to ~10% at 2% delta-alpha;
  `scripts/check_kernel.py`). Text: paragraph on the construction; sentence on why there is no doughnut
  hole in 2D (quarter-period phase shift of a line-source far field). FD solver gained point sources
  (explosion, body force), a fixed-dt option and optional z sponges.
* Lecture 18: Figs. 4 and 5 reduced to 0.62 and 0.7 textwidth to improve placement at the end of the lecture.
* References section added at the end of lectures.tex (hand-written list, author-year in text); in-text citations made consistent (Dziewonski & Gilbert year corrected 1975 -> 1971; Ishii & Tromp 1999, Lau et al. 2017, Backus & Gilbert 1961, Gilbert 1971, Woodhouse & Dziewonski 1989, Chandrasekhar 1969, Newton 1687, Truesdell & Noll 1965, Hörmander 1971, Abramson et al. 1997 added). DAVID TO CHECK the bibliographic details, especially Maitra & Al-Attar 2024 page numbers.
* Lecture 20: derivation of the adjoint source for the cross-correlation delay time added (four equations,
  labels renumbered positionally) in place of "details left to the interested reader", since the former
  problem-set question on this has gone. Assumes s_obs(t - tau_bar) = s(t).
* Lecture 20 (later): Bozdağ et al. (2016) figure replaced by three figures from the open-access GLAD-M35
  paper (Cui et al. 2024; PDF supplied by David as figures/ggae270.pdf, crops made at 400 dpi): Fig. 1
  earthquake distribution (in the numerical-simulation section), Fig. 3 map views at seven depths, Fig. 4
  uncertainty maps (with a new paragraph on Hessian-based uncertainty), Fig. 5 comparison of six joint P/S
  models at 2800 km (with a paragraph on how each model was made). Kernel figures are now Figs. 6–7, Liu &
  Tromp Fig. 8. Kernel computation box padded by 0.5 on every side and trimmed for plotting (David).
  Licence statements dropped from captions ("From Cui et al. (2024)"). Six references added (Cui, Koelemeijer,
  Lu, Simmons, Thrastarson, van Herwaarden) — DAVID TO CHECK details. "10^2–10^3 earthquakes" -> "10^3".
* Lecture 20 reordered (David): banana-doughnut section now precedes a new closing section "Global waveform
  tomographic models" (GLAD-M35 maps, uncertainty, model comparison). Figures: 1 earthquakes, 2 Tape, 3–4
  kernels (scripts L8F3.py, L8F4.py), 5 Liu & Tromp, 6 maps, 7 uncertainty, 8 comparison. Master now loads
  placeins [section] and relaxed float fractions so figures stay in their sections; paper figures at 0.66 width.
* Lecture 20: paper-figure crops redone with wider bounds and automatic white-margin trimming (David saw
  clipping); Fig. 9 added (Cui et al. Fig. 14, P-wave cross-sections beneath hot spots in four models).
