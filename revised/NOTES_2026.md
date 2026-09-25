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
* Preface added to lectures.tex (TiS availability, errors to da380@cam.ac.uk, problem sets to be appended to the
  document, solutions on the TiS). TO DO at tidy-up: append problem sets 1 and 2 to the end of lectures.tex
  (David, 25 Sept; supersedes "lectures only").
* Lecture 21 (lecture9.tex): factor 2 dropped from the centrifugal scaling; "final problem set" -> "second";
  new Fig. 1 (PREM g and p from pepseis.prem, forward reference to Lecture 24 for the density) and Fig. 2
  (yspec synthetics with/without self-gravitation, Mw 8, Delta = 90, 0.2–50 mHz, six hours, no attenuation, stacked band-passes, 10-min start taper; inputs and outputs in
  pepcode/data/yspec, generated by scripts/run_yspec.py with YSPEC=<binary>; yspec built from ~/dev/YSpec).
  Reference Al-Attar & Woodhouse (2008) added.
* Lecture 20 Fig. 1 reduced to 0.78 textwidth so it sits on the first page.
* Lecture 21 Fig. 2 now uses acceleration seismograms (yspec output 2), broadband and 1–2 mHz band-passed
  columns (David). Lecture title blocks compacted in lectures.tex (\@maketitle redefined: no empty author line,
  smaller skips) so that Lecture 20's Fig. 1 fits on its first page; that figure is a non-floating \captionof.
  Preface trimmed (solutions "available on the TiS", no editorial remark).
* Functional-derivative notation audit (David): Lecture 13 now introduces $\delta\mathcal S=\langle D\mathcal S(\bphi)|\delta\bphi\rangle$
  and the term "functional derivative" where the first variation is defined, so that the back-reference in
  Lecture 19 is correct; Lectures 19–21 use $D(\cdot)$ with the bra-ket pairing consistently; Lecture 21's
  $\delta\mathcal L/\delta\varphi_i$ is a locally defined shorthand (also used in Set 2 Q4); Lecture 24's
  $\delta\gamma/\delta\beta$ is a ratio of perturbations, not a functional derivative. No other changes needed.
* mode_lab_2 (Alex Myhill's SEM normal-mode code, ~/dev/mode_lab_2, installed editable into the pepcode venv)
  wrapped in `pepseis/modes.py`: PremModes (toroidal/spheroidal/radial solves on a p=5, 100 km mesh, <1 s per
  degree; eigenfunctions sampled on radius), toroidal kernels K_mu, K_rho from the lecture's toroidal weak
  form (verified: int K_mu mu dr = omega/2, int K_rho rho dr = -omega/2), and access to the reviewed PREM
  catalogues (toroidal 150k modes, spheroidal 62k accepted). Draft figures (not yet in the text, awaiting
  David): L12F1/L12F2 dispersion diagrams regenerated; L12F4 PREM density from prem.py; L12F6 eigenfunctions
  (0T2, 1T2, 0T10, 0T40; 0S2, 1S2, 0S10, 0S0); L12F7 toroidal kernels. (An L11F3 yspec-synthetic spectrum was drafted and dropped: David does not
  want comparisons with 1D synthetics here, 3D effects being too large; data/yspec/long.* kept.)
* Lecture 22 (lecture10.tex): Fig. 1 regenerated (PREM density + 1/epsilon from a new Clairaut solver,
  pepseis/clairaut.py: surface 1/eps = 299.9, centre 415) and now cited, with a sentence on the hydrostatic
  (1/299.9) versus observed (1/298.3) flattening. Lecture titles now left-aligned (lectures.tex).
* Lecture 22: David confirmed the removal of the old summary items (ii)/(vi). Uniqueness footnote now cites
  Lichtenstein (1918, 1933) for the general theorem and mentions the quasi-linear-elliptic symmetry viewpoint
  (David: Gårding?). Two Lichtenstein references added — DAVID TO CHECK the 1918 details.
* Lecture 22: figures renumbered (Fig. 1 TikZ cartoon of coincident level surfaces, cited in the hydrostatic
  section; Fig. 2 PREM density + Clairaut). New descriptive section "Rapidly rotating bodies" (David: cartoons,
  Haumea, Jacobi/bifurcation, Poincaré–Darwin fission, qualitative): Fig. 3 Maclaurin/Jacobi bifurcation diagram
  from pepseis/ellipsoids.py (index symbols; verified against Chandrasekhar: bifurcation e=0.8127,
  c/a=0.5827, Omega^2/piG rho=0.3742; pear point a2/a1=0.4322, 0.2840; Maclaurin max 0.4493) with section
  shapes; Fig. 4 Maclaurin, Jacobi (a2/a1=0.73) and Haumea (Ortiz et al. 2017 axes) rendered. Text notes that
  Haumea's axes match the Jacobi ellipsoid but Omega^2/piG rho = 0.50 exceeds the homogeneous maximum
  (differentiated interior). References Ortiz et al. 2017 and Poincaré 1885 added. DAVID TO CHECK the
  history (Liapounov/Cartan instability, Darwin fission) and the Haumea interpretation.
* Lecture 22: fission paragraph rebalanced at David's request (giant-impact isotope problem, angular-momentum
  objection to fission needs an ad hoc removal mechanism, renewed interest); Fig. 3 gained schematic pear and
  dumbbell sketches.
* Lecture 23 (lecture11.tex): David: the old Figs 1–2 were synthetic ("not data"); replaced with real
  records. obspy added as a dependency; `scripts/fetch_modes_data.py` downloads ten days of LH? data after the
  2011 Tohoku earthquake (IU ANMO/KIP/CTAO, II BFO; EarthScope FDSN), removes the response to acceleration and
  caches MiniSEED in `pepcode/data/seismograms/` (45 MB; KIP has a glitch mid-record and is unused).
  `pepseis/seismograms.py` gives access, R/T rotation and Hann-windowed spectra. Fig. 1 (L11F1.py): BFO
  vertical, first 12 h unfiltered + 6 h–10 d band-passed 0.3–1 mHz (gravest modes ring for ~4 days).
  Fig. 2 (L11F2.py): amplitude spectra 0.25–3 mHz, vertical (10 d) with PREM 0S_l and 0S_0 marked, transverse
  (first 3 d, horizontals noisier) with 0T_l marked; the catalogue lines sit on the observed peaks, and 0S2/0S3
  splitting is visible (caption points forward to Lecture 24). No 1D synthetic comparison (David). Completeness
  of the modes left as an assumption, no remark added (David: physics not maths).
* Lecture 23: eqs. (14) and (15) (Coriolis and centrifugal self-adjointness chains) split over two lines with `multline` (David).
* Lecture 24 dispersion diagrams (David, 2026-09-25): extended to 50 mHz and degree 500, radial modes no
  longer distinguished from spheroidal ones. Toroidal from the reviewed catalogue (reaches 50.9 mHz at l=500).
  The reviewed spheroidal catalogue stops at l=200 (~20 mHz on the fundamental), so `pepseis.modes.
  spheroidal_extension` solves l=201–500 with mode_lab_2 (order 5, 50 km mesh) and applies the catalogue's
  own acceptance tests (`diagnose_spectrum`: p/h-refinement Rayleigh error < 1%, ocean energy fraction < 0.5,
  numerical-zero states dropped first); validated at l=200 (26/26 states, agreement to 6e-6). Results cached
  in `pepcode/data/modes/spheroidal_l201-500_50mHz.npz` (~13 min to compute). Captions updated.
* Lecture 23 stability: David's decision — lecture stays qualitative, plus one sentence that PREM calculations return no
  negative squared eigenfrequencies (checked with mode_lab_2: spheroidal l≤30 to 20 mHz and toroidal l≤30, none negative;
  the many zero-frequency fluid-core states are deliberately NOT mentioned — David: undertones/stratification are out of scope).
  The homogeneous-planet compression argument is now Set 2 Q8 (binding energy, hydrostatic pressure, virial relation,
  second-order energy → κ > (4/3) p̄ = (16π/45)Gρ²a² ≈ 92 GPa, sound speed > sqrt(4ga/15) ≈ 4.1 km/s; γ > 4/3 analogue).
  This covers the "binding energy" item on David's Set 2 wish list.
* Lecture 20 (lecture8.tex), David 2026-09-25: time reversal ADDED after all, purely as the device converting the terminal
  value problem into an initial value one (eqs. for u†(x,t)=u'(x,T−t), h†, initial conditions; kernels restated with explicit
  time arguments pairing u(x,t) with u†(x,T−t); u† now called the adjoint wavefield, footnote on naming). Reverses the earlier
  "no time reversal" consistency note. Fig. 3 gained a fourth column (instantaneous integrand −2ρα ∇·u ∇·u†) and the caption
  recalls the integral expression; `KernelExperiment.integrand(i)` added.
* Lecture 13 (lecture1.tex): sentence after eq. (1) that the freedom of choice of reference body can be put to use
  (aspherical model rewritten on a spherical reference body). David first asked for an Al-Attar & Crawford (2016)
  citation and borrowed figure, then withdrew both; no reference, no figure.
* Lecture 24 dispersion diagrams regenerated to 50 mHz / l ≤ 500 (L12F1.pdf, L12F2.pdf now in the text; 2294 accepted
  extension states cached). Both figures now included as PDFs.
* Lecture 24 (lecture12.tex), 2026-09-25 afternoon, David's decisions applied:
  - Symbols: toroidal test function W' → w(r); bulk sound speed γ → v_φ. Five-moduli footnote cites Love's A,C,F,L,N
    (Lecture 14) and states isotropy is assumed. Fig. 5 (S-model depth slices) now credited to Ritsema et al. (2011);
    reference added.
  - New figures and text: L12F6 (representative eigenfunctions; paragraph on nodes/overtones, surface confinement —
    NB no surface waves in the course, David), L12F8 ("curious" modes: Slichter 1S1 on a 100 km mesh [50 km mesh fails
    the degree-one translation test], CMB Stoneley 2S16, inner-core 11S2, l=25 tsunami mode zoomed to the top 16 km; all
    picked by kinetic-energy partition, see scripts/L12F8.py), L12F7 (kernels: toroidal analytic from eq. 18, spheroidal
    by element-wise re-assembly + Rayleigh's principle in `pepseis.modes.spheroidal_kernels`, validated against direct
    re-solves to <1% and against PREM's tabulated Q and dispersion-corrected frequencies to 0.1% via
    `modal_q_and_dispersion`), L12F9 (0S2/0S3 singlets from 20-day Tohoku records at CTAO, NNA, BFO; Hann window,
    8x zero padding; PREM lines are dispersion-corrected — the elastic catalogue is 0.5% high at these periods).
    Kernel paragraph: Rayleigh's principle (footnote), toroidal K_mu, K_rho written out, spheroidal only described.
  - Closing section expanded: full coupling = expansion in the spherical basis (word "Galerkin" avoided), self-coupling,
    splitting functions and their use/limits for density, Akbarashrafi et al. (2018) cited (reference added), Stoneley/CMB
    modes flagged as where the problem is worst (David: the paper understates this).
  - Figure files no longer match figure numbers (L12F6 is Fig. 3, L12F8 Fig. 4, L12F3 Fig. 5, L12F9 Fig. 6, L12F7 Fig. 7,
    L12F4 Fig. 8, L12F5 Fig. 9). Cosmetic; renumber files at the tidy-up stage if wanted.
  - Possible refinement: Lecture 23 Fig. 2 mode lines use elastic (1 s) catalogue frequencies, 0.3–0.5% high; a
    dispersion correction per mode is now available (needs toroidal Q too). Not done.
  - Data: `pepcode/data/seismograms/tohoku20_*.mseed` (8 stations, 20 days; MAJO and TUC have glitches, unused).
* Lecture 19 (lecture7.tex): Fig. 1 (Woodhouse & Dziewonski 1989 image) reduced to 0.5\textwidth (David).
* Set 2 homogeneous-sphere toroidal question: David — probably no space; discuss at the Set 2 stage.
* Lecture 24 curious-modes figure: tsunami panel dropped (David: looked off — the ocean is a single 3 km element and the
  sampled V ramped linearly; not investigated further) and replaced by the ICB Stoneley mode 4S8 (accepted catalogue
  mode, 94% of kinetic energy within 300 km of the ICB). The OC-only states near the ICB found by the raw solver
  (2.5–3.2 mHz at l=2 etc.) are NOT accepted catalogue modes — spurious fluid-core states; do not use.
* Lecture 19 Fig. 1 made non-floating (\captionof) so it sits on the first page, at 0.42\textwidth.
* Lecture 24 closing section (David): added the direct solution method (Al-Attar, Woodhouse & Deuss 2012; reference
  added): forced frequency-domain system in the spherical basis (eqs. 38–39 in the current numbering), iterative solution
  described only as diagonally dominant away from the block of modes near the target frequency and "solved efficiently by
  modern iterative methods" (David: the coupling matrix is NOT sparse in general; the IDSM preconditioner inverts the
  target block directly and is diagonal elsewhere, but students don't know what a preconditioner is — keep it vague). Eigenfunction paragraph + Figs. 3–4 moved to the end of the spherical-model section (after spheroidal
  modes are introduced). Figure order now: 1–2 dispersion, 3 eigenfunctions, 4 curious modes, 5 Benioff, 6 splitting,
  7 kernels, 8 density, 9 tomography. References currently spill one line onto a final page (146); will change once the
  problem sets are appended.
  Also (David): solutions are computed on a contour shifted into the complex plane (Bromwich-like), not at real frequencies;
  the time-domain result comes from an inverse FFT times a growing exponential over a targeted interval. Text now says only
  "along a line shifted slightly into the complex plane ... amounts to computing a seismogram of finite length".
* Lecture 24 final paragraph (David): link to adjoint kernels for spectra built on the DSM (Adourian, Dursun, Lau &
  Al-Attar 2024, GJI 238, 257–271 — David remembered it as 2025; reference added) and to waveform-inversion methods moving
  to free oscillations.
* Lecture 24 summary list: item (vi) added — know the idea of full coupling (expansion in the spherical basis) and why first-order theory fails for density; details of the calculations, splitting functions and adjoint spectra NOT needed (David).
* Problem Set 2 extended to ten questions (2026-09-25): new Q1 guided Herglotz–Wiechert (tau(q), dT/dDelta = q, Abel
  inversion with the supplied identity, cosh^-1 form, discussion of LVZs and triplications); Q7 (point source) gained
  part (a) on the point-source approximation from the Taylor expansion of the eigenfunction strain (L << v/f; 100 km
  fault, 5 km/s -> fine below ~10 mHz); new Q10 Coriolis splitting (omega_m = omega_k + m Omega beta from eq. 28 of
  Lecture 24, consistency with Q8, toroidal beta = 1/[l(l+1)] derived with C_lm, spheroidal beta quoted:
  0S2 0.40, 0S3 0.19 computed from mode_lab_2 eigenfunctions with beta = int rho(2UV+V^2) r^2 dr / int rho(U^2+l(l+1)V^2)
  r^2 dr in the lecture's unnormalised-B convention; predicted spacings 4.6 / 2.2 microHz match the observed 0S2/0S3
  singlets of Lecture 24 Fig. 6). Solutions written for all three. David: NO numerical Herglotz–Wiechert application
  (a PREM check was coded and worked to 1e-5 away from discontinuities, but was dropped — "can't work" in practice /
  not needed; script S2F1.py deleted). Adjoint question kept as is. Homogeneous-sphere toroidal question not added.
* Problem sets incorporated into lectures.tex via `\problemset{n}{problemn_body}` (same title treatment as lectures,
  label prefix Pn, TOC entries) after Lecture 24 and before the References. `problem1.tex`/`problem2.tex` are now thin
  standalone wrappers that \input the `problemN_body.tex` files (single source of truth); solutions stay standalone.
* Problem sets trimmed and marked (David, 2026-09-25): Set 1 Q9 (rays in a depth-varying half-space) and Set 2 Q6 (Poisson
  boundary perturbation) DROPPED — both kept with solutions in `revised/spare_questions.tex` as possible exam questions
  (David: the half-space one is "not a bad exam question to have in the pocket"). The Poisson boundary-condition
  derivation is now a non-examinable appendix to Lecture 22 (lecture10.tex). Nine questions per set.
  Difficulty marked at PART level with a dagger (`\hardq`/`\hardp` macros defined by \providecommand at the top of each
  body file; note at the head of each sheet): Set 1 — Q5(b),(c) and Q9(b),(c); Set 2 — Q1(c), Q4 (whole question), Q9(b).
  Solutions carry no daggers. Set 2 Q1(a) gained a hint about Leibniz's rule (the turning-point limit depends on q).
  Set 2 cross-reference in Q9 now points to question 7 (rotating eigenvalue).
* FULL REVIEW PASS (2026-09-25, evening): seven read-only reviewers (two lectures each + problem sets) checked text,
  figures and maths; five fixers then applied all editorial [A] and figure/layout [C] items plus the safe [B]
  clarifications. Everything rebuilt: master 159 pp, 0 errors, 0 undefined refs; standalone problem/solution docs 0
  errors/0 overfull. Highlights of what changed: L24 "toroidal modes vanish at the CMB" corrected to vanishing shear
  traction; overtone/node statement qualified for spheroidal modes; 0S0 density-kernel claim dropped (only 0S2 changes
  sign); complex-contour sentence reworded; L20 Fig. 1 caption cross-ref fixed, receiver delta made a surface delta,
  "iff" -> "requires", CC delay sign-convention remark added, model-dimension symbols unified; L21 gravity integrals use
  y_i, yspec comparison wording made precise (switch 0 = all gravitational effects off); L22 mantle-flow speed 1e-9 m/s
  and 1e6 Pa; L18/19 many wording fixes, null-space vector renamed m_null, stale problem-set pointer fixed, Bayesian
  cross-ref fixed; L19 toy-tomography figures regenerated (thinner rays, localised singular vector j=278, central spike
  block, colourbars); L22/L24 figures regenerated at printable font sizes; references: order, van Heijst, Bozdağ/
  Dahlen & Tromp/Clairaut/Poincaré/Reid now cited in text. REVERTED one fixer change: the auxiliary strain-energy
  functions stay U(x,C) and V(x,J) (David's earlier decision), not W-tilde/W-hat.
  LEFT FOR DAVID (not changed): (1) L18 Oldham/Gutenberg attribution of the P shadow zone; (2) L23 "complete basis for
  an arbitrary vector field" vs "field satisfying the boundary conditions" (lines ~197 vs ~213) — which wording;
  (3) L24: perturbation parameter s next to eigenfunction bold s (rename to epsilon?), "co-ordinates" spelling;
  (4) Lichtenstein 1918 volume; (5) L11F2 (Lecture 23 Fig. 2) mode lines are elastic 1 s-reference frequencies, 0.3–0.5%
  high (dispersion correction available). Standalone solution1.pdf has a \newpage before the S1F1 figure so figure and
  part (c) share the last page.
* Follow-ups after the review (David's decisions): L24 node counting — spheroidal modes "no such simple rule" with a
  footnote on Woodhouse's minor-vector generalisation (Woodhouse 1988 reference added; Al-Attar & Woodhouse 2008); L24
  perturbation parameter renamed s -> epsilon (eqs 21–27); NOT renamed elsewhere because epsilon is the ellipticity in L22
  and the Levi-Civita symbol in L21/L23 — s stays there. L18 Oldham/Gutenberg: Oldham (1906) argued for the core from late
  arrivals at large distances, Gutenberg (1914) explained the shadow zone and found the CMB depth (Gutenberg 1914 reference
  added). L23 completeness now stated without the boundary-condition qualification ("complete basis for vector fields";
  David: H1-completeness needs no reference to boundary conditions, just say complete). Set 2 Q4 footnote kept and the
  solution gained a paragraph on well-posedness of pure-traction static problems (net force/torque, rigid motions fixed by
  constraints, kernels unaffected). L23 Fig. 2 mode lines now dispersion-corrected (`fundamental_frequencies_corrected`,
  `toroidal_q_and_dispersion`: 0T2 Q 249 vs PREM 250.5; 0T10 Q 173 vs 184 — spline-derivative kernel near the crust —
  corrected frequencies within 0.25%, invisible at plot scale); cached in pepcode/data/modes/fundamentals_corrected_l20.npz.
