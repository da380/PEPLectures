# PEP seismology lectures (Lectures 13–24): notes for maintaining the 2026 version

## Layout and build
* `lectures.tex` is the single master (preface, contents, `\lecture{N}{Title}{file}` for each lecture,
  `\problemset{N}{file}` for the two problem sets, references last). `lectureN.tex` (N = 1..12, i.e. Lecture N+12)
  and `problemN_body.tex` are bodies only: no preamble, no title. `./build.sh` runs pdflatex three times and
  produces `lectures.pdf`.
* `problemN.tex` are thin standalone wrappers around the bodies (10 pt); `solutionN.tex` are standalone (10 pt).
  Compile with `pdflatex` twice. Solutions are not part of the master document.
* `spare_questions.tex` holds questions removed from the sets (with solutions), kept as possible exam questions.
* Last year's sources and the latexdiff comparisons were removed when the revision was folded into the root
  (September 2026); they remain in the git history (last pre-revision commit `d38c2f9`).
* Labels are per lecture (`eq:N`, `fig:N`), made unique by the master's per-lecture prefix; never cross-reference
  an equation of another lecture by number, say "in Lecture 20" instead.
* `figures/` holds only files the documents reference: generated PDFs, TikZ output, and borrowed PNGs
  (Loma Prieta record, Benioff 1961 spectrum, GLAD-M35 crops, tomographic-model comparisons, Woodhouse &
  Dziewonski 1989). TikZ sources in `figsrc/` (compile with pdflatex, copy the PDF to `figures/`).
* `pepcode/` is a poetry package (`pepseis`) with every figure script (`scripts/L<n>F<m>.py`, `S1F1.py`,
  `make_figures.py`); see `pepcode/README.md`. Figure file names follow the *file* number (`L12F7` = lecture12.tex),
  not the printed figure number.

## House style
* Sentence-case titles; section headings `\section*`/`\subsection*`; non-examinable material marked
  "(non-examinable)" in the heading. Summary section "What you need to know and be able to do" as an
  `enumerate` with `\item[(i)]` labels.
* British spelling. "the Earth" for the planet, "earth model" for a model. P-wave, S-wave. "data are".
  "travel time" (noun, unhyphenated everywhere), "core--mantle boundary", "source--receiver", "first-order" as
  adjective. "Part IB" for the earlier course.
* En-dashes for name pairs and ranges (Piola--Kirchhoff, Cauchy--Green, 1859--1944, `13--24`); no spaced dashes
  as parentheses; footnote markers after punctuation; `Fig.~`, `eq.~(\ref{})`, `Lecture~20`; `e.g.\ `, `i.e.\ `;
  ``quotes''; expand acronyms at first use (PREM, CMB, ICB, LLSVP, MCMC).
* Displayed equations are punctuated as part of the sentence. Label every displayed equation; never hard-code
  numbers. Long equations: `multline` rather than an overfull line.
* Notation: `\dd` for differentials, `\ddns` for ordinary derivatives, `\ii`, `\ee`, `\bm{}` for bold Greek,
  `\bphi` for the motion, `\|\cdot\|` for norms, `\mathrm{}` for roman labels and units (`5\,\mathrm{km\,s^{-1}}`),
  bra-kets `\langle \mathbf{w}|P|\mathbf{u}\rangle` with italic operators. Surface radius of the model is `b`.
  The perturbation parameter is `s` (Lectures 14, 22, 23), except `\epsilon` in Lecture 24 where `\mathbf{s}`
  is the eigenfunction; `\epsilon` is also the ellipticity (Lecture 22) and the Levi-Civita symbol.
* Author's voice, first-person asides and footnotes are kept. Bold only where a term is first defined.

## Content decisions worth knowing (2026)
* Transverse isotropy is written in Love's `A, C, F, L, N` (Lecture 14, Problem Set 1); the elastic tensor is
  the second derivative of the elastic part of the strain energy. The auxiliary strain-energy functions are
  `U(x, C)` and `V(x, J)` (kept despite the clash with the stretch tensors).
* Lecture 17: caustics kept simple (envelope of the ray family); Euler's theorem proved in an appendix.
* Lecture 19: pixel basis stated explicitly; SVD/null space via the eigen-decomposition of `A^T A`; regularisation
  via the discrepancy principle only (no L-curve); Bayesian methods presented positively; function-space
  (pygeoinf) example at the end, not examinable beyond the idea.
* Lecture 20: adjoint derived with terminal conditions; time reversal introduced only as the device that turns
  the terminal-value problem into an initial-value one, with explicit time arguments in the kernels; the
  cross-correlation delay has the opposite sign convention to Lecture 19's delay time (remarked on).
* Lecture 22: "Rapidly rotating bodies" is descriptive and non-examinable (Maclaurin, Jacobi, pear-shaped
  branch, fission viewed sympathetically); boundary-perturbation derivation is an appendix.
* Lecture 23: completeness of the modes is stated as a fact ("complete basis for vector fields"), no
  qualification; stability left qualitative plus the statement that PREM calculations give no negative squared
  eigenfrequencies (undertones and stratification deliberately not discussed); Figs. 1–2 are real data.
* Lecture 24: full coupling = expansion in the spherical basis (the word Galerkin avoided); direct solution
  method described without preconditioner details, solved along a contour shifted into the complex plane
  (not at real frequencies); coupling matrix is not sparse; splitting functions and their limits for density
  (Akbarashrafi et al. 2018); adjoint kernels for spectra (Adourian et al. 2024). Students need the idea of full
  coupling only (summary item vi). No surface waves anywhere in the course.
* Problem sets: nine questions each; challenging parts carry a dagger (`\hardq` for a whole question, `\hardp`
  for a part); everything unmarked should be attempted. Set 2 Q1 is a guided Herglotz–Wiechert question
  (no numerical application), Q10 Coriolis splitting (toroidal `beta = 1/[l(l+1)]`; spheroidal
  `beta(0S2) = 0.40`, `beta(0S3) = 0.19` from PREM eigenfunctions, matching the observed 4.6 and 2.2 microHz
  singlet spacings).

## Data and computations behind the figures
* Seismograms: 2011 Tohoku (Mw 9.1, 2011-03-11 05:46:24), EarthScope FDSN via obspy, response removed to
  acceleration, cached as MiniSEED in `pepcode/data/seismograms/` (II.BFO 10 days for Lecture 23; CTAO, NNA, BFO
  20 days for the splitting figure). Spectra use a Hann window and zero padding.
* Normal modes: Alex Myhill's `mode_lab_2` (private, installed editable from `~/dev/mode_lab_2`), PREM at 1 s.
  Catalogue frequencies are elastic; where lines are compared with data they are corrected for physical
  dispersion using the modal Q from the kernels (`pepseis.modes`). The spheroidal catalogue is extended to
  degree 500 with the catalogue's own acceptance tests (cached in `pepcode/data/modes/`).
* Sensitivity kernels: toroidal analytic from the lecture's weak form; spheroidal by element-wise re-assembly
  and Rayleigh's principle, validated against direct re-solves (<1 %) and PREM's tabulated Q (0.1 %).
* Synthetic seismograms (Lecture 21): David's `yspec`, gravitation switches 2 (full) and 0 (none); parameter files
  and outputs in `pepcode/data/yspec/`.
* Ellipsoid sequences (Lecture 22) reproduce Chandrasekhar's classical numbers; Clairaut solution gives
  1/epsilon = 299.9 at the surface.

## Open items
* Reference list: Lichtenstein (1918) has no volume; David to check the Lichtenstein and Maitra & Al-Attar details.
* The 2025 worked-example sets and their verification scripts were deleted (not to be shared). Spherical Bessel
  functions are not part of the course.
