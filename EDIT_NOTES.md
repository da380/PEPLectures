# Revision of the seismology lectures (12–23) and problem sets: notes for review

## What is where

| Path | Contents |
|---|---|
| `revised/*.tex` | Proposed new versions of the 12 lectures, 2 problem sets and 2 solution sets (same file names as the originals; the originals in the repo root are untouched). Compiled PDFs alongside. |
| `revised/CHANGES_<file>.md` | Per-file change log: **A** corrections and rewordings (before → after, with reasons), **B** things left alone / questions for you, **C** the style pass. |
| `diffs/<file>_diff.pdf` | `latexdiff` mark-up of original → revised (blue underlined = added, red struck = removed; changed equations are shown whole). |
| `diffs/all_diffs.pdf` | All sixteen diffs concatenated in course order (172 pages). |
| `make_diffs.sh` | Regenerates the diffs (`./make_diffs.sh` or `./make_diffs.sh lecture5`). |
| `STYLE.md` | The house-style conventions that were applied, written so that any one of them can be reverted in bulk. |
| `examples/` | Worked-example sets (`examples1.tex` … `examples12.tex`) and `EXAMPLES_PLAN.md`; see the second half of this file. |
| `python/` | Verification scripts for the worked examples, figure scripts, and `FIGURE_PLAN.md` for regenerating lecture figures. |

To accept a revised file: `cp revised/lecture5.tex .` and recompile. Baseline note: `lecture1.tex` and
`lecture2.tex` already carried uncommitted edits in the working tree when this work started; those
working-tree versions (not `HEAD`) were taken as the starting point.

## Systematic changes (see `STYLE.md` for the full list)

* Sentence-case titles ("Finite elasticity", "Free oscillations", "Outline and motivation"); `\title` before `\begin{document}`; redundant `\usepackage{amsmath}` removed (it is in `mycommands.sty`).
* Every displayed equation labelled positionally (`eq:1`, `eq:2`, …) and **all hard-coded equation and figure numbers replaced by `\ref`s**. The hard-coded numbers were checked one by one; four were wrong (listed below). Figure labels are `fig:1`, `fig:2`, … in every file.
* Displayed equations punctuated as part of the sentence (roughly 150 missing commas/full stops added across the set).
* Notation: `\|\cdot\|` for norms; `\dd`/`\ddns` for roman differentials throughout (was mixed with `dt`, `dS`, `d\sigma`, `\frac{d}{dt}`); `\ii`, `\ee`; `\bm{}` for bold Greek and `\bphi` for the motion; roman `\mathrm{obs}`, `\mathrm{in}`, `\mathrm{diag}`; `\dots`/`\cdots` instead of `...` and `\cdot\cdot\cdot`; units as `2\,\mathrm{Hz}`, `5\,\mathrm{km\,s^{-1}}`; bold `\mathbf{x}` in all function arguments (Lecture 16 had many `u_i(x,t)`).
* Text: British spelling; ``quotes''; en-dashes for name pairs (Piola--Kirchhoff, Cauchy--Green, Euler--Lagrange, Sturm--Liouville) and ranges; footnote markers after punctuation; `e.g.\ `; "data are"; consistent compounds (wavefront, wavefield, waveform, travel time / travel-time curve, length scale, underdetermined, squared eigenfrequency, frequency domain); `Lecture~12`, `Fig.~`, `eq.~`.
* Summary sections all use `enumerate` with `(i)`, `(ii)`, …; non-examinable headings are "… (non-examinable)".

## Substantive corrections (please check these in particular)

Wrong cross-references in the originals (now fixed by the labels):
* Lecture 16: "eq.(36) has the trivial solution $T=\sigma$" referred to Euler's-theorem identity; it now points to $\mathrm{d}T/\mathrm{d}\sigma=1$.
* Lecture 17: "described in Lecture 17" → Lecture 16 (the ray-tracing ODEs); "Within the second problem set … horizontally stratified earth model" → first problem set (it is Problem 7 of Set 1). Lecture 16 item (iv) "final example sheet" → first problem set. Lecture 18 "Example sheet 3" → first problem set. Solutions 1: "Lecture 12" → 13 (chain rule for $\mathbf{C}$) and "Lecture 14" → 15 (perturbation theory).
* Lecture 18: the four-line `align` for $\partial J/\partial m_p$ carried four equation numbers, so "eq.(31)" pointed at the wrong display; the block is now one equation and the reference is to the regularised least-squares solution.
* Lecture 23: "Using eq.(8)" (the $\mathbf{A}_{lm}$ orthogonality) → the $\mathbf{C}_{lm}$ orthogonality relation, which is the one used.

Mathematical slips corrected:
* Lecture 14, polynomial example: "$m>1$" → "$m\ge 1$" (collinear points are fitted by a degree-one polynomial); "unique if and only if $m+1=n$" → "$m+1\le n$" (an overdetermined but consistent Vandermonde system has a unique solution).
* Lecture 16: Fourier-integral display was missing $\omega$ in the exponent and the sentence around it was a run-on; "travel time depends linearly on $\mathbf{p}$" → "on $\mathbf{x}$"; "divided by a factor of $1/(\mathrm{i}\omega)$" → "multiplied by a further factor of".
* Lecture 20: $\gamma_i=R_{ij}\overline{\gamma}_i$ → $R_{ij}\overline{\gamma}_j$; $\overline{v}=\sqrt{\overline{\mu}/\rho}$ → $\sqrt{\overline{\mu}/\overline{\rho}}$; "perturbation to the referential gravitational potential" → "acceleration"; $\partial\zeta_0/\partial x_i\partial x_j$ → second derivative; the dummy index in $R_{kj}\dot R_{kl}=\epsilon_{jkl}\Omega_k$ renamed ($\epsilon_{jml}\Omega_m$).
* Lecture 22: "pointwise stability holds in an isotropic medium if and only if P- and S-waves have positive phase speeds with the P-wave faster" replaced by "if and only if the bulk and shear moduli are both positive, and this in turn implies … $\alpha>\beta$" (positive wave speeds is the weaker strong-ellipticity condition; e.g. $\lambda=-0.9\mu$ gives $\alpha>\beta>0$ but $\kappa<0$). Convolution kernel $\tilde h(\omega'-\omega)$ → $\tilde h(\omega-\omega')$. Orthogonality stated for $\omega_1^2\ne\omega_2^2$. Symbol clashes removed ($k\to c_0$ for the stability constant, $A_{ij}\to B_{ij}$ for the infinitesimal rotation).
* Lecture 15: the isotropic eigenvector displays used an undefined $\hat{\mathbf{a}}$; now $\bm{\Gamma}\hat{\mathbf p}=\frac{\lambda+2\mu}{\rho}\hat{\mathbf p}$ and $\bm{\Gamma}\mathbf a=\frac{\mu}{\rho}\mathbf a$ for $\mathbf a\perp\hat{\mathbf p}$. The scalar in $\bm\Gamma(\lambda\mathbf p)=\lambda^2\bm\Gamma(\mathbf p)$ renamed $\gamma$ to avoid the Lamé parameter.
* Lecture 14: $n_l$ → $\hat n_l$ in the general-fault stress glut (the text calls the normal $\hat{\mathbf n}$).
* Lecture 18: "null space not empty" → "non-trivial"; summary item (ii) had the standard form as $\mathbf{Am}=\mathbf d+\mathbf e$, now $\mathbf d=\mathbf{Am}+\mathbf e$.
* Solutions 1, Q3: the identity $\partial C_{kl}/\partial F_{ij}$ had its indices swapped (now $\delta_{jk}F_{il}+\delta_{jl}F_{ik}$, as in Lecture 13), and the second-derivative line had a dangling index $q$ (now $\partial C_{pl}$).
* Solutions 2: Q2(b) $(\mathbf I+\mu\mathbf A^T\mathbf A)\mathbf m=2\mu\mathbf A^T\mathbf d$ → $\mu\mathbf A^T\mathbf d$; Q3 the shear-modulus kernel was a factor 2 too large (now $K_\mu=\tfrac12(\partial_ju_i+\partial_iu_j)(\partial_ju'_i+\partial_iu'_j)$); Q4 the first display had the wrong overall sign ($+\tfrac12 G\int\dots$); Q5 the azimuth was written $\phi$ (clashing with the potential) → $\varphi$.
* Problem Set 2, Q3: the standard deviations were $\sigma_i$ while the load is $\sigma$; the standard deviations are now $\varepsilon_i$ (also in the solutions). Q7: bra-ket operators written $P,W,H$ as in Lecture 22.
* Lecture 17: the surface radius was $a$ here but $b$ in Lectures 21 and 23; now $b$ throughout. The PREM paragraph now refers to the figure and expands the acronym.
* Lecture 16: $H_k(\mathbf x,\mathbf p)=\tfrac12\|\mathbf p\|^2c_k(\mathbf x,\hat{\mathbf p})^2$ is now written with independent phase-space arguments (the original had $\mathbf p(\mathbf x)$ inside the definition of a function on phase space).

## Decisions needed from you

1. **Lecture 21 summary.** Item (ii) mentioned "the two compatibility conditions" and item (vi) "minimum stress fields and their relation to slow viscous flow"; neither topic is in the lecture. I reworded (ii) and removed (vi). If that material was meant to be in the lecture (or a problem set), it needs to be restored instead.
2. **Lecture 23 summary item (ii)** pointed to a toroidal-mode/spherical-Bessel-function question "in the second problem set" that does not exist. It now points to the worked examples for Lecture 23 (Example 1 there is exactly this calculation). Alternatively add such a question to Problem Set 2.
3. **Lecture 13**: "the 1B dynamics course" — should this read "Part IB"?
4. **Lecture 14**: the wording "To rectify things, we will not produce a complete dynamical theory…" reads oddly; left as is.
5. **Lecture 21 figure** and **Lecture 17 Fig. 3** are never cited in the text; no citations were added because that would be new text.
6. The en-dash convention for name pairs (Piola--Kirchhoff etc.) and "data are" are the two style choices most likely to be a matter of taste; both are one `sed` to undo.
7. Lecture 23 uses $W'(r)$ for the toroidal test function (a prime that is not a derivative) alongside $\mathrm{d}W'/\mathrm{d}r$; a different symbol would be clearer. Left as is.

More minor observations are in section B of each `revised/CHANGES_*.md`.

## Worked examples (task 2)

Twelve sets, `examples/examples1.tex` … `examples12.tex` (one per lecture, 7–11 pages each, 112 pages
in `examples/all_examples.pdf`). Each set has five or six boxed problems with full solutions; the
specification is in `examples/EXAMPLES_PLAN.md`. Every numerical or closed-form claim is checked by
`python/verify_examplesN.py` (all checks pass; run e.g. `python3 python/verify_examples12.py`), and
eight sets carry a generated figure from `python/exN_*.py` (in `examples/figures_ex/`), including the
polar-decomposition picture for the $\mathbf F$ of Example 2 of Lecture 12.

Points the authors of the sets flagged for you to look at by eye (everything else was derived and
checked numerically):
* Lecture 12, Ex. 5 closing remark: a frame-indifferent $W$ quadratic in $\mathbf F$ must be linear in $\mathbf C$ and so cannot be stress-free at $\mathbf F=\mathbf I$ (argument is correct).
* Lecture 13, Ex. 3: the symmetries of the transversely isotropic energy are the rotations with $\mathbf Q\hat{\bm\nu}=\pm\hat{\bm\nu}$ (not only rotations about the axis); Ex. 2 introduces the Cauchy stress via Nanson's formula.
* Lecture 14, Ex. 4 uses the Hanks–Kanamori $M_w$ formula and the statement that far-field displacement is proportional to $\dot M_0$; Ex. 5 quotes Vandermonde condition numbers.
* Lecture 15, Ex. 5 introduces "energy velocity" (normal to the slowness surface).
* Lecture 16, Ex. 4(b) quotes the heterogeneous transport equation $\nabla\cdot(\rho\alpha^2A_0^2\nabla T)=0$ without derivation.
* Lecture 17, Ex. 2 quotes the Herglotz–Wiechert formula; Ex. 5 gives the shadow zone of a two-layer model as 114°–142° (the 103° edge of the real Earth needs the mantle gradient).
* Lecture 18, Ex. 4 defines a resolution matrix in one line (this is also Problem 1 of Set 2).
* Lecture 19, Ex. 4: with the lecture's definition $C(\tau)=\int s^{\rm obs}(t-\tau)s(t)\,\mathrm dt$ the adjoint traction is $-w\dot s/\int\dot s^2$; sign conventions are stated and checked numerically.
* Lecture 20, Ex. 1 quotes $10^{18}$–$10^{19}$ J for the largest earthquakes.
* Lecture 21, Ex. 4 has a short historical footnote (Newton 1/230, Clairaut 1743) and quotes PREM gravity values from memory (9.82 and 10.68 m s$^{-2}$).
* Lecture 22, Ex. 4 uses an observed $_0S_2$ singlet spacing of about 4.6 μHz and $Q\approx510$.
* Lecture 23, Ex. 1 is the toroidal-mode/spherical-Bessel calculation now cited in the lecture summary ($_0T_2$ = 42.3 min for a homogeneous sphere with $\beta=6.3$ km s$^{-1}$, observed 44.2 min); Ex. 3 states the Coriolis matrix element as $\langle km'|W|km\rangle=-2\mathrm i\Omega\beta_k m\,\delta_{mm'}$ (sign derived from the lecture's conventions, verified numerically) so that $\delta\omega_m=+m\Omega\beta_k$; Ex. 4 uses illustrative mineral-physics derivatives; Ex. 5 mentions the Slichter mode $_1S_1$.

## Figures in Python

`python/FIGURE_PLAN.md` classifies the 28 lecture figures: eight can be regenerated with short
scripts (polar decomposition, slowness surfaces, PREM profiles, Clairaut ellipticity, ray fans
with a caustic, Earth ray paths), two need a small mode/wave solver, and the rest are
observational or hand-drawn and should stay. The polar-decomposition, ray-fan, shadow-zone and
toroidal-mode scripts written for the examples are the starting points for the first group.
