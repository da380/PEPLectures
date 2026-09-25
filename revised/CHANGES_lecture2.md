# Changes to `lecture2.tex` (Lecture 13: Constitutive theory)

## A. Corrections and rewordings

**Preamble**
* `\documentclass [a4paper,12pt]{paper}` → `\documentclass[a4paper,12pt]{paper}`; `\title`/`\author` moved above `\begin{document}` (house skeleton).

**Proof of the polar decomposition theorem**
* Subsection title `... theorem - NON-EXAMINABLE` → `... theorem (non-examinable)`.
* Footnote on the square root: `\text{diag}` → `\mathrm{diag}` (twice); footnote marker moved after the comma (`defined,\footnote{...} and we set`); inside the footnote, `\mathbf{Q}^T$ which is well-defined` → `\mathbf{Q}^T$, which is well-defined` (non-restrictive clause).
* "For the second identity, we merely note that" → "For the second form of the decomposition, we merely note that" (the "first identity" earlier in the paragraph is $\mathbf{R}^T\mathbf{R}=\mathbf{I}$, so "second identity" was ambiguous; the sentence is about $\mathbf{F}=\mathbf{V}\mathbf{R}$).

**Implications for the strain energy function**
* "only through the symmetric matrix $\mathbf{U}$, or equivalently on $\mathbf{C}$" → "... or equivalently through $\mathbf{C}$" (parallel construction).
* Footnote on the left Cauchy--Green tensor: bold extended to cover "left Cauchy--Green deformation tensor"; "tensor which is useful" → "tensor, which is useful"; footnote marker moved after the full stop.
* Footnote on Noether: "For those that know" → "For those who know"; marker moved after the full stop.
* Hard-coded references `eq.(1)`, `eq. (9)`, `eq. (15)` → `eq.~(\ref{eq:1})`, `eq.~(\ref{eq:9})`, `eq.~(\ref{eq:15})`. All pointed at the correct equations; no mapping corrections were needed.

**The equilibrium equations**
* `eq.(19)`, `eq.(20)` → `eq.~(\ref{eq:19})`, `eq.~(\ref{eq:20})` (both already correct).

**Linearised equations of motion**
* "Working always to first-order in $s$" → "to first order in $s$" (adverbial use, per style guide).
* "linearisation of the \textbf{finite-strain} $\mathbf{E} = \dots$" → "linearisation of the \textbf{finite strain tensor} $\mathbf{E} = \dots$" (missing noun).
* "in terms of $\mathbf{E}$ and not $\mathbf{C}$ but there is little gained" → "in terms of $\mathbf{E}$ rather than $\mathbf{C}$, but there is little gained".
* `Eqs.(36) and (37)` → `Eqs.~(\ref{eq:36}) and (\ref{eq:37})` (correct as they stood).
* "core-mantle boundary" → "core--mantle boundary" (name-pair en-dash).

**Material symmetries / Isotropic materials**
* `eq. (38)` → `eq.~(\ref{eq:38})`; `eq.(43)` → `eq.~(\ref{eq:43})`; `eq.(39)` → `eq.~(\ref{eq:39})` (all correct as they stood).
* "To simplify notations" → "To simplify notation".
* "with $\mathbf{C} = \mathbf{F}^{T}\mathbf{F}$. And hence for an isotropic material" → "with $\mathbf{C} = \mathbf{F}^{T}\mathbf{F}$, and hence for an isotropic material".
* "the strain energy depends only on $\mathbf{C}$ through its eigenvalues" → "depends on $\mathbf{C}$ only through its eigenvalues" (what is meant is that the dependence on $\mathbf{C}$ is only via its eigenvalues).
* Isotropic homogeneous equation of motion (eq. 43): right-hand side `= 0` → `= \mathbf{0}` (vector equation; matches `\mathbf{T}^0(\mathbf{x}) = \mathbf{0}` used earlier in the file). **Mathematical notation correction -- please check.**
* `(e.g. isotropy, ...)` → `(e.g.\ isotropy, ...)`.

**Transversely isotropic materials / Elastic fluids**
* `\hat{\boldsymbol{\nu}}` → `\hat{\bm{\nu}}` (three occurrences).
* "$\gamma, \xi$, and $\zeta$" → "$\gamma$, $\xi$ and $\zeta$".
* Noll footnote: marker moved after the full stop; "(1925-2017)" → "(1925--2017)".

**Punctuation of displayed equations**: every displayed equation was checked; all already ended in the comma/full stop the sentence requires (eqs. 8, 16 and 32 sit mid-sentence with no punctuation, correctly), so no changes were needed here.

## B. Left unchanged / questions for the author

* "This is the form of the equations of elasticity covered in the 1B dynamics course." -- should this read "Part IB dynamics course" (or "the Part IB Dynamics and Relativity course")? Left as is.
* "reduce the number of independent components from 45 to 21" -- correct (hyperelastic symmetry alone gives 45), left.
* "Let $\bphi^0$ be an equilibrium configuration of a body." (equilibrium section) repeats the definition given two sentences earlier; harmless, left.
* "for some given vector $\mathbf{u}^{0}$" -- strictly a vector field; left as written.
* "In the first problem set you will show that this tensor has additional symmetries" -- fine, but the phrase "in the first problem set" also appears earlier for the angular-momentum result; both are presumably in the first set. Not changed.
* In the transversely isotropic elastic tensor (eq. 45) the coefficients $8\gamma$, $4\xi$, $-\zeta$ follow one particular convention; I did not check them against a source and left them untouched.
* The footnote on Noll states there is no proper subgroup intermediate between $\mathrm{SO}(3)$ and $\mathrm{SL}(3)$ -- this is Noll's result for the unimodular group; left as is.

## C. Style pass

Applied: labels `eq:1`--`eq:48` (48 numbered equations, `split` counted as one) and `Fig.~\ref{fig:1}` already present; all hard-coded `eq.(n)`/`Eqs.(n)` replaced by `eq.~(\ref{eq:n})`; Piola--Kirchhoff (7), Cauchy--Green (2), core--mantle, year range en-dashes; `\boldsymbol` → `\bm`; `\text{diag}` → `\mathrm{diag}`; footnote markers after punctuation; `e.g.\ `; subsection non-examinable tag; document skeleton. No place where the style guide could not be applied cleanly. Line breaks and indentation of untouched lines preserved (each `\label` is inserted on its own line before `\end{equation}`, as in `revised/lecture1.tex`).

## D. Second pass (25 September 2026)

* Renumbered: now "Lecture 14: Constitutive theory" (title set in `lectures.tex`; this file is body only).
* "the 1B dynamics course" → "the Part IB dynamics course".
* Fig. 1 regenerated in matplotlib from an explicit deformation gradient (`figsrc/L2F1.py` → `figures/L2F1.pdf`).
* $U$ (stretch) and $U(\mathbf{x},\mathbf{C})$ (auxiliary strain energy) deliberately left as they are.
* Repeated sentence in the equilibrium section cut.
* Eq. (45): transversely isotropic elastic tensor now in Love's $A,C,F,L,N$ notation, with the component
  meaning for a vertical axis and the isotropic limit stated. Problem Set 1 Q4 and Solutions 1 Q4 updated to match.
