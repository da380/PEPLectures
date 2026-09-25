# Changes to `lecture3.tex` (Lecture 14: Seismic sources)

Equation numbering: 28 numbered displayed equations in the original, now labelled `eq:1`--`eq:28`
in order. All five hard-coded references in the original ("Eq. (3)", "eq. (9)", "eq. (10)" twice,
"eq. (20)", "eq. (12)") pointed at the correct equations; no remapping was needed. Figure labels
`fig:1`--`fig:3` were already in place.

## A. Corrections and rewordings

**Outline and motivation**
* `\textbf{kinematic  source approximation}` → single space.

**The elastic rebound theory for earthquakes**
* "slow build up, and then rapid release" → "slow build-up" (noun form).
* "whose cause remains uncertain\footnote{...}. There are also" → "uncertain.\footnote{...} There are also" (footnote marker after the full stop, house style).
* "the process must seemingly be different than for shallow earthquakes" → "different from that for shallow earthquakes".
* "including nuclear weapon tests, and explosions generated as part of oil prospecting, and small ... earthquakes" → "including nuclear weapon tests, explosions generated as part of oil prospecting, and small ... earthquakes" (stray "and" in a three-item list).
* "Harry Fielding Reid (1859-1944)" → "(1859--1944)".
* "as summarised in fig.2" → `Fig.~\ref{fig:2}`.
* `"theory"` → ``theory''.
* "earthquakes are a so-called \textbf{critical phenomena}" → "critical phenomenon".

**Introduction of the stress glut**
* "the time-interval of interest" → "the time interval of interest".
* Added the comma that the sentence grammar requires before eqs. (3), (5), (10) and after "is defined in the expected manner" (eq. (4)); eq. (7) now ends with a comma ("... = 0," followed by "and similarly for the boundary conditions"); eq. (10) now ends with a comma before "where".
* "Eq. (3) implies" → `Eq.~(\ref{eq:3}) implies`.
* "does not alter the net linear nor angular momentum" → "does not alter the net linear or angular momentum".
* "the resulting Euler-Lagrange equations remains unchanged" → "the resulting Euler--Lagrange equations remain unchanged," (subject--verb agreement).
* "Using eq. (9), we see that eq. (10) implies" → `eq.~(\ref{eq:9})`, `eq.~(\ref{eq:10})`.
* Eq. (11) ended with a full stop followed by "and can use this relation": full stop → comma, and "and can use" → "and we can use".
* Removed double spaces in "the  linearised equations", "the  elastic constitutive relation", "is  consistent with".

**Linearised equations of motion incorporating a stress glut**
* `u_{i}(\mathbf{x}, 0) = 0, \quad` → `= 0,\quad` (house style).
* "Motivated by eq.(10)" → `eq.~(\ref{eq:10})`.
* Trailing double space removed after "linearised equations of motion".

**The stress glut for an idealised fault**
* Fig. 3 caption: "a zone of finite-width which" → "a zone of finite width which"; likewise "across the finite-width of a fault zone" → "across the finite width of a fault zone" (noun, no hyphen).
* "wavelength substantially larger that the width" → "larger than the width".
* "fault zones widths of metres" → "fault zone widths of metres".
* "\textbf{fault plane }, $\Sigma$,  passes" → "\textbf{fault plane}, $\Sigma$, passes".
* "and hence it can be decomposed in the form" → "and hence the displacement vector can be decomposed in the form" (the "it" grammatically referred to the discontinuity, but it is the displacement field that is decomposed).
* "the second has a jump discontinuity the fault plane" → "across the fault plane".
* "within eq. (20)" → `eq.~(\ref{eq:20})`; "Specialising eq. (12)" → `eq.~(\ref{eq:12})`.
* Eq. (22) (and its repeat as eq. (26)): `A_{ijkl}w_{k}n_{l}\delta_{\Sigma}` → `A_{ijkl}w_{k}\hat{n}_{l}\delta_{\Sigma}`, to match the text, which calls the unit normal $\hat{\mathbf{n}}$, and the $\hat{n}_{j}$ used in the boundary conditions.
* Eq. (23): `\int_{\Sigma}f\,dS` → `\int_{\Sigma}f\dd S`.
* Double spaces removed in "than  something" and "surface,  $\Sigma$".

**Formulating an inverse problem**
* Eq. (25) now ends with a comma before "while the stress glut ...".
* Model-parameter item 4: "Geometric parameters describing the shape of the earth model including any internal discontinuities and for the fault plane." → "Geometric parameters describing the shape of the earth model, including any internal discontinuities, and the geometry of the fault surface."
* "to solve the equations though this generally requires" → "to solve the equations, though this generally requires".
* "$i=1, ..., n$" → `$i=1,\dots,n$` (twice); "Suppose  we", "seismometers,  with" double spaces removed.
* "a polynomial, $p(x)$ of degree $m$ such that" → "a polynomial, $p(x)$, of degree $m$ such that".
* **Maths correction.** "(e.g. the $y_i$ all lie on a line, and $m>1$)" → "$m\ge 1$": a degree-one polynomial already interpolates collinear points, so $m=1$ must be included.
* **Maths correction.** "the polynomial is unique if and only if $m+1 = n$, while if $m+1 > n$ there are an infinite number of possible answers" → "the polynomial is unique if and only if $m+1 \le n$, while if $m+1 > n$ there are infinitely many possible answers". When $m+1<n$ and a solution exists it is unique, because the $n\times(m+1)$ Vandermonde matrix for distinct $x_i$ has full column rank; the original statement wrongly excluded this case (it also contradicted the preceding sentence, which had just discussed existence for $m+1<n$).
* "with  fit to the data need not be exact, but in a statistical sense" → "with the fit to the data no longer required to be exact, but only adequate in a statistical sense".
* `e.g.`/`i.e.` followed by a word given the control space (`e.g.\ `, `i.e.\ `) in four places.

**What you need to know and be able to do**
* `itemize` → `enumerate` (existing `\item[(i)]` labels kept).
* "Reid's elastic rebound theory for Earthquakes" → "for earthquakes".

## B. Left unchanged / questions for the author

* Eq. (27), `p(x_i) = y_i.`, does not say "for $i=1,\dots,n$"; the preceding sentence makes this clear, so I left it, but `,\quad i=1,\dots,n` could be added if desired.
* Eq. (12) defines the stress glut as $\overline{T}_{ij}=F_{ik}\overline{S}_{kj}$, i.e. a first Piola--Kirchhoff-type (not necessarily symmetric) object, while the later text says "within this linearised theory the tensors $\overline{S}_{ij}$ and $\overline{T}_{ij}$ agree to first order". This is correct as stated ($\mathbf{F}=\mathbf{I}+O(s)$ and $\overline{S}=O(s)$), so nothing was changed; noted only because the phrase "either being called the stress glut" might be worth a half-sentence of explanation.
* "The exceptions to this picture are (i) ... and (ii) deep earthquakes ... . There are also man-made earthquakes" -- the man-made examples are arguably a third exception; left the author's structure.
* The sentence "In doing so, we see immediately that our equations of motion imply that earthquakes cannot happen. To rectify things, we will not produce a complete dynamical theory ..." reads slightly oddly ("To rectify things, we will not ..."); left as it is the author's voice, but "To rectify things we will not attempt a complete dynamical theory ..." might read better.
* Item 3 of the well-posedness list, "The output depends continuously on the input", and the later remark "these conditions only make sense sequentially" are fine; no change.

## C. Style pass

Applied: `\documentclass[...]` spacing; `\boldsymbol{\varphi}` → `\bphi` (4); `||\cdot||` → `\|\cdot\|` (3); `dS` → `\dd S`; `n_{l}` → `\hat{n}_{l}` (2); en-dashes for Piola--Kirchhoff (2), Cauchy--Green, Euler--Lagrange, (1859--1944); `Lecture~12`, `Lecture~13` (3); `Fig.~\ref{}`/`eq.~(\ref{})` throughout; straight quotes → ``''; `...` → `\dots` (2); `,\quad`; `e.g.\ `/`i.e.\ `; footnote marker after punctuation; `\label{eq:N}` on all 28 displayed equations; summary list `itemize` → `enumerate`. No place where the style guide could not be applied cleanly. Untouched lines (including their trailing whitespace) are byte-identical to the original; the diff touches only edited lines and the 28 `\begin{equation}` lines.

## D. Second pass (25 September 2026)

* Renumbered: "Lecture 15: Seismic sources" (body-only file; title in `lectures.tex`).
* Outline: "To rectify things, we will not produce ..." → "We will not attempt to rectify this by producing a complete dynamical theory ..., as no such theory yet exists. Instead, we simply modify ...".
* Fig. 1 replaced by a TikZ cross-section of a normal fault (the previous image carried a Pearson copyright notice) and cited in the first paragraph; Fig. 2 redrawn in TikZ; Fig. 3 (photograph) now cited where the finite fault zone is discussed.
* "it is readily shown that the total energy ... is conserved" → "it can be shown (and you will do so in the first problem set)".
