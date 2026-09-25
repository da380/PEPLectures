# Changes to lecture12.tex (Lecture 23: Mode splitting and coupling)

Equation count in the original: 35 numbered displays (30 `equation` + 3 lines of the
vector-spherical-harmonic orthogonality `align` + 1 each from the two `align`s that use
`\nonumber`). The coordinating agent's note said 40; the three hard-coded references
nevertheless resolve consistently with a count of 35, as checked against the original PDF.

## A. Corrections and rewordings

### Outline and motivation
* "you are not expected to the details" → "you are not expected to know the details" (missing word).
* "radial and laterally variations of density" → "radial and lateral variations of density".

### Vector spherical harmonics
* "the eigenspace of $\nabla_{1}^{2}$ with degree l" → "$l$"; "for suitably regular function $f$" → "functions $f$".
* "The following orthogonality relations can be established" → added a colon before the display.
* `\hat{\boldsymbol{\theta}}`, `\hat{\boldsymbol{\varphi}}` → `\hat{\bm{\theta}}`, `\hat{\bm{\varphi}}` (house style).
* "unit-sphere" → "unit sphere"; "scale-factors" → "scale factors".

### Spheroidal and toroidal modes
* eq. 12 given its terminal comma ("$\ldots=0,$ where $\omega$ is ...").
* Footnote: "five elastic modulii" → "moduli"; footnote marker moved after the comma.
* "the possible eigenfunctions take two possible types" → "the eigenfunctions are of two possible types".
* Fig. 1 caption: "eigenfrequencies at different degree" → "degrees".
* "depends on $l$ but not m" → "$m$".
* **Reference correction:** "Using eq.(8) it is then easy to obtain" → `eq.~(\ref{eq:10})`. Eq. 8 is the
  $\mathbf{A}_{lm}$ orthogonality relation (first line of the `align`), whereas the kinetic-energy calculation
  with $\mathbf{s}=W\mathbf{C}_{lm}$ and $\mathbf{w}=W'\mathbf{C}_{lm}$ uses the $\mathbf{C}_{lm}$ relation, which
  is eq. 10. (`eq.(13)` → `eq:13` and `eq.(22)` → `eq:22` were correct and simply converted.)
* eq. 15 (test function) given a terminal full stop; eq. 18 "$=0 ,$" → "$=0,$".
* eq. 17: plain `( ... )` around the $r\,\ddns W/\ddns r - W$ factors replaced by `\left( ... \right)`, matching eq. 18
  (no change to content).
* "ODE-eigenvalue problem of Sturm-Lioville type" → "ODE eigenvalue problem of Sturm--Liouville type" (spelling).
* "$n=0,1,2,...$" → `n=0,1,2,\dots`.
* "there exist a countable infinity of toroidal modes" → "there exists a countable infinity".
* Fig. 2 caption: "As for Fig 1" → `Fig.~\ref{fig:1}`.
* "two families of normal modes: toroidals and spheroidal" → "toroidal and spheroidal".
* Footnote after "$q=2$ toroidal" moved after the full stop.
* "the triplet $(q, l, n)$" → "$(q,l,n)$" (consistency with the other occurrences).

### Splitting of a multiplet due to rotation and lateral heterogeneity
* "Runge-Lenz" → "Runge--Lenz"; "Zeeman-effect for a hydrogen atom" → "Zeeman effect for the hydrogen atom".
* "what rotation and laterally heterogeneities do" → "lateral heterogeneities".
* eq. 19: terminal "." → "," since the sentence continues with "where we note ...".
* "identified by the triplet (q,l,n)" → "$(q,l,n)$"; "where s is a perturbation parameter" → "$s$".
* `\cdot\cdot\cdot` → `\cdots` (7 occurrences); `s~\delta\mathbf{s}` → `s\,\delta\mathbf{s}`.
* eqs. 23, 24, 26, 28 given terminal commas (each is followed by "while"/"where").
* "A similar expansion to first-order then leads" → "to first order" (adverbial use).
* "in closed-form for small $l$" → "in closed form"; "left hand side" → "left-hand side" (twice in the file).

### The diagonal sum rule and spherically symmetric density variations
* "we consider a Hermitian eigenvalue problem" → "we considered" (referring back).
* eq. 29 given a terminal full stop.
* Footnote "That are near identical to those arising in the quantum theory of angular momentum" →
  "These are nearly identical to those arising in the quantum theory of angular momentum." (was a sentence
  fragment without a full stop); marker moved after the comma. "zero-trace" → "zero trace".
* "played an key role in early studies of Earth structure using observations of free oscillation" →
  "played a key role ... of free oscillations".
* "Great" footnote: comma splice "earthquakes, it does not suggest" → "earthquakes; it does not suggest"; stray
  space before the closing brace removed.
* "in a homogeneous solids compressed under its own weight" → "in a homogeneous solid compressed under its own weight".
* "An example of this can be seen in Fig 3" → `Fig.~\ref{fig:3}`.
* "It is notable, firstly, that these predictions are in reasonable agreement with an observed eigenfrequency. But
  also that while ... there seems to be at least two in the data." → the two sentences joined with ", but also that
  while ...", and "there seems to be" → "there seem to be" (the second sentence was a fragment).
* "to publish in (1961) the first" → "to publish in 1961 the first".
* "plots of travel-time against epicentral angle" → "travel time" (noun); "data base" → "database".
* "shown in Fig. 4" → `Fig.~\ref{fig:4}`.
* "comprised of an iron-nickel alloy as is suggested observations of meteorites" → "iron--nickel alloy as is
  suggested by observations of meteorites".
* "how planets form and differentiation" → "form and differentiate".
* Fig. 4 caption: Unicode apostrophe in "Earth’s" → "Earth's".

### Lateral density variations and mode coupling
* "Lectures 18 and 19" → `Lectures~18 and~19`; "In Fig. 5 images ... are shown" → "In `Fig.~\ref{fig:5}`, images ...".
* "in the early 90s" → "in the early 1990s"; "in the late 90s" → "in the late 1990s".
* "the LLSVPs seem are associated with" → "seem to be associated with".
* eq. 31: terminal "." → "," (followed by "where $\delta T$ ..."); eq. 32: "." → "," (followed by "as required");
  eq. 33: ";" → "," (followed by "from which it is clear ...").
* "we write $X$ for a compositional variables" → "variable".
* "then we cannot determined from the tomographic students discussed whether LLSVPs are be light or dense" →
  "then we cannot determine from the tomographic studies discussed whether LLSVPs are light or dense".
* "and it free oscillations that have the potential" → "and it is free oscillations that have the potential".
* "undertaken in the late 90s by Miaki Ishii and Jeroen Tromp then both at Harvard and Miaki a PhD student" →
  "undertaken in the late 1990s by Miaki Ishii and Jeroen Tromp, then both at Harvard (Ishii being a PhD student
  at the time)".
* "a slight elaboration the first-order splitting theory" → "a slight elaboration of the first-order splitting theory".
* Straight quotes `"slab graveyard"` → ``slab graveyard''; "...accumulated at the base of the mantle through
  subduction." → "subduction?" (it is a question, parallel to the following "Or are they ...?").
* "Earth Scientists" → "Earth scientists"; double spaces in "interact  with", "only  been", "inversions.  Nevertheless,  " removed.
* eqs. 34 and 35: "$\rangle ,$" / "$=0 ,$" → "$\rangle,$" / "$=0,$"; `i\omega` → `\ii\omega` in eq. 35.
* "for all (k,m)" → "for all $(k,m)$".

### What you need to know and be able to do
* `itemize` → `enumerate` (house style).
* **Item (ii) pointed to a non-existent problem.** "How to solve the toroidal mode eigenvalue problem in a homogeneous
  earth model using spherical Bessel functions see the second problem set." → "... using spherical Bessel functions
  (see the worked examples accompanying this lecture). Note that you might be asked to perform such calculations in
  the exam, but would be provided with the necessary information about special functions." There is no such
  problem in Problem Set 2; the author should confirm that the worked examples do contain this calculation.

## B. Left unchanged / questions for the author
* **Summary item (ii)** (see above): the reference now points to "the worked examples accompanying this lecture" as
  instructed. Please confirm that such a worked example exists (or is planned); otherwise the item should be dropped
  or the reference changed.
* "Explicit within this theory is the assumption that on the addition of laterally varying structure distinct
  multiplets do not interact with each other." Arguably "Implicit" is meant (the assumption enters through the
  ansatz eq. 24 rather than being stated); left as written since the ansatz does restrict $\mathbf{s}$ explicitly.
* The test function for toroidal modes is written $W'(r)$, with the prime denoting a different function rather than a
  derivative, while $\ddns W'/\ddns r$ also appears. This is the author's notation and was kept, but a different
  symbol (e.g. $\tilde W$) would avoid possible confusion.
* Footnote on "Great" (1960 Great Chilean earthquake) sits mid-phrase; it could not be moved after punctuation
  without changing the sentence and was left in place.
* PREM is not expanded at its first use in this file (Fig. 1 caption); it is expanded in Lecture 17
  (revised/lecture6.tex), so this was left as is.
* "Benioff et al. (1961)", "Backus and Gilbert ... 1961", "Gilbert in 1971": dates not checked against the literature.
* The sentence "In the case of a homogeneous planet closed-form solutions can be obtained using spherical Bessel
  functions, though such a solution is of little practical use" sits slightly uneasily with summary item (ii), which
  says students may be examined on exactly this calculation. Left as written.
* eq. 27: the Coriolis term is written $\ii\omega[s\langle\mathbf{w}|W|\mathbf{s}\rangle]$ with square brackets around a
  single term (parallel to the other bracketed expansions); left as is.

## C. Style pass
Applied: sentence-case section title; `\label{eq:1}`--`\label{eq:35}` on every numbered display (the two
`\nonumber` `align`s carry one label each); `fig:fig1..5` → `fig:1..5`; all hard-coded `eq.(N)`/`Fig N`/`Fig. N`
replaced by `\ref`s; `\ii` for the imaginary unit; `\dd S` in the three vector-harmonic orthogonality integrals;
`\cdots`, `\dots`; thin space `s\,\delta\mathbf{s}`; `\bm` for bold Greek; `\left(\right)` tall brackets in eq. 17;
en-dashes for Sturm--Liouville, Runge--Lenz, iron--nickel; ``quotes''; no Unicode apostrophes; `Lectures~18 and~19`;
`Fig.~`, `eq.~`; footnote markers after punctuation (except the mid-phrase "Great" footnote); variables in prose set
in maths mode; `enumerate` in the summary; terminal punctuation on all displays. No place where the style guide could
not be applied cleanly. Line breaks of untouched text preserved; `latexdiff` runs cleanly (45 kB output); two
`pdflatex` passes give no errors, no undefined or multiply-defined labels, and the same single 0.8 pt overfull box
as the original.
