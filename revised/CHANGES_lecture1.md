# Changes to lecture1.tex ("Lecture 12: Finite elasticity")

## A. Corrections and rewordings

**Preamble / title**
- `\title{Lecture 12: Finite Elasticity}` → `\title{Lecture 12: Finite elasticity}` (sentence case); `\title`/`\author` moved above `\begin{document}` per the STYLE.md skeleton.

**Outline and motivation**
- "a good though not complete model" → "a good, though not complete, model" (parenthetical commas).

**Velocity, deformation gradient, and Jacobian**
- "which is accurate to first-order in $\delta\mathbf{x}$" → "accurate to first order in" (adverbial use: no hyphen, per STYLE.md §4).
- Footnote after "vector field" moved after the full stop; stray space before the closing brace of the footnote removed.

**Conservation of mass**
- "such that eq.~(1) holds at some time $t=t_{0}$ then at this time" → "... $t=t_{0}$, then at this time" (comma after the conditional clause).

**Potential energy**
- Footnote after "dissipative processes" moved after the full stop; inside it "i.e., those forces" → "i.e.\ those forces" (STYLE.md §4).
- (ii) "Motions that do not alter its shape nor volume" → "do not alter its shape or volume" (grammar: "not ... nor" → "not ... or").
- (ii) "Again, necessity of this condition can also be proven." → "Again, the necessity of this condition ...".
- (ii) "the second argument on the right hand side" → "right-hand side".
- (ii) "it is useful to again re-define the potential energy density" → "it is useful once again to redefine the potential energy density" (split infinitive; "redefine").
- "this potential energy density is  called" → single space.

**Hamilton's principle**
- "deformation gradient on the right hand side" → "right-hand side".
- First-variation definition: `h~\delta\bphi` → `h\,\delta\bphi`.
- Footnote after "independently equal to zero" moved after the full stop.
- "the equations of motion  are defined" → single space.

**Equation references**
- No hard-coded numbers were present; all ten `eq.(\ref{eq:N})` changed to `eq.~(\ref{eq:N})`. Existing labels (`eq:1`–`eq:28`, `eq:29a`, `eq:29b`, `eq:30`–`eq:34`, `eq:35a`–`eq:35c`, `eq:36`–`eq:38`) kept unchanged, as instructed. All references resolve to the same printed numbers as in the original PDF.

## B. Left unchanged / questions for the author

- The footnote in "Hamilton's principle" ("For those with less experience with variational principles ...") is attached mid-sentence after "in the usual manner", where there is no punctuation to place it after. Left in place; it could be moved to the end of the sentence (after eq. (24)) if you prefer.
- Existing labels `eq:30`, `eq:31` are on the 31st and 32nd numbered displays (because of `eq:29a`/`eq:29b`); left as they are since the instructions say not to renumber.
- "rigid body motion(s)" (used adjectivally several times) left unhyphenated; STYLE.md does not list it.
- The `\the\year` in `\author` is kept.

## C. Style pass

Applied: title/author relocated to the preamble and title put in sentence case; `eq.(` → `eq.~(` (10×); Piola-Kirchhoff → Piola--Kirchhoff (3×), Euler-Lagrange → Euler--Lagrange (2×); three footnote markers moved after punctuation; `~` → `\,` thin space in `h\,\delta\bphi`; double spaces removed (2×); `i.e.,` → `i.e.\ `. Already conformant and untouched: `\|...\|` norms, `\dd`/`\ddns`, `\bphi`, `\hat{\mathbf{n}}`, `Fig.~\ref{fig:1}`, enumerate `(i)`–`(iv)` in the summary. Compiles with no errors/warnings (8 pages, as the original); latexdiff runs cleanly.

## D. Second pass (25 September 2026)

* Renumbered for the 2026 syllabus: title is now "Lecture 13: Finite elasticity" (all lectures shifted by one; every `Lecture~N` cross-reference in the set updated).
* Fig. 1 replaced by a TikZ redrawing of the hand-drawn sketch (`figsrc/L1F1.tex` → `figures/L1F1.pdf`); same content, vector output.
* Text re-read in full: no further corrections needed. The mid-sentence footnote in "Hamilton's principle" is left where it is (it sits immediately before a displayed equation, where there is no punctuation to follow).
