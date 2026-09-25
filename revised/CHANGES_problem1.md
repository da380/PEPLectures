# Changes to problem1.tex ("Seismology Problem Set 1")

## A. Corrections and rewordings

**Preamble / title**
- `\title{Problem Set 1}` → `\title{Seismology Problem Set 1}` (to match the solutions file); `\usepackage{amsmath}` removed (loaded by `mycommands.sty`).

**Q1**
- Equation of motion display ended with `;` → `,` (the sentence continues "for an elastic body").
- `$\bm{\varphi}$` → `$\bphi$` (house macro for the motion).
- `\frac{d}{dt}` → `\frac{\ddns}{\ddns t}`; "Within Lecture 12" → `Lecture~12`.

**Q2**
- "The Lagrangian density for the linearised motion for an elastic body" → "... linearised motion of an elastic body".
- "the standard Euler Lagrange equations" → "Euler--Lagrange".
- "show the following equality holds" → "show that the following equality holds".

**Q3**
- "The elastic tensor, $A_{ijkl}$ is defined by" → "The elastic tensor, $A_{ijkl}$, is defined by" (matching comma).
- "right hand side" → "right-hand side"; Cauchy-Green → Cauchy--Green. `\mathbf{C}=\mathbf{F}^{T}\mathbf{F}` and the symmetries checked; correct.

**Q4**
- `\hat{\boldsymbol{\nu}}` → `\hat{\bm{\nu}}`.

**Q5**
- `i\hbar` → `\ii\hbar` (2×); `-i(E t-\ldots)` → `-\ii(\ldots)` (2×).
- Eikonal display `H(\mathbf{x},\nabla\varphi)=E` → `=E,` (followed by "where").

**Q6**
- "The travel-time of a ray" → "The travel time of a ray" (noun).
- `\frac{d x_{i}}{d\gamma}` → `\frac{\ddns x_{i}}{\ddns\gamma}`.
- "where the $\gamma\mapsto \mathbf{x}(\gamma)$ is the ray path" → "where $\gamma\mapsto\mathbf{x}(\gamma)$ is the ray path" (stray "the").
- Euler-Lagrange → Euler--Lagrange. Hint kept.

**Q7**
- "lies in the $x-z$ plane" → "lies in the $(x,z)$-plane"; "is traveling horizontally" → "is travelling horizontally".
- `\frac{d x_{i}}{dz}` → `\frac{\ddns x_{i}}{\ddns z}`; `2q\alpha~\dd z` → `2q\alpha\,\dd z`.
- Snell's-law display `\frac{\sin\theta}{\alpha}=q` → `=q,` (followed by "where").

**Equation references**
- None present; all displays are `equation*`/`align*` and remain unnumbered and unlabelled, as instructed.

## B. Left unchanged / questions for the author

- "Hamiltonian ray tracing equations" (Q6) left unhyphenated; STYLE.md lists "travel-time curve" but not "ray-tracing".
- "co-ordinate" (Q7) left as is (also used in the solutions); STYLE.md gives no ruling.
- The plane-wave display in Q5 (`\psi=a\exp[\ldots]`) is followed by "and obtain ..." with no punctuation; grammatically fine, left as is.
- "Schrödinger" is entered as a UTF-8 character; compiles fine, so left.

## C. Style pass

Applied: amsmath removed; title aligned with solutions; `\frac{d}{dt}`/`d x_i/d\gamma`/`d x_i/dz` → `\ddns`; `i` → `\ii`; `\boldsymbol` → `\bm`; `\bm{\varphi}` → `\bphi`; hyphen → en-dash in Piola--Kirchhoff (2×), Cauchy--Green, Euler--Lagrange (2×); `~` → `\,` thin space; `Lecture~12`; `$(x,z)$-plane`; British "travelling"; trailing commas added to two displays that run into "where". Compiles with no errors or warnings (3 pages, as the original); latexdiff runs cleanly.
