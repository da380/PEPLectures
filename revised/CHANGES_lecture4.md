# Changes to lecture4.tex (Lecture 15: Plane wave propagation)

Equation numbering: 44 numbered displayed equations in the original (two `align`
blocks count 2 each: eqs 39--40 and 42--43). All hard-coded references in the
original -- eq.(1), (2), (10), (13) x2, (14), (16) x2, (22), (38) and fig.1 --
pointed at the correct equations/figure, so no mapping corrections were needed.
They are now `eq.~(\ref{eq:N})` / `Fig.~\ref{fig:1}`; the compiled PDF resolves
each to the same number as the original.

## A. Corrections and rewordings

### Preamble
* Removed the leading blank line and the space in `\documentclass [a4paper,12pt]`.

### Derivation of the Christoffel equation
* "the terms ``homogeneous'' and ``isotropic'' are sometimes confused they mean"
  -> "are sometimes confused, they mean" (missing comma).
* Plane wave ansatz: `u_{i}(x,t)` -> `u_{i}(\mathbf{x},t)` (position argument is a vector).
* Terminal punctuation added to displayed equations that lacked it: the ansatz
  (eq:2, comma), `c=1/\|\mathbf{p}\|` (eq:4, full stop), the divergence identity
  (eq:6, comma before "because"), the combined second-derivative result (eq:9,
  full stop), the Christoffel-matrix definition (eq:11, full stop), the
  Christoffel equation (eq:12, comma), the determinant condition (eq:13, full stop).
* Christoffel equation: right-hand side `0` -> `\mathbf{0}` (the equation is
  vector-valued; the next sentence already writes the trivial solution as
  `\mathbf{a}=\mathbf{0}`).
* "Using this definition, then eq.(10) can be written" -> "Using this
  definition, eq.~(\ref{eq:10}) can be written" (stray "then").

### Reduction to a symmetric eigenvalue problem
* Homogeneity relation: `\bm{\Gamma}(\lambda\mathbf{p})=\lambda^{2}\bm{\Gamma}(\mathbf{p})`
  "for any $\lambda$" -> `\bm{\Gamma}(\gamma\mathbf{p})=\gamma^{2}\bm{\Gamma}(\mathbf{p})`
  "for any real $\gamma$" ($\lambda$ is the first Lame parameter later in the
  same file; STYLE.md section 3 forbids reusing a symbol).
* Footnote marker moved after the full stop: "for any real $\gamma$ and
  $\mathbf{p}$.\footnote{...}" (the footnote text itself is unchanged, including
  the "third different way" remark).
* "propagating either forward and backwards along" -> "either forwards or backwards along".
* Removed a double space before "What does this mean physically?".
* "in a fixed propagation direction all we need do is vary" -> "in a fixed
  propagation direction, all we need do is vary" (comma after the long
  introductory clause).
* Footnote on convexity of the inner sheet moved after the full stop
  ("interested.\footnote{...}"); text unchanged.
* Figure: `[h]` -> `[ht]`, label `fig:slowness` -> `fig:1`, "In fig.1" -> "In Fig.~\ref{fig:1}".

### The Christoffel equation in an isotropic material
* "that you have likely have seen before" -> "that you have likely seen before".
* "are the components of a matrix that projects a vector parallel to
  $\hat{\mathbf{p}}$ while" -> "of the matrix that projects ... , while" (comma).
* **Eigenvector displays (eqs 24 and 26).** The original wrote both as
  `\bm{\Gamma}(\hat{\mathbf{p}})\hat{\mathbf{a}}=\ldots\hat{\mathbf{a}}` although
  the surrounding text says "take the polarisation vector parallel to the
  propagation direction" and then "$\mathbf{a}\propto\hat{\mathbf{p}}$ is an
  eigenvector", and $\hat{\mathbf{a}}$ is never defined. Rewritten as
  - "if we take the polarisation vector equal to the propagation direction we obtain"
    `\bm{\Gamma}(\hat{\mathbf{p}})\hat{\mathbf{p}}=\frac{\lambda+2\mu}{\rho}\hat{\mathbf{p}}.`
  - for $\mathbf{a}\perp\hat{\mathbf{p}}$:
    `\bm{\Gamma}(\hat{\mathbf{p}})\mathbf{a}=\frac{\mu}{\rho}\mathbf{a},`
  No change to the mathematics; this only makes the displays consistent with the
  prose and with eq:23.
* "(i.e., $\mathbf{a}\perp\hat{\mathbf{p}}$)" -> "(i.e.\ $\mathbf{a}\perp\hat{\mathbf{p}}$)".
* "The first are P-waves which have phase speed" -> "P-waves, which have"; "The
  next are S-waves which have phase speed $\beta$ which is again independent of"
  -> "S-waves, which have phase speed $\beta$, again independent of" (removes the
  double "which"). Double space before "lie" removed.
* `\mu>0 \quad \kappa=\lambda+\frac{2}{3}\mu>0` -> `\mu>0,\quad \kappa=\lambda+\tfrac{2}{3}\mu>0,`
  (punctuation inside and after the display).
* "while the outer one for S-waves is two-fold degenerate" -> "while the outer
  one, for S-waves, is two-fold degenerate".

### Perturbation theory for slightly anisotropic materials
* Terminal punctuation added: the restated eigenvalue problem (eq:29, full
  stop), the perturbed elastic tensor (eq:30, full stop), the perturbed
  eigenvalue problem (eq:33, comma before "where").
* "Here s is a perturbation parameter" -> "Here $s$ is a perturbation parameter".
* "while the first order term is given by" -> "first-order term" (adjective).
* "from first-principles" -> "from first principles" (noun phrase, no hyphen).
* Thin spaces made consistent in the first-order relations:
  `2\alpha\delta\alpha\hat{\mathbf{p}}` -> `2\alpha\,\delta\alpha\,\hat{\mathbf{p}}`;
  `2\beta\delta\beta` -> `2\beta\,\delta\beta` in eqs 41, 43 and 44 (eq 42
  already had the thin space; the two lines of the align block now match).
* "where we have used the symmetry of the Christoffel equation along with" ->
  "the symmetry of the Christoffel matrix" (it is the matrix that is symmetric,
  as shown in eq:17).
* Double space in "however,  now" removed.
* "It therefore seems reasonable to look for solutions" -> "look for solutions
  of the form" (the sentence leads into a display of the assumed form).
* "solve the above eigenvalue problem in closed-form" -> "in closed form".
* "note the first-order term $\delta\mathbf{a}$ will in general" -> "note that
  the first-order term".
* "This phenomena occurs for waves in the Earth" -> "This phenomenon occurs".

### Summary
* `itemize` -> `enumerate` (labels `\item[(i)]` etc. unchanged).

## B. Left unchanged / questions for the author
* "In a previous lecture we obtained the linearised equations of motion" and
  "We recall from the last lecture that in an isotropic material ..." -- left as
  is; if a specific reference is wanted it should presumably be `Lecture~13`
  (linearised equations) and `Lecture~14` (isotropic elastic tensor), but I did
  not want to guess the lecture numbers.
* Eq:3 keeps `\text{constant}`; STYLE.md only forbids `\text{diag}`. Could be
  `\mathrm{constant}` if preferred.
* Eqs 39--40 and 42--43 are `align` blocks whose lines are numbered separately
  (as in the original). They are never referred to, so nothing depends on this,
  but the first pair could equally be a single `split` equation.
* The physical remark "in a fluid there is one positive phase speed and two that
  vanish" and "It can be shown that in any real isotropic solid $\mu>0$,
  $\kappa>0$" are stated without proof, as in the original; left unchanged.
* The two overfull `\hbox` warnings (the "anisotropic perturbation" paragraph
  and item (iv) of the summary) are present in the original compile too and are
  a few points at most; not touched.

## C. Style pass
Applied: `||`->`\|`, `u_i(x,t)`->`u_{i}(\mathbf{x},t)`, `\label{eq:N}` on all 44
displayed equations, `eq.~(\ref{...})`/`Eq.~`/`Fig.~\ref{fig:1}`, `[ht]` figure
placement, footnote markers after punctuation, `i.e.\ `, `,\quad` inside the
inequality display, `\tfrac{2}{3}` inline, thin spaces `\,` between products of
perturbation quantities, summary list as `enumerate`. No `\dd`, `\ee`, `\ii`,
units, dashes or `\boldsymbol` occurrences were present in this file. Lame
accent kept (`Lamé`), P-wave/S-wave capitalisation already correct. No place
where the style guide could not be applied.

## D. Second pass (25 September 2026)

* Renumbered: "Lecture 16: Plane wave propagation" (body-only file).
* Fig. 1 regenerated in matplotlib (`figsrc/L4F1.py`): transversely isotropic medium and olivine, labelled axes in s/km; caption rewritten accordingly.
