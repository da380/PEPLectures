# Changes to lecture5.tex ("Lecture 16: Ray theory")

## A. Corrections and rewordings

**Preamble / title**
- Leading blank line removed; `\documentclass [a4paper,12pt]` → `\documentclass[a4paper,12pt]`; "Outline and Motivation" → "Outline and motivation".

**Outline and motivation**
- "the behaviour of quantum mechanical systems are studied" → "is studied" (subject is "behaviour").
- "length-scale of the heterogeneity" → "length scale" (noun form, STYLE.md §4).

**Relative length scales of heterogeneity**
- "contained within the half-plane $p_{i}^{0}x_{i}\ge0$. Within the homogeneous half plane" → "half-space" in both places (the regions are three-dimensional).
- Fig. 1 caption: "inter-\naction" → "interaction" (hard hyphenation from PDF paste); "plane p-wave" → "plane P-wave"; "length-scale" → "length scale" (2×); "wave field" → "wavefield".
- Fig. 2 caption: "plane p-wave" → "plane P-wave"; "length-scale" → "length scale"; "wave front" → "wavefront".

**The ray series ansatz**
- List item 1: "as they propagates" → "as they propagate".
- "the wavefronts of eq.(5) are, by definition level surfaces" → "are, by definition, the level surfaces".
- "within homogeneous regions where the travel time depends linearly on $\mathbf{p}$" → "depends linearly on $\mathbf{x}$" (in a homogeneous region $T=p_{i}^{0}x_{i}$, which is linear in position; $\mathbf{p}$ is the constant gradient).
- Fourier paragraph rewritten minimally. Before: "To proceed, it is useful to decompose displacement vector fields into an integral over harmonic plane waves [eq. 8, with `e^{\ii[t-T(x)]}` -- the $\omega$ missing from the exponent], It is clear that eq.(5) is the temporal convolution ...". After: "To proceed, it is useful to write the waveform $f$ as a Fourier integral, so that eq.~(\ref{eq:5}) becomes [eq. 8 with `\ee^{\ii\omega[t-T(\mathbf{x})]}`], where $\tilde{f}$ denotes the Fourier transform of $f$. It is clear that eq.~(\ref{eq:5}) is the temporal convolution ..." (missing $\omega$ in the exponent restored; the run-on sentence split; $\tilde{f}$ defined before it is used).
- "with the waveform f" → "with the waveform $f$".
- "is equivalent to the eq.(9)" → "is equivalent to eq.~(\ref{eq:9})".
- "each successive term within the series is divided by a factor of $1/(\ii\omega)$" → "is multiplied by a further factor of $1/(\ii\omega)$" (dividing by $1/(\ii\omega)$ would be multiplying by $\ii\omega$, i.e. differentiation, contradicting the following clause).

**Substitution into the equations of motion**
- Divergence `multline` (eq. 16) now ends with "," instead of "." since the sentence continues "and substituting into ...".
- Eq. 17: "$=0 .$" → "$=0.$".
- Eq. 20: `(A_{ijkl}\frac{\partial a_{k}^{n-2}}{\partial x_{l}})` → `\left( ... \right)`; terminal full stop added.
- "gives the orientation of $a^{0}(x)$" → "$\mathbf{a}^{0}(\mathbf{x})$"; "all the $a^{n}$" → "all the $\mathbf{a}^{n}$" (consistent with the bold vector notation used from eq. 22 onwards).
- "it is reasonable to ask where exactly we have made an approximation?" → "... approximation." (indirect question).
- "\textbf{the ray series  is  not convergent ...}" double spaces removed; "with  an approximate" likewise.
- Footnote marker moved after the full stop: "...need not be small}\footnote{...}." → "...need not be small}.\footnote{...}"; "Lars Hormander" → "Lars H\"{o}rmander".

**The eikonal equation**
- "we must first reduce eq.(18) into a single partial differential equation" → "reduce eq.~(\ref{eq:18}) to a single ...".
- Eq. 21 now ends with "," (sentence continues "for all ..."); "for all $x\in\mathbb{R}^{3}$ and vectors $\mathbf{a}$ and $\mathbf{p}$" → "for all $\mathbf{x}\in\mathbb{R}^{3}$ and $\mathbf{p}$" ($\mathbf{a}$ does not appear in the definition).
- "where $A_{0}(\mathbf{x})$ a scalar amplitude" → "is a scalar amplitude"; "normalised eigenvector  which" single space.
- "and can therefore, always write" → "and can therefore always write".
- Eq. 25: third factor had `||p(\mathbf{x})||` (non-bold p) → `\|\mathbf{p}(\mathbf{x})\|`.

**The method of characteristics**
- Eq. 27 (ray Hamiltonian): `H_{k}(\mathbf{x},\mathbf{p})=\frac{1}{2}||\mathbf{p}(\mathbf{x})||^{2}c_{k}[\mathbf{x},\hat{\mathbf{p}}(\mathbf{x})]^{2}` → `\frac{1}{2}\|\mathbf{p}\|^{2}c_{k}(\mathbf{x},\hat{\mathbf{p}})^{2}`. Reason: $H_{k}$ is defined as a function on the six-dimensional phase space of independent $(\mathbf{x},\mathbf{p})$, as the following clause says; writing $\mathbf{p}(\mathbf{x})$ inside the definition presupposes the field being solved for. Eq. 28 (which evaluates $H_{k}$ on $\mathbf{p}(\mathbf{x})$) is unchanged. **Please check.**
- "to write the eikonal equation eq.(26) in the alternate form" → "the eikonal equation, eq.~(\ref{eq:26}), in the alternative form".
- "corresponds to a \textbf{three-dimensional level surfaces}" → "\textbf{three-dimensional level surface}".
- "The \textbf{initial data} we have for the problem  is that  on the plane" spacing tidied.
- Eqs. 30, 32: `d\sigma` → `\ddns\sigma` (five occurrences); eq. 30 terminal comma added.
- "As we have $H_{k}(\mathbf{x}^{0},\mathbf{p}^{0})=\frac{1}{2}$ each point on the initial surface" → "at each point".
- "Using the method's of \textbf{symplectic geometry}" → "methods".
- "We first make use of eq.(7) and (30)" → "eqs.~(\ref{eq:7}) and (\ref{eq:30})"; eq. 34 terminal full stop added.
- Footnote: "the only case you might be asked to reproduce this derivation" → "the only case in which you might be asked ...".
- **Reference correction:** "and so eq.(36) has the trivial solution $T[\mathbf{x}(\sigma)]=\sigma$" → `eq.~(\ref{eq:37})`. Eq. 36 is Euler's identity $p_{i}\partial H/\partial p_{i}=2H_{k}$; the equation whose solution is $T=\sigma$ is $\frac{\ddns}{\ddns\sigma}T[\mathbf{x}(\sigma)]=1$, the 37th display.
- "In order for the travel time $T$ ... to be a single-valued function of $\mathbf{x}$, then it is necessary that" → "..., it is necessary that" (stray "then").
- "calculate a travel-time through the above method" → "travel time" (noun).
- "at the the behaviour" → "at the behaviour"; "scalar amplitude $A^{0}(\mathbf{x})$ introduced in eq.(24)" → "$A_{0}(\mathbf{x})$" (matches eq. 24).
- "asymptotic techniques  allow", "specified on  a curve" double spaces removed.

**What you need to know and be able to do**
- `itemize` → `enumerate` (labels `(i)`--`(v)` unchanged).
- Item (iv): "within the final example sheet" → "within the first problem set" (the ray-tracing problems are Problems 6 and 7 of Problem Set 1).

**Equation and figure references**
- 38 displayed equations labelled `eq:1`--`eq:38` in order (the two-line `align` is `eq:14`/`eq:15`; the two `multline`s are `eq:16`/`eq:17`). All hard-coded `eq.(N)` / `Eq.(N)` replaced by `eq.~(\ref{eq:N})` / `Eq.~(\ref{eq:N})`; "eq.(7) and eq.(13)" → "eqs.~(\ref{eq:7}) and (\ref{eq:13})". Every original number was correct except eq.(36) → eq:37 (see above). Figure labels `fig:fig1..3` → `fig:1..3` (they are not referenced in the text).
- Terminal punctuation added to displayed equations that lacked it where the sentence requires it: eqs. 1, 2, 3, 4, 5, 7, 12, 20, 21, 25, 26, 30, 34 (commas or full stops as appropriate). Eqs. 9 and 38, which the sentence runs straight on from ("with the waveform $f$", "along the curve ..."), left without.

## B. Left unchanged / questions for the author

- Eq. 27 change described above (dropping the $(\mathbf{x})$ arguments on $\mathbf{p}$ inside the definition of $H_{k}$) is the only edit that touches mathematical notation beyond typography; revert if you prefer the original.
- Eqs. 32, 34, 36 use $H$ without the subscript $k$ (and eq. 36 mixes $H$ on the left with $H_{k}$ on the right). Left as written since the text says arguments are suppressed for clarity, but you may wish to make the subscript consistent.
- "for all real $\gamma$" in the homogeneity statement (eq. 35): for $\gamma<0$ this relies on $c_{k}(\mathbf{x},-\hat{\mathbf{p}})=c_{k}(\mathbf{x},\hat{\mathbf{p}})$, which holds by the symmetry of the Christoffel operator. Left unchanged.
- "Let $(\mathbf{x}^{0},\mathbf{p}^{0})$ be a point on the initial plane $S^{0}$" -- $S^{0}$ was introduced as a "two-dimensional level surface" and is later called "the initial surface". It is indeed an affine plane in phase space, so left as is.
- Eq. 6 uses `\text{constant}`; left (readable and compiles; `\mathrm{constant}` would be equivalent).
- The statement in the Fourier paragraph is now that eq. 5 is the convolution of the singular wave with $f$; the sentence "from now on we set $\tilde{f}(\omega)=1$ or equivalently $f(t)=\delta(t)$" follows naturally. No content added beyond the definition of $\tilde{f}$, which the original used without introducing.
- "right hand side" in "where on the right hand side we have suppressed arguments" left unhyphenated (adverbial position, not attributive).

## C. Style pass

Applied: `\documentclass[...]` spacing; section title case; straight `"..."` → ``...'' (2×); `x` → `\mathbf{x}` in all position arguments (`u_{i}(\mathbf{x},t)`, `T(\mathbf{x})`, `a_{i}(\mathbf{x})`, `a_{i}^{n}(\mathbf{x})`, `u_{i}(\mathbf{x},\omega)`, `\Gamma_{ik}(\mathbf{x},\mathbf{p})`, `\mathbf{x}\in\mathbb{R}^{3}`); `e^{` → `\ee^{` (8×); `||...||` → `\|...\|` (6×); `d\sigma` → `\ddns\sigma` (5×); `\left( \right)` in eq. 20; "length-scale" → "length scale" (5×), "wave field"/"wave front" → "wavefield"/"wavefront", "travel-time" (noun) → "travel time", "half-plane" → "half-space"; "p-wave" → "P-wave" (2×); footnote marker after punctuation (1×); double spaces removed (8×); `itemize` → `enumerate` in the summary; labels `eq:1`--`eq:38`, `fig:1`--`fig:3`; all references via `\ref` with `~`. Compiles twice with no errors, no undefined/multiply-defined labels and no overfull boxes (10 pages, as the original); pdftotext confirms each reference resolves to the same number as in the original PDF apart from the intended (36) → (37); latexdiff runs cleanly (38 kB output). Line breaks and paragraphing of untouched text preserved.

## D. Second pass (25 September 2026)

* Renumbered: "Lecture 17: Ray theory" (body-only file).
* Caustics: now defined as the envelope of the ray family bounding the region of multipathing, with the amplitude infinite on the caustic.
* $H$ → $H_k$ in the conservation calculation and Euler's theorem line; "right hand side" hyphenated.
* Appendix (non-examinable) added with the proof of Euler's theorem for homogeneous functions; the footnote now points to it.
* Fig. 3 cited in the text; new Fig. 4 (rays through a smooth random medium) added with a sentence on the generic formation of caustics.
* Figs. 1–4 regenerated from new Python codes (`figsrc/fd2d_psv.py`, `figsrc/raytrace2d.py`); captions updated.
