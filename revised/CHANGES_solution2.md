# Changes to solution2.tex ("Seismology Problem Set 2 -- Solutions")

## A. Corrections and rewordings

**Preamble / title**
- `\title{Seismology Problem Set 2 - Solutions}` → `... 2 -- Solutions}` (en-dash). `\usepackage{amsmath}` deleted. 10pt kept. All displays stay unnumbered, so no `\label`s were added.

**Q1**
- "From Lecture 18" → "From Lecture~18"; `_{in}`, `_{out}` → `\mathrm{in}`, `\mathrm{out}` (3×).
- Comma added after the first display (sentence continues "with $\mathbf{d}$ the data vector").

**Q2**
- (a) "unbounded in these null space directions" → "null-space directions" (compound adjective).
- (b) **Maths correction.** After "Rearranging this gives", the display read
  `(\mathbf{I}+\mu\mathbf{A}^{T}\mathbf{A})\mathbf{m}=2\mu\mathbf{A}^{T}\mathbf{d}` → `(\mathbf{I}+\mu\mathbf{A}^{T}\mathbf{A})\mathbf{m}=\mu\mathbf{A}^{T}\mathbf{d}`.
  Reason: from $2\mathbf{m}+2\mu\mathbf{A}^{T}(\mathbf{A}\mathbf{m}-\mathbf{d})=\mathbf{0}$ both terms carry the factor 2, which cancels; the stray 2 was inconsistent with the following (correct) line $\mathbf{m}=(\mathbf{A}^{T}\mathbf{A}+\mu^{-1}\mathbf{I})^{-1}\mathbf{A}^{T}\mathbf{d}$.
- (b) `( \mathbf{A}^{T}\mathbf{A} + \frac{1}{\mu}\mathbf{I})^{-1}` → `\left( ... \right)^{-1}` (tall brackets).
- (b) `\mathbf{m}_{reg}` → `\mathbf{m}_{\mathrm{reg}}`; "In standard  regularisation" → single space.

**Q3**
- $\sigma_{i}$ → $\varepsilon_{i}$ in the definition of $\mathbf{h}'$, to match the renamed standard deviations in the problem (where $\sigma$ is the surface load).
- `\mathbf{u}_{i}^{obs}` → `\mathbf{u}_{i}^{\mathrm{obs}}`.
- "variation w.r.t to the multipliers" → "variation with respect to the multipliers"; "the sensitivity kernel w.r.t $\sigma$" → "with respect to $\sigma$".
- **Maths correction.** The shear-modulus term. With $\delta A_{ijkl}=\delta\mu(\delta_{ik}\delta_{jl}+\delta_{il}\delta_{jk})$ one finds $\delta A_{ijkl}\,\partial_{l}u_{k}\,\partial_{j}u'_{i}=\delta\mu\,\partial_{j}u'_{i}(\partial_{j}u_{i}+\partial_{i}u_{j})$, which is exactly one half of $(\partial_{j}u_{i}+\partial_{i}u_{j})(\partial_{j}u'_{i}+\partial_{i}u'_{j})$ (expand the product and relabel $i\leftrightarrow j$ in two of the four terms). A factor $\frac{1}{2}$ has therefore been inserted in front of the $\delta\mu$ integral in the $\delta J$ display and in
  `K_{\mu}=\left(\ldots\right)\left(\ldots\right)` → `K_{\mu}=\frac{1}{2}\left(\ldots\right)\left(\ldots\right)`.
  $K_{\lambda}$ and $K_{\sigma}=-g_{i}u_{i}'$ are unchanged (checked).

**Q4**
- **Sign correction.** "The part of the Lagrangian that depends explicitly on the motion is" `-\frac{1}{2}G\int_{M}\frac{\rho\rho'}{\{\ldots\}^{1/2}}\dd^{3}\mathbf{x}'` → `+\frac{1}{2}G\int_{M}\ldots` (leading minus removed). Reason: $\mathcal{L}\supset-\frac{1}{2}\rho\zeta$ with $\zeta=-G\int\rho'/|\bphi-\bphi'|\,\dd^{3}\mathbf{x}'$, so the term is $+\frac{1}{2}G\rho\int\rho'/|\bphi-\bphi'|$. Differentiating $|\bphi-\bphi'|^{-1}$ produces the minus sign, so all subsequent displays (with $-\frac{1}{2}G$, $-G$ and the $3/2$ power) are correct as they stood and were left unchanged; the final result $\delta\mathcal{L}/\delta\varphi_{i}=\rho\gamma_{i}$ follows.
- `\boldsymbol{\varphi}+\delta\boldsymbol{\varphi}` → `\bphi+\delta\bphi`; "un-primed functions" → "unprimed functions".

**Q5**
- Azimuthal coordinate: `r=b+sh_{1}(\theta,\phi)` → `r=b+sh_{1}(\theta,\varphi)` — the problem uses $\varphi$, and $\phi$ is the gravitational potential in this question, so `\phi` was a clash.
- `\frac{d\phi_{0}}{dr}` → `\frac{\ddns\phi_{0}}{\ddns r}` (2×), `\frac{d}{dr}` → `\frac{\ddns}{\ddns r}`; `\cdot\cdot\cdot` → `\cdots` (4×).
- The first-order expansion of the jump `[\phi_{0}+s\frac{d\phi_{0}}{dr}h_{1}+\ldots]_{-}^{+}` given tall brackets `\left[ \ldots \right]_{-}^{+}` (contains a fraction).

**Q6**
- "that on the source locations and source time is non-linear" → "on the source location and source time" (single source).
- "find the best fitting source parameters" → "best-fitting"; "it is likely to be over-determined" → "overdetermined" (STYLE.md §4).
- "(in a spherically symmetric earth at least)" → "(in a spherically symmetric earth model at least)" (STYLE.md capitalisation rule).

**Q7**
- "Note that that the Coriolis operator" → "Note that the Coriolis operator" (duplicated word).
- `\mathbf{P}`, `\mathbf{W}`, `\mathbf{H}` → plain `P`, `W`, `H` in all bra-kets (6×); `i\omega`, `\frac{1}{2}i` → `\ii\omega`, `\frac{1}{2}\ii`.
- `|(\mathbf{s}\cdot\boldsymbol{\Omega})|^2-||\mathbf{s}||^2||\boldsymbol{\Omega}||^2` → `|\mathbf{s}\cdot\bm{\Omega}|^{2}-\|\mathbf{s}\|^{2}\|\bm{\Omega}\|^{2}`; "Cauchy-Schwarz" → "Cauchy--Schwarz".
- Comma added after the centrifugal-energy display (sentence continues "and using the identity ...").

## B. Left unchanged / questions for the author

- Q2(a): the answer describes the compatible region as a "hyper-cylinder (or more precisely, an ellipsoidal cylinder)". Fine for $m>n$ with $\mathbf{A}$ of full row rank; if $\mathbf{A}$ is rank-deficient the cross-section is a degenerate ellipsoid, but this level of detail seems unnecessary here.
- Q3: after the integration by parts the condition is written with the sum over the observation points hidden inside $\mathbf{h}'$, which contains delta functions on the surface; the surface integral $\int_{\partial M}h_{i}'\delta u_{i}\,\dd S$ is therefore the usual formal shorthand. Left as is.
- Q3: the adjoint boundary condition $\hat{n}_{j}(A_{ijkl}\partial_{l}u'_{k})=-h'_{i}$ and "$w_{i}'=u_{i}'$ on $\partial M$" are consistent with the stated Lagrangian; no change.
- Q5: $\phi_{0}$ depends on $r$ only, but the display for $(\hat{\mathbf{r}}\cdot\nabla)^{2}\phi_{0}$ uses partial derivatives $\partial^{2}\phi_{0}/\partial r^{2}$, $\partial\phi_{0}/\partial r$ while the spherically symmetric Poisson equation two lines earlier uses $\ddns/\ddns r$ (now with `\ddns`). Harmless mixed notation; left as written.
- Q5: the statement that $\hat{\mathbf{n}}_{1}\cdot\nabla\phi_{0}$ "is zero" is slightly loose — it is the *jump* of this term that vanishes (as the sentence goes on to explain). Not reworded, since the following clause makes the meaning clear.
- Q6: the time integral is evaluated as $[\cos[\omega_{k}(t-t')]/\omega_{k}^{2}]_{t_{s}}^{t}$; checked, correct (the derivative of $\cos[\omega_{k}(t-t')]$ with respect to $t'$ is $+\omega_{k}\sin[\omega_{k}(t-t')]$).
- Q7: the sign of the centrifugal contribution to the potential energy form is asserted rather than derived ("This latter part can be written explicitly as ..."); it is consistent with the conclusion that centrifugal forces are destabilising, so left unchanged.
- Q7: "Coriolis operator is anti-Hermitian, and hence $\langle\mathbf{s}|W|\mathbf{s}\rangle$ is an imaginary number" — correct; the phrase "the contribution from the Coriolis term always being non-negative" refers to $-\frac{1}{4}\langle\mathbf{s}|W|\mathbf{s}\rangle^{2}\ge0$. Left as is.

## C. Style pass

Applied: `amsmath` removed; title en-dash; `_{in}`/`_{out}`/`_{reg}`/`^{obs}` → `\mathrm`; `\boldsymbol` → `\bphi`/`\bm{\Omega}`; `||` → `\|`; `i` → `\ii` (2×); `d/dr` → `\ddns` (3×); `\cdot\cdot\cdot` → `\cdots` (4×); `\mathbf{P/W/H}` → `P/W/H`; `w.r.t` → "with respect to" (2×); `Lecture 18` → `Lecture~18`; Cauchy--Schwarz; over-determined → overdetermined; `\left( \right)` and `\left[ \right]` on two displays containing fractions; two terminal commas added to displays; trailing spaces at line ends removed; line breaks otherwise preserved. Not applied: I tried `\left[ \right]` on the long "condition now reads" display in Q3 but it increased a pre-existing overfull box (4.0pt in the original) to 13pt, so the plain `[ ]` was kept there; the 4.0pt overfull hbox at that equation is inherited from the original. Compiles with no errors or undefined references (5 pages, as the original); latexdiff runs cleanly.
