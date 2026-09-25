# Changes to solution1.tex ("Seismology Problem Set 1 -- Solutions")

## A. Corrections and rewordings

**Preamble / title**
- Unicode en dash in `\title{Seismology Problem Set 1 – Solutions}` → `-- Solutions`; `\usepackage{amsmath}` removed (loaded by `mycommands.sty`).

**Q1**
- `\frac{d}{dt}` → `\frac{\ddns}{\ddns t}` (5× across Q1(a), Q1(b) and Q2); Piola-Kirchhoff → Piola--Kirchhoff.
- "(i.e. the conservation of linear momentum)" → `i.e.\ the`.
- "(the sufficiency is clear, for necessity, act $\epsilon_{ipq}$ on the identity and use ...)" → "(the sufficiency is clear; for necessity, contract the identity with $\epsilon_{ipq}$ and use ...)" (comma splice; "act ... on" → "contract ... with").
- "right hand side" → "right-hand side" (3× in file).

**Q2**
- "Euler Lagrange" → "Euler--Lagrange".
- `\frac{\partial}{\partial x_{j}}(A_{ijkl}\frac{\partial u_{k}}{\partial x_{l}})` → `\left( \ldots \right)`; trailing comma added to that display (followed by "which is what we obtained ...").
- `E~\dd^{3}\mathbf{x}` → `E\,\dd^{3}\mathbf{x}`.

**Q3** (three mathematical/reference corrections)
- "we repeat the calculation from Lecture 12" → `Lecture~13` (the chain-rule calculation `\partial W/\partial F_{ij}` is in Lecture 13, not 12).
- Identity `\frac{\partial C_{kl}}{\partial F_{ij}}=\delta_{ik}F_{jl}+\delta_{il}F_{jk}` → `\frac{\partial C_{kl}}{\partial F_{ij}}=\delta_{jk}F_{il}+\delta_{jl}F_{ik}` (indices were swapped: since $C_{kl}=F_{mk}F_{ml}$, $\partial C_{kl}/\partial F_{ij}=\delta_{jk}F_{il}+\delta_{jl}F_{ik}$; this is the form that gives the quoted result $2F_{im}\,\partial U/\partial C_{mj}$ and matches Lecture 13).
- Last line of the second-derivative `align*`: `4F_{im}F_{kp}\frac{\partial^{2}U}{\partial C_{mj}\partial C_{pq}}` → `\ldots\partial C_{pl}}` (free index must be $l$, not the dummy $q$; then at $\mathbf{F}=\mathbf{I}$ it gives $4\,\partial^{2}U/\partial C_{ij}\partial C_{kl}$ as stated).
- "symmetric second order tensors" → "second-order tensors"; "and hence the elastic tensor can be identified with a $6\times6$ matrix, and hence has at most" → "..., and so has at most" (repeated "hence").
- Pre-stressed elastic tensor display `+\delta_{ik}\sigma_{lj}` → `+\delta_{ik}\sigma_{lj},` (followed by "where").

**Q4**
- "derived in Lecture 14" → `Lecture~15` (perturbation theory for plane waves is in Lecture 15).
- `\boldsymbol{\Gamma}` → `\bm{\Gamma}` (4×). The Christoffel matrix and the intermediate `-\frac{\zeta}{\rho}[\cos^{2}\theta+2\cos^{2}\theta+\cos^{2}\theta]` were checked against the given $A_{ijkl}$ and are correct (sum $=4\cos^{2}\theta$); unchanged.

**Q5**
- `i\hbar` → `\ii\hbar`; `-i(E t-\ldots)` → `-\ii(\ldots)` (2×); `e^{-i(Et-\varphi)/\hbar}` → `\ee^{-\ii(Et-\varphi)/\hbar}` (4×) with thin spaces `a\,\ee^{\ldots}`; `\frac{i}{\hbar}\nabla\varphi a` → `\frac{\ii}{\hbar}\nabla\varphi\,a`.
- `||\mathbf{p}||`, `||\nabla\varphi||` → `\|\ldots\|`; `\frac{d\mathbf{x}}{d\sigma}`, `\frac{d\mathbf{p}}{d\sigma}`, `\frac{d}{d\sigma}` → `\ddns`; `d\sigma'` → `\dd\sigma'`.

**Q6**
- "from first-principles" → "from first principles"; Euler-Lagrange → Euler--Lagrange (3× in Q6/Q7).
- "and where we write primes for derivatives with respect to $\gamma$ as we have already used dots for those with respect to $\sigma$" → "and where primes denote derivatives with respect to $\gamma$ (dots will be used below for derivatives with respect to $\sigma$)" (dots had not yet been introduced at that point).
- "the Hamiltonian form of the ray equations are obtained by starting from the p-wave Hamiltonian" → "... is obtained ... P-wave Hamiltonian"; "the p-wave speed" → "the P-wave speed".
- `ds=||\mathbf{x}'||\dd\gamma=\alpha~\dd\sigma` → `\dd s=\|\mathbf{x}'\|\dd\gamma=\alpha\,\dd\sigma`; all `||` → `\|` (22 pairs in file); `\frac{dx_{i}}{d\gamma}` → `\ddns`; `||\mathbf{p}||^2` → `\|\mathbf{p}\|^{2}`.
- Trailing commas added to the Hamiltonian display and to the `\dot{\mathbf{x}}=\alpha^{2}\mathbf{p}` display (both followed by "where").

**Q7**
- `ds`, `dz` → `\dd s`, `\dd z`; `\frac{d\mathbf{x}}{dz}`, `\frac{d}{dz}`, `\frac{dx}{dz}`, `d\mathbf{x}/dz` → `\ddns`; `2q\alpha~dz` → `2q\alpha\,\dd z`.
- "for $T,$ where" → "for $T$, where".
- "the upwards going half of the ray" → "the upward-going half of the ray".
- "Considering the initial direction of the ray lies in the $x-z$ plane, it follows that" → "Since the initial direction of the ray lies in the $(x,z)$-plane, it follows that".

**Equation references**
- None present; all displays are `equation*`/`align*` and remain unnumbered and unlabelled, as instructed.

## B. Left unchanged / questions for the author

- Q5: "The eikonal equation restricts $(\mathbf{x},\mathbf{p})$ to lie on a three-dimensional level surface on which the Hamiltonian is equal to the energy $E$." The level set $H=E$ in six-dimensional phase space is five-dimensional; the three-dimensional object is the Lagrangian submanifold $\mathbf{p}=\nabla\varphi(\mathbf{x})$ that lies inside it and is swept out by the characteristics. Left as written, but "level surface" may be worth rewording.
- Q3: "the equilibrium has been assumed to be stress free" left unhyphenated (predicative use); the problem sheet has attributive "stress-free".
- Q7: "Because $\alpha$ increases monotonically, this equation has a unique solution for each initial ray orientation" implicitly assumes $\alpha$ is unbounded (or at least exceeds $1/q$) at depth; left as is.
- Q4: `\cos^{2}\theta\hat{\nu}_{i}\hat{\nu}_{k}` etc. could take a thin space before `\hat{\nu}`; not required by STYLE.md, left.
- Line beginning " Considering the quasi P-wave" has a leading space in the source; harmless, left untouched.
- "co-ordinates" (Q7) left as is; STYLE.md gives no ruling.

## C. Style pass

Applied: amsmath removed; Unicode en dash in title → `--`; `d/dt`, `d/d\sigma`, `d/dz`, `dx_i/d\gamma` → `\ddns`; `ds`, `dz`, `d\sigma'` → `\dd`; `||` → `\|` (22×); `i` → `\ii`, `e^{}` → `\ee^{}`; `\boldsymbol` → `\bm` (4×); `~` → `\,` thin spaces (3×); `\left(\right)` for the tall bracket in Q2; hyphen → en-dash in Piola--Kirchhoff, Euler--Lagrange (4×); `Lecture~13`, `Lecture~15`; `i.e.\ `; "right-hand side" (3×); P-wave capitalised; `$(x,z)$-plane`; trailing commas on five displays that run into "where"/"which". Compiles with no errors or warnings (7 pages, as the original); latexdiff runs cleanly.
