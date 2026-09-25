# Changes to lecture9.tex ("Lecture 20: Self-gravitation and rotation")

## A. Corrections and rewordings

**Equation numbering and references**
- The original has 56 numbered displays (52 `equation` environments plus the 4-line `align` in the scaling section), not 54 as noted in the brief. Every display now carries `\label{eq:N}` (`eq:1`–`eq:56`).
- All eight hard-coded references were checked by counting and were correct: `eq.(4)`, `eq.(6)` → `eq.~(\ref{eq:4})`, `eq.~(\ref{eq:6})` (closed-form potential; definition of $\zeta$); `eq.(5)` → `eq:5` (closed-form acceleration); `eq.(18)` → `eq:18` ($\mathcal{V}_{g}=\frac{1}{2}\int\varrho\phi$); `eq.(6)`, `eq.(8)` → `eq:6`, `eq:8` (definition of $\zeta$; conservation of mass); `eq.(29)` (3×) → `eq:29` (decomposition of the motion); `eq.(10)` (2×) → `eq:10` (referential acceleration integral); `eq.(45)` → `eq:45` (explicit $\gamma_{i}^{1}$). The revised PDF prints the same numbers as the original.

**Outline and motivation**
- "how that travel time and waveform observations" → "how travel time and waveform observations".
- "long period free-oscillations of the Earth" → "long-period free oscillations of the Earth" (hyphenation).

**A review of Newtonian gravitation**
- "satisfies Poisson equation" → "satisfies Poisson's equation"; "The above Poisson equation is subject to" → "Poisson's equation is subject to" ("Poisson's equation" used consistently throughout, also in "Using Poisson's equation for $\phi$" and "using Poisson's equation $\nabla^{2}\zeta_{0}=4\pi G\rho$").
- Eq. (1) ended with a full stop followed by "in $\mathbb{R}^{3}$, where" → comma. Eq. (2) had no terminal punctuation → full stop. Eq. (3) → comma before "where"; `, \quad` → `,\quad`. Eq. (9) → comma before "which". Eq. (10) ended with a comma before a section break → full stop.
- "in the direction of outward normal" → "in the direction of the outward normal".
- "we similar find" → "we similarly find".
- "where  from Lecture 12 we know that  conservation" → single spaces, `Lecture~12` (2×).

**Gravitational binding energy**
- "The total binding energy associated with a perturbation $\delta\varrho$ to body's density by is, therefore, given by" → "The change in the binding energy associated with a perturbation $\delta\varrho$ to the body's density is, therefore, given by" (the displayed quantity is $\delta\mathcal{V}_{g}$, an increment, not the total).
- Eq. (11) "." + "to first-order accuracy" → ","; eq. (13) → comma added; eq. (14) ". where" → ", where"; eq. (15) → comma before "and so"; eq. (18) "$\dd^{3}\mathbf{y}. ,$" → "$\dd^{3}\mathbf{y},$"; eq. (19) ". where" → ", where".
- "Using the Poisson equation for $\phi$ we can write the above expression in the form" → "Using Poisson's equation for $\phi$, we can write the integrand of the above expression in the form" (the display that follows is an identity for the integrand, not for $\mathcal{V}_{g}$ itself).

**The equations of motion relative to inertial space**
- "the mathematical statement of Hamilton's principle is that vanishing of the functional derivative" → "is the vanishing of".
- "Here the the left hand side" → "Here the left-hand side"; "right hand side", "left hand side" → "right-hand side", "left-hand side" (as in the revised Lecture 12).
- Eq. (23) ". where" → ", where"; eqs. (26), (27) given terminal commas; `, \quad` → `,\quad` in eq. (26).
- "first Piola Kirchhoff stress tensor" → "first Piola--Kirchhoff"; "Euler-Lagrange" → "Euler--Lagrange" (2×).
- "action at a distance property of Newtonian gravitation" → "action-at-a-distance property" (2×, here and in the linearised section).

**Use of a co-rotating reference frame**
- MATHS (clarity): eq. (32) `R_{kj}\frac{\ddns R_{kl}}{\ddns t}=\epsilon_{jkl}\Omega_{k}` → `=\epsilon_{jml}\Omega_{m}`. The index $k$ was a dummy on both sides with different meanings; the statement was not wrong but was confusing. Eq. (33) already uses a free choice of dummies and is unchanged.
- MATHS (error): eq. (39) `\gamma_{i}=R_{ij}\overline{\gamma}_{i}` → `\gamma_{i}=R_{ij}\overline{\gamma}_{j}` (index $i$ appeared three times).
- "Here $\mathbf{\Phi}(t)$ denotes  translation of the reference frame" → "denotes the translation of the reference frame".
- Straight quotes `"fixed stars"` → ``fixed stars''.
- "carries  no net linear nor angular momentum\footnote{...}." → "carries no net linear or angular momentum.\footnote{...}" (grammar; footnote marker moved after the full stop).
- Eqs. (30), (31), (32), (34), (36), (37), (38), (39) given terminal punctuation (comma or full stop as the sentence requires); "$) .$" → "$).$" in eq. (33); `, \quad` → `,\quad` in eq. (35).
- "Piola-Kirchhoff" → "Piola--Kirchhoff".

**Linearised equations of motion**
- "Strictly, we should  use the term" → single space; "undergoing a steadily rotation" → "a steady rotation"; "and  hence the strain energy to be" → "and hence take the strain energy to be".
- Eq. (41) → comma before "with $s$"; eq. (42) ended "," before a new sentence → "."; "At zeroth-order the equations" → "At zeroth order, the equations".
- "a balance between the centrifugal, gravitational forces on the body and a resulting equilibrium stress" → "between the centrifugal and gravitational forces on the body and the resulting equilibrium stress".
- ERROR: "where $\gamma_{i}^{1}$ denotes the first-order perturbation to the referential gravitational potential" → "acceleration" ($\gamma_{i}$ is the referential gravitational acceleration, eq. (9)).
- Eq. (46) "$\hat{n}_{j} .$" → "$\hat{n}_{j}.$".
- "but  \textbf{integro partial differential equations}" → "but \textbf{integro-differential equations}".

**The neglect of self-gravitation and rotation**
- "$||\mathbf{u}||\sim U$, a length-scale $L$ and time-scale $T$" → "$\|\mathbf{u}\|\sim U$, a length scale $L$ and a time scale $T$"; "order of magnitude scaling estimates" → "order-of-magnitude scaling estimates"; "characteristics length-scale" → "characteristic length scale"; "small length-scales" → "small length scales"; "whose time-scale" → "whose time scale".
- Scaling `align`: line endings were ", ; ," → all commas. The final line keeps a comma (not a full stop) because the sentence continues with "where $\overline{\rho}$ and $\overline{\mu}$ denote ...".
- "elastic modulii" → "elastic moduli".
- MATHS (notation): `\frac{\partial\zeta_{0}}{\partial x_{i}\partial x_{j}}` → `\frac{\partial^{2}\zeta_{0}}{\partial x_{i}\partial x_{j}}` (second derivative).
- Eq. (53) `\rho\gamma_{i}^{1}\sim4\pi G\overline{\rho}^{2}U` given a full stop; eq. (56) `=\Omega T` given a full stop.
- MATHS (error): "$\overline{v}=\sqrt{\frac{\overline{\mu}}{\rho}}$" → "$\overline{v}=\sqrt{\overline{\mu}/\overline{\rho}}$" (the typical density is $\overline{\rho}$, as in the ratio it is used to simplify).
- "We should then ask whether gravitational forces will ever be of importance within the Earth?" → full stop (indirect question).
- Units: "6371 km" → `$6371\,\mathrm{km}$`; `6.67\times10^{-11} \, \text{m}^{3}\text{kg}^{-1}\text{s}^{-2}` → `6.67\times10^{-11}\,\mathrm{m^{3}\,kg^{-1}\,s^{-2}}`; `5000 \, \text{kg m}^{-3}` → `5000\,\mathrm{kg\,m^{-3}}`; `8000 \, \text{ms}^{-1}` → `8000\,\mathrm{m\,s^{-1}}`; comma added before "we find".
- "rotational effects play a relatively small (but not negligible) role the problem" → "role in the problem"; "e.g. those" → "e.g.\ those".

**What you need to know and be able to do**
- `itemize` → `enumerate` (labels `(i)`–`(iv)` unchanged); "terms with in an equation" → "within"; "you should insure you are happy with it" → "ensure".

## B. Left unchanged / questions for the author

- In "Use of a co-rotating reference frame", material frame indifference is written $W(\mathbf{x},\mathbf{F})=W(\mathbf{x},\overline{\mathbf{F}})$, dropping the $t$ argument that $W(\mathbf{x},t,\mathbf{F})$ carries elsewhere in the lecture. Harmless, left as is.
- Eq. (35) sets $\ddns\bm{\Phi}/\ddns t=\mathbf{0}$; the acceleration $\ddns^{2}\bm{\Phi}/\ddns t^{2}$ then vanishes as used in eq. (40), but the constant velocity $\ddns\bm{\Phi}/\ddns t$ is set to zero rather than merely constant. This is the standard choice and is unchanged.
- The inline `\frac{\delta\mathcal{L}}{\delta\varphi_{i}}` etc. in the prose after eq. (23) are left as inline fractions (author's choice).
- "anti-symmetric matrix" left hyphenated (STYLE.md does not list this word).
- "up till this point" left (acceptable British usage).
- The estimate "gravitational/elastic $\sim 3$" was checked: $4\pi G\overline{\rho}L^{2}/\overline{v}^{2}\approx 4\pi(6.67\times10^{-11})(5000)(6.371\times10^{6})^{2}/(8000)^{2}\approx 2.7$, consistent with "$\sim3$". Unchanged.
- One overfull `\hbox` (4 pt, 6 pt in the original) in the prose paragraph after eq. (14) ("where we have used the continuity of the potentials ...") is inherited from the original and left as is.

## C. Style pass

Applied: labels `eq:1`–`eq:56` and `eq.~(\ref{eq:N})` throughout; `\dd\mathbf{x}`/`\dd\mathbf{x}'` → `\dd^{3}\mathbf{x}`/`\dd^{3}\mathbf{x}'` (9×); `||...||` → `\|...\|` (3×); `\mathbf{\Phi}`, `\mathbf{\Omega}` → `\bm{\Phi}`, `\bm{\Omega}` (5×); `\overline{\boldsymbol{\gamma}}` → `\overline{\bm{\gamma}}`; `\partial^{2}` for the second derivative; straight quotes → ``''; Piola--Kirchhoff (2×), Euler--Lagrange (2×); `Lecture~12` (2×); `,\quad` between parallel relations (3×); `\mathrm` units with thin spaces; `e.g.\ `; footnote marker after punctuation; double spaces removed (6×); terminal punctuation on every display; summary list `itemize` → `enumerate`. `\text{gravitational forces}` etc. inside the ratio displays are prose words and were kept as `\text`. Compiles with no errors, no undefined or multiply-defined labels (9 pages, as the original); `latexdiff` runs cleanly.
