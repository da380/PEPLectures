# Changes to lecture11.tex ("Lecture 22: Free oscillations")

## A. Corrections and rewordings

**Preamble / title**
- `\title{Lecture 22: Free Oscillations}` → `\title{Lecture 22: Free oscillations}`; `\section*{Outline and Motivation}` → `Outline and motivation` (sentence case).

**Outline and motivation**
- "It is due to this fact that the observations of the  eigenfrequencies can be made from sufficiently long time-series" → "that observations of the eigenfrequencies can be made from sufficiently long time series"; stray double spaces removed.

**Linearised equations of motion**
- "Within Lecture 20" → "Within Lecture~20".
- After eq. (1): "$\Omega_{i}$ is the angular velocity, $A_{ijkl}$ the elastic tensor, $\gamma_{i}^{1}$ the perturbation to the gravitational acceleration which is given" → "$\Omega_{i}$ the angular velocity, $A_{ijkl}$ the elastic tensor, and $\gamma_{i}^{1}$ the perturbation to the gravitational acceleration, which is given" (parallel list; non-restrictive "which").
- Eqs. (4) and (5): `(A_{ijkl}\frac{\partial u_{k}}{\partial x_{l}}-\overline{S}_{ij})` → `\left( ... \right)` (tall brackets).
- Eq. (4) ended "." before "where" → ",".
- "satisfies eq. (1) and (3)" / "must satisfy eq.(1) and (3)" → `eqs.~(\ref{eq:1}) and (\ref{eq:3})` (twice); "eq.(5)" → `eq.~(\ref{eq:5})`; "Eq.(5)" → `Eq.~(\ref{eq:5})`.
- Footnote on complex conjugate moved after the full stop.
- Footnote on "sesquilinear forms": 'The term "sesqui" is used instead of "bilinear" to indicate that the form is linear in one argument and conjugate linear in the other.' → "The prefix ``sesqui'' (one and a half) is used instead of ``bi'' to indicate that the form is linear in one argument and conjugate-linear in the other." (straight quotes removed; wording corrected); footnote marker moved after the full stop.
- Eq. (7) ended "." before "where time-derivatives" → ","; "time-derivatives" → "time derivatives".
- Eq. (8) (Coriolis form) had no terminal punctuation before "which we call" → comma added.
- "eq.(12)" → `eq.~(\ref{eq:12})`.

**Properties of the sesquilinear forms**
- Eq. (16): the two integrals were written `\int` without a domain while every other integral in the file is `\int_{M}` → `\int_{M}` added to both (consistency; no change of meaning). The equation ended "." before "where we have included arguments" → ",".
- "Proof of this identity not difficult but nor is it enlightening." → "The proof of this identity is not difficult, but nor is it enlightening."

**Eigenfunctions and eigenvalues in a non-rotating earth model**
- "Lecture 20" → "Lecture~20"; "low frequency seismology" → "low-frequency seismology".
- Eq. (18): `e^{\ii\omega t}` → `\ee^{\ii\omega t}`; terminal comma added before "with $\mathbf{s}$ the spatial part".
- "eq. (17)" → `eq.~(\ref{eq:17})`; "eq. (21)" → `eq.~(\ref{eq:21})`.
- Eq. (23): "=0 ." → "=0.".
- **Maths:** "if we suppose that $\omega_{1}\ne\omega_{2}$, then we must have [orthogonality]" → "$\omega_{1}^{2}\ne\omega_{2}^{2}$". The derivation (eq. (23)) only gives orthogonality when the *squared* eigenfrequencies differ, and the following sentence says exactly that ("eigenfunctions with different squared eigenfrequencies are orthogonal"); $\omega_{1}\ne\omega_{2}$ is not sufficient if $\omega_{1}=-\omega_{2}$.
- Eq. (24) had no terminal punctuation before "this meaning that" → comma added.
- "``mass weighted'' inner product" → "``mass-weighted''".
- Footnote after "stated without proof" moved after the colon; "in-viscid" → "inviscid".
- "in a finite-interval" → "in a finite interval"; double spaces after full stops in the enumerate removed.
- "we let $k=0,1,\dots$, label the different squared eigenfrequencies" → "we let $k=0,1,\dots$ label ..." (spurious comma).
- Eq. (26): "|km\rangle ," → "|km\rangle,".

**Eigenfunction expansions**
- Fig. 1 caption: "A seven day long vertical component seismogram following the a magnitude 8.2 earthquake" → "A seven-day-long vertical-component seismogram following the magnitude 8.2 earthquake"; "Note also that the overall amplitude" → "Note that the overall amplitude" (nothing precedes it for "also" to refer to). Label `fig:fig1` → `fig:1`.
- Fig. 2 caption: "corresponding to Fig. 1" → `Fig.~\ref{fig:1}`; "at higher-frequencies" → "at higher frequencies". Label `fig:fig2` → `fig:2`.
- "for any test-function" → "test function"; "studied in the frequency-domain, where" → "in the frequency domain, where".
- Eq. (28) (Fourier convention): `e^{-\ii\omega t}`, `e^{\ii\omega t}` → `\ee^{...}`; `\dd \omega` → `\dd\omega`; ended "." before "noting that" → ","; "the initial conditions that ... vanish before time $t=0$ has been built into" → "have been built into".
- Eq. (30): "|k'm'\rangle ," → "|k'm'\rangle,".
- "eq. (29)" → `eq.~(\ref{eq:29})`; eq. (31) ended "." before "which can be simplified" → ",".
- Eq. (32): "\rangle ." → "\rangle.".
- "eq.(25)" → `eq.~(\ref{eq:25})`; eq. (33) (coefficient formula) had no terminal punctuation → full stop added.
- Eq. (34): "|km\rangle ," → "|km\rangle,"; "into the time-domain can be readily done" → "into the time domain".
- "eq. (34)" → `eq.~(\ref{eq:34})` (three times); "for each k" → "for each $k$"; "squared-eigenfrequencies are all positive" → "squared eigenfrequencies".
- "the observed time-series" / "such a time-series" / "length of the time-series" / "longer time-series provide" → "time series" (noun); "in the frequency-domain." → "in the frequency domain."
- Eq. (35): `\mathbf{u}_{obs}` → `\mathbf{u}_{\mathrm{obs}}`; terminal comma added before "where we neglect spatial arguments".
- **Maths:** Eq. (36): `\tilde{\mathbf{u}}_{obs}(\omega)\propto\int\tilde{\mathbf{u}}(\omega')\tilde{h}(\omega'-\omega)d\omega'` → `\tilde{\mathbf{u}}_{\mathrm{obs}}(\omega)\propto\int\tilde{\mathbf{u}}(\omega')\tilde{h}(\omega-\omega')\dd\omega'`. The convolution theorem gives the kernel $\tilde{h}(\omega-\omega')$, not $\tilde{h}(\omega'-\omega)$ (the latter is a cross-correlation; for real $h$ it equals $\tilde{h}(\omega-\omega')^{*}$, so the qualitative conclusion is unaffected). Also `d\omega'` → `\dd\omega'`.
- "with the right hand side" → "right-hand side".

**Gravitational stability**
- "If, instead, a squared eigenfrequency, $\omega_{k}^{2}$, was negative" → "were negative" (subjunctive).
- "a double-root of the denominator" → "a double root"; "can be shown to results in a term" → "result in".
- Eq. (37) ended "." before "we can set the test function" → ","; eq. (38) ended "." before "where we have made use" → ",".
- **Notation:** pointwise stability constant renamed $k$ → $c_{0}$ in eq. (41) and in "for some constant $c_{0}>0$", since $k$ is the mode index throughout the lecture. Eq. (41) also given a terminal comma, and a thin space inserted: `c_{0}\,e_{ij}^{*}e_{ij}`.
- Inline definition of the linearised strain: `\frac{1}{2}(\frac{\partial u_{i}}{\partial x_{j}}+...)` → `\left( ... \right)`.
- **Notation:** rigid-body displacement `u_{i}=a_{i}+A_{ij}x_{j}` with "$A_{ij}$ an anti-symmetric matrix" → `B_{ij}` (both places), because $A_{ijkl}$ is the elastic tensor in this file (STYLE.md §3).
- **PHYSICS CORRECTION (please check):** "It can be shown that the pointwise stability condition holds in an isotropic medium if and only if P- and S-waves have positive phase speeds, with the P-wave speed being faster." → "It can be shown that the pointwise stability condition holds in an isotropic medium if and only if the bulk and shear moduli are both positive, and this in turn implies that P- and S-waves have positive phase speeds with $\alpha>\beta$." Reason: positive wave speeds ($\lambda+2\mu>0$, $\mu>0$) is the weaker strong-ellipticity condition and does not imply pointwise stability ($\kappa>0$, $\mu>0$); the implication runs only one way, as now stated.
- "squared-eigenfrequencies are all non-negative" → "squared eigenfrequencies"; "eigenspace at zero-frequency" → "at zero frequency" (twice).
- "with $S_{ij}$ symmetric, we see that these trivial modes" → "with $\overline{S}_{ij}$ symmetric" (the stress glut is written $\overline{S}_{ij}$ everywhere else, including in eq. (44) immediately above).
- "With gravitation  added" → single space; footnote on the gravitational binding energy moved after the full stop.

**What you need to know and be able to do**
- `itemize` → `enumerate` (labels `\item[(i)]` etc. unchanged).
- "bra-ket like notation" → "bra-ket-like notation"; "How to derive reality of the squared-eigenfrequencies" → "the reality of the squared eigenfrequencies"; "in the frequency-domain. You do not need" → "in the frequency domain."

**Equation references / labels**
- All 44 displayed equations are `equation` environments (no `align`), so `\label{eq:N}` was attached to the Nth one. The hard-coded references eq.(1), (3), (5), (12), (17), (21), (25), (29), (34) were all checked against this count and all point at the intended equations; no renumbering corrections were needed. The compiled PDF prints the same numbers as the original.

## B. Left unchanged / questions for the author

- "then the field $\mathbf{u}$ must satisfy eqs. (1) and (3), this being the weak form of the problem." — "this" refers back to eq. (5) holding for all $\mathbf{w}$ rather than to the sentence's last clause. Left as is; a possible tidier version is "..., and so eq. (5) holding for all $\mathbf{w}$ is the weak form of the problem."
- Fig. 2 caption says "The amplitude spectra" (plural) for the spectrum of the single seismogram in Fig. 1. Left; change to "spectrum" if the figure shows only one panel.
- Eq. (23) is written $(\omega_{1}^{2}-\omega_{2}^{2*})\langle\mathbf{s}_{1}|P|\mathbf{s}_{2}\rangle^{*}=0$; direct substitution of the conjugate of eq. (22) into eq. (21) gives $(\omega_{2}^{2*}-\omega_{1}^{2})\langle\ldots\rangle^{*}=0$. Equivalent, so left alone.
- "rigid body motion(s)" left unhyphenated, as in the Lecture 12 revision.
- In the eigenvalue-problem discussion the eigenfrequency is called $\omega$ although only $\omega^{2}$ is determined; this is the standard usage and was left.
- Sign conventions in eqs. (1)–(2) (in particular $-\rho\gamma_{i}^{1}$ and the kernel of $\gamma_{i}^{1}$) were taken on trust from Lecture 20 and not re-derived.

## C. Style pass

Applied: sentence-case title/heading; `\label{eq:1}`–`\label{eq:44}` and `fig:1`, `fig:2`; all references via `\ref` with `~`; `Lecture~20`; `\ee^{\ii\omega t}`; `\dd\omega`, `\dd\omega'`; `\mathrm{obs}`; `\left(\right)` around tall brackets; footnote markers after punctuation; ``quotes''; terminal punctuation on every display; "squared eigenfrequency", "frequency domain"/"time domain" (noun) vs "frequency-domain"/"time-domain" (adjective), "time series" (noun), "right-hand side", "low-frequency", "double root", "inviscid"; symbol clashes $k\to c_{0}$, $A_{ij}\to B_{ij}$; summary list `enumerate`. No place where the style guide could not be applied. Overfull boxes: 6 in the revised file versus 8 in the original, none new.
