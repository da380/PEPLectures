# Changes to lecture8.tex ("Lecture 19: Waveform tomography")

## A. Corrections and rewordings

**Preamble**
- Removed `\usepackage{amsmath} % Ensure amsmath is loaded ...` (already loaded by `mycommands.sty`).

**Numerical wavefield simulations**
- "the time for solving the necessary forward problem is of order hours" → "the necessary forward problems" (one simulation per earthquake, so plural).
- "depends  on" → single space.

**Waveform measurements**
- "at the surface location, $\mathbf{x}_{r}$ following an earthquake" → "at the surface location, $\mathbf{x}_{r}$, following an earthquake" (closing comma; matches the later section).
- Fig. 1 caption: "Data is shown in black" → "Data are shown in black"; "right hand side" → "right-hand side"; "Tape et al. (2010)" → "Tape et al.\ (2010)".
- "define the least squares misfit" → "least-squares misfit" (adjectival; consistent with the later use).
- "closer their true values" → "closer to their true values".
- "In fact, least-squares misfit between seismograms is rarely used" → "In fact, the least-squares misfit ...".
- "illustrate wave form inversion" → "illustrate waveform inversion".

**Functional derivatives and sensitivity kernels**
- "To simplify notations" → "To simplify notation".
- "functional derivative of J with respect to" → "$J$".
- "supposing that ... is the only model parameter we can concretely write" → comma added after "parameter".
- Eq. (3): added terminal comma (followed by "where ..."); `K_{\beta}\delta\beta` → `K_{\beta}\,\delta\beta`.

**Gradient-based optimisation**
- Section title "Gradient based optimisation" → "Gradient-based optimisation"; likewise "Gradient based optimisation depends" and item (iii) of the summary (the file already used "gradient-based" elsewhere).
- `"grid search"` → ``grid search''.
- "Lecture 18" → `Lecture~18`.
- "had n-dimensions", "m points", "n is of order", → "$n$ dimensions", "$m$ points", "$n$ is of order".
- "Markov Chain Monte Carlo" → "Markov chain Monte Carlo".
- "eq. (5)" → `eq.~(\ref{eq:5})`.
- "But  this", "there do exist  approximate" → single spaces.
- Footnote after "uncertainty quantification" moved after the full stop.
- Eq. (6): `\lambda~DJ` → `\lambda\,DJ` and terminal comma added (followed by "where"); eq. (7) likewise `\lambda\,DJ`.
- "step-length" → "step length" (noun).
- "Having now arrived at a new model, $\mathbf{m}_{1}$, with a lower misfit we can" → comma after "misfit".
- Footnote after "found" moved after the full stop; inside it "decreased sufficiently much from a statistical perspective" → "decreased sufficiently from a statistical perspective".
- "The above method is the simplest example of gradient-based optimisation known as \textbf{the method of steepest descent}" → "... optimisation, and is known as ...".
- "More sophisticated algorithms for choosing the descent direction tend to be used in practice which have superior convergence properties." → "More sophisticated algorithms for choosing the descent direction, which have superior convergence properties, tend to be used in practice." (dangling relative clause).

**How not to calculate the gradient of $J$**
- "depends on being able to both solve the forward problem, and to calculate" → "being able both to solve the forward problem and to calculate".
- "The problem is, as noted before, the number of model parameters, $m$, is typically large and so" → "... that the number ..., is typically large, and so".
- "needed in iterative solution" → "needed in the iterative solution".

**How you should calculate the gradient of $J$**
- "method of \textbf{Lagrange multipliers} which many of you" → comma before "which"; "For those interested a non-examinable outline" → comma after "interested"; double space before "We begin" removed; "a \textbf{Lagrangian} for the problem which takes" → comma before "which".
- "right hand side" → "right-hand side" (twice).
- "features a second Lagrangian multiplier field" → "Lagrange multiplier field".
- Eq. (10): removed the trailing comma ("the equality [eq] holds so long as" takes no comma).
- `eq. \eqref{eq:11}`, `\eqref{eq:12}`, `\eqref{eq:13}` → `eq.~(\ref{eq:N})`.
- "eq. (1)", "eq. (20)", "eq.(11)", "eq. (10)" → `eq.~(\ref{eq:1})`, `(\ref{eq:20})`, `(\ref{eq:11})`, `(\ref{eq:10})`.
- "Starting with the time-derivatives" → "time derivatives".
- "Making use of the hyperelastic symmetry, $A_{ijkl}=A_{klij}$ we can" → closing comma after the symmetry relation.
- "for an arbitrary function its boundary value and that of its normal derivative can be chosen separately" → "for an arbitrary function, its boundary value ... separately".
- "the second term implies that $w'$ is the restriction of $u'$ to $\partial M$" → bold vectors `$\mathbf{w}'$`, `$\mathbf{u}'$` (consistent with the rest of the section).
- "with this problem driven by a traction applied at the observation point and being equal to the difference" → "... at the observation point that is equal to the difference".
- "It has also shown that the second Lagrange multiplier field" → "We have also shown that ..." (the subject "it" had no clear referent).
- "Integrating by parts once, this expression can be simplified" → "... can be simplified to".
- "From this solution determine J" → "$J$".
- Footnote after "\textbf{adjoint traction}" moved after the full stop.
- "this being done so the forward field" → "so that the forward field".
- "$10^{4}-10^{6}$ this is an enormous improvement" → "$10^{4}$--$10^{6}$, this is ...".
- "It is worth  commenting" → single space.
- Fig. 2 caption: "Bozdag et al. (2016)" → "Bozda\u{g} et al.\ (2016)"; "starting model S362ANI which was obtained using ray theoretic methods" → "S362ANI, which was obtained using ray-theoretic methods"; "While the difference in the models might not seem huge, they are still significant" → "While the differences between the models might not seem huge, they are still significant" (number agreement); "a range of hot spots including that under Yellowstone" → "hot spots, including that"; "as methods continue to improve the difference between ray theoretic and waveform" → "as methods continue to improve, the difference between ray-theoretic and waveform".

**Banana doughnut kernels**
- Section title "Banana Doughnut kernels" → "Banana doughnut kernels".
- Receiver location `\mathbf{x}_{s}` → `\mathbf{x}_{r}` throughout this section (prose and eqs. (35), (36)) for consistency with the rest of the lecture, where the receiver is `\mathbf{x}_{r}`.
- Eq. (34): terminal full stop → comma (followed by "with $\hat{\bm{\nu}}$ ...").
- `\hat{\boldsymbol{\nu}}` → `\hat{\bm{\nu}}`.
- "zeros everything out" → "zeroes everything out".
- "when $s^{obs}$ and s are maximally aligned" → "$s^{\mathrm{obs}}$ and $s$".
- Eq. (37): `\frac{dC}{d\tau}` → `\frac{\ddns C}{\ddns\tau}`; terminal comma added (followed by "and this acts to define ...").
- "functional derivatives of $\overline{\tau}$ and here we need only" → comma before "and here".
- Eq. (38): terminal comma added (followed by "for the first-order perturbation ..."); thin spaces `K_{\alpha}\,\delta\alpha+K_{\beta}\,\delta\beta`.
- "the ray theoretic result" → "the ray-theoretic result"; eq. (39): `\int_{\text{ray}}` → `\int_{\mathrm{ray}}`.
- "with c the appropriate phase speed" → "$c$".
- Fig. 3 caption: "Cross sections" → "Cross-sections".
- "i.e. the region" → "i.e.\ the region".
- "The form of these  sensitivity kernels deserves  some comment", "the true kernels  vanish", "named  them" → single spaces.
- "When this was first observed from numerical calculations this result seemed" → comma after "calculations".
- "the kernels in Fig.3" → `Fig.~\ref{fig:3}`.
- "(1942-2007)" → "(1942--2007)".

**What you need to know and be able to do**
- `itemize` → `enumerate` (labels `(i)`–`(iv)` unchanged).
- Item (iii): "Know how gradient based optimisation works" → "How gradient-based optimisation works" (parallel with the other items).
- Item (iv): "you may be asked perform parts" → "asked to perform parts".

**Appendix**
- `\section*{The Lagrange multiplier theorem - NON-EXAMINABLE}` → `\section*{The Lagrange multiplier theorem (non-examinable)}`.
- Eq. (40): terminal full stop → comma (followed by "with $\mathbf{m}$ the model vector").
- "Here $a$ is some potentially non-linear vector-valued function" → `$\mathbf{a}$`; "partial derivatives of $a$" → `$\mathbf{a}$`.
- "in which $\mathbf{m}$ act as parameters" → "in which $\mathbf{m}$ acts as a parameter" (agreement).
- "between $J(\mathbf{u})$ and $\hat{J}(\mathbf{m})$ but here" → comma before "but".
- `"brute-force"` → ``brute-force''.
- "expand to first-order to obtain", "Expanding this out to first-order" → "to first order" (adverbial, STYLE.md §4).
- `\cdot\cdot\cdot` → `\cdots` (six places).
- Terminal punctuation added to displays: eq. (44) comma (followed by "where"), eq. (45) comma ("where it is understood"), eq. (46) full stop, eq. (48) full stop, eq. (51) comma ("with $\mathbf{u}'$ ..."), eq. (54) comma ("for any"), eq. (55) comma ("and as"), eq. (58) comma ("for any").
- "an n-dimensional vector", "an m-dimensional one", "also n-dimensional vectors" → `$n$-dimensional`, `$m$-dimensional`, `$n$-dimensional`.
- "Assuming that the linear operator, $D_{\mathbf{u}}\mathbf{a}(\mathbf{u},\mathbf{m})$ is invertible" → matching comma after the operator.
- "eq. (42)", "eq. (41)" (×4), "eq.(44)", "eq. (53)", "eq. (52)", "eq. (56)" → `eq.~(\ref{eq:N})` with the same N.

**Equation references and labels**
- Every displayed equation now carries `\label{eq:N}` with N its position in the file (`eq:1`–`eq:59`; the pre-existing `eq:11`–`eq:13` coincide with their positional numbers). Multi-line `align` blocks meant as one equation keep `\nonumber` on all but the last line; the kernel `align` is labelled `eq:31`–`eq:33`.
- All hard-coded numbers in the original were checked against the original PDF and were correct; no mapping corrections were needed. The compiled revised PDF prints the same numbers as the original for every reference.

## B. Left unchanged / questions for the author

- Eq. (9) (the Lagrangian) and eq. (16) use the boundary multiplier $\mathbf{w}'$ only through the traction term; the statement that eq. (13) "gives the boundary conditions for $\mathbf{u}$" is correct as written (free-surface traction condition). No change made.
- The step-length symbol $\lambda$ in eqs. (6)–(7) is not the Lamé parameter; left as is, since no Lamé parameter appears in this lecture.
- "a factor of $\frac{3}{m+1}$" left as an inline `\frac`; could be written `$3/(m+1)$` if preferred.
- "for eq. (41) to admit unique solutions" (appendix): eq. (41) is `a[û(m),m]=0`; arguably the reference should be to eq. (40) `a(u,m)=0`, which is the equation being solved. Left pointing at eq. (41) as in the original since the two are the same equation with the solution substituted.
- "eq. (56) for the Lagrange multiplier ... involves the adjoint" (last sentence): eq. (56) is the linear equation, eq. (57) its solution; either fits the sentence. Left as eq. (56).
- "counter-intuitive" left hyphenated; "Southern California" left capitalised as in the original.
- The inner product in the appendix (`⟨δm | [D_m a]† u'⟩`) is written with the arguments in the opposite order to the main text (`⟨DJ | δm⟩`); harmless for a real inner product, left unchanged.

## C. Style pass

Applied: `||`→`\|`; `^{obs}`→`^{\mathrm{obs}}`; `dt`→`\dd t`, `ds`→`\dd s`; `\frac{dC}{d\tau}`→`\ddns`; `\text{ray}`→`\mathrm{ray}`; `\boldsymbol`→`\bm`; `\cdot\cdot\cdot`→`\cdots`; `~`→`\,` before `DJ`; straight quotes→``''; hyphen→en-dash in ranges and life dates; `Lecture~18`; `Fig.~\ref`, `eq.~(\ref{...})` (`\eqref` replaced by `\ref` for uniformity); footnote markers moved after punctuation; `e.g.`/`i.e.`/`et al.` control spaces; "ray theoretic"→"ray-theoretic" (adjective); summary list `enumerate`. Positional `\label{eq:N}` on all 59 displays; `fig:1`–`fig:3` already conformed. No place where the style guide could not be applied.
