# Changes to lecture7.tex ("Lecture 18: Delay time tomography")

## A. Corrections and rewordings

**Equation numbering and references**
- The four-line `align` for $\partial J/\partial m_{p}$ numbered every line in the original (4 numbers). `\nonumber` added to the first three lines so the block carries one number; the file now has 40 numbered displays (original 43). All displays labelled `eq:1`--`eq:40` by position.
- "eq.(4)" → `eq.~(\ref{eq:4})` ($\delta T=-\int\delta v/v^{2}\,\dd s$) and "eq.(22)" → `eq.~(\ref{eq:22})` (modified misfit with $\lambda$): both verified correct by counting.
- "What we have obtained in eq.(31)" → `eq.~(\ref{eq:31})`, the regularised least squares solution $\mathbf{m}=(\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{A}+\lambda\mathbf{B})^{-1}\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{d}$. In the ORIGINAL numbering this formula was eq. (34) and "(31)" pointed at the bare matrix display $\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{A}+\lambda\mathbf{B}$; the sentence is clearly about the solution formula. (With the `align` counted as one equation, as the author evidently intended, the solution is indeed the 31st display, so the printed number is unchanged.)
- "Fig. 1" → `Fig.~\ref{fig:1}`. Figure labels `fig:1`, `fig:2` already present and in order.

**Outline and motivation**
- Section title "Outline and Motivation" → "Outline and motivation".

**Laterally varying structure from seismology**
- "The result of such studies are tomographic models" → "The results of such studies are" (agreement).

**Delay time measurements**
- Display $\delta T = T - T_{0}$: terminal "." → "," (followed by "where").

**Relating delay times to lateral variations**
- "if the Earth's velocity structure was equal to" → "were equal to" (counterfactual subjunctive).
- "which is a scalar-function in the earth model" → "scalar function".
- "in Lecture 12" → `Lecture~12`.
- "which you have established in Example sheet 3" → "which you established in the first problem set" (Problem 6 of Problem Set 1).
- "arc-length" → "arc length" (noun).

**Setting up the inverse problem**
- "Our data in such a study comprises delay times" → "comprise" (data plural).
- "Near identical methods" → "Nearly identical methods".
- "assume that (i) velocity variations ... and (ii) that the length scale" → "and (ii) the length scale" (redundant second "that").
- Display $\delta T_{i} = \dots + e_{i}.$ → "," (followed by "where").
- "(piece-wise) continuous" → "(piecewise)"; "under-determined problem" → "underdetermined problem".
- "There are a range of methods" → "There is a range of methods".
- "the $i$-th delay time" → "the $i$th delay time".

**Non-uniqueness and the null space**
- "Having reduced delay time tomography into a set of linear algebraic equations" → "to a set"; comma added after the following display $\mathbf{d}=\mathbf{A}\mathbf{m}+\mathbf{e}$ (sentence continues "we can start...").
- Comma added after display $\mathbf{A}\mathbf{m}_{0}=\mathbf{0}$ (followed by "then").
- "If the null space is not empty" → "If the null space is non-trivial (i.e.\ contains non-zero vectors)"; "non-empty null space" → "non-trivial null space" (twice). A null space always contains $\mathbf{0}$, so "non-empty" was mathematically inaccurate.

**Regularised least squares solutions**
- "a simple least-squares method" → "least squares method" (consistent with "least squares misfit/solution" elsewhere in the file).
- Comma added after the weighted least squares misfit display (followed by "where the factor of one half...").
- `\text{diag}(\sigma_{1}^{2},...,\sigma_{n}^{2})` → `\mathrm{diag}(\sigma_{1}^{2},\dots,\sigma_{n}^{2})`.
- "right hand side measures how well the data is fit" → "right-hand side ... the data are fit"; "left hand side" → "left-hand side".
- "a means of trading-off between" → "trading off between".
- "Reverting back to matrix vector notation" → "Reverting to matrix-vector notation" (redundant "back"; hyphen as in "matrix-vector notation" earlier in the file).
- Comma added after display $(\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{A}+\lambda\mathbf{B})\mathbf{m}=\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{d}$ (followed by "which will have").
- Eigenvalue display "$=\mu\mathbf{m}.$ with" → ",".
- "the unique solution of $\nabla J(\mathbf{m})=0$" → "$=\mathbf{0}$" (bold zero vector, consistent with the earlier display $\nabla J(\mathbf{m})=\mathbf{0}$).
- "must be a minimum of J" → "$J$"; `$|\!|\mathbf{m}|\!|$` → `$\|\mathbf{m}\|$`.
- "There  are a range of different schemes" → "There is a range of different schemes".

**A Bayesian perspective**
- "given that the data, $\mathbf{d}$, has been observed" → "have been observed".
- "Markov Chain Monte Carlo" → "Markov chain Monte Carlo".
- Comma added after the prior-PDF display (followed by "where we recall") and after the factored posterior display (followed by "where we have defined").
- "determined in closed-form" → "closed form"; "re-frames" → "reframes".
- "justifies use of the regularised least squares" → "justifies the use of regularised least squares".

**What you need to know and be able to do**
- `itemize` → `enumerate` (labels `\item[(i)]` etc. unchanged).
- Item (ii): "reduced into the standard form $\mathbf{A}\mathbf{m}=\mathbf{d}+\mathbf{e}$" → "reduced to the standard form $\mathbf{d}=\mathbf{A}\mathbf{m}+\mathbf{e}$" (the form actually used throughout the lecture; the original had the error on the wrong side).
- "from first-principles" → "from first principles".

## B. Left unchanged / questions for the author

- "source-receiver pair(s)" left hyphenated (compound adjective, not a name pair); STYLE.md does not list it.
- Bayes' theorem display uses `\dd \mathbf{m}` for the measure on model space. STYLE.md lists `\dd\mathbf{x}` under "Not" (in favour of `\dd^{3}\mathbf{x}`), but for an $m$-dimensional model space there is no natural analogue, so left as is.
- Displays $\nabla J(\mathbf{m})=\mathbf{0}$ ("... the equation [display] has a unique solution") and $\mathbf{A}^{T}\mathbf{C}^{-1}\mathbf{A}+\lambda\mathbf{B}$ ("... symmetric matrix [display] is invertible") deliberately carry no terminal punctuation, since each is the grammatical subject of the clause that follows.
- "end-points", "ray theoretic approximations", "matrix-vector notation" left as written.
- The statement that the second term $\lambda\mathbf{m}^{T}\mathbf{B}\mathbf{m}$ "is positive by definition" assumes the basis functions $\varphi_{j}$ are linearly independent (so that the Gram matrix $\mathbf{B}$ is positive definite). This is implicit and reasonable; not changed.
- "Many of the basic ideas, however, carry over to more complicated tomographic techniques including waveform tomography that we will discuss next time." reads slightly awkwardly (a comma before "including" and "which" for "that" would help), but was left as the author's voice.
- The single 1.28pt overfull `\hbox` (the paragraph containing "$n$-by-$m$ matrix") is present in the original as well.

## C. Style pass

Applied: `\nonumber` on multi-line `align`, `\label{eq:N}` on every display, `eq.~(\ref{...})`/`Fig.~\ref{...}`, `Lecture~12`, `\dots` for `...`, `\mathrm{diag}`, `\mathrm{ray}` (was `\text{ray}`) in integral subscripts, `\|\cdot\|` norm, `$10^{5}$--$10^{8}$` range, `e.g.\ `/`i.e.\ `, hyphen removal in underdetermined/piecewise/closed form/first principles/reframes/trading off, "right-hand/left-hand side", "data are/comprise", `enumerate` in the summary, and removal of doubled spaces inside sentences. Line breaks and paragraphing are otherwise unchanged. No place where the style guide could not be applied cleanly.

## D. Second pass (25 September 2026)

* Renumbered: "Lecture 19: Delay time tomography" (body-only file).
* Fig. 2 regenerated with `figsrc/L7F2.py` and its caption rewritten to describe the inversion actually used.
* Four further figures from the toy problem added (null space, damping and trade-off curve, checkerboard test, posterior standard deviation), each with linking text; new short paragraph on synthetic resolution tests.
* Closing paragraph of the Bayesian section rewritten at David's request.
* Toy problem introduced explicitly (24x24 pixel basis); sparse form of A for block parameterisations; SVD paragraph and Fig. 4; discrepancy principle; spike test (Fig. 6); Bayesian figure with correlated prior (Fig. 7); function-space subsection with Fig. 8 (pygeoinf).
