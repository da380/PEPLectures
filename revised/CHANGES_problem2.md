# Changes to problem2.tex ("Seismology Problem Set 2")

## A. Corrections and rewordings

**Preamble / title**
- `\title{Problem Set 2}` → `\title{Seismology Problem Set 2}` (matches the solutions file). `\usepackage{amsmath}` deleted (loaded by `mycommands.sty`). 10pt kept. All displays stay unnumbered (`equation*`/`gather*`), so no `\label`s were added.

**Q1**
- "Suppose that the data takes the form" → "the data take the form" (data plural, STYLE.md §4).
- `\mathbf{m}_{in}`, `\mathbf{m}_{out}` → `\mathbf{m}_{\mathrm{in}}`, `\mathbf{m}_{\mathrm{out}}` (4×).

**Q2**
- "and $\mathbf{e}$ and $n$-dimensional vector" → "and $\mathbf{e}$ an $n$-dimensional vector" (typo).
- "We further suppose  that the data has been processed" → "We further suppose that the data have been processed" (double space; data plural).
- (a) "geometric region of the model-space" → "of the model space" (consistent with "the model space" in (b)).
- (b) "must lie on the boundary-set defined by" → "on the boundary set defined by" (no hyphen needed).
- (b) `$ \|\mathbf{d}\|^{2} > ...$` → `$\|\mathbf{d}\|^{2} > ...$` (stray space).

**Q3**
- Symbol clash resolved: $\sigma$ is the applied surface load, so the standard deviations in the misfit are renamed $\sigma_{i}$ → $\varepsilon_{i}$ in the display for $J$ and in the text ("where $\varepsilon_{i}$ is the standard deviation for the $i$th observation"). The same change is made in `solution2.tex`.
- `\mathbf{u}_{i}^{obs}` → `\mathbf{u}_{i}^{\mathrm{obs}}` (2×); `$i=1,...,m$` → `$i=1,\dots,m$`; `||...||^{2}` → `\|...\|^{2}`; "the $i$-th observation" → "the $i$th observation"; "Lecture 19" → "Lecture~19".
- Comma added after the display for $J$ (the sentence continues with "where ...").
- In the `gather*` block the first line ended `=0;` → `=0,` (comma between parallel relations).

**Q4**
- `\boldsymbol{\varphi}` → `\bphi` in the argument list of $\mathcal{L}$; `||\mathbf{v}||^{2}` → `\|\mathbf{v}\|^{2}`.
- Comma added after the display $\delta\mathcal{L}/\delta\varphi_{i}=\rho\gamma_{i}$ ("where $\gamma_{i}$ ...").

**Q5**
- `\cdot\cdot\cdot` → `\cdots` (3×); `s \rho_{1}`, `s \phi_{1}` → `s\rho_{1}`, `s\phi_{1}` (stray spaces).
- "Similarly writing the potential in the form $\phi=\ldots$ show that" → "... in the form $\phi=\ldots$, show that" (comma after the participial clause).

**Q6**
- "Lecture 22" → "Lecture~22". "(which you don't need to know for this course)" kept as instructed.

**Q7**
- `\mathbf{P}`, `\mathbf{W}`, `\mathbf{H}` inside the bra-kets → plain italic `P`, `W`, `H` (the convention of Lecture~22 and STYLE.md §3) — 6 occurrences.
- "$\omega$ the eigenvalue" → "$\omega$ the eigenfrequency" ($\omega$ is the frequency; the eigenvalue of the quadratic problem is $\omega$ only in the generalised sense, and the next sentence already says "the eigenfrequency can be written").
- Comma added after the first display (the sentence continues with "where $\mathbf{s}$ is ...").

## B. Left unchanged / questions for the author

- Q3: the "equations of motion" are written as the static equilibrium equation $-\partial T_{ij}/\partial x_{j}=0$ (no inertia, no body force). This is consistent with a static loading problem and with the solution, so left as is; you might prefer "equilibrium equations".
- Q3: the boundary condition $T_{ij}\hat{n}_{j}=\sigma g_{i}$ has the load acting in the direction of $\mathbf{g}$ (i.e. pushing into the body); sign convention left untouched.
- Q5: "spherical polar co-ordinates" left with the hyphen (both spellings are British; STYLE.md does not list it).
- Q5 uses $\varphi$ for the azimuthal angle while Q4 uses $\varphi_{i}$/$\bphi$ for the motion. The two questions are independent so no clash arises; noted only because the solutions file had to be brought into line with $\varphi$ (see `CHANGES_solution2.md`).
- Q6: "the ten source parameters $(M_{ij},\mathbf{x}_{s},t_{s})$" — six independent moment-tensor components plus three plus one; correct as stated.
- Q7: the problem asks the student to "show that the eigenfrequency can be written ..." after having called $\omega$ "the eigenvalue" — I have made the wording consistent ("eigenfrequency") rather than the other way round; revert if you prefer "eigenvalue" throughout.

## C. Style pass

Applied: `amsmath` removed; `_{in}`/`_{out}`/`^{obs}` → `\mathrm`; `||` → `\|` (3×); `...` → `\dots`; `\cdot\cdot\cdot` → `\cdots` (3×); `\boldsymbol{\varphi}` → `\bphi`; `\mathbf{P/W/H}` → `P/W/H`; `Lecture N` → `Lecture~N` (2×); `$i$-th` → `$i$th`; trailing spaces at line ends and a double space after `\item` removed; terminal punctuation added to three displays. Line breaks and paragraphing otherwise preserved. No equation labels/references exist in this file (all displays unnumbered), so the labelling rule was not applicable. Compiles with no errors or warnings (3 pages, as the original); latexdiff runs cleanly. `Lamé` kept as the original UTF-8 character (compiles as before).
